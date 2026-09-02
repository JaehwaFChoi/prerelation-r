"""parity_driver.py -- cross-language parity harness, Python side and report.

Runs the Python reference implementation over every committed golden
fixture, reads the values the R implementation wrote for the same
fixtures, and prints a table of quantity / Python value / R value /
absolute difference / verdict.

The contract (tests/golden/README.md) is agreement to 1e-12 on the
component quantities and *exact* agreement on the permutation p-values.
Bit-identity is not the target: the independence baseline is a double sum
and the two languages accumulate it in a different order, so a difference
of one unit in the last place is expected and passes.

Both sides read the same fixture CSV and the same permutation index matrix
from disk, so neither is comparing against numbers it generated itself.

Usage: python parity_driver.py <golden_dir> <r_out_csv> [--tol 1e-12]
"""

import argparse
import csv
import json
import os
import sys

import numpy as np

from prerelation.core import direction, prereq_index
from prerelation.reference import pi_envelope
from prerelation.core import _baseline_mean  # noqa: F401  (definition check)

DELTA = 0.05
TOP_Q = 0.8

FIXTURES = ["product", "min", "independent", "equivalence",
            "partial_equivalence", "ecpe_slice"]

# Quantities compared at the tolerance, and quantities compared exactly.
EXACT = {"n", "n_interior", "perm_p", "n_tail_band"}


def read_fixture(path):
    """Read a fixture CSV. Explicit float conversion, as on the R side."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    if header != ["x", "y"]:
        raise ValueError(f"{path}: expected header x,y, got {header}")
    x = np.array([float(r[0]) for r in body])
    y = np.array([float(r[1]) for r in body])
    return x, y


def read_perm_indices(path):
    """Read the committed permutation matrix. Zero-based, as written."""
    P = np.loadtxt(path, delimiter=",", dtype=int)
    n = P.shape[1]
    for r in range(P.shape[0]):
        if sorted(P[r]) != list(range(n)):
            raise ValueError(f"{path} row {r} is not a permutation of 0..{n - 1}")
    return P


def components(x, y, P):
    n = x.size
    res = prereq_index(x, y)
    rev = prereq_index(y, x)
    dl = direction(x, y)

    v = float(np.mean(np.maximum(y - x, 0.0)))
    v0 = float(np.mean(np.maximum(y[None, :] - x[:, None], 0.0)))

    u = np.clip(y / np.maximum(x, 1e-9), 0.0, 1.0)
    ceil_mask = u >= 1.0 - DELTA
    x_top = x >= np.quantile(x, TOP_Q)
    p1_top = float(np.mean(ceil_mask[x_top])) if x_top.sum() > 0 else 1.0

    obs = res["PI"]
    cnt = sum(prereq_index(x, y[P[r]])["PI"] >= obs for r in range(P.shape[0]))
    perm_p = (cnt + 1) / (P.shape[0] + 1)
    env = pi_envelope(x, y)

    return {
        "n": float(n),
        "v": v,
        "v0": v0,
        "A1": res["A1"],
        "mass_ceiling_band": float(np.mean(ceil_mask)),
        "mass_interior": float(np.mean(~ceil_mask)),
        "n_interior": float(np.sum(~ceil_mask)),
        "p1_top": p1_top,
        "q": res["q"],
        "ell": res["ell"],
        "A2": res["A2"],
        "PI": res["PI"],
        "PI_reverse": rev["PI"],
        "Delta": dl[0],
        "perm_p": perm_p,
        # reference class and envelope (E-6a golden keys)
        "n_tail_band": float(env["n_tail"]),
        "D_star": env["D_star"],
        "sup_q": env["sup_q"],
        "inf_q": env["inf_q"],
        "PI_hi": env["PI_hi"],
    }


def read_r_output(path):
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out.setdefault(row["fixture"], {})[row["quantity"]] = float(row["value"])
    return out


def ulps(a, b):
    """Distance in units in the last place, for reporting only."""
    if a == b:
        return 0
    if np.isnan(a) or np.isnan(b):
        return float("nan")
    ia = np.abs(np.frombuffer(np.float64(a).tobytes(), dtype=np.int64)[0])
    ib = np.abs(np.frombuffer(np.float64(b).tobytes(), dtype=np.int64)[0])
    return int(abs(ia - ib))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("golden")
    ap.add_argument("r_out")
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--expected", default=None,
                    help="expected.json to cross-check the Python side against")
    ap.add_argument("--label", default="R", help="name of the other language")
    ap.add_argument("--table", default=None,
                    help="write the per-quantity table as CSV to this path")
    args = ap.parse_args()

    r_values = read_r_output(args.r_out)
    expected = None
    if args.expected:
        with open(args.expected) as fh:
            expected = json.load(fh)

    rows = []
    n_fail = 0
    n_bit_identical = 0
    worst = (0.0, None, None)
    drift = []

    for name in FIXTURES:
        x, y = read_fixture(os.path.join(args.golden, f"fixture_{name}.csv"))
        P = read_perm_indices(
            os.path.join(args.golden, f"perm_indices_n{x.size}.csv"))
        py = components(x, y, P)

        if name not in r_values:
            print(f"FAIL: R output has no rows for fixture {name}")
            n_fail += 1
            continue

        if expected is not None and name in expected:
            for key, val in expected[name].items():
                if key in py and float(val) != py[key]:
                    drift.append((name, key, float(val), py[key]))

        for key, pv in py.items():
            if key not in r_values[name]:
                rows.append((name, key, pv, float("nan"), float("nan"),
                             "FAIL (missing in R output)"))
                n_fail += 1
                continue
            rv = r_values[name][key]
            diff = abs(pv - rv)
            if key in EXACT:
                ok = (pv == rv)
                verdict = "PASS (exact)" if ok else "FAIL (exact required)"
            else:
                ok = diff <= args.tol
                verdict = "PASS" if ok else "FAIL"
            if not ok:
                n_fail += 1
            if pv == rv:
                n_bit_identical += 1
            if diff > worst[0]:
                worst = (diff, name, key)
            rows.append((name, key, pv, rv, diff, verdict))

    width = max(len(r[1]) for r in rows) if rows else 10
    print()
    print(f"{'fixture':<20}{'quantity':<{width + 2}}"
          f"{'python':>24}{args.label:>24}{'abs diff':>12}  verdict")
    print("-" * (20 + width + 2 + 24 + 24 + 12 + 10))
    for name, key, pv, rv, diff, verdict in rows:
        print(f"{name:<20}{key:<{width + 2}}{pv:>24.17g}{rv:>24.17g}"
              f"{diff:>12.3e}  {verdict}")

    print("-" * (20 + width + 2 + 24 + 24 + 12 + 10))
    if args.table:
        with open(args.table, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["fixture", "quantity", "python", args.label,
                        "abs_diff", "bit_identical", "tolerance", "verdict"])
            for name, key, pv, rv, diff, verdict in rows:
                w.writerow([name, key, repr(pv), repr(rv), f"{diff:.3e}",
                            int(pv == rv),
                            "exact" if key in EXACT else f"{args.tol:g}", verdict])
    print(f"quantities compared      : {len(rows)}")
    print(f"bit-identical            : {n_bit_identical}")
    print(f"failures                 : {n_fail}")
    if worst[1] is not None:
        d, fx, ky = worst
        pv = next(r[2] for r in rows if r[0] == fx and r[1] == ky)
        rv = next(r[3] for r in rows if r[0] == fx and r[1] == ky)
        print(f"worst absolute difference: {d:.3e}  ({fx}/{ky}, {ulps(pv, rv)} ulp)")
    print(f"tolerance                : {args.tol:g} "
          f"(exact for {', '.join(sorted(EXACT))})")

    if drift:
        print()
        print("DRIFT against expected.json (Python side no longer reproduces "
              "the committed reference):")
        for name, key, exp, got in drift:
            print(f"  {name}/{key}: expected {exp!r} got {got!r}")
        n_fail += len(drift)
    elif expected is not None:
        print("expected.json cross-check: the Python side reproduces every "
              "committed value exactly")

    print()
    print("RESULT:", "PASS" if n_fail == 0 else f"FAIL ({n_fail})")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
