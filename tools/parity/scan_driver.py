"""scan_driver.py -- parity harness for the pairwise scan.

Two modes.

  --setup   writes theta.csv (persons x attributes) and perm_scan.csv (a
            pool of permutations) with %.17g / integer text, then stops.
            Both languages read these files back, so neither is comparing
            against numbers it generated in memory.

  --compare runs the Python scan with the shared permutations, reads what
            the R implementation wrote for the same inputs, and prints the
            comparison.

The permutation pool holds n_perm + n_pairs rows; the ordered pair at
pair_position p uses rows p .. p + n_perm - 1. That is a fixed reference
set rather than an independent draw per pair -- adequate here, because the
object under test is agreement between two implementations and not the
sampling behaviour of the p-values, and it keeps one file on disk instead
of twelve.

Design floor: K = k(k-1) ordered pairs at level alpha need
n_perm >= K/alpha - 1. With k = 4 and alpha = 0.05 that is 239.
"""

import argparse
import csv
import json
import os
import sys

import numpy as np

from prerelation.core import prereq_index
from prerelation.scan import bh_fdr, find_cycles, transitive_reduction

N = 200
K_ATTR = 4
N_PERM = 239
ALPHA = 0.05
NAMES = ["a1", "a2", "a3", "a4"]


def ceiling(x):
    return np.clip(0.15 + 0.85 * x ** 0.8, 0.0, 1.0)


def make_theta():
    """A three-step chain a1 -> a2 -> a3 with a4 free of the others.

    The chain gives the scan direct edges, an indirect pair for the
    transitive reduction to remove, and an attribute that should stay
    isolated.
    """
    rng = np.random.default_rng(np.random.SeedSequence([20260828, 3]))
    a1 = rng.uniform(0.02, 0.98, N)
    a2 = np.clip(ceiling(a1) * rng.uniform(0.0, 1.0, N), 0.0, 1.0)
    a3 = np.clip(ceiling(a2) * rng.uniform(0.0, 1.0, N), 0.0, 1.0)
    a4 = rng.uniform(0.02, 0.98, N)
    return np.column_stack([a1, a2, a3, a4])


def setup(outdir):
    theta = make_theta()
    with open(os.path.join(outdir, "theta.csv"), "w") as fh:
        fh.write(",".join(NAMES) + "\n")
        for row in theta:
            fh.write(",".join("%.17g" % v for v in row) + "\n")

    n_pairs = K_ATTR * (K_ATTR - 1)
    rng = np.random.default_rng(20260828)
    pool = np.stack([rng.permutation(N) for _ in range(N_PERM + n_pairs)])
    with open(os.path.join(outdir, "perm_scan.csv"), "w") as fh:
        for row in pool:
            fh.write(",".join(str(int(i)) for i in row) + "\n")
    print(f"setup: theta {theta.shape}, permutation pool {pool.shape}, "
          f"n_perm={N_PERM}, pairs={n_pairs}")
    print(f"design floor K/alpha - 1 = {n_pairs}/{ALPHA} - 1 = "
          f"{n_pairs / ALPHA - 1:g}; n_perm = {N_PERM}")


def read_theta(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    names = rows[0]
    data = np.array([[float(v) for v in r] for r in rows[1:]])
    return names, data


def read_pool(path):
    return np.loadtxt(path, delimiter=",", dtype=int)


def python_scan(theta, names, pool):
    """The package scan, with the shared permutation pool substituted for
    the internal generator. Everything except the source of the
    permutations is prerelation.scan's own logic, re-expressed here because
    the Python scan has no indices hook."""
    n, k = theta.shape
    pi_matrix = np.full((k, k), np.nan)
    comps = {}
    for i in range(k):
        for j in range(k):
            if i != j:
                res = prereq_index(theta[:, i], theta[:, j])
                pi_matrix[i, j] = res["PI"]
                comps[(i, j)] = res

    records = []
    pair_position = 0
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            x, y = theta[:, i], theta[:, j]
            res = comps[(i, j)]
            obs = res["PI"]
            P = pool[pair_position:pair_position + N_PERM]
            cnt = sum(prereq_index(x, y[P[r]])["PI"] >= obs for r in range(N_PERM))
            records.append({
                "source": names[i], "target": names[j],
                "pi": res["PI"], "pi_reverse": pi_matrix[j, i],
                "delta": res["PI"] - pi_matrix[j, i],
                "A1": res["A1"], "A2": res["A2"], "q": res["q"], "ell": res["ell"],
                "p_value": (cnt + 1) / (N_PERM + 1),
                "n": n, "n_perm": N_PERM,
            })
            pair_position += 1

    p_adj = bh_fdr([r["p_value"] for r in records])
    edges = []
    for rec, pa in zip(records, p_adj):
        rec["p_adj"] = float(pa)
        keep = bool(pa <= ALPHA and rec["pi"] >= 0.0 and rec["delta"] > 0)
        rec["edge"] = keep
        if keep:
            edges.append([rec["source"], rec["target"]])
    cycles = find_cycles(names, [tuple(e) for e in edges])
    reduced = None if cycles else [list(e) for e in transitive_reduction(
        names, [tuple(e) for e in edges])]
    return records, edges, cycles, reduced


NUMERIC = ["pi", "pi_reverse", "delta", "A1", "A2", "q", "ell", "p_value", "p_adj"]


def compare(outdir, tol):
    names, theta = read_theta(os.path.join(outdir, "theta.csv"))
    pool = read_pool(os.path.join(outdir, "perm_scan.csv"))
    py_records, py_edges, py_cycles, py_reduced = python_scan(theta, names, pool)

    with open(os.path.join(outdir, "scan_r_out.json")) as fh:
        r_out = json.load(fh)
    r_records = r_out["records"]

    print()
    print(f"{'pair':<10}{'quantity':<12}{'python':>24}{'R':>24}"
          f"{'abs diff':>12}  verdict")
    print("-" * 92)
    n_fail = 0
    n_bit = 0
    worst = (0.0, None)
    if len(py_records) != len(r_records):
        print(f"FAIL: {len(py_records)} python records vs {len(r_records)} R")
        return 1
    for pr, rr in zip(py_records, r_records):
        pair = f"{pr['source']}->{pr['target']}"
        if (pr["source"], pr["target"]) != (rr["source"], rr["target"]):
            print(f"FAIL: record order differs at {pair}")
            n_fail += 1
            continue
        for key in NUMERIC:
            pv, rv = float(pr[key]), float(rr[key])
            diff = abs(pv - rv)
            exact = key in ("p_value", "p_adj")
            ok = (pv == rv) if exact else diff <= tol
            if pv == rv:
                n_bit += 1
            if diff > worst[0]:
                worst = (diff, f"{pair}/{key}")
            if not ok:
                n_fail += 1
            print(f"{pair:<10}{key:<12}{pv:>24.17g}{rv:>24.17g}{diff:>12.3e}"
                  f"  {'PASS' + (' (exact)' if exact else '') if ok else 'FAIL'}")
        if bool(pr["edge"]) != bool(rr["edge"]):
            print(f"FAIL: edge flag differs on {pair}")
            n_fail += 1

    def norm(x):
        return [list(e) for e in x] if x is not None else None

    checks = [
        ("edge set", norm(py_edges), norm(r_out["edges"])),
        ("cycles", norm(py_cycles), norm(r_out["cycles"])),
        ("transitive reduction", norm(py_reduced), norm(r_out["reduced_edges"])),
    ]
    print("-" * 92)
    for label, a, b in checks:
        ok = a == b
        if not ok:
            n_fail += 1
        print(f"{label:<24} python={a}")
        print(f"{'':<24} R     ={b}   {'PASS' if ok else 'FAIL'}")

    print("-" * 92)
    print(f"numeric quantities compared: {len(py_records) * len(NUMERIC)}")
    print(f"bit-identical              : {n_bit}")
    print(f"failures                   : {n_fail}")
    print(f"worst absolute difference  : {worst[0]:.3e} ({worst[1]})")
    print()
    print("RESULT:", "PASS" if n_fail == 0 else f"FAIL ({n_fail})")
    return 0 if n_fail == 0 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--setup", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()
    if args.setup:
        setup(args.outdir)
        return 0
    if args.compare:
        return compare(args.outdir, args.tol)
    ap.error("choose --setup or --compare")


if __name__ == "__main__":
    sys.exit(main())
