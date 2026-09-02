# prerelation (R)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22135973.svg)](https://doi.org/10.5281/zenodo.22135973)

An R implementation of the **prerelation coefficient** — a coefficient for
prerequisite relations between traits reported on a common anchored scale —
held in parity with the Python reference implementation
[`prerelation`](https://github.com/JaehwaFChoi/prerelation) and with the
JavaScript implementation
[`prerelation-js`](https://github.com/JaehwaFChoi/prerelation-js).

```
Pi(X -> Y) = A1 * A2        in [0, 1]

A1   the corner {Y > X} is empty, relative to what independence would give
A2   below the ceiling Y varies as a free component should (q), and the
     censoring thins out at high x (ell)

Delta = Pi(X -> Y) - Pi(Y -> X)
```

The product structure is what lets a single number separate the four
extremes. Independence is annihilated by `A1`; exact equivalence is
annihilated by `A2`.

## Scope of this implementation

This package is **version 0.1.0** while the Python and JavaScript packages
are at 0.2.0. The version numbers differ on purpose: this port covers the
coefficient and the scan, and does not cover the rest of the reference
package. What it contains:

| | |
|---|---|
| **present** | `prereq_index` (with all five components), `direction`, `perm_pvalue`, the pairwise scan with Benjamini-Hochberg control, cycle detection, transitive reduction, `condense` (strongly connected components to a quotient order), the golden-vector readers; the admissible reference class (`admissibility`), the interior component at a declared reference (`interior_q`), the family member (`prereq_index_family`) and the exact upper envelope (`pi_envelope`), with `uniform_reference`, `beta_reference` (base `pbeta`), `point_mass_reference` and `attaining_reference` |
| **absent** | `ceiling_fit` and the postulate correction, `pv_correct` (plausible-value handling — there is no R counterpart for the grid posteriors it consumes), `run_study` |

The absent pieces are absent by decision, not by oversight. Use the Python
package for them.

## How to read the two scales

**Delta — the prerelation direction coefficient.**
`Delta = Pi(X -> Y) - Pi(Y -> X)` lies in `[-1, +1]` and is antisymmetric.
Its **sign is the direction** — positive means X is the prerequisite side —
and its **magnitude is the strength of the asymmetry**.

*Read Delta together with the Pi pair:* `Delta = 0` by itself does not
distinguish "no relation" from "equivalent skills" — both put the two
directions on an equal footing.

**Pi — per-direction strength.** `Pi` lies in `[0, 1]`: **0 means no
prerequisite relation, 1 means a perfect prerequisite relation**. It is a
continuous quantity, read like a correlation magnitude; the package defines
no thresholds and no cutoffs.

**The reading ladder.** Permutation p (is there a relation at all) ->
Delta (which direction, how asymmetric) -> the Pi pair (per-direction
strength).

## Relation to the correlation coefficient

Pearson r answers a symmetric question: do X and Y move together?
Prerequisite-ness is asymmetric: does progress in Y require X first? A high
r cannot separate X -> Y from Y -> X, nor either from "both reflect one
shared ability", and two nearly identical skills correlate almost perfectly
while neither is a prerequisite for the other. Pi scores the one-sided
ceiling footprint instead, which is why the equivalence case splits the two
apart: r is close to 1 while `Delta = 0` and both `Pi = 0`. That is a
property of the definitions. The two coefficients answer different
questions; Pi complements r rather than replacing it.

Anchored scales are an interpretability requirement, not a claim about the
measurement precision of any scoring model: the ratio `Y / X` and the
corner moment `(Y - X)_+` only carry the reading "how much of the ceiling
granted by X is used by Y" when both endpoints are substantive anchors. On
an unanchored scale Pi carries no prerequisite interpretation.

## Install

Base R only. No `Imports`, no `Suggests`, nothing from CRAN.

```r
# from a local clone
R CMD INSTALL .

# or, from R
install.packages("remotes")
remotes::install_github("JaehwaFChoi/prerelation-r")
```

## Quick start

```r
library(prerelation)

x <- c(...)   # trait values in [0, 1]
y <- c(...)   # same length, same scale

prereq_index(x, y)   # list(PI, A1, A2, q, ell)
direction(x, y)      # list(delta_stat, forward, reverse)
perm_pvalue(x, y, n_perm = 999, seed = 20260827)

theta <- matrix(...)  # rows are persons, columns attributes
res <- prereq_scan(theta, attr_names = c("A", "B", "C"), n_perm = 199)
res$edges                # pairs surviving BH-FDR control
res$reduced_edges        # transitive reduction, or NULL when cyclic
res$equivalence_classes  # mutually dominating attributes, condensed
```

`prereq_index` returns the same five components as the Python and
JavaScript implementations (`PI`, `A1`, `A2`, `q`, `ell`), so results
transfer between the three without renaming.

### Design floor on permutation replicates

With `k` attributes there are `K = k (k - 1)` ordered pairs, and the
smallest attainable permutation p-value is `1 / (n_perm + 1)`. For any pair
to survive Benjamini-Hochberg control at level `alpha`, the replicate count
must satisfy

```
n_perm >= K / alpha - 1
```

(`K = 6`, `alpha = 0.05` needs `n_perm >= 119`; `K = 56` needs
`n_perm >= 1119`). Below the floor the scan cannot return any edge,
regardless of the data.

## What the scan recovers

The edge set — and its transitive reduction, and the quotient order the
equivalence classes condense to — is a **dominance preorder** over the
attributes: which attributes act as ceilings on which others. It is not a
direct-prerequisite acyclic graph. Indirect dominance produces edges of its
own, and siblings under a common ceiling can be linked to each other even
though neither is a prerequisite for the other. Directed cycles, and the
merged nodes they condense into, are expected behaviour of a pairwise
index. A disagreement between the recovered order and an expert-specified
prerequisite graph is a difference between two concepts, not by itself an
error in either.

## API mapping across the three implementations

Names differ where a Python or JavaScript name would mask a base R
function. **The numerical definition behind each row is the same in all
three implementations**, and the parity harness below is what holds them
there.

| R | Python | JavaScript | numerical definition |
|---|---|---|---|
| `prereq_index(x, y)` | `prereq_index(x, y)` | `prereqIndex(x, y)` | unchanged |
| `direction(x, y)` | `direction(x, y)` | `direction(x, y)` | unchanged |
| `perm_pvalue(x, y, ...)` | `perm_pvalue(x, y, ...)` | `permPvalue(x, y, ...)` | unchanged |
| `prereq_scan(theta, ...)` | `scan(theta, ...)` | `scan(theta, ...)` | unchanged |
| `attr_names =` | `names =` | `names:` | unchanged |
| `baseline_mean(x, y)` | `_baseline_mean(x, y)` | (internal) | unchanged |
| `bh_fdr(p)` | `bh_fdr(p)` | `bhFdr(p)` | unchanged |
| `find_cycles(...)` | `find_cycles(...)` | `findCycles(...)` | unchanged |
| `transitive_reduction(...)` | `transitive_reduction(...)` | `transitiveReduction(...)` | unchanged |
| `condense(nodes, edges)` | — (not in 0.2.0) | `condense(nodes, edges)` | unchanged |
| `DELTA`, `TOP_Q`, `MIN_INTERIOR` | same names | same names | frozen constants, identical values |

`scan` and `names` are base R functions, so `prereq_scan` and `attr_names`
are used here to avoid masking them. The renaming is confined to the
surface: the arguments, the returned components and every number are the
same.

## Parity with the reference implementation

The Python package is the reference implementation and its committed
oracle fixes the definition. This package is checked against the golden
vectors of the published release **prerelation 0.2.0** (concept DOI
[10.5281/zenodo.22132819](https://doi.org/10.5281/zenodo.22132819)) and
against **prerelation-js 0.2.0** (concept DOI
[10.5281/zenodo.22133624](https://doi.org/10.5281/zenodo.22133624)).

`tools/parity/` holds the harness. Both languages read the same fixture
CSV and the same permutation index matrix from disk, so neither side is
comparing against numbers it produced itself. Measured on the current tree:

| layer | what is compared | result |
|---|---|---|
| closed-form components (`v`, `v0`, `A1`, band masses, `q`, `ell`, `A2`, `PI`, `PI_reverse`, `Delta`) over the six golden fixtures | 90 quantities | 0 failures, 74 bit-identical, worst absolute difference 1.665e-16 |
| pairwise scan, 4 attributes, n = 200, 12 ordered pairs, `n_perm = 239` | 108 numeric quantities | 0 failures, 90 bit-identical, worst absolute difference 2.220e-16; p-values and Benjamini-Hochberg adjusted p-values agree exactly on all 12 pairs; edge set, cycle report and transitive reduction identical |
| `condense` against the JavaScript reference over six graphs | canonical rendering | byte-identical |

The agreement criterion is `1e-12` absolute on the component quantities,
not bit-identity: the independence baseline is a double sum and the two
languages accumulate it in a different order, so a difference of a few
units in the last place is expected and passes. The summation is
deliberately not rewritten to chase the last bit.

**On the permutation p-values.** The strong evidence is the scan layer,
where the twelve ordered pairs carry interior p-values and a perturbed
permutation pool fails immediately. On the six golden fixtures the
permutation p-values are saturated at the ends of their range, so their
exact agreement is a *consistency* check — it shows the shared index
matrices were applied correctly — and is not by itself evidence about
which indices were consumed.

Note on random numbers: a seeded run in one language does not reproduce
the seeded run of another. Cross-implementation p-value equality is
defined only through shared index matrices read from disk.

## Repository layout

```
DESCRIPTION, NAMESPACE, LICENSE
R/            core.R, scan.R, golden.R
tests/        run_tests.R          (plain Rscript runner, no test framework)
tools/parity/ parity_r.R, parity_driver.py, scan_r.R, scan_driver.py,
              condense_r.R, condense_js.mjs, graphs.json,
              korean_gate.py, destruction.sh
```

`tests/run_tests.R` takes the golden directory as an argument or from
`PRERELATION_GOLDEN`. The golden fixtures are deliberately not copied into
this package: their canonical home is the reference repository, and a
second copy could drift. Without them the runner reports a skip and says
plainly that a skipped check is not a pass.

## Citation

Cite the reference implementation for the method and this package for the
R port.

| | concept DOI (latest version) |
|---|---|
| `prerelation-r` (this package) | [10.5281/zenodo.22135973](https://doi.org/10.5281/zenodo.22135973) |
| `prerelation` (Python, reference) | [10.5281/zenodo.22132819](https://doi.org/10.5281/zenodo.22132819) |
| `prerelation-js` (JavaScript) | [10.5281/zenodo.22133624](https://doi.org/10.5281/zenodo.22133624) |

Machine-readable metadata is in `CITATION.cff`.

## License

MIT. See `LICENSE.md`; `LICENSE` is the two-line form R's
`License: MIT + file LICENSE` field requires.
