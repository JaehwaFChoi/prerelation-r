# parity_r.R -- R side of the cross-language parity harness.
#
# Computes, for every committed golden fixture, the twenty quantities that
# expected.json records (fifteen core, five reference-class), and writes them as full-precision text (%.17g) so
# that the Python driver reads back exactly the float64 R produced. The
# fixture CSV and the permutation index matrix are read from the same files
# the Python side reads, so both languages consume identical bytes.
#
# Usage: Rscript parity_r.R <golden_dir> <out_csv>

args <- commandArgs(trailingOnly = TRUE)
golden <- args[1]
out_csv <- args[2]

if (nzchar(Sys.getenv("R_PRERELATION_LIB"))) .libPaths(Sys.getenv("R_PRERELATION_LIB"))
library(prerelation)

FIXTURES <- c("product", "min", "independent", "equivalence",
              "partial_equivalence", "ecpe_slice")
DELTA_ <- 0.05
TOP_Q_ <- 0.8

quantile_type7 <- prerelation:::.quantile_type7

components <- function(x, y, P) {
  n <- length(x)
  res <- prereq_index(x, y)
  rev <- prereq_index(y, x)
  dl <- direction(x, y)

  v <- mean(pmax(y - x, 0))
  v0 <- baseline_mean(x, y)

  # The auxiliary masses are recomputed here from the definition in
  # core.R, with the same 1e-9 denominator floor the coefficient uses.
  u <- pmin(pmax(y / pmax(x, 1e-9), 0), 1)
  ceil_mask <- u >= 1 - DELTA_
  x_top <- x >= quantile_type7(x, TOP_Q_)
  p1_top <- if (sum(x_top) > 0) mean(ceil_mask[x_top]) else 1

  pp <- perm_pvalue(x, y, perm_indices = P)
  env <- pi_envelope(x, y)

  list(n = n, v = v, v0 = v0, A1 = res$A1,
       mass_ceiling_band = mean(ceil_mask),
       mass_interior = mean(!ceil_mask),
       n_interior = sum(!ceil_mask),
       p1_top = p1_top, q = res$q, ell = res$ell, A2 = res$A2,
       PI = res$PI, PI_reverse = rev$PI, Delta = dl$delta_stat,
       perm_p = pp$p_value,
       # reference class and envelope (E-6a golden keys)
       n_tail_band = env$n_tail, D_star = env$D_star, sup_q = env$sup_q,
       inf_q = env$inf_q, PI_hi = env$PI_hi)
}

lines <- c("fixture,quantity,value")
for (name in FIXTURES) {
  f <- read_fixture(file.path(golden, sprintf("fixture_%s.csv", name)))
  n <- length(f$x)
  P <- read_perm_indices(file.path(golden, sprintf("perm_indices_n%d.csv", n)))
  comp <- components(f$x, f$y, P)
  for (q in names(comp)) {
    lines <- c(lines, sprintf("%s,%s,%.17g", name, q, comp[[q]]))
  }
  cat(sprintf("computed %-20s n=%d PI=%.17g\n", name, n, comp$PI))
}

writeLines(lines, out_csv)
cat(sprintf("wrote %s (%d rows)\n", out_csv, length(lines) - 1L))
