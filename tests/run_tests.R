# run_tests.R -- the test runner for prerelation.
#
# Plain Rscript, no testthat: the package must install and test with zero
# CRAN access, so the runner is written out here.
#
#   Rscript tests/run_tests.R [golden_dir]
#
# The golden-vector check needs the committed fixtures from the Python
# package (tests/golden). They are not copied into this package -- their
# canonical home is the reference implementation. Pass the directory as an
# argument or set PRERELATION_GOLDEN. When it is absent the golden check
# reports SKIPPED and the run says so in its summary; a skipped check and
# a passing check are not allowed to look the same.

args <- commandArgs(trailingOnly = TRUE)
if (nzchar(Sys.getenv("R_PRERELATION_LIB"))) .libPaths(Sys.getenv("R_PRERELATION_LIB"))
library(prerelation)

.n_pass <- 0L
.n_fail <- 0L
.n_skip <- 0L

ok <- function(label, condition) {
  if (isTRUE(condition)) {
    .n_pass <<- .n_pass + 1L
    cat(sprintf("  PASS  %s\n", label))
  } else {
    .n_fail <<- .n_fail + 1L
    cat(sprintf("  FAIL  %s\n", label))
  }
}

near <- function(label, got, want, tol = 1e-12) {
  ok(sprintf("%s (got %.17g, want %.17g)", label, got, want),
     is.finite(got) && abs(got - want) <= tol)
}

skip <- function(label, why) {
  .n_skip <<- .n_skip + 1L
  cat(sprintf("  SKIP  %s -- %s\n", label, why))
}

errors_with <- function(label, expr, pattern) {
  res <- tryCatch({ force(expr); NULL },
                  error = function(e) conditionMessage(e))
  ok(sprintf("%s (errors, message matches '%s')", label, pattern),
     !is.null(res) && grepl(pattern, res, fixed = FALSE))
}


cat("\n[1] Benjamini-Hochberg, hand-checked\n")
# m = 5, p sorted ascending. Adjusted value at rank i is p_(i) * m / i,
# then monotonised from the largest downward:
#   rank 1: 0.001 * 5 / 1 = 0.005
#   rank 2: 0.008 * 5 / 2 = 0.020
#   rank 3: 0.039 * 5 / 3 = 0.065  -> pulled down to 0.042 by rank 5
#   rank 4: 0.041 * 5 / 4 = 0.05125 -> pulled down to 0.042
#   rank 5: 0.042 * 5 / 5 = 0.042
# so the adjusted vector is (0.005, 0.020, 0.042, 0.042, 0.042). Rank 3 is
# the case worth having: its own step-up value exceeds the value of a
# larger p, and the monotonisation is what stops the adjusted sequence
# from decreasing.
p_sorted <- c(0.001, 0.008, 0.039, 0.041, 0.042)
adj <- bh_fdr(p_sorted)
near("rank 1", adj[1], 0.005)
near("rank 2", adj[2], 0.020)
near("rank 3 (monotonised down from 0.065)", adj[3], 0.042)
near("rank 4 (monotonised down from 0.05125)", adj[4], 0.042)
near("rank 5", adj[5], 0.042)

# The same p-values out of order: the answer must follow the value, not
# the position.
perm <- c(5L, 1L, 3L, 2L, 4L)
adj2 <- bh_fdr(p_sorted[perm])
ok("order restoration", isTRUE(all.equal(adj2, adj[perm], tolerance = 0)))
ok("empty input", length(bh_fdr(numeric(0))) == 0L)
# The largest adjusted value is p_max * m / m = p_max <= 1, and
# monotonisation only lowers the others, so the clip at 1 is defensive.
ok("bounded by one", max(bh_fdr(c(0.9, 0.95, 0.99))) <= 1)


cat("\n[2] Cycles, transitive reduction, condensation\n")
nodes <- c("a", "b", "c", "d")
chain <- rbind(c("a", "b"), c("b", "c"), c("a", "c"))
ok("chain has no cycle", length(find_cycles(nodes, chain)) == 0L)
red <- transitive_reduction(nodes, chain)
ok("indirect edge a->c removed",
   nrow(red) == 2L && all(red[, 1] == c("a", "b")) && all(red[, 2] == c("b", "c")))

cyc <- rbind(c("a", "b"), c("b", "a"), c("b", "c"))
ok("two-cycle detected", length(find_cycles(nodes, cyc)) == 1L)
errors_with("reduction refuses a cyclic graph",
            transitive_reduction(nodes, cyc), "acyclic")

cond <- condense(nodes, cyc)
ok("a and b collapse into one class",
   identical(cond$classes[[1L]], c("a", "b")))
ok("c and d are singletons",
   identical(cond$classes[[2L]], "c") && identical(cond$classes[[3L]], "d"))
ok("quotient edge class1 -> class2",
   nrow(cond$quotient_edges) == 1L && all(cond$quotient_edges[1, ] == c(1L, 2L)))
ok("condensation is defined on a cyclic graph", nrow(cond$hasse_edges) == 1L)

# Duplicate parallel class edges are emitted once.
cross <- rbind(c("a", "b"), c("b", "a"), c("c", "d"), c("d", "c"),
               c("b", "c"), c("a", "d"))
ok("parallel quotient edges deduplicated",
   nrow(condense(nodes, cross)$quotient_edges) == 1L)


cat("\n[3] Coefficient behaviour at the extremes\n")
set.seed(11)
n <- 400
x <- runif(n, 0.02, 0.98)
prod_y <- x * runif(n)
r_prod <- prereq_index(x, prod_y)
ok("product form gives a positive PI", r_prod$PI > 0.4)
ok("PI is in [0, 1]", r_prod$PI >= 0 && r_prod$PI <= 1)
r_rev <- prereq_index(prod_y, x)
ok("reverse direction is exactly zero", r_rev$PI == 0)
d <- direction(x, prod_y)
near("Delta equals the difference of the pair",
     d$delta_stat, r_prod$PI - r_rev$PI, tol = 0)

r_eq <- prereq_index(x, x)
ok("equivalence: A1 = 1 but A2 = 0", r_eq$A1 == 1 && r_eq$A2 == 0)
ok("equivalence: ell = 0 is what kills PI", r_eq$ell == 0 && r_eq$PI == 0)

r_ind <- prereq_index(x, runif(n))
ok("independence: A1 = 0", r_ind$A1 == 0 && r_ind$PI == 0)

# The MIN_INTERIOR guard: fewer than max(10, 0.05 n) interior points
# forces q = 0. n = 40 needs 10.
xs <- runif(40, 0.5, 0.6)
ok("interior guard forces q = 0 under equivalence",
   prereq_index(xs, xs)$q == 0)

ok("dense and sorted baselines agree",
   abs(baseline_mean(x, prod_y, dense_max_n = 10000) -
       baseline_mean(x, prod_y, dense_max_n = 1)) < 1e-12)

errors_with("unequal lengths rejected", prereq_index(1:3, 1:4), "equal length")
errors_with("empty input rejected", prereq_index(numeric(0), numeric(0)),
            "non-empty")


cat("\n[4] Permutation indices are read as data, one-based\n")
P0 <- rbind(c(0L, 1L, 2L), c(2L, 0L, 1L))
tmp <- tempfile(fileext = ".csv")
writeLines(apply(P0, 1L, paste, collapse = ","), tmp)
P1 <- read_perm_indices(tmp)
ok("zero-based file becomes one-based matrix",
   identical(P1, matrix(c(1L, 2L, 3L, 3L, 1L, 2L), nrow = 2L, byrow = TRUE)))
errors_with("zero-based indices rejected downstream",
            prerelation:::.check_perm_indices(P0, 3L), "one-based")

# The repr-style CSV trap: a wrapped float literal must stop the read where
# the bad byte is, not produce a character column that fails later.
bad <- tempfile(fileext = ".csv")
writeLines(c("x,y", "np.float64(0.5),0.25"), bad)
errors_with("wrapped float literal rejected at read time",
            read_fixture(bad), "not numeric")


cat("\n[5] Golden vectors\n")
golden <- if (length(args) >= 1L) args[1L] else Sys.getenv("PRERELATION_GOLDEN")
if (!nzchar(golden) || !dir.exists(golden)) {
  skip("golden fixtures", "no golden directory supplied")
} else {
  # Provenance: the committed expected.json of prerelation (Python), the
  # reference implementation. Contract: 1e-12 on components, exact on the
  # permutation p-value.
  expected <- list(
    product = list(PI = 0.6169018554399411, A1 = 0.6939863400902884,
                   A2 = 0.8889250692739594, q = 0.9515254262650833,
                   ell = 0.9342105263157895, perm_p = 0.005),
    min = list(PI = 0.26202393315496875, A1 = 0.34354474546234226,
               A2 = 0.7627068573042418, q = 0.9502577238544653,
               ell = 0.8026315789473684, perm_p = 0.005),
    independent = list(PI = 0.0, A1 = 0.0, A2 = 0.8768248536237636,
                       q = 0.9519812696486577, ell = 0.9210526315789473,
                       perm_p = 1.0),
    equivalence = list(PI = 0.0, A1 = 1.0, A2 = 0.0, q = 0.0, ell = 0.0,
                       perm_p = 1.0),
    partial_equivalence = list(PI = 0.13863815253871326,
                               A1 = 0.6575998915009345,
                               A2 = 0.21082447599296214,
                               q = 0.26704433625775204,
                               ell = 0.7894736842105263, perm_p = 0.005),
    ecpe_slice = list(PI = 0.20408619298974254, A1 = 0.4385386661965045,
                      A2 = 0.4653778759346472, q = 0.4653778759346472,
                      ell = 1.0, perm_p = 0.005)
  )
  for (name in names(expected)) {
    f <- read_fixture(file.path(golden, sprintf("fixture_%s.csv", name)))
    res <- prereq_index(f$x, f$y)
    for (key in c("PI", "A1", "A2", "q", "ell")) {
      near(sprintf("%s/%s", name, key), res[[key]], expected[[name]][[key]])
    }
    P <- read_perm_indices(file.path(golden,
                                     sprintf("perm_indices_n%d.csv", length(f$x))))
    pp <- perm_pvalue(f$x, f$y, perm_indices = P)
    ok(sprintf("%s/perm_p exact (%.17g)", name, pp$p_value),
       pp$p_value == expected[[name]]$perm_p)
  }
}

cat("\n[6] Reference class and envelope (base R, pbeta)\n")
# Admissibility is the class B, a ceiling condition -- not dominance.
a_u <- admissibility(uniform_reference())
ok("Uniform is admissible (boundary point)", isTRUE(a_u$admissible))
near("Uniform tail mass is delta", a_u$tail_mass, 0.05)
a_21 <- admissibility(beta_reference(2, 1))
ok("Beta(2,1) rejected (tail 0.0975)", !isTRUE(a_21$admissible) && abs(a_21$tail_mass - 0.0975) < 5e-5)
a_82 <- admissibility(beta_reference(8, 2))
ok("Beta(8,2) rejected (tail 0.0712)", !isTRUE(a_82$admissible) && abs(a_82$tail_mass - 0.0712) < 5e-5)
F210 <- beta_reference(2, 10)
ok("Beta(2,10) fails dominance at the floor: F(0.01) < 0.01", F210(0.01) < 0.01)
a_210 <- admissibility(F210)
ok("Beta(2,10) ACCEPTED by B (tail 1.03e-12)", isTRUE(a_210$admissible) && a_210$tail_mass < 1e-11)
ok("point mass at 0 is admissible", isTRUE(admissibility(point_mass_reference(0))$admissible))

# E-6a case 3c: the only tail-band point sits below its index, so the
# positive part in D* is what keeps sup_q at 1 rather than above it.
m3 <- 100L
t3 <- c(seq(0.005, 0.90, length.out = m3 - 1L), 0.96)
x3 <- rep(0.5, m3)
y3 <- x3 * t3 * (1 - 0.05)
e3 <- pi_envelope(x3, y3)
ok("case 3c: n_tail == 1", e3$n_tail == 1L)
ok("case 3c: D_star == 0", e3$D_star == 0)
ok("case 3c: sup_q == 1", e3$sup_q == 1)

if (nzchar(golden) && dir.exists(golden)) {
  # Provenance: the 30 reference-class keys added to expected.json by E-6a
  # (n_tail_band, D_star, sup_q, inf_q, PI_hi per fixture). Tolerance 1e-12.
  expected_ref <- list(
    product = list(n_tail_band = 21L, D_star = 0.024720479594442346,
      sup_q = 0.9752795204055577, inf_q = 0.003676470588235294, PI_hi = 0.6323023317121026),
    min = list(n_tail_band = 11L, D_star = 0.025073797152681543,
      sup_q = 0.9749262028473185, inf_q = 0.005952380952380952, PI_hi = 0.268826016135615),
    independent = list(n_tail_band = 15L, D_star = 0.02444612128830015,
      sup_q = 0.9755538787116999, inf_q = 0.005154639175257732, PI_hi = 0.0),
    equivalence = list(n_tail_band = 0L, D_star = 0.0,
      sup_q = 0.0, inf_q = 0.0, PI_hi = 0.0),
    partial_equivalence = list(n_tail_band = 28L, D_star = 0.10006445443418799,
      sup_q = 0.899935545565812, inf_q = 0.005555555555555556, PI_hi = 0.4672085661488781),
    ecpe_slice = list(n_tail_band = 21L, D_star = 0.11872185619096265,
      sup_q = 0.8812781438090374, inf_q = 0.008403361344537815, PI_hi = 0.3864745417341465)
  )
  for (name in names(expected_ref)) {
    f <- read_fixture(file.path(golden, sprintf("fixture_%s.csv", name)))
    env <- pi_envelope(f$x, f$y)
    ok(sprintf("%s/n_tail_band (%d)", name, env$n_tail),
       env$n_tail == expected_ref[[name]]$n_tail_band)
    for (key in c("D_star", "sup_q", "inf_q", "PI_hi")) {
      near(sprintf("%s/%s", name, key), env[[key]], expected_ref[[name]][[key]])
    }
    ok(sprintf("%s/attained", name), isTRUE(env$attained))
    core <- prereq_index(f$x, f$y)
    ok(sprintf("%s/interior_q at Uniform == core q (bitwise)", name),
       interior_q(f$x, f$y) == core$q)
    ok(sprintf("%s/family member at Uniform == PI (bitwise)", name),
       prereq_index_family(f$x, f$y)$PI == core$PI)
    if (env$m > 0L) {
      Fs <- attaining_reference(f$x, f$y)
      ok(sprintf("%s/F0* admissible", name), isTRUE(admissibility(Fs)$admissible))
      near(sprintf("%s/q(F0*) attains sup_q", name), interior_q(f$x, f$y, Fs), env$sup_q)
      near(sprintf("%s/q(point mass) == inf_q", name),
           interior_q(f$x, f$y, point_mass_reference(0)), env$inf_q)
      for (F in list(F210, beta_reference(1, 2), point_mass_reference(0))) {
        ok(sprintf("%s/q(%s) <= sup_q", name, attr(F, "reference_name")),
           interior_q(f$x, f$y, F) <= env$sup_q + 1e-12)
      }
    }
  }
} else {
  skip("reference-class golden keys", "no golden directory supplied")
}

cat(sprintf("\npassed %d, failed %d, skipped %d\n", .n_pass, .n_fail, .n_skip))
if (.n_skip > 0L) cat("NOTE: skipped checks are not passes.\n")
if (.n_fail > 0L) quit(status = 1L)
