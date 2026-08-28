# prerelation/core.R -- the prerelation coefficient and its direction statistic.
#
# Pi(X -> Y) = A1 * A2
#
#   A1  = max(0, 1 - v / v0)   emptiness of the corner {Y > X}, measured
#                              against the independence baseline v0
#   A2  = q * ell              conditional freedom, censoring-aware
#   q                          uniformity of the interior of
#                              W = min(Y / X, 1) below the ceiling band
#   ell                        legitimacy of the ceiling: censoring must
#                              thin out at high x, which is what separates
#                              a genuine ceiling from equivalence
#
#   Delta = Pi(X -> Y) - Pi(Y -> X).
#
# Anchored scales are an interpretability requirement, not a claim about the
# measurement precision of any scoring model. On an unanchored scale Pi
# carries no prerequisite interpretation.
#
# Correctness standard
# --------------------
# This file is a transcription of the Python oracle
# `tests/oracle/prereq_index_v2.py` in the `prerelation` package, made with
# the source open. Every branch of the definition is transcribed literally.
# The following are part of the *definition*, not implementation detail:
#
#   * the independence baseline is a V-statistic -- the double sum runs over
#     all n^2 ordered pairs, the diagonal i = j included;
#   * v0 <= 1e-9 forces A1 = 0;
#   * fewer than max(10, 0.05 n) interior points forces q = 0;
#   * an empty top-x stratum forces p1_top = 1 (hence ell = 0);
#   * the ratio u = y / x is clipped to [0, 1] and its denominator floored
#     at 1e-9;
#   * the ceiling band is u >= 1 - delta.
#
# Only the summation order of the baseline is allowed to differ. The
# defaults delta = 0.05, TOP_Q = 0.8 and MIN_INTERIOR = max(10, 0.05 n) are
# fixed and are not adjustable by an implementation.

DELTA <- 0.05          # ceiling band width
TOP_Q <- 0.8           # top-x quantile for the legitimacy check
MIN_INTERIOR <- 10     # minimum interior points before freedom is credited

.EPS_DEN <- 1e-9       # floor of the ratio denominator
.EPS_V0 <- 1e-9        # guard on the independence baseline

# Above this sample size the baseline is accumulated by sorting instead of
# forming the n x n matrix. The two paths differ only in summation order;
# below the threshold the dense expression of the oracle is used.
DENSE_MAX_N <- 3000


# Local replacements for the two base-package helpers this port would
# otherwise borrow. The package declares no Imports at all (D17: zero CRAN
# dependencies), so nothing outside base is referenced.

.set_names <- function(v, nm) {
  names(v) <- nm
  v
}

# Type-7 quantile, written out rather than delegated to stats::quantile.
# This is NumPy's "linear" method, including its lerp branch at gamma >= 0.5
# (numpy.lib.function_base._lerp): the branch changes the last bit of the
# result, and the threshold it produces is compared against the data with
# `>=`, so transcribing it removes a way for the top-x stratum to differ
# between the two implementations on tied data.
.quantile_type7 <- function(x, p) {
  xs <- sort(x)
  n <- length(xs)
  if (n == 1L) return(xs[1L])
  h <- (n - 1) * p
  lo <- floor(h)
  gamma <- h - lo
  a <- xs[lo + 1L]
  b <- if (lo + 1L < n) xs[lo + 2L] else xs[n]
  d <- b - a
  if (gamma >= 0.5) b - d * (1 - gamma) else a + d * gamma
}


.as_pair <- function(x, y) {
  x <- as.double(x)
  y <- as.double(y)
  if (!is.null(dim(x)) || !is.null(dim(y))) {
    stop("x and y must be one-dimensional")
  }
  if (length(x) != length(y)) {
    stop(sprintf("x and y must have equal length, got %d and %d",
                 length(x), length(y)))
  }
  if (length(x) == 0L) {
    stop("x and y must be non-empty")
  }
  list(x = x, y = y)
}


#' Independence baseline v0 = mean over all n^2 ordered pairs of (y_j - x_i)_+.
#'
#' The diagonal is included: this is a V-statistic, not a U-statistic.
baseline_mean <- function(x, y, dense_max_n = DENSE_MAX_N) {
  n <- length(x)
  if (n <= dense_max_n) {
    # Dense expression of the oracle. outer(y, x, "-")[j, i] = y_j - x_i;
    # the mean over the whole matrix is the same quantity whichever way
    # round the matrix is laid out.
    return(mean(pmax(outer(y, x, "-"), 0)))
  }

  # Sort-and-accumulate: for each x_i, sum_j (y_j - x_i)_+ equals the sum of
  # the y above x_i minus x_i times how many there are. O(n log n) instead
  # of O(n^2), which keeps large scans and permutation loops feasible.
  ys <- sort(y)
  tail_sum <- c(rev(cumsum(rev(ys))), 0)  # tail_sum[m] = sum(ys[m:n])
  k <- findInterval(x, ys)                # count of y_j <= x_i
  inner <- tail_sum[k + 1L] - x * (n - k)
  mean(inner) / n
}


#' Prerelation coefficient of the ordered pair (x -> y).
#'
#' @param x,y numeric vectors of length n; trait values on a common anchored
#'   scale, x the candidate prerequisite. Values are expected in [0, 1].
#' @param delta ceiling band width; the default 0.05 is the fixed convention.
#' @param dense_max_n sample size up to which the independence baseline is
#'   formed as a dense n x n matrix.
#' @return list with elements PI, A1, A2, q, ell.
#'
#' The statistic is deliberately not invariant to monotone rescaling of
#' either axis, and it is not symmetric: prereq_index(y, x) answers a
#' different question.
prereq_index <- function(x, y, delta = DELTA, dense_max_n = DENSE_MAX_N) {
  pair <- .as_pair(x, y)
  x <- pair$x
  y <- pair$y
  n <- length(x)

  # A1: corner emptiness relative to the independence baseline.
  v <- mean(pmax(y - x, 0))
  v0 <- baseline_mean(x, y, dense_max_n)
  a1 <- if (v0 > .EPS_V0) max(0, 1 - v / v0) else 0

  # A2: conditional freedom with a censoring-aware benchmark.
  u <- pmin(pmax(y / pmax(x, .EPS_DEN), 0), 1)
  ceil_mask <- u >= 1 - delta
  interior <- u[!ceil_mask]

  if (length(interior) < max(MIN_INTERIOR, 0.05 * n)) {
    q <- 0
  } else {
    tt <- sort(interior / (1 - delta))
    ecdf_vals <- seq_len(length(tt)) / length(tt)
    q <- 1 - max(abs(ecdf_vals - tt))
  }

  x_top <- x >= .quantile_type7(x, TOP_Q)
  p1_top <- if (sum(x_top) > 0) mean(ceil_mask[x_top]) else 1
  ell <- 1 - max(0, p1_top - delta) / (1 - delta)

  a2 <- q * ell
  list(PI = a1 * a2, A1 = a1, A2 = a2, q = q, ell = ell)
}


#' Directional contrast of the pair.
#'
#' @return list with elements delta_stat = pi_xy - pi_yx, pi_xy, pi_yx. The
#'   element order follows the oracle so that parity tests can compare
#'   element by element.
direction <- function(x, y, delta = DELTA, dense_max_n = DENSE_MAX_N) {
  pi_xy <- prereq_index(x, y, delta, dense_max_n)$PI
  pi_yx <- prereq_index(y, x, delta, dense_max_n)$PI
  list(delta_stat = pi_xy - pi_yx, pi_xy = pi_xy, pi_yx = pi_yx)
}


#' Permutation test of independence for the forward statistic.
#'
#' Under the null of independence the joint law of the sample is invariant
#' under permutations of the y-labels, so the test is exact for the full
#' group; the Monte-Carlo version inherits validity because the observed
#' configuration is counted in the reference set (the add-one rule).
#'
#' @param perm_indices optional integer matrix of ONE-BASED permutation
#'   indices, one permutation of 1..n per row. Cross-language parity is
#'   only defined through this argument: R's generator cannot reproduce the
#'   NumPy stream that the Python package consumes, so the shared
#'   permutations are read as data (see read_perm_indices). When it is
#'   supplied, n_perm is the row count of the matrix and seed is ignored.
#' @param n_perm,seed used only when perm_indices is NULL, in which case the
#'   reference set comes from R's own generator and is NOT comparable to the
#'   Python or JavaScript implementations.
#' @return list with elements observed and p_value.
#'
#' The count uses >=, which makes it sensitive to exact ties; ties occur in
#' degenerate configurations where many replicates give an identical value
#' (typically 0).
perm_pvalue <- function(x, y, n_perm = 1000, seed = 0, delta = DELTA,
                        dense_max_n = DENSE_MAX_N, perm_indices = NULL) {
  pair <- .as_pair(x, y)
  x <- pair$x
  y <- pair$y
  n <- length(x)
  obs <- prereq_index(x, y, delta, dense_max_n)$PI

  if (is.null(perm_indices)) {
    set.seed(seed)
    cnt <- 0L
    for (r in seq_len(n_perm)) {
      yp <- y[sample.int(n)]
      if (prereq_index(x, yp, delta, dense_max_n)$PI >= obs) cnt <- cnt + 1L
    }
    return(list(observed = obs, p_value = (cnt + 1) / (n_perm + 1)))
  }

  P <- .check_perm_indices(perm_indices, n)
  n_perm <- nrow(P)
  cnt <- 0L
  for (r in seq_len(n_perm)) {
    yp <- y[P[r, ]]
    if (prereq_index(x, yp, delta, dense_max_n)$PI >= obs) cnt <- cnt + 1L
  }
  list(observed = obs, p_value = (cnt + 1) / (n_perm + 1))
}


.check_perm_indices <- function(P, n) {
  P <- as.matrix(P)
  if (!is.numeric(P)) {
    stop("perm_indices must be numeric; a character matrix usually means the ",
         "index file was read with a non-numeric column")
  }
  if (ncol(P) != n) {
    stop(sprintf("perm_indices has %d columns but the sample has %d rows",
                 ncol(P), n))
  }
  if (any(is.na(P)) || min(P) < 1L || max(P) > n) {
    stop("perm_indices must hold one-based indices in 1..n; ",
         "the committed matrices are zero-based and need read_perm_indices")
  }
  storage.mode(P) <- "integer"
  P
}
