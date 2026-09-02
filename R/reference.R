# prerelation/reference.R -- the admissible reference class and the exact
# upper envelope. Port of prerelation/reference.py (the Python reference
# implementation), checked against the same golden vectors. Base R only.
#
# The interior component q compares the rescaled interior of u = y / x
# with a reference law on [0, 1]; the definition fixes that reference at
# Uniform(0, 1). These functions generalise the reference to an exogenously
# DECLARED law F0 and supply:
#
#   * the admissible reference class
#         B = { F0 on [0,1] : F0(t) >= t for all t in [1 - delta, 1] },
#     references placing no more mass near the ceiling than Uniform does.
#     Uniform meets the condition with equality (the boundary point of B).
#     Membership is strictly weaker than first-order stochastic dominance by
#     Uniform: Beta(2, 10) fails dominance at the floor and belongs to B.
#
#   * the exact upper envelope of q over B:
#         sup_{F0 in B} q(F0) = 1 - D*,
#         D* = max { (t_(i) - i/m)_+ : t_(i) >= 1 - delta }   (0 if empty),
#     attained (for distinct t) by F0*(t) = max(ECDF_m(t), t 1{t >= 1-delta}),
#     itself a member of B; PI_hi = A1 * ell * (1 - D*). With tied t values
#     1 - D* remains a valid upper bound for every member of B but need not
#     be attained; `attained` reports which case holds.
#
#   * the vacuous lower end: the point mass at 0 is admissible and gives
#     q = 1/m exactly, a function of the interior sample size alone.
#
#   * the family member PI(F0) = A1 * (q(F0) * ell), composed from the
#     definition's own components (A1 and ell do not depend on F0).
#
# Both uses of delta are on the rescaled scale t: the interior is
# u < 1 - delta, then t = u / (1 - delta), and the admissibility threshold
# t >= 1 - delta is applied to t (equivalently u >= (1 - delta)^2).
#
# A reference is an R function F0(t) taking a numeric vector in [0, 1] and
# returning distribution-function values. The reference must be declared
# before the data are seen and never fitted from the same data.

.REF_TOL <- 1e-12


# --------------------------------------------------------------------------
# reference constructors
# --------------------------------------------------------------------------

#' The package default F0(t) = t: the boundary point of B.
uniform_reference <- function() {
  F0 <- function(t) pmin(pmax(as.numeric(t), 0), 1)
  attr(F0, "reference_name") <- "Uniform(0,1)"
  F0
}

#' Beta(a, b) distribution function as a reference (pbeta is base R).
beta_reference <- function(a, b) {
  if (!(a > 0 && b > 0)) stop("Beta parameters must be positive")
  F0 <- function(t) pbeta(pmin(pmax(as.numeric(t), 0), 1), a, b)
  attr(F0, "reference_name") <- sprintf("Beta(%g,%g)", a, b)
  F0
}

#' Degenerate law at `at`: F0(t) = 1{t >= at}. The point mass at 0 attains
#' the (vacuous) lower end q = 1/m.
point_mass_reference <- function(at = 0) {
  F0 <- function(t) as.numeric(as.numeric(t) >= at)
  attr(F0, "reference_name") <- sprintf("PointMass(%g)", at)
  F0
}


# --------------------------------------------------------------------------
# the rescaled interior, exactly as core forms it
# --------------------------------------------------------------------------

.rescaled_interior <- function(x, y, delta) {
  pair <- .as_pair(x, y)
  x <- pair$x
  y <- pair$y
  u <- pmin(pmax(y / pmax(x, .EPS_DEN), 0), 1)
  ceil_mask <- u >= 1 - delta
  interior <- u[!ceil_mask]
  list(t = sort(interior / (1 - delta)), n = length(x))
}

.guard_fires <- function(m, n) m < max(MIN_INTERIOR, 0.05 * n)

# ECDF of a sorted vector evaluated at v: #{t_i <= v} / m
.ecdf_sorted <- function(t, v) findInterval(v, t) / length(t)

#' The member of B that attains the supremum for this pair (distinct t):
#' F0*(t) = max(ECDF_m(t), t 1{t >= 1 - delta}).
attaining_reference <- function(x, y, delta = DELTA) {
  ri <- .rescaled_interior(x, y, delta)
  t_sorted <- ri$t
  if (length(t_sorted) == 0L) stop("no interior points")
  lo <- 1 - delta
  F0 <- function(v) {
    v <- as.numeric(v)
    pmax(.ecdf_sorted(t_sorted, v), ifelse(v >= lo, v, 0))
  }
  attr(F0, "reference_name") <- "F0star"
  F0
}


# --------------------------------------------------------------------------
# admissibility
# --------------------------------------------------------------------------

#' Does the declared reference belong to the admissible class B?
#'
#' Checks the defining pointwise condition F0(t) >= t on `n_grid` equally
#' spaced points of [1 - delta, 1] (both endpoints), NOT first-order
#' stochastic dominance: nothing is checked below 1 - delta. The tail mass
#' 1 - F0(1 - delta) is a consequence of admissibility (the condition at
#' the single point t = 1 - delta), not an equivalent restatement, so both
#' are reported.
#'
#' @return list(admissible, tail_mass, worst_slack, worst_t)
admissibility <- function(F0, delta = DELTA, n_grid = 2001L, tol = .REF_TOL) {
  if (!is.function(F0)) stop("F0 must be a function")
  lo <- 1 - delta
  t <- seq(lo, 1, length.out = n_grid)
  Fv <- as.numeric(F0(t))
  if (length(Fv) != length(t)) stop("F0 must return one value per input point")
  slack <- Fv - t
  k <- which.min(slack)
  list(admissible = slack[k] >= -tol,
       tail_mass = 1 - as.numeric(F0(lo))[1L],
       worst_slack = slack[k],
       worst_t = t[k])
}


# --------------------------------------------------------------------------
# q at a declared reference, the family member, and the envelope
# --------------------------------------------------------------------------

#' The interior component q(F0) = 1 - max_i |i/m - F0(t_(i))|. With F0 = NULL
#' (Uniform) this equals prereq_index(x, y)$q bit for bit. Returns 0 when
#' the interior guard fires, as core does.
interior_q <- function(x, y, F0 = NULL, delta = DELTA) {
  ri <- .rescaled_interior(x, y, delta)
  t <- ri$t
  m <- length(t)
  if (.guard_fires(m, ri$n)) return(0)
  Fm <- seq_len(m) / m
  s <- if (is.null(F0)) t else as.numeric(F0(t))
  1 - max(abs(Fm - s))
}

#' The family member PI(F0) = A1 * (q(F0) * ell), composed from the same core
#' components; at Uniform it equals prereq_index(x, y)$PI bit for bit.
prereq_index_family <- function(x, y, F0 = NULL, delta = DELTA,
                                dense_max_n = DENSE_MAX_N) {
  res <- prereq_index(x, y, delta, dense_max_n)
  q <- interior_q(x, y, F0, delta)
  a2 <- q * res$ell
  nm <- if (is.null(F0)) "Uniform(0,1)" else {
    r <- attr(F0, "reference_name"); if (is.null(r)) "F0" else r
  }
  list(PI = res$A1 * a2, A1 = res$A1, A2 = a2, q = q, ell = res$ell,
       reference = nm)
}

#' The exact upper envelope of the coefficient over B.
#'
#' @return list(PI_hi, sup_q, inf_q, D_star, n_tail, m, attained, A1, ell,
#'   q, PI). When the interior guard fires, q is 0 for every reference by
#'   definition, so sup_q = inf_q = PI_hi = 0.
pi_envelope <- function(x, y, delta = DELTA, dense_max_n = DENSE_MAX_N) {
  res <- prereq_index(x, y, delta, dense_max_n)
  ri <- .rescaled_interior(x, y, delta)
  t <- ri$t
  m <- length(t)
  out <- list(A1 = res$A1, ell = res$ell, q = res$q, PI = res$PI, m = m)

  if (.guard_fires(m, ri$n)) {
    return(c(out, list(PI_hi = 0, sup_q = 0, inf_q = 0, D_star = 0,
                       n_tail = 0L, attained = TRUE)))
  }

  Fm <- seq_len(m) / m
  tail <- t >= 1 - delta
  n_tail <- sum(tail)
  d_star <- if (n_tail == 0L) 0 else max(pmax(t[tail] - Fm[tail], 0))
  sup_q <- 1 - d_star

  # Direct evaluation at the attaining reference as a FUNCTION (ties get one
  # common ECDF value), so attainment is measured rather than assumed.
  s_star <- pmax(.ecdf_sorted(t, t), ifelse(tail, t, 0))
  q_star <- 1 - max(abs(Fm - s_star))

  c(out, list(PI_hi = res$A1 * res$ell * sup_q,
              sup_q = sup_q,
              inf_q = 1 / m,
              D_star = d_star,
              n_tail = as.integer(n_tail),
              attained = abs(q_star - sup_q) <= .REF_TOL))
}
