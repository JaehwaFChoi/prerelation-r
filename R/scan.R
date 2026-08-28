# prerelation/scan.R -- pairwise screening of a whole attribute set.
#
# `scan` takes an (n_persons x n_attributes) matrix of trait values on a
# common anchored scale and returns
#
#   * a tidy record per ordered pair with Pi, the reverse Pi, Delta and a
#     permutation p-value;
#   * the edge set surviving Benjamini-Hochberg control of the false
#     discovery rate;
#   * a cycle report;
#   * the transitive reduction of the edge set when it is acyclic
#     (Aho, Garey and Ullman, 1972: for a directed acyclic graph the
#     transitive reduction is unique and is a subgraph of the original);
#   * the equivalence-class condensation with its quotient order (see
#     `condense`), which is always defined even when the raw edge set has
#     cycles -- the Hasse diagram of the scan is drawn on the quotient.
#
# What the scan recovers is a *dominance preorder* over the attributes --
# the ordering induced by which attributes act as ceilings on which others
# -- not a direct-prerequisite DAG. Indirect dominance produces edges of
# its own (removed only along chains by the transitive reduction), and
# siblings that share a common ceiling can be linked to each other even
# though neither is a prerequisite for the other. A disagreement between
# the recovered order and an expert-specified prerequisite graph is
# therefore a difference between two concepts, not by itself an error in
# either. Mutually dominating attributes (directed cycles) are expected
# behavior of a pairwise index and are read as equivalence classes of the
# preorder.
#
# Provenance: `bh_fdr`, `find_cycles`, `transitive_reduction` and `scan`
# are transcribed from prerelation/scan.py (Python, the reference
# implementation); `condense` is transcribed from prerelation-js
# src/scan.mjs, which is where it lives -- the Python package does not
# have it.


.as_edges <- function(edges) {
  if (is.null(edges) || length(edges) == 0L) {
    return(matrix(character(0), nrow = 0L, ncol = 2L))
  }
  if (is.list(edges) && !is.data.frame(edges)) {
    E <- do.call(rbind, lapply(edges, function(e) as.character(e[1:2])))
    return(matrix(as.character(E), ncol = 2L))
  }
  E <- as.matrix(edges)
  if (ncol(E) != 2L) stop("edges must have exactly two columns (source, target)")
  matrix(as.character(E), ncol = 2L)
}


.adjacency <- function(nodes, E, unique_successors = FALSE) {
  adj <- .set_names(rep(list(character(0)), length(nodes)), nodes)
  if (nrow(E) > 0L) {
    for (i in seq_len(nrow(E))) {
      u <- E[i, 1L]
      if (!(u %in% nodes)) stop(sprintf("edge source '%s' is not in nodes", u))
      if (!(E[i, 2L] %in% nodes)) {
        stop(sprintf("edge target '%s' is not in nodes", E[i, 2L]))
      }
      adj[[u]] <- c(adj[[u]], E[i, 2L])
    }
  }
  if (unique_successors) adj <- lapply(adj, unique)
  adj
}


#' Benjamini-Hochberg adjusted p-values (step-up, monotonised).
bh_fdr <- function(pvalues) {
  p <- as.double(pvalues)
  m <- length(p)
  if (m == 0L) return(p)
  ord <- order(p, method = "radix")   # stable, as NumPy's mergesort argsort
  ranked <- p[ord] * m / seq_len(m)
  ranked <- rev(cummin(rev(ranked)))  # monotonise from the largest p downward
  out <- numeric(m)
  out[ord] <- pmin(pmax(ranked, 0), 1)
  out
}


#' Return the directed cycles found by depth-first search (may be empty).
find_cycles <- function(nodes, edges) {
  nodes <- as.character(nodes)
  adj <- .adjacency(nodes, .as_edges(edges))
  colour <- .set_names(rep(0L, length(nodes)), nodes)  # 0 white 1 grey 2 black
  stack <- character(0)
  cycles <- list()

  visit <- function(u) {
    colour[[u]] <<- 1L
    stack <<- c(stack, u)
    for (v in adj[[u]]) {
      if (colour[[v]] == 0L) {
        visit(v)
      } else if (colour[[v]] == 1L) {
        pos <- match(v, stack)
        cycles[[length(cycles) + 1L]] <<- c(stack[pos:length(stack)], v)
      }
    }
    stack <<- stack[-length(stack)]
    colour[[u]] <<- 2L
  }

  for (u in nodes) if (colour[[u]] == 0L) visit(u)
  cycles
}


# Is dst reachable from src without traversing the single edge (skip_u, skip_v)?
.reachable <- function(adj, src, dst, skip_u, skip_v) {
  stack <- src
  seen <- src
  while (length(stack) > 0L) {
    u <- stack[length(stack)]
    stack <- stack[-length(stack)]
    for (v in adj[[u]]) {
      if (u == skip_u && v == skip_v) next
      if (v == dst) return(TRUE)
      if (!(v %in% seen)) {
        seen <- c(seen, v)
        stack <- c(stack, v)
      }
    }
  }
  FALSE
}


#' Transitive reduction of a directed acyclic graph.
#'
#' An edge is dropped when the same ordering is already implied by a path of
#' length two or more. Errors if the graph has a cycle, where the reduction
#' is not unique.
transitive_reduction <- function(nodes, edges) {
  nodes <- as.character(nodes)
  E <- .as_edges(edges)
  cycles <- find_cycles(nodes, E)
  if (length(cycles) > 0L) {
    stop(sprintf(paste0("transitive reduction is only defined for acyclic ",
                        "graphs; found cycle %s"),
                 paste(cycles[[1L]], collapse = " -> ")))
  }
  adj <- .adjacency(nodes, E, unique_successors = TRUE)
  keep <- logical(nrow(E))
  if (nrow(E) > 0L) {
    for (i in seq_len(nrow(E))) {
      keep[i] <- !.reachable(adj, E[i, 1L], E[i, 2L], E[i, 1L], E[i, 2L])
    }
  }
  E[keep, , drop = FALSE]
}


#' Condensation of a directed graph into its equivalence classes.
#'
#' Tarjan's strongly connected components, then the quotient graph over the
#' classes, then the transitive reduction of that quotient. Mutually
#' dominating attributes collapse into one class, so the condensation is
#' defined even when the raw edge set has cycles.
#'
#' @return list with elements
#'   classes       list of character vectors, one per equivalence class;
#'   class_of      named integer vector, node -> class id;
#'   quotient_edges  two-column integer matrix of class ids;
#'   hasse_edges     transitive reduction of quotient_edges.
#'
#' Class ids are ONE-BASED here, following R. The JavaScript reference
#' emits the same classes in the same order with zero-based ids; a parity
#' check must add one to the JavaScript ids (the same off-by-one that the
#' shared permutation matrices carry).
condense <- function(nodes, edges) {
  nodes <- as.character(nodes)
  E <- .as_edges(edges)
  adj <- .adjacency(nodes, E)

  idx <- .set_names(rep(NA_integer_, length(nodes)), nodes)
  low <- .set_names(rep(NA_integer_, length(nodes)), nodes)
  on_stack <- .set_names(rep(FALSE, length(nodes)), nodes)
  tstack <- character(0)
  components <- list()
  counter <- 0L

  for (start in nodes) {
    if (!is.na(idx[[start]])) next
    work_node <- start
    work_ptr <- 0L
    while (length(work_node) > 0L) {
      top <- length(work_node)
      u <- work_node[top]
      if (work_ptr[top] == 0L) {
        idx[[u]] <- counter
        low[[u]] <- counter
        counter <- counter + 1L
        tstack <- c(tstack, u)
        on_stack[[u]] <- TRUE
      }
      nb <- adj[[u]]
      advanced <- FALSE
      while (work_ptr[top] < length(nb)) {
        v <- nb[work_ptr[top] + 1L]
        work_ptr[top] <- work_ptr[top] + 1L
        if (is.na(idx[[v]])) {
          work_node <- c(work_node, v)
          work_ptr <- c(work_ptr, 0L)
          advanced <- TRUE
          break
        } else if (on_stack[[v]]) {
          low[[u]] <- min(low[[u]], idx[[v]])
        }
      }
      if (advanced) next
      if (low[[u]] == idx[[u]]) {
        comp <- character(0)
        repeat {
          w <- tstack[length(tstack)]
          tstack <- tstack[-length(tstack)]
          on_stack[[w]] <- FALSE
          comp <- c(comp, w)
          if (identical(w, u)) break
        }
        components[[length(components) + 1L]] <- comp
      }
      work_node <- work_node[-top]
      work_ptr <- work_ptr[-top]
      if (length(work_node) > 0L) {
        parent <- work_node[length(work_node)]
        low[[parent]] <- min(low[[parent]], low[[u]])
      }
    }
  }

  # Members within a class, and the classes themselves, are ordered by the
  # position of the node in `nodes`.
  classes <- lapply(components, function(comp) comp[order(match(comp, nodes))])
  classes <- classes[order(vapply(classes, function(cl) match(cl[1L], nodes), 1L))]

  class_of <- .set_names(rep(NA_integer_, length(nodes)), nodes)
  for (ci in seq_along(classes)) {
    for (u in classes[[ci]]) class_of[[u]] <- ci
  }

  quotient_edges <- matrix(integer(0), nrow = 0L, ncol = 2L)
  if (nrow(E) > 0L) {
    seen <- character(0)
    for (i in seq_len(nrow(E))) {
      cu <- class_of[[E[i, 1L]]]
      cv <- class_of[[E[i, 2L]]]
      if (cu == cv) next
      key <- paste0(cu, "->", cv)
      if (!(key %in% seen)) {
        seen <- c(seen, key)
        quotient_edges <- rbind(quotient_edges, c(cu, cv))
      }
    }
  }

  class_ids <- as.character(seq_along(classes))
  hasse_chr <- transitive_reduction(class_ids, matrix(as.character(quotient_edges),
                                                     ncol = 2L))
  hasse_edges <- matrix(as.integer(hasse_chr), ncol = 2L)

  list(classes = classes, class_of = class_of,
       quotient_edges = quotient_edges, hasse_edges = hasse_edges)
}


#' Screen every ordered pair of attributes for ceiling dominance.
#'
#' The edge set (and its transitive reduction) is read as a *dominance
#' preorder* over the attributes, not as a recovered direct-prerequisite
#' DAG.
#'
#' Design floor on permutation replicates. With k attributes there are
#' K = k (k - 1) ordered pairs, and the smallest attainable permutation
#' p-value is 1 / (n_perm + 1). For any pair to survive Benjamini-Hochberg
#' control at level alpha the replicate count must satisfy
#' n_perm >= K / alpha - 1 (K = 6, alpha = 0.05 requires n_perm >= 119;
#' K = 56 requires n_perm >= 1119). Below the floor the scan cannot return
#' any edge, regardless of the data.
#'
#' @param theta_matrix numeric matrix, persons x attributes, on a common
#'   anchored scale.
#' @param alpha target false discovery rate for the BH step-up procedure,
#'   applied jointly to all ordered pairs.
#' @param attr_names attribute labels; defaults to A1, A2, ... (the Python
#'   and JavaScript implementations spell this argument `names`, which is
#'   avoided here because it shadows base::names).
#' @param n_perm permutation replicates per ordered pair, used only when
#'   indices_provider is NULL. R's generator does not reproduce the NumPy
#'   or JavaScript streams, so a run without indices_provider is
#'   reproducible within R but is NOT comparable across implementations.
#' @param indices_provider parity hook: a function(pair_position, n)
#'   returning a matrix of ONE-BASED permutation indices for the ordered
#'   pair at pair_position (zero-based, row-major, diagonal skipped).
#'   n_perm is then that matrix's row count.
#' @param min_pi additional floor on Pi for an edge to be kept.
#' @param require_positive_delta keep an edge only when the forward
#'   direction dominates the reverse.
prereq_scan <- function(theta_matrix, alpha = 0.05, attr_names = NULL, n_perm = 999,
                 seed = 0, delta = DELTA, min_pi = 0,
                 require_positive_delta = TRUE, indices_provider = NULL,
                 dense_max_n = DENSE_MAX_N) {
  theta <- as.matrix(theta_matrix)
  storage.mode(theta) <- "double"
  if (length(dim(theta)) != 2L) {
    stop("theta_matrix must be two-dimensional (persons x attributes)")
  }
  n <- nrow(theta)
  k <- ncol(theta)
  if (k < 2L) stop("need at least two attributes to scan")
  if (is.null(attr_names)) attr_names <- paste0("A", seq_len(k))
  attr_names <- as.character(attr_names)
  if (length(attr_names) != k) stop("attr_names must have one entry per attribute")

  pi_matrix <- matrix(NA_real_, k, k)
  comp <- vector("list", k * k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i != j) {
        res <- prereq_index(theta[, i], theta[, j], delta, dense_max_n)
        pi_matrix[i, j] <- res$PI
        comp[[(i - 1L) * k + j]] <- res
      }
    }
  }

  rows <- list()
  pair_position <- 0L
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i == j) next
      x <- theta[, i]
      y <- theta[, j]
      res <- comp[[(i - 1L) * k + j]]
      obs <- res$PI

      if (!is.null(indices_provider)) {
        P <- .check_perm_indices(indices_provider(pair_position, n), n)
        this_n_perm <- nrow(P)
        cnt <- 0L
        for (r in seq_len(this_n_perm)) {
          if (prereq_index(x, y[P[r, ]], delta, dense_max_n)$PI >= obs) {
            cnt <- cnt + 1L
          }
        }
      } else {
        this_n_perm <- n_perm
        set.seed(seed + pair_position)
        cnt <- 0L
        for (r in seq_len(this_n_perm)) {
          if (prereq_index(x, y[sample.int(n)], delta, dense_max_n)$PI >= obs) {
            cnt <- cnt + 1L
          }
        }
      }
      p_value <- (cnt + 1) / (this_n_perm + 1)

      rows[[length(rows) + 1L]] <- data.frame(
        source = attr_names[i], target = attr_names[j],
        pi = obs, pi_reverse = pi_matrix[j, i],
        delta = obs - pi_matrix[j, i],
        A1 = res$A1, A2 = res$A2, q = res$q, ell = res$ell,
        p_value = p_value, n = n, n_perm = this_n_perm,
        stringsAsFactors = FALSE)
      pair_position <- pair_position + 1L
    }
  }

  records <- do.call(rbind, rows)
  records$p_adj <- bh_fdr(records$p_value)
  keep <- records$p_adj <= alpha & records$pi >= min_pi
  if (require_positive_delta) keep <- keep & records$delta > 0
  records$edge <- keep

  edges <- as.matrix(records[keep, c("source", "target"), drop = FALSE])
  dimnames(edges) <- NULL
  edges <- matrix(as.character(edges), ncol = 2L)

  cycles <- find_cycles(attr_names, edges)
  reduced <- if (length(cycles) > 0L) NULL else transitive_reduction(attr_names, edges)
  quotient <- condense(attr_names, edges)

  list(records = records, edges = edges, reduced_edges = reduced,
       cycles = cycles, equivalence_classes = quotient$classes,
       quotient = quotient, names = attr_names, alpha = alpha,
       meta = list(n = n, n_perm = n_perm, seed = seed, delta = delta))
}
