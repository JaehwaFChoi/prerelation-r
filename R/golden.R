# prerelation/golden.R -- readers for the shared parity data.
#
# The `prerelation` Python package carries tests/golden: six fixtures, the
# committed permutation index matrices and expected.json. Those files are
# the specification this implementation is checked against, so reading them
# is part of the package rather than of a throwaway script.
#
# Both readers are deliberately strict about type. A CSV writer that emits
# a wrapped float literal (for example a repr-style writer producing
# "np.float64(0.123)") makes read.csv return a character column, and the
# failure then surfaces far downstream as "non-numeric argument to binary
# operator". Here the conversion is explicit and a non-numeric field stops
# the read at the point where the bad byte was found.

.split_numeric_line <- function(line, path, lineno) {
  parts <- strsplit(line, ",", fixed = TRUE)[[1L]]
  vals <- suppressWarnings(as.numeric(parts))
  if (any(is.na(vals))) {
    bad <- parts[is.na(vals)][1L]
    stop(sprintf("%s line %d: field '%s' is not numeric", path, lineno, bad))
  }
  vals
}


#' Read a golden fixture CSV (header "x,y", one pair per row).
#'
#' @return list with numeric elements x and y.
read_fixture <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  header <- strsplit(lines[1L], ",", fixed = TRUE)[[1L]]
  if (!identical(header, c("x", "y"))) {
    stop(sprintf("%s: expected header 'x,y', got '%s'", path, lines[1L]))
  }
  body <- lines[-1L]
  m <- length(body)
  x <- numeric(m)
  y <- numeric(m)
  for (i in seq_len(m)) {
    vals <- .split_numeric_line(body[i], path, i + 1L)
    if (length(vals) != 2L) {
      stop(sprintf("%s line %d: expected 2 fields, got %d", path, i + 1L,
                   length(vals)))
    }
    x[i] <- vals[1L]
    y[i] <- vals[2L]
  }
  list(x = x, y = y)
}


#' Read a committed permutation index matrix.
#'
#' The committed matrices hold ZERO-BASED indices, because they were written
#' from Python. R is one-indexed, so this reader adds one; everything
#' downstream (perm_pvalue, prereq_scan's indices_provider) expects
#' one-based indices. This is the single place the offset is applied.
#'
#' @param path CSV file, one permutation per row, no header.
#' @param zero_based whether the file holds zero-based indices (TRUE for
#'   every matrix committed in the Python package).
#' @return integer matrix, n_perm x n, one-based.
read_perm_indices <- function(path, zero_based = TRUE) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  rows <- lapply(seq_along(lines), function(i) {
    v <- .split_numeric_line(lines[i], path, i)
    if (any(v != round(v))) {
      stop(sprintf("%s line %d: permutation indices must be integers", path, i))
    }
    as.integer(v)
  })
  widths <- vapply(rows, length, 1L)
  if (length(unique(widths)) != 1L) {
    stop(sprintf("%s: rows have differing lengths (%s)", path,
                 paste(unique(widths), collapse = ", ")))
  }
  P <- matrix(unlist(rows, use.names = FALSE), nrow = length(rows),
              ncol = widths[1L], byrow = TRUE)
  if (zero_based) P <- P + 1L
  n <- ncol(P)
  for (i in seq_len(nrow(P))) {
    if (!identical(sort(P[i, ]), seq_len(n))) {
      stop(sprintf("%s line %d: row is not a permutation of 1..%d", path, i, n))
    }
  }
  P
}
