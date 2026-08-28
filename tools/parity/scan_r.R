# scan_r.R -- R side of the scan parity harness.
#
# Reads the same theta.csv and perm_scan.csv the Python driver wrote, runs
# prereq_scan with the shared permutations supplied through the
# indices_provider hook, and writes the result as JSON for the driver to
# compare. JSON is emitted by hand: the package declares no dependencies,
# and this script keeps to the same rule.
#
# Usage: Rscript scan_r.R <dir>

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

if (nzchar(Sys.getenv("R_PRERELATION_LIB"))) .libPaths(Sys.getenv("R_PRERELATION_LIB"))
library(prerelation)

N_PERM <- 239
ALPHA <- 0.05

lines <- readLines(file.path(dir, "theta.csv"), warn = FALSE)
lines <- lines[nzchar(lines)]
attr_names <- strsplit(lines[1L], ",", fixed = TRUE)[[1L]]
body <- lines[-1L]
theta <- matrix(NA_real_, nrow = length(body), ncol = length(attr_names))
for (i in seq_along(body)) {
  v <- suppressWarnings(as.numeric(strsplit(body[i], ",", fixed = TRUE)[[1L]]))
  if (any(is.na(v))) stop(sprintf("theta.csv line %d is not numeric", i + 1L))
  theta[i, ] <- v
}

# The pool is written zero-based by the Python side, exactly like the
# committed golden matrices; read_perm_indices applies the +1.
pool <- read_perm_indices(file.path(dir, "perm_scan.csv"))
cat(sprintf("theta %d x %d, permutation pool %d x %d\n",
            nrow(theta), ncol(theta), nrow(pool), ncol(pool)))

provider <- function(pair_position, n) {
  pool[(pair_position + 1L):(pair_position + N_PERM), , drop = FALSE]
}

res <- prereq_scan(theta, alpha = ALPHA, attr_names = attr_names,
                   indices_provider = provider)

# --- minimal JSON writer -----------------------------------------------
jnum <- function(v) {
  if (length(v) == 0L) return("null")
  sprintf("%.17g", v)
}
jstr <- function(s) paste0("\"", s, "\"")
jarr <- function(items) paste0("[", paste(items, collapse = ","), "]")
jpair <- function(m) {
  if (is.null(m)) return("null")
  if (nrow(m) == 0L) return("[]")
  jarr(apply(m, 1L, function(r) jarr(vapply(r, jstr, ""))))
}

rec_json <- character(nrow(res$records))
for (i in seq_len(nrow(res$records))) {
  r <- res$records[i, ]
  fields <- c(
    paste0(jstr("source"), ":", jstr(r$source)),
    paste0(jstr("target"), ":", jstr(r$target)),
    paste0(jstr("pi"), ":", jnum(r$pi)),
    paste0(jstr("pi_reverse"), ":", jnum(r$pi_reverse)),
    paste0(jstr("delta"), ":", jnum(r$delta)),
    paste0(jstr("A1"), ":", jnum(r$A1)),
    paste0(jstr("A2"), ":", jnum(r$A2)),
    paste0(jstr("q"), ":", jnum(r$q)),
    paste0(jstr("ell"), ":", jnum(r$ell)),
    paste0(jstr("p_value"), ":", jnum(r$p_value)),
    paste0(jstr("p_adj"), ":", jnum(r$p_adj)),
    paste0(jstr("edge"), ":", if (r$edge) "true" else "false")
  )
  rec_json[i] <- paste0("{", paste(fields, collapse = ","), "}")
}

cycles_json <- if (length(res$cycles) == 0L) "[]" else
  jarr(vapply(res$cycles, function(cy) jarr(vapply(cy, jstr, "")), ""))

out <- paste0("{",
  jstr("records"), ":", jarr(rec_json), ",",
  jstr("edges"), ":", jpair(res$edges), ",",
  jstr("cycles"), ":", cycles_json, ",",
  jstr("reduced_edges"), ":", jpair(res$reduced_edges),
  "}")
writeLines(out, file.path(dir, "scan_r_out.json"))

cat(sprintf("edges=%d cycles=%d reduced=%s\n", nrow(res$edges),
            length(res$cycles),
            if (is.null(res$reduced_edges)) "NA" else nrow(res$reduced_edges)))
cat("equivalence classes:",
    paste(vapply(res$equivalence_classes, paste, "", collapse = "+"),
          collapse = " | "), "\n")
