# condense_r.R -- runs condense on the fixed graph set and prints a
# canonical text rendering. Class ids are one-based on the R side.
# Usage: Rscript condense_r.R <graphs.json>
args <- commandArgs(trailingOnly = TRUE)
if (nzchar(Sys.getenv("R_PRERELATION_LIB"))) .libPaths(Sys.getenv("R_PRERELATION_LIB"))
library(prerelation)

# Tiny reader for the fixed-shape graph file (no JSON dependency).
txt <- paste(readLines(args[1], warn = FALSE), collapse = "")
blocks <- regmatches(txt, gregexpr('"[A-Za-z_]+"[[:space:]]*:[[:space:]]*\\{[^}]*\\}', txt))[[1]]
for (b in blocks) {
  name <- sub('^"([A-Za-z_]+)".*$', "\\1", b)
  nodes_s <- sub('.*"nodes"[^\\[]*\\[([^]]*)\\].*', "\\1", b)
  nodes <- gsub('[" ]', "", strsplit(nodes_s, ",", fixed = TRUE)[[1]])
  edges_s <- sub('.*"edges"[^\\[]*\\[(.*)\\][^]]*$', "\\1", b)
  pairs <- regmatches(edges_s, gregexpr('\\[[^]]*\\]', edges_s))[[1]]
  if (length(pairs) == 0) {
    E <- matrix(character(0), nrow = 0, ncol = 2)
  } else {
    E <- do.call(rbind, lapply(pairs, function(p) {
      gsub('[]["  ]', "", strsplit(gsub('[][]', "", p), ",", fixed = TRUE)[[1]])
    }))
  }
  r <- condense(nodes, E)
  cat(sprintf("graph=%s\n", name))
  cat(sprintf("  classes=%s\n", paste(vapply(r$classes, paste, "", collapse = "+"), collapse = "|")))
  cat(sprintf("  class_of=%s\n", paste(sprintf("%s:%d", names(r$class_of), r$class_of), collapse = ",")))
  qe <- if (nrow(r$quotient_edges) == 0) "-" else paste(apply(r$quotient_edges, 1, function(z) sprintf("%d>%d", z[1], z[2])), collapse = ",")
  he <- if (nrow(r$hasse_edges) == 0) "-" else paste(apply(r$hasse_edges, 1, function(z) sprintf("%d>%d", z[1], z[2])), collapse = ",")
  cat(sprintf("  quotient=%s\n  hasse=%s\n", qe, he))
}
