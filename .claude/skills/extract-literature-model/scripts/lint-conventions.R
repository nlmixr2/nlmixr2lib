#!/usr/bin/env Rscript

# Pretty-print nlmixr2lib::checkModelConventions() output so the result
# is paste-ready for a PR body. Used by the extract-literature-model
# skill at Phase 6 step 3.
#
# Usage:
#   Rscript lint-conventions.R <model-name>
#   Rscript lint-conventions.R --all                (audit all registered models)
#
# Exit codes:
#   0  clean (no warnings or errors)
#   1  warnings only
#   2  errors present
#   3  argument / environment problem

# Load the worktree's nlmixr2lib package when run from a package directory
# (DESCRIPTION present at cwd or its first parent) so that lint reflects the
# in-development source, not a possibly-stale installed copy. Fall back to
# the installed package when not in a development tree.
.loadNlmixr2lib <- function() {
  cwd <- getwd()
  desc <- NULL
  for (d in c(cwd, dirname(cwd), dirname(dirname(cwd)))) {
    p <- file.path(d, "DESCRIPTION")
    if (file.exists(p)) {
      desc <- p
      break
    }
  }
  is_pkg <- !is.null(desc) && grepl(
    "^Package:\\s*nlmixr2lib",
    readLines(desc, n = 1L, warn = FALSE)
  )
  if (is_pkg && requireNamespace("pkgload", quietly = TRUE)) {
    suppressPackageStartupMessages(
      pkgload::load_all(dirname(desc), quiet = TRUE)
    )
  } else if (is_pkg && requireNamespace("devtools", quietly = TRUE)) {
    suppressPackageStartupMessages(
      devtools::load_all(dirname(desc), quiet = TRUE)
    )
  } else {
    suppressPackageStartupMessages(
      library(nlmixr2lib, warn.conflicts = FALSE)
    )
  }
  invisible(NULL)
}
.loadNlmixr2lib()

formatRow <- function(row) {
  prefix <- switch(
    row$severity,
    error = "[ERROR]",
    warning = "[WARN] ",
    info = "[info] ",
    paste0("[", row$severity, "]")
  )
  bits <- c(
    sprintf("%s %s: %s", prefix, row$category, row$name),
    sprintf("        %s", row$message)
  )
  if (!is.na(row$suggestion) && nzchar(row$suggestion)) {
    bits <- c(bits, sprintf("        suggestion: %s", row$suggestion))
  }
  paste(bits, collapse = "\n")
}

lintOne <- function(model) {
  res <- tryCatch(
    suppressWarnings(checkModelConventions(model, verbose = FALSE)),
    error = function(e) {
      message(sprintf("[ERROR] %s could not be parsed: %s",
                      model, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(res)) return(2L)
  if (nrow(res) == 0L) {
    cat(sprintf("OK: %s — no convention issues.\n", model))
    return(0L)
  }
  cat(sprintf("\n=== %s — %d issue(s) ===\n", model, nrow(res)))
  for (i in seq_len(nrow(res))) {
    cat(formatRow(res[i, ]), "\n\n", sep = "")
  }
  if (any(res$severity == "error")) {
    2L
  } else if (any(res$severity == "warning")) {
    1L
  } else {
    0L
  }
}

main <- function(argv) {
  if (length(argv) == 0L) {
    message("Usage: lint-conventions.R <model-name> | --all")
    quit(status = 3L)
  }
  if (identical(argv[1], "--all")) {
    models <- nlmixr2lib::modellib()$name
    worst <- 0L
    for (m in models) {
      worst <- max(worst, lintOne(m))
    }
    quit(status = worst)
  }
  worst <- 0L
  for (m in argv) {
    worst <- max(worst, lintOne(m))
  }
  quit(status = worst)
}

if (!interactive()) {
  argv <- commandArgs(trailingOnly = TRUE)
  main(argv)
}
