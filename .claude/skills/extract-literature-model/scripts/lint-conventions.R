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

# Double quotes in a trailing `ini()` comment. checkModelConventions() cannot
# see this: it is a property of the FILE TEXT, not of the parsed model object.
# rxode2 promotes a trailing comment on an ini() line that carries no label()
# into `label("<comment>")`, and an embedded double quote terminates that
# generated string early, so the model fails to re-parse and its vignette dies
# at render.
#
# The package already has an enumerating test for this
# (tests/testthat/test-checkModelConventions.R, "no ini() line carries a quoted
# trailing comment rxode2 would promote into a broken label"), but it only runs
# over the whole library once the file is written and the suite is run. Two
# consecutive consolidation merges had to repair the defect after the fact
# (2026-09-05, 15 lines; 2026-09-09, 18 lines), so the same check runs here,
# on the one file the extraction just wrote.
.lintIniQuotes <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(character())
  }
  lines <- readLines(path, warn = FALSE)
  inIni <- FALSE
  bad <- character()
  for (i in seq_along(lines)) {
    ln <- lines[[i]]
    if (grepl("^\\s*ini\\(\\{", ln)) {
      inIni <- TRUE
      next
    }
    if (grepl("^\\s*model\\(\\{", ln)) inIni <- FALSE
    if (!inIni) next
    # A standalone comment is not attached to a parameter; a line that already
    # calls label() is not promoted.
    if (grepl("^\\s*#", ln)) next
    if (grepl("label(", ln, fixed = TRUE)) next
    # First `#` that is not itself inside a string literal.
    chars <- strsplit(ln, "", fixed = TRUE)[[1]]
    nq <- 0L
    hash <- 0L
    for (k in seq_along(chars)) {
      if (chars[[k]] == '"') {
        nq <- nq + 1L
      } else if (chars[[k]] == "#" && nq %% 2L == 0L) {
        hash <- k
        break
      }
    }
    if (hash == 0L) next
    if (grepl('"', substring(ln, hash), fixed = TRUE)) {
      bad <- c(bad, sprintf("%s:%d: %s", basename(path), i, trimws(ln)))
    }
  }
  bad
}

.modelFilePath <- function(model) {
  db <- tryCatch(nlmixr2lib::modeldb, error = function(e) NULL)
  if (is.null(db) || !("filename" %in% names(db))) return(NULL)
  row <- db[db$name == model, , drop = FALSE]
  if (nrow(row) != 1L) return(NULL)
  for (root in c(file.path(getwd(), "inst", "modeldb"),
                 system.file("modeldb", package = "nlmixr2lib"))) {
    p <- file.path(root, row$filename[[1]])
    if (nzchar(root) && file.exists(p)) return(p)
  }
  NULL
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
  quoted <- .lintIniQuotes(.modelFilePath(model))
  if (length(quoted)) {
    cat(sprintf("\n=== %s — %d ini() comment(s) with a double quote ===\n",
                model, length(quoted)))
    for (q in quoted) {
      cat(sprintf("  [error] %s\n", q))
    }
    cat(paste0("        rxode2 promotes a trailing comment on an ini() line ",
               "with no label() into\n        label(\"<comment>\"); an embedded ",
               "double quote terminates that string early\n        and the model ",
               "will not re-parse. Use single quotes.\n\n"))
  }
  if (nrow(res) == 0L && length(quoted) == 0L) {
    cat(sprintf("OK: %s — no convention issues.\n", model))
    return(0L)
  }
  if (nrow(res) == 0L) return(2L)
  cat(sprintf("\n=== %s — %d issue(s) ===\n", model, nrow(res)))
  for (i in seq_len(nrow(res))) {
    cat(formatRow(res[i, ]), "\n\n", sep = "")
  }
  if (length(quoted) || any(res$severity == "error")) {
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
