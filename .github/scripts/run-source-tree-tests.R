#!/usr/bin/env Rscript
#
# Runs, from the checkout, every test file that reads the package source tree,
# and fails on any failure, error or skip.
#
# WHY
#
# Tests that read _pkgdown.yml, vignettes/articles/, .claude/, R/ or
# inst/references/ through the checkout skip under R CMD check and covr, which
# run the tests against an installed copy. Those tests gate the naming
# registers, the article links and the locale-stable navbar, and while they
# skipped in every CI job they went red only when someone ran them by hand.
#
# WHICH FILES
#
# Every tests/testthat/test-*.R that calls .skipIfNoSourceTree() or
# .hasSourceTree(), the gates in tests/testthat/helper-sourceTree.R; the
# selection reads parsed calls, so a mention in a comment does not count. With
# NLMIXR2LIB_REQUIRE_SOURCE_TREE set, a gate errors instead of skipping, and
# this script also fails on any other skip (a missing Suggests package, say),
# so the job cannot pass by not running a test. A test in that directory
# checks that nothing reaches the source tree except through a gate.
#
# Run from the package root: Rscript .github/scripts/run-source-tree-tests.R

testDir <- file.path("tests", "testthat")
helpers <- new.env()
sys.source(file.path(testDir, "helper-collate.R"), envir = helpers)
sys.source(file.path(testDir, "helper-sourceTree.R"), envir = helpers)
files <- basename(helpers$.sourceTreeTestFiles(testDir))
if (length(files) == 0L) {
  stop("no test file calls a source-tree gate; run this from the package root", call. = FALSE)
}
cat("Source-tree test files:\n", paste0("  ", files, "\n"), sep = "")

Sys.setenv(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "true")
stems <- sub("^test-(.*)[.]R$", "\\1", files)
filter <- paste0("^(", paste(gsub(".", "[.]", stems, fixed = TRUE), collapse = "|"), ")$")
res <- as.data.frame(testthat::test_local(".", filter = filter, stop_on_failure = FALSE))

bad <- res[res$failed > 0L | res$error | res$skipped, c("file", "test", "failed", "error", "skipped")]
ran <- unique(res$file)
cat(sprintf(
  "\n%d tests in %d of %d files: %d failed, %d errored, %d skipped\n",
  nrow(res),
  length(ran),
  length(files),
  sum(res$failed > 0L),
  sum(res$error),
  sum(res$skipped)
))
if (!setequal(ran, files)) {
  cat("Files that reported no tests:", setdiff(files, ran), "\n")
  quit(status = 1L)
}
if (nrow(bad) > 0L) {
  print(bad, row.names = FALSE)
  quit(status = 1L)
}
