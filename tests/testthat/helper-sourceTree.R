# Some tests read the package SOURCE tree: files that R CMD build leaves out
# (_pkgdown.yml, vignettes/articles/ and .claude/ are in .Rbuildignore) or that
# a test can only reach through the checkout (R/, inst/references/). Under
# R CMD check and covr the tests run against an installed copy, where
# test_path("..", "..") is not the source tree, so these tests skip there.
#
# The source-tree-gates CI job runs every test file that calls
# .skipIfNoSourceTree() from the checkout, with NLMIXR2LIB_REQUIRE_SOURCE_TREE
# set; a missing path is then an error rather than a skip, so the gate cannot
# pass by not running. Use this helper for any new test that reads the source
# tree, so that the job picks it up.

# The package root when the tests run from a source checkout.
.sourceTreeRoot <- function() {
  normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
}

# Skip the calling test when `path` (a file or directory) is absent, unless the
# source tree is required, in which case fail. `what` names the path in the
# message.
.skipIfNoSourceTree <- function(path, what = basename(path)) {
  if (file.exists(path)) {
    return(invisible(path))
  }
  if (nzchar(Sys.getenv("NLMIXR2LIB_REQUIRE_SOURCE_TREE"))) {
    stop(
      what,
      " not found at ",
      path,
      ", and NLMIXR2LIB_REQUIRE_SOURCE_TREE is set, so this test must not skip",
      call. = FALSE
    )
  }
  testthat::skip(paste(what, "not present (installed package)"))
}
