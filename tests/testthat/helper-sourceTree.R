# Some tests read the package SOURCE tree: files that R CMD build leaves out
# (_pkgdown.yml, vignettes/articles/ and .claude/ are in .Rbuildignore) or that
# a test can only reach through the checkout (R/, inst/references/). Under
# R CMD check and covr the tests run against an installed copy, where
# test_path("..", "..") is not the source tree, so these tests skip there.
# (Installed paths such as system.file("modeldb") are not source-tree reads:
# they exist in the installed package and those tests run everywhere.)
#
# The source-tree-gates CI job runs every test file that calls one of the
# gates below, from the checkout, with NLMIXR2LIB_REQUIRE_SOURCE_TREE set; a
# missing path is then an error rather than a skip, so the gate cannot pass by
# not running. A test reaches the source tree only through .sourceTreeRoot(),
# behind one of the gates, in a test file (not a helper file), so the job
# selects every such test; .ungatedSourceTreeAccess() enforces that.
#
# The call walkers .collateCallNodes() and .collateCalleeName() are in
# helper-collate.R.

.sourceTreeGates <- c(".skipIfNoSourceTree", ".hasSourceTree")

# The package root when the tests run from a source checkout.
.sourceTreeRoot <- function() {
  normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
}

# TRUE when `path` (a file or directory) exists. FALSE when it does not,
# unless the source tree is required, in which case an error. For a test that
# does more when the source tree is present but also runs without it. `what`
# names the path in the message.
.hasSourceTree <- function(path, what = basename(path)) {
  if (file.exists(path)) {
    return(TRUE)
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
  FALSE
}

# Skip the calling test when `path` is absent, unless the source tree is
# required, in which case fail.
.skipIfNoSourceTree <- function(path, what = basename(path)) {
  if (!.hasSourceTree(path, what)) {
    testthat::skip(paste(what, "not present (installed package)"))
  }
  invisible(path)
}

# Names of every function called in the R file at `path`, from its parse tree,
# so a mention in a comment or a string does not count.
.calleesIn <- function(path) {
  calls <- .collateCallNodes(parse(path, keep.source = FALSE))
  unname(vapply(calls, .collateCalleeName, character(1)))
}

# The test files in `dir` that read the source tree: those calling a gate.
.sourceTreeTestFiles <- function(dir) {
  files <- list.files(dir, pattern = "^test-.*[.]R$", full.names = TRUE)
  files[vapply(files, function(f) any(.sourceTreeGates %in% .calleesIn(f)), logical(1))]
}

# One message per way a file in `dir` reaches the source tree other than
# through .sourceTreeRoot() behind a gate in a test file: test_path() going up
# from tests/testthat, a "../.." path string, .sourceTreeRoot() in a top-level
# expression (a test_that() block, say) that calls no gate, and any helper file
# touching the source tree, whose callers the CI job could not select.
#
# The two files that define and test the scanner are not scanned: their
# fixtures and messages are made of exactly the patterns it looks for.
.ungatedSourceTreeAccess <- function(dir) {
  out <- character()
  for (f in list.files(dir, pattern = "[.]R$", full.names = TRUE)) {
    file <- basename(f)
    if (file %in% c("helper-sourceTree.R", "test-helper-sourceTree.R")) {
      next
    }
    exprs <- parse(f, keep.source = TRUE)
    starts <- utils::getSrcLocation(exprs, "line")
    for (i in seq_along(exprs)) {
      where <- paste0(file, ":", starts[[i]])
      calls <- .collateCallNodes(exprs[[i]])
      callee <- unname(vapply(calls, .collateCalleeName, character(1)))
      upward <- vapply(
        calls[callee == "test_path"],
        function(cl) length(cl) >= 2L && identical(cl[[2]], ".."),
        logical(1)
      )
      if (any(upward)) {
        out <- c(out, paste0(where, ": test_path(\"..\", ...) reaches the source tree; use .sourceTreeRoot()"))
      }
      if (startsWith(file, "helper-") && any(c(".sourceTreeRoot", .sourceTreeGates) %in% callee)) {
        out <- c(out, paste0(where, ": a helper file reads the source tree; do it in the test file"))
      } else if (".sourceTreeRoot" %in% callee && !any(.sourceTreeGates %in% callee)) {
        out <- c(out, paste0(where, ": .sourceTreeRoot() without .skipIfNoSourceTree() or .hasSourceTree()"))
      }
    }
    pd <- utils::getParseData(exprs)
    for (j in which(pd$token == "STR_CONST" & grepl("../..", pd$text, fixed = TRUE))) {
      out <- c(
        out,
        paste0(file, ":", pd$line1[[j]], ": a \"../..\" path reaches the source tree; use .sourceTreeRoot()")
      )
    }
  }
  out
}
