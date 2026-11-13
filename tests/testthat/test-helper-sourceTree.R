test_that(".skipIfNoSourceTree() and .hasSourceTree() pass an existing path through", {
  dir <- withr::local_tempdir()
  file <- file.path(dir, "present.txt")
  writeLines("", file)
  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "true")
  expect_identical(.skipIfNoSourceTree(file), file)
  expect_identical(.skipIfNoSourceTree(dir), dir)
  expect_true(.hasSourceTree(file))
})

test_that("a missing path skips or is FALSE unless the source tree is required", {
  missing <- file.path(withr::local_tempdir(), "absent")
  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "")
  # testthat signals a skip as a condition of class "skip", and prefixes its
  # message with "Reason: ".
  cnd <- tryCatch(.skipIfNoSourceTree(missing, "the thing"), condition = identity)
  expect_s3_class(cnd, "skip")
  expect_identical(sub("^Reason: ", "", conditionMessage(cnd)), "the thing not present (installed package)")
  expect_false(.hasSourceTree(missing))

  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "true")
  expected <- paste0(
    "the thing not found at ",
    missing,
    ", and NLMIXR2LIB_REQUIRE_SOURCE_TREE is set, so this test must not skip"
  )
  err <- expect_error(.skipIfNoSourceTree(missing, "the thing"))
  expect_identical(conditionMessage(err), expected)
  err <- expect_error(.hasSourceTree(missing, "the thing"))
  expect_identical(conditionMessage(err), expected)
})

test_that("the source-tree scanner and file selection find what they claim to", {
  dir <- withr::local_tempdir()
  writeLines(
    c(
      "test_that('gated', {",
      "  root <- .sourceTreeRoot()",
      "  .skipIfNoSourceTree(file.path(root, '_pkgdown.yml'))",
      "})",
      "test_that('conditional', {",
      "  if (.hasSourceTree(.sourceTreeRoot())) NULL",
      "})"
    ),
    file.path(dir, "test-good.R")
  )
  writeLines(
    c(
      "# .skipIfNoSourceTree() in a comment does not select this file",
      "test_that('ungated', {",
      "  readLines(file.path(.sourceTreeRoot(), '_pkgdown.yml'))",
      "})",
      "test_that('upward', {",
      "  testthat::test_path('..', '..', 'R')",
      "})",
      "x <- '../../vignettes'"
    ),
    file.path(dir, "test-bad.R")
  )
  writeLines("f <- function() .skipIfNoSourceTree(.sourceTreeRoot())", file.path(dir, "helper-bad.R"))

  expect_identical(basename(.sourceTreeTestFiles(dir)), "test-good.R")
  expect_identical(
    .ungatedSourceTreeAccess(dir),
    c(
      "helper-bad.R:1: a helper file reads the source tree; do it in the test file",
      "test-bad.R:2: .sourceTreeRoot() without .skipIfNoSourceTree() or .hasSourceTree()",
      "test-bad.R:5: test_path(\"..\", ...) reaches the source tree; use .sourceTreeRoot()",
      "test-bad.R:8: a \"../..\" path reaches the source tree; use .sourceTreeRoot()"
    )
  )
})

test_that("every test reaches the source tree only through a gate", {
  # Runs everywhere: under R CMD check the test files are installed with the
  # tests, so the scan needs no source tree.
  issues <- .ungatedSourceTreeAccess(testthat::test_path())
  expect_identical(issues, character(), info = paste0("\n", paste(issues, collapse = "\n")))
})
