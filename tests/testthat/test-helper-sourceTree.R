test_that(".skipIfNoSourceTree() passes an existing path through", {
  dir <- withr::local_tempdir()
  file <- file.path(dir, "present.txt")
  writeLines("", file)
  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "true")
  expect_identical(.skipIfNoSourceTree(file), file)
  expect_identical(.skipIfNoSourceTree(dir), dir)
})

test_that(".skipIfNoSourceTree() skips on a missing path unless the source tree is required", {
  missing <- file.path(withr::local_tempdir(), "absent")
  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "")
  # testthat signals a skip as a condition of class "skip", and prefixes its
  # message with "Reason: ".
  cnd <- tryCatch(.skipIfNoSourceTree(missing, "the thing"), condition = identity)
  expect_s3_class(cnd, "skip")
  expect_identical(sub("^Reason: ", "", conditionMessage(cnd)), "the thing not present (installed package)")

  withr::local_envvar(NLMIXR2LIB_REQUIRE_SOURCE_TREE = "true")
  err <- expect_error(.skipIfNoSourceTree(missing, "the thing"))
  expect_identical(
    conditionMessage(err),
    paste0(
      "the thing not found at ",
      missing,
      ", and NLMIXR2LIB_REQUIRE_SOURCE_TREE is set, so this test must not skip"
    )
  )
})
