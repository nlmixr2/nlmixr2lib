test_that("modeldb is correctly generated", {
  expect_true("PK_3cmt" %in% modeldb$name)
})

# Spies installed with local_mocked_bindings(). Each records into an
# environment and delegates to the real function stored there by the test.
.modeldbSpy <- new.env(parent = emptyenv())

.spyCheckModelConventions <- function(model, verbose = TRUE) {
  .modeldbSpy$checkedNames <- c(.modeldbSpy$checkedNames, model$modelName)
  res <- .modeldbSpy$checkModelConventions(model, verbose = verbose)
  .modeldbSpy$issues <- rbind(.modeldbSpy$issues, res)
  res
}

test_that("the build checks a model under its own name, so its specific-scope covariates are allowed", {
  # CTROUGH is registered `Scope: specific`, and Chen_2021_lorlatinib_teae_grade3
  # is one of the models its entry lists. The model used to be checked under
  # the name rxode2 gave it, "eval(parsedFile)", which no entry lists.
  src <- system.file("modeldb", "specificDrugs", "Chen_2021_lorlatinib_teae_grade3.R", package = "nlmixr2lib")
  skip_if_not(nzchar(src), "Chen_2021_lorlatinib_teae_grade3 is not in the library")
  dir <- withr::local_tempdir("modeldb-fixture")
  lines <- readLines(src)
  writeLines(lines, file.path(dir, "Chen_2021_lorlatinib_teae_grade3.R"))
  # The same model under a name the CTROUGH entry does not list.
  writeLines(
    sub("^Chen_2021_lorlatinib_teae_grade3 <- function", "Fixture_2026_ctrough <- function", lines),
    file.path(dir, "Fixture_2026_ctrough.R")
  )
  .modeldbSpy$checkedNames <- character()
  .modeldbSpy$issues <- NULL
  .modeldbSpy$checkModelConventions <- nlmixr2lib::checkModelConventions
  testthat::local_mocked_bindings(checkModelConventions = .spyCheckModelConventions, .package = "nlmixr2lib")

  suppressMessages(nlmixr2lib:::.addDirToModelDbCached(dir, list()))
  expect_identical(.modeldbSpy$checkedNames, c("Chen_2021_lorlatinib_teae_grade3", "Fixture_2026_ctrough"))
  scope <- .modeldbSpy$issues[grepl("scoped 'specific'", .modeldbSpy$issues$message, fixed = TRUE), , drop = FALSE]
  expect_false("Chen_2021_lorlatinib_teae_grade3" %in% scope$model)
  expect_identical(scope$model[scope$name == "CTROUGH"], "Fixture_2026_ctrough")
})
