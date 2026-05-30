test_that("updateOncologyXenograftSimeoni2004", {
  suppressMessages(
    origModel <- rxode2::rxode(readModelDb("oncology_xenograft_simeoni_2004"))
  )
  suppressMessages(
    newModel <- updateOncologyXenograftSimeoni2004(origModel, ncmt = 5)
  )
  expect_s3_class(newModel, class = "rxUi")

  # the new states are added
  expect_false("damaged_cells5" %in% origModel$state)
  expect_true("damaged_cells5" %in% newModel$state)

  # the new tumor line is added
  expect_no_match(
    deparse(as.function(newModel), width.cutoff = 500),
    "tumor_vol <- cycling_cells \\+ damaged_cells1 \\+ damaged_cells2 \\+ damaged_cells3 \\+ damaged_cells4 \\+ damaged_cells5", # nolint: line_length_linter.
    all = FALSE
  )
  expect_match(
    deparse(as.function(newModel), width.cutoff = 500),
    "tumor_vol <- cycling_cells \\+ damaged_cells1 \\+ damaged_cells2 \\+ damaged_cells3 \\+ damaged_cells4 \\+ damaged_cells5", # nolint: line_length_linter.
    all = FALSE
  )
})
