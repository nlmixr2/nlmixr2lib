test_that("updateOncologyXenograftSimeoni2004", {
  suppressMessages(
    origModel <- rxode2::rxode(readModelDb("oncology_xenograft_simeoni_2004"))
  )
  suppressMessages(
    newModel <- updateOncologyXenograftSimeoni2004(origModel, ncmt = 5)
  )
  expect_s3_class(newModel, class = "rxUi")

  # the new states are added
  expect_false("damagedCells5" %in% origModel$state)
  expect_true("damagedCells5" %in% newModel$state)

  # the new tumor line is added
  expect_no_match(
    deparse(as.function(newModel), width.cutoff = 500),
    "tumorVol <- cyclingCells \\+ damagedCells1 \\+ damagedCells2 \\+ damagedCells3 \\+ damagedCells4 \\+ damagedCells5",
    all = FALSE
  )
  expect_match(
    deparse(as.function(newModel), width.cutoff = 500),
    "tumorVol <- cyclingCells \\+ damagedCells1 \\+ damagedCells2 \\+ damagedCells3 \\+ damagedCells4 \\+ damagedCells5",
    all = FALSE
  )
})
