test_that("updateOncologyXenograftSimeoni2004", {
  newModel <- updateOncologyXenograftSimeoni2004(readModelDb("oncology_xenograft_simeoni_2004"), ncmt = 5)
  expect_s3_class(newModel, class = "rxUi")
  expect_true("damagedCells5" %in% newModel$state)
})
