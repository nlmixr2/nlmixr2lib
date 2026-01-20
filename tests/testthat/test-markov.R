test_that("createMarkovModelDataset.default", {
  expect_equal(
    createMarkovModelDataset.default(colPrev = c(1, 2, 1), colCur = c(1, 2, 3)),
    data.frame(
      prevX1 = c(TRUE, FALSE, TRUE),
      curX1 = c(TRUE, FALSE, FALSE),
      prevX2 = c(FALSE, TRUE, FALSE),
      curX2 = c(FALSE, TRUE, FALSE),
      prevX3 = c(FALSE, FALSE, FALSE),
      curX3 = c(FALSE, FALSE, TRUE)
    )
  )
})
