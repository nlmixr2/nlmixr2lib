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

test_that("simMarkovId", {
  d <-
    data.frame(
      prAtoA = c(0.2, 0.1),
      prAtoB = c(0.8, 0.9),
      prBtoA = c(0.2, 0.1),
      prBtoB = c(0.8, 0.9)
    )
  probCols <-
    list(
      A = c(A = "prAtoA", B = "prAtoB"),
      B = c(A = "prBtoA", B = "prBtoB")
    )

  expect_equal(
    withr::with_seed(
      5,
      simMarkovId(data = d, initialState = "A", prCols = probCols)
    ),
    data.frame(prev = c("A", "B"), cur = c("B", "B"))
  )

  # A probability column that is not in the data is an error
  probCols <-
    list(
      A = c(A = "prAtoX", B = "prAtoB"),
      B = c(A = "prBtoA", B = "prBtoB")
    )
  expect_error(
    simMarkovId(data = d, initialState = "A", prCols = probCols),
    regexp = "Names must include the elements.*but is missing elements.*prAtoX"
  )
  # A name in the probability column list that is not within the outer list is
  # an error
  probCols <-
    list(
      A = c(A = "prAtoA", B = "prAtoB"),
      B = c(X = "prBtoA", B = "prBtoB")
    )
  expect_error(
    simMarkovId(data = d, initialState = "A", prCols = probCols),
    regexp = "Names must be a subset of.*but has additional elements.*'X'"
  )
})