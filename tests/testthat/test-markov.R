test_that("createMarkovTransitionMatrix", {
  expect_message(
    transition <- createMarkovTransitionMatrix(colPrev = c(1, 2, 1), colCur = c(2, 1, 1)),
    regexp = "The following appear only to be collecting state(s): 2",
    fixed = TRUE
  )
  expect_equal(
    transition,
    matrix(c(0.5, 0.5, 1, 0), nrow = 2, byrow = TRUE, dimnames = list(c("1", "2"), c("1", "2")))
  )
  expect_error(
    createMarkovTransitionMatrix(colPrev = c(1, 2, NA), colCur = c(2, 1, 1)),
    regexp = "`colPrev` cannot be `NA`"
  )
  expect_error(
    createMarkovTransitionMatrix(colPrev = c(1, 2, 1), colCur = c(2, 1, NA)),
    regexp = "`colCur` cannot be `NA`"
  )
  expect_error(
    createMarkovTransitionMatrix(colPrev = c(1, 1, 1), colCur = c(1, 1, 1)),
    regexp = "Only one state detected, cannot create a nontrivial Markov model"
  )
  expect_message(
    createMarkovTransitionMatrix(colPrev = c(1, 1, 2, 3, 2, 3), colCur = c(2, 3, 3, 2, 2, 3)),
    regexp = "The following state(s) appear only to be initial states: 1",
    fixed = TRUE
  )
  expect_message(
    createMarkovTransitionMatrix(colPrev = c(1, 2, 3, 1, 2), colCur = c(3, 3, 3, 2, 1)),
    regexp = "The following appear only to be collecting state(s): 3",
    fixed = TRUE
  )
})

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

test_that("createMarkovModelDataset.factor", {
  expect_equal(
    createMarkovModelDataset(factor(c("A", "B", "A")), factor(c("B", "A", "B"))),
    data.frame(
      prevA = c(TRUE, FALSE, TRUE),
      curA = c(FALSE, TRUE, FALSE),
      prevB = c(FALSE, TRUE, FALSE),
      curB = c(TRUE, FALSE, TRUE)
    )
  )
  expect_error(
    createMarkovModelDataset(factor(c("A", "B", "A")), factor(c("B", "A", "C"))),
    regexp = "Assertion on 'colCur' failed: Must have levels: A,B"
  )
})

test_that("createMarkovModelDataset.data.frame", {
  expect_equal(
    createMarkovModelDataset(
      data.frame(prev = c(1, 2, 1), cur = c(1, 2, 3)),
      colPrev = "prev", colCur = "cur"
    ),
    data.frame(
      prev = c(1, 2, 1),
      cur = c(1, 2, 3),
      prevX1 = c(TRUE, FALSE, TRUE),
      curX1 = c(TRUE, FALSE, FALSE),
      prevX2 = c(FALSE, TRUE, FALSE),
      curX2 = c(FALSE, TRUE, FALSE),
      prevX3 = c(FALSE, FALSE, FALSE),
      curX3 = c(FALSE, FALSE, TRUE)
    )
  )
  expect_error(
    createMarkovModelDataset(
      data.frame(prev = c(1, 2, 1), cur = c(1, 2, 3)),
      colPrev = "prev", colCur = "foo"
    ),
    regexp = "Names must include the elements {'prev','foo'}, but is missing elements {'foo'}",
    fixed = TRUE
  )
})

test_that("simMarkov", {
  d <-
    data.frame(
      id = c(1, 2),
      prAtoA = rep(c(0.2, 0.1), 2),
      prAtoB = rep(c(0.8, 0.9), 2),
      prBtoA = rep(c(0.2, 0.1), 2),
      prBtoB = rep(c(0.8, 0.9), 2)
    )
  simresult <- simMarkov(ui = d, initialState = "A", states = c("A", "B"))
  expect_named(
    simresult,
    c("id", "prAtoA", "prAtoB", "prBtoA", "prBtoB", "current", "previous")
  )
  expect_error(
    simMarkov(ui = d[, 2:5], initialState = "A", states = c("A", "B")),
    regexp = "Could not find an ID column in `ui`; expected a column named 'id' (case-insensitive).",
    fixed = TRUE
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