test_that("addEta named parameter", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addEta(model, eta = "lka"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("etalka ~ 0.1")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[2]],
    str2lang("ka <- exp(lka + etalka)")
  )
})

test_that("addEta mu-ref parameter", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addEta(model, eta = "ka"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("etalka ~ 0.1")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[2]],
    str2lang("ka <- exp(lka + etalka)")
  )
})

test_that("addEta named parameter", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addEta(model, eta = "lka"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("etalka ~ 0.1")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[2]],
    str2lang("ka <- exp(lka + etalka)")
  )
})

test_that("addEta non-existent parameter", {
  model <- readModelDb("PK_1cmt")
  expect_error(
    addEta(model, eta = "foo"),
    "cannot find parameter 'etafoo'"
  )
})
