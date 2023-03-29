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

test_that("addEta multiple parameter, mu-ref and not", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addEta(model, eta = c("lvc", "ka")))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("etalvc ~ 0.1")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[11]],
    str2lang("etalka ~ 0.1")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[2]],
    str2lang("ka <- exp(lka + etalka)")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[4]],
    str2lang("vc <- exp(lvc + etalvc)")
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
  suppressMessages(expect_error(
    addEta(model, eta = "foo")
  ))
})


test_that("compiled ui object", {
  model <- readModelDb("PK_1cmt")
  model <- rxode2::rxode2(model)
  expect_true(inherits(model, "rxUi"))
  suppressMessages(modelUpdate <- addEta(model, eta = "lka"))
  expect_true(inherits(modelUpdate, "rxUi"))
})

test_that("addEta() correctly adds IIV when there is a covariate (#27)", {
  model <- function() {
    ini({
      lka <- 0.45
      lcl <- 1
      lvc <- 3.45
      propSd <- c(0, 0.5)
      allo_cl <- 0.75
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl + allo_cl * WT)
      vc <- exp(lvc)
      cp <- linCmt()
      cp ~ prop(propSd)
    })
  }
  # Update the model detecting the correct parameter for cl
  suppressMessages(
    newEtaRemap <- addEta(model, "cl")
  )
  # Update the model where the correct parameter cor cl is given
  suppressMessages(
    newEta <- addEta(model, "lcl")
  )
  expect_equal(newEtaRemap, newEta, ignore_function_env = TRUE)
})
