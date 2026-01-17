test_that("addResErr with each expected residual distribution and combinations", {

  model <- readModelDb("PK_1cmt")

  suppressMessages(modelUpdateAdd <- as.function(addResErr(model, reserr = "addSd")))
  suppressMessages(modelUpdateProp <- as.function(addResErr(model, reserr = "propSd")))
  suppressMessages(modelUpdateLnorm <- as.function(addResErr(model, reserr = "lnormSd")))

  suppressMessages(modelUpdateAll <- as.function(addResErr(model, reserr = c("addSd", "propSd", "lnormSd"))))

  expect_equal(
    functionBody(
      modelUpdateAdd
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcAddSd)")
  )

  expect_equal(
    functionBody(
      modelUpdateProp
    )[[4]][[2]][[6]],
    str2lang("Cc ~ prop(CcPropSd)")
  )
  expect_equal(
    functionBody(
      modelUpdateLnorm
    )[[4]][[2]][[6]],
    str2lang("Cc ~ lnorm(CcLnormSd)")
  )
  expect_equal(
    functionBody(
      modelUpdateAll
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcAddSd) + prop(CcPropSd) + lnorm(CcLnormSd)")
  )
  # Check that initial conditions were set correctly
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAdd)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcAddSd = 1)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateProp)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcPropSd = 0.5)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateLnorm)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcLnormSd = 0.5)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAll)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcAddSd = 1, CcPropSd = 0.5, CcLnormSd = 0.5)
  ))
})

test_that("addResErr with named numeric values sets the reserr", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdateAdd <- as.function(addResErr(model, reserr = c(addSd = 10))))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAdd)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcAddSd = 10)
  ))
})

test_that("addResErr with des model, changing to additive error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- as.function(addResErr(model, reserr = "addSd")))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[8]],
    str2lang("CcAddSd <- c(0, 1)")
  )
  # residual error model is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[5]][[2]][[9]],
    str2lang("Cc ~ add(CcAddSd)")
  )
})

test_that("addResErr with linCmt model, changing to additive error", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- as.function(addResErr(model, reserr = "addSd")))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("CcAddSd <- c(0, 1)")
  )
  # residual error model is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcAddSd)")
  )
})

test_that("addResErr with des model, changing to additive and proportional error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- as.function(addResErr(model, reserr = c("addSd", "propSd"))))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[8]],
    str2lang("CcAddSd <- c(0, 1)")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("CcPropSd <- c(0, 0.5)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[5]][[2]][[9]],
    str2lang("Cc ~ add(CcAddSd) + prop(CcPropSd)")
  )
})

test_that("addResErr bad error type", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(expect_error(
    addResErr(model, reserr = "foo"),
    regexp = "Must be a subset of"
  ))
  suppressMessages(expect_error(
    addResErr(model, reserr = 5),
    regexp = "Must have names"
  ))
})

test_that("addResErr with zero endpoints", {
  model <- function() {
    ini({
      a <- 1
      b <- 1
    })
    model({
      e <- exp(a)
      f <- exp(b)
    })
  }
  expect_error(
    addResErr(model, reserr = "addSd"),
    regexp = "there must be at least one prediction in the model({}) block",
    fixed = TRUE
  )
})

test_that("addResErr with multiple endpoints", {
  model <- function() {
    ini({
      a <- 1
      b <- 1
      c <- 1
      d <- 1
    })
    model({
      e <- exp(a)
      f <- exp(b)
      e ~ add(c)
      f ~ add(d)
    })
  }
  expect_error(
    addResErr(model, reserr = "addSd"),
    regexp = "multiple endpoints detected, choose one: e, f",
    fixed = TRUE
  )
  expect_error(
    addResErr(model, reserr = "addSd", endpoint = "g"),
    regexp = "requested to add/change residual error for 'g' but not defined as a modeled endpoint",
    fixed = TRUE
  )
  expect_error(
    addResErr(model, reserr = c("addSd" = factor("A")), endpoint = "f"),
    regexp = "reserr must be a character string or a named numeric vector",
    fixed = TRUE
  )
  addSdE <-
    suppressMessages(
      addResErr(model, reserr = "addSd", endpoint = "e")
    )
  expect_equal(
    functionBody(as.function(addSdE))[[2]][[2]][[5]],
    str2lang("eAddSd <- c(0, 1)")
  )
  expect_equal(
    functionBody(as.function(addSdE))[[3]][[2]][[4]],
    str2lang("e ~ add(eAddSd)")
  )
  # The other endpoint is untouched
  expect_equal(
    functionBody(as.function(addSdE))[[3]][[2]][[5]],
    str2lang("f ~ add(d)")
  )
})
