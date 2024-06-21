test_that("addResErr with each expected residual distribution and combinations", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdateAdd <- addResErr(model, reserr = "addSd"))
  suppressMessages(modelUpdateProp <- addResErr(model, reserr = "propSd"))
  suppressMessages(modelUpdateLnorm <- addResErr(model, reserr = "lnormSd"))
  suppressMessages(modelUpdateAll <- addResErr(model, reserr = c("addSd", "propSd", "lnormSd")))

  expect_equal(
    functionBody(
      modelUpdateAdd
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcaddSd)")
  )
  expect_equal(
    functionBody(
      modelUpdateProp
    )[[4]][[2]][[6]],
    str2lang("Cc ~ prop(CcpropSd)")
  )
  expect_equal(
    functionBody(
      modelUpdateLnorm
    )[[4]][[2]][[6]],
    str2lang("Cc ~ lnorm(CclnormSd)")
  )
  expect_equal(
    functionBody(
      modelUpdateAll
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcaddSd) + prop(CcpropSd) + lnorm(CclnormSd)")
  )
  # Check that initial conditions were set correctly
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAdd)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcaddSd = 1)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateProp)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcpropSd = 0.5)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateLnorm)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CclnormSd = 0.5)
  ))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAll)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcaddSd = 1, CcpropSd = 0.5, CclnormSd = 0.5)
  ))
})

test_that("addResErr with named numeric values sets the reserr", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdateAdd <- addResErr(model, reserr = c(addSd = 10)))
  suppressMessages(expect_equal(
    nlmixr2est::fixef(rxode2::rxode(modelUpdateAdd)),
    c(lka = 0.45, lcl = 1, lvc = 3.45, CcaddSd = 10)
  ))
})

test_that("addResErr with des model, changing to additive error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "addSd"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[8]],
    str2lang("CcaddSd <- c(0, 1)")
  )
  # residual error model is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[5]][[2]][[9]],
    str2lang("Cc ~ add(CcaddSd)")
  )
})

test_that("addResErr with linCmt model, changing to additive error", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "addSd"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("CcaddSd <- c(0, 1)")
  )
  # residual error model is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[6]],
    str2lang("Cc ~ add(CcaddSd)")
  )
})

test_that("addResErr with des model, changing to additive and proportional error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = c("addSd", "propSd")))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[8]],
    str2lang("CcaddSd <- c(0, 1)")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("CcpropSd <- c(0, 0.5)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[5]][[2]][[9]],
    str2lang("Cc ~ add(CcaddSd) + prop(CcpropSd)")
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
