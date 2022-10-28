test_that("addResErr with des model, changing to additive error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "addSd"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("addSd <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("cp ~ add(addSd)")
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
    str2lang("addSd <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[6]],
    str2lang("cp ~ add(addSd)")
  )
})

test_that("addResErr with des model, changing to additive and proportinoal error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "addSd+propSd"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("propSd <- c(0, 0.5)")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("addSd <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("cp ~ add(addSd) + prop(propSd)")
  )
})

test_that("addResErr bad error type", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(expect_error(
    addResErr(model, reserr = "foo"),
    regexp = "unknown residual error"
  ))
  suppressMessages(expect_error(
    addResErr(model, reserr = 5),
    regexp = "reserr must be character"
  ))
})
