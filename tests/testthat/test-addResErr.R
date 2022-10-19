test_that("addResErr with des model, changing to additive error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "add"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("add.err <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("cp ~ add(add.err)")
  )
})

test_that("addResErr with linCmt model, changing to additive error", {
  model <- readModelDb("PK_1cmt")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "add"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("add.err <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[6]],
    str2lang("cp ~ add(add.err)")
  )
})

test_that("addResErr with des model, changing to additive and proportinoal error", {
  model <- readModelDb("PK_1cmt_des")
  suppressMessages(modelUpdate <- addResErr(model, reserr = "add+prop"))
  # initial conditions are added
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[8]],
    str2lang("prop.err <- c(0, 0.5)")
  )
  expect_equal(
    functionBody(
      modelUpdate
    )[[3]][[2]][[10]],
    str2lang("add.err <- c(0, 1)")
  )
  # eta is added
  expect_equal(
    functionBody(
      modelUpdate
    )[[4]][[2]][[9]],
    str2lang("cp ~ add(add.err) + prop(prop.err)")
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
