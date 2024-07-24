test_that("addDepot adds depot", {
  suppressMessages(
    model <- readModelDb("PK_1cmt_des") |> removeDepot()
  )
  suppressMessages(
    modelUpdate <- addDepot(model, central = "central", depot = "depot")
  )
  # check for change in lka in ini block
  expect_false("lka" %in% model$iniDf$name)
  expect_true("lka" %in% modelUpdate$iniDf$name)

  suppressMessages(
    model <- readModelDb("PK_1cmt_des") |> removeDepot()
  )
  suppressMessages(
    modelUpdate <- addDepot(model, central = "central", depot = "depot")
  )
  # check for change in lka in ini block
  expect_false("lka" %in% model$iniDf$name)
  expect_true("lka" %in% modelUpdate$iniDf$name)

  # check for labels for lka
  expect_true("First order absorption rate (ka)" %in% modelUpdate$iniDf$label)

  # check for ka in model block
  suppressMessages(kaLine <- rxode2::modelExtract(modelUpdate, "ka", lines = TRUE))
  expect_true(grepl("\\s*^ka", kaLine))

  # check for ODE for depot
  suppressMessages(depotLine <- rxode2::modelExtract(modelUpdate, "d/dt(depot)", lines = TRUE))
  expect_true(grepl("\\s*ka\\s*\\*\\s*depot", depotLine))
})

test_that("addDepot adds other than default agruments", {
  suppressMessages(
    model <- readModelDb("PK_1cmt_des") |> removeDepot()
  )
  suppressMessages(
    modelUpdate <- addDepot(model, central = "central", depot = "depot", ka = "ktr")
  )
  mv <- rxode2::rxModelVars(modelUpdate)
  expect_true("ktr" %in% mv$lhs)
})
