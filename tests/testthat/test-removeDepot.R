# test_that("removeDepot throws error in model with no depot compartment",{
#
#   expect_error(removeDepot(readModelDb("PK_2cmt_no_depot"), "central", "depot"), "'depot' needs to be in the model")
# })

test_that("removeDepot throws error in model with invalid central compartment", {
  expect_error(removeDepot(readModelDb("PK_1cmt_des"), "cent", "depot"), "'cent' compartment is not in the model")
})

test_that("removeDepot removes depot compartment", {
  modelTest <- readModelDb("PK_1cmt_des")
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeDepot(modelTest, central = "central", depot = "depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("depot" %in% mv$state, FALSE)
  expect_error(removeDepot(temp, "central", "depot"),
               "'depot' compartment is not in the model")
})

test_that("removeDepot removes ka in model block", {
  modelTest <- readModelDb("PK_1cmt_des")
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeDepot(modelTest, central = "central", depot = "depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("ka" %in% mv$lhs, FALSE)
})

test_that("removeDepot modifies ODE for central compartment", {
  modelTest <- readModelDb("PK_1cmt_des")
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeDepot(modelTest, central = "central", depot = "depot"))
  suppressMessages(centralLine <- rxode2::modelExtract(modelUpdate, "d/dt(central)", lines = TRUE))
  expect_equal(grepl(".*<-\\s*", centralLine), TRUE)
})
test_that("removeDepot removes lka in ini block", {
  modelTest <- readModelDb("PK_2cmt_des")
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeDepot(modelTest, central = "central", depot = "depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("lka" %in% mv$params, FALSE)
})
