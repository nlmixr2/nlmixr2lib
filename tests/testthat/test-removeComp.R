test_that("removeComp throws error in model with no peripheral compartment",{
  
  expect_error(removeComp(readModelDb("PK_1cmt_des")), "'peripheral need to be in the model")
})
test_that("removeComp throws error in model with invalid central compartment",{
  
  expect_error(removeComp(readModelDb("PK_1cmt_des"),central="cent"), "'cent' needs to be in the model")
})
test_that("removeComp removes peripheral compartments", {
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- removeComp(modelTest, central="central",depot="depot", peripheralComp ="peripheral"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("peripheral1" %in% mv$state,FALSE)
})

test_that("removeComp removes k12 in model block",{
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- removeComp(modelTest, central="central",depot="depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("k12" %in% mv$lhs,FALSE)
})

test_that("removeComp removes k21 in model block",{
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- removeComp(modelTest, central="central",depot="depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("k21" %in% mv$lhs,FALSE)
})

test_that("removeComp removes vp in model block",{
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- removeDepot(modelTest, central="central",depot="depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("vp" %in% mv$lhs,FALSE)
})

test_that("removeComp removes q in model block",{
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- removeDepot(modelTest, central="central",depot="depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("q" %in% mv$lhs,FALSE)
})





