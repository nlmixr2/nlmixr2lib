test_that("addComp adds a peripheral compartment", {
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- addComp(modelTest, 1))
  #check for lvp in ini block
  temp <- rxode2::assertRxUi(modelUpdate)
  temp2<-temp$iniDf
  expect_equal("lvp" %in% temp2$name,TRUE)
  
  #check for lq in ini block
  temp <- rxode2::assertRxUi(modelUpdate)
  temp2<-temp$iniDf
  expect_equal("lq" %in% temp2$name,TRUE)
  
  #check for k12 in model block
  suppressMessages(kLine <- rxode2::modelExtract(modelUpdate,"k12",lines = TRUE))
  expect_equal(grepl("\\s*^k12",kLine),TRUE)
  
  #check for ODE for peripheral1
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("peripheral1" %in% mv$state,TRUE)
  
})

# Test if the function throws an error when invalid input for 'central' is provided
test_that("addComp throws an error with invalid 'central'", {
  expect_error(addComp(readModelDb("PK_2cmt_des"), 3, "cent", "depot"), "'cent' needs to be in the model")
})

test_that("addComp throws an error with invalid 'numPeripheral'", {
  expect_error(addComp(readModelDb("PK_2cmt_des"), -1), "Assertion on 'numPeripheral' failed: Element 1 is not >= 1")
})

test_that("addComp removes existing peripheral compartments before adding new", {
  modelTest <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- addComp(modelTest, 1))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("peripheral2" %in% mv$state,FALSE)
  
})



