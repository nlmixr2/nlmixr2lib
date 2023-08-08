test_that("removeTransit throws error in model with no transit compartment",{
  
  expect_error(removeTransit(readModelDb("PK_2cmt_des"), "central", "depot",transitComp ="transit"), "Assertion on 'transit' failed: Must be of type 'integerish', not 'character'")
})
test_that("removeTransit throws error in model with invalid central compartment",{
  
  expect_error(removeTransit(readModelDb("PK_2cmt_des"), "cent", "depot"), "Assertion on 'transit' failed: Must be of type 'integerish', not 'character'")
})
test_that("removeTransit removes transit compartment", {
  modelTest <- readModelDb("indirect_0cpt_transitEx")
  suppressMessages(modelUpdate <- removeTransit(modelTest, central="central",depot="depot", transitComp ="transit"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("transit1" %in% mv$state,FALSE)
})
test_that("removeTransit removes ktr1 in model block",{
  modelTest <- readModelDb("indirect_0cpt_transitEx")
  suppressMessages(modelUpdate <- removeTransit(modelTest, central="central",depot="depot", transitComp ="transit"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("ktr1" %in% mv$lhs,FALSE)
})
test_that("removeTransit modifies ODE for central compartment", {
  modelTest <- readModelDb("indirect_0cpt_transitEx")
  suppressMessages(modelUpdate <- removeTransit(modelTest, central="central",depot="depot"))
  suppressMessages(centralLine <- rxode2::modelExtract(modelUpdate,"d/dt(central)",lines = TRUE))
  expect_equal(grepl(".*<-\\s*",centralLine),TRUE)
})
test_that("removeTransit removes lktr in ini block",{
  modelTest <- readModelDb("indirect_0cpt_transitEx")
  suppressMessages(modelUpdate <- removeTransit(modelTest, central="central",depot="depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("lktr1" %in% mv$params,FALSE)
})

