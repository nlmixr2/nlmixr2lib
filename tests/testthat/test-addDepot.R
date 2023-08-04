test_that("addDepot adds depot", {
  model <- readModelDb("PK_2cmt_des")
  suppressMessages(modelUpdate <- addDepot(model, central="central",depot="depot"))
  #check for lka in ini block
  suppressMessages(temp <- modelUpdate)
  temp2<-temp$iniDf
  expect_equal("lka" %in% temp2$name,TRUE)
  
  #check for lfdepot in ini block
  suppressMessages(temp <- modelUpdate)
  temp2<-temp$iniDf
  expect_equal("lfdepot" %in% temp2$name,TRUE)
  
  #check for labels for lka
  suppressMessages(temp <- modelUpdate)
  temp2<-temp$iniDf
  expect_equal("First order absorption rate (ka)" %in% temp2$label,TRUE)
  
  #check for labels for lfdepot
  suppressMessages(temp <- modelUpdate)
  temp2<-temp$iniDf
  expect_equal("Proportional residual error (fraction)" %in% temp2$label,TRUE)
  
  #check for ka in model block
  suppressMessages(kaLine <- rxode2::modelExtract(modelUpdate,"ka",lines = TRUE))
  expect_equal(grepl("\\s*^ka",kaLine),TRUE)
  
  #check for ODE for depot
  suppressMessages(depotLine <- rxode2::modelExtract(modelUpdate,"d/dt(depot)",lines = TRUE))
  expect_equal(grepl("\\s*ka\\s*\\*\\s*depot",depotLine),TRUE)
  
})


test_that("addDepot adds other than default agruments",{
  model <- readModelDb("PK_2cmt_mm")
  
  # suppressMessages(modelUpdate <- addDepot(model, central="cent",depot="depot",absRate = "ka"))
  # temp <- rxode2::assertRxUi(modelUpdate)
  # mv <- rxode2::rxModelVars(temp)
  # expect_equal("cent" %in% mv$state,TRUE)
  # 
  # suppressMessages(modelUpdate <- addDepot(model, central="central",depot="depo",absRate = "ka"))
  # temp <- rxode2::assertRxUi(modelUpdate)
  # mv <- rxode2::rxModelVars(temp)
  # expect_equal("depo" %in% mv$state,TRUE)
  
  suppressMessages(modelUpdate <- addDepot(model, central="central",depot="depot",absRate = "ktr"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("ktr" %in% mv$lhs,TRUE)
})







































