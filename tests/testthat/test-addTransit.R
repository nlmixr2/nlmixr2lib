test_that("addTransit adds transit compartment", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(3)
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- addTransit(modelTest, 1))
  # check for lktr1 in ini block
  temp <- rxode2::assertRxUi(modelUpdate)
  temp2 <- temp$iniDf
  expect_equal("lktr" %in% temp2$name, TRUE)

  # check for ktr in model block
  suppressMessages(kaLine <- rxode2::modelExtract(modelUpdate, "ktr", lines = TRUE))
  expect_equal(grepl("\\s*^ktr", kaLine), TRUE)
})


# Test if the function throws an error when invalid input for 'central' is provided
test_that("addTransit throws an error with invalid 'central'", {
  expect_error(addTransit(readModelDb("PK_2cmt_des"), 3, "cent", "depot"), "'cent' needs to be in the model")
})

# Test if the function throws an error when invalid input for 'transit' is provided
test_that("addTransit throws an error with invalid 'transit'", {
  expect_error(addTransit(readModelDb("PK_2cmt_des"), -1), "Assertion on 'transit' failed: Element 1 is not >= 1")
})

# Test if the function adds transit compartments correctly when 'depot' is present
test_that("addTransit adds transit compartments correctly with 'depot'", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(3)
  modelTest <- rxode2::assertRxUi(modelTest)
  modelUpdate <- addTransit(modelTest, 3)
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("transit1" %in% mv$state, TRUE)
  expect_equal("transit2" %in% mv$state, TRUE)
  expect_equal("transit3" %in% mv$state, TRUE)
})

# Test if the function adds transit compartments  when 'depot' is not present
test_that("addTransit does not add transit compartments  without 'depot'", {
  expect_warning(addTransit(readModelDb("PK_2cmt_no_depot"), 1), "'depot' added to model for transit model")
})
