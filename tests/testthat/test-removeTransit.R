test_that("removeTransit throws error in model with no transit compartment", {
  expect_error(removeTransit(readModelDb("PK_2cmt_des"), "central", transit = "transit"), "Assertion on 'ntransit' failed: Must be of type 'integerish', not 'character'")
})
test_that("removeTransit removes transit compartment", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(1)
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeTransit(modelTest, central = "central", depot = "depot", transit = "transit"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("transit1" %in% mv$state, FALSE)
})
test_that("removeTransit removes ktr in model block", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(1)
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeTransit(modelTest, central = "central", depot = "depot", transit = "transit"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("ktr" %in% mv$lhs, FALSE)
})
test_that("removeTransit modifies ODE for central compartment", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(1)
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeTransit(modelTest, central = "central", depot = "depot"))
  suppressMessages(centralLine <- rxode2::modelExtract(modelUpdate, "d/dt(central)", lines = TRUE))
  expect_equal(grepl(".*<-\\s*", centralLine), TRUE)
})
test_that("removeTransit removes lktr in ini block", {
  modelTest <- readModelDb("PK_1cmt_des") |> addTransit(1)
  modelTest <- rxode2::assertRxUi(modelTest)
  suppressMessages(modelUpdate <- removeTransit(modelTest, central = "central", depot = "depot"))
  temp <- rxode2::assertRxUi(modelUpdate)
  mv <- rxode2::rxModelVars(temp)
  expect_equal("lktr1" %in% mv$params, FALSE)
})


test_that("extreme model cases", {

  f <- function() {
    model({
      ktr <- exp(lktr)
      d/dt(depot) <- -ktr * depot
      d/dt(transit1) <- ktr * depot - ktr * transit1
      d/dt(transit2) <- ktr * transit1 - ktr * transit2
      d/dt(transit3) <- ktr * transit2 - ktr * transit3
      d/dt(transit4) <- ktr * transit3 - ka * transit4
      d/dt(central) <- ka * transit4 - kel * central
      Cc <- central/vc
    })
  }

  f <- rxode2::rxode2(f)

  expect_error(f %>% removeTransit(), NA)
  expect_warning(f %>% removeTransit(), NA)
  tmp <- f %>% removeTransit()
  expect_true(length(tmp$iniDf$name) == 0)

  f <- function() {
    ini({
      lktr <- 1
    })
    model({
      ktr <- exp(lktr)
      d/dt(depot) <- -ktr * depot
      d/dt(transit1) <- ktr * depot - ktr * transit1
      d/dt(transit2) <- ktr * transit1 - ktr * transit2
      d/dt(transit3) <- ktr * transit2 - ktr * transit3
      d/dt(transit4) <- ktr * transit3 - ka * transit4
      d/dt(central) <- ka * transit4 - kel * central
      Cc <- central/vc
    })
  }

  f <- rxode2::rxode2(f)

  expect_error(f %>% removeTransit(), NA)
  expect_warning(f %>% removeTransit(), NA)
  tmp <- f %>% removeTransit()
  expect_true(length(tmp$iniDf$name) == 0)

  f <- function() {
    ini({
      lktr ~ 1
    })
    model({
      ktr <- exp(lktr)
      d/dt(depot) <- -ktr * depot
      d/dt(transit1) <- ktr * depot - ktr * transit1
      d/dt(transit2) <- ktr * transit1 - ktr * transit2
      d/dt(transit3) <- ktr * transit2 - ktr * transit3
      d/dt(transit4) <- ktr * transit3 - ka * transit4
      d/dt(central) <- ka * transit4 - kel * central
      Cc <- central/vc
    })
  }

  f <- rxode2::rxode2(f)

  expect_error(f %>% removeTransit(), NA)
  expect_warning(f %>% removeTransit(), NA)
  tmp <- f %>% removeTransit()
  expect_true(length(tmp$iniDf$name) == 0)

  f <- rxode2::rxode2(f)

  expect_error(f %>% removeTransit(), NA)

  tmp <- f %>% removeTransit()

  expect_equal(tmp$iniDf$name, "ka")
  expect_equal(tmp$iniDf$ntheta, 1)

  f <- function() {
    ini({
      ka ~ 1
      lktr ~ 1
    })
    model({
      ktr <- exp(lktr)
      d/dt(depot) <- -ktr * depot
      d/dt(transit1) <- ktr * depot - ktr * transit1
      d/dt(transit2) <- ktr * transit1 - ktr * transit2
      d/dt(transit3) <- ktr * transit2 - ktr * transit3
      d/dt(transit4) <- ktr * transit3 - ka * transit4
      d/dt(central) <- ka * transit4 - kel * central
      Cc <- central/vc
    })
  }

  f <- rxode2::rxode2(f)

  expect_error(f %>% removeTransit(), NA)

  tmp <- f %>% removeTransit()

  expect_equal(tmp$iniDf$name, "ka")
  expect_equal(tmp$iniDf$neta1, 1)


})
