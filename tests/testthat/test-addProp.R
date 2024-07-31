f <- function() {
  description <- "A two compartment model with a direct effect , no endpoints and no thetas"
  model({
    d/dt(central) <- -kel * central - k12 * central + k21 *
      peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central/vc
  })

}

test_that("addCmpProp removes description", {
  res <- rxode2::rxode2(f) |> addCmtProp("f", "central")
  expect_null(res$meta$description)
})

res0 <- readModelDb("PK_1cmt_des")

test_that("control --nothing exists", {
  expect_equal(rxode2::modelExtract(res0, "f(depot)"),
               character(0))
})

test_that("addCmtProp test for F", {

  res <- res0 |> addCmtProp("f", "depot")

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
               "f(depot) <- fDepot")


  res <- res0 |> addCmtProp("f", depot)

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
               "f(depot) <- fDepot")

  res <- res0 |> addBioavailability("depot")

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
               "f(depot) <- fDepot")

  res <- res0 |> addBioavailability(depot)

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
               "f(depot) <- fDepot")

  expect_error(res0 |> addBioavailability(matt))

  expect_error(addBioavailability())

})

test_that("addCmtProp test for Dur", {

  res <- res0 |> addCmtProp("dur", "depot")

  expect_equal(rxode2::modelExtract(res, "dur(depot)"),
               "dur(depot) <- durDepot")


  res <- res0 |> addCmtProp("dur", depot)

  expect_equal(rxode2::modelExtract(res, "dur(depot)"),
               "dur(depot) <- durDepot")

  res <- res0 |> addDur("depot")

  expect_equal(rxode2::modelExtract(res, "dur(depot)"),
               "dur(depot) <- durDepot")

  res <- res0 |> addDur(depot)

  expect_equal(rxode2::modelExtract(res, "dur(depot)"),
               "dur(depot) <- durDepot")

  expect_error(res0 |> addDur(matt))

  expect_error(addDur())

})


test_that("addCmtProp test for Rate", {

  res <- res0 |> addCmtProp("rate", "depot")

  expect_equal(rxode2::modelExtract(res, "rate(depot)"),
               "rate(depot) <- rateDepot")


  res <- res0 |> addCmtProp("rate", depot)

  expect_equal(rxode2::modelExtract(res, "rate(depot)"),
               "rate(depot) <- rateDepot")

  res <- res0 |> addRate("depot")

  expect_equal(rxode2::modelExtract(res, "rate(depot)"),
               "rate(depot) <- rateDepot")

  res <- res0 |> addRate(depot)

  expect_equal(rxode2::modelExtract(res, "rate(depot)"),
               "rate(depot) <- rateDepot")

  expect_error(res0 |> addRate(matt))

  expect_error(addRate())

})


test_that("addCmtProp test for lag", {

  res <- res0 |> addCmtProp("lag", "depot")

  expect_equal(rxode2::modelExtract(res, "lag(depot)"),
               "lag(depot) <- lagDepot")


  res <- res0 |> addCmtProp("lag", depot)

  expect_equal(rxode2::modelExtract(res, "lag(depot)"),
               "lag(depot) <- lagDepot")

  res <- res0 |> addLag("depot")

  expect_equal(rxode2::modelExtract(res, "lag(depot)"),
               "lag(depot) <- lagDepot")

  res <- res0 |> addLag(depot)

  expect_equal(rxode2::modelExtract(res, "lag(depot)"),
               "lag(depot) <- lagDepot")

  expect_error(res0 |> addLag(matt))

  expect_error(addLag())

})


test_that("addCmtProp test for ini", {

  res <- res0 |> addCmtProp("ini", "depot")

  expect_equal(rxode2::modelExtract(res, "depot(0)"),
               "depot(0) <- iniDepot")


  res <- res0 |> addCmtProp("ini", depot)

  expect_equal(rxode2::modelExtract(res, "depot(0)"),
               "depot(0) <- iniDepot")

  res <- res0 |> addIni("depot")

  expect_equal(rxode2::modelExtract(res, "depot(0)"),
               "depot(0) <- iniDepot")

  res <- res0 |> addIni(depot)

  expect_equal(rxode2::modelExtract(res, "depot(0)"),
               "depot(0) <- iniDepot")

  expect_error(res0 |> addIni(matt))

  expect_error(addIni())

})
