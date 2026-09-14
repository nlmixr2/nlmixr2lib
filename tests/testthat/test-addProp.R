f <- function() {
  description <- "A two compartment model with a direct effect , no endpoints and no thetas"
  model({
    d / dt(central) <- -kel * central - k12 * central + k21 *
      peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc
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

test_that("addLogitBioavailability for a single compartment", {

  res <- res0 |> addLogitBioavailability("depot")

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
    "f(depot) <- fDepot")

  res <- res0 |> addLogitBioavailability(depot)

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
    "f(depot) <- fDepot")

  expect_equal(rxode2::modelExtract(res, "fDepot"),
    "fDepot <- expit(lgfDepot, 0, 1)")

  expect_error(res0 |> addLogitBioavailability(matt))

  expect_error(addLogitBioavailability())

})

test_that("addLogitBioavailability keeps f bounded and mu-referenced", {

  res <- res0 |> addLogitBioavailability(depot)

  .id <- res$iniDf
  .w <- which(.id$name == "lgfDepot")

  expect_equal(length(.w), 1L)
  # est is on the logit scale of the default f=0.8
  expect_equal(.id$est[.w], logit(0.8))
  expect_false(.id$fix[.w])

  .cur <- res$muRefCurEval
  .w2 <- which(.cur$parameter == "lgfDepot")

  expect_equal(length(.w2), 1L)
  expect_equal(.cur$curEval[.w2], "expit")
  expect_equal(.cur$low[.w2], 0)
  expect_equal(.cur$hi[.w2], 1)

  res2 <- res0 |> addLogitBioavailability(depot, f = 0.25)

  expect_equal(res2$iniDf$est[which(res2$iniDf$name == "lgfDepot")],
    logit(0.25))

})

test_that("addLogitBioavailability splits the dose between two paths", {

  res <- res0 |>
    addDepot(depot = "depot2", ka = "ka2") |>
    addLogitBioavailability(depot, depot2)

  expect_equal(rxode2::modelExtract(res, "f(depot)"),
    "f(depot) <- fDepot")

  expect_equal(rxode2::modelExtract(res, "f(depot2)"),
    "f(depot2) <- 1 - fDepot")

  # a single fraction drives both paths
  expect_equal(sum(res$iniDf$name == "lgfDepot"), 1L)
  expect_equal(sum(grepl("^lgfDepot", res$iniDf$name)), 1L)

  expect_error(res0 |> addLogitBioavailability(depot, depot))

  expect_error(res0 |> addLogitBioavailability(depot, matt))

})

test_that("addLogitBioavailability dose split solves with the expected ratio", {

  res <- res0 |>
    addDepot(depot = "depot2", ka = "ka2") |>
    addLogitBioavailability(depot, depot2, f = 0.7)

  .mkEv <- function(.cmt) {
    data.frame(ID = 1, time = 0, evid = 1, cmt = .cmt, amt = 100)
  }
  .p <- c(lka = log(1.2), lka2 = log(1.2), lcl = log(0.1),
    lvc = log(10), lgfDepot = logit(0.7))

  .sA <- as.data.frame(rxode2::rxSolve(res, .mkEv("depot"), params = .p))
  .sB <- as.data.frame(rxode2::rxSolve(res, .mkEv("depot2"), params = .p))

  .tEnd <- max(.sA$time)
  .cA <- .sA$Cc[.sA$time == .tEnd]
  .cB <- .sB$Cc[.sB$time == .tEnd]

  # with matched ka the late concentration ratio is F1/(1-F1) = 0.7/0.3
  expect_equal(.cA / .cB, 0.7 / 0.3, tolerance = 1e-5)

})

test_that("addLogitBioavailability rejects f outside (0,1)", {

  expect_error(res0 |> addLogitBioavailability(depot, f = 0),
    regexp = "f must be in")

  expect_error(res0 |> addLogitBioavailability(depot, f = 1),
    regexp = "f must be in")

  expect_error(res0 |> addLogitBioavailability(depot, f = 1.5),
    regexp = "f must be in")

  expect_error(res0 |> addLogitBioavailability(depot, f = -0.2),
    regexp = "f must be in")

})

test_that("addLogitBioavailability removes the description", {

  res <- res0 |> addLogitBioavailability(depot)

  expect_null(res$meta$description)

})

test_that("addLogitBioavailability refuses to double-define f", {

  # addBioavailability() leaves fDepot <- exp(lfDepot); applying the
  # logit form on top would leave two fDepot definitions and the
  # solved model would silently reflect neither
  expect_error(
    res0 |> addBioavailability(depot) |> addLogitBioavailability(depot),
    regexp = "bioavailability already present for compartment 'depot'"
  )

  expect_error(
    res0 |>
      addDepot(depot = "depot2", ka = "ka2") |>
      addLogitBioavailability(depot) |>
      addLogitBioavailability(depot, depot2),
    regexp = "bioavailability already present for compartment 'depot'"
  )

  # the guard reads the model's own state properties, so it also sees
  # an f() written directly in a library seed model
  .da <- rxode2::rxode2(readModelDb("PK_double_sim_01"))

  expect_error(.da |> addLogitBioavailability("depot1"),
    regexp = "bioavailability already present for compartment 'depot1'"
  )

  # and it sees the f() this function itself added to both paths
  .split <- res0 |>
    addDepot(depot = "depot2", ka = "ka2") |>
    addLogitBioavailability(depot, depot2)

  expect_equal(.split$props$cmtProp$Property, c("f", "f"))

  expect_error(.split |> addLogitBioavailability(depot2),
    regexp = "bioavailability already present for compartment 'depot2'"
  )

})

test_that("addLogitBioavailability treats a 'NULL' string cmt2 as a compartment name", {

  # a quoted "NULL" is not the NULL default; it must fail as an
  # unknown compartment instead of silently becoming a single-path call
  expect_error(res0 |> addLogitBioavailability(depot, "NULL"),
    regexp = "not in the model"
  )

})

test_that("addLogitBioavailability does not inherit template backTransform", {

  res <- res0 |> addLogitBioavailability(depot)

  .id <- res$iniDf
  .w <- which(.id$name == "lgfDepot")

  expect_true(is.na(.id$backTransform[.w]))
  expect_true(is.na(.id$condition[.w]))
  expect_true(is.na(.id$prior[.w]))

})
