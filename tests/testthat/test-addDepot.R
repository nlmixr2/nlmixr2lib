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

  # Add a depot when there are no parameters in the model
  modelNoDepot <-
    suppressMessages(
      readModelDb("PK_1cmt_des") |>
        rxode2::model(-Cc~.) |>
        rxode2::ini(-lka, -lcl, -lvc) |>
        removeDepot()
    )
  expect_s3_class(modelNoDepot, "rxUi")
  expect_equal(modelNoDepot$state, "central")
  modelDepot <-
    suppressMessages(
      addDepot(modelNoDepot)
    )
  expect_s3_class(modelDepot, "rxUi")
  expect_equal(modelDepot$state, c("depot", "central"))
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

# Helper: collapse expression list to deparsed character vector for order checks
.deparseLines <- function(lstExpr) {
  vapply(lstExpr, function(e) paste(deparse(e), collapse = " "), character(1))
}

# Helper: TRUE if `needle` lines appear in `haystack` in the same relative order
.isSubseq <- function(haystack, needle) {
  i <- 1L
  for (h in haystack) {
    if (i > length(needle)) break
    if (identical(h, needle[[i]])) i <- i + 1L
  }
  i > length(needle)
}

# Issue #77 / #78 — addDepot edge cases ----

test_that("addDepot works when d/dt(central) is the only / last model line (#77)", {
  # #77 "central at end": residual error precedes d/dt(central)
  m <- function() {
    ini({
      lcl <- 1; label("CL")
      lvc <- 3.45; label("V")
      propSd <- 0.5; label("prop err")
    })
    model({
      central ~ prop(propSd)
      d/dt(central) <- -exp(lcl)/exp(lvc) * central
    })
  }
  res <- addDepot(m)
  expect_s3_class(res, "rxUi")
  expect_true("lka" %in% res$iniDf$name)
  lines <- .deparseLines(res$lstExpr)
  # d/dt(depot) appears exactly once
  expect_equal(sum(grepl("^d/dt\\(depot\\)", lines)), 1L)
  # endpoint kept in its original position
  expect_true(grepl("^central ~ prop", lines[[1]]))
  # central ODE is last and has the absorption term added
  expect_true(grepl("^d/dt\\(central\\)", lines[[length(lines)]]))
  expect_true(grepl("\\+ ka \\* depot", lines[[length(lines)]]))
})

test_that("addDepot works when d/dt(central) is the first model line (#77)", {
  m <- function() {
    ini({
      lcl <- 1; label("CL")
      lvc <- 3.45; label("V")
      propSd <- 0.5; label("prop err")
    })
    model({
      d/dt(central) <- -exp(lcl)/exp(lvc) * central
      central ~ prop(propSd)
    })
  }
  res <- addDepot(m)
  expect_s3_class(res, "rxUi")
  expect_true("lka" %in% res$iniDf$name)
  lines <- .deparseLines(res$lstExpr)
  expect_equal(sum(grepl("^d/dt\\(depot\\)", lines)), 1L)
  # endpoint kept as the last line
  expect_true(grepl("^central ~ prop", lines[[length(lines)]]))
})

test_that("addDepot preserves the relative order of interleaved residual-error and assignment lines", {
  # User-stated constraint from issue #77 review: multiple residual-error
  # blocks interleaved with assignments must keep their source-order
  # semantics. Using a non-state variable name so rxode2 accepts the
  # repeated assignment form.
  m <- function() {
    ini({
      lcl <- 1; lvc <- 3.45; addSd <- 0.1; propSd <- 0.5
    })
    model({
      obs <- 1
      obs ~ add(addSd)
      obs <- 2
      d/dt(central) <- -exp(lcl)/exp(lvc) * central
      obs ~ prop(propSd)
    })
  }
  res <- addDepot(m)
  expect_s3_class(res, "rxUi")
  lines <- .deparseLines(res$lstExpr)
  expect_equal(sum(grepl("^d/dt\\(depot\\)", lines)), 1L)
  # Every input line (except the modified d/dt(central)) appears verbatim and
  # in its original relative order.
  original <- c(
    "obs <- 1",
    "obs ~ add(addSd)",
    "obs <- 2",
    "obs ~ prop(propSd)"
  )
  normalized <- gsub("\\s+", " ", trimws(lines))
  expect_true(.isSubseq(normalized, original))
  # Two new lines sit immediately before the modified central ODE.
  wCentral <- which(grepl("^d/dt\\(central\\)", normalized))
  expect_length(wCentral, 1L)
  expect_equal(normalized[[wCentral - 2L]], "ka <- exp(lka)")
  expect_equal(normalized[[wCentral - 1L]], "d/dt(depot) <- -ka * depot")
})

test_that("addDepot works when a transit ODE is the first model line (#78)", {
  m <- function() {
    ini({
      lktr <- 0; lcl <- 1; lvc <- 3.45; propSd <- 0.5
    })
    model({
      d/dt(transit1) <- -exp(lktr) * transit1
      d/dt(central) <- exp(lktr) * transit1 - exp(lcl)/exp(lvc) * central
      Cc <- central / exp(lvc)
      Cc ~ prop(propSd)
    })
  }
  res <- addDepot(m)
  expect_s3_class(res, "rxUi")
  expect_true("lka" %in% res$iniDf$name)
  lines <- .deparseLines(res$lstExpr)
  normalized <- gsub("\\s+", " ", trimws(lines))
  # transit ODE still at position 1
  expect_true(grepl("^d/dt\\(transit1\\)", normalized[[1]]))
  # endpoint still last
  expect_true(grepl("^Cc ~ prop", normalized[[length(normalized)]]))
  expect_equal(sum(grepl("^d/dt\\(depot\\)", normalized)), 1L)
})

test_that("addDepot works when a transit ODE sits above d/dt(central) (#78)", {
  m <- function() {
    ini({
      lktr <- 0; lcl <- 1; lvc <- 3.45; propSd <- 0.5
    })
    model({
      kel <- exp(lcl) / exp(lvc)
      d/dt(transit1) <- -exp(lktr) * transit1
      d/dt(central) <- exp(lktr) * transit1 - kel * central
      Cc <- central / exp(lvc)
      Cc ~ prop(propSd)
    })
  }
  res <- addDepot(m)
  expect_s3_class(res, "rxUi")
  lines <- .deparseLines(res$lstExpr)
  normalized <- gsub("\\s+", " ", trimws(lines))
  # Original lines keep their relative order.
  original <- c(
    "kel <- exp(lcl)/exp(lvc)",
    "d/dt(transit1) <- -exp(lktr) * transit1",
    "Cc <- central/exp(lvc)",
    "Cc ~ prop(propSd)"
  )
  expect_true(.isSubseq(normalized, original))
  # transit1 ODE appears above the modified central ODE
  wTransit <- which(grepl("^d/dt\\(transit1\\)", normalized))
  wCentral <- which(grepl("^d/dt\\(central\\)", normalized))
  expect_true(wTransit < wCentral)
  expect_equal(sum(grepl("^d/dt\\(depot\\)", normalized)), 1L)
})
