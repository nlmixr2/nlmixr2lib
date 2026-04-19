test_that("checkModelConventions returns a data.frame with the expected columns", {
  res <- suppressWarnings(checkModelConventions("PK_1cmt", verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("model", "category", "severity", "name", "message", "suggestion"),
    ignore.order = TRUE
  )
  expect_true(all(res$severity %in% c("error", "warning", "info")))
})

test_that("conventional parameter names produce no naming issues", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("Absorption rate (ka, 1/day)")
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      etalcl ~ 0.09
      etalvc ~ 0.09
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl + etalcl)
      vc <- exp(lvc + etalvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  naming <- res[res$category == "parameter_naming", ]
  expect_equal(nrow(naming), 0)
})

test_that("etacl on lcl is flagged with etalcl suggestion", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a (1/day)")
      lcl <- 1;   label("b (L/day)")
      lvc <- 1;   label("c (L)")
      etacl ~ 0.09
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl + etacl)
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  naming <- res[res$category == "parameter_naming" & res$name == "etacl", ]
  expect_equal(nrow(naming), 1)
  expect_match(naming$message, "etalcl")
})

test_that("deprecated residual error names are flagged", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1;   label("b")
      lvc <- 1;   label("c")
      prop.err <- 0.1; label("d")
      add.err  <- 0.01; label("e")
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ add(add.err) + prop(prop.err)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  expect_true("prop.err" %in% dep$name)
  expect_true("add.err" %in% dep$name)
})

test_that("missing units$concentration is flagged as an error when residual error exists", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1;   label("b")
      lvc <- 1;   label("c")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  u <- res[res$category == "units" & res$name == "concentration", ]
  expect_equal(nrow(u), 1)
  expect_equal(u$severity, "error")
})

test_that("a covariate used but not in covariateData is an error", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("a"); lvc <- 1; label("b")
      propSd <- 0.1; label("d")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^0.75
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  cov <- res[res$category == "covariates" & res$name == "WT", ]
  expect_true(any(cov$severity == "error"))
})

test_that("a covariate alias without declared source_name produces a warning", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      ADA = list(description = "Anti-drug antibody", units = "(binary)",
                 type = "binary")
    )
    ini({
      lcl <- 1; label("a"); lvc <- 1; label("b")
      e_ada_cl <- 0.5; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      cl <- exp(lcl) * exp(e_ada_cl * ADA)
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  cov <- res[res$category == "covariates" & res$name == "ADA" &
               grepl("alias", res$message), ]
  expect_equal(nrow(cov), 1)
  expect_equal(cov$severity, "warning")
})

test_that("declared alias via source_name produces no alias warning", {
  good <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      ADA_POS = list(description = "ADA-positive", units = "(binary)",
                     type = "binary", source_name = "ADA")
    )
    ini({
      lcl <- 1; label("a"); lvc <- 1; label("b")
      e_ada_cl <- 0.5; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      cl <- exp(lcl) * exp(e_ada_cl * ADA_POS)
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  cov <- res[res$category == "covariates", ]
  expect_equal(nrow(cov), 0)
})

test_that("non-canonical compartment is flagged as a warning", {
  bad <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1; label("b"); lvc <- 1; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      kel <- cl / vc
      d/dt(absorption) <- -ka * absorption
      d/dt(central) <- ka * absorption - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  cmt <- res[res$category == "compartments" & res$name == "absorption", ]
  expect_equal(nrow(cmt), 1)
  expect_equal(cmt$severity, "warning")
})

test_that("warning is emitted when issues exist and suppressed when clean", {
  good <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a (1/day)")
      lcl <- 1; label("b (L/day)")
      lvc <- 1; label("c (L)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  expect_no_warning(checkModelConventions(good, verbose = FALSE))

  bad <- function() {
    ini({
      lcl <- 1
      propSd <- 0.1
    })
    model({
      cl <- exp(lcl)
      d/dt(central) <- -cl * central
      Cc <- central
      Cc ~ prop(propSd)
    })
  }
  expect_warning(
    checkModelConventions(bad, verbose = FALSE),
    "convention issue"
  )
})

test_that("iterating with no argument returns stacked data.frame with model column", {
  res <- suppressWarnings(checkModelConventions(verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_true("model" %in% names(res))
  expect_gt(length(unique(res$model)), 1)
})

test_that("canonical covariates are parsed from inst/references/covariate-columns.md", {
  canon <- nlmixr2lib:::.loadCanonicalCovariates()
  expect_true(length(canon) > 20)
  expect_true(all(c("WT", "SEXF", "ADA_POS", "RACE_BLACK",
                    "RACE_BLACK_OTH", "CREAT", "ALB", "hsCRP") %in%
                    names(canon)))
  expect_true("SEXM" %in% canon$SEXF$aliases)
  expect_true("ADA" %in% canon$ADA_POS$aliases)
  expect_true("BLACK_OTH" %in% canon$RACE_BLACK_OTH$aliases)
  # ALB has no source aliases in the register ("none; ALB is the universal...")
  expect_equal(length(canon$ALB$aliases), 0)
})

test_that("covariate alias map resolves document-order last-writes-win", {
  map <- nlmixr2lib:::.nlmixr2libCovariateAliasMap()
  expect_equal(map[["BLACK_OTH"]], "RACE_BLACK_OTH")
  expect_equal(map[["SEXM"]], "SEXF")
  expect_equal(map[["ADA"]], "ADA_POS")
})

test_that("an rxUi object is accepted directly", {
  good <- function() {
    description <- "A"; reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a (1/day)")
      lcl <- 1; label("b (L/day)")
      lvc <- 1; label("c (L)")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  ui <- nlmixr2est::nlmixr(good)
  res <- suppressWarnings(checkModelConventions(ui, verbose = FALSE))
  expect_s3_class(res, "data.frame")
})
