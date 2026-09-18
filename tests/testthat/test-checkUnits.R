# Model fixtures below use the rxode2/nlmixr2 DSL inside ini() and model()
# blocks: `param <- value; label("...")` on one line and `d/dt(cmt)` without
# spaces around `/`. Both are intentional.
# nolint start: semicolon_linter, infix_spaces_linter

skip_if_not_installed("units")

# A one-compartment oral model with log-scale parameters and an allometric
# clearance; `Cc` is written as an mg/L concentration.
.unitFixture1cmt <- function() {
  ini({
    lka <- 0.1
    lcl <- 1
    lvc <- 1
    etalcl ~ 0.1
    propSd <- 0.1
    addSd <- 1
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^0.75
    vc <- exp(lvc)
    kel <- cl / vc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

.unitRow <- function(res, name) {
  res[res$name == name, , drop = FALSE]
}

# The same model without an eta, so rxSolve() is deterministic.
.unitFixture1cmtFixed <- function() {
  ini({
    lka <- 0.1
    lcl <- 1
    lvc <- 1
    propSd <- 0.1
    addSd <- 1
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl) * (WT / 70)^0.75
    vc <- exp(lvc)
    kel <- cl / vc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}

# ---- units-package wrapper ----------------------------------------------------

test_that("unit strings the library writes parse to canonical spellings", {
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("mg/L")), "mg/L")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("ng/mL")), "ng/mL")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("L/hr")), "L/h")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("mL/kg/hour")), "mL/kg/h")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("1/h")), "1/h")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("IU/mL")), "IU/mL")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("nM")), "nM")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("mcg")), "ug")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("mcg/mL")), "ug/mL")
  micro <- intToUtf8(0xB5)
  middot <- intToUtf8(0xB7)
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse(paste0(micro, "g/mL"))), "ug/mL")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse(paste0("mg", middot, "h/L"))), "mg*h/L")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("L/wk")), "L/week")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("kg^2")), "kg^2")
  # a trailing descriptor is not a factor; grouping parentheses are kept
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("umol/L (uM)")), "umol/L")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("nmol (convert mg via mw)")), "nmol")
  expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse("1/(nM*h)")), "1/(nM*h)")
  for (s in c("unitless", "fraction", "", "1", "none", "n/a", "Unitless")) {
    expect_equal(nlmixr2lib:::.unitDeparse(nlmixr2lib:::.unitParse(s)), "unitless", info = s)
  }
})

test_that("placeholders, fractional exponents and unparseable strings fail parsing", {
  for (s in c("dose_unit", "conc_unit/vol_unit", "time_unit", "mg per L", "AUC")) {
    expect_error(nlmixr2lib:::.unitParse(s), class = "nlmixr2libUnitError", info = s)
  }
  # prose never reaches udunits, which can take minutes over a sentence
  t0 <- Sys.time()
  expect_error(
    nlmixr2lib:::.unitParse("(none; static concentration-response model driven by the substrate covariate)"),
    "not a unit",
    class = "nlmixr2libUnitError"
  )
  expect_error(nlmixr2lib:::.unitParse("pmol/min/pmol recombinant CYP1A2"), "not a unit", class = "nlmixr2libUnitError")
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  expect_true(nlmixr2lib:::.unitLooksLikeUnit("mL/min/1.73 m^2"))
  expect_error(nlmixr2lib:::.unitParse("kg^0.75"), "fractional exponent", class = "nlmixr2libUnitError")
  expect_error(nlmixr2lib:::.unitParse("kg**2"), "fractional exponent", class = "nlmixr2libUnitError")
  expect_error(nlmixr2lib:::.unitParse("mg-h/L"), class = "nlmixr2libUnitError")
})

test_that("convertibility and conversion factors come from udunits", {
  expect_true(nlmixr2lib:::.unitConvertible("mg/L", "ng/mL"))
  expect_equal(nlmixr2lib:::.unitFactor("mg/L", "ng/mL"), 1000)
  expect_equal(nlmixr2lib:::.unitFactor("L/h", "mL/min"), 1000 / 60)
  expect_equal(nlmixr2lib:::.unitFactor("%", "unitless"), 0.01)
  expect_false(nlmixr2lib:::.unitConvertible("mg", "nmol"))
  expect_false(nlmixr2lib:::.unitConvertible("mg/kg", "mg"))
  expect_true(nlmixr2lib:::.unitConvertible("mL/min/1.73m^2", "mL/min/m^2"))
  expect_true(nlmixr2lib:::.unitIsDimensionless(nlmixr2lib:::.unitParse("%")))
  expect_false(nlmixr2lib:::.unitIsDimensionless(nlmixr2lib:::.unitParse("mg")))
})

test_that("custom units install once and stay installed", {
  expect_true(nlmixr2lib:::.unitInstallCustom())
  expect_true(nlmixr2lib:::.unitInstallCustom())
  expect_true(units::ud_are_convertible("IU", "IU"))
  expect_false(units::ud_are_convertible("IU", "IU/mL"))
  expect_true(units::ud_are_convertible("uM", "nmol/L"))
  expect_equal(units::ud_convert(1, "uM", "nmol/L"), 1000)
})

test_that("checkUnits() errors with an installation hint when units is absent", {
  testthat::with_mocked_bindings(
    .unitsAvailable = function() FALSE,
    {
      expect_error(checkUnits(.unitFixture1cmt, time = "h"), "install.packages", class = "nlmixr2libUnitError")
      expect_error(addUnits(.unitFixture1cmt, time = "h"), "install.packages", class = "nlmixr2libUnitError")
    }
  )
})

# ---- checkUnits() -----------------------------------------------------------------

test_that("a consistent one-compartment model resolves every symbol without issues", {
  res <- checkUnits(.unitFixture1cmt, time = "h", depot = "mg", Cc = "mg/L")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("name", "type", "unit", "source", "transformOf", "line", "conversion", "issue"))
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  expect_equal(nrow(attr(res, "conversions")), 0L)
  expected <- c(
    time = "h",
    depot = "mg",
    central = "mg",
    Cc = "mg/L",
    lka = "unitless",
    lcl = "unitless",
    lvc = "unitless",
    propSd = "unitless",
    addSd = "mg/L",
    etalcl = "unitless",
    WT = NA_character_,
    ka = "1/h",
    cl = "L/h",
    vc = "L",
    kel = "1/h",
    "d/dt(depot)" = "mg/h",
    "d/dt(central)" = "mg/h"
  )
  expect_equal(stats::setNames(res$unit, res$name), expected)
  expect_equal(.unitRow(res, "lcl")$transformOf, "cl")
  expect_equal(.unitRow(res, "etalcl")$transformOf, "cl")
  expect_equal(.unitRow(res, "lvc")$transformOf, "vc")
  expect_equal(.unitRow(res, "lka")$transformOf, "ka")
  expect_equal(.unitRow(res, "WT")$source, "compatible")
  expect_equal(.unitRow(res, "cl")$source, "inferred")
  expect_equal(.unitRow(res, "addSd")$source, "boundary")
  expect_equal(.unitRow(res, "depot")$source, "declared")
  expect_equal(.unitRow(res, "Cc")$line, 7L)
  expect_equal(.unitRow(res, "kel")$line, 4L)
  expect_equal(.unitRow(res, "d/dt(central)")$type, "property")
})

test_that("an ng/mL output needs a factor of 1000 that is reported, not an issue", {
  res <- checkUnits(.unitFixture1cmt, time = "h", depot = "mg", Cc = "ng/mL")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "Cc")$conversion, "mg/L * 1000 = ng/mL")
  expect_equal(.unitRow(res, "addSd")$unit, "ng/mL")
  conv <- attr(res, "conversions")
  expect_equal(nrow(conv), 1L)
  expect_equal(conv$target, "Cc")
  expect_equal(conv$factor, 1000)
  expect_equal(conv$line, 7L)
  expect_equal(conv$termIndex, 1L)
})

test_that("a conversion the model already carries is recognised, any other literal is not", {
  withConv <- function() {
    ini({ lka <- 0.1; lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc * 1000
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(withConv, time = "h", depot = "mg", Cc = "ng/mL")
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  otherLiteral <- function() {
    ini({ lka <- 0.1; lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc * 500
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(otherLiteral, time = "h", depot = "mg", Cc = "ng/mL")
  expect_equal(.unitRow(res, "Cc")$conversion, "mg/L * 1000 = ng/mL")
  # a literal on a line whose units already agree is a model constant
  halfLife <- function() {
    ini({ lthalf <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      thalf <- exp(lthalf); vc <- exp(lvc)
      kel <- 0.693 / thalf
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(halfLife, time = "h", central = "mg", Cc = "mg/L")
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  expect_equal(.unitRow(res, "thalf")$unit, "h")
})

test_that("a dose per body weight uses mL/kg defaults and needs the mg/mL factor", {
  res <- checkUnits(.unitFixture1cmt, time = "h", depot = "mg/kg", Cc = "mg/L")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "vc")$unit, "mL/kg")
  expect_equal(.unitRow(res, "cl")$unit, "mL/kg/h")
  expect_equal(.unitRow(res, "central")$unit, "mg/kg")
  expect_equal(.unitRow(res, "d/dt(central)")$unit, "mg/kg/h")
  expect_equal(.unitRow(res, "Cc")$conversion, "mg/mL * 1000 = mg/L")
  expect_equal(nrow(attr(res, "conversions")), 1L)
})

test_that("an overridden clearance unit puts the conversion on the line that uses it", {
  res <- checkUnits(.unitFixture1cmt, time = "h", depot = "mg", Cc = "mg/L", cl = "mL/min")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "cl")$unit, "mL/min")
  expect_equal(.unitRow(res, "cl")$source, "declared")
  expect_equal(.unitRow(res, "kel")$unit, "1/h")
  expect_equal(.unitRow(res, "kel")$conversion, "mL/L/min * 0.06 = 1/h")
  expect_equal(attr(res, "conversions")$factor, 0.06)
})

test_that("default units follow the model's own time unit", {
  res <- checkUnits(.unitFixture1cmt, time = "day", depot = "mg", Cc = "mg/L", WT = "kg")
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  expect_equal(.unitRow(res, "time")$unit, "d")
  expect_equal(.unitRow(res, "cl")$unit, "L/d")
  expect_equal(.unitRow(res, "kel")$unit, "1/d")
  expect_equal(.unitRow(res, "WT")$unit, "kg")
  expect_equal(.unitRow(res, "WT")$source, "declared")
})

test_that("additive and multiplicative etas take the unit their use implies", {
  mult <- function() {
    ini({ tvc <- 1; tcl <- 1; etavc ~ 0.1; etacl ~ 0.1; propSd <- 0.1 })
    model({
      vc <- tvc * exp(etavc)
      cl <- tcl + etacl
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(mult, time = "h", central = "mg", Cc = "mg/L")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "tvc")$unit, "L")
  expect_true(is.na(.unitRow(res, "tvc")$transformOf))
  expect_equal(.unitRow(res, "etavc")$unit, "unitless")
  expect_equal(.unitRow(res, "etavc")$transformOf, "vc")
  expect_equal(.unitRow(res, "tcl")$unit, "L/h")
  expect_equal(.unitRow(res, "etacl")$unit, "L/h")
})

test_that("a variable assigned in both branches of an if must agree", {
  same <- function() {
    ini({ lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      cl <- exp(lcl); vc <- exp(lvc)
      if (WT > 70) { kel <- cl / vc } else { kel <- cl / vc }
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(same, time = "h", central = "mg", Cc = "mg/L")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "kel")$unit, "1/h")
  differ <- function() {
    ini({ lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      cl <- exp(lcl); vc <- exp(lvc)
      if (WT > 70) { kel <- cl / vc } else { kel <- cl }
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(differ, time = "h", central = "mg", Cc = "mg/L")
  expect_match(.unitRow(res, "kel")$issue, "'cl' \\(line 3\\) has units L/h but kel is 1/h")
  expect_error(addUnits(differ, time = "h", central = "mg", Cc = "mg/L"), "kel", class = "nlmixr2libUnitError")
})

test_that("transcendental functions give a unitless result and validate their argument's arithmetic", {
  m <- function() {
    ini({ lcl <- 1; lvc <- 1; propSd <- 0.1; hill <- 1; te50 <- 1 })
    model({
      cl <- exp(lcl); vc <- exp(lvc)
      x <- log(cl)
      y <- exp(-cl / vc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc * exp(-cl / vc * t)
      Cc ~ prop(propSd)
      E <- Cc^hill / (te50^hill + Cc^hill)
      # once cl and vc are known, a sum of them inside exp() is caught
      z <- exp(cl + vc)
    })
  }
  res <- checkUnits(m, time = "h", central = "mg", Cc = "mg/L")
  expect_true(is.na(.unitRow(res, "x")$issue))
  expect_equal(.unitRow(res, "x")$unit, "unitless")
  expect_true(is.na(.unitRow(res, "y")$issue))
  # `exp(cl + vc)` is a modelling error (a clearance plus a volume). The
  # checker reports it either on that line or, when the sum is what first
  # named cl's unit, on the line where the contradiction surfaces.
  iss <- res$issue[!is.na(res$issue)]
  expect_gte(length(iss), 1L)
  expect_true(any(grepl(
    "adds L/h to L, which are not convertible|has units mg but d/dt\\(central\\) is mg/h",
    iss
  )))
  expect_true(is.na(.unitRow(res, "Cc")$issue))
  # a dimensioned base to a symbolic power is unknown, not an issue
  expect_true(is.na(.unitRow(res, "E")$issue))
  expect_true(is.na(.unitRow(res, "E")$unit))
  expect_true(is.na(.unitRow(res, "te50")$unit))
})

test_that("comparisons follow the sum rules", {
  m <- function() {
    ini({ tcl <- 1; lvc <- 1; propSd <- 0.1; tthr <- 10; tfrac <- 0.5 })
    model({
      vc <- exp(lvc)
      cl <- tcl * ifelse(t < 8, 1, tfrac)
      if (AGE <= tthr) {
        cl <- cl * 0.5
      }
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(m, time = "h", central = "mg", Cc = "mg/L", AGE = "year")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "tthr")$unit, "year")
  expect_equal(.unitRow(res, "tfrac")$unit, "unitless")
  bad <- function() {
    ini({ tcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      vc <- exp(lvc)
      cl <- tcl
      if (AGE <= WT) {
        cl <- cl * 0.5
      }
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(bad, time = "h", central = "mg", Cc = "mg/L", AGE = "year", WT = "kg")
  expect_match(res$issue[res$type == "line"], "compares year with kg")
})

test_that("covariate normalisation is unitless and mixed sums are issues", {
  m <- function() {
    ini({ tcl <- 1; lvc <- 1; e_age_cl <- 0.01; propSd <- 0.1 })
    model({
      cl <- tcl * (WT / 70)^0.75 * (1 + e_age_cl * (AGE - 40))
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(m, time = "h", central = "mg", Cc = "mg/L", WT = "kg", AGE = "year")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "tcl")$unit, "L/h")
  # the literal 1 makes the covariate effect unitless, so its slope is per year
  expect_equal(.unitRow(res, "e_age_cl")$unit, "1/year")
  res <- checkUnits(m, time = "h", central = "mg", Cc = "mg/L")
  expect_equal(.unitRow(res, "e_age_cl")$unit, "unitless")
  expect_equal(.unitRow(res, "tcl")$unit, "L/h")
  # an Emax-type sum names the unit of its unknown side
  emax <- function() {
    ini({ te0 <- 1; temax <- 1; tec50 <- 1; lvc <- 1; tkel <- 0.1; propSd <- 0.1; addSd <- 1 })
    model({
      vc <- exp(lvc)
      d/dt(central) <- -tkel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
      E <- te0 * (1 + temax * Cc / (tec50 + Cc))
      E ~ add(addSd)
    })
  }
  res <- checkUnits(emax, time = "h", central = "mg", Cc = "ng/mL", E = "mmHg")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "tec50")$unit, "ng/mL")
  expect_equal(.unitRow(res, "temax")$unit, "unitless")
  expect_equal(.unitRow(res, "te0")$unit, "mmHg")
  raw <- function() {
    ini({ tcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      cl <- tcl * WT^0.75
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(raw, time = "h", central = "mg", Cc = "mg/L", WT = "kg")
  expect_match(.unitRow(res, "cl")$issue, "raises a quantity in kg to the power 0.75")
  mixed <- function() {
    ini({ tcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      cl <- tcl * (AGE + WT)
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(mixed, time = "h", central = "mg", Cc = "mg/L", WT = "kg", AGE = "year")
  expect_match(.unitRow(res, "cl")$issue, "adds year to kg, which are not convertible")
  # without declared units the covariates are compatible with anything
  res <- checkUnits(mixed, time = "h", central = "mg", Cc = "mg/L")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "tcl")$unit, "L/h")
})

test_that("transit absorption, bioavailability and residual structures map to their units", {
  m <- function() {
    ini({
      lka <- 0.1; lcl <- 1; lvc <- 1; lmtt <- 1; lfdepot <- 0.5
      propSd <- 0.1; addSd <- 1; pSd <- 0.1; pExp <- 1; sdE <- 0.1; lnSd <- 0.1; addSdE <- 1
    })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc); mtt <- exp(lmtt)
      f(depot) <- exp(lfdepot)
      d/dt(depot) <- transit(3, mtt) - ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc
      Cc ~ add(addSd) + prop(propSd)
      Cp <- Cc
      Cp ~ pow(pSd, pExp)
      E <- Cc / (1 + Cc)
      E ~ logitNorm(sdE, 0, 1)
      G <- Cc
      G ~ lnorm(lnSd)
      H <- E
      H ~ add(addSdE)
    })
  }
  res <- checkUnits(m, time = "h", depot = "mg", Cc = "mg/L", Cp = "mg/L", G = "mg/L", E = "unitless", H = "unitless")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "mtt")$unit, "h")
  expect_equal(.unitRow(res, "lmtt")$transformOf, "mtt")
  expect_equal(.unitRow(res, "f(depot)")$unit, "unitless")
  expect_equal(.unitRow(res, "lfdepot")$unit, "unitless")
  expect_equal(.unitRow(res, "lfdepot")$transformOf, "f(depot)")
  expect_equal(.unitRow(res, "addSd")$unit, "mg/L")
  expect_equal(.unitRow(res, "propSd")$unit, "unitless")
  expect_equal(.unitRow(res, "pExp")$unit, "unitless")
  expect_true(is.na(.unitRow(res, "pSd")$unit))
  expect_equal(.unitRow(res, "pSd")$source, "unresolved")
  expect_equal(.unitRow(res, "sdE")$unit, "unitless")
  expect_equal(.unitRow(res, "lnSd")$unit, "unitless")
  expect_equal(.unitRow(res, "G")$unit, "mg/L")
  expect_equal(.unitRow(res, "addSdE")$unit, "unitless")
  expect_equal(.unitRow(res, "E")$issue, NA_character_)
})

test_that("a dimensionless time unit is accepted", {
  # an exposure-response model with no time axis declares time = "none"
  er <- function() {
    units <- list(time = "none", dosing = "mg", concentration = "unitless")
    ini({ tint <- 0; tslope <- 0.1; propSd <- 0.1 })
    model({
      lp <- tint + tslope * AUC
      p <- expit(lp)
      p ~ prop(propSd)
    })
  }
  res <- checkUnits(er, AUC = "mg*h/L")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "time")$unit, "unitless")
  expect_equal(.unitRow(res, "p")$unit, "unitless")
  # the argument of expit() is not judged, so the linear predictor and its
  # coefficients stay unresolved rather than being forced unitless
  expect_equal(.unitRow(res, "lp")$source, "unresolved")
  expect_true(is.na(.unitRow(res, "tslope")$unit))
  expect_true(is.na(.unitRow(res, "tint")$unit))
  expect_equal(nlmixr2lib:::.unitToUdunits("L/unitless"), "L/1")
})

test_that("an undosed turnover state is unresolved until declared", {
  pd <- function() {
    ini({ tkin <- 1; tkout <- 1; addSd <- 1 })
    model({
      kin <- tkin; kout <- tkout
      d/dt(R) <- kin - kout * R
      E <- R
      E ~ add(addSd)
    })
  }
  res <- checkUnits(pd, time = "h")
  expect_true(all(is.na(res$issue)))
  expect_true(is.na(.unitRow(res, "R")$unit))
  expect_true(is.na(.unitRow(res, "kin")$unit))
  expect_equal(.unitRow(res, "kin")$source, "unresolved")
  res <- checkUnits(pd, time = "h", R = "mmHg")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "kin")$unit, "mmHg/h")
  expect_equal(.unitRow(res, "tkin")$unit, "mmHg/h")
  expect_equal(.unitRow(res, "kout")$unit, "1/h")
  expect_equal(.unitRow(res, "E")$unit, "mmHg")
  expect_equal(.unitRow(res, "addSd")$unit, "mmHg")
})

test_that("a linCmt() model keeps its output at the boundary and leaves parameters unresolved", {
  lc <- function() {
    ini({ lcl <- 1; lvc <- 1; lka <- 1; propSd <- 0.1 })
    model({
      cl <- exp(lcl); vc <- exp(lvc); ka <- exp(lka)
      linCmt() ~ prop(propSd)
    })
  }
  res <- checkUnits(lc, time = "h", depot = "mg", linCmt = "ng/mL")
  expect_true(all(is.na(res$issue)))
  expect_equal(.unitRow(res, "rxLinCmt")$unit, "ng/mL")
  expect_equal(.unitRow(res, "rxLinCmt")$source, "declared")
  expect_true(is.na(.unitRow(res, "cl")$unit))
  expect_equal(.unitRow(res, "propSd")$unit, "unitless")
})

test_that("legacy metadata is the default and arguments override it", {
  legacy <- function() {
    units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
    ini({ lka <- 0.1; lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc * 1000
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(legacy)
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  expect_equal(.unitRow(res, "depot")$unit, "mg")
  expect_equal(.unitRow(res, "Cc")$unit, "ng/mL")
  expect_equal(.unitRow(res, "time")$unit, "h")
  res <- checkUnits(legacy, Cc = "ug/mL")
  expect_equal(.unitRow(res, "Cc")$unit, "ug/mL")
  expect_equal(.unitRow(res, "Cc")$conversion, "mg/L * 1 = ug/mL")
  placeholder <- function() {
    units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
    ini({ lka <- 0.1; lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(placeholder)
  expect_true(all(is.na(res$issue)))
  expect_true(is.na(.unitRow(res, "time")$unit))
  expect_equal(length(attr(res, "notes")), 3L)
  expect_match(attr(res, "notes")[1], "units\\$time = 'time_unit' is not a parseable unit")
  # a free-form metadata key that names nothing in the model is noted, not an error
  extraKey <- function() {
    units <- list(time = "h", dosing = "mg", concentration = "mg/L", weight = "kg")
    ini({ lka <- 0.1; lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl); vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(extraKey)
  expect_true(all(is.na(res$issue)))
  expect_equal(attr(res, "notes"), "metadata units$weight names nothing in the model and was ignored")
  expect_false("weight" %in% res$name)
})

test_that("covariate units come from covariateData and go back into it", {
  withCov <- function() {
    covariateData <- list(
      WT = list(description = "Body weight", units = "kg", type = "continuous"),
      SEX = list(description = "Sex", units = "(binary)", type = "binary")
    )
    ini({ tcl <- 1; lvc <- 1; e_age_cl <- 0.01; propSd <- 0.1 })
    model({
      cl <- tcl * (WT / 70)^0.75 * (1 + 0.2 * SEX) * (1 + e_age_cl * (AGE - 40))
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- checkUnits(withCov, time = "h", central = "mg", Cc = "mg/L")
  expect_equal(.unitRow(res, "WT")$unit, "kg")
  expect_equal(.unitRow(res, "WT")$source, "declared")
  expect_equal(.unitRow(res, "SEX")$source, "compatible")
  expect_equal(.unitRow(res, "AGE")$source, "compatible")
  out <- addUnits(withCov, time = "h", central = "mg", Cc = "mg/L", AGE = "year")
  expect_equal(out$meta$covariateData$WT$units, "kg")
  expect_equal(out$meta$covariateData$SEX$units, "(binary)")
  expect_equal(out$meta$covariateData$AGE, list(description = NA_character_, units = "year", type = "continuous"))
  expect_equal(out$meta$units$WT, "kg")
  expect_equal(out$meta$units$AGE, "year")
  expect_equal(out$meta$units$e_age_cl, "1/year")
  expect_null(out$meta$units$SEX)
})

test_that("bad input is rejected at the boundary", {
  expect_error(checkUnits(.unitFixture1cmt, time = "h", foo = "mg"), "foo", class = "nlmixr2libUnitError")
  expect_error(checkUnits(.unitFixture1cmt, time = "h", Cc = "mg per L"), "Cc", class = "nlmixr2libUnitError")
  expect_error(checkUnits(.unitFixture1cmt, time = "h", Cc = "kg^0.75"), "fractional", class = "nlmixr2libUnitError")
  expect_error(checkUnits(.unitFixture1cmt, "mg"), "must be named", class = "nlmixr2libUnitError")
  expect_error(checkUnits(.unitFixture1cmt, time = "h", time = "day"), "Duplicated", class = "nlmixr2libUnitError")
  expect_error(checkUnits(.unitFixture1cmt, units = list(1)), "names")
  res <- checkUnits(.unitFixture1cmt, time = "h", depot = "mg", Cc = "mg/L", cl = "mg")
  expect_match(.unitRow(res, "kel")$issue, "has units mg/L but kel is 1/h")
})

# ---- addUnits() --------------------------------------------------------------------

test_that("addUnits() leaves a consistent model's solution unchanged and writes its units", {
  ev <- rxode2::et(amt = 100, cmt = "depot") |> rxode2::et(seq(0, 24, by = 4))
  before <- rxode2::rxSolve(.unitFixture1cmtFixed, ev, params = c(WT = 70), returnType = "data.frame")
  out <- addUnits(.unitFixture1cmtFixed, time = "h", depot = "mg", Cc = "mg/L")
  expect_s3_class(out, "rxUi")
  after <- rxode2::rxSolve(out, ev, params = c(WT = 70), returnType = "data.frame")
  expect_equal(after$Cc, before$Cc)
  expect_equal(
    out$meta$units,
    list(
      time = "h",
      depot = "mg",
      central = "mg",
      Cc = "mg/L",
      lka = "unitless",
      lcl = "unitless",
      lvc = "unitless",
      propSd = "unitless",
      addSd = "mg/L",
      ka = "1/h",
      cl = "L/h",
      vc = "L",
      kel = "1/h"
    )
  )
  expect_equal(out$meta$dosing, "depot")
  expect_null(out$meta$unitConversions)
  withEta <- addUnits(.unitFixture1cmt, time = "h", depot = "mg", Cc = "mg/L")
  expect_equal(withEta$meta$units$etalcl, "unitless")
})

test_that("addUnits() inserts the factor the arithmetic needs and records it", {
  ev <- rxode2::et(amt = 100, cmt = "depot") |> rxode2::et(seq(0, 24, by = 4))
  before <- rxode2::rxSolve(.unitFixture1cmtFixed, ev, params = c(WT = 70), returnType = "data.frame")
  expect_message(
    out <- addUnits(.unitFixture1cmtFixed, time = "h", depot = "mg", Cc = "ng/mL"),
    "Inserted 1 unit conversion"
  )
  after <- rxode2::rxSolve(out, ev, params = c(WT = 70), returnType = "data.frame")
  # the ODE system is unchanged; the tolerance covers the solver's own noise
  expect_equal(after$Cc, before$Cc * 1000, tolerance = 1e-5)
  expect_equal(rxode2::modelExtract(out, "Cc"), "Cc <- central/vc * 1000")
  expect_equal(out$meta$unitConversions, list(Cc = "mg/L * 1000 = ng/mL"))
  expect_equal(out$meta$units$Cc, "ng/mL")
  expect_equal(out$meta$units$addSd, "ng/mL")
  # idempotent: the second pass finds the conversion already in place
  again <- checkUnits(out)
  expect_true(all(is.na(again$conversion)))
  expect_true(all(is.na(again$issue)))
  out2 <- addUnits(out)
  expect_equal(rxode2::modelExtract(out2, "Cc"), "Cc <- central/vc * 1000")
  expect_equal(out2$meta$unitConversions, list(Cc = "mg/L * 1000 = ng/mL"))
})

test_that("addUnits() applies a conversion inside an if branch", {
  branchy <- function() {
    ini({ lcl <- 1; lvc <- 1; propSd <- 0.1 })
    model({
      cl <- exp(lcl); vc <- exp(lvc)
      if (WT > 70) { kel <- cl / vc } else { kel <- cl / vc }
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  out <- addUnits(branchy, time = "h", central = "mg", Cc = "mg/L", cl = "mL/min")
  text <- paste(rxode2::modelExtract(out), collapse = "\n")
  expect_equal(lengths(regmatches(text, gregexpr("kel <- cl/vc * 0.06", text, fixed = TRUE))), 2L)
  expect_equal(out$meta$unitConversions, list(kel = "mL/L/min * 0.06 = 1/h"))
  ev <- rxode2::et(amt = 100, cmt = "central") |> rxode2::et(c(0, 4))
  before <- rxode2::rxSolve(branchy, ev, params = c(WT = 80), returnType = "data.frame")
  after <- rxode2::rxSolve(out, ev, params = c(WT = 80), returnType = "data.frame")
  expect_equal(after$kel, before$kel * 0.06)
})

test_that("addUnits() metadata survives further piping and a legacy library model converts", {
  out <- addUnits(.unitFixture1cmtFixed, time = "h", depot = "mg", Cc = "ng/mL")
  piped <- rxode2::model(out, Cc <- central / vc * 1000)
  expect_equal(piped$meta$units$Cc, "ng/mL")
  expect_equal(piped$meta$dosing, "depot")
  expect_equal(piped$meta$unitConversions, list(Cc = "mg/L * 1000 = ng/mL"))
  van <- rxode2::rxode2(readModelDb("Zhao_2014_vancomycin"))
  res <- checkUnits(van)
  expect_true(all(is.na(res$issue)))
  expect_true(all(is.na(res$conversion)))
  expect_equal(.unitRow(res, "vc")$unit, "L")
  expect_equal(.unitRow(res, "cl")$unit, "L/h")
  expect_equal(.unitRow(res, "WT")$unit, "kg")
  out <- addUnits(van)
  expect_equal(out$meta$units$central, "mg")
  expect_equal(out$meta$units$Cc, "mg/L")
  expect_null(out$meta$units$dosing)
  expect_equal(out$meta$dosing, "central")
  expect_equal(out$meta$description, van$meta$description)
})

# nolint end
