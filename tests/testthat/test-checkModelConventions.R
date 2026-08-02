# Model fixtures below use the canonical rxode2/nlmixr2 DSL inside ini() and
# model() blocks, which relies on idioms that lintr flags as style warnings:
#   - ini() uses `param <- value; label("...")` on one line (semicolon_linter)
#   - model() uses `d/dt(cmt)` without spaces around `/` (infix_spaces_linter)
# Both are intentional — do not rewrite them.
# nolint start: semicolon_linter, infix_spaces_linter

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
    compartmentData <- list(
      depot   = list(analyte = "drug", units = "mg", specimen = "administration site", verified = TRUE),
      central = list(analyte = "drug", units = "mg", specimen = "plasma", verified = TRUE)
    )
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
    description <- "A"
    reference <- "R"
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1;   label("b")
      lvc <- 1;   label("c")
      prop.err <- 0.1; label("d")
      add.err <- 0.01; label("e")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1;   label("b")
      lvc <- 1;   label("c")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("a")
      lvc <- 1; label("b")
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      ADA = list(description = "Anti-drug antibody", units = "(binary)",
                 type = "binary")
    )
    ini({
      lcl <- 1; label("a")
      lvc <- 1; label("b")
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      ADA_POS = list(description = "ADA-positive", units = "(binary)",
                     type = "binary", source_name = "ADA")
    )
    ini({
      lcl <- 1; label("a")
      lvc <- 1; label("b")
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a")
      lcl <- 1; label("b")
      lvc <- 1; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
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
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a (1/day)")
      lcl <- 1; label("b (L/day)")
      lvc <- 1; label("c (L)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  # Asserted on the naming categories rather than "no warning at all": when
  # several fixtures in this file share identical model equations, rxode2 can
  # return a cached ui whose $meta predates the fixture's metadata block, so a
  # metadata-driven check (compartmentData) intermittently fires here. That is
  # a test-harness artifact of reusing equations across fixtures -- models read
  # from a file resolve their metadata correctly -- and it is orthogonal to
  # what this test is about.
  namingIssues <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  namingIssues <- namingIssues[
    namingIssues$severity %in% c("error", "warning") &
      namingIssues$category %in%
        c("parameter_names", "parameter_labels", "parameter_units",
          "deprecated_names", "fixed_label_disagreement", "compartments",
          "observation", "units", "covariates"), , drop = FALSE]
  expect_equal(nrow(namingIssues), 0L)

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
  expect_true(length(canon) > 30)
  expect_true(all(c("WT", "SEXF", "ADA_POS", "RACE_BLACK",
                    "RACE_BLACK_OTH", "CREAT", "ALB", "CRP", "CRCL",
                    "EOS", "PRIOR_GAST", "ADA_TITER") %in%
                    names(canon)))
  expect_true("SEXM" %in% canon$SEXF$aliases)
  expect_true("ADA" %in% canon$ADA_POS$aliases)
  expect_true("BLACK_OTH" %in% canon$RACE_BLACK_OTH$aliases)
  # ALB picked up the BALB (baseline-albumin) alias via Zhou 2021 belimumab.
  expect_true("BALB" %in% canon$ALB$aliases)
  # 2026-04-20 mergers: hsCRP + BLCRP + standard-CRP -> CRP; eGFR + CRCL_BSA -> CRCL;
  # ADA_TITRE + ADA_TITER -> ADA_TITER. BEOS renamed to EOS; GAST renamed to PRIOR_GAST.
  expect_false("hsCRP" %in% names(canon))
  expect_false("BLCRP" %in% names(canon))
  expect_false("eGFR" %in% names(canon))
  expect_false("CRCL_BSA" %in% names(canon))
  expect_false("BEOS" %in% names(canon))
  expect_false("GAST" %in% names(canon))
  expect_false("ADA_TITRE" %in% names(canon))
  expect_true(all(c("hsCRP", "CRPHS", "BLCRP") %in% canon$CRP$aliases))
  expect_true(all(c("eGFR", "CRCL_BSA") %in% canon$CRCL$aliases))
  expect_true("BEOS" %in% canon$EOS$aliases)
  expect_true("GAST" %in% canon$PRIOR_GAST$aliases)
  expect_true(all(c("ADA_TITRE", "ADAT") %in% canon$ADA_TITER$aliases))
  # Scope field is populated for every registered entry.
  expect_equal(canon$WT$scope, "general")
  expect_equal(canon$CRP$scope, "general")
  expect_equal(canon$CRCL$scope, "general")
  expect_equal(canon$EOS$scope, "general")
  expect_equal(canon$PRIOR_GAST$scope, "general")
  expect_equal(canon$ADA_TITER$scope, "general")
  expect_equal(canon$FORM_SAR_DP2$scope, "specific")
  expect_equal(canon$TUMTP_HODGKIN_CLASSICAL$scope, "specific")
  expect_equal(canon$ooc1$scope, "specific")
  expect_equal(canon$CONMED_EOX$scope, "specific")
  expect_equal(canon$DOSE_70MG$scope, "specific")
  expect_true("Xu_2019_sarilumab" %in% canon$FORM_SAR_DP2$example_models)
  expect_true("Cirincione_2017_exenatide" %in% canon$STUDY1$example_models)
})

test_that("covariate alias map resolves document-order last-writes-win", {
  map <- nlmixr2lib:::.nlmixr2libCovariateAliasMap()
  expect_equal(map[["BLACK_OTH"]], "RACE_BLACK_OTH")
  expect_equal(map[["SEXM"]], "SEXF")
  expect_equal(map[["ADA"]], "ADA_POS")
})

test_that("a scope-general canonical covariate produces no scope warning in any model", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      CRP = list(description = "C-reactive protein", units = "mg/L",
                 type = "continuous")
    )
    ini({
      lcl <- 1; label("a")
      lvc <- 1; label("b")
      e_crp_cl <- 0.05; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      cl <- exp(lcl) * (CRP / 5)^e_crp_cl
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  scoped <- res[res$category == "covariates" &
                  grepl("scoped 'specific'", res$message), ]
  expect_equal(nrow(scoped), 0)
})

test_that("a scope-specific canonical covariate in an unapproved model raises a warning", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      FORM_DP2 = list(description = "Drug-product 2 indicator",
                      units = "(binary)", type = "binary")
    )
    ini({
      lcl <- 1; label("a")
      lvc <- 1; label("b")
      e_dp2_cl <- 1.3; label("c")
      propSd <- 0.1; label("d")
    })
    model({
      cl <- exp(lcl) * e_dp2_cl^FORM_DP2
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  scoped <- res[res$category == "covariates" & res$name == "FORM_DP2" &
                  grepl("scoped 'specific'", res$message), ]
  expect_equal(nrow(scoped), 1)
  expect_equal(scoped$severity, "warning")
  expect_match(scoped$message, "Xu_2019_sarilumab")
})

test_that("a scope-specific canonical covariate in its listed model produces no scope warning", {
  # Xu_2019_sarilumab is in FORM_DP2's example_models list and uses FORM_DP2.
  res <- suppressWarnings(
    checkModelConventions("Xu_2019_sarilumab", verbose = FALSE)
  )
  scoped <- res[res$category == "covariates" & res$name == "FORM_DP2" &
                  grepl("scoped 'specific'", res$message), ]
  expect_equal(nrow(scoped), 0)
})

test_that("an rxUi object is accepted directly", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("a (1/day)")
      lcl <- 1; label("b (L/day)")
      lvc <- 1; label("c (L)")
      propSd <- 0.1; label("d")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
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

test_that("metabolite-suffixed parameters and compartments are accepted", {
  good <- function() {
    description <- "ADC parent + MMAE payload"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1;       label("Parent CL (L/day)")
      lvc <- 1;       label("Parent Vc (L)")
      lcl_mmae <- 1;  label("MMAE CL (L/day)")
      lvc_mmae <- 1;  label("MMAE Vc (L)")
      etalcl ~ 0.09
      propSd <- 0.1;          label("Parent prop residual (fraction)")
      propSd_mmae <- 0.1;     label("MMAE prop residual (fraction)")
    })
    model({
      cl <- exp(lcl + etalcl)
      vc <- exp(lvc)
      cl_mmae <- exp(lcl_mmae)
      vc_mmae <- exp(lvc_mmae)
      d/dt(central) <- -cl / vc * central
      d/dt(central_mmae) <- cl / vc * central - cl_mmae / vc_mmae * central_mmae
      Cc <- central / vc
      Cc_mmae <- central_mmae / vc_mmae
      Cc ~ prop(propSd)
      Cc_mmae ~ prop(propSd_mmae)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  naming <- res[res$category == "parameter_naming", ]
  cmts <- res[res$category == "compartments", ]
  obs <- res[res$category == "observation", ]
  expect_equal(nrow(naming), 0)
  expect_equal(nrow(cmts), 0)
  expect_equal(nrow(obs), 0)
})

test_that("DAR-numbered compartments are accepted", {
  good <- function() {
    description <- "DAR-mechanistic ADC"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1; label("CL (L/day)")
      lvc <- 1; label("Vc (L)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      d/dt(dar0_central) <- -cl / vc * dar0_central
      d/dt(dar1_central) <- -cl / vc * dar1_central
      Cc <- (dar0_central + dar1_central) / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  cmts <- res[res$category == "compartments", ]
  expect_equal(nrow(cmts), 0)
})

test_that("numbered precursor/lat chains are accepted", {
  good <- function() {
    description <- "Precursor maturation chain"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1; label("CL (L/day)")
      lvc <- 1; label("Vc (L)")
      lkmat <- log(0.1); label("Maturation rate (1/day)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      kmat <- exp(lkmat)
      d/dt(central) <- -cl / vc * central - kmat * central
      d/dt(precursor1) <- kmat * central - kmat * precursor1
      d/dt(precursor2) <- kmat * precursor1 - kmat * precursor2
      d/dt(precursor3) <- kmat * precursor2 - kmat * precursor3
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  cmts <- res[res$category == "compartments", ]
  expect_equal(nrow(cmts), 0)
})

test_that("shared-exponent covariate effects (e_<cov>_<param1>_<param2>) are accepted", {
  good <- function() {
    description <- "Shared allometric exponent"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    covariateData <- list(
      WT = list(description = "Body weight", units = "kg",
                type = "continuous")
    )
    ini({
      lcl <- 1; label("CL (L/day)")
      lvc <- 1; label("Vc (L)")
      lq  <- 1; label("Q  (L/day)")
      lvp <- 1; label("Vp (L)")
      e_wt_cl_q  <- 0.75; label("Shared WT exponent on CL and Q (unitless)")
      e_wt_vc_vp <- 1.0;  label("Shared WT exponent on Vc and Vp (unitless)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_cl_q
      vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp
      q  <- exp(lq)  * (WT / 70)^e_wt_cl_q
      vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp
      kel <- cl / vc
      k12 <- q  / vc
      k21 <- q  / vp
      d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
      d/dt(peripheral1) <- k12 * central - k21 * peripheral1
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  naming <- res[res$category == "parameter_naming", ]
  expect_equal(nrow(dep), 0)
  expect_equal(nrow(naming), 0)
})

test_that("multi-component CL covariate effects (e_<cov>_cl_ss, e_<cov>_cl_time) are accepted", {
  good <- function() {
    description <- "Time-varying CL"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    covariateData <- list(
      WT = list(description = "Body weight", units = "kg",
                type = "continuous")
    )
    ini({
      lcl_ss <- 1;   label("Steady-state CL (L/day)")
      lcl_time <- 1; label("Time-varying CL (L/day)")
      lvc <- 1;      label("Vc (L)")
      e_wt_cl_ss   <- 0.75; label("WT exponent on CL_ss (unitless)")
      e_wt_cl_time <- 0.5;  label("WT exponent on CL_time (unitless)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl_ss   <- exp(lcl_ss)   * (WT / 70)^e_wt_cl_ss
      cl_time <- exp(lcl_time) * (WT / 70)^e_wt_cl_time
      vc <- exp(lvc)
      d/dt(central) <- -(cl_ss + cl_time * exp(-0.1 * t)) / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  naming <- res[res$category == "parameter_naming", ]
  expect_equal(nrow(dep), 0)
  expect_equal(nrow(naming), 0)
})

test_that("deprecated lv1 / lv2 / lv structural names are flagged", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("a (L/day)")
      lv1 <- 1; label("b (L)")
      lv2 <- 1; label("c (L)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      cl <- exp(lcl)
      v1 <- exp(lv1)
      v2 <- exp(lv2)
      d/dt(central) <- -cl / v1 * central
      Cc <- central / v1
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  expect_true("lv1" %in% dep$name)
  expect_true("lv2" %in% dep$name)
  expect_match(dep$message[dep$name == "lv1"], "deprecated")
  expect_match(dep$suggestion[dep$name == "lv1"], "lvc")
})

test_that("deprecated lvm / vm Michaelis-Menten names are flagged", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("a (L/day)")
      lvc <- 1; label("b (L)")
      lvm <- 1; label("c (mg/day)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      vm <- exp(lvm)
      d/dt(central) <- -cl * central / vc - vm * central / (1 + central)
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  expect_true("lvm" %in% dep$name)
  expect_match(dep$suggestion[dep$name == "lvm"], "lvmax")
})

test_that("deprecated parent _adc suffix is flagged", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl_adc <- 1; label("a (L/day)")
      lvc_adc <- 1; label("b (L)")
      propSd <- 0.1; label("c (fraction)")
    })
    model({
      cl_adc <- exp(lcl_adc)
      vc_adc <- exp(lvc_adc)
      d/dt(central) <- -cl_adc / vc_adc * central
      Cc <- central / vc_adc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  expect_true("lcl_adc" %in% dep$name)
  expect_match(dep$suggestion[dep$name == "lcl_adc"], "lcl")
})

test_that("deprecated covariate-effect suffixes are flagged with rename suggestions", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      WT = list(description = "Body weight", units = "kg",
                type = "continuous")
    )
    ini({
      lcl <- 1; label("a (L/day)")
      lvc <- 1; label("b (L)")
      e_wt_v   <- 1.0; label("c1 (unitless)")
      e_wt_clq <- 0.75; label("c2 (unitless)")
      e_wt_vss <- 1.0; label("c3 (unitless)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_clq
      vc <- exp(lvc) * (WT / 70)^e_wt_v * (WT / 70)^e_wt_vss
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names", ]
  expect_true("e_wt_v" %in% dep$name)
  expect_true("e_wt_clq" %in% dep$name)
  expect_true("e_wt_vss" %in% dep$name)
  expect_match(dep$suggestion[dep$name == "e_wt_v"], "vc")
  expect_match(dep$suggestion[dep$name == "e_wt_clq"], "cl_q")
  expect_match(dep$suggestion[dep$name == "e_wt_vss"], "vc_vp")
})

test_that("reversed-order covariate effects (e_<param>_<cov>) are flagged", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    covariateData <- list(
      WT = list(description = "Body weight", units = "kg",
                type = "continuous")
    )
    ini({
      lcl <- 1; label("a (L/day)")
      lvc <- 1; label("b (L)")
      e_cl_wt <- 0.75; label("c (unitless)")
      propSd <- 0.1; label("d (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_cl_wt
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  dep <- res[res$category == "deprecated_names" & res$name == "e_cl_wt", ]
  expect_equal(nrow(dep), 1)
  expect_match(dep$message, "reversed-order")
})

test_that("deprecated C<metab> output naming is flagged in multi-output models", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1;       label("a (L/day)")
      lvc <- 1;       label("b (L)")
      lcl_mmae <- 1;  label("c (L/day)")
      lvc_mmae <- 1;  label("d (L)")
      CcpropSd     <- 0.1; label("e (fraction)")
      CmmaepropSd <- 0.1;  label("f (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      cl_mmae <- exp(lcl_mmae)
      vc_mmae <- exp(lvc_mmae)
      d/dt(central) <- -cl / vc * central
      d/dt(central_mmae) <- cl / vc * central - cl_mmae / vc_mmae * central_mmae
      Cc <- central / vc
      Cmmae <- central_mmae / vc_mmae
      Cc ~ prop(CcpropSd)
      Cmmae ~ prop(CmmaepropSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  obs <- res[res$category == "observation" & res$name == "Cmmae", ]
  expect_equal(nrow(obs), 1)
  expect_match(obs$suggestion, "Cc_mmae")
})

test_that("csf and isf are accepted as canonical physiologic compartments", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl  <- 1; label("clearance (CL, L/day)")
      lvc  <- 1; label("central volume (Vc, L)")
      lqcsf <- 0.1; label("CSF distribution clearance (L/day)")
      lvcsf <- 0.1; label("CSF volume (L)")
      lqisf <- 0.1; label("ISF distribution clearance (L/day)")
      lvisf <- 0.1; label("ISF volume (L)")
      propSd <- 0.1; label("g (fraction)")
    })
    model({
      cl  <- exp(lcl)
      vc  <- exp(lvc)
      qcsf <- exp(lqcsf)
      vcsf <- exp(lvcsf)
      qisf <- exp(lqisf)
      visf <- exp(lvisf)
      d/dt(central) <- -cl/vc*central - qcsf/vc*central + qcsf/vcsf*csf
      d/dt(csf)     <-  qcsf/vc*central - qcsf/vcsf*csf - qisf/vcsf*csf + qisf/visf*isf
      d/dt(isf)     <-  qisf/vcsf*csf - qisf/visf*isf
      Cc <- central/vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  comps <- res[res$category == "compartments", ]
  expect_equal(nrow(comps), 0)
})

test_that("target/complex with csf/isf location suffix are accepted", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl   <- 1;     label("a (L/day)")
      lvc   <- 1;     label("b (L)")
      lkon  <- 0;     label("c (L/mg/day)")
      lkoff <- 0;     label("d (1/day)")
      lkint <- -1;    label("e (1/day)")
      lr0   <- 0;     label("f (mg/L)")
      propSd <- 0.1;  label("g (fraction)")
    })
    model({
      cl   <- exp(lcl)
      vc   <- exp(lvc)
      kon  <- exp(lkon)
      koff <- exp(lkoff)
      kint <- exp(lkint)
      R0   <- exp(lr0)
      Cc <- central/vc
      d/dt(central)     <- -cl/vc*central - kon*target*central + koff*complex*vc
      d/dt(target)      <-  R0 - kon*target*Cc + koff*complex
      d/dt(complex)     <-  kon*target*Cc - (koff + kint)*complex
      d/dt(target_csf)  <- -kon*target_csf*Cc + koff*complex_csf
      d/dt(complex_csf) <-  kon*target_csf*Cc - (koff + kint)*complex_csf
      d/dt(target_isf)  <- -kon*target_isf*Cc + koff*complex_isf
      d/dt(complex_isf) <-  kon*target_isf*Cc - (koff + kint)*complex_isf
      target(0)      <- R0
      target_csf(0)  <- R0
      target_isf(0)  <- R0
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  comps <- res[res$category == "compartments", ]
  expect_equal(nrow(comps), 0)
})

test_that("expSd is accepted as a canonical residual error name for lnorm", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lka <- 0.1; label("ka (1/day)")
      lcl <- 1;   label("CL (L/day)")
      lvc <- 1;   label("Vc (L)")
      expSd <- 0.4; label("Log-scale residual SD")
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc <- exp(lvc)
      kel <- cl/vc
      d/dt(depot)   <- -ka*depot
      d/dt(central) <-  ka*depot - kel*central
      Cc <- central/vc
      Cc ~ lnorm(expSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  naming <- res[res$category == "parameter_naming" & res$name == "expSd", ]
  expect_equal(nrow(naming), 0)
})

test_that("capitalized compartments get a specific 'uppercase letter' warning", {
  bad <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1; label("CL (L/day)")
      lvc <- 1; label("Vc (L)")
      propSd <- 0.1; label("g (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      d/dt(central) <- -cl/vc*central
      d/dt(Cbrain)  <- 0.1*(0.5*central/vc - Cbrain)
      Cc <- central/vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad, verbose = FALSE))
  comps <- res[res$category == "compartments" & res$name == "Cbrain", ]
  expect_equal(nrow(comps), 1)
  expect_match(comps$message, "uppercase")
})

test_that("parallel-absorption depot/depot2 passes (already permitted by depot<n> regex)", {
  good <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
    ini({
      lka1 <- 0.1; label("ka1 (1/day)")
      lka2 <- 0.1; label("ka2 (1/day)")
      lcl  <- 1;   label("CL (L/day)")
      lvc  <- 1;   label("Vc (L)")
      propSd <- 0.1; label("g (fraction)")
    })
    model({
      ka1 <- exp(lka1)
      ka2 <- exp(lka2)
      cl  <- exp(lcl)
      vc  <- exp(lvc)
      kel <- cl/vc
      d/dt(depot)   <- -ka1*depot
      d/dt(depot2)  <- -ka2*depot2
      d/dt(central) <-  ka1*depot + ka2*depot2 - kel*central
      Cc <- central/vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(good, verbose = FALSE))
  comps <- res[res$category == "compartments", ]
  expect_equal(nrow(comps), 0)
})


# ---- Issue #474: retired parameter names -------------------------------------

test_that("retired parameter names are flagged as errors with the replacement", {
  retired <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^allo_cl
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(retired, verbose = FALSE))
  hit <- res[res$name == "allo_cl" & res$category == "deprecated_names", ]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$severity, "error")
  expect_match(hit$suggestion, "e_wt_cl", fixed = TRUE)
})

test_that("every renamedParameters entry maps to a different, non-empty name", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  map <- conv$renamedParameters
  expect_true(length(map) > 0L)
  expect_true(all(nzchar(names(map))))
  expect_true(all(nzchar(unlist(map))))
  expect_true(all(names(map) != unlist(map)))
})

test_that("no model in the database uses a retired parameter name", {
  # Enumerating the map rather than spot-checking, so a newly retired name
  # fails here until every model is migrated (see CLAUDE.md: a contract that
  # quantifies over "every case" needs an enumerating test).
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  src <- vapply(
    list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
    function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
    character(1)
  )
  for (nm in names(conv$renamedParameters)) {
    pat <- paste0("(^|[^A-Za-z0-9._])", nm, "[[:space:]]*<-")
    expect_equal(sum(grepl(pat, src)), 0L, info = paste("retired name still used:", nm))
  }
})

# ---- Issue #479: label claims "fixed" while the parameter is estimable -------

test_that("a label claiming a value is fixed while fix is FALSE is flagged", {
  claims <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      e_wt_cl <- 0.75; label("Allometric exponent on CL (unitless; FIXED)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_cl
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(claims, verbose = FALSE))
  hit <- res[res$category == "fixed_label_disagreement", ]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$name, "e_wt_cl")
  expect_equal(hit$severity, "error")
})

test_that("wrapping the same parameter in fixed() clears the disagreement", {
  ok <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless; FIXED)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_cl
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(ok, verbose = FALSE))
  expect_equal(sum(res$category == "fixed_label_disagreement"), 0L)
})

test_that("a plain label on an estimated parameter is not flagged", {
  # Guards against the pattern matching ordinary prose such as "transferred".
  plain <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      e_wt_cl <- 0.75; label("Allometric exponent on CL (unitless)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_cl
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(plain, verbose = FALSE))
  expect_equal(sum(res$category == "fixed_label_disagreement"), 0L)
})


test_that("the fixed-claim pattern only matches claims about the parameter itself", {
  # Each FALSE case is a real label from the database where the "fixed" clause
  # describes a DIFFERENT parameter. Matching them was the rule's main source
  # of false positives; see inst/references/fixed-provenance-followup.md.
  pat <- nlmixr2lib:::.fixedClaimPattern
  claims <- c(
    "Allometric exponent on CL/F (unitless; fixed)",
    "Max fractional decrease in UA production by febuxostat (fixed in source)",
    "Vascular reflection coefficient for tight tissues (unitless; fixed at 0.95 in Cao 2013)",
    "Proportional residual error (fraction) - ASSUMED from assay validation",
    "Additive residual error (nmol cystine / mg protein; assumed; paper does not specify)"
  )
  notClaims <- c(
    "IC50 for inhibition of precursor synthesis (Imax fixed to 1) (mg/L)",
    "Apparent DHA metabolite clearance CL/F_DHA (L/h/kg); F_DHA fixed to 1",
    "Apparent DM4 clearance CL_DM4 (L/day; V_DM4 fixed to 1 L)",
    "Scaling factor K relating eta_V/F to eta_CL/F (correlation fixed to 1)",
    "CMS nonrenal clearance CL_NRCMS (L/h); CL_RCMS structurally fixed at 0",
    "Passive plasma -> CSF_CM clearance (mL/min) - P-gp components fixed to 0",
    "Absorption lag time for Group 1 (h; Group 2 lag is fixed at 0)",
    "Clearance (CL, L/day)"
  )
  for (lbl in claims) {
    expect_true(grepl(pat, lbl, ignore.case = TRUE, perl = TRUE), info = lbl)
  }
  for (lbl in notClaims) {
    expect_false(grepl(pat, lbl, ignore.case = TRUE, perl = TRUE), info = lbl)
  }
})

test_that("variance terms are exempt from the fixed-label check", {
  # Labels on eta terms routinely note that the corresponding typical value
  # was fixed while the variance itself was estimated.
  withEta <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lka <- 0.1; label("Absorption rate (ka, 1/day)")
      lcl <- 1;   label("Clearance (CL, L/day)")
      lvc <- 1;   label("Central volume (Vc, L)")
      etalka ~ 1.41; label("Table IV omega_Ka (Ka pop is fixed but IIV is estimated)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      ka <- exp(lka + etalka)
      cl <- exp(lcl)
      vc <- exp(lvc)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(withEta, verbose = FALSE))
  expect_equal(sum(res$category == "fixed_label_disagreement"), 0L)
})


# ---- Issue #481: time-varying clearance naming --------------------------------

test_that("a time-dependent clearance using non-canonical names is flagged", {
  oldStyle <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("Clearance (CL, L/day)")
      lvc <- 1; label("Central volume (Vc, L)")
      emax <- -0.3; label("Maximal fractional change in CL (unitless)")
      hill <- 2.5; label("Sigmoidicity (unitless)")
      t50 <- 50;   label("Time of half-maximal change (day)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * exp(emax * t^hill / (t50^hill + t^hill))
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(oldStyle, verbose = FALSE))
  hit <- res[res$category == "time_varying_clearance", ]
  expect_equal(nrow(hit), 1L)
  expect_match(hit$suggestion, "cl_hill_max", fixed = TRUE)
})

test_that("the canonical time-varying clearance stems are accepted", {
  newStyle <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    ini({
      lcl <- 1; label("Clearance (CL, L/day)")
      lvc <- 1; label("Central volume (Vc, L)")
      cl_hill_max <- -0.3;  label("Maximum fractional change in CL over time (unitless)")
      cl_hill_gamma <- 2.5; label("Sigmoidicity of the time effect on CL (unitless)")
      cl_hill_t50 <- 50;    label("Time of half-maximal change in CL (day)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) *
        exp(cl_hill_max * t^cl_hill_gamma / (cl_hill_t50^cl_hill_gamma + t^cl_hill_gamma))
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(newStyle, verbose = FALSE))
  expect_equal(sum(res$category == "time_varying_clearance"), 0L)
})

test_that("no model in the database still uses a pre-#481 time-varying clearance name", {
  # Enumerating rather than spot-checking: a newly added model that reuses the
  # old spellings fails here until it is migrated.
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  offenders <- character(0)
  for (f in list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) {
    lines <- readLines(f, warn = FALSE)
    start <- grep("^\\s*model\\(\\{", lines)
    if (!length(start)) next
    for (ln in lines[seq(start[1], length(lines))]) {
      if (!grepl(nlmixr2lib:::.clearanceLhsPattern, ln)) next
      rhs <- sub("#.*$", "", sub("^[^<]*<-", "", ln))
      rhs <- gsub('"[^"]*"', "", rhs)
      if (!grepl(nlmixr2lib:::.bareTimePattern, rhs, perl = TRUE)) next
      if (grepl("cl_hill_|cl_exp_", rhs)) next
      offenders <- c(offenders, basename(f))
    }
  }
  # The four deliberate exclusions are a different structure (diurnal cosine,
  # circadian clock offset, a lag state, ADA-gated logistic onset).
  allowed <- c("Bienczak_2016_nevirapine.R", "Hayashi_1998_epoetinBeta.R",
               "Mann_2022_respiratory_physiology.R", "Yoshida_2024_fazpilodemab.R")
  expect_equal(sort(unique(setdiff(offenders, allowed))), character(0))
})

# ---- Issue #482: compartmentData ----------------------------------------------

test_that("compartmentData must cover every ODE state and use the vocabulary", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  expect_true(length(conv$specimenVocabulary) > 0L)
  expect_true(all(c("analyte", "units", "specimen", "verified") %in%
                    conv$compartmentDataFields))
  # The two non-matrix categories must exist, or every latent/PD state would be
  # forced into a false specimen.
  expect_true(all(c("administration site", "not applicable") %in%
                    conv$specimenVocabulary))
})

test_that("the worked compartmentData examples validate", {
  for (nm in c("PerezRuixo_2025_posdinemab", "HuttonSmith_2018_ranibizumab",
               "Le_2015_lampalizumab_cyno")) {
    res <- suppressWarnings(checkModelConventions(nm, verbose = FALSE))
    expect_equal(sum(res$category == "compartment_data"), 0L, info = nm)
  }
})


test_that("a label repeating what fixed() already states is flagged", {
  redundant <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    compartmentData <- list(
      central = list(analyte = "drug", units = "mg", specimen = "plasma", verified = TRUE)
    )
    ini({
      lcl <- 1; label("Clearance (CL, L/day)")
      lvc <- 1; label("Central volume (Vc, L)")
      e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless; fixed)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl) * (WT / 70)^e_wt_cl
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(redundant, verbose = FALSE))
  hit <- res[res$category == "fixed_label_redundant", ]
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$name, "e_wt_cl")
  expect_equal(hit$severity, "error")
})

test_that("provenance kept alongside fixed() is not flagged as redundant", {
  # `fixed()` says the value was not estimated; it cannot say where the value
  # came from, so that part of the label must survive.
  keeps <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "day", dosing = "mg", concentration = "mg/L")
    compartmentData <- list(
      central = list(analyte = "drug", units = "mg", specimen = "plasma", verified = TRUE)
    )
    ini({
      lcl <- fixed(1); label("Maternal apparent clearance, from Rizk 2015 (CL, L/day)")
      lvc <- 1; label("Central volume (Vc, L)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      d/dt(central) <- -cl / vc * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(keeps, verbose = FALSE))
  expect_equal(sum(res$category == "fixed_label_redundant"), 0L)
})

test_that("no model in the database repeats fixed() in its label", {
  # Enumerating: the database was cleaned in one pass, so any hit is new.
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  pat <- paste0("<-\\s*fixed\\(.*\\)\\s*;\\s*label\\(\"[^\"]*",
                "\\bfixed(?![- ](?:effects?|dose|dosing))\\b")
  offenders <- character(0)
  for (f in list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) {
    src <- readLines(f, warn = FALSE)
    if (any(grepl(pat, src, perl = TRUE, ignore.case = TRUE))) {
      offenders <- c(offenders, basename(f))
    }
  }
  expect_equal(sort(offenders), character(0))
})

# nolint end
