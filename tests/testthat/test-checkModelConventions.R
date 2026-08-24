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

test_that("every blessed chain prefix accepts numbered members and the register documents the same regex", {
  # `compartmentRegex` is the ONLY thing that clears the "Compartment 'x1' is
  # not a canonical name" warning for a numbered chain -- a compartment-names.md
  # entry alone does nothing. Ratifying a chain canonical is therefore two edits
  # (the regex in R, the documentation in the register) and they silently drift.
  # This test enumerates the prefixes out of the regex itself, so a newly
  # blessed prefix is covered the moment it is added, and it pins the R constant
  # to the string the register advertises.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  prefixes <- strsplit(
    sub("\\)\\[0-9\\]\\+\\$$", "", sub("^\\^\\(", "", conv$compartmentRegex)),
    "|", fixed = TRUE
  )[[1]]
  expect_true(length(prefixes) >= 9L)
  expect_true(all(grepl("^[a-z_]+$", prefixes)))
  # Every ratified chain prefix accepts its numbered members ...
  for (p in prefixes) {
    expect_true(nlmixr2lib:::.matchesCompartment(paste0(p, "1"), conv), info = p)
    expect_true(nlmixr2lib:::.matchesCompartment(paste0(p, "12"), conv), info = p)
  }
  # ... and an unblessed stem is still rejected, so the check can go red.
  expect_false(nlmixr2lib:::.matchesCompartment("necrotic1", conv))
  expect_false(nlmixr2lib:::.matchesCompartment("granuloma3", conv))

  # The register's documented regex must be byte-identical to the R constant.
  md <- readLines(nlmixr2lib:::.compartmentNamesPath(), warn = FALSE)
  bullet <- grep("^- `compartmentRegex = ", md, value = TRUE)
  expect_length(bullet, 1L)
  documented <- sub("^- `compartmentRegex = \"(.*?)\"`.*$", "\\1", bullet,
                    perl = TRUE)
  expect_equal(documented, conv$compartmentRegex)
  # Each prefix is also named in that bullet's prose, so the reader-facing
  # description cannot fall behind the pattern.
  for (p in prefixes) {
    expect_true(grepl(paste0("`", p, "1`"), bullet, fixed = TRUE), info = p)
  }
})

test_that("the caseum<n> TB granuloma chain is accepted", {
  # Ratified 2026-08-24 with Karakitsios_2025_bedaquiline_mouse_pbpk.R; before
  # ratification the mouse and human models needed the
  # `paper_specific_compartments` escape hatch for these six states.
  good <- function() {
    description <- "Caseous-granuloma catenary chain"
    reference <- "R"
    units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
    ini({
      lcl <- 1; label("CL (L/h)")
      lvc <- 1; label("Vc (L)")
      lkcas <- log(0.01); label("Caseum diffusion rate (1/h)")
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl)
      vc <- exp(lvc)
      kcas <- exp(lkcas)
      d/dt(central) <- -cl / vc * central
      d/dt(lesion) <- kcas * (central / vc - lesion)
      d/dt(caseum1) <- kcas * (lesion - caseum1) - kcas * caseum1 + kcas * caseum2
      d/dt(caseum2) <- kcas * caseum1 - kcas * caseum2
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


test_that("the extraction skill never teaches a name the rules reject", {
  # The checks below only fire once a model exists. This one fires on the
  # INSTRUCTIONS, so a retired name cannot survive in the guidance that new
  # models are written from -- which is how `allo_cl`, `logitfr`, `Km` and
  # `Vmax` were still being taught after the models themselves were migrated.
  skillDir <- testthat::test_path("..", "..", ".claude", "skills",
                                  "extract-literature-model")
  skip_if(!dir.exists(skillDir), "extraction skill not present (installed package)")
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  docs <- list.files(skillDir, pattern = "[.]md$", recursive = TRUE, full.names = TRUE)
  expect_true(length(docs) > 0L)
  offenders <- character(0)
  for (f in docs) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (nm in names(conv$renamedParameters)) {
      if (grepl(paste0("(^|[^A-Za-z0-9._])", nm, "([^A-Za-z0-9._]|$)"), src)) {
        offenders <- c(offenders, paste0(basename(f), ": ", nm))
      }
    }
  }
  expect_equal(sort(unique(offenders)), character(0))
})

test_that("the extraction skill never shows fixed() repeated in a label", {
  skillDir <- testthat::test_path("..", "..", ".claude", "skills",
                                  "extract-literature-model")
  skip_if(!dir.exists(skillDir), "extraction skill not present (installed package)")
  pat <- paste0("<-\\s*fixed\\(.*\\)\\s*;\\s*label\\(\"[^\"]*",
                "\\bfixed(?![- ](?:effects?|dose|dosing))\\b")
  offenders <- character(0)
  for (f in list.files(skillDir, pattern = "[.]md$", recursive = TRUE, full.names = TRUE)) {
    src <- readLines(f, warn = FALSE)
    # skip the "Not this" column of the guidance table, which shows the
    # anti-pattern on purpose
    src <- src[!grepl("Not this|^\\|", src)]
    if (any(grepl(pat, src, perl = TRUE, ignore.case = TRUE))) {
      offenders <- c(offenders, basename(f))
    }
  }
  expect_equal(sort(unique(offenders)), character(0))
})


# ---- Unit spellings ----------------------------------------------------------

test_that("non-canonical time and dose spellings are flagged", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  chk <- function(tu, du) {
    ui <- list(meta = list(units = list(time = tu, dosing = du, concentration = "mg/L")))
    nlmixr2lib:::.checkUnitSpellings(ui, conv)
  }
  expect_equal(nrow(chk("hour", "mg")), 1L)
  expect_equal(nrow(chk("hr", "mg")), 1L)
  expect_equal(nrow(chk("h", "microgram")), 1L)
  expect_equal(chk("hour", "mg")$severity, "error")
})

test_that("canonical spellings and distinct units are accepted", {
  # "min" and "h" are BOTH canonical. Normalising between them would misstate
  # every value in the model, so the check must never conflate them.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  chk <- function(tu, du) {
    ui <- list(meta = list(units = list(time = tu, dosing = du, concentration = "mg/L")))
    nlmixr2lib:::.checkUnitSpellings(ui, conv)
  }
  for (tu in c("h", "min", "day", "week", "month", "year", "s")) {
    expect_equal(nrow(chk(tu, "mg")), 0L, info = tu)
  }
  for (du in c("mg", "ug", "ng", "g", "nmol", "umol", "IU")) {
    expect_equal(nrow(chk("h", du)), 0L, info = du)
  }
})

test_that("generic models' placeholder units are exempt", {
  # PK_1cmt and friends are dimensionless by design; Beal_2001_iv1cmt_bql
  # expresses time in half-lives. Those are correct, not unnormalised.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  chk <- function(tu, du) {
    ui <- list(meta = list(units = list(time = tu, dosing = du, concentration = "c")))
    nlmixr2lib:::.checkUnitSpellings(ui, conv)
  }
  expect_equal(nrow(chk("time_unit", "dose_unit")), 0L)
  expect_equal(nrow(chk("half_life", "dose_unit")), 0L)
})

test_that("no model in the database uses a non-canonical unit spelling", {
  # Enumerating: a newly added model with time = "hour" fails here.
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  bad <- character(0)
  for (f in list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) {
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (fld in c("time", "dosing")) {
      map <- if (fld == "time") conv$timeUnitSpellings else conv$doseUnitSpellings
      m <- regmatches(src, regexpr(sprintf('%s\\s*=\\s*"[^"]*"', fld), src))
      if (!length(m)) next
      val <- sub('.*"([^"]*)".*', "\\1", m)
      if (tolower(val) %in% names(map)) bad <- c(bad, paste0(basename(f), ": ", fld, "=", val))
    }
  }
  expect_equal(sort(unique(bad)), character(0))
})

test_that("no ini() line carries a quoted trailing comment rxode2 would promote into a broken label", {
  # Enumerating: a newly added model whose ini() line ends in a comment that
  # quotes the source paper fails here.
  #
  # rxode2 promotes a trailing comment on an ini() line into `label("<comment>")`
  # when the function is parsed with source refs kept. A double quote inside
  # that comment terminates the generated string early and the re-parse fails
  # inside rxode2:::.rxReplaceCommentWithLabel(), with no indication of which
  # model or line is responsible. ANY embedded double quote breaks it, not just
  # an unbalanced one, because the generated label is itself double-quoted.
  #
  # This is invisible to buildModelDb(), which resolves without source refs and
  # so never runs the promotion -- the build goes green while the test suite
  # fails. Moein_2024_apitolisib_human.R hit this on an eta line quoting the
  # Online Resource text. Fix by using single quotes inside the comment.
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  bad <- character(0)
  for (f in list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)) {
    lines <- readLines(f, warn = FALSE)
    inIni <- FALSE
    for (i in seq_along(lines)) {
      ln <- lines[[i]]
      if (grepl("^\\s*ini\\(\\{", ln)) {
        inIni <- TRUE
        next
      }
      if (grepl("^\\s*model\\(\\{", ln)) inIni <- FALSE
      if (!inIni) next
      # A standalone comment line is not attached to a parameter, and a line
      # that already calls label() is not promoted.
      if (grepl("^\\s*#", ln)) next
      if (grepl("label(", ln, fixed = TRUE)) next
      # Locate the first # that is not itself inside a string literal.
      chars <- strsplit(ln, "", fixed = TRUE)[[1]]
      nq <- 0L
      hash <- 0L
      for (k in seq_along(chars)) {
        if (chars[[k]] == '"') {
          nq <- nq + 1L
        } else if (chars[[k]] == "#" && nq %% 2L == 0L) {
          hash <- k
          break
        }
      }
      if (hash == 0L) next
      if (grepl('"', substring(ln, hash), fixed = TRUE)) {
        bad <- c(bad, paste0(basename(f), ":", i))
      }
    }
  }
  expect_equal(sort(bad), character(0))
})

# Within-subject random-effect levels. `etaiov_<param>_<occ>` (inter-occasion)
# and `etabvv_<param>_<visit>` (between-visit) both pair with the log-scale
# fixed effect `l<param>` rather than with a same-named fixed effect, so
# .isIOVEtaSuffix must accept both prefixes. Abdelgawad_2024_linezolid fits
# five-occasion IOV on ka/mtt alongside a two-visit BVV on vmax; before the
# bvv_ prefix was recognised, every etabvv_ slot raised a spurious "no
# matching fixed-effect parameter" warning.
test_that("iov_ and bvv_ level etas pair with their l<param> fixed effect", {
  levels_model <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "h", dosing = "mg", concentration = "mg/L")
    compartmentData <- list(
      central = list(analyte = "drug", units = "mg", specimen = "plasma", verified = TRUE)
    )
    covariateData <- list(
      OCC = list(
        description = "Occasion index", units = "(count)", type = "categorical",
        reference_category = NULL, notes = "n", source_name = "OCC"
      )
    )
    ini({
      lvmax <- 5; label("Maximum elimination rate (Vmax, mg/h)")
      lkm <- 3;   label("Michaelis-Menten constant (km, mg/L)")
      lvc <- 1;   label("Central volume (Vc, L)")
      etalvmax ~ 0.01
      etabvv_vmax_1 ~ 0.04
      etabvv_vmax_2 ~ fix(0.04)
      etaiov_vc_1 ~ 0.09
      etaiov_vc_2 ~ fix(0.09)
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      oc1 <- (OCC == 1)
      oc2 <- (OCC == 2)
      bvv <- oc1 * etabvv_vmax_1 + oc2 * etabvv_vmax_2
      iov <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
      vmax <- exp(lvmax + etalvmax + bvv)
      km <- exp(lkm)
      vc <- exp(lvc + iov)
      Cc <- central / vc
      d/dt(central) <- -vmax * Cc / (km + Cc)
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(levels_model, verbose = FALSE))
  naming <- res[res$category == "parameter_naming", ]
  expect_equal(nrow(naming), 0)
})

test_that("a bvv_ eta whose parameter has no fixed effect is still flagged", {
  bad_bvv <- function() {
    description <- "A"
    reference <- "R"
    units <- list(time = "h", dosing = "mg", concentration = "mg/L")
    compartmentData <- list(
      central = list(analyte = "drug", units = "mg", specimen = "plasma", verified = TRUE)
    )
    ini({
      lcl <- 1; label("Clearance (CL, L/h)")
      lvc <- 1; label("Central volume (Vc, L)")
      etabvv_nosuch_1 ~ 0.04
      propSd <- 0.1; label("Proportional residual error (fraction)")
    })
    model({
      cl <- exp(lcl + etabvv_nosuch_1)
      vc <- exp(lvc)
      kel <- cl / vc
      d/dt(central) <- -kel * central
      Cc <- central / vc
      Cc ~ prop(propSd)
    })
  }
  res <- suppressWarnings(checkModelConventions(bad_bvv, verbose = FALSE))
  expect_true(any(grepl("etabvv_nosuch_1", res$name, fixed = TRUE)))
})


# Reference-register duplicates. The registers resolve in document order,
# last one wins, so two blocks sharing a name AND a Type silently discard the
# earlier one. `col`, `mic` and `cloca` each hit this. A name repeated under
# DIFFERENT Types is deliberate and must not be flagged.

test_that("no register entry is duplicated at the same Type", {
  # Enumerating: a merge that leaves two same-Type blocks fails here.
  root <- system.file("references", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "references not installed")
  issues <- nlmixr2lib:::.referenceDuplicateIssues(
    Sys.glob(file.path(root, "*.md")))
  expect_equal(nrow(issues), 0L)
})

test_that("two blocks sharing a name and a Type are an error", {
  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "### foo (**canonical foo**)",
    "- **Type:** compartment",
    "- **Role:** first.",
    "",
    "### foo (**canonical foo again**)",
    "- **Type:** compartment",
    "- **Role:** second, silently discarded."
  ), tmp)
  issues <- nlmixr2lib:::.referenceDuplicateIssues(tmp)
  expect_equal(nrow(issues), 1L)
  expect_equal(issues$severity, "error")
  expect_equal(issues$category, "reference_duplicate")
  expect_equal(issues$name, "foo")
  expect_true(grepl("compartment", issues$message, fixed = TRUE))
})

test_that("the same name under two different Types is NOT flagged", {
  # col / complex / dap / lzd / mer / mero / plasma / van are each both a bare
  # compartment and a metabolite-suffix, on purpose.
  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "### col (**canonical colistin bare drug-state compartment**)",
    "- **Type:** compartment",
    "- **Role:** bare state.",
    "",
    "### col (**canonical colistin metabolite suffix**)",
    "- **Type:** metabolite-suffix",
    "- **Role:** suffix."
  ), tmp)
  expect_equal(nrow(nlmixr2lib:::.referenceDuplicateIssues(tmp)), 0L)
})

test_that("register blocks are parsed with and without the bold parenthetical", {
  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "## A section heading is not an entry",
    "### bare",
    "- **Type:** paper-named-param",
    "",
    "### withParen (**canonical thing**)",
    "- **Type:** compartment",
    "",
    "### Files using ALB",
    "(prose heading, no Type line)"
  ), tmp)
  b <- nlmixr2lib:::.referenceRegisterBlocks(tmp)
  expect_true(all(c("bare", "withParen") %in% b$name))
  expect_equal(b$type[b$name == "bare"], "paper-named-param")
  expect_equal(b$type[b$name == "withParen"], "compartment")
  # A prose heading has no Type; it must not collide with other prose headings
  # on the first word alone.
  expect_false("Files" %in% b$name)
})

# nolint end


# Fraction-metabolised pathway family (`fm_<pathway>`). `fm` was registered
# but its pathway suffixes were not, so nothing consumed the fact that
# `fm_other` / `fm_others` and `fm_h4` / `fm_H4` were each in the tree twice
# under two spellings. The register heading is the enumeration; these tests
# check that the heading is actually load-bearing.

# `.checkFmFamily()` needs only the bound names, so stub the two fields it
# reads rather than paying for a full nlmixr2() parse per case.
.fmStubUi <- function(ini = character(0), model = character(0)) {
  list(
    iniDf = if (length(ini)) {
      data.frame(name = ini, stringsAsFactors = FALSE)
    } else {
      NULL
    },
    lstExpr = lapply(model, function(s) parse(text = s)[[1]])
  )
}

test_that("every registered fm_<pathway> member is accepted", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  registered <- grep("^fm_", conv$paperNamedParams, value = TRUE)
  # Guard against a vacuous gate: if a register edit broke the heading, the
  # parse would yield nothing and every test below would pass for free.
  expect_gte(length(registered), 7L)
  expect_true(all(c("fm_cyp3a4", "fm_cyp3a5", "fm_h4", "fm_ko516_frac",
                    "fm_m3g", "fm_m6g", "fm_other") %in% registered))
  issues <- nlmixr2lib:::.checkFmFamily(.fmStubUi(ini = registered), conv)
  expect_equal(nrow(issues), 0L)
})

test_that("the two spellings this family exists to stop are not registered", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  registered <- grep("^fm_", conv$paperNamedParams, value = TRUE)
  # Plural residual route, and the capitalised clopidogrel H4 pathway.
  expect_false("fm_others" %in% registered)
  expect_false("fm_H4" %in% registered)
})

test_that("a plural residual route is an error naming the singular", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  issues <- nlmixr2lib:::.checkFmFamily(.fmStubUi(ini = "fm_others"), conv)
  expect_equal(nrow(issues), 1L)
  expect_equal(issues$severity, "error")
  expect_equal(issues$category, "fm_family")
  expect_equal(issues$name, "fm_others")
  expect_true(grepl("Rename to 'fm_other'", issues$suggestion, fixed = TRUE))
})

test_that("a capitalised pathway token is an error naming the lowercase form", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  issues <- nlmixr2lib:::.checkFmFamily(.fmStubUi(ini = "fm_H4"), conv)
  expect_equal(nrow(issues), 1L)
  expect_equal(issues$severity, "error")
  expect_equal(issues$name, "fm_H4")
  expect_true(grepl("Rename to 'fm_h4'", issues$suggestion, fixed = TRUE))
})

test_that("an unregistered pathway is an error pointing at the register", {
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  issues <- nlmixr2lib:::.checkFmFamily(.fmStubUi(ini = "fm_ugt1a1"), conv)
  expect_equal(nrow(issues), 1L)
  expect_equal(issues$severity, "error")
  expect_equal(issues$name, "fm_ugt1a1")
  expect_true(grepl("parameter-names.md", issues$suggestion, fixed = TRUE))
  # No near-miss exists, so it must not invent one.
  expect_false(grepl("Rename to", issues$suggestion, fixed = TRUE))
})

test_that("model() locals are in scope, not just ini() parameters", {
  # Mitra_2026_ziftomenib.R binds `fm_ko516_frac` inside model(), never in
  # ini(); an ini-only check would not see it.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  expect_equal(nrow(nlmixr2lib:::.checkFmFamily(
    .fmStubUi(model = "fm_ko516_frac <- 0.5"), conv)), 0L)
  bad <- nlmixr2lib:::.checkFmFamily(
    .fmStubUi(model = "fm_ko516frac <- 0.5"), conv)
  expect_equal(nrow(bad), 1L)
  expect_equal(bad$name, "fm_ko516frac")
})

test_that("a name bound inside an if() branch is still in scope", {
  # rxode2 permits branching in model(); a flat scan of the top-level
  # statements would let a respelling hide one level down.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  issues <- nlmixr2lib:::.checkFmFamily(
    .fmStubUi(model = "if (COV > 0) { fm_others <- 0.1 } else { fm_others <- 0.2 }"),
    conv)
  expect_equal(nrow(issues), 1L)
  expect_equal(issues$name, "fm_others")
  expect_true(grepl("Rename to 'fm_other'", issues$suggestion, fixed = TRUE))
})

test_that("source-paper notation in prose is out of scope", {
  # Both of these are correct as written and must never be flagged:
  #   * Svensson_2013_bedaquiline.R quotes the paper's `fm_M2` inside a
  #     label() string, beside the `CL_M2` / `V_M2` symbols it belongs with.
  #   * Robarge_2017_efavirenz.R's `fm_range` is FAT MASS in kg -- the
  #     canonical `FM` covariate -- sitting in population metadata next to
  #     `ffm_range`, not a fraction metabolised at all.
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  ui <- .fmStubUi(
    ini = "lcl_m2",
    model = c('lab <- "Apparent M2 metabolite clearance CL_M2/(F*fm_M2)"',
              'population <- list(ffm_range = "35.6-75.1", fm_range = "6.3-42.7")')
  )
  expect_equal(nrow(nlmixr2lib:::.checkFmFamily(ui, conv)), 0L)
})

test_that("no model in the database binds an unregistered fm_ parameter", {
  # Enumerating: a newly added model using a new or respelled pathway fails
  # here. The candidate set is derived from the tree rather than listed, so a
  # new model file is picked up without editing this test.
  root <- system.file("modeldb", package = "nlmixr2lib")
  skip_if(!nzchar(root) || !dir.exists(root), "modeldb sources not installed")
  conv <- nlmixr2lib:::.nlmixr2libConventions()
  files <- list.files(root, pattern = "[.]R$", recursive = TRUE,
                      full.names = TRUE)
  candidates <- Filter(function(f) {
    any(grepl("\\bfm_[A-Za-z0-9_]+", readLines(f, warn = FALSE), perl = TRUE))
  }, files)
  # If this drops to zero the test has stopped testing anything.
  expect_gte(length(candidates), 1L)
  bad <- character(0)
  for (f in candidates) {
    nm <- sub("[.]R$", "", basename(f))
    ui <- suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr(readModelDb(nm))))
    issues <- nlmixr2lib:::.checkFmFamily(ui, conv)
    if (nrow(issues)) bad <- c(bad, paste0(nm, ": ", issues$name))
  }
  expect_equal(sort(bad), character(0))
})
