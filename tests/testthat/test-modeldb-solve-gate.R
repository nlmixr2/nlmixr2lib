# A shipped model that simulates identically zero is the worst kind of defect
# this registry can carry: there is no error, no warning and no NA, so every
# downstream check that divides one simulated number by another -- an exposure
# ratio, a fold-change, an accumulation index -- returns 0/0 and is duly
# recorded as a match. Eleven models reached main that way, and the vignette of
# one of them reported four-decimal agreement with the source paper for
# quantities it had never actually computed.
#
# WHAT GOES WRONG
#
# rxode2's `rxSolve()` defaulted to `useLinCmt = TRUE` until the conversion was
# made opt-in in September 2026 (rxode2 issue 1389); the probes below pass the
# flag explicitly, so the gate is unaffected by the default. Under conversion
# `.odeToLinDetect()` pattern-matches a depot/central/peripheral ODE system and
# rewrites it as the analytical `linCmt()` solution. The rewrite keeps the
# compartment topology and the transfer-rate parameters and DISCARDS every
# right-hand-side term that is not proportional to a state -- a `transit()`
# absorption rate, a zero-order input, an endogenous production term. Those
# terms are invisible to both of the converter's guards:
# `.odeToLinDetectTopology()` skips each term for which `is.na(term$state)`,
# and `.odeToLinMassBalanced()` only balances state-proportional coefficients.
#
# So a model written as
#
#   d/dt(depot)   <- transit(nn, mtt) - ka * depot
#   d/dt(central) <- ka * depot - kel * central
#   f(depot)      <- 0
#   Cc            <- central / vc
#   Cc ~ add(addSd)
#
# is converted to `Cc <- linCmt(ka, kel, vc)` with `f(depot) <- 0` retained.
# The transit input is gone, the only remaining drug input is the depot bolus,
# and `f(depot) <- 0` zeroes that. Every concentration is exactly 0.
#
# The same conversion also has a quieter failure mode: without `f(depot) <- 0`
# the model still solves, but as plain first-order absorption with the transit
# chain deleted. That is a wrong model rather than a dead one, which is why the
# second test below compares the two solve paths instead of only looking for
# zeros.
#
# This is an rxode2 defect, not a model-file idiom -- `?odeToLin` documents
# that detection "requires linear ODE right-hand sides" and that conversion
# "preserv[es] all other model lines", and both claims are violated. It is
# recorded in refharvest/reports/maint-001-zero_solving_model_gate_and_sweep.md
# together with a minimal reproducible pair and the suggested upstream patch.
# Until rxode2 ships a fix the affected models are quarantined below.
#
# WHY A TEST AND NOT A RULE
#
# The idiom is not what distinguishes the broken models from the sound ones:
# of 46 models carrying `f(depot) <- 0` together with `transit()` or
# `podo()`/`tad()`, 11 were dead and 34 were fine, because only some of them
# match the converter's topology. No reviewer is going to run
# `rxode2:::.odeToLinDetect()` in their head. The gate has to be mechanical.

# EMPTIED 2026-09-16. This list quarantined models defeated by rxode2's
# ODE-to-linCmt conversion silently dropping exogenous input terms
# (https://github.com/nlmixr2/rxode2/issues/1370). The entry condition it set
# for its own removal was "the same commit that raises the rxode2 minimum in
# DESCRIPTION to the first release that refuses to convert an ODE carrying an
# exogenous input term" -- DESCRIPTION now requires rxode2 (>= 5.1.8), and that
# release fixes it.
#
# Verified rather than assumed: all eleven former entries were re-probed under
# rxode2 5.1.8 with this file's own probeSolveModel(), and every one responds
# to a dose again (doseEffect 5.4 to 130.4, none zero):
#
#   Bienczak_2016_efavirenz 73.0   Chigutsa_2011_rifampicin 53.9
#   Chigutsa_2012_ofloxacin 59.6   Dong_2014_mycophenolic_acid 130.4
#   Marques_2025_salbutamol 5.4    ResendizGalvan_2025_cycloserine 93.1
#   Sloan_2017_rifampicin 82.2     Smythe_2013_gatifloxacin 91.1
#   Tikiso_2021_abacavir 60.0      Vinnard_2017_rifampicin 73.3
#   Wilkins_2008_rifampicin 61.5
#
# No model file was rewritten to achieve this -- the fix was entirely upstream,
# which is why the models were left faithful to their papers rather than
# reshaped to dodge the converter. Keep the vector (empty) and the third test:
# they are what will catch the next model that stops responding to a dose.
knownLinCmtDropModels <- character(0)

# Models the probe cannot drive, SPLIT BY CAUSE. The test's own message says
# "either the model is broken or helper-solveProbe.R needs to learn how to
# drive it -- work out which before adding it". Both causes turned up on
# 2026-09-17, and they need different follow-up, so they get different lists.

# (1) The probe's limitation, NOT a model defect.
#
# Kim_2024_islatravir: the probe's covariate derivation pairs an oral tablet
# amount with the implant route. DOSE_ISLATRAVIR_MG takes 0.1 (the first
# `DOSE_ISLATRAVIR_MG - c` constant in model(), from `is_0p1`) while ROUTE_ORAL
# takes its reference category 0 (implant). The model sets kimpl and
# kce_implant to NA for an implant load outside its five studied loads (48, 52,
# 54, 56, 62 mg), so every ODE state is NaN. That is the model refusing an
# unsupported regimen, as designed. Its vignette solves the studied regimens.
#
# Hopkins_2024_amisulpride: the formulation is chosen by ten one-hot
# FORM_SEP4199_* indicators with no reference level, so the probe's reference
# value of 0 for each selects no formulation. Every bioavailability fraction is
# then 0 and the dose changes nothing. The vignette sets exactly one indicator
# per arm.
#
# Tan_2025_adomeglivant_ffh: the probe reads FPG = 3.5 from `FPG - 3.5`, which
# is the HOMA glucose floor in `bcf0 <- INS_BL / (ki * (FPG - 3.5))`, not a
# reference value. bcf0 is then infinite and Hba1c NaN from the first step. The
# model requires INS_BL < 5 * (FPG - 3.5), as its covariateData says; at the
# cohort reference FPG = 8.2 (the `FPG - 8.2` in kout) every state is finite.
probeUnsupportedModels <- c(
  "Kim_2024_islatravir",
  "Hopkins_2024_amisulpride",
  "Tan_2025_adomeglivant_ffh"
)

# (2) Genuinely BROKEN models, quarantined so this gate can go green while the
# defects are recorded rather than hidden. Both predate the gate, neither has a
# vignette -- which is precisely why nothing caught them. Do NOT "fix" either by
# guessing the algebra; both need their source reparameterisation checked.
#
#   indirect_1cpt_stim_kin_r0rmaxcrmax
#     IC50 <- crmax*(rbase*(1+emax)-rmax)/(rmax-rbase) evaluates to -0.9838 at
#     the shipped ini values (rbase 1.2214, rmax 2.4596, emax 0.5034, crmax
#     1.9542). A negative IC50 makes Cc + IC50 cross zero as Cc rises, so the
#     solver hits a singularity: "could not solve the system". The sibling
#     indirect_1cpt_stim_kout_r0rmaxcrmax uses the mirrored form
#     crmax*(rmax*(1+emax)-rbase)/(rbase-rmax) and solves, which is what
#     suggests a sign transposition here rather than a modelling choice.
#
#   indirect_prec_1cpt_stim_r0rmaxcrmax
#     Same negative-EC50 expression, AND line 42 defines kin in terms of
#     itself: kin <- (kout*(kin+kpout)*rbase)/kpin.
#
#   tgi_sat_genVonBertalanffy, tgi_no_sat_powerLaw
#     Both write `lgamma <- 0.95` and `gamma <- exp(lgamma)`, so the growth
#     exponent is exp(0.95) = 2.59, though both label gamma "proliferative
#     cells as a fraction of the full tumor volume". Every other lgamma in the
#     registry is written `log(x)`. With gamma > 1, `kg * tumor_size^gamma`
#     outgrows any linear loss and the tumour size reaches infinity in finite
#     time: NaN from t = 0.5 in genVonBertalanffy, "could not solve the
#     system" in powerLaw. Both solve at gamma = 0.95 or 2/3.
knownBrokenModels <- c(
  "indirect_1cpt_stim_kin_r0rmaxcrmax",
  "indirect_prec_1cpt_stim_r0rmaxcrmax",
  "tgi_sat_genVonBertalanffy",
  "tgi_no_sat_powerLaw"
)

# Solving every model in the registry takes roughly an hour. That is the right
# cost for a release gate and the wrong cost for an edit-run-edit loop, so the
# full sweep runs when NLMIXR2LIB_SOLVE_GATE=full and the at-risk subset runs
# otherwise. The subset is NOT a hand-written list: it is computed by
# linCmtRiskCandidates() below, which walks every model file in the registry,
# so a model added next week is screened automatically.
#
# STATING THE CAP RATHER THAN LETTING IT READ AS FULL COVERAGE. On the default
# path both gates below run over the screened candidates only. The screen's
# three preconditions are specific to this defect class -- it looks for an
# exogenous input term, which is the thing the conversion drops -- so what the
# default path gives up is (a) a model that solves to zero for some reason
# unrelated to the conversion, and (b) a model the conversion changes for some
# reason other than a dropped input term. Under `full` both gates widen to the
# entire registry, which covers every model rxode2 would convert (2912 of 2912
# screened, 802 actually convertible as of 2026-09-13) rather than only the
# 289 the cheap screen flags.
.solveGateScope <- function() {
  if (identical(tolower(Sys.getenv("NLMIXR2LIB_SOLVE_GATE", "")), "full")) "full" else "screened"
}

# ---- Probing in bounded subprocesses ----------------------------------------
#
# Every probe used to run in the test process. On GitHub's ubuntu R-devel and
# R-oldrel-1 legs that process grew from 1.9 GB to 15.5 GB plus 2.7 GB of swap
# in the 22 minutes after this file started -- about 1 GB a minute -- and the
# hosted runner was shut down ("The runner has received a shutdown signal",
# exit 143) on every push to main and every pull request since the 2026-09-17
# consolidation. The same loop is flat on R 4.6.1 with the identical rxode2
# build (240 -> 330 MB over 40 models, warm or cold compile cache, 1 to 16
# solver threads) and in an R 4.5.3 container, so the growth is specific to
# the hosted runner. Probing each block of models in a fresh R process bounds
# the damage to one block whatever leaks, at the cost of one package load per
# block; setting NLMIXR2LIB_SOLVE_GATE_TRACE=1 prints each child's and the
# parent's memory after every block. Both solve paths are probed once per model and
# shared by the tests below, which is no more work than before.
.probeBlockSize <- 20L

.probeInChild <- function(nms) {
  dev <- pkgload::is_dev_package("nlmixr2lib")
  root <- if (dev) normalizePath(testthat::test_path("..", "..")) else NA_character_
  helper <- normalizePath(testthat::test_path("helper-solveProbe.R"))
  trace <- nzchar(Sys.getenv("NLMIXR2LIB_SOLVE_GATE_TRACE")) && .Platform$OS.type == "unix"
  res <- callr::r(
    function(nms, dev, root, helper, trace) {
      if (dev) {
        pkgload::load_all(root, quiet = TRUE)
      } else {
        library(nlmixr2lib)
      }
      probes <- new.env()
      sys.source(helper, envir = probes)
      db <- nlmixr2lib::modeldb
      out <- lapply(nms, function(n) {
        list(
          lin = probes$probeSolveModel(n, db, useLinCmt = TRUE),
          ode = probes$probeSolveModel(n, db, useLinCmt = FALSE)
        )
      })
      names(out) <- nms
      if (trace) {
        rss <- as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
        message(sprintf(
          "solve-gate block of %d models: child RSS %.0f MB, %d DLLs loaded",
          length(nms),
          rss,
          length(getLoadedDLLs())
        ))
      }
      out
    },
    args = list(nms = nms, dev = dev, root = root, helper = helper, trace = trace),
    show = trace
  )
  if (trace) {
    # The parent's own footprint, to tell a child-side leak from one in the
    # process that collects the results.
    rss <- as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
    message(sprintf("solve-gate parent after block: RSS %.0f MB, %d DLLs loaded", rss, length(getLoadedDLLs())))
  }
  res
}

# Probe results keyed by model name, each model probed at most once however
# many tests ask for it: the tests below cover overlapping sets of models. An
# environment rather than a variable so the memo survives across test_that()
# blocks without superassignment.
.gateCache <- new.env(parent = emptyenv())
.gateCache$res <- list()
.gateProbe <- function(names_) {
  todo <- setdiff(names_, names(.gateCache$res))
  if (length(todo)) {
    blocks <- split(todo, ceiling(seq_along(todo) / .probeBlockSize))
    .gateCache$res <- c(.gateCache$res, do.call(c, unname(lapply(blocks, .probeInChild))))
  }
  .gateCache$res[names_]
}

test_that("every model in the registry responds to a dose", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  scope <- .solveGateScope()
  names_ <- if (scope == "full") db$name else intersect(db$name, linCmtRiskCandidates())
  expect_gt(length(names_), 0)

  res <- .gateProbe(names_)
  dead <- character(0)
  brokenByConversion <- character(0)
  unsupported <- character(0)
  for (n in names_) {
    r <- res[[n]]$lin
    if (identical(r$status, "skip")) {
      next
    }
    if (!identical(r$status, "ok")) {
      # An outright solve failure has two very different causes, and they must
      # not be conflated. If the same model solves with the conversion
      # disabled, the conversion is the problem -- it renumbers compartments,
      # so an event table that names an endpoint compartment stops resolving
      # and rxode2 raises "'dvid'->'cmt' or 'cmt' on observation record". If it
      # fails both ways, the probe cannot drive the model.
      ode <- res[[n]]$ode
      if (identical(ode$status, "ok")) {
        brokenByConversion <- c(brokenByConversion, n)
      } else {
        unsupported <- c(unsupported, n)
      }
      next
    }
    # A dose that changes nothing anywhere in the model -- not one state, not
    # one computed quantity -- did not arrive. This is a difference and not a
    # ratio on purpose: a ratio of two zeros passes, which is how the
    # quarantined models accumulated in the first place.
    if (is.na(r$doseEffect) || r$doseEffect <= 0) dead <- c(dead, n)
  }

  # A model the probe cannot drive, or that is known to be broken, can also
  # look dead (Hopkins_2024_amisulpride, whose probe selects no formulation).
  deadExempt <- c(knownLinCmtDropModels, probeUnsupportedModels, knownBrokenModels)
  expect_equal(
    setdiff(dead, deadExempt),
    character(0),
    info = paste0(
      "Model(s) that a dose does not reach under default rxSolve() arguments: ",
      paste(setdiff(dead, deadExempt), collapse = ", "),
      ". Try rxSolve(..., useLinCmt = FALSE); if that fixes it the model has ",
      "hit the ODE-to-linCmt conversion defect described at the top of this ",
      "file."
    )
  )
  expect_equal(
    setdiff(brokenByConversion, knownLinCmtDropModels),
    character(0),
    info = paste0(
      "Model(s) that fail to solve under default rxSolve() arguments but ",
      "solve with useLinCmt = FALSE: ",
      paste(setdiff(brokenByConversion, knownLinCmtDropModels), collapse = ", "),
      ". The conversion has renumbered the compartments out from under the ",
      "event table."
    )
  )
  expect_equal(
    setdiff(unsupported, c(probeUnsupportedModels, knownBrokenModels)),
    character(0),
    info = paste0(
      "Model(s) the probe could not solve either way: ",
      paste(setdiff(unsupported, c(probeUnsupportedModels, knownBrokenModels)), collapse = ", "),
      ". Either the model is broken or helper-solveProbe.R needs to learn how ",
      "to drive it -- work out which, then add it to knownBrokenModels or to ",
      "probeUnsupportedModels accordingly. Do not put a broken model in ",
      "probeUnsupportedModels: that claims the model is fine."
    )
  )
})

test_that("the non-finite check reads ODE states only, not timing outputs", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtStomachDosedModel))
  solved <- data.frame(
    time = c(0, 1),
    stomach = c(1, NaN),
    central = c(0, Inf),
    Cc = c(0, NA),
    tad = c(NA, 1)
  )
  expect_identical(probeNonFiniteStates(ui, solved), c("stomach", "central"))
  solved$stomach <- c(1, 0.5)
  solved$central <- c(0, 1)
  # `Cc` and `tad` are still NA, and are not states.
  expect_identical(probeNonFiniteStates(ui, solved), character(0))
  # A state missing from the output (a linCmt() conversion drops them) is not
  # looked for.
  expect_identical(probeNonFiniteStates(ui, solved[, c("time", "Cc")]), character(0))
})

test_that("the probe doses each declared compartment whose dose history the model reads", {
  # .dtTwoDepotModel reads tad(depot), tad0(depot2) and tafd(depot3).
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtTwoDepotModel))
  db <- data.frame(name = "m", dosing = "depot,depot2")
  expect_identical(probeDoseCmts("m", db, ui), c("depot", "depot2"))
  # An undeclared compartment is never dosed, so a wrong `dosing` still leaves
  # its dose-history reads NA.
  db$dosing <- "depot"
  expect_identical(probeDoseCmts("m", db, ui), "depot")
  # The first declared compartment is always dosed, read or not.
  db$dosing <- "central,depot2"
  expect_identical(probeDoseCmts("m", db, ui), c("central", "depot2"))
  db$dosing <- NA_character_
  expect_identical(probeDoseCmts("m", db, ui), NA_character_)
})

test_that("the probe event table holds one dose record per dosed compartment", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtTwoDepotModel))
  ev <- as.data.frame(probeEvents(ui, c("depot", "depot2"), 100, c(0, 1)))
  doses <- ev[ev$evid == 1L, ]
  expect_identical(doses$cmt, c("depot", "depot2"))
  expect_equal(doses$amt, c(100, 100))
  expect_equal(doses$time, c(0, 0))
})

test_that("the probe is unchanged for a model reading only its first dosing compartment", {
  # Most of the registry: a second dosing compartment is declared (or none),
  # and every dose-history read names the first one. The probe must dose just
  # that compartment, as it did before it learned to dose more than one, so
  # these models are probed exactly as before.
  times <- c(0, 1, 2)
  # tad(depot), with central declared but not read.
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtFallbackDepotModel))
  db <- data.frame(name = "m", dosing = "depot,central")
  expect_identical(probeDoseCmts("m", db, ui), "depot")
  expect_identical(
    as.data.frame(probeEvents(ui, probeDoseCmts("m", db, ui), 100, times)),
    as.data.frame(rxode2::et(rxode2::et(amt = 100, cmt = "depot"), times))
  )
  # podo(stomach) and tad(stomach), stomach the only dosing compartment.
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtStomachDosedModel))
  db <- data.frame(name = "m", dosing = "stomach")
  expect_identical(probeDoseCmts("m", db, ui), "stomach")
  # Bare podo() and tad() name no compartment, so nothing is added.
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.dtBareTimingModel))
  db <- data.frame(name = "m", dosing = "central,stomach")
  expect_identical(probeDoseCmts("m", db, ui), "central")
})

test_that("doseTimingCandidates() finds the models that read a compartment's dose history", {
  skip_on_cran()
  cand <- doseTimingCandidates()
  # podo(stomach) and tad(stomach)
  expect_true("Bulitta_2009_cefuroxime_axetil" %in% cand)
  # tafd(bact_susceptible), the bacterial inoculum
  expect_true("Thorsted_2023_neutrophil_mouse" %in% cand)
  # Bare first-order absorption: no dose-history call at all
  expect_false("PK_1cmt_des" %in% cand)
})

test_that("the probe gives a within-subject grouping variable two levels", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbOccModel))
  expect_identical(.probeGroupingVars(ui), "OCC")
  ev <- probeEvents(ui, "central", 100, c(0, 1, 2))
  # Rows: observation at 0, dose at 0, observations at 1 and 2.
  expect_identical(ev$OCC, c(1L, 1L, 1L, 2L))
  # Without the column rxSolve cannot find OCC; with one level it asks for the
  # level's THETA. With two it solves.
  r <- .probeSolveUi(ui, "central", useLinCmt = FALSE)
  expect_identical(r$status, "ok")
  expect_gt(r$doseEffect, 0)
  # Subject-level etas only: no extra column.
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbLagModel))
  expect_identical(.probeGroupingVars(ui), character(0))
  expect_false("OCC" %in% names(as.data.frame(probeEvents(ui, "depot", 100, c(0, 1)))))
})

test_that("the `first` observation form tags every observation with the first endpoint", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbTwoEndpointModel))
  ev <- as.data.frame(probeEvents(ui, "central", 100, c(0, 1), "first"))
  expect_identical(ev$cmt[ev$evid == 0L], c("Cc", "Cc"))
  ev <- as.data.frame(probeEvents(ui, "central", 100, c(0, 1), "cmt"))
  expect_identical(sort(ev$cmt[ev$evid == 0L]), c("Cc", "Cc", "effect", "effect"))
})

test_that("covariates the unit table would put out of domain get typical values", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbNeonateModel))
  cov <- probeCovariateValues(ui)
  expect_identical(
    cov$values[c("PNA", "GA", "HCT", "T_INFUSION_LAG")],
    c(PNA = 0.46, GA = 38, HCT = 0.45, T_INFUSION_LAG = 0.17)
  )
  expect_identical(
    unname(cov$why[c("PNA", "HCT")]),
    c("typical PNA in months", "typical HCT in fraction")
  )
  r <- .probeSolveUi(ui, "central", useLinCmt = FALSE)
  expect_identical(r$status, "ok")
  expect_gt(r$doseEffect, 0)
})

test_that(".probeUnitKey() reads the scale from the leading word of a units string", {
  expect_identical(.probeUnitKey("days in this paper (canonical PNA is months)"), "days")
  expect_identical(.probeUnitKey("months"), "months")
  expect_identical(.probeUnitKey("week"), "weeks")
  expect_identical(.probeUnitKey("h"), "h")
  expect_identical(.probeUnitKey("hours"), "h")
  expect_identical(.probeUnitKey("% (volume fraction times 100)"), "percent")
  expect_identical(.probeUnitKey("percent"), "percent")
  expect_identical(.probeUnitKey("(fraction)"), "fraction")
  expect_identical(.probeUnitKey("fraction (0-1)"), "fraction")
  expect_identical(.probeUnitKey("L/L (volume fraction)"), "fraction")
  expect_identical(.probeUnitKey("pmol/L"), "pmol")
  expect_identical(.probeUnitKey("mU/L (paper internal units)"), "mu")
  expect_identical(.probeUnitKey("uU/mL"), "mu")
  expect_identical(.probeUnitKey("kg"), NA_character_)
  expect_identical(.probeUnitKey(NULL), NA_character_)
})

test_that("a normalising constant held in a named variable is found", {
  body <- quote({
    ref_wt <- 88
    neg <- -2
    twice <- 1
    twice <- 2
    cl <- exp(lcl) * (1 + e_wt * (WT - ref_wt)) * (AGE - twice)
  })
  # `twice` has no single value, so it is left out.
  expect_identical(.probeConstantAssignments(body), c(ref_wt = 88, neg = -2))
  expect_identical(.probeNormalisingConstant(body, "WT"), 88)
  expect_identical(.probeNormalisingConstant(body, "AGE"), NA_real_)
  expect_identical(
    .probeNormalisingConstant(
      quote({
    ref <- 70
    x <- ref / WT
  }),
      "WT"
    ),
    70
  )
  # A literal still wins without any assignment.
  expect_identical(
    .probeNormalisingConstant(
      quote({
    x <- (WT / 70)^0.75
  }),
      "WT"
    ),
    70
  )
})

test_that("the probe lengthens the grid when an absorption lag outlasts it", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbLagModel))
  r <- .probeSolveUi(ui, "depot", useLinCmt = FALSE)
  expect_identical(r$status, "ok")
  expect_identical(r$horizon, .probeLongHorizon)
  expect_identical(max(r$dosed$time), 480)
  expect_gt(r$doseEffect, 0)
  # A 1 h lag responds on the default grid, which is kept.
  short <- rxode2::rxUiDecompress(rxode2::as.rxUi(rxode2::ini(.pbLagModel, ltlag = log(1))))
  r <- .probeSolveUi(short, "depot", useLinCmt = FALSE)
  expect_identical(r$horizon, 1)
  expect_identical(max(r$dosed$time), 48)
})

test_that("two paths on different grids or observation forms count as fully divergent", {
  a <- list(dosed = data.frame(time = c(0, 1), Cc = c(1, 2)), cols = "Cc", horizon = 1, eventForm = "cmt")
  expect_identical(probeMaxRelDiff(a, a), 0)
  expect_identical(probeMaxRelDiff(a, utils::modifyList(a, list(horizon = 10))), 1)
  expect_identical(probeMaxRelDiff(a, utils::modifyList(a, list(eventForm = "first"))), 1)
})

test_that("the probe falls back to dop853 when liblsoda cannot solve a model", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  db <- nlmixr2lib::modeldb
  # liblsoda stops with "could not solve the system".
  r <- probeSolveModel("Svensson_2016_rifampicin", db, useLinCmt = FALSE)
  expect_identical(r$status, "ok")
  expect_identical(r$method, "dop853")
  expect_identical(r$nonFiniteStates, character(0))
  # liblsoda returns every state NaN, with no warning.
  r <- probeSolveModel("Ho_2025_betalactam_pbp78_CRAB", db, useLinCmt = FALSE)
  expect_identical(r$method, "dop853")
  expect_identical(r$nonFiniteStates, character(0))
  # A model liblsoda solves keeps the default.
  expect_identical(probeSolveModel("PK_1cmt", db, useLinCmt = FALSE)$method, "default")
})

test_that("a model mixing an algebraic endpoint with ODE-state endpoints probes with the `first` form", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  r <- probeSolveModel("AlvarezJimenez_2016_scopolamine", nlmixr2lib::modeldb, useLinCmt = FALSE)
  expect_identical(r$status, "ok")
  expect_identical(r$eventForm, "first")
  expect_gt(r$doseEffect, 0)
})

test_that("a linCmt divergence that is ODE solver error disappears at tight tolerances", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  db <- nlmixr2lib::modeldb
  lin <- probeSolveModel("Matcha_2025_amikacin", db, useLinCmt = TRUE)
  ode <- probeSolveModel("Matcha_2025_amikacin", db, useLinCmt = FALSE)
  expect_gt(probeMaxRelDiff(lin, ode), 1e-4)
  expect_lt(probeTightRelDiff("Matcha_2025_amikacin", db), 1e-4)
})

test_that("a dosed probe solve leaves every ODE state finite", {
  # A state goes NaN when a derivative does, and the usual cause in this
  # registry is a dose-history function -- `tad(cmt)`, `podo(cmt)`,
  # `tafd(cmt)` -- reading a compartment the probe never doses. The probe
  # doses the registry's first `dosing` compartment and any other declared one
  # whose dose history the model reads (probeDoseCmts()), so a model whose
  # `dosing` metadata leaves out the compartment its timing reads (or that has
  # none, and gets the depot/central fallback) solves with NaN states. The solve does not
  # error: every state is NA from the first step after the dose. The first test
  # above passes such a model, because `doseEffect` ignores NA and the time-0
  # row still shows the dose.
  #
  # Timing outputs themselves (`tad()`, a time since inoculation) are NA until
  # their first event by design, so only `d/dt()` states are checked; see
  # probeNonFiniteStates().
  #
  # On the default path this covers the linCmt-risk screen plus every model
  # that reads a named compartment's dose history (doseTimingCandidates()),
  # which is the class this check exists for; under `full`, every model.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  names_ <- if (.solveGateScope() == "full") {
    db$name
  } else {
    intersect(db$name, union(linCmtRiskCandidates(), doseTimingCandidates()))
  }
  expect_gt(length(names_), 0)

  res <- .gateProbe(names_)
  exempt <- c(probeUnsupportedModels, knownBrokenModels)
  nonFinite <- character(0)
  for (n in setdiff(names_, exempt)) {
    states <- unique(c(res[[n]]$lin$nonFiniteStates, res[[n]]$ode$nonFiniteStates))
    if (length(states)) {
      nonFinite <- c(nonFinite, sprintf("%s (%s)", n, paste(states, collapse = ", ")))
    }
  }

  expect_equal(
    nonFinite,
    character(0),
    info = paste0(
      "Model(s) whose dosed probe solve leaves an ODE state NaN or infinite: ",
      paste(nonFinite, collapse = "; "),
      ". Check that every compartment the model's tad()/podo()/tafd()/... ",
      "calls name is in its `dosing` metadata (checkModelConventions() ",
      "reports this as dose_timing_compartment), and that the first entry is ",
      "the usual route. If the model is genuinely broken, add it to ",
      "knownBrokenModels with the evidence."
    )
  )
})

test_that("rxode2's linCmt() optimisation never changes a model's solution", {
  # `useLinCmt = TRUE` is an optimisation: it must return the same numbers as
  # the ODE solve, or it is not an optimisation but a different model. This
  # catches both failure modes of the conversion -- the dead model and the
  # silently-wrong one that drops a transit chain but still produces plausible
  # concentrations. The quiet one is the larger exposure, because a plausible
  # profile is invisible to every other check in this package.
  #
  # Under `full` this runs over the WHOLE registry rather than the screened
  # subset. That is deliberate: the screen looks for an exogenous input term,
  # so it covers the models this defect is known to hit, but the invariant
  # being asserted -- the optimisation does not change the answer -- belongs to
  # every model rxode2 is willing to convert, and enumerating the registry is
  # the only way to say that without re-implementing rxode2's detector here.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  names_ <- if (.solveGateScope() == "full") {
    db$name
  } else {
    intersect(db$name, linCmtRiskCandidates())
  }
  expect_gt(length(names_), 0)

  res <- .gateProbe(names_)
  divergent <- character(0)
  for (n in names_) {
    lin <- res[[n]]$lin
    ode <- res[[n]]$ode
    if (!identical(lin$status, "ok") || !identical(ode$status, "ok")) {
      next
    }
    if (probeMaxRelDiff(lin, ode) > 1e-4) divergent <- c(divergent, n)
  }
  # ODE solver error alone can exceed 1e-4; a conversion that changes the model
  # does not shrink with the tolerance. So what is left is re-solved at tight
  # tolerances before it is reported (see .probeTightSolver). Quarantined,
  # unsupported and broken models are not re-solved: none of them is expected
  # to agree.
  candidates <- setdiff(divergent, c(knownLinCmtDropModels, probeUnsupportedModels, knownBrokenModels))
  tight <- vapply(candidates, probeTightRelDiff, numeric(1), db = db)
  newDivergent <- candidates[is.na(tight) | tight > 1e-4]

  expect_equal(
    newDivergent,
    character(0),
    info = paste0(
      "Model(s) that solve differently with and without rxode2's linCmt ",
      "conversion, even at tight ODE tolerances: ",
      paste(newDivergent, collapse = ", "),
      ". The conversion has changed the model; see the comment above ",
      "knownLinCmtDropModels for the causes found so far."
    )
  )
})

test_that("the linCmt-drop quarantine contains no model that already solves", {
  # Without this, a quarantine entry outlives the bug it documents and the
  # gates above silently stop covering that model. When rxode2 ships the fix,
  # this test goes red until the list is emptied -- which is the point.
  #
  # The invariant checked here is "the conversion still changes the outcome",
  # not zero-ness, because the defect has three distinct symptoms and the
  # quarantine has to cover all of them: a model that dies, a model that loses
  # its transit chain but still returns plausible concentrations, and a model
  # that fails to solve at all because the conversion renumbered its
  # compartments.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  #
  # The quarantined models are probed in child processes through .gateProbe(),
  # like every other probe in this file, and share its memo.
  db <- nlmixr2lib::modeldb
  expect_equal(
    setdiff(knownLinCmtDropModels, db$name),
    character(0),
    info = "Quarantined but not in the registry; remove from knownLinCmtDropModels."
  )
  res <- .gateProbe(intersect(knownLinCmtDropModels, db$name))
  for (nm in names(res)) {
    lin <- res[[nm]]$lin
    ode <- res[[nm]]$ode
    stillBroken <- identical(ode$status, "ok") &&
      (!identical(lin$status, "ok") || probeMaxRelDiff(lin, ode) > 1e-4)
    expect_true(
      stillBroken,
      # nolint next: line_length_linter.
      label = paste0(
        nm,
        " now solves the same with and without rxode2's linCmt conversion,",
        " so it must be removed from knownLinCmtDropModels"
      )
    )
  }
})
