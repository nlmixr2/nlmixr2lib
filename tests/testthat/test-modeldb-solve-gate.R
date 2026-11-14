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

# Models whose solution rxode2's ODE-to-linCmt conversion changes, by cause.
# Each needs a fix in rxode2, not in the model file; the models stay faithful
# to their papers. The conversion is opt-in (`useLinCmt = TRUE`), so a solve
# with rxSolve()'s defaults is unaffected. The input-dropping defect this list
# was first made for (https://github.com/nlmixr2/rxode2/issues/1370) is fixed
# in rxode2 5.1.8, the DESCRIPTION minimum. Every entry was confirmed with
# probeRecheck(): at atol = rtol = 1e-12, and, where the two paths needed
# different observation forms, in the form the converted model accepts. The
# quarantine test below goes red for any entry the conversion no longer
# changes.
#
# (1) No observation form resolves after the conversion (4 models). Several
#     endpoints, some algebraic: the conversion renumbers the compartments,
#     and every observation form fails with "'dvid'->'cmt' or 'cmt' on
#     observation record or on a undefined compartment".
#
# (2) Renumbered, and the values drift too (7 Chen 2024 rat models). The
#     converted model accepts only the `first` observation form; solved in
#     that form on both paths, the values still differ by 5.7e-4 to 4.9e-2.
#
# (3) States the model still reads (2 models). Porchet 1987 names its tissue
#     states `other` and `brain` and reads them in `Cother <- other / ...`
#     and `Cbrain <- brain / ...`. The conversion replaces them with linCmt's
#     internal peripheral1 and peripheral2, and the solve stops with "The
#     following parameter(s) are required for solving: peripheral2,
#     peripheral1".
#
# (4) Time-dependent parameters (27 models). A clearance or absorption rate
#     that depends on `t`, `time` or `tad()` (time-varying clearance of the
#     monoclonal antibodies, Desai 2016's Weibull absorption) is held constant
#     between records by the analytical solution, which only approximates
#     the ODE: 1.1e-4 to 0.21 apart at tight tolerances.
#
# (5) `kel <- cl / vc` with estimated micro-constants (9 models). A
#     two-compartment model that derives `kel` from cl and vc but estimates
#     k12 and k21 directly converts to a ONE-compartment solution: the
#     peripheral compartment is dropped and Cc follows `dose / vc *
#     exp(-kel * t)` (Valadez 2025 cefepime: Cc 44% off). Writing
#     `-cl / vc * central` inline, estimating kel directly, or deriving k12
#     and k21 from q avoids it.
knownLinCmtDropModels <- c(
  # (1) no observation form resolves
  "Dimelow_2018_avibactam",
  "Dimelow_2018_ceftazidime",
  "Ji_2025_eravacycline",
  "vanSteeg_2007_atenolol_noniso_rat",
  # (2) renumbered, and the values drift
  "Chen_2024_chlorogenicAcid_rat",
  "Chen_2024_hydroxysafflorYellowA_rat",
  "Chen_2024_kaempferol3OGlucoside_rat",
  "Chen_2024_kaempferol3ORutinoside_rat",
  "Chen_2024_nAcetylglutamine_rat",
  "Chen_2024_pCoumaricAcid_rat",
  "Chen_2024_rutin_rat",
  # (3) states the model still reads
  "Porchet_1987_nicotine_human_pbpk",
  "Porchet_1987_nicotine_rabbit_pbpk",
  # (4) time-dependent parameters
  "Abegesah_2025_durvalumab",
  "Bajaj_2017_nivolumab",
  "Bonate_2004_apomine",
  "Chang_2023_rifapentine_mouse",
  "Chen_2024_adebrelimab",
  "Cloesmeijer_2020_clonidine",
  "Desai_2016_isavuconazole",
  "Gibiansky_2014_obinutuzumab",
  "Gibiansky_2021_rituximab",
  "Gupta_2006_peginterferon_alfa_2b",
  "Huppe_2023_fosfomycin",
  "Miao_2023_teclistamab",
  "OHanlon_2024_milrinone",
  "Papathanasiou_2025_belantamab",
  "PK_2cmt_tdcl_des",
  "Prasanchaimontri_2024_favipiravir",
  "Rozman_2017_rituximab",
  "Sheng_2026_sugemalimab",
  "Wada_2023_sparsentan",
  "Wang_2024_sugemalimab",
  "Wang_2025_serplulimab",
  "Yamada_2025_zolbetuximab",
  "Yang_2021_cemiplimab",
  "Zeng_2026_atezolizumab",
  "Zhang_2019b_nivolumab",
  "Zhao_2025_nivolumab",
  "Zhao_2026_durvalumab",
  # (5) kel = cl / vc with estimated k12 and k21
  "Carral_2023_levofloxacin",
  "Cleary_2026_methylprednisolone",
  "Ekhart_2008_carboplatin",
  "Mohanan_2017_fludarabine",
  "Nath_2007_melphalan",
  "Stott_2018_amphotericinB_deoxycholate",
  "Stott_2022_amphotericinB_liposomal",
  "Valadez_2025_cefepime",
  "Yang_2025_rifampicin"
)

# Models the conversion renumbers without changing their values (32). Each has
# several endpoints, some algebraic. The converted model rejects the ODE
# model's per-endpoint observation records, so the two paths need different
# observation forms, and probeMaxRelDiff() scores that as a full divergence.
# Solved on both paths in the one form the converted model accepts (`first`),
# they agree to within 5e-12. A user who observes with `cmt = <endpoint>` and
# opts into `useLinCmt = TRUE` still gets the error, which is the rxode2 issue
# these entries record; the solution itself is unchanged.
knownLinCmtFormModels <- c(
  "Ahmed_2015_topiramate",
  "Beijer_2026_cloxacillin",
  "Broeker_2018_tigecycline",
  "Bukkems_2022_factorviii_plasminPeak",
  "Bukkems_2022_factorviii_thrombinPeak",
  "Bukkems_2022_factorviii_thrombinPotential",
  "ButraguenoLaiseca_2025_teicoplanin",
  "Diao_2016_daclizumab_treg",
  "Eyler_2014_ertapenem",
  "Gastonguay_2005_efaproxiral",
  "Jager_2020_flucloxacillin",
  "Jin_2025_benralizumab",
  "Kim_2025_evogliptin",
  "Komatsu_2026_cefazolin",
  "Lee_2024_eculizumab",
  "Mo_2024_shikimicAcid_pig_iv",
  "Mukap_2025_amoxicillin",
  "Papachristos_2020_bevacizumab_pkpd",
  "PerezGuille_2018_dexmedetomidine",
  "Renaud_2020_glenzocimab",
  "Valke_2024_factorviii_plasminPeak",
  "Valke_2024_factorviii_thrombinPeak",
  "Valke_2024_factorviii_thrombinPotential",
  "vanHasselt_2014_cefazolin_empirical",
  "vanHasselt_2014_cefazolin_semiphysiological",
  "vanSteeg_2007_isoprenaline_rat",
  "Wald_2022_brexanolone",
  "Wang_2024_5hydroxysaxagliptin_rat",
  "Wang_2024_saxagliptin_rat",
  "Wittau_2015_meropenem",
  "Yan_2019_benralizumab",
  "Zhang_2012_bivalirudin"
)

# Models the probe solves only with its fallbacks (see helper-solveProbe.R).
# Each list is pinned: a model newly needing a fallback fails the gate below,
# because under rxSolve()'s defaults it would fail or give NaN states, and an
# entry that no longer needs its fallback fails until it is removed.
#
# liblsoda, rxSolve()'s default solver, cannot solve these (6), while dop853
# can. Ho 2025 and the four Mahadevan 2026 time-kill models hold an undosed
# drug at exactly 0 inside `c^hill` with hill < 1: liblsoda returned NaN
# states without a warning before rxode2 PR 1542, and fails the subject since.
# Their vignettes pass `method`. Svensson 2016 fails under liblsoda outright.
knownLsodaFailModels <- c(
  "Ho_2025_betalactam_pbp78_CRAB",
  "Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP28",
  "Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP61",
  "Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP67",
  "Mahadevan_2026_polymyxinB_meropenem_fosfomycin_KP00521",
  "Svensson_2016_rifampicin"
)

# An absorption lag longer than the probe's default grid (3): the dose changes
# nothing until 85.7 h (Jeon 2013), 235 h (Khwarg 2024) and 299 h (Park 2025),
# past the 48 h grid, so the probe retries on a grid ten times longer.
knownLongLagModels <- c(
  "Jeon_2013_interferonAlfa2a",
  "Khwarg_2024_donepezil_im",
  "Park_2025_donepezil_im"
)

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
#
# Wang_2024_evorpacept_qsp: rxSolve() cannot solve it at its default tolerances
# ("could not solve the system" on both paths). Its amounts span about ten
# orders of magnitude in nmol, from whole-body antibody pools down to
# per-synapse receptor densities, and its vignette solves with atol = 1e-18,
# rtol = 1e-8 and maxsteps = 5e6. With those settings the probe's dose into
# q_V_C_aCD47 solves on both paths and every state stays finite.
probeUnsupportedModels <- c(
  "Kim_2024_islatravir",
  "Hopkins_2024_amisulpride",
  "Tan_2025_adomeglivant_ffh",
  "Wang_2024_evorpacept_qsp"
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
knownBrokenModels <- c(
  "indirect_1cpt_stim_kin_r0rmaxcrmax",
  "indirect_prec_1cpt_stim_r0rmaxcrmax"
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

.probeInChild <- function(nms, task = "pair", forms = NULL) {
  dev <- pkgload::is_dev_package("nlmixr2lib")
  root <- if (dev) normalizePath(testthat::test_path("..", "..")) else NA_character_
  helper <- normalizePath(testthat::test_path("helper-solveProbe.R"))
  trace <- nzchar(Sys.getenv("NLMIXR2LIB_SOLVE_GATE_TRACE")) && .Platform$OS.type == "unix"
  res <- callr::r(
    function(nms, dev, root, helper, trace, task, forms) {
      if (dev) {
        pkgload::load_all(root, quiet = TRUE)
      } else {
        library(nlmixr2lib)
      }
      probes <- new.env()
      sys.source(helper, envir = probes)
      db <- nlmixr2lib::modeldb
      out <- lapply(nms, function(n) {
        if (identical(task, "recheck")) {
          return(probes$probeRecheck(n, db, forms[[n]]))
        }
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
    args = list(nms = nms, dev = dev, root = root, helper = helper, trace = trace, task = task, forms = forms),
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

# probeRecheck() of each model, in child processes in blocks like
# .gateProbe(). `forms` names the common observation form per model (NULL to
# let each path find its own). Returns a named numeric vector.
.gateRecheck <- function(names_, forms) {
  if (!length(names_)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  blocks <- split(names_, ceiling(seq_along(names_) / .probeBlockSize))
  out <- unlist(lapply(blocks, function(b) .probeInChild(b, task = "recheck", forms = forms[b])))
  stats::setNames(as.numeric(out), names_)
}

# Which fallback a probe pair needed: dop853 for either path, or the longer grid.
.usedSolverFallback <- function(pair) {
  any(c(pair$lin$method, pair$ode$method) %in% "dop853")
}
.usedLongGrid <- function(pair) {
  any(c(pair$lin$horizon, pair$ode$horizon) %in% .probeLongHorizon)
}

test_that("every model in the registry responds to a dose", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  scope <- .solveGateScope()
  names_ <- if (scope == "full") db$name else intersect(db$name, linCmtRiskCandidates())
  # A floor on the screened scope, so an empty screen cannot pass unnoticed.
  expect_gte(length(names_), 452L)

  res <- .gateProbe(names_)
  dead <- character(0)
  deadOde <- character(0)
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
    if (is.na(r$doseEffect) || r$doseEffect <= 0) {
      dead <- c(dead, n)
    }
    o <- res[[n]]$ode
    if (identical(o$status, "ok") && !isTRUE(o$doseEffect > 0)) deadOde <- c(deadOde, n)
  }

  # A model the probe cannot drive, or that is known to be broken, can look
  # dead on both paths (Hopkins_2024_amisulpride, whose probe selects no
  # formulation). A quarantined model is exempt only on the converted path:
  # its ODE path must still respond to a dose.
  quarantined <- c(knownLinCmtDropModels, knownLinCmtFormModels)
  notDriven <- c(probeUnsupportedModels, knownBrokenModels)
  deadFail <- union(
    setdiff(dead, c(quarantined, notDriven)),
    setdiff(intersect(deadOde, quarantined), notDriven)
  )
  expect_equal(
    deadFail,
    character(0),
    info = paste0(
      "Model(s) that a dose does not reach under default rxSolve() arguments: ",
      paste(deadFail, collapse = ", "),
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

# A single-compartment probe event table, written independently of
# probeEvents(): one dose record into `doseCmt`, then the observations in the
# requested form. probeEvents() must build exactly this table when it is given
# one compartment.
.singleCmtEvents <- function(ui, doseCmt, amt, times, form = "cmt") {
  ev <- rxode2::et(amt = amt, cmt = doseCmt)
  eps <- ui$predDf$cond
  if (length(eps) > 1L && !identical(form, "plain")) {
    if (identical(form, "dvid")) {
      for (i in seq_along(eps)) {
        ev <- rxode2::et(ev, times, dvid = i)
      }
    } else {
      for (ep in eps) {
        ev <- rxode2::et(ev, times, cmt = ep)
      }
    }
  } else {
    ev <- rxode2::et(ev, times)
  }
  ev
}

test_that("the probe is unchanged for a model reading only its first dosing compartment", {
  # Most of the registry: a second dosing compartment is declared (or none),
  # and every dose-history read names the first one. The probe doses just that
  # compartment, with the single-compartment event table, in every form.
  cases <- list(
    # tad(depot), with central declared but not read.
    list(model = .dtFallbackDepotModel, dosing = "depot,central", first = "depot"),
    # podo(stomach) and tad(stomach), stomach the only dosing compartment.
    list(model = .dtStomachDosedModel, dosing = "stomach", first = "stomach"),
    # Bare podo() and tad() name no compartment, so nothing is added.
    list(model = .dtBareTimingModel, dosing = "central,stomach", first = "central")
  )
  times <- c(0, 1, 2)
  for (cs in cases) {
    ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(cs$model))
    db <- data.frame(name = "m", dosing = cs$dosing)
    expect_identical(probeDoseCmts("m", db, ui), cs$first)
    for (form in .probeEventForms) {
      expect_identical(
        as.data.frame(probeEvents(ui, probeDoseCmts("m", db, ui), 100, times, form)),
        as.data.frame(.singleCmtEvents(ui, cs$first, 100, times, form)),
        info = paste(cs$first, form)
      )
    }
  }
})

test_that("the probe stops when it cannot read a model() block for its dose-history reads", {
  # Reading nothing would quietly dose only the first compartment.
  ui <- new.env()
  ui$lstExpr <- NULL
  db <- data.frame(name = "m", dosing = "depot,depot2")
  expect_error(probeDoseCmts("m", db, ui), "cannot read the model\\(\\) block of m")
})

test_that("a model file without a model() block stops the dose-history screen", {
  path <- withr::local_tempfile(fileext = ".R")
  writeLines("notTheFileName <- function() 1", path)
  expect_error(.probeReadsDoseHistory(path), "no model\\(\\) block in ")
})

test_that("doseTimingCandidates() finds the models that read a compartment's dose history", {
  skip_on_cran()
  cand <- doseTimingCandidates()
  # A floor rather than an exact count, so a newly added model that reads a
  # compartment's dose history does not fail this test.
  expect_gte(length(cand), 126L)
  # podo(stomach) and tad(stomach)
  expect_true("Bulitta_2009_cefuroxime_axetil" %in% cand)
  # tafd(bact_susceptible), the bacterial inoculum
  expect_true("Thorsted_2023_neutrophil_mouse" %in% cand)
  # Bare first-order absorption: no dose-history call at all
  expect_false("PK_1cmt_des" %in% cand)
})

test_that("linCmtRiskCandidates() screens every model file and stops on an unreadable one", {
  skip_on_cran()
  # A floor rather than an exact count, so a newly added model the screen
  # selects does not fail this test.
  expect_gte(length(linCmtRiskCandidates()), 452L)
  path <- withr::local_tempfile(fileext = ".R")
  writeLines("notTheFileName <- function() 1", path)
  expect_error(.probeScreenOneFile(path), "no model\\(\\) block in ")
})

test_that("a reference value in covariateData beats a constant read from the equations", {
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(.pbMetadataReferenceModel))
  cov <- probeCovariateValues(ui)
  # The model also contains `FPG - 3.5`, which the equation rule would read.
  expect_identical(cov$values[["FPG"]], 8.2)
  expect_identical(cov$why[["FPG"]], "covariateData$reference_value")
  # Without metadata the equation rule still applies.
  expect_identical(cov$values[["WT"]], 70)
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

test_that("probeCommonForm() picks the converted path's form only when the paths differ", {
  expect_identical(probeCommonForm(list(eventForm = "first"), list(eventForm = "cmt")), "first")
  expect_null(probeCommonForm(list(eventForm = "cmt"), list(eventForm = "cmt")))
})

test_that("a form-only linCmt mismatch agrees once both paths use one form", {
  skip_on_cran()
  skip_if_not_installed("rxode2")
  db <- nlmixr2lib::modeldb
  # Jin_2025_benralizumab: the converted model accepts only the `first` form.
  lin <- probeSolveModel("Jin_2025_benralizumab", db, useLinCmt = TRUE)
  ode <- probeSolveModel("Jin_2025_benralizumab", db, useLinCmt = FALSE)
  expect_identical(c(lin$eventForm, ode$eventForm), c("first", "cmt"))
  expect_identical(probeMaxRelDiff(lin, ode), 1)
  expect_lt(probeRecheck("Jin_2025_benralizumab", db, probeCommonForm(lin, ode)), 1e-4)
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
  expect_lt(probeRecheck("Matcha_2025_amikacin", db), 1e-4)
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
  dtc <- doseTimingCandidates()
  # The dose-history models are what this test exists for; checking only that
  # the scope is non-empty would pass on the linCmt screen alone.
  expect_gte(length(intersect(db$name, dtc)), 126L)
  names_ <- if (.solveGateScope() == "full") {
    db$name
  } else {
    intersect(db$name, union(linCmtRiskCandidates(), dtc))
  }

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
  expect_gte(length(names_), 452L)

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
  # What is left is re-solved with probeRecheck() before it is reported: at
  # tight tolerances, since ODE solver error alone can exceed 1e-4, and, when
  # the two paths needed different observation forms, in the one the converted
  # model accepts, since probeMaxRelDiff() scores a form mismatch as a full
  # divergence without comparing values. Models on a list are not re-solved
  # here; the quarantine test re-checks them.
  listed <- c(knownLinCmtDropModels, knownLinCmtFormModels, probeUnsupportedModels, knownBrokenModels)
  candidates <- setdiff(divergent, listed)
  forms <- stats::setNames(lapply(candidates, function(n) probeCommonForm(res[[n]]$lin, res[[n]]$ode)), candidates)
  rc <- .gateRecheck(candidates, forms)
  formOnly <- !vapply(forms, is.null, logical(1)) & !is.na(rc) & rc <= 1e-4
  newDivergent <- candidates[is.na(rc) | rc > 1e-4]
  newFormOnly <- candidates[formOnly]

  expect_equal(
    newDivergent,
    character(0),
    info = paste0(
      "Model(s) that solve differently with and without rxode2's linCmt ",
      "conversion, even at tight ODE tolerances and in a common observation ",
      "form: ",
      paste(newDivergent, collapse = ", "),
      ". The conversion has changed the model; see the comment above ",
      "knownLinCmtDropModels for the causes found so far."
    )
  )
  expect_equal(
    newFormOnly,
    character(0),
    info = paste0(
      "Model(s) whose converted form rejects the ODE model's observation ",
      "records but agrees in a common form: ",
      paste(newFormOnly, collapse = ", "),
      ". Add them to knownLinCmtFormModels."
    )
  )
})

test_that("the linCmt quarantine contains no model the conversion no longer changes", {
  # Without this, a quarantine entry outlives the bug it documents and the
  # gates above silently stop covering that model. When rxode2 ships a fix,
  # this test goes red until the list is emptied -- which is the point.
  #
  # Each entry is re-checked as it was admitted: with probeRecheck(), at
  # tight tolerances and in the common observation form, so an entry that
  # only ever diverged through ODE solver error or a form mismatch does not
  # stay quarantined as a change of values. Probes and rechecks run in child
  # processes, like every other probe in this file.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  listed <- c(knownLinCmtDropModels, knownLinCmtFormModels)
  expect_equal(
    setdiff(listed, db$name),
    character(0),
    info = "Quarantined but not in the registry; remove from knownLinCmtDropModels or knownLinCmtFormModels."
  )
  listed <- intersect(listed, db$name)
  res <- .gateProbe(listed)
  bothOk <- listed[vapply(
    listed,
    function(n) identical(res[[n]]$lin$status, "ok") && identical(res[[n]]$ode$status, "ok"),
    logical(1)
  )]
  forms <- stats::setNames(lapply(bothOk, function(n) probeCommonForm(res[[n]]$lin, res[[n]]$ode)), bothOk)
  rc <- .gateRecheck(bothOk, forms)

  for (nm in intersect(knownLinCmtDropModels, listed)) {
    linFails <- !identical(res[[nm]]$lin$status, "ok")
    stillChanged <- identical(res[[nm]]$ode$status, "ok") &&
      (linFails || is.na(rc[nm]) || rc[nm] > 1e-4)
    expect_true(
      stillChanged,
      label = paste0(
        nm,
        " now solves the same with and without rxode2's linCmt conversion,",
        " so it must be removed from knownLinCmtDropModels"
      )
    )
  }
  for (nm in intersect(knownLinCmtFormModels, listed)) {
    stillFormOnly <- nm %in% bothOk && !is.null(forms[[nm]]) && !is.na(rc[nm]) && rc[nm] <= 1e-4
    expect_true(
      stillFormOnly,
      label = paste0(
        nm,
        " no longer needs a different observation form under the conversion, or its values now",
        " differ in the common form; remove it from knownLinCmtFormModels or move it to knownLinCmtDropModels"
      )
    )
  }
})

test_that("only the known models need the probe's solver or grid fallback", {
  # A model that needs dop853 fails, or gives NaN states, under rxSolve()'s
  # default solver; one that needs the longer grid shows no dose effect on the
  # default one. The probe hides both, so they are pinned here.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  names_ <- if (.solveGateScope() == "full") {
    db$name
  } else {
    intersect(db$name, union(linCmtRiskCandidates(), doseTimingCandidates()))
  }
  res <- .gateProbe(names_)
  solver <- names_[vapply(names_, function(n) .usedSolverFallback(res[[n]]), logical(1))]
  grid <- names_[vapply(names_, function(n) .usedLongGrid(res[[n]]), logical(1))]
  expect_equal(
    setdiff(solver, knownLsodaFailModels),
    character(0),
    info = "Model(s) the probe solved only with dop853; add them to knownLsodaFailModels with the evidence."
  )
  expect_equal(
    setdiff(grid, knownLongLagModels),
    character(0),
    info = "Model(s) whose dose shows only on the longer grid; add them to knownLongLagModels with the evidence."
  )
})

test_that("every model on the fallback lists still needs its fallback", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  res <- .gateProbe(intersect(c(knownLsodaFailModels, knownLongLagModels), db$name))
  for (nm in knownLsodaFailModels) {
    expect_true(
      .usedSolverFallback(res[[nm]]),
      label = paste(nm, "no longer needs dop853; remove it from knownLsodaFailModels")
    )
  }
  for (nm in knownLongLagModels) {
    expect_true(
      .usedLongGrid(res[[nm]]),
      label = paste(nm, "no longer needs the longer grid; remove it from knownLongLagModels")
    )
  }
})
