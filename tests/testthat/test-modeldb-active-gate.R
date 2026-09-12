# A time-varying on/off gate covariate (`*_ACTIVE`) that does not change the
# solved system is a silent, invisible correctness defect. This file is the
# mechanical gate for a class of bug that shipped in three models undetected:
#
#   For some model shapes rxode2 does not integrate the `d/dt()` the file
#   declares: it solves the linear-compartment system with its analytic kernel
#   driven by variables named `cl` (or `CL`) and `vc` (or `v`), and the
#   explicit right-hand side is discarded. So the widely-copied idiom
#
#       cl       <- <off-therapy arm>
#       cl_total <- cl + RRT_HEMODIAL_ACTIVE * cl_hemodialysis
#       kel      <- cl_total / vc
#       d/dt(central) <- -kel * central
#
#   reports a correct-looking `cl_total` and `kel` while the simulated AMOUNTS
#   decay at the off-therapy rate in BOTH gate states. The dialysis arm does
#   nothing at all. The fix is to assign the gated SUM to `cl` itself.
#
# Nothing else in this package catches it. `checkModelConventions()` does not;
# `modeldb$linCmt` reads FALSE for an affected model; and a validation vignette
# that simulates only one gate state cannot see it. Confirmed inert and repaired
# 2026-09-12: Veinstein_2013_gentamicin, Dohmann_2025_piperacillin,
# Eyler_2014_ertapenem.
#
# The assertion here is deliberately weak in VALUE and strong in DIRECTION: it
# does not check that any number matches the paper (that is each model's
# vignette's job), only that flipping the gate from 0 to 1 actually moves the
# ODE states. That is precisely the property the defect destroys.

# Physiologically plausible covariate values used only to put each model in a
# solvable state. They are NOT paper-specific and nothing here asserts accuracy
# against a source -- they exist so the solve produces finite numbers. Values
# that would make a term degenerate are avoided on purpose: BFR and DFR must
# differ and be non-zero or the Michaels equation in Liesenfeld 2013 evaluates
# 0/0.
activeGateCovariateValues <- c(
  WT = 70, HT = 170, BSA = 1.9, BMI = 25, AGE = 55, PAGE = 55, PNA = 30,
  SEXF = 0, CRCL = 60, CREAT = 1.2, ALB = 35, TBILI = 10, HCT = 30, PLT = 200,
  CRP = 50, GGT = 40, AST = 30, PT_SEC = 13, ANURIA = 0,
  BFR = 250, DFR = 500, QEFF = 2000, RRT_CRRT_EFFLUENT_FLOW = 2000,
  URINE_FLOW = 50, URINE_VOL_24H = 500, URINE_VOL_INTERVAL = 200,
  FILT_SA = 1.5, FILT_SA_MED = 1, FILT_SA_LARGE = 0,
  T_HEMODIAL_INIT = 24, T_POST_HEMODIAL = 12,
  RRT_HEMODIAL_STATUS = 1, RRT_CRRT_STATUS = 1
)

# Indicator families default to their reference level; anything else defaults to
# 1, which is a safe multiplier and a safe ratio denominator. A covariate that
# needs a specific value for the model to solve at all belongs in the map above,
# not here -- an unsolvable model is reported as a failure, not skipped.
activeGateCovariateValue <- function(nm) {
  if (nm %in% names(activeGateCovariateValues)) {
    return(unname(activeGateCovariateValues[[nm]]))
  }
  indicatorPrefixes <-
    "^(DIS_|STUDY_|RACE_|REGION_|RENALIMP_|VASCACC_|SNP_|CONMED_|APACHE_)"
  if (grepl(indicatorPrefixes, nm)) return(0)
  if (grepl("ACTIVE$", nm)) return(0)
  1
}

# Pull the model function plus the two metadata lists this test needs off a
# model file, by walking the parsed body rather than by grepping, so a name
# appearing in a comment or a prose note cannot produce a false positive.
activeGateModelPieces <- function(path) {
  env <- new.env(parent = globalenv())
  fname <- sub("[.]R$", "", basename(path))
  sys.source(path, envir = env)
  if (!exists(fname, envir = env, inherits = FALSE)) {
    stop("`", fname, "` not defined by ", path)
  }
  fn <- get(fname, envir = env)
  covariateData <- NULL
  dosing <- NULL
  for (st in as.list(body(fn))[-1]) {
    isAssign <- is.call(st) && length(st) >= 3 &&
      (identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))
    if (!isAssign || !is.name(st[[2]])) next
    target <- as.character(st[[2]])
    if (target == "covariateData") covariateData <- names(eval(st[[3]]))
    if (target == "dosing") dosing <- eval(st[[3]])
  }
  list(fn = fn, name = fname, covariateData = covariateData, dosing = dosing)
}

activeGateNames <- function(covariateData) {
  grep("[A-Z0-9]_ACTIVE$", covariateData, value = TRUE)
}

# Two warnings are expected while building and solving these models and say
# nothing about the gate:
#   * "non-mu referenced" -- emitted for the inter-occasion variability idiom
#     the affected models already ship with;
#   * "No omega parameters in the model" -- emitted by zeroRe() for the
#     deterministic models that carry no IIV to zero out.
# Both are muffled BY MESSAGE PATTERN so that every OTHER warning still
# propagates; a new warning here would be worth seeing.
activeGateExpectedWarnings <- c("non-mu referenced", "No omega parameters in the model")

activeGateMuffle <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      isExpected <- vapply(activeGateExpectedWarnings, grepl, logical(1),
                           x = msg, fixed = TRUE)
      if (any(isExpected)) invokeRestart("muffleWarning")
    }
  )
}

# Largest relative change in any ODE state when `gate` flips from 0 to 1, with
# every other covariate held fixed. Returns NA when the solve cannot be built,
# which the caller treats as a failure rather than a skip.
activeGateStateResponse <- function(pieces, gate) {
  # Build from the model FUNCTION, not from calling it: `ini()` / `model()` are
  # only resolvable once rxode2 has the function in hand.
  #
  ui <- rxode2::zeroRe(rxode2::rxUiDecompress(rxode2::rxode2(pieces$fn)))
  states <- ui$state
  doseCmt <- intersect(as.character(pieces$dosing), states)
  if (!length(doseCmt)) {
    doseCmt <- if ("depot" %in% states) "depot" else states[[1]]
  }
  ev <- rxode2::et(seq(0, 24, by = 1))
  for (cmt in doseCmt) {
    ev <- rxode2::et(ev, amt = 100, cmt = cmt, time = 0)
  }
  dat <- as.data.frame(ev)
  # Multi-endpoint models need an explicit dvid on observation rows or rxode2
  # cannot resolve which endpoint an un-compartmented observation belongs to.
  if (NROW(ui$predDf) > 1L) {
    dat$dvid <- ifelse(dat$evid == 0, 1L, NA_integer_)
  }
  for (cv in ui$allCovs) dat[[cv]] <- activeGateCovariateValue(cv)
  solveAt <- function(value) {
    dat[[gate]] <- value
    as.matrix(
      rxode2::rxSolve(ui, dat, returnType = "data.frame",
                      addDosing = FALSE)[, states, drop = FALSE]
    )
  }
  off <- solveAt(0)
  on <- solveAt(1)
  if (anyNA(off) || anyNA(on)) return(NA_real_)
  max(abs(off - on) / pmax(abs(off), abs(on), 1e-12))
}

activeGateModelFiles <- function() {
  root <- system.file("modeldb", package = "nlmixr2lib")
  if (!nzchar(root)) skip("nlmixr2lib modeldb directory not found")
  list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
}

test_that("every *_ACTIVE gate covariate changes the solved system", {
  skip_if_not_installed("rxode2")

  files <- activeGateModelFiles()
  expect_gt(length(files), 0)

  # Enumerate first, so the set this test covers is visible and a model that
  # stops carrying a gate cannot silently shrink the coverage to nothing.
  cases <- list()
  for (path in files) {
    pieces <- try(activeGateModelPieces(path), silent = TRUE)
    if (inherits(pieces, "try-error")) next
    for (gate in activeGateNames(pieces$covariateData)) {
      cases[[length(cases) + 1L]] <- list(pieces = pieces, gate = gate)
    }
  }

  # An empty (or quietly shrunken) enumeration would make every assertion below
  # vacuous, so it is a failure, not a pass. 37 model/gate pairs across 35
  # models were covered when this gate was added; this is a floor guarding
  # against a vacuous run, not an exact inventory.
  expect_gte(length(cases), 35)

  inert <- character(0)
  unsolvable <- character(0)
  for (case in cases) {
    label <- paste0(case$pieces$name, " / ", case$gate)
    response <- try(activeGateMuffle(activeGateStateResponse(case$pieces, case$gate)),
                    silent = TRUE)
    if (inherits(response, "try-error")) {
      unsolvable <- c(unsolvable, paste0(label, ": ", conditionMessage(attr(response, "condition"))))
    } else if (is.na(response)) {
      unsolvable <- c(unsolvable, paste0(label, ": solve produced NA/NaN"))
    } else if (response < 1e-10) {
      inert <- c(inert, label)
    }
  }

  expect_equal(
    inert, character(0),
    info = paste0(
      "Gate covariate(s) that do not change any ODE state: ",
      paste(inert, collapse = "; "),
      ". The usual cause is that the gated total clearance is assigned to a ",
      "separate variable (`cl_total`) while `cl` and `vc` are both defined, so ",
      "rxode2 solves the system analytically from `cl`/`vc` and discards the ",
      "explicit d/dt(). Assign the gated SUM to `cl` itself -- see ",
      "Berthaud_2025_cefazolin.R or Lee_2024_gentamicin_teigen.R."
    )
  )

  # A model this test cannot solve is NOT covered by it, so fail loudly rather
  # than degrade to a silent pass. Add a value to activeGateCovariateValues (or a
  # `dosing` field to the model) until the solve succeeds.
  expect_equal(
    unsolvable, character(0),
    info = paste0(
      "Gate covariate(s) whose model could not be solved by this gate: ",
      paste(unsolvable, collapse = "; ")
    )
  )
})
