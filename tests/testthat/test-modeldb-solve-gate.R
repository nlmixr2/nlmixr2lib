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
# rxode2's `rxSolve()` defaults to `useLinCmt = TRUE`. Under that default
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

# Models that cannot survive rxode2's ODE-to-linCmt conversion. Every entry
# here is a CORRECT model file defeated by the upstream defect described
# above; none of them is a transcription error, and none should be "fixed" by
# rewriting the model to dodge the converter. Remove entries -- and this
# comment -- in the same commit that raises the rxode2 minimum in DESCRIPTION
# to the first release that refuses to convert an ODE carrying an exogenous
# input term.
knownLinCmtDropModels <- c(
  "Bienczak_2016_efavirenz",
  "Chigutsa_2011_rifampicin",
  "Chigutsa_2012_ofloxacin",
  "Dong_2014_mycophenolic_acid",
  "Marques_2025_salbutamol",
  "ResendizGalvan_2025_cycloserine",
  "Sloan_2017_rifampicin",
  "Smythe_2013_gatifloxacin",
  "Tikiso_2021_abacavir",
  "Vinnard_2017_rifampicin",
  "Wilkins_2008_rifampicin"
)

# Models the probe cannot drive at all -- they fail to solve with the linCmt
# conversion both enabled and disabled, so whatever is wrong is not the defect
# this file is about. Keep this list as close to empty as possible. An entry
# here is an admission that the gate does not cover that model, so never add
# one without first establishing that the model itself is sound.
probeUnsupportedModels <- character(0)

# Solving every model in the registry takes roughly an hour. That is the right
# cost for a release gate and the wrong cost for an edit-run-edit loop, so the
# full sweep runs when NLMIXR2LIB_SOLVE_GATE=full and the at-risk subset runs
# otherwise. The subset is NOT a hand-written list: it is computed by
# linCmtRiskCandidates() below, which walks every model file in the registry,
# so a model added next week is screened automatically. What the default path
# gives up is the ability to catch a model that solves to zero for some reason
# unrelated to the linCmt conversion; the full sweep covers that.
.solveGateScope <- function() {
  if (identical(tolower(Sys.getenv("NLMIXR2LIB_SOLVE_GATE", "")), "full")) "full" else "screened"
}

test_that("every model in the registry responds to a dose", {
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  scope <- .solveGateScope()
  names_ <- if (scope == "full") db$name else intersect(db$name, linCmtRiskCandidates())
  expect_gt(length(names_), 0)

  dead <- character(0)
  brokenByConversion <- character(0)
  unsupported <- character(0)
  for (n in names_) {
    r <- probeSolveModel(n, db)
    if (identical(r$status, "skip")) next
    if (!identical(r$status, "ok")) {
      # An outright solve failure has two very different causes, and they must
      # not be conflated. If the same model solves with the conversion
      # disabled, the conversion is the problem -- it renumbers compartments,
      # so an event table that names an endpoint compartment stops resolving
      # and rxode2 raises "'dvid'->'cmt' or 'cmt' on observation record". If it
      # fails both ways, the probe cannot drive the model.
      ode <- probeSolveModel(n, db, useLinCmt = FALSE)
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

  expect_equal(
    setdiff(dead, knownLinCmtDropModels), character(0),
    info = paste0(
      "Model(s) that a dose does not reach under default rxSolve() arguments: ",
      paste(setdiff(dead, knownLinCmtDropModels), collapse = ", "),
      ". Try rxSolve(..., useLinCmt = FALSE); if that fixes it the model has ",
      "hit the ODE-to-linCmt conversion defect described at the top of this ",
      "file."
    )
  )
  expect_equal(
    setdiff(brokenByConversion, knownLinCmtDropModels), character(0),
    info = paste0(
      "Model(s) that fail to solve under default rxSolve() arguments but ",
      "solve with useLinCmt = FALSE: ",
      paste(setdiff(brokenByConversion, knownLinCmtDropModels), collapse = ", "),
      ". The conversion has renumbered the compartments out from under the ",
      "event table."
    )
  )
  expect_equal(
    setdiff(unsupported, probeUnsupportedModels), character(0),
    info = paste0(
      "Model(s) the probe could not solve either way: ",
      paste(setdiff(unsupported, probeUnsupportedModels), collapse = ", "),
      ". Either the model is broken or helper-solveProbe.R needs to learn how ",
      "to drive it -- work out which before adding it to probeUnsupportedModels."
    )
  )
})

test_that("rxode2's linCmt() optimisation never changes a model's solution", {
  # `useLinCmt = TRUE` is an optimisation: it must return the same numbers as
  # the ODE solve, or it is not an optimisation but a different model. This
  # catches both failure modes of the conversion -- the dead model and the
  # silently-wrong one that drops a transit chain but still produces plausible
  # concentrations.
  skip_on_cran()
  skip_if_not_installed("rxode2")

  db <- nlmixr2lib::modeldb
  names_ <- intersect(db$name, linCmtRiskCandidates())
  expect_gt(length(names_), 0)

  divergent <- character(0)
  for (n in names_) {
    lin <- probeSolveModel(n, db, useLinCmt = TRUE)
    ode <- probeSolveModel(n, db, useLinCmt = FALSE)
    if (!identical(lin$status, "ok") || !identical(ode$status, "ok")) next
    if (probeMaxRelDiff(lin, ode) > 1e-4) divergent <- c(divergent, n)
  }

  expect_equal(
    setdiff(divergent, knownLinCmtDropModels), character(0),
    info = paste0(
      "Model(s) that solve differently with and without rxode2's linCmt ",
      "conversion: ", paste(setdiff(divergent, knownLinCmtDropModels), collapse = ", "),
      ". The conversion has dropped a term the analytical solution cannot ",
      "represent; see the comment at the top of this file."
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

  db <- nlmixr2lib::modeldb
  for (nm in knownLinCmtDropModels) {
    expect_true(
      nm %in% db$name,
      info = paste0(nm, " is quarantined but is not in the registry; ",
                    "remove it from knownLinCmtDropModels.")
    )
    if (!nm %in% db$name) next
    lin <- probeSolveModel(nm, db, useLinCmt = TRUE)
    ode <- probeSolveModel(nm, db, useLinCmt = FALSE)
    stillBroken <- identical(ode$status, "ok") &&
      (!identical(lin$status, "ok") || probeMaxRelDiff(lin, ode) > 1e-4)
    expect_true(
      stillBroken,
      # nolint next: line_length_linter.
      label = paste0(nm, " now solves the same with and without rxode2's linCmt conversion, so it must be removed from knownLinCmtDropModels")
    )
  }
})
