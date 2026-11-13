# Machinery for solving every model in the registry at typical values.
#
# Used by test-modeldb-solve-gate.R.  It exists as a helper rather than
# inline in the test file so that the covariate-derivation rules -- the part
# most likely to need extending as new models land -- are documented in one
# place and can be exercised on their own.
#
# THE PROBE MUST NOT BE DEGENERATE.  Setting every covariate to 0 is the
# obvious shortcut and it is wrong: nearly every model in this registry
# carries an allometric term of the form `(WT / 70)^0.75`, so a zero
# covariate collapses the denominator and manufactures BOTH false zeros and
# absurd non-zeros (a probe that did exactly this reported max(Cc) = 1.35e18
# for one model and a solve failure for another; see
# refharvest/reports/oasweep_PMC12648368.md).  Every covariate value below is
# therefore derived, and the rule that produced it is recorded so a surprising
# probe result can be traced back to the value that caused it.

# Fallback typical values by unit string, used ONLY when neither the model's
# covariateData nor its model() body supplies a reference.  These are
# order-of-magnitude physiological placeholders, not claims about any
# particular paper: the gate they feed asks whether a dose reaches the system
# at all, which is insensitive to the exact value but very sensitive to a
# zero in a denominator.  Nothing here may be 0.
.probeUnitDefaults <- c(
  "kg" = 70,
  "g" = 70000,
  "lb" = 154,
  "kg/m^2" = 25,
  "cm" = 170,
  "m" = 1.7,
  "m^2" = 1.73,
  "years" = 40,
  "year" = 40,
  "y" = 40,
  "months" = 480,
  "month" = 480,
  "weeks" = 2080,
  "week" = 2080,
  "days" = 14600,
  "day" = 14600,
  "h" = 350000,
  "mL/min" = 100,
  "mL/min/1.73 m^2" = 90,
  "L/h" = 6,
  "mL/h" = 6000,
  "g/L" = 40,
  "g/dL" = 4,
  "mg/dL" = 1,
  "mg/L" = 10,
  "umol/L" = 80,
  "mmol/L" = 5,
  "nmol/L" = 100,
  "nM" = 100,
  "uM" = 1,
  "U/L" = 30,
  "IU/L" = 30,
  "U" = 30,
  "ng/mL" = 10,
  "ug/mL" = 1,
  "pg/mL" = 100,
  "mg" = 100,
  "ug" = 100,
  "g " = 1,
  "nmol" = 1,
  "umol" = 1,
  "mmol" = 1,
  "%" = 50,
  "fraction" = 0.5,
  "mg/day" = 100,
  "mg/kg" = 1,
  "mg/m^2" = 100,
  "(count)" = 1,
  "(binary)" = 0,
  "(categorical)" = 0
)

# Read the leading number out of a reference_category, which is written either
# as a bare number (`reference_category = 0`) or as an annotated string
# (`reference_category = "0 (male)"`).  Returns NA when there is no leading
# number, which is the signal to fall through to the next rule.
.probeLeadingNumber <- function(x) {
  if (is.null(x) || length(x) != 1L) {
    return(NA_real_)
  }
  if (is.numeric(x)) {
    return(if (is.finite(x)) as.numeric(x) else NA_real_)
  }
  if (!is.character(x)) {
    return(NA_real_)
  }
  m <- regmatches(x, regexpr("^\\s*[-+]?[0-9]*\\.?[0-9]+", x))
  if (!length(m)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(m))
}

# Variables a model() block sets to a plain number exactly once, such as
# `ref_wt <- 88`, as a named numeric vector. A name assigned more than once has
# no single value and is left out.
.probeConstantAssignments <- function(expr) {
  if (!is.call(expr) || !identical(expr[[1]], as.name("{"))) {
    return(numeric(0))
  }
  nms <- character(0)
  vals <- numeric(0)
  everyName <- character(0)
  for (i in seq_along(expr)[-1L]) {
    st <- expr[[i]]
    if (!is.call(st) || length(st) != 3L || !(identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))) {
      next
    }
    if (!is.name(st[[2]])) {
      next
    }
    nm <- as.character(st[[2]])
    everyName <- c(everyName, nm)
    rhs <- st[[3]]
    neg <- is.call(rhs) && length(rhs) == 2L && identical(rhs[[1]], as.name("-"))
    if (neg) {
      rhs <- rhs[[2]]
    }
    if (is.numeric(rhs) && length(rhs) == 1L && is.finite(rhs)) {
      nms <- c(nms, nm)
      vals <- c(vals, if (neg) -as.numeric(rhs) else as.numeric(rhs))
    }
  }
  once <- names(which(table(everyName) == 1L))
  keep <- nms %in% once
  stats::setNames(vals[keep], nms[keep])
}

# A number written literally, or a name `consts` gives a value; NA otherwise.
.probeResolveConstant <- function(x, consts) {
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) {
    return(as.numeric(x))
  }
  if (is.name(x) && as.character(x) %in% names(consts)) {
    return(consts[[as.character(x)]])
  }
  NA_real_
}

# Walk a model() block looking for the constant a covariate is normalised
# against.  Population-PK models almost always write their reference value
# into the equation -- `(WT / 70)^0.75`, `CRCL / 68`, `AGE - 40` -- or into a
# named constant used there (`ref_wt <- 88; ...(WT - ref_wt)`), so the author's
# own typical value is recoverable from the source rather than guessed at.
# Returns NA when no such constant appears.
# nolint start: cyclocomp_linter. One branch per expression shape; the
# branch count is the rule set, not accidental complexity.
.probeNormalisingWalk <- function(e, covName, hit, consts = numeric(0)) {
  if (is.finite(hit$value) || !is.call(e)) {
    return(invisible(NULL))
  }
  if (
    length(e) == 3L &&
      (identical(e[[1]], as.name("/")) || identical(e[[1]], as.name("-")) || identical(e[[1]], as.name("+")))
  ) {
    lhs <- e[[2]]
    rhs <- e[[3]]
    rv <- .probeResolveConstant(rhs, consts)
    lv <- .probeResolveConstant(lhs, consts)
    if (is.name(lhs) && as.character(lhs) == covName && is.finite(rv) && rv != 0) {
      hit$value <- rv
    } else if (
      identical(e[[1]], as.name("/")) && is.name(rhs) && as.character(rhs) == covName && is.finite(lv) && lv != 0
    ) {
      # `70 / WT` -- the reference sits on the other side
      hit$value <- lv
    }
  }
  for (i in seq_along(e)) {
    if (!is.null(e[[i]])) .probeNormalisingWalk(e[[i]], covName, hit, consts)
  }
  invisible(NULL)
}
# nolint end

.probeNormalisingConstant <- function(expr, covName) {
  # An environment rather than a closure variable: reference semantics mean the
  # recursion can record its hit without a superassignment, and the assignment
  # target is named right where the assignment happens.
  hit <- new.env(parent = emptyenv())
  hit$value <- NA_real_
  .probeNormalisingWalk(expr, covName, hit, .probeConstantAssignments(expr))
  hit$value
}

# Typical values for covariates that the unit table would put out of domain.
# The unit table reads every duration as an adult's age (480 months, 14600
# days, 350000 h), which suits AGE. For a postnatal or gestational age, a time
# since injury or an infusion lag, it lands far outside the range the model was
# fitted on, and a hematocrit given as a fraction falls through to the type
# default of 1, which leaves no plasma. Each value is keyed by the leading word
# of the declared units (see .probeUnitKey()), so a covariate gets the same
# typical value on whichever scale a model carries it.
.probeNameDefaults <- list(
  # Two weeks after birth.
  PNA = c(h = 336, days = 14, weeks = 2, months = 0.46, years = 0.038),
  # Term.
  PAGE = c(days = 280, weeks = 40, months = 9.2, years = 0.77),
  GA = c(days = 266, weeks = 38),
  DAY_AFTER_INJURY = c(days = 7),
  # Ten minutes.
  T_INFUSION_LAG = c(h = 0.17),
  HCT = c(percent = 45, fraction = 0.45),
  INS_BL = c(pmol = 60, mu = 10)
)

# The scale a free-text units string names, from its leading word: "days in
# this paper (canonical PNA is months)" is days, "% (volume fraction times
# 100)" is percent, "(fraction)" and "L/L" are fractions, "mU/L" and "uU/mL"
# are both micro-units. NA when nothing matches.
.probeUnitKey <- function(u) {
  if (!is.character(u) || length(u) != 1L) {
    return(NA_character_)
  }
  u <- tolower(trimws(u))
  keys <- c(
    percent = "^(%|percent)",
    fraction = "^[(]?fraction|^l/l",
    pmol = "^pmol",
    mu = "^(mu|uu)/",
    h = "^(h|hr|hours?)\\b",
    days = "^days?\\b",
    weeks = "^weeks?\\b",
    months = "^months?\\b",
    years = "^years?\\b"
  )
  hit <- names(keys)[vapply(keys, grepl, logical(1), x = u)]
  if (length(hit)) hit[[1]] else NA_character_
}

# Pull the model() block out of an rxUi so the normalising-constant search has
# something to walk.
.probeModelExpr <- function(ui) {
  lst <- try(ui$lstExpr, silent = TRUE)
  if (inherits(lst, "try-error") || is.null(lst)) {
    return(NULL)
  }
  as.call(c(list(as.name("{")), lst))
}

#' Derive a probe value for every covariate a model requires
#'
#' @param ui A decompressed rxUi.
#' @return A list with `values` (named numeric, one per covariate in
#'   `ui$allCovs`) and `why` (named character, the rule that produced each).
#' @noRd
probeCovariateValues <- function(ui) {
  covs <- ui$allCovs
  if (is.null(covs) || !length(covs)) {
    return(list(values = numeric(0), why = character(0)))
  }
  meta <- as.list(ui$meta)
  cd <- meta$covariateData
  body <- .probeModelExpr(ui)
  values <- stats::setNames(rep(NA_real_, length(covs)), covs)
  why <- stats::setNames(rep(NA_character_, length(covs)), covs)
  for (nm in covs) {
    entry <- if (is.list(cd) && !is.null(cd[[nm]]) && is.list(cd[[nm]])) cd[[nm]] else list()

    # 1. An explicit reference category.  This is the right answer for every
    #    binary and categorical covariate: it is the value the paper's typical
    #    subject has.
    v <- .probeLeadingNumber(entry$reference_category)
    if (is.finite(v)) {
      values[[nm]] <- v
      why[[nm]] <- paste0("covariateData$reference_category (", as.character(entry$reference_category)[1], ")")
      next
    }

    # 2. A numeric reference the metadata states. It comes before the
    #    constant inferred from the equations because a model can say which
    #    value is its reference; an equation can hold other constants of the
    #    same covariate (Tan_2025's `FPG - 3.5` is the HOMA glucose floor, not
    #    its reference).
    for (field in c("reference_value", "reference", "median", "typical_value")) {
      v <- .probeLeadingNumber(entry[[field]])
      if (is.finite(v)) break
    }
    if (is.finite(v)) {
      values[[nm]] <- v
      why[[nm]] <- paste0("covariateData$", field)
      next
    }

    # 3. The constant the model itself normalises against.  For a continuous
    #    covariate this is usually the population reference the authors
    #    fitted to, read out of their own equation.
    if (!is.null(body)) {
      v <- .probeNormalisingConstant(body, nm)
      if (is.finite(v)) {
        values[[nm]] <- v
        why[[nm]] <- paste0("normalising constant in model() (", format(v), ")")
        next
      }
    }

    # 4. A typical value for a covariate the unit table would get wrong.
    byScale <- .probeNameDefaults[[nm]]
    scale <- .probeUnitKey(entry$units)
    if (!is.null(byScale) && !is.na(scale) && scale %in% names(byScale)) {
      values[[nm]] <- unname(byScale[[scale]])
      why[[nm]] <- paste0("typical ", nm, " in ", scale)
      next
    }

    # 5. A physiological placeholder chosen from the declared units.
    u <- entry$units
    if (is.character(u) && length(u) == 1L && u %in% names(.probeUnitDefaults)) {
      values[[nm]] <- unname(.probeUnitDefaults[[u]])
      why[[nm]] <- paste0("unit default for '", u, "'")
      next
    }

    # 6. Last resort, by declared type.  1 rather than 0 for anything
    #    continuous, because a continuous covariate is far more likely to sit
    #    in a denominator than a binary one.
    ty <- entry$type
    ty <- if (is.character(ty) && length(ty) == 1L) ty else "unknown"
    values[[nm]] <- if (ty %in% c("binary", "categorical")) 0 else 1
    why[[nm]] <- paste0("type default for '", ty, "' (no reference available)")
  }
  list(values = values, why = why)
}

# Time grid, in the model's own time units.  A monoclonal antibody on a
# 0-48 grid looks indistinguishable from a dead model if the units are days.
probeTimeGrid <- function(ui, horizon = 1) {
  u <- as.list(ui$meta)$units
  tu <- if (is.list(u) && is.character(u$time) && length(u$time) == 1L) tolower(u$time) else ""
  end <- switch(
    tu,
    "min" = 720,
    "minute" = 720,
    "minutes" = 720,
    "h" = 48,
    "hr" = 48,
    "hour" = 48,
    "hours" = 48,
    "day" = 56,
    "days" = 56,
    "d" = 56,
    "week" = 24,
    "weeks" = 24,
    "month" = 12,
    "months" = 12,
    48
  )
  seq(0, end * horizon, length.out = 97L)
}

# How much longer the grid is made when a dose changes nothing on the default
# grid: an absorption lag can exceed it (Khwarg_2024_donepezil_im's IM depot
# lag is 235 h against a 48 h grid), which delays every effect past the last
# observation. A model that is genuinely dead stays dead on the longer grid.
.probeLongHorizon <- 10

# The compartment a model expects to be dosed into, taken from the registry
# rather than guessed.  Models that declare none are not probeable.
probeDoseCmt <- function(name, db) {
  d <- db$dosing[db$name == name]
  if (!length(d) || is.na(d) || !nzchar(d)) {
    return(NA_character_)
  }
  strsplit(d, ",")[[1]][1]
}

# Every compartment the probe doses: the registry's first dosing compartment,
# plus any other declared dosing compartment whose dose history the model
# reads. A split dose such as `f(depot) <- frel; f(depot2) <- 1 - frel` with
# `tad(depot2)` in a rate needs both records, or that tad() is NA and the
# derivatives using it are NaN. Compartments the model does not declare are
# never added, so a wrong `dosing` still shows up as non-finite states.
probeDoseCmts <- function(name, db, ui) {
  first <- probeDoseCmt(name, db)
  if (is.na(first)) {
    return(first)
  }
  declared <- strsplit(db$dosing[db$name == name], ",")[[1]]
  body <- .probeModelExpr(ui)
  if (is.null(body)) {
    stop("cannot read the model() block of ", name, " to find its dose-history reads", call. = FALSE)
  }
  calls <- nlmixr2lib:::.doseTimingCalls(body)
  read <- vapply(calls, `[[`, character(1), "cmt")
  unique(c(first, intersect(declared, read)))
}

# Build the event table.  Observation records need a compartment or dvid once
# a model has more than one endpoint, otherwise rxode2 raises
# "'dvid'->'cmt' or 'cmt' on observation record" -- which is an event-table
# error, not a model defect, and must not be reported as one.
#
# `form` selects how multi-endpoint observations are tagged.  No single form
# works for every model in the registry, so probeSolveModel() tries them in
# turn: naming the endpoint compartment is the most readable and works for the
# large majority, `dvid` indices cover models whose endpoint names are not
# compartment names, and untagged observations cover more. `first` tags every
# observation with the first endpoint, as the vignettes of models mixing an
# algebraic endpoint with ODE-state endpoints do (AlvarezJimenez_2016's
# `Cc <- 1e6 * linCmt()` with PD states); rxSolve still returns every
# endpoint's column.
#
# A model with a random effect on a grouping variable other than the subject
# (`eta ~ v | OCC`) gets that column too: every record at level 1 except the
# last, at level 2. See .probeGroupingVars().
probeEvents <- function(ui, doseCmt, amt, times, form = "cmt") {
  ev <- rxode2::et(amt = amt, cmt = doseCmt[[1]])
  for (cm in doseCmt[-1]) {
    ev <- rxode2::et(ev, amt = amt, cmt = cm)
  }
  eps <- ui$predDf$cond
  if (length(eps) > 1L && !identical(form, "plain")) {
    if (identical(form, "dvid")) {
      for (i in seq_along(eps)) {
        ev <- rxode2::et(ev, times, dvid = i)
      }
    } else if (identical(form, "first")) {
      ev <- rxode2::et(ev, times, cmt = eps[[1]])
    } else {
      for (ep in eps) {
        ev <- rxode2::et(ev, times, cmt = ep)
      }
    }
  } else {
    ev <- rxode2::et(ev, times)
  }
  grp <- .probeGroupingVars(ui)
  if (length(grp)) {
    ev <- as.data.frame(ev)
    if (!"id" %in% names(ev)) {
      ev$id <- 1L
    }
    for (g in grp) {
      ev[[g]] <- c(rep(1L, nrow(ev) - 1L), 2L)
    }
  }
  ev
}

.probeEventForms <- c("cmt", "dvid", "plain", "first")

# Grouping variables of random effects that vary within a subject, such as
# `etaiov ~ 0.1 | OCC`. rxSolve needs the column in the data even with the
# etas zeroed. It also needs two levels of it: with one, it stops with "The
# following parameter(s) are required for solving: THETA[1]". The zeroed etas
# make the level irrelevant to every value solved. A grouping variable that is
# also a covariate is supplied as one already.
.probeGroupingVars <- function(ui) {
  df <- ui$iniDf
  cond <- df$condition[!is.na(df$neta1)]
  setdiff(unique(cond[!is.na(cond)]), c("id", ui$allCovs))
}

# Every column worth comparing between two solves of the same model: the ODE
# states and the model's own computed quantities, minus the stochastic
# columns, which differ run to run and say nothing about whether the dose
# arrived.
.probeValueCols <- function(df) {
  drop <- c("id", "time", "evid", "amt", "cmt", "dvid", "ii", "addl", "dur", "rate", "ss", "sim", "ipredSim")
  setdiff(names(df), drop)
}

#' Largest relative disagreement between two solves of the same model
#'
#' Scaled per column so a 1e-4 threshold means the same thing for a
#' concentration in ng/mL and an amount in mg. Columns present in only one of
#' the two solves are ignored: the point is to compare the numbers the two
#' paths both claim to produce. Two paths that needed a different time grid
#' (see .probeSolveUi()) or a different observation form (see probeEvents())
#' disagree completely, and the result is 1: on the first the dose reached the
#' model at different times, and on the second the event table the ODE model
#' accepts no longer resolves after the conversion renumbered its
#' compartments.
#'
#' @param a,b Results of `probeSolveModel()`, both with `status == "ok"`.
#' @return A single non-negative number; 0 when the two solves agree.
#' @noRd
probeMaxRelDiff <- function(a, b) {
  if (!identical(a$horizon, b$horizon) || !identical(a$eventForm, b$eventForm)) {
    return(1)
  }
  cols <- intersect(a$cols, b$cols)
  if (!length(cols) || nrow(a$dosed) != nrow(b$dosed)) {
    return(0)
  }
  rel <- vapply(
    cols,
    function(cc) {
      x <- a$dosed[[cc]]
      y <- b$dosed[[cc]]
      sc <- suppressWarnings(max(abs(c(x, y)), na.rm = TRUE))
      if (!is.finite(sc) || sc == 0) {
        return(0)
      }
      d <- suppressWarnings(max(abs(x - y), na.rm = TRUE))
      if (!is.finite(d)) {
        return(0)
      }
      d / sc
    },
    numeric(1)
  )
  if (!length(rel) || all(is.na(rel))) {
    return(0)
  }
  max(rel, na.rm = TRUE)
}

# ODE solver error alone can put the two paths more than 1e-4 apart:
# Matcha_2025_amikacin's peripheral amount differs by 8e-3 at rxSolve()'s
# default tolerances and by 6e-5 at these. A conversion that changes the model
# does not shrink with the tolerance, so a divergent model is re-solved at
# these before it is reported.
.probeTightSolver <- list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e6)

#' Re-solve a model's two paths for a closer comparison
#'
#' Both paths are solved at tight tolerances (see .probeTightSolver), and, when
#' `form` is given, only in that observation form. Comparing the two paths in
#' one common form separates a conversion that changes the values from one
#' that only changes which observation records the converted model accepts.
#'
#' @param name Registry model name.
#' @param db The model registry.
#' @param form One of .probeEventForms, or NULL to let each path find its own.
#' @return probeMaxRelDiff() of the two solves, or NA when either path fails.
#' @noRd
probeRecheck <- function(name, db, form = NULL) {
  forms <- if (is.null(form)) .probeEventForms else form
  lin <- probeSolveModel(name, db, useLinCmt = TRUE, solverArgs = .probeTightSolver, forms = forms)
  ode <- probeSolveModel(name, db, useLinCmt = FALSE, solverArgs = .probeTightSolver, forms = forms)
  if (!identical(lin$status, "ok") || !identical(ode$status, "ok")) {
    return(NA_real_)
  }
  probeMaxRelDiff(lin, ode)
}

# The observation form to re-solve a divergent pair in: the one the converted
# model accepted, when the two paths needed different forms; NULL otherwise.
probeCommonForm <- function(lin, ode) {
  if (!identical(lin$eventForm, ode$eventForm)) lin$eventForm else NULL
}

# Screening for models at risk of the ODE-to-linCmt conversion defect ------
#
# rxode2 converts an ODE system to `linCmt()` only when it matches a
# depot/central/peripheral topology, and the conversion only loses information
# when some right-hand side carries a term that is not proportional to a
# state. Both conditions are visible in the source, so the candidate set can be
# computed from the model files in a couple of seconds -- no compiling, no
# solving. The screen deliberately errs towards including models the converter
# would reject (it checks the cheap structural preconditions, not the full
# topology), because a false positive costs two solves and a false negative
# costs a silently broken model.

.probeModelBlockOfFile <- function(path) {
  env <- new.env(parent = globalenv())
  fname <- sub("[.]R$", "", basename(path))
  sys.source(path, envir = env)
  if (!exists(fname, envir = env, inherits = FALSE)) {
    return(NULL)
  }
  b <- body(get(fname, envir = env))
  for (i in seq_along(b)) {
    node <- b[[i]]
    if (is.call(node) && identical(node[[1]], as.name("model"))) return(node[[2]])
  }
  NULL
}

# `d/dt(x)` parses as `/`(d, dt(x))
.probeIsDdt <- function(target) {
  is.call(target) &&
    identical(target[[1]], as.name("/")) &&
    identical(target[[2]], as.name("d")) &&
    is.call(target[[3]]) &&
    identical(target[[3]][[1]], as.name("dt"))
}

# Split an expression into signed additive terms, the same decomposition
# rxode2's `.collectAddTerms()` performs.
.probeAddTerms <- function(e, sign = 1L, acc = list()) {
  if (is.call(e) && length(e) == 3L && identical(e[[1]], as.name("+"))) {
    return(.probeAddTerms(e[[3]], sign, .probeAddTerms(e[[2]], sign, acc)))
  }
  if (is.call(e) && length(e) == 3L && identical(e[[1]], as.name("-"))) {
    return(.probeAddTerms(e[[3]], -sign, .probeAddTerms(e[[2]], sign, acc)))
  }
  if (is.call(e) && length(e) == 2L && identical(e[[1]], as.name("-"))) {
    return(.probeAddTerms(e[[2]], -sign, acc))
  }
  if (is.call(e) && length(e) == 2L && identical(e[[1]], as.name("("))) {
    return(.probeAddTerms(e[[2]], sign, acc))
  }
  c(acc, list(list(sign = sign, expr = e)))
}

# nolint start: cyclocomp_linter. Screens each model file for
# linCmt-convertible topology; the branch count is the rule set.
.probeScreenOneFile <- function(path) {
  mb <- .probeModelBlockOfFile(path)
  if (is.null(mb)) {
    stop("no model() block in ", path, call. = FALSE)
  }
  odes <- list()
  assigns <- list()
  for (i in seq_along(mb)) {
    st <- mb[[i]]
    isAssign <- is.call(st) &&
      length(st) >= 3 &&
      (identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))
    if (!isAssign) {
      next
    }
    if (.probeIsDdt(st[[2]])) {
      odes[[length(odes) + 1L]] <- list(cmt = as.character(st[[2]][[3]][[2]]), rhs = st[[3]])
    } else if (is.name(st[[2]])) {
      assigns[[length(assigns) + 1L]] <- st[[3]]
    }
  }
  # The converter handles at most four compartments.
  if (!length(odes) || length(odes) > 4L) {
    return(FALSE)
  }
  states <- vapply(odes, function(o) o$cmt, character(1))
  # It also needs an output line of the form `var <- <state> / <expr>`.
  hasOut <- FALSE
  for (r in assigns) {
    if (
      is.call(r) &&
        length(r) == 3L &&
        identical(r[[1]], as.name("/")) &&
        is.name(r[[2]]) &&
        as.character(r[[2]]) %in% states
    ) {
      hasOut <- TRUE
      break
    }
  }
  if (!hasOut) {
    return(FALSE)
  }
  # Finally, something on a right-hand side that references no state at all:
  # an exogenous input the analytical solution cannot carry.
  for (o in odes) {
    for (tm in .probeAddTerms(o$rhs)) {
      if (is.numeric(tm$expr) && length(tm$expr) == 1L && tm$expr == 0) {
        next
      }
      if (!any(all.vars(tm$expr) %in% states)) return(TRUE)
    }
  }
  FALSE
}

#' Registry models structurally at risk of losing a term to `linCmt()`
#'
#' Enumerates every model file shipped in `inst/modeldb`; nothing is
#' hand-listed, so a model added later is screened without touching this code.
#'
#' @return Character vector of model names.
#' @noRd
# nolint end
linCmtRiskCandidates <- function() {
  files <- .probeModelFiles()
  keep <- vapply(files, .probeScreenOneFile, logical(1))
  unname(sub("[.]R$", "", basename(files[keep])))
}

# Every model file shipped in the installed modeldb/ directory. Finding none is
# an error: an empty list would leave the gate's screened scope empty, and a
# test that only checks the scope is non-empty could still pass on the other
# screen.
.probeModelFiles <- function() {
  root <- system.file("modeldb", package = "nlmixr2lib")
  files <- if (nzchar(root)) list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE) else character(0)
  if (!length(files)) {
    stop("no model files found under the installed modeldb/ directory", call. = FALSE)
  }
  files
}

# Does this model file's model() block read the dose history of a named
# compartment (`tad(depot2)`, `podo(stomach)`, ...)? Uses the convention
# checker's own walker, so the screen and the rule agree on what counts. A
# file that fails to load, or has no model() block, is an error rather than a
# model that reads nothing: dropped silently, it would leave the gate's
# default scope.
.probeReadsDoseHistory <- function(path) {
  mb <- .probeModelBlockOfFile(path)
  if (is.null(mb)) {
    stop("no model() block in ", path, call. = FALSE)
  }
  length(nlmixr2lib:::.doseTimingCalls(mb)) > 0L
}

#' Registry models that read a named compartment's dose history
#'
#' These are the models a wrong dosing compartment breaks: the probe doses only
#' the registry's `dosing` compartments (see probeDoseCmts()), and a
#' dose-history function of any other compartment is NA, which makes the
#' derivatives that use it NaN.
#' Enumerated from the model files like `linCmtRiskCandidates()`.
#'
#' @return Character vector of model names.
#' @noRd
doseTimingCandidates <- function() {
  files <- .probeModelFiles()
  keep <- vapply(files, .probeReadsDoseHistory, logical(1))
  unname(sub("[.]R$", "", basename(files[keep])))
}

# One solve of an already-prepared model. `amt = 0` gives the no-dose control,
# which keeps the two solves identical in every other respect -- same event
# grid, same covariates, same solver settings. `solverArgs` are passed to
# rxSolve() as they are (`method`, `atol`, `rtol`, ...).
.probeSolveOnce <- function(ui, doseCmt, amt, times, covValues, useLinCmt, form, solverArgs = list()) {
  ev <- probeEvents(ui, doseCmt, amt, times, form)
  args <- list(object = ui, events = ev, returnType = "data.frame", addDosing = FALSE, useLinCmt = useLinCmt)
  if (length(covValues)) {
    args$params <- covValues
  }
  args <- c(args, solverArgs)
  suppressWarnings(suppressMessages(do.call(rxode2::rxSolve, args)))
}

# Find an observation-record form this model accepts, and return the dosed
# solve made with it.  Returns a `try-error` carrying the FIRST form's message
# when none works, because that message names the real obstacle.
.probeSolveWithAnyForm <- function(
  ui,
  doseCmt,
  amt,
  times,
  covValues,
  useLinCmt,
  solverArgs = list(),
  forms = .probeEventForms
) {
  first <- NULL
  for (form in forms) {
    r <- try(.probeSolveOnce(ui, doseCmt, amt, times, covValues, useLinCmt, form, solverArgs), silent = TRUE)
    if (!inherits(r, "try-error")) {
      return(list(result = r, form = form))
    }
    if (is.null(first)) {
      first <- r
    }
    # Only the multi-endpoint tagging varies between forms; for a
    # single-endpoint model all three build the same table, so retrying is
    # pointless. NROW() because a model with no residual-error line has a NULL
    # predDf.
    if (NROW(ui$predDf) <= 1L) break
  }
  list(result = first, form = NA_character_)
}

# Solvers tried in turn, as extra rxSolve() arguments. rxode2's default
# liblsoda returns every state NaN, with no warning, for the stiff time-kill
# systems of Ho_2025_betalactam_pbp78_CRAB and the four
# Mahadevan_2026_polymyxinB_meropenem_fosfomycin models, and cannot solve
# Svensson_2016_rifampicin at all; dop853 solves all three, and their vignettes
# already pass `method`. A model that is genuinely broken stays non-finite, or
# unsolvable, under both.
.probeSolverMethods <- list(default = list(), dop853 = list(method = "dop853"))

# The dosed solve, with the first solver that gives finite states. When none
# does, the default solver's attempt is returned, so a failure is reported as
# rxSolve() gives it under its defaults.
.probeSolveDosed <- function(
  ui,
  doseCmt,
  amt,
  times,
  covValues,
  useLinCmt,
  solverArgs = list(),
  forms = .probeEventForms
) {
  first <- NULL
  for (m in names(.probeSolverMethods)) {
    args <- utils::modifyList(.probeSolverMethods[[m]], solverArgs)
    attempt <- .probeSolveWithAnyForm(ui, doseCmt, amt, times, covValues, useLinCmt, args, forms)
    attempt$method <- m
    attempt$solverArgs <- args
    if (!inherits(attempt$result, "try-error") && !length(probeNonFiniteStates(ui, attempt$result))) {
      return(attempt)
    }
    if (is.null(first)) {
      first <- attempt
    }
  }
  first
}

#' ODE states that a solve leaves NaN or infinite
#'
#' Only `d/dt()` states count. A timing output such as `tad()`, a time since
#' inoculation or a trigger time is NA before its first event by design, but a
#' state is an amount: a non-finite one means a derivative was non-finite,
#' which newer rxode2 builds refuse outright by failing the subject.
#' States a `linCmt()` conversion removed from the output are not looked for.
#'
#' @param ui The rxUi that was solved.
#' @param solved The solve, as a data frame.
#' @return Character vector of state names, empty when every state is finite.
#' @noRd
probeNonFiniteStates <- function(ui, solved) {
  states <- intersect(ui$state, names(solved))
  bad <- colSums(!is.finite(as.matrix(solved[, states, drop = FALSE]))) > 0
  states[bad]
}

#' Solve one registry model at typical values, with and without a dose
#'
#' The no-dose solve is the control.  Asking "is the solve identically zero"
#' misses a model whose PK pathway died but whose PD baseline is non-zero;
#' asking "did the dose change anything" cannot.  It is also a difference
#' rather than a ratio, so it cannot pass by dividing one zero by another.
#'
#' @param name Registry model name.
#' @param db The model registry (`modeldb`), passed in so the caller can
#'   subset it.
#' @param useLinCmt Passed through to `rxSolve()`; `FALSE` disables rxode2's
#'   ODE-to-`linCmt()` optimisation.
#' @param solverArgs Extra `rxSolve()` arguments for every solve, such as
#'   tighter `atol` and `rtol`.
#' @return A list with `status` (one of "ok", "skip", "error"), and on success
#'   `doseEffect` (the largest absolute difference the dose made to any
#'   model quantity), `nonFiniteStates` (see `probeNonFiniteStates()`), the
#'   `dosed` data frame, the solver (`method`) and grid (`horizon`) used, and
#'   the covariate provenance.
#' @noRd
probeSolveModel <- function(name, db, useLinCmt = TRUE, amt = 100, solverArgs = list(), forms = .probeEventForms) {
  out <- list(name = name, status = "error", message = NA_character_)
  ui <- try(rxode2::rxUiDecompress(rxode2::as.rxUi(nlmixr2lib::readModelDb(name))), silent = TRUE)
  if (inherits(ui, "try-error")) {
    out$message <- paste("could not load:", conditionMessage(attr(ui, "condition")))
    return(out)
  }
  doseCmt <- probeDoseCmts(name, db, ui)
  if (is.na(doseCmt[[1]])) {
    out$status <- "skip"
    out$message <- "no dosing compartment declared in the registry"
    return(out)
  }
  c(list(name = name), .probeSolveUi(ui, doseCmt, useLinCmt, amt, solverArgs, forms))
}

# probeSolveModel() for a model already loaded, so fixture models can be
# probed without a registry entry. Probes on the default grid, and again on a
# grid .probeLongHorizon times longer when the dose changed nothing there.
.probeSolveUi <- function(ui, doseCmt, useLinCmt = TRUE, amt = 100, solverArgs = list(), forms = .probeEventForms) {
  out <- list(status = "error", message = NA_character_)
  ui <- try(.probeZeroRe(ui), silent = TRUE)
  if (inherits(ui, "try-error")) {
    out$message <- paste("zeroRe() failed:", conditionMessage(attr(ui, "condition")))
    return(out)
  }
  cov <- probeCovariateValues(ui)
  res <- .probeSolveGrid(ui, doseCmt, amt, probeTimeGrid(ui), cov, useLinCmt, solverArgs, forms)
  res$horizon <- 1
  if (identical(res$status, "ok") && !isTRUE(res$doseEffect > 0)) {
    longer <- .probeSolveGrid(ui, doseCmt, amt, probeTimeGrid(ui, .probeLongHorizon), cov, useLinCmt, solverArgs, forms)
    if (identical(longer$status, "ok") && isTRUE(longer$doseEffect > 0)) {
      res <- longer
      res$horizon <- .probeLongHorizon
    }
  }
  res
}

# zeroRe() for the kinds of random effect the model has. zeroRe() warns "No
# omega parameters in the model" (or sigma) for a kind that is absent.
.probeZeroRe <- function(ui) {
  df <- ui$iniDf
  which <- c(if (any(!is.na(df$neta1))) "omega", if (any(!is.na(df$err))) "sigma")
  if (!length(which)) {
    return(ui)
  }
  rxode2::zeroRe(ui, which = which)
}

# The dosed and no-dose solves on one time grid, and what the dose changed.
.probeSolveGrid <- function(ui, doseCmt, amt, times, cov, useLinCmt, solverArgs, forms = .probeEventForms) {
  out <- list(status = "error", message = NA_character_, covWhy = cov$why)
  attempt <- .probeSolveDosed(ui, doseCmt, amt, times, cov$values, useLinCmt, solverArgs, forms)
  dosed <- attempt$result
  if (inherits(dosed, "try-error")) {
    out$message <- paste("solve failed:", conditionMessage(attr(dosed, "condition")))
    return(out)
  }
  out$eventForm <- attempt$form
  out$method <- attempt$method
  undosed <- try(
    .probeSolveOnce(ui, doseCmt, 0, times, cov$values, useLinCmt, attempt$form, attempt$solverArgs),
    silent = TRUE
  )
  if (inherits(undosed, "try-error")) {
    out$message <- paste("control solve failed:", conditionMessage(attr(undosed, "condition")))
    return(out)
  }
  cols <- intersect(.probeValueCols(dosed), .probeValueCols(undosed))
  cols <- cols[vapply(cols, function(cc) is.numeric(dosed[[cc]]), logical(1))]
  if (!length(cols) || nrow(dosed) != nrow(undosed)) {
    out$status <- "skip"
    out$message <- "no comparable numeric output columns"
    return(out)
  }
  diffs <- vapply(
    cols,
    function(cc) {
      d <- abs(dosed[[cc]] - undosed[[cc]])
      if (all(is.na(d))) {
        return(NA_real_)
      }
      max(d, na.rm = TRUE)
    },
    numeric(1)
  )
  out$status <- "ok"
  out$doseCmt <- doseCmt
  out$doseEffect <- if (all(is.na(diffs))) NA_real_ else max(diffs, na.rm = TRUE)
  out$nonFiniteStates <- probeNonFiniteStates(ui, dosed)
  out$dosed <- dosed
  out$cols <- cols
  out$covWhy <- cov$why
  out$covValues <- cov$values
  out
}
