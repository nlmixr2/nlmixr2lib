# Unit tracking for rxode2/nlmixr2 models: one walker over `units` objects, one
# solver for the single unknown in a term, and a fixed-point loop.
#
# A symbol's unit is a `units` object made by `units::set_units(1, "<unit>")`.
# The walker (`.unitWalk`) turns an expression into such an object, or NULL
# when it depends on a symbol whose unit is not known yet; the solver
# (`.unitSolve`) works a single unknown symbol out of a term whose target unit
# is known. Both are pure functions over a named list `u` of units objects, so
# each can be tested on a quoted expression alone. `checkUnits()` runs them to
# a fixed point and validates every line; `addUnits()` applies the result.

# ---- units helpers ----------------------------------------------------------

.unitsAvailable <- function() {
  requireNamespace("units", quietly = TRUE)
}

.unitError <- function(message) {
  structure(class = c("nlmixr2libUnitError", "error", "condition"), list(message = message, call = NULL))
}

# Units udunits lacks but pharmacometric models use; "" installs a base unit.
# Installing is process-global and a repeated bare install errors, hence the
# convertible-to-itself guard.
.unitCustom <- c(IU = "", U = "", M = "mol/L", cells = "", CFU = "", copies = "", molecules = "")

.unitInstall <- function() {
  if (!.unitsAvailable()) {
    stop(.unitError(paste0(
      "checkUnits() and addUnits() need the 'units' package and its udunits-2 ",
      "system library: install.packages(\"units\")"
    )))
  }
  for (nm in names(.unitCustom)) {
    if (units::ud_are_convertible(nm, nm)) {
      next
    }
    if (nzchar(.unitCustom[[nm]])) units::install_unit(nm, .unitCustom[[nm]]) else units::install_unit(nm)
  }
  invisible(TRUE)
}

.unitlessSpellings <- c("", "1", "unitless", "dimensionless", "fraction", "none", "n/a")

# A unit string as a units object with value 1. "unitless" is this package's
# spelling of the dimensionless unit and becomes udunits' "1" here; installing
# it as a unit instead makes udunits name every bare 1 with it and misjudge
# convertibility.
.unitOf <- function(s) {
  .unitInstall()
  checkmate::assertString(s)
  s <- trimws(s)
  s <- gsub(intToUtf8(0xB5), "u", s, fixed = TRUE)
  s <- gsub(intToUtf8(0x3BC), "u", s, fixed = TRUE)
  if (tolower(s) %in% .unitlessSpellings) {
    s <- "unitless"
  }
  aliases <- .nlmixr2libConventionsStatic$unitSpellingAliases
  for (from in names(aliases)) {
    s <- gsub(paste0("(?<![A-Za-z])", from, "(?![A-Za-z])"), aliases[[from]], s, perl = TRUE)
  }
  # units only has integer powers and reads kg^0.75 as the scalar 75
  if (grepl("\\^\\s*[-+]?[0-9]*\\.[0-9]", s) || grepl("**", s, fixed = TRUE)) {
    stop(.unitError(sprintf("'%s': the units library cannot represent a fractional exponent", s)))
  }
  s <- gsub("(?<![A-Za-z])unitless(?![A-Za-z])", "1", s, perl = TRUE)
  tryCatch(
    withCallingHandlers(
      units::set_units(1, s, mode = "standard"),
      warning = function(w) stop(conditionMessage(w))
    ),
    error = function(e) stop(.unitError(sprintf("'%s' is not a unit: %s", s, conditionMessage(e))))
  )
}

.unitDimless <- function() {
  units::set_units(1, "1", mode = "standard")
}

.unitIsDimless <- function(x) {
  units::ud_are_convertible(units::deparse_unit(x), "1")
}

# The unit as units prints it, with identical tokens cancelled (`L*h/h` is
# `L`); "unitless" for the dimensionless unit. checkUnits() runs with units'
# own simplification off because that would also fold `mg/kg` into a bare
# number, and a dose per body weight is its own dimension here.
.unitStr <- function(x) {
  s <- units(x)
  num <- s$numerator
  den <- s$denominator
  for (tok in intersect(num, den)) {
    n <- min(sum(num == tok), sum(den == tok))
    num <- num[-which(num == tok)[seq_len(n)]]
    den <- den[-which(den == tok)[seq_len(n)]]
  }
  if (length(num) == 0 && length(den) == 0) {
    return("unitless")
  }
  tokens <- c(num, if (length(den) > 0) paste0(den, "-1"))
  as.character(units(units::set_units(1, paste(tokens, collapse = " "), mode = "standard")))
}

# The same unit with the scale dropped.
.unitOne <- function(x) {
  x[] <- 1
  x
}

.unitConv <- function(a, b) {
  units::ud_are_convertible(units::deparse_unit(a), units::deparse_unit(b))
}

# The number one unit of `a` (scale included) is in `b` (scale included).
.unitFactor <- function(a, b) {
  as.numeric(units::set_units(a, units(b), mode = "standard")) / as.numeric(b)
}

# ---- the walker ---------------------------------------------------------------

.unitTranscendental <- c(
  "exp",
  "expm1",
  "log",
  "log10",
  "log2",
  "log1p",
  "sqrt",
  "sin",
  "cos",
  "tan",
  "asin",
  "acos",
  "atan",
  "sinh",
  "cosh",
  "tanh",
  "expit",
  "logit",
  "probit",
  "probitInv",
  "gammafn",
  "lgamma",
  "digamma",
  "gamma",
  "beta",
  "lbeta",
  "factorial",
  "lfactorial",
  "pnorm",
  "qnorm",
  "erf",
  "erfc"
)
.unitPassThrough <- c("abs", "floor", "ceiling", "round", "trunc", "signif", "max", "min", "ifelse")
.unitLogical <- c("<", ">", "<=", ">=", "==", "!=", "&&", "||", "&", "|", "!")

.isLiteral <- function(e) {
  is.numeric(e) && length(e) == 1L
}

.uw <- function(unit = NULL, issues = character()) {
  list(unit = unit, issues = issues)
}

# The unit of an expression given the units `u` of its symbols: a list with the
# units object (NULL when unknown) and the issues found along the way. Rules:
# a literal is a bare number, except that one added to or subtracted from a
# quantity is a value in that quantity's unit and one dividing a quantity is a
# reference value (`WT / 70` is unitless); transcendental functions,
# comparisons and integer powers of a unitless base are unitless; a
# dimensioned base to a symbolic power (`Cc^hill`) is unknown.
.unitWalk <- function(e, u) {
  if (.isLiteral(e)) {
    return(.uw(.unitDimless()))
  }
  if (is.name(e)) {
    return(.uw(u[[as.character(e)]]))
  }
  if (!is.call(e) || !is.name(e[[1]])) {
    return(.uw())
  }
  fn <- as.character(e[[1]])
  a <- as.list(e)[-1]
  if (fn == "(" || (fn %in% c("+", "-") && length(a) == 1L)) {
    return(.unitWalk(a[[1]], u))
  }
  parts <- lapply(a, .unitWalk, u = u)
  issues <- unlist(lapply(parts, `[[`, "issues"))
  x <- parts[[1]]$unit
  y <- if (length(parts) > 1L) parts[[2]]$unit
  known <- function(i) !.isLiteral(a[[i]]) && !is.null(parts[[i]]$unit)
  unit <- switch(
    fn,
    "*" = if (!is.null(x) && !is.null(y)) x * y,
    "/" = if (.isLiteral(a[[2]])) {
      .unitDimless()
    } else if (!is.null(x) && !is.null(y)) {
      x / y
    },
    "+" = ,
    "-" = {
      if (.isLiteral(a[[1]]) || .isLiteral(a[[2]])) {
        other <- if (.isLiteral(a[[1]])) y else x
        if (is.null(other)) .unitDimless() else .unitOne(other)
      } else if (!is.null(x) && !is.null(y)) {
        if (!.unitConv(x, y)) {
          issues <- c(issues, sprintf("'%s' adds %s to %s", deparse1(e), .unitStr(x), .unitStr(y)))
        }
        .unitOne(x)
      } else if (!is.null(x)) {
        .unitOne(x)
      } else if (!is.null(y)) {
        .unitOne(y)
      }
    },
    "^" = {
      k <- a[[2]]
      if (is.null(x)) {
        NULL
      } else if (.isLiteral(k) && abs(k - round(k)) < 1e-12) {
        x^as.integer(round(k))
      } else if (.unitIsDimless(x)) {
        .unitDimless()
      } else if (.isLiteral(k)) {
        issues <- c(
          issues,
          sprintf(
            "'%s' raises %s to the power %s; normalise by a reference value first",
            deparse1(e),
            .unitStr(x),
            format(k)
          )
        )
        NULL
      }
    },
    {
      if (fn %in% c(.unitLogical, .unitTranscendental)) {
        .unitDimless()
      } else if (fn %in% .unitPassThrough) {
        # the first known operand's unit; with none known, a bare number among
        # the operands makes the result unitless (`ifelse(t < 8, 1, tfrac)`)
        idx <- if (fn == "ifelse") seq_along(a)[-1] else seq_along(a)
        knownIdx <- Filter(known, idx)
        if (length(knownIdx) > 0) {
          .unitOne(parts[[knownIdx[[1]]]]$unit)
        } else if (any(vapply(a[idx], .isLiteral, logical(1)))) {
          .unitDimless()
        }
      }
    }
  )
  .uw(unit, issues)
}

# The ini() parameters a back-transform names: `exp(lcl + etalcl)`, or a
# product of such a call with other factors, is a parameter on its estimation
# scale, so the assigned variable keeps its own unit and these are unitless.
.unitBackTransform <- function(e, iniNames) {
  e <- .stripParens(e)
  if (!is.call(e) || !is.name(e[[1]])) {
    return(character())
  }
  fn <- as.character(e[[1]])
  if (fn %in% c("*", "/") && length(e) == 3L) {
    return(c(.unitBackTransform(e[[2]], iniNames), if (fn == "*") .unitBackTransform(e[[3]], iniNames)))
  }
  if (!fn %in% c("exp", "expit", "probitInv") || length(e) != 2L) {
    return(character())
  }
  syms <- .modelLineSymbols(list(e[[2]]))
  if (length(syms) > 0 && all(syms %in% iniNames)) syms else character()
}

# ---- constraints ----------------------------------------------------------------

# One record per assignment in the model block, with the index path of the
# assignment inside `lstExpr[[line]]` so an `if` branch can be edited in place.
.unitConstraints <- function(lstExpr) {
  acc <- new.env(parent = emptyenv())
  acc$rows <- list()
  for (i in seq_along(lstExpr)) {
    .unitConstraintsFrom(lstExpr[[i]], i, integer(), acc)
  }
  acc$rows
}

.unitConstraintsFrom <- function(e, line, path, acc) {
  if (!is.call(e) || !is.name(e[[1]])) {
    return(invisible(NULL))
  }
  fn <- as.character(e[[1]])
  if (fn %in% c("{", "if")) {
    if (fn == "if") {
      acc$rows[[length(acc$rows) + 1L]] <- list(
        target = sprintf("line %d", line),
        kind = "check",
        line = line,
        path = c(path, 2L),
        rhs = e[[2]]
      )
    }
    for (j in seq_along(e)[-(1:(if (fn == "if") 2L else 1L))]) {
      .unitConstraintsFrom(e[[j]], line, c(path, j), acc)
    }
    return(invisible(NULL))
  }
  if (!fn %in% c("<-", "=")) {
    return(invisible(NULL))
  }
  lhs <- e[[2]]
  row <- NULL
  if (is.name(lhs)) {
    row <- list(target = as.character(lhs), kind = "sym")
  } else if (is.call(lhs) && is.name(lhs[[1]])) {
    lfn <- as.character(lhs[[1]])
    if (lfn == "/" && length(lhs) == 3L && is.call(lhs[[3]]) && identical(lhs[[3]][[1]], as.name("dt"))) {
      cmt <- as.character(lhs[[3]][[2]])
      row <- list(target = sprintf("d/dt(%s)", cmt), kind = "ddt", cmt = cmt)
    } else if (lfn %in% c("f", "F", "rate", "dur", "alag", "lag") && length(lhs) == 2L) {
      kind <- switch(lfn, f = "f", F = "f", rate = "rate", "time")
      row <- list(target = deparse1(lhs), kind = kind, cmt = as.character(lhs[[2]]))
    }
  }
  if (!is.null(row)) {
    acc$rows[[length(acc$rows) + 1L]] <- c(row, list(line = line, path = path, rhs = e[[3]]))
  }
  invisible(NULL)
}

# The unit a constraint's target must have, or NULL when not yet known.
.unitTarget <- function(cs, u) {
  switch(
    cs$kind,
    sym = u[[cs$target]],
    ddt = ,
    rate = if (!is.null(u[[cs$cmt]]) && !is.null(u$time)) u[[cs$cmt]] / u$time,
    f = ,
    check = .unitDimless(),
    time = u$time
  )
}

# Top-level additive terms of an expression, each with its sign.
.unitTerms <- function(e) {
  e <- .stripParens(e)
  if (is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("+", "-") && length(e) == 3L) {
    rhs <- .unitTerms(e[[3]])
    if (identical(e[[1]], as.name("-"))) {
      rhs <- lapply(rhs, function(t) call("-", t))
    }
    return(c(.unitTerms(e[[2]]), rhs))
  }
  list(e)
}

# ---- the solver ------------------------------------------------------------------

# Multiplicative leaves of a term with their integer powers. `x / 70` stays one
# leaf so the walker's reference-value rule applies to it.
.unitLeaves <- function(e, power = 1L) {
  e <- .stripParens(e)
  if (is.call(e) && is.name(e[[1]]) && length(e) == 3L) {
    fn <- as.character(e[[1]])
    if (fn == "*") {
      return(c(.unitLeaves(e[[2]], power), .unitLeaves(e[[3]], power)))
    }
    if (fn == "/" && !.isLiteral(e[[3]])) {
      return(c(.unitLeaves(e[[2]], power), .unitLeaves(e[[3]], -power)))
    }
    if (fn == "^" && .isLiteral(e[[3]]) && abs(e[[3]] - round(e[[3]])) < 1e-12) {
      return(.unitLeaves(e[[2]], power * as.integer(round(e[[3]]))))
    }
  }
  if (is.call(e) && is.name(e[[1]]) && length(e) == 2L && as.character(e[[1]]) %in% c("-", "+")) {
    return(.unitLeaves(e[[2]], power))
  }
  list(list(expr = e, power = power))
}

# The product of the bare numeric literals in a term (`central / vc * 1000`
# gives 1000): a conversion the model already carries.
.unitLiteralFactor <- function(term) {
  prod(vapply(.unitLeaves(term), function(lf) if (.isLiteral(lf$expr)) as.numeric(lf$expr)^lf$power else 1, numeric(1)))
}

# The unit of the one unknown symbol in a term whose target unit is known:
# `list(name, unit)`, or NULL when the term has no single solvable unknown.
.unitSolve <- function(term, target, u, defaults) {
  unknown <- list()
  rest <- .unitDimless()
  for (lf in .unitLeaves(term)) {
    r <- .unitWalk(lf$expr, u)$unit
    if (is.null(r)) {
      if (!is.name(lf$expr)) {
        return(NULL)
      }
      unknown[[length(unknown) + 1L]] <- lf
    } else {
      rest <- rest * r^lf$power
    }
  }
  if (length(unknown) != 1L || abs(unknown[[1]]$power) != 1L) {
    return(NULL)
  }
  q <- target / rest
  if (unknown[[1]]$power == -1L) {
    q <- .unitDimless() / q
  }
  list(name = as.character(unknown[[1]]$expr), unit = .unitSnap(q, defaults))
}

# A volume or clearance is written in the model's default spelling; anything
# else keeps what the equations gave, with the scale dropped.
.unitSnap <- function(q, defaults) {
  for (cand in defaults) {
    if (.unitConv(q, cand)) {
      return(cand)
    }
  }
  .unitOne(q)
}

# The volume and clearance defaults: L and L/<time>, or mL/kg and mL/kg/<time>
# for a dose per body weight.
.unitDefaults <- function(time, dose) {
  weight <- !is.null(dose) && "kg" %in% units(dose)$denominator
  scheme <- .nlmixr2libConventionsStatic$unitDefaults[[if (weight) "weight" else "mass"]]
  out <- list(.unitOf(scheme[["volume"]]))
  if (!is.null(time)) {
    out <- c(list(.unitOf(sub("<time>", .unitStr(time), scheme[["clearance"]], fixed = TRUE))), out)
  }
  out
}

# ---- inference and validation ---------------------------------------------------------

# What a residual-error parameter represents: the endpoint's unit for an
# additive SD, unitless for a fraction or a transformed-scale SD, NA when it
# cannot be written (a power-law SD carries var^(1-c)).
.unitResidual <- c(
  add = "endpoint",
  dnorm = "endpoint",
  prop = "unitless",
  propT = "unitless",
  propF = "unitless",
  pow = NA,
  powT = NA,
  powF = NA,
  pow2 = "unitless",
  powT2 = "unitless",
  powF2 = "unitless",
  lnorm = "unitless",
  logn = "unitless",
  dlnorm = "unitless",
  logitNorm = "unitless",
  probitNorm = "unitless",
  boxCox = "unitless",
  yeoJohnson = "unitless",
  t = "unitless",
  dt = "unitless",
  cauchy = "unitless",
  dcauchy = "unitless",
  dpois = "unitless",
  dbinom = "unitless",
  dbeta = "unitless",
  dgeom = "unitless",
  dnbinom = "unitless",
  dexp = "unitless"
)

# Forward and backward inference to a fixed point. `u` only ever gains
# entries, so a pass that adds none ends the loop.
.unitInfer <- function(cs, u, iniNames, residuals, defaults) {
  repeat {
    n <- length(u)
    for (i in seq_len(nrow(residuals))) {
      nm <- residuals$name[i]
      kind <- unname(.unitResidual[residuals$err[i]])
      if (!is.null(u[[nm]]) || is.na(kind)) {
        next
      }
      if (kind == "unitless") {
        u[[nm]] <- .unitDimless()
      } else if (!is.null(u[[residuals$endpoint[i]]])) {
        u[[nm]] <- .unitOne(u[[residuals$endpoint[i]]])
      }
    }
    for (c in cs) {
      if (c$kind == "check") {
        next
      }
      target <- .unitTarget(c, u)
      for (term in .unitTerms(c$rhs)) {
        bt <- .unitBackTransform(term, iniNames)
        for (p in setdiff(bt, names(u))) {
          u[[p]] <- .unitDimless()
        }
        if (length(bt) > 0 && is.null(target)) {
          next # the assigned variable is inferred from its own use
        }
        r <- .unitWalk(term, u)$unit
        if (!is.null(target) && is.null(r)) {
          s <- .unitSolve(term, target, u, defaults)
          if (!is.null(s)) {
            u[[s$name]] <- s$unit
          }
        } else if (is.null(target) && !is.null(r)) {
          if (c$kind == "sym") {
            u[[c$target]] <- .unitSnap(r, defaults)
          } else if (c$kind == "ddt" && !is.null(u$time)) {
            u[[c$cmt]] <- .unitSnap(r * u$time, defaults)
          }
          target <- .unitTarget(c, u)
        }
      }
    }
    if (length(u) == n) {
      break
    }
  }
  u
}

# Every term against its target: an issue when not convertible, a conversion
# when convertible but off by a factor the model does not already carry.
.unitValidate <- function(cs, u, iniNames) {
  issues <- list()
  conv <- list()
  for (c in cs) {
    target <- .unitTarget(c, u)
    terms <- if (c$kind == "check") list(c$rhs) else .unitTerms(c$rhs)
    for (i in seq_along(terms)) {
      w <- .unitWalk(terms[[i]], u)
      if (length(w$issues) > 0) {
        issues[[c$target]] <- c(issues[[c$target]], w$issues)
      }
      if (c$kind == "check" || is.null(target) || is.null(w$unit)) {
        next
      }
      if (.unitIsDimless(w$unit) && length(.unitBackTransform(terms[[i]], iniNames)) > 0) {
        next # a parameter on its estimation scale, not a unitless value
      }
      if (!.unitConv(w$unit, target)) {
        issues[[c$target]] <- c(
          issues[[c$target]],
          sprintf(
            "'%s' (line %d) has units %s but %s is %s",
            deparse1(terms[[i]]),
            c$line,
            .unitStr(w$unit),
            c$target,
            .unitStr(target)
          )
        )
        next
      }
      k <- .unitFactor(w$unit, target)
      if (abs(k - 1) < 1e-8 || abs(.unitLiteralFactor(terms[[i]]) / k - 1) < 1e-6) {
        next
      }
      conv[[length(conv) + 1L]] <- data.frame(
        line = c$line,
        path = paste(c$path, collapse = ","),
        target = c$target,
        termIndex = i,
        from = .unitStr(w$unit),
        to = .unitStr(target),
        factor = k,
        stringsAsFactors = FALSE
      )
    }
  }
  conversions <- if (length(conv) > 0) {
    do.call(rbind, conv)
  } else {
    data.frame(
      line = integer(),
      path = character(),
      target = character(),
      termIndex = integer(),
      from = character(),
      to = character(),
      factor = numeric(),
      stringsAsFactors = FALSE
    )
  }
  list(issues = issues, conversions = conversions)
}

# ---- inputs ---------------------------------------------------------------------------

# The compartments doses go into when the metadata does not say.
.unitDosedCompartments <- function(ui) {
  meta <- as.list(ui$meta)
  states <- ui$state
  if (is.character(meta$dosing) && length(meta$dosing) > 0) {
    return(intersect(meta$dosing, states))
  }
  props <- ui$props$cmtProp
  if (is.data.frame(props) && nrow(props) > 0) {
    return(intersect(unique(props$Compartment), states))
  }
  utils::head(intersect(c("depot", "central"), states), 1)
}

.unitParses <- function(s) {
  tryCatch(
    {
      .unitOf(s)
      TRUE
    },
    error = function(e) FALSE
  )
}

# Unit strings from the model's `units` metadata. Legacy `dosing` and
# `concentration` keys map onto the dosed compartments and a single endpoint;
# a value that does not parse or a key that names nothing in the model is
# skipped with a note.
.unitDeclaredMeta <- function(legacy, known, dosed, endpoints) {
  decl <- list()
  notes <- character()
  for (nm in names(legacy)) {
    val <- legacy[[nm]]
    if (!is.character(val) || length(val) != 1L || is.na(val)) {
      next
    }
    keys <- switch(
      nm,
      dosing = dosed,
      concentration = if (length(endpoints) == 1L) endpoints else character(),
      if (nm %in% known) nm else character()
    )
    if (length(keys) == 0) {
      notes <- c(notes, sprintf("units$%s ignored: it names nothing in the model", nm))
    } else if (!.unitParses(val)) {
      notes <- c(notes, sprintf("units$%s = '%s' ignored: not a unit", nm, val))
    } else {
      for (k in keys) {
        decl[[k]] <- val
      }
    }
  }
  list(decl = decl, notes = notes)
}

# Declared unit strings, `...` over `units` over `covariateData` over the
# model's `units` metadata. An argument that does not parse or names nothing
# in the model is an error.
.unitDeclared <- function(ui, dots, unitsArg, known, dosed, endpoints) {
  meta <- as.list(ui$meta)
  fromMeta <- .unitDeclaredMeta(meta$units, known, dosed, endpoints)
  decl <- fromMeta$decl
  notes <- fromMeta$notes
  for (nm in names(meta$covariateData)) {
    v <- if (is.list(meta$covariateData[[nm]])) meta$covariateData[[nm]]$units
    if (is.character(v) && length(v) == 1L && !is.na(v) && is.null(decl[[nm]]) && nm %in% known && .unitParses(v)) {
      decl[[nm]] <- v
    }
  }
  if (!is.null(unitsArg)) {
    checkmate::assertList(unitsArg, types = "character", names = "unique")
    decl[names(unitsArg)] <- unitsArg
  }
  if (length(dots) > 0) {
    if (is.null(names(dots)) || any(!nzchar(names(dots))) || anyDuplicated(names(dots)) > 0) {
      stop(.unitError("units passed through `...` must have unique names, e.g. `depot = \"mg\"`"))
    }
    for (nm in names(dots)) {
      checkmate::assertString(dots[[nm]], .var.name = nm)
      decl[[nm]] <- dots[[nm]]
    }
  }
  if (!is.null(decl$linCmt)) {
    decl$rxLinCmt <- decl$linCmt
    decl$linCmt <- NULL
  }
  list(decl = decl, notes = notes)
}

# "mg/L * 1000 = ng/mL", or a division when the reciprocal is a clean integer.
.unitConversionText <- function(from, factor, to) {
  op <- .unitConversionOperator(factor)
  sprintf("%s %s %s = %s", from, op$op, format(op$value, digits = 12, scientific = FALSE), to)
}

.unitConversionOperator <- function(factor) {
  factor <- signif(factor, 12)
  inv <- signif(1 / factor, 12)
  if (factor < 1 && abs(inv - round(inv)) < 1e-8) list(op = "/", value = round(inv)) else list(op = "*", value = factor)
}

# ---- checkUnits() ----------------------------------------------------------------------

#' Check and infer the physical units of an rxode2 / nlmixr2 model
#'
#' `checkUnits()` walks the `model({})` block and works out a unit for every
#' symbol from the units declared at the model's boundary: `time`, the dose
#' unit of each dosed compartment (named by the compartment, `depot = "mg"`),
#' the unit of each output variable (`Cc = "ng/mL"`), and any per-parameter
#' override (`cl = "mL/min"`). Units flow forward through assignments and
#' backward from a line whose target unit is known to the one symbol in it
#' that is not, until nothing more can be learned. A solved volume or
#' clearance is written in the model's default spelling (`L` and `L/<time>`,
#' or `mL/kg` and `mL/kg/<time>` for a dose per body weight); anything else
#' keeps what the equations gave. `addUnits()` applies the result.
#'
#' Units in the model's `units` metadata are the defaults, arguments override
#' them, and a legacy `list(time =, dosing =, concentration =)` block is read
#' as the dose unit of the dosed compartments and the unit of a single
#' endpoint. Covariate units come from `covariateData`; an undeclared covariate
#' is unknown like any other symbol.
#'
#' The rules: `+` and `-` need convertible operands; a bare number added to or
#' subtracted from a quantity is a value in that quantity's unit (`AGE - 40`),
#' one dividing a quantity is a reference value (`WT / 70` is unitless), and
#' any other bare number is unitless. Transcendental functions and comparisons
#' are unitless; the arithmetic inside their arguments is still checked.
#' `exp()` (or `expit()`, `probitInv()`) of nothing but `ini()` parameters is
#' a back-transform: the assigned variable keeps its own unit and the
#' parameters are unitless. A dimensioned quantity raised to a symbolic power
#' (`Cc^hill`) is unknown; raised to a fractional literal power it is an issue.
#' A variable assigned more than once must agree with itself. Residual-error
#' parameters take the unit their distribution implies: `add()` the endpoint's
#' unit, `prop()` and the log, logit and probit families unitless.
#'
#' The `units` package and its udunits-2 system library do the parsing,
#' convertibility and conversion factors; both functions error with an
#' installation hint when it is absent. On first use a few units udunits
#' lacks are installed for the session (`IU`, `U`, `M` for molar, `cells`,
#' `CFU`, `copies`, `molecules`), which is process-global state.
#'
#' @param ui A model: an `rxUi`, a model function, or anything
#'   [rxode2::assertRxUi()] accepts.
#' @param ... Named unit strings, one per symbol: `time`, a dosed compartment,
#'   an endpoint variable (`linCmt` for a `linCmt()` model), an `ini()`
#'   parameter, a covariate, or a variable assigned in `model({})`.
#' @param units The same as `...`, as a named list. `...` takes precedence.
#' @return A `data.frame` with one row per symbol and the columns `name`,
#'   `unit` (as the `units` package spells it, `"unitless"` for a
#'   dimensionless quantity, `NA` when unresolved), `issue` (`NA` when the
#'   symbol is consistent), and `conversion` (the constant the symbol's line
#'   needs, as `"mg/L * 1000 = ng/mL"`, or `NA`). The attribute
#'   `"conversions"` is the data frame [addUnits()] applies and `"notes"` lists
#'   metadata that was ignored.
#' @family units
#' @export
#' @author Bill Denney
#' @examples
#' if (requireNamespace("units", quietly = TRUE)) {
#'   checkUnits(readModelDb("PK_1cmt_des"), time = "h", depot = "mg", Cc = "ng/mL")
#' }
checkUnits <- function(ui, ..., units = NULL) {
  .unitInstall()
  oldSimplify <- units::units_options("simplify")
  units::units_options(simplify = FALSE)
  on.exit(units::units_options(simplify = oldSimplify), add = TRUE)
  .ui <- rxode2::assertRxUi(ui)
  states <- .ui$state %||% character()
  predDf <- .ui$predDf
  endpoints <- if (is.data.frame(predDf)) unique(predDf$var) else character()
  ini <- .ui$iniDf
  if (!is.data.frame(ini)) {
    ini <- data.frame(
      name = character(),
      neta1 = numeric(),
      neta2 = numeric(),
      condition = character(),
      err = character()
    )
  }
  thetas <- ini$name[is.na(ini$neta1) & is.na(ini$err)]
  etas <- ini$name[!is.na(ini$neta1) & ini$neta1 == ini$neta2]
  residuals <- data.frame(name = ini$name, endpoint = ini$condition, err = ini$err, stringsAsFactors = FALSE)[
    !is.na(ini$err),
    ,
    drop = FALSE
  ]
  covs <- .ui$allCovs %||% character()
  cs <- .unitConstraints(.ui$lstExpr)
  isSym <- vapply(cs, function(c) c$kind == "sym", logical(1))
  inter <- setdiff(unique(vapply(cs[isSym], `[[`, character(1), "target")), c(states, endpoints))
  props <- unique(vapply(
    cs[!isSym & vapply(cs, function(c) c$kind != "check", logical(1))],
    `[[`,
    character(1),
    "target"
  ))
  known <- c("time", states, endpoints, ini$name, covs, inter)
  dosed <- .unitDosedCompartments(.ui)
  d <- .unitDeclared(.ui, list(...), units, c(known, "linCmt"), dosed, endpoints)
  bad <- setdiff(names(d$decl), known)
  if (length(bad) > 0) {
    stop(.unitError(sprintf("units given for names that are not in the model: %s", paste(bad, collapse = ", "))))
  }
  u <- list()
  for (nm in names(d$decl)) {
    u[[nm]] <- tryCatch(.unitOf(d$decl[[nm]]), nlmixr2libUnitError = function(e) {
      stop(.unitError(paste0(nm, ": ", conditionMessage(e))))
    })
  }
  reserved <- rxode2::rxReservedKeywords[["Reserved Name"]]
  for (nm in setdiff(reserved[!is.na(reserved)], c("time", "podo", "tlast"))) {
    u[[nm]] <- .unitDimless()
  }
  if (!is.null(u$time)) {
    u$t <- u$time
    u$tlast <- u$time
  }
  dosedDecl <- intersect(names(d$decl), states)
  dose <- if (length(dosedDecl) > 0) u[[dosedDecl[[1]]]]
  if (!is.null(dose)) {
    u$podo <- dose
  }
  u <- .unitInfer(cs, u, ini$name, residuals, .unitDefaults(u$time, dose))
  for (e in setdiff(etas, names(u))) {
    u[[e]] <- .unitDimless()
  }
  v <- .unitValidate(cs, u, ini$name)
  order <- unique(c("time", dosedDecl, states, endpoints, thetas, residuals$name, etas, covs, inter, props))
  order <- c(order, setdiff(names(v$issues), order))
  unitOf <- function(nm) {
    if (nm %in% props) .unitTarget(Filter(function(c) identical(c$target, nm), cs)[[1]], u) else u[[nm]]
  }
  out <- data.frame(
    name = order,
    unit = vapply(order, function(nm) if (is.null(unitOf(nm))) NA_character_ else .unitStr(unitOf(nm)), character(1)),
    issue = vapply(
      order,
      function(nm) if (is.null(v$issues[[nm]])) NA_character_ else paste(v$issues[[nm]], collapse = "; "),
      character(1)
    ),
    conversion = vapply(
      order,
      function(nm) {
        rows <- v$conversions[v$conversions$target == nm, , drop = FALSE]
        if (nrow(rows) == 0) {
          NA_character_
        } else {
          paste(unique(mapply(.unitConversionText, rows$from, rows$factor, rows$to)), collapse = "; ")
        }
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "conversions") <- v$conversions
  attr(out, "notes") <- d$notes
  out
}

# ---- addUnits() ------------------------------------------------------------------------

#' Add units to a model, inserting the conversions its arithmetic needs
#'
#' Runs [checkUnits()] and, when it reports no issue, writes the units into the
#' model: every conversion the arithmetic needs is inserted as a bare constant
#' (`Cc <- central/vc * 1000`) and explained in the `unitConversions`
#' metadata, the `units` metadata lists the unit of every resolved symbol
#' (`"unitless"` for dimensionless ones), `dosing` names the dosed
#' compartments, and a covariate given a unit gets it in `covariateData`.
#' `ini()` estimates are never rescaled: declaring `cl = "mL/min"` for a model
#' written in hours inserts the factor the equations need and leaves the
#' estimate to be read in the new unit.
#'
#' @inheritParams checkUnits
#' @return The model as an `rxUi` with the units metadata and any inserted
#'   conversion constants. An issue reported by [checkUnits()] is an error.
#' @family units
#' @export
#' @author Bill Denney
#' @examples
#' if (requireNamespace("units", quietly = TRUE)) {
#'   res <- addUnits(readModelDb("PK_1cmt_des"), time = "h", depot = "mg", Cc = "ng/mL")
#'   res$meta$units
#'   res$meta$unitConversions
#' }
addUnits <- function(ui, ..., units = NULL) {
  .ui <- rxode2::assertRxUi(ui)
  res <- checkUnits(.ui, ..., units = units)
  issues <- res[!is.na(res$issue), , drop = FALSE]
  if (nrow(issues) > 0) {
    stop(.unitError(paste0(
      "checkUnits() found unit issues; fix them before adding units:\n",
      paste(sprintf("  %s: %s", issues$name, issues$issue), collapse = "\n")
    )))
  }
  conv <- attr(res, "conversions")
  lst <- .ui$lstExpr
  for (i in seq_len(nrow(conv))) {
    path <- as.integer(strsplit(conv$path[i], ",", fixed = TRUE)[[1]])
    lst[[conv$line[i]]] <- .unitInsert(lst[[conv$line[i]]], path, conv$termIndex[i], conv$factor[i])
  }
  .ui <- rxode2::rxUiDecompress(.ui)
  meta <- .ui$meta
  resolved <- res[!is.na(res$unit) & !grepl("[()]", res$name), , drop = FALSE]
  assign("units", stats::setNames(as.list(resolved$unit), resolved$name), envir = meta)
  declared <- c(names(list(...)), names(units))
  dosing <- intersect(unique(c(.unitDosedCompartments(.ui), declared)), .ui$state)
  if (length(dosing) > 0) {
    assign("dosing", dosing, envir = meta)
  }
  if (nrow(conv) > 0) {
    existing <- if (exists("unitConversions", envir = meta)) get("unitConversions", envir = meta) else list()
    for (i in seq_len(nrow(conv))) {
      existing[[conv$target[i]]] <- .unitConversionText(conv$from[i], conv$factor[i], conv$to[i])
    }
    assign("unitConversions", existing, envir = meta)
  }
  covs <- intersect(res$name[!is.na(res$unit)], .ui$allCovs)
  if (length(covs) > 0) {
    covData <- if (exists("covariateData", envir = meta)) get("covariateData", envir = meta) else list()
    for (nm in covs) {
      entry <- covData[[nm]]
      if (!is.list(entry)) {
        entry <- list(description = NA_character_, units = NA_character_, type = "continuous")
      }
      entry$units <- res$unit[res$name == nm]
      covData[[nm]] <- entry
    }
    assign("covariateData", covData, envir = meta)
  }
  rxode2::model(.ui) <- lst
  if (nrow(conv) > 0) {
    lines <- unique(sprintf("%s: %s", conv$target, mapply(.unitConversionText, conv$from, conv$factor, conv$to)))
    cli::cli_inform(c(
      "i" = "Inserted {length(lines)} unit conversion{?s}:",
      stats::setNames(lines, rep("*", length(lines)))
    ))
  }
  rxode2::rxUiCompress(.ui)
}

# Multiply the `termIndex`-th additive term of the assignment at `path` inside
# `e` by `factor`, editing the expression rather than its text.
.unitInsert <- function(e, path, termIndex, factor) {
  if (length(path) > 0) {
    e[[path[[1]]]] <- .unitInsert(e[[path[[1]]]], path[-1], termIndex, factor)
    return(e)
  }
  acc <- new.env(parent = emptyenv())
  acc$n <- 0L
  e[[3]] <- .unitScaleTerm(e[[3]], termIndex, factor, acc)
  e
}

.unitScaleTerm <- function(e, termIndex, factor, acc) {
  inner <- .stripParens(e)
  if (is.call(inner) && is.name(inner[[1]]) && as.character(inner[[1]]) %in% c("+", "-") && length(inner) == 3L) {
    inner[[2]] <- .unitScaleTerm(inner[[2]], termIndex, factor, acc)
    inner[[3]] <- .unitScaleTerm(inner[[3]], termIndex, factor, acc)
    return(inner)
  }
  acc$n <- acc$n + 1L
  if (acc$n != termIndex) {
    return(e)
  }
  op <- .unitConversionOperator(factor)
  call(op$op, e, op$value)
}
