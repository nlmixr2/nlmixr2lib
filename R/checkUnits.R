# Unit tracking for rxode2/nlmixr2 models.
#
# `checkUnits()` walks the model block's expressions, infers a physical unit
# for every symbol from the units declared at the model's boundary (time, the
# dose unit of each dosed compartment, the unit of each output variable, and
# any per-parameter overrides), and reports the result as a data frame.
# `addUnits()` applies that result: it inserts the conversion constants the
# arithmetic needs and writes the units into the model metadata.
#
# The `units` package (and its udunits-2 system library) does the parsing,
# convertibility tests, conversion factors and unit arithmetic. Everything that
# touches it is in the `.unit*` helpers at the top of this file so the rest of
# the code only sees units objects.

# ---- units-package wrapper -------------------------------------------------

.unitsAvailable <- function() {
  requireNamespace("units", quietly = TRUE)
}

.assertUnitsInstalled <- function() {
  if (!.unitsAvailable()) {
    .unitStop(paste0(
      "checkUnits() and addUnits() need the 'units' package, which in turn needs the ",
      "udunits-2 system library. Install it with `install.packages(\"units\")`."
    ))
  }
  invisible(TRUE)
}

.unitStop <- function(message) {
  stop(structure(
    class = c("nlmixr2libUnitError", "error", "condition"),
    list(message = message, call = NULL)
  ))
}

# Units that udunits does not know but pharmacometric models use. A definition
# of "" installs an independent base unit. Installing is process-global and
# permanent for the session; a second bare install of the same symbol errors,
# hence the convertible-to-itself guard. "unitless" is deliberately not
# installed: udunits then names every bare `1` with it, and a string such as
# "h unitless unitless-1" is judged convertible to "unitless h-1". The
# dimensionless unit is written "unitless" everywhere in this package and
# becomes udunits' "1" only in .unitToUdunits().
.unitCustomDefinitions <- c(
  IU = "",
  U = "",
  M = "mol/L",
  cells = "",
  CFU = "",
  copies = "",
  molecules = ""
)

.unitInstallCustom <- function() {
  .assertUnitsInstalled()
  for (sym in names(.unitCustomDefinitions)) {
    if (units::ud_are_convertible(sym, sym)) {
      next
    }
    def <- .unitCustomDefinitions[[sym]]
    if (nzchar(def)) {
      units::install_unit(sym, def)
    } else {
      units::install_unit(sym)
    }
  }
  invisible(TRUE)
}

# The spellings that mean "no unit". Anything here becomes "unitless".
.unitlessSpellings <- c("", "1", "unitless", "dimensionless", "fraction", "none", "n/a", "na")

# Canonical spelling of a unit string before it is parsed. Applies the alias
# table from the conventions, folds the micro sign, and turns the middle dot
# into a multiplication. This is a regex on a unit STRING, never on model code.
.unitNormalizeSpelling <- function(s) {
  checkmate::assertString(s)
  s <- trimws(s)
  # MICRO SIGN, GREEK SMALL LETTER MU, MIDDLE DOT; built at run time so the
  # source stays ASCII
  s <- gsub(intToUtf8(0xB5), "u", s, fixed = TRUE)
  s <- gsub(intToUtf8(0x3BC), "u", s, fixed = TRUE)
  s <- gsub(intToUtf8(0xB7), "*", s, fixed = TRUE)
  # a trailing parenthetical after a space is a descriptor ("umol/L (uM)",
  # "nmol (convert mg via ...)"), not a factor; grouping parentheses follow an
  # operator directly ("1/(nM*h)") and stay
  s <- sub("\\s+\\(.*\\)\\s*$", "", s)
  if (tolower(s) %in% .unitlessSpellings) {
    return("unitless")
  }
  aliases <- .nlmixr2libConventionsStatic$unitSpellingAliases
  for (from in names(aliases)) {
    s <- gsub(paste0("(?<![A-Za-z])", from, "(?![A-Za-z])"), aliases[[from]], s, perl = TRUE)
  }
  s
}

# Widest string that can be a unit. udunits takes tens of seconds over a
# sentence such as "(none; static concentration-response model driven by ...)"
# that the library's legacy metadata carries, so anything longer, or holding a
# separator, is refused before it reaches the parser.
.unitMaxChars <- 40L

.unitLooksLikeUnit <- function(s) {
  nchar(s) <= .unitMaxChars && !grepl("[;,]", s) && !grepl("[[:alpha:]]{3,} [[:alpha:]]{3,}", s)
}

# Parse a unit string into a units object, failing loudly. `units` silently
# reads `kg^0.75` as the scalar 75 and falls back to an opaque symbol with a
# warning when it cannot parse, so both are rejected here.
.unitParse <- function(s) {
  .unitInstallCustom()
  norm <- .unitNormalizeSpelling(s)
  if (!.unitLooksLikeUnit(norm)) {
    .unitStop(sprintf("'%s' is not a unit: it is too long or contains prose.", s))
  }
  if (grepl("\\^\\s*[-+]?[0-9]*\\.[0-9]", norm) || grepl("**", norm, fixed = TRUE)) {
    .unitStop(sprintf(
      "Unit '%s' has a fractional exponent, which the units library cannot represent.",
      s
    ))
  }
  withCallingHandlers(
    tryCatch(
      units::as_units(.unitToUdunits(norm)),
      error = function(e) {
        .unitStop(sprintf("Unit '%s' is not recognized: %s", s, conditionMessage(e)))
      }
    ),
    warning = function(w) {
      .unitStop(sprintf("Unit '%s' could not be parsed: %s", s, conditionMessage(w)))
    }
  )
}

# The last step before a string reaches udunits, which spells the
# dimensionless unit "1".
.unitToUdunits <- function(s) {
  if (identical(s, "unitless")) "1" else s
}

# A units object with value 1 in the given unit.
.unitOne <- function(s) {
  u <- .unitParse(s)
  units::set_units(1, .unitStr(u), mode = "standard")
}

# A unit string or units object as a udunits-parseable string, built from the
# cancelled symbols so nothing but real unit tokens reaches udunits.
.unitStr <- function(x) {
  if (!inherits(x, "units")) {
    return(.unitToUdunits(.unitNormalizeSpelling(x)))
  }
  s <- .unitSymbols(x)
  if (length(s$num) == 0 && length(s$den) == 0) {
    return("1")
  }
  paste(c(s$num, if (length(s$den) > 0) paste0(s$den, "-1")), collapse = " ")
}

# Time units, kept last in a denominator so a clearance reads "mL/kg/h".
.unitTimeSymbols <- c("s", "min", "h", "d", "day", "week", "month", "year", "years", "yr")

# Numerator and denominator symbols with identical tokens cancelled
# (`L*h/h` is `L`), without the "unitless" token udunits reports for a bare 1,
# and ordered with time last. Only identical tokens cancel: `mg/kg` stays
# `mg/kg`, which is why checkUnits() runs with units simplification off.
.unitSymbols <- function(u) {
  s <- units(u)
  num <- s$numerator[!s$numerator %in% c("unitless", "1")]
  den <- s$denominator[!s$denominator %in% c("unitless", "1")]
  for (tok in unique(num)) {
    n <- min(sum(num == tok), sum(den == tok))
    if (n > 0) {
      num <- num[-which(num == tok)[seq_len(n)]]
      den <- den[-which(den == tok)[seq_len(n)]]
    }
  }
  timeLast <- function(x) {
    x <- sort(x)
    c(x[!x %in% .unitTimeSymbols], x[x %in% .unitTimeSymbols])
  }
  list(num = timeLast(num), den = timeLast(den))
}

.unitIdentical <- function(a, b) {
  identical(.unitSymbols(a), .unitSymbols(b))
}

.unitConvertible <- function(a, b) {
  sa <- .unitStr(a)
  sb <- .unitStr(b)
  if (isTRUE(units::ud_are_convertible(sa, sb))) {
    return(TRUE)
  }
  # ud_are_convertible() false-negatives on units with a numeric literal
  # (mL/min/1.73m^2); set_units() gets those right.
  ok <- tryCatch(
    {
      units::set_units(units::set_units(1, sa, mode = "standard"), sb, mode = "standard")
      TRUE
    },
    error = function(e) FALSE
  )
  isTRUE(ok)
}

# Numeric factor that turns 1 `a` into `b`.
.unitFactor <- function(a, b) {
  units::ud_convert(1, .unitStr(a), .unitStr(b))
}

# The value of quantity `q` (a units object whose value carries the scale
# accumulated by unit arithmetic) expressed in `target`.
.unitScaleTo <- function(q, target) {
  as.numeric(units::set_units(q, .unitStr(target), mode = "standard"))
}

.unitIsDimensionless <- function(u) {
  .unitConvertible(u, "unitless")
}

.unitMul <- function(a, b) {
  a * b
}

.unitDiv <- function(a, b) {
  a / b
}

.unitPow <- function(u, k) {
  checkmate::assertIntegerish(k, len = 1)
  u^as.integer(k)
}

# Canonical "mg/L", "L/h", "1/h", "mL/kg/h", "unitless" spelling of a units
# object for metadata. `units::deparse_unit()` gives "mg L-1", which is not
# what a model file should read.
.unitDeparse <- function(u) {
  s <- .unitSymbols(u)
  if (length(s$num) == 0 && length(s$den) == 0) {
    return("unitless")
  }
  collapse <- function(x) {
    if (length(x) == 0) {
      return(character())
    }
    vapply(
      unique(x),
      function(nm) {
        n <- sum(x == nm)
        if (n == 1) nm else sprintf("%s^%d", nm, n)
      },
      character(1)
    )
  }
  num <- collapse(s$num)
  den <- collapse(s$den)
  if (length(num) == 0) {
    if (length(den) == 1) {
      return(paste0("1/", den))
    }
    return(sprintf("1/(%s)", paste(den, collapse = "*")))
  }
  out <- paste(num, collapse = "*")
  if (length(den) > 0) {
    out <- paste(c(out, den), collapse = "/")
  }
  out
}

# ---- symbol table -------------------------------------------------------------

# Functions whose result has no unit. The arithmetic inside their argument is
# still validated, and a dimensioned argument is an issue.
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
  "trigamma",
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

# Back-transforms: applied to nothing but ini() parameters and covariate terms
# they name a parameter on its estimation scale, so the assigned variable keeps
# its own unit rather than becoming unitless.
.unitBackTransform <- c("exp", "expit", "probitInv", "expm1")

# Functions whose result carries the unit of their (first) argument.
.unitPassThrough <- c("abs", "floor", "ceiling", "round", "trunc", "signif", "max", "min", "ifelse")

.unitComparison <- c("<", ">", "<=", ">=", "==", "!=", "&&", "||", "&", "|", "!")

# The residual-error structures rxode2 writes into `iniDf$err`, mapped to what
# the parameter represents in its distribution: "endpoint" means the endpoint
# variable's unit, "unitless" a fraction or a transformed-scale quantity, and
# NA a unit that cannot be written (a power-law SD carries var^(1-c)).
.unitResidualErrorMap <- c(
  add = "endpoint",
  dnorm = "endpoint",
  prop = "unitless",
  propT = "unitless",
  propF = "unitless",
  pow = NA_character_,
  powT = NA_character_,
  powF = NA_character_,
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

# Result of walking one expression. `q` is NULL when the unit is unknown; `any`
# marks a covariate without a declared unit (compatible with anything);
# `factor` is the product of bare numeric literals in multiplicative position;
# `bt` marks a back-transform of ini() parameters; `btSyms` are those
# parameters.
.uq <- function(q = NULL, factor = 1, any = FALSE, bt = FALSE, btSyms = character()) {
  list(q = q, factor = factor, any = any, bt = bt, btSyms = btSyms)
}

.uqKnown <- function(r) {
  !is.null(r$q) && !r$any
}

.uqUnknown <- function(r) {
  is.null(r$q) && !r$any
}

# The state of one checkUnits() run: an environment so the recursive walkers
# can record without superassignment.
.unitState <- function() {
  st <- new.env(parent = emptyenv())
  st$u <- list() # symbol -> units object (value 1)
  st$source <- character() # symbol -> how the unit was determined
  st$type <- character() # symbol -> symbol class
  st$transformOf <- character() # symbol -> back-transformed variable
  st$anySyms <- character() # covariates without a declared unit
  st$defLine <- list() # symbol -> integer line indices
  st$issues <- list() # name -> character messages
  st$conversions <- list() # list of data.frame rows
  st$notes <- character() # human-readable notes for the caller
  st$progress <- FALSE # set by .unitSet(); read by the inference loop
  st$solving <- FALSE # TRUE while inferring, so a sum may solve its unknown side
  st
}

.unitSet <- function(st, sym, u, source, type = NULL) {
  st$u[[sym]] <- u
  st$source[sym] <- source
  if (!is.null(type)) {
    st$type[sym] <- type
  }
  st$progress <- TRUE
  invisible(NULL)
}

.unitGet <- function(st, sym) {
  st$u[[sym]]
}

.unitKnown <- function(st, sym) {
  !is.null(st$u[[sym]])
}

.unitIssue <- function(st, name, message) {
  st$issues[[name]] <- unique(c(st$issues[[name]], message))
  invisible(NULL)
}

# A literal number inside an expression.
.isNumericLiteral <- function(e) {
  is.numeric(e) && length(e) == 1L && !is.na(e)
}

# The symbols an expression reads, excluding function names.
.unitSymbolsIn <- function(e) {
  .modelLineSymbols(list(e))
}

# TRUE when every symbol in `e` is an ini() parameter, a covariate, or a
# reserved constant: the shape of a back-transform argument.
.backTransformEligible <- function(e, st) {
  syms <- .unitSymbolsIn(e)
  if (length(syms) == 0) {
    return(FALSE)
  }
  ok <- st$type[syms]
  all(!is.na(ok) & ok %in% c("parameter", "eta", "covariate", "reserved"))
}

# TRUE when every symbol in `e` is a covariate or a reserved constant: the
# shape of a covariate normalisation such as `WT / 70`.
.covariateOnly <- function(e, st) {
  syms <- .unitSymbolsIn(e)
  if (length(syms) == 0) {
    return(FALSE)
  }
  ok <- st$type[syms]
  all(!is.na(ok) & ok %in% c("covariate", "reserved"))
}

# ---- expression walker ----------------------------------------------------------

# Unit of one expression. `where` names the line for issue messages.
.unitOfExpr <- function(e, st, where) {
  if (.isNumericLiteral(e)) {
    return(.uq(q = st$unitless, factor = as.numeric(e)))
  }
  if (is.name(e)) {
    sym <- as.character(e)
    if (sym %in% st$anySyms) {
      return(.uq(any = TRUE))
    }
    if (.unitKnown(st, sym)) {
      return(.uq(q = .unitGet(st, sym)))
    }
    return(.uq())
  }
  if (!is.call(e)) {
    return(.uq())
  }
  if (!is.name(e[[1]])) {
    return(.uq())
  }
  fn <- as.character(e[[1]])
  args <- as.list(e)[-1]
  if (fn == "(") {
    return(.unitOfExpr(args[[1]], st, where))
  }
  if (fn %in% c("-", "+") && length(args) == 1) {
    return(.unitOfExpr(args[[1]], st, where))
  }
  if (fn %in% c("*", "/")) {
    return(.unitOfProduct(fn, args[[1]], args[[2]], st, where))
  }
  if (fn %in% c("+", "-")) {
    return(.unitOfSum(fn, args[[1]], args[[2]], st, where))
  }
  if (fn == "^") {
    return(.unitOfPower(args[[1]], args[[2]], st, where))
  }
  if (fn %in% .unitComparison) {
    return(.unitOfComparison(fn, args, st, where))
  }
  if (fn %in% .unitTranscendental) {
    return(.unitOfTranscendental(fn, args, st, where))
  }
  if (fn %in% .unitPassThrough) {
    return(.unitOfPassThrough(fn, args, st, where))
  }
  if (fn %in% c("transit", "linCmt")) {
    if (fn == "transit" && length(args) >= 2) {
      mtt <- .unitOfExpr(args[[2]], st, where)
      if (.uqUnknown(mtt) && is.name(args[[2]]) && !is.null(st$time)) {
        .unitSet(st, as.character(args[[2]]), st$time, "inferred")
      } else if (.uqKnown(mtt) && !is.null(st$time) && !.unitConvertible(mtt$q, st$time)) {
        .unitIssue(
          st,
          where,
          sprintf(
            "transit() mean transit time has units %s, not a time.",
            .unitDeparse(mtt$q)
          )
        )
      }
    }
    return(.uq(any = TRUE))
  }
  # An unknown function: validate its arguments, result unknown.
  for (a in args) {
    .unitOfExpr(a, st, where)
  }
  .uq()
}

.unitOfProduct <- function(fn, lhs, rhs, st, where) {
  a <- .unitOfExpr(lhs, st, where)
  b <- .unitOfExpr(rhs, st, where)
  bt <- a$bt || b$bt
  btSyms <- c(a$btSyms, b$btSyms)
  # A literal divisor of a covariate term is that covariate's reference value
  # (`WT / 70`), so the ratio has no unit.
  if (fn == "/" && .isNumericLiteral(rhs) && !.isNumericLiteral(lhs) && .covariateOnly(lhs, st)) {
    if (a$any) {
      return(.uq(any = TRUE))
    }
    if (.uqKnown(a)) {
      return(.uq(q = st$unitless, factor = a$factor))
    }
    return(.uq())
  }
  if (a$any) {
    a <- .uq(q = st$unitless)
  }
  if (b$any) {
    b <- .uq(q = st$unitless)
  }
  if (!.uqKnown(a) || !.uqKnown(b)) {
    return(.uq(bt = bt, btSyms = btSyms))
  }
  if (fn == "*") {
    return(.uq(q = .unitMul(a$q, b$q), factor = a$factor * b$factor, bt = bt, btSyms = btSyms))
  }
  .uq(q = .unitDiv(a$q, b$q), factor = a$factor / b$factor, bt = bt, btSyms = btSyms)
}

.unitOfSum <- function(fn, lhs, rhs, st, where) {
  a <- .unitOfExpr(lhs, st, where)
  b <- .unitOfExpr(rhs, st, where)
  # A literal added to or subtracted from a known quantity is a value in that
  # quantity's unit (`AGE - 40`, `t - 2`); added to an unknown one it is a bare
  # number, so the unknown side must be unitless (`1 + e_age_cl * (AGE - 40)`).
  if (.isNumericLiteral(lhs) || .isNumericLiteral(rhs)) {
    other <- if (.isNumericLiteral(lhs)) rhs else lhs
    r <- if (.isNumericLiteral(lhs)) b else a
    if (.uqUnknown(r) && st$solving && .unitSolveTerm(other, st$unitless, st, where)) {
      r <- .unitOfExpr(other, st, where)
    }
    return(.uq(q = r$q, any = r$any))
  }
  if (a$any && b$any) {
    return(.uq(any = TRUE))
  }
  if (a$any) {
    return(.uq(q = b$q))
  }
  if (b$any) {
    return(.uq(q = a$q))
  }
  # A known side names the unit of an unknown side (`ec50 + Cc`).
  if (st$solving && .uqKnown(a) && .uqUnknown(b)) {
    .unitSolveTerm(rhs, a$q, st, where)
    return(.uq(q = a$q))
  }
  if (st$solving && .uqKnown(b) && .uqUnknown(a)) {
    .unitSolveTerm(lhs, b$q, st, where)
    return(.uq(q = b$q))
  }
  if (.uqKnown(a) && .uqKnown(b)) {
    if (!.unitConvertible(a$q, b$q)) {
      .unitIssue(
        st,
        where,
        sprintf(
          "'%s' adds %s to %s, which are not convertible.",
          deparse1(call(fn, lhs, rhs)),
          .unitDeparse(a$q),
          .unitDeparse(b$q)
        )
      )
    } else if (!.unitIdentical(a$q, b$q) || abs(.unitScaleTo(a$q, b$q) - 1) > 1e-8) {
      .unitIssue(
        st,
        where,
        sprintf(
          "'%s' mixes %s and %s in a sum; write both terms in the same unit.",
          deparse1(call(fn, lhs, rhs)),
          .unitDeparse(a$q),
          .unitDeparse(b$q)
        )
      )
    }
    return(.uq(q = a$q))
  }
  if (.uqKnown(a)) {
    return(.uq(q = a$q))
  }
  if (.uqKnown(b)) {
    return(.uq(q = b$q))
  }
  .uq()
}

.unitOfPower <- function(base, expo, st, where) {
  a <- .unitOfExpr(base, st, where)
  if (.isNumericLiteral(expo)) {
    k <- as.numeric(expo)
    if (.isNumericLiteral(base)) {
      # `10^x` is a literal power: a bare number.
      return(.uq(q = st$unitless, factor = as.numeric(base)^k))
    }
    if (a$any) {
      return(.uq(any = TRUE))
    }
    if (!.uqKnown(a)) {
      return(.uq())
    }
    if (abs(k - round(k)) < 1e-12) {
      return(.uq(q = .unitPow(a$q, round(k)), factor = a$factor^k))
    }
    if (.unitIsDimensionless(a$q)) {
      return(.uq(q = st$unitless, factor = a$factor^k))
    }
    .unitIssue(
      st,
      where,
      sprintf(
        "'%s' raises a quantity in %s to the power %s; normalise it by a reference value first.",
        deparse1(base),
        .unitDeparse(a$q),
        format(k)
      )
    )
    return(.uq(q = st$unitless))
  }
  # A symbolic exponent must be unitless; `10^ltheta` is a back-transform.
  b <- .unitOfExpr(expo, st, where)
  if (.uqKnown(b) && !.unitIsDimensionless(b$q)) {
    .unitIssue(
      st,
      where,
      sprintf(
        "the exponent '%s' has units %s; exponents must be unitless.",
        deparse1(expo),
        .unitDeparse(b$q)
      )
    )
  }
  if (.isNumericLiteral(base) && .backTransformEligible(expo, st)) {
    return(.uq(q = st$unitless, bt = TRUE, btSyms = .unitSymbolsIn(expo)))
  }
  # a dimensioned base raised to a symbolic power (`Cc^hill`) has no writable
  # unit, so it is unknown rather than an issue
  if (a$any) {
    return(.uq(any = TRUE))
  }
  if (.uqKnown(a) && .unitIsDimensionless(a$q)) {
    return(.uq(q = st$unitless))
  }
  .uq()
}

# A comparison is unitless. Its two sides follow the sum rules: a literal is
# a value in the other side's unit (`t < 8`), a known side names the unit of
# an unknown side (`PNA <= pna_threshold`), and two known sides must be
# convertible.
.unitOfComparison <- function(fn, args, st, where) {
  rs <- lapply(args, .unitOfExpr, st = st, where = where)
  if (length(rs) == 2 && !.isNumericLiteral(args[[1]]) && !.isNumericLiteral(args[[2]])) {
    a <- rs[[1]]
    b <- rs[[2]]
    if (.uqKnown(a) && .uqKnown(b) && !.unitConvertible(a$q, b$q)) {
      .unitIssue(
        st,
        where,
        sprintf(
          "'%s' compares %s with %s, which are not convertible.",
          deparse1(as.call(c(list(as.name(fn)), args))),
          .unitDeparse(a$q),
          .unitDeparse(b$q)
        )
      )
    } else if (st$solving && .uqKnown(a) && .uqUnknown(b)) {
      .unitSolveTerm(args[[2]], a$q, st, where)
    } else if (st$solving && .uqKnown(b) && .uqUnknown(a)) {
      .unitSolveTerm(args[[1]], b$q, st, where)
    }
  }
  .uq(q = st$unitless)
}

.unitOfTranscendental <- function(fn, args, st, where) {
  if (length(args) == 0) {
    return(.uq(q = st$unitless))
  }
  arg <- args[[1]]
  # the arithmetic inside is validated; the argument's own unit is not judged
  # (`log(CFU/mL)` is an ordinary endpoint)
  .unitOfExpr(arg, st, where)
  for (extra in args[-1]) {
    .unitOfExpr(extra, st, where)
  }
  if (fn %in% .unitBackTransform && .backTransformEligible(arg, st)) {
    return(.uq(q = st$unitless, bt = TRUE, btSyms = .unitSymbolsIn(arg)))
  }
  .uq(q = st$unitless)
}

.unitOfPassThrough <- function(fn, args, st, where) {
  if (fn == "ifelse") {
    .unitOfExpr(args[[1]], st, where)
    args <- args[-1]
  }
  rs <- lapply(args, .unitOfExpr, st = st, where = where)
  if (length(rs) == 1) {
    return(.uq(q = rs[[1]]$q, factor = rs[[1]]$factor, any = rs[[1]]$any))
  }
  # The operands follow the sum rules: a literal is a value in the unit of a
  # known operand (`max(0, Cc - thr)`); with no known operand a non-zero
  # literal is a bare number and the unknown operands must be unitless
  # (`ifelse(t < 8, 1, tfrac)`); zero is the same in any unit.
  lit <- vapply(args, .isNumericLiteral, logical(1))
  knownIdx <- which(!lit & vapply(rs, .uqKnown, logical(1)))
  unknownIdx <- which(!lit & vapply(rs, .uqUnknown, logical(1)))
  if (length(knownIdx) >= 1) {
    ref <- rs[[knownIdx[[1]]]]$q
    for (i in knownIdx[-1]) {
      if (!.unitConvertible(ref, rs[[i]]$q)) {
        .unitIssue(st, where, sprintf("%s() mixes %s and %s.", fn, .unitDeparse(ref), .unitDeparse(rs[[i]]$q)))
      }
    }
    if (st$solving) {
      for (i in unknownIdx) {
        .unitSolveTerm(args[[i]], ref, st, where)
      }
    }
    return(.uq(q = ref))
  }
  nonZeroLit <- lit & vapply(args, function(a) .isNumericLiteral(a) && a != 0, logical(1))
  if (any(nonZeroLit)) {
    if (st$solving) {
      for (i in unknownIdx) {
        .unitSolveTerm(args[[i]], st$unitless, st, where)
      }
    }
    return(.uq(q = st$unitless))
  }
  if (any(vapply(rs, function(r) r$any, logical(1)))) {
    return(.uq(any = TRUE))
  }
  .uq()
}

# ---- constraints -----------------------------------------------------------------

# One constraint per assignment in the model block: `target` is the assigned
# symbol or property (`d/dt(central)`), `kind` says how the target's unit
# relates to a symbol, `rhs` is the expression, `path` indexes the assignment
# call inside `lstExpr[[line]]` (so an `if` branch can be edited in place).
.unitConstraints <- function(lstExpr) {
  acc <- new.env(parent = emptyenv())
  acc$rows <- list()
  for (i in seq_along(lstExpr)) {
    .unitConstraintsFrom(lstExpr[[i]], i, integer(), acc)
  }
  acc$rows
}

.unitConstraintsFrom <- function(e, line, path, acc) {
  if (!is.call(e)) {
    return(invisible(NULL))
  }
  head <- e[[1]]
  if (!is.name(head)) {
    return(invisible(NULL))
  }
  fn <- as.character(head)
  if (fn == "{") {
    for (j in seq_along(e)[-1]) {
      .unitConstraintsFrom(e[[j]], line, c(path, j), acc)
    }
    return(invisible(NULL))
  }
  if (fn == "if") {
    acc$rows[[length(acc$rows) + 1]] <- list(
      line = line,
      path = c(path, 2L),
      target = NA_character_,
      kind = "check",
      rhs = e[[2]]
    )
    .unitConstraintsFrom(e[[3]], line, c(path, 3L), acc)
    if (length(e) >= 4) {
      .unitConstraintsFrom(e[[4]], line, c(path, 4L), acc)
    }
    return(invisible(NULL))
  }
  if (!fn %in% c("<-", "=")) {
    return(invisible(NULL))
  }
  lhs <- e[[2]]
  rhs <- e[[3]]
  if (is.name(lhs)) {
    acc$rows[[length(acc$rows) + 1]] <- list(
      line = line,
      path = path,
      target = as.character(lhs),
      kind = "assign",
      rhs = rhs
    )
    return(invisible(NULL))
  }
  if (is.call(lhs) && is.name(lhs[[1]])) {
    lfn <- as.character(lhs[[1]])
    if (lfn == "/" && length(lhs) == 3 && is.call(lhs[[3]]) && identical(as.character(lhs[[3]][[1]]), "dt")) {
      cmt <- as.character(lhs[[3]][[2]])
      acc$rows[[length(acc$rows) + 1]] <- list(
        line = line,
        path = path,
        target = sprintf("d/dt(%s)", cmt),
        kind = "ddt",
        rhs = rhs,
        cmt = cmt
      )
      return(invisible(NULL))
    }
    if (lfn %in% c("f", "F", "rate", "dur", "alag", "lag") && length(lhs) == 2 && is.name(lhs[[2]])) {
      cmt <- as.character(lhs[[2]])
      kind <- switch(lfn, f = "f", F = "f", rate = "rate", dur = "time", alag = "time", lag = "time")
      acc$rows[[length(acc$rows) + 1]] <- list(
        line = line,
        path = path,
        target = sprintf("%s(%s)", lfn, cmt),
        kind = kind,
        rhs = rhs,
        cmt = cmt
      )
    }
  }
  invisible(NULL)
}

# Top-level additive terms of an expression, each with its sign.
.unitSplitTerms <- function(e) {
  e <- .stripParens(e)
  if (is.call(e) && is.name(e[[1]])) {
    fn <- as.character(e[[1]])
    if (fn %in% c("+", "-") && length(e) == 3) {
      lhs <- .unitSplitTerms(e[[2]])
      rhs <- .unitSplitTerms(e[[3]])
      if (fn == "-") {
        rhs <- lapply(rhs, function(t) call("-", t))
      }
      return(c(lhs, rhs))
    }
  }
  list(e)
}

# Multiplicative leaves of a term with their integer powers.
.unitTermLeaves <- function(e, power = 1L) {
  e <- .stripParens(e)
  if (is.call(e) && is.name(e[[1]])) {
    fn <- as.character(e[[1]])
    if (fn %in% c("-", "+") && length(e) == 2) {
      return(.unitTermLeaves(e[[2]], power))
    }
    if (fn == "*" && length(e) == 3) {
      return(c(.unitTermLeaves(e[[2]], power), .unitTermLeaves(e[[3]], power)))
    }
    if (fn == "/" && length(e) == 3) {
      return(c(.unitTermLeaves(e[[2]], power), .unitTermLeaves(e[[3]], -power)))
    }
    if (fn == "^" && length(e) == 3 && .isNumericLiteral(e[[3]])) {
      k <- as.numeric(e[[3]])
      if (abs(k - round(k)) < 1e-12) {
        return(.unitTermLeaves(e[[2]], power * as.integer(round(k))))
      }
    }
  }
  list(list(expr = e, power = power))
}

# The unit a constraint's target must have, or NULL when not yet known.
.unitTargetUnit <- function(cs, st) {
  switch(
    cs$kind,
    check = st$unitless,
    assign = .unitGet(st, cs$target),
    ddt = ,
    rate = {
      su <- .unitGet(st, cs$cmt)
      if (is.null(su) || is.null(st$time)) NULL else .unitDiv(su, st$time)
    },
    f = st$unitless,
    time = st$time
  )
}

# Give the target of a constraint a unit found by forward inference.
.unitSetTarget <- function(cs, q, st) {
  if (cs$kind == "assign") {
    .unitSet(st, cs$target, .unitRepresentative(q, st), "inferred")
    return(TRUE)
  }
  if (cs$kind == "ddt" && !is.null(st$time)) {
    .unitSet(st, cs$cmt, .unitRepresentative(.unitMul(q, st$time), st), "inferred", "state")
    return(TRUE)
  }
  FALSE
}

# Choose the unit written for an inferred quantity. A volume or a clearance
# is written in the model's default spelling (L and L/h, or mL/kg and
# mL/kg/h for a dose per body weight). Any other inferred unit is kept as the
# equations gave it (`1/year` for a covariate slope, `h` for a half-life)
# unless it mixes prefixes of one dimension (`mg*mL/ng` from an mg dose and
# an ng/mL concentration), in which case the model's representative unit for
# that dimension is used and the arithmetic then needs a conversion.
.unitRepresentative <- function(q, st) {
  s <- .unitSymbols(q)
  if (length(s$num) == 0 && length(s$den) == 0) {
    return(st$unitless)
  }
  for (cand in st$defaultCandidates) {
    if (.unitConvertible(q, cand)) {
      return(cand)
    }
  }
  if (!.unitMixesPrefixes(s)) {
    return(units::set_units(1, .unitStr(q), mode = "standard"))
  }
  for (cand in st$candidates) {
    if (.unitConvertible(q, cand)) {
      return(cand)
    }
  }
  units::set_units(1, .unitStr(q), mode = "standard")
}

# TRUE when a numerator token and a denominator token are convertible to each
# other without being identical (mg over ng): the quotient is a bare number
# that the model's arithmetic has to carry.
.unitMixesPrefixes <- function(s) {
  for (a in unique(s$num)) {
    for (b in unique(s$den)) {
      if (a != b && units::ud_are_convertible(a, b)) {
        return(TRUE)
      }
    }
  }
  FALSE
}

# Representative units by dimension, built from the model's time unit, the
# dose unit of the first dosed compartment, and the first endpoint's unit.
.unitCandidates <- function(st, doseUnit, endpointUnit) {
  defaults <- .nlmixr2libConventionsStatic$unitDefaults
  weightBased <- !is.null(doseUnit) && "kg" %in% .unitSymbols(doseUnit)$den
  scheme <- if (weightBased) defaults$weight else defaults$mass
  out <- list()
  if (!is.null(st$time)) {
    tm <- .unitDeparse(st$time)
    st$defaultCandidates <- list(
      .unitOne(sub("<time>", tm, scheme[["clearance"]], fixed = TRUE)),
      .unitOne(scheme[["volume"]])
    )
    out <- c(out, list(.unitOne(paste0("1/", tm)), st$time))
  } else {
    st$defaultCandidates <- list(.unitOne(scheme[["volume"]]))
  }
  if (!is.null(doseUnit)) {
    out <- c(out, list(doseUnit))
    if (!is.null(st$time)) {
      out <- c(out, list(.unitDiv(doseUnit, st$time)))
    }
  }
  if (!is.null(endpointUnit)) {
    out <- c(out, list(endpointUnit))
    if (!is.null(st$time)) {
      out <- c(out, list(.unitMul(endpointUnit, st$time)))
    }
  }
  out
}

# Solve a constraint's single unknown symbol from the target unit.
.unitSolveTerm <- function(term, target, st, where) {
  leaves <- .unitTermLeaves(term)
  unknown <- list()
  rest <- st$unitless
  for (lf in leaves) {
    r <- .unitOfExpr(lf$expr, st, where)
    if (.uqUnknown(r)) {
      if (!is.name(lf$expr)) {
        return(FALSE)
      }
      unknown[[length(unknown) + 1]] <- lf
      next
    }
    q <- if (r$any) st$unitless else r$q
    rest <- .unitMul(rest, .unitPow(q, lf$power))
  }
  if (length(unknown) != 1L || abs(unknown[[1]]$power) != 1L) {
    return(FALSE)
  }
  sym <- as.character(unknown[[1]]$expr)
  q <- .unitDiv(target, rest)
  if (unknown[[1]]$power == -1L) {
    q <- .unitDiv(st$unitless, q)
  }
  type <- st$type[sym]
  if (!is.na(type) && type %in% c("parameter", "eta", "state")) {
    .unitSet(st, sym, .unitRepresentative(q, st), "inferred")
  } else {
    .unitSet(st, sym, .unitRepresentative(q, st), "inferred", "intermediate")
  }
  TRUE
}

# Residual-error parameters take their unit from the endpoint they describe.
.unitResidualPass <- function(st) {
  for (i in seq_len(nrow(st$residuals))) {
    row <- st$residuals[i, ]
    if (.unitKnown(st, row$name)) {
      next
    }
    kind <- .unitResidualErrorMap[row$err]
    if (is.na(row$err) || !row$err %in% names(.unitResidualErrorMap) || is.na(kind)) {
      next
    }
    if (kind == "unitless") {
      .unitSet(st, row$name, st$unitless, "boundary")
    } else if (.unitKnown(st, row$endpoint)) {
      .unitSet(st, row$name, .unitGet(st, row$endpoint), "boundary")
    }
  }
  invisible(NULL)
}

# Record the parameters inside a back-transform as unitless.
.unitRecordBackTransform <- function(r, cs, st) {
  if (!r$bt || cs$kind == "check") {
    return(invisible(NULL))
  }
  for (sym in unique(r$btSyms)) {
    if (!.unitKnown(st, sym) && st$type[sym] %in% c("parameter", "eta")) {
      .unitSet(st, sym, st$unitless, "inferred")
      st$transformOf[sym] <- cs$target
    }
  }
  invisible(NULL)
}

.unitWhere <- function(cs) {
  if (is.na(cs$target)) {
    return(sprintf("line %d", cs$line))
  }
  cs$target
}

# Forward and backward inference to a fixed point.
.unitInfer <- function(constraints, st) {
  st$solving <- TRUE
  on.exit(st$solving <- FALSE, add = TRUE)
  repeat {
    st$progress <- FALSE
    .unitResidualPass(st)
    for (cs in constraints) {
      where <- .unitWhere(cs)
      if (cs$kind == "check") {
        # an `if` condition solves nothing by itself, but a comparison inside
        # it names the unit of an unknown side
        .unitOfExpr(cs$rhs, st, where)
        next
      }
      target <- .unitTargetUnit(cs, st)
      terms <- .unitSplitTerms(cs$rhs)
      for (term in terms) {
        r <- .unitOfExpr(term, st, where)
        .unitRecordBackTransform(r, cs, st)
        if (!is.null(target)) {
          if (.uqUnknown(r)) {
            .unitSolveTerm(term, target, st, where)
          }
        } else if (.uqKnown(r) && !r$bt) {
          if (.unitSetTarget(cs, r$q, st)) {
            target <- .unitTargetUnit(cs, st)
          }
        }
      }
    }
    if (!st$progress) {
      break
    }
  }
  invisible(NULL)
}

# Validate every constraint against the final units and record the
# conversions the arithmetic needs.
.unitValidate <- function(constraints, st) {
  st$issues <- list()
  for (cs in constraints) {
    where <- .unitWhere(cs)
    if (cs$kind == "check") {
      .unitOfExpr(cs$rhs, st, where)
      next
    }
    target <- .unitTargetUnit(cs, st)
    terms <- .unitSplitTerms(cs$rhs)
    for (i in seq_along(terms)) {
      r <- .unitOfExpr(terms[[i]], st, where)
      if (is.null(target) || !.uqKnown(r) || r$bt) {
        next
      }
      if (!.unitConvertible(r$q, target)) {
        .unitIssue(
          st,
          where,
          sprintf(
            "'%s' (line %d) has units %s but %s is %s.",
            deparse1(terms[[i]]),
            cs$line,
            .unitDeparse(r$q),
            cs$target,
            .unitDeparse(target)
          )
        )
        next
      }
      k <- .unitScaleTo(r$q, target)
      if (.unitIdentical(r$q, target) && abs(k - 1) < 1e-8) {
        next
      }
      if (abs(r$factor / k - 1) < 1e-6) {
        next
      }
      st$conversions[[length(st$conversions) + 1]] <- data.frame(
        line = cs$line,
        path = paste(cs$path, collapse = ","),
        target = cs$target,
        termIndex = i,
        from = .unitDeparse(r$q),
        to = .unitDeparse(target),
        factor = k,
        stringsAsFactors = FALSE
      )
    }
  }
  invisible(NULL)
}

# ---- input resolution ----------------------------------------------------------

# TRUE when a unit string parses.
.unitParses <- function(s) {
  tryCatch(
    {
      .unitParse(s)
      TRUE
    },
    error = function(e) FALSE
  )
}

# Units from the model's own `units` metadata. Legacy `dosing` and
# `concentration` keys are mapped onto the dosed compartments and the single
# endpoint; values that do not parse, and keys that name nothing in the model,
# are dropped with a note (the metadata is free-form; arguments are not).
.unitDeclarationsFromMeta <- function(legacy, st, known) {
  declared <- list()
  if (!is.list(legacy)) {
    return(declared)
  }
  for (nm in names(legacy)) {
    val <- legacy[[nm]]
    if (!is.character(val) || length(val) != 1 || is.na(val)) {
      next
    }
    if (!nm %in% c("dosing", "concentration", known)) {
      st$notes <- c(st$notes, sprintf("metadata units$%s names nothing in the model and was ignored", nm))
      next
    }
    if (!.unitParses(val)) {
      st$notes <- c(st$notes, sprintf("metadata units$%s = '%s' is not a parseable unit and was ignored", nm, val))
      next
    }
    if (nm == "dosing") {
      for (cmt in st$dosedDefault) {
        declared[[cmt]] <- val
      }
    } else if (nm == "concentration") {
      if (length(st$endpoints) == 1) {
        declared[[st$endpoints]] <- val
      } else {
        st$notes <- c(
          st$notes,
          "metadata units$concentration was ignored because the model has more than one endpoint"
        )
      }
    } else {
      declared[[nm]] <- val
    }
  }
  declared
}

# Covariate units from `covariateData`; "(binary)"-style descriptors are not
# units.
.unitDeclarationsFromCovariateData <- function(cov) {
  declared <- list()
  if (!is.list(cov)) {
    return(declared)
  }
  for (nm in names(cov)) {
    entry <- cov[[nm]]
    u <- if (is.list(entry)) entry$units else NULL
    if (is.character(u) && length(u) == 1 && !is.na(u) && !grepl("^\\(", u) && .unitParses(u)) {
      declared[[nm]] <- u
    }
  }
  declared
}

# Units from the `units` list argument and from `...`, validated.
.unitDeclarationsFromArgs <- function(dots, unitsArg) {
  declared <- list()
  if (!is.null(unitsArg)) {
    checkmate::assertList(unitsArg, types = "character", names = "unique")
    for (nm in names(unitsArg)) {
      declared[[nm]] <- unitsArg[[nm]]
    }
  }
  if (length(dots) == 0) {
    return(declared)
  }
  if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
    .unitStop("Every unit passed through `...` must be named, e.g. `depot = \"mg\"`.")
  }
  if (anyDuplicated(names(dots))) {
    dup <- unique(names(dots)[duplicated(names(dots))])
    .unitStop(sprintf("Duplicated unit names: %s.", paste(dup, collapse = ", ")))
  }
  for (nm in names(dots)) {
    checkmate::assertString(dots[[nm]], .var.name = nm)
    declared[[nm]] <- dots[[nm]]
  }
  declared
}

# Merge the caller's units over the model's own metadata: `...` over the
# `units` argument over `covariateData` over the `units` metadata.
.unitDeclarations <- function(ui, dots, unitsArg, st, known) {
  meta <- as.list(ui$meta)
  declared <- .unitDeclarationsFromMeta(meta$units, st, known)
  fromCov <- .unitDeclarationsFromCovariateData(meta$covariateData)
  for (nm in setdiff(names(fromCov), names(declared))) {
    declared[[nm]] <- fromCov[[nm]]
  }
  fromArgs <- .unitDeclarationsFromArgs(dots, unitsArg)
  for (nm in names(fromArgs)) {
    declared[[nm]] <- fromArgs[[nm]]
  }
  if ("linCmt" %in% names(declared)) {
    declared[["rxLinCmt"]] <- declared[["linCmt"]]
    declared[["linCmt"]] <- NULL
  }
  declared
}

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
  if ("depot" %in% states) {
    return("depot")
  }
  if ("central" %in% states) {
    return("central")
  }
  character()
}

# ---- checkUnits() ----------------------------------------------------------------

#' Check and infer the physical units of an rxode2 / nlmixr2 model
#'
#' `checkUnits()` walks the `model({})` block and works out a unit for every
#' symbol from the units declared at the model's boundary: `time`, the dose
#' unit of each dosed compartment (named by the compartment, `depot = "mg"`),
#' the unit of each output variable (`Cc = "ng/mL"`), and any per-parameter
#' override (`cl = "mL/min"`). Parameters that the equations determine only up
#' to a convertible unit are written in the model's own representative unit
#' for that dimension: `L` and `L/<time>` for a mass or molar dose, `mL/kg`
#' and `mL/kg/<time>` for a per-body-weight dose, `1/<time>` for rate
#' constants. `addUnits()` applies the result to the model.
#'
#' Units declared in the model's `units` metadata are the defaults; a legacy
#' `list(time =, dosing =, concentration =)` block is read as the dose unit of
#' the dosed compartments and the unit of a single endpoint. Arguments override
#' the metadata. A covariate without a unit (in `covariateData` or the
#' arguments) is taken as compatible with whatever it is combined with.
#'
#' The rules the walker applies:
#' * `+`, `-`, `*`, `/` and comparisons validate their operands; a number added
#'   to or subtracted from a quantity, or dividing a covariate, is a value in
#'   that quantity's unit (`AGE - 40`, `WT / 70`), and any other bare number is
#'   unitless.
#' * `exp()`, `log()`, `sqrt()` and the other transcendental functions give a
#'   unitless result; the arithmetic inside their argument is validated but
#'   the argument's own unit is not judged. `exp()` (and `expit()`,
#'   `probitInv()`, `10^x`) of nothing but `ini()` parameters and covariates is a
#'   back-transform: the assigned variable keeps its own unit and the
#'   parameters are unitless with `transformOf` naming that variable. A
#'   dimensioned quantity raised to a symbolic power (`Cc^hill`) is unknown;
#'   raised to a literal fractional power it is an issue.
#' * A variable assigned more than once (in both branches of an `if`, say) must
#'   have the same unit each time.
#' * Residual-error parameters take the unit their distribution implies:
#'   `add()` the endpoint's unit, `prop()` and the log/logit/probit families
#'   unitless.
#'
#' The `units` package, with its udunits-2 system library, does the parsing and
#' conversions. It is a suggested dependency; both functions error with an
#' installation hint when it is absent. On first use a few units udunits lacks
#' are installed for the session (`unitless`, `IU`, `U`, `M` for molar, `mcg`,
#' `cells`, `CFU`, `copies`, `molecules`); this is process-global state.
#'
#' @param ui A model: an `rxUi`, a model function, or anything
#'   [rxode2::assertRxUi()] accepts.
#' @param ... Named unit strings, one per symbol: `time`, a dosed compartment,
#'   an endpoint variable, an `ini()` parameter, a covariate, or a variable
#'   assigned in `model({})`.
#' @param units The same as `...`, as a named list. `...` takes precedence.
#' @return A `data.frame` with one row per symbol and the columns `name`,
#'   `type` (`time`, `state`, `parameter`, `eta`, `residual`, `covariate`,
#'   `intermediate`, `endpoint`, `property`, `reserved`), `unit` (a canonical
#'   string, `"unitless"`, or `NA` when unresolved), `source` (`declared`,
#'   `boundary`, `inferred`, `default`, `unresolved`), `transformOf`, `line`
#'   (the defining model line), `conversion` (the constant the line needs, as
#'   `"mg/L * 1000 = ng/mL"`, or `NA`), and `issue` (`NA` when the symbol is
#'   consistent). The attribute `"conversions"` is a data frame that
#'   [addUnits()] applies, and `"notes"` lists metadata that was ignored.
#' @family units
#' @export
#' @author Bill Denney
#' @examples
#' if (requireNamespace("units", quietly = TRUE)) {
#'   mod <- readModelDb("PK_1cmt_des")
#'   checkUnits(mod, time = "h", depot = "mg", Cc = "ng/mL")
#' }
checkUnits <- function(ui, ..., units = NULL) {
  .assertUnitsInstalled()
  .unitInstallCustom()
  # Unit arithmetic must not fold `mg/kg` into a dimensionless 1e-6: a dose
  # per body weight is its own dimension here. .unitSymbols() cancels only
  # identical tokens.
  oldSimplify <- units::units_options("simplify")
  units::units_options(simplify = FALSE)
  on.exit(units::units_options(simplify = oldSimplify), add = TRUE)
  .ui <- rxode2::assertRxUi(ui)
  dots <- list(...)
  st <- .unitState()
  st$unitless <- .unitOne("unitless")
  states <- .ui$state
  if (is.null(states)) {
    states <- character()
  }
  predDf <- .ui$predDf
  endpoints <- if (is.data.frame(predDf)) unique(predDf$var) else character()
  st$endpoints <- endpoints
  st$dosedDefault <- .unitDosedCompartments(.ui)
  ini <- .ui$iniDf
  if (!is.data.frame(ini)) {
    ini <- data.frame(name = character(), neta1 = numeric(), condition = character(), err = character())
  }
  thetas <- ini$name[is.na(ini$neta1) & is.na(ini$err)]
  etas <- ini$name[!is.na(ini$neta1) & ini$neta1 == ini$neta2]
  residuals <- ini[!is.na(ini$err), c("name", "condition", "err"), drop = FALSE]
  names(residuals) <- c("name", "endpoint", "err")
  st$residuals <- residuals
  covs <- .ui$allCovs
  if (is.null(covs)) {
    covs <- character()
  }
  lstExpr <- .ui$lstExpr
  constraints <- .unitConstraints(lstExpr)
  intermediates <- unique(vapply(
    Filter(function(cs) cs$kind == "assign", constraints),
    function(cs) cs$target,
    character(1)
  ))
  intermediates <- setdiff(intermediates, c(states, endpoints))
  reserved <- rxode2::rxReservedKeywords[["Reserved Name"]]
  reserved <- reserved[!is.na(reserved)]
  timeNames <- c("t", "time", "tlast")
  constants <- setdiff(reserved, c(timeNames, "podo"))

  # symbol classes
  for (s in states) {
    st$type[s] <- "state"
  }
  for (s in endpoints) {
    st$type[s] <- "endpoint"
  }
  for (s in thetas) {
    st$type[s] <- "parameter"
  }
  for (s in etas) {
    st$type[s] <- "eta"
  }
  for (s in residuals$name) {
    st$type[s] <- "residual"
  }
  for (s in covs) {
    st$type[s] <- "covariate"
  }
  for (s in intermediates) {
    st$type[s] <- "intermediate"
  }
  for (s in c(timeNames, "podo", constants)) {
    st$type[s] <- "reserved"
  }
  for (cs in constraints) {
    if (cs$kind != "check") {
      st$defLine[[cs$target]] <- c(st$defLine[[cs$target]], cs$line)
    }
  }

  known <- c("time", "linCmt", names(st$type))
  declared <- .unitDeclarations(.ui, dots, units, st, known)
  bad <- setdiff(names(declared), known)
  if (length(bad) > 0) {
    .unitStop(sprintf(
      "Unit(s) given for names that are not in the model: %s.",
      paste(bad, collapse = ", ")
    ))
  }
  for (nm in names(declared)) {
    u <- tryCatch(.unitOne(declared[[nm]]), nlmixr2libUnitError = function(e) {
      .unitStop(sprintf("%s: %s", nm, conditionMessage(e)))
    })
    if (nm == "time") {
      st$time <- u
      next
    }
    .unitSet(st, nm, u, "declared")
  }
  if (!is.null(st$time)) {
    for (s in timeNames) {
      .unitSet(st, s, st$time, "boundary")
    }
  }
  for (s in constants) {
    .unitSet(st, s, st$unitless, "boundary")
  }
  dosed <- intersect(names(declared), states)
  st$dosed <- dosed
  doseUnit <- if (length(dosed) > 0) .unitGet(st, dosed[[1]]) else NULL
  if (!is.null(doseUnit)) {
    .unitSet(st, "podo", doseUnit, "boundary")
  }
  endpointUnit <- NULL
  for (s in endpoints) {
    if (.unitKnown(st, s)) {
      endpointUnit <- .unitGet(st, s)
      break
    }
  }
  st$candidates <- .unitCandidates(st, doseUnit, endpointUnit)
  st$anySyms <- setdiff(covs, names(declared))

  .unitInfer(constraints, st)
  for (s in etas) {
    if (!.unitKnown(st, s)) {
      .unitSet(st, s, st$unitless, "default")
    }
  }
  .unitValidate(constraints, st)
  .unitResultFrame(st, constraints, states, endpoints, thetas, etas, residuals$name, covs, intermediates, lstExpr)
}

# Assemble the result data frame.
.unitResultFrame <- function(
  st,
  constraints,
  states,
  endpoints,
  thetas,
  etas,
  residualNames,
  covs,
  intermediates,
  lstExpr
) {
  usedReserved <- intersect(
    names(st$type)[st$type == "reserved"],
    unique(unlist(lapply(lstExpr, .unitSymbolsIn)))
  )
  props <- vapply(
    Filter(function(cs) cs$kind %in% c("ddt", "rate", "f", "time"), constraints),
    function(cs) cs$target,
    character(1)
  )
  order <- c(
    "time",
    unique(c(st$dosed, states)),
    endpoints,
    thetas,
    residualNames,
    etas,
    covs,
    intermediates,
    unique(props),
    usedReserved
  )
  order <- unique(order)
  convByTarget <- list()
  for (cv in st$conversions) {
    convByTarget[[cv$target]] <- c(
      convByTarget[[cv$target]],
      .unitConversionText(cv$from, cv$factor, cv$to)
    )
  }
  rows <- lapply(order, function(nm) {
    if (nm == "time") {
      u <- st$time
      type <- "time"
      source <- if (is.null(u)) "unresolved" else "declared"
    } else if (nm %in% names(st$type)) {
      type <- st$type[[nm]]
      if (type == "reserved") {
        u <- .unitGet(st, nm)
      } else if (grepl("^d/dt\\(|^(f|F|rate|dur|alag|lag)\\(", nm)) {
        u <- NULL
      } else {
        u <- .unitGet(st, nm)
      }
      source <- if (!is.null(u)) {
        st$source[[nm]]
      } else if (nm %in% st$anySyms) {
        "compatible"
      } else {
        "unresolved"
      }
    } else {
      type <- "property"
      cs <- Filter(function(c) identical(c$target, nm), constraints)[[1]]
      u <- .unitTargetUnit(cs, st)
      source <- if (is.null(u)) "unresolved" else "boundary"
    }
    line <- st$defLine[[nm]]
    data.frame(
      name = nm,
      type = type,
      unit = if (is.null(u)) NA_character_ else .unitDeparse(u),
      source = source,
      transformOf = if (nm %in% names(st$transformOf)) st$transformOf[[nm]] else NA_character_,
      line = if (is.null(line)) NA_integer_ else as.integer(line[[1]]),
      conversion = if (is.null(convByTarget[[nm]])) NA_character_ else paste(convByTarget[[nm]], collapse = "; "),
      issue = if (is.null(st$issues[[nm]])) NA_character_ else paste(st$issues[[nm]], collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # issues recorded against a line rather than a symbol (an `if` condition)
  lineIssues <- st$issues[grepl("^line [0-9]+$", names(st$issues))]
  if (length(lineIssues) > 0) {
    extra <- data.frame(
      name = names(lineIssues),
      type = "line",
      unit = NA_character_,
      source = NA_character_,
      transformOf = NA_character_,
      line = as.integer(sub("^line ", "", names(lineIssues))),
      conversion = NA_character_,
      issue = vapply(lineIssues, paste, character(1), collapse = "; "),
      stringsAsFactors = FALSE
    )
    out <- rbind(out, extra)
  }
  rownames(out) <- NULL
  conv <- if (length(st$conversions) == 0) {
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
  } else {
    do.call(rbind, st$conversions)
  }
  attr(out, "conversions") <- conv
  attr(out, "notes") <- st$notes
  out
}

# "mg/L * 1000 = ng/mL" or "1/day / 24 = 1/h".
.unitConversionText <- function(from, factor, to) {
  op <- .unitConversionOperator(factor)
  sprintf("%s %s %s = %s", from, op$op, format(op$value, digits = 12, scientific = FALSE), to)
}

# Write a factor below 1 as a division when its reciprocal is a clean number.
.unitConversionOperator <- function(factor) {
  factor <- signif(factor, 12)
  if (factor < 1) {
    inv <- signif(1 / factor, 12)
    if (abs(inv - round(inv)) < 1e-8) {
      return(list(op = "/", value = round(inv)))
    }
  }
  list(op = "*", value = factor)
}

# ---- addUnits() --------------------------------------------------------------------

#' Add units to a model, inserting the conversions its arithmetic needs
#'
#' Runs [checkUnits()] and, when it reports no issue, writes the units into the
#' model: every conversion the arithmetic needs is inserted as a bare constant
#' (`Cc <- central/vc * 1000`) and explained in the `unitConversions`
#' metadata, the `units` metadata lists the unit of every resolved symbol
#' (`"unitless"` for dimensionless ones), `dosing` names the dosed
#' compartments, and a covariate given a unit gets it in `covariateData`.
#' `ini()` estimates are never rescaled: declaring `cl = "mL/min"` for a model
#' written in hours inserts the `/ 60` the equations need and leaves the
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
#'   mod <- readModelDb("PK_1cmt_des")
#'   res <- addUnits(mod, time = "h", depot = "mg", Cc = "ng/mL")
#'   res$meta$units
#'   res$meta$unitConversions
#' }
addUnits <- function(ui, ..., units = NULL) {
  .ui <- rxode2::assertRxUi(ui)
  res <- checkUnits(.ui, ..., units = units)
  issues <- res[!is.na(res$issue), , drop = FALSE]
  if (nrow(issues) > 0) {
    .unitStop(paste0(
      "checkUnits() found unit issues; fix them before adding units:\n",
      paste(sprintf("  %s: %s", issues$name, issues$issue), collapse = "\n")
    ))
  }
  conv <- attr(res, "conversions")
  lst <- .ui$lstExpr
  for (i in seq_len(nrow(conv))) {
    path <- as.integer(strsplit(conv$path[i], ",", fixed = TRUE)[[1]])
    lst[[conv$line[i]]] <- .unitInsertConversion(lst[[conv$line[i]]], path, conv$termIndex[i], conv$factor[i])
  }
  resolved <- res[!is.na(res$unit) & res$type != "property" & res$type != "reserved", , drop = FALSE]
  unitsMeta <- stats::setNames(as.list(resolved$unit), resolved$name)
  .ui <- rxode2::rxUiDecompress(.ui)
  meta <- .ui$meta
  assign("units", unitsMeta, envir = meta)
  existingDosing <- if (exists("dosing", envir = meta)) get("dosing", envir = meta) else character()
  dosing <- unique(c(existingDosing, res$name[res$type == "state" & res$source == "declared"]))
  dosing <- intersect(dosing, .ui$state)
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
  covRows <- res[res$type == "covariate" & !is.na(res$unit), , drop = FALSE]
  if (nrow(covRows) > 0) {
    covData <- if (exists("covariateData", envir = meta)) get("covariateData", envir = meta) else list()
    for (i in seq_len(nrow(covRows))) {
      nm <- covRows$name[i]
      entry <- covData[[nm]]
      if (!is.list(entry)) {
        entry <- list(description = NA_character_, units = covRows$unit[i], type = "continuous")
      } else {
        entry$units <- covRows$unit[i]
      }
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

# Multiply the `termIndex`-th additive term of the assignment at `path`
# inside `e` by `factor`, editing the expression rather than its text.
.unitInsertConversion <- function(e, path, termIndex, factor) {
  if (length(path) == 0) {
    e[[3]] <- .unitScaleTerm(e[[3]], termIndex, factor)
    return(e)
  }
  e[[path[[1]]]] <- .unitInsertConversion(e[[path[[1]]]], path[-1], termIndex, factor)
  e
}

.unitScaleTerm <- function(rhs, termIndex, factor) {
  acc <- new.env(parent = emptyenv())
  acc$n <- 0L
  .unitScaleWalk(rhs, termIndex, factor, acc)
}

.unitScaleWalk <- function(e, termIndex, factor, acc) {
  inner <- .stripParens(e)
  if (is.call(inner) && is.name(inner[[1]]) && as.character(inner[[1]]) %in% c("+", "-") && length(inner) == 3) {
    inner[[2]] <- .unitScaleWalk(inner[[2]], termIndex, factor, acc)
    inner[[3]] <- .unitScaleWalk(inner[[3]], termIndex, factor, acc)
    return(inner)
  }
  acc$n <- acc$n + 1L
  if (acc$n != termIndex) {
    return(e)
  }
  op <- .unitConversionOperator(factor)
  if (op$op == "/") {
    return(call("/", e, op$value))
  }
  call("*", e, op$value)
}
