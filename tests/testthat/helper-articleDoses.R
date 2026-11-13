# The compartments an article doses, read from its code without running it.
#
# An article's R chunks are extracted with knitr::purl() and parsed, and
# magrittr pipes are rewritten into nested calls, the form the native pipe
# already parses to. A call is a dose record when it has an `evid` of 1 or 4
# (a dose, or a reset and dose), or else an `amt` that is not a literal 0 or
# NA. Its compartment is its own `cmt` argument or, when it has none, the
# `cmt` of the call piped into it (`crossing(cmt = c("depot", "depot2")) |>
# mutate(amt = dose, evid = 1)`).
#
# A `cmt` argument names compartments through the string literals in it
# (`"depot"`, `ifelse(oral, "depot", "central")`) and through every variable
# it reads, which stands for every string ever assigned to that variable:
# `cmt <- if (iv) "venous" else "stomach"` anywhere in the article, a
# `for (cmt in c(...))` loop, or a `function(cmt = "depot")` default. Where
# `cmt` and `evid` (or `amt`) are literal vectors of one length, only the
# dose positions count. A numeric `cmt` is returned as "#<n>".
#
# What is not seen: a compartment passed positionally into a helper function,
# or built some other way. The result is a lower bound.

# The R code of an article, parsed, with magrittr pipes nested.
.articleExprs <- function(rmd) {
  out <- tempfile(fileext = ".R")
  on.exit(unlink(out))
  # purl() evaluates chunk options, and one that reads an article variable
  # (`eval = have_data`) prints an error it then recovers from.
  utils::capture.output(
    invisible(suppressMessages(knitr::purl(rmd, output = out, quiet = TRUE, documentation = 0L))),
    type = "message"
  )
  lapply(parse(out, keep.source = FALSE), .unpipe)
}

# `lhs %>% f(x)` becomes `f(lhs, x)`; `lhs %>% f` becomes `f(lhs)`.
.unpipe <- function(e) {
  if (!is.call(e)) {
    return(e)
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) e[[i]] <- .unpipe(e[[i]])
  }
  if (identical(e[[1]], as.name("%>%")) && length(e) == 3L) {
    rhs <- e[[3]]
    if (is.call(rhs)) {
      return(as.call(c(rhs[[1]], e[[2]], as.list(rhs)[-1])))
    }
    if (is.name(rhs)) {
      return(as.call(list(rhs, e[[2]])))
    }
  }
  e
}

# Every string literal in an expression.
.stringsIn <- function(e) {
  if (is.character(e)) {
    return(e)
  }
  if (!is.call(e)) {
    return(character(0))
  }
  out <- character(0)
  # `e[[i]]` is tested in place: bound to a variable, an empty argument such
  # as the one in `x[, 1]` errors when the variable is read.
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]]) || is.character(e[[i]])) out <- c(out, .stringsIn(e[[i]]))
  }
  out
}

# The compartments an expression can evaluate to: its string literals, and
# its whole-number literals as "#<n>". The conditions of `if`, ifelse() and
# if_else(), comparisons, sequences and pasted strings are skipped:
# `if (route == "iv") "central" else "depot"` names central and depot, not
# "iv".
.cmtLiteralsIn <- function(e) {
  if (is.character(e)) {
    return(e)
  }
  if (is.numeric(e)) {
    return(if (length(e) == 1L && !is.na(e) && e >= 1 && e == round(e)) paste0("#", e) else character(0))
  }
  skipCalls <- c("==", "!=", "<", ">", "<=", ">=", "%in%", ":", "seq", "seq_len", "paste", "paste0", "sprintf")
  if (!is.call(e) || .isCallTo(e, skipCalls)) {
    return(character(0))
  }
  skip <- if (identical(e[[1]], as.name("if")) || .isCallTo(e, c("ifelse", "if_else"))) 2L else integer(0)
  out <- character(0)
  for (i in setdiff(seq_along(e)[-1], skip)) {
    if (is.call(e[[i]]) || is.character(e[[i]]) || is.numeric(e[[i]])) out <- c(out, .cmtLiteralsIn(e[[i]]))
  }
  out
}

# Every symbol an expression reads, including the names of the functions it
# calls (an article's own helper functions are variables too).
.symbolsRead <- function(e) {
  if (is.name(e)) {
    return(as.character(e))
  }
  if (!is.call(e)) {
    return(character(0))
  }
  out <- if (is.name(e[[1]])) as.character(e[[1]]) else character(0)
  # `x$name` and `x@name` read `x`, not a variable called `name`.
  idx <- seq_along(e)[-1]
  if (identical(e[[1]], as.name("$")) || identical(e[[1]], as.name("@"))) {
    idx <- 2L
  }
  for (i in idx) {
    if (is.name(e[[i]]) && nzchar(as.character(e[[i]]))) {
      out <- c(out, as.character(e[[i]]))
    } else if (is.call(e[[i]])) {
      out <- c(out, .symbolsRead(e[[i]]))
    }
  }
  out
}

# Adds `strings` to `vars[[name]]`.
.addVar <- function(vars, name, strings) {
  if (length(strings)) {
    vars[[name]] <- unique(c(vars[[name]], strings))
  }
  vars
}

# The compartments a value assigned to a variable can name: its strings, and
# its numbers only when it is a literal or a choice between literals (`6L`,
# `c(1, 2)`, `ifelse(iv, 6L, 1L)`), since numbers elsewhere are rarely
# compartments.
.assignedCmts <- function(rhs) {
  choice <- is.numeric(rhs) ||
    .isCallTo(rhs, c("c", "ifelse", "if_else", "switch")) ||
    is.call(rhs) && identical(rhs[[1]], as.name("if"))
  found <- .cmtLiteralsIn(rhs)
  if (choice) found else found[!startsWith(found, "#")]
}

# The compartments assigned to each variable anywhere in an expression (see
# .assignedCmts()), in the order written: by `<-` or `=` (with the strings of
# the variables the value reads), as a `for` loop's values, or as a function
# argument's default.
.stringVarsIn <- function(e, vars = list()) {
  if (!is.call(e)) {
    return(vars)
  }
  head <- e[[1]]
  isAssign <- identical(head, as.name("<-")) || identical(head, as.name("="))
  if (isAssign && length(e) == 3L && is.name(e[[2]])) {
    # A row or column taken from a table (`row <- arms[i, ]`) can name what
    # the table names.
    inherited <- as.character(unlist(vars[intersect(.symbolsRead(e[[3]]), names(vars))], use.names = FALSE))
    inherited <- inherited[!startsWith(inherited, "#")]
    vars <- .addVar(vars, as.character(e[[2]]), c(.assignedCmts(e[[3]]), inherited))
  } else if (identical(head, as.name("for"))) {
    vars <- .addVar(vars, as.character(e[[2]]), .assignedCmts(e[[3]]))
  } else if (identical(head, as.name("function"))) {
    fm <- e[[2]]
    for (nm in names(fm)) {
      if (is.call(fm[[nm]]) || is.character(fm[[nm]]) || is.numeric(fm[[nm]])) {
        vars <- .addVar(vars, nm, .assignedCmts(fm[[nm]]))
      }
    }
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) vars <- .stringVarsIn(e[[i]], vars)
  }
  vars
}

.isLiteralScalar <- function(x) {
  is.atomic(x) && length(x) == 1L
}

.isCallTo <- function(x, fns) {
  is.call(x) &&
    (is.name(x[[1]]) &&
      as.character(x[[1]]) %in% fns ||
      is.call(x[[1]]) && identical(x[[1]][[1]], as.name("::")) && as.character(x[[1]][[3]]) %in% fns)
}

# One element of a column: a literal scalar, `rep(<literal>, n)`, or an
# optional element `if (cond) <literal> else numeric(0)`; NULL when it is
# anything else.
.segmentValue <- function(x) {
  if (.isLiteralScalar(x)) {
    return(x)
  }
  if (.isCallTo(x, c("rep", "rep_len")) && length(x) >= 2L && .isLiteralScalar(x[[2]])) {
    return(x[[2]])
  }
  isOptional <- is.call(x) &&
    identical(x[[1]], as.name("if")) &&
    length(x) == 4L &&
    .isCallTo(x[[4]], c("numeric", "integer", "character", "logical")) &&
    length(x[[4]]) == 2L
  if (isOptional && .isLiteralScalar(x[[3]])) {
    return(x[[3]])
  }
  NULL
}

# One element of a `c()` column as a list of segments: `rep(x, k)` with a
# small literal `k` is k copies of `x`, so `c(rep("depot", 3), "central")`
# lines up with `c(1, 1, 1, 0)`.
.expandRep <- function(x) {
  isRep <- .isCallTo(x, "rep") && length(x) == 3L && .argNames(as.list(x))[3] %in% c("", "times")
  k <- if (isRep) x[[3]]
  if (is.numeric(k) && length(k) == 1L && k >= 1 && k <= 50 && k == round(k)) {
    return(rep(list(x[[2]]), k))
  }
  list(x)
}

# A record column as a list of segments (`values`, with the expressions they
# came from in `exprs`), with a key that says how the segments line up with
# another column's:
# - a literal or `rep(<literal>, n)`: one segment, key "1";
# - `c(...)`, or `rep(c(...), n)`: one segment per argument, key "c<n>";
# - `ifelse(cond, a, b)`, `if_else()` or `if (cond) a else b`: two segments,
#   keyed by the condition.
# A segment that is not literal is NULL. Returns NULL for any other form.
.segments <- function(x) {
  v <- .segmentValue(x)
  if (!is.null(v)) {
    return(list(values = list(v), exprs = list(x), key = "1"))
  }
  if (.isCallTo(x, c("rep", "rep_len")) && length(x) >= 2L && .isCallTo(x[[2]], "c")) {
    x <- x[[2]]
  }
  if (.isCallTo(x, "c") && length(x) > 1L) {
    ex <- unlist(lapply(as.list(x)[-1], .expandRep), recursive = FALSE)
    return(list(values = lapply(ex, .segmentValue), exprs = ex, key = paste0("c", length(ex))))
  }
  isIf <- is.call(x) && identical(x[[1]], as.name("if")) && length(x) == 4L
  if (isIf || .isCallTo(x, c("ifelse", "if_else")) && length(x) == 4L) {
    ex <- list(x[[3]], x[[4]])
    return(list(values = lapply(ex, .segmentValue), exprs = ex, key = paste(deparse(x[[2]]), collapse = "")))
  }
  NULL
}

# Can an `evid` segment be a dose (1 or 4)? A segment that is not literal can.
.evidSegmentIsDose <- function(v) {
  is.null(v) || isTRUE(v %in% c(1, 4))
}

# Can an `amt` segment be a dose (not 0 or NA)? A segment that is not literal can.
.amtSegmentIsDose <- function(v) {
  is.null(v) || !(is.na(v) || isTRUE(v == 0))
}

# Which segments of a record can be doses: the `evid` or `amt` column's
# segments (see .segments()) with a logical `mask`; NULL when the record has
# neither column in a form .segments() reads. Where both columns are read, a
# segment is a dose only if both allow it; when their segments do not line
# up, and neither is a single segment, `evid` decides.
.doseSegments <- function(e) {
  ev <- if ("evid" %in% names(e)) .segments(e[["evid"]])
  am <- if ("amt" %in% names(e)) .segments(e[["amt"]])
  if (is.null(am) && .isEtDose(e)) {
    am <- .segments(e[["dose"]])
  }
  if (!is.null(ev)) {
    ev$mask <- vapply(ev$values, .evidSegmentIsDose, logical(1))
  }
  if (!is.null(am)) {
    am$mask <- vapply(am$values, .amtSegmentIsDose, logical(1))
  }
  if (is.null(ev) || is.null(am)) {
    return(if (is.null(ev)) am else ev)
  }
  if (identical(ev$key, am$key) || length(am$mask) == 1L) {
    ev$mask <- ev$mask & am$mask
    return(ev)
  }
  if (length(ev$mask) == 1L) {
    am$mask <- am$mask & ev$mask
    return(am)
  }
  ev
}

# Does an event-table call give its amount as `dose =`, which et() and
# add.dosing() take for `amt`?
.isEtDose <- function(e) {
  "dose" %in% names(e) && !("amt" %in% names(e)) && .isCallTo(e, c("et", "add.dosing"))
}

# Is a call a dose record? It has an `evid` or `amt` column (or an event
# table's `dose`), and no segment of it is ruled out as a dose by literal
# values.
.isDoseCall <- function(e) {
  nms <- names(e)
  if (is.null(nms) || !any(c("evid", "amt") %in% nms) && !.isEtDose(e)) {
    return(FALSE)
  }
  ds <- .doseSegments(e)
  is.null(ds) || any(ds$mask)
}

# The `cmt` argument of a call, or of the first call piped into it. The walk
# stops at an event table (`et()`): its records carry their own `amt`, so a
# `cmt` there without one belongs to observations.
.recordCmt <- function(e) {
  if ("cmt" %in% names(e)) {
    return(e[["cmt"]])
  }
  while (length(e) >= 2L && is.call(e[[2]])) {
    e <- e[[2]]
    if (.isCallTo(e, "et")) {
      return(NULL)
    }
    if ("cmt" %in% names(e)) {
      return(e[["cmt"]])
    }
  }
  NULL
}

# The compartments one dose record names, with how each was read: "literal"
# (a string or number in the record's `cmt`) or "variable" (a string assigned
# to a variable the `cmt` reads).
.recordCmts <- function(e, vars) {
  cm <- .recordCmt(e)
  if (is.null(cm)) {
    return(NULL)
  }
  seg <- .segments(cm)
  ds <- .doseSegments(e)
  parts <- list(cm)
  lead <- .leadingDoses(ds)
  if (!is.null(seg) && !is.null(ds) && identical(seg$key, ds$key)) {
    parts <- seg$exprs[ds$mask]
  } else if (!is.null(seg) && startsWith(seg$key, "c") && lead > 0L && length(seg$exprs) > lead) {
    # Dose rows first, then observation rows whose segments the columns
    # split differently: `amt = c(d, rep(NA, 2 * n))` with
    # `cmt = c("depot", rep("a", n), rep("b", n))`.
    parts <- seg$exprs[seq_len(lead)]
  } else if (!is.null(seg) && all(!vapply(seg$values, is.null, logical(1)))) {
    parts <- seg$exprs
  }
  out <- lapply(parts, .cmtsOfExpr, vars = vars)
  do.call(rbind, c(list(data.frame(cmt = character(0), via = character(0))), out))
}

# How many leading segments of a `c()` record are doses, when the record is
# those doses followed only by non-doses; 0 otherwise.
.leadingDoses <- function(ds) {
  if (is.null(ds) || !startsWith(ds$key, "c")) {
    return(0L)
  }
  k <- match(FALSE, ds$mask) - 1L
  if (is.na(k) || k == 0L || any(ds$mask[-seq_len(k)])) 0L else k
}

# The compartments one `cmt` expression names: its literals ("literal"), and
# the values of the variables it reads ("variable").
.cmtsOfExpr <- function(x, vars) {
  lit <- .cmtLiteralsIn(x)
  viaVar <- unlist(vars[intersect(.symbolsRead(x), names(vars))], use.names = FALSE)
  data.frame(
    cmt = c(lit, viaVar),
    via = rep(c("literal", "variable"), c(length(lit), length(viaVar)))
  )
}

# Is an argument a tribble column header (`~name`)?
.isColumnHeader <- function(a) {
  is.call(a) && identical(a[[1]], as.name("~")) && length(a) == 2L && is.name(a[[2]])
}

# The rows of a `tribble(~a, ~b, 1, 2, ...)` call, each as a
# `data.frame(a = 1, b = 2)` call; an empty list for any other call.
.tribbleRows <- function(e) {
  if (!.isCallTo(e, c("tribble", "frame_data"))) {
    return(list())
  }
  args <- as.list(e)[-1]
  nCol <- match(FALSE, vapply(args, .isColumnHeader, logical(1))) - 1L
  if (is.na(nCol) || nCol == 0L || (length(args) - nCol) %% nCol != 0L) {
    return(list())
  }
  cols <- vapply(args[seq_len(nCol)], all.vars, character(1))
  vals <- args[-seq_len(nCol)]
  rows <- list()
  for (r in seq_len(length(vals) / nCol)) {
    cells <- vals[(r - 1L) * nCol + seq_len(nCol)]
    names(cells) <- cols
    rows[[r]] <- as.call(c(as.name("data.frame"), cells))
  }
  rows
}

# Every dose record in an expression: its compartments, how each was read,
# and the record's code.
.doseRecordsIn <- function(e, vars, acc = list()) {
  if (!is.call(e)) {
    return(acc)
  }
  for (row in .tribbleRows(e)) {
    acc <- .doseRecordsIn(row, vars, acc)
  }
  if (.isDoseCall(e)) {
    rc <- .recordCmts(e, vars)
    if (NROW(rc)) {
      rc$record <- paste(deparse(e, width.cutoff = 500L, nlines = 1L), collapse = "")
      acc[[length(acc) + 1L]] <- rc
    }
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) acc <- .doseRecordsIn(e[[i]], vars, acc)
  }
  acc
}

# The names of a call's arguments, "" for an unnamed one.
.argNames <- function(args) {
  nms <- names(args)
  if (is.null(nms)) rep("", length(args)) else nms
}

# The argument a call passes for parameter `param` of a function with
# arguments `formals`: by exact name, or else by position among the
# parameters not named in the call and before any `...`. NULL if not passed.
.argFor <- function(e, formals, param) {
  args <- as.list(e)[-1]
  nms <- .argNames(args)
  if (param %in% nms) {
    return(args[[param]])
  }
  free <- setdiff(names(formals), nms)
  if ("..." %in% free) {
    free <- free[seq_len(match("...", free) - 1L)]
  }
  k <- match(param, free)
  positional <- args[!nzchar(nms)]
  if (!is.na(k) && k <= length(positional)) positional[[k]]
}

# The variable a dose record's `cmt` names directly (`cmt = x`, `cmt = x[i]`
# or `cmt = x[[i]]`), or NULL.
.recordCmtVariable <- function(r) {
  cm <- .recordCmt(r)
  if (is.call(cm) && (identical(cm[[1]], as.name("[")) || identical(cm[[1]], as.name("[[")))) {
    cm <- cm[[2]]
  }
  if (is.name(cm)) as.character(cm)
}

# The column a `V$col <- x` or `V[["col"]] <- x` statement assigns, as
# list(var, col), or NULL.
.columnAssignment <- function(st) {
  isAssign <- is.call(st) && (identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))
  if (!isAssign || length(st) != 3L || !is.call(st[[2]])) {
    return(NULL)
  }
  tgt <- st[[2]]
  ok <- length(tgt) == 3L &&
    is.name(tgt[[2]]) &&
    (identical(tgt[[1]], as.name("$")) &&
      is.name(tgt[[3]]) ||
      identical(tgt[[1]], as.name("[[")) && is.character(tgt[[3]]))
  if (ok) list(var = as.character(tgt[[2]]), col = as.character(tgt[[3]]))
}

# Dose records assembled column by column in a sequence of statements:
# `d <- data.frame(time = 0, amt = dose); d$cmt <- "depot"`. For each
# variable given an `amt`, `evid` or `cmt` column this way, the columns of the
# call that created it and those assigned after are put together into one
# data.frame() call, emitted when the variable is assigned anew or the
# statements end. Returns a list of list(at, record), `at` being the
# statement that assigned the last of those columns.
.assembleColumns <- function(stmts) {
  cols <- list()
  at <- integer(0)
  recs <- list()
  for (k in seq_along(stmts)) {
    st <- stmts[[k]]
    ca <- .columnAssignment(st)
    lhs <- .assignedName(st)
    if (!is.null(ca)) {
      if (!is.list(cols[[ca$var]])) {
        cols[[ca$var]] <- list()
      }
      cols[[ca$var]][[ca$col]] <- st[[3]]
      if (ca$col %in% c("amt", "evid", "cmt")) at[[ca$var]] <- k
    } else if (!is.null(lhs)) {
      if (lhs %in% names(at)) {
        recs[[length(recs) + 1L]] <- list(at = at[[lhs]], record = as.call(c(as.name("data.frame"), cols[[lhs]])))
        at <- at[names(at) != lhs]
      }
      named <- if (is.call(st[[3]])) as.list(st[[3]])[-1] else list()
      cols[lhs] <- list(named[nzchar(.argNames(named))])
    }
  }
  for (v in names(at)) {
    recs[[length(recs) + 1L]] <- list(at = at[[v]], record = as.call(c(as.name("data.frame"), cols[[v]])))
  }
  recs
}

# Column-by-column dose records in the `{` blocks inside one statement.
.nestedColumnRecords <- function(e, acc = list()) {
  if (!is.call(e)) {
    return(acc)
  }
  if (identical(e[[1]], as.name("{"))) {
    for (r in .assembleColumns(as.list(e)[-1])) {
      acc[[length(acc) + 1L]] <- r$record
    }
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) acc <- .nestedColumnRecords(e[[i]], acc)
  }
  acc
}

# The article's own functions that build dose records from an argument: a
# top-level `F <- function(...)` with a dose record whose `cmt` is one of F's
# arguments (or an element of it). Returns, per function, its argument list
# and those arguments' names.
.doseFunctions <- function(exprs, vars) {
  out <- list()
  for (e in exprs) {
    lhs <- .assignedName(e)
    fn <- if (!is.null(lhs)) e[[3]]
    if (!is.call(fn) || !identical(fn[[1]], as.name("function"))) {
      next
    }
    recs <- .doseCallsIn(fn[[3]])
    used <- unique(unlist(lapply(recs, .recordCmtVariable)))
    params <- intersect(names(fn[[2]]), used)
    if (length(params)) out[[lhs]] <- list(formals = fn[[2]], params = params)
  }
  out
}

# Every dose record call in an expression.
.doseCallsIn <- function(e, acc = list()) {
  if (!is.call(e)) {
    return(acc)
  }
  if (.isDoseCall(e) && !is.null(.recordCmt(e))) {
    acc[[length(acc) + 1L]] <- e
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) acc <- .doseCallsIn(e[[i]], acc)
  }
  acc
}

# Dose records made by calling the article's own dose functions: a call
# `F(..., cmt = "depot")` doses what it passes for F's dose argument. Adds
# .recordCmts() rows, with the call as the record, to `acc`.
.doseFunctionCallsIn <- function(e, doseFns, vars, acc = list()) {
  if (!is.call(e)) {
    return(acc)
  }
  f <- if (is.name(e[[1]])) as.character(e[[1]])
  if (!is.null(f) && f %in% names(doseFns)) {
    for (p in doseFns[[f]]$params) {
      arg <- .argFor(e, doseFns[[f]]$formals, p)
      rc <- if (!is.null(arg)) .recordCmts(call("data.frame", amt = quote(dose), cmt = arg), vars)
      if (NROW(rc)) {
        rc$record <- paste(deparse(e, width.cutoff = 500L, nlines = 1L), collapse = "")
        acc[[length(acc) + 1L]] <- rc
      }
    }
  }
  for (i in seq_along(e)[-1]) {
    if (is.call(e[[i]])) acc <- .doseFunctionCallsIn(e[[i]], doseFns, vars, acc)
  }
  acc
}

# The models an expression refers to: by a model name written as a string,
# or through a variable that holds one (`mod_rat <- readModelDb("...")`).
.modelsReferencedIn <- function(e, models, modelVars) {
  viaVars <- unlist(modelVars[intersect(.symbolsRead(e), names(modelVars))], use.names = FALSE)
  intersect(unique(c(.stringsIn(e), viaVars)), models)
}

# Does an expression call `fn` anywhere?
.hasCallTo <- function(e, fn) {
  if (!is.call(e)) {
    return(FALSE)
  }
  if (.isCallTo(e, fn)) {
    return(TRUE)
  }
  for (i in seq_along(e)) {
    if (is.call(e[[i]]) && .hasCallTo(e[[i]], fn)) {
      return(TRUE)
    }
  }
  FALSE
}

# Does a top-level expression define a model in the article itself
# (`ref <- rxode2::rxode2({ ... })`, or a model function: no arguments, with
# `ini()` and `model()` blocks)? Doses after it, until a registry model is
# referred to again, are for that model.
.definesInlineModel <- function(e) {
  isAssign <- is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))
  if (!isAssign || length(e) != 3L) {
    return(FALSE)
  }
  rhs <- e[[3]]
  if (
    .isCallTo(rhs, c("rxode2", "rxode", "RxODE")) &&
      length(rhs) >= 2L &&
      is.call(rhs[[2]]) &&
      identical(rhs[[2]][[1]], as.name("{"))
  ) {
    return(TRUE)
  }
  is.call(rhs) &&
    identical(rhs[[1]], as.name("function")) &&
    !length(rhs[[2]]) &&
    .hasCallTo(rhs[[3]], "ini") &&
    .hasCallTo(rhs[[3]], "model")
}

# Updates the variable-to-model map with a top-level assignment whose
# right-hand side refers to models.
.recordModelVar <- function(e, refs, modelVars) {
  isAssign <- is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))
  if (isAssign && length(e) == 3L && is.name(e[[2]]) && length(refs)) {
    modelVars[[as.character(e[[2]])]] <- refs
  }
  modelVars
}

# The variable a top-level assignment binds, or NULL.
.assignedName <- function(e) {
  isAssign <- is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))
  if (isAssign && length(e) == 3L && is.name(e[[2]])) as.character(e[[2]])
}

# For each top-level expression: the registry models it refers to (`refs`),
# whether it defines or reads a model the article defines itself (`inline`),
# the variable it assigns (`lhs`) and the symbols it reads (`symbols`).
.articleModelRefs <- function(exprs, models) {
  modelVars <- list()
  inlineVars <- character(0)
  out <- vector("list", length(exprs))
  for (i in seq_along(exprs)) {
    e <- exprs[[i]]
    lhs <- .assignedName(e)
    # An inline model's own variables (`ntr <- mtt * ktr - 1`) can share a
    # name with an article variable that holds registry models.
    defines <- .definesInlineModel(e)
    refs <- if (defines) character(0) else .modelsReferencedIn(e, models, modelVars)
    modelVars <- .recordModelVar(e, refs, modelVars)
    symbols <- unique(.symbolsRead(if (is.null(lhs)) e else e[[3]]))
    out[[i]] <- list(
      refs = refs,
      lhs = lhs,
      symbols = symbols,
      inline = defines || any(symbols %in% inlineVars),
      defines = defines
    )
    if (defines) inlineVars <- c(inlineVars, lhs)
  }
  out
}

# The models referred to by the expressions after `i` that read the variable
# expression `i` assigns, directly or through variables assigned from it
# (`cohort <- function() design_events(...)`). A reading expression that
# refers to no registry model but reads an inline model is a solve of that
# model and adds none. NULL when no reading expression refers to a model.
.modelUses <- function(info, i) {
  reach <- info[[i]]$lhs
  if (is.null(reach) || i == length(info)) {
    return(NULL)
  }
  found <- NULL
  for (j in seq(i + 1L, length(info))) {
    if (!any(reach %in% info[[j]]$symbols)) {
      next
    }
    if (info[[j]]$inline) {
      found <- c(found, character(0))
    } else if (length(info[[j]]$refs)) {
      found <- unique(c(found, info[[j]]$refs))
    } else {
      reach <- c(reach, info[[j]]$lhs)
    }
  }
  found
}

# Which models the dose records of expression `i` are for, and how that was
# decided; `current` is the models referred to most recently before it
# (character(0) after an inline model, NULL before any).
.doseOwners <- function(info, i, current, models) {
  if (length(info[[i]]$refs)) {
    return(list(owner = "referenced", models = info[[i]]$refs))
  }
  used <- .modelUses(info, i)
  if (!is.null(used)) {
    return(list(owner = "used", models = used))
  }
  if (!is.null(current)) {
    return(list(owner = "recent", models = current))
  }
  list(owner = "all", models = models)
}

#' Compartments an article doses, and the models it doses them in
#'
#' Each dose record is credited to the models its top-level expression refers
#' to ("referenced"); or, when that expression assigns a variable, to the
#' models referred to by the later expressions that read the variable
#' ("used"; none, when those solve a model the article defines itself); or
#' else to the models the article referred to most recently before
#' it ("recent"; none, after a model the article defines itself); or, before
#' any, to every model the article covers ("all").
#'
#' @param rmd Path to the article (.Rmd).
#' @param models The registry models the article covers.
#' @return A data frame with one row per credited compartment and dose
#'   record: `model`, `cmt` (a numeric `cmt` is "#<n>"), `via` ("literal" or
#'   "variable"), `owner` ("referenced", "used", "recent" or "all") and
#'   `record`, the record's code.
#' @noRd
articleDoseCmts <- function(rmd, models) {
  exprs <- .articleExprs(rmd)
  vars <- list()
  for (e in exprs) {
    vars <- .stringVarsIn(e, vars)
  }
  info <- .articleModelRefs(exprs, models)
  doseFns <- .doseFunctions(exprs, vars)
  topColumns <- .assembleColumns(exprs)
  current <- NULL
  rows <- list()
  for (i in seq_along(exprs)) {
    rec <- .doseRecordsIn(exprs[[i]], vars)
    for (b in .nestedColumnRecords(exprs[[i]])) {
      rec <- .doseRecordsIn(b, vars, rec)
    }
    for (tc in topColumns) {
      if (tc$at == i) rec <- .doseRecordsIn(tc$record, vars, rec)
    }
    rec <- .doseFunctionCallsIn(exprs[[i]], doseFns, vars, rec)
    if (length(rec)) {
      rec <- unique(do.call(rbind, rec))
      own <- .doseOwners(info, i, current, models)
      idx <- expand.grid(m = seq_along(own$models), r = seq_len(nrow(rec)))
      if (nrow(idx)) {
        rows[[length(rows) + 1L]] <- data.frame(
          model = own$models[idx$m],
          rec[idx$r, , drop = FALSE],
          owner = own$owner,
          row.names = NULL
        )
      }
    }
    if (length(info[[i]]$refs)) {
      current <- info[[i]]$refs
    }
    if (info[[i]]$defines) current <- character(0)
  }
  if (!length(rows)) {
    return(data.frame(
      model = character(0),
      cmt = character(0),
      via = character(0),
      record = character(0),
      owner = character(0)
    ))
  }
  unique(do.call(rbind, rows))
}

# The ODE states a model file declares with `d/dt()`, in file order.
.declaredStates <- function(file) {
  mb <- .probeModelBlockOfFile(file)
  st <- character(0)
  for (i in seq_along(mb)[-1]) {
    x <- mb[[i]]
    if (is.call(x) && length(x) == 3L && .probeIsDdt(x[[2]])) st <- c(st, as.character(x[[2]][[3]][[2]]))
  }
  unique(st)
}

# The model state each credited dose names, NA where it names none: a name
# must be one of the model's `d/dt()` states, and "#<n>" is the model's n-th
# state in rxode2's order.
.articleDoseStates <- function(d, files) {
  out <- rep(NA_character_, nrow(d))
  for (m in unique(d$model)) {
    i <- d$model == m
    states <- .declaredStates(files[basename(files) == paste0(m, ".R")])
    out[i & d$cmt %in% states] <- d$cmt[i & d$cmt %in% states]
    num <- i & !is.na(d$cmt) & startsWith(d$cmt, "#")
    if (any(num)) {
      order <- suppressWarnings(suppressMessages(rxode2::rxModelVars(readModelDb(m))$state))
      k <- as.integer(substring(d$cmt[num], 2L))
      # A number past the last state names none.
      out[num] <- order[match(k, seq_along(order))]
    }
  }
  out
}

# Is `state` one of the comma-separated `dosing` entries?
.inDosing <- function(state, dosing) {
  !is.na(dosing) && state %in% strsplit(dosing, ",", fixed = TRUE)[[1]]
}

#' Article doses missing from their models' `dosing`
#'
#' Reads every article in `articleDir` that covers a registry model (see
#' articleDoseCmts()) and returns each model state it doses that the model's
#' `dosing` lacks.
#'
#' @param articleDir The articles directory.
#' @param db The model registry.
#' @param exempt A data frame with columns `model` and `state`: doses that
#'   set or displace an endogenous state rather than give something.
#' @return A list: `gaps`, a data frame with one row per missing model and
#'   state (`article`, `model`, `state`, and one `record` that doses it), and
#'   `exempted`, the `exempt` rows the articles do dose.
#' @noRd
articleDoseGaps <- function(articleDir, db, exempt) {
  files <- list.files(system.file("modeldb", package = "nlmixr2lib"), "[.]R$", recursive = TRUE, full.names = TRUE)
  article <- ifelse(is.na(db$vignette), db$name, db$vignette)
  rmds <- file.path(articleDir, paste0(unique(article), ".Rmd"))
  rows <- list()
  for (rmd in rmds[file.exists(rmds)]) {
    a <- sub("[.]Rmd$", "", basename(rmd))
    d <- articleDoseCmts(rmd, db$name[article == a])
    if (nrow(d)) {
      d$state <- .articleDoseStates(d, files)
    }
    d <- d[!is.na(d$state), ]
    missing <- !mapply(.inDosing, d$state, db$dosing[match(d$model, db$name)], USE.NAMES = FALSE)
    if (any(missing)) rows[[a]] <- data.frame(article = a, d[missing, c("model", "state", "record")])
  }
  all <- do.call(
    rbind,
    c(list(data.frame(article = character(0), model = character(0), state = character(0), record = character(0))), rows)
  )
  all <- all[!duplicated(all[, c("model", "state")]), ]
  isExempt <- paste(all$model, all$state) %in% paste(exempt$model, exempt$state)
  list(gaps = all[!isExempt, ], exempted = exempt[paste(exempt$model, exempt$state) %in% paste(all$model, all$state), ])
}
