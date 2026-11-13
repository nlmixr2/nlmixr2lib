# Support for the per-level random-effect draw gate (test-articleNesting.R).
#
# An `eta ~ variance | G` eta is drawn per level of G. rxode2 decides from the
# simulation data whether G sits below or above the subject:
#
#   * some subject has two or more levels of G  -> per subject and level
#     (what an occasion or a catheter needs);
#   * every subject has exactly one level       -> G is taken as a level ABOVE
#     the subject, and every subject on a level shares one draw;
#   * the whole table has a single level        -> rxSolve() fails with
#     "parameter(s) are required for solving: THETA[1]".
#
# For an occasion-type G the shared draw removes the IOV spread without any
# error, so an article that simulates such a model has to assert the number of
# distinct draws per level. This file finds, statically, the articles that
# simulate one of those models without doing so.

# Grouping columns that name a level ABOVE the subject on purpose (a study),
# where one draw per level is the model. Each entry is (model, group) with the
# evidence for classifying it that way; anything not listed is below-subject.
.nestAboveSubjectGroups <- data.frame(
  model = c("Qi_2024_vosoritide", "Rosenborg_2025_salmeterol"),
  group = c("SIDN", "SIDN"),
  evidence = c(
    paste(
      "Qi 2024 Sect. 3.2: SIDN is the study identity number; StudyCL and StudyV",
      "are an additional hierarchical level of effect above the subject"
    ),
    "Rosenborg 2025: SIDN maps to the STUD column, studies 1 to 3, one study per subject"
  ),
  stringsAsFactors = FALSE
)

# Depth-first walk over a language object, calling `fn` on every node. A `fn`
# that returns FALSE stops the walk below that node. Empty arguments (the
# `x[, 1]` form) are skipped.
.nestWalk <- function(expr, fn) {
  if (is.expression(expr)) {
    for (i in seq_along(expr)) {
      .nestWalk(expr[[i]], fn)
    }
    return(invisible(NULL))
  }
  if (identical(fn(expr), FALSE)) {
    return(invisible(NULL))
  }
  if (is.call(expr)) {
    kids <- as.list(expr)
    for (i in seq_along(kids)) {
      # An empty argument (`x[, 1]`) is an unbound symbol that cannot be
      # passed on, so it is compared inside its list slot and skipped.
      if (!identical(kids[i], list(quote(expr = )))) {
        .nestWalk(kids[[i]], fn)
      }
    }
  }
  invisible(NULL)
}

# The name a call is made through, with any `pkg::` prefix dropped, or NA.
.nestCallName <- function(e) {
  if (!is.call(e)) {
    return(NA_character_)
  }
  head <- e[[1L]]
  if (is.name(head)) {
    return(as.character(head))
  }
  if (is.call(head) && length(head) == 3L && is.name(head[[1L]]) && as.character(head[[1L]]) %in% c("::", ":::")) {
    return(as.character(head[[3L]]))
  }
  NA_character_
}

# The `eta ~ variance | G` etas declared in a model file's ini() block, as one
# row per (eta, G). Parsing rather than grepping means a comment or a string
# cannot produce a row, and the `y ~ add(a) | cmt` endpoint form in model() is
# not mistaken for one. A file without both a `~` and a `|` anywhere in its text
# cannot declare such an eta however its lines are split, so only those are
# parsed. A group that is not a plain column name (`| interaction(a, b)`) is an
# error: the gate cannot say which column it would have to be asserted on.
.nestedEtas <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  empty <- data.frame(eta = character(), group = character(), stringsAsFactors = FALSE)
  if (!grepl("~", txt, fixed = TRUE) || !grepl("|", txt, fixed = TRUE)) {
    return(empty)
  }
  exprs <- parse(text = txt, keep.source = FALSE)
  acc <- new.env(parent = emptyenv())
  acc$rows <- list()
  .nestWalk(exprs, function(e) {
    nm <- .nestCallName(e)
    if (identical(nm, "model")) {
      return(FALSE)
    }
    if (identical(nm, "ini")) {
      .nestWalk(e, function(s) {
        if (!is.call(s) || !identical(s[[1L]], as.name("~")) || length(s) != 3L) {
          return(invisible(NULL))
        }
        rhs <- s[[3L]]
        if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
          if (!is.name(rhs[[3L]])) {
            stop(errorCondition(
              sprintf(
                "%s: the grouping of `%s` is `%s`, not a plain column name; extend .nestedEtas() to resolve it",
                basename(path),
                paste(deparse(s[[2L]]), collapse = " "),
                paste(deparse(rhs[[3L]]), collapse = " ")
              ),
              class = "nestUnsupportedGroup"
            ))
          }
          acc$rows[[length(acc$rows) + 1L]] <- data.frame(
            eta = all.vars(s[[2L]]),
            group = as.character(rhs[[3L]]),
            stringsAsFactors = FALSE
          )
        }
        invisible(NULL)
      })
      return(FALSE)
    }
    invisible(NULL)
  })
  if (length(acc$rows) == 0L) {
    return(empty)
  }
  unique(do.call(rbind, acc$rows))
}

# Every model under `model_root` that has a `| G` eta, as (model, eta, group,
# above_subject). The registry is not consulted: it keeps no ini() structure.
.nestedModels <- function(model_root, above = .nestAboveSubjectGroups) {
  files <- list.files(model_root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  rows <- lapply(files, function(f) {
    e <- .nestedEtas(f)
    if (nrow(e) == 0L) {
      return(NULL)
    }
    e$model <- sub("[.]R$", "", basename(f))
    e
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    return(data.frame(
      model = character(),
      eta = character(),
      group = character(),
      above_subject = logical(),
      stringsAsFactors = FALSE
    ))
  }
  key <- paste(rows$model, rows$group)
  rows$above_subject <- key %in% paste(above$model, above$group)
  rows[, c("model", "eta", "group", "above_subject")]
}

# The R chunks of an Rmd, parsed, as list(chunks, unparseable). A chunk is
# skipped when it does not run (`eval = FALSE` in the header, or `#| eval: false`
# / `#| eval = FALSE` in its body). `unparseable` holds the line of each running
# chunk that does not parse, so the caller cannot mistake it for an empty chunk.
.nestArticleExprs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  is_open <- grepl("^\\s*```\\{[rR][ ,}]", lines)
  is_close <- grepl("^\\s*```\\s*$", lines)
  chunks <- list()
  unparseable <- integer()
  inchunk <- FALSE
  skip <- FALSE
  start <- 0L
  buf <- character()
  for (i in seq_along(lines)) {
    if (!inchunk && is_open[[i]]) {
      inchunk <- TRUE
      skip <- grepl("eval\\s*=\\s*(FALSE|F)\\b", lines[[i]])
      start <- i
      buf <- character()
    } else if (inchunk && is_close[[i]]) {
      inchunk <- FALSE
      if (!skip) {
        parsed <- tryCatch(parse(text = buf, keep.source = FALSE), error = function(e) NULL)
        if (is.null(parsed)) {
          unparseable <- c(unparseable, start)
        } else {
          chunks[[length(chunks) + 1L]] <- parsed
        }
      }
    } else if (inchunk) {
      if (grepl("^#\\|\\s*eval\\s*[:=]\\s*(false|FALSE|False|F)\\s*$", lines[[i]])) {
        skip <- TRUE
      }
      buf <- c(buf, lines[[i]])
    }
  }
  list(chunks = chunks, unparseable = unparseable)
}

# Every symbol and string constant in a language object, or only the string
# constants.
.nestTokens <- function(expr, strings_only = FALSE) {
  acc <- new.env(parent = emptyenv())
  acc$tok <- character()
  .nestWalk(expr, function(e) {
    if (is.name(e) && !strings_only) {
      acc$tok <- c(acc$tok, as.character(e))
    } else if (is.character(e)) {
      acc$tok <- c(acc$tok, e[!is.na(e)])
    }
    invisible(NULL)
  })
  unique(acc$tok)
}

# The tokens of every element of a list of language objects; empty arguments
# are skipped.
.nestTokensList <- function(lst, strings_only = FALSE) {
  out <- character()
  for (i in seq_along(lst)) {
    if (!identical(lst[i], list(quote(expr = )))) {
      out <- c(out, .nestTokens(lst[[i]], strings_only))
    }
  }
  unique(out)
}

# Follow the symbols in `tokens` back through the assignments that define them,
# up to four levels, and return every token met on the way. `assign` maps a
# symbol to the tokens of every right-hand side assigned to it.
.nestExpandTokens <- function(tokens, assign) {
  seen <- character()
  todo <- tokens
  for (depth in seq_len(4L)) {
    nxt <- intersect(setdiff(todo, seen), names(assign))
    seen <- union(seen, todo)
    todo <- unique(unlist(assign[nxt], use.names = FALSE))
    if (length(todo) == 0L) {
      break
    }
    tokens <- union(tokens, todo)
  }
  tokens
}

# Functions that count distinct values.
.nestCountWords <- c("n_distinct", "unique", "uniqueN", "duplicated", "anyDuplicated")

# Calls whose first argument (the data, the columns, the vector iterated over)
# or whose `.by` / `by` argument names what a count inside them is taken over,
# and calls whose every argument does.
.nestScopeFirstArg <- c(
  "summarise",
  "summarize",
  "mutate",
  "reframe",
  "transmute",
  "across",
  "with",
  "within",
  "lapply",
  "sapply",
  "vapply",
  "map",
  "map_dbl",
  "map_int",
  "map_lgl",
  "map_chr",
  "map_df",
  "%>%",
  "%<>%",
  "%T>%"
)
.nestScopeAllArgs <- c("tapply", "aggregate", "by", "ave", "mapply", "Map", "split", "outer")

# One token vector per count call (n_distinct, unique, ...) in `expr`: the
# tokens of the call's own arguments, plus those of the data and grouping
# arguments of the calls around it. `anc` is the chain of enclosing calls.
.nestScopeWalk <- function(expr, anc, acc) {
  if (is.expression(expr)) {
    for (i in seq_along(expr)) {
      .nestScopeWalk(expr[[i]], anc, acc)
    }
    return(invisible(NULL))
  }
  if (!is.call(expr)) {
    return(invisible(NULL))
  }
  nm <- .nestCallName(expr)
  if (!is.na(nm) && nm %in% .nestCountWords) {
    scope <- .nestTokensList(as.list(expr)[-1L])
    for (a in anc) {
      an <- .nestCallName(a)
      if (is.na(an)) {
        next
      }
      args <- as.list(a)[-1L]
      if (an %in% .nestScopeFirstArg) {
        by_args <- if (is.null(names(args))) list() else args[names(args) %in% c(".by", "by")]
        scope <- c(scope, .nestTokensList(args[seq_len(min(1L, length(args)))]), .nestTokensList(by_args))
      }
      if (an %in% .nestScopeAllArgs) {
        scope <- c(scope, .nestTokensList(args))
      }
    }
    acc$scopes[[length(acc$scopes) + 1L]] <- unique(scope)
  }
  kids <- as.list(expr)
  for (i in seq_along(kids)) {
    if (!identical(kids[i], list(quote(expr = )))) {
      .nestScopeWalk(kids[[i]], c(anc, list(expr)), acc)
    }
  }
  invisible(NULL)
}

# The count-call scopes of a list of expressions.
.nestCountScopes <- function(exprs) {
  acc <- new.env(parent = emptyenv())
  acc$scopes <- list()
  for (x in exprs) {
    .nestScopeWalk(x, list(), acc)
  }
  acc$scopes
}

# Functions that join pieces into a string; the literal pieces of a model name
# built with one are fragments of the name.
.nestConcatenators <- c("paste", "paste0", "sprintf", "glue", "file.path", "str_c", "str_glue")

# Calls through which the value of a model name passes unchanged: the vector or
# table is built from their arguments (all of them, or the first), and a
# mutate()-like call adds columns to its first argument.
.nestPassAll <- c("c", "list", "tibble", "data.frame", "data_frame", "tribble", "expand.grid", "bind_rows")
.nestPassFirst <- c(
  "setNames",
  "unique",
  "rev",
  "sort",
  "rep",
  "unlist",
  "as.character",
  "as.vector",
  "head",
  "tail",
  "filter",
  "arrange",
  "select",
  "slice",
  "distinct",
  "rename",
  "left_join",
  "inner_join",
  "pull",
  "unname",
  "return",
  "(",
  "["
)
.nestPassBranches <- c("ifelse", "if_else", "if", "switch")
.nestPassMutate <- c("mutate", "transmute")

# The loader functions that take a model name.
.nestLoaders <- c("readModelDb", "modellib")

# Calls that apply a function over a vector; the loader may be passed as the
# function, and a function literal's parameters take the vector's values.
.nestIterators <- c(
  "lapply",
  "sapply",
  "vapply",
  "map",
  "map_dbl",
  "map_int",
  "map_lgl",
  "map_chr",
  "map_df",
  "Map",
  "mapply",
  "walk",
  "Filter"
)

# The model-name argument of a readModelDb() / modellib() call: the argument named
# `name`, else the first unnamed one; NULL when the call passes none.
.nestLoaderArg <- function(e) {
  args <- as.list(e)[-1L]
  if (length(args) == 0L) {
    return(NULL)
  }
  nms <- names(args)
  if (is.null(nms)) {
    return(args[[1L]])
  }
  if ("name" %in% nms) {
    return(args[["name"]])
  }
  unnamed <- which(!nzchar(nms))
  if (length(unnamed) == 0L) NULL else args[[unnamed[[1L]]]]
}

# Record that `key` takes the value of `expr`: its tokens and the expression
# itself.
.nestBind <- function(acc, key, expr) {
  acc$assign[[key]] <- c(acc$assign[[key]], .nestTokens(expr))
  acc$assign_expr[[key]] <- c(acc$assign_expr[[key]], list(expr))
  invisible(NULL)
}

# The symbol an assignment target writes to: `x`, `x$a`, `x[i]`, `x[[i]]`.
.nestAssignTarget <- function(lhs) {
  while (is.call(lhs) && length(lhs) >= 2L && .nestCallName(lhs) %in% c("$", "[", "[[", "@")) {
    lhs <- lhs[[2L]]
  }
  if (is.name(lhs)) as.character(lhs) else NA_character_
}

# Whether `x` is the loader itself passed as a function (`readModelDb`,
# `nlmixr2lib::readModelDb`).
.nestIsLoaderRef <- function(x) {
  if (is.name(x)) {
    return(as.character(x) %in% .nestLoaders)
  }
  is.call(x) && .nestCallName(x) %in% c("::", ":::") && as.character(x[[3L]]) %in% .nestLoaders
}

# Bind what a call node defines: an assignment (and the parameters of a function
# it assigns), a loop variable, and the parameters of a user-defined function
# called here.
.nestCollectBindings <- function(e, nm, args, acc) {
  if (nm %in% c("<-", "=") && length(e) == 3L) {
    key <- .nestAssignTarget(e[[2L]])
    if (!is.na(key)) {
      .nestBind(acc, key, e[[3L]])
      rhs <- e[[3L]]
      if (is.call(rhs) && identical(.nestCallName(rhs), "function")) {
        acc$formals[[key]] <- names(rhs[[2L]])
      }
    }
  }
  if (nm == "for" && length(e) >= 3L && is.name(e[[2L]])) {
    .nestBind(acc, as.character(e[[2L]]), e[[3L]])
  }
  if (nm %in% names(acc$formals)) {
    .nestBindCall(acc, acc$formals[[nm]], args)
  }
  invisible(NULL)
}

# An iterator call (lapply, map, ...): the parameters of a function literal, the
# `.x` of a formula and the parameters of a user-defined function given to it
# take the values of the vector iterated over, and a loader given to it is a
# loader call on that vector.
.nestCollectIterator <- function(args, acc) {
  is_fn <- vapply(args, function(a) .nestIsLoaderRef(a) || identical(.nestCallName(a), "function"), logical(1L))
  data_args <- args[!is_fn]
  if (!is.null(names(data_args))) {
    data_args <- data_args[
      !(names(data_args) %in% c("FUN", "f", ".f", "FUN.VALUE", ".progress", "SIMPLIFY", "USE.NAMES"))
    ]
  }
  for (a in args) {
    if (identical(.nestCallName(a), "function")) {
      for (p in names(a[[2L]])) {
        for (d in data_args) {
          .nestBind(acc, p, d)
        }
      }
    }
    if (identical(.nestCallName(a), "~")) {
      for (d in data_args) {
        for (p in c(".x", ".", "..1")) {
          .nestBind(acc, p, d)
        }
      }
    }
    if (.nestIsLoaderRef(a) && length(data_args) > 0L) {
      acc$loader_args[[length(acc$loader_args) + 1L]] <- data_args[[1L]]
    }
    if (is.name(a) && as.character(a) %in% names(acc$formals) && length(data_args) > 0L) {
      .nestBindCall(acc, acc$formals[[as.character(a)]], data_args)
    }
  }
  invisible(NULL)
}

# Collect from one call node: bindings, assertions, loader arguments and
# simulation calls.
.nestCollect <- function(e, acc) {
  nm <- .nestCallName(e)
  if (is.na(nm)) {
    return(invisible(NULL))
  }
  args <- as.list(e)[-1L]
  .nestCollectBindings(e, nm, args, acc)
  if (nm %in% .nestIterators) {
    .nestCollectIterator(args, acc)
  }
  if (nm %in% c("stopifnot", "assert_that")) {
    acc$asserts <- c(acc$asserts, args)
  } else if (grepl("^expect_", nm)) {
    acc$asserts[[length(acc$asserts) + 1L]] <- e
  }
  if (nm %in% .nestLoaders) {
    arg <- .nestLoaderArg(e)
    if (!is.null(arg) && !is.character(arg)) {
      acc$loader_args[[length(acc$loader_args) + 1L]] <- arg
    }
  }
  if (nm %in% c("rxSolve", "simulate", "solve") && length(e) >= 2L) {
    if (!identical(.nestCallName(e[[2L]]), "zeroRe")) {
      acc$simulates <- TRUE
    }
  }
  invisible(NULL)
}

# Bind the parameters of a user-defined function to the arguments of one call
# to it, by name and then by position.
.nestBindCall <- function(acc, formals, args) {
  an <- names(args)
  if (is.null(an)) {
    an <- rep("", length(args))
  }
  free <- setdiff(formals, an)
  for (i in seq_along(args)) {
    if (identical(args[i], list(quote(expr = )))) {
      next
    }
    target <- if (nzchar(an[[i]]) && an[[i]] %in% formals) {
      an[[i]]
    } else if (length(free) > 0L) {
      free[[1L]]
    } else {
      NA_character_
    }
    if (!is.na(target)) {
      free <- setdiff(free, target)
      .nestBind(acc, target, args[[i]])
    }
  }
  invisible(NULL)
}

# A regular expression for each literal in `x`.
.nestEscape <- function(x) {
  gsub("([^A-Za-z0-9_ ])", "\\\\\\1", x, perl = TRUE)
}

# The pattern a concatenator call can produce: its literal pieces, escaped, with
# `.*` for every piece that is not a string constant.
.nestConcatPattern <- function(x) {
  nm <- .nestCallName(x)
  args <- as.list(x)[-1L]
  an <- names(args)
  if (is.null(an)) {
    an <- rep("", length(args))
  }
  if (nm %in% c("sprintf", "glue", "str_glue")) {
    fmt <- if (length(args) > 0L) args[[1L]] else NULL
    if (!is.character(fmt) || length(fmt) != 1L) {
      return(".*")
    }
    marked <- gsub("%[-+ 0#]*[0-9*]*(?:[.][0-9*]+)?[a-zA-Z]|\\{[^}]*\\}", "\001", fmt, perl = TRUE)
    pieces <- strsplit(marked, "\001", fixed = TRUE)[[1L]]
    if (endsWith(marked, "\001")) {
      pieces <- c(pieces, "")
    }
    return(paste(.nestEscape(pieces), collapse = ".*"))
  }
  sep <- if (nm %in% c("paste0", "str_c")) {
    ""
  } else if (nm == "file.path") {
    "/"
  } else {
    " "
  }
  if ("sep" %in% an && is.character(args[[which(an == "sep")[[1L]]]])) {
    sep <- args[[which(an == "sep")[[1L]]]]
  }
  keep <- !(an %in% c("sep", "collapse"))
  pieces <- vapply(
    args[keep],
    function(a) if (is.character(a) && length(a) == 1L) .nestEscape(a) else ".*",
    character(1L)
  )
  paste(pieces, collapse = .nestEscape(sep))
}

.nestMergeNames <- function(parts) {
  list(
    pat = unique(unlist(lapply(parts, `[[`, "pat"), use.names = FALSE)),
    unknown = any(vapply(parts, `[[`, logical(1L), "unknown")),
    why = unique(unlist(lapply(parts, `[[`, "why"), use.names = FALSE))
  )
}

.nestUnknownName <- function(why) {
  list(pat = character(), unknown = TRUE, why = why)
}

# What a model-name expression can evaluate to, as list(pat, unknown): the
# patterns its values match (a literal gives its own text, a built name a
# pattern with `.*`) and whether any value comes from somewhere the gate cannot
# see (a symbol with no value bound in the article, a package object, a string
# transformer, a function it does not follow). Symbols are followed through the
# values bound to them; `x$col` and `x[["col"]]` select the named argument of a
# table built in the article.
.nestNameValue <- function(x, acc, depth = 0L, seen = character(), col = NA_character_) {
  if (depth > 10L) {
    return(.nestUnknownName("too deep"))
  }
  if (is.character(x)) {
    return(list(pat = .nestEscape(x[!is.na(x)]), unknown = FALSE, why = character()))
  }
  if (is.numeric(x) || is.logical(x) || is.null(x)) {
    return(list(pat = character(), unknown = FALSE, why = character()))
  }
  if (is.name(x)) {
    return(.nestNameValueSymbol(as.character(x), acc, depth, seen, col))
  }
  if (!is.call(x)) {
    return(.nestUnknownName("not an expression"))
  }
  nm <- .nestCallName(x)
  if (is.na(nm)) {
    return(.nestUnknownName("computed function"))
  }
  .nestNameValueCall(x, nm, as.list(x)[-1L], acc, depth, seen, col)
}

# The values of each expression in `exprs`, merged.
.nestNameValueAll <- function(exprs, acc, depth, seen, col) {
  .nestMergeNames(lapply(exprs, .nestNameValue, acc = acc, depth = depth + 1L, seen = seen, col = col))
}

# The value of a symbol: the values bound to it. A symbol reused for unrelated
# values is resolved by the values that resolve.
.nestNameValueSymbol <- function(key, acc, depth, seen, col) {
  if (key %in% seen) {
    return(list(pat = character(), unknown = FALSE, why = character()))
  }
  rhs <- acc$assign_expr[[key]]
  if (length(rhs) == 0L) {
    return(.nestUnknownName(paste0("`", key, "` has no value bound in the article")))
  }
  parts <- lapply(rhs, .nestNameValue, acc = acc, depth = depth + 1L, seen = c(seen, key), col = col)
  resolved <- parts[!vapply(parts, `[[`, logical(1L), "unknown")]
  .nestMergeNames(if (length(resolved) > 0L) resolved else parts)
}

# The values of a tribble() column: the data cells at that column's positions.
.nestNameValueTribble <- function(args, col, acc, depth, seen) {
  is_head <- vapply(args, function(a) identical(.nestCallName(a), "~"), logical(1L))
  heads <- vapply(args[is_head], function(a) paste(deparse(a[[2L]]), collapse = ""), character(1L))
  if (!(col %in% heads)) {
    return(NULL)
  }
  data <- args[!is_head]
  .nestNameValueAll(
    data[seq(which(heads == col)[[1L]], length(data), by = length(heads))],
    acc,
    depth,
    seen,
    NA_character_
  )
}

# The value of a call: a literal-building call gives the patterns of its pieces,
# a call that passes its argument through gives that argument's value, and any
# other call is unseen. `col` names the column being selected from a table.
.nestNameValueCall <- function(x, nm, args, acc, depth, seen, col) {
  an <- names(args)
  if (is.null(an)) {
    an <- rep("", length(args))
  }
  named <- !is.na(col) && col %in% an
  if (nm == "~") {
    return(list(pat = character(), unknown = FALSE, why = character()))
  }
  if (nm %in% .nestConcatenators) {
    return(list(pat = .nestConcatPattern(x), unknown = FALSE, why = character()))
  }
  if (nm %in% c("$", "[[") && length(args) == 2L) {
    pick <- if (nm == "$") {
      as.character(args[[2L]])
    } else if (is.character(args[[2L]])) {
      args[[2L]]
    } else {
      NA_character_
    }
    return(.nestNameValue(args[[1L]], acc, depth + 1L, seen, pick))
  }
  if (nm == "%>%" && length(args) == 2L && is.call(args[[2L]])) {
    piped <- as.call(c(args[[2L]][[1L]], list(args[[1L]]), as.list(args[[2L]])[-1L]))
    return(.nestNameValue(piped, acc, depth + 1L, seen, col))
  }
  if (nm == "{" && length(args) >= 1L) {
    return(.nestNameValue(args[[length(args)]], acc, depth + 1L, seen, col))
  }
  if (nm %in% .nestPassBranches && length(args) >= 2L) {
    return(.nestNameValueAll(args[-1L], acc, depth, seen, col))
  }
  if (nm == "tribble" && !is.na(col)) {
    value <- .nestNameValueTribble(args, col, acc, depth, seen)
    if (!is.null(value)) {
      return(value)
    }
  }
  if (nm %in% .nestPassAll) {
    return(.nestNameValueAll(if (named) args[an == col] else args, acc, depth, seen, NA_character_))
  }
  if (nm %in% .nestPassMutate && length(args) >= 1L) {
    return(.nestNameValueAll(
      if (named) args[an == col] else args[1L],
      acc,
      depth,
      seen,
      if (named) NA_character_ else col
    ))
  }
  if (nm %in% .nestPassFirst && length(args) >= 1L) {
    return(.nestNameValue(args[[1L]], acc, depth + 1L, seen, col))
  }
  if (nm %in% names(acc$formals)) {
    defs <- Filter(function(r) identical(.nestCallName(r), "function"), acc$assign_expr[[nm]])
    if (length(defs) > 0L) {
      return(.nestMergeNames(lapply(defs, function(f) .nestNameValue(f[[3L]], acc, depth + 1L, seen, col))))
    }
  }
  .nestUnknownName(paste0("value of `", nm, "()`"))
}

# Why a model name given to a loader cannot be matched against `models`, as a
# character vector that is empty when it can: some value behind it is unseen
# (see .nestNameValue()), or a name built from pieces can be one of `models`
# (`paste0("Smith_2020_", drug)` cannot when no model starts with `Smith_2020_`).
# Literal names are matched elsewhere.
.nestDynamicName <- function(arg, acc, models) {
  value <- .nestNameValue(arg, acc)
  if (value$unknown) {
    return(value$why)
  }
  built <- value$pat[grepl(".*", value$pat, fixed = TRUE)]
  hit <- unique(unlist(lapply(built, function(p) models[grepl(paste0("^", p, "$"), models, perl = TRUE)])))
  if (length(hit) > 0L) {
    return(paste0("a name built from pieces can be ", paste(hit, collapse = ", ")))
  }
  character()
}

.nestNoFindings <- function() {
  list(
    models = character(),
    etas = character(),
    simulates = FALSE,
    uncovered = character(),
    dynamic = character(),
    unparseable = integer()
  )
}

# The article's findings:
#   models, etas  the below-subject models it names and their `| G` etas;
#   simulates     it calls rxSolve / simulate / solve with a first argument that
#                 is not zeroRe();
#   uncovered     etas with no per-level distinct-draw assertion. A count call
#                 (n_distinct, unique, ...) covers an eta when its own arguments,
#                 together with the data and grouping arguments of the calls
#                 around it, hold the eta and its grouping column once their
#                 symbols are followed back through assignments. The count is
#                 looked for in each stopifnot() argument or expect_*() call and
#                 in the assignments those reach, so a summary computed earlier
#                 and asserted later counts, while a `unique()` that only builds
#                 the event table does not;
#   dynamic       the reasons a readModelDb() / modellib() call, or a loader
#                 passed to lapply() and the like, has a model name that cannot
#                 be matched (see .nestDynamicName()); empty when none;
#   unparseable   lines of running R chunks that do not parse.
# Most articles neither name a below-subject model nor pass a non-literal name
# to a loader; those are not parsed.
.nestArticleFindings <- function(path, nested) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  below <- nested[!nested$above_subject, , drop = FALSE]
  names_model <- any(vapply(unique(below$model), grepl, logical(1L), x = text, fixed = TRUE))
  non_literal_loader <- grepl("(?:readModelDb|modellib)\\s*(?:\\(\\s*[^\"'\\s)]|[,)])", text, perl = TRUE)
  if (!names_model && !non_literal_loader) {
    return(.nestNoFindings())
  }
  parsed <- .nestArticleExprs(path)
  acc <- new.env(parent = emptyenv())
  acc$strings <- character()
  acc$assign <- list()
  acc$assign_expr <- list()
  acc$formals <- list()
  acc$asserts <- list()
  acc$loader_args <- list()
  acc$simulates <- FALSE
  for (ch in parsed$chunks) {
    .nestWalk(ch, function(e) {
      if (is.character(e)) {
        acc$strings <- c(acc$strings, e)
      } else if (is.name(e)) {
        acc$strings <- c(acc$strings, as.character(e))
      }
      .nestCollect(e, acc)
    })
  }
  acc$assign <- lapply(acc$assign, unique)
  used <- nested[nested$model %in% acc$strings & !nested$above_subject, , drop = FALSE]
  scopes <- unlist(
    lapply(acc$asserts, function(a) {
      tok <- .nestExpandTokens(.nestTokens(a), acc$assign)
      reached <- unlist(acc$assign_expr[intersect(tok, names(acc$assign_expr))], recursive = FALSE)
      .nestCountScopes(c(list(a), reached))
    }),
    recursive = FALSE
  )
  scope_tokens <- lapply(scopes, .nestExpandTokens, assign = acc$assign)
  covered <- vapply(
    seq_len(nrow(used)),
    function(i) {
      any(vapply(scope_tokens, function(tok) all(c(used$eta[[i]], used$group[[i]]) %in% tok), logical(1L)))
    },
    logical(1L)
  )
  dynamic <- unique(unlist(lapply(acc$loader_args, .nestDynamicName, acc = acc, models = unique(below$model))))
  list(
    models = unique(used$model),
    etas = unique(used$eta),
    simulates = acc$simulates,
    uncovered = unique(used$eta[!covered]),
    dynamic = dynamic,
    unparseable = parsed$unparseable
  )
}

.nestViolationRow <- function(path, models, etas, reason) {
  data.frame(
    article = basename(path),
    models = models,
    etas = etas,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

# One row per (article, reason): an article that simulates a below-subject
# `| G` model without a per-level distinct-draw assertion on every one of that
# model's etas, one that simulates a model it cannot identify, and one with a
# running R chunk that does not parse.
.nestDrawViolations <- function(model_root, article_paths, above = .nestAboveSubjectGroups) {
  nested <- .nestedModels(model_root, above)
  rows <- lapply(article_paths, function(p) {
    f <- .nestArticleFindings(p, nested)
    out <- list()
    if (length(f$models) > 0L && f$simulates && length(f$uncovered) > 0L) {
      out[[1L]] <- .nestViolationRow(
        p,
        paste(f$models, collapse = ", "),
        paste(f$uncovered, collapse = ", "),
        "no per-level distinct-draw assertion"
      )
    }
    if (f$simulates && length(f$dynamic) > 0L) {
      out[[length(out) + 1L]] <- .nestViolationRow(
        p,
        "",
        "",
        paste0(
          "unverifiable: a model name given to readModelDb() or modellib() cannot be resolved (",
          paste(utils::head(f$dynamic, 3L), collapse = "; "),
          ")"
        )
      )
    }
    if (length(f$unparseable) > 0L) {
      out[[length(out) + 1L]] <- .nestViolationRow(
        p,
        "",
        "",
        paste0("unverifiable: R chunk at line ", paste(f$unparseable, collapse = ", "), " does not parse")
      )
    }
    do.call(rbind, out)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    return(data.frame(
      article = character(),
      models = character(),
      etas = character(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows
}
