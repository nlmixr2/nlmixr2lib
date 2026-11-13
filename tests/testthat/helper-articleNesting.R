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

# Known limits. The gate is a static reading of an article's R chunks; it is
# written to flag or report an article it cannot read, never to pass one silently,
# but it is a heuristic:
#   * Code that is not in a running chunk is invisible: source(), eval(parse()),
#     child documents, knit parameters and chunks with eval = FALSE.
#   * A count is recognised as n_distinct / unique / uniqueN / duplicated /
#     anyDuplicated, called or passed by name to an iterator. table(),
#     nrow(distinct()), lengths(), dplyr::count() and tally() are not, so an
#     article asserting with them is flagged.
#   * A count restricted to one level of the group (`x[group == 1]`) covers the
#     eta on purpose: when each subject has one level the draws are shared and
#     one level already shows a single distinct value, so the gate does not
#     require every level.
#   * Assignments are followed up to four levels, by symbol name. A parameter is
#     looked up in its own function and an assignment inside a function is
#     visible only inside it; a symbol reused at the top level is not separated.
#   * Model names are evaluated by an interpreter of a short list of calls (c,
#     list, tibble, paste*, sprintf, `$`, ifelse, ...). Any other call is an
#     unseen value and is reported; the five articles in
#     .nestDynamicLoaderExceptions were excused by hand.
#   * A simulation is a call to rxSolve / simulate / solve / predict, a solver
#     passed to an iterator or do.call(), or nlmixr2(est = "rxSolve"); a helper
#     that wraps the solve is found through the solver call inside it. A solver
#     aliased to a symbol (`s <- rxSolve; s(mod, ev)`), stored in a list, or passed
#     to a function outside .nestApplyCalls is not found.
#   * The iterators that hand a function the values it iterates over are
#     .nestIterators (lapply, sapply, vapply, map, map_dbl, map_int, map_lgl,
#     map_chr, map_df, Map, mapply, walk, Filter); .nestApplyCalls adds do.call,
#     match.fun, get, pmap, invoke_map, Reduce and Position for finding a solver
#     or loader passed by name. tapply, aggregate, by, ave and map2 are in
#     neither, so a model name that reaches a loader through them is an unbound
#     parameter and is reported.
#   * The gate checks that a count is present, not what its predicate tests:
#     `stopifnot(is.data.frame(k))` over a per-level summary counts, and a tuple
#     count such as `n_distinct(id, OCC, etaiov)` counts although it cannot show
#     whether the draws differ between levels.
#   * A count is found wherever it is written: one inside a function the article
#     never calls, or under `if (FALSE)`, counts.
#
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

# Articles whose loader names cannot be followed by the gate but come from one
# family of models. Each row gives the name prefix of that family and the
# evidence; the article is reported again if a below-subject model starts with the
# prefix, and the test pins that each article still needs its row.
.nestDynamicLoaderExceptions <- data.frame(
  article = c(
    "Lee_2023_tylosin.Rmd",
    "Li_2025_mosunetuzumab.Rmd",
    "Shemesh_2020_atezolizumab_tmb.Rmd",
    "Sten_2023_galnacSirna.Rmd",
    "Yadav_2017_imipenem_aminoglycoside_pseudomonas.Rmd"
  ),
  prefix = c(
    "Lee_2023_tylosin_",
    "Li_2025_mosunetuzumab_",
    "Shemesh_2020_",
    "Sten_2023_",
    "Yadav_2017_"
  ),
  evidence = c(
    paste(
      "every name is one of the four Lee_2023_tylosin_* strings of the literal vector `mods`,",
      "passed through the data frames `fits` and `chk` to helper functions"
    ),
    paste(
      "names are string literals, or the `model` column of the literal tribbles `crs_fits` and",
      "`eff_check`, passed to helper functions"
    ),
    "names are the `Model` column of the literal tribbles of coefficients, passed to helper functions",
    "names are string literals and the `model` column of `tc`, assigned by a literal case_when() of Sten_2023_* names",
    paste(
      "names are the literal vector `yadav_models` and the `model` column of the literal tribble",
      "`panels`, passed through rowwise() and mutate()"
    )
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

# Calls that simulate: a model solved or predicted from.
.nestSimCalls <- c("rxSolve", "simulate", "solve", "predict")

# Calls whose first argument is the data a count inside them is taken over, so
# a grouping call in that argument groups the count.
.nestScopeData <- c(
  "summarise",
  "summarize",
  "mutate",
  "reframe",
  "transmute",
  "with",
  "within",
  "%>%",
  "%<>%",
  "%T>%"
)

# Calls whose first argument is the set of columns or elements iterated over, so
# it names the variable counted inside the function they apply.
.nestScopeColumns <- c(
  "across",
  "lapply",
  "sapply",
  "vapply",
  "map",
  "map_dbl",
  "map_int",
  "map_lgl",
  "map_chr",
  "map_df"
)

# Calls that group or split by their own arguments, and the group column of a
# count inside them.
.nestGroupCalls <- c(
  "group_by",
  "group_by_at",
  "group_by_all",
  "group_by_if",
  "split",
  "tapply",
  "aggregate",
  "ave",
  "by"
)

# Whether `expr` contains a call that simulates.
.nestContainsSolve <- function(expr) {
  found <- new.env(parent = emptyenv())
  found$hit <- FALSE
  .nestWalk(expr, function(e) {
    if (isTRUE(.nestCallName(e) %in% .nestSimCalls)) {
      found$hit <- TRUE
    }
    invisible(NULL)
  })
  found$hit
}

# The argument of `args` named `name`, else the `pos`-th unnamed one, as a list.
.nestArgNamedOrPos <- function(args, an, name, pos) {
  if (name %in% an) {
    return(args[an == name])
  }
  unnamed <- which(!nzchar(an))
  if (length(unnamed) >= pos) args[unnamed[[pos]]] else list()
}

# What a grouping or splitting call is applied to (`value`) and what it groups
# by (`group`), each as a list of expressions: group_by(data, ...) and
# split(x, f) group by every argument after the first, tapply(X, INDEX),
# by(data, INDICES) and aggregate(x, by) by their second, and aggregate(y ~ g)
# by the right of its formula.
.nestGroupSplit <- function(e) {
  nm <- .nestCallName(e)
  args <- as.list(e)[-1L]
  an <- names(args)
  if (is.null(an)) {
    an <- rep("", length(args))
  }
  if (nm == "tapply") {
    return(list(value = .nestArgNamedOrPos(args, an, "X", 1L), group = .nestArgNamedOrPos(args, an, "INDEX", 2L)))
  }
  if (nm == "by") {
    return(list(value = .nestArgNamedOrPos(args, an, "data", 1L), group = .nestArgNamedOrPos(args, an, "INDICES", 2L)))
  }
  if (nm == "aggregate") {
    first <- .nestArgNamedOrPos(args, an, "x", 1L)
    if (length(first) == 1L && identical(.nestCallName(first[[1L]]), "~") && length(first[[1L]]) == 3L) {
      return(list(value = list(first[[1L]][[2L]]), group = list(first[[1L]][[3L]])))
    }
    return(list(value = first, group = .nestArgNamedOrPos(args, an, "by", 2L)))
  }
  rest <- args[-1L]
  rest_names <- an[-1L]
  list(value = args[seq_len(min(1L, length(args)))], group = rest[!(rest_names %in% c("FUN", ".add", ".drop", "drop"))])
}

# The tokens that group in `expr`: those in the grouping arguments of a grouping
# or splitting call (never its data argument, so a filter inside it does not
# group), and in a `.by` / `by` argument. With `follow`, the values of the
# symbols it holds are searched too (up to three assignments back, never through
# a simulation: the events a simulation was given do not group a count over its
# output). `ungroup()` and the extraction of a column (`$`, `[[`) end the
# grouping, so what is inside them is not searched.
.nestGroupTokens <- function(expr, acc, depth = 0L, seen = character(), follow = TRUE) {
  found <- new.env(parent = emptyenv())
  found$tok <- character()
  .nestWalk(expr, function(e) {
    nm <- .nestCallName(e)
    if (isTRUE(nm %in% c("ungroup", "$", "[[", "@"))) {
      return(FALSE)
    }
    if (isTRUE(nm %in% c("%>%", "%<>%", "%T>%")) && length(e) == 3L && is.call(e[[3L]])) {
      piped <- as.call(c(e[[3L]][[1L]], list(e[[2L]]), as.list(e[[3L]])[-1L]))
      found$tok <- c(found$tok, .nestGroupTokens(piped, acc, depth, seen, follow))
      return(FALSE)
    }
    if (isTRUE(nm %in% .nestGroupCalls)) {
      found$tok <- c(found$tok, .nestTokensList(.nestGroupSplit(e)$group))
      return(FALSE)
    }
    if (!is.na(nm) && !is.null(names(e))) {
      found$tok <- c(found$tok, .nestTokensList(as.list(e)[names(e) %in% c(".by", "by")]))
    }
    if (follow && is.name(e) && depth < 3L && !(as.character(e) %in% seen)) {
      key <- as.character(e)
      for (r in acc$assign_expr[[key]]) {
        if (!.nestContainsSolve(r)) {
          found$tok <- c(found$tok, .nestGroupTokens(r, acc, depth + 1L, c(seen, key)))
        }
      }
    }
    invisible(NULL)
  })
  unique(found$tok)
}

# The two token sets of one count call: those that name what is counted (its own
# arguments and the columns iterated over around it) and those that group it (its
# own arguments, the grouping arguments of the grouping calls in the data it is
# summarised over, `.by` arguments, and the grouping arguments of tapply and the
# like). A symbol among the count's own arguments is not followed for its
# grouping: `g$etaiov` is a column, not grouped data. `anc` is the chain of
# enclosing calls.
.nestCountScope <- function(own_args, anc, acc) {
  own <- .nestTokensList(own_args)
  eta <- own
  group <- c(own, .nestGroupTokens(as.call(c(as.name("c"), own_args)), acc, follow = FALSE))
  for (a in anc) {
    an <- .nestCallName(a)
    if (is.na(an)) {
      next
    }
    args <- as.list(a)[-1L]
    first <- args[seq_len(min(1L, length(args)))]
    if (an %in% .nestScopeColumns) {
      eta <- c(eta, .nestTokensList(first))
    }
    if (an %in% .nestScopeData) {
      by_args <- if (is.null(names(args))) list() else args[names(args) %in% c(".by", "by")]
      group <- c(group, .nestGroupTokens(as.call(c(as.name("c"), first)), acc), .nestTokensList(by_args))
    }
    if (an %in% .nestGroupCalls) {
      gs <- .nestGroupSplit(a)
      eta <- c(eta, .nestTokensList(gs$value))
      group <- c(group, .nestTokensList(gs$group))
    }
  }
  list(eta = unique(eta), group = unique(group))
}

# One scope per count in `expr`: a call to a counting function, or a counting
# function passed by name to an iterator (`sapply(by_occasion, n_distinct)`).
.nestScopeWalk <- function(expr, anc, acc, out) {
  if (is.expression(expr)) {
    for (i in seq_along(expr)) {
      .nestScopeWalk(expr[[i]], anc, acc, out)
    }
    return(invisible(NULL))
  }
  if (!is.call(expr)) {
    return(invisible(NULL))
  }
  nm <- .nestCallName(expr)
  args <- as.list(expr)[-1L]
  if (!is.na(nm) && nm %in% .nestCountWords) {
    out$scopes[[length(out$scopes) + 1L]] <- .nestCountScope(args, anc, acc)
  }
  if (!is.na(nm) && nm %in% .nestIterators) {
    by_name <- vapply(args, function(a) is.name(a) && as.character(a) %in% .nestCountWords, logical(1L))
    data_args <- args[!by_name & !vapply(args, function(a) identical(.nestCallName(a), "function"), logical(1L))]
    if (any(by_name) && length(data_args) > 0L) {
      out$scopes[[length(out$scopes) + 1L]] <- .nestCountScope(data_args[1L], c(anc, list(expr)), acc)
    }
  }
  kids <- as.list(expr)
  for (i in seq_along(kids)) {
    if (!identical(kids[i], list(quote(expr = )))) {
      .nestScopeWalk(kids[[i]], c(anc, list(expr)), acc, out)
    }
  }
  invisible(NULL)
}

# The count scopes of a list of expressions.
.nestCountScopes <- function(exprs, acc) {
  out <- new.env(parent = emptyenv())
  out$scopes <- list()
  for (x in exprs) {
    .nestScopeWalk(x, list(), acc, out)
  }
  out$scopes
}

# Functions that join pieces into a string; the literal pieces of a model name
# built with one are fragments of the name.
.nestConcatenators <- c("paste", "paste0", "sprintf", "glue", "file.path", "str_c", "str_glue")

# Calls through which the value of a model name passes unchanged: the vector or
# table is built from their arguments (all of them, or the first), and a
# mutate()-like call adds columns to its first argument.
.nestPassAll <- c(
  "c",
  "list",
  "tibble",
  "data.frame",
  "data_frame",
  "tribble",
  "expand.grid",
  "expand_grid",
  "crossing",
  "bind_rows",
  "cbind",
  "rbind"
)
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
  "setdiff",
  "intersect",
  "union",
  "factor",
  "as.factor",
  "split",
  "first",
  "last",
  "ungroup",
  "group_by",
  "as.data.frame",
  "as_tibble",
  "suppressMessages",
  "suppressWarnings",
  "invisible",
  "identity",
  "(",
  "["
)
.nestPassBranches <- c("ifelse", "if_else", "if", "switch")

# Calls whose value is never a model name: numbers, logicals, functions, models
# and simulations.
.nestNoString <- c(
  "numeric",
  "integer",
  "double",
  "logical",
  "character",
  "as.numeric",
  "as.integer",
  "as.double",
  "as.logical",
  "seq_len",
  "seq_along",
  "seq",
  "nrow",
  "ncol",
  "length",
  "sum",
  "mean",
  "median",
  "sd",
  "var",
  "exp",
  "log",
  "log10",
  "sqrt",
  "abs",
  "round",
  "floor",
  "ceiling",
  "cumsum",
  "diff",
  "*",
  "+",
  "-",
  "/",
  "^",
  "%%",
  "%/%",
  ":",
  "pmax",
  "pmin",
  "n",
  "rnorm",
  "runif",
  "rbinom",
  "rlnorm",
  "rpois",
  "which",
  "any",
  "all",
  "is.na",
  "!",
  "==",
  "!=",
  "<",
  ">",
  "<=",
  ">=",
  "&",
  "&&",
  "|",
  "||",
  "function",
  "readModelDb",
  "rxode",
  "rxode2",
  "rxSolve",
  "zeroRe",
  "nlmixr2",
  "ini",
  "model"
)
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
# itself, with the functions around the binding. The tokens of a simulation are not kept: the events and model it was
# given say nothing about the draws in its output.
.nestBind <- function(acc, key, expr) {
  acc$assign[[key]] <- c(acc$assign[[key]], if (.nestContainsSolve(expr)) character() else .nestTokens(expr))
  acc$assign_expr[[key]] <- c(acc$assign_expr[[key]], list(expr))
  acc$assign_ctx[[key]] <- c(acc$assign_ctx[[key]], list(acc$ctx))
  invisible(NULL)
}

# The key of a function literal or formula: its text.
.nestFnKey <- function(fn) {
  paste(deparse(fn), collapse = "\n")
}

# Record that parameter `param` of the function with key `fk` takes the value of
# `expr` at one call. A parameter is looked up in the function it belongs to, never
# among the article's other variables of the same name.
.nestBindScoped <- function(acc, fk, param, expr) {
  acc$scoped[[fk]][[param]] <- c(acc$scoped[[fk]][[param]], list(list(expr = expr, ctx = acc$ctx)))
  invisible(NULL)
}

# The symbol an assignment target writes to: `x`, `x$a`, `x[i]`, `x[[i]]`.
.nestAssignTarget <- function(lhs) {
  while (is.call(lhs) && length(lhs) >= 2L && .nestCallName(lhs) %in% c("$", "[", "[[", "@")) {
    lhs <- lhs[[2L]]
  }
  if (is.name(lhs)) as.character(lhs) else NA_character_
}

# Whether `x` is a loader passed as a function (`readModelDb`,
# `nlmixr2lib::readModelDb`, `"readModelDb"`, or an alias of one).
.nestIsLoaderRef <- function(x, loaders) {
  if (is.name(x)) {
    return(as.character(x) %in% loaders)
  }
  if (is.character(x) && length(x) == 1L) {
    return(isTRUE(x %in% loaders))
  }
  is.call(x) && .nestCallName(x) %in% c("::", ":::") && as.character(x[[3L]]) %in% loaders
}

# Whether `rhs` makes its target an alias of a loader: the loader itself, or a
# lookup of one by name (`get("readModelDb")`, `match.fun(...)`).
.nestIsLoaderAlias <- function(rhs, loaders) {
  if (.nestIsLoaderRef(rhs, loaders)) {
    return(TRUE)
  }
  lookup <- c("get", "get0", "match.fun", "getFromNamespace", "getExportedValue", "eval", "parse")
  is.call(rhs) && isTRUE(.nestCallName(rhs) %in% lookup) && any(loaders %in% .nestTokens(rhs))
}

# Bind what a call node defines: an assignment (and the parameters of a function
# it assigns, and the alias of a loader), a loop variable, and the parameters of a
# user-defined function called here.
.nestCollectBindings <- function(e, nm, args, acc) {
  if (nm %in% c("<-", "=") && length(e) == 3L) {
    key <- .nestAssignTarget(e[[2L]])
    if (!is.na(key)) {
      rhs <- e[[3L]]
      .nestBind(acc, key, rhs)
      if (is.call(rhs) && identical(.nestCallName(rhs), "function")) {
        acc$formals[[key]] <- names(rhs[[2L]])
        acc$fn_key[[key]] <- .nestFnKey(rhs)
      }
      if (.nestIsLoaderAlias(rhs, acc$loaders)) {
        acc$loaders <- union(acc$loaders, key)
      }
    }
  }
  if (nm == "for" && length(e) >= 3L && is.name(e[[2L]])) {
    .nestBind(acc, as.character(e[[2L]]), e[[3L]])
  }
  if (nm %in% names(acc$formals)) {
    .nestBindCall(acc, acc$fn_key[[nm]], acc$formals[[nm]], args)
  }
  invisible(NULL)
}

# An iterator call (lapply, map, ...): the parameters of a function literal, the
# `.x` of a formula and the parameters of a user-defined function given to it
# take the values of the vector iterated over, and a loader given to it is a
# loader call on that vector. A formula is the function argument, not data.
.nestCollectIterator <- function(args, acc) {
  is_fn <- vapply(
    args,
    function(a) {
      .nestIsLoaderRef(a, acc$loaders) || isTRUE(.nestCallName(a) %in% c("function", "~"))
    },
    logical(1L)
  )
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
          .nestBindScoped(acc, .nestFnKey(a), p, d)
        }
      }
    }
    if (identical(.nestCallName(a), "~")) {
      for (d in data_args) {
        for (p in .nestFormulaParams) {
          .nestBindScoped(acc, .nestFnKey(a), p, d)
        }
      }
    }
    if (.nestIsLoaderRef(a, acc$loaders) && length(data_args) > 0L) {
      .nestAddLoaderArg(acc, data_args[[1L]])
    }
    if (is.name(a) && as.character(a) %in% names(acc$formals) && length(data_args) > 0L) {
      .nestBindCall(acc, acc$fn_key[[as.character(a)]], acc$formals[[as.character(a)]], data_args)
    }
  }
  invisible(NULL)
}

# Record a model-name argument given to a loader, with the functions around the
# call.
.nestAddLoaderArg <- function(acc, expr) {
  acc$loader_args[[length(acc$loader_args) + 1L]] <- list(expr = expr, ctx = acc$ctx)
  invisible(NULL)
}

# The parameters a formula `~ f(.x)` takes.
.nestFormulaParams <- c(".x", ".", "..1", ".y", "..2")

# Walk an expression collecting from every node, with the chain of functions
# (literals and formulas) the node is inside, so a parameter is seen only from
# within its function.
.nestCollectWalk <- function(expr, ctx, acc) {
  if (is.expression(expr)) {
    for (i in seq_along(expr)) {
      .nestCollectWalk(expr[[i]], ctx, acc)
    }
    return(invisible(NULL))
  }
  if (is.character(expr)) {
    acc$strings <- c(acc$strings, expr)
  } else if (is.name(expr)) {
    acc$strings <- c(acc$strings, as.character(expr))
  }
  if (!is.call(expr)) {
    return(invisible(NULL))
  }
  acc$ctx <- ctx
  .nestCollect(expr, acc)
  inner <- ctx
  head <- .nestCallName(expr)
  if (identical(head, "function") || identical(head, "~")) {
    fk <- .nestFnKey(expr)
    acc$fn_formals[[fk]] <- if (identical(head, "function")) names(expr[[2L]]) else .nestFormulaParams
    inner <- c(ctx, fk)
    if (identical(head, "function")) {
      acc$ctx <- inner
      .nestBindDefaults(acc, fk, expr[[2L]])
      # A default is code that runs: its strings and calls are read like the body's.
      defaults <- as.list(expr[[2L]])
      for (i in seq_along(defaults)) {
        if (!identical(unname(defaults[i]), list(quote(expr = )))) {
          .nestCollectWalk(defaults[[i]], inner, acc)
        }
      }
    }
  }
  kids <- as.list(expr)
  for (i in seq_along(kids)) {
    if (!identical(kids[i], list(quote(expr = )))) {
      .nestCollectWalk(kids[[i]], inner, acc)
    }
  }
  invisible(NULL)
}

# A parameter's default is one of its values: bind each default of a function's
# parameters (in the function's own scope).
.nestBindDefaults <- function(acc, fk, formals) {
  defaults <- as.list(formals)
  for (i in seq_along(defaults)) {
    if (!identical(unname(defaults[i]), list(quote(expr = )))) {
      .nestBindScoped(acc, fk, names(defaults)[[i]], defaults[[i]])
    }
  }
  invisible(NULL)
}

# A function literal called where it is written, `(\(x) f(x))(v)`: its
# parameters take the arguments.
.nestCollectLiteralCall <- function(head, args, acc) {
  if (is.call(head) && identical(.nestCallName(head), "(") && length(head) == 2L) {
    head <- head[[2L]]
  }
  if (is.call(head) && identical(.nestCallName(head), "function")) {
    .nestBindCall(acc, .nestFnKey(head), names(head[[2L]]), args)
  }
  invisible(NULL)
}

# `lhs %>% rhs`: collected as the call with `lhs` inserted as the first argument,
# so a loader or an iterator on the right of a magrittr pipe sees its data.
.nestCollectPipe <- function(lhs, rhs, acc) {
  piped <- if (.nestIsLoaderRef(rhs, acc$loaders)) {
    as.call(list(rhs, lhs))
  } else if (is.call(rhs)) {
    as.call(c(rhs[[1L]], list(lhs), as.list(rhs)[-1L]))
  } else {
    NULL
  }
  if (!is.null(piped)) {
    .nestCollect(piped, acc)
  }
  invisible(NULL)
}

# Whether a zeroRe() call zeroes the omegas: its `which` is absent (the default
# zeroes both) or holds "omega". `zeroRe(m, which = "sigma")` keeps the omegas.
.nestZeroesOmega <- function(x) {
  if (!identical(.nestCallName(x), "zeroRe")) {
    return(FALSE)
  }
  args <- as.list(x)[-1L]
  an <- names(args)
  if (is.null(an)) {
    an <- rep("", length(args))
  }
  which <- .nestArgNamedOrPos(args[-1L], an[-1L], "which", 1L)
  length(which) == 0L || "omega" %in% .nestTokens(which[[1L]], strings_only = TRUE)
}

# Whether `x` is a solver passed as a function: `rxSolve`, `"rxSolve"`,
# `rxode2::rxSolve`.
.nestIsSimRef <- function(x) {
  if (is.name(x)) {
    return(as.character(x) %in% .nestSimCalls)
  }
  if (is.character(x) && length(x) == 1L) {
    return(isTRUE(x %in% .nestSimCalls))
  }
  is.call(x) && isTRUE(.nestCallName(x) %in% c("::", ":::")) && as.character(x[[3L]]) %in% .nestSimCalls
}

# Calls that apply or look up a function given as an argument.
.nestApplyCalls <- c(.nestIterators, "do.call", "match.fun", "get", "pmap", "invoke_map", "Reduce", "Position")

# Whether a call simulates: rxSolve, simulate, solve or predict on something other
# than a model with its omegas zeroed, a solver passed to an iterator or do.call(),
# or an nlmixr2() fit with est = "rxSolve".
.nestIsSimulation <- function(e, nm) {
  if (nm %in% .nestSimCalls) {
    return(length(e) < 2L || !.nestZeroesOmega(e[[2L]]))
  }
  args <- as.list(e)[-1L]
  if (nm %in% .nestApplyCalls && any(vapply(args, .nestIsSimRef, logical(1L)))) {
    return(TRUE)
  }
  nm %in% c("nlmixr2", "nlmixr") && any(vapply(args, function(a) identical(a, "rxSolve"), logical(1L)))
}

# Collect from one call node: bindings, assertions, loader arguments and
# simulation calls.
.nestCollect <- function(e, acc) {
  if (!is.call(e)) {
    return(invisible(NULL))
  }
  nm <- .nestCallName(e)
  args <- as.list(e)[-1L]
  if (is.na(nm)) {
    return(.nestCollectLiteralCall(e[[1L]], args, acc))
  }
  .nestCollectBindings(e, nm, args, acc)
  if (nm %in% c("%>%", "%<>%", "%T>%") && length(args) == 2L) {
    .nestCollectPipe(args[[1L]], args[[2L]], acc)
  }
  if (nm %in% .nestIterators) {
    .nestCollectIterator(args, acc)
  }
  if (nm == "do.call" && length(args) >= 2L && .nestIsLoaderRef(args[[1L]], acc$loaders)) {
    .nestAddLoaderArg(acc, args[[2L]])
  }
  if (nm %in% c("stopifnot", "assert_that")) {
    acc$asserts <- c(acc$asserts, args)
  } else if (grepl("^expect_", nm)) {
    acc$asserts[[length(acc$asserts) + 1L]] <- e
  }
  if (nm %in% acc$loaders) {
    arg <- .nestLoaderArg(e)
    if (!is.null(arg) && !is.character(arg)) {
      .nestAddLoaderArg(acc, arg)
    }
  }
  if (.nestIsSimulation(e, nm)) {
    acc$simulates <- TRUE
  }
  invisible(NULL)
}

# Bind the parameters of a user-defined function to the arguments of one call
# to it, by name and then by position.
.nestBindCall <- function(acc, fk, formals, args) {
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
      .nestBindScoped(acc, fk, target, args[[i]])
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
    why = unique(unlist(lapply(parts, `[[`, "why"), use.names = FALSE)),
    hint = unique(unlist(lapply(parts, `[[`, "hint"), use.names = FALSE))
  )
}

# An unseen value, with the string constants found inside it (`hint`).
.nestUnknownName <- function(why, hint = character()) {
  list(pat = character(), unknown = TRUE, why = why, hint = hint)
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
    return(list(pat = .nestEscape(x[!is.na(x)]), unknown = FALSE, why = character(), hint = character()))
  }
  if (is.numeric(x) || is.logical(x) || is.null(x)) {
    return(list(pat = character(), unknown = FALSE, why = character(), hint = character()))
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

# The value of an expression evaluated in a context of enclosing functions.
.nestNameValueIn <- function(expr, ctx, acc, depth, seen, col) {
  old <- acc$ctx
  acc$ctx <- ctx
  on.exit(acc$ctx <- old)
  .nestNameValue(expr, acc, depth + 1L, seen, col)
}

# The value of a symbol: the values bound to it. A parameter of an enclosing
# function takes only the values given to that function; any other symbol takes
# the values assigned to it. One value that cannot be seen makes the symbol unseen,
# whatever else is bound to it.
.nestNameValueSymbol <- function(key, acc, depth, seen, col) {
  if (key %in% seen) {
    return(list(pat = character(), unknown = FALSE, why = character(), hint = character()))
  }
  bindings <- NULL
  owner <- Find(function(fk) key %in% acc$fn_formals[[fk]], rev(acc$ctx))
  if (!is.null(owner)) {
    bindings <- acc$scoped[[owner]][[key]]
    if (length(bindings) == 0L) {
      return(.nestUnknownName(paste0("parameter `", key, "` is never given a value the article shows")))
    }
  } else {
    # An assignment made inside a function is visible only inside that function.
    here <- acc$ctx
    visible <- vapply(
      acc$assign_ctx[[key]],
      function(ctx) length(ctx) <= length(here) && all(ctx == here[seq_along(ctx)]),
      logical(1L)
    )
    rhs <- acc$assign_expr[[key]][visible]
    if (length(rhs) == 0L) {
      return(.nestUnknownName(paste0("`", key, "` has no value bound in the article")))
    }
    bindings <- Map(function(e, ctx) list(expr = e, ctx = ctx), rhs, acc$assign_ctx[[key]][visible])
  }
  .nestMergeNames(lapply(
    bindings,
    function(b) .nestNameValueIn(b$expr, b$ctx, acc, depth, c(seen, key), col)
  ))
}

# The value of `names(x)`: the argument names of a table or vector built in the
# article, when `x` is one; otherwise unseen.
.nestNameValueNames <- function(x, acc, depth, seen) {
  parts <- list(x)
  if (is.name(x)) {
    key <- as.character(x)
    if (key %in% seen || length(acc$assign_expr[[key]]) == 0L) {
      return(.nestUnknownName(paste0("`names(", key, ")` has no value bound in the article")))
    }
    parts <- acc$assign_expr[[key]]
  }
  built <- vapply(parts, function(p) is.call(p) && isTRUE(.nestCallName(p) %in% .nestPassAll), logical(1L))
  if (!all(built)) {
    return(.nestUnknownName("value of `names()`"))
  }
  an <- unlist(lapply(parts, function(p) names(as.list(p)[-1L])))
  if (length(an) == 0L || any(!nzchar(an))) {
    return(.nestUnknownName("value of `names()`"))
  }
  list(pat = .nestEscape(unique(an)), unknown = FALSE, why = character(), hint = character())
}

# The value of an iterator call (lapply, vapply, map, ...): the value of the
# function applied, a loader giving models (never names), or unseen.
.nestNameValueIterated <- function(args, acc, depth, seen) {
  fns <- args[vapply(args, function(a) is.call(a) && identical(.nestCallName(a), "function"), logical(1L))]
  if (length(fns) > 0L) {
    return(.nestMergeNames(lapply(fns, function(f) {
      .nestNameValueIn(f[[3L]], c(acc$ctx, .nestFnKey(f)), acc, depth, seen, NA_character_)
    })))
  }
  if (any(vapply(args, .nestIsLoaderRef, logical(1L), loaders = acc$loaders))) {
    return(list(pat = character(), unknown = FALSE, why = character(), hint = character()))
  }
  users <- args[vapply(args, function(a) is.name(a) && as.character(a) %in% names(acc$formals), logical(1L))]
  if (length(users) > 0L) {
    return(.nestNameValueAll(users, acc, depth, seen, NA_character_))
  }
  .nestUnknownName("value of an iterator call")
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

# The value of a call that passes its arguments through: a table or vector built
# from them (the named argument when a column is selected), a mutate()-like call,
# or a call that returns its first argument; NULL for any other call.
.nestNameValuePass <- function(nm, args, an, named, col, acc, depth, seen) {
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
  NULL
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
  if (nm == "~" || nm %in% .nestNoString) {
    return(list(pat = character(), unknown = FALSE, why = character(), hint = character()))
  }
  if (nm == "names" && length(args) == 1L) {
    return(.nestNameValueNames(args[[1L]], acc, depth, seen))
  }
  if (nm %in% .nestConcatenators) {
    return(list(pat = .nestConcatPattern(x), unknown = FALSE, why = character(), hint = character()))
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
  passed <- .nestNameValuePass(nm, args, an, named, col, acc, depth, seen)
  if (!is.null(passed)) {
    return(passed)
  }
  if (nm %in% .nestIterators) {
    return(.nestNameValueIterated(args, acc, depth, seen))
  }
  if (nm %in% names(acc$formals)) {
    defs <- Filter(function(r) identical(.nestCallName(r), "function"), acc$assign_expr[[nm]])
    if (length(defs) > 0L) {
      return(.nestMergeNames(lapply(defs, function(f) {
        .nestNameValueIn(f[[3L]], c(acc$ctx, .nestFnKey(f)), acc, depth, seen, col)
      })))
    }
  }
  .nestUnknownName(paste0("value of `", nm, "()`"), .nestTokens(x, strings_only = TRUE))
}

# Why a model name given to a loader cannot be matched against `models`:
# list(reasons, firm). `reasons` is empty when the name can be matched: some value
# behind it is unseen (see .nestNameValue()), or a name built from pieces can be
# one of `models` (`paste0("Smith_2020_", drug)` cannot when no model starts with
# `Smith_2020_`). `firm` are the reasons no exception can excuse: a built pattern
# that can be a below-subject model, or a string inside an unseen value that is
# part of one. Literal names are matched elsewhere.
.nestDynamicName <- function(arg, acc, models) {
  value <- .nestNameValueIn(arg$expr, arg$ctx, acc, 0L, character(), NA_character_)
  built <- value$pat[grepl(".*", value$pat, fixed = TRUE)]
  firm <- character()
  hit <- unique(unlist(lapply(built, function(p) models[grepl(paste0("^", p, "$"), models, perl = TRUE)])))
  if (length(hit) > 0L) {
    firm <- paste0("a name built from pieces can be ", paste(hit, collapse = ", "))
  }
  hints <- value$hint[nchar(value$hint) >= 5L]
  inside <- unique(unlist(lapply(hints, function(h) models[grepl(h, models, fixed = TRUE)])))
  if (length(inside) > 0L) {
    firm <- c(firm, paste0("a string inside an unseen value is part of ", paste(inside, collapse = ", ")))
  }
  list(reasons = c(if (value$unknown) value$why, firm), firm = firm)
}

.nestNoFindings <- function() {
  list(
    models = character(),
    etas = character(),
    simulates = FALSE,
    uncovered = character(),
    dynamic = character(),
    dynamic_firm = character(),
    unparseable = integer()
  )
}

# The article's findings:
#   models, etas  the below-subject models it names and their `| G` etas;
#   simulates     it calls rxSolve / simulate / solve / predict on something other
#                 than a zeroRe()'d model, or nlmixr2(est = "rxSolve");
#   uncovered     etas with no per-level distinct-draw assertion. A count call
#                 (n_distinct, unique, ...) covers an eta when the eta is in the
#                 count's own arguments (or the columns iterated over around it,
#                 followed back through assignments) and its grouping column is in
#                 the count's own arguments, in a grouping or splitting call around
#                 it (group_by, split, tapply, `.by`), or in such a call in the
#                 value of a symbol it holds. A column that only filters the data,
#                 or arrives through the events a simulation was given, does not
#                 group the count. The count is looked for in each stopifnot()
#                 argument or expect_*() call and in the assignments those reach,
#                 so a summary computed earlier and asserted later counts, while a
#                 `unique()` that only builds the event table does not;
#   dynamic       the reasons a readModelDb() / modellib() call (also through a
#                 pipe, an alias, do.call() or an iterator) has a model name that
#                 cannot be matched (see .nestDynamicName()); empty when none;
#   unparseable   lines of running R chunks that do not parse.
# Most articles neither name a below-subject model nor pass a non-literal name
# to a loader; those are not parsed.
.nestArticleFindings <- function(path, nested) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  below <- nested[!nested$above_subject, , drop = FALSE]
  names_model <- any(vapply(unique(below$model), grepl, logical(1L), x = text, fixed = TRUE))
  non_literal_loader <- grepl(
    paste0(
      "(?:\\|>|%>%|%<>%|%T>%)\\s*(?:\\w+::)?(?:readModelDb|modellib)\\b|",
      "(?:readModelDb|modellib)(?!\\s*\\(\\s*[\"')])"
    ),
    text,
    perl = TRUE
  )
  if (!names_model && !non_literal_loader) {
    return(.nestNoFindings())
  }
  parsed <- .nestArticleExprs(path)
  acc <- new.env(parent = emptyenv())
  acc$strings <- character()
  acc$assign <- list()
  acc$assign_expr <- list()
  acc$assign_ctx <- list()
  acc$formals <- list()
  acc$fn_key <- list()
  acc$fn_formals <- list()
  acc$scoped <- list()
  acc$ctx <- character()
  acc$asserts <- list()
  acc$loader_args <- list()
  acc$loaders <- .nestLoaders
  acc$simulates <- FALSE
  for (ch in parsed$chunks) {
    .nestCollectWalk(ch, character(), acc)
  }
  acc$assign <- lapply(acc$assign, unique)
  used <- nested[nested$model %in% acc$strings & !nested$above_subject, , drop = FALSE]
  scopes <- unlist(
    lapply(acc$asserts, function(a) {
      tok <- .nestExpandTokens(.nestTokens(a), acc$assign)
      reached <- unlist(acc$assign_expr[intersect(tok, names(acc$assign_expr))], recursive = FALSE)
      .nestCountScopes(c(list(a), reached), acc)
    }),
    recursive = FALSE
  )
  eta_tokens <- lapply(scopes, function(sc) .nestExpandTokens(sc$eta, acc$assign))
  covered <- vapply(
    seq_len(nrow(used)),
    function(i) {
      any(vapply(
        seq_along(scopes),
        function(k) used$eta[[i]] %in% eta_tokens[[k]] && used$group[[i]] %in% scopes[[k]]$group,
        logical(1L)
      ))
    },
    logical(1L)
  )
  names_found <- lapply(acc$loader_args, .nestDynamicName, acc = acc, models = unique(below$model))
  dynamic <- unique(unlist(lapply(names_found, `[[`, "reasons")))
  dynamic_firm <- unique(unlist(lapply(names_found, `[[`, "firm")))
  list(
    models = unique(used$model),
    etas = unique(used$eta),
    simulates = acc$simulates,
    uncovered = unique(used$eta[!covered]),
    dynamic = dynamic,
    dynamic_firm = dynamic_firm,
    unparseable = parsed$unparseable
  )
}

# Whether an article's unresolved loader names are excused: it has an entry in
# `exceptions` and no below-subject model starts with the entry's prefix. A
# built name or a string that can be a below-subject model is never excused.
.nestExcepted <- function(path, nested, exceptions) {
  prefix <- exceptions$prefix[exceptions$article == basename(path)]
  length(prefix) == 1L && !any(startsWith(nested$model[!nested$above_subject], prefix))
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
.nestDrawViolations <- function(
  model_root,
  article_paths,
  above = .nestAboveSubjectGroups,
  exceptions = .nestDynamicLoaderExceptions
) {
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
    excused <- length(f$dynamic_firm) == 0L && .nestExcepted(p, nested, exceptions)
    if (f$simulates && length(f$dynamic) > 0L && !excused) {
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
