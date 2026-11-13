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
      acc$tok <- c(acc$tok, e)
    }
    invisible(NULL)
  })
  unique(acc$tok)
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

# Functions that build a string; a model name passed through one of them is
# not a literal the gate can match.
.nestStringBuilders <- c(
  "paste",
  "paste0",
  "sprintf",
  "glue",
  "file.path",
  "gsub",
  "sub",
  "toupper",
  "tolower",
  "substr",
  "substring",
  "strsplit",
  "format",
  "formatC",
  "str_c",
  "str_glue",
  "str_replace",
  "str_sub"
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

# Whether a model name given to a loader is built at run time in a way that
# could produce one of `models`: the argument, followed back through up to four
# assignments, calls a string function, and either holds no literal of five or
# more characters (the name comes from data) or holds one that is part of one of
# the model names. `paste0("Smith_2020_", drug)` is not such a call when no model
# of interest starts with `Smith_2020_`. Format placeholders (`%s`, `{x}`) split
# a literal into fragments. `assign` and `assign_str` map a symbol to the
# tokens, and to the string constants, of every right-hand side assigned to it.
.nestDynamicName <- function(arg, assign, assign_str, models) {
  tok <- .nestExpandTokens(.nestTokens(arg), assign)
  if (!any(.nestStringBuilders %in% tok)) {
    return(FALSE)
  }
  lits <- c(
    .nestTokens(arg, strings_only = TRUE),
    unlist(assign_str[intersect(tok, names(assign_str))], use.names = FALSE)
  )
  frags <- unlist(strsplit(lits, "%[-+0-9.]*[a-zA-Z]|\\{[^}]*\\}"))
  frags <- frags[nchar(frags) >= 5L]
  if (length(frags) == 0L) {
    return(TRUE)
  }
  any(vapply(frags, function(f) any(grepl(f, models, fixed = TRUE)), logical(1L)))
}

.nestNoFindings <- function() {
  list(
    models = character(),
    etas = character(),
    simulates = FALSE,
    uncovered = character(),
    dynamic = FALSE,
    unparseable = integer()
  )
}

# The article's findings:
#   models, etas  the below-subject models it names and their `| G` etas;
#   simulates     it calls rxSolve / simulate / solve with a first argument that
#                 is not zeroRe();
#   uncovered     etas with no per-level distinct-draw assertion. An assertion
#                 is one stopifnot() / assert_that() argument, or one expect_*()
#                 call, whose tokens (followed back through up to four
#                 assignments) include a distinct-value count, the eta and the
#                 eta's grouping column, so a count of subjects or a finiteness
#                 check on the eta alone does not cover it;
#   dynamic       a readModelDb() / modellib() call whose name argument is built
#                 by a string function in a way that could name one of the
#                 below-subject models (see .nestDynamicName());
#   unparseable   lines of running R chunks that do not parse.
# Most articles neither name a below-subject model nor pass a non-literal name
# to a loader; those are not parsed.
.nestArticleFindings <- function(path, nested) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  below <- nested[!nested$above_subject, , drop = FALSE]
  names_model <- any(vapply(unique(below$model), grepl, logical(1L), x = text, fixed = TRUE))
  non_literal_loader <- grepl("(?:readModelDb|modellib)\\s*\\(\\s*[^\"'\\s)]", text, perl = TRUE)
  if (!names_model && !non_literal_loader) {
    return(.nestNoFindings())
  }
  parsed <- .nestArticleExprs(path)
  acc <- new.env(parent = emptyenv())
  acc$strings <- character()
  acc$assign <- list()
  acc$assign_str <- list()
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
      nm <- .nestCallName(e)
      if (is.na(nm)) {
        return(invisible(NULL))
      }
      if (nm %in% c("<-", "=") && length(e) == 3L && is.name(e[[2L]])) {
        key <- as.character(e[[2L]])
        acc$assign[[key]] <- c(acc$assign[[key]], .nestTokens(e[[3L]]))
        acc$assign_str[[key]] <- c(acc$assign_str[[key]], .nestTokens(e[[3L]], strings_only = TRUE))
      }
      if (nm %in% c("stopifnot", "assert_that")) {
        acc$asserts <- c(acc$asserts, as.list(e)[-1L])
      } else if (grepl("^expect_", nm)) {
        acc$asserts[[length(acc$asserts) + 1L]] <- e
      }
      if (nm %in% c("readModelDb", "modellib")) {
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
    })
  }
  acc$assign <- lapply(acc$assign, unique)
  acc$assign_str <- lapply(acc$assign_str, unique)
  used <- nested[nested$model %in% acc$strings & !nested$above_subject, , drop = FALSE]
  count_words <- c("n_distinct", "unique", "uniqueN", "duplicated", "anyDuplicated")
  assert_tokens <- lapply(acc$asserts, function(a) .nestExpandTokens(.nestTokens(a), acc$assign))
  has_count <- vapply(assert_tokens, function(tok) any(count_words %in% tok), logical(1L))
  covered <- vapply(
    seq_len(nrow(used)),
    function(i) {
      any(vapply(
        assert_tokens[has_count],
        function(tok) all(c(used$eta[[i]], used$group[[i]]) %in% tok),
        logical(1L)
      ))
    },
    logical(1L)
  )
  dynamic <- any(vapply(
    acc$loader_args,
    .nestDynamicName,
    logical(1L),
    assign = acc$assign,
    assign_str = acc$assign_str,
    models = unique(below$model)
  ))
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
    if (f$simulates && f$dynamic) {
      out[[length(out) + 1L]] <- .nestViolationRow(
        p,
        "",
        "",
        "unverifiable: the model name given to readModelDb() or modellib() is built at run time"
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
