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

# Depth-first walk over a language object, calling `fn` on every node. Empty
# arguments (the `x[, 1]` form) are skipped.
.nestWalk <- function(expr, fn) {
  if (is.expression(expr)) {
    for (i in seq_along(expr)) {
      .nestWalk(expr[[i]], fn)
    }
    return(invisible(NULL))
  }
  fn(expr)
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
# not mistaken for one.
.nestedEtas <- function(path) {
  lines <- readLines(path, warn = FALSE)
  empty <- data.frame(eta = character(), group = character(), stringsAsFactors = FALSE)
  if (!any(grepl("~[^#]*\\|", lines))) {
    return(empty)
  }
  exprs <- parse(text = lines, keep.source = FALSE)
  acc <- new.env(parent = emptyenv())
  acc$rows <- list()
  .nestWalk(exprs, function(e) {
    if (identical(.nestCallName(e), "ini")) {
      .nestWalk(e, function(s) {
        if (!is.call(s) || !identical(s[[1L]], as.name("~")) || length(s) != 3L) {
          return(invisible(NULL))
        }
        rhs <- s[[3L]]
        if (is.call(rhs) && identical(rhs[[1L]], as.name("|")) && is.name(rhs[[3L]])) {
          acc$rows[[length(acc$rows) + 1L]] <- data.frame(
            eta = all.vars(s[[2L]]),
            group = as.character(rhs[[3L]]),
            stringsAsFactors = FALSE
          )
        }
        invisible(NULL)
      })
    }
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

# The parsed R chunks of an Rmd, skipping `eval = FALSE` chunks and any chunk
# that does not parse.
.nestArticleExprs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  inchunk <- FALSE
  skip <- FALSE
  buf <- character()
  for (ln in lines) {
    if (!inchunk && grepl("^```\\{[rR][ ,}]", ln)) {
      inchunk <- TRUE
      skip <- grepl("eval\\s*=\\s*(FALSE|F)\\b", ln)
      buf <- character()
    } else if (inchunk && grepl("^```\\s*$", ln)) {
      inchunk <- FALSE
      if (!skip) {
        parsed <- tryCatch(parse(text = buf, keep.source = FALSE), error = function(e) NULL)
        if (!is.null(parsed)) {
          out[[length(out) + 1L]] <- parsed
        }
      }
    } else if (inchunk) {
      buf <- c(buf, ln)
    }
  }
  out
}

# Every symbol and string constant in a language object.
.nestTokens <- function(expr) {
  acc <- new.env(parent = emptyenv())
  acc$tok <- character()
  .nestWalk(expr, function(e) {
    if (is.name(e)) {
      acc$tok <- c(acc$tok, as.character(e))
    } else if (is.character(e)) {
      acc$tok <- c(acc$tok, e)
    }
  })
  unique(acc$tok)
}

# Follow the symbols in `tokens` back through the assignments that define them,
# up to four levels, and return every token met on the way.
.nestExpandTokens <- function(tokens, assign) {
  seen <- character()
  todo <- tokens
  for (depth in seq_len(4L)) {
    nxt <- intersect(setdiff(todo, seen), names(assign))
    seen <- union(seen, todo)
    todo <- unique(unlist(lapply(nxt, function(k) {
      unlist(lapply(assign[[k]], .nestTokens))
    })))
    if (length(todo) == 0L) {
      break
    }
    tokens <- union(tokens, todo)
  }
  tokens
}

# The article's findings: which below-subject models it names, whether it
# simulates (rxSolve / simulate / solve whose first argument is not zeroRe()),
# and which of the below-subject etas are covered by a per-level distinct-draw
# assertion. An assertion is a stopifnot() / expect_*() / assert_that() call
# whose arguments, followed back through up to four assignments, mention both
# the eta and a distinct-value count.
.nestArticleFindings <- function(path, nested) {
  # Most articles name no below-subject model; skip parsing those.
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  below <- nested$model[!nested$above_subject]
  if (!any(vapply(unique(below), grepl, logical(1L), x = text, fixed = TRUE))) {
    return(list(models = character(), etas = character(), simulates = FALSE, uncovered = character()))
  }
  chunks <- .nestArticleExprs(path)
  acc <- new.env(parent = emptyenv())
  acc$strings <- character()
  acc$assign <- list()
  acc$asserts <- list()
  acc$simulates <- FALSE
  for (ch in chunks) {
    .nestWalk(ch, function(e) {
      if (is.character(e)) {
        acc$strings <- c(acc$strings, e)
      }
      nm <- .nestCallName(e)
      if (is.na(nm)) {
        return(invisible(NULL))
      }
      if (nm %in% c("<-", "=") && length(e) == 3L && is.name(e[[2L]])) {
        key <- as.character(e[[2L]])
        acc$assign[[key]] <- c(acc$assign[[key]], list(e[[3L]]))
      }
      if (nm %in% c("stopifnot", "assert_that") || grepl("^expect_", nm)) {
        acc$asserts[[length(acc$asserts) + 1L]] <- e
      }
      if (nm %in% c("rxSolve", "simulate", "solve") && length(e) >= 2L) {
        first <- e[[2L]]
        if (!identical(.nestCallName(first), "zeroRe")) {
          acc$simulates <- TRUE
        }
      }
    })
  }
  used <- nested[nested$model %in% acc$strings & !nested$above_subject, , drop = FALSE]
  count_words <- c("n_distinct", "unique", "uniqueN", "duplicated", "anyDuplicated")
  covered <- character()
  for (a in acc$asserts) {
    tok <- .nestExpandTokens(.nestTokens(a), acc$assign)
    if (any(count_words %in% tok)) {
      covered <- union(covered, intersect(tok, used$eta))
    }
  }
  list(
    models = unique(used$model),
    etas = unique(used$eta),
    simulates = acc$simulates,
    uncovered = setdiff(unique(used$eta), covered)
  )
}

# One row per article that simulates a below-subject `| G` model without a
# per-level distinct-draw assertion on every one of that model's etas.
.nestDrawViolations <- function(model_root, article_paths, above = .nestAboveSubjectGroups) {
  nested <- .nestedModels(model_root, above)
  rows <- lapply(article_paths, function(p) {
    f <- .nestArticleFindings(p, nested)
    if (length(f$models) == 0L || !f$simulates || length(f$uncovered) == 0L) {
      return(NULL)
    }
    data.frame(
      article = basename(p),
      models = paste(f$models, collapse = ", "),
      etas = paste(f$uncovered, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows)) {
    return(data.frame(article = character(), models = character(), etas = character()))
  }
  rows
}
