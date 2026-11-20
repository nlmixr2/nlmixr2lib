# Read the R code an Rmd article runs, without evaluating it and without knitr
# (knitr is only in Suggests). Used by the `nesting_assertion` rule of
# checkModelConventions(); any check that has to know what an article does reads
# it through .articleExprs().

# Spellings of FALSE that a YAML 1.1 parser (the one knitr uses for `#| eval:`)
# reads as a logical.
.yamlFalse <- c("false", "False", "FALSE", "no", "No", "NO", "n", "N", "off", "Off", "OFF")

# Does a chunk not run? `eval = FALSE` (or `F`) in the header, or `#| eval: <v>` /
# `#| eval: !expr FALSE` / `#| eval = FALSE` in the body. An `eval` that is
# computed (`eval = have_data`) is taken to run.
.chunkHeaderSkips <- function(header) {
  grepl("(?<![A-Za-z0-9_.])eval\\s*=\\s*(FALSE|F)\\b", header, perl = TRUE)
}

.chunkBodyOptionSkips <- function(line) {
  yaml <- regmatches(line, regexec("^#\\|\\s*eval\\s*:\\s*(!expr\\s+)?(\\S+)\\s*$", line))[[1L]]
  if (length(yaml) == 3L) {
    # `!expr FALSE` is an R expression, which knitr evaluates; any other value is
    # a YAML scalar.
    return(if (nzchar(yaml[[2L]])) yaml[[3L]] %in% c("FALSE", "F") else yaml[[3L]] %in% .yamlFalse)
  }
  grepl("^#\\|\\s*eval\\s*=\\s*(FALSE|F)\\s*$", line)
}

# Does the call `e` have an argument named `eval`?
.hasEvalArg <- function(e) {
  "eval" %in% names(e)
}

# Is `e` a call `knitr::opts_chunk$set(...)` or `opts_chunk$set(...)` that sets
# the chunk option `eval`, as a named argument or inside a `list()` argument? It
# then decides which of the later chunks run, which the reader cannot follow.
.setsEvalOption <- function(e) {
  head <- e[[1L]]
  if (!(is.call(head) && length(head) == 3L && identical(head[[1L]], as.name("$")) && identical(head[[3L]], as.name("set")))) {
    return(FALSE)
  }
  target <- head[[2L]]
  is_opts <- (is.name(target) && identical(as.character(target), "opts_chunk")) ||
    (is.call(target) &&
      length(target) == 3L &&
      is.name(target[[1L]]) &&
      as.character(target[[1L]]) %in% c("::", ":::") &&
      identical(as.character(target[[3L]]), "opts_chunk"))
  if (!is_opts) {
    return(FALSE)
  }
  args <- as.list(e)[-1L]
  lists <- Filter(function(a) is.call(a) && identical(.calledName(a), "list"), args)
  .hasEvalArg(e) || any(vapply(lists, .hasEvalArg, logical(1L)))
}

# The running R chunks of an Rmd, parsed, as list(chunks, unreadable).
#
# `chunks` holds one list of expressions per chunk; magrittr `%>%` calls are
# rewritten into the nested call the native pipe already parses to (.unpipe()),
# so a caller sees every call the same way. A chunk is fenced by three or more
# backticks and closes on a line of at least as many; text inside a fenced block
# that is not an R chunk is not read.
#
# `unreadable` holds, for each running chunk the reader cannot take at face
# value, the line it opens on, named by the reason: it does not parse, it is
# never closed, it pipes into a form `.unpipe()` does not rewrite, or it sets the
# knitr chunk option `eval` for the chunks after it (so which chunks run is
# unknown). A caller must treat such a chunk as unread code, not as an empty one.
# A chunk that parses is still listed in `chunks`, with `%>%` left where it was
# not rewritten.
.articleExprs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  chunks <- list()
  unreadable <- integer()
  fence <- 0L # backticks of the open fence, 0 when none
  is_r <- FALSE
  skip <- FALSE
  start <- 0L
  buf <- character()
  for (i in seq_along(lines)) {
    if (fence == 0L) {
      r <- regmatches(lines[[i]], regexec("^\\s*(`{3,})\\{[rR][ ,}]", lines[[i]]))[[1L]]
      any_fence <- regmatches(lines[[i]], regexec("^\\s*(`{3,})", lines[[i]]))[[1L]]
      if (length(r) == 2L) {
        fence <- nchar(r[[2L]])
        is_r <- TRUE
        skip <- .chunkHeaderSkips(lines[[i]])
        start <- i
        buf <- character()
      } else if (length(any_fence) == 2L) {
        fence <- nchar(any_fence[[2L]])
        is_r <- FALSE
      }
      next
    }
    closes <- regmatches(lines[[i]], regexec("^\\s*(`{3,})\\s*$", lines[[i]]))[[1L]]
    if (length(closes) == 2L && nchar(closes[[2L]]) >= fence) {
      if (is_r && !skip) {
        parsed <- tryCatch(parse(text = buf, keep.source = FALSE), error = function(e) NULL)
        if (is.null(parsed)) {
          unreadable <- c(unreadable, "does not parse" = start)
        } else {
          chunk <- lapply(as.list(parsed), .unpipe)
          chunks[[length(chunks) + 1L]] <- chunk
          unreadable <- c(unreadable, .chunkUnreadable(chunk, start))
        }
      }
      fence <- 0L
      is_r <- FALSE
    } else if (is_r) {
      if (.chunkBodyOptionSkips(lines[[i]])) {
        skip <- TRUE
      }
      buf <- c(buf, lines[[i]])
    }
  }
  if (fence > 0L && is_r && !skip) {
    unreadable <- c(unreadable, "is never closed" = start)
  }
  list(chunks = chunks, unreadable = unreadable)
}

# The reasons a parsed, unpiped chunk that opens on line `start` cannot be taken
# at face value, as a named integer (empty when it can).
.chunkUnreadable <- function(chunk, start) {
  acc <- new.env(parent = emptyenv())
  acc$pipe <- FALSE
  acc$eval <- FALSE
  .walkCalls(chunk, function(e) {
    if (identical(.calledName(e), "%>%")) {
      acc$pipe <- TRUE
    }
    if (.setsEvalOption(e)) {
      acc$eval <- TRUE
    }
  })
  out <- integer()
  if (acc$pipe) {
    out <- c(out, "pipes into a `{`, `(` or `pkg::f` form that the reader does not rewrite" = start)
  }
  if (acc$eval) {
    out <- c(out, "sets the knitr chunk option `eval` for the chunks after it, so which chunks run is unknown" = start)
  }
  out
}

# How the right-hand side of a magrittr `lhs %>% rhs` takes `lhs`: "call" for
# `f(x)` (and `pkg::f(x)`), "name" for `f`, NA for a form magrittr evaluates
# differently (a `{` block, a parenthesised function, `pkg::f` without
# parentheses, or anything else), which the reader does not rewrite.
.pipeRhsKind <- function(rhs) {
  if (is.name(rhs)) {
    return(if (identical(as.character(rhs), ".")) NA_character_ else "name")
  }
  if (is.call(rhs)) {
    head <- rhs[[1L]]
    if (is.name(head) && as.character(head) %in% c("{", "(", "::", ":::", "function")) {
      return(NA_character_)
    }
    return("call")
  }
  NA_character_
}

# Is the symbol `.` one of the direct arguments of `call`? magrittr then does not
# add the left-hand side as the first argument. `call[[i]]` is tested in place:
# bound to a variable, an empty argument (`x[, 1]`) errors when the variable is read.
.hasDirectDot <- function(call) {
  for (i in seq_along(call)[-1L]) {
    if (is.name(call[[i]]) && identical(as.character(call[[i]]), ".")) {
      return(TRUE)
    }
  }
  FALSE
}

# `lhs %>% f(x)` becomes `f(lhs, x)`, `lhs %>% f(x, .)` becomes `f(x, lhs)` (a
# direct `.` takes the place of the lhs, which is then not added as the first
# argument) and `lhs %>% f` becomes `f(lhs)`, at every depth; `lhs %>% pkg::f(x)`
# keeps the qualified function. A `%>%` whose right-hand side is another form
# (see .pipeRhsKind()) is left as it is.
.unpipe <- function(e) {
  if (!is.call(e)) {
    return(e)
  }
  # `e[[i]]` is tested in place: bound to a variable, an empty argument such as
  # the one in `x[, 1]` errors when the variable is read.
  for (i in seq_along(e)[-1L]) {
    if (is.call(e[[i]])) {
      e[[i]] <- .unpipe(e[[i]])
    }
  }
  if (identical(e[[1L]], as.name("%>%")) && length(e) == 3L) {
    rhs <- e[[3L]]
    kind <- .pipeRhsKind(rhs)
    if (identical(kind, "name")) {
      return(as.call(list(rhs, e[[2L]])))
    }
    if (identical(kind, "call")) {
      if (.hasDirectDot(rhs)) {
        for (i in seq_along(rhs)[-1L]) {
          if (is.name(rhs[[i]]) && identical(as.character(rhs[[i]]), ".")) {
            rhs[[i]] <- e[[2L]]
          }
        }
        return(rhs)
      }
      return(as.call(c(list(rhs[[1L]]), list(e[[2L]]), as.list(rhs)[-1L])))
    }
  }
  e
}

# Call `fn` on every call node of a language object, depth first. The object may
# be a call, an expression vector, a list of those (the chunks of
# .articleExprs()) or a function's formals, whose defaults hold calls. Empty
# arguments (`x[, 1]`) are unbound symbols that cannot be passed on, so they are
# skipped in their list slot.
.walkCalls <- function(expr, fn) {
  if (is.call(expr)) {
    fn(expr)
  }
  if (is.call(expr) || is.expression(expr) || is.list(expr) || is.pairlist(expr)) {
    kids <- as.list(expr)
    for (i in seq_along(kids)) {
      if (!identical(kids[i], list(quote(expr = )))) {
        .walkCalls(kids[[i]], fn)
      }
    }
  }
  invisible(NULL)
}

# The function a call is made through, with any `pkg::` or `pkg:::` prefix
# dropped, or NA when `e` is not a call or its head is not a name.
.calledName <- function(e) {
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
