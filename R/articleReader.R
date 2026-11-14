# Read the R code an Rmd article runs, without evaluating it and without knitr
# (knitr is only in Suggests). Used by the `nesting_assertion` rule of
# checkModelConventions(); any check that has to know what an article does reads
# it through .articleExprs().

# Spellings of FALSE that a YAML 1.1 parser (the one knitr uses for `#| eval:`)
# reads as a logical.
.yamlFalse <- c("false", "False", "FALSE", "no", "No", "NO", "n", "N", "off", "Off", "OFF")

# Does a chunk not run? `eval = FALSE` (or `F`) in the header, or `#| eval: <v>` /
# `#| eval = FALSE` in the body. An `eval` that is computed (`eval = have_data`)
# is taken to run.
.chunkHeaderSkips <- function(header) {
  grepl("(?<![A-Za-z0-9_.])eval\\s*=\\s*(FALSE|F)\\b", header, perl = TRUE)
}

.chunkBodyOptionSkips <- function(line) {
  yaml <- regmatches(line, regexec("^#\\|\\s*eval\\s*:\\s*(\\S+)\\s*$", line))[[1L]]
  if (length(yaml) == 2L) {
    return(yaml[[2L]] %in% .yamlFalse)
  }
  grepl("^#\\|\\s*eval\\s*=\\s*(FALSE|F)\\s*$", line)
}

# The running R chunks of an Rmd, parsed, as list(chunks, unparseable).
#
# `chunks` holds one list of expressions per chunk; magrittr `%>%` calls are
# rewritten into the nested call the native pipe already parses to (.unpipe()),
# so a caller sees every call the same way. A chunk is fenced by three or more
# backticks and closes on a line of at least as many; text inside a fenced block
# that is not an R chunk is not read. `unparseable` holds, for each running chunk
# that does not parse or is never closed, the line it opens on, named by the
# reason; a caller must treat it as unread code, not as an empty chunk.
.articleExprs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  chunks <- list()
  unparseable <- integer()
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
          unparseable <- c(unparseable, "does not parse" = start)
        } else {
          chunks[[length(chunks) + 1L]] <- lapply(as.list(parsed), .unpipe)
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
    unparseable <- c(unparseable, "is never closed" = start)
  }
  list(chunks = chunks, unparseable = unparseable)
}

# `lhs %>% f(x)` becomes `f(lhs, x)` and `lhs %>% f` becomes `f(lhs)`, at every
# depth; `lhs %>% pkg::f(x)` keeps the qualified function.
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
    if (is.call(rhs)) {
      return(as.call(c(list(rhs[[1L]]), list(e[[2L]]), as.list(rhs)[-1L])))
    }
    if (is.name(rhs)) {
      return(as.call(list(rhs, e[[2L]])))
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
