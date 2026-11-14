# .articleExprs() and the call walkers: what an article runs, read without
# evaluating it and without knitr.

.readerArticle <- function(lines, env = parent.frame()) {
  path <- tempfile(fileext = ".Rmd")
  withr::defer(unlink(path), envir = env)
  writeLines(lines, path)
  path
}

# The calls of every running chunk, deparsed, in order.
.readerCalls <- function(read) {
  acc <- new.env(parent = emptyenv())
  acc$calls <- character()
  .walkCalls(read$chunks, function(e) acc$calls <- c(acc$calls, paste(deparse(e), collapse = " ")))
  acc$calls
}

test_that("the running chunks come back parsed, one list per chunk, with nothing unreadable", {
  path <- .readerArticle(c(
    "# Title",
    "",
    "```{r setup, include = FALSE}",
    "library(a)",
    "x <- 1",
    "```",
    "",
    "Some prose with `inline` code.",
    "",
    "```{r second}",
    "y <- f(x)",
    "```"
  ))
  read <- .articleExprs(path)
  expect_named(read, c("chunks", "unreadable"))
  expect_length(read$chunks, 2L)
  expect_length(read$chunks[[1L]], 2L)
  expect_identical(read$chunks[[2L]][[1L]], quote(y <- f(x)))
  expect_length(read$unreadable, 0L)
})

test_that("a chunk that does not run is not read, however eval is written", {
  for (skip in c("eval = FALSE", "eval=F")) {
    path <- .readerArticle(c(sprintf("```{r a, %s}", skip), "hidden()", "```", "```{r b}", "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = skip)
  }
  for (opt in c(paste("#| eval:", .yamlFalse), "#| eval = FALSE")) {
    path <- .readerArticle(c("```{r a}", opt, "hidden()", "```", "```{r b}", "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = opt)
  }
  # A computed eval, and an eval that is not false, run.
  for (header in c("```{r a, eval = have_data}", "```{r a, eval = TRUE}", "```{r a, cache.eval = FALSE}")) {
    path <- .readerArticle(c(header, "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = header)
  }
  for (opt in c("#| eval: true", "#| eval: yes", "#| eval: y", "#| eval: on", "#| echo: false")) {
    path <- .readerArticle(c("```{r a}", opt, "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = opt)
  }
})

test_that("an eval given as a YAML !expr tag is read as the R expression it holds", {
  for (opt in c("#| eval: !expr FALSE", "#| eval: !expr F", "#|   eval:   !expr   FALSE  ")) {
    path <- .readerArticle(c("```{r a}", opt, "hidden()", "```", "```{r b}", "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = opt)
  }
  # TRUE, a computed expression, and the YAML word "false" after the tag run.
  for (opt in c("#| eval: !expr TRUE", "#| eval: !expr have_data", "#| eval: !expr false")) {
    path <- .readerArticle(c("```{r a}", opt, "shown()", "```"))
    expect_identical(.readerCalls(.articleExprs(path)), "shown()", label = opt)
  }
})

test_that("the YAML spellings of FALSE are the ones the yaml package reads as FALSE", {
  skip_if_not_installed("yaml")
  spellings <- c(
    "false",
    "False",
    "FALSE",
    "no",
    "No",
    "NO",
    "n",
    "N",
    "off",
    "Off",
    "OFF",
    "true",
    "True",
    "TRUE",
    "yes",
    "Yes",
    "YES",
    "y",
    "Y",
    "on",
    "On",
    "ON",
    "f",
    "F",
    "0",
    "nope",
    "fALSE"
  )
  isFalse <- vapply(spellings, function(sp) identical(yaml::yaml.load(paste0("eval: ", sp))$eval, FALSE), logical(1L))
  expect_setequal(spellings[isFalse], .yamlFalse)
})

test_that("chunks fenced with more than three backticks are read and close on a line of as many", {
  path <- .readerArticle(c(
    "````{r four}",
    "a()",
    "x <- '",
    "```",
    "'",
    "b()",
    "````",
    "```{r three}",
    "c()",
    "```"
  ))
  read <- .articleExprs(path)
  expect_length(read$unreadable, 0L)
  # The string literal that holds a three-backtick line is one call on its own.
  calls <- .readerCalls(read)
  expect_length(calls, 4L)
  expect_identical(calls[c(1L, 3L, 4L)], c("a()", "b()", "c()"))
})

test_that("text in a fenced block that is not an R chunk is not read", {
  path <- .readerArticle(c(
    "````markdown",
    "```{r}",
    "shown_in_the_text()",
    "```",
    "````",
    "```python",
    "not_r()",
    "```",
    "```{r real}",
    "real()",
    "```"
  ))
  expect_identical(.readerCalls(.articleExprs(path)), "real()")
})

test_that("a chunk that does not parse, or is never closed, is reported by its opening line and reason", {
  path <- .readerArticle(c("```{r ok}", "ok()", "```", "", "```{r bad}", "x <- (1 +", "```", "```{r open}", "never()"))
  read <- .articleExprs(path)
  expect_identical(.readerCalls(read), "ok()")
  expect_identical(read$unreadable, c("does not parse" = 5L, "is never closed" = 8L))
  # A chunk that does not run is not code, closed or not.
  path <- .readerArticle(c("```{r bad, eval = FALSE}", "x <- (1 +", "```", "```{r open, eval = FALSE}", "never()"))
  read <- .articleExprs(path)
  expect_length(read$chunks, 0L)
  expect_length(read$unreadable, 0L)
})

test_that("an article with no chunks reads as empty, not as an error", {
  read <- .articleExprs(.readerArticle(c("# Only prose", "", "text")))
  expect_identical(read$chunks, list())
  expect_length(read$unreadable, 0L)
})

test_that(".unpipe() rewrites magrittr pipes at every depth and leaves everything else alone", {
  expect_identical(.unpipe(quote(x %>% f(a, b = 2))), quote(f(x, a, b = 2)))
  expect_identical(.unpipe(quote(x %>% f)), quote(f(x)))
  expect_identical(.unpipe(quote(x %>% pkg::f(a))), quote(pkg::f(x, a)))
  expect_identical(.unpipe(quote(x %>% pkg:::f(a))), quote(pkg:::f(x, a)))
  expect_identical(.unpipe(quote(x %>% f() %>% g(1))), quote(g(f(x), 1)))
  expect_identical(.unpipe(quote(h(x %>% f(a), y))), quote(h(f(x, a), y)))
  expect_identical(.unpipe(quote(if (a) x %>% f(1) else y)), quote(if (a) f(x, 1) else y))
  # The native pipe is already a nested call; an empty argument and non-calls survive.
  expect_identical(.unpipe(quote(x |> f(a))), quote(f(x, a)))
  expect_identical(.unpipe(quote(m[, 1])), quote(m[, 1]))
  expect_identical(.unpipe(quote(x)), quote(x))
  expect_identical(.unpipe(1L), 1L)
  expect_identical(.unpipe(quote(x %in% y)), quote(x %in% y))
})

test_that(".unpipe() puts the lhs where a direct `.` is, and does not add it as the first argument then", {
  expect_identical(.unpipe(quote(x %>% f(y, .))), quote(f(y, x)))
  expect_identical(.unpipe(quote(x %>% f(., y, .))), quote(f(x, y, x)))
  expect_identical(.unpipe(quote(x %>% f(.))), quote(f(x)))
  expect_identical(.unpipe(quote(x %>% .$a)), quote(x$a))
  expect_identical(.unpipe(quote(x %>% pkg::f(y, .))), quote(pkg::f(y, x)))
  # A `.` that is not a direct argument does not stop the lhs being added.
  expect_identical(.unpipe(quote(x %>% f(g(.)))), quote(f(x, g(.))))
  expect_identical(.unpipe(quote(x %>% f(y = .z))), quote(f(x, y = .z)))
  expect_identical(.unpipe(quote(x %>% f(y, .) %>% h(2))), quote(h(f(y, x), 2)))
  # An empty argument next to the dot.
  expect_identical(.unpipe(quote(x %>% m[, .])), quote(m[, x]))
})

test_that(".unpipe() gives the value magrittr gives, for every form it rewrites", {
  # magrittr is not a dependency of the package: it is the reference when present.
  skip_if_not_installed("magrittr")
  `%>%` <- getExportedValue("magrittr", "%>%")
  env <- new.env(parent = environment())
  env$f <- function(...) list(...)
  env$g <- function(v) v * 10
  env$h <- function(v, w = 1) v + w
  env$x <- 3
  # magrittr binds `.` to the lhs, also inside a nested call; the rewrite leaves
  # that `.` in place, so it is bound here to the same value.
  env$. <- env$x
  env$y <- 7
  env$lst <- list(a = 5)
  env$m <- matrix(1:4, 2)
  forms <- list(
    quote(x %>% f(y)),
    quote(x %>% f(y, .)),
    quote(x %>% f(., y, .)),
    quote(x %>% f(g(.))),
    quote(x %>% g),
    quote(x %>% h(1, w = .)),
    quote(x %>% h(1, w = .)),
    quote(lst %>% .$a),
    quote(x %>% f(y) %>% length()),
    quote(x %>% h(y, .) %>% g()),
    quote(h(x %>% g, y)),
    quote(m %>% .[1, 2])
  )
  for (form in forms) {
    expect_identical(eval(.unpipe(form), env), eval(form, env), label = paste(deparse(form), collapse = " "))
  }
})

test_that("a chunk that pipes into a form the reader does not rewrite is reported, and still returned", {
  for (code in c("y <- x %>% pkg::f", "y <- x %>% {f(.)}", "y <- x %>% (function(d) f(d))", "z <- x %>% h() %>% {k(.)}")) {
    path <- .readerArticle(c("```{r ok}", "ok()", "```", "```{r p}", code, "```"))
    read <- .articleExprs(path)
    expect_identical(
      read$unreadable,
      c("pipes into a `{`, `(` or `pkg::f` form that the reader does not rewrite" = 4L),
      label = code
    )
    expect_length(read$chunks, 2L)
  }
  # The forms it does rewrite raise nothing.
  path <- .readerArticle(c("```{r p}", "y <- x %>% f(a, .) %>% g %>% pkg::h(1)", "```"))
  expect_length(.articleExprs(path)$unreadable, 0L)
})

test_that("a chunk that sets the knitr option eval makes which chunks run unknown, and is reported", {
  sets <- c(
    "knitr::opts_chunk$set(eval = FALSE)",
    "opts_chunk$set(eval = TRUE)",
    "knitr::opts_chunk$set(collapse = TRUE, eval = have_data)",
    "knitr:::opts_chunk$set(eval = FALSE)",
    "knitr::opts_chunk$set(list(echo = TRUE, eval = FALSE))"
  )
  for (code in sets) {
    path <- .readerArticle(c("```{r setup, include = FALSE}", code, "```", "```{r b}", "shown()", "```"))
    read <- .articleExprs(path)
    expect_identical(
      read$unreadable,
      c("sets the knitr chunk option `eval` for the chunks after it, so which chunks run is unknown" = 1L),
      label = code
    )
    expect_identical(.readerCalls(read)[[length(.readerCalls(read))]], "shown()", label = code)
  }
  # Other options, other functions and other objects named `set` do not.
  quiet <- c(
    "knitr::opts_chunk$set(echo = FALSE, collapse = TRUE)",
    "knitr::opts_chunk$set(list(echo = FALSE))",
    "knitr::opts_knit$set(eval = FALSE)",
    "other$set(eval = FALSE)",
    "set(eval = FALSE)",
    "opts_chunk$get(\"eval\")"
  )
  for (code in quiet) {
    path <- .readerArticle(c("```{r setup}", code, "```"))
    expect_length(.articleExprs(path)$unreadable, 0L)
  }
  # Inside a chunk that does not run, it sets nothing.
  path <- .readerArticle(c("```{r setup, eval = FALSE}", "knitr::opts_chunk$set(eval = FALSE)", "```"))
  expect_length(.articleExprs(path)$unreadable, 0L)
})

test_that(".walkCalls() visits every call of a call, an expression, a list and a function's defaults", {
  seen <- function(obj) {
    acc <- new.env(parent = emptyenv())
    acc$calls <- character()
    .walkCalls(obj, function(e) acc$calls <- c(acc$calls, paste(deparse(e), collapse = " ")))
    acc$calls
  }
  expect_identical(seen(quote(f(g(1), 2))), c("f(g(1), 2)", "g(1)"))
  expect_identical(seen(expression(a(), b(c()))), c("a()", "b(c())", "c()"))
  expect_identical(seen(list(quote(a()), list(quote(b())))), c("a()", "b()"))
  # A default argument of a function the article defines is code too.
  fn <- quote(function(cmt = h("depot"), n) n + 1)
  expect_identical(seen(fn), c(paste(deparse(fn), collapse = " "), "h(\"depot\")", "n + 1"))
  # An empty argument cannot be passed on and is skipped.
  expect_identical(seen(quote(m[, 1])), "m[, 1]")
  expect_identical(seen(NULL), character())
  expect_identical(seen(quote(x)), character())
})

test_that(".calledName() drops a package prefix and is NA for what is not a named call", {
  expect_identical(.calledName(quote(f(x))), "f")
  expect_identical(.calledName(quote(pkg::f(x))), "f")
  expect_identical(.calledName(quote(pkg:::f(x))), "f")
  expect_identical(.calledName(quote(x$f(1))), NA_character_)
  expect_identical(.calledName(quote((function(x) x)(1))), NA_character_)
  expect_identical(.calledName(quote(x)), NA_character_)
  expect_identical(.calledName(1), NA_character_)
  expect_identical(.calledName(NULL), NA_character_)
})
