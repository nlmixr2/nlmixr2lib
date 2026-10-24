#!/usr/bin/env Rscript
# equals-count-probe.R <paper.pdf | paper.txt | paper.md>
# equals-count-probe.R --self-test
#
# Count the `=` signs in a paper's text layer, but only after checking that the
# text layer encodes mathematical glyphs at all. A small `=` count is used as
# corroboration that a platform PBPK / QSP paper prints few or no equations; it
# is never the gate (the Litjens volume grep is).
#
# The count fails in both directions, so it is reported only with a verdict:
#
# * False zero. Some publisher math fonts do not map `=` to U+003D, so
#   `pdftotext` emits no `=` while the display equations sit in the vector
#   layer. The same font turns operators into digits: superscript minus
#   becomes `2` (`min21` or `min 2 1` for min^-1), and `=` and `<=` become `5`
#   (`n 5 10`, `BLQ 5 10 ng/mL`). Salerno 2017 solithromycin
#   (doi:10.1002/psp4.12252) scores 0 `=` while printing three equations,
#   including a complete ODE; Li 2018 (doi:10.1002/psp4.12260) keeps `<` and
#   `>` but prints every `=` as `5`.
# * False non-zero. Other fonts render the solidus as `=`, so unit strings
#   (`lg=L`, `ml=hr`, `cigarettes=day`) inflate the count of a paper with
#   almost no math. El-Masri 2018 arsenic (PMC6108830) scores 45, 25 of them
#   unit strings.
# * Both at once. Some Springer/Adis fonts print `=` as the vulgar fraction
#   one-quarter (U+00BC) and `/` as `=`, with `(` `)` as `eth` / `thorn`
#   and the decimal point as `:` (`CL=F (L=h) <1/4> 2:1`, Altcheh 2023,
#   PMC10259795). Every real equation is missing from the count and every unit
#   is in it. ASM journals (Antimicrob Agents Chemother) print `=` as U+2AFD
#   and `+/-` as U+2AFE (Lakota 2018, PMC5971566).
#
# Verdicts (one `EQ_PROBE` line, then the evidence lines):
#   OK          the text layer carries math glyphs; the count stands.
#   LOWER_BOUND operators are mangled somewhere but `=` encodes; the count is
#               a lower bound.
#   UNUSABLE    no math glyph survives, or `=` is absent while operators are
#               mangled, or most `=` hits are mangled solidi, or `=` is
#               printed as another glyph. The count is not a finding; read
#               the equation pages off a raster
#               (`pdftoppm -r 160 -f <p> -l <p> -png`). Short documents
#               (corrections, perspectives) with no mathematics at all also
#               land here, and the raster read settles them quickly.
#
# Exit status: 0 OK, 1 LOWER_BOUND or UNUSABLE, 2 usage / environment error.
# `--self-test` runs the classifier on synthetic known-answer text and exits 0
# only if every verdict is the expected one.
#
# Usage:  Rscript equals-count-probe.R paper.pdf

# Operators other than `=` that a pharmacometric paper nearly always prints
# somewhere (mean +/- SD, P < 0.05, x 10^n, a minus sign). Their absence means
# the text layer does not encode math, whatever the `=` count says.
glyphs <- c(
  "+/-" = "\u00b1",
  "<" = "<",
  ">" = ">",
  "<=" = "\u2264",
  ">=" = "\u2265",
  "x" = "\u00d7",
  "minus" = "\u2212",
  "approx" = "\u2248"
)
# Operators mangled into digits: a unit fused to `21` / ` 2 1` (min^-1), and a
# statistic, count or limit followed by a bare `5` where `=` or `<=` belongs.
minus_re <- "\\b(min|h|hr|hrs|s|sec|d|day|days|wk|L|mL|kg|g|cm)(21|\\s2\\s1)\\b"
op5_re <- "\\b(BLQ|LLOQ|LOQ|LOD|P|p|n|N|r|R2|R)\\s+5\\s+[0-9.]"
# The solidus rendered as `=`: two short unit-like tokens fused by `=`.
slash_re <- "\\b[A-Za-z\u00b5]{1,4}=[A-Za-z]{1,11}\\b"
# `=` rendered as a free-standing one-quarter between two tokens, or as the
# ASM double-solidus operator. Two hits are required so a single literal
# one-quarter in prose does not condemn a paper.
eqsub_re <- "(?<=\\S)\\s\u00bc\\s(?=\\S)|\u2afd"

count_matches <- function(pattern, x) {
  m <- gregexpr(pattern, x, perl = TRUE)
  sum(vapply(m, function(v) sum(v > 0L), integer(1L)))
}

# URLs carry `=`, `<` and `>` that are not mathematics. pdftotext -layout
# splits a bracketed URL across lines and columns, so the opening `<` before a
# scheme and a closing `>` glued to a URL-ish character are dropped on their
# own as well as whole URLs.
strip_urls <- function(lines) {
  x <- gsub("<[^<>[:space:]]*>", " ", lines, perl = TRUE)
  x <- gsub("https?://\\S+|www\\.\\S+|doi:\\S+", " ", x, perl = TRUE)
  x <- gsub("<(?=\\s*$|\\s{2,}|https?:|www\\.)", " ", x, perl = TRUE)
  gsub("(?<=[/.A-Za-z0-9_%-])>(?=[\\s.,;)]|$)", " ", x, perl = TRUE)
}

probe_lines <- function(lines) {
  body <- strip_urls(lines)
  n_eq <- count_matches("=", body)
  glyph_n <- vapply(glyphs, function(g) count_matches(g, body), integer(1L))
  n_glyph <- sum(glyph_n)
  n_minus <- count_matches(minus_re, body)
  n_op5 <- count_matches(op5_re, body)
  n_mangled <- n_minus + n_op5
  n_slash <- count_matches(slash_re, body)
  n_eqsub <- count_matches(eqsub_re, body)
  verdict <- if (n_eqsub >= 2L) {
    "UNUSABLE"
  } else if (n_eq == 0L && n_glyph == 0L) {
    "UNUSABLE"
  } else if (n_eq == 0L && n_mangled > 0L) {
    "UNUSABLE"
  } else if (n_eq > 0L && n_slash * 2L >= n_eq) {
    "UNUSABLE"
  } else if (n_mangled > 0L) {
    "LOWER_BOUND"
  } else {
    "OK"
  }
  list(
    body = body,
    eq = n_eq,
    glyph_n = glyph_n,
    glyphs = n_glyph,
    minus = n_minus,
    op5 = n_op5,
    slash = n_slash,
    eqsub = n_eqsub,
    verdict = verdict
  )
}

hit_lines <- function(pattern, x, label, max_show = 5L) {
  idx <- grep(pattern, x, perl = TRUE)
  if (length(idx) == 0L) {
    character(0L)
  } else {
    shown <- utils::head(idx, max_show)
    c(
      sprintf("  %s: %d line(s)", label, length(idx)),
      sprintf("    %d: %s", shown, substr(trimws(x[shown]), 1L, 110L))
    )
  }
}

self_test <- function() {
  cases <- list(
    list(
      name = "no math glyph at all (Salerno-shaped false zero)",
      text = c("Clearance was 12 L/h.", "Volume was 40 L."),
      expect = "UNUSABLE"
    ),
    list(
      name = "operators mangled into digits, no `=`",
      text = c("kinact (min21) 0.022", "samples BLQ 5 10 ng/mL", "P < 0.05"),
      expect = "UNUSABLE"
    ),
    list(
      name = "`=` rendered as 5 while `<` survives (Li 2018-shaped)",
      text = c("food effect (n 5 10)", "AUCR <1.25"),
      expect = "UNUSABLE"
    ),
    list(
      name = "`=` as one-quarter, `/` as `=` (Altcheh-shaped)",
      text = c(
        "CL=F (L=h) \u00bc 2:1 (L=h) x",
        "V \u00bc 1:75 (L)",
        "P < 0.05",
        "CL = 2"
      ),
      expect = "UNUSABLE"
    ),
    list(
      name = "`=` as U+2AFD, +/- as U+2AFE (ASM / Lakota-shaped)",
      text = c("subjects (n \u2afd 42)", "CL \u2afd 3.1 \u2afe 0.4 L/h"),
      expect = "UNUSABLE"
    ),
    list(
      name = "solidus rendered as `=` (El-Masri-shaped false non-zero)",
      text = c("As conc. (lg=L)", "Rice (g=d)", "Milk (ml=d)", "CL = 2"),
      expect = "UNUSABLE"
    ),
    list(
      name = "`=` encodes but superscript minus is mangled",
      text = c("CL = 2.1 L/h", "V = 30 L", "ka (h 2 1) 1.2", "P < 0.05"),
      expect = "LOWER_BOUND"
    ),
    list(
      name = "clean text layer",
      text = c("CL = 2.1 L/h", "age 45 \u00b1 12 y", "P < 0.05"),
      expect = "OK"
    ),
    list(
      name = "clean text layer with no `=` is a true zero",
      text = c("age 45 \u00b1 12 y", "P < 0.05"),
      expect = "OK"
    ),
    list(
      name = "one literal one-quarter in prose is not a mangled `=`",
      text = c("CL = 2.1 L/h", "about \u00bc of patients", "P < 0.05"),
      expect = "OK"
    ),
    list(
      name = "URL brackets and query strings are not math",
      text = c(
        "See <https://example.org/a.php?x=1> and <http://",
        "www.example.org/b.pdf> (2010).",
        "SB Suite. <      35. Shi, J.",
        "Clearance was 12 L/h."
      ),
      expect = "UNUSABLE"
    )
  )
  got <- vapply(cases, function(k) probe_lines(k$text)$verdict, character(1L))
  want <- vapply(cases, function(k) k$expect, character(1L))
  nm <- vapply(cases, function(k) k$name, character(1L))
  cat(
    sprintf(
      "%s  %s (want %s, got %s)\n",
      ifelse(got == want, "PASS", "FAIL"),
      nm,
      want,
      got
    ),
    sep = ""
  )
  all(got == want)
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--self-test")) {
  quit(status = if (self_test()) 0L else 1L)
}
if (length(args) != 1L) {
  cat("usage: equals-count-probe.R <paper.pdf | paper.txt | paper.md>\n")
  cat("       equals-count-probe.R --self-test\n")
  quit(status = 2L)
}
f <- args[[1L]]
if (!file.exists(f)) {
  cat(sprintf("equals-count-probe: no such file: %s\n", f))
  quit(status = 2L)
}

if (grepl("\\.pdf$", f, ignore.case = TRUE)) {
  if (!nzchar(Sys.which("pdftotext"))) {
    cat("equals-count-probe: pdftotext (poppler-utils) is not on PATH\n")
    quit(status = 2L)
  }
  lines <- suppressWarnings(
    system2("pdftotext", c("-layout", shQuote(f), "-"), stdout = TRUE)
  )
} else {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
}
if (length(lines) == 0L || !any(nzchar(trimws(lines)))) {
  cat(sprintf("equals-count-probe: no text layer in %s\n", f))
  quit(status = 2L)
}

res <- probe_lines(lines)
cat(sprintf(
  paste0(
    "EQ_PROBE file=%s eq=%d glyphs=%d mangled_minus=%d mangled_eq_leq=%d ",
    "solidus_as_eq=%d eq_as_other_glyph=%d verdict=%s\n"
  ),
  basename(f),
  res$eq,
  res$glyphs,
  res$minus,
  res$op5,
  res$slash,
  res$eqsub,
  res$verdict
))
cat(sprintf(
  "  glyph inventory: %s\n",
  paste(sprintf("%s:%d", names(res$glyph_n), res$glyph_n), collapse = " ")
))
evidence <- c(
  hit_lines(minus_re, res$body, "mangled superscript minus"),
  hit_lines(op5_re, res$body, "`=` or `<=` rendered as 5"),
  hit_lines(slash_re, res$body, "solidus rendered as ="),
  hit_lines(eqsub_re, res$body, "`=` rendered as another glyph")
)
if (length(evidence) > 0L) {
  cat(evidence, sep = "\n")
}
if (res$verdict == "UNUSABLE") {
  cat(paste0(
    "  The `=` count is not a finding. Read the equation pages off a raster ",
    "(pdftoppm -r 160 -f <p> -l <p> -png) before reporting how many ",
    "equations the paper prints.\n"
  ))
}
quit(status = if (res$verdict == "OK") 0L else 1L)
