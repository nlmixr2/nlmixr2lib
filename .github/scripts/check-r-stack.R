#!/usr/bin/env Rscript
#
# Fails the job when the installed nlmixr2 stack is not mutually compatible.
#
# WHY
#
# rxode2 hands its dependents a positional vector of C function pointers --
# `rxode2:::.rxode2ptrs()` -- which `inst/include/rxode2ptr.h` consumes as
# `VECTOR_ELT(p, 0) ... VECTOR_ELT(p, N)` at FIXED indices and with NO length
# check. The vector is append-only and it grows: 92 entries in rxode2 5.1.6, 99
# in 5.1.7. A package compiled against the 99-entry header but loaded against a
# 92-entry rxode2 reads seven elements past the end of the list, converts
# whatever it finds there into a function pointer, and calls it. Nothing errors:
# the process crashes, or computes nonsense, somewhere far from the cause.
#
# CI can produce exactly that pairing. The nlmixr2 stack is served both by
# CRAN/RSPM and by nlmixr2.r-universe.dev; resolution has differed per platform
# within a single run (rxode2 5.1.7 on one job and 5.1.8 on another, nlmixr2est
# 7.0.2 against 7.1.0); and r-universe rebuilds continuously under a FIXED
# version string while the dependency cache is keyed on version -- so a restored
# "5.1.8" can be an older build than the "5.1.8" its companions were compiled
# against.
#
# This script turns that class of fault from a silent crash into a red job with
# a readable message, and prints the resolved stack either way so a human can
# see at a glance what was installed and from where.

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

STACK <- c(
  "rxode2", "rxode2ll", "lotri", "nlmixr2est", "nlmixr2data",
  "PreciseSums", "dparser", "n1qn1", "lbfgsb3c"
)

fail <- function(...) {
  msg <- paste0(...)
  cat("\n::error title=Incompatible R package stack::", msg, "\n", sep = "")
  cat("\n", msg, "\n", sep = "")
  quit(status = 1L, save = "no")
}

# ---- 1. What is actually installed, and where did it come from? -------------

cat("== Installed nlmixr2 stack ==\n")
rows <- list()
for (p in STACK) {
  if (!nzchar(system.file(package = p))) {
    rows[[p]] <- data.frame(package = p, version = "<not installed>",
                            repository = "-", built = "-",
                            stringsAsFactors = FALSE)
    next
  }
  d <- utils::packageDescription(p)
  # r-universe builds record RemoteUrl; CRAN/RSPM record Repository.
  src <- d$RemoteUrl %||% (d$Repository %||% "<local/unknown>")
  if (!is.null(d$RemoteSha) && !is.na(d$RemoteSha) && nzchar(d$RemoteSha)) {
    src <- paste0(src, "@", substr(d$RemoteSha, 1, 8))
  }
  rows[[p]] <- data.frame(
    package = p,
    version = as.character(utils::packageVersion(p)),
    repository = src,
    built = (strsplit(d$Built %||% "", ";")[[1]][1]) %||% "-",
    # The build timestamp is the only field that tells two r-universe builds
    # of the same version apart: pak records RemoteSha as the version string
    # for standard-repository installs, so "@5.1.8" above is not a commit.
    packaged = (strsplit(d$Packaged %||% "", ";")[[1]][1]) %||% "-",
    stringsAsFactors = FALSE
  )
}
tbl <- do.call(rbind, rows)
print(tbl, row.names = FALSE)
cat("\n")

# ---- 1b. Is an r-universe build stale against what r-universe serves now? ---
#
# setup-r-dependencies keys its cache on resolved versions, and r-universe
# rebuilds under a fixed version string, so a restored cache can be days older
# than the build pak resolved. The 2026-09-24 pull-request run tested the
# release legs against a pre-issue-1381 rxode2 while the freshly compiled legs
# had the current one. The workflows now key the cache on the r-universe build
# time; this check reports the comparison so a stale build is visible in the
# log, and annotates the job when one slips through.
runiverse_created <- function(pkg) {
  tryCatch({
    con <- url(sprintf("https://nlmixr2.r-universe.dev/api/packages/%s", pkg), open = "rb")
    on.exit(close(con), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    m <- regmatches(txt, regexpr("\"_created\": *\"[^\"]+\"", txt))
    if (!length(m)) return(NA_character_)
    sub("\"_created\": *\"([^\"]+)\"", "\\1", m)
  }, error = function(e) NA_character_)
}
for (p in tbl$package) {
  if (!grepl("nlmixr2/|r-universe", tbl$repository[tbl$package == p])) next
  have <- tbl$packaged[tbl$package == p]
  have_t <- suppressWarnings(as.POSIXct(have, tz = "UTC"))
  now <- runiverse_created(p)
  now_t <- suppressWarnings(as.POSIXct(now, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))
  if (is.na(have_t) || is.na(now_t)) {
    cat(sprintf("%s: installed build %s; r-universe build time unavailable (%s)\n", p, have, now %||% "no response"))
    next
  }
  lag <- as.numeric(difftime(now_t, have_t, units = "hours"))
  cat(sprintf("%s: installed build %s; r-universe current %s (%+.1f h)\n", p, have, now, lag))
  if (lag > 1) {
    cat(sprintf("::warning title=Stale r-universe build::%s installed here was built %s but nlmixr2.r-universe.dev now serves a build from %s; the dependency cache restored an older build. Check that the cache key includes the r-universe build time.\n", p, have, now))
  }
}
cat("\n")

# ---- 2. rxode2 satisfies this package's declared floor ----------------------

desc <- read.dcf("DESCRIPTION")
deps <- paste(c(
  if ("Imports" %in% colnames(desc)) desc[, "Imports"] else NULL,
  if ("Depends" %in% colnames(desc)) desc[, "Depends"] else NULL
), collapse = ",")
m <- regmatches(deps, regexpr("rxode2 *\\(>= *([0-9.]+)\\)", deps))
if (length(m) == 1L) {
  floor_ver <- gsub(".*>= *([0-9.]+).*", "\\1", m)
  have <- utils::packageVersion("rxode2")
  cat(sprintf("DESCRIPTION requires rxode2 (>= %s); installed %s\n",
              floor_ver, have))
  if (have < package_version(floor_ver)) {
    fail(sprintf(
      "rxode2 %s is older than the DESCRIPTION floor of %s. The job would ",
      have, floor_ver),
      "test against a different rxode2 than this package supports. ",
      "Check that nlmixr2.r-universe.dev is reachable and is being ",
      "consulted -- RSPM's snapshot lags it.")
  }
} else {
  cat("DESCRIPTION declares no rxode2 version floor.\n")
}

# ---- 3. rxode2's pointer vector is long enough for its own header -----------
#
# The header a dependent compiles against is rxode2's own, so if the installed
# rxode2 does not supply enough entries for the header shipped beside it, every
# dependent built against that header is already reading out of bounds.

hdr <- system.file("include", "rxode2ptr.h", package = "rxode2")
if (nzchar(hdr)) {
  lines <- readLines(hdr, warn = FALSE)
  idx <- regmatches(lines, gregexpr("VECTOR_ELT\\(p, *([0-9]+)\\)", lines))
  idx <- suppressWarnings(as.integer(gsub("\\D", "", unlist(idx))))
  idx <- idx[!is.na(idx)]
  if (length(idx)) {
    need <- max(idx) + 1L
    have_n <- length(rxode2:::.rxode2ptrs())
    cat(sprintf("rxode2 C pointer table: header consumes %d entries, ",
                need))
    cat(sprintf("installed rxode2 supplies %d\n", have_n))
    if (have_n < need) {
      fail(sprintf(
        "rxode2's pointer table has %d entries but its own header indexes %d. ",
        have_n, need),
        "Any package compiled against this header reads past the end of the ",
        "list and calls a garbage function pointer -- the crash is silent and ",
        "arbitrarily far from here. The installed rxode2 binary and its ",
        "headers are from different builds; clear the dependency cache ",
        "(bump cache-version) and reinstall.")
    }
  }
}

# ---- 4. Exercise the linkage for real --------------------------------------
#
# The checks above cannot see what version of the header a DEPENDENT was
# compiled against, so they cannot rule out the reverse skew. Loading the
# dependents runs their .onLoad pointer handoff, and a solve calls through the
# table. If the pointers are wrong this is where it dies -- in a five-second
# job that says so, instead of forty minutes into a vignette render.

cat("\nExercising the rxode2 C interface...\n")
suppressMessages({
  library(rxode2)
  if (nzchar(system.file(package = "nlmixr2est"))) library(nlmixr2est)
})

mod <- rxode2::rxode2({
  d/dt(central) <- -kel * central
  cp <- central / vc
})
ev <- rxode2::et(amt = 100, cmt = "central") |> rxode2::et(seq(0, 24, by = 4))
res <- rxode2::rxSolve(mod, ev, params = c(kel = 0.1, vc = 10),
                       returnType = "data.frame")

stopifnot(
  nrow(res) > 0L,
  all(is.finite(res$cp)),
  # C = 10 at t=0 for a 100 mg bolus into vc = 10, decaying thereafter.
  abs(res$cp[1] - 10) < 1e-8,
  res$cp[nrow(res)] < res$cp[1]
)
cat("  solve returned", nrow(res), "rows; Cc(0) =", res$cp[1],
    "as expected\n")
cat("\nStack is self-consistent.\n")
