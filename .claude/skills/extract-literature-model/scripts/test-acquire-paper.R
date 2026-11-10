#!/usr/bin/env Rscript

# Offline known-answer tests for acquire-paper.R.
#
# Usage (from the repository root; needs curl, pdftotext and jsonlite):
#   Rscript .claude/skills/extract-literature-model/scripts/test-acquire-paper.R
#
# Exit status 0 when every check passes, 1 otherwise (each failed check is
# listed on stderr).
#
# The ladder's source lookups (CrossRef, Unpaywall, Europe PMC, NCBI,
# publisher patterns) are stubbed, so no network is touched: CrossRef is
# made to offer `file://` URLs of local fixtures, and the real `curl`,
# `pdftotext`, validity and title checks then run on them. `quit()` is
# replaced by a condition carrying the status so `main()` can be driven
# in-process.

results <- new.env(parent = emptyenv())
results$failures <- character()
results$checks <- 0L

check <- function(cond, what) {
  results$checks <- results$checks + 1L
  if (!isTRUE(cond)) {
    results$failures <- c(results$failures, what)
  }
  invisible(NULL)
}

scriptDir <- function() {
  fileArg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fileArg) != 1L) {
    stop("run this file with Rscript so its location can be found")
  }
  dirname(normalizePath(sub("^--file=", "", fileArg)))
}

# A PDF above the script's 10 kB floor whose first page carries `title`.
makePdf <- function(path, title) {
  grDevices::pdf(path, compress = FALSE)
  graphics::plot.new()
  graphics::text(0.5, 0.9, title)
  for (i in seq_len(120L)) {
    graphics::text(
      0.5, 0.9 - i * 0.007,
      sprintf("filler line %d of the fixture body text", i),
      cex = 0.3
    )
  }
  invisible(grDevices::dev.off())
  stopifnot(file.info(path)$size >= 10000L)
  path
}

makeHtml <- function(path, padBytes = 0L) {
  body <- c(
    "<!DOCTYPE html><html><head><title>Just a moment...</title></head>",
    "<body><h1>Checking your browser before accessing the site.</h1>",
    strrep("<!-- pad -->", ceiling(padBytes / 12)),
    "</body></html>"
  )
  writeLines(body, path)
  path
}

fileUrl <- function(path) paste0("file://", normalizePath(path))

loadAcquirer <- function(scriptPath, candidates) {
  env <- new.env(parent = globalenv())
  # The script runs main() at load time unless the session is interactive.
  env$interactive <- function() TRUE
  sys.source(scriptPath, envir = env)
  env$crossrefCandidates <- function(doi, agent) candidates
  env$unpaywallCandidate <- function(doi, email, agent) character()
  env$europePmcId <- function(doi, agent) NULL
  env$ncbiPmcId <- function(doi, email, agent) NULL
  env$publisherCandidates <- function(doi) character()
  env$quit <- function(status = 0L, ...) {
    stop(structure(
      class = c("quitCalled", "condition"),
      list(message = "quit", call = NULL, status = status)
    ))
  }
  env
}

# Returns list(status = <exit status, or NA when main() raised an error>,
#              error = <error message, or NULL>).
runMain <- function(env, argv) {
  suppressMessages(
    tryCatch(
      {
        env$main(argv)
        list(status = NA_integer_, error = "main() returned without quit()")
      },
      quitCalled = function(q) list(status = q$status, error = NULL),
      error = function(e) {
        list(status = NA_integer_, error = conditionMessage(e))
      }
    )
  )
}

md5 <- function(path) unname(tools::md5sum(path))

readLog <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

newWorkDir <- function(root, name) {
  dir <- file.path(root, name)
  dir.create(dir)
  dir
}

argvFor <- function(dir, title = "propofol") {
  c(
    "--doi", "10.0000/fixture",
    "--out", file.path(dir, "lead.pdf"),
    "--title", title,
    "--log", file.path(dir, "acquire-log.json")
  )
}

runScenario <- function(name, fn) {
  tryCatch(
    fn(),
    error = function(e) {
      check(FALSE, sprintf("%s: scenario crashed: %s", name, conditionMessage(e)))
    }
  )
}

scenarioLadderWalk <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "ladder-walk")
  env <- loadAcquirer(
    scriptPath,
    c(fileUrl(fx$html), fileUrl(fx$wrong), fileUrl(fx$good))
  )
  res <- runMain(env, argvFor(dir))
  out <- file.path(dir, "lead.pdf")
  check(is.null(res$error), sprintf("walk: main raised: %s", res$error))
  check(identical(res$status, 0L), "walk: exit status is 0")
  check(
    file.exists(out) && md5(out) == md5(fx$good),
    "walk: --out holds the title-matching PDF"
  )
  check(
    identical(sort(list.files(dir)), c("acquire-log.json", "lead.pdf")),
    "walk: no download leftovers beside --out"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (!file.exists(logPath)) {
    check(FALSE, "walk: log file written")
    return(invisible(NULL))
  }
  lg <- readLog(logPath)
  check(identical(lg$result, "ok"), "walk: log result is ok")
  check(length(lg$attempts) == 3L, "walk: log records all three attempts")
  if (length(lg$attempts) == 3L) {
    reasons <- vapply(lg$attempts, function(a) a$reason, character(1L))
    check(
      grepl("not a valid PDF", reasons[1L], fixed = TRUE),
      "walk: attempt 1 (HTML challenge page) rejected as not a PDF"
    )
    check(
      grepl("title did not match", reasons[2L], fixed = TRUE),
      "walk: attempt 2 (wrong-title PDF) rejected on title"
    )
    check(
      grepl("title matched", reasons[3L], fixed = TRUE),
      "walk: attempt 3 (right PDF) accepted"
    )
    sources <- vapply(lg$attempts, function(a) a$source, character(1L))
    check(
      identical(sources, rep("crossref", 3L)),
      "walk: attempts are labelled with their ladder source, in order"
    )
    urls <- vapply(lg$attempts, function(a) a$url, character(1L))
    tried <- vapply(
      c(fx$html, fx$wrong, fx$good), fileUrl, character(1L), USE.NAMES = FALSE
    )
    check(
      identical(urls, tried),
      "walk: attempts are recorded in the order tried"
    )
  }
}

scenarioAllRejected <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "all-rejected")
  env <- loadAcquirer(
    scriptPath,
    c(fileUrl(fx$html), fileUrl(fx$bigHtml), fileUrl(fx$wrong))
  )
  res <- runMain(env, argvFor(dir))
  check(is.null(res$error), sprintf("rejected: main raised: %s", res$error))
  check(identical(res$status, 2L), "rejected: exit status is 2")
  check(
    !file.exists(file.path(dir, "lead.pdf")),
    "rejected: no file left at --out"
  )
  check(
    identical(list.files(dir), "acquire-log.json"),
    "rejected: no download leftovers in the output directory"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (file.exists(logPath)) {
    lg <- readLog(logPath)
    check(
      identical(lg$result, "all-sources-failed"),
      "rejected: log result is all-sources-failed"
    )
    check(length(lg$attempts) == 3L, "rejected: log records every attempt")
  } else {
    check(FALSE, "rejected: log file written")
  }
}

scenarioCurlFailure <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "curl-failure")
  missing <- file.path(root, "does-not-exist.pdf")
  env <- loadAcquirer(scriptPath, paste0("file://", missing))
  res <- runMain(env, argvFor(dir))
  check(is.null(res$error), sprintf("curl-fail: main raised: %s", res$error))
  check(identical(res$status, 2L), "curl-fail: exit status is 2")
  check(
    !file.exists(file.path(dir, "lead.pdf")),
    "curl-fail: no file left at --out"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (file.exists(logPath)) {
    lg <- readLog(logPath)
    check(length(lg$attempts) == 1L, "curl-fail: one attempt recorded")
    if (length(lg$attempts) == 1L) {
      check(
        isFALSE(lg$attempts[[1L]]$downloaded),
        "curl-fail: attempt recorded as not downloaded"
      )
      check(
        grepl("curl failed", lg$attempts[[1L]]$reason, fixed = TRUE),
        "curl-fail: reason names the curl failure"
      )
    }
  } else {
    check(FALSE, "curl-fail: log file written")
  }
}

scenarioErrorInsideAttempt <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "error-inside")
  env <- loadAcquirer(scriptPath, fileUrl(fx$good))
  env$titleMatches <- function(pdfPath, expectedFragment) stop("boom")
  res <- runMain(env, argvFor(dir))
  check(
    identical(res$error, "boom"),
    sprintf("error: the underlying error propagates (got %s)", res$error)
  )
  check(
    length(list.files(dir)) == 0L,
    "error: an error mid-attempt leaves nothing in the output directory"
  )
}

scenarioPreexistingFile <- function(root, fx, scriptPath) {
  # A failed run must not destroy a file the caller already had at --out.
  dir <- newWorkDir(root, "preexisting-kept")
  out <- file.path(dir, "lead.pdf")
  file.copy(fx$wrong, out)
  before <- md5(out)
  env <- loadAcquirer(scriptPath, c(fileUrl(fx$html), fileUrl(fx$bigHtml)))
  res <- runMain(env, argvFor(dir))
  check(is.null(res$error), sprintf("preexisting: main raised: %s", res$error))
  check(identical(res$status, 2L), "preexisting: exit status is 2")
  check(
    file.exists(out) && md5(out) == before,
    "preexisting: a failed run leaves the caller's file untouched"
  )

  # An accepted download replaces it.
  dir2 <- newWorkDir(root, "preexisting-replaced")
  out2 <- file.path(dir2, "lead.pdf")
  file.copy(fx$wrong, out2)
  env2 <- loadAcquirer(scriptPath, fileUrl(fx$good))
  res2 <- runMain(env2, argvFor(dir2))
  check(is.null(res2$error), sprintf("replace: main raised: %s", res2$error))
  check(identical(res2$status, 0L), "replace: exit status is 0")
  check(
    file.exists(out2) && md5(out2) == md5(fx$good),
    "replace: an accepted PDF replaces the mismatched file"
  )
}

scenarioAlreadyOnDisk <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "already-on-disk")
  out <- file.path(dir, "lead.pdf")
  file.copy(fx$good, out)
  env <- loadAcquirer(scriptPath, fileUrl(fx$html))
  res <- runMain(env, argvFor(dir))
  check(is.null(res$error), sprintf("on-disk: main raised: %s", res$error))
  check(identical(res$status, 0L), "on-disk: exit status is 0")
  check(md5(out) == md5(fx$good), "on-disk: the existing PDF is untouched")
  lg <- readLog(file.path(dir, "acquire-log.json"))
  check(
    identical(lg$result, "already-on-disk") && length(lg$attempts) == 0L,
    "on-disk: no download is attempted"
  )
}

main <- function() {
  for (tool in c("curl", "pdftotext")) {
    if (!nzchar(Sys.which(tool))) {
      stop(sprintf("%s is required on PATH to run these tests", tool))
    }
  }
  scriptPath <- file.path(scriptDir(), "acquire-paper.R")
  stopifnot(file.exists(scriptPath))

  root <- tempfile("acquire-paper-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  fixtures <- newWorkDir(root, "fixtures")
  fx <- list(
    good = makePdf(
      file.path(fixtures, "good.pdf"),
      "Population pharmacokinetics of propofol in adults"
    ),
    wrong = makePdf(
      file.path(fixtures, "wrong.pdf"),
      "Gentamicin dosing in neonates"
    ),
    html = makeHtml(file.path(fixtures, "challenge.html")),
    bigHtml = makeHtml(file.path(fixtures, "challenge-big.html"), 12000L)
  )
  stopifnot(file.info(fx$html)$size < 10000L, file.info(fx$bigHtml)$size >= 10000L)

  runScenario("ladder walk", function() scenarioLadderWalk(root, fx, scriptPath))
  runScenario("all rejected", function() scenarioAllRejected(root, fx, scriptPath))
  runScenario("curl failure", function() scenarioCurlFailure(root, fx, scriptPath))
  runScenario("error inside attempt", function() {
    scenarioErrorInsideAttempt(root, fx, scriptPath)
  })
  runScenario("pre-existing file", function() {
    scenarioPreexistingFile(root, fx, scriptPath)
  })
  runScenario("already on disk", function() {
    scenarioAlreadyOnDisk(root, fx, scriptPath)
  })

  nFail <- length(results$failures)
  message(sprintf("%d checks, %d failed", results$checks, nFail))
  if (nFail > 0L) {
    message(paste0("FAIL: ", results$failures, collapse = "\n"))
    return(1L)
  }
  0L
}

quit(status = main())
