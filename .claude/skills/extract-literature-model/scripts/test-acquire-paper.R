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
#              error = <error message, or NULL>,
#              messages = <character vector of stderr messages>).
runMain <- function(env, argv) {
  seen <- new.env(parent = emptyenv())
  seen$messages <- character()
  res <- withCallingHandlers(
    tryCatch(
      {
        env$main(argv)
        list(status = NA_integer_, error = "main() returned without quit()")
      },
      quitCalled = function(q) list(status = q$status, error = NULL),
      error = function(e) {
        list(status = NA_integer_, error = conditionMessage(e))
      }
    ),
    message = function(m) {
      seen$messages <- c(seen$messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  res$messages <- seen$messages
  res
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
    !any(grepl("left as it was", res$messages, fixed = TRUE)),
    "rejected: no kept-file warning when --out did not exist"
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
    check(
      is.null(lg$out_left_in_place),
      "rejected: a fresh --out is not reported as left in place"
    )
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
  check(is.null(res$error), sprintf("error: main raised: %s", res$error))
  check(identical(res$status, 1L), "error: an unexpected error exits 1")
  check(
    any(grepl("ERROR: boom", res$messages, fixed = TRUE)),
    "error: the underlying message is reported on stderr"
  )
  check(
    identical(list.files(dir), "acquire-log.json"),
    "error: only the log is left in the output directory"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (file.exists(logPath)) {
    lg <- readLog(logPath)
    check(
      identical(lg$result, "error") && identical(lg$error, "boom"),
      "error: the log records the error"
    )
  } else {
    check(FALSE, "error: log file written")
  }
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
  check(
    any(
      grepl(out, res$messages, fixed = TRUE) &
        grepl("did not validate", res$messages, fixed = TRUE) &
        grepl("left as it was", res$messages, fixed = TRUE)
    ),
    "preexisting: stderr names --out as unvalidated and left in place"
  )
  logPath <- file.path(dir, "acquire-log.json")
  check(
    file.exists(logPath) && isTRUE(readLog(logPath)$out_left_in_place),
    "preexisting: the log flags --out as left in place"
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

scenarioTitleCheckSkipped <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "title-skipped")
  env <- loadAcquirer(scriptPath, fileUrl(fx$wrong))
  env$titleMatches <- function(pdfPath, expectedFragment) NA
  res <- runMain(env, argvFor(dir))
  out <- file.path(dir, "lead.pdf")
  check(is.null(res$error), sprintf("skipped: main raised: %s", res$error))
  check(identical(res$status, 0L), "skipped: exit status is 0")
  check(
    file.exists(out) && md5(out) == md5(fx$wrong),
    "skipped: a valid PDF is accepted when the title cannot be checked"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (file.exists(logPath)) {
    lg <- readLog(logPath)
    check(
      length(lg$attempts) == 1L &&
        grepl("title check skipped", lg$attempts[[1L]]$reason, fixed = TRUE),
      "skipped: the log says the title check was skipped"
    )
  } else {
    check(FALSE, "skipped: log file written")
  }
}

scenarioInstallFailure <- function(root, fx, scriptPath) {
  # A directory at --out cannot be replaced by file.rename().
  dir <- newWorkDir(root, "install-failure")
  out <- file.path(dir, "lead.pdf")
  dir.create(out)
  env <- loadAcquirer(scriptPath, fileUrl(fx$good))
  res <- runMain(env, argvFor(dir))
  check(is.null(res$error), sprintf("install: main raised: %s", res$error))
  check(identical(res$status, 1L), "install: exit status is 1")
  check(
    any(grepl("ERROR", res$messages, fixed = TRUE)),
    "install: the failure is reported on stderr"
  )
  check(dir.exists(out), "install: the directory at --out is untouched")
  check(
    identical(sort(list.files(dir)), c("acquire-log.json", "lead.pdf")),
    "install: no staging file left behind"
  )
  logPath <- file.path(dir, "acquire-log.json")
  if (file.exists(logPath)) {
    lg <- readLog(logPath)
    check(
      identical(lg$result, "error") && length(lg$attempts) == 1L,
      "install: the log records the failed attempt and the error result"
    )
    if (length(lg$attempts) == 1L) {
      check(
        grepl("could not be moved onto --out", lg$attempts[[1L]]$reason, fixed = TRUE),
        "install: the attempt's reason names the failed move"
      )
    }
  } else {
    check(FALSE, "install: log file written")
  }
}

scenarioArgumentErrors <- function(root, fx, scriptPath) {
  dir <- newWorkDir(root, "arg-errors")
  out <- file.path(dir, "lead.pdf")
  cases <- list(
    "missing --doi" = c("--out", out),
    "missing --out" = c("--doi", "10.0000/fixture"),
    "unknown option" = c("--doi", "10.0000/fixture", "--out", out, "--bogus", "1"),
    "option without a value" = c("--out", out, "--doi")
  )
  for (nm in names(cases)) {
    env <- loadAcquirer(scriptPath, fileUrl(fx$good))
    res <- runMain(env, cases[[nm]])
    check(
      identical(res$status, 3L),
      sprintf("args: %s exits 3 (got %s)", nm, format(res$status))
    )
    check(
      any(grepl("ERROR", res$messages, fixed = TRUE)),
      sprintf("args: %s is reported on stderr", nm)
    )
  }
  env <- loadAcquirer(scriptPath, fileUrl(fx$good))
  env$Sys.which <- function(names) stats::setNames(rep("", length(names)), names)
  res <- runMain(env, argvFor(dir))
  check(identical(res$status, 3L), "args: curl absent from PATH exits 3")
  check(
    length(list.files(dir)) == 0L,
    "args: argument errors write neither a PDF nor a log"
  )
}

scenarioStaleStagingFile <- function(root, fx, scriptPath) {
  # A staging file left by an interrupted earlier run is overwritten by the
  # next download and never survives the run, whatever the outcome.
  dirOk <- newWorkDir(root, "stale-part-ok")
  writeLines("stale", file.path(dirOk, "lead.pdf.part"))
  env <- loadAcquirer(scriptPath, fileUrl(fx$good))
  res <- runMain(env, argvFor(dirOk))
  check(identical(res$status, 0L), "stale: a run over a stale .part succeeds")
  check(
    identical(sort(list.files(dirOk)), c("acquire-log.json", "lead.pdf")),
    "stale: the stale .part is gone after a successful run"
  )
  check(
    md5(file.path(dirOk, "lead.pdf")) == md5(fx$good),
    "stale: --out holds the downloaded PDF, not the stale bytes"
  )

  dirFail <- newWorkDir(root, "stale-part-fail")
  writeLines("stale", file.path(dirFail, "lead.pdf.part"))
  env2 <- loadAcquirer(scriptPath, fileUrl(fx$html))
  res2 <- runMain(env2, argvFor(dirFail))
  check(identical(res2$status, 2L), "stale: a rejected run over a stale .part exits 2")
  check(
    identical(list.files(dirFail), "acquire-log.json"),
    "stale: the stale .part is gone after a failed run"
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
  runScenario("title check skipped", function() {
    scenarioTitleCheckSkipped(root, fx, scriptPath)
  })
  runScenario("install failure", function() {
    scenarioInstallFailure(root, fx, scriptPath)
  })
  runScenario("argument errors", function() {
    scenarioArgumentErrors(root, fx, scriptPath)
  })
  runScenario("stale staging file", function() {
    scenarioStaleStagingFile(root, fx, scriptPath)
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
