#!/usr/bin/env Rscript

# Acquire an open-access PDF for a DOI, walking the 5-source ladder
# documented in `references/oa-acquisition.md`. Used by the
# extract-literature-model skill at Phase 1 Step 0.
#
# Usage:
#   Rscript acquire-paper.R --doi <DOI> --out <path/to/lead.pdf>
#                           [--title <expected-title-fragment>]
#                           [--retries 1]   (accepted for compatibility,
#                                            currently ignored)
#                           [--email <mailto-for-unpaywall>]
#                           [--log <path/to/acquire-log.json>]
#
# Exit codes:
#   0  success; PDF at --out
#   1  unexpected error (e.g. an accepted download could not be moved onto
#      --out); the message goes to stderr and, with --log, into the log as
#      result "error"
#   2  all 5 sources tried, none produced a valid PDF (sidecar to operator)
#   3  argument / environment error (missing curl, bad DOI, etc.)
#
# Downloads are staged in `<out>.part` and moved onto --out only once they
# pass the PDF and title checks, so --out holds an accepted PDF only after
# exit 0; a file already at --out is replaced by an accepted PDF and
# otherwise left as it was. When exit 2 leaves such a file in place, stderr
# says so and the log carries `out_left_in_place: true`: it did not validate,
# so a caller must not use it.

suppressPackageStartupMessages({
  library(jsonlite, warn.conflicts = FALSE)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) b else a
}

knownOptions <- c("--doi", "--out", "--title", "--retries", "--email", "--log")

parseArgs <- function(argv) {
  defaults <- list(
    doi = NULL,
    out = NULL,
    title = NULL,
    retries = 1L,
    email = "wdenney@humanpredictions.com",
    log = NULL
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[i]
    val <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (key %in% knownOptions && is.na(val)) {
      stop(sprintf("%s requires a value", key), call. = FALSE)
    }
    if (key == "--doi") {
      defaults$doi <- val
    } else if (key == "--out") {
      defaults$out <- val
    } else if (key == "--title") {
      defaults$title <- val
    } else if (key == "--retries") {
      if (!grepl("^[0-9]+$", val)) {
        stop(
          sprintf("--retries must be a non-negative integer, got '%s'", val),
          call. = FALSE
        )
      }
      defaults$retries <- as.integer(val)
    } else if (key == "--email") {
      defaults$email <- val
    } else if (key == "--log") {
      defaults$log <- val
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
    i <- i + 2L
  }
  if (is.null(defaults$doi) || !nzchar(defaults$doi)) {
    stop("--doi is required", call. = FALSE)
  }
  if (is.null(defaults$out) || !nzchar(defaults$out)) {
    stop("--out is required", call. = FALSE)
  }
  defaults
}

isValidPdf <- function(path, minBytes = 10000L) {
  if (!file.exists(path)) return(FALSE)
  if (file.info(path)$size < minBytes) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con))
  # Compared as raw bytes: rawToChar() errors on an embedded NUL, and a short
  # read never equals the four magic bytes.
  identical(readBin(con, "raw", 4L), charToRaw("%PDF"))
}

curlDownload <- function(url, dest, agent = NULL) {
  args <- c("-sS", "-L", "--max-time", "60", "-o", shQuote(dest))
  if (!is.null(agent)) args <- c(args, "-A", shQuote(agent))
  args <- c(args, shQuote(url))
  status <- suppressWarnings(
    system2("curl", args, stdout = FALSE, stderr = FALSE)
  )
  status == 0L
}

curlJson <- function(url, agent = NULL) {
  args <- c("-sS", "-L", "--max-time", "30")
  if (!is.null(agent)) args <- c(args, "-A", shQuote(agent))
  args <- c(args, shQuote(url))
  out <- suppressWarnings(
    system2("curl", args, stdout = TRUE, stderr = FALSE)
  )
  if (length(out) == 0L) return(NULL)
  tryCatch(
    jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = TRUE),
    error = function(e) NULL
  )
}

titleMatches <- function(pdfPath, expectedFragment) {
  if (is.null(expectedFragment) || !nzchar(expectedFragment)) return(TRUE)
  if (!nzchar(Sys.which("pdftotext"))) return(NA)
  out <- suppressWarnings(
    system2(
      "pdftotext",
      c("-l", "1", shQuote(pdfPath), "-"),
      stdout = TRUE,
      stderr = FALSE
    )
  )
  body <- tolower(paste(out, collapse = " "))
  grepl(tolower(expectedFragment), body, fixed = TRUE)
}

# Downloads `url` beside `dest` and moves it onto `dest` only when it is a
# valid PDF whose first page matches `title`, so a rejected or partial
# download is never visible at `dest` and an existing `dest` survives a
# failed attempt. Returns list(log = <log with this attempt appended>,
# ok = <accepted>, error = <message when an accepted download could not be
# moved onto `dest`, else NULL>).
attemptAndRecord <- function(log, label, url, dest, agent, title) {
  part <- paste0(dest, ".part")
  # A staging file left by an interrupted run must not be mistaken for this
  # download when curl exits 0 without writing anything.
  unlink(part)
  on.exit(unlink(part), add = TRUE)
  ok <- curlDownload(url, part, agent = agent)
  outcome <- list(source = label, url = url, downloaded = ok)
  accepted <- FALSE
  error <- NULL
  if (!ok) {
    outcome$reason <- "curl failed (network / timeout / non-200)"
  } else if (!isValidPdf(part)) {
    outcome$reason <- sprintf(
      "not a valid PDF (size=%s, not '%%PDF' head)",
      if (file.exists(part)) file.info(part)$size else "missing"
    )
  } else {
    tcheck <- titleMatches(part, title)
    if (isTRUE(tcheck)) {
      outcome$reason <- "valid PDF, title matched"
      accepted <- TRUE
    } else if (is.na(tcheck)) {
      outcome$reason <- "valid PDF, title check skipped (pdftotext absent)"
      accepted <- TRUE
    } else {
      outcome$reason <- "valid PDF but title did not match --title expectation"
    }
  }
  if (accepted && !suppressWarnings(file.rename(part, dest))) {
    accepted <- FALSE
    error <- sprintf("%s passed the checks but could not be moved onto %s", part, dest)
    outcome$reason <- sprintf("%s, but could not be moved onto --out", outcome$reason)
  }
  log$attempts <- append(log$attempts, list(outcome))
  list(log = log, ok = accepted, error = error)
}

writeLog <- function(log, path) {
  if (is.null(path)) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(
    jsonlite::toJSON(log, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    path
  )
}

exitOk <- function(log, logPath, out) {
  log$result <- "ok"
  # The log is part of the contract: if it cannot be written, main()'s handler
  # exits 1 although the PDF is already at --out.
  writeLog(log, logPath)
  message("OK: ", out)
  quit(status = 0L)
}

exitError <- function(log, logPath, msg) {
  message("ERROR: ", msg)
  log$result <- "error"
  log$error <- msg
  tryCatch(
    writeLog(log, logPath),
    error = function(e) {
      message("WARNING: could not write the log: ", conditionMessage(e))
    }
  )
  quit(status = 1L)
}

exitUsage <- function(msg) {
  message("ERROR: ", msg)
  quit(status = 3L)
}

# One ladder attempt for main(): exits 0 on an accepted PDF and 1 when an
# accepted PDF cannot be installed at --out; otherwise returns the log with
# the rejected attempt appended.
tryLadderSource <- function(log, label, url, args, agent) {
  step <- attemptAndRecord(log, label, url, args$out, agent, args$title)
  if (step$ok) exitOk(step$log, args$log, args$out)
  if (!is.null(step$error)) exitError(step$log, args$log, step$error)
  step$log
}

crossrefCandidates <- function(doi, agent) {
  cr <- curlJson(sprintf("https://api.crossref.org/works/%s", doi), agent)
  if (is.null(cr) || is.null(cr$message$link)) return(character())
  links <- cr$message$link
  if (!is.data.frame(links)) return(character())
  hits <- links$URL[
    grepl("pdf", links[["content-type"]] %||% "", ignore.case = TRUE) |
      grepl("\\.pdf$", links$URL, ignore.case = TRUE)
  ]
  unique(hits)
}

unpaywallCandidate <- function(doi, email, agent) {
  url <- sprintf(
    "https://api.unpaywall.org/v2/%s?email=%s",
    doi, email
  )
  up <- curlJson(url, agent)
  if (is.null(up) || is.null(up$best_oa_location$url_for_pdf)) {
    return(character())
  }
  up$best_oa_location$url_for_pdf
}

europePmcId <- function(doi, agent) {
  url <- sprintf(
    paste0(
      "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
      "?query=DOI:%s&format=json"
    ),
    doi
  )
  pmc <- curlJson(url, agent)
  if (is.null(pmc) || is.null(pmc$resultList$result)) return(NULL)
  res <- pmc$resultList$result
  if (!is.data.frame(res) || is.null(res$pmcid)) return(NULL)
  cand <- res$pmcid[nzchar(res$pmcid %||% "")]
  if (length(cand) == 0L) NULL else cand[1L]
}

ncbiPmcId <- function(doi, email, agent) {
  url <- sprintf(
    paste0(
      "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
      "?ids=%s&format=json&tool=literature-acquisition&email=%s"
    ),
    doi, email
  )
  idconv <- curlJson(url, agent)
  if (is.null(idconv) || is.null(idconv$records)) return(NULL)
  rec <- idconv$records
  if (!is.data.frame(rec) || is.null(rec$pmcid)) return(NULL)
  cand <- rec$pmcid[nzchar(rec$pmcid %||% "")]
  if (length(cand) == 0L) NULL else cand[1L]
}

publisherCandidates <- function(doi) {
  prefix <- sub("/.*$", "", doi)
  tail <- sub("^[^/]+/", "", doi)
  hits <- character()
  if (grepl("^10\\.1186$", prefix)) {
    hits <- c(
      hits,
      sprintf("https://link.springer.com/content/pdf/%s.pdf", doi)
    )
  }
  if (grepl("^10\\.3389$", prefix)) {
    hits <- c(
      hits,
      sprintf("https://www.frontiersin.org/articles/%s/pdf", doi)
    )
  }
  if (grepl("^10\\.1038$", prefix)) {
    id_seg <- sub("^.*/", "", doi)
    hits <- c(
      hits,
      sprintf("https://www.nature.com/articles/%s.pdf", id_seg)
    )
  }
  if (grepl("^10\\.2147$", prefix)) {
    hits <- c(
      hits,
      sprintf("https://www.dovepress.com/getfile.php?fileID=%s", tail)
    )
  }
  hits
}

main <- function(argv) {
  args <- tryCatch(
    parseArgs(argv),
    error = function(e) exitUsage(conditionMessage(e))
  )
  if (!nzchar(Sys.which("curl"))) {
    exitUsage("curl not on PATH; cannot acquire PDFs.")
  }
  doi <- args$doi
  out <- args$out
  agent <- sprintf("literature-acquisition (mailto:%s)", args$email)
  log <- list(doi = doi, out = out, attempts = list())

  # Anything unexpected exits 1 with the log as far as it got. `log` is read
  # from this frame when the handler runs.
  tryCatch(
    {
      if (isValidPdf(out) && isTRUE(titleMatches(out, args$title))) {
        log$result <- "already-on-disk"
        writeLog(log, args$log)
        message("OK: ", out, " already on disk and valid.")
        quit(status = 0L)
      }

      for (url in crossrefCandidates(doi, agent)) {
        log <- tryLadderSource(log, "crossref", url, args, agent)
      }

      for (url in unpaywallCandidate(doi, args$email, agent)) {
        log <- tryLadderSource(log, "unpaywall", url, args, agent)
      }

      pmcid <- europePmcId(doi, agent)
      if (!is.null(pmcid)) {
        url <- sprintf("https://europepmc.org/articles/%s?pdf=render", pmcid)
        log <- tryLadderSource(log, "europepmc", url, args, agent)
      }

      if (is.null(pmcid)) {
        pmcid <- ncbiPmcId(doi, args$email, agent)
      }
      if (!is.null(pmcid)) {
        url <- sprintf(
          "https://www.ncbi.nlm.nih.gov/pmc/articles/%s/pdf/",
          pmcid
        )
        log <- tryLadderSource(log, "ncbi-pmc", url, args, agent)
      }

      for (url in publisherCandidates(doi)) {
        log <- tryLadderSource(log, "publisher-pattern", url, args, agent)
      }

      log$result <- "all-sources-failed"
      # Downloads never touch --out unless accepted, so a file here predates
      # the run and was not accepted by the initial check.
      if (file.exists(out)) {
        log$out_left_in_place <- TRUE
        message(sprintf(
          paste0(
            "WARNING: %s was already present and did not validate (not a ",
            "valid PDF, or its title did not match --title or could not be ",
            "checked); it was left as it was. Do not use it."
          ),
          out
        ))
      }
      writeLog(log, args$log)
      message("FAIL: all 5 sources tried; see log for details.")
      quit(status = 2L)
    },
    error = function(e) exitError(log, args$log, conditionMessage(e))
  )
}

if (!interactive()) {
  argv <- commandArgs(trailingOnly = TRUE)
  main(argv)
}
