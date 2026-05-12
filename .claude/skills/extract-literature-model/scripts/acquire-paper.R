#!/usr/bin/env Rscript

# Acquire an open-access PDF for a DOI, walking the 5-source ladder
# documented in `references/oa-acquisition.md`. Used by the
# extract-literature-model skill at Phase 1 Step 0.
#
# Usage:
#   Rscript acquire-paper.R --doi <DOI> --out <path/to/lead.pdf>
#                           [--title <expected-title-fragment>]
#                           [--retries 1]
#                           [--email <mailto-for-unpaywall>]
#                           [--log <path/to/acquire-log.json>]
#
# Exit codes:
#   0  success; PDF at --out
#   2  all 5 sources tried, none produced a valid PDF (sidecar to operator)
#   3  argument / environment error (missing curl, bad DOI, etc.)

suppressPackageStartupMessages({
  library(jsonlite, warn.conflicts = FALSE)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) b else a
}

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
    if (key == "--doi") {
      defaults$doi <- val
    } else if (key == "--out") {
      defaults$out <- val
    } else if (key == "--title") {
      defaults$title <- val
    } else if (key == "--retries") {
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
  head <- readBin(con, "raw", 4L)
  identical(rawToChar(head), "%PDF")
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

attemptAndRecord <- function(label, url, dest, agent, log, title) {
  ok <- curlDownload(url, dest, agent = agent)
  outcome <- list(source = label, url = url, downloaded = ok)
  if (!ok) {
    outcome$reason <- "curl failed (network / timeout / non-200)"
    log$attempts <<- append(log$attempts, list(outcome))
    return(FALSE)
  }
  if (!isValidPdf(dest)) {
    outcome$reason <- sprintf(
      "not a valid PDF (size=%s, not '%%PDF' head)",
      if (file.exists(dest)) file.info(dest)$size else "missing"
    )
    log$attempts <<- append(log$attempts, list(outcome))
    unlink(dest)
    return(FALSE)
  }
  tcheck <- titleMatches(dest, title)
  if (isTRUE(tcheck)) {
    outcome$reason <- "valid PDF, title matched"
    log$attempts <<- append(log$attempts, list(outcome))
    TRUE
  } else if (is.na(tcheck)) {
    outcome$reason <- "valid PDF, title check skipped (pdftotext absent)"
    log$attempts <<- append(log$attempts, list(outcome))
    TRUE
  } else {
    outcome$reason <- "valid PDF but title did not match --title expectation"
    log$attempts <<- append(log$attempts, list(outcome))
    unlink(dest)
    FALSE
  }
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
  writeLog(log, logPath)
  message("OK: ", out)
  quit(status = 0L)
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
  args <- parseArgs(argv)
  if (!nzchar(Sys.which("curl"))) {
    message("ERROR: curl not on PATH; cannot acquire PDFs.")
    quit(status = 3L)
  }
  doi <- args$doi
  out <- args$out
  agent <- sprintf("literature-acquisition (mailto:%s)", args$email)
  log <- list(doi = doi, out = out, attempts = list())

  if (isValidPdf(out) && isTRUE(titleMatches(out, args$title))) {
    log$result <- "already-on-disk"
    writeLog(log, args$log)
    message("OK: ", out, " already on disk and valid.")
    quit(status = 0L)
  }

  for (url in crossrefCandidates(doi, agent)) {
    if (attemptAndRecord("crossref", url, out, agent, log, args$title)) {
      exitOk(log, args$log, out)
    }
  }

  for (url in unpaywallCandidate(doi, args$email, agent)) {
    if (attemptAndRecord("unpaywall", url, out, agent, log, args$title)) {
      exitOk(log, args$log, out)
    }
  }

  pmcid <- europePmcId(doi, agent)
  if (!is.null(pmcid)) {
    url <- sprintf("https://europepmc.org/articles/%s?pdf=render", pmcid)
    if (attemptAndRecord("europepmc", url, out, agent, log, args$title)) {
      exitOk(log, args$log, out)
    }
  }

  if (is.null(pmcid)) {
    pmcid <- ncbiPmcId(doi, args$email, agent)
  }
  if (!is.null(pmcid)) {
    url <- sprintf(
      "https://www.ncbi.nlm.nih.gov/pmc/articles/%s/pdf/",
      pmcid
    )
    if (attemptAndRecord("ncbi-pmc", url, out, agent, log, args$title)) {
      exitOk(log, args$log, out)
    }
  }

  for (url in publisherCandidates(doi)) {
    matched <- attemptAndRecord(
      "publisher-pattern", url, out, agent, log, args$title
    )
    if (matched) {
      exitOk(log, args$log, out)
    }
  }

  log$result <- "all-sources-failed"
  writeLog(log, args$log)
  message("FAIL: all 5 sources tried; see log for details.")
  quit(status = 2L)
}

if (!interactive()) {
  argv <- commandArgs(trailingOnly = TRUE)
  main(argv)
}
