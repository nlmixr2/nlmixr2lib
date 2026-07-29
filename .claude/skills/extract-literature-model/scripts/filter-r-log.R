#!/usr/bin/env Rscript
# filter-r-log.R <logfile>
#
# Compact a verbose R log (devtools::check / R CMD check / vignette render)
# down to the lines that matter: the status summary, ERROR / WARNING / NOTE
# blocks (with their indented detail), render tracebacks, and the final result
# line. Drops the hundreds of "* checking ... OK" and pandoc-verbose lines so
# the agent reads tens of lines instead of thousands.
#
# Token economics: a full `devtools::check()` log is commonly 1,000-5,000 lines;
# the agent re-reads it (cached) every turn of the build-fix loop. Reading the
# filtered view instead is the single cheapest per-task token win. Open the raw
# log ONLY if this summary is insufficient.
#
# Usage:  Rscript filter-r-log.R /tmp/vignette-render-Foo_2024_drug.log

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  cat("usage: filter-r-log.R <logfile>\n"); quit(status = 2L)
}
f <- args[[1L]]
if (!file.exists(f)) {
  cat(sprintf("filter-r-log: no such file: %s\n", f)); quit(status = 2L)
}

lines <- readLines(f, warn = FALSE)
n <- length(lines)
if (n == 0L) { cat("(empty log)\n"); quit(status = 0L) }

# Lines worth keeping. "* checking ... OK" is pure noise and dropped below.
keep_re <- paste0(
  "(ERROR|WARNING|NOTE|Error|Warning|Quitting|Status:|FAIL|failed|Failed|",
  "cannot|unable|Execution halted|halt|exit=|html_bytes|RENDER_GATE|",
  "Requesting an AUC|traceback|Backtrace|\\bin chunk\\b|",
  "\\* checking .*\\.\\.\\. (WARNING|ERROR|NOTE))"
)
ok_noise_re <- "^\\* checking .*\\.\\.\\. OK\\s*$"

keep <- logical(n)
context_after <- 6L  # keep the indented detail block under a status hit
i <- 1L
while (i <= n) {
  ln <- lines[[i]]
  if (!grepl(ok_noise_re, ln) && grepl(keep_re, ln, perl = TRUE)) {
    keep[[i]] <- TRUE
    j <- i + 1L; c <- 0L
    while (j <= n && c < context_after &&
           (grepl("^\\s", lines[[j]]) ||
            nchar(trimws(lines[[j]])) == 0L ||
            grepl(keep_re, lines[[j]], perl = TRUE))) {
      keep[[j]] <- TRUE; j <- j + 1L; c <- c + 1L
    }
  }
  i <- i + 1L
}
# Always keep the final lines (overall status / pandoc result / RENDER_GATE).
keep[seq.int(max(1L, n - 7L), n)] <- TRUE

sel <- which(keep)
out <- lines[sel]

cap <- 250L
if (length(out) > cap) {
  out <- c(out[seq_len(cap)],
           sprintf("...[filtered view truncated at %d of %d kept lines; raw log: %s]",
                   cap, length(out), f))
}

cat(sprintf("# filter-r-log: %d raw lines -> %d kept (raw log: %s)\n", n, length(sel), f))
cat(out, sep = "\n")
cat("\n")
