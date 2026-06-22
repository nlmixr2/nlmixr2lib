#!/usr/bin/env Rscript
# lookup-canonical.R <term> [covariate|parameter|compartment]
#
# Return the canonical-register entries matching <term> (matched against the
# heading name, source aliases, and description) WITHOUT loading the whole
# register into the agent's context.
#
# Token economics: inst/references/covariate-columns.md is ~1.1 MB (~284k
# tokens). A single accidental whole-file Read blows the per-task budget.
# This script reads the register on disk and returns only the matching H2/H3
# entries (a few lines), so the 284k-token file never enters context. Use it
# for every canonical name check (covariate / parameter / compartment).
#
# Run from the repo root (the agent's worktree CWD). Exit 0 = match(es) printed;
# exit 0 with a "no match" note = likely a NEW canonical (follow the
# stop-and-ask rule in SKILL.md).
#
# Usage:  Rscript .../scripts/lookup-canonical.R WT covariate
#         Rscript .../scripts/lookup-canonical.R vc parameter

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  cat("usage: lookup-canonical.R <term> [covariate|parameter|compartment]\n")
  quit(status = 2L)
}
term <- args[[1L]]
kind <- if (length(args) >= 2L) tolower(args[[2L]]) else "all"

skill <- ".claude/skills/extract-literature-model"
cand <- list(
  covariate   = c("inst/references/covariate-columns.md",
                  file.path(skill, "references/covariate-columns.md")),
  parameter   = file.path(skill, "references/parameter-names.md"),
  compartment = file.path(skill, "references/compartment-names.md")
)
if (kind != "all") {
  if (!kind %in% names(cand)) {
    cat(sprintf("unknown kind '%s' (use covariate|parameter|compartment)\n", kind))
    quit(status = 2L)
  }
  cand <- cand[kind]
}

resolve <- function(paths) { for (p in paths) if (file.exists(p)) return(p); NA_character_ }
# whole-word, case-insensitive match (so "WT" hits "### WT" but not "BWT"/"WTKG").
# \Q...\E quotes the term literally, so any regex metacharacters are harmless.
wb <- paste0("\\b\\Q", term, "\\E\\b")
wmatch <- function(x) grepl(wb, x, ignore.case = TRUE, perl = TRUE)

entries <- list()
for (k in names(cand)) {
  f <- resolve(cand[[k]])
  if (is.na(f)) next
  lines <- readLines(f, warn = FALSE)
  head_idx <- grep("^#{2,4}\\s", lines)
  if (length(head_idx) == 0L) next
  starts <- head_idx; ends <- c(head_idx[-1L] - 1L, length(lines))
  for (b in seq_along(starts)) {
    blk <- lines[starts[[b]]:ends[[b]]]
    heading <- sub("^#+\\s*", "", lines[[starts[[b]]]])
    hmatch <- wmatch(heading)
    bmatch <- any(wmatch(blk))
    if (hmatch || bmatch)
      entries[[length(entries) + 1L]] <-
        list(k = k, f = f, heading = heading, blk = blk, hmatch = hmatch)
  }
}

# Heading matches are what "is there a canonical named X?" wants; if any exist,
# drop body-only (alias/description) mentions. Cap entries + per-entry lines.
hm <- Filter(function(e) e$hmatch, entries)
use <- if (length(hm) > 0L) hm else entries
cap_entries <- 8L
dropped <- max(0L, length(use) - cap_entries)
use <- utils::head(use, cap_entries)

if (length(use) == 0L) {
  cat(sprintf("No %s register entry matches '%s'. If the concept is genuinely absent, this is a NEW canonical -> follow the stop-and-ask rule before adding it.\n",
              if (kind == "all") "" else kind, term))
} else {
  for (e in use) {
    tag <- if (e$hmatch) "" else " [body/alias match]"
    cat(sprintf("=== [%s] %s%s  (%s)\n", e$k, e$heading, tag, e$f))
    blk2 <- if (length(e$blk) > 40L) c(e$blk[seq_len(40L)], "...[entry truncated]") else e$blk
    cat(blk2, sep = "\n"); cat("\n\n")
  }
  if (dropped > 0L)
    cat(sprintf("...[%d more matching entries omitted; refine the term]\n", dropped))
}
