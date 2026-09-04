# Structural checks on the naming registers themselves.
#
# checkModelConventions() validates the names a MODEL uses. This validates the
# REGISTERS that define them: a different failure surface, and one that had
# never been checked. An audit on 2026-09-01 found 8 duplicate canonicals, 12
# canonicals no model used, and 86 citations of model files that did not exist
# (42 distinct, 22 of which had never existed at any commit, i.e. filenames an
# extraction invented rather than the file it actually wrote).
#
# Ported from _scripts/audit_naming_registers.py in the ingestion repo so the
# invariants run in package CI instead of by hand.

.registerFiles <- function() {
  c("covariate-columns.md", "parameter-names.md", "compartment-names.md")
}

# Canonicals that legitimately carry no `Example models:` line. The structural
# PK parameters are universal: every model uses them, so listing examples would
# be noise rather than provenance. Kept as an explicit list (operator ruling
# 2026-09-01) so that a NEW example-less entry still fails the check.
.exampleExempt <- c(
  "lka", "lcl", "lvc", "lvp", "lvp2", "lq", "lq2", "lkel", "lvss", "lf", "lalag",
  "ka", "cl", "vc", "vp", "vp2", "q", "q2", "kel", "vss", "f", "alag",
  "propSd", "addSd", "lnSd", "logitSd", "probitSd"
)

# The `- **Type:**` vocabulary. An entry whose Type is absent drops silently out
# of the canonical name list that checkModelConventions() builds from these
# files, so the model-facing check stops recognising it -- that is how `depot`
# and `central` were briefly lost during a merge, with every other register
# check still green. Kept as a closed set so a typo is caught rather than
# silently minting a new type; a genuinely new type is one edit here.
.knownTypes <- c(
  "compartment", "continuous", "binary", "categorical", "count",
  "metabolite-suffix", "paper-named-param", "log-transformed-pk", "bare-pk"
)

.allMatches <- function(x, pat) {
  m <- gregexpr(pat, x, perl = TRUE)
  hits <- regmatches(x, m)[[1]]
  if (!length(hits)) {
    return(character())
  }
  sub(pat, "\\1", hits, perl = TRUE)
}

.parseRegister <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  sec <- NA_character_
  cur <- NULL
  for (i in seq_along(lines)) {
    ln <- lines[[i]]
    if (grepl("^## [^#]", ln)) {
      if (!is.null(cur)) out <- c(out, list(cur))
      cur <- NULL
      sec <- trimws(sub("^## ", "", ln))
    } else if (startsWith(ln, "### ")) {
      if (!is.null(cur)) out <- c(out, list(cur))
      head <- trimws(sub(" \\(.*$", "", sub("^### ", "", ln)))
      nms <- trimws(strsplit(head, ",")[[1]])
      nms <- nms[nzchar(nms)]
      if (!length(nms)) nms <- ""
      # `###` is also used for policy notes and for PATTERNS such as
      # `<tissue>_slab<n>`. Those are not canonicals: they have no example
      # models and no source token, so they would fake up two issue classes.
      pseudo <- grepl("[*<`]", head) || grepl(" -- ", head) ||
        grepl("[[:space:]]", nms[[1]])
      # A DEPRECATED entry is a tombstone pointing at the canonical that
      # replaced it. It has no example models and no current use BY DESIGN,
      # so it is exempt from those two checks but still has to parse and
      # still has to resolve its cross-references.
      deprecated <- grepl("DEPRECATED", ln, fixed = TRUE)
      cur <- list(name = nms[[1]], names = nms, section = sec, line = i,
                  pseudo = pseudo, deprecated = deprecated,
                  examples = character(), xrefs = character(),
                  hasExampleField = FALSE, inEx = FALSE,
                  hasTypeField = FALSE, type = NA_character_)
    } else if (!is.null(cur)) {
      st <- trimws(ln)
      if (startsWith(st, "- **Type:**")) {
        cur$hasTypeField <- TRUE
        # Trailing parentheticals qualify the type ("metabolite-suffix
        # (deprecated)"); the leading token is the type itself.
        cur$type <- trimws(sub("\\s*\\(.*$", "",
                              sub("^- \\*\\*Type:\\*\\*\\s*", "", st)))
      }
      if (grepl("Example models", st)) {
        cur$hasExampleField <- TRUE
        cur$inEx <- TRUE
      } else if (startsWith(st, "- **")) {
        cur$inEx <- FALSE
      }
      if (isTRUE(cur$inEx)) {
        # The list often continues onto following lines; reading only the header
        # line reports a populated entry as example-less.
        cur$examples <- c(cur$examples, .allMatches(ln, "`([A-Za-z0-9_.-]+[.]R)`"))
      }
      cur$xrefs <- c(cur$xrefs, .allMatches(ln, "\\[\\[([^]]+)\\]\\]"))
    }
  }
  if (!is.null(cur)) out <- c(out, list(cur))
  out
}

#' Check the naming registers for structural defects
#'
#' Validates `inst/references/*.md`, the registers that define canonical
#' covariate, parameter and compartment names, rather than the models that use
#' them. Complements [checkModelConventions()], which checks the reverse
#' direction.
#'
#' Checks performed:
#' \itemize{
#'   \item duplicate canonical entries (globally for covariates; per section
#'     for compartments, where a token may legitimately be both a bare
#'     compartment and a suffix);
#'   \item citations of model files that do not exist on disk;
#'   \item canonicals that no model uses;
#'   \item `[[cross-references]]` with no matching entry;
#'   \item entries with no `Example models:` line, outside an explicit
#'     exemption list of universal parameters;
#'   \item entries with no `Type:` line, or whose `Type:` is outside the
#'     known vocabulary. A missing `Type:` is the quietest defect of the set:
#'     the entry drops out of the canonical name list
#'     [checkModelConventions()] builds, so the model-facing check stops
#'     recognising the name while every other register check stays green.
#' }
#'
#' @param root Package root to check. `NULL` (default) uses the installed
#'   package.
#' @return A data.frame with one row per issue and columns `register`,
#'   `check`, `name`, `line`, `detail`. Zero rows means clean.
#' @export
checkNamingRegisters <- function(root = NULL) {
  if (is.null(root)) {
    modeldbDir <- system.file("modeldb", package = "nlmixr2lib")
    refDir <- system.file("references", package = "nlmixr2lib")
  } else {
    modeldbDir <- file.path(root, "inst", "modeldb")
    refDir <- file.path(root, "inst", "references")
  }
  empty <- data.frame(register = character(), check = character(),
                      name = character(), line = integer(),
                      detail = character(), stringsAsFactors = FALSE)
  if (!dir.exists(modeldbDir) || !dir.exists(refDir)) {
    return(empty)
  }
  modelFiles <- list.files(modeldbDir, pattern = "[.]R$", recursive = TRUE)
  onDisk <- basename(modelFiles)
  srcTokens <- unique(unlist(lapply(file.path(modeldbDir, modelFiles), function(f) {
    .allMatches(paste(readLines(f, warn = FALSE), collapse = "\n"),
                "([A-Za-z_][A-Za-z0-9_.]*)")
  })))

  # Accumulate into an explicit environment rather than using `<<-` from the
  # inner add() closure: environments have reference semantics, so a plain
  # `<-` into `acc` mutates the shared accumulator, and the assignment target
  # is named at the assignment site. Matches `.parseCovariateColumns()`;
  # lintr's assignment_linter rejects `<<-`.
  acc <- new.env(parent = emptyenv())
  acc$iss <- list()
  add <- function(register, check, name, line, detail) {
    acc$iss[[length(acc$iss) + 1L]] <- data.frame(
      register = register, check = check, name = name, line = line,
      detail = detail, stringsAsFactors = FALSE)
  }

  for (f in .registerFiles()) {
    path <- file.path(refDir, f)
    if (!file.exists(path)) next
    ents <- .parseRegister(path)
    real <- Filter(function(e) !isTRUE(e$pseudo), ents)
    if (!length(real)) next
    canon <- unique(unlist(lapply(ents, `[[`, "names")))

    # Covariate canonicals are globally unique. A COMPARTMENT token may
    # legitimately be both a bare compartment and a suffix (`igg`, `plasma`,
    # `anc`), so those are checked for uniqueness per `##` section only; a
    # global check reports 10 false positives on compartment-names.md.
    glob <- identical(f, "covariate-columns.md")
    keys <- unlist(lapply(real, function(e) {
      if (glob) e$names else paste0(e$section, "\r", e$names)
    }))
    lns <- unlist(lapply(real, function(e) rep(e$line, length(e$names))))
    nms <- unlist(lapply(real, `[[`, "names"))
    for (k in unique(keys[duplicated(keys)])) {
      sel <- which(keys == k)
      add(f, "duplicate-canonical", nms[sel][[1]], lns[sel][[1]],
          paste("also at line", paste(lns[sel][-1], collapse = ", ")))
    }

    for (e in real) {
      miss <- setdiff(e$examples, onDisk)
      if (length(miss)) {
        add(f, "orphan-example-model", e$name, e$line,
            paste("cites absent file(s):", paste(miss, collapse = ", ")))
      }
      tok <- e$name
      # Suffix sections register the SUFFIX (`dox`), which appears in model
      # source only inside a compound token (`auc_dox`), so match that too.
      used <- tok %in% srcTokens ||
        any(endsWith(srcTokens, paste0("_", tok))) ||
        any(startsWith(srcTokens, paste0(tok, "_")))
      if (nchar(tok) >= 3L && !used && !isTRUE(e$deprecated)) {
        add(f, "registered-but-unused", tok, e$line, "no model on disk uses it")
      }
      if (!isTRUE(e$hasExampleField) && !(tok %in% .exampleExempt) &&
          !isTRUE(e$deprecated)) {
        add(f, "no-example-model", tok, e$line, "no `Example models:` line")
      }
      # A DEPRECATED tombstone is exempt: it is a pointer to the replacement
      # canonical, not a declaration of a name that is still in use.
      if (!isTRUE(e$hasTypeField) && !isTRUE(e$deprecated)) {
        add(f, "no-type", tok, e$line, "no `Type:` line")
      } else if (!is.na(e$type) && nzchar(e$type) &&
                 !(e$type %in% .knownTypes)) {
        add(f, "unknown-type", tok, e$line,
            paste0("Type `", e$type, "` is not one of: ",
                   paste(.knownTypes, collapse = ", ")))
      }
      for (x in unique(e$xrefs)) {
        t <- gsub('^["\']+|["\']+$', "", trimws(x))
        # Prose placeholders such as [[<name>]] or [[...]] are not links.
        if (!nzchar(t) || grepl("[<.]", t)) next
        if (!(t %in% canon)) {
          add(f, "broken-xref", e$name, e$line,
              paste0("[[", t, "]] has no entry"))
        }
      }
    }
  }
  if (!length(acc$iss)) {
    return(empty)
  }
  res <- do.call(rbind, acc$iss)
  rownames(res) <- NULL
  res
}
