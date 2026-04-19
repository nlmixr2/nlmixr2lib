#' nlmixr2lib convention standards
#'
#' Internal register of canonical parameter, compartment, covariate, and
#' residual-error names used by [checkModelConventions()]. The canonical
#' covariate list is parsed at runtime from
#' `inst/references/covariate-columns.md` (installed as
#' `system.file("references", "covariate-columns.md", package = "nlmixr2lib")`),
#' so that register remains the single authoritative source. Remaining fields
#' mirror the `extract-literature-model` skill's `naming-conventions.md` and
#' `vignettes/create-model-library.Rmd`.
#'
#' @keywords internal
#' @noRd
.nlmixr2libConventionsStatic <- list(
  pkParams = c(
    "lka", "lcl", "lvc", "lvp", "lvp2", "lq", "lq2", "lfdepot"
  ),
  pkBareParams = c(
    "ka", "cl", "vc", "vp", "vp2", "q", "q2", "kel",
    "k12", "k21", "k13", "k31", "fdepot"
  ),
  compartments = c(
    "depot", "central", "peripheral1", "peripheral2", "effect",
    "target", "complex", "total_target"
  ),
  compartmentRegex = "^(transit|effect)[0-9]+$",
  observationVar = "Cc",
  residualError = c("propSd", "addSd"),
  transformPrefixes = c("l", "logit", "probit"),
  covEffectPattern = "^e_[A-Za-z0-9]+_[A-Za-z0-9]+$",
  requiredUnits = c("time", "dosing", "concentration"),
  requiredMetadata = c("description", "reference", "units"),
  deprecatedResidualError = c(
    "prop.err", "add.err", "propErr", "addErr",
    "err.prop", "err.add"
  ),
  deprecatedIivPrefixes = c("iiv_", "IIV_", "bsv_", "BSV_")
)

.covariateRegisterCache <- new.env(parent = emptyenv())

.covariateColumnsPath <- function() {
  p <- system.file("references", "covariate-columns.md",
                   package = "nlmixr2lib")
  if (nzchar(p)) return(p)
  stop("Could not locate inst/references/covariate-columns.md in the ",
       "nlmixr2lib package. Check the installation.", call. = FALSE)
}

#' Parse the canonical covariate register from covariate-columns.md.
#'
#' Walks the Markdown register and extracts one entry per H3 heading
#' (`### NAME (**...**)`). For each entry, captures the `Units`, `Type`, and
#' the names of any `Source aliases` given as backticked identifiers. Aliases
#' whose backticked content is not a bare R identifier (e.g. `DVID = "study1"`)
#' are skipped.
#'
#' @param path Path to the markdown file.
#' @return A named list keyed by canonical name.
#' @keywords internal
#' @noRd
.parseCovariateColumns <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- list()
  current <- NULL
  state <- "idle"

  flush <- function() {
    if (is.null(current)) return(invisible())
    for (nm in current$names) {
      entries[[nm]] <<- list(
        units = current$units %||% "",
        type = current$type %||% "",
        aliases = current$aliases %||% character()
      )
    }
  }

  aliasRegex <- "^\\s*-\\s*`([^`]+)`"
  identRegex <- "^[A-Za-z_][A-Za-z0-9_]*$"

  for (line in lines) {
    if (startsWith(line, "## ") && !startsWith(line, "### ")) {
      flush()
      current <- NULL
      state <- "idle"
      next
    }
    if (startsWith(line, "### ")) {
      flush()
      heading <- sub("^###\\s+", "", line)
      heading <- sub("\\s*\\(\\*\\*.*\\*\\*\\)\\s*$", "", heading)
      nms <- trimws(strsplit(heading, ",")[[1]])
      nms <- nms[grepl(identRegex, nms)]
      current <- list(names = nms, aliases = character())
      state <- "header"
      next
    }
    if (is.null(current)) next

    m <- regmatches(line, regexec("^- \\*\\*Units:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$units <- trimws(m[[2]])
      state <- "header"
      next
    }
    m <- regmatches(line, regexec("^- \\*\\*Type:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$type <- trimws(m[[2]])
      state <- "header"
      next
    }
    if (grepl("^- \\*\\*Source aliases:\\*\\*", line)) {
      state <- "aliases"
      after <- sub("^- \\*\\*Source aliases:\\*\\*\\s*", "", line)
      # "none", "none known", "none;"- style declarations have no aliases.
      if (grepl("^none\\b", after, ignore.case = TRUE)) next
      # Capture inline aliases up to the first em-dash prose separator.
      after <- strsplit(after, "\\s+\u2014\\s+", perl = TRUE)[[1]][1]
      inline <- regmatches(after, gregexpr("`([^`]+)`", after))[[1]]
      for (tok in inline) {
        inner <- gsub("`", "", tok)
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
      }
      next
    }

    if (state == "aliases") {
      m <- regmatches(line, regexec(aliasRegex, line))[[1]]
      if (length(m) == 2) {
        inner <- m[[2]]
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
        next
      }
      if (grepl("^- \\*\\*", line)) {
        state <- "header"
      }
    }
  }
  flush()
  entries
}

.loadCanonicalCovariates <- function(force = FALSE) {
  if (!force && !is.null(.covariateRegisterCache$canonical)) {
    return(.covariateRegisterCache$canonical)
  }
  entries <- .parseCovariateColumns(.covariateColumnsPath())
  .covariateRegisterCache$canonical <- entries
  entries
}

.nlmixr2libConventions <- function() {
  out <- .nlmixr2libConventionsStatic
  out$canonicalCovariates <- .loadCanonicalCovariates()
  out
}

.nlmixr2libCovariateAliasMap <- function() {
  conv <- .loadCanonicalCovariates()
  out <- character()
  for (canon in names(conv)) {
    for (a in conv[[canon]]$aliases) {
      out[a] <- canon
    }
  }
  out
}
