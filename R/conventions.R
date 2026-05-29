#' nlmixr2lib convention standards
#'
#' Internal register of canonical parameter, compartment, covariate, and
#' residual-error names used by [checkModelConventions()]. The canonical
#' covariate, parameter, and compartment / metabolite-suffix lists are parsed
#' at runtime from `inst/references/covariate-columns.md`,
#' `inst/references/parameter-names.md`, and
#' `inst/references/compartment-names.md`
#' (installed as `system.file("references", "<file>.md", package = "nlmixr2lib")`),
#' so the markdown registers remain the single authoritative sources.
#' The static list below retains only structural regex constants and the
#' deprecation lists. Remaining fields mirror the `extract-literature-model`
#' skill's `naming-conventions.md` and `vignettes/create-model-library.Rmd`.
#'
#' @keywords internal
#' @noRd
.nlmixr2libConventionsStatic <- list(
  # Numbered-chain compartment pattern. Bare numbered chains (transit /
  # effect / precursor / lat / depot) and metabolite-suffixed
  # compartments are validated separately via .matchesCompartment() so
  # that the registered metabolite list can be honored at runtime; this
  # static regex covers only the numbered-chain patterns. `depot[0-9]+`
  # accommodates parallel-absorption models with two or more depots.
  compartmentRegex = "^(transit|effect|precursor|lat|depot)[0-9]+$",
  # Membrane-limited PBPK sub-compartment pattern: paper-prefix +
  # spelled-out organ name. Recognises the recurring `<sub>_<organ>`
  # shape used in Shah 2012 mAb PBPK and Parhiz 2024 mRNA-LNP
  # extractions (bc / eu / eb / fr / is / int / mrna / luc prefixes).
  pbpkSubCompartmentRegex = "^(bc|eu|eb|fr|is|int|mrna|luc)_(liver|lung|kidney|spleen|heart|muscle|skin|adipose|bone|brain|small_intestine|large_intestine|pancreas|thymus|portal|remainder|other|hepatic|fat|rapidly_perfused|slowly_perfused|venous|arterial|urine|gut)$",
  # DAR-numbered ADC isoform compartments (`dar0_central`,
  # `dar4_peripheral1`, ...).
  darCompartmentRegex = "^dar[0-9]+_(central|peripheral[0-9]?)$",
  # Target species in physiologic body-fluid or named peripheral
  # compartments (e.g., target_csf, target_isf, target_peripheral,
  # target_peripheral1, complex_csf, complex_isf, complex_peripheral).
  targetLocationRegex = "^(target|complex)_(csf|isf|peripheral[0-9]?)$",
  observationVar = "Cc",
  # propSd and addSd are the canonical proportional and additive
  # residual-error SDs; expSd is the log-scale residual SD used with
  # `~ lnorm(...)`.
  residualError = c("propSd", "addSd", "expSd"),
  transformPrefixes = c("l", "logit", "probit"),
  # Covariate-effect names match e_<cov>(_<continuation>)+_<param>.
  # The pattern accepts up to 6 underscore-separated tokens after the
  # leading "e_" to accommodate compound covariates (RACE_BLACK,
  # ADA_POSITIVE, FORM_CHO_PHASE2). Semantic interpretation of the
  # trailing tokens is handled by .classifyCovEffect().
  covEffectPattern = "^e_[A-Za-z0-9]+(_[A-Za-z0-9]+){1,5}$",
  # Suffixes allowed for multi-component CL parameters. `_ss` denotes
  # the steady-state arm; `_time` the time-varying decay arm; `_renal`
  # the glomerular-filtration / tubular-secretion arm; `_nonren` the
  # non-renal (hepatic / metabolic / extra-renal) arm (e.g. Jonckheere
  # 2019 cefepime: CL_total = CL_renal + CL_nonren).
  clComponents = c("ss", "time", "renal", "nonren"),
  requiredUnits = c("time", "dosing", "concentration"),
  requiredMetadata = c("description", "reference", "units"),
  deprecatedResidualError = c(
    "prop.err", "add.err", "propErr", "addErr",
    "err.prop", "err.add"
  ),
  deprecatedIivPrefixes = c("iiv_", "IIV_", "bsv_", "BSV_"),
  # Bare volume names that should be replaced with vc / vp / vp2.
  deprecatedVolumeNames = c("v", "v1", "v2", "v3", "lv", "lv1", "lv2", "lv3"),
  # Deprecated Michaelis-Menten Vmax names.
  deprecatedVmaxNames = c("vm", "lvm"),
  # Deprecated parent-suffix marker. A model that names a parent-side
  # parameter `<base>_adc` should drop the `_adc` suffix; the parent
  # uses the canonical name unsuffixed.
  deprecatedParentSuffix = "_adc"
)

# The following canonical-name lists are NOT carried on
# .nlmixr2libConventionsStatic above; they are parsed at runtime from
# the markdown registers in inst/references and merged into the
# returned object by .nlmixr2libConventions(): pkParams, pkBareParams,
# and paperNamedParams come from parameter-names.md; compartments and
# registeredMetabolites come from compartment-names.md; and
# canonicalCovariates comes from covariate-columns.md.

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
#' (`### NAME (**...**)`). For each entry, captures the `Units`, `Type`,
#' `Scope`, `Source aliases`, and `Example models` fields. Aliases whose
#' backticked content is not a bare R identifier (e.g. `DVID = "study1"`)
#' are skipped. Example-model tokens are accepted as backticked file names
#' ending in `.R`; the `.R` suffix is stripped so the value matches the
#' bare model function name used throughout the rest of the package.
#'
#' @param path Path to the markdown file.
#' @return A named list keyed by canonical name. Each entry is a list with
#'   `units`, `type`, `scope` (one of `"general"` / `"specific"` / `NA`),
#'   `aliases` (character vector of alias names), and `example_models`
#'   (character vector of model function names).
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
        scope = current$scope %||% NA_character_,
        aliases = current$aliases %||% character(),
        example_models = current$example_models %||% character()
      )
    }
  }

  aliasRegex <- "^\\s*-\\s*`([^`]+)`"
  identRegex <- "^[A-Za-z_][A-Za-z0-9_]*$"
  modelFileRegex <- "^[A-Za-z_][A-Za-z0-9_-]*\\.R$"

  extractBacktickedModels <- function(text) {
    toks <- regmatches(text, gregexpr("`([^`]+)`", text))[[1]]
    models <- character()
    for (tok in toks) {
      inner <- gsub("`", "", tok)
      if (grepl(modelFileRegex, inner)) {
        models <- c(models, sub("\\.R$", "", inner))
      }
    }
    models
  }

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
      current <- list(names = nms, aliases = character(),
                      example_models = character())
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
    m <- regmatches(line, regexec("^- \\*\\*Scope:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      scope_raw <- tolower(trimws(sub("\\.$", "", m[[2]])))
      if (scope_raw %in% c("general", "specific")) {
        current$scope <- scope_raw
      }
      state <- "header"
      next
    }
    if (grepl("^- \\*\\*Source aliases:\\*\\*", line)) {
      state <- "aliases"
      after <- sub("^- \\*\\*Source aliases:\\*\\*\\s*", "", line)
      # "none", "none known", "none;"-style declarations have no aliases.
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
    if (grepl("^- \\*\\*Example models:\\*\\*", line)) {
      state <- "example_models"
      after <- sub("^- \\*\\*Example models:\\*\\*\\s*", "", line)
      current$example_models <- c(current$example_models,
                                  extractBacktickedModels(after))
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

    if (state == "example_models") {
      # Continuation bullet lines in a multi-line Example-models list.
      if (grepl("^\\s+-\\s", line)) {
        current$example_models <- c(current$example_models,
                                    extractBacktickedModels(line))
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

.parameterNamesPath <- function() {
  p <- system.file("references", "parameter-names.md",
                   package = "nlmixr2lib")
  if (nzchar(p)) return(p)
  stop("Could not locate inst/references/parameter-names.md in the ",
       "nlmixr2lib package. Check the installation.", call. = FALSE)
}

.compartmentNamesPath <- function() {
  p <- system.file("references", "compartment-names.md",
                   package = "nlmixr2lib")
  if (nzchar(p)) return(p)
  stop("Could not locate inst/references/compartment-names.md in the ",
       "nlmixr2lib package. Check the installation.", call. = FALSE)
}

#' Parse a canonical-name markdown register from inst/references/.
#'
#' Walks the markdown file and extracts one entry per H3 heading
#' (`### NAME (**...**)`). For each entry, captures `Type`, `Source aliases`,
#' and `Example models`. A single canonical name may appear under more than
#' one H3 entry (with distinct Type values); each appearance is returned as
#' a separate list element so the caller can route by Type without losing
#' duplicates.
#'
#' @param path Path to the markdown file.
#' @return A list of entries; each entry is a list with `name`, `type`,
#'   `aliases` (character vector), and `example_models` (character vector).
#' @keywords internal
#' @noRd
.parseTypedNamesMd <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- list()
  current <- NULL
  state <- "idle"

  identRegex <- "^[A-Za-z_0-9][A-Za-z0-9_]*$"
  modelFileRegex <- "^[A-Za-z_][A-Za-z0-9_-]*\\.R$"

  extractBacktickedModels <- function(text) {
    toks <- regmatches(text, gregexpr("`([^`]+)`", text))[[1]]
    models <- character()
    for (tok in toks) {
      inner <- gsub("`", "", tok)
      if (grepl(modelFileRegex, inner)) {
        models <- c(models, sub("\\.R$", "", inner))
      }
    }
    models
  }

  flush <- function() {
    if (is.null(current)) return(invisible())
    if (is.null(current$type) || !nzchar(current$type)) return(invisible())
    for (nm in current$names) {
      entries[[length(entries) + 1]] <<- list(
        name = nm,
        type = current$type,
        aliases = current$aliases %||% character(),
        example_models = current$example_models %||% character()
      )
    }
  }

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
      current <- list(names = nms, aliases = character(),
                      example_models = character())
      state <- "header"
      next
    }
    if (is.null(current)) next

    m <- regmatches(line, regexec("^- \\*\\*Type:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$type <- trimws(m[[2]])
      state <- "header"
      next
    }
    if (grepl("^- \\*\\*Source aliases:\\*\\*", line)) {
      state <- "aliases"
      after <- sub("^- \\*\\*Source aliases:\\*\\*\\s*", "", line)
      if (grepl("^none\\b", after, ignore.case = TRUE)) next
      after <- strsplit(after, "\\s+—\\s+", perl = TRUE)[[1]][1]
      inline <- regmatches(after, gregexpr("`([^`]+)`", after))[[1]]
      for (tok in inline) {
        inner <- gsub("`", "", tok)
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
      }
      next
    }
    if (grepl("^- \\*\\*Example models:\\*\\*", line)) {
      state <- "example_models"
      after <- sub("^- \\*\\*Example models:\\*\\*\\s*", "", line)
      current$example_models <- c(current$example_models,
                                  extractBacktickedModels(after))
      next
    }

    if (state == "aliases") {
      m <- regmatches(line, regexec("^\\s*-\\s*`([^`]+)`", line))[[1]]
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

    if (state == "example_models") {
      if (grepl("^\\s+-\\s", line)) {
        current$example_models <- c(current$example_models,
                                    extractBacktickedModels(line))
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

.namesByType <- function(entries, type) {
  out <- character()
  for (e in entries) {
    if (identical(e$type, type)) out <- c(out, e$name)
  }
  out
}

.loadCanonicalParameters <- function(force = FALSE) {
  if (!force && !is.null(.covariateRegisterCache$parameters)) {
    return(.covariateRegisterCache$parameters)
  }
  entries <- .parseTypedNamesMd(.parameterNamesPath())
  .covariateRegisterCache$parameters <- entries
  entries
}

.loadCanonicalCompartments <- function(force = FALSE) {
  if (!force && !is.null(.covariateRegisterCache$compartments)) {
    return(.covariateRegisterCache$compartments)
  }
  entries <- .parseTypedNamesMd(.compartmentNamesPath())
  .covariateRegisterCache$compartments <- entries
  entries
}

.nlmixr2libConventions <- function() {
  out <- .nlmixr2libConventionsStatic
  out$canonicalCovariates <- .loadCanonicalCovariates()
  paramEntries <- .loadCanonicalParameters()
  out$pkParams <- .namesByType(paramEntries, "log-transformed-pk")
  out$pkBareParams <- .namesByType(paramEntries, "bare-pk")
  out$paperNamedParams <- .namesByType(paramEntries, "paper-named-param")
  compEntries <- .loadCanonicalCompartments()
  out$compartments <- .namesByType(compEntries, "compartment")
  out$registeredMetabolites <- .namesByType(compEntries, "metabolite-suffix")
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

# Return TRUE when `name` is a canonical log-transformed PK parameter or
# a metabolite-suffixed PK parameter (`l<base>_<metab>`).
.isPkParam <- function(name, conv) {
  if (name %in% conv$pkParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkParams) return(TRUE)
    }
  }
  FALSE
}

# Return TRUE when `name` is a canonical bare PK parameter or a
# metabolite-suffixed bare PK parameter (`<base>_<metab>`).
.isPkBareParam <- function(name, conv) {
  if (name %in% conv$pkBareParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkBareParams) return(TRUE)
    }
  }
  FALSE
}

# Compartment name validator. Recognizes:
#   - canonical names from conv$compartments
#   - numbered chains via conv$compartmentRegex (transit/effect/precursor/lat/depot)
#   - DAR-numbered ADC isoforms via conv$darCompartmentRegex
#   - target species in physiologic compartments via conv$targetLocationRegex
#   - metabolite-suffixed compartments: <canonical>_<metab>
.matchesCompartment <- function(name, conv) {
  if (name %in% conv$compartments) return(TRUE)
  if (grepl(conv$compartmentRegex, name)) return(TRUE)
  if (grepl(conv$darCompartmentRegex, name)) return(TRUE)
  if (grepl(conv$targetLocationRegex, name)) return(TRUE)
  if (!is.null(conv$pbpkSubCompartmentRegex) &&
      grepl(conv$pbpkSubCompartmentRegex, name)) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$compartments) return(TRUE)
      # Compositions of a numbered-chain prefix with a metabolite
      # suffix are canonical: e.g., `transit1_m3g`, `precursor2_dxd`,
      # `lat1_complex`. Used by formation-delay transit chains feeding
      # a metabolite central compartment (deHoogd 2017 morphine model:
      # transit1_m3g..transit5_m3g and transit1_m6g..transit2_m6g).
      if (grepl(conv$compartmentRegex, base)) return(TRUE)
      # Recursive: strip one metabolite suffix and re-check whether
      # the base matches another canonical pattern (compartment list,
      # chain regex, DAR regex, target-location regex, PBPK sub-
      # compartment regex, or another metabolite-suffixed canonical).
      # Lets multi-suffix compartments pass, e.g.
      # `liver_endo_asn1` -> strip `_asn1` -> `liver_endo` -> strip
      # `_endo` -> `liver` (canonical). Used by Ayyar 2024 givosiran
      # parent + metabolite PBPK extraction.
      if (.matchesCompartment(base, conv)) return(TRUE)
    }
  }
  FALSE
}

# TRUE when `name` ends with `_<metab>` for any registered metabolite.
.endsWithMetabolite <- function(name, conv) {
  for (metab in conv$registeredMetabolites) {
    if (endsWith(name, paste0("_", metab))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<component>` for any registered CL component.
.endsWithClComponent <- function(name, conv) {
  for (comp in conv$clComponents) {
    if (endsWith(name, paste0("_", comp))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<param>` for any bare PK parameter.
# Used to detect shared-exponent covariate effects like e_wt_cl_q.
.endsWithBarePkParam <- function(name, conv) {
  for (p in conv$pkBareParams) {
    if (endsWith(name, paste0("_", p))) return(TRUE)
  }
  FALSE
}

# Classify a covariate-effect name by its trailing suffix. Returns one of:
#   "two_token"  - matches e_<cov>_<param> with no third-token suffix
#   "metabolite" - matches e_<cov>_<param>_<metab>
#   "shared"     - matches e_<cov>_<param>_<param2> (shared exponent)
#   "component"  - matches e_<cov>_<param>_<component> (multi-CL arm)
#   "unknown"    - has a third-token suffix that doesn't match any
#                  registered category
.classifyCovEffect <- function(name, conv) {
  if (!startsWith(name, "e_")) return(NA_character_)
  if (!grepl(conv$covEffectPattern, name)) return(NA_character_)
  if (.endsWithMetabolite(name, conv)) return("metabolite")
  if (.endsWithClComponent(name, conv)) return("component")
  if (.endsWithBarePkParam(name, conv)) return("shared")
  # Strip the leading `e_` and check whether the rest is a single
  # `<cov>_<param>` pair (no third-token suffix). If yes, two_token;
  # otherwise the name has an unrecognized trailing suffix.
  rest <- substr(name, 3, nchar(name))
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  if (length(parts) == 2) return("two_token")
  "unknown"
}
