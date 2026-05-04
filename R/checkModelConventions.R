#' Check a model against nlmixr2lib conventions
#'
#' Parses a model and reports deviations from the nlmixr2lib conventions
#' documented in `vignettes/create-model-library.Rmd` and the
#' `extract-literature-model` skill references (especially
#' `naming-conventions.md` and `inst/references/covariate-columns.md`). The checker inspects:
#' file-level metadata (description, reference, units, covariateData);
#' parameter names (log-prefix PK params, `eta`-prefix IIV, `propSd`/`addSd`
#' residual error); parameter labels; covariates (canonical register, units,
#' declared aliases); compartment names; the observation variable (`Cc`); and
#' a syntactic dosing-vs-concentration unit cross-check.
#'
#' When any issue of severity `"error"` or `"warning"` is found, a `warning()`
#' is emitted whose message directs the user to call this function for the
#' full report.
#'
#' @param model Model to check. Accepted forms:
#'   - character scalar: resolved via [readModelDb()].
#'   - function: evaluated via `nlmixr2est::nlmixr()`.
#'   - `rxUi` object: used directly.
#'   - missing / `NULL`: iterate over every model in `modeldb`.
#' @param verbose Logical; if `TRUE` (default) print the per-category report
#'   to the console via `cli`.
#' @return Invisibly, a `data.frame` with one row per issue and columns
#'   `model`, `category`, `severity`, `name`, `message`, `suggestion`.
#' @export
#' @examples
#' \dontrun{
#' checkModelConventions("PK_1cmt_des")
#' checkModelConventions()            # check every model in modeldb
#' }
checkModelConventions <- function(model, verbose = TRUE) {
  checkmate::assertLogical(verbose, any.missing = FALSE, len = 1)
  if (missing(model) || is.null(model)) {
    names_all <- get0("modeldb", envir = asNamespace("nlmixr2lib"))$name
    if (is.null(names_all)) {
      stop("modeldb is not available in this package namespace",
           call. = FALSE)
    }
    issues <- do.call(rbind, lapply(names_all, function(nm) {
      .checkOneModel(nm, verbose = verbose)
    }))
    .emitSummaryWarning(issues)
    return(invisible(issues))
  }
  issues <- .checkOneModel(model, verbose = verbose)
  .emitSummaryWarning(issues)
  invisible(issues)
}

.checkOneModel <- function(model, verbose) {
  resolved <- .resolveModel(model)
  ui <- resolved$ui
  model_name <- resolved$name
  conv <- .nlmixr2libConventions()

  checks <- list(
    .checkFileMetadata(ui, conv),
    .checkParameterNames(ui, conv),
    .checkParameterLabels(ui, conv),
    .checkParameterUnits(ui, conv),
    .checkCovariates(ui, conv, model_name),
    .checkCompartments(ui, conv),
    .checkObservation(ui, conv),
    .checkUnits(ui, conv),
    .checkDeprecatedNames(ui, conv)
  )
  issues <- do.call(rbind, checks)
  if (is.null(issues) || nrow(issues) == 0) {
    issues <- .emptyIssues()
  }
  issues <- data.frame(
    model = rep(model_name, nrow(issues)),
    issues,
    stringsAsFactors = FALSE
  )
  if (verbose) .printReport(model_name, issues)
  issues
}

.resolveModel <- function(model) {
  if (is.character(model)) {
    checkmate::assertCharacter(model, len = 1, any.missing = FALSE,
                               min.chars = 1)
    fun <- readModelDb(model)
    ui <- nlmixr2est::nlmixr(fun)
    return(list(ui = ui, name = model))
  }
  if (inherits(model, "rxUi")) {
    nm <- tryCatch(model$modelName, error = function(e) NA_character_)
    nm <- .pickModelName(nm, fallback = "<rxUi>")
    return(list(ui = model, name = nm))
  }
  if (is.function(model)) {
    ui <- nlmixr2est::nlmixr(model)
    # rxode2 sets modelName to the formal parameter name (e.g. "model") when a
    # bare function is passed in, which is uninformative. Prefer the name
    # stashed in the function's closure by readModelDb(); fall back to the
    # explicit `<function>` placeholder so downstream scope checks can
    # recognize that the true identity is unknown rather than trusting the
    # misleading formal-parameter string.
    closure_nm <- tryCatch(get0("name", envir = environment(model),
                                inherits = FALSE),
                           error = function(e) NULL)
    nm <- .pickModelName(closure_nm, fallback = "<function>")
    return(list(ui = ui, name = nm))
  }
  stop("`model` must be a character name, function, rxUi, or missing.",
       call. = FALSE)
}

# rxode2 may populate `ui$modelName` with the full call chain that produced the
# UI (e.g., `c("readModelDb", "Valenzuela_2025_nipocalimab")` when a model was
# loaded via `rxode2(readModelDb("..."))`). Collapse to a single character so
# downstream `is.null(nm) || is.na(nm)` guards don't blow up on length > 1.
.pickModelName <- function(nm, fallback) {
  if (is.null(nm)) return(fallback)
  nm <- nm[!is.na(nm) & nzchar(nm)]
  if (length(nm) == 0) return(fallback)
  # Prefer the last (innermost) call-chain entry -- that's the actual model
  # function name, not the wrapper that returned it.
  nm[[length(nm)]]
}

.emptyIssues <- function() {
  data.frame(
    category = character(),
    severity = character(),
    name = character(),
    message = character(),
    suggestion = character(),
    stringsAsFactors = FALSE
  )
}

.issue <- function(category, severity, name, message, suggestion = NA_character_) {
  data.frame(
    category = category,
    severity = severity,
    name = name,
    message = message,
    suggestion = suggestion,
    stringsAsFactors = FALSE
  )
}

.emitSummaryWarning <- function(issues) {
  if (is.null(issues) || nrow(issues) == 0) return(invisible())
  bad <- issues[issues$severity %in% c("error", "warning"), , drop = FALSE]
  if (nrow(bad) == 0) return(invisible())
  models <- unique(bad$model)
  if (length(models) == 1) {
    msg <- sprintf(
      "%d convention issue(s) found for model '%s'. Run `checkModelConventions('%s')` for the full report.",
      nrow(bad), models, models
    )
  } else {
    msg <- sprintf(
      paste0(
        "%d convention issue(s) found across %d models. ",
        "Run `checkModelConventions(<name>)` for the full report per model."
      ),
      nrow(bad), length(models)
    )
  }
  warning(msg, call. = FALSE)
}

.printReport <- function(model_name, issues) {
  cli::cli_h1("Convention check: {.val {model_name}}")
  if (nrow(issues) == 0) {
    cli::cli_alert_success("No convention issues found.")
    return(invisible())
  }
  counts <- table(issues$severity)
  parts <- vapply(names(counts), function(s) {
    sprintf("%d %s", counts[[s]], s)
  }, character(1))
  cli::cli_alert_info("{paste(parts, collapse = ', ')}")
  for (cat in unique(issues$category)) {
    sub <- issues[issues$category == cat, , drop = FALSE]
    cli::cli_h2(cat)
    for (i in seq_len(nrow(sub))) {
      sev <- sub$severity[i]
      msg <- sprintf("[%s] %s: %s", sev, sub$name[i], sub$message[i])
      if (sev == "error") {
        cli::cli_alert_danger(msg)
      } else if (sev == "warning") {
        cli::cli_alert_warning(msg)
      } else {
        cli::cli_alert_info(msg)
      }
      if (!is.na(sub$suggestion[i]) && nzchar(sub$suggestion[i])) {
        cli::cli_text("    {.emph suggestion:} {sub$suggestion[i]}")
      }
    }
  }
}

.isEndogenousOrTemplate <- function(ui) {
  pred <- ui$predDf
  is.null(pred) || nrow(pred) == 0
}

.checkFileMetadata <- function(ui, conv) {
  issues <- .emptyIssues()
  meta <- as.list(ui$meta)
  endo <- .isEndogenousOrTemplate(ui)
  for (fld in conv$requiredMetadata) {
    val <- meta[[fld]]
    missing <- is.null(val) || (is.character(val) && !nzchar(val))
    if (!missing) next
    sev <- "error"
    if (fld == "reference" && endo) sev <- "info"
    if (fld == "units" && endo) sev <- "info"
    issues <- rbind(issues, .issue(
      "file_metadata", sev, fld,
      sprintf("Metadata field '%s' missing or empty.", fld),
      sprintf("Add `%s <- ...` at the top of the model function body.", fld)
    ))
  }
  uses_covariates <- length(ui$covariates %||% character()) > 0
  if (uses_covariates && is.null(meta$covariateData)) {
    issues <- rbind(issues, .issue(
      "file_metadata", "error", "covariateData",
      "Model references covariates but `covariateData` metadata is missing.",
      "Add a `covariateData <- list(<COV> = list(description=..., units=..., type=...))` block."
    ))
  }
  if (is.null(meta$population) && !endo) {
    issues <- rbind(issues, .issue(
      "file_metadata", "info", "population",
      "Optional `population` metadata block not present.",
      "Consider adding a `population` list (n_subjects, age_range, weight_range, ...) per the skill template."
    ))
  }
  issues
}

.classifyParam <- function(name, conv) {
  if (.isPkParam(name, conv)) return("canonical_pk")
  if (grepl(conv$covEffectPattern, name) && startsWith(name, "e_")) {
    return("cov_effect")
  }
  if (grepl("^l[A-Za-z]", name)) return("log_transformed")
  if (grepl("^logit[A-Za-z0-9]", name)) return("logit_transformed")
  if (grepl("^probit[A-Za-z0-9]", name)) return("probit_transformed")
  if (.isPkBareParam(name, conv)) return("bare_pk")
  "other"
}

.checkParameterNames <- function(ui, conv) {
  issues <- .emptyIssues()
  ini <- ui$iniDf
  if (is.null(ini) || nrow(ini) == 0) return(issues)

  fixed <- ini[is.na(ini$neta1) & is.na(ini$err), , drop = FALSE]
  # Only diagonal rows represent distinct IIV parameters; off-diagonal rows
  # carry the covariance of a block (e.g., name "(etalcl,etalvc)") and are
  # fully specified by the two diagonal entries.
  iiv <- ini[!is.na(ini$neta1) & !is.na(ini$neta2) &
               ini$neta1 == ini$neta2, , drop = FALSE]
  reserr <- ini[!is.na(ini$err), , drop = FALSE]

  for (nm in fixed$name) {
    cls <- .classifyParam(nm, conv)
    if (cls == "bare_pk") {
      issues <- rbind(issues, .issue(
        "parameter_naming", "warning", nm,
        sprintf("Fixed-effect PK parameter '%s' should be log-transformed (named 'l%s').", nm, nm),
        sprintf("Rename to 'l%s' in ini() and back-transform in model() via `%s <- exp(l%s)`.", nm, nm, nm)
      ))
    }
  }

  for (i in seq_len(nrow(iiv))) {
    nm <- iiv$name[i]
    if (!grepl("^eta", nm)) {
      issues <- rbind(issues, .issue(
        "parameter_naming", "warning", nm,
        sprintf("IIV parameter '%s' does not start with 'eta'.", nm),
        sprintf("Rename to 'eta<transformed-param>' (e.g., 'etalcl' for IIV on lcl).")
      ))
      next
    }
    suffix <- sub("^eta", "", nm)
    if (!(suffix %in% fixed$name)) {
      canonical <- paste0("l", suffix)
      if (canonical %in% fixed$name) {
        issues <- rbind(issues, .issue(
          "parameter_naming", "warning", nm,
          sprintf("IIV '%s' should include the log prefix of its parameter ('%s'), i.e., 'eta%s'.",
                  nm, canonical, canonical),
          sprintf("Rename '%s' to 'eta%s' in ini() and model() to match the transformed parameter name.",
                  nm, canonical)
        ))
      } else {
        issues <- rbind(issues, .issue(
          "parameter_naming", "warning", nm,
          sprintf("IIV '%s' has no matching fixed-effect parameter '%s'.", nm, suffix),
          "Ensure every eta<x> pairs with a fixed-effect parameter named x."
        ))
      }
    }
  }

  obs_vars <- unique(ui$predDf$cond %||% character())
  canonical_reserr <- .canonicalResidualErrorNames(obs_vars, conv)
  for (nm in reserr$name) {
    if (!(nm %in% canonical_reserr) && !.matchesDeprecatedReserr(nm, conv)) {
      issues <- rbind(issues, .issue(
        "parameter_naming", "warning", nm,
        sprintf(
          paste0(
            "Residual-error parameter '%s' does not match canonical ",
            "propSd/addSd (or '<output>propSd'/'<output>addSd')."
          ),
          nm
        ),
        sprintf("Rename to one of: %s.", paste(canonical_reserr, collapse = ", "))
      ))
    }
  }

  issues
}

.canonicalResidualErrorNames <- function(obs_vars, conv) {
  out <- conv$residualError
  for (v in obs_vars) {
    if (is.na(v)) next
    if (v == conv$observationVar) {
      # Parent observation Cc: canonical bare propSd / addSd. The
      # output-prefixed CcpropSd / CcaddSd form is no longer accepted
      # because every observation uses propSd_<X> / addSd_<X> with the
      # parent special-cased to the bare suffix-free form.
      next
    }
    # Every non-parent output uses parameter-name-then-output-suffix:
    # propSd_<x> / addSd_<x>. For metabolite outputs Cc_<metab> the
    # suffix is the metabolite name; for non-PK paper-named outputs
    # (Rtot, freeIgE, totalIgE, das28, tumorSize, ...) the suffix is
    # the output name itself, preserving its case.
    suffix <- v
    if (startsWith(v, "Cc_")) {
      candidate <- substr(v, 4, nchar(v))
      if (candidate %in% conv$registeredMetabolites) {
        suffix <- candidate
      }
    }
    out <- c(
      out,
      paste0("propSd_", suffix),
      paste0("addSd_",  suffix)
    )
  }
  unique(out)
}

.matchesDeprecatedReserr <- function(nm, conv) {
  nm %in% conv$deprecatedResidualError
}

.checkParameterLabels <- function(ui, conv) {
  issues <- .emptyIssues()
  ini <- ui$iniDf
  if (is.null(ini) || nrow(ini) == 0) return(issues)
  fixed_or_err <- ini[is.na(ini$neta1), , drop = FALSE]
  for (i in seq_len(nrow(fixed_or_err))) {
    nm <- fixed_or_err$name[i]
    lbl <- fixed_or_err$label[i]
    if (is.na(lbl) || !nzchar(lbl)) {
      issues <- rbind(issues, .issue(
        "parameter_labels", "warning", nm,
        sprintf("Parameter '%s' has no label.", nm),
        sprintf("Add `label(\"<description with units>\")` to the ini() entry for %s.", nm)
      ))
    }
  }
  issues
}

.checkParameterUnits <- function(ui, conv) {
  issues <- .emptyIssues()
  ini <- ui$iniDf
  if (is.null(ini) || nrow(ini) == 0) return(issues)
  fixed <- ini[is.na(ini$neta1) & is.na(ini$err), , drop = FALSE]
  for (i in seq_len(nrow(fixed))) {
    nm <- fixed$name[i]
    lbl <- fixed$label[i]
    if (is.na(lbl) || !nzchar(lbl)) next
    has_units_hint <- grepl("\\([^)]*[A-Za-z/0-9%][^)]*\\)", lbl)
    needs_units <- nm %in% conv$pkParams ||
      grepl("^l[a-z]", nm) ||
      nm %in% conv$pkBareParams
    if (needs_units && !has_units_hint) {
      issues <- rbind(issues, .issue(
        "parameter_units", "info", nm,
        sprintf("Label for '%s' does not appear to include a unit hint in parentheses.", nm),
        "Include units in the label, e.g. `label(\"Clearance (CL, L/day)\")`."
      ))
    }
  }
  issues
}

.checkCovariates <- function(ui, conv, model_name) {
  issues <- .emptyIssues()
  covs_all <- ui$covariates %||% character()
  depends <- as.list(ui$meta)$depends %||% character()
  # `depends` names are model inputs supplied by an upstream model (e.g., Cc
  # threaded into a PD model); they are not epidemiological covariates.
  covs <- setdiff(covs_all, depends)
  covData <- as.list(ui$meta)$covariateData
  covDataNames <- names(covData %||% list())

  alias_map <- .nlmixr2libCovariateAliasMap()
  canonical <- names(conv$canonicalCovariates)

  for (nm in covs) {
    if (!(nm %in% covDataNames)) {
      issues <- rbind(issues, .issue(
        "covariates", "error", nm,
        sprintf("Covariate '%s' is used in model() but has no entry in covariateData.", nm),
        sprintf("Add `%s = list(description=..., units=..., type=...)` to covariateData.", nm)
      ))
    }
    # Resolve the canonical entry for `nm` (may be nm itself or via alias map).
    canon_name <- if (nm %in% canonical) {
      nm
    } else if (nm %in% names(alias_map)) {
      alias_map[[nm]]
    } else {
      NA_character_
    }
    if (!(nm %in% canonical)) {
      if (nm %in% names(alias_map)) {
        canon <- alias_map[[nm]]
        entry <- covData[[nm]]
        declared <- !is.null(entry) && is.list(entry) &&
          !is.null(entry$source_name) && nzchar(entry$source_name)
        if (!declared) {
          issues <- rbind(issues, .issue(
            "covariates", "warning", nm,
            sprintf(
              paste0(
                "Covariate '%s' is an alias of canonical '%s'; ",
                "alias mapping is not declared via source_name."
              ),
              nm, canon
            ),
            sprintf(
              paste0(
                "Either rename to '%s' in the model and data, or ",
                "document the mapping by adding `source_name = \"%s\"` ",
                "under covariateData[[\"%s\"]] using the canonical name."
              ),
              canon, nm, canon
            )
          ))
        }
      } else {
        issues <- rbind(issues, .issue(
          "covariates", "warning", nm,
          sprintf("Covariate '%s' is not in the canonical register.", nm),
          "Rename to a canonical name (see `inst/references/covariate-columns.md`) or register a new canonical entry."
        ))
      }
    }
    # Scope check: a canonical covariate marked `scope: specific` is only
    # permitted in the models listed under its `Example models` field.
    # This catches authors who unknowingly reuse a study-specific name
    # (e.g., STUDY1, FORM_DP2, TUMTP_CHL) for a different purpose.
    # When model_name is a placeholder (pre-registration), the warning is an
    # informative todo: once the author registers under a real name they can
    # add it to the `Example models` list for covariates that are legitimate.
    if (!is.na(canon_name)) {
      entry <- conv$canonicalCovariates[[canon_name]]
      if (identical(entry$scope, "specific") &&
          !(model_name %in% entry$example_models)) {
        allowed <- if (length(entry$example_models) == 0) {
          "(none)"
        } else {
          paste(entry$example_models, collapse = ", ")
        }
        issues <- rbind(issues, .issue(
          "covariates", "warning", nm,
          sprintf(
            paste0(
              "Covariate '%s' is canonical but scoped 'specific' to ",
              "model(s) %s; using it in '%s' is not permitted."
            ),
            nm, allowed, model_name
          ),
          sprintf(
            paste0(
              "Rename to a different canonical name, promote '%s' to ",
              "`Scope: general` in inst/references/covariate-columns.md, ",
              "or add '%s' to that entry's Example models list."
            ),
            canon_name, model_name
          )
        ))
      }
    }
  }

  for (nm in covDataNames) {
    entry <- covData[[nm]]
    if (is.character(entry)) {
      issues <- rbind(issues, .issue(
        "covariates", "warning", nm,
        sprintf("covariateData[['%s']] is a bare string (old style); missing structured units/type.", nm),
        "Convert to a list: `list(description=..., units=..., type=<continuous|binary|categorical|count>)`."
      ))
      next
    }
    if (!is.list(entry)) next
    if (is.null(entry$description) || !nzchar(entry$description %||% "")) {
      issues <- rbind(issues, .issue(
        "covariates", "error", nm,
        sprintf("covariateData[['%s']] is missing `description`.", nm),
        "Add a one-line `description` field."
      ))
    }
    if (is.null(entry$units) || !nzchar(entry$units %||% "")) {
      issues <- rbind(issues, .issue(
        "covariates", "error", nm,
        sprintf("covariateData[['%s']] is missing `units`.", nm),
        "Add a `units` field (use \"(binary)\" or \"(categorical)\" where appropriate)."
      ))
    }
    if (!(nm %in% covs)) {
      issues <- rbind(issues, .issue(
        "covariates", "warning", nm,
        sprintf("covariateData[['%s']] has an entry but is not referenced in model().", nm),
        "Remove the unused entry or confirm the covariate is still used."
      ))
    }
  }

  issues
}

.checkCompartments <- function(ui, conv) {
  issues <- .emptyIssues()
  cmts <- ui$props$cmt %||% character()
  for (cm in cmts) {
    if (.matchesCompartment(cm, conv)) next
    issues <- rbind(issues, .issue(
      "compartments", "warning", cm,
      sprintf("Compartment '%s' is not a canonical name.", cm),
      sprintf(
        paste0(
          "Use one of: %s; numbered chains (transit<n>, effect<n>, ",
          "precursor<n>, lat<n>); DAR-numbered (dar<n>_central, ",
          "dar<n>_peripheral<m>); or metabolite-suffixed ",
          "(<canonical>_<metab>, e.g. central_mmae). Registered ",
          "metabolites: %s. For new payloads or therapeutic-area ",
          "compartments, register the new metabolite/compartment in ",
          "R/conventions.R first."
        ),
        paste(conv$compartments, collapse = ", "),
        paste(conv$registeredMetabolites, collapse = ", ")
      )
    ))
  }
  issues
}

.checkObservation <- function(ui, conv) {
  issues <- .emptyIssues()
  pred <- ui$predDf
  if (is.null(pred) || nrow(pred) == 0) return(issues)
  obs_vars <- unique(pred$cond)
  if (length(obs_vars) == 1 && obs_vars != conv$observationVar) {
    issues <- rbind(issues, .issue(
      "observation", "warning", obs_vars,
      sprintf("Single-output observation variable '%s' should be named '%s'.",
              obs_vars, conv$observationVar),
      sprintf("Rename observation and its residual-error assignment to '%s'.",
              conv$observationVar)
    ))
    return(issues)
  }
  # Multi-output: flag deprecated `C<metab>` style (e.g. Cmmae, Cdxd, Cdar0)
  # and suggest the canonical `Cc_<metab>` form for PK metabolite outputs.
  for (v in obs_vars) {
    if (is.na(v) || v == conv$observationVar) next
    if (startsWith(v, "Cc_")) next  # already canonical metabolite output
    if (!startsWith(v, "C")) next   # non-PK output (tumorSize, freeIgE, ...)
    rest <- substr(v, 2, nchar(v))
    rest_lc <- tolower(rest)
    if (rest_lc %in% conv$registeredMetabolites) {
      issues <- rbind(issues, .issue(
        "observation", "warning", v,
        sprintf(
          paste0(
            "Multi-output observation '%s' uses the deprecated 'C<metab>' ",
            "form for a PK metabolite output."
          ),
          v
        ),
        sprintf("Rename to 'Cc_%s' to match the canonical metabolite-suffix convention.",
                rest_lc)
      ))
    }
  }
  issues
}

.checkUnits <- function(ui, conv) {
  issues <- .emptyIssues()
  units <- as.list(ui$meta)$units
  if (is.null(units)) return(issues)
  endo <- .isEndogenousOrTemplate(ui)
  for (fld in conv$requiredUnits) {
    if (is.null(units[[fld]]) || !nzchar(units[[fld]] %||% "")) {
      sev <- if (endo) "info" else "error"
      issues <- rbind(issues, .issue(
        "units", sev, fld,
        sprintf("units$%s is missing.", fld),
        sprintf("Add `%s = \"<unit>\"` to the `units <- list(...)` metadata block.", fld)
      ))
    }
  }
  dosing <- units$dosing
  conc <- units$concentration
  if (!is.null(dosing) && !is.null(conc) && nzchar(dosing) && nzchar(conc)) {
    if (!grepl("/", conc)) {
      issues <- rbind(issues, .issue(
        "units", "warning", "concentration",
        sprintf("units$concentration='%s' does not contain '/' (mass/volume).", conc),
        "Concentration units usually look like 'mg/L', 'ug/mL', 'ng/mL'."
      ))
    } else {
      conc_num <- trimws(sub("/.*$", "", conc))
      if (!.unitsCompatible(dosing, conc_num) &&
          !.unitsSameDimension(dosing, conc_num)) {
        issues <- rbind(issues, .issue(
          "units", "warning", "dosing_concentration",
          sprintf(
            paste0(
              "units$dosing ('%s') and units$concentration numerator ",
              "('%s') appear dimensionally incompatible."
            ),
            dosing, conc_num
          ),
          paste0(
            "Confirm the dosing unit dimension matches the concentration ",
            "numerator dimension (both mass, or both molar, etc.)."
          )
        ))
      } else if (!.unitsCompatible(dosing, conc_num)) {
        issues <- rbind(issues, .issue(
          "units", "info", "dosing_concentration",
          sprintf(
            paste0(
              "units$dosing ('%s') and units$concentration numerator ",
              "('%s') differ in magnitude; ensure scaling is applied in model()."
            ),
            dosing, conc_num
          ),
          paste0(
            "When dosing is mg but concentration is ug/mL (= mg/L), no ",
            "conversion is needed if volume is in L. Verify the relationship."
          )
        ))
      }
    }
  }
  issues
}

.unitsCompatible <- function(a, b) {
  if (identical(a, b)) return(TRUE)
  identical(.normUnit(a), .normUnit(b))
}

.normUnit <- function(x) {
  x <- gsub("\u00b5", "u", x)
  x <- gsub("\u03bc", "u", x)
  tolower(trimws(x))
}

.unitsSameDimension <- function(a, b) {
  mass <- c("kg", "g", "mg", "ug", "ng", "pg")
  molar <- c("mol", "mmol", "umol", "nmol", "pmol")
  iu <- c("iu", "miu", "uiu", "niu")
  a <- .normUnit(a)
  b <- .normUnit(b)
  (a %in% mass && b %in% mass) ||
    (a %in% molar && b %in% molar) ||
    (a %in% iu && b %in% iu)
}

.checkDeprecatedNames <- function(ui, conv) {
  issues <- .emptyIssues()
  ini <- ui$iniDf
  if (is.null(ini) || nrow(ini) == 0) return(issues)
  for (nm in ini$name) {
    if (nm %in% conv$deprecatedResidualError) {
      issues <- rbind(issues, .issue(
        "deprecated_names", "warning", nm,
        sprintf("'%s' is a deprecated residual-error name.", nm),
        "Rename to 'propSd' (proportional) or 'addSd' (additive), with output prefix for multi-output models."
      ))
    }
    for (pfx in conv$deprecatedIivPrefixes) {
      if (startsWith(nm, pfx)) {
        issues <- rbind(issues, .issue(
          "deprecated_names", "warning", nm,
          sprintf("'%s' uses deprecated IIV prefix '%s'.", nm, pfx),
          "Use `eta<transformed-param>` (e.g., 'etalcl')."
        ))
        break
      }
    }
    issues <- rbind(issues, .checkDeprecatedVolumeOrVmaxName(nm, conv))
    issues <- rbind(issues, .checkDeprecatedAdcSuffix(nm, conv))
    issues <- rbind(issues, .checkDeprecatedCovEffectSuffix(nm, conv))
  }
  issues
}

# Flag bare-volume and Vmax names that have canonical replacements:
#   Change: v / v1 / lv / lv1 -> vc / lvc
#   Change: v2 / lv2 / v3 / lv3 -> ambiguous (vp / vp2); ask the source paper
#   Change: vm / lvm -> vmax / lvmax
.checkDeprecatedVolumeOrVmaxName <- function(nm, conv) {
  issues <- .emptyIssues()
  if (nm %in% c("v", "v1", "lv", "lv1")) {
    canonical <- if (startsWith(nm, "l")) "lvc" else "vc"
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf("'%s' is a deprecated central-volume name.", nm),
      sprintf("Rename to '%s' (canonical central volume).", canonical)
    ))
  } else if (nm %in% c("v2", "lv2")) {
    canonical <- if (startsWith(nm, "l")) "lvp (likely) or lvc if v1 is depot" else "vp (likely) or vc"
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf("'%s' is a deprecated numbered-volume name.", nm),
      sprintf("Verify against the source paper and rename to %s.", canonical)
    ))
  } else if (nm %in% c("v3", "lv3")) {
    canonical <- if (startsWith(nm, "l")) "lvp2" else "vp2"
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf("'%s' is a deprecated numbered-volume name.", nm),
      sprintf("Rename to '%s' (second peripheral volume).", canonical)
    ))
  } else if (nm %in% c("vm", "lvm")) {
    canonical <- if (startsWith(nm, "l")) "lvmax" else "vmax"
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf("'%s' is a deprecated Michaelis-Menten Vmax name.", nm),
      sprintf("Rename to '%s'.", canonical)
    ))
  }
  issues
}

# Flag the parent `_adc` suffix: a parent-side ADC parameter should drop
# the `_adc` since the parent uses canonical names. Detect both
# parameter names (lcl_adc) and covariate-effect names (e_wt_cl_adc).
.checkDeprecatedAdcSuffix <- function(nm, conv) {
  issues <- .emptyIssues()
  if (endsWith(nm, "_adc")) {
    suggested <- substr(nm, 1, nchar(nm) - 4)
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf(
        paste0(
          "'%s' uses the deprecated parent-ADC suffix '_adc'. ",
          "Parent-drug parameters should use the canonical name ",
          "without a suffix; metabolite parameters carry the suffix."
        ),
        nm
      ),
      sprintf("Drop the '_adc' suffix and rename to '%s'.", suggested)
    ))
  }
  issues
}

# Flag covariate-effect names whose parameter token is a deprecated
# numbered/synthesized form. Examples:
#   Change: e_wt_v / _v1 -> e_wt_vc;  _v2 -> e_wt_vp (verify)
#   Change: e_wt_clq -> e_wt_cl_q;  e_wt_vcvp -> e_wt_vc_vp
#   Change: e_wt_clinf -> e_wt_cl_ss;  e_wt_clt -> e_wt_cl_time
#   Change: e_wt_vss -> e_wt_vc_vp;  e_wt_clss -> e_wt_cl_ss
.checkDeprecatedCovEffectSuffix <- function(nm, conv) {
  issues <- .emptyIssues()
  if (!startsWith(nm, "e_")) return(issues)
  if (!grepl(conv$covEffectPattern, nm)) return(issues)
  # Reverse-order: e_<param>_<cov>(_<extra>)*. Detect when the FIRST
  # token is a bare PK parameter (cl, vc, vp, q, ka, ...) and the SECOND
  # token (or the joined remainder for compound covariates) maps to a
  # known canonical covariate. The canonical order is e_<cov>_<param>.
  rest <- substr(nm, 3, nchar(nm))
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    first <- parts[1]
    canonical_covs <- names(conv$canonicalCovariates %||% list())
    canonical_covs_lc <- tolower(canonical_covs)
    alias_map <- .nlmixr2libCovariateAliasMap()
    canonical_aliases_lc <- tolower(names(alias_map))
    is_known_cov <- function(tok) {
      tolower(tok) %in% c(canonical_covs_lc, canonical_aliases_lc)
    }
    # Try the second token, then second+third, then second+third+fourth
    # — covers compound covariates like RACE_BLACK, ADA_POSITIVE,
    # FORM_CHO_PHASE2.
    second_or_more <- character()
    for (j in 2:length(parts)) {
      second_or_more <- c(
        second_or_more,
        paste(parts[2:j], collapse = "_")
      )
    }
    if (first %in% conv$pkBareParams &&
        any(vapply(second_or_more, is_known_cov, logical(1)))) {
      issues <- rbind(issues, .issue(
        "deprecated_names", "warning", nm,
        sprintf(
          paste0(
            "'%s' looks like a reversed-order covariate effect ",
            "(e_<param>_<cov>); the canonical order is e_<cov>_<param>."
          ),
          nm
        ),
        sprintf("Rename to e_<cov>_<param>; verify the covariate identity in the source.")
      ))
      return(issues)
    }
  }
  # Trailing token deprecations (parameter token that's a deprecated form).
  trailing_token <- parts[length(parts)]
  trailing_map <- list(
    v   = "vc",
    v1  = "vc",
    v2  = "vp (verify against source)",
    v3  = "vp2",
    vm  = "vmax",
    clq = "cl_q (split shared exponent)",
    vcvp = "vc_vp (split shared exponent)",
    vss = "vc_vp (Vss = Vc + Vp; split shared exponent)",
    clinf = "cl_ss",
    clss  = "cl_ss",
    clt   = "cl_time"
  )
  if (trailing_token %in% names(trailing_map)) {
    canonical <- trailing_map[[trailing_token]]
    issues <- rbind(issues, .issue(
      "deprecated_names", "warning", nm,
      sprintf(
        paste0(
          "'%s' uses a deprecated parameter-token '%s' in a ",
          "covariate-effect name."
        ),
        nm, trailing_token
      ),
      sprintf("Rename the parameter portion to '%s'.", canonical)
    ))
  }
  issues
}

`%||%` <- function(a, b) if (is.null(a)) b else a
