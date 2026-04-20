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
    if (is.null(nm) || is.na(nm)) nm <- "<rxUi>"
    return(list(ui = model, name = nm))
  }
  if (is.function(model)) {
    ui <- nlmixr2est::nlmixr(model)
    nm <- tryCatch(ui$modelName, error = function(e) NA_character_)
    if (is.null(nm) || is.na(nm)) nm <- "<function>"
    return(list(ui = ui, name = nm))
  }
  stop("`model` must be a character name, function, rxUi, or missing.",
       call. = FALSE)
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
  if (name %in% conv$pkParams) return("canonical_pk")
  if (grepl(conv$covEffectPattern, name)) return("cov_effect")
  if (grepl("^l[A-Za-z]", name)) return("log_transformed")
  if (grepl("^logit[A-Za-z0-9]", name)) return("logit_transformed")
  if (grepl("^probit[A-Za-z0-9]", name)) return("probit_transformed")
  if (name %in% conv$pkBareParams) return("bare_pk")
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
    if (!is.na(v) && v != conv$observationVar) {
      out <- c(out, paste0(v, "propSd"), paste0(v, "addSd"))
    } else if (!is.na(v) && v == conv$observationVar) {
      out <- c(out, "CcpropSd", "CcaddSd")
    }
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
  allowed <- conv$compartments
  for (cm in cmts) {
    if (cm %in% allowed) next
    if (grepl(conv$compartmentRegex, cm)) next
    issues <- rbind(issues, .issue(
      "compartments", "warning", cm,
      sprintf("Compartment '%s' is not a canonical name.", cm),
      sprintf("Use one of: %s (or transit<n>). For new therapeutic-area compartments, open a GitHub issue first.",
              paste(allowed, collapse = ", "))
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
  }
  issues
}

`%||%` <- function(a, b) if (is.null(a)) b else a
