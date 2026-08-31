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
    .checkUnitSpellings(ui, conv),
    .checkCovariates(ui, conv, model_name),
    .checkCompartments(ui, conv),
    .checkObservation(ui, conv),
    .checkUnits(ui, conv),
    .checkDeprecatedNames(ui, conv),
    .checkFixedLabelAgreement(ui, conv),
    .checkTimeVaryingClearanceNames(ui, conv),
    .checkCompartmentData(ui, conv),
    .checkFmFamily(ui, conv)
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
  # `depends` lists names that are provided by an upstream model when
  # this model is composed downstream (e.g. PD templates inheriting Cc
  # from a PK model); subtract them so they don't trigger a spurious
  # covariateData-missing error.
  depends <- meta$depends %||% character()
  covs_for_check <- setdiff(ui$covariates %||% character(), depends)
  uses_covariates <- length(covs_for_check) > 0
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

# Recognize a shared-eta suffix of the form `l<p1>_<p2>(_<p3>...)` where
# every "_"-separated token corresponds to a fixed-effect parameter
# (matched either as the bare token `<pi>` or with the canonical log
# prefix `l<pi>`). Used to accept names like `etalkdeg_kint` (a single
# eta shared between `lkdeg` and `lkint`) without flagging them as
# missing a structural pair.
.isSharedEtaSuffix <- function(suffix, fixed_names) {
  if (!startsWith(suffix, "l")) return(FALSE)
  body <- substr(suffix, 2, nchar(suffix))
  if (!grepl("_", body, fixed = TRUE)) return(FALSE)
  parts <- strsplit(body, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  matches <- vapply(parts, function(p) {
    p %in% fixed_names || paste0("l", p) %in% fixed_names
  }, logical(1))
  all(matches)
}

# Recognize a within-subject-level variability eta suffix of the form
# `iov_<param>_<occ>` or `bvv_<param>_<visit>` where `<param>` corresponds
# to an existing fixed-effect parameter (matched as the bare token
# `<param>`, with the canonical log prefix `l<param>`, or with the
# typical-value prefixes `ltv<param>` / `lv<param>`) and the trailing
# integer is the occasion / visit index. Used to accept names like
# `etaiov_k13_1` (a single IOV eta on `k13` for occasion 1) and
# `etabvv_vmax_1` (a between-visit eta on `vmax` for visit 1) without
# flagging them as missing a structural pair.
#
# `iov_` and `bvv_` are separate prefixes because a source paper can fit
# BOTH levels simultaneously on different parameters: Abdelgawad 2024
# linezolid carries five-occasion IOV on `ka` and `mtt` alongside a
# two-visit BVV on `vmax`, and collapsing them would lose the published
# random-effect hierarchy.
.isIOVEtaSuffix <- function(suffix, fixed_names) {
  prefix <- NULL
  for (p in c("iov_", "bvv_")) {
    if (startsWith(suffix, p)) prefix <- p
  }
  if (is.null(prefix)) return(FALSE)
  rest <- substr(suffix, nchar(prefix) + 1L, nchar(suffix))
  m <- regmatches(rest, regexec("^(.+)_([0-9]+)$", rest))[[1]]
  if (length(m) != 3L) return(FALSE)
  param <- m[[2L]]
  param %in% fixed_names ||
    paste0("l",   param) %in% fixed_names ||
    paste0("ltv", param) %in% fixed_names ||
    paste0("lv",  param) %in% fixed_names
}

# Accept paper-mechanistic stratified-typical-value etas that don't
# have a 1-to-1 matching `lX` parameter because the underlying
# typical value is split by lesion, population, study, age bracket,
# or biomarker. Two acceptance routes:
#
#   1. Any fixed-effect parameter shares a prefix with the eta's
#      suffix. Example: etalbase_les1 (suffix = "lbase_les1") accepts
#      because the ini block has lbase_les1, lbase_les2, ... typical
#      values. etalcl_form_m3g (suffix = "lcl_form_m3g") accepts
#      because the ini block has lcl_form_m3g_le10 / lcl_form_m3g_gt10
#      age-bracketed typical values. eta_study_pmax_f accepts because
#      the ini block has e.g. lkp_f, dmax_f, pmax_f, etc., one of
#      which starts with "study" stripped. Implementation: an ini
#      param `p` matches the eta suffix `s` when
#      startsWith(p, s) || startsWith(s, p), with at least one strict
#      prefix-extension (i.e. the names are not identical because that
#      case is already handled by the strict pairing rule above).
#
#   2. The eta suffix matches a known paper-mechanistic stratification
#      pattern (eta_study_*, eta<X>_les<n>, eta<X>_pop<n>, eta<X>_l<n>,
#      eta<X>_px<n>) AND the eta's underlying parameter root matches
#      an existing ini parameter (with optional l/ltv/lv log-prefix).
#      Example: etalbase_les1 strips the `_les1` suffix to `lbase`
#      and accepts when `lbase` is in the ini block. etalS0_mtd_pop1
#      strips `_pop1` to `lS0_mtd` and accepts when `lS0_mtd_pop2`
#      etc. exists.
.isPaperMechanisticEtaSuffix <- function(suffix, fixed_names) {
  # Route 1: prefix-extension match.
  for (p in fixed_names) {
    if (nchar(p) > 0L &&
        (startsWith(p, suffix) || startsWith(suffix, p)) &&
        p != suffix) {
      return(TRUE)
    }
  }
  # Route 2: known stratification suffix patterns. Strip the trailing
  # stratification token (_les<n>, _pop<n>, _l<n>, _px<n>) and look
  # for any parameter that starts with the base. This is more lenient
  # than Route 1 because it lets etalbase_les1 pair with lbase_suv /
  # lbase_sld (same `lbase` root, different per-paper sub-strata).
  base <- sub("_(les|pop|px|l|sld|suv)[0-9]*$", "", suffix)
  if (base != suffix && nchar(base) > 0L) {
    for (p in fixed_names) {
      if (startsWith(p, base)) return(TRUE)
    }
  }
  # Study-stratified pattern: eta_study_<param>_<stratum>.
  if (startsWith(suffix, "study_")) {
    rest <- substr(suffix, 7L, nchar(suffix))
    # Accept if the trailing token matches a stratum-grouped ini
    # param. We don't enforce the exact stratification family so the
    # check is lenient.
    for (p in fixed_names) {
      if (grepl(rest, p, fixed = TRUE)) return(TRUE)
    }
  }
  # IOV with a paper-named occasion grouping (etaiov_<group>_<n>):
  # if the eta is etaiov-prefixed and the group name appears as part
  # of any ini parameter, accept. Catches etaiov_bio_1 / etaiov_bio_2
  # where the paper's bioavailability typical value is implied by a
  # named ini parameter (lf, fdepot, etc.) rather than a literal
  # `bio` ini name.
  if (startsWith(suffix, "iov_")) {
    return(TRUE)
  }
  # Paper-named etas where the suffix is a generic paper-mechanistic
  # parameter name that the source paper uses as an additive shift
  # (etalogit, etap1..p5, etaibase, etafrel, etalmrt_pooled,
  # etaclge_px<n>, etalec50). These don't have a 1-to-1 fixed-effect
  # parameter because the underlying typical value is part of a
  # paper-mechanistic structural equation (logit transform of a
  # mixture probability, transit-rate fractions, indirect-response
  # baseline, ...). Accept any eta whose suffix:
  #   (a) starts with a letter and contains no further `_`-separated
  #       semantic tokens that would conflict with a structural-PK
  #       canonical, OR
  #   (b) matches a `<root>_px<n>` / `<root>_<biomarker>` / `_pooled`
  #       paper-mechanistic stratification suffix.
  if (grepl("_(pooled|px[0-9]+|vegf|svegfr2|svegfr3|skit|hb|f|pmax|dmax|lkp|lkdrug|ic50|drug_slope)$", suffix)) {
    return(TRUE)
  }
  FALSE
}

# Recognize an inter-study / inter-arm-variability (IAV) eta of the form
# `eta_study_<param>` where `<param>` corresponds to an existing
# fixed-effect parameter. The eta name passed in has had its leading
# "eta" stripped, so `suffix` is "_study_<param>" (note the leading
# underscore -- this is why the `.isPaperMechanisticEtaSuffix` "study_"
# branch, which expects no leading underscore, does not catch it). Used
# by model-based meta-analysis (MBMA) extractions that encode
# study-arm-level random effects in place of between-subject variability
# (Wang 2018 daclatasvir/asunaprevir; Boucher 2018 naproxen). Acceptance
# routes for `<param>`:
#   (a) `<param>` (or `l<param>`) is itself a fixed-effect parameter
#       -- covers eta_study_lcl (lcl), eta_study_e0 (e0), eta_study_ld_asv;
#   (b) some fixed-effect name extends `<param>` (startsWith);
#   (c) token-bridge: a fixed effect shares `<param>`'s first AND last
#       underscore-separated tokens -- covers a model()-internal combined
#       typical value built from sibling anchors, e.g. logit_fk_asv
#       pairing with logit_fk_cap_asv / logit_fk_sol_asv.
.isStudyEtaSuffix <- function(suffix, fixed_names) {
  if (!startsWith(suffix, "_study_")) return(FALSE)
  param <- substr(suffix, 8L, nchar(suffix))   # strip leading "_study_"
  if (!nzchar(param)) return(FALSE)
  if (param %in% fixed_names) return(TRUE)
  if (paste0("l", param) %in% fixed_names) return(TRUE)
  for (p in fixed_names) {
    if (nchar(p) > 0L && startsWith(p, param) && p != param) return(TRUE)
  }
  ptok <- strsplit(param, "_", fixed = TRUE)[[1]]
  if (length(ptok) >= 2L) {
    first <- ptok[[1L]]
    last <- ptok[[length(ptok)]]
    for (p in fixed_names) {
      ftok <- strsplit(p, "_", fixed = TRUE)[[1]]
      if (length(ftok) >= 2L &&
          ftok[[1L]] == first &&
          ftok[[length(ftok)]] == last) return(TRUE)
    }
  }
  FALSE
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

  # Optional paper-specific exception fields (analogous to the
  # `paper_specific_compartments` mechanism for compartment names):
  #   paper_specific_etas <- c("etalogit", "etap1", ...) -- IIV names
  #     whose typical-value parameter is a paper-mechanistic structural
  #     equation rather than a 1-to-1 `lX` ini parameter; the
  #     "no matching fixed-effect parameter" check is skipped.
  #   paper_specific_residual_sds <- c("propSd_vact_l1", ...) -- residual
  #     SD names with paper-specific multi-token output suffixes that
  #     the canonical propSd_<output> matcher does not recognise.
  meta <- as.list(ui$meta)
  paper_specific_etas <- meta$paper_specific_etas %||% character()
  paper_specific_reserr <- meta$paper_specific_residual_sds %||% character()

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
      } else if (.isSharedEtaSuffix(suffix, fixed$name)) {
        # Shared-eta on multiple structural parameters: etal<p1>_<p2>...
        # is accepted when every "_"-separated token matches an existing
        # fixed-effect parameter (with an optional leading "l" prefix).
        # No issue is emitted.
      } else if (.isIOVEtaSuffix(suffix, fixed$name)) {
        # Inter-occasion variability eta of the form
        # etaiov_<param>_<occ> where <param> is an existing fixed-effect
        # parameter (with optional l/ltv/lv log-prefix) and <occ> is a
        # positive integer occasion index. No issue is emitted.
      } else if (.isStudyEtaSuffix(suffix, fixed$name)) {
        # Inter-study / inter-arm-variability eta of the form
        # eta_study_<param> (MBMA study-arm-level random effect) where
        # <param> is an existing fixed-effect parameter (bare, log-
        # prefixed, prefix-extended, or token-bridged to a combined
        # typical value). No issue is emitted.
      } else if (.isPaperMechanisticEtaSuffix(suffix, fixed$name)) {
        # Paper-mechanistic stratified-typical-value etas where the
        # underlying typical value is split by lesion, population,
        # study, age bracket, or biomarker (e.g. etalbase_les1..les5
        # pairing with lbase_les1_*, etalcl_form_m3g pairing with
        # lcl_form_m3g_le10 / _gt10 age brackets, eta_study_<param>_<f|hb>
        # pairing with paper-stratified ini params, etc.). No issue
        # is emitted.
      } else if (nm %in% paper_specific_etas) {
        # Paper-specific eta declared via the `paper_specific_etas`
        # metadata field. The author has explicitly documented that
        # the underlying typical-value parameter is a paper-mechanistic
        # structural equation rather than a 1-to-1 `lX` ini parameter.
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
    if (nm %in% paper_specific_reserr) next
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
      paste0("addSd_",  suffix),
      paste0("expSd_",  suffix),
      paste0("powExp_", suffix)
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

# Issue: the library wrote the same time unit three ways in the machine-readable
# `units` block -- "hour" (643 models), "h" (208) and "hr" (28) -- so a consumer
# parsing units$time could not canonicalise without a spelling table of its own.
# This checks SPELLING, not dimension: "min" and "h" are both canonical and are
# never conflated, because normalising between them would misstate every value.
.checkUnitSpellings <- function(ui, conv) {
  issues <- .emptyIssues()
  units <- as.list(ui$meta)$units
  if (is.null(units)) return(issues)
  spellings <- list(time = conv$timeUnitSpellings, dosing = conv$doseUnitSpellings)
  for (fld in names(spellings)) {
    val <- units[[fld]]
    if (is.null(val) || !is.character(val) || length(val) != 1 || is.na(val)) next
    if (val %in% conv$placeholderUnits) next
    map <- spellings[[fld]]
    key <- tolower(val)
    if (key %in% names(map)) {
      issues <- rbind(issues, .issue(
        "unit_spelling", "error", sprintf("units$%s", fld),
        sprintf("units$%s = '%s' is a non-canonical spelling.", fld, val),
        sprintf("Use '%s'. Same unit, one spelling -- consumers cannot canonicalise otherwise.",
                map[[key]])))
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

  # Covariates the source model documents but intentionally does NOT use
  # in model() (e.g. FREM-screened demographics whose effects were not
  # clinically meaningful and have no published point estimate) live in a
  # `covariatesDataExcluded` metadata list (same shape as covariateData).
  # They are documentation only: the checker does not require them to
  # appear in model() and does not flag them as unused. The one guard --
  # a name listed there must NOT also be referenced in model() (that
  # would be a mis-filing; it belongs in covariateData instead).
  covExcluded <- as.list(ui$meta)$covariatesDataExcluded
  for (nm in names(covExcluded %||% list())) {
    if (nm %in% covs_all) {
      issues <- rbind(issues, .issue(
        "covariates", "warning", nm,
        sprintf("Covariate '%s' is listed in covariatesDataExcluded but is referenced in model().", nm),
        "Move it to covariateData (it is actually used), or stop referencing it in model()."
      ))
    }
  }

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
    # Check that the covariateData entry corresponds to a name actually
    # used inside model(). Use the unfiltered ui$covariates list so
    # that names also declared in `depends` (e.g. operational data-table
    # identifiers like TYPE / DILmer / DILcol / STUDY_DD that the
    # model body reads but treats as upstream inputs rather than
    # epidemiological effect modifiers) are correctly recognised when
    # they appear in BOTH depends and covariateData.
    if (!(nm %in% covs_all)) {
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
  # Model files may declare paper-mechanistic compartments that are not
  # in the canonical register via a `paper_specific_compartments`
  # metadata field. The validator subtracts these from the warning set
  # so the author can explicitly document a model's per-paper named
  # states (similar in spirit to the `depends` mechanism for upstream
  # covariate inputs). Two declaration forms are accepted:
  #
  #   paper_specific_compartments <- c("name1", "name2", ...)
  #     a literal vector of compartment names to skip
  #
  #   paper_specific_compartment_pattern <- "^bact_"
  #     a single regex pattern matched against compartment names; if
  #     any element matches, the name is skipped
  meta <- as.list(ui$meta)
  paper_specific <- meta$paper_specific_compartments %||% character()
  paper_specific_re <- meta$paper_specific_compartment_pattern %||% character()
  for (cm in cmts) {
    if (.matchesCompartment(cm, conv)) next
    if (length(paper_specific) > 0L && cm %in% paper_specific) next
    if (length(paper_specific_re) > 0L &&
        any(sapply(paper_specific_re,
                   function(p) grepl(p, cm)))) next
    # Tailor the message: a capital-prefixed name is almost never
    # canonical because the convention is lowercase compartment names
    # (observation variables like Cc are the exception). Surface that
    # specifically so the fix is obvious.
    if (grepl("^[A-Z]", cm)) {
      msg <- sprintf(
        paste0(
          "Compartment '%s' starts with an uppercase letter; the ",
          "convention is lowercase compartment names (observation ",
          "variables like Cc are the only canonical capitalized form)."
        ),
        cm
      )
      sug <- paste0(
        "Rename the state to lowercase. If a capital-prefixed ",
        "concentration alias is needed for output mapping, derive it ",
        "as a non-state assignment in model() (e.g., ",
        "`d/dt(brain) <- ...; Cbrain <- brain/vbrain`)."
      )
    } else {
      msg <- sprintf("Compartment '%s' is not a canonical name.", cm)
      sug <- sprintf(
        paste0(
          "Use one of: %s; numbered chains (transit<n>, effect<n>, ",
          "precursor<n>, lat<n>, depot<n>); DAR-numbered ",
          "(dar<n>_central, dar<n>_peripheral<m>); target species in a ",
          "physiologic compartment (target_csf, target_isf, ",
          "complex_csf, complex_isf); or metabolite-suffixed ",
          "(<canonical>_<metab>, e.g. central_mmae). Registered ",
          "metabolites: %s. For new payloads or therapeutic-area ",
          "compartments, register the new metabolite/compartment in ",
          "R/conventions.R first."
        ),
        paste(conv$compartments, collapse = ", "),
        paste(conv$registeredMetabolites, collapse = ", ")
      )
    }
    issues <- rbind(issues, .issue(
      "compartments", "warning", cm, msg, sug
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
    obs <- obs_vars
    # Cc is canonical for drug-concentration outputs; per the 2026-05-28
    # naming-audit operator clarification, single-output PD models may
    # use any registered output-state name (tumor_size, das28, ANC via
    # circ_anc, etc.). Treat the observation as canonical when it
    # matches the compartment register (which now includes the PD-output
    # canonicals registered in D19) or one of the canonical observation
    # variants (Cc itself and the Cc_<metab> metabolite outputs). Also
    # accept derived `C<canonical-compartment>` aliases (e.g.
    # Cbrain_csf <- ... where brain_csf is the underlying canonical
    # compartment); the C-prefix denotes the concentration-derived
    # output state corresponding to a registered compartment amount.
    is_canon_pd <- .matchesCompartment(obs, conv) ||
      startsWith(obs, "Cc_") ||
      (startsWith(obs, "C") && nchar(obs) > 1 &&
       .matchesCompartment(substr(obs, 2, nchar(obs)), conv))
    if (!is_canon_pd) {
      issues <- rbind(issues, .issue(
        "observation", "warning", obs,
        sprintf(
          paste0(
            "Single-output observation variable '%s' is not canonical: ",
            "use 'Cc' for drug-concentration outputs, or a registered ",
            "PD-output canonical compartment / state name otherwise."
          ),
          obs
        ),
        sprintf(
          paste0(
            "Rename to '%s' for plasma-drug-concentration outputs, ",
            "or register the PD output name as a canonical compartment ",
            "in R/conventions.R if it is a recurring paper-mechanistic ",
            "endpoint."
          ),
          conv$observationVar
        )
      ))
    }
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
  # The dosing<->concentration dimensional checks below are only
  # meaningful for classical PK models where `dosing` is a simple
  # administered amount (mass / molar / IU) and `concentration` is a
  # plasma-drug concentration. They are skipped when dosing is a rate or
  # concentration (contains "/"), is declared not-applicable, or is
  # otherwise not a clean amount token -- i.e. PD / MBMA / in-vitro /
  # endogenous models. They are further RELAXED for PD / QSP / bacterial-
  # kill endpoints (2026-06-28 naming-audit decision):
  #   * a `concentration` whose numerator is not a recognized
  #     pharmacological amount is a biological / PD endpoint (CFU, cells,
  #     mmHg, ordinal grade, receptor occupancy, probability, ...) and a
  #     mass/volume check does not apply -> skipped (no issue);
  #   * a `concentration` numerator that IS a pharmacological amount but of
  #     a different dimension than the dosing amount (mass dose vs molar
  #     concentration, mg vs IU, ...) is reported as `info`, not a warning,
  #     because the molecular-weight / potency conversion is expected and
  #     handled in model().
  # Descriptive parentheticals on either token (e.g.
  # "nmol (convert mg via dose_nmol = ...)") are stripped before the
  # comparison so they do not spuriously defeat a same-dimension match.
  # The relaxation only removes or downgrades issues; it never introduces a
  # new warning.
  dose_is_amount <- !endo && !is.null(dosing) && nzchar(dosing) &&
    .unitsRecognizedAmount(dosing)
  if (dose_is_amount && !is.null(conc) && nzchar(conc)) {
    dose_core <- trimws(sub("\\(.*$", "", dosing))
    if (!grepl("/", conc)) {
      # No mass/volume slash. Warn only when the concentration field is
      # itself a bare pharmacological amount (an amount written where a
      # concentration belongs); a non-amount PD endpoint (mmHg, ordinal
      # grade, ...) is a legitimate non-plasma-concentration output.
      if (.unitsRecognizedAmount(conc)) {
        issues <- rbind(issues, .issue(
          "units", "warning", "concentration",
          sprintf("units$concentration='%s' does not contain '/' (mass/volume).", conc),
          "Concentration units usually look like 'mg/L', 'ug/mL', 'ng/mL'."
        ))
      }
    } else {
      conc_num <- trimws(sub("\\(.*$", "", trimws(sub("/.*$", "", conc))))
      if (!.unitsRecognizedAmount(conc_num)) {
        # Non-pharmacological-amount numerator (CFU, cells, mU, ...):
        # PD / QSP / bacterial-kill endpoint -> mass/volume check N/A, skip.
      } else if (.unitsCompatible(dose_core, conc_num)) {
        # Identical pharmacological dimension and token -> nothing to flag.
      } else if (.unitsSameDimension(dose_core, conc_num)) {
        issues <- rbind(issues, .issue(
          "units", "info", "dosing_concentration",
          sprintf(
            paste0(
              "units$dosing ('%s') and units$concentration numerator ",
              "('%s') differ in magnitude; ensure scaling is applied in model()."
            ),
            dose_core, conc_num
          ),
          paste0(
            "When dosing is mg but concentration is ug/mL (= mg/L), no ",
            "conversion is needed if volume is in L. Verify the relationship."
          )
        ))
      } else {
        # Both are recognized amounts but of different pharmacological
        # dimensions (mass vs molar vs IU). The conversion (molecular
        # weight / potency) is expected and performed in model() ->
        # informational, not a warning.
        issues <- rbind(issues, .issue(
          "units", "info", "dosing_concentration",
          sprintf(
            paste0(
              "units$dosing ('%s') and units$concentration numerator ",
              "('%s') differ in dimension (mass vs molar vs IU); ensure the ",
              "molecular-weight / potency conversion is applied in model()."
            ),
            dose_core, conc_num
          ),
          paste0(
            "Confirm the dosing-to-concentration conversion (e.g. mg -> nmol ",
            "via molecular weight) is performed in model()."
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

# TRUE when `x` is a single recognized administered-amount unit token
# (mass / molar / IU / mEq), ignoring case and any trailing parenthetical
# descriptor. Used to gate the dosing<->concentration dimensional checks
# to classical PK models: when dosing is a rate or concentration (so the
# token contains "/"), is "n/a" / "not applicable", or is any other
# non-amount string, this returns FALSE and the checks are skipped.
.unitsRecognizedAmount <- function(x) {
  core <- .normUnit(trimws(sub("\\(.*$", "", x)))
  amounts <- c(
    "kg", "g", "mg", "ug", "ng", "pg",
    "mol", "mmol", "umol", "nmol", "pmol",
    "iu", "miu", "uiu", "niu",
    "meq", "mu", "u"
  )
  core %in% amounts
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
    if (!is.null(conv$renamedParameters) && nm %in% names(conv$renamedParameters)) {
      issues <- rbind(issues, .issue(
        "deprecated_names", "error", nm,
        sprintf("'%s' was renamed; this spelling defeats name-based discovery.", nm),
        sprintf("Use '%s'.", conv$renamedParameters[[nm]])
      ))
    }
    issues <- rbind(issues, .checkDeprecatedVolumeOrVmaxName(nm, conv))
    issues <- rbind(issues, .checkDeprecatedAdcSuffix(nm, conv))
    issues <- rbind(issues, .checkDeprecatedCovEffectSuffix(nm, conv))
  }
  issues
}

# Flag bare-volume and Vmax names that have canonical replacements.
# Change deprecated central-volume names "v", "v1", "lv", "lv1" to
# "vc" or "lvc"; numbered "v2", "lv2", "v3", "lv3" are ambiguous
# (could be vp or vp2) and must be verified against the source paper;
# Michaelis-Menten names "vm" and "lvm" become "vmax" and "lvmax".
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
# numbered or synthesized form. Examples of expected renames: the
# deprecated tokens "v" and "v1" become "vc" (e.g. an effect named
# e_wt_v should become e_wt_vc) and "v2" becomes "vp" (verify against
# the source paper); the run-together tokens "clq" and "vcvp" become
# "cl_q" and "vc_vp" with an underscore separator; "clinf" and "clt"
# become "cl_ss" and "cl_time" respectively; and "vss" and "clss"
# become "vc_vp" and "cl_ss".
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


# Severity is "error". The backlog this rule found when it was added (15
# parameters in 10 models beyond the ones issue #479 enumerated by hand) has
# been cleared, so any hit is a new mistake. Keep it that way: if this ever
# needs to be demoted to "warning" to land a change, the change is
# reintroducing the class the rule exists to prevent.
# History and the false-positive classes are in
# inst/references/fixed-provenance-followup.md.
#
# Issue #479: `iniDf$fix` is the only machine-readable signal that a value was
# not estimated from the study's own data. Models that SAY in the label that a
# parameter was fixed, assumed, borrowed from another paper, or set to a
# theoretical value, while leaving it estimable, cause a downstream consumer to
# treat an assumption as evidence -- e.g. counting a fixed 0.75 allometric
# exponent as data supporting the 0.75 convention.
# A label claims THIS parameter was not estimated only in a self-referential
# construction. Matching the bare word "fixed" anywhere in the label produces
# mostly false positives, because labels routinely describe a DIFFERENT
# parameter as fixed: "the typical-value Ka is fixed but eta on Ka is
# estimated", "corr(CL,Vp) = 0 (fixed)", "equivalent under fixed Km",
# "F_HIV- = 1 fixed". Anchor on trailing qualifiers and on verb forms that
# take the parameter itself as their subject.
.fixedClaimPattern <- paste0(
  # The claim must open its own clause. When another name precedes it in the
  # same clause the sentence is about THAT parameter, not this one -- e.g.
  # "IC50 ... (Imax fixed to 1)", "CL/F_DHA (L/h/kg); F_DHA fixed to 1",
  # "CL_DM4 (L/day; V_DM4 fixed to 1 L)", "(correlation fixed to 1)",
  # "P-gp components fixed to 0". Anchoring on the clause start keeps the
  # self-referential "(unitless; fixed)" and "(fixed in source)" forms and
  # drops the rest.
  "(?:^|[;(),\u2014-])\\s*",
  # "fixed effect" is the mixed-models term for a population parameter, not a
  # claim that the value was held constant (Lehr 2010: "fixed effect, no IIV
  # per the source paper", reported with RSE 6.9%).
  "(?:fixed(?!\\s+effects?\\b)|assumed|not published|not paper-derived",
  "|literature value|taken from|transferred from)\\b"
)

.estimatedDisclaimerPattern <- "\\b(?:is|was|were|are) estimated\\b|\\bestimated\\b[^.;]*\\bper\\b"

# Issue #481: a clearance that depends on time must say so in its parameter
# names. The old spellings cannot simply be banned by name -- `emax`, `imax`,
# `gamma`, `hill`, `t50` and `kdes` are all legitimate PD names elsewhere in
# the library -- so the check is structural: find a clearance symbol whose
# definition references `t` or `time`, then require the canonical stem to
# appear in that same expression.
.clearanceLhsPattern <- "^\\s*(cl[a-z0-9_]*|[a-z0-9_]*_cl|td_cl)\\s*<-"
.bareTimePattern <- "(?<![A-Za-z0-9_.])(t|time)(?![A-Za-z0-9_])"
# The canonical stems that make a time-dependent clearance expression
# acceptable. Kept as ONE constant because the registry-wide enumerating test in
# test-checkModelConventions.R applies the same rule to the model files on disk:
# when this lived as a duplicated literal in both places, ratifying `tclchange`
# updated the checker and left the test flagging the very model that founded the
# canonical. `tclchange` is the breakpoint symbol of the piecewise-constant
# (step / NONMEM MTIME) form, which is the only one of the three whose early arm
# keeps the plain `cl` name -- so there is no `cl_*` stem on the left of its
# switch and the breakpoint is what identifies the structure.
.timeVaryingClearanceAcceptPattern <- "cl_time_max|cl_t50|cl_time_hill|cl_exp_|tclchange"

# The `ini({})` block declares parameters and their labels; only `model({})`
# contains the equations. Scanning the whole function makes label prose such as
# "Sigmoidicity exponent of time on CL" look like a reference to `t`.
.modelBlockLines <- function(ui) {
  lines <- tryCatch(deparse(ui$fun), error = function(e) character(0))
  if (!length(lines)) return(character(0))
  start <- grep("^\\s*model\\(\\{", lines)
  if (!length(start)) return(character(0))
  lines[seq(start[1], length(lines))]
}

.checkTimeVaryingClearanceNames <- function(ui, conv) {
  issues <- .emptyIssues()
  lines <- .modelBlockLines(ui)
  if (!length(lines)) return(issues)
  for (ln in lines) {
    if (!grepl(.clearanceLhsPattern, ln)) next
    rhs <- sub("^[^<]*<-", "", ln)
    rhs <- sub("#.*$", "", rhs)
    # Strip string literals: a label such as "...exponent of time on CL..."
    # is prose about the parameter, not a reference to the time variable.
    rhs <- gsub("\"[^\"]*\"", "", rhs)
    if (!grepl(.bareTimePattern, rhs, perl = TRUE)) next
    if (grepl(.timeVaryingClearanceAcceptPattern, rhs)) next
    nm <- trimws(sub("\\s*<-.*$", "", ln))
    issues <- rbind(issues, .issue(
      "time_varying_clearance", "warning", nm,
      sprintf("'%s' makes clearance depend on time but uses none of the canonical names.", nm),
      paste("Use cl_time_max / cl_t50 / cl_time_hill for a sigmoidal-in-time",
            "clearance, or cl_exp_inf / cl_exp_component / cl_exp_kdes for an",
            "exponential decay, so the structure can be found by name (issue #481).")
    ))
  }
  issues
}

# Issue #482: what molecule a compartment holds, in what units, and in what
# biological matrix. Inferring the matrix from the state name works only
# because most models happen to use `central` / `peripheral1`; where they do
# not, inference fails silently or gets it wrong.
.checkCompartmentData <- function(ui, conv) {
  issues <- .emptyIssues()
  meta <- as.list(ui$meta)
  states <- tryCatch(ui$state, error = function(e) character(0))
  if (!length(states)) return(issues)
  cd <- meta$compartmentData
  if (is.null(cd)) {
    # Warning rather than error while the database is being backfilled; see
    # inst/references/compartment-data-followup.md for the remaining models.
    return(rbind(issues, .issue(
      "compartment_data", "warning", "compartmentData",
      "Model has ODE states but no `compartmentData` metadata.",
      paste("Add compartmentData <- list(<state> = list(analyte=, units=,",
            "specimen=, verified=)) for every d/dt() state (issue #482).")
    )))
  }
  missing <- setdiff(states, names(cd))
  if (length(missing)) {
    issues <- rbind(issues, .issue(
      "compartment_data", "error", paste(missing, collapse = ", "),
      "ODE state(s) have no `compartmentData` entry.",
      "Add one entry per d/dt() state."))
  }
  extra <- setdiff(names(cd), states)
  if (length(extra)) {
    issues <- rbind(issues, .issue(
      "compartment_data", "error", paste(extra, collapse = ", "),
      "`compartmentData` names a compartment that is not an ODE state.",
      "Remove the entry, or correct the state name."))
  }
  for (nm in intersect(names(cd), states)) {
    entry <- cd[[nm]]
    absent <- setdiff(conv$compartmentDataFields, names(entry))
    if (length(absent)) {
      issues <- rbind(issues, .issue(
        "compartment_data", "error", nm,
        sprintf("`compartmentData$%s` is missing: %s.", nm, paste(absent, collapse = ", ")),
        paste("Every entry needs", paste(conv$compartmentDataFields, collapse = ", "), ".")))
      next
    }
    if (!entry$specimen %in% conv$specimenVocabulary) {
      issues <- rbind(issues, .issue(
        "compartment_data", "error", nm,
        sprintf("`%s` is not in the specimen vocabulary.", entry$specimen),
        paste("Use one of:", paste(conv$specimenVocabulary, collapse = ", "))))
    }
  }
  issues
}

# "fixed" in a label whose value is already `fixed(...)`. Excludes compound
# terms where "fixed" modifies something else and is not a claim about this
# parameter's estimation status: "fixed effect" is the mixed-models term,
# "fixed dose" / "fixed-dose regimen" describes a dosing strategy.
.redundantFixedPattern <- "\\bfixed(?!\\s*[- ]\\s*(?:effects?|dose|dosing))\\b"

.checkFixedLabelAgreement <- function(ui, conv) {
  issues <- .emptyIssues()
  ini <- ui$iniDf
  if (is.null(ini) || nrow(ini) == 0) return(issues)
  if (!all(c("name", "label", "fix") %in% names(ini))) return(issues)
  for (i in seq_len(nrow(ini))) {
    lbl <- ini$label[[i]]
    if (is.na(lbl) || !nzchar(lbl)) next
    if (isTRUE(ini$fix[[i]])) {
      # Mirror case: the value IS fixed, so saying so again in the label is
      # redundant. `fixed()` is the machine-readable statement; the label
      # should carry only what `fixed()` cannot express -- where the value
      # came from ("from Rizk 2015"), or that the encoder assumed it.
      if (grepl(.redundantFixedPattern, lbl, ignore.case = TRUE, perl = TRUE)) {
        issues <- rbind(issues, .issue(
          "fixed_label_redundant", "error", ini$name[[i]],
          sprintf("Label of '%s' says the value is fixed, which `fixed()` already states.",
                  ini$name[[i]]),
          paste("Drop the word from the label. Keep any provenance around it --",
                "\"fixed from Rizk 2015\" becomes \"from Rizk 2015\".")))
      }
      next
    }
    # Variance terms: their labels almost always discuss the fixed status of
    # the corresponding typical value, not of the variance itself.
    if (grepl("^eta", ini$name[[i]])) next
    if (grepl(.estimatedDisclaimerPattern, lbl, ignore.case = TRUE)) next
    if (grepl(.fixedClaimPattern, lbl, ignore.case = TRUE, perl = TRUE)) {
      issues <- rbind(issues, .issue(
        "fixed_label_disagreement", "error", ini$name[[i]],
        sprintf("Label of '%s' says the value was fixed/assumed/borrowed, but fix = FALSE.",
                ini$name[[i]]),
        paste("Wrap the value in `fixed(...)` so `iniDf$fix` matches the label,",
              "or reword the label if it was in fact estimated.")
      ))
    }
  }
  issues
}

# Fraction-metabolised pathway family (`fm_<pathway>`) -----------------------
#
# `fm` has been registered as a canonical for a long time, but the `<pathway>`
# suffixes were not, so nothing consumed the fact that two models had spelled
# the same pathway differently. `fm_other` and `fm_others` both shipped, as did
# `fm_h4` and `fm_H4`; the drift was found by hand in 2026-08, not by the
# build. Registering the family in parameter-names.md is what makes this check
# enumerating: a pathway is legal because it appears in that heading, so a new
# spelling fails here until someone adds it deliberately.
#
# Scope is every name BOUND inside `ini()` or `model()` -- parameters, not
# prose. That distinction is load-bearing, and it is why this reads the parsed
# `ui` rather than grepping the source:
#   * Source-paper notation quoted in a `label()` string is prose.
#     `Svensson_2013_bedaquiline.R` labels its apparent metabolite parameters
#     "...CL_M2/(F*fm_M2)...", matching the `CL_M2` / `V_M2` symbols beside it.
#   * `fm` also abbreviates FAT MASS, the canonical `FM` covariate.
#     `Robarge_2017_efavirenz.R` carries `fm_range` in its `population`
#     metadata next to `ffm_range`.
# Both are correct as written; neither is a fraction metabolised, and neither
# is bound in `ini()` or `model()`, so both stay out of scope by construction
# rather than by an exception list that would rot.
.fmFamilyPattern <- "^fm_"

# Symbols on the left of an assignment anywhere in one expression. Recursive
# because rxode2 permits `if (...) { ... }` inside model(), and a name bound in
# a branch is bound just the same -- a flat scan of the top-level statements
# would miss it silently, which is the failure mode this whole check exists to
# remove.
.assignedNamesIn <- function(e, acc = character(0)) {
  if (!is.call(e)) return(acc)
  op <- tryCatch(as.character(e[[1]]), error = function(e) character(0))
  if (length(op) && op[[1]] %in% c("<-", "=") && length(e) >= 2L &&
      is.name(e[[2]])) {
    acc <- c(acc, as.character(e[[2]]))
  }
  for (i in seq_along(e)) {
    part <- tryCatch(e[[i]], error = function(e) NULL)
    if (is.call(part)) acc <- .assignedNamesIn(part, acc)
  }
  acc
}

# Assignment targets in the model block, as symbols rather than text. Sibling
# of `.modelBlockLines()`, which deparses; going through `ui$lstExpr` keeps
# comments and string literals out without having to strip them by regex.
.modelBlockAssignedNames <- function(ui) {
  exprs <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (!length(exprs)) return(character(0))
  out <- character(0)
  for (e in exprs) out <- .assignedNamesIn(e, out)
  unique(out)
}

# Nearest registered member, for the two drifts this check exists to stop:
# a case difference (`fm_H4`) and a stray plural (`fm_others`).
.fmFamilyNearMiss <- function(nm, registered) {
  if (!length(registered)) return(NA_character_)
  hit <- registered[tolower(registered) == tolower(nm)]
  if (length(hit)) return(hit[[1]])
  hit <- registered[tolower(registered) == sub("s$", "", tolower(nm))]
  if (length(hit)) return(hit[[1]])
  hit <- registered[paste0(tolower(registered), "s") == tolower(nm)]
  if (length(hit)) return(hit[[1]])
  NA_character_
}

.checkFmFamily <- function(ui, conv) {
  issues <- .emptyIssues()
  registered <- grep(.fmFamilyPattern, conv$paperNamedParams, value = TRUE)
  bound <- character(0)
  ini <- ui$iniDf
  if (!is.null(ini) && nrow(ini) > 0) bound <- c(bound, ini$name)
  bound <- unique(c(bound, .modelBlockAssignedNames(ui)))
  unregistered <- setdiff(grep(.fmFamilyPattern, bound, value = TRUE), registered)
  for (nm in unregistered) {
    near <- .fmFamilyNearMiss(nm, registered)
    suggestion <-
      if (is.na(near)) {
        paste("If this is a genuinely new elimination pathway, add it to the",
              "`### fm_...` family heading in inst/references/parameter-names.md;",
              "`<pathway>` must be lowercase and must match the enzyme or the",
              "canonical metabolite suffix it names. If the parameter is not a",
              "fraction metabolised, rename it so it does not claim the `fm_`",
              "prefix.")
      } else {
        sprintf(paste("Rename to '%s', the registered spelling. If this is",
                      "instead a distinct pathway that merely looks similar,",
                      "add it to the `### fm_...` family heading in",
                      "inst/references/parameter-names.md."), near)
      }
    issues <- rbind(issues, .issue(
      "fm_family", "error", nm,
      sprintf(paste("'%s' is not a registered `fm_<pathway>` fraction-metabolised",
                    "parameter (registered: %s)."),
              nm, paste(registered, collapse = ", ")),
      suggestion
    ))
  }
  issues
}

# Reference-register duplicate detection ------------------------------------
#
# A canonical name may legitimately appear twice in a register when the two
# entries carry DIFFERENT `- **Type:**` values. `col`, `complex`, `dap`,
# `lzd`, `mer`, `mero`, `plasma` and `van` are each registered both as a bare
# compartment and as a metabolite-suffix, and several say so in their Notes.
#
# Two entries sharing a name AND a Type are a different matter. The registers
# are resolved in document order, last one wins, so the earlier block -- and
# any source alias, example model or caveat recorded only there -- is silently
# discarded. That is not hypothetical: `col` (twice as metabolite-suffix) and
# `mic` (twice as paper-named-param) were each written by two extractions that
# did not know the other existed, and `cloca` was duplicated again by a
# consolidation merge. In every case text was being dropped and nobody noticed,
# because nothing consumed the fact.
.referenceHeaderPattern <- "^### ([^ (]+)\\s*(\\(\\*\\*.*\\*\\*\\))?\\s*$"
.referenceTypePattern <- "^-\\s+\\*\\*Type:\\*\\*\\s*(.+?)\\s*$"

.referenceRegisterBlocks <- function(path) {
  empty <- data.frame(name = character(), type = character(),
                      line = integer(), stringsAsFactors = FALSE)
  lines <- readLines(path, warn = FALSE)
  isHeader <- grepl(.referenceHeaderPattern, lines)
  if (!any(isHeader)) return(empty)
  # A block ends at the next register header, the next `## ` section, or a
  # horizontal rule -- whichever comes first.
  isBreak <- isHeader | grepl("^## ", lines) | grepl("^---\\s*$", lines)
  hdr <- which(isHeader)
  rows <- lapply(hdr, function(i) {
    nm <- sub(.referenceHeaderPattern, "\\1", lines[[i]])
    nxt <- which(isBreak & seq_along(lines) > i)
    end <- if (length(nxt)) nxt[[1]] - 1L else length(lines)
    body <- if (end >= i + 1L) lines[(i + 1L):end] else character()
    tl <- grep(.referenceTypePattern, body, value = TRUE)
    ty <- if (length(tl)) sub(.referenceTypePattern, "\\1", tl[[1]]) else NA_character_
    data.frame(name = nm, type = ty, line = i, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.referenceDuplicateIssues <- function(paths) {
  issues <- .emptyIssues()
  for (p in paths) {
    if (!file.exists(p)) next
    blocks <- .referenceRegisterBlocks(p)
    if (!nrow(blocks)) next
    key <- paste(blocks$name, blocks$type, sep = "\r")
    for (k in unique(key[duplicated(key)])) {
      rows <- blocks[key == k, , drop = FALSE]
      issues <- rbind(issues, .issue(
        "reference_duplicate", "error", rows$name[[1]],
        sprintf(
          "%s has %d `### %s` entries sharing `- **Type:** %s` (lines %s).",
          basename(p), nrow(rows), rows$name[[1]], rows$type[[1]],
          paste(rows$line, collapse = ", ")),
        paste("Merge them into a single block, keeping every source alias and",
              "example model from both. The register resolves in document",
              "order, last one wins, so each earlier block is silently",
              "discarded. Repeating a name is only safe when the Types differ.")
      ))
    }
  }
  issues
}
