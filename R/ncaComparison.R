#' Map PKNCA parameter codes to friendly human-readable labels
#'
#' Vectorized recoder that maps the PKNCA parameter codes (`PPTESTCD`
#' values such as `"cmax"`, `"aucinf.obs"`, `"half.life"`) to the friendly
#' pharmacometric labels readers expect in tables and prose (`"Cmax"`,
#' `"AUC0-∞ (obs)"`, `"t½"`). Codes that are not in the mapping
#' table are returned unchanged with a warning, so unfamiliar parameters
#' never get silently dropped.
#'
#' @param code Character vector of PKNCA parameter codes. `NA` values are
#'   passed through unchanged.
#' @param units Optional named character vector keyed by PKNCA code. For
#'   each code with a matching entry, the friendly label is suffixed with
#'   the unit string in parentheses (e.g. `"Cmax (ng/mL)"`).
#' @return Character vector the same length as `code`.
#' @export
#' @author Bill Denney
#' @examples
#' ncaParamLabel(c("cmax", "tmax", "auclast", "aucinf.obs", "half.life"))
#' ncaParamLabel(
#'   c("cmax", "auclast"),
#'   units = c(cmax = "ng/mL", auclast = "ng*h/mL")
#' )
ncaParamLabel <- function(code, units = NULL) {
  checkmate::assertCharacter(code)
  if (!is.null(units)) {
    checkmate::assertCharacter(units, names = "named", any.missing = FALSE)
  }
  tbl <- .ncaParamLabelTable()
  out <- unname(tbl[code])
  unknown <- which(!is.na(code) & is.na(out))
  if (length(unknown) > 0L) {
    warning(
      "ncaParamLabel(): unknown PKNCA code(s) returned as-is: ",
      paste(shQuote(unique(code[unknown])), collapse = ", "),
      call. = FALSE
    )
    out[unknown] <- code[unknown]
  }
  if (!is.null(units)) {
    has_units <- !is.na(code) & code %in% names(units)
    out[has_units] <- paste0(out[has_units], " (", units[code[has_units]], ")")
  }
  out
}

.ncaParamLabelTable <- function() {
  c(
    cmax           = "Cmax",
    cmin           = "Cmin",
    tmax           = "Tmax",
    tlast          = "Tlast",
    clast.obs      = "Clast",
    aucinf.obs     = "AUC0-∞ (obs)",
    aucinf.pred    = "AUC0-∞ (pred)",
    auclast        = "AUClast",
    aucall         = "AUCall",
    half.life      = "t½",
    lambda.z       = "λz",
    lambda.z.n.points = "λz n points",
    cl.obs         = "CL/F",
    cl.pred        = "CL/F",
    vss.obs        = "Vss/F",
    vss.pred       = "Vss/F",
    vz.obs         = "Vz/F",
    vz.pred        = "Vz/F",
    cav            = "Cavg",
    ctau           = "Cτ",
    ctrough        = "Ctrough",
    accumulation   = "Rac",
    swing          = "Swing",
    swing_over     = "Swing",
    fluctuation    = "Fluctuation",
    fluctuation_pct = "Fluctuation (%)",
    ptr            = "PTR",
    aumcinf.obs    = "AUMC0-∞ (obs)",
    aumcinf.pred   = "AUMC0-∞ (pred)",
    aumclast       = "AUMC0-t",
    mrtinf.obs     = "MRT0-∞ (obs)",
    mrtinf.pred    = "MRT0-∞ (pred)",
    mrtlast        = "MRT0-t",
    mrt.iv.obs     = "MRTiv (obs)",
    mrt.iv.pred    = "MRTiv (pred)",
    thalf.eff.obs  = "t½,eff",
    f              = "F",
    n.samples      = "N samples",
    ae             = "Ae (amount excreted)",
    fe             = "Fe (fraction excreted)",
    clr.last       = "CLr (last)",
    clr.obs        = "CLr (obs)",
    clr.pred       = "CLr (pred)"
  )
}

#' Build a side-by-side NCA comparison table (simulated vs. reference)
#'
#' Combines model-predicted NCA results with a reference set (e.g. paper
#' values, prior estimates, alternative model runs, observed data) into a
#' single tidy comparison frame ready for `knitr::kable()` or interactive
#' review. Each row is one NCA parameter \eqn{\times} optional grouping
#' level; columns show the reference value, simulated value, and percent
#' difference. Rows whose discrepancy exceeds `tolerance_pct` are flagged
#' with a trailing asterisk and the function attaches a footnote string
#' to the result.
#'
#' Designed for any nlmixr2 user comparing predicted NCA to a reference
#' \U2014 not solely a vignette utility. Accepts a `PKNCAresults` object,
#' the `$result` data.frame, or a tidy long frame with columns `PPTESTCD`
#' and `PPORRES` for `simulated`; the same shapes plus a wide
#' `data.frame` (one row per group, one column per PKNCA parameter code)
#' for `reference`. Multiple rows per group in `simulated` (the typical
#' per-subject case from [PKNCA::pk.nca()]) are aggregated to the group
#' level using [stats::median()] before joining; pass a pre-aggregated
#' frame if a different summary statistic is required.
#'
#' @param simulated Model-predicted NCA. One of: a `PKNCAresults` object,
#'   a `data.frame` with columns `PPTESTCD` and `PPORRES` plus any
#'   grouping columns named in `by`, or a wide `data.frame` (one row per
#'   group, one column per PKNCA code).
#' @param reference Reference NCA. Same shapes as `simulated`. A wide
#'   `data.frame` is the typical form when transcribing values from a
#'   paper table.
#' @param by Optional character vector of grouping column name(s) (e.g.
#'   `"treatment"`, `c("cohort","weight_band")`). When `NULL`, the inputs
#'   are treated as a single ungrouped comparison.
#' @param params Optional character vector of PKNCA codes restricting
#'   which parameters appear in the output. Defaults to the intersection
#'   of codes present in both inputs.
#' @param tolerance_pct Numeric scalar; rows with
#'   `|% diff| > tolerance_pct` get a trailing `*` and the function
#'   attaches a footnote string via the `"footnote"` attribute on the
#'   result. Pass `Inf` to disable flagging.
#' @param units Optional named character vector keyed by PKNCA code,
#'   passed through to [ncaParamLabel()].
#' @param label_first_column Header for the parameter column. Defaults to
#'   `"NCA parameter"`.
#' @return A `data.frame` whose first column is the friendly parameter
#'   label (under `label_first_column`); when `by` is non-NULL, the
#'   grouping column(s) follow; the final three columns are `Reference`,
#'   `Simulated`, and `% diff`. Rows are sorted by canonical PKNCA
#'   parameter order then by group. Carries a `"footnote"` attribute when
#'   any row exceeds the tolerance.
#' @export
#' @author Bill Denney
#' @examples
#' simulated <- data.frame(
#'   treatment = rep(c("50 mg", "100 mg"), each = 3),
#'   PPTESTCD  = rep(c("cmax", "tmax", "auclast"), 2),
#'   PPORRES   = c(15.2, 2.0, 96.0, 29.1, 2.1, 191.0)
#' )
#' reference <- data.frame(
#'   treatment = c("50 mg", "100 mg"),
#'   cmax      = c(14.8, 28.5),
#'   tmax      = c(2.0,  2.1),
#'   auclast   = c(95.0, 190.0)
#' )
#' tbl <- ncaComparisonTable(
#'   simulated, reference,
#'   by = "treatment",
#'   units = c(cmax = "ug/mL", auclast = "ug*h/mL", tmax = "h")
#' )
#' tbl
#' attr(tbl, "footnote")
ncaComparisonTable <- function(simulated, reference,
                               by = NULL,
                               params = NULL,
                               tolerance_pct = 20,
                               units = NULL,
                               label_first_column = "NCA parameter") {
  checkmate::assertNumeric(
    tolerance_pct, len = 1L, lower = 0, any.missing = FALSE
  )
  checkmate::assertCharacter(label_first_column, len = 1L, any.missing = FALSE)
  if (!is.null(by)) {
    checkmate::assertCharacter(by, any.missing = FALSE, min.len = 1L)
  }
  if (!is.null(params)) {
    checkmate::assertCharacter(params, any.missing = FALSE, min.len = 1L)
  }

  sim_long <- .toLongNca(simulated, by, value_name = "Simulated")
  ref_long <- .toLongNca(reference, by, value_name = "Reference")

  group_order <- if (!is.null(by)) {
    setNames(
      lapply(by, function(b) unique(ref_long[[b]])),
      by
    )
  } else NULL

  if (!is.null(params)) {
    sim_long <- sim_long[sim_long$PPTESTCD %in% params, , drop = FALSE]
    ref_long <- ref_long[ref_long$PPTESTCD %in% params, , drop = FALSE]
  }

  group_cols <- c(by, "PPTESTCD")
  sim_long <- .aggregateMedian(sim_long, group_cols, value_col = "Simulated")
  ref_long <- .aggregateMedian(ref_long, group_cols, value_col = "Reference")

  joined <- merge(
    ref_long, sim_long,
    by = group_cols, all = FALSE, sort = FALSE
  )
  if (nrow(joined) == 0L) {
    stop(
      "No NCA parameters overlap between `simulated` and `reference`.",
      " Check `by`, `params`, and column-name conventions.",
      call. = FALSE
    )
  }

  joined$`% diff` <- ifelse(
    is.na(joined$Reference) | joined$Reference == 0,
    NA_real_,
    (joined$Simulated - joined$Reference) / joined$Reference * 100
  )
  flagged <- !is.na(joined$`% diff`) & abs(joined$`% diff`) > tolerance_pct

  ref_fmt <- .fmtNcaValue(joined$Reference)
  sim_fmt <- .fmtNcaValue(joined$Simulated)
  pct_fmt <- .fmtNcaPct(joined$`% diff`, flagged)

  canon <- names(.ncaParamLabelTable())
  sort_idx <- match(joined$PPTESTCD, canon)
  sort_idx[is.na(sort_idx)] <- length(canon) + 1L
  if (!is.null(by)) {
    group_idx <- vapply(seq_len(nrow(joined)), function(i) {
      sum(vapply(seq_along(by), function(j) {
        match(joined[[by[j]]][i], group_order[[by[j]]]) *
          10L^(length(by) - j)
      }, numeric(1L)))
    }, numeric(1L))
    ord <- order(sort_idx, joined$PPTESTCD, group_idx)
  } else {
    ord <- order(sort_idx, joined$PPTESTCD)
  }

  out <- data.frame(
    .label = ncaParamLabel(joined$PPTESTCD, units = units),
    stringsAsFactors = FALSE
  )
  names(out) <- label_first_column
  if (!is.null(by)) {
    for (bcol in by) out[[bcol]] <- joined[[bcol]]
  }
  out$Reference <- ref_fmt
  out$Simulated <- sim_fmt
  out$`% diff`  <- pct_fmt
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  if (any(flagged)) {
    attr(out, "footnote") <- sprintf(
      "* differs from reference by more than ±%g%%.",
      tolerance_pct
    )
  }
  out
}

.toLongNca <- function(x, by, value_name) {
  if (inherits(x, "PKNCAresults")) {
    x <- as.data.frame(x$result)
  }
  if (!is.data.frame(x)) {
    stop(
      "`", value_name,
      "` must be a data.frame, tibble, or PKNCAresults.",
      call. = FALSE
    )
  }
  x <- as.data.frame(x)
  has_long_cols <- all(c("PPTESTCD", "PPORRES") %in% names(x))
  if (has_long_cols) {
    keep <- intersect(c(by, "PPTESTCD", "PPORRES"), names(x))
    out <- x[, keep, drop = FALSE]
    names(out)[match("PPORRES", names(out))] <- value_name
    return(out)
  }
  by_present <- intersect(by, names(x))
  miss_by <- setdiff(by, names(x))
  if (length(miss_by) > 0L) {
    stop(
      "`", value_name, "` is missing grouping column(s): ",
      paste(miss_by, collapse = ", "),
      call. = FALSE
    )
  }
  value_cols <- setdiff(names(x), by_present)
  if (length(value_cols) == 0L) {
    stop(
      "`", value_name, "` (wide form) has no parameter columns.",
      call. = FALSE
    )
  }
  parts <- lapply(value_cols, function(col) {
    df <- x[, c(by_present, col), drop = FALSE]
    names(df)[length(names(df))] <- value_name
    df$PPTESTCD <- col
    df[, c(by_present, "PPTESTCD", value_name), drop = FALSE]
  })
  do.call(rbind, parts)
}

.aggregateMedian <- function(df, group_cols, value_col) {
  if (nrow(df) == 0L) return(df)
  if (anyDuplicated(df[, group_cols, drop = FALSE]) == 0L) {
    return(df)
  }
  agg_fmla <- stats::as.formula(
    paste0("`", value_col, "` ~ ",
           paste0("`", group_cols, "`", collapse = " + "))
  )
  out <- stats::aggregate(
    agg_fmla, data = df,
    FUN = function(v) stats::median(v, na.rm = TRUE)
  )
  out
}

.fmtNcaValue <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) "—"
    else format(signif(v, 3), scientific = FALSE, trim = TRUE)
  }, character(1L))
}

.fmtNcaPct <- function(x, flagged) {
  out <- vapply(x, function(v) {
    if (is.na(v)) "—"
    else sprintf("%+0.1f%%", v)
  }, character(1L))
  out[flagged] <- paste0(out[flagged], "*")
  out
}
