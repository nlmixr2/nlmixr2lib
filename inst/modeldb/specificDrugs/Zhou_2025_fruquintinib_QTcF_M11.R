Zhou_2025_fruquintinib_QTcF_M11 <- function() {
  description <- paste(
    "Concentration-QTc (C-QTc) linear mixed-effects PD model relating the",
    "change from baseline in the Fridericia-corrected QT interval",
    "(DeltaQTcF, ms) to the plasma concentration of M11, the major",
    "metabolite of fruquintinib, in patients with previously treated",
    "metastatic colorectal cancer enrolled in the ECG substudy of the",
    "phase 3 FRESCO-2 trial (NCT04322539). This is the SUPPORTIVE",
    "analysis of Zhou 2025: the authors selected the population-based",
    "correction (QTcP) for the primary analysis because a graphical",
    "assessment showed Fridericia's formula did not adequately correct",
    "for heart rate in this cohort, and repeated the analysis on QTcF",
    "for cross-comparison with conventional thorough-QT analyses. The",
    "development trajectory was the same as for DeltaQTcP and the",
    "selected model has the same structure, differing only in the",
    "coefficient values and in the baseline centering reference",
    "(409.5 ms for QTcF vs 419.3 ms for QTcP). The model is:",
    "DeltaQTcF = e0 + etae0 + e_on_treatment_e0 * ON_TREATMENT",
    "+ (slope + etaslope) * CP_FRUQUINTINIB_M11_NGML",
    "+ sum_j(e_ntime<j>_e0 * [NTIME == j]) + e_day21_e0 * DAY21",
    "+ e_qtc_bl_e0 * (QTC_BL - 409.5),",
    "with independent (diagonal) between-subject variability on the",
    "intercept and on the M11 concentration slope, and additive residual",
    "error. The placebo-corrected contrast is mean DeltaDeltaQTcF =",
    "e_on_treatment_e0 + slope * CP_FRUQUINTINIB_M11_NGML.",
    "PD-only model: M11 plasma concentration is supplied as a",
    "time-varying covariate (ng/mL); placebo records carry 0. The source",
    "publication does not fit a population PK model, so users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "M11 concentration trajectory (no fruquintinib or M11 popPK model",
    "exists in the nlmixr2lib registry). Companion model files:",
    "Zhou_2025_fruquintinib_QTcP_M11.R (the paper's final model) and",
    "Zhou_2025_fruquintinib_QTcP_parent.R (the additional analysis driven",
    "by parent fruquintinib concentration).",
    sep = " "
  )

  reference <- paste(
    "Zhou X, Toms A, Morton D, Wang X, Taylor A, Dasari A, Gupta N,",
    "Chien C. Assessment of the Effects of Fruquintinib on Cardiac Safety",
    "in Patients with Metastatic Colorectal Cancer.",
    "The Journal of Clinical Pharmacology 2025;65(11):1411-1419.",
    "doi:10.1002/jcph.70051.",
    "Parameter values from Supplementary Table 4 of the same article.",
    sep = " "
  )

  vignette <- "Zhou_2025_fruquintinib"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external M11 plasma-concentration covariate)",
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF = QT / RR^(1/3) with RR in seconds, in ms; driving covariate CP_FRUQUINTINIB_M11_NGML is in ng/mL)"
  )

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Randomized treatment-arm indicator: 1 = fruquintinib 5 mg once daily, 0 = matching placebo. Time-fixed per subject (parallel-group randomized design).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = paste(
        "Source paper writes this as TRT_i, 'the fixed effect associated",
        "with treatment TRT_i received by patient i (TRT_i = 0 for placebo",
        "and 1 for active drug)' (Zhou 2025 Methods 'C-QTc Model'), and is",
        "explicit that the coefficient is empirical: 'the theta_1 term has",
        "no physiological interpretation but allows flexibility in the",
        "event of model specification'. It does NOT vanish at zero",
        "concentration: the placebo-corrected contrast at concentration C",
        "is e_on_treatment_e0 + slope * C.",
        "All placebo records were assigned a fruquintinib and M11",
        "concentration of 0 at each nominal time point (Zhou 2025 Methods",
        "'Overview of Data')."
      ),
      source_name        = "Trt (paper notation TRT_i)"
    ),
    CP_FRUQUINTINIB_M11_NGML = list(
      description        = "Instantaneous total plasma concentration of M11, the major metabolite of fruquintinib, at the time of each time-matched ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. The slope is reported directly in ms",
        "per ng/mL (Zhou 2025 Supplementary Table 4 units column), so no",
        "in-model unit rescaling is required.",
        "Observed M11 plasma concentration drawn at the same NOMINAL time",
        "point as the triplicate Holter-extracted 12-lead ECG; pairs whose",
        "ECG and PK sampling were more than 30 minutes apart were",
        "excluded. Placebo records carry 0. 326 of 936 (34.8%) M11",
        "concentrations in the fruquintinib arm were below the limit of",
        "quantification and were set to zero and retained.",
        "Reference values observed: the geometric mean steady-state Cmax",
        "of M11 after fruquintinib 5 mg QD is 77 ng/mL and twice that is",
        "154 ng/mL; these are the two scenarios tabulated in Zhou 2025",
        "Table 4 for this model. The QTcF model predicts the upper bound",
        "of the 90% CI of mean DeltaDeltaQTcF not to exceed 10 ms at M11",
        "concentrations up to 177 ng/mL, 2.3-fold the observed",
        "steady-state geometric mean Cmax (Zhou 2025 Results 'Supportive",
        "Analyses').",
        "Set to 0 for placebo records and outside the drug-exposure",
        "window; the concentration-slope term then collapses to 0."
      ),
      source_name        = "CONC2 (supplement model-development logs); 'M11 Conc.' (Supplementary Table 4)"
    ),
    NTIME = list(
      description        = "Protocol-scheduled (nominal) time in whole hours after the most recent dose at which the time-matched ECG / PK pair was collected. Takes the values 0, 1, 2, 3 and 4 h in this study; 0 h is the reference level.",
      units              = "h",
      type               = "categorical",
      reference_category = "0 (pre-dose nominal time point)",
      notes              = paste(
        "Decomposed inside model() into mutually-exclusive binary",
        "indicators nt1 .. nt4 multiplied by the per-timepoint fixed",
        "effects e_ntime1_e0 .. e_ntime4_e0. NTIME = 0 sets every",
        "indicator to zero and recovers the reference intercept.",
        "Source paper: 'theta_3 is the fixed effect associated with",
        "nominal time (one value estimated nonreference per nominal time",
        "point)' (Zhou 2025 Methods 'C-QTc Model'); the estimated rows are",
        "labelled 'NTime=1' through 'NTime=4' with units 'msec' in",
        "Supplementary Table 4, whose footnote defines 'NTime, nominal",
        "time in hours'.",
        "The paper interprets these coefficients as drug-independent",
        "diurnal variation, so they are common to both arms and cancel",
        "out of the placebo-corrected DeltaDeltaQTc contrast.",
        "NOMINAL, not actual, time: populate this column from the protocol",
        "schedule, not from rxode2's tad(). Because the model()",
        "decomposition tests exact equality, supply whole-number hours."
      ),
      source_name        = "NTIME / NTime"
    ),
    DAY21 = list(
      description        = "Cycle 1 day 21 visit indicator: 1 = the observation was collected at the cycle 1 day 21 (steady-state) ECG visit, 0 = the observation was collected at the cycle 1 day 1 (first-dose) visit. Time-varying within subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cycle 1 day 1 visit)",
      notes              = paste(
        "Source paper: 'theta_4 is the fixed effect associated with cycle",
        "1 day 21 (taking cycle 1 day 1 as the reference visit)' (Zhou",
        "2025 Methods 'C-QTc Model'); tabulated as the 'Cycle 1 Day 21",
        "Visit' row of Supplementary Table 4.",
        "Fruquintinib was given on days 1-21 of a 28-day cycle, so the",
        "cycle 1 day 1 visit captures the first dose and the cycle 1 day",
        "21 visit captures steady state after 21 consecutive daily doses.",
        "Unlike the parent-concentration model",
        "(Zhou_2025_fruquintinib_QTcP_parent.R), which drops this term,",
        "the QTcF/M11 model retains it (Supplementary Table 4 reports a",
        "significant estimate, P = 0.0214)."
      ),
      source_name        = "visit (paper notation Visit_k, k = cycle 1 day 1 or 21)"
    ),
    QTC_BL = list(
      description        = "Subject's pre-dose (cycle 1 day 1) baseline Fridericia-corrected QT interval (QTcF), treated as a per-subject time-fixed covariate. Enters the linear-mixed-effects intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 409.5). Set QTC_BL = 409.5 ms for the typical (mean-baseline) subject -- the centered term then collapses to 0.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject and MEAN-centered (Zhou 2025 Methods",
        "'C-QTc Model': 'QTc0 is the overall mean baseline QTc'), not",
        "median-centered as in Darpo 2014.",
        "Centering reference: reported directly by the source -- 'The mean",
        "baseline QTcP and QTcF in the analysis dataset were 419.3 and",
        "409.5 ms, respectively' (Zhou 2025 Results 'Summary of Data Used",
        "for Analysis'). This file uses qtc_bl_ref = 409.5 ms, the QTcF",
        "value; it is a reported value, not an assumption.",
        "Correction method: Fridericia, QTcF = QT / RR^(1/3) with RR in",
        "seconds, computed per replicate (Zhou 2025 Methods 'QTcF and",
        "QTcP'). The authors note this correction was inadequate in this",
        "cohort -- the baseline QTcF-RR slope was 0.0493 (90% CI 0.0393,",
        "0.0592), statistically significant and much larger in magnitude",
        "than the QTcP-RR slope of -3.73e-05 (90% CI -0.0102, 0.0101) --",
        "which is why QTcP was chosen for the primary analysis and this",
        "model is supportive (Zhou 2025 Discussion and Supplementary",
        "Figure 1).",
        "For the population-based-corrected companion model the baseline",
        "is the QTcP baseline with a 419.3 ms centering reference; see",
        "Zhou_2025_fruquintinib_QTcP_M11.R."
      ),
      source_name        = "baseline QTcF (paper notation QTc0,i)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 205L,
    n_studies        = 1L,
    n_observations   = 1456L,
    age_range        = NA_character_,
    weight_range     = NA_character_,
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = paste(
      "Adults with previously treated metastatic colorectal cancer (mCRC)",
      "enrolled in the phase 3 FRESCO-2 trial (NCT04322539) and included",
      "in the cardiovascular-safety ECG subset. The study excluded",
      "patients with a QTcF > 480 ms, patients with any other factor that",
      "could prolong the QTc interval or increase arrhythmic risk, and",
      "patients receiving concomitant QTc-prolonging medication (Zhou 2025",
      "Methods 'Overview of Data')."
    ),
    dose_range       = paste(
      "Fruquintinib 5 mg orally once daily on days 1-21 of a 28-day cycle",
      "plus best supportive care, or matching placebo plus best supportive",
      "care. Time-matched ECG / PK pairs at nominal times 0, 1, 2, 3 and",
      "4 h after dose at the cycle 1 day 1 and cycle 1 day 21 visits."
    ),
    regions          = NA_character_,
    notes            = paste(
      "The QTcF analysis uses the same 205-patient C-QTc analysis set as",
      "the primary QTcP analysis (137 fruquintinib, 68 placebo; 1456",
      "time-matched pairs). Mean baseline QTcF 409.5 ms.",
      "The QTcF model is the SUPPORTIVE analysis. Zhou 2025 Results",
      "'Supportive Analyses': 'The development trajectory for the C-QTcF",
      "model was similar to that of DeltaQTcP, and the selected model had",
      "only M11 concentrations with the same structure as that of the",
      "model developed using DeltaQTcP.'",
      "Analysis software: R 4.2.0 with the nlme package (version 3.1-157)",
      "and the contrast package (version 0.24.2)."
    )
  )

  ini({
    # ==================================================================
    # Linear mixed-effects C-QTc model, Fridericia-corrected endpoint
    # (Zhou 2025 Supplementary Table 4, 'Final Parameter Estimates for
    # the C-QTc Model with DeltaQTcF'). Structure identical to the QTcP
    # final model (Zhou 2025 Results 'Supportive Analyses': "the selected
    # model had only M11 concentrations with the same structure as that
    # of the model developed using DeltaQTcP"), so the between-subject
    # covariance matrix is likewise DIAGONAL.
    #
    # Every parameter is on the LINEAR DeltaQTcF (ms) scale and every
    # random effect is ADDITIVE on that scale, matching the source nlme
    # fit. See Zhou_2025_fruquintinib_QTcP_M11.R for why neither the
    # intercept nor the slope is log-transformed.
    # ==================================================================

    e0 <- 2.37
    label("Linear-mixed-effects intercept on DeltaQTcF, placebo / cycle 1 day 1 / 0 h nominal time / mean-baseline reference (ms)")
    # Zhou 2025 Supplementary Table 4 'Intercept' = 2.37 (SE 1.31;
    # RSE 55.3%; 95% CI -0.182, 4.93; P = 0.0695).

    slope <- 0.0477
    label("Linear M11-concentration slope on DeltaQTcF (ms per ng/mL)")
    # Zhou 2025 Supplementary Table 4 'M11 Conc. Slope' = 0.0477
    # (SE 0.0163; RSE 34.2%; 95% CI 0.0158, 0.0795; P = 0.0035).
    # Reported directly in msec per ng/mL, so no in-model rescaling.

    # ------------------------------------------------------------------
    # Covariate effects (Zhou 2025 Supplementary Table 4). All estimated.
    # ------------------------------------------------------------------

    e_on_treatment_e0 <- -3.34
    label("Effect of active-treatment arm on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'Treatment' = -3.34 (SE 1.32;
    # RSE 39.5%; 95% CI -5.94, -0.743; P = 0.0122). Empirical
    # flexibility term; with `slope` it defines the placebo-corrected
    # contrast mean DeltaDeltaQTcF = -3.34 + 0.0477 * C.

    e_ntime1_e0 <- 0.964
    label("Effect of the 1 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'NTime=1' = 0.964 (SE 0.826;
    # RSE 85.7%; 95% CI -0.651, 2.58; P = 0.2432).

    e_ntime2_e0 <- 2.73
    label("Effect of the 2 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'NTime=2' = 2.73 (SE 0.826;
    # RSE 30.3%; 95% CI 1.12, 4.35; P = 0.0010).

    e_ntime3_e0 <- 2.69
    label("Effect of the 3 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'NTime=3' = 2.69 (SE 0.826;
    # RSE 30.7%; 95% CI 1.07, 4.31; P = 0.0012).

    e_ntime4_e0 <- 1.42
    label("Effect of the 4 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'NTime=4' = 1.42 (SE 0.833;
    # RSE 58.7%; 95% CI -0.210, 3.05; P = 0.0887).

    e_day21_e0 <- -1.50
    label("Effect of the cycle 1 day 21 visit on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 4 'Cycle 1 Day 21 Visit' = -1.50
    # (SE 0.652; RSE 43.5%; 95% CI -2.78, -0.226; P = 0.0214).

    e_qtc_bl_e0 <- -0.118
    label("Effect of mean-centered baseline QTcF on the intercept (additive, ms per ms)")
    # Zhou 2025 Supplementary Table 4 'Baseline QTcF' = -0.118
    # (SE 0.0286; RSE 24.2%; 95% CI -0.174, -0.0621; P < 0.0001).

    # ==================================================================
    # Between-subject variability. Supplementary Table 4 reports SDs on
    # the linear DeltaQTcF scale; ini() takes variances, so each is
    # squared. Diagonal (independent) covariance structure, matching the
    # QTcP final model.
    #
    # These two values (8.02 and 0.116) are the ones misquoted in the
    # main text's QTcP paragraph (Zhou 2025 Results 'C-QTcP Model and
    # Evaluation'); they belong here, to the QTcF model, not to Table 2.
    # See Zhou_2025_fruquintinib_QTcP_M11.R and the vignette Errata.
    # ==================================================================

    etae0 ~ 64.3204
    # Zhou 2025 Supplementary Table 4 'BSV SD for Intercept' = 8.02 msec
    # (SE 0.473; RSE 5.90%; 95% CI 7.14, 9.00). Variance = 8.02^2
    # = 64.3204.

    etaslope ~ 0.013456
    # Zhou 2025 Supplementary Table 4 'BSV SD of Conc. Slope' = 0.116
    # msec per ng/mL (SE 0.013; RSE 11.2%; 95% CI 0.0938, 0.145).
    # Variance = 0.116^2 = 0.013456.

    # ==================================================================
    # Residual error: additive on the DeltaQTcF scale (Zhou 2025 Methods
    # 'C-QTc Model': residuals normal with mean 0 and variance R).
    # ==================================================================
    addSd <- 7.54
    label("Additive residual error standard deviation on DeltaQTcF (ms)")
    # Zhou 2025 Supplementary Table 4 'Residual Error SD' = 7.54 msec
    # (SE 0.157; RSE 2.08%; 95% CI 7.24, 7.86).
  })

  model({
    # ==================================================================
    # 1. Nominal-time decomposition (reference level NTIME = 0).
    # ==================================================================
    nt1 <- (NTIME == 1)
    nt2 <- (NTIME == 2)
    nt3 <- (NTIME == 3)
    nt4 <- (NTIME == 4)
    ntime_effect <-
      e_ntime1_e0 * nt1 + e_ntime2_e0 * nt2 +
      e_ntime3_e0 * nt3 + e_ntime4_e0 * nt4

    # ==================================================================
    # 2. Mean-centered baseline QTcF. The centering reference is
    #    REPORTED: mean baseline QTcF = 409.5 ms (Zhou 2025 Results
    #    'Summary of Data Used for Analysis').
    # ==================================================================
    qtc_bl_ref      <- 409.5
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    # ==================================================================
    # 3. Per-subject intercept and slope; both random effects additive
    #    on the linear DeltaQTcF scale.
    # ==================================================================
    e0_i <-
      e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      ntime_effect +
      e_day21_e0 * DAY21 +
      e_qtc_bl_e0 * qtc_bl_centered
    slope_i <- slope + etaslope

    # ==================================================================
    # 4. Prediction: change from baseline in the Fridericia-corrected QT
    #    interval (DeltaQTcF), in ms.
    # ==================================================================
    QTcF <- e0_i + slope_i * CP_FRUQUINTINIB_M11_NGML

    QTcF ~ add(addSd)
  })
}
