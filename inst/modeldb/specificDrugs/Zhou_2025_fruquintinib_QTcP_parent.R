Zhou_2025_fruquintinib_QTcP_parent <- function() {
  description <- paste(
    "Concentration-QTc (C-QTc) linear mixed-effects PD model relating the",
    "change from baseline in the population-based corrected QT interval",
    "(DeltaQTcP, ms) to the plasma concentration of PARENT fruquintinib",
    "in patients with previously treated metastatic colorectal cancer",
    "enrolled in the ECG substudy of the phase 3 FRESCO-2 trial",
    "(NCT04322539). This is the ADDITIONAL analysis of Zhou 2025,",
    "generated because 34.8% of M11 concentrations in the analysis set",
    "were below the limit of quantification; it showed NO relationship",
    "between DeltaQTcP and fruquintinib concentration, with a slope of",
    "0.00778 ms per ng/mL (95% CI -0.00247, 0.018; P = 0.1377). It is a",
    "genuine negative result and is packaged so that the paper's full",
    "analysis set is reproducible.",
    "Two structural differences from the two M11-driven companion models",
    "are load-bearing and follow the authors' own model selection:",
    "(1) the cycle 1 day 21 VISIT term is ABSENT -- the authors removed",
    "it because its 95% CI included 0 and the likelihood-ratio-test",
    "p-value was 0.1318 (Supplementary Table 1, model fit.f.2);",
    "(2) the between-subject variability uses a FULL (correlated)",
    "covariance matrix, with an estimated intercept-slope correlation of",
    "-0.340, rather than the diagonal matrix of the M11 models. The",
    "model is:",
    "DeltaQTcP = e0 + etae0 + e_on_treatment_e0 * ON_TREATMENT",
    "+ (slope + etaslope) * CP_FRUQUINTINIB_NGML",
    "+ sum_j(e_ntime<j>_e0 * [NTIME == j])",
    "+ e_qtc_bl_e0 * (QTC_BL - 419.3),",
    "with additive residual error. The placebo-corrected contrast is",
    "mean DeltaDeltaQTcP = e_on_treatment_e0 + slope *",
    "CP_FRUQUINTINIB_NGML.",
    "PD-only model: fruquintinib plasma concentration is supplied as a",
    "time-varying covariate (ng/mL); placebo records carry 0. The source",
    "publication does not fit a population PK model, so users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "concentration trajectory (no fruquintinib popPK model exists in the",
    "nlmixr2lib registry). Companion model files:",
    "Zhou_2025_fruquintinib_QTcP_M11.R (the paper's final model) and",
    "Zhou_2025_fruquintinib_QTcF_M11.R (the supportive QTcF analysis).",
    sep = " "
  )

  reference <- paste(
    "Zhou X, Toms A, Morton D, Wang X, Taylor A, Dasari A, Gupta N,",
    "Chien C. Assessment of the Effects of Fruquintinib on Cardiac Safety",
    "in Patients with Metastatic Colorectal Cancer.",
    "The Journal of Clinical Pharmacology 2025;65(11):1411-1419.",
    "doi:10.1002/jcph.70051.",
    "Parameter values from Supplementary Table 5 of the same article.",
    sep = " "
  )

  vignette <- "Zhou_2025_fruquintinib"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external fruquintinib plasma-concentration covariate)",
    concentration = "(observation QTcP is the CHANGE FROM BASELINE in the population-based corrected QT interval, DeltaQTcP, in ms; driving covariate CP_FRUQUINTINIB_NGML is in ng/mL)"
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
        "event of model specification'.",
        "In this model the treatment offset carries essentially all of the",
        "apparent arm difference (-5.24 ms), while the concentration slope",
        "is not statistically significant -- which is exactly the paper's",
        "conclusion that DeltaQTcP has no relationship with parent",
        "fruquintinib concentration.",
        "All placebo records were assigned a fruquintinib concentration of",
        "0 at each nominal time point (Zhou 2025 Methods 'Overview of",
        "Data')."
      ),
      source_name        = "Trt (paper notation TRT_i)"
    ),
    CP_FRUQUINTINIB_NGML = list(
      description        = "Instantaneous total plasma concentration of parent fruquintinib at the time of each time-matched ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. The slope is reported directly in ms",
        "per ng/mL (Zhou 2025 Supplementary Table 5 units column), so no",
        "in-model unit rescaling is required.",
        "Observed fruquintinib plasma concentration drawn at the same",
        "NOMINAL time point as the triplicate Holter-extracted 12-lead",
        "ECG; pairs whose ECG and PK sampling were more than 30 minutes",
        "apart were excluded. Placebo records carry 0. Only 16 of 936",
        "(1.7%) fruquintinib concentrations in the active arm were below",
        "the limit of quantification (versus 34.8% for M11), which is why",
        "the authors ran this parent-driven analysis as a cross-check on",
        "the M11-driven final model (Zhou 2025 Results 'Summary of Data",
        "Used for Analysis' and Discussion).",
        "Reference values observed: the geometric mean steady-state Cmax",
        "of fruquintinib after 5 mg QD is approximately 290 ng/mL total",
        "(unbound approximately 13.6 ng/mL at 95.3% protein binding), and",
        "twice the geometric mean steady-state Cmax is 580 ng/mL -- the",
        "concentration at which the predicted upper bound of the 90% CI",
        "of mean DeltaDeltaQTcP was 3.96 ms (Zhou 2025 Discussion).",
        "Set to 0 for placebo records and outside the drug-exposure",
        "window; the concentration-slope term then collapses to 0."
      ),
      source_name        = "CONC1 (supplement model-development logs); 'Fruquintinib Conc.' (Supplementary Table 5)"
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
        "point)' (Zhou 2025 Methods 'C-QTc Model'); tabulated as 'NTime=1'",
        "through 'NTime=4' in Supplementary Table 5, whose footnote",
        "defines 'NTime, nominal time in hours'.",
        "NOMINAL, not actual, time: populate this column from the protocol",
        "schedule, not from rxode2's tad(). Because the model()",
        "decomposition tests exact equality, supply whole-number hours."
      ),
      source_name        = "NTIME / NTime"
    ),
    QTC_BL = list(
      description        = "Subject's pre-dose (cycle 1 day 1) baseline population-based corrected QT interval (QTcP), treated as a per-subject time-fixed covariate. Enters the linear-mixed-effects intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 419.3). Set QTC_BL = 419.3 ms for the typical (mean-baseline) subject -- the centered term then collapses to 0.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject and MEAN-centered (Zhou 2025 Methods",
        "'C-QTc Model': 'QTc0 is the overall mean baseline QTc').",
        "Centering reference: reported directly by the source -- 'The mean",
        "baseline QTcP and QTcF in the analysis dataset were 419.3 and",
        "409.5 ms, respectively' (Zhou 2025 Results 'Summary of Data Used",
        "for Analysis'). This model uses the QTcP value, 419.3 ms.",
        "Correction method: QTcP is the POPULATION-BASED correction,",
        "QTcP = QT / RR^beta with beta estimated from a log-log regression",
        "of QT on RR fitted to the cycle 1 day 1 pre-dose individual",
        "replicates of all patients (Zhou 2025 Methods 'QTcF and QTcP')."
      ),
      source_name        = "baseline QTcP (paper notation QTc0,i)"
    )
  )

  covariatesDataExcluded <- list(
    DAY21 = list(
      description = "Cycle 1 day 21 visit indicator (1 = cycle 1 day 21 steady-state ECG visit, 0 = cycle 1 day 1 first-dose visit).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested but NOT retained in this model. Zhou 2025 Supplementary",
        "Table 1 documents the two-step development log for the",
        "fruquintinib-only C-QTcP model: fit.f.1 included the visit term",
        "(AIC 10754.45) but 'the visit coefficient 95% CI included 0; the",
        "LRT p-value was 0.1318 when visit was dropped from the model';",
        "fit.f.2 dropped it (AIC 10754.72, dAIC +0.27) and is flagged",
        "'Best model with only fruquintinib concentration'. Supplementary",
        "Table 5, the parameter table for this model, accordingly has no",
        "visit row.",
        "The two M11-driven companion models",
        "(Zhou_2025_fruquintinib_QTcP_M11.R and",
        "Zhou_2025_fruquintinib_QTcF_M11.R) DO retain the visit term."
      )
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
      "Same 205-patient C-QTc analysis set as the primary QTcP analysis",
      "(137 fruquintinib, 68 placebo; 1456 time-matched pairs). Mean",
      "baseline QTcP 419.3 ms.",
      "This is the paper's NEGATIVE additional analysis: 'Additional C-QTc",
      "analysis including only fruquintinib concentrations showed no",
      "relationship between DeltaQTcP and fruquintinib concentrations with",
      "a slope (95% CI) of 0.00778 (-0.00247, 0.018) ms/ng/mL",
      "(P = .1377)' (Zhou 2025 Results 'Supportive Analyses').",
      "Model selection: this family had the highest AIC of the three",
      "(best 10754.45 / 10754.72 versus 10728.75 for fruquintinib + M11",
      "and 10727.07 for M11 alone) -- Zhou 2025 Supplementary Tables 1-3.",
      "Analysis software: R 4.2.0 with the nlme package (version 3.1-157)",
      "and the contrast package (version 0.24.2)."
    )
  )

  ini({
    # ==================================================================
    # Linear mixed-effects C-QTc model driven by PARENT fruquintinib
    # concentration (Zhou 2025 Supplementary Table 5, 'Parameter
    # Estimates of the C-QTc Model with Only Fruquintinib
    # Concentration'), corresponding to model fit.f.2 of Supplementary
    # Table 1 -- the visit term dropped, full BSV covariance matrix.
    #
    # Every parameter is on the LINEAR DeltaQTcP (ms) scale and every
    # random effect is ADDITIVE on that scale, matching the source nlme
    # fit. Neither the intercept nor the slope is log-transformed; see
    # Zhou_2025_fruquintinib_QTcP_M11.R for the reasoning.
    # ==================================================================

    e0 <- 2.33
    label("Linear-mixed-effects intercept on DeltaQTcP, placebo / 0 h nominal time / mean-baseline reference (ms)")
    # Zhou 2025 Supplementary Table 5 'Intercept' = 2.33 (SE 1.31;
    # RSE 56.2%; 95% CI -0.233, 4.89; P = 0.0756).

    slope <- 0.00778
    label("Linear fruquintinib-concentration slope on DeltaQTcP (ms per ng/mL; NOT statistically significant)")
    # Zhou 2025 Supplementary Table 5 'Fruquintinib Conc. Slope' =
    # 0.00778 (SE 0.00524; RSE 67.4%; 95% CI -0.00247, 0.018;
    # P = 0.1377). Also quoted verbatim in the main text (Results
    # 'Supportive Analyses'). The 95% CI spans zero: this is the paper's
    # negative finding, encoded as reported.

    # ------------------------------------------------------------------
    # Covariate effects (Zhou 2025 Supplementary Table 5). All estimated.
    # NOTE: there is NO cycle 1 day 21 visit term in this model -- the
    # authors dropped it (Supplementary Table 1, model fit.f.2). See
    # covariatesDataExcluded above.
    # ------------------------------------------------------------------

    e_on_treatment_e0 <- -5.24
    label("Effect of active-treatment arm on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 5 'Treatment' = -5.24 (SE 1.46;
    # RSE 27.9%; 95% CI -8.12, -2.37; P = 0.0004).

    e_ntime1_e0 <- -0.254
    label("Effect of the 1 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 5 'NTime=1' = -0.254 (SE 0.851;
    # RSE 335%; 95% CI -1.92, 1.41; P = 0.7655).

    e_ntime2_e0 <- 1.6
    label("Effect of the 2 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 5 'NTime=2' = 1.6 (SE 0.846;
    # RSE 52.9%; 95% CI -0.0524, 3.26; P = 0.0584).

    e_ntime3_e0 <- 1.6
    label("Effect of the 3 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 5 'NTime=3' = 1.6 (SE 0.846;
    # RSE 52.9%; 95% CI -0.0611, 3.25; P = 0.0598). Printed with the same
    # point estimate and SE as NTime=2 but a distinct 95% CI and P-value,
    # so the two coefficients are separately estimated values that agree
    # to the two significant figures the table prints.

    e_ntime4_e0 <- 0.353
    label("Effect of the 4 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Supplementary Table 5 'NTime=4' = 0.353 (SE 0.854;
    # RSE 242%; 95% CI -1.32, 2.02; P = 0.6797).

    e_qtc_bl_e0 <- -0.105
    label("Effect of mean-centered baseline QTcP on the intercept (additive, ms per ms)")
    # Zhou 2025 Supplementary Table 5 'Baseline QTcP' = -0.105
    # (SE 0.0309; RSE 29.4%; 95% CI -0.166, -0.0442; P = 0.0008).

    # ==================================================================
    # Between-subject variability -- FULL (correlated) covariance matrix,
    # unlike the diagonal structure of the two M11-driven models. Zhou
    # 2025 Supplementary Table 1 describes both fruquintinib-only
    # candidates as carrying a "full BSV variance-covariance matrix", and
    # Supplementary Table 5 reports the correlation explicitly.
    #
    # Supplementary Table 5 reports SDs and a CORRELATION; ini() takes
    # variances and a COVARIANCE, so:
    #   var(etae0)     = 8.62^2   = 74.3044
    #   var(etaslope)  = 0.0400^2 = 0.0016
    #   cov            = -0.340 * 8.62 * 0.0400 = -0.117232
    # The resulting 2x2 matrix is positive definite (determinant
    # 74.3044 * 0.0016 - 0.117232^2 = 0.1051), so no nudge is needed.
    # ==================================================================

    etae0 + etaslope ~ c(
      74.3044,
      -0.117232, 0.0016
    )
    # Zhou 2025 Supplementary Table 5:
    #   'BSV SD for Intercept'            = 8.62 msec (SE 0.591;
    #                                        RSE 6.86%; 95% CI 7.54, 9.85)
    #   'BSV SD of Conc. Slope'           = 0.0400 msec per ng/mL
    #                                        (SE 0.00485; RSE 12.1%;
    #                                        95% CI 0.0316, 0.0505)
    #   'Slope:Intercept BSV Correlation' = -0.340 (SE 0.141; RSE 41.5%;
    #                                        95% CI -0.551, -0.0872)

    # ==================================================================
    # Residual error: additive on the DeltaQTcP scale.
    # ==================================================================
    addSd <- 8.02
    label("Additive residual error standard deviation on DeltaQTcP (ms)")
    # Zhou 2025 Supplementary Table 5 'Residual Error SD' = 8.02 msec
    # (SE 0.167; RSE 2.08%; 95% CI 7.70, 8.35).
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
    # 2. Mean-centered baseline QTcP (reported reference 419.3 ms).
    # ==================================================================
    qtc_bl_ref      <- 419.3
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    # ==================================================================
    # 3. Per-subject intercept and slope. NO visit term (dropped by the
    #    authors; see covariatesDataExcluded). Both random effects are
    #    additive on the linear DeltaQTcP scale and correlated
    #    (rho = -0.340).
    # ==================================================================
    e0_i <-
      e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      ntime_effect +
      e_qtc_bl_e0 * qtc_bl_centered
    slope_i <- slope + etaslope

    # ==================================================================
    # 4. Prediction: change from baseline in the population-based
    #    corrected QT interval (DeltaQTcP), in ms.
    # ==================================================================
    QTcP <- e0_i + slope_i * CP_FRUQUINTINIB_NGML

    QTcP ~ add(addSd)
  })
}
