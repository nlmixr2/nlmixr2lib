Zhou_2025_fruquintinib_QTcP_M11 <- function() {
  description <- paste(
    "Concentration-QTc (C-QTc) linear mixed-effects PD model relating the",
    "change from baseline in the population-based corrected QT interval",
    "(DeltaQTcP, ms) to the plasma concentration of M11, the major",
    "metabolite of fruquintinib, in 205 patients with previously treated",
    "metastatic colorectal cancer (137 fruquintinib 5 mg once daily, 68",
    "matching placebo) enrolled in the ECG substudy of the phase 3",
    "FRESCO-2 trial (NCT04322539). This is the FINAL model of Zhou 2025:",
    "it had the lowest AIC (10727.07) of the three model families the",
    "authors fitted (fruquintinib only, M11 only, fruquintinib + M11).",
    "The model is:",
    "DeltaQTcP = e0 + etae0 + e_on_treatment_e0 * ON_TREATMENT",
    "+ (slope + etaslope) * CP_FRUQUINTINIB_M11_NGML",
    "+ sum_j(e_ntime<j>_e0 * [NTIME == j]) + e_day21_e0 * DAY21",
    "+ e_qtc_bl_e0 * (QTC_BL - 419.3),",
    "with independent (diagonal) between-subject variability on the",
    "intercept and on the M11 concentration slope, and additive residual",
    "error. The four nominal-time coefficients (1, 2, 3 and 4 h after",
    "dose, with the 0 h nominal time as the reference) carry the",
    "drug-independent diurnal variation of DeltaQTcP. The treatment term",
    "e_on_treatment_e0 is an empirical flexibility term with no",
    "physiological interpretation (Zhou 2025 Methods, citing Garnett",
    "2018); together with the concentration slope it defines the",
    "placebo-corrected contrast, mean DeltaDeltaQTcP =",
    "e_on_treatment_e0 + slope * CP_FRUQUINTINIB_M11_NGML, because every",
    "other term is common to the active and placebo arms and cancels.",
    "PD-only model: M11 plasma concentration is supplied as a",
    "time-varying covariate (ng/mL); placebo records carry 0. The source",
    "publication does not fit a population PK model, so users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "M11 concentration trajectory (no fruquintinib or M11 popPK model",
    "exists in the nlmixr2lib registry). Companion model files:",
    "Zhou_2025_fruquintinib_QTcF_M11.R (the supportive analysis on the",
    "Fridericia-corrected endpoint) and",
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
    sep = " "
  )

  vignette <- "Zhou_2025_fruquintinib"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external M11 plasma-concentration covariate)",
    concentration = "(observation QTcP is the CHANGE FROM BASELINE in the population-based corrected QT interval, DeltaQTcP, in ms; driving covariate CP_FRUQUINTINIB_M11_NGML is in ng/mL)"
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
        "and 1 for active drug)' (Zhou 2025 Methods 'C-QTc Model').",
        "The paper is explicit that the coefficient is empirical: 'the",
        "theta_1 term has no physiological interpretation but allows",
        "flexibility in the event of model specification' and 'is an",
        "empirical term allowing the linear mixed effects model to be",
        "flexible, for example, in cases where DeltaQTc is nonlinear'",
        "(Zhou 2025 Methods and Results 'C-QTcP Model and Evaluation').",
        "Because it is an arm-level offset rather than an exposure effect,",
        "it does NOT vanish at zero concentration: the placebo-corrected",
        "contrast at concentration C is e_on_treatment_e0 + slope * C.",
        "Trial design: FRESCO-2 randomized patients 2:1 to fruquintinib",
        "5 mg QD on days 1-21 of a 28-day cycle plus best supportive care,",
        "or matching placebo plus best supportive care. In the C-QTc",
        "analysis set 137 patients were in the fruquintinib arm and 68 in",
        "the placebo arm (Zhou 2025 Results 'Summary of Data Used for",
        "Analysis'). All placebo records were assigned a fruquintinib and",
        "M11 concentration of 0 at each nominal time point (Zhou 2025",
        "Methods 'Overview of Data')."
      ),
      source_name        = "Trt (paper notation TRT_i)"
    ),
    CP_FRUQUINTINIB_M11_NGML = list(
      description        = "Instantaneous total plasma concentration of M11, the major metabolite of fruquintinib, at the time of each time-matched ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear",
        "concentration-DeltaQTcP expression DeltaQTcP = ... +",
        "(slope + etaslope) * CP_FRUQUINTINIB_M11_NGML + ... , with the",
        "slope reported directly in ms per ng/mL (Zhou 2025 Table 2 units",
        "column), so no in-model unit rescaling is required.",
        "In Zhou 2025 this was the observed M11 plasma concentration drawn",
        "at the same NOMINAL time point as the triplicate Holter-extracted",
        "12-lead ECG; pairs whose ECG and PK sampling were performed more",
        "than 30 minutes apart were excluded (Zhou 2025 Methods 'Overview",
        "of Data' and Results 'Summary of Data Used for Analysis').",
        "Placebo records were assigned 0 at every nominal time point.",
        "Of the 936 fruquintinib-arm concentration pairs, 326 (34.8%) M11",
        "concentrations were below the limit of quantification and were",
        "set to zero and retained in the analysis; the high BLQ fraction",
        "is why the authors also report the parent-driven model (see",
        "Zhou_2025_fruquintinib_QTcP_parent.R).",
        "Reference values observed: the geometric mean steady-state Cmax",
        "of M11 after fruquintinib 5 mg QD is 77 ng/mL, and twice that",
        "value is 154 ng/mL; these are the two scenarios tabulated in",
        "Zhou 2025 Table 3. The model predicts the upper bound of the 90%",
        "CI of mean DeltaDeltaQTcP to exceed 10 ms at 262 ng/mL, 3.4-fold",
        "the observed steady-state geometric mean Cmax (Zhou 2025 Results",
        "'Final Model Predictions').",
        "Set to 0 for placebo records and outside the drug-exposure",
        "window; the concentration-slope term then collapses to 0."
      ),
      source_name        = "CONC2 (supplement model-development logs); 'M11 conc.' (Table 2)"
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
        "labelled 'NTime = 1' through 'NTime = 4' with units 'ms' in",
        "Table 2, and the Table 2 footnote defines 'NTime, nominal time in",
        "hours'. The supplement's model-development logs write the term as",
        "NTIME and gloss it as 'nominal time point after the most recent",
        "dose' (Zhou 2025 Supplementary Tables 1-3 footnotes).",
        "The paper interprets these coefficients as drug-independent:",
        "'The nominal time coefficients reflect the drug-independent",
        "diurnal variation of DeltaQTc over time' (Zhou 2025 Results",
        "'C-QTcP Model and Evaluation'). They are therefore common to",
        "both treatment arms and cancel out of the placebo-corrected",
        "DeltaDeltaQTc contrast.",
        "NOMINAL, not actual, time: a sample drawn at 0.4 h actual time",
        "after dose belongs to the 1 h nominal time point if that is the",
        "protocol slot it was collected in. Populate this column from the",
        "protocol schedule, not from rxode2's tad(). Because the model()",
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
        "2025 Methods 'C-QTc Model'); tabulated as the 'Cycle 1 day 21",
        "visit' row of Table 2.",
        "Design: fruquintinib was given on days 1-21 of a 28-day cycle, so",
        "the cycle 1 day 1 visit captures the first dose and the cycle 1",
        "day 21 visit captures steady state after 21 consecutive daily",
        "doses (Zhou 2025 Methods 'Overview of Data' and Discussion).",
        "The visit coefficient is drug-independent in the same sense as",
        "the nominal-time coefficients -- it is estimated across both arms",
        "and cancels out of the placebo-corrected DeltaDeltaQTc contrast.",
        "The companion parent-concentration model",
        "(Zhou_2025_fruquintinib_QTcP_parent.R) DROPS this term: the",
        "authors removed it from the fruquintinib-only model because its",
        "95% CI included 0 and the likelihood-ratio-test p-value was",
        "0.1318 (Zhou 2025 Supplementary Table 1, model fit.f.2)."
      ),
      source_name        = "visit (paper notation Visit_k, k = cycle 1 day 1 or 21)"
    ),
    QTC_BL = list(
      description        = "Subject's pre-dose (cycle 1 day 1) baseline population-based corrected QT interval (QTcP), treated as a per-subject time-fixed covariate. Enters the linear-mixed-effects intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 419.3). Set QTC_BL = 419.3 ms for the typical (mean-baseline) subject -- the centered term then collapses to 0.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Source paper: 'theta_5 is the fixed",
        "effect associated with baseline QTc, QTc0,i is the baseline QTc",
        "of patient i, QTc0 is the overall mean baseline QTc' (Zhou 2025",
        "Methods 'C-QTc Model'), i.e. the covariate is MEAN-centered, not",
        "median-centered as in Darpo 2014.",
        "Centering reference: unlike Darpo 2014 (where the centering value",
        "was not quoted and a rounded standard had to be assumed), Zhou",
        "2025 reports the value directly: 'The mean baseline QTcP and QTcF",
        "in the analysis dataset were 419.3 and 409.5 ms, respectively'",
        "(Zhou 2025 Results 'Summary of Data Used for Analysis'). This",
        "file therefore uses qtc_bl_ref = 419.3 ms -- a reported value,",
        "not an assumption.",
        "Correction method: QTcP is the POPULATION-BASED correction,",
        "derived by fitting log(QT) = alpha + beta * log(RR) with no",
        "random effects to the cycle 1 day 1 pre-dose individual",
        "replicates of all patients and applying QTcP = QT / RR^beta",
        "(RR in ms). The authors selected QTcP over QTcF because a",
        "graphical assessment showed Fridericia's formula did not",
        "adequately correct for heart rate in this cohort: the baseline",
        "QTcP-RR slope was -3.73e-05 with a 90% CI (-0.0102, 0.0101)",
        "spanning 0, whereas the baseline QTcF-RR slope was 0.0493",
        "(90% CI 0.0393, 0.0592) (Zhou 2025 Methods 'QTcF and QTcP',",
        "Supplementary Figure 1, and Discussion).",
        "For the Fridericia-corrected companion model the baseline is the",
        "QTcF baseline with a 409.5 ms centering reference; see",
        "Zhou_2025_fruquintinib_QTcF_M11.R."
      ),
      source_name        = "baseline QTcP (paper notation QTc0,i)"
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
      "in the cardiovascular-safety ECG subset. Patients had progressed on",
      "or were intolerant to standard therapies including",
      "fluoropyrimidine-, oxaliplatin- and irinotecan-based chemotherapy,",
      "anti-VEGF agents and anti-EGFR agents. The study excluded patients",
      "with a QTcF > 480 ms, patients with any other factor that could",
      "prolong the QTc interval or increase arrhythmic risk, and patients",
      "receiving concomitant QTc-prolonging medication; concomitant",
      "medications known to be associated with QTc prolongation were",
      "prohibited (Zhou 2025 Methods 'Overview of Data')."
    ),
    dose_range       = paste(
      "Fruquintinib 5 mg orally once daily on days 1-21 of a 28-day cycle",
      "plus best supportive care, or matching placebo plus best supportive",
      "care. 5 mg QD is the only dose in this analysis; the maximum",
      "tolerated dose is 6 mg QD (one dose-limiting toxicity in six",
      "patients with advanced solid tumors), so a supratherapeutic dose",
      "could not be given and no positive control (e.g. moxifloxacin) was",
      "included, consistent with ICH E14 Q&A 6.1 approaches for anticancer",
      "agents (Zhou 2025 Introduction and Discussion).",
      "Time-matched ECG / PK pairs were collected at nominal times 0, 1, 2,",
      "3 and 4 h after dose at the cycle 1 day 1 and cycle 1 day 21 visits",
      "(Zhou 2025 Table 5 and Supplementary Tables 6-7; the cycle 1 day 1",
      "visit contributes only the 1-4 h nominal times because the 0 h",
      "reading is the baseline)."
    ),
    regions          = NA_character_,
    notes            = paste(
      "C-QTc analysis set: 1456 time-matched concentration-DeltaQTcP",
      "pairs from 205 patients (137 fruquintinib, 68 placebo), drawn from",
      "1954 sets of Holter ECG measurements in 243 patients (163",
      "fruquintinib, 80 placebo). Exclusions were for missing baseline",
      "QTc, ECG measurements without time-matched concentrations, and ECG",
      "and PK sampling performed more than 30 minutes apart (Zhou 2025",
      "Results 'Summary of Data Used for Analysis').",
      "Of the 936 fruquintinib-arm concentration pairs, 16 (1.7%)",
      "fruquintinib and 326 (34.8%) M11 concentrations were below the",
      "limit of quantification; BLQ values were set to zero and included.",
      "Mean baseline QTcP 419.3 ms; mean baseline QTcF 409.5 ms.",
      "ECG: 12-lead ECGs extracted in triplicate from continuous Holter",
      "recordings and processed centrally by eResearch Technology, Inc.",
      "(Zhou 2025 Methods 'QTcF and QTcP').",
      "Analysis software: R 4.2.0 with the nlme package (version 3.1-157)",
      "for the mixed-effects fits and the contrast package (version",
      "0.24.2) for the mean DeltaDeltaQTc predictions and their 90% CIs",
      "(Zhou 2025 Methods 'Cardiac Safety Analyses' and 'Categorical, by-",
      "Time Point, and Central Tendency Analyses'). Model selection was by",
      "AIC across the three model families: fruquintinib only (best AIC",
      "10754.45), fruquintinib + M11 (best AIC 10728.75), and M11 only",
      "(best AIC 10727.07, the final model encoded here) -- Zhou 2025",
      "Supplementary Tables 1-3."
    )
  )

  ini({
    # ==================================================================
    # Linear mixed-effects C-QTc model (Zhou 2025 Methods 'C-QTc Model'
    # and Table 2, 'Final C-QTc Model Parameter Estimates - QTcP Model
    # with M11 Concentration'). The published equation is:
    #   DeltaQTc_ijk = theta0 + eta0_i + theta1 * TRT_i
    #                  + (theta2 + eta2_i) * C_ijk
    #                  + theta3(NTime_j) + theta4 * Visit_k
    #                  + theta5 * (QTc0_i - mean(QTc0)) + eps_ijk
    # for patient i at nominal time j at visit k. Random effects are
    # normal with mean 0; the final model uses a DIAGONAL between-subject
    # covariance matrix (Zhou 2025 Supplementary Table 3, model fit.m.1a
    # 'Intercept with BSV + CONC2 with BSV + trt + visit + baseline QTcP
    # + NTIME; diagonal BSV variance-covariance matrix', flagged
    # 'Final model', AIC 10727.07).
    #
    # NOTE ON SCALE. Every parameter below is on the LINEAR DeltaQTcP
    # (ms) scale, and every random effect is ADDITIVE on that scale,
    # exactly as the source nlme fit specifies. Neither the intercept nor
    # the slope is log-transformed: the intercept of a change-from-
    # baseline endpoint is sign-free (cf. Darpo_2014_racSotalol_QTcF.R),
    # and log-transforming the slope would misrepresent its random
    # effect, whose SD (0.0953) is nearly three times the typical value
    # (0.0339) so that a substantial fraction of subjects have a
    # NEGATIVE individual slope under the published normal BSV.
    # ==================================================================

    e0 <- 2.82
    label("Linear-mixed-effects intercept on DeltaQTcP, placebo / cycle 1 day 1 / 0 h nominal time / mean-baseline reference (ms)")
    # Zhou 2025 Table 2 'Intercept' = 2.82 (SE 1.35; RSE 47.9%;
    # 95% CI 0.176, 5.47; P = 0.0372).

    slope <- 0.0339
    label("Linear M11-concentration slope on DeltaQTcP (ms per ng/mL)")
    # Zhou 2025 Table 2 'M11 conc. slope' = 0.0339 (SE 0.0147; RSE 43.4%;
    # 95% CI 0.00516, 0.0625; P = 0.0212). The paper highlights this as
    # "a statistically significant slope of 0.0339 (95% CI 0.00516,
    # 0.0625) ms per ng/mL (P = .0212)" (Results 'C-QTcP Model and
    # Evaluation'). Reported directly in ms per ng/mL, so the covariate
    # CP_FRUQUINTINIB_M11_NGML is consumed in ng/mL with no rescaling.

    # ------------------------------------------------------------------
    # Covariate effects (Zhou 2025 Table 2). All are estimated point
    # estimates with reported SE / RSE / 95% CI / P-value; none is fixed.
    # ------------------------------------------------------------------

    e_on_treatment_e0 <- -5.19
    label("Effect of active-treatment arm on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'Treatment' = -5.19 (SE 1.36; RSE 26.2%;
    # 95% CI -7.86, -2.52; P = 0.0002). Empirical flexibility term with
    # no physiological interpretation (Zhou 2025 Methods 'C-QTc Model'
    # and Results). Together with `slope` it defines the placebo-
    # corrected contrast: mean DeltaDeltaQTcP = -5.19 + 0.0339 * C.

    e_ntime1_e0 <- -0.25
    label("Effect of the 1 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'NTime = 1' = -0.25 (SE 0.867; RSE 347%;
    # 95% CI -1.95, 1.45; P = 0.7733). Not statistically significant;
    # retained per the Garnett 2018 C-QTc white-paper recommendation
    # ("The nonsignificant coefficients and the visit term were kept in
    # the model, which aligns with recommendations from the scientific
    # white paper on concentration-QTc modeling", Zhou 2025 Results).

    e_ntime2_e0 <- 1.78
    label("Effect of the 2 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'NTime = 2' = 1.78 (SE 0.867; RSE 48.7%;
    # 95% CI 0.0873, 3.48; P = 0.0399).

    e_ntime3_e0 <- 1.84
    label("Effect of the 3 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'NTime = 3' = 1.84 (SE 0.868; RSE 47.2%;
    # 95% CI 0.144, 3.54; P = 0.0340).

    e_ntime4_e0 <- 0.570
    label("Effect of the 4 h nominal time point on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'NTime = 4' = 0.570 (SE 0.874; RSE 153%;
    # 95% CI -1.14, 2.28; P = 0.5150). Not statistically significant;
    # retained per the same white-paper recommendation as NTime = 1.

    e_day21_e0 <- -1.39
    label("Effect of the cycle 1 day 21 visit on the intercept (additive, ms)")
    # Zhou 2025 Table 2 'Cycle 1 day 21 visit' = -1.39 (SE 0.676;
    # RSE 48.6%; 95% CI -2.71, -0.0638; P = 0.0406).

    e_qtc_bl_e0 <- -0.102
    label("Effect of mean-centered baseline QTcP on the intercept (additive, ms per ms)")
    # Zhou 2025 Table 2 'Baseline QTcP' = -0.102 (SE 0.0306; RSE 30.0%;
    # 95% CI -0.162, -0.0415; P = 0.0011). The negative coefficient is
    # the usual regression-to-the-mean shrinkage of a change-from-
    # baseline endpoint at high baseline.

    # ==================================================================
    # Between-subject variability. Zhou 2025 Table 2 reports the BSV as
    # STANDARD DEVIATIONS on the linear DeltaQTcP scale; nlmixr2's ini()
    # takes VARIANCES, so each is squared here. The final model uses a
    # DIAGONAL covariance matrix -- the intercept and slope random
    # effects are independent (Zhou 2025 Supplementary Table 3, model
    # fit.m.1a; the full-covariance model fit.m.1 estimated a
    # correlation of only 0.018 and had a higher AIC by 1.99). The two
    # etas are therefore declared separately rather than as a block.
    #
    # WARNING ABOUT THE MAIN-TEXT PROSE. Zhou 2025 Results 'C-QTcP Model
    # and Evaluation' states "The standard deviations of the between-
    # subject variability (BSV) on the intercept and the M11
    # concentration slope were 8.02 ms and 0.116 ms per ng/mL." Those
    # two numbers are NOT the QTcP model's values -- they are exactly
    # the QTcF model's values from Supplementary Table 4 (8.02 and
    # 0.116), i.e. a copy-paste error in the prose. Table 2, the table
    # the sentence points at, reports 8.23 and 0.0953 with internally
    # consistent SEs, RSEs and 95% CIs (7.33-9.24 and 0.0751-0.121).
    # This file uses the Table 2 values; see the vignette Errata.
    # ==================================================================

    etae0 ~ 67.7329
    # Zhou 2025 Table 2 'BSV SD for intercept' = 8.23 ms (SE 0.486;
    # RSE 5.91%; 95% CI 7.33, 9.24). Variance = 8.23^2 = 67.7329.

    etaslope ~ 0.00908209
    # Zhou 2025 Table 2 'BSV SD of the M11 conc. slope' = 0.0953 ms per
    # ng/mL (SE 0.0118; RSE 12.4%; 95% CI 0.0751, 0.121).
    # Variance = 0.0953^2 = 0.00908209.

    # ==================================================================
    # Residual error. Zhou 2025 Methods 'C-QTc Model': "Residuals were
    # assumed to be normally distributed with a mean of 0 and a variance
    # R", i.e. additive on the DeltaQTcP scale. Table 2 reports the SD.
    # ==================================================================
    addSd <- 7.93
    label("Additive residual error standard deviation on DeltaQTcP (ms)")
    # Zhou 2025 Table 2 'Residual error SD' = 7.93 ms (SE 0.165;
    # RSE 2.08%; 95% CI 7.61, 8.26).
  })

  model({
    # ==================================================================
    # 1. Nominal-time decomposition. NTIME is the canonical integer
    #    nominal-time column (protocol-scheduled hours after the most
    #    recent dose); the reference level NTIME = 0 sets every
    #    indicator to zero. Zhou 2025 Methods 'C-QTc Model': "theta_3 is
    #    the fixed effect associated with nominal time (one value
    #    estimated nonreference per nominal time point)".
    # ==================================================================
    nt1 <- (NTIME == 1)
    nt2 <- (NTIME == 2)
    nt3 <- (NTIME == 3)
    nt4 <- (NTIME == 4)
    ntime_effect <-
      e_ntime1_e0 * nt1 + e_ntime2_e0 * nt2 +
      e_ntime3_e0 * nt3 + e_ntime4_e0 * nt4

    # ==================================================================
    # 2. Mean-centered baseline QTcP. Unlike Darpo 2014, the centering
    #    reference here is REPORTED, not assumed: "The mean baseline
    #    QTcP and QTcF in the analysis dataset were 419.3 and 409.5 ms,
    #    respectively" (Zhou 2025 Results 'Summary of Data Used for
    #    Analysis').
    # ==================================================================
    qtc_bl_ref      <- 419.3
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    # ==================================================================
    # 3. Per-subject intercept and slope. Both random effects are
    #    ADDITIVE on the linear DeltaQTcP scale, as the source nlme fit
    #    specifies.
    # ==================================================================
    e0_i <-
      e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      ntime_effect +
      e_day21_e0 * DAY21 +
      e_qtc_bl_e0 * qtc_bl_centered
    slope_i <- slope + etaslope

    # ==================================================================
    # 4. Prediction. The observation variable is the CHANGE FROM
    #    BASELINE in the population-based corrected QT interval
    #    (DeltaQTcP), in ms, named QTcP per the canonical
    #    heart-rate-corrected QT endpoint register.
    # ==================================================================
    QTcP <- e0_i + slope_i * CP_FRUQUINTINIB_M11_NGML

    QTcP ~ add(addSd)
  })
}
