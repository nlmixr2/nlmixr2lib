Elhefnawy_2023_recurrent_ischemic_stroke_nondm <- function() {
  description <- paste(
    "Parametric time-to-event model for recurrent ischemic stroke (IS) after a",
    "first (index) IS, in the NON-DIABETIC stratum of the National Neurology",
    "Registry of Malaysia (4,204 patients without diabetes mellitus, 138",
    "recurrences, up to 7.37 years of follow-up), developed in NONMEM 7.5.",
    "The hazard is Gompertz, h(t) = theta_x * exp(theta_y * t), with the scale",
    "switching at 0.5 year (theta1 -> theta3) and the shape switching",
    "independently at 3 years (theta2 -> theta4). Three covariates act",
    "log-linearly on the hazard: hyperlipidemia and hypertension recorded",
    "before the index stroke raise it, and an antihyperlipidemic prescribed for",
    "secondary prevention lowers it. Time is in years and there is no drug",
    "exposure term, so this is a disease-progression / event-risk model rather",
    "than a PK-PD model. The cumulative hazard is evaluated in closed form (the",
    "piecewise Gompertz integral is elementary) rather than by an ODE, which",
    "avoids integrating across the two hazard discontinuities; the model",
    "exposes the instantaneous hazard `hazard_ris`, the cumulative hazard",
    "`cumhaz_ris` and the recurrence-free survivor probability `sur`. Sister",
    "model file for the diabetic stratum of the same paper:",
    "modellib('Elhefnawy_2023_recurrent_ischemic_stroke_dm').",
    "IMPORTANT -- the published baseline-hazard scale does not reproduce the",
    "paper's own event count. As printed, the model predicts a 32.6 percent",
    "chance of recurrence within 7.37 years for a covariate-free non-diabetic",
    "patient, against the 138 of 4,204 = 3.28 percent actually reported in",
    "Elhefnawy 2023 Table 1: an 11.8-fold overprediction on the",
    "cumulative-hazard scale, rising to 13.8-fold once the cohort covariate",
    "distribution is applied. No reading of the piecewise structure removes it",
    "-- theta1 = 0.253 per year held flat over just the first 6 months already",
    "contributes 0.127, 3.8 times the entire 7.37-year observed cumulative",
    "hazard of 0.0334. The covariate block, by contrast, reproduces the paper's",
    "own adjusted hazard ratios and half-life column exactly. The values are",
    "transcribed as published and have NOT been tuned; see the validation",
    "vignette's Errata section for the full arithmetic. The same defect, at the",
    "same magnitude, is present in this group's companion pooled-cohort",
    "publication (Elhefnawy 2023, Front Neurol 14:1118711, PMC10176964).",
    sep = " "
  )
  reference <- paste(
    "Elhefnawy M, Noor Harun S, Leykhim T, Tangiisuran B, Zainal H, Looi I,",
    "Sidek N, Abdul Aziz Z, Sheikh Ghadzi SM.",
    "A parametric time-to-event modelling of recurrent ischemic stroke after",
    "index stroke among patients with and without diabetes mellitus:",
    "implementation of temporal validation of the model.",
    "Cureus. 2023;15(12):e50794. doi:10.7759/cureus.50794.",
    "Open Access under CC BY 4.0.",
    "The survival / hazard relationship is Equation 1 and the Gompertz hazard",
    "is Table 2 row 4; the piecewise scale / shape switching is defined in the",
    "Table 2 footnote and restated in Methods, 'Base Model Development'; the",
    "covariate form is the Methods equation h = h0 * exp(beta_t*t +",
    "beta_cov*(cov)) (Equation 5, printed inline on page 2); all parameter",
    "estimates are Table 3, non-DM columns. The EuropePMC open-access deposit",
    "for PMC10796130 contains no supplementary files.",
    "Sister model file from the same paper:",
    "modellib('Elhefnawy_2023_recurrent_ischemic_stroke_dm').",
    sep = " "
  )
  vignette <- "Elhefnawy_2023_recurrent_ischemic_stroke_diabetes"

  units <- list(
    time          = "year",
    dosing        = "n/a (no dosing events; secondary-prevention therapy enters only through the binary CONMED_LIPIDLOWER covariate)",
    concentration = "n/a (the model outputs are a hazard in 1/year, a unitless cumulative hazard and a unitless recurrence-free survivor probability, not a drug concentration)"
  )

  covariateData <- list(
    DIS_HYPERLIP = list(
      description        = "1 = the patient carried a hyperlipidemia (HPLD) diagnosis before the index ischemic stroke; 0 = no hyperlipidemia. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hyperlipidemia before the index stroke)",
      notes              = "Ascertained from the National Neurology Registry concurrent-disease fields (Elhefnawy 2023 Methods, 'Data collection'). Prevalence in the 4,204-patient non-diabetic stratum is (63 + 865) / 4,204 = 22.07 percent (Table 1). The strongest predictor in this stratum: aHR = exp(1.03) = 2.80 (95 percent CI 2.00-3.90). The same covariate is retained, more weakly, in the diabetic stratum (aHR 1.88) -- see modellib('Elhefnawy_2023_recurrent_ischemic_stroke_dm').",
      source_name        = "HPLD"
    ),
    DIS_HYPERT = list(
      description        = "1 = the patient carried a hypertension (HTN) diagnosis before the index ischemic stroke; 0 = no hypertension. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hypertension before the index stroke)",
      notes              = "Ascertained the same way as the other concurrent-disease flags. Prevalence in the non-diabetic stratum is (108 + 2,355) / 4,204 = 58.59 percent (Table 1). aHR = exp(0.789) = 2.20 (95 percent CI 1.53-3.14). Retained in this stratum only: Elhefnawy 2023 Discussion reads the absence of an HTN term in the diabetic stratum as collinearity with diabetes, and hypertension prevalence there is 87.1 percent, leaving little contrast. Table 1 additionally stratifies hypertension duration at 5 years; no duration effect was retained.",
      source_name        = "HTN"
    ),
    CONMED_LIPIDLOWER = list(
      description        = "1 = the patient was prescribed an antihyperlipidemic (lipid-lowering) medication for secondary prevention or concurrent-disease control around the index ischemic stroke admission; 0 = no antihyperlipidemic prescribed. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no antihyperlipidemic prescribed)",
      notes              = "Elhefnawy 2023 Methods, 'Data collection': the registry recorded 'medications for secondary prevention or concurrent disease control', and Table 1 tabulates the class as 'Antihyperlipidemic' under 'Received medication for concurrent disease control and/or secondary prevention'. The paper does not enumerate which agents were pooled into the class, so the class composition is unknown here; the Discussion argues the effect specifically in terms of statins, but the covariate as fitted is the broader antihyperlipidemic indicator. Prevalence in the non-diabetic stratum is (120 + 3,682) / 4,204 = 90.44 percent (Table 1). This is the only protective term retained in either stratum: aHR = exp(-0.835) = 0.434 (the paper prints the interval descending, '0.433(0.65-0.28)'), i.e. Results, 'receiving antihyperlipidemic medication among IS patients with non-DM decreased recurrent IS by 43.3%'. Because the indicator is fixed at the index admission it encodes prescription rather than adherence, and the model carries no time-varying exposure term.",
      source_name        = "Antihyperlipidemic"
    )
  )

  # Covariates that Elhefnawy 2023 collected and screened for the non-diabetic
  # stratum but did not retain in the final model. Listed here for provenance
  # only; none is referenced in model(). The paper reports no point estimate
  # for any of them -- the covariate screen itself is reported only by
  # reference to a previously published preprint (Methods, 'Covariate Model
  # Development and Model Evaluation'), and the EuropePMC deposit for
  # PMC10796130 carries no supplementary file, so no dOFV is available here
  # either.
  #
  # IHD is the one worth calling out: it is the strongest predictor in the
  # diabetic stratum (aHR 2.40) and absent here, which is the asymmetry the
  # paper's Discussion attributes to diabetes-associated angiopathy making the
  # shared atherosclerotic pathophysiology more prominent.
  #
  # These names are documentation labels for the paper's screen, not registered
  # canonical covariate columns: none is used in model(), so none is added to
  # inst/references/covariate-columns.md.
  covariatesDataExcluded <- list(
    DIS_IHD = list(
      description = "Ischemic heart disease (IHD) before the index stroke; (25 + 382) / 4,204 = 9.68 percent of the non-diabetic stratum.",
      units = "(binary)", type = "binary",
      notes = "Screened but not retained in the non-diabetic stratum, while it IS the strongest predictor in the diabetic stratum (aHR 2.40). Prevalence here is also lower (9.68 vs 13.51 percent). Registered canonical DIS_IHD is used for it in the sister DM model file."
    ),
    DIS_HYPERURICEMIA = list(
      description = "Hyperuricemia before the index stroke; (6 + 97) / 4,204 = 2.45 percent.",
      units = "(binary)", type = "binary",
      notes = "Named among the investigated concurrent diseases in Methods, 'Data collection' and tabulated in Table 1; no effect retained in Table 3."
    ),
    DIS_AF = list(
      description = "Atrial fibrillation before the index stroke; (5 + 172) / 4,204 = 4.21 percent.",
      units = "(binary)", type = "binary",
      notes = "Named among the investigated concurrent diseases and tabulated in Table 1; no effect retained."
    ),
    AGE = list(
      description = "Age at the index ischemic stroke; median 62.9 years across the whole study, dichotomised at 60 years in Table 1.",
      units = "year", type = "continuous",
      notes = "Screened as a demographic covariate (Methods, 'Data collection'); not retained. Table 1: 57.97 percent of recurrent and 60.45 percent of non-recurrent non-diabetic patients were older than 60 years."
    ),
    SEXF = list(
      description = "Female sex; (53 + 1,607) / 4,204 = 39.49 percent of the non-diabetic stratum.",
      units = "(binary)", type = "binary",
      notes = "Screened as a demographic covariate; not retained. Table 1 shows near-identical proportions in the recurrent (38.40 percent) and non-recurrent (39.52 percent) non-diabetic groups."
    ),
    SMOKER = list(
      description = "Current smoker at the index stroke; (89 + 1,917) / 4,204 = 47.72 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1 with a sizeable unadjusted imbalance (64.49 percent of recurrent vs 47.14 percent of non-recurrent non-diabetic patients); no effect retained in Table 3."
    ),
    FAMHX_STROKE = list(
      description = "Family history of stroke; (11 + 247) / 4,204 = 6.14 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    NIHSS = list(
      description = "National Institutes of Health Stroke Scale severity of the index stroke, dichotomised by the paper into minor vs moderate/severe.",
      units = "(score)", type = "continuous",
      notes = "Tabulated in Table 1 (46.37 percent minor among recurrent non-diabetic patients vs 46.06 percent among non-recurrent); no NIHSS term appears in Table 3."
    ),
    RACE = list(
      description = "Ethnicity, recorded by the registry as Malay / Chinese / Indian / Others.",
      units = "(categorical)", type = "categorical",
      notes = "Tabulated in Table 1 with a large unadjusted imbalance (54.34 percent Malay among recurrent non-diabetic patients vs 19.23 percent among non-recurrent); no ethnicity effect appears in Table 3."
    ),
    CONMED_ANTIPLATELET = list(
      description = "Antiplatelet (APLT) prescribed for secondary prevention; (118 + 3,635) / 4,204 = 89.27 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained in Table 3 for either stratum. Contrast the companion pooled-cohort publication (Elhefnawy 2023, Front Neurol 14:1118711), where antiplatelet at discharge IS the single protective covariate retained (aHR 0.59)."
    ),
    CONMED_ANTIDIABETIC = list(
      description = "Antidiabetic (ADM) prescribed; (15 + 301) / 4,204 = 7.52 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1 even for the non-diabetic stratum (these are patients without a registry diabetes diagnosis who nonetheless received an antidiabetic); no effect retained."
    ),
    CONMED_ACEI = list(
      description = "Angiotensin-converting-enzyme inhibitor prescribed; (37 + 1,154) / 4,204 = 28.33 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_BETABLOCKER = list(
      description = "Beta-blocker prescribed; (15 + 343) / 4,204 = 8.52 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_CCB = list(
      description = "Calcium-channel blocker prescribed; (21 + 736) / 4,204 = 18.01 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_DIURETIC = list(
      description = "Diuretic prescribed; (7 + 170) / 4,204 = 4.21 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4204L,
    n_events       = 138L,
    n_studies      = 1L,
    age_range      = "adults; median 62.9 years at the index stroke across the whole study population. In the non-diabetic stratum 57.97 percent of the patients who recurred and 60.45 percent of those who did not were older than 60 years (Table 1).",
    sex_female_pct = 39.49,
    race_ethnicity = "Multiethnic Malaysian registry cohort. Elhefnawy 2023 Table 1 reports Malay / Chinese / Indian / Others separately for the recurrent (54.34 / 1.44 / 0 / 44.20 percent) and non-recurrent (19.23 / 2.90 / 0.74 / 77.15 percent) non-diabetic groups; ethnicity was screened but not retained.",
    disease_state  = "Adults WITHOUT diabetes mellitus who had a first (index) ischemic stroke diagnosed by WHO criteria and confirmed by brain CT or MRI. Diabetes status was assigned from physician diagnosis, antidiabetic medication history, the patient's electronic record, or antidiabetic medications prescribed at discharge; patients with none of these formed this stratum. The endpoint is a subsequent ischemic stroke recorded in the registry after the index event.",
    dose_range     = "n/a (no drug exposure is modelled; secondary prevention enters only as the binary antihyperlipidemic indicator)",
    regions        = "Malaysia -- National Neurology Registry (NNEUR), a hospital-based registry covering 13 states; index strokes registered August 2009 to December 2016.",
    notes          = paste(
      "138 of 4,204 non-diabetic patients (3.28 percent) had a recurrent",
      "ischemic stroke within the maximum 7.37 years of follow-up.",
      "Baseline concurrent-disease prevalences in this stratum, computed from",
      "the Table 1 recurrent and non-recurrent counts: hypertension 58.59",
      "percent, hyperlipidemia 22.07 percent, ischemic heart disease 9.68",
      "percent, atrial fibrillation 4.21 percent, hyperuricemia 2.45 percent,",
      "current smoking 47.72 percent, family history of stroke 6.14 percent.",
      "Secondary-prevention prescribing: antihyperlipidemic 90.44 percent,",
      "antiplatelet 89.27 percent, ACE inhibitor 28.33 percent,",
      "calcium-channel blocker 18.01 percent, beta-blocker 8.52 percent,",
      "antidiabetic 7.52 percent, diuretic 4.21 percent.",
      "Unlike the diabetic stratum, this model received no temporal",
      "validation: the separate 2017-2020 NNEUR cohort of 1,262 patients was",
      "restricted to patients with diabetes (Methods, 'Temporal validation of",
      "the final TTE model among IS patients with DM'; Figure 4). The internal",
      "Kaplan-Meier visual predictive check in Figure 3 panel ii is the only",
      "predictive-performance evidence reported for this stratum.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # All values are Elhefnawy 2023 Table 3, "Estimated parameters of
    # the final model for recurrent IS after index IS among patients
    # with and without DM", non-DM columns. The hazard form is Table 2
    # row 4 (the 4-parameter interval Gompertz), h(t) = h0 *
    # exp(theta_y * t), and the Table 2 footnote fixes the switching:
    # "h0 equals theta1 if time < 0.5 years, theta3 if time >= 0.5;
    # theta_y equals theta2 if time < 3 years, theta4 if time >= 3
    # years." Covariates enter through the Methods equation printed
    # inline on page 2, h = h0 * exp(beta_t * t + beta_cov * (cov)).
    #
    # The two baseline-hazard scales are strictly positive rates and are
    # carried on the log scale; the two Gompertz shapes are exponents
    # that could in principle take either sign and are carried linearly.
    #
    # INTERNAL CONSISTENCY OF THE COVARIATE AND SHAPE BLOCKS -- all
    # reproduce the paper's own derived columns:
    #   exp( 1.030) = 2.8011  vs the published aHR 2.801
    #   exp( 0.789) = 2.2012  vs the published aHR 2.201
    #   exp(-0.835) = 0.4339  vs the published aHR 0.433
    #   ln(2)/1.70   = 0.4077 year = 4.89 months
    #                         vs Table 3 "0.40 (4.88 months)"
    #   ln(2)/0.213  = 3.2542 year
    #                         vs Table 3 "3.24years"
    # So the transcription below is right; the discrepancy documented
    # next is the publication's, not this file's.
    #
    # KNOWN DISCREPANCY IN THE BASELINE-HAZARD SCALE -- NOT TUNED.
    # The scale is too large to be reconciled with the paper's own event
    # count. Elhefnawy 2023 Table 1 reports 138 recurrences among 4,204
    # non-diabetic patients over at most 7.37 years, a cumulative hazard
    # of -log(1 - 138/4204) = 0.0334. The model as printed gives 0.3952
    # for a covariate-free patient at 7.37 years -- 11.8-fold too high
    # -- and 0.4592 (13.8-fold) after applying the cohort-average
    # covariate multiplier of 1.162. Every alternative reading also
    # fails:
    #   as published                       H(7.37) =  0.3952  ( 11.8x)
    #   constant scale, no Gompertz shape             0.1389  (  4.2x)
    #   decaying Gompertz (shapes negated)            0.0884  (  2.6x)
    #   the two scales swapped                       27.5214  (824.6x)
    #   Table 1 count                                 0.0334  (  1.0x)
    # The constant-scale row is the load-bearing one: theta1 = 0.253 per
    # year held flat over just the first 6 months contributes 0.1265 by
    # itself, already 3.8 times the entire 7.37-year observed cumulative
    # hazard, so the defect is in the scale and not in a misreading of
    # the piecewise structure. A frailty term cannot rescue it either: a
    # zero-mean random effect on the log hazard leaves at least half the
    # cohort at or above the typical hazard, which bounds the population
    # survivor function well below the published Kaplan-Meier curve.
    # The published values are transcribed verbatim regardless; the
    # vignette's Errata section carries the full arithmetic, asserted
    # numerically.
    # ==================================================================

    lh0_early_ris <- log(0.253)
    label("Log Gompertz baseline-hazard scale over the first 6 months after the index stroke, non-diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 1, non-DM column: theta1 (<6 months) = 0.253, RSE 24.38 percent (sampling importance resampling). Also Abstract and Results: "the index IS attack was predicted to contribute to the hazard of recurrent IS by 0.356 and 0.253 within the first six months after the index IS among patients with and without DM, respectively". See the scale discrepancy note above.

    lh0_late_ris <- log(0.0018)
    label("Log Gompertz baseline-hazard scale from 6 months after the index stroke onward, non-diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 2, non-DM column: theta3 (>=6 months) = 0.0018, RSE 23.09 percent. Also Results: "Even after six months of index IS, the baseline hazard of recurrent IS was not equal to zero among both groups (0.0023, 0.0018)."

    shape_early_ris <- 1.7
    label("Gompertz shape (hazard exponent) during the first 3 years after the index stroke, non-diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 3, non-DM column: alpha (<3) = theta2 = 1.7, RSE 6.37 percent, "Shape parameter in the first three years after index IS". Positive, so the hazard rises within each of the first two intervals.

    shape_late_ris <- 0.213
    label("Gompertz shape (hazard exponent) from 3 years after the index stroke onward, non-diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 4, non-DM column: alpha (>=3) = theta4 = 0.213, RSE 33.07 percent, "Shape parameter after three years of index IS". Printed as POSITIVE. Note the tension with the Results sentence "and then exponentially reduced afterwards": with theta4 = +0.213 the hazard still rises after 3 years, 8-fold more slowly, and the "reduction" is the downward jump at t = 3 that the shape switch produces (the hazard falls from 0.2952 to 0.00341 per year, an 86.6-fold drop). Encoded with the printed sign; see the vignette Errata.

    e_dis_hyperlip_ris <- 1.03
    label("Log-hazard coefficient for pre-index hyperlipidemia; aHR = exp(1.03) = 2.80")
    # Elhefnawy 2023 Table 3 row 6, non-DM column: theta5 = 1.03, RSE 16.51 percent, aHR 2.801 (95 percent CI 2.00-3.90). Results: "the recurrent IS rate among patients with HPLD and HTN was more than two times higher than that in those with no-HPLD or HTN prior to index IS (HR, 2.80; 95% CI (2.00-3.90))".

    e_dis_hypert_ris <- 0.789
    label("Log-hazard coefficient for pre-index hypertension; aHR = exp(0.789) = 2.20")
    # Elhefnawy 2023 Table 3 row 7, non-DM column: theta6 = 0.789, RSE 23.17 percent, aHR 2.201 (95 percent CI 1.53-3.14).

    e_conmed_lipidlower_ris <- -0.835
    label("Log-hazard coefficient for an antihyperlipidemic prescribed for secondary prevention; aHR = exp(-0.835) = 0.434")
    # Elhefnawy 2023 Table 3 row 8, non-DM column: theta7 = -0.835, RSE 24.80 percent, aHR 0.433 (the paper prints the interval descending, "0.433(0.65-0.28)"). Results: "receiving antihyperlipidemic medication among IS patients with non-DM decreased recurrent IS by 43.3%". This is the only protective covariate in either stratum, and its absence from the diabetic model is the paper's headline conclusion.

    # Interval boundaries. These are structural design choices, not
    # estimated quantities: Elhefnawy 2023 Table 2 counts the interval
    # Gompertz model as having 4 parameters (theta1-theta4), so the two
    # breakpoints carry no degrees of freedom and are fixed here.
    tbreak_scale_ris <- fixed(0.5)
    label("Time after the index stroke at which the baseline-hazard scale switches from theta1 to theta3 (year)")
    # Elhefnawy 2023 Table 2 footnote: "h0 equals theta1 if time < 0.5 years, theta3 if time >= 0.5". Table 3 labels the same split "(<6months)" / "(>=6months)".

    tbreak_shape_ris <- fixed(3)
    label("Time after the index stroke at which the Gompertz shape switches from theta2 to theta4 (year)")
    # Elhefnawy 2023 Table 2 footnote: "theta_y equals theta2 if time < 3 years, theta4 if time >= 3 years". Table 3 labels the same split "alpha (<3)" / "alpha (>=3)".

    # Between-subject variability. Elhefnawy 2023 reports NO random
    # effect for either stratum: Methods describes only the structural
    # hazard and the covariate search, and Table 3 has no variance, no
    # CV percent, no shrinkage and no omega row. This is encoded
    # faithfully as a typical-value model with no eta parameters -- no
    # variance is invented.
    #
    # The source fits this model with the parametric survival (event-
    # density) likelihood under LAPLACE (Methods: "the LAPACE (ADVAN=6
    # TOL=9 NSIG=3) estimation method"), so there is no observation-
    # error model to translate. This placeholder additive residual is
    # attached to the survivor-probability output so the nlmixr2
    # likelihood machinery accepts the model for forward simulation. It
    # is NOT from the source. Same device as
    # Lindauer_2017_lacosamide_dropout.R.
    addSd <- 0.001
    label("Placeholder additive residual error on the survivor-probability output sur (unitless); not from the source")
  })

  model({
    # --- Baseline-hazard scales, back-transformed.
    h0_early_ris <- exp(lh0_early_ris)
    h0_late_ris  <- exp(lh0_late_ris)

    # --- Covariate multiplier. Elhefnawy 2023 Methods, the equation
    # --- printed inline on page 2 and referenced as Equation 5:
    # ---   h = h0 * exp(beta_t * t + beta_cov * (cov))
    # --- All three retained covariates are 0/1 indicators, so exp() of
    # --- each coefficient is the adjusted hazard ratio the paper
    # --- tabulates in Table 3 and plots in Figure 2 (right panel).
    cov_ris <-
      exp(e_dis_hyperlip_ris      * DIS_HYPERLIP +
          e_dis_hypert_ris        * DIS_HYPERT +
          e_conmed_lipidlower_ris * CONMED_LIPIDLOWER)

    # --- Instantaneous baseline hazard, Elhefnawy 2023 Table 2 row 4
    # --- combined with the Table 2 footnote. The scale switches at
    # --- tbreak_scale_ris (0.5 year) and the shape switches
    # --- independently at tbreak_shape_ris (3 years), which gives three
    # --- regimes rather than two. Both switches are downward jumps in
    # --- the hazard (140.6-fold at 0.5 year, 86.6-fold at 3 years).
    if (t < tbreak_scale_ris) {
      h0_ris <- h0_early_ris * exp(shape_early_ris * t)
    } else if (t < tbreak_shape_ris) {
      h0_ris <- h0_late_ris * exp(shape_early_ris * t)
    } else {
      h0_ris <- h0_late_ris * exp(shape_late_ris * t)
    }
    hazard_ris <- h0_ris * cov_ris

    # --- Cumulative baseline hazard in closed form. Equation 1 defines
    # --- S(t) = exp(-integral of h over 0..t); the piecewise Gompertz
    # --- integrates elementally, so the integral is evaluated exactly
    # --- rather than by an ODE. Doing it this way keeps the two hazard
    # --- discontinuities out of the solver, where they would otherwise
    # --- be stepped over and silently smeared.
    # ---   segment 1, 0 .. 0.5 yr : h0_early * exp(shape_early * t)
    # ---   segment 2, 0.5 .. 3 yr : h0_late  * exp(shape_early * t)
    # ---   segment 3, 3 yr ..     : h0_late  * exp(shape_late  * t)
    seg1_full_ris <-
      h0_early_ris / shape_early_ris *
      (exp(shape_early_ris * tbreak_scale_ris) - 1)
    seg2_full_ris <-
      h0_late_ris / shape_early_ris *
      (exp(shape_early_ris * tbreak_shape_ris) -
         exp(shape_early_ris * tbreak_scale_ris))

    if (t < tbreak_scale_ris) {
      cumhaz0_ris <-
        h0_early_ris / shape_early_ris * (exp(shape_early_ris * t) - 1)
    } else if (t < tbreak_shape_ris) {
      cumhaz0_ris <-
        seg1_full_ris +
        h0_late_ris / shape_early_ris *
        (exp(shape_early_ris * t) - exp(shape_early_ris * tbreak_scale_ris))
    } else {
      cumhaz0_ris <-
        seg1_full_ris + seg2_full_ris +
        h0_late_ris / shape_late_ris *
        (exp(shape_late_ris * t) - exp(shape_late_ris * tbreak_shape_ris))
    }

    # --- The covariate multiplier is time-constant, so it factors
    # --- straight out of the integral (proportional hazards).
    cumhaz_ris <- cumhaz0_ris * cov_ris

    # --- Elhefnawy 2023 Equation 1: the probability of not experiencing
    # --- a recurrent ischemic stroke within [0, t].
    sur <- exp(-cumhaz_ris)

    sur ~ add(addSd)
  })
}
