Elhefnawy_2023_recurrent_ischemic_stroke_dm <- function() {
  description <- paste(
    "Parametric time-to-event model for recurrent ischemic stroke (IS) after a",
    "first (index) IS, in the DIABETIC stratum of the National Neurology",
    "Registry of Malaysia (3,493 patients with diabetes mellitus, 195",
    "recurrences, up to 7.37 years of follow-up), developed in NONMEM 7.5.",
    "The hazard is Gompertz, h(t) = theta_x * exp(theta_y * t), with the scale",
    "switching at 0.5 year (theta1 -> theta3) and the shape switching",
    "independently at 3 years (theta2 -> theta4). Two baseline comorbidity",
    "covariates act log-linearly on the hazard, both raising it: ischemic heart",
    "disease and hyperlipidemia recorded before the index stroke. Time is in",
    "years and there is no drug exposure term, so this is a disease-progression",
    "/ event-risk model rather than a PK-PD model. The cumulative hazard is",
    "evaluated in closed form (the piecewise Gompertz integral is elementary)",
    "rather than by an ODE, which avoids integrating across the two hazard",
    "discontinuities; the model exposes the instantaneous hazard `hazard_ris`,",
    "the cumulative hazard `cumhaz_ris` and the recurrence-free survivor",
    "probability `sur`. Sister model file for the non-diabetic stratum of the",
    "same paper: modellib('Elhefnawy_2023_recurrent_ischemic_stroke_nondm').",
    "IMPORTANT -- the published baseline-hazard scale does not reproduce the",
    "paper's own event count. As printed, the model predicts a 37.6 percent",
    "chance of recurrence within 7.37 years for a covariate-free diabetic",
    "patient, against the 195 of 3,493 = 5.58 percent actually reported in",
    "Elhefnawy 2023 Table 1: an 8.2-fold overprediction on the cumulative-hazard",
    "scale, rising to 12.5-fold once the cohort covariate distribution is",
    "applied. No reading of the piecewise structure removes it -- theta1 = 0.356",
    "per year held flat over just the first 6 months already contributes 0.178,",
    "3.1 times the entire 7.37-year observed cumulative hazard of 0.0574. The",
    "covariate block, by contrast, reproduces the paper's own adjusted hazard",
    "ratios and half-life column exactly. The values are transcribed as",
    "published and have NOT been tuned; see the validation vignette's Errata",
    "section for the full arithmetic. The same defect, at the same magnitude,",
    "is present in this group's companion pooled-cohort publication",
    "(Elhefnawy 2023, Front Neurol 14:1118711, PMC10176964).",
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
    "estimates are Table 3, DM columns. The EuropePMC open-access deposit for",
    "PMC10796130 contains no supplementary files.",
    "Sister model file from the same paper:",
    "modellib('Elhefnawy_2023_recurrent_ischemic_stroke_nondm').",
    sep = " "
  )
  vignette <- "Elhefnawy_2023_recurrent_ischemic_stroke_diabetes"

  units <- list(
    time          = "year",
    dosing        = "n/a (no dosing events; no drug exposure term enters this stratum's final model)",
    concentration = "n/a (the model outputs are a hazard in 1/year, a unitless cumulative hazard and a unitless recurrence-free survivor probability, not a drug concentration)"
  )

  covariateData <- list(
    DIS_IHD = list(
      description        = "1 = the patient carried an ischemic heart disease (IHD) diagnosis before the index ischemic stroke; 0 = no IHD. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ischemic heart disease before the index stroke)",
      notes              = "Ascertained from the National Neurology Registry concurrent-disease fields (Elhefnawy 2023 Methods, 'Data collection': 'Patients' demographic data and concurrent disease data, including HPLD, hypertension (HTN), IHD, hyperuricemia, atrial fibrillation (AF)... were investigated'). Prevalence in the 3,493-patient diabetic stratum is (52 + 420) / 3,493 = 13.51 percent (Table 1). Retained in the diabetic stratum only: aHR = exp(0.876) = 2.40 (95 percent CI 1.79-3.20), the strongest single predictor in this model. Elhefnawy 2023 Discussion attributes the effect to the atherosclerotic pathophysiology shared by IHD and ischemic stroke, made more prominent by diabetes-associated angiopathy.",
      source_name        = "IHD"
    ),
    DIS_HYPERLIP = list(
      description        = "1 = the patient carried a hyperlipidemia (HPLD) diagnosis before the index ischemic stroke; 0 = no hyperlipidemia. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hyperlipidemia before the index stroke)",
      notes              = "Ascertained the same way as the other concurrent-disease flags. Prevalence in the diabetic stratum is (96 + 1,004) / 3,493 = 31.49 percent (Table 1). aHR = exp(0.633) = 1.88 (95 percent CI 1.44-2.45). The same covariate is retained, more strongly, in the non-diabetic stratum (aHR 2.80) -- see modellib('Elhefnawy_2023_recurrent_ischemic_stroke_nondm').",
      source_name        = "HPLD"
    )
  )

  # Covariates that Elhefnawy 2023 collected and screened for the diabetic
  # stratum but did not retain in the final model. Listed here for provenance
  # only; none is referenced in model(). The paper reports no point estimate
  # for any of them -- the covariate screen itself is reported only by
  # reference to a previously published preprint (Methods, 'Covariate Model
  # Development and Model Evaluation': "The covariate model development and
  # model evaluation were described comprehensively in a previously published
  # preprint [23]"), and the EuropePMC deposit for PMC10796130 carries no
  # supplementary file, so no dOFV is available here either.
  #
  # Two of these names are worth calling out. HTN is absent from the diabetic
  # stratum although it is retained in the non-diabetic one; Elhefnawy 2023
  # Discussion reads this as collinearity ("The co-existence of HTN could
  # explain the non-significance of HTN among DM patients"), and hypertension
  # prevalence in this stratum is 87.1 percent, leaving little contrast.
  # Antihyperlipidemic therapy is likewise absent here but retained (and
  # protective) in the non-diabetic stratum -- the paper's headline negative
  # finding, stated in the Abstract Conclusion and the Conclusions section.
  #
  # These names are documentation labels for the paper's screen, not registered
  # canonical covariate columns: none is used in model(), so none is added to
  # inst/references/covariate-columns.md.
  covariatesDataExcluded <- list(
    DIS_HYPERT = list(
      description = "Hypertension (HTN) before the index stroke; (180 + 2,863) / 3,493 = 87.12 percent of the diabetic stratum.",
      units = "(binary)", type = "binary",
      notes = "Screened but not retained in the diabetic stratum, while it IS retained in the non-diabetic stratum (aHR 2.20). Elhefnawy 2023 Discussion: 'The co-existence of HTN could explain the non-significance of HTN among DM patients in this study.' Table 1 additionally stratifies hypertension duration at 5 years; no duration effect was retained either. Registered canonical DIS_HYPERT is used for it in the sister non-DM model file."
    ),
    CONMED_LIPIDLOWER = list(
      description = "Antihyperlipidemic medication prescribed at discharge for secondary prevention; (167 + 2,926) / 3,493 = 88.55 percent of the diabetic stratum.",
      units = "(binary)", type = "binary",
      notes = "The paper's headline negative result. Retained and protective in the non-diabetic stratum (aHR 0.433) but NOT significant among patients with diabetes: Abstract Conclusion, 'receiving medications for secondary prevention failed to demonstrate a significant association with reducing IS recurrence among IS patients with DM'. Elhefnawy 2023 Discussion notes agreement with Zhang et al. and raises statin-driven worsening of insulin resistance as a candidate mechanism. Registered canonical CONMED_LIPIDLOWER is used for it in the sister non-DM model file."
    ),
    DIS_HYPERURICEMIA = list(
      description = "Hyperuricemia before the index stroke; (10 + 121) / 3,493 = 3.75 percent.",
      units = "(binary)", type = "binary",
      notes = "Named among the investigated concurrent diseases in Methods, 'Data collection' and tabulated in Table 1; no effect retained in Table 3."
    ),
    DIS_AF = list(
      description = "Atrial fibrillation before the index stroke; (4 + 87) / 3,493 = 2.61 percent.",
      units = "(binary)", type = "binary",
      notes = "Named among the investigated concurrent diseases and tabulated in Table 1; no effect retained. Table 1 shows the direction opposite to the usual clinical expectation in this stratum (2.05 percent of recurrent vs 2.63 percent of non-recurrent patients)."
    ),
    AGE = list(
      description = "Age at the index ischemic stroke; median 62.9 years across the whole study, dichotomised at 60 years in Table 1.",
      units = "year", type = "continuous",
      notes = "Screened as a demographic covariate (Methods, 'Data collection'); not retained. Table 1: 52.82 percent of recurrent and 60.06 percent of non-recurrent diabetic patients were older than 60 years."
    ),
    SEXF = list(
      description = "Female sex; (101 + 1,647) / 3,493 = 50.04 percent of the diabetic stratum.",
      units = "(binary)", type = "binary",
      notes = "Screened as a demographic covariate; not retained. Table 1 shows near-identical proportions in the recurrent (51.79 percent) and non-recurrent (49.93 percent) diabetic groups."
    ),
    SMOKER = list(
      description = "Current smoker at the index stroke; (113 + 1,630) / 3,493 = 49.90 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1 with a sizeable unadjusted imbalance (57.94 percent of recurrent vs 49.42 percent of non-recurrent diabetic patients); no effect retained in Table 3."
    ),
    FAMHX_STROKE = list(
      description = "Family history of stroke; (16 + 152) / 3,493 = 4.81 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    NIHSS = list(
      description = "National Institutes of Health Stroke Scale severity of the index stroke, dichotomised by the paper into minor vs moderate/severe.",
      units = "(score)", type = "continuous",
      notes = "Tabulated in Table 1 (41.53 percent minor among recurrent diabetic patients vs 46.54 percent among non-recurrent); no NIHSS term appears in Table 3."
    ),
    DUR_DIAB = list(
      description = "Duration of diabetes before the index stroke, banded by the paper into <1, 1-5, 6-10 and >10 years.",
      units = "year", type = "categorical",
      notes = "Tabulated in Table 1 for the diabetic stratum only; no duration effect appears in Table 3. This is the one screened covariate with no counterpart in the non-diabetic stratum."
    ),
    RACE = list(
      description = "Ethnicity, recorded by the registry as Malay / Chinese / Indian / Others.",
      units = "(categorical)", type = "categorical",
      notes = "Tabulated in Table 1 with a large unadjusted imbalance (41.02 percent Malay among recurrent diabetic patients vs 21.13 percent among non-recurrent); no ethnicity effect appears in Table 3. Elhefnawy 2023 Results attributes the large 'Others' stratum to the East Malaysian hospitals contributing most of the registry data."
    ),
    CONMED_ANTIPLATELET = list(
      description = "Antiplatelet (APLT) prescribed at discharge; (167 + 2,978) / 3,493 = 90.04 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1 among the secondary-prevention medications; no effect retained in Table 3 for either stratum. Contrast the companion pooled-cohort publication (Elhefnawy 2023, Front Neurol 14:1118711), where antiplatelet at discharge IS the single protective covariate retained (aHR 0.59)."
    ),
    CONMED_ANTIDIABETIC = list(
      description = "Antidiabetic (ADM) prescribed at discharge; (117 + 2,005) / 3,493 = 60.75 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained in Table 3."
    ),
    CONMED_ACEI = list(
      description = "Angiotensin-converting-enzyme inhibitor prescribed at discharge; (61 + 1,144) / 3,493 = 34.50 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_BETABLOCKER = list(
      description = "Beta-blocker prescribed at discharge; (24 + 407) / 3,493 = 12.34 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_CCB = list(
      description = "Calcium-channel blocker prescribed at discharge; (59 + 784) / 3,493 = 24.13 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    ),
    CONMED_DIURETIC = list(
      description = "Diuretic prescribed at discharge; (22 + 255) / 3,493 = 7.93 percent.",
      units = "(binary)", type = "binary",
      notes = "Tabulated in Table 1; no effect retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3493L,
    n_events       = 195L,
    n_studies      = 1L,
    age_range      = "adults; median 62.9 years at the index stroke across the whole study population. In the diabetic stratum 52.82 percent of the patients who recurred and 60.06 percent of those who did not were older than 60 years (Table 1).",
    sex_female_pct = 50.04,
    race_ethnicity = "Multiethnic Malaysian registry cohort. Elhefnawy 2023 Table 1 reports Malay / Chinese / Indian / Others separately for the recurrent (41.02 / 3.07 / 1.53 / 54.35 percent) and non-recurrent (21.13 / 2.63 / 1.51 / 74.71 percent) diabetic groups; ethnicity was screened but not retained.",
    disease_state  = "Adults with diabetes mellitus and a first (index) ischemic stroke diagnosed by WHO criteria and confirmed by brain CT or MRI. Diabetes was identified from physician diagnosis, antidiabetic medication history, the patient's electronic record, or antidiabetic medications prescribed at discharge. The endpoint is a subsequent ischemic stroke recorded in the registry after the index event.",
    dose_range     = "n/a (no drug exposure is modelled; secondary-prevention prescribing was screened as binary indicators and none was retained in this stratum)",
    regions        = "Malaysia -- National Neurology Registry (NNEUR), a hospital-based registry covering 13 states; index strokes registered August 2009 to December 2016.",
    notes          = paste(
      "195 of 3,493 diabetic patients (5.58 percent) had a recurrent ischemic",
      "stroke within the maximum 7.37 years of follow-up. Note that the paper's",
      "Abstract reports this proportion as 5.82 percent and its Results section",
      "as 5.55 percent; 195 / 3,493 is 5.58 percent, so both printed",
      "percentages are slightly off and the counts are used here.",
      "Baseline concurrent-disease prevalences in the diabetic stratum,",
      "computed from the Table 1 recurrent and non-recurrent counts:",
      "hypertension 87.12 percent, hyperlipidemia 31.49 percent, ischemic",
      "heart disease 13.51 percent, hyperuricemia 3.75 percent, atrial",
      "fibrillation 2.61 percent, current smoking 49.90 percent, family",
      "history of stroke 4.81 percent. Secondary-prevention prescribing at",
      "discharge: antiplatelet 90.04 percent, antihyperlipidemic 88.55",
      "percent, antidiabetic 60.75 percent, ACE inhibitor 34.50 percent,",
      "calcium-channel blocker 24.13 percent, beta-blocker 12.34 percent,",
      "diuretic 7.93 percent.",
      "The diabetic stratum is the one this paper additionally validated",
      "temporally, on a separate NNEUR cohort of 1,262 diabetic patients with",
      "index strokes registered January 2017 to December 2020",
      "(Methods, 'Temporal validation'; Figure 4). No parameter estimates are",
      "reported for that validation cohort, so it is described here but not",
      "encoded.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # All values are Elhefnawy 2023 Table 3, "Estimated parameters of
    # the final model for recurrent IS after index IS among patients
    # with and without DM", DM columns. The hazard form is Table 2 row 4
    # (the 4-parameter interval Gompertz), h(t) = h0 * exp(theta_y * t),
    # and the Table 2 footnote fixes the switching: "h0 equals theta1 if
    # time < 0.5 years, theta3 if time >= 0.5; theta_y equals theta2 if
    # time < 3 years, theta4 if time >= 3 years." Covariates enter
    # through the Methods equation printed inline on page 2,
    # h = h0 * exp(beta_t * t + beta_cov * (cov)).
    #
    # The two baseline-hazard scales are strictly positive rates and are
    # carried on the log scale; the two Gompertz shapes are exponents
    # that could in principle take either sign and are carried linearly.
    #
    # INTERNAL CONSISTENCY OF THE COVARIATE AND SHAPE BLOCKS -- all
    # reproduce the paper's own derived columns:
    #   exp(0.876) = 2.4013  vs the published aHR 2.40
    #   exp(0.633) = 1.8833  vs the published aHR 1.88
    #   ln(2)/1.58  = 0.4387 year = 5.26 months
    #                        vs Table 3 "0.43 (5.25 months)"
    #   ln(2)/0.242 = 2.8642 year
    #                        vs Table 3 "2.85 years"
    # So the transcription below is right; the discrepancy documented
    # next is the publication's, not this file's.
    #
    # KNOWN DISCREPANCY IN THE BASELINE-HAZARD SCALE -- NOT TUNED.
    # The scale is too large to be reconciled with the paper's own event
    # count. Elhefnawy 2023 Table 1 reports 195 recurrences among 3,493
    # diabetic patients over at most 7.37 years, a cumulative hazard of
    # -log(1 - 195/3493) = 0.0574. The model as printed gives 0.4714 for
    # a covariate-free patient at 7.37 years -- 8.2-fold too high -- and
    # 0.7166 (12.5-fold) after applying the cohort-average covariate
    # multiplier of 1.520. Every alternative reading also fails:
    #   as published                       H(7.37) = 0.4714  ( 8.2x)
    #   constant scale, no Gompertz shape            0.1938  ( 3.4x)
    #   decaying Gompertz (shapes negated)           0.1267  ( 2.2x)
    #   the two scales swapped                      31.0030  (540x)
    #   Table 1 count                                0.0574  ( 1.0x)
    # The constant-scale row is the load-bearing one: theta1 = 0.356/year
    # held flat over just the first 6 months contributes 0.178 by itself,
    # already 3.1 times the entire 7.37-year observed cumulative hazard,
    # so the defect is in the scale and not in a misreading of the
    # piecewise structure. A frailty term cannot rescue it either: a
    # zero-mean random effect on the log hazard leaves at least half the
    # cohort at or above the typical hazard, which bounds the population
    # survivor function well below the published Kaplan-Meier curve.
    # The paper's own clinical-calculator scenarios do not reproduce
    # either, and are internally inconsistent (its scenario-1 recurrence
    # probability FALLS from 8.396 percent at 1 year to 4.917 percent at
    # 4 years, which is impossible for a cumulative probability).
    # The published values are transcribed verbatim regardless; the
    # vignette's Errata section carries the full arithmetic, asserted
    # numerically.
    # ==================================================================

    lh0_early_ris <- log(0.356)
    label("Log Gompertz baseline-hazard scale over the first 6 months after the index stroke, diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 1, DM column: theta1 (<6 months) = 0.356, RSE 13.69 percent (sampling importance resampling). Also Abstract and Results: "the index IS attack was predicted to contribute to the hazard of recurrent IS by 0.356 ... within the first six months after the index IS among patients with ... DM". See the scale discrepancy note above.

    lh0_late_ris <- log(0.0023)
    label("Log Gompertz baseline-hazard scale from 6 months after the index stroke onward, diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 2, DM column: theta3 (>=6 months) = 0.0023, RSE 17.01 percent. Also Results: "Even after six months of index IS, the baseline hazard of recurrent IS was not equal to zero among both groups (0.0023, 0.0018)."

    shape_early_ris <- 1.58
    label("Gompertz shape (hazard exponent) during the first 3 years after the index stroke, diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 3, DM column: alpha (<3) = theta2 = 1.58, RSE 5.98 percent, "Shape parameter in the first three years after index IS". Positive, so the hazard rises within each of the first two intervals; Results: "the recurrent hazard increased exponentially during the first three years after the index IS".

    shape_late_ris <- 0.242
    label("Gompertz shape (hazard exponent) from 3 years after the index stroke onward, diabetic stratum (1/year)")
    # Elhefnawy 2023 Table 3 row 4, DM column: alpha (>=3) = theta4 = 0.242, RSE 22.36 percent, "Shape parameter after three years of index IS". Printed as POSITIVE. Note the tension with the Results sentence "and then exponentially reduced afterwards": with theta4 = +0.242 the hazard still rises after 3 years, 6.5-fold more slowly, and the "reduction" is the downward jump at t = 3 that the shape switch produces (the hazard falls from 0.2632 to 0.00475 per year, a 55-fold drop). Encoded with the printed sign; see the vignette Errata.

    e_dis_ihd_ris <- 0.876
    label("Log-hazard coefficient for pre-index ischemic heart disease; aHR = exp(0.876) = 2.40")
    # Elhefnawy 2023 Table 3 row 5, DM column: theta5 = 0.876, RSE 16.88 percent, aHR 2.40 (95 percent CI 1.79-3.20). Results: "the recurrent IS rate among DM patients with IHD was 2.40 times higher than that in patients with no-IHD prior to index IS".

    e_dis_hyperlip_ris <- 0.633
    label("Log-hazard coefficient for pre-index hyperlipidemia; aHR = exp(0.633) = 1.88")
    # Elhefnawy 2023 Table 3 row 6, DM column: theta6 = 0.633, RSE 21.38 percent, aHR 1.88 (95 percent CI 1.44-2.45). Results: "the presence of HPLD prior to index IS increased the recurrent IS rate among DM by 88%".

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
    # variance is invented. (The companion pooled-cohort publication,
    # Front Neurol 14:1118711, does state that between-subject
    # variability was estimated but likewise never reports its
    # magnitude; that model file carries fixed(0) etas for it. Here
    # there is nothing to carry.)
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
    # --- Both retained covariates are 0/1 indicators, so exp() of each
    # --- coefficient is the adjusted hazard ratio the paper tabulates
    # --- in Table 3 and plots in Figure 2 (left panel).
    cov_ris <-
      exp(e_dis_ihd_ris      * DIS_IHD +
          e_dis_hyperlip_ris * DIS_HYPERLIP)

    # --- Instantaneous baseline hazard, Elhefnawy 2023 Table 2 row 4
    # --- combined with the Table 2 footnote. The scale switches at
    # --- tbreak_scale_ris (0.5 year) and the shape switches
    # --- independently at tbreak_shape_ris (3 years), which gives three
    # --- regimes rather than two. Both switches are downward jumps in
    # --- the hazard (154.8-fold at 0.5 year, 55.4-fold at 3 years).
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
