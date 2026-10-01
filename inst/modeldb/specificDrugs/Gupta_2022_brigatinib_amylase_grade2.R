Gupta_2022_brigatinib_amylase_grade2 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-SAFETY model for grade >= 2 ",
    "increased amylase as a function of the brigatinib daily AUC ",
    "time-averaged from the first dose to the first occurrence of the ",
    "event (or to the end of treatment), in adults with ",
    "ALK-inhibitor-naive ALK-positive advanced non-small cell lung ",
    "cancer treated first-line with brigatinib 180 mg orally once daily ",
    "after a 7-day lead-in at 90 mg (Gupta 2022, ALTA-1L brigatinib ",
    "arm, n = 123). The probability is expit(-2.225 + 0.03922 * ",
    "AUC_BRIG_EVT), with the exposure in ug*h/mL/day entering linearly ",
    "and uncentred. A weak, non-significant trend; with the alternative ",
    "early-exposure metric of days 8-14 the same endpoint becomes ",
    "significant (Gupta_2022_brigatinib_amylase_grade2_d8_14). The ",
    "slope is the printed odds ratio 1.04 per 1 ug*h/mL/day; the paper ",
    "prints no intercept, which was recovered from the vector-drawn ",
    "fitted curve of Supporting Figure S7g with the slope held at the ",
    "printed value. There is no PK layer and no ODE: the exposure ",
    "metric is supplied as a data column, derived in the source ",
    "analysis from the individual Bayesian CL/F of the population PK ",
    "model packaged as Gupta_2021_brigatinib. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Twenty companion exposure-response models in the ",
    "Gupta_2022_brigatinib_* family."
  )
  reference <- paste(
    "Gupta N, Reckamp KL, Camidge DR, Kleijn HJ, Ouerdani A, Bellanti F,",
    "Maringwa J, Hanley MJ, Wang S, Zhang P, Venkatakrishnan K. (2022).",
    "Population pharmacokinetic and exposure-response analyses from ALTA-1L:",
    "Model-based analyses supporting the brigatinib dose in ALK-positive NSCLC.",
    "Clin Transl Sci 15:1143-1154.",
    "doi:10.1111/cts.13231.",
    "Individual brigatinib exposures derive from the previously published",
    "population PK model applied without modification by Bayesian",
    "re-estimation; see modellib('Gupta_2021_brigatinib').",
    "Supporting Figures S6 and S7 are Supporting Information files s002 and",
    "s005; Methods S1 is file s007."
  )
  vignette <- "Gupta_2022_brigatinib_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_BRIG_EVT covariate column)",
    concentration = "prob_amylase_increase_grade2 (probability of grade >= 2 amylase increase, 0-1; also logit_amylase_increase_grade2)"
  )

  covariateData <- list(
    AUC_BRIG_EVT = list(
      description = "Individual brigatinib daily AUC time-averaged from the first dose to the first occurrence of the adverse event, or to the end of treatment for a patient without the event. Supplied as data: this model has no PK layer, and the source analysis derived it from the individual Bayesian CL/F of the Gupta 2021 population PK model and the actual dosing history (Gupta 2022 Methods, Exposure-safety analyses).",
      units = "ug*h/mL/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cumulative dose to the event (or to the end of treatment) divided",
        "by the individual CL/F and by the days elapsed (Gupta 2022",
        "Methods S1), so the averaging window differs per endpoint and per",
        "patient: the same patient carries a different value in each",
        "adverse-event model. Includes the 7-day 90 mg lead-in, so an early",
        "event is attached to a LOWER exposure than the 180 mg",
        "maintenance AUC. Enters LINEARLY and uncentred. Gupta 2022",
        "Figures S6 and S7 plot the observed values over roughly 6 to 70",
        "ug*h/mL/day.",
        "In this model the exposure enters with slope log(1.04) per ug*h/mL/day."
      ),
      source_name = "Time-Averaged Exposure (ug.h/mL/day)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 123L,
    n_studies = 1L,
    n_observations = "123 binary grade >= 2 amylase increase records, one per patient (events from the first dose to 30 days after the last dose, CTCAE v4.03)",
    age_range = "27-85 years (median 57; Gupta 2022 Table 1, PK-evaluable ALTA-1L brigatinib arm, n = 123)",
    weight_range = "43-111 kg (median 67)",
    sex_female_pct = 48.8,
    race_ethnicity = c(White = 55.3, Asian = 43.1, Other = 1.6),
    disease_state = "ALK-inhibitor-naive ALK-positive advanced non-small cell lung cancer, first-line treatment in the phase 3 ALTA-1L trial (NCT02737501), second interim analysis (data cut-off 28 June 2019). ECOG performance status 0 / 1 / 2 in 39.8 / 56.1 / 4.1 percent; prior chemotherapy 28.5 percent; brain metastases at baseline 29.3 percent; baseline albumin median 41 g/L (range 24-48).",
    dose_range = "Brigatinib 180 mg orally once daily after a 7-day lead-in at 90 mg once daily (28-day cycles); the crizotinib comparator arm (250 mg twice daily) contributes no exposure and is not modelled.",
    regions = "Multiregional (Gupta 2022 Discussion: ALTA-1L was a multiregional clinical trial); per-region counts not reported.",
    notes = paste0(
      "Of 137 intent-to-treat brigatinib patients, 13 had no ",
      "quantifiable brigatinib concentration and one was not dosed, ",
      "leaving 123 for the population PK re-estimation (1069 samples) ",
      "and the exposure-response analyses of PFS, ORR (both by ",
      "blinded independent review) and safety. Observed incidence ",
      "24/123 (19.5 percent); Supporting Figure S7g quartile counts ",
      "5/31, 6/31, 4/30, 9/31 (events/patients)."
    )
  )

  ini({
    # ==================================================================
    # Gupta 2022 Supporting Figure S7g: logit(p) = logit_ref + e_auc_logit *
    # AUC_BRIG_EVT.
    #
    # Recovered from Supporting Figure S7g, whose fitted curve (annotated
    # 'Exposure effect P = 0.0943') is drawn as vector paths in the
    # publisher PDF: every Bezier node of the curve lies exactly on it, so
    # the node coordinates are read without pixel measurement.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -2.225; label("Logit of the probability of grade >= 2 amylase increase at zero brigatinib exposure (unitless logit)") # figure-derived: Supporting Figure S7g fitted curve, intercept recovered with the slope held at the printed odds ratio; predicts 23.5 of 24 observed events

    # ----- Exposure slope on the logit -----
    e_auc_logit <- log(1.04); label("Log-odds of grade >= 2 amylase increase per 1 ug*h/mL/day increase in AUC_BRIG_EVT (unitless logit)") # odds ratio 1.04 (95% CI 0.99-1.09), p = 0.094 (Gupta 2022 Results, Exposure-safety analysis); Supporting Figure S7g panel annotates P = 0.0943

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_amylase_increase_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Exposure enters linearly and uncentred, so logit_ref is the logit at
    # zero exposure, an extrapolation below the observed range.
    logit_amylase_increase_grade2 <- logit_ref + e_auc_logit * AUC_BRIG_EVT

    prob_amylase_increase_grade2 <- expit(logit_amylase_increase_grade2)

    prob_amylase_increase_grade2 ~ add(addSd_prob_amylase_increase_grade2)
  })
}
