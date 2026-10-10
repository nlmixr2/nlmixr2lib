Gupta_2022_brigatinib_bradycardia_grade2 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-SAFETY model for grade >= 2 ",
    "bradycardia as a function of the brigatinib daily AUC ",
    "time-averaged from the first dose to the first occurrence of the ",
    "event (or to the end of treatment), in adults with ",
    "ALK-inhibitor-naive ALK-positive advanced non-small cell lung ",
    "cancer treated first-line with brigatinib 180 mg orally once daily ",
    "after a 7-day lead-in at 90 mg (Gupta 2022, ALTA-1L brigatinib ",
    "arm, n = 123). The probability is expit(-2.807 + -0.05964 * ",
    "AUC_BRIG_EVT), with the exposure in ug*h/mL/day entering linearly ",
    "and uncentred. The relationship is NOT statistically significant ",
    "(P = 0.522 annotated on Supporting Figure S7c); the model is ",
    "packaged so that the null result is reproducible rather than ",
    "merely asserted. The paper prints no coefficient for this ",
    "regression; intercept and slope were both recovered from the ",
    "vector-drawn fitted curve of Supporting Figure S7c. There is no PK ",
    "layer and no ODE: the exposure metric is supplied as a data ",
    "column, derived in the source analysis from the individual ",
    "Bayesian CL/F of the population PK model packaged as ",
    "Gupta_2021_brigatinib. No between-subject random effect and no ",
    "residual error are estimated (Bernoulli likelihood). Twenty ",
    "companion exposure-response models in the Gupta_2022_brigatinib_* ",
    "family."
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
    concentration = "prob_bradycardia_grade2 (probability of grade >= 2 bradycardia, 0-1; also logit_bradycardia_grade2)"
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
        "In this model the exposure enters with slope -0.05964 per ug*h/mL/day (figure-derived)."
      ),
      source_name = "Time-Averaged Exposure (ug.h/mL/day)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 123L,
    n_studies = 1L,
    n_observations = "123 binary grade >= 2 bradycardia records, one per patient (events from the first dose to 30 days after the last dose, CTCAE v4.03)",
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
      "2/123 (1.6 percent); Supporting Figure S7c quartile counts ",
      "1/31, 0/31, 1/30, 0/31 (events/patients)."
    )
  )

  ini({
    # ==================================================================
    # Gupta 2022 Supporting Figure S7c: logit(p) = logit_ref + e_auc_logit *
    # AUC_BRIG_EVT.
    #
    # Recovered from Supporting Figure S7c, whose fitted curve (annotated
    # 'Exposure effect P = 0.522') is drawn as vector paths in the publisher
    # PDF: every Bezier node of the curve lies exactly on it, so the node
    # coordinates are read without pixel measurement.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -2.807; label("Logit of the probability of grade >= 2 bradycardia at zero brigatinib exposure (unitless logit)") # figure-derived: Supporting Figure S7c fitted curve, intercept recovered; predicts 2.5 of 2 observed events

    # ----- Exposure slope on the logit -----
    e_auc_logit <- -0.05964; label("Log-odds of grade >= 2 bradycardia per 1 ug*h/mL/day increase in AUC_BRIG_EVT (unitless logit)") # figure-derived: Supporting Figure S7c fitted curve (no coefficient printed); panel annotates P = 0.522

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_bradycardia_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Exposure enters linearly and uncentred, so logit_ref is the logit at
    # zero exposure, an extrapolation below the observed range.
    logit_bradycardia_grade2 <- logit_ref + e_auc_logit * AUC_BRIG_EVT

    prob_bradycardia_grade2 <- expit(logit_bradycardia_grade2)

    prob_bradycardia_grade2 ~ add(addSd_prob_bradycardia_grade2)
  })
}
