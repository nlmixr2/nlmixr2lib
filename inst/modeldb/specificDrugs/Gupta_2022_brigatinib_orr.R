Gupta_2022_brigatinib_orr <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-EFFICACY model for confirmed ",
    "objective response (complete or partial response) by blinded ",
    "independent review committee (BIRC) as a function of the ",
    "brigatinib daily AUC time-averaged over the last scan interval ",
    "before the best confirmed response, in adults with ",
    "ALK-inhibitor-naive ALK-positive advanced non-small cell lung ",
    "cancer treated first-line with brigatinib 180 mg orally once daily ",
    "after a 7-day lead-in at 90 mg (Gupta 2022, ALTA-1L brigatinib ",
    "arm, n = 123). The probability is expit(2.173 + -0.03046 * ",
    "AUC_BRIG_SCAN), with the exposure in ug*h/mL/day entering linearly ",
    "and uncentred. The relationship is NOT statistically significant ",
    "(p = 0.108); the paper concludes that brigatinib exposure is not a ",
    "predictor of ORR at 180 mg once daily. The slope is the printed ",
    "odds ratio 0.97 per 1 ug*h/mL/day; the paper prints no intercept, ",
    "which was recovered from the vector-drawn fitted curve of Figure ",
    "2b with the slope held at the printed value. There is no PK layer ",
    "and no ODE: the exposure metric is supplied as a data column, ",
    "derived in the source analysis from the individual Bayesian CL/F ",
    "of the population PK model packaged as Gupta_2021_brigatinib. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Twenty companion exposure-response models ",
    "in the Gupta_2022_brigatinib_* family."
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
    dosing = "n/a (no dose events; exposure enters as the AUC_BRIG_SCAN covariate column)",
    concentration = "prob_orr_central (probability of a BIRC-confirmed objective response, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    AUC_BRIG_SCAN = list(
      description = "Individual brigatinib daily AUC averaged over the disease-assessment scan interval -- in the static models, the interval between the last two scans preceding the event, best confirmed response or censoring. Supplied as data: this model has no PK layer, and the source analysis derived it from the individual Bayesian CL/F of the Gupta 2021 population PK model and the actual dosing history (Gupta 2022 Methods S1).",
      units = "ug*h/mL/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Daily AUC on day i is the increment of the cumulative AUC,",
        "AUCD(i) = CAUC(i + 1) - CAUC(i), with CAUC(i) = cumulative dose to",
        "day i / individual CL/F (Gupta 2022 Methods S1); the metric is the",
        "mean of AUCD over the scan interval, so dose reductions and",
        "interruptions inside the interval lower it. It is NOT a",
        "steady-state AUC(0-24) at the labelled dose, although for an",
        "uninterrupted 180 mg once-daily course the two coincide",
        "(180 mg / CL/F). Enters LINEARLY and uncentred.",
        "Gupta 2022 Figure 2a gives the PFS analysis-set median (range)",
        "as 15.2 (0, 77.0) and the quartile medians as 9.7, 13.0, 17.1 and",
        "28.3 ug*h/mL/day.",
        "In this model the exposure enters with slope log(0.97) per ug*h/mL/day."
      ),
      source_name = "Time-averaged AUC between last two scans (ug.h/mL/day)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 123L,
    n_studies = 1L,
    n_observations = "123 binary best-confirmed-response records, one per patient",
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
      "blinded independent review) and safety. Observed ",
      "BIRC-confirmed ORR in the analysis set 101/123 (82.1 percent; ",
      "Figure 2b quartile counts 27/31, 25/31, 25/30, 24/31)."
    )
  )

  ini({
    # ==================================================================
    # Gupta 2022 Figure 2b: logit(p) = logit_ref + e_auc_logit *
    # AUC_BRIG_SCAN.
    #
    # Recovered from Figure 2b, whose fitted curve (annotated 'Exposure
    # effect P = 0.108') is drawn as vector paths in the publisher PDF:
    # every Bezier node of the curve lies exactly on it, so the node
    # coordinates are read without pixel measurement.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- 2.173; label("Logit of the probability of a BIRC-confirmed objective response at zero brigatinib exposure (unitless logit)") # figure-derived: Figure 2b fitted curve, intercept recovered with the slope held at the printed odds ratio; predicts 101.4 of 101 observed events

    # ----- Exposure slope on the logit -----
    e_auc_logit <- log(0.97); label("Log-odds of a BIRC-confirmed objective response per 1 ug*h/mL/day increase in AUC_BRIG_SCAN (unitless logit)") # odds ratio 0.97 (95% CI 0.93-1.01), p = 0.108 (Gupta 2022 Results, Exposure-efficacy analysis); Figure 2b panel annotates P = 0.108

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Exposure enters linearly and uncentred, so logit_ref is the logit at
    # zero exposure, an extrapolation below the observed range.
    logit_orr_central <- logit_ref + e_auc_logit * AUC_BRIG_SCAN

    prob_orr_central <- expit(logit_orr_central)

    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
