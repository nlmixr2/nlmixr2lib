Gupta_2022_brigatinib_iorr <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-EFFICACY model for confirmed ",
    "intracranial objective response by blinded independent review ",
    "committee in patients with intracranial central nervous system ",
    "metastases at baseline as a function of the brigatinib daily AUC ",
    "time-averaged over the last scan interval before the best ",
    "confirmed response, in adults with ALK-inhibitor-naive ",
    "ALK-positive advanced non-small cell lung cancer treated ",
    "first-line with brigatinib 180 mg orally once daily after a 7-day ",
    "lead-in at 90 mg (Gupta 2022, ALTA-1L brigatinib arm, n = 42). The ",
    "probability is expit(-0.913 + 0.1222 * AUC_BRIG_SCAN), with the ",
    "exposure in ug*h/mL/day entering linearly and uncentred. The ",
    "relationship is statistically significant (p = 0.049) and is the ",
    "paper's only significant efficacy result: the predicted ",
    "probability of intracranial response is 0.58, 0.83 and 0.99 at the ",
    "5th, 50th and 95th percentiles of steady-state exposure after 180 ",
    "mg once daily. Covariate analysis retained no covariate. The slope ",
    "is the printed odds ratio 1.13 per 1 ug*h/mL/day; the paper prints ",
    "no intercept, which was recovered from the vector-drawn fitted ",
    "curve of Figure 2c with the slope held at the printed value. There ",
    "is no PK layer and no ODE: the exposure metric is supplied as a ",
    "data column, derived in the source analysis from the individual ",
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
    dosing = "n/a (no dose events; exposure enters as the AUC_BRIG_SCAN covariate column)",
    concentration = "prob_icorr (probability of a BIRC-confirmed intracranial objective response, 0-1; also logit_icorr)"
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
        "In this model the exposure enters with slope log(1.13) per ug*h/mL/day."
      ),
      source_name = "Time-averaged AUC between last two scans (ug.h/mL/day)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 1L,
    n_observations = "42 binary best-confirmed-intracranial-response records, one per patient with baseline CNS metastases",
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
      "blinded independent review) and safety. The iORR analysis set ",
      "is the 42 of 47 patients with intracranial CNS metastases by ",
      "BIRC at baseline who had brigatinib concentration data; 31/42 ",
      "responded (Figure 2c tertile counts 8/14, 11/14, 12/14)."
    )
  )

  ini({
    # ==================================================================
    # Gupta 2022 Figure 2c: logit(p) = logit_ref + e_auc_logit *
    # AUC_BRIG_SCAN.
    #
    # Recovered from Figure 2c, whose fitted curve (annotated 'Exposure
    # effect P = 0.0486') is drawn as vector paths in the publisher PDF:
    # every Bezier node of the curve lies exactly on it, so the node
    # coordinates are read without pixel measurement.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -0.913; label("Logit of the probability of a BIRC-confirmed intracranial objective response at zero brigatinib exposure (unitless logit)") # figure-derived: Figure 2c fitted curve, intercept recovered with the slope held at the printed odds ratio; predicts 31.3 of 31 observed events

    # ----- Exposure slope on the logit -----
    e_auc_logit <- log(1.13); label("Log-odds of a BIRC-confirmed intracranial objective response per 1 ug*h/mL/day increase in AUC_BRIG_SCAN (unitless logit)") # odds ratio 1.13 (95% CI 1.01-1.28), p = 0.049 (Gupta 2022 Results, Exposure-efficacy analysis); Figure 2c panel annotates P = 0.0486

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_icorr <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Exposure enters linearly and uncentred, so logit_ref is the logit at
    # zero exposure, an extrapolation below the observed range.
    logit_icorr <- logit_ref + e_auc_logit * AUC_BRIG_SCAN

    prob_icorr <- expit(logit_icorr)

    prob_icorr ~ add(addSd_prob_icorr)
  })
}
