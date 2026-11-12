Li_2025_patritumabDeruxtecan_discontinuation <- function() {
  description <- paste(
    "MARGINAL discrete-time (per 21-day cycle) logit-hazard exposure-safety",
    "model for an adverse event leading to treatment discontinuation in",
    "adults with breast cancer, non-small-cell lung cancer or colorectal",
    "cancer given patritumab deruxtecan (HER3-DXd) in the J101, U102, U201",
    "and U202 studies (Li 2025, n = 715). The source is a Bayesian",
    "discrete-time survival model in which the probability of a first event",
    "in cycle j, given none earlier, is expit(alpha0 + alpha_j + f(C_ij)),",
    "with cycle random effects alpha_j and the cycle-varying DXd",
    "(unconjugated payload) Cavg; no covariate was influential. Li 2025",
    "prints no coefficients and no estimate of the cycle-effect variance.",
    "The baseline per-cycle logit hazard and the log-linear exposure slope",
    "here are back-solved from the population-simulated 6-month rates of",
    "Table 4 at 3.2, 4.8, 5.6 and 6.4 mg/kg Q3W, taking the exposure as",
    "proportional to dose about the Table S5 median at 5.6 mg/kg (2.1",
    "ng/mL) and every cycle effect as zero. The per-cycle hazard h is",
    "integrated as a continuous rate -log(1 - h) / 21 per day into cumhaz,",
    "so survival exp(-cumhaz) equals the discrete-time product over whole",
    "cycles. There is no PK layer. One of nine Li_2025_patritumabDeruxtecan",
    "exposure-response models."
  )
  reference <- paste(
    "Li L, Lee M, Joshi R, Xu Y, Garcia R, Baron K, Anderson E, Waterhouse T,",
    "Sternberg D, Shrestha P, Vaddady P, Abutarif M, Garimella T.",
    "Exposure-Response Relationships in Patients with Non-Small-Cell Lung",
    "Cancer and Other Solid Tumors Treated with Patritumab Deruxtecan",
    "(HER3-DXd). Clin Pharmacol Ther. 2025;118(2):408-417.",
    "doi:10.1002/cpt.3674."
  )
  vignette <- "Li_2025_patritumabDeruxtecan_exposure_response"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; exposure enters as the time-varying CAV covariate column)",
    concentration = "prob_treatment_discontinuation (cumulative probability of a first AE leading to treatment discontinuation by time t, 0-1)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "treatment-discontinuation hazard",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CAV = list(
      description = "Average plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, over the current dosing cycle; TIME-VARYING, one value per cycle.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd). Supply one record per cycle carrying that cycle's Cavg (Li 2025 Supplementary methods). Centred at 2.1 ng/mL inside model(), the Table S5 median Cycle 1 DXd Cavg at 5.6 mg/kg Q3W. Back-solved against dose-proportional exposure, so the model is most faithful when supplied a regimen's typical (median) exposure per cycle.",
      source_name = "DXd Cavg (ng/mL)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 715L,
    n_studies = 4L,
    age_range = "29-83 years (median 61)",
    weight_range = "32.4-137 kg (median 60.0)",
    sex_female_pct = 69.2,
    race_ethnicity = c(Asian = 54.7, White = 37.5, Other = 5.7, Black = 2.1),
    disease_state = "HER3-expressing advanced breast cancer (J101), metastatic or unresectable NSCLC (U102, U201; EGFR-mutated in 446 of 493) and advanced colorectal cancer (U202).",
    dose_range = "Patritumab deruxtecan 1.6-8 mg/kg IV Q3W, up-titration (3.2 / 4.8 / 6.4 mg/kg) or 4.2 mg/kg Q2W for Cycles 1-3 then 6.4 mg/kg Q3W.",
    regions = "Global; country of origin Japan 34.4%, other Asian 13.3%, non-Asian 52.3%.",
    notes = "Exposure-safety population of 715 patients (Li 2025 Table S1); 69 (9.7%) had an AE leading to discontinuation, median time to first event 94 days (Table S4). The Table 4 predictions used to back-solve this model were simulated by the authors in 1,000 virtual patients resampled from the NSCLC patients."
  )

  ini({
    # Least-squares fit of 1 - (1 - h)^(182.625 / 21) to the Table 4
    # population-simulated 6-month rates 4.7 / 7.4 / 9.0 / 10 % at
    # 3.2 / 4.8 / 5.6 / 6.4 mg/kg Q3W, with h = expit(logit_ref +
    # e_cav_logit * log(dose / 5.6)) (fitted 4.8 / 7.4 / 8.8 / 10.1).
    # Hold-out: the up-titration regimen predicts 9.2 % against the printed
    # 9.3 (6.4, 13).
    logit_ref <- -4.5444; label("Logit of the per-cycle treatment-discontinuation hazard at the reference DXd Cavg (unitless logit)") # back-solved from Table 4 discontinuation row; h = 0.0105 per cycle at 2.1 ng/mL
    e_cav_logit <- 1.1383; label("Change in the per-cycle hazard logit per unit natural log of the cycle DXd Cavg (unitless logit)") # back-solved from Table 4 discontinuation row 4.7/7.4/9.0/10 across 3.2-6.4 mg/kg

    addSd_prob_treatment_discontinuation <- fixed(0.001); label("Placeholder additive residual SD on the cumulative event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # 2.1 ng/mL: Table S5 median Cycle 1 DXd Cavg at 5.6 mg/kg Q3W.
    logit_haz_cycle <- logit_ref + e_cav_logit * log(CAV / 2.1)
    haz_cycle <- expit(logit_haz_cycle)
    # Discrete per-cycle hazard spread evenly over the 21-day Q3W cycle, so
    # each whole cycle multiplies survival by exactly (1 - haz_cycle).
    d/dt(cumhaz) <- -log(1 - haz_cycle) / 21
    sur <- exp(-cumhaz)
    prob_treatment_discontinuation <- 1 - sur
    prob_treatment_discontinuation ~ add(addSd_prob_treatment_discontinuation)
  })
}
