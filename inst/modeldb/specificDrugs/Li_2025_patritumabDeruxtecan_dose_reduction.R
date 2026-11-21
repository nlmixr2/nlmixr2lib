Li_2025_patritumabDeruxtecan_dose_reduction <- function() {
  description <- paste(
    "MARGINAL discrete-time (per 21-day cycle) logit-hazard exposure-safety",
    "model for an adverse event leading to dose reduction in adults with",
    "breast cancer, non-small-cell lung cancer or colorectal cancer given",
    "patritumab deruxtecan (HER3-DXd) in the J101, U102, U201 and U202",
    "studies (Li 2025, n = 715). The source is a Bayesian discrete-time",
    "survival model in which the probability of a first event in cycle j,",
    "given none earlier, is expit(alpha0 + alpha_j + f(C_ij) + x'gamma),",
    "with cycle random effects alpha_j, the cycle-varying DXd (unconjugated",
    "payload) Cmax and an influential age effect. Li 2025 prints no",
    "coefficients and no estimate of the cycle-effect variance. The baseline",
    "per-cycle logit hazard and the log-linear exposure slope here are",
    "back-solved from the population-simulated 6-month rates of Table 4 at",
    "3.2, 4.8, 5.6 and 6.4 mg/kg Q3W, taking the exposure as proportional to",
    "dose about the Table S5 median at 5.6 mg/kg (24.2 ng/mL) and every",
    "cycle effect as zero; age is folded into the intercept. This back-solve",
    "FAILS its hold-out: it predicts 24% at 6 months for the up-titration",
    "regimen against the printed 20%, because dose reductions are",
    "front-loaded (median onset day 23) and the unreported cycle effects",
    "matter. Use with that caveat. The per-cycle hazard h is integrated as a",
    "continuous rate -log(1 - h) / 21 per day into cumhaz. There is no PK",
    "layer. One of nine Li_2025_patritumabDeruxtecan exposure-response",
    "models."
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
    dosing = "n/a (no dose events; exposure enters as the time-varying CMAX covariate column)",
    concentration = "prob_dose_reduction (cumulative probability of a first AE leading to dose reduction by time t, 0-1)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "dose-reduction hazard",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CMAX = list(
      description = "Maximum plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, in the current dosing cycle; TIME-VARYING, one value per cycle.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd). Supply one record per cycle carrying that cycle's Cmax (Li 2025 Supplementary methods). Centred at 24.2 ng/mL inside model(), the Table S5 median Cycle 1 DXd Cmax at 5.6 mg/kg Q3W. Back-solved against dose-proportional exposure, so the model is most faithful when supplied a regimen's typical (median) exposure per cycle.",
      source_name = "DXd Cmax (ng/mL)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Influential covariate of the source dose-reduction model (Li 2025 Table 2), but no coefficient, forest plot or covariate-stratified prediction is printed, so it cannot be recovered; it is averaged into the marginal intercept."
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
    notes = "Exposure-safety population of 715 patients (Li 2025 Table S1); 163 (22.8%) had an AE leading to dose reduction, median time to first event 23 days (Table S4). The Table 4 predictions used to back-solve this model were simulated by the authors in 1,000 virtual patients resampled from the NSCLC patients."
  )

  ini({
    # Least-squares fit of 1 - (1 - h)^(182.625 / 21) to the Table 4
    # population-simulated 6-month rates 13 / 20 / 23 / 26 % at
    # 3.2 / 4.8 / 5.6 / 6.4 mg/kg Q3W, with h = expit(logit_ref +
    # e_cmax_logit * log(dose / 5.6)) (fitted 13.1 / 19.8 / 23.0 / 26.1).
    # Hold-out FAILS: the up-titration regimen predicts 24.0 % against the
    # printed 20 (16, 24).
    logit_ref <- -3.4901; label("Logit of the per-cycle dose-reduction hazard at the reference DXd Cmax (unitless logit)") # back-solved from Table 4 dose-reduction row; h = 0.0296 per cycle at 24.2 ng/mL
    e_cmax_logit <- 1.1172; label("Change in the per-cycle hazard logit per unit natural log of the cycle DXd Cmax (unitless logit)") # back-solved from Table 4 dose-reduction row 13/20/23/26 across 3.2-6.4 mg/kg

    addSd_prob_dose_reduction <- fixed(0.001); label("Placeholder additive residual SD on the cumulative event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # 24.2 ng/mL: Table S5 median Cycle 1 DXd Cmax at 5.6 mg/kg Q3W.
    logit_haz_cycle <- logit_ref + e_cmax_logit * log(CMAX / 24.2)
    haz_cycle <- expit(logit_haz_cycle)
    # Discrete per-cycle hazard spread evenly over the 21-day Q3W cycle, so
    # each whole cycle multiplies survival by exactly (1 - haz_cycle).
    d/dt(cumhaz) <- -log(1 - haz_cycle) / 21
    sur <- exp(-cumhaz)
    prob_dose_reduction <- 1 - sur
    prob_dose_reduction ~ add(addSd_prob_dose_reduction)
  })
}
