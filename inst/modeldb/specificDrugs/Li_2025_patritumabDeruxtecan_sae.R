Li_2025_patritumabDeruxtecan_sae <- function() {
  description <- paste(
    "MARGINAL (population-averaged) logistic exposure-safety model for any",
    "serious adverse event in adults with breast cancer, non-small-cell lung",
    "cancer or colorectal cancer given patritumab deruxtecan (HER3-DXd) in the",
    "J101, U102, U201 and U202 studies (Li 2025, n = 715). The source is a",
    "Bayesian logistic regression on the Cycle 1 DXd (unconjugated payload)",
    "Cavg with influential age and Japan country-of-origin effects, but Li",
    "2025 prints no coefficients and no forest plot for this endpoint. The",
    "intercept and the log-linear exposure slope here are back-solved from",
    "the population-simulated SAE rates of Table 4 at 3.2, 4.8, 5.6 and",
    "6.4 mg/kg Q3W, taking the exposure as proportional to dose about the",
    "Table S5 median at 5.6 mg/kg (2.1 ng/mL). The covariate effects are",
    "therefore folded into the intercept: the model returns a",
    "population-average probability at a given exposure, not the",
    "conditional probability for a patient with given covariates. There is",
    "no PK layer. One of nine Li_2025_patritumabDeruxtecan exposure-response",
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
    time = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_sae (population-average probability of any serious adverse event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cycle 1 average plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, over the dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd). Cycle 1 landmark (median SAE onset day 51.5, within the first 3 cycles; Li 2025 Methods and Table S4). Centred at 2.1 ng/mL inside model(), the Table S5 median Cycle 1 DXd Cavg at 5.6 mg/kg Q3W. Because the slope was back-solved against dose-proportional exposure, the model is most faithful when supplied a regimen's typical (median) exposure; supplying individual exposures does not recover the source's conditional model.",
      source_name = "Cycle 1 DXd Cavg (ng/mL)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Influential covariate of the source SAE model (Li 2025 Table 2), but no coefficient, forest plot or covariate-stratified prediction is printed, so it cannot be recovered; it is averaged into the marginal intercept."
    ),
    REGION_JAPAN = list(
      description = "Country of origin Japan indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Influential covariate of the source SAE model (Li 2025 Table 2, 'Country of origin: Japan'), not recoverable for the same reason as AGE; averaged into the marginal intercept."
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
    notes = "Exposure-safety population of 715 patients (Li 2025 Table S1); 272 (38.0%) had an SAE (Table S4). The Table 4 predictions used to back-solve this model were simulated by the authors in 1,000 virtual patients resampled from the NSCLC patients."
  )

  ini({
    # Li 2025 prints neither coefficients nor a forest plot for SAE. The two
    # values below are a least-squares logit fit to the Table 4
    # population-simulated rates 28 / 38 / 42 / 46 % at 3.2 / 4.8 / 5.6 /
    # 6.4 mg/kg Q3W with exposure 2.1 * dose / 5.6 ng/mL (residuals at most
    # 0.007 logit). The up-titration rate 28 % (same Cycle 1 dose as
    # 3.2 mg/kg) is reproduced exactly and is a hold-out check.
    logit_ref <- -0.3161; label("Logit of the population-average SAE probability at the reference DXd Cavg (unitless logit)") # back-solved from Table 4 SAE row; 42.2% at 2.1 ng/mL vs printed 42 (37, 47) at 5.6 mg/kg
    e_cav_logit <- 1.1253; label("Change in the logit per unit natural log of Cycle 1 DXd Cavg (unitless logit)") # back-solved from Table 4 SAE row 28/38/42/46 across 3.2-6.4 mg/kg

    addSd_prob_sae <- fixed(0.001); label("Placeholder additive residual SD on the event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # 2.1 ng/mL: Table S5 median Cycle 1 DXd Cavg at 5.6 mg/kg Q3W.
    logit_sae <- logit_ref + e_cav_logit * log(CAV / 2.1)
    prob_sae <- expit(logit_sae)
    prob_sae ~ add(addSd_prob_sae)
  })
}
