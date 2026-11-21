Li_2025_patritumabDeruxtecan_nausea_vomiting <- function() {
  description <- paste(
    "MARGINAL (population-averaged) logistic exposure-safety model for grade",
    ">= 2 nausea or vomiting in adults with breast cancer, non-small-cell",
    "lung cancer or colorectal cancer given patritumab deruxtecan (HER3-DXd)",
    "in the J101, U102, U201 and U202 studies (Li 2025, n = 715). The source",
    "is a Bayesian logistic regression on the Cycle 1 anti-HER3-ac-DXd",
    "(conjugated payload) Cavg with no influential covariate, but Li 2025",
    "prints no coefficients for this endpoint. The intercept and the",
    "log-linear exposure slope here are back-solved from the",
    "population-simulated rates of Table 4 at 3.2, 4.8, 5.6 and 6.4 mg/kg",
    "Q3W, taking the exposure as proportional to dose about the Table S5",
    "median at 5.6 mg/kg (418 ng/mL). With no covariates in the source model",
    "the marginal and conditional forms differ only through the spread of",
    "exposure within a regimen. There is no PK layer. One of nine",
    "Li_2025_patritumabDeruxtecan exposure-response models."
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
    concentration = "prob_nausea_vomiting_grade2 (population-average probability of grade >= 2 nausea or vomiting, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cycle 1 average plasma concentration of anti-HER3-ac-DXd (antibody-conjugated DXd) over the dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Conjugated payload (anti-HER3-ac-DXd), not unconjugated DXd -- the conjugate fitted nausea/vomiting better than DXd (Li 2025 Results). Cycle 1 landmark (median onset day 4; Table S4). Centred at 418 ng/mL inside model(), the Table S5 median Cycle 1 anti-HER3-ac-DXd Cavg at 5.6 mg/kg Q3W. Back-solved against dose-proportional exposure, so the model is most faithful when supplied a regimen's typical (median) exposure.",
      source_name = "Cycle 1 anti-HER3-ac-DXd Cavg (ng/mL)"
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
    notes = "Exposure-safety population of 715 patients (Li 2025 Table S1); 234 (32.7%) had grade >= 2 nausea or vomiting (Table S4). The Table 4 predictions used to back-solve this model were simulated by the authors in 1,000 virtual patients resampled from the NSCLC patients."
  )

  ini({
    # Least-squares logit fit to the Table 4 population-simulated rates
    # 26 / 31 / 33 / 35 % at 3.2 / 4.8 / 5.6 / 6.4 mg/kg Q3W with exposure
    # 418 * dose / 5.6 ng/mL (residuals at most 0.004 logit). The
    # up-titration rate 26 % (same Cycle 1 dose as 3.2 mg/kg) is reproduced
    # exactly and is a hold-out check.
    logit_ref <- -0.7045; label("Logit of the population-average grade >= 2 nausea or vomiting probability at the reference anti-HER3-ac-DXd Cavg (unitless logit)") # back-solved from Table 4 nausea/vomiting row; 33.1% at 418 ng/mL vs printed 33 (28, 39) at 5.6 mg/kg
    e_cav_logit <- 0.6124; label("Change in the logit per unit natural log of Cycle 1 anti-HER3-ac-DXd Cavg (unitless logit)") # back-solved from Table 4 nausea/vomiting row 26/31/33/35 across 3.2-6.4 mg/kg

    addSd_prob_nausea_vomiting_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # 418 ng/mL: Table S5 median Cycle 1 anti-HER3-ac-DXd Cavg at 5.6 mg/kg Q3W.
    logit_nausea_vomiting_grade2 <- logit_ref + e_cav_logit * log(CAV / 418)
    prob_nausea_vomiting_grade2 <- expit(logit_nausea_vomiting_grade2)
    prob_nausea_vomiting_grade2 ~ add(addSd_prob_nausea_vomiting_grade2)
  })
}
