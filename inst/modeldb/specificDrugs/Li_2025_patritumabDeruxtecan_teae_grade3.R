Li_2025_patritumabDeruxtecan_teae_grade3 <- function() {
  description <- paste(
    "Bayesian logistic-regression exposure-safety model for any grade >= 3",
    "treatment-emergent adverse event in adults with breast cancer,",
    "non-small-cell lung cancer or colorectal cancer given patritumab",
    "deruxtecan (HER3-DXd) in the J101, U102, U201 and U202 studies (Li 2025,",
    "n = 715; 1.6-8 mg/kg IV Q3W, up-titration or an alternative Q2W/Q3W",
    "regimen). logit(P) is log-linear in the Cycle 1 DXd (unconjugated",
    "payload) Cavg and linear in age, with additive shifts for race (Asian,",
    "Black, Other vs White). Li 2025 prints no coefficient table; every value",
    "is back-solved from the typical-value forest plot (Figure S4).",
    "Covariate-by-exposure interaction terms of the source model are not",
    "recoverable and are omitted. There is no PK layer: exposure is a data",
    "column from the companion population PK model (see",
    "Lee_2023_patritumab). One of nine Li_2025_patritumabDeruxtecan",
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
    time = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_teae_grade3 (probability of any grade >= 3 TEAE, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cycle 1 average plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, over the dosing interval, per subject.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd), not the conjugate. Cycle 1 landmark; median onset of the event was day 15 (Li 2025 Table S4). Empirical-Bayes predictions from the population PK model 'previously developed' (uncited in Li 2025; the same group's model is Lee_2023_patritumab). Centred at 2.22 ng/mL inside model(), the reference exposure implied by Figure S4; Table S5 reports a median of 2.1 ng/mL at 5.6 mg/kg Q3W. Figure S4 5th/25th/75th/95th percentiles: 0.6 / 1.4 / 3.2 / 5.2 ng/mL.",
      source_name = "Cycle 1 DXd Cavg (ng/mL)"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Influential but not significant (Li 2025 Table 2). Centred at the 61-year safety-population median inside model() (Results; the Figure S4 rows imply 61.6 years). Figure S4 5th/25th/75th/95th percentiles: 41 / 53 / 68 / 77 years.",
      source_name = "Baseline age (y)"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_BLACK and RACE_OTHER are also 0).",
      notes = "Influential but not significant. Figure S4 'Race: Asian' 73.5% vs reference 67.7%.",
      source_name = "Race: Asian"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator (1 = Black, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_OTHER are also 0).",
      notes = "Figure S4 'Race: Black' 67.6% vs reference 67.7%.",
      source_name = "Race: Black"
    ),
    RACE_OTHER = list(
      description = "Other race indicator (1 = race other than White, Asian or Black, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_BLACK are also 0).",
      notes = "Figure S4 'Race: Other' 66.3% vs reference 67.7%.",
      source_name = "Race: Other"
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
    notes = "Exposure-safety population of 715 patients (J101 n = 182, U102 n = 219, U201 n = 274, U202 n = 40; Li 2025 Table S1); 482 (67.4%) had a grade >= 3 TEAE (Table S4)."
  )

  ini({
    # Li 2025 prints no coefficient table. Values below are reconstructed
    # from the typical-value forest plot of Figure S4 under the
    # logit-log-linear exposure form of Supplementary methods (residuals at
    # most 0.018 logit across the four exposure rows; the linear form leaves
    # 0.37 logit). Age is linear on the logit (residuals at most 0.012
    # logit). Race shifts are read directly from the printed probabilities.
    logit_ref <- log(0.677 / (1 - 0.677)); label("Logit of grade >= 3 TEAE probability for the reference patient (unitless logit)") # Figure S4 'Reference patient' 67.7 [58.6, 74.8]
    e_cav_logit <- 1.2608; label("Change in the logit per unit natural log of Cycle 1 DXd Cavg (unitless logit)") # back-solved from Figure S4 exposure rows 0.6/1.4/3.2/5.2 ng/mL -> 28.4/54.4/76.9/85.8 (reference 2.22 ng/mL implied)
    e_age_logit <- 0.01192; label("Change in the logit per year of age (unitless logit)") # back-solved from Figure S4 rows 41/53/68/77 y -> 62.2/65.3/69.5/71.5 (reference 61 y)
    e_race_asian_logit <- log(0.735 / 0.265) - log(0.677 / 0.323); label("Shift in the logit for Asian vs White race (unitless logit)") # Figure S4 'Race: Asian' 73.5 [63.7, 81.9]
    e_race_black_logit <- log(0.676 / 0.324) - log(0.677 / 0.323); label("Shift in the logit for Black vs White race (unitless logit)") # Figure S4 'Race: Black' 67.6 [55.3, 78.5]
    e_race_other_logit <- log(0.663 / 0.337) - log(0.677 / 0.323); label("Shift in the logit for Other vs White race (unitless logit)") # Figure S4 'Race: Other' 66.3 [52.1, 74.8]

    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Reference values are those implied by the Figure S4 rows (back-solved).
    logit_teae_grade3 <- logit_ref +
      e_cav_logit * log(CAV / 2.22) +
      e_age_logit * (AGE - 61) +
      e_race_asian_logit * RACE_ASIAN +
      e_race_black_logit * RACE_BLACK +
      e_race_other_logit * RACE_OTHER
    prob_teae_grade3 <- expit(logit_teae_grade3)
    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
