Li_2025_mosunetuzumab_crs_groupb_astct <- function() {
  description <- paste(
    "Linear logistic-regression exposure-safety model for the occurrence",
    "of at least one GRADE >= 2 CYTOKINE RELEASE SYNDROME event (ASTCT",
    "grading) versus the maximum model-predicted CD20 receptor occupancy",
    "over the first two cycles (ROMAX, percent) in adults with",
    "relapsed/refractory B-cell NHL on intravenous mosunetuzumab with",
    "Cycle 1 STEP-UP dosing (Li 2025, study GO29781 Group B, 0.4/1/2.8",
    "to 1/2/60/30 mg; Table S12). This is the base CRS model behind the",
    "paper's dose-regimen predictions (Figure 4: median 11.8 percent at",
    "1/2/60/30 mg). One of four companion fits (Group A or B x ASTCT or",
    "Lee 2014 grading)."
  )
  reference <- paste(
    "Li C-C, Bender B, Wilkins J, Li F, Turner DC, Wang B, Deng R,",
    "Vadhavkar S, Li Z, Kwan A, Huang H, Peng K, Penuel E, Huw L-Y, Chanu P,",
    "Li C, Yin S, Wei MC. A Novel Step-Up Dosage Regimen for Enhancing the",
    "Benefit-to-Risk Ratio of Mosunetuzumab in Relapsed or Refractory",
    "Follicular Lymphoma. Clin Pharmacol Ther. 2025;117(2):465-474.",
    "doi:10.1002/cpt.3445."
  )
  vignette <- "Li_2025_mosunetuzumab"
  units <- list(
    time = "n/a (static landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the ROMAX data column)",
    concentration = "prob_crs_grade2 (probability of grade >= 2 cytokine release syndrome, ASTCT grading, 0-1)"
  )

  covariateData <- list(
    ROMAX = list(
      description = "Maximum model-predicted mosunetuzumab CD20 receptor occupancy over days 0-42",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maximum model-predicted CD20 receptor occupancy (percent,",
        "0-100) over days 0-42 (the first two 21-day cycles), Li 2025's",
        "primary CRS exposure metric ('model-estimated maximum CD20 RO",
        "(ROmax%) within the first two cycles (42 days)'). RO% is Bender",
        "2024's competitive equilibrium-binding occupancy of",
        "mosunetuzumab against residual rituximab and obinutuzumab: the",
        "`RO` output of Bender_2024_mosunetuzumab, maximised over days",
        "0-42 of a solve under the patient's regimen. Table S11: median",
        "0.809 percent in Group A, 7.5 in B1-B10, 36 in B11 (range",
        "0.0294-64.3). Enters linearly and uncentred, per 1 percentage",
        "point."
      ),
      source_name = "ROmax0-42"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 387,
    n_studies = 1,
    age_range = "19-96 years (median 63; whole safety population, Table S11)",
    weight_range = "37.1-163 kg (median 77.9)",
    sex_female_pct = 35.3,
    race_ethnicity = c(
      White = 76.2,
      Asian = 17.6,
      Black = 2.75,
      `American Indian or Alaska Native` = 0.458,
      Unknown = 2.97
    ),
    disease_state = "Relapsed/refractory B-cell non-Hodgkin lymphoma (FL 148, DLBCL 142, transformed FL 46, other histologies 51: MCL 34, Richter's syndrome 10, MZL 2, SLL, FL3B, DLBCL/MCL, melanoma and transformed MZL 1 each)",
    dose_range = "Cycle 1 step-up dosing from 0.4/1.0/2.8 to 1/2/60 and 1/2/60/30 mg IV q3w (Group B)",
    regions = "Multinational",
    notes = "Group B (step-up dosing) of the phase I/II GO29781 study (NCT02500407), ROmax-evaluable subset (387 patients; Li 2025 Methods); safety clinical cutoff 27 August 2021. Demographics are those of the whole 439-patient safety population (Table S11; Group B is 407 of them). ROmax0-42 median 7.5 percent in B1-B10 and 36 percent in B11 (1/2/60/30 mg). Observed ASTCT grade >= 2 CRS: B7 (1/2/13.5 mg) 1/73, B11 33/198."
  )

  ini({
    logite0 <- -2.64; label("Logit of the probability of grade >= 2 CRS at ROMAX = 0 percent, beta0 (unitless logit)") # Li 2025 Table S12, Intercept, Group B, ASTCT = -2.64 (SE 0.279)
    e_romax_logit <- 0.0196; label("Log-odds of grade >= 2 CRS per 1 percentage-point increase in ROMAX (unitless logit)") # Li 2025 Table S12, Slope per 1% ROmax, Group B, ASTCT = 0.0196 (SE 0.00804; OR 1.02, 95% CI 1.00-1.04; p = 0.0150)

    # No between-subject variability and no residual error: the source is a
    # Bernoulli-likelihood logistic regression. The placeholder additive term
    # exists only because rxode2 requires an observation declaration.
    addSd_prob_crs_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Li 2025 Methods, linear logistic regression:
    #   logit(p) = beta0 + beta_ROmax * ROmax
    # with ROmax in percent (0-100), untransformed and uncentred.
    logit_crs <- logite0 + e_romax_logit * ROMAX

    prob_crs_grade2 <- expit(logit_crs)

    # Deterministic probability of the event. Downstream callers can sample
    # binary outcomes with rbinom(n, 1, prob_crs_grade2) on the rxSolve output.
    prob_crs_grade2 ~ add(addSd_prob_crs_grade2)
  })
}
