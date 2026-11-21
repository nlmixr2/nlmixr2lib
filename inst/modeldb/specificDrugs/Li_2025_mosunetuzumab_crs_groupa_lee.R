Li_2025_mosunetuzumab_crs_groupa_lee <- function() {
  description <- paste(
    "Linear logistic-regression exposure-safety model for the occurrence",
    "of at least one GRADE >= 2 CYTOKINE RELEASE SYNDROME event (Lee",
    "2014 grading) versus the maximum model-predicted CD20 receptor",
    "occupancy over the first two cycles (ROMAX, percent) in adults with",
    "relapsed/refractory B-cell NHL on FIXED-DOSE intravenous",
    "mosunetuzumab without step-up dosing (Li 2025, study GO29781 Group",
    "A, 0.2-2.8 mg q3w; Table S12). The slope is not statistically",
    "significant (p = 0.677). One of four companion fits (Group A or B x",
    "ASTCT or Lee 2014 grading)."
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
    concentration = "prob_crs_grade2 (probability of grade >= 2 cytokine release syndrome, Lee 2014 grading, 0-1)"
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
    n_subjects = 31,
    n_studies = 1,
    age_range = "30-84 years (median 64; Group A safety population, n = 32)",
    weight_range = "45.9-111 kg (median 76.2)",
    sex_female_pct = 40.6,
    race_ethnicity = c(White = 61.3, Asian = 32.3, Black = 3.23, Unknown = 3.23),
    disease_state = "Relapsed/refractory B-cell non-Hodgkin lymphoma (DLBCL 40.6, FL 28.1, transformed FL 15.6, MCL 6.25, MZL, RS and SLL 3.12 percent each)",
    dose_range = "Fixed dosing without step-up, 0.2-2.8 mg IV q3w (Group A, ROmax-evaluable patients)",
    regions = "Multinational",
    notes = "Group A (fixed-dose escalation) of the phase I/II GO29781 study (NCT02500407), ROmax-evaluable subset (31 of 32; Li 2025 Methods). Demographics from Table S11 (cohorts A2-A8). ROmax0-42 median 0.809 percent (range 0.032-3.43)."
  )

  ini({
    logite0 <- -2.14; label("Logit of the probability of grade >= 2 CRS at ROMAX = 0 percent, beta0 (unitless logit)") # Li 2025 Table S12, Intercept, Group A, Lee = -2.14 (SE 0.796)
    e_romax_logit <- 0.194; label("Log-odds of grade >= 2 CRS per 1 percentage-point increase in ROMAX (unitless logit)") # Li 2025 Table S12, Slope per 1% ROmax, Group A, Lee = 0.194 (SE 0.467; OR 1.21, 95% CI 0.436-3.04; p = 0.677)

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
