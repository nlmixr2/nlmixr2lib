Li_2025_mosunetuzumab_crs_stepup1 <- function() {
  description <- paste(
    "Linear logistic-regression exposure-safety model for GRADE >= 2",
    "CYTOKINE RELEASE SYNDROME (ASTCT grading) AFTER THE FIRST (CYCLE 1",
    "DAY 1, 0.4-1 mg) DOSE of intravenous mosunetuzumab versus the",
    "maximum model-predicted CD20 receptor occupancy within that dosing",
    "window (ROMAX, percent) in adults with relapsed/refractory B-cell",
    "NHL (Li 2025, study GO29781, N = 418; Figure 3a). One of three",
    "per-dosing-window fits whose curves shift right with each dose, the",
    "paper's evidence that step-up dosing desensitises CRS. Li 2025",
    "prints no coefficients: both were digitised from the Figure 3a",
    "fitted curve."
  )
  reference <- paste(
    "Li C-C, Bender B, Wilkins J, Li F, Turner DC, Wang B, Deng R,",
    "Vadhavkar S, Li Z, Kwan A, Huang H, Peng K, Penuel E, Huw L-Y, Chanu P,",
    "Li C, Yin S, Wei MC. A Novel Step-Up Dosage Regimen for Enhancing the",
    "Benefit-to-Risk Ratio of Mosunetuzumab in Relapsed or Refractory",
    "Follicular Lymphoma. Clin Pharmacol Ther. 2025;117(2):465-474.",
    "doi:10.1002/cpt.3445. Coefficients digitised from Figure 3a."
  )
  vignette <- "Li_2025_mosunetuzumab"
  units <- list(
    time = "n/a (static landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the ROMAX data column)",
    concentration = "prob_crs_stepup1_grade2 (probability of grade >= 2 cytokine release syndrome, ASTCT grading, 0-1)"
  )

  covariateData <- list(
    ROMAX = list(
      description = "Maximum model-predicted mosunetuzumab CD20 receptor occupancy within the Cycle 1 Day 1 dosing window",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maximum model-predicted CD20 receptor occupancy (percent,",
        "0-100) WITHIN THE CYCLE 1 DAY 1 DOSING WINDOW (from the first",
        "dose to the Day 8 dose); the Figure 3a curve spans 0.03-2.4",
        "percent (Li 2025 Figure 3 caption: 'model-predicted maximum",
        "CD20 receptor occupancy (RO%) for each dosing window'), not",
        "over days 0-42. RO% is Bender 2024's competitive",
        "equilibrium-binding occupancy against residual rituximab and",
        "obinutuzumab (the `RO` output of Bender_2024_mosunetuzumab).",
        "Enters linearly and uncentred, per 1 percentage point."
      ),
      source_name = "ROmax (first dose)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 418,
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
    disease_state = "Relapsed/refractory B-cell non-Hodgkin lymphoma (DLBCL 35.8, FL 37.1, transformed FL 13.0, MCL 8.88, Richter's syndrome 3.19 percent, other rare histologies)",
    dose_range = "Group A fixed dosing 0.2-2.8 mg IV q3w (n = 31) and Group B Cycle 1 step-up dosing 0.4/1.0/2.8 to 1/2/60 mg IV (n = 387); Day 1 dose 0.4-1 mg, Day 8 dose 1-2 mg, Day 15 dose 2.8-60 mg",
    regions = "Multinational",
    notes = "ROmax-evaluable patients of the phase I/II GO29781 study (NCT02500407; Li 2025 Methods and Figure 3 caption, N = 418). Demographics are those of the whole 439-patient safety population (Table S11)."
  )

  ini({
    # Li 2025 prints no coefficients for the per-dosing-window models. Both
    # were recovered from the fitted curve in Figure 3a (ASTCT grading):
    # the curve was pixel-traced on the calibrated log10(ROmax) x linear
    # probability axes and a logistic linear in untransformed ROmax was
    # fitted by least squares (RMSE < 0.001 in probability). The curve spans
    # ROmax 0.031-2.1 percent, probability 0.017 to 0.044.
    logite0 <- -4.03; label("Logit of the probability of grade >= 2 CRS at ROMAX = 0 percent, beta0 (unitless logit)") # digitised from Figure 3a (pixel-traced fitted curve, least-squares logistic fit; Li 2025 prints no coefficients for the per-dose-window models)
    e_romax_logit <- 0.461; label("Log-odds of grade >= 2 CRS per 1 percentage-point increase in ROMAX (unitless logit)") # digitised from Figure 3a (pixel-traced fitted curve, least-squares logistic fit; Li 2025 prints no coefficients for the per-dose-window models)

    # No between-subject variability and no residual error: the source is a
    # Bernoulli-likelihood logistic regression. The placeholder additive term
    # exists only because rxode2 requires an observation declaration.
    addSd_prob_crs_stepup1_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Li 2025 Methods, linear logistic regression:
    #   logit(p) = beta0 + beta_ROmax * ROmax
    # with ROmax in percent (0-100), untransformed and uncentred.
    logit_crs <- logite0 + e_romax_logit * ROMAX

    prob_crs_stepup1_grade2 <- expit(logit_crs)

    # Deterministic probability of the event. Downstream callers can sample
    # binary outcomes with rbinom(n, 1, prob_crs_stepup1_grade2) on the rxSolve output.
    prob_crs_stepup1_grade2 ~ add(addSd_prob_crs_stepup1_grade2)
  })
}
