Yarur_2025_ustekinumab_ehi50_cl <- function() {
  description <- paste0(
    "Logistic exposure-response model for an ENDOSCOPIC HEALING INDEX ",
    "(EHI) ABOVE 50 POINTS -- a serum-biomarker surrogate of ",
    "endoscopically active disease -- versus ustekinumab CLEARANCE and ",
    "prior biologic exposure in patients with inflammatory bowel disease ",
    "on ustekinumab maintenance (Yarur 2025, 83 patients, 303 dose ",
    "cycles with an EHI). Clearance is the per-cycle empirical-Bayes ",
    "estimate from modellib('Yarur_2025_ustekinumab'), supplied as the ",
    "CL_INDIV column. Higher clearance raises the probability of active ",
    "disease (+18.0 log-odds per L/day, P < 0.001); prior biologic ",
    "exposure is not significant (P = 0.659). Typical-value model: the ",
    "between-subject variance of the mixed-effects logistic fit was not ",
    "reported. Reproduces Figure 3B. One of five ",
    "Yarur_2025_ustekinumab_* exposure-response models."
  )
  reference <- paste(
    "Yarur AJ, Dervieux T, Ungaro R, Spencer EA, Bruss A, Nunez L,",
    "Berens B, Vermeire S, Wang Z, Panetta JC, Dreesen E, Dubinsky MC.",
    "Ustekinumab Drug Clearance Is Better Associated with Disease Control",
    "than Serum Trough Concentrations in a Prospective Cohort of",
    "Inflammatory Bowel Disease. Pharmaceutics. 2025;17(2):187.",
    "doi:10.3390/pharmaceutics17020187.",
    "Coefficients from Supplementary Table S1, row 'EHI>50 with prior",
    "biologics as covariate and UST CL as regressor'.",
    "Clearance supplied by modellib('Yarur_2025_ustekinumab').",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; clearance enters as the CL_INDIV covariate column)",
    concentration = "prob_ehi50 (probability of an endoscopic healing index above 50 points, 0-1)"
  )

  covariateData <- list(
    CL_INDIV = list(
      description = "Individual empirical-Bayes ustekinumab clearance for the dose cycle",
      units = "L/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters the logit linearly and UNCENTERED, so logite0 is the ",
        "logit at zero clearance. Estimated per dose cycle by Bayesian ",
        "estimation (Monolix) from the trough concentration, albumin and ",
        "ATU status, with every parameter of ",
        "modellib('Yarur_2025_ustekinumab') fixed as the prior (Methods ",
        "2.3). For a typical-value simulation use ",
        "0.161 x (ALB_gdL / 4.0)^-0.887. Study median 0.157 L/day, IQR ",
        "0.120-0.190 (Table 4); Figure 3B spans 0.05-0.5 L/day."
      ),
      source_name = "UST CL"
    ),
    PRIOR_BIO = list(
      description = "Prior exposure to another biologic before starting ustekinumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (biologic-naive; ustekinumab was the first biologic)",
      notes = "89% (74/83) were biologic-experienced (Table 4). Figure 3B labels the two levels 'naive' and 'experienced'.",
      source_name = "prior biologics"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 1L,
    n_observations = "303 dose cycles with an EHI score; 54 (17.8%) above 50 points",
    age_median = "42 years (IQR 28-57)",
    sex_female_pct = 52,
    disease_state = "Crohn's disease (78, 94%) or ulcerative colitis (5) on ustekinumab maintenance for at least 16 weeks",
    dose_range = "Weight-based IV loading dose, then 90 mg SC every 4 weeks (77 cycles) or every 8 weeks",
    regions = "United States (Medical College of Wisconsin; Mount Sinai, New York)",
    notes = "Prospective two-centre cohort; baseline characteristics in Yarur 2025 Table 4. Median EHI 33 points (IQR 23-43)."
  )

  ini({
    # Supplementary Table S1: logit(P) = theta_pop + cov1 * theta_cov1 +
    # cov2 * theta_cov2, estimates printed as estimate +/- SE.
    logite0 <- -6.29; label("Logit of EHI above 50 points for a biologic-naive patient at zero clearance (unitless logit)")  # Table S1 'theta_pop = -6.29 +/- 1.45'
    e_priorbio_ehi50 <- 0.639; label("Log-odds shift of EHI above 50 points for biologic-experienced vs naive patients (unitless logit)")  # Table S1 'theta_prior biologics = +0.639 +/- 1.45 (p = 0.659)'
    e_cl_ehi50 <- 18.0; label("Log-odds of EHI above 50 points per L/day of ustekinumab clearance (day/L)")  # Table S1 'theta_CL = 18.0 +/- 3.2 (p < 0.001)'

    addSd_prob_ehi50 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_ehi50 <- logite0 + e_priorbio_ehi50 * PRIOR_BIO + e_cl_ehi50 * CL_INDIV
    prob_ehi50 <- expit(logit_ehi50)
    prob_ehi50 ~ add(addSd_prob_ehi50)
  })
}
