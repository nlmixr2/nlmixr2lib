Yarur_2025_ustekinumab_clinbiorem_cl <- function() {
  description <- paste0(
    "Logistic exposure-response model for combined CLINICAL AND ",
    "BIOCHEMICAL REMISSION versus ustekinumab CLEARANCE and prior ",
    "biologic exposure in patients with inflammatory bowel disease on ",
    "ustekinumab maintenance (Yarur 2025, 83 patients, 312 dose cycles). ",
    "Remission is Harvey-Bradshaw Index < 5 (Crohn's disease) or partial ",
    "Mayo score < 2 (ulcerative colitis) together with CRP < 3 mg/L. ",
    "Clearance is the per-cycle empirical-Bayes estimate from ",
    "modellib('Yarur_2025_ustekinumab'), supplied as the CL_INDIV column. ",
    "Higher clearance lowers the probability of remission (-13.5 log-odds ",
    "per L/day, P < 0.001). Typical-value model: the between-subject ",
    "variance of the mixed-effects logistic fit was not reported. ",
    "Reproduces Figure 2B. One of five Yarur_2025_ustekinumab_* ",
    "exposure-response models."
  )
  reference <- paste(
    "Yarur AJ, Dervieux T, Ungaro R, Spencer EA, Bruss A, Nunez L,",
    "Berens B, Vermeire S, Wang Z, Panetta JC, Dreesen E, Dubinsky MC.",
    "Ustekinumab Drug Clearance Is Better Associated with Disease Control",
    "than Serum Trough Concentrations in a Prospective Cohort of",
    "Inflammatory Bowel Disease. Pharmaceutics. 2025;17(2):187.",
    "doi:10.3390/pharmaceutics17020187.",
    "Coefficients from Supplementary Table S1, row 'Clinical & biochemical",
    "remission with prior biologics as covariate and UST CL as regressor'.",
    "Clearance supplied by modellib('Yarur_2025_ustekinumab').",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; clearance enters as the CL_INDIV covariate column)",
    concentration = "prob_clinbiorem (probability of clinical and biochemical remission, 0-1)"
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
        "MAP/MCMC estimation (Monolix) from the trough concentration, ",
        "albumin and ATU status, with every parameter of ",
        "modellib('Yarur_2025_ustekinumab') fixed as the prior (Methods ",
        "2.3). For a typical-value simulation use ",
        "0.161 x (ALB_gdL / 4.0)^-0.887. Study median 0.157 L/day, IQR ",
        "0.120-0.190 (Table 4); Figure 2B spans 0.05-0.5 L/day."
      ),
      source_name = "UST CL"
    ),
    PRIOR_BIO = list(
      description = "Prior exposure to another biologic before starting ustekinumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (biologic-naive; ustekinumab was the first biologic)",
      notes = "89% (74/83) were biologic-experienced (Table 4). Figure 2B labels the two levels 'naive' and 'experienced'.",
      source_name = "prior biologics"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 1L,
    n_observations = "312 dose cycles (median 3 per patient, range 1-6); 84 (26%) in clinical and biochemical remission",
    age_median = "42 years (IQR 28-57)",
    sex_female_pct = 52,
    disease_state = "Crohn's disease (78, 94%) or ulcerative colitis (5) on ustekinumab maintenance for at least 16 weeks",
    dose_range = "Weight-based IV loading dose, then 90 mg SC every 4 weeks (77 cycles) or every 8 weeks",
    regions = "United States (Medical College of Wisconsin; Mount Sinai, New York)",
    notes = "Prospective two-centre cohort; baseline characteristics in Yarur 2025 Table 4. Median follow-up 279 days (IQR 210-351)."
  )

  ini({
    # Supplementary Table S1: logit(P) = theta_pop + cov1 * theta_cov1 +
    # cov2 * theta_cov2, estimates printed as estimate +/- SE.
    logite0 <- 2.29; label("Logit of clinical and biochemical remission for a biologic-naive patient at zero clearance (unitless logit)")  # Table S1 'theta_pop = +2.29 +/- 0.88'
    e_priorbio_clinbiorem <- -2.02; label("Log-odds shift of remission for biologic-experienced vs naive patients (unitless logit)")  # Table S1 'theta_prior biologics = -2.02 +/- 0.81 (p = 0.013)'
    e_cl_clinbiorem <- -13.5; label("Log-odds of remission per L/day of ustekinumab clearance (day/L)")  # Table S1 'theta_CL = -13.5 +/- 3.7 (p < 0.001)'

    addSd_prob_clinbiorem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_clinbiorem <- logite0 + e_priorbio_clinbiorem * PRIOR_BIO + e_cl_clinbiorem * CL_INDIV
    prob_clinbiorem <- expit(logit_clinbiorem)
    prob_clinbiorem ~ add(addSd_prob_clinbiorem)
  })
}
