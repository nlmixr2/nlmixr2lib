Yarur_2025_ustekinumab_clinbiorem_time <- function() {
  description <- paste0(
    "Logistic model for combined CLINICAL AND BIOCHEMICAL REMISSION ",
    "versus follow-up TIME and prior biologic exposure in patients with ",
    "inflammatory bowel disease on ustekinumab maintenance (Yarur 2025, ",
    "83 patients, 312 dose cycles). Remission is Harvey-Bradshaw ",
    "Index < 5 (Crohn's disease) or partial Mayo score < 2 (ulcerative ",
    "colitis) together with CRP < 3 mg/L. No exposure term: prior ",
    "biologic exposure lowers the probability of remission (P = 0.004) ",
    "while time has little effect (P = 0.158). Time is the model time in ",
    "days from study baseline. Typical-value model: the between-subject ",
    "variance of the mixed-effects logistic fit was not reported. ",
    "Reproduces the 66% versus 14% remission probabilities at 279 days ",
    "quoted in Results 3.3. One of five Yarur_2025_ustekinumab_* ",
    "exposure-response models."
  )
  reference <- paste(
    "Yarur AJ, Dervieux T, Ungaro R, Spencer EA, Bruss A, Nunez L,",
    "Berens B, Vermeire S, Wang Z, Panetta JC, Dreesen E, Dubinsky MC.",
    "Ustekinumab Drug Clearance Is Better Associated with Disease Control",
    "than Serum Trough Concentrations in a Prospective Cohort of",
    "Inflammatory Bowel Disease. Pharmaceutics. 2025;17(2):187.",
    "doi:10.3390/pharmaceutics17020187.",
    "Coefficients from Supplementary Table S1, row 'Clinical and",
    "biochemical remission with prior biologics as covariate and time as",
    "regressor'.",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events)",
    concentration = "prob_clinbiorem (probability of clinical and biochemical remission, 0-1)"
  )

  covariateData <- list(
    PRIOR_BIO = list(
      description = "Prior exposure to another biologic before starting ustekinumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (biologic-naive; ustekinumab was the first biologic)",
      notes = "89% (74/83) were biologic-experienced (Table 4).",
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
    logite0 <- 0.006; label("Logit of clinical and biochemical remission for a biologic-naive patient at study baseline (unitless logit)")  # Table S1 'theta_pop = +0.006 +/- 0.85'
    e_priorbio_clinbiorem <- -2.48; label("Log-odds shift of remission for biologic-experienced vs naive patients (unitless logit)")  # Table S1 'theta_prior biologics = -2.48 +/- 0.87 (p = 0.004)'
    e_time_clinbiorem <- 0.0024; label("Log-odds of remission per day of follow-up (1/day)")  # Table S1 'theta_time = +0.0024 +/- 0.0017 (p = 0.158)'

    addSd_prob_clinbiorem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # t is days from study baseline: Results 3.3 evaluates the model at
    # 279 days, the median follow-up from baseline (Table 4)
    logit_clinbiorem <- logite0 + e_priorbio_clinbiorem * PRIOR_BIO + e_time_clinbiorem * t
    prob_clinbiorem <- expit(logit_clinbiorem)
    prob_clinbiorem ~ add(addSd_prob_clinbiorem)
  })
}
