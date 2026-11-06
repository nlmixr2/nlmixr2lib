Yarur_2025_ustekinumab_clinbiorem_ctrough <- function() {
  description <- paste0(
    "Logistic exposure-response model for combined CLINICAL AND ",
    "BIOCHEMICAL REMISSION versus ustekinumab TROUGH CONCENTRATION and ",
    "prior biologic exposure in patients with inflammatory bowel disease ",
    "on ustekinumab maintenance (Yarur 2025, 83 patients, 312 dose ",
    "cycles). Remission is Harvey-Bradshaw Index < 5 (Crohn's disease) ",
    "or partial Mayo score < 2 (ulcerative colitis) together with ",
    "CRP < 3 mg/L. The trough is the HMSA-measured serum concentration ",
    "supplied as the CTROUGH column. The concentration slope is not ",
    "significant (+0.024 log-odds per ug/mL, P = 0.533); prior biologic ",
    "exposure drives the outcome. Typical-value model: the ",
    "between-subject variance of the mixed-effects logistic fit was not ",
    "reported. Reproduces Figure 2A. One of five ",
    "Yarur_2025_ustekinumab_* exposure-response models."
  )
  reference <- paste(
    "Yarur AJ, Dervieux T, Ungaro R, Spencer EA, Bruss A, Nunez L,",
    "Berens B, Vermeire S, Wang Z, Panetta JC, Dreesen E, Dubinsky MC.",
    "Ustekinumab Drug Clearance Is Better Associated with Disease Control",
    "than Serum Trough Concentrations in a Prospective Cohort of",
    "Inflammatory Bowel Disease. Pharmaceutics. 2025;17(2):187.",
    "doi:10.3390/pharmaceutics17020187.",
    "Coefficients from Supplementary Table S1, row 'Clinical & biochemical",
    "remission with prior biologics as covariate and UST concentration as",
    "regressor'.",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_clinbiorem (probability of clinical and biochemical remission, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Ustekinumab serum trough concentration for the dose cycle (HMSA)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Measured immediately before the SC injection with the ",
        "drug-tolerant homogeneous mobility shift assay (LLOQ 0.9, ULOQ ",
        "25 ug/mL; Methods 2.2). HMSA concentrations run about 1.5-fold ",
        "higher than solid-phase ELISA (Results 3.1), so ELISA troughs ",
        "are not interchangeable here. Enters the logit linearly and ",
        "UNCENTERED. Study median 5.0 ug/mL, IQR 3.0-8.7 (Table 4); ",
        "Figure 2A spans 1-25 ug/mL. Can be simulated with ",
        "modellib('Yarur_2025_ustekinumab')."
      ),
      source_name = "UST concentrations"
    ),
    PRIOR_BIO = list(
      description = "Prior exposure to another biologic before starting ustekinumab",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (biologic-naive; ustekinumab was the first biologic)",
      notes = "89% (74/83) were biologic-experienced (Table 4). Figure 2A labels the two levels 'naive' and 'experienced'.",
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
    logite0 <- 0.22; label("Logit of clinical and biochemical remission for a biologic-naive patient at zero trough (unitless logit)")  # Table S1 'theta_pop = 0.22 +/- 0.84'
    e_priorbio_clinbiorem <- -2.36; label("Log-odds shift of remission for biologic-experienced vs naive patients (unitless logit)")  # Table S1 'theta_prior biologics = -2.36 +/- 0.85 (p = 0.005)'
    e_ctrough_clinbiorem <- 0.024; label("Log-odds of remission per ug/mL of ustekinumab trough (mL/ug)")  # Table S1 'theta_conc = +0.024 +/- 0.045 (p = 0.533)'

    addSd_prob_clinbiorem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_clinbiorem <- logite0 + e_priorbio_clinbiorem * PRIOR_BIO + e_ctrough_clinbiorem * CTROUGH
    prob_clinbiorem <- expit(logit_clinbiorem)
    prob_clinbiorem ~ add(addSd_prob_clinbiorem)
  })
}
