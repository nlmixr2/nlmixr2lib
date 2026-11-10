Yarur_2025_ustekinumab_ehi50_ctrough <- function() {
  description <- paste0(
    "Logistic exposure-response model for an ENDOSCOPIC HEALING INDEX ",
    "(EHI) ABOVE 50 POINTS -- a serum-biomarker surrogate of ",
    "endoscopically active disease -- versus ustekinumab TROUGH ",
    "CONCENTRATION and prior biologic exposure in patients with ",
    "inflammatory bowel disease on ustekinumab maintenance (Yarur 2025, ",
    "83 patients, 303 dose cycles with an EHI). The trough is the ",
    "HMSA-measured serum concentration supplied as the CTROUGH column. ",
    "Coefficients are transcribed from Supplementary Table S1 ",
    "(-0.072 log-odds per ug/mL); the published Figure 3A curves are ",
    "steeper and order the prior-biologic groups the other way, so this ",
    "model does NOT reproduce Figure 3A -- see the vignette. ",
    "Typical-value model: the between-subject variance of the ",
    "mixed-effects logistic fit was not reported. One of five ",
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
    "biologics as covariate and UST concentration as regressor'.",
    sep = " "
  )
  vignette <- "Yarur_2025_ustekinumab"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_ehi50 (probability of an endoscopic healing index above 50 points, 0-1)"
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
        "higher than solid-phase ELISA (Results 3.1). Enters the logit ",
        "linearly and UNCENTERED. Study median 5.0 ug/mL, IQR 3.0-8.7 ",
        "(Table 4); Figure 3A spans 1-25 ug/mL. Can be simulated with ",
        "modellib('Yarur_2025_ustekinumab')."
      ),
      source_name = "UST concentrations"
    ),
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
    logite0 <- -1.85; label("Logit of EHI above 50 points for a biologic-naive patient at zero trough (unitless logit)")  # Table S1 'theta_pop = -1.85 +/- 1.55'
    e_priorbio_ehi50 <- -0.52; label("Log-odds shift of EHI above 50 points for biologic-experienced vs naive patients (unitless logit)")  # Table S1 'theta_prior biologics = -0.52 +/- 1.50 (p = 0.731)'
    e_ctrough_ehi50 <- -0.072; label("Log-odds of EHI above 50 points per ug/mL of ustekinumab trough (mL/ug)")  # Table S1 'theta_conc = -0.072 +/- 0.045 (p = 0.001)'

    addSd_prob_ehi50 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_ehi50 <- logite0 + e_priorbio_ehi50 * PRIOR_BIO + e_ctrough_ehi50 * CTROUGH
    prob_ehi50 <- expit(logit_ehi50)
    prob_ehi50 ~ add(addSd_prob_ehi50)
  })
}
