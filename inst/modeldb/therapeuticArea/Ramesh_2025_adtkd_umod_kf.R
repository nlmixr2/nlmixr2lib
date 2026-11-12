Ramesh_2025_adtkd_umod_kf <- function() {
  description <- "Parametric Weibull time-to-event model for the age at kidney failure (eGFR <= 10 mL/min/1.73 m^2, dialysis initiation) in adults with autosomal dominant tubulointerstitial kidney disease caused by UMOD mutations (ADTKD-UMOD; Ramesh 2025). The hazard h(age) = (p / Te) * (age / Te)^(p - 1) has a characteristic age Te (survival exp(-1) = 0.368) that is shorter in males (exponential effect) and longer with older age at the first clinic visit (power effect). No drug input and no random effects; the model time is the subject's age in years, and the outputs are the hazard, the cumulative hazard and the kidney-failure-free survival probability from birth."
  reference <- paste(
    "Ramesh SS, Rogge M, Kidd KO, Williams AH, Yoon DY, Roignot J, Blakeslee K,",
    "Bleyer AJ, Kim S. (2025). Quantifying clinical and genetic factors influencing",
    "rate and severity of autosomal dominant tubulointerstitial kidney disease",
    "progression. J Pharmacokinet Pharmacodyn 52(4):41.",
    "doi:10.1007/s10928-025-09989-0.",
    "Parameter estimates: Table 2B (ADTKD-UMOD column); hazard form: Equation 2",
    "and Supplementary Material S2 (Monolix weibull_model_singleEvent).",
    sep = " "
  )
  vignette <- "Ramesh_2025_adtkd"
  units <- list(
    time = "year (the subject's age; the model time IS age)",
    dosing = "n/a (natural-history model, no drug input)",
    concentration = "probability (the output `sur` is the probability of being free of kidney failure at that age)"
  )

  covariateData <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = "Female is the reference; males have Te multiplied by exp(-0.07), about 3.8 years earlier (Discussion). The paper's beta_Te_Male applies to (1 - SEXF).",
      source_name = "SEX (Male vs Female)"
    ),
    AGE_BASE = list(
      description = "Age at the subject's first clinic visit (AFCV), time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE_BASE / 35.2138)^0.59 on Te. The paper does not print the centring value for the time-to-event models; 35.2138 is the observation-weighted mean of the same covariate printed in Supplementary S2 for the eGFR models, and it reproduces the Figure 4 model-predicted medians by age stratum (see the vignette). Distinct from the subject's current age, which is the model time.",
      source_name = "BL_VA (AFCV in the paper)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 371L,
    n_studies = 1L,
    n_observations = "105 kidney-failure events among 371 UMOD patients (Table 1); event ages right-censored at the end of follow-up for the rest. Fitted to a random 4:1 training split.",
    age_range = "18-80 years at assessment (whole data set); age at kidney failure median 52 years (95% interval 31-64)",
    age_median = "36.77 years at the first clinic visit (95% interval 18.34-58.93)",
    sex_female_pct = 58,
    race_ethnicity = "92% 'American White' and 89% not Hispanic or Latinx across the whole data set.",
    disease_state = "Autosomal dominant tubulointerstitial kidney disease with a pathogenic UMOD (uromodulin) variant.",
    dose_range = "Not applicable -- natural-history model, no drug.",
    regions = "United States (Wake Forest University School of Medicine inherited-kidney-disease registry, referrals since 1996)",
    notes = "Table 1 ADTKD-UMOD column: 215 female / 156 male. Covariates explored in the time-to-event models: UMOD sub-mutation family, sex and age at first clinic visit (Methods 'Covariates'); sex and age at first clinic visit were retained (Table 2B)."
  )

  ini({
    # Ramesh 2025 Table 2B, ADTKD-UMOD 'Final estimates (%RSE)' column. The
    # Supplementary S2 project for this model lists p_pop = 20 and Te_pop =
    # 64.29 with sub-mutation covariates: those are initial values of an
    # earlier covariate model, not the final model of Table 2B.
    #
    # The Weibull scale is carried as the canonical rate lambda = 1 / Te;
    # model() rebuilds Te so the paper's covariate exponents apply to Te
    # with their published signs.
    llam_haz <- log(1 / 54.32); label("Weibull rate lambda = 1/Te for a female at the reference age at first clinic visit (1/year)") # Table 2B Te = 54.32 years (RSE 1.86%)
    lalfa_haz <- log(8.69); label("Weibull shape p (unitless)") # Table 2B p = 8.69 (RSE 7.50%)
    e_male_te <- -0.07; label("Exponential effect of male sex on Te (unitless)") # Table 2B beta_Te_Male = -0.07 (RSE 33.8%)
    e_age_base_te <- 0.59; label("Power exponent of age at first clinic visit on Te (unitless)") # Table 2B beta_Te_AFCV = 0.59 (RSE 6.05%)
    # No interindividual variability ('no-variability' in Supplementary S2)
    # and no residual error: the time-to-event likelihood is carried by the
    # event / censoring density.
  })
  model({
    age <- t

    te <- exp(-llam_haz) *
      exp(e_male_te * (1 - SEXF)) *
      (AGE_BASE / 35.2138)^e_age_base_te
    lam <- 1 / te
    alfa <- exp(lalfa_haz)

    # Ramesh 2025 Equation 2, h = (p / Te) * (t / Te)^(p - 1), written with
    # lambda = 1 / Te. The cumulative hazard has the closed form (t / Te)^p,
    # so no ODE state is needed.
    hazard <- lam * alfa * (lam * age)^(alfa - 1)
    cumhaz <- (lam * age)^alfa
    sur <- exp(-cumhaz)
  })
}
