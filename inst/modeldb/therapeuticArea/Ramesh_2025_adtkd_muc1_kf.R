Ramesh_2025_adtkd_muc1_kf <- function() {
  description <- "Parametric Weibull time-to-event model for the age at kidney failure (eGFR <= 10 mL/min/1.73 m^2, dialysis initiation) in adults with autosomal dominant tubulointerstitial kidney disease caused by MUC1 mutations (ADTKD-MUC1; Ramesh 2025). The hazard h(age) = (p / Te) * (age / Te)^(p - 1) has a characteristic age Te (survival exp(-1) = 0.368) that is longer with older age at the first clinic visit (power effect); sex was not retained. No drug input and no random effects; the model time is the subject's age in years, and the outputs are the hazard, the cumulative hazard and the kidney-failure-free survival probability from birth."
  reference <- paste(
    "Ramesh SS, Rogge M, Kidd KO, Williams AH, Yoon DY, Roignot J, Blakeslee K,",
    "Bleyer AJ, Kim S. (2025). Quantifying clinical and genetic factors influencing",
    "rate and severity of autosomal dominant tubulointerstitial kidney disease",
    "progression. J Pharmacokinet Pharmacodyn 52(4):41.",
    "doi:10.1007/s10928-025-09989-0.",
    "Parameter estimates: Table 2B (ADTKD-MUC1 column); hazard form: Equation 2",
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
    AGE_BASE = list(
      description = "Age at the subject's first clinic visit (AFCV), time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE_BASE / 35.2138)^0.81 on Te. The paper does not print the centring value for the time-to-event models; 35.2138 is the observation-weighted mean of the same covariate printed in Supplementary S2 for the eGFR models, and it reproduces the Figure 5 model-predicted medians by age stratum (see the vignette). Distinct from the subject's current age, which is the model time.",
      source_name = "BL_VA (AFCV in the paper)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the time-to-event covariate analysis but not retained (Discussion: 'For MUC1, sex was not a significant covariate influencing the time to KF')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 233L,
    n_studies = 1L,
    n_observations = "145 kidney-failure events among 233 MUC1 patients (Table 1); event ages right-censored at the end of follow-up for the rest. Fitted to a random 4:1 training split.",
    age_range = "18-80 years at assessment (whole data set); age at kidney failure median 49 years (95% interval 27-74)",
    age_median = "37.68 years at the first clinic visit (95% interval 20.56-64.99)",
    sex_female_pct = 61,
    race_ethnicity = "92% 'American White' and 89% not Hispanic or Latinx across the whole data set.",
    disease_state = "Autosomal dominant tubulointerstitial kidney disease with a pathogenic MUC1 (mucin-1) variant.",
    dose_range = "Not applicable -- natural-history model, no drug.",
    regions = "United States (Wake Forest University School of Medicine inherited-kidney-disease registry, referrals since 1996)",
    notes = "Table 1 ADTKD-MUC1 column: 142 female / 91 male. Covariates explored in the time-to-event models: sex and age at first clinic visit (Methods 'Covariates'); only age at first clinic visit was retained (Table 2B; Discussion: 'sex was not a significant covariate')."
  )

  ini({
    # Ramesh 2025 Table 2B, ADTKD-MUC1 'Final estimates (%RSE)' column. The
    # Supplementary S2 project for this model lists p_pop = 4 and Te_pop =
    # 58.82 with no covariate: those are initial values of the base model,
    # not the final model of Table 2B.
    #
    # The Weibull scale is carried as the canonical rate lambda = 1 / Te;
    # model() rebuilds Te so the paper's covariate exponents apply to Te
    # with their published signs.
    llam_haz <- log(1 / 51.54); label("Weibull rate lambda = 1/Te at the reference age at first clinic visit (1/year)") # Table 2B Te = 51.54 years (RSE 1.97%)
    lalfa_haz <- log(6.82); label("Weibull shape p (unitless)") # Table 2B p = 6.82 (RSE 8.81%)
    e_age_base_te <- 0.81; label("Power exponent of age at first clinic visit on Te (unitless)") # Table 2B beta_Te_AFCV = 0.81 (RSE 7.01%)
    # No interindividual variability ('no-variability' in Supplementary S2)
    # and no residual error: the time-to-event likelihood is carried by the
    # event / censoring density.
  })
  model({
    age <- t

    te <- exp(-llam_haz) * (AGE_BASE / 35.2138)^e_age_base_te
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
