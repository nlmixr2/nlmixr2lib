Huh_2025_ritlecitinib_herpes_zoster <- function() {
  description <- paste(
    "Poisson exposure-response model for the first occurrence of a",
    "herpes zoster infection adverse event during oral ritlecitinib",
    "treatment, pooled across",
    "healthy volunteers and patients with alopecia areata or vitiligo.",
    "The mean count over the interval from the first dose is",
    "lambda = t * exp(theta0 + theta_Cave * sqrt(CAV) + covariates), with",
    "t the time since the first dose in days and CAV the time-weighted",
    "average ritlecitinib concentration (cumulative AUC divided by t).",
    "Age shifts the log baseline rate linearly around 34 years. The",
    "exposure slope was not statistically significant (95% CI including",
    "zero; only 22 events). With only the",
    "first event per subject and a log(t) offset, the model is",
    "equivalent to a constant-hazard (exponential) time-to-first-event",
    "model with hazard exp(theta0 + theta_Cave * sqrt(CAV) + covariates)",
    "per day. No random effects."
  )

  reference <- paste(
    "Huh Y, Qiu R, Prybylski J, Wojciechowski J, Wang Y, Purohit VS.",
    "Comprehensive Safety Exposure-Response Analysis to Support",
    "Ritlecitinib Dose Selection. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(6):1108-1118. doi:10.1002/psp4.70030.",
    sep = " "
  )

  vignette <- "Huh_2025_ritlecitinib_safety"

  # There is no PK component: exposure is supplied as CAV, the time-weighted
  # average concentration from the first dose, derived in Huh 2025 from the
  # empirical Bayes estimates of the Wojciechowski 2023 final population PK
  # model (modellib('Wojciechowski_2023_ritlecitinib_final')).
  depends <- "CAV"

  units <- list(
    time = "day",
    dosing = "(none; exposure enters through the CAV covariate)",
    concentration = "ng/mL (CAV covariate)"
  )

  covariateData <- list(
    CAV = list(
      description = "Time-weighted average ritlecitinib plasma concentration from the first dose to the record time: cumulative AUC(0-t) divided by t.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Huh 2025 Methods ('Exposure-Response Modeling of Incidence of Infections, Herpes Zoster, and Rash') and Figure S1: 'the cumulative AUC divided by the time of the adverse event with respect to the time of the first dose', derived from the empirical Bayes estimates of the Wojciechowski 2023 final population PK model.",
        "Time-varying: on an event record it is the value at the time of the event; on a censored record the value at the end of observation. 0 for placebo.",
        "Enters the log rate through its square root (see the slope parameter); geometric mean steady-state values are 27 ng/mL at 30 mg QD and 257 ng/mL at 200 mg QD in alopecia areata patients (Results)."
      ),
      source_name = "CAVE"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on the log baseline rate centered at 34 years, log(rate) + 0.0396 * (AGE - 34) (Huh 2025 Table 2 row 'Effect of age (referenced to 34 years) on theta0'; Supplement 2 covariate AGEcov). The Results print the relationship as theta0 = population theta0 x (Age/34)^0.0396, but applied to the negative log-scale intercept that form would LOWER the rate with age and give a 65-versus-34-year risk ratio of 0.75, against the reported 3.62 and the text's 'elderly patients of 65 years old have higher' incidence; the linear log-scale form gives 3.41 and 0.47 for 65 and 15 years (reported 3.62 and 0.500).",
      source_name = "AGEcov"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1744L,
    n_events = 22L,
    n_studies = 9L,
    disease_state = "Healthy volunteers (121), alopecia areata patients (1259; adolescents and adults) and vitiligo patients (364)",
    dose_range = "Healthy volunteers 5-800 mg single doses and up to 400 mg/day multiple doses; alopecia areata and vitiligo 10-50 mg QD with or without 100-200 mg QD loading doses (Table 1)",
    regions = "Global",
    notes = "Huh 2025 Results: 1744 subjects (121 healthy volunteers, 1259 alopecia areata and 364 vitiligo patients) from phase 1 studies B7981001, B7981003, B7981008, B7981022 and B7981036, phase 2/3 alopecia areata studies B7931005, B7981015 and B7981032 and phase 2b vitiligo study B7981019 (Table 1); 22 subjects had at least one herpes zoster event (MedDRA preferred terms herpes zoster and herpes zoster infection neurological, Table S1). Only the first event per subject was analysed. Fitted with glm(family = poisson) in R (Supplement 2)."
  )

  ini({
    lbase <- -11.3; label("Log baseline herpes zoster rate at CAV = 0 for a 34-year-old subject (log events/day)") # Table 2 Herpes Zoster ER model: Intercept (theta 0) = -11.3 (95% CI -12.2, -10.5)
    e_age_lbase <- 0.0396; label("Log-scale effect per year of age above 34 years on the baseline rate (1/year)") # Table 2: Effect of age (referenced to 34 years) on theta 0 = 0.0396 (95% CI 0.0106, 0.0686)
    slope <- 0.0873; label("Slope of the log rate on the square root of CAV (per sqrt(ng/mL))") # Table 2: Slope (theta Cave; unit/ng/mL) = 0.0873 (95% CI -0.014, 0.189); acts on sqrt(CAV), see vignette
  })
  model({
    # Huh 2025 Equation 3: log(lambda) = theta0 + log(t) + theta_Cave * f(Cave)
    # + covariates, with f = square root (selected by AIC, Supplement 2).
    aerate <- exp(lbase + e_age_lbase * (AGE - 34) + slope * sqrt(CAV))

    # Expected number of first events over the interval from the first dose
    # (time-interval offset t in days).
    lambda <- time * aerate
    herpes_zoster_count ~ pois(lambda)
  })
}
