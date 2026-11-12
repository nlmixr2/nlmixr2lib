Huh_2025_ritlecitinib_infections <- function() {
  description <- paste(
    "Poisson exposure-response model for the first occurrence of an",
    "infection adverse event (moderate, severe, or leading to",
    "discontinuation) during oral ritlecitinib treatment, pooled across",
    "healthy volunteers and patients with alopecia areata or vitiligo.",
    "The mean count over the interval from the first dose is",
    "lambda = t * exp(theta0 + theta_Cave * sqrt(CAV) + covariates), with",
    "t the time since the first dose in days and CAV the time-weighted",
    "average ritlecitinib concentration (cumulative AUC divided by t).",
    "Geographical location (Asia versus North America/Europe/Rest of",
    "World) and female sex shift the log baseline rate. With only the",
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
    REGION_ASIA = list(
      description = "Asia geographical-location indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (North America, Europe or Rest of World)",
      notes = "Huh 2025 Table 2 row 'Effect of Asia (geographical location) on theta0' (Supplement 2 covariate LOCATIONAf: Asia versus North America/Europe/Rest of World). A study-site location, not self-reported race (race was screened separately as RACEf and not retained).",
      source_name = "LOCATIONAf"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Huh 2025 Table 2 row 'Effect of female sex on theta0' (Supplement 2 covariate SEXf).",
      source_name = "SEXf"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1744L,
    n_events = 241L,
    n_studies = 9L,
    disease_state = "Healthy volunteers (121), alopecia areata patients (1259; adolescents and adults) and vitiligo patients (364)",
    dose_range = "Healthy volunteers 5-800 mg single doses and up to 400 mg/day multiple doses; alopecia areata and vitiligo 10-50 mg QD with or without 100-200 mg QD loading doses (Table 1)",
    regions = "Global; Asia versus North America/Europe/Rest of World retained as a covariate",
    notes = "Huh 2025 Results: 1744 subjects (121 healthy volunteers, 1259 alopecia areata and 364 vitiligo patients) from phase 1 studies B7981001, B7981003, B7981008, B7981022 and B7981036, phase 2/3 alopecia areata studies B7931005, B7981015 and B7981032 and phase 2b vitiligo study B7981019 (Table 1); 241 subjects had at least one infection event. Infection = moderate or severe infections and those leading to discontinuation (MedDRA preferred terms in Table S1). Only the first event per subject was analysed. Fitted with glm(family = poisson) in R (Supplement 2)."
  )

  ini({
    lbase <- -9.00; label("Log baseline infection rate at CAV = 0 for a male subject outside Asia (log events/day)") # Table 2 Infections ER model: Intercept (theta 0) = -9.00 (95% CI -9.34, -8.67)
    e_region_asia_lbase <- -1.31; label("Log-scale effect of Asia geographical location on the baseline rate (unitless)") # Table 2: Effect of Asia (geographical location) on theta 0 = -1.31 (95% CI -1.90, -0.733)
    e_sexf_lbase <- 0.748; label("Log-scale effect of female sex on the baseline rate (unitless)") # Table 2: Effect of female sex on theta 0 = 0.748 (95% CI 0.449, 1.05)
    slope <- 0.0989; label("Slope of the log rate on the square root of CAV (per sqrt(ng/mL))") # Table 2: Slope (theta Cave; unit/ng/mL) = 0.0989 (95% CI 0.0718, 0.126); acts on sqrt(CAV), see vignette
  })
  model({
    # Huh 2025 Equation 3: log(lambda) = theta0 + log(t) + theta_Cave * f(Cave)
    # + covariates, with f = square root (selected by AIC, Supplement 2).
    aerate <- exp(lbase + e_region_asia_lbase * REGION_ASIA +
      e_sexf_lbase * SEXF + slope * sqrt(CAV))

    # Expected number of first events over the interval from the first dose
    # (time-interval offset t in days).
    lambda <- time * aerate
    infection_count ~ pois(lambda)
  })
}
