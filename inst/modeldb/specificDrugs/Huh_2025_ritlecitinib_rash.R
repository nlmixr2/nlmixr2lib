Huh_2025_ritlecitinib_rash <- function() {
  description <- paste(
    "Poisson exposure-response model for the first occurrence of a rash",
    "adverse event (rash/dermatitis, drug eruption, urticaria, acne and",
    "related skin terms) during oral ritlecitinib treatment, pooled across",
    "healthy volunteers and patients with alopecia areata or vitiligo.",
    "The mean count over the interval from the first dose is",
    "lambda = t * exp(theta0 + theta_Cave * sqrt(CAV) + covariates), with",
    "t the time since the first dose in days and CAV the time-weighted",
    "average ritlecitinib concentration (cumulative AUC divided by t).",
    "Patient type (alopecia areata or vitiligo versus healthy",
    "volunteer), age (linear around 34 years), geographical location",
    "(Asia versus North America/Europe/Rest of World) and female sex",
    "shift the log baseline rate. With only the",
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
    DIS_ALOPECIA_AREATA = list(
      description = "Alopecia areata patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer, when DIS_VITILIGO is also 0)",
      notes = "Huh 2025 Table 2 row 'Effect of AA patients on theta0' (Supplement 2 covariate PTSTf, patient type HV / AA / vitiligo with healthy volunteers as the reference). Mutually exclusive with DIS_VITILIGO.",
      source_name = "PTSTf"
    ),
    DIS_VITILIGO = list(
      description = "Vitiligo patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer, when DIS_ALOPECIA_AREATA is also 0)",
      notes = "Huh 2025 Table 2 row 'Effect of Vitiligo patients on theta0' (Supplement 2 covariate PTSTf). Vitiligo patients came from phase 2b study B7981019 (Table 1). Mutually exclusive with DIS_ALOPECIA_AREATA.",
      source_name = "PTSTf"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on the log baseline rate centered at 34 years, log(rate) - 0.0185 * (AGE - 34) (Huh 2025 Table 2 row 'Effect of age (referenced to 34 years) on theta0'; Supplement 2 covariate AGEcov). The same linear log-scale form as the herpes zoster model, whose printed power form is not usable on a negative log-scale intercept (see Huh_2025_ritlecitinib_herpes_zoster). Gives a 15-versus-34-year risk ratio of 1.42 (reported 1.56, Table 3).",
      source_name = "AGEcov"
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
    n_events = 434L,
    n_studies = 9L,
    disease_state = "Healthy volunteers (121), alopecia areata patients (1259; adolescents and adults) and vitiligo patients (364)",
    dose_range = "Healthy volunteers 5-800 mg single doses and up to 400 mg/day multiple doses; alopecia areata and vitiligo 10-50 mg QD with or without 100-200 mg QD loading doses (Table 1)",
    regions = "Global; Asia versus North America/Europe/Rest of World retained as a covariate",
    notes = "Huh 2025 Results: 1744 subjects (121 healthy volunteers, 1259 alopecia areata and 364 vitiligo patients) from phase 1 studies B7981001, B7981003, B7981008, B7981022 and B7981036, phase 2/3 alopecia areata studies B7931005, B7981015 and B7981032 and phase 2b vitiligo study B7981019 (Table 1); 434 subjects had at least one rash event (rash/dermatitis, drug eruption, urticaria, acne and other skin terms; MedDRA preferred terms in Table S1). Only the first event per subject was analysed. Fitted with glm(family = poisson) in R (Supplement 2)."
  )

  ini({
    lbase <- -6.02; label("Log baseline rash rate at CAV = 0 for a 34-year-old male healthy volunteer outside Asia (log events/day)") # Table 2 Rash ER model: Intercept (theta 0) = -6.02 (95% CI -6.41, -5.63)
    e_alopecia_areata_lbase <- -2.86; label("Log-scale effect of alopecia areata (versus healthy volunteer) on the baseline rate (unitless)") # Table 2: Effect of AA patients on theta 0 = -2.86 (95% CI -3.27, -2.46)
    e_vitiligo_lbase <- -2.73; label("Log-scale effect of vitiligo (versus healthy volunteer) on the baseline rate (unitless)") # Table 2: Effect of Vitiligo patients on theta 0 = -2.73 (95% CI -3.21, -2.25)
    e_region_asia_lbase <- 0.632; label("Log-scale effect of Asia geographical location on the baseline rate (unitless)") # Table 2: Effect of Asia (geographical location) on theta 0 = 0.632 (95% CI 0.408, 0.856)
    e_age_lbase <- -0.0185; label("Log-scale effect per year of age above 34 years on the baseline rate (1/year)") # Table 2: Effect of age (referenced to 34 years) on theta 0 = -0.0185 (95% CI -0.0258, -0.0112)
    e_sexf_lbase <- 0.521; label("Log-scale effect of female sex on the baseline rate (unitless)") # Table 2: Effect of female sex on theta 0 = 0.521 (95% CI 0.307, 0.734)
    slope <- 0.142; label("Slope of the log rate on the square root of CAV (per sqrt(ng/mL))") # Table 2: Slope (theta Cave; unit/ng/mL) = 0.142 (95% CI 0.129, 0.155); acts on sqrt(CAV), see vignette
  })
  model({
    # Huh 2025 Equation 3: log(lambda) = theta0 + log(t) + theta_Cave * f(Cave)
    # + covariates, with f = square root (selected by AIC, Supplement 2).
    aerate <- exp(lbase +
      e_alopecia_areata_lbase * DIS_ALOPECIA_AREATA +
      e_vitiligo_lbase * DIS_VITILIGO +
      e_region_asia_lbase * REGION_ASIA +
      e_age_lbase * (AGE - 34) +
      e_sexf_lbase * SEXF +
      slope * sqrt(CAV))

    # Expected number of first events over the interval from the first dose
    # (time-interval offset t in days).
    lambda <- time * aerate
    rash_count ~ pois(lambda)
  })
}
