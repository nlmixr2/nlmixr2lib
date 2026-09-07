Liu_2024_saf189s_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for ",
    "PROGRESSION-FREE SURVIVAL in Chinese adults with ALK-positive or ",
    "ROS1-positive advanced non-small cell lung cancer treated with the ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor SAF-189s ",
    "(Liu 2024, n = 244, the phase II analysis set of study SAF001, ",
    "NCT04237805, 80-210 mg orally once daily). The model returns the ",
    "RELATIVE hazard exp(0.0071 * CTROUGH) where CTROUGH is the ",
    "individual SAF-189s steady-state trough concentration in ng/mL at ",
    "the patient's most prevalent dose, entering LINEARLY rather than ",
    "on a log scale, so the printed hazard ratio 1.007 is per ng/mL. ",
    "The relationship is NOT statistically significant (95% CI ",
    "1.000-1.015, p = 0.059), matching the Kaplan-Meier analysis in ",
    "which PFS above and below the median exposure did not separate; ",
    "the model is packaged so that this null result is reproducible. ",
    "NO BASELINE HAZARD IS ENCODED: a Cox regression is semiparametric ",
    "and leaves h0(t) unspecified, so the absolute survivor function is ",
    "a quantity the fit never produced rather than an unreported ",
    "parameter. There is no PK layer and no ODE: the exposure metric is ",
    "supplied as a data column, derived in the source analysis from the ",
    "individual post hoc parameters of the companion population PK ",
    "model, packaged as Liu_2024_saf189s. No between-subject random ",
    "effect and no residual error are estimated. Companion ",
    "duration-of-response model in Liu_2024_saf189s_dor; six companion ",
    "models in the Liu_2024_saf189s_* family."
  )
  reference <- paste(
    "Liu Y, Tan Y, Hu L, Li J, Yang J, Diao L, Yang J.",
    "Population pharmacokinetics and exposure-response analyses of",
    "SAF-189s in Chinese patients with ALK+/ROS1+ non-small cell lung",
    "cancer.",
    "Front Pharmacol. 2024;15:1418549.",
    "doi:10.3389/fphar.2024.1418549.",
    "Individual SAF-189s trough concentrations derive from the companion",
    "population pharmacokinetic model reported in the same paper; see",
    "modellib('Liu_2024_saf189s').",
    sep = " "
  )
  vignette <- "Liu_2024_saf189s"
  units <- list(
    time          = "n/a (semiparametric Cox relative hazard; the baseline hazard and therefore the time scale are unspecified by the fit)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "hr (relative hazard of progression or death, unitless; also lhr, the log relative hazard)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual SAF-189s steady-state trough plasma concentration (Cmin,ss), per subject. Supplied as data: this model has no PK layer, and the source analysis used the individual post hoc parameters of the companion population PK model together with the patient's most prevalent dose level (Liu 2024 Methods, E-R analysis).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration, at STEADY STATE, and",
        "computed at the MOST PREVALENT dose the patient received. The",
        "companion exposure-SAFETY models use AUCss computed at the",
        "FIRST dose instead, so the two exposure columns are not",
        "interchangeable even after a unit conversion.",
        "Enters LINEARLY, not on a log scale -- this is the one place",
        "where Liu 2024 departs from the log-transformed exposure it",
        "uses in every logistic fit. Table 4 prints an Estimate of",
        "0.0071 alongside a hazard ratio of 1.007, and",
        "exp(0.0071) = 1.00713, which confirms the coefficient is per",
        "ng/mL on the untransformed scale. A log-scale reading would",
        "give exp(0.0071 * log(C)) and could not reproduce the printed",
        "hazard ratio.",
        "Units are load-bearing in a different way than for the logistic",
        "models: because there is no intercept to absorb a rescaling,",
        "supplying ug/mL instead of ng/mL would change the relative",
        "hazard by a factor of exp(0.0071 * 999 * C) -- that is, it",
        "would silently destroy the result rather than shift it.",
        "For calibration, Liu 2024 Table 3 gives the analysis set's",
        "Cmin,ss quartile boundaries as 18, 60.7, 79, 108 and",
        "182 ng/mL. Over that full observed span the relative hazard",
        "moves only from exp(0.0071*18) = 1.14 to",
        "exp(0.0071*182) = 3.64 relative to a hypothetical",
        "zero-exposure patient, and by a factor of 3.2 across the",
        "observed range -- with a confidence interval whose lower bound",
        "is 1.000."
      ),
      source_name        = "C min,ss (steady-state minimum concentration) estimated using the most prevalent dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 244L,
    n_studies      = 1L,
    n_observations = "244 right-censored progression-free-survival times (one per patient)",
    age_range      = "SAF001 phase II: median 54.1 years, range 20.0-84.0 (Liu 2024 Supplementary Table 1)",
    weight_range   = "median 63.2 kg, range 37.3-92.5 across the PK-evaluable cohort",
    sex_female_pct = 53.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer; progression assessed by an independent review committee under RECIST 1.1",
    dose_range     = "SAF-189s 80, 120, 160 or 210 mg orally once daily in 21-day cycles (phase IIa) or 160 mg once daily (phase IIb)",
    regions        = "China",
    notes          = paste0(
      "Same 244-patient phase II analysis set as the companion ORR ",
      "model. Liu 2024 states explicitly that the PFS data were ",
      "IMMATURE at the time of analysis ('the PFS and DOR data were ",
      "immature for patients participating in the phase II study at ",
      "the time of this analysis'), which is the main reason to treat ",
      "the borderline p = 0.059 as uninformative rather than as a ",
      "near-miss. No median PFS is reported for this pooled set; the ",
      "Introduction quotes a median PFS of 16.5 months for ROS1-naive ",
      "patients in phase IIa from a separate report."
    )
  )

  ini({
    # ==================================================================
    # Liu 2024 Table 4, "Cox regression models for PFS and DOR in all
    # patients", PFS row: n = 244, Parameter C min,ss, Estimate 0.0071,
    # SE 0.0038, Hazard ratio 1.007 (95% CI 1.000-1.015), p = 0.059.
    #
    # The coefficient below is the natural log of the printed hazard
    # ratio, which is the scale a Cox model estimates on, and it is the
    # printed Estimate itself: exp(0.0071) = 1.00713, rounding to the
    # printed 1.007. Because the paper prints the Estimate directly,
    # nothing here is back-solved or digitised.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is left completely unspecified by the method, so it is not
    # an unreported parameter but a quantity the fit never produced.
    # This model therefore returns the RELATIVE hazard hr only and
    # deliberately does not offer a survivor function. Liu 2024
    # Figure 8A gives the paper's Kaplan-Meier PFS curves stratified at
    # the median exposure as an empirical description of the absolute
    # time course; the log-rank test across that split was not
    # significant. See the vignette Assumptions and deviations for why
    # a calibrated parametric baseline was rejected. This matches the
    # relative-hazard-only pattern of
    # Rayner_2013_oseltamivir_shedding.R.
    #
    # The regression was fitted with standard survival software rather
    # than NONMEM; the paper reports a coefficient, its standard error
    # and a hazard ratio with no variance components, so there is no
    # IIV and no residual error to encode and no observation endpoint
    # is declared.
    # ==================================================================

    e_ctrough_haz <- 0.0071 ; label("Log hazard ratio for progression or death per 1 ng/mL increase in SAF-189s steady-state trough concentration (log scale; HR 1.007)")  # Liu 2024 Table 4, PFS row: Estimate 0.0071, SE 0.0038, HR 1.007 (95% CI 1.000-1.015), p = 0.059
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor and relative hazard against a hypothetical
    # zero-exposure patient. Under the proportional-hazards assumption
    # hr is constant in time, so the subject hazard is
    # h(t) = h0(t) * hr for any baseline h0(t) the user supplies. Only
    # RATIOS of hr between two exposures are meaningful, since the
    # CTROUGH = 0 reference is outside the observed data.
    #
    # The coefficient is POSITIVE, so the fitted hazard of progression
    # rises with exposure -- the same counter-intuitive direction as the
    # companion ORR model's negative slope, and null for the same
    # reason (the confidence interval's lower bound is 1.000). Liu 2024
    # attributes the flat efficacy exposure-response to exposures
    # sitting on the plateau of the response curve over 80-210 mg,
    # with confounding by dose reduction and by the treatment-resistant
    # subgroup the likely source of the sign.
    # ------------------------------------------------------------------
    lhr <- e_ctrough_haz * CTROUGH

    hr <- exp(lhr)
  })
}
