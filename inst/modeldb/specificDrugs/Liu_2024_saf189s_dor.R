Liu_2024_saf189s_dor <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for DURATION OF ",
    "RESPONSE in Chinese adults with ALK-positive or ROS1-positive ",
    "advanced non-small cell lung cancer treated with the ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor SAF-189s ",
    "(Liu 2024, n = 192, the responder subset of the phase II analysis ",
    "set of study SAF001, NCT04237805, 80-210 mg orally once daily). ",
    "The model returns the RELATIVE hazard exp(0.0059 * CTROUGH) where ",
    "CTROUGH is the individual SAF-189s steady-state trough ",
    "concentration in ng/mL at the patient's most prevalent dose, ",
    "entering LINEARLY rather than on a log scale, so the printed ",
    "hazard ratio 1.006 is per ng/mL. The relationship is clearly NOT ",
    "statistically significant (95% CI 0.996-1.016, p = 0.251) and is ",
    "the weakest of the paper's four exposure-efficacy analyses; the ",
    "model is packaged so that this null result is reproducible. NO ",
    "BASELINE HAZARD IS ENCODED: a Cox regression is semiparametric and ",
    "leaves h0(t) unspecified, so the absolute survivor function is a ",
    "quantity the fit never produced rather than an unreported ",
    "parameter. There is no PK layer and no ODE: the exposure metric is ",
    "supplied as a data column, derived in the source analysis from the ",
    "individual post hoc parameters of the companion population PK ",
    "model, packaged as Liu_2024_saf189s. No between-subject random ",
    "effect and no residual error are estimated. Companion ",
    "progression-free-survival model in Liu_2024_saf189s_pfs; six ",
    "companion models in the Liu_2024_saf189s_* family."
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
    concentration = "hr (relative hazard of loss of response, unitless; also lhr, the log relative hazard)"
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
        "Enters LINEARLY, not on a log scale, exactly as in the",
        "companion PFS model. Table 4 prints an Estimate of 0.0059",
        "alongside a hazard ratio of 1.006, and",
        "exp(0.0059) = 1.00592, which confirms the coefficient is per",
        "ng/mL on the untransformed scale.",
        "Because there is no intercept to absorb a rescaling, supplying",
        "ug/mL instead of ng/mL would silently destroy the result",
        "rather than shift it.",
        "For calibration, Liu 2024 Table 3 gives the phase II analysis",
        "set's Cmin,ss quartile boundaries as 18, 60.7, 79, 108 and",
        "182 ng/mL; the DOR subset is the 192 of those 244 patients who",
        "responded, and no separate exposure distribution is reported",
        "for it."
      ),
      source_name        = "C min,ss (steady-state minimum concentration) estimated using the most prevalent dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 192L,
    n_studies      = 1L,
    n_observations = "192 right-censored duration-of-response times (one per responding patient)",
    age_range      = "SAF001 phase II: median 54.1 years, range 20.0-84.0 (Liu 2024 Supplementary Table 1)",
    weight_range   = "median 63.2 kg, range 37.3-92.5 across the PK-evaluable cohort",
    sex_female_pct = 53.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer, restricted to patients achieving an independent-review-committee-assessed complete or partial response; duration measured from first response to progression or death",
    dose_range     = "SAF-189s 80, 120, 160 or 210 mg orally once daily in 21-day cycles (phase IIa) or 160 mg once daily (phase IIb)",
    regions        = "China",
    notes          = paste0(
      "The DOR analysis set is the 192 responders among the 244 phase ",
      "II patients used for the companion ORR and PFS models (Liu 2024 ",
      "Results, Exposure-response analysis; Table 4). That is a ",
      "RESPONDER-CONDITIONED set, so any exposure effect it estimates ",
      "is conditional on having responded and is not comparable with ",
      "the unconditioned PFS estimate. Liu 2024 states that the DOR ",
      "data were IMMATURE at the time of analysis."
    )
  )

  ini({
    # ==================================================================
    # Liu 2024 Table 4, "Cox regression models for PFS and DOR in all
    # patients", DOR row: n = 192, Parameter C min,ss, Estimate 0.0059,
    # SE 0.0052, Hazard ratio 1.006 (95% CI 0.996-1.016), p = 0.251.
    #
    # The coefficient below is the natural log of the printed hazard
    # ratio, which is the scale a Cox model estimates on, and it is the
    # printed Estimate itself: exp(0.0059) = 1.00592, rounding to the
    # printed 1.006. Because the paper prints the Estimate directly,
    # nothing here is back-solved or digitised. The estimate is roughly
    # one standard error from zero (0.0059 against SE 0.0052), so it is
    # statistically indistinguishable from no effect.
    #
    # NO BASELINE HAZARD IS ENCODED, for the same semiparametric reason
    # as in the companion PFS model: h0(t) is left completely
    # unspecified by a Cox fit, so it is not an unreported parameter but
    # a quantity the method never produces. This model returns the
    # RELATIVE hazard hr only and deliberately does not offer a
    # survivor function. Liu 2024 Figure 8B gives the paper's
    # Kaplan-Meier DOR curves stratified at the median exposure as an
    # empirical description of the absolute time course; the log-rank
    # test across that split was not significant. This matches the
    # relative-hazard-only pattern of
    # Rayner_2013_oseltamivir_shedding.R.
    #
    # No IIV and no residual error: the paper reports a coefficient, its
    # standard error and a hazard ratio with no variance components, so
    # no observation endpoint is declared.
    # ==================================================================

    e_ctrough_haz <- 0.0059 ; label("Log hazard ratio for loss of response per 1 ng/mL increase in SAF-189s steady-state trough concentration (log scale; HR 1.006)")  # Liu 2024 Table 4, DOR row: Estimate 0.0059, SE 0.0052, HR 1.006 (95% CI 0.996-1.016), p = 0.251
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor and relative hazard against a hypothetical
    # zero-exposure patient. Under the proportional-hazards assumption
    # hr is constant in time, so the subject hazard is
    # h(t) = h0(t) * hr for any baseline h0(t) the user supplies. Only
    # RATIOS of hr between two exposures are meaningful, since the
    # CTROUGH = 0 reference is outside the observed data.
    # ------------------------------------------------------------------
    lhr <- e_ctrough_haz * CTROUGH

    hr <- exp(lhr)
  })
}
