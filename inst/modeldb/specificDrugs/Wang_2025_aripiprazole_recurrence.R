Wang_2025_aripiprazole_recurrence <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for the time to recurrence",
    "of any mood episode in adults with bipolar I disorder maintained on",
    "aripiprazole once monthly (AOM) 400 mg or withdrawn to placebo in the",
    "randomised double-blind withdrawal phase of trial NCT01567527 (Wang 2025,",
    "final model). The exposure driver is CTROUGH, the model-predicted",
    "aripiprazole plasma concentration 672 h (28 days) after the first injection",
    "of the double-blind phase ('Ctau'), generated from the companion population",
    "PK model modellib('Wang_2022_aripiprazole'); for placebo subjects it is the",
    "residual concentration still washing out from the last open-label AOM dose.",
    "CTROUGH enters the log hazard linearly and uncentred: each 1 ng/mL increase",
    "multiplies the hazard of recurrence by exp(-0.00346) = 0.99655, a 0.34",
    "percent reduction. No other covariate was retained.",
    "Because the fit is a semiparametric Cox regression, the baseline hazard is",
    "unspecified by construction and is not reported; the model therefore",
    "returns the RELATIVE hazard hr against a hypothetical subject with",
    "CTROUGH = 0, not an absolute survivor function. Multiply hr by a",
    "user-supplied baseline hazard h0(t) to obtain a subject hazard.",
    "The model is algebraic and deterministic (no ODE state, no drug input, no",
    "IIV, no residual error).",
    "The separately fitted categorical model that dichotomises CTROUGH at",
    "95 ng/mL is modellib('Wang_2025_aripiprazole_recurrence_cutoff')."
  )
  reference <- paste(
    "Wang X, Phillips L, Harlin M, Bell Lynum KS, Larsen F, Such P,",
    "Madera-McDonough J, Yildirim M, Procyshyn RM, Chepke C, Passarell J.",
    "Exposure-Response Analysis for Aripiprazole Once-Monthly in Patients",
    "Diagnosed With Bipolar I Disorder.",
    "Clin Pharmacol Drug Dev. 2025;14(11):869-878. doi:10.1002/cpdd.1580.",
    "The coefficient is from the Results section 'Exposure-Response Model'; the",
    "model equation is in Appendix S4 of the Supporting Information.",
    "The upstream population PK model that generates the CTROUGH covariate",
    "(reproduced as Appendix S2 Table S2) is Wang X, Raoufinia A, Bihorel S, et",
    "al. Clin Pharmacol Drug Dev. 2022;11(2):150-164, doi:10.1002/cpdd.1022; see",
    "modellib('Wang_2022_aripiprazole').",
    sep = " "
  )
  vignette <- "Wang_2025_aripiprazole_bipolar"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method; the source analysis measured time in days from randomisation)",
    dosing = "n/a (exposure-response model; aripiprazole exposure enters as the CTROUGH covariate, not as a dose record)",
    concentration = "hr (hazard of mood-episode recurrence relative to CTROUGH = 0, unitless); not a drug concentration"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted aripiprazole plasma concentration 672 h (28 days) after the first injection of the double-blind withdrawal phase (Wang 2025 'Ctau'), i.e. the end-of-interval trough of the first double-blind dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "(a) TOTAL aripiprazole plasma concentration (parent drug only;",
        "dehydro-aripiprazole is not modelled).",
        "(b) Landmark: 672 h after the first double-blind (Phase D) injection,",
        "which follows at least 12 weeks of open-label AOM in Phase C, so for",
        "the AOM arm it represents the steady-state monthly trough; one",
        "time-fixed value per subject is applied over the whole follow-up.",
        "(c) Empirical-Bayes prediction: the Wang 2022 population PK model",
        "(modellib('Wang_2022_aripiprazole'), all parameters fixed) was applied",
        "to each subject's recorded dosing history, with a $MIX step assigning",
        "an unknown CYP2D6 phenotype, and the concentration was obtained by",
        "numerical integration from the individual Bayesian parameter",
        "estimates (Wang 2025 Methods and Appendix S3).",
        "Placebo records are NOT zero: subjects randomised to placebo still",
        "carried aripiprazole from their last open-label AOM dose, about 1344 h",
        "(56 days) earlier at the landmark. Reported means (SD): 170.5 (71.6)",
        "ng/mL in the AOM arm and 78.7 (36.1) ng/mL in the placebo arm (Wang",
        "2025 Results and Table 1)."
      ),
      source_name = "Ctau"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the Wang 2025 full-model backward elimination (P < 0.01 to retain) and not retained; no estimate reported."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the Wang 2025 full-model backward elimination and not retained; no estimate reported."
    ),
    SCORE_MADRS = list(
      description = "Montgomery-Asberg Depression Rating Scale total score at randomisation (baseline of Phase D)",
      units = "(SCORE_MADRS units, 0-60 score)",
      type = "continuous",
      notes = "Screened in the Wang 2025 full-model backward elimination and not retained; no estimate reported. Mean (SD) 2.7 (3.4) in the exposure-response set (Table S5)."
    ),
    SCORE_YMRS = list(
      description = "Young Mania Rating Scale total score at randomisation (baseline of Phase D)",
      units = "(0-60 score)",
      type = "continuous",
      notes = "Screened in the Wang 2025 full-model backward elimination and not retained; no estimate reported. Mean (SD) 2.8 (3.3) in the exposure-response set (Table S5). Named after the SCORE_<instrument> pattern of SCORE_MADRS; not a registered canonical, recorded here as documentation only."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the Wang 2025 full-model backward elimination and not retained; no estimate reported."
    ),
    RACE = list(
      description = "Race",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened in the Wang 2025 full-model backward elimination and not retained; no estimate reported, and the paper does not state how race was coded for the regression."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 265L,
    n_studies = 1L,
    n_events = 103L,
    age_range = "18-65 years by inclusion criterion; mean (SD) 40.1 (11.0) years in the exposure-response set (Wang 2025 Table S5)",
    weight_range = "Body mass index mean (SD) 30.2 (7.2) kg/m^2 in the exposure-response set (Table S5); body weight mean (SD) 87.17 (22.5) kg in the 420-patient population PK set (Table S3)",
    sex_female_pct = 57.4,
    race_ethnicity = c(
      White = 54.3,
      `Black or African American` = 27.9,
      Asian = 14.0,
      `American Indian or Alaska Native` = 0.8,
      Other = 3.0
    ),
    disease_state = paste(
      "Adults with bipolar I disorder (DSM-IV-TR, confirmed by the Mini",
      "International Neuropsychiatric Interview) who were in a current manic",
      "episode (YMRS >= 20) at trial entry and had reached stability (YMRS <= 12,",
      "MADRS <= 12, outpatient, no active suicidality) for at least 8 weeks on",
      "open-label AOM before randomisation. Baseline (randomisation) mean (SD)",
      "MADRS 2.7 (3.4) and YMRS 2.8 (3.3)."
    ),
    dose_range = paste(
      "Phase C (open-label, 12-28 weeks): AOM 400 mg gluteal IM every 4 weeks",
      "with oral aripiprazole 10 or 15 mg/day for the first 14 days; a single",
      "decrease to AOM 300 mg was allowed for tolerability (68 patients, 16.2",
      "percent). Phase D (double-blind, 52 weeks): randomised 1:1 to continue",
      "AOM 400 (or 300) mg every 4 weeks (n = 132) or switch to placebo",
      "injections (n = 133)."
    ),
    regions = "103 sites in the USA, Canada, Poland, Romania, Taiwan, South Korea and Japan (Appendix S1, Table S1)",
    notes = paste(
      "Event: recurrence of any mood episode during Phase D, defined by any of",
      "hospitalisation for a mood episode; YMRS >= 15, MADRS >= 15 and/or",
      "CGI-BP-S > 4; a serious adverse event of worsening disease;",
      "discontinuation for lack of efficacy or worsening disease; clinical",
      "worsening needing an added mood stabiliser, antidepressant,",
      "antipsychotic or above-limit benzodiazepine; or active suicidality",
      "(Wang 2025 Methods, 'Efficacy Data'). 103 events: 35 of 132 (26.5",
      "percent) on AOM and 68 of 133 (51.1 percent) on placebo; about 65",
      "percent of the events occurred in the first 120 days (Table 1, Figure 1).",
      "The proportional-hazards assumption was tested with Schoenfeld",
      "residuals and met (P = .8173). A log-linear CTROUGH term was also tried;",
      "the linear term was the more significant (P = .0170). Fitted with the R",
      "survival package (R 4.1.3)."
    )
  )

  ini({
    # ==================================================================
    # Wang 2025 final exposure-response model (Results, 'Exposure-
    # Response Model'; equation in Appendix S4 step 5):
    #
    #   S_i(t) = S0(t) ^ exp(beta * Ctau_i)
    #
    # equivalently h_i(t) = h0(t) * exp(beta * Ctau_i), with Ctau in
    # ng/mL, linear and uncentred. The coefficient is printed on the
    # log-hazard scale, which is the scale a Cox model estimates on, so
    # nothing here is back-solved or digitised.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # S0(t) is left completely unspecified by the method, so it is not an
    # unreported parameter but a quantity the fit never produced. This
    # model therefore returns the RELATIVE hazard hr only and
    # deliberately does not offer a survivor function, matching the
    # relative-hazard-only pattern of Rayner_2013_oseltamivir_
    # symptomalleviation.R and Chan_2025_atezolizumab_pfs.R.
    #
    # The regression was fitted with the R survival package rather than
    # NONMEM; the paper reports a coefficient, its RSE and a hazard ratio
    # with 95% CI, and no variance components, so there is no IIV and no
    # residual error to encode and no observation endpoint is declared.
    # ==================================================================
    e_ctrough_haz <- -0.00346; label("Log hazard ratio for recurrence of any mood episode per 1 ng/mL increase in aripiprazole Ctau (log scale)")  # Wang 2025 Results, 'Exposure-Response Model': coefficient for the effect of Ctau = -0.00346 (RSE 44%), HR 0.99655 (95% CI 0.9936-0.9995), P = .0230
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor and relative hazard against a hypothetical subject
    # with CTROUGH = 0. Under the proportional-hazards assumption hr is
    # constant in time, so the subject hazard is h(t) = h0(t) * hr for
    # any baseline h0(t) the user supplies; only RATIOS of hr between two
    # CTROUGH values are interpretable without a baseline.
    # ------------------------------------------------------------------
    lhr <- e_ctrough_haz * CTROUGH
    hr <- exp(lhr)
  })
}
