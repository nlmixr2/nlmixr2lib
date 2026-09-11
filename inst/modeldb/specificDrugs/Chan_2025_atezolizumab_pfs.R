Chan_2025_atezolizumab_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for ",
    "PROGRESSION-FREE SURVIVAL in adults with locally advanced or ",
    "metastatic non-small cell lung cancer receiving SUBCUTANEOUS ",
    "atezolizumab 1875 mg every 3 weeks (Chan 2025, n = 246, cohort 5 of ",
    "the phase III portion of IMscin001, NCT03735121, clinical cut-off ",
    "26 April 2022). The model returns the RELATIVE hazard ",
    "exp(-0.002 * (CTROUGH / 10) - 0.025 * AGE - 0.114 * (WT / 10) - ",
    "0.062 * ALB + 0.532 * LMET + 0.046 * NLR + 0.670 * HEPIMP_MILD + ",
    "0.409 * HEPIMP_MOD). THE EXPOSURE TERM IS NOT STATISTICALLY ",
    "SIGNIFICANT and is essentially zero (-0.002 per 10 ug/mL, 95 ",
    "percent CI -0.052 to 0.047); that flat exposure-response is the ",
    "paper's headline result, not a transcription gap. Six of the seven ",
    "baseline covariates ARE significant, so the model is in practice a ",
    "baseline-prognostic model for advanced NSCLC with a null exposure ",
    "term carried alongside it. NO BASELINE HAZARD IS ENCODED: a Cox ",
    "regression is semiparametric and leaves h0(t) unspecified, so the ",
    "absolute survivor function is a quantity the fit never produced ",
    "rather than an unreported parameter; only RATIOS of hr between two ",
    "covariate settings are meaningful. There is no PK layer and no ODE: ",
    "the exposure metric is supplied as a data column, derived in the ",
    "source analysis from the individual empirical-Bayes predictions of ",
    "the companion population PK model packaged as ",
    "Chan_2025_atezolizumab. No between-subject random effect and no ",
    "residual error are estimated. Five companion exposure-response ",
    "models in the Chan_2025_atezolizumab_* family."
  )
  reference <- paste(
    "Chan P, Liu SN, Gosselin N, Sauve Z, Marchand M, Lin A,",
    "Herraez-Baranda L, Zanghi J, Shearer-Kang E, Liu X, Wu B, Chanu P.",
    "Population pharmacokinetics and exposure-response of subcutaneous",
    "atezolizumab in patients with non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):726-737.",
    "doi:10.1002/psp4.13310.",
    "Coefficients are transcribed from Supporting Information Table S10B.",
    "Individual Cycle-1 Ctrough values derive from the companion",
    "population PK model reported in the same paper; see",
    "modellib('Chan_2025_atezolizumab').",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(
    time          = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH data column)",
    concentration = "hr (relative hazard of progression or death, unitless; also lhr, the linear predictor)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Model-predicted Cycle-1 trough serum atezolizumab concentration",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 10, not centred: Chan 2025",
        "Table S10B labels the row 'C_trough (10 ug/mL)', so the",
        "coefficient -0.002 is the log-hazard-ratio per 10 ug/mL.",
        "Derived from the companion population PK model, NOT observed.",
        "Chan 2025 Methods selected 'the exposure metric with the lowest",
        "p-value' per endpoint from AUC0-21d and Ctrough for the efficacy",
        "endpoints, and Ctrough won for progression-free survival",
        "(AUC0-21d won for overall survival). The retained coefficient is",
        "nonetheless effectively zero with a confidence interval",
        "straddling it symmetrically. Cohort 5 distribution (Chan 2025",
        "Table S5B): geometric mean 97.2 ug/mL (geoCV 42.3 percent),",
        "median 99.2, range 20.0-218."
      ),
      source_name        = "Ctrough"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY per year and is not centred: Chan 2025",
        "Table S10B labels the row 'Age (year)'. Coefficient -0.025",
        "(95 percent CI -0.043 to -0.008, p < 0.01), i.e. older patients",
        "had a LOWER hazard of progression or death in this cohort."
      ),
      source_name        = "Age"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 10, not centred: Chan 2025",
        "Table S10B labels the row 'Body weight (10 kg)', so the",
        "coefficient -0.114 is the log-hazard-ratio per 10 kg (95 percent",
        "CI -0.223 to -0.006, p < 0.05). Cohort 5 median 67.8 kg",
        "[30.0, 117] (Chan 2025 Table 1). Note the direction: heavier",
        "patients had a lower hazard, the opposite of what a",
        "weight-driven-exposure argument would predict, which is",
        "consistent with weight acting here as a proxy for nutritional",
        "and performance status rather than through drug exposure."
      ),
      source_name        = "Body weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY per g/L and is not centred: Chan 2025",
        "Table S10B labels the row 'Albumin (g/L)'. Coefficient -0.062",
        "(95 percent CI -0.096 to -0.027, p < 0.001), the strongest",
        "single predictor in the model; higher albumin means a lower",
        "hazard, the well-established prognostic direction in advanced",
        "NSCLC. Reported in g/L (SI), consistent with Chan 2025 Table 1",
        "(cohort 5 median 40.0 g/L) and with the population PK model's",
        "40 g/L covariate reference. NOTE that the companion overall",
        "survival table (Table S10C) labels the same covariate 'g/dL';",
        "that label is a typographical slip -- both Cox models were",
        "fitted on the same merged analysis dataset, whose albumin column",
        "is in g/L, and the two coefficients (-0.062 here, -0.082 there)",
        "are of the same order, which they could not be if they were on",
        "scales differing by a factor of ten. See the vignette Errata."
      ),
      source_name        = "Albumin"
    ),
    LMET = list(
      description        = "Baseline presence of liver metastases (1 = present, 0 = absent)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no liver metastases at baseline)",
      notes              = paste(
        "Chan 2025 Table S10B labels the row 'Presence of metastasis in",
        "the liver'. Coefficient 0.532 (95 percent CI 0.189 to 0.875,",
        "p < 0.01), i.e. a hazard ratio of exp(0.532) = 1.70 for patients",
        "with liver metastases. Chan 2025 Methods lists 'liver metastasis'",
        "among the baseline characteristics merged into the",
        "exposure-response dataset."
      ),
      source_name        = "Presence of metastasis in the liver"
    ),
    NLR = list(
      description        = "Baseline neutrophil-to-lymphocyte ratio",
      units              = "ratio (unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY per unit ratio and is not centred: Chan 2025",
        "Table S10B labels the row 'Neutrophil-to-lymphocyte ratio' with",
        "no rescaling parenthetical, unlike the '(10 kg)' and",
        "'(10 ug/mL)' rows in the same table. Coefficient 0.046",
        "(95 percent CI 0.026 to 0.065, p < 0.001); a higher ratio means",
        "a higher hazard, the established inflammatory-prognostic",
        "direction."
      ),
      source_name        = "Neutrophil-to-lymphocyte ratio"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (1 = mild, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; mutually exclusive with HEPIMP_MOD)",
      notes              = paste(
        "One of two non-reference levels of a three-level hepatic",
        "function factor (normal [reference], mild, moderate).",
        "Coefficient 0.670 (95 percent CI 0.196 to 1.145, p < 0.01), a",
        "hazard ratio of 1.95. Chan 2025 Methods cites a standard",
        "hepatic-impairment classification (reference 21 of the paper)",
        "for the categorisation."
      ),
      source_name        = "Hepatic impairment-mild"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = moderate, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; mutually exclusive with HEPIMP_MILD)",
      notes              = paste(
        "The second non-reference level of the hepatic function factor.",
        "Coefficient 0.409 with a 95 percent CI of -1.939 to 2.038, i.e.",
        "wholly uninformative -- the moderate-impairment stratum is very",
        "small. It is retained because dropping it would silently move",
        "those patients into the normal-function reference and bias the",
        "mild-impairment coefficient a downstream user reproduces.",
        "Interpret the point estimate with caution."
      ),
      source_name        = "Hepatic impairment-moderate"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 246,
    n_studies      = 1,
    disease_state  = "locally advanced or metastatic non-small cell lung cancer",
    dose_range     = "subcutaneous atezolizumab 1875 mg every 3 weeks in the thigh, as a ready-to-use co-formulation with recombinant human hyaluronidase PH20",
    notes          = paste(
      "Cohort 5 of the phase III (part 2) portion of IMscin001",
      "(NCT03735121). Exposure-response was assessed on the subcutaneous",
      "arm ONLY -- Chan 2025 Methods notes that intravenous",
      "exposure-response has been evaluated extensively elsewhere and was",
      "therefore excluded. Time-to-event endpoints 'were first explored",
      "using Kaplan-Meier estimation stratified by quartiles of exposure",
      "metrics, then by using Cox regression with exposure metrics as",
      "continuous variables'. Baseline medians (Chan 2025 Table 1): body",
      "weight 67.8 kg, tumor burden 79.5 mm, albumin 40.0 g/L, hemoglobin",
      "123 g/L; 29.3 percent female. Cycle-1 exposure metrics were used",
      "deliberately 'to minimize the potential effect of",
      "response-dependent time-varying clearance'."
    )
  )

  ini({
    # ==================================================================
    # Chan 2025 Supporting Information Table S10B, "Cox Proportional
    # Hazard Model for PFS and Atezolizumab Exposure with Relevant
    # Covariates in Cohort 5 (Atezolizumab SC 1875 mg Q3W)". Every entry
    # is printed as a Covariate Effect with a 95 percent CI on the
    # log-hazard scale, which is the scale a Cox model estimates on, so
    # nothing here is back-solved or digitised.
    #
    #   Ctrough (10 ug/mL)                  -0.002  (-0.052,  0.047)
    #   Age (year)                          -0.025  (-0.043, -0.008)  **
    #   Body weight (10 kg)                 -0.114  (-0.223, -0.006)  *
    #   Albumin (g/L)                       -0.062  (-0.096, -0.027)  ***
    #   Presence of metastasis in the liver  0.532  ( 0.189,  0.875)  **
    #   Neutrophil-to-lymphocyte ratio       0.046  ( 0.026,  0.065)  ***
    #   Hepatic impairment-mild              0.670  ( 0.196,  1.145)  **
    #   Hepatic impairment-moderate          0.409  (-1.939,  2.038)
    #     (* p < 0.05, ** p < 0.01, *** p < 0.001)
    #
    # Table S10B's own footnote fixes the sign convention: "positive
    # covariate effect corresponds to an increase in the risk for disease
    # progression or death; negative covariate effect corresponds to a
    # decrease".
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is left completely unspecified by the method, so it is not an
    # unreported parameter but a quantity the fit never produced. This
    # model therefore returns the RELATIVE hazard hr only and
    # deliberately does not offer a survivor function, matching the
    # relative-hazard-only pattern of Liu_2024_saf189s_pfs.R and
    # Rayner_2013_oseltamivir_shedding.R.
    #
    # The regression was fitted in R 4.1.1 rather than NONMEM; the paper
    # reports coefficients and confidence intervals with no variance
    # components, so there is no IIV and no residual error to encode and
    # no observation endpoint is declared.
    # ==================================================================
    e_ctrough_haz <- -0.002; label("Log hazard ratio for progression or death per 10 ug/mL increase in Cycle-1 trough atezolizumab concentration (log scale)")  # Chan 2025 Table S10B, Ctrough row: -0.002 (95% CI -0.052 to 0.047), not significant
    e_age_haz <- -0.025; label("Log hazard ratio for progression or death per year of baseline age (log scale)")                                                # Chan 2025 Table S10B, Age row: -0.025 (95% CI -0.043 to -0.008), p < 0.01
    e_wt_haz <- -0.114; label("Log hazard ratio for progression or death per 10 kg of baseline body weight (log scale)")                                        # Chan 2025 Table S10B, Body weight row: -0.114 (95% CI -0.223 to -0.006), p < 0.05
    e_alb_haz <- -0.062; label("Log hazard ratio for progression or death per g/L of baseline serum albumin (log scale)")                                       # Chan 2025 Table S10B, Albumin row: -0.062 (95% CI -0.096 to -0.027), p < 0.001
    e_lmet_haz <- 0.532; label("Log hazard ratio for progression or death for baseline liver metastases versus none (log scale; HR 1.70)")                      # Chan 2025 Table S10B, liver-metastasis row: 0.532 (95% CI 0.189 to 0.875), p < 0.01
    e_nlr_haz <- 0.046; label("Log hazard ratio for progression or death per unit of baseline neutrophil-to-lymphocyte ratio (log scale)")                      # Chan 2025 Table S10B, NLR row: 0.046 (95% CI 0.026 to 0.065), p < 0.001
    e_hepimp_mild_haz <- 0.670; label("Log hazard ratio for progression or death for mild hepatic impairment versus normal (log scale; HR 1.95)")               # Chan 2025 Table S10B, Hepatic impairment-mild row: 0.670 (95% CI 0.196 to 1.145), p < 0.01
    e_hepimp_mod_haz <- 0.409; label("Log hazard ratio for progression or death for moderate hepatic impairment versus normal (log scale; CI spans zero)")      # Chan 2025 Table S10B, Hepatic impairment-moderate row: 0.409 (95% CI -1.939 to 2.038), not significant
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor and relative hazard against a hypothetical patient
    # with every covariate at zero. Under the proportional-hazards
    # assumption hr is constant in time, so the subject hazard is
    # h(t) = h0(t) * hr for any baseline h0(t) the user supplies. Only
    # RATIOS of hr between two covariate settings are meaningful, since
    # the all-zero reference (age 0, weight 0, albumin 0) is not a
    # patient.
    #
    # The divisors 10 convert the raw canonical data columns (ug/mL for
    # CTROUGH, kg for WT) into the per-10 units in which Table S10B
    # reports those two coefficients. AGE, ALB and NLR are per-unit in
    # the source and so are used undivided. HEPIMP_MILD and HEPIMP_MOD
    # are mutually exclusive; setting both to 0 selects normal hepatic
    # function.
    # ------------------------------------------------------------------
    lhr <- e_ctrough_haz * (CTROUGH / 10) +
      e_age_haz * AGE +
      e_wt_haz * (WT / 10) +
      e_alb_haz * ALB +
      e_lmet_haz * LMET +
      e_nlr_haz * NLR +
      e_hepimp_mild_haz * HEPIMP_MILD +
      e_hepimp_mod_haz * HEPIMP_MOD

    hr <- exp(lhr)
  })
}
