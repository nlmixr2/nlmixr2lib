Chan_2025_atezolizumab_os <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for OVERALL ",
    "SURVIVAL in adults with locally advanced or metastatic non-small ",
    "cell lung cancer receiving SUBCUTANEOUS atezolizumab 1875 mg every ",
    "3 weeks (Chan 2025, n = 246, cohort 5 of the phase III portion of ",
    "IMscin001, NCT03735121, clinical cut-off 26 April 2022). The model ",
    "returns the RELATIVE hazard ",
    "exp(0.134 * (AUC_ATEZO / 1000) - 0.082 * ALB + 0.124 * (LDH / 100) ",
    "+ 0.081 * (CRP / 10) + 0.040 * NLR + 0.934 * HEPIMP_MILD + 0.414 * ",
    "HEPIMP_MOD). THE EXPOSURE TERM IS NOT STATISTICALLY SIGNIFICANT ",
    "(0.134 per 1000 ug*day/mL, 95 percent CI -0.110 to 0.377); that ",
    "flat exposure-response is the paper's headline result, not a ",
    "transcription gap. Note that this endpoint retained AUC0-21d ",
    "whereas the companion progression-free-survival model retained ",
    "Ctrough: Chan 2025 selected, per endpoint, whichever of the two ",
    "efficacy exposure metrics gave the lowest p-value. All six baseline ",
    "covariates that carry a significance marker ARE significant, so the ",
    "model is in practice a baseline-prognostic model for advanced NSCLC ",
    "with a null exposure term carried alongside it. NO BASELINE HAZARD ",
    "IS ENCODED: a Cox regression is semiparametric and leaves h0(t) ",
    "unspecified, so the absolute survivor function is a quantity the ",
    "fit never produced rather than an unreported parameter; only RATIOS ",
    "of hr between two covariate settings are meaningful. There is no PK ",
    "layer and no ODE: the exposure metric is supplied as a data column, ",
    "derived in the source analysis from the individual empirical-Bayes ",
    "predictions of the companion population PK model packaged as ",
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
    "Coefficients are transcribed from Supporting Information Table S10C.",
    "Individual Cycle-1 AUC0-21d values derive from the companion",
    "population PK model reported in the same paper; see",
    "modellib('Chan_2025_atezolizumab').",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(
    time          = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_ATEZO data column)",
    concentration = "hr (relative hazard of death, unitless; also lhr, the linear predictor)"
  )

  covariateData <- list(
    AUC_ATEZO = list(
      description        = "Model-predicted Cycle-1 atezolizumab AUC from day 0 to day 21",
      units              = "ug*day/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 1000, not centred: Chan 2025",
        "Table S10C labels the row 'AUC0-21d (1000 ug*day/mL)', so the",
        "coefficient 0.134 is the log-hazard-ratio per 1000 ug*day/mL.",
        "Derived from the companion population PK model, NOT observed.",
        "Chan 2025 Methods selected 'the exposure metric with the lowest",
        "p-value' per endpoint from AUC0-21d and Ctrough for the efficacy",
        "endpoints, and AUC0-21d won for overall survival (Ctrough won",
        "for progression-free survival). The confidence interval spans",
        "zero. Cohort 5 distribution (Chan 2025 Table S5B): geometric",
        "mean 2907 ug*day/mL (geoCV 35.9 percent), median 2974, range",
        "666-6572."
      ),
      source_name        = "AUC0-21d"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY per g/L and is not centred. Coefficient -0.082",
        "(95 percent CI -0.135 to -0.030, p < 0.01); higher albumin means",
        "a lower hazard of death, the well-established prognostic",
        "direction in advanced NSCLC. UNITS ERRATUM: Chan 2025",
        "Table S10C labels this row 'Albumin (g/dL)', but the value is",
        "encoded here per g/L. Three independent lines of evidence fix",
        "the scale as g/L. (1) The companion progression-free-survival",
        "table (Table S10B) labels the identical covariate from the",
        "identical merged analysis dataset 'Albumin (g/L)'. (2) That",
        "dataset's albumin column is in g/L throughout the paper: Chan",
        "2025 Table 1 reports cohort medians of 39.0-41.2 g/L and the",
        "population PK model normalizes albumin to a 40 g/L reference.",
        "(3) The two Cox coefficients, -0.062 (PFS) and -0.082 (OS), are",
        "of the same order; they could not be if they were on scales",
        "differing by a factor of ten. The sibling rows in Table S10C are",
        "also SI throughout ('LDH (100 U/L)', 'CRP (10 mg/L)'), so the",
        "'g/dL' label is an isolated typographical slip. See the vignette",
        "Errata."
      ),
      source_name        = "Albumin"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 100, not centred: Chan 2025",
        "Table S10C labels the row 'LDH (100 U/L)', so the coefficient",
        "0.124 is the log-hazard-ratio per 100 U/L (95 percent CI 0.040",
        "to 0.209, p < 0.01). Higher lactate dehydrogenase means a higher",
        "hazard, the established prognostic direction. Chan 2025 Methods",
        "lists lactate dehydrogenase among the baseline characteristics",
        "merged into the exposure-response dataset. Retained for overall",
        "survival but not for progression-free survival."
      ),
      source_name        = "LDH"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 10, not centred: Chan 2025",
        "Table S10C labels the row 'CRP (10 mg/L)', so the coefficient",
        "0.081 is the log-hazard-ratio per 10 mg/L (95 percent CI 0.030",
        "to 0.131, p < 0.01). The Chan 2025 Discussion singles CRP out:",
        "'baseline CRP level has a statistically significant impact on",
        "most of the efficacy endpoints'. It also appears in the",
        "companion SAE and AEG35 safety models. Reported in mg/L (SI)."
      ),
      source_name        = "CRP"
    ),
    NLR = list(
      description        = "Baseline neutrophil-to-lymphocyte ratio",
      units              = "ratio (unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY per unit ratio and is not centred: Chan 2025",
        "Table S10C labels the row 'Neutrophil-to-Lymphocyte ratio' with",
        "no rescaling parenthetical, unlike the '(100 U/L)' and",
        "'(10 mg/L)' rows in the same table. Coefficient 0.040",
        "(95 percent CI 0.016 to 0.064, p < 0.01); a higher ratio means a",
        "higher hazard, consistent with the companion",
        "progression-free-survival model's 0.046."
      ),
      source_name        = "Neutrophil-to-Lymphocyte ratio"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (1 = mild, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; mutually exclusive with HEPIMP_MOD)",
      notes              = paste(
        "One of two non-reference levels of a three-level hepatic",
        "function factor (normal [reference], mild, moderate).",
        "Coefficient 0.934 (95 percent CI 0.341 to 1.527, p < 0.01), a",
        "hazard ratio of 2.54 -- the largest effect in the model, and",
        "larger than the same covariate's 0.670 in the companion",
        "progression-free-survival model. Chan 2025 Methods cites a",
        "standard hepatic-impairment classification (reference 21 of the",
        "paper) for the categorisation."
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
        "Coefficient 0.414 with a 95 percent CI of -1.597 to 2.426, i.e.",
        "wholly uninformative -- the moderate-impairment stratum is very",
        "small, and the interval is nearly identical to the one in the",
        "companion progression-free-survival model. It is retained",
        "because dropping it would silently move those patients into the",
        "normal-function reference and bias the mild-impairment",
        "coefficient a downstream user reproduces. Interpret the point",
        "estimate with caution."
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
    # Chan 2025 Supporting Information Table S10C, "Cox Proportional
    # Hazard Model for OS and Atezolizumab Exposure with Relevant
    # Covariates in Cohort 5 (Atezolizumab SC 1875 mg Q3W)". Every entry
    # is printed as a Coefficient with a 95 percent CI on the log-hazard
    # scale, which is the scale a Cox model estimates on, so nothing here
    # is back-solved or digitised.
    #
    #   AUC0-21d (1000 ug*day/mL)       0.134  (-0.110, 0.377)
    #   Albumin (g/dL -> see erratum)  -0.082  (-0.135, -0.030)  **
    #   LDH (100 U/L)                   0.124  ( 0.040,  0.209)  **
    #   CRP (10 mg/L)                   0.081  ( 0.030,  0.131)  **
    #   Neutrophil-to-Lymphocyte ratio  0.040  ( 0.016,  0.064)  **
    #   Hepatic impairment-mild         0.934  ( 0.341,  1.527)  **
    #   Hepatic impairment-moderate     0.414  (-1.597,  2.426)
    #     (* p < 0.05, ** p < 0.01, *** p < 0.001)
    #
    # Table S10C's own footnote fixes the sign convention: "positive
    # covariate effect corresponds to an increase in the risk of death;
    # negative covariate effect corresponds to a decrease".
    #
    # The albumin row's printed unit is 'g/dL'; it is encoded here per
    # g/L. See the ALB covariateData note for the three-way argument and
    # the vignette Errata.
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
    e_auc_atezo_haz <- 0.134; label("Log hazard ratio for death per 1000 ug*day/mL increase in Cycle-1 AUC0-21d (log scale)")                   # Chan 2025 Table S10C, AUC0-21d row: 0.134 (95% CI -0.110 to 0.377), not significant
    e_alb_haz <- -0.082; label("Log hazard ratio for death per g/L of baseline serum albumin (log scale)")                                      # Chan 2025 Table S10C, Albumin row: -0.082 (95% CI -0.135 to -0.030), p < 0.01; printed unit "g/dL" is a typographical slip, see the ALB covariateData note
    e_ldh_haz <- 0.124; label("Log hazard ratio for death per 100 U/L of baseline lactate dehydrogenase (log scale)")                           # Chan 2025 Table S10C, LDH row: 0.124 (95% CI 0.040 to 0.209), p < 0.01
    e_crp_haz <- 0.081; label("Log hazard ratio for death per 10 mg/L of baseline C-reactive protein (log scale)")                              # Chan 2025 Table S10C, CRP row: 0.081 (95% CI 0.030 to 0.131), p < 0.01
    e_nlr_haz <- 0.040; label("Log hazard ratio for death per unit of baseline neutrophil-to-lymphocyte ratio (log scale)")                     # Chan 2025 Table S10C, NLR row: 0.040 (95% CI 0.016 to 0.064), p < 0.01
    e_hepimp_mild_haz <- 0.934; label("Log hazard ratio for death for mild hepatic impairment versus normal (log scale; HR 2.54)")              # Chan 2025 Table S10C, Hepatic impairment-mild row: 0.934 (95% CI 0.341 to 1.527), p < 0.01
    e_hepimp_mod_haz <- 0.414; label("Log hazard ratio for death for moderate hepatic impairment versus normal (log scale; CI spans zero)")     # Chan 2025 Table S10C, Hepatic impairment-moderate row: 0.414 (95% CI -1.597 to 2.426), not significant
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor and relative hazard against a hypothetical patient
    # with every covariate at zero. Under the proportional-hazards
    # assumption hr is constant in time, so the subject hazard is
    # h(t) = h0(t) * hr for any baseline h0(t) the user supplies. Only
    # RATIOS of hr between two covariate settings are meaningful, since
    # the all-zero reference (albumin 0, LDH 0) is not a patient.
    #
    # The divisors 1000, 100 and 10 convert the raw canonical data
    # columns (ug*day/mL, U/L, mg/L) into the per-1000, per-100 and
    # per-10 units in which Table S10C reports those coefficients. ALB
    # and NLR are per-unit in the source and so are used undivided.
    # HEPIMP_MILD and HEPIMP_MOD are mutually exclusive; setting both to
    # 0 selects normal hepatic function.
    # ------------------------------------------------------------------
    lhr <- e_auc_atezo_haz * (AUC_ATEZO / 1000) +
      e_alb_haz * ALB +
      e_ldh_haz * (LDH / 100) +
      e_crp_haz * (CRP / 10) +
      e_nlr_haz * NLR +
      e_hepimp_mild_haz * HEPIMP_MILD +
      e_hepimp_mod_haz * HEPIMP_MOD

    hr <- exp(lhr)
  })
}
