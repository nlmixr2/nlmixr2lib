DAgate_2021_bph_aurs_mbma <- function() {
  description <- paste(
    "MBMA (individual patient data). Parametric time-to-event model for the",
    "time to first acute urinary retention or benign prostatic hyperplasia",
    "(BPH)-related surgery (AUR/S) in men with moderate or severe lower",
    "urinary tract symptoms due to BPH, pooled from six phase III/IV",
    "dutasteride trials (10,238 patients, up to 4 years). The hazard is",
    "exponential (constant in time): h = lambda * exp(sum of beta_k *",
    "(x_k - median_k)) * HR_treatment, with baseline hazard lambda =",
    "7.78e-5 /day (2.84 percent per year) for a placebo patient at the pooled",
    "median baseline covariates. Four baseline covariates act log-linearly",
    "on the hazard, centred at their pooled medians: International Prostate",
    "Symptom Score (HR 1.04 per point, median 16), prostate-specific antigen",
    "(HR 1.08 per ng/mL, median 3.4), prostate volume (HR 1.01 per mL, median",
    "48.5) and maximum urinary flow rate Qmax (HR 0.91 per mL/s, median 10.2).",
    "Dutasteride monotherapy (HR 0.432) and tamsulosin-dutasteride combination",
    "therapy (HR 0.336) reduce the hazard; tamsulosin monotherapy was not",
    "different from placebo and its hazard ratio is fixed to 1. The model",
    "exposes the hazard `hazard` (1/day), the cumulative hazard `cumhaz` and",
    "the AUR/S-free survival probability `sur`. Sister IPSS drug-disease",
    "model fitted to the same pooled data: DAgate_2020_bph_ipss_mbma. Time",
    "in days.",
    sep = " "
  )
  reference <- paste(
    "D'Agate S, Chavan C, Manyak M, Palacios-Moreno JM, Oelke M, Michel MC,",
    "Roehrborn CG, Della Pasqua O.",
    "Model-based meta-analysis of the time to first acute urinary retention",
    "or benign prostatic hyperplasia-related surgery in patients with",
    "moderate or severe symptoms.",
    "Br J Clin Pharmacol. 2021;87:2777-2789. doi:10.1111/bcp.14682.",
    "Hazard form from Methods 2.3 (Equations 2 and 3, the exponential",
    "density lambda * exp(-lambda * t)); covariate centring at the median",
    "from Results 3.3; parameter values from Table 3; centring medians from",
    "Table 2. The Supporting Information (BCP-87-2777-s001) contains the",
    "survival-function definitions and Figures S1-S3 only; the control",
    "stream the Methods refer to is not part of the deposit.",
    sep = " "
  )
  vignette <- "DAgate_2021_bph_aurs_mbma"

  units <- list(
    time = "day",
    dosing = "(no dose events; each treatment is the single fixed regimen studied -- tamsulosin 0.4 mg, dutasteride 0.5 mg, or both, once daily -- and is identified by its arm indicator)",
    concentration = "(no concentration; the outputs are a hazard in 1/day, a unitless cumulative hazard and a unitless AUR/S-free survival probability)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "none",
      units = "(unitless cumulative hazard of first AUR/S)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    IPSS_BL = list(
      description = "Observed International Prostate Symptom Score at the baseline visit (last day of the placebo run-in).",
      units = "IPSS points (0-35)",
      type = "continuous",
      reference_category = "16 points (pooled-population median, Table 2)",
      notes = paste(
        "Log-linear effect on the hazard, exp(log(1.04) * (IPSS_BL - 16)):",
        "each extra point raises the instantaneous risk of AUR/S by 4",
        "percent (Table 3, Results 3.3). Centring at the median is stated in",
        "Results 3.3 ('starting from the median value of the covariate",
        "factor'); the median is Table 2. Missing baseline values were",
        "imputed with the study-population median (Methods 2.2)."
      ),
      source_name = "B_VAR"
    ),
    PSA_BL = list(
      description = "Baseline serum prostate-specific antigen concentration.",
      units = "ng/mL",
      type = "continuous",
      reference_category = "3.4 ng/mL (pooled-population median, Table 2)",
      notes = paste(
        "Log-linear effect on the hazard, exp(log(1.08) * (PSA_BL - 3.4)):",
        "+8 percent instantaneous risk per ng/mL (Table 3). Measured,",
        "unadjusted PSA; dutasteride roughly halves PSA on treatment, but",
        "the model uses the pre-treatment value only."
      ),
      source_name = "B_PSA"
    ),
    PROSTATE_VOL_BL = list(
      description = "Baseline prostate volume (transrectal ultrasound).",
      units = "mL",
      type = "continuous",
      reference_category = "48.5 mL (pooled-population median, Table 2)",
      notes = paste(
        "Log-linear effect on the hazard,",
        "exp(log(1.01) * (PROSTATE_VOL_BL - 48.5)): +1 percent per mL",
        "(Table 3). The published Figures 4-6 imply an unrounded hazard",
        "ratio near 1.0094 per mL, inside the Table 3 95 percent CI",
        "1.007-1.012; the printed 1.01 is used (see the vignette). Missing",
        "values (< 4 percent) were imputed with the study median."
      ),
      source_name = "B_PV"
    ),
    QMAX_BL = list(
      description = "Baseline maximum urinary flow rate (Qmax) from uroflowmetry.",
      units = "mL/s",
      type = "continuous",
      reference_category = "10.2 mL/s (pooled-population median, Table 2)",
      notes = paste(
        "Log-linear effect on the hazard, exp(log(0.91) * (QMAX_BL - 10.2)):",
        "-9 percent instantaneous risk per mL/s (Table 3). 10.5 percent of",
        "baseline Qmax values were missing and imputed with the study median",
        "(Methods 2.2)."
      ),
      source_name = "Qmax"
    ),
    TRT_TAMSULOSIN = list(
      description = "1 = the patient was randomised to tamsulosin 0.4 mg once-daily monotherapy, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all TRT_ indicators are 0)",
      notes = "Table 3 hazard ratio 1 with no uncertainty: the tamsulosin effect was not statistically different from placebo and was fixed to 1 (Results 3.3). Carried with a fixed zero log-hazard ratio so the arm structure is explicit. Tamsulosin-monotherapy patients come from CombAT only.",
      source_name = "TAM"
    ),
    TRT_DUTASTERIDE = list(
      description = "1 = the patient was randomised to dutasteride 0.5 mg once-daily monotherapy, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all TRT_ indicators are 0)",
      notes = "Hazard ratio 0.432 versus placebo (Table 3; a 56.8 percent lower instantaneous risk of AUR/S at any time).",
      source_name = "DUT"
    ),
    TRT_TAMSULOSIN_DUTASTERIDE = list(
      description = "1 = the patient was randomised to tamsulosin 0.4 mg + dutasteride 0.5 mg once-daily combination therapy, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all TRT_ indicators are 0)",
      notes = "Hazard ratio 0.336 versus placebo (Table 3; a 66.4 percent lower instantaneous risk). The combination has its own estimated hazard ratio; it is not the product of the monotherapy ratios.",
      source_name = "CT"
    )
  )

  covariatesDataExcluded <- list(
    TRT_WATCHFUL_WAITING = list(
      description = "1 = the patient was randomised to watchful waiting with protocol-defined initiation of tamsulosin (CONDUCT), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      notes = "Watchful waiting was handled as an active arm (Methods 2.3), but 'reliable estimates of the effect of WW could not be obtained due to the small sample size and very low number of events' (Results 3.3), and Table 3 has no watchful-waiting row. A watchful-waiting patient therefore carries the placebo hazard in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 9832L,
    n_studies = 6L,
    age_range = "47-94 years (median 66, mean 66.2)",
    weight_range = "37-179 kg (median 82, mean 83.2)",
    sex_female_pct = 0,
    race_ethnicity = "white 9268, black 229, Hispanic 276, Asian 374 (Table 2)",
    disease_state = paste(
      "Men with moderate or severe LUTS due to BPH at risk of progression;",
      "baseline IPSS median 16 (range 1-35), PSA median 3.4 ng/mL (0.6-23.2),",
      "prostate volume median 48.5 mL (16.6-296.9), Qmax median 10.2 mL/s",
      "(2.2-36.2) (Table 2)."
    ),
    dose_range = paste(
      "Placebo (n = 2158, 2 y); watchful waiting with protocol-defined",
      "tamsulosin initiation (n = 373); tamsulosin 0.4 mg once daily",
      "(n = 1611); dutasteride 0.5 mg once daily (n = 3790); tamsulosin +",
      "dutasteride combination (n = 2143); up to 4 years (Table 2)."
    ),
    regions = "multinational",
    notes = paste(
      "Studies ARIA3001, ARIA3002, ARI40002, CombAT, CONDUCT and ARIB3003",
      "(Table 1). Of 10,238 pooled patients, 402 (229 tamsulosin step-up",
      "patients in CONDUCT and 173 ARI40002 patients who switched from",
      "combination therapy to dutasteride or did not complete the second",
      "phase) were excluded, giving the 9832-patient final data set",
      "(Methods 2.2, Figure 1). NONMEM 7.3, Laplacian estimation. No",
      "variance components are estimated in a time-to-event model (Methods",
      "2.3)."
    )
  )

  ini({
    # Baseline hazard. Table 3 prints 'Baseline hazard rate lambda (per
    # 100 000 subjects) [day-1]' = 7.78, i.e. 7.78e-5 events per subject per
    # day; 7.78e-5 * 365.25 = 2.84 percent per year, the incidence Results 3.3
    # quotes for the baseline hazard.
    llam_haz <- log(7.78e-5); label("Baseline hazard of first AUR/S for a placebo patient at the median covariates (1/day)") # Table 3: 7.78 per 100 000 subjects per day (RSE 6%, 95% CI 6.86-8.70)

    # Covariate and treatment effects, entered as log hazard ratios. Table 3
    # prints the hazard ratios; Methods 2.3 Equation 2 defines HR = exp(beta)
    # per one-unit increase and Equation 3 h_j(t) = exp(alpha_j) * h0(t).
    e_ipss_bl_haz <- log(1.04); label("Log hazard ratio per IPSS point above the median of 16 (unitless)") # Table 3: baseline IPSS HR 1.04 (RSE 19%, 95% CI 1.022-1.050)
    e_psa_bl_haz <- log(1.08); label("Log hazard ratio per ng/mL of baseline PSA above the median of 3.4 ng/mL (unitless)") # Table 3: baseline PSA HR 1.08 (RSE 24%, 95% CI 1.041-1.120)
    e_prostate_vol_bl_haz <- log(1.01); label("Log hazard ratio per mL of baseline prostate volume above the median of 48.5 mL (unitless)") # Table 3: baseline prostate volume HR 1.01 (RSE 14%, 95% CI 1.007-1.012)
    e_qmax_bl_haz <- log(0.91); label("Log hazard ratio per mL/s of baseline maximum urinary flow above the median of 10.2 mL/s (unitless)") # Table 3: baseline maximum urinary flow HR 0.91 (RSE 14%, 95% CI 0.889-0.935)
    e_trt_tamsulosin_haz <- fixed(log(1)); label("Log hazard ratio for tamsulosin monotherapy versus placebo (unitless)") # Table 3: tamsulosin HR 1 with no RSE or CI; Results 3.3: not significant, 'therefore fixed to 1'
    e_trt_dutasteride_haz <- log(0.432); label("Log hazard ratio for dutasteride monotherapy versus placebo (unitless)") # Table 3: dutasteride HR 0.432 (RSE 7%, 95% CI 0.352-0.512)
    e_trt_tamsulosin_dutasteride_haz <- log(0.336); label("Log hazard ratio for tamsulosin-dutasteride combination therapy versus placebo (unitless)") # Table 3: combination HR 0.336 (RSE 7%, 95% CI 0.249-0.423)

    # No between-subject variability: Methods 2.3, 'no variance component
    # (interindividual, interoccasion or residual variability) is obtained
    # for model parameters'. The likelihood of a time-to-event model is the
    # event density itself, so there is no residual error either. The tiny
    # additive residual below is a placeholder on the survival-probability
    # output so the nlmixr2 likelihood machinery accepts the model for
    # forward simulation; it is not from the source.
    addSd <- fixed(0.001); label("Placeholder additive residual on the survival probability sur (unitless); not from the source")
  })

  model({
    # Log-linear covariate model centred at the pooled medians of Table 2
    # (Results 3.3: 'starting from the median value of the covariate factor
    # a difference of 1 unit in the covariate value will correspond to a
    # percentage increase or reduction in the baseline hazard').
    lp_cov <- e_ipss_bl_haz * (IPSS_BL - 16) +
      e_psa_bl_haz * (PSA_BL - 3.4) +
      e_prostate_vol_bl_haz * (PROSTATE_VOL_BL - 48.5) +
      e_qmax_bl_haz * (QMAX_BL - 10.2)

    # Treatment enters proportionally on the baseline hazard (Equation 3).
    # The arm indicators are mutually exclusive; all zero = placebo (and
    # watchful waiting, which has no estimated effect).
    lp_trt <- e_trt_tamsulosin_haz * TRT_TAMSULOSIN +
      e_trt_dutasteride_haz * TRT_DUTASTERIDE +
      e_trt_tamsulosin_dutasteride_haz * TRT_TAMSULOSIN_DUTASTERIDE

    # Exponential (time-constant) hazard, Methods 2.3 density lambda *
    # exp(-lambda * t), and the survival function S(t) = exp(-H(t)).
    hazard <- exp(llam_haz + lp_cov + lp_trt)
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
