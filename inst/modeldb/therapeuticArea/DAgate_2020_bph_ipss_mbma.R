DAgate_2020_bph_ipss_mbma <- function() {
  description <- paste(
    "MBMA (individual patient data). Longitudinal drug-disease model of the",
    "International Prostate Symptom Score (IPSS) in men with moderate or",
    "severe lower urinary tract symptoms due to benign prostatic hyperplasia",
    "(BPH), pooled from six phase III/IV dutasteride trials (10,238 patients,",
    "up to 4 years). IPSS follows",
    "dIPSS/dt = DISP * (1 - IPSS/35) - (PLACEBO + TREATMENT) * IPSS with",
    "IPSS(0) = observed baseline IPSS: a zero-order symptom-worsening term",
    "(DISP, 0.347 points/month) that saturates at the 35-point scale maximum,",
    "a placebo effect that decays exponentially (maximum 0.061 /month,",
    "half-life 7.26 months), and a constant first-order treatment effect for",
    "tamsulosin (0.015 /month), dutasteride (0.016 /month), watchful waiting",
    "(0.018 /month) or tamsulosin-dutasteride combination therapy",
    "(0.032 /month). Covariates: baseline IPSS on DISP, duration of BPH",
    "symptoms on the placebo magnitude, alcohol use on the placebo half-life.",
    "As in the authors' control stream, the between-subject variability on",
    "DISP and on the placebo magnitude applies only to placebo-arm patients.",
    "Time in months."
  )
  reference <- paste(
    "D'Agate S, Wilson T, Adalig B, Manyak M, Palacios-Moreno JM, Chavan C,",
    "Oelke M, Roehrborn C, Della Pasqua O.",
    "Model-based meta-analysis of individual International Prostate Symptom",
    "Score trajectories in patients with benign prostatic hyperplasia with",
    "moderate or severe symptoms.",
    "Br J Clin Pharmacol. 2020;86(8):1585-1599. doi:10.1111/bcp.14268.",
    "Parameter values from Tables 2 and 3; structural equation, covariate",
    "functional forms and the placebo random-effect block from the NONMEM",
    "control stream in the Supporting Information appendix.",
    sep = " "
  )
  vignette <- "DAgate_2020_bph_ipss_mbma"

  units <- list(
    time = "month",
    dosing = "(no dose events; each treatment is the single fixed regimen studied -- tamsulosin 0.4 mg, dutasteride 0.5 mg, or both, once daily -- and is identified by its arm indicator)",
    concentration = "(IPSS points; the modelled quantity is a 0-35 symptom score, not a drug concentration, so the dosing-versus-concentration dimensional check is not applicable)"
  )

  compartmentData <- list(
    ipss = list(
      analyte = "none",
      units = "IPSS points (0-35 clinical score)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    IPSS_BL = list(
      description = "Observed International Prostate Symptom Score at the baseline (randomisation) visit.",
      units = "IPSS points (0-35)",
      type = "continuous",
      reference_category = "16 points (centring value of the covariate effect on DISP)",
      notes = paste(
        "Used twice: as the initial condition ipss(0) = IPSS_BL",
        "(control stream A_0(1) = obsIPSS) and as a linear covariate on the",
        "disease progression rate, DISP * (1 + 0.027 * (IPSS_BL - 16)).",
        "Screening IPSS is a different quantity and was excluded from the",
        "analysis (Discussion 4.1). The control stream imputes a missing",
        "baseline (-99) as 16, the median."
      ),
      source_name = "IPSS0"
    ),
    T_SYMPT_BPH = list(
      description = "Duration of BPH symptoms at baseline (time since symptom onset).",
      units = "years",
      type = "continuous",
      reference_category = "4 years (centring value; Table 1 median 4 y)",
      notes = paste(
        "Linear effect on the magnitude of the placebo effect,",
        "DELTA_placebo * (1 - 0.025 * (T_SYMPT_BPH - 4)). Distinct from the",
        "time since BPH DIAGNOSIS (BPHTIM, Table 1 median 2.3 y), which was",
        "screened but not retained. Table 1 reports the duration in years",
        "(median 4, range 0-54.8) and the control stream centres it at 4.00,",
        "so the Supporting Information 'months' label for this column is",
        "taken as a slip. A missing value (-99) sets the factor to 1."
      ),
      source_name = "BPHDUR"
    ),
    ALCOHOL_USE = list(
      description = "1 = alcohol user at baseline, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not an alcohol user)",
      notes = paste(
        "Linear categorical effect on the placebo half-life,",
        "T1/2 * (1 - 0.135 * ALCOHOL_USE). The control stream codes ALCU as",
        "0, 1 or 2 and applies the effect to both 1 and 2, so any non-zero",
        "ALCU maps to ALCOHOL_USE = 1; the meaning of the 1 / 2 split is not",
        "published. Table 1: 6198 users / 3992 non-users."
      ),
      source_name = "ALCU"
    ),
    TRT_TAMSULOSIN = list(
      description = "1 = the patient was randomised to tamsulosin 0.4 mg once-daily monotherapy, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all four TRT_ indicators are 0)",
      notes = "Control stream ARM = 3. Selects DELTA_tamsulosin and its IIV and the tamsulosin residual SD. The four TRT_ indicators are mutually exclusive; all zero = placebo (ARM = 0).",
      source_name = "ARM"
    ),
    TRT_DUTASTERIDE = list(
      description = "1 = the patient was randomised to dutasteride 0.5 mg once-daily monotherapy, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all four TRT_ indicators are 0)",
      notes = "Control stream ARM = 4. Selects DELTA_dutasteride and its IIV and the dutasteride residual SD.",
      source_name = "ARM"
    ),
    TRT_TAMSULOSIN_DUTASTERIDE = list(
      description = "1 = the patient was randomised to tamsulosin 0.4 mg + dutasteride 0.5 mg once-daily combination therapy (free or fixed-dose combination), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all four TRT_ indicators are 0)",
      notes = "Control stream ARM = 5 ('FDC'). The combination carries its own estimated effect DELTA_FDC; it is NOT the sum of the two monotherapy effects (the Supporting Information reports that additive-combination parameterisations were tested and rejected).",
      source_name = "ARM"
    ),
    TRT_WATCHFUL_WAITING = list(
      description = "1 = the patient was randomised to watchful waiting with protocol-defined initiation of tamsulosin (CONDUCT), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo when all four TRT_ indicators are 0)",
      notes = "Control stream ARM = 1 or 2 (the two are pooled). The arm's effect DELTA_WW is a constant from time 0 regardless of when tamsulosin was stepped up; the step-up time is not a model input.",
      source_name = "ARM"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index at baseline.",
      units = "kg/m^2",
      type = "continuous",
      notes = "The control stream carries BMI effects on DISP and on the placebo half-life, centred at 26.89 kg/m^2, but both THETAs are '(0) FIX', so BMI has no effect in the final model and is absent from Table 2."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10238L,
    n_studies = 6L,
    age_range = "47-94 years (median 66, mean 66.2)",
    weight_range = "37-179 kg (median 82, mean 83.2)",
    sex_female_pct = 0,
    race_ethnicity = "white 9268, black 229, Hispanic 276, Asian 374 (Table 1)",
    disease_state = paste(
      "Men with moderate or severe LUTS due to BPH; baseline IPSS median 16",
      "(range 1-35), prostate volume median 48.6 cm^3, PSA median 3.4 ng/mL,",
      "BPH symptom duration median 4 y."
    ),
    dose_range = paste(
      "Placebo (n = 2158, 2 y); watchful waiting with protocol-defined",
      "tamsulosin initiation (n = 373); tamsulosin 0.4 mg once daily",
      "(n = 1611); dutasteride 0.5 mg once daily (n = 3790); tamsulosin +",
      "dutasteride combination (n = 2143); up to 4 years."
    ),
    regions = "multinational",
    notes = paste(
      "Studies ARIA3001, ARIA3002, ARI40002, CombAT, CONDUCT and ARIB3003;",
      "140,733 IPSS observations (Methods 2.1, Table 1, Figure 1). NONMEM 7.3",
      "FOCE-I. Placebo and disease-progression parameters (Table 2) were",
      "estimated first and then fixed while the treatment effects (Table 3)",
      "were estimated."
    )
  )

  ini({
    # ================================================================
    # Disease progression and placebo effect (Table 2). The authors
    # estimated these in a first step and then held them fixed while
    # estimating the treatment effects; in the published final control
    # stream every one of these THETAs carries FIX, so fixed() throughout.
    # The control stream runs in days (THETA(1) 0.0116 /day, THETA(2)
    # 0.00202 /day, THETA(3) 218 days); Table 2 reports the same values
    # per 30-day month, which is the time unit used here.
    # ================================================================
    ldisp <- fixed(log(0.347))
    label("Disease (symptom) progression rate DISP at IPSS_BL = 16 (IPSS points/month)")
    # Table 2 'DISP, disease progression rate' 0.347 (printed 'months-1'); control stream THETA(1) 0.0116 /day x 30 = 0.348.

    ldelta_pbo <- fixed(log(0.061))
    label("Maximum placebo effect DELTA_placebo at 4 y symptom duration (1/month)")
    # Table 2 'DELTA placebo, placebo effect' 0.061 months-1; control stream THETA(2) 0.00202 /day x 30 = 0.0606.

    lthalf_pbo <- fixed(log(7.263))
    label("Half-life of the placebo effect in non-drinkers (month)")
    # Table 2 'T1/2, placebo t1/2' 7.263 months; control stream THETA(3) 218 days / 30 = 7.27.

    e_t_sympt_bph_delta_pbo <- fixed(-0.025)
    label("Linear effect of BPH symptom duration on the placebo magnitude (1/year)")
    # Table 2 'BPHDUR on magnitude of placebo effect' -0.025; control stream THETA(4) -0.0245, centred at 4.00.

    e_ipss_bl_disp <- fixed(0.027)
    label("Linear effect of baseline IPSS on DISP (1/IPSS point)")
    # Table 2 'IPSSb on disease progression rate' 0.027; control stream THETA(6) 0.0267, centred at 16.00.

    e_alcohol_use_thalf_pbo <- fixed(-0.135)
    label("Fractional change in the placebo half-life for alcohol users (unitless)")
    # Table 2 'Alcohol user status at baseline on placebo t1/2' -0.135; control stream THETA(7) -0.135.

    # ================================================================
    # Treatment effects (Table 3), estimated with the placebo model fixed.
    # Each is a constant first-order rate added to the placebo effect.
    # ================================================================
    ldelta_tamsulosin <- log(0.015)
    label("Tamsulosin monotherapy effect DELTA_tamsulosin (1/month)")
    # Table 3 'DELTA tamsulosin, effect of tamsulosin' 0.015 months-1.

    ldelta_dutasteride <- log(0.016)
    label("Dutasteride monotherapy effect DELTA_dutasteride (1/month)")
    # Table 3 'DELTA dutasteride, effect of dutasteride' 0.016 months-1.

    ldelta_ww <- log(0.018)
    label("Watchful-waiting effect DELTA_WW (1/month)")
    # Table 3 'DELTA WW, effect of watchful waiting' 0.018 months-1.

    ldelta_fdc <- log(0.032)
    label("Tamsulosin-dutasteride combination effect DELTA_FDC (1/month)")
    # Table 3 'DELTA FDC, effect of combination therapy' 0.032 months-1.

    # ================================================================
    # Between-subject variability (variances, log-normal).
    # Placebo block: diagonal from Table 2. Table 2 prints variances only,
    # but the Supporting Information states that the final model carried
    # eta correlations between DISP and the placebo magnitude and
    # half-life, and the control stream fixes an $OMEGA BLOCK(3)
    # (0.982; -0.977, 3.29; 0.524, -3.01, 2.96). Its correlations
    # (-0.544, 0.307, -0.965) are applied to the Table 2 variances:
    # cov = r * sqrt(var_i * var_j).
    # ================================================================
    etaldisp + etaldelta_pbo + etalthalf_pbo ~ fixed(c(
      0.997,
      -0.9117, 2.822,
      0.5097, -2.6909, 2.758
    ))
    # Table 2 'IIV on disease progression rate' 0.997, 'IIV on placebo effect' 2.822, 'IIV on placebo t1/2' 2.758; off-diagonals from control stream $OMEGA BLOCK(3) correlations

    etaldelta_tamsulosin ~ 1.923  # Table 3 'IIV on the effect of tamsulosin' 1.923
    etaldelta_dutasteride ~ 1.734  # Table 3 'IIV on the effect of dutasteride' 1.734
    etaldelta_ww ~ 1.812  # Table 3 'IIV on the effect of watchful waiting' 1.812
    etaldelta_fdc ~ 1.449  # Table 3 'IIV on the effect of combination therapy' 1.449

    # ================================================================
    # Residual error: additive, one SD per treatment arm (Tables 2 and 3,
    # reported as SD; control stream $SIGMA carries variances).
    # ================================================================
    addSd_pbo <- fixed(3.224)
    label("Additive residual SD, placebo arm (IPSS points)")
    # Table 2 'Additive RUV' 3.224; control stream SIGMA(1,1) 10.3 FIX = 3.21^2.

    addSd_tamsulosin <- 3.613
    label("Additive residual SD, tamsulosin arm (IPSS points)")
    # Table 3 tamsulosin 'Additive RUV' 3.613.

    addSd_dutasteride <- 3.506
    label("Additive residual SD, dutasteride arm (IPSS points)")
    # Table 3 dutasteride 'Additive RUV' 3.506.

    addSd_ww <- 2.778
    label("Additive residual SD, watchful-waiting arm (IPSS points)")
    # Table 3 watchful waiting 'Additive RUV' 2.778.

    addSd_fdc <- 3.141
    label("Additive residual SD, combination-therapy arm (IPSS points)")
    # Table 3 combination therapy 'Additive RUV' 3.141.
  })
  model({
    # Placebo-arm indicator (control stream ARM = 0): all four TRT_ = 0.
    pbo <- 1 - TRT_TAMSULOSIN - TRT_DUTASTERIDE - TRT_TAMSULOSIN_DUTASTERIDE - TRT_WATCHFUL_WAITING

    # Control stream: DISP and DELTp carry their etas only when ARM = 0;
    # T12p carries its eta in every arm.
    disp <- exp(ldisp + etaldisp * pbo) * (1 + e_ipss_bl_disp * (IPSS_BL - 16))
    delta_pbo <- exp(ldelta_pbo + etaldelta_pbo * pbo) * (1 + e_t_sympt_bph_delta_pbo * (T_SYMPT_BPH - 4))
    thalf_pbo <- exp(lthalf_pbo + etalthalf_pbo) * (1 + e_alcohol_use_thalf_pbo * ALCOHOL_USE)

    delta_tamsulosin <- exp(ldelta_tamsulosin + etaldelta_tamsulosin)
    delta_dutasteride <- exp(ldelta_dutasteride + etaldelta_dutasteride)
    delta_ww <- exp(ldelta_ww + etaldelta_ww)
    delta_fdc <- exp(ldelta_fdc + etaldelta_fdc)

    # Equation 2: placebo effect decays from the start of treatment (t = 0).
    eff_pbo <- delta_pbo * exp(-t * log(2) / thalf_pbo)
    # Equation 3: constant treatment effect selected by the arm.
    eff_trt <- delta_tamsulosin * TRT_TAMSULOSIN + delta_dutasteride * TRT_DUTASTERIDE +
      delta_fdc * TRT_TAMSULOSIN_DUTASTERIDE + delta_ww * TRT_WATCHFUL_WAITING

    # Equation 1 (control stream DADT(1)); 35 is the IPSS scale maximum.
    ipss(0) <- IPSS_BL
    d/dt(ipss) <- disp * (1 - ipss / 35) - (eff_pbo + eff_trt) * ipss

    ruvAdd <- addSd_pbo * pbo + addSd_tamsulosin * TRT_TAMSULOSIN + addSd_dutasteride * TRT_DUTASTERIDE +
      addSd_fdc * TRT_TAMSULOSIN_DUTASTERIDE + addSd_ww * TRT_WATCHFUL_WAITING
    ipss ~ add(ruvAdd)
  })
}
