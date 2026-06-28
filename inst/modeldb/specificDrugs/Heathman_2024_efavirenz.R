Heathman_2024_efavirenz <- function() {
  description <- paste0(
    "Population pharmacokinetic model for efavirenz (EFV) and its 8-hydroxy- ",
    "and 7-hydroxy-metabolites in 135 healthy volunteers receiving a single ",
    "600 mg dose followed by 17 days of 600 mg/day (4594 plasma concentration ",
    "samples). Each of EFV, 8-OH EFV, and 7-OH EFV is a 2-compartment model ",
    "with the metabolite central volume fixed equal to that of EFV (lack of ",
    "identifiability). EFV absorption is sequential zero- (D1 = 1.74 h) plus ",
    "first-order (KA = 0.165/h). Two independent enzyme-turnover models drive ",
    "CYP2B6 and CYP2A6 autoinduction: R(t) = kout * (1 + Emax * Cc_EFV / ",
    "(EC50 + Cc_EFV)); dE/dt = R(t) - kout * E; E(0) = 1. CYP2B6 modulates the ",
    "EFV-to-8-OH formation arm (CL-EFV,2B6 = 3.64 L/h; Emax-2B6 = 15.5; ",
    "EC50-2B6 = 32000 nM = 10.10 mg/L) and the 8-OH-onward CYP2B6 arm ",
    "(CL-8OH,2B6 = 0.758 L/h); CYP2A6 modulates the EFV-to-7-OH formation arm ",
    "(CL-EFV,2A6 = 0.0947 L/h; Emax-2A6 = 4.22; EC50-2A6 = 12500 nM = 3.95 ",
    "mg/L). UGT2B7 elimination arms (CL-EFV,UGT = 0.0504 L/h; CL-8OH,UGT = ",
    "5.44 L/h) and the 7-OH total clearance (CL-7OH = 3.39 L/h) are not ",
    "autoinduced. CYP2B6 phenotype reduces the EFV-to-8-OH formation arm by ",
    "9.72% (IM) / 9.06% (PM, encoded as canonical CYP2B6_SM) and reduces the ",
    "CYP2B6 Emax by 53.5% (IM) / 93.2% (PM); NM (extensive metabolizer) is the ",
    "reference. PM subjects show essentially no autoinduction (effective ",
    "Emax-2B6 ~ 1.05) and accumulate over 2-3 weeks. IIV is exponential on PK ",
    "parameters (estimated for most; 15.1% CV fixed for several); residual ",
    "variability is proportional for all three analytes (25.8% / 28.0% / 29.9% ",
    "CV)."
  )
  reference <- paste(
    "Heathman MA, Metzger IF, Lu J, Gufford BT, Desta Z, Aruldhas BW.",
    "(2024). The Effect of CYP2B6 Genotype on the Clearance and Autoinduction",
    "of Efavirenz in Healthy Subjects and the Subsequent Impact on Efavirenz",
    "Exposure. Conference poster, Metrum Research Group and Indiana University",
    "School of Medicine. doi:10.70534/pgia9927",
    sep = " "
  )
  vignette <- "Heathman_2024_efavirenz"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L (= ug/mL) for plasma EFV (Cc), 8-OH EFV (Cc_8oh), and 7-OH EFV (Cc_7oh)"
  )

  covariateData <- list(
    CYP2B6_IM = list(
      description        = "1 = CYP2B6 intermediate-metabolizer phenotype, 0 = otherwise. Reference category (both CYP2B6_IM and CYP2B6_SM equal to 0) is the CYP2B6 normal (extensive) metabolizer phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal / extensive metabolizer (CYP2B6_IM = 0 and CYP2B6_SM = 0)",
      notes              = paste0(
        "Time-fixed (germline genotype). Heathman 2024 Methods 'Methods' ",
        "paragraph: 'CYP2B6 genotype was obtained and phenotype classified as ",
        "normal (NM), intermediate (IM), or poor metabolizer (PM)'. The ",
        "Heathman cohort frequency for NM / IM / PM is not reported in the ",
        "poster. IM reduces the EFV-to-8-OH formation arm (CL-EFV,2B6) by ",
        "9.72% (95% CI -2.85 to 20.8) and reduces Emax-2B6 by 53.5% (95% CI ",
        "37.1-65.6); both reductions are applied multiplicatively on the ",
        "log scale via `e_cyp2b6_im_cl_2b6` and `e_cyp2b6_im_emax_2b6`."
      ),
      source_name        = "CYP2B6 phenotype (intermediate metabolizer level)"
    ),
    CYP2B6_SM = list(
      description        = "1 = CYP2B6 slow- / poor-metabolizer phenotype (Heathman's 'PM'); 0 = otherwise. Reference category is the CYP2B6 normal (extensive) metabolizer phenotype.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Normal / extensive metabolizer (CYP2B6_IM = 0 and CYP2B6_SM = 0)",
      notes              = paste0(
        "Time-fixed (germline genotype). Heathman 2024's 'poor metabolizer' ",
        "(PM) corresponds to the canonical 'slow metabolizer' (SM; see ",
        "Robarge 2017, Bienczak 2016 precedents). PM reduces the EFV-to-8-OH ",
        "formation arm (CL-EFV,2B6) by 9.06% (95% CI -12.7 to 26.6) and ",
        "reduces Emax-2B6 by 93.2% (95% CI 87.6-96.3) -- i.e., PM subjects ",
        "show essentially no CYP2B6 autoinduction (effective Emax-2B6 ~ 1.05), ",
        "which drives the 2-3 week EFV accumulation kinetics reported in ",
        "Heathman 2024 Results paragraph 2 ('concentrations in PM subjects ",
        "continue to accumulate for 2 to 3 weeks')."
      ),
      source_name        = "CYP2B6 phenotype (poor metabolizer level; Heathman PM = canonical SM)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 135L,
    n_studies      = NA_integer_,
    age_range      = "(not reported in the conference poster)",
    weight_range   = "(not reported in the conference poster)",
    sex_female_pct = NA_real_,
    race_ethnicity = "(not reported in the conference poster)",
    disease_state  = "healthy volunteers",
    dose_range     = "single 600 mg oral dose of efavirenz, followed by 600 mg once-daily for 17 days (Heathman 2024 Methods 'Methods' paragraph)",
    regions        = "(not reported in the conference poster)",
    n_observations = "4594 plasma concentration samples across the three analytes, collected up to 144 h post-dose (Heathman 2024 Methods 'Methods' paragraph)",
    notes          = paste0(
      "Heathman 2024 is a Metrum Research Group / Indiana University School ",
      "of Medicine conference poster (PAGE-style; doi:10.70534/pgia9927). The ",
      "poster reports full final-model point estimates (Table 1, structural ",
      "parameters + autoinduction parameters + CYP2B6 phenotype effects), a ",
      "complete IIV table with %CV and shrinkage (visible in the original PDF; ",
      "some entries marked FIXED at 15.1% CV), and proportional residual ",
      "%CV for each of the three analytes (25.8% / 28.0% / 29.9%). Baseline ",
      "demographics (age, weight, sex, race, ethnicity) and study counts are ",
      "not provided in the poster text. The model was fit in NONMEM using the ",
      "SAEM estimation method (Heathman 2024 Methods 'Methods' paragraph)."
    )
  )

  ini({
    # ============================================================================
    # Structural PK -- EFV parent (2-compartment with sequential zero+first-order
    # absorption). Heathman 2024 Table 1 'Structural Parameters' block.
    # ============================================================================
    lka     <- log(0.165)
    label("First-order absorption rate ka (1/h)")
    # Heathman 2024 Table 1: KA = 0.165 1/h (95% CI 0.144-0.190).

    ld1     <- log(1.74)
    label("Zero-order absorption duration D1 (h); sequential zero-order input precedes first-order ka")
    # Heathman 2024 Table 1: D1 = 1.74 h (95% CI 1.59-1.91). The mean absorption
    # time is D1/2 + 1/ka = 0.87 + 6.06 = 6.93 h.

    lcl_2b6 <- log(3.64)
    label("EFV CYP2B6-mediated clearance arm CL-EFV,2B6 (L/h) at baseline E_2B6 = 1 and NM phenotype")
    # Heathman 2024 Table 1: CL-EFV,2B6 = 3.64 L/h (95% CI 3.33-3.98). This is
    # the formation rate of 8-OH EFV from EFV via CYP2B6.

    lcl_2a6 <- log(0.0947)
    label("EFV CYP2A6-mediated clearance arm CL-EFV,2A6 (L/h) at baseline E_2A6 = 1")
    # Heathman 2024 Table 1: CL-EFV,2A6 = 0.0947 L/h (95% CI 0.0788-0.114). This
    # is the formation rate of 7-OH EFV from EFV via CYP2A6.

    lcl_ugt <- log(0.0504)
    label("EFV UGT-mediated clearance arm CL-EFV,UGT (L/h); constant in time (not autoinduced)")
    # Heathman 2024 Table 1: CL-EFV,UGT = 0.0504 L/h (95% CI 0.0453-0.0560).

    lvc     <- log(3.99)
    label("EFV central volume of distribution Vc-EFV (L)")
    # Heathman 2024 Table 1: VC-EFV = 3.99 L (95% CI 3.44-4.63). Vc of both
    # metabolites is fixed equal to this value per Heathman 2024 Results
    # paragraph 1 ('the central volume of both metabolites were set equal to
    # that of EFV, due to lack of identifiability').

    lvp     <- log(520)
    label("EFV peripheral volume Vp-EFV (L)")
    # Heathman 2024 Table 1: VP-EFV = 520 L (95% CI 476-567).

    lq      <- log(28.3)
    label("EFV inter-compartmental clearance Q-EFV (L/h)")
    # Heathman 2024 Table 1: Q-EFV = 28.3 L/h (95% CI 24.7-32.4).

    # ============================================================================
    # Structural PK -- 8-OH EFV metabolite (2-compartment; Vc fixed equal to Vc-EFV).
    # ============================================================================
    lcl_2b6_8oh <- log(0.758)
    label("8-OH EFV CYP2B6-mediated clearance arm CL-8OH,2B6 (L/h) at baseline E_2B6 = 1")
    # Heathman 2024 Table 1: CL-8OH,2B6 = 0.758 L/h (95% CI 0.700-0.819).
    # Further CYP2B6-mediated metabolism of the 8-OH EFV metabolite; scaled
    # by enzyme_2b6(t).

    lcl_ugt_8oh <- log(5.44)
    label("8-OH EFV UGT-mediated clearance arm CL-8OH,UGT (L/h); constant in time")
    # Heathman 2024 Table 1: CL-8OH,UGT = 5.44 L/h (95% CI 4.97-5.94).

    lvp_8oh <- log(133)
    label("8-OH EFV peripheral volume Vp-8OH (L)")
    # Heathman 2024 Table 1: VP-8OH = 133 L (95% CI 120-147).

    lq_8oh  <- log(5.62)
    label("8-OH EFV inter-compartmental clearance Q-8OH (L/h)")
    # Heathman 2024 Table 1: Q-8OH = 5.62 L/h (95% CI 5.29-5.97).

    # ============================================================================
    # Structural PK -- 7-OH EFV metabolite (2-compartment; Vc fixed equal to Vc-EFV).
    # ============================================================================
    lcl_7oh <- log(3.39)
    label("7-OH EFV total clearance CL-7OH (L/h); constant in time")
    # Heathman 2024 Table 1: CL-7OH = 3.39 L/h (95% CI 2.96-3.87). The poster
    # reports a single lumped CL-7OH (no CYP2B6 or UGT decomposition), so this
    # arm is constant in time.

    lvp_7oh <- log(29.4)
    label("7-OH EFV peripheral volume Vp-7OH (L)")
    # Heathman 2024 Table 1: VP-7OH = 29.4 L (95% CI 22.6-38.2).

    lq_7oh  <- log(2.25)
    label("7-OH EFV inter-compartmental clearance Q-7OH (L/h)")
    # Heathman 2024 Table 1: Q-7OH = 2.25 L/h (95% CI 2.08-2.43).

    # ============================================================================
    # Autoinduction -- CYP2B6 enzyme turnover (Heathman 2024 Figure 1 ODE box).
    # R_2B6(t) = R_0 * (1 + Emax * Cc_EFV / (EC50 + Cc_EFV));
    # dE_2B6/dt = R_2B6(t) - kout_2B6 * E_2B6;  E_2B6(0) = 1; R_0 = kout_2B6 * E_0.
    # ============================================================================
    lemax_2b6 <- log(15.5)
    label("CYP2B6 autoinduction Emax (unitless multiplier on baseline enzyme synthesis rate)")
    # Heathman 2024 Table 1: EMAX-2B6 = 15.5 (95% CI 12.7-19.0). The
    # peak steady-state induction factor at saturating EFV is 1 + Emax = ~16.5.

    lkout_2b6 <- log(0.00684)
    label("CYP2B6 enzyme turnover first-order elimination rate (1/h)")
    # Heathman 2024 Table 1: KOUT-2B6 = 0.00684 1/h (95% CI 0.00630-0.00744);
    # turnover half-life ln(2)/kout = ~101 h.

    lec50_2b6 <- log(32000 * 315.68e-6)
    label("CYP2B6 autoinduction EC50 (mg/L = ug/mL); 32000 nM from Heathman 2024 Table 1 converted via EFV MW = 315.68 g/mol")
    # Heathman 2024 Table 1: EC50-2B6 = 32000 nM (95% CI 29100-35100).
    # Conversion: 32000 nM * 315.68 g/mol * 1e-6 mg-nmol/(L-g) = 10.10 mg/L.
    # EFV molecular formula C14H9ClF3NO2 -> MW = 315.68 g/mol.

    # ============================================================================
    # Autoinduction -- CYP2A6 enzyme turnover.
    # ============================================================================
    lemax_2a6 <- log(4.22)
    label("CYP2A6 autoinduction Emax (unitless multiplier on baseline enzyme synthesis rate)")
    # Heathman 2024 Table 1: EMAX-2A6 = 4.22 (95% CI 3.33-5.33).

    lkout_2a6 <- log(0.00982)
    label("CYP2A6 enzyme turnover first-order elimination rate (1/h)")
    # Heathman 2024 Table 1: KOUT-2A6 = 0.00982 1/h (95% CI 0.00895-0.0108);
    # turnover half-life ln(2)/kout = ~71 h.

    lec50_2a6 <- log(12500 * 315.68e-6)
    label("CYP2A6 autoinduction EC50 (mg/L); 12500 nM from Heathman 2024 Table 1 converted via EFV MW = 315.68 g/mol")
    # Heathman 2024 Table 1: EC50-2A6 = 12500 nM (95% CI 11400-13600).
    # Conversion: 12500 nM * 315.68 g/mol * 1e-6 = 3.946 mg/L.

    # ============================================================================
    # CYP2B6 phenotype effects (Heathman 2024 Table 1 'Effect of Phenotype' block).
    # Heathman's NM / IM / PM map to canonical (CYP2B6_IM = 0 and CYP2B6_SM = 0) /
    # CYP2B6_IM = 1 / CYP2B6_SM = 1. Encoded as multiplicative log-scale effects.
    # ============================================================================

    # Effect of CYP2B6 phenotype on CL-EFV,2B6 (the EFV-to-8-OH formation arm).
    e_cyp2b6_im_cl_2b6 <- log(1 - 0.0972)
    label("Log-ratio of CL-EFV,2B6 for CYP2B6 IM vs NM (unitless)")
    # Heathman 2024 Table 1: IM 'CL-EFV,2B6' % reduction = 9.72 (95% CI -2.85
    # to 20.8). The CI crosses zero; the point estimate is retained. Factor
    # = 1 - 0.0972 = 0.9028; log(0.9028) = -0.1023.

    e_cyp2b6_sm_cl_2b6 <- log(1 - 0.0906)
    label("Log-ratio of CL-EFV,2B6 for CYP2B6 PM (= canonical SM) vs NM (unitless)")
    # Heathman 2024 Table 1: PM 'CL-EFV,2B6' % reduction = 9.06 (95% CI -12.7
    # to 26.6). The CI crosses zero; the point estimate is retained. Factor
    # = 1 - 0.0906 = 0.9094; log(0.9094) = -0.0950.

    # Effect of CYP2B6 phenotype on EMAX-2B6 (the CYP2B6 autoinduction maximum).
    e_cyp2b6_im_emax_2b6 <- log(1 - 0.535)
    label("Log-ratio of EMAX-2B6 for CYP2B6 IM vs NM (unitless)")
    # Heathman 2024 Table 1: IM 'EMAX-2B6' % reduction = 53.5 (95% CI 37.1-65.6).
    # Factor = 1 - 0.535 = 0.465; log(0.465) = -0.7657.

    e_cyp2b6_sm_emax_2b6 <- log(1 - 0.932)
    label("Log-ratio of EMAX-2B6 for CYP2B6 PM (= canonical SM) vs NM (unitless)")
    # Heathman 2024 Table 1: PM 'EMAX-2B6' % reduction = 93.2 (95% CI 87.6-96.3).
    # Factor = 1 - 0.932 = 0.068; log(0.068) = -2.6882. The PM-effective
    # Emax-2B6 = 15.5 * 0.068 = ~1.05, i.e., PM subjects show essentially no
    # CYP2B6 autoinduction.

    # ============================================================================
    # IIV -- exponential on PK parameters; omega^2 = log(1 + CV^2).
    # Heathman 2024 IIV table (visible in the original PDF; full %CV column).
    # Entries marked '15.1 FIXED' in the source are wrapped in fixed() here.
    # The paper does not report IIV correlations, so each IIV is on the diagonal.
    # ============================================================================
    etalka         ~ log(1 + 0.893^2)
    # Heathman 2024 IIV table: IIV KA = 89.3% CV (95% CI 71.2-107). omega^2 =
    # log(1 + 0.893^2) = 0.580.

    etald1         ~ log(1 + 0.541^2)
    # Heathman 2024 IIV table: IIV D1 = 54.1% CV (95% CI 44.3-63.0). omega^2 =
    # log(1 + 0.541^2) = 0.255.

    etalcl_2b6     ~ log(1 + 0.356^2)
    # Heathman 2024 IIV table: IIV CL-EFV,2B6 = 35.6% CV (95% CI 27.1-42.7).

    etalcl_2a6     ~ log(1 + 1.37^2)
    # Heathman 2024 IIV table: IIV CL-EFV,2A6 = 137% CV (95% CI 109-168).

    etalcl_ugt     ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV CL-EFV,UGT = 15.1% CV FIXED.

    etalvc         ~ log(1 + 0.886^2)
    # Heathman 2024 IIV table: IIV VC-EFV = 88.6% CV (95% CI 66.3-110).

    etalvp         ~ log(1 + 0.537^2)
    # Heathman 2024 IIV table: IIV VP-EFV = 53.7% CV (95% CI 47.6-59.5).

    etalq          ~ log(1 + 0.865^2)
    # Heathman 2024 IIV table: IIV Q-EFV = 86.5% CV (95% CI 59.7-112).

    etalcl_2b6_8oh ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV CL-8OH,2B6 = 15.1% CV FIXED.

    etalcl_ugt_8oh ~ log(1 + 0.494^2)
    # Heathman 2024 IIV table: IIV CL-8OH,UGT = 49.4% CV (95% CI 41.8-56.3).

    etalvp_8oh     ~ log(1 + 0.369^2)
    # Heathman 2024 IIV table: IIV VP-8OH = 36.9% CV (95% CI 20.3-48.9).

    etalq_8oh      ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV Q-8OH = 15.1% CV FIXED.

    etalcl_7oh     ~ log(1 + 0.740^2)
    # Heathman 2024 IIV table: IIV CL-7OH = 74.0% CV (95% CI 59.8-87.3).

    etalvp_7oh     ~ log(1 + 0.937^2)
    # Heathman 2024 IIV table: IIV VP-7OH = 93.7% CV (95% CI 48.1-137).

    etalq_7oh      ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV Q-7OH = 15.1% CV FIXED.

    etalemax_2b6   ~ log(1 + 0.839^2)
    # Heathman 2024 IIV table: IIV EMAX-2B6 = 83.9% CV (95% CI 53.3-112).

    etalkout_2b6   ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV KOUT-2B6 = 15.1% CV FIXED.

    etalec50_2b6   ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV EC50-2B6 = 15.1% CV FIXED.

    etalemax_2a6   ~ log(1 + 1.71^2)
    # Heathman 2024 IIV table: IIV EMAX-2A6 = 171% CV (95% CI 137-209).

    etalkout_2a6   ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV KOUT-2A6 = 15.1% CV FIXED.

    etalec50_2a6   ~ fixed(log(1 + 0.151^2))
    # Heathman 2024 IIV table: IIV EC50-2A6 = 15.1% CV FIXED.

    # ============================================================================
    # Residual error -- proportional model for each of the three analytes.
    # Heathman 2024 Results paragraph 1: 'Residual variability was described using
    # a proportional model for all three analytes'.
    # ============================================================================
    propSd      <- 0.258
    label("Proportional residual error for plasma EFV (Cc; fraction)")
    # Heathman 2024 RUV table: RUV EFV = 25.8% CV (95% CI 25.2-26.4).

    propSd_8oh <- 0.280
    label("Proportional residual error for plasma 8-OH EFV (Cc_8oh; fraction)")
    # Heathman 2024 RUV table: RUV 8-OH EFV = 28.0% CV (95% CI 27.3-28.6).

    propSd_7oh <- 0.299
    label("Proportional residual error for plasma 7-OH EFV (Cc_7oh; fraction)")
    # Heathman 2024 RUV table: RUV 7-OH EFV = 29.9% CV (95% CI 29.2-30.5).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters.
    #    The CYP2B6 phenotype factor enters as a log-additive multiplicative
    #    effect on CL-EFV,2B6 and EMAX-2B6 (NM is the reference with both
    #    CYP2B6_IM and CYP2B6_SM = 0).
    # -----------------------------------------------------------------------
    ka         <- exp(lka + etalka)
    d1         <- exp(ld1 + etald1)
    vc         <- exp(lvc + etalvc)
    vp         <- exp(lvp + etalvp)
    q          <- exp(lq  + etalq)

    cl_2b6     <- exp(lcl_2b6 + etalcl_2b6
                      + e_cyp2b6_im_cl_2b6 * CYP2B6_IM
                      + e_cyp2b6_sm_cl_2b6 * CYP2B6_SM)
    cl_2a6     <- exp(lcl_2a6 + etalcl_2a6)
    cl_ugt     <- exp(lcl_ugt + etalcl_ugt)

    # 8-OH EFV metabolite parameters; Vc-8OH set equal to Vc-EFV per the paper.
    vc_8oh     <- vc
    vp_8oh     <- exp(lvp_8oh + etalvp_8oh)
    q_8oh      <- exp(lq_8oh  + etalq_8oh)
    cl_2b6_8oh <- exp(lcl_2b6_8oh + etalcl_2b6_8oh)
    cl_ugt_8oh <- exp(lcl_ugt_8oh + etalcl_ugt_8oh)

    # 7-OH EFV metabolite parameters; Vc-7OH set equal to Vc-EFV per the paper.
    vc_7oh     <- vc
    vp_7oh     <- exp(lvp_7oh + etalvp_7oh)
    q_7oh      <- exp(lq_7oh  + etalq_7oh)
    cl_7oh     <- exp(lcl_7oh + etalcl_7oh)

    # Autoinduction parameters. CYP2B6 phenotype enters EMAX-2B6
    # multiplicatively (log-additive on the log scale).
    emax_2b6   <- exp(lemax_2b6 + etalemax_2b6
                      + e_cyp2b6_im_emax_2b6 * CYP2B6_IM
                      + e_cyp2b6_sm_emax_2b6 * CYP2B6_SM)
    kout_2b6   <- exp(lkout_2b6 + etalkout_2b6)
    ec50_2b6   <- exp(lec50_2b6 + etalec50_2b6)

    emax_2a6   <- exp(lemax_2a6 + etalemax_2a6)
    kout_2a6   <- exp(lkout_2a6 + etalkout_2a6)
    ec50_2a6   <- exp(lec50_2a6 + etalec50_2a6)

    # -----------------------------------------------------------------------
    # 2. EFV plasma concentration (mg/L = ug/mL).
    # -----------------------------------------------------------------------
    Cc <- central / vc

    # -----------------------------------------------------------------------
    # 3. Enzyme synthesis rates (Heathman 2024 Figure 1 ODE box).
    #    R_iso(t) = R_0 * (1 + Emax * Cc_EFV / (EC50 + Cc_EFV));
    #    R_0 = kout_iso * E_0 with E_0 = 1, so the baseline synthesis rate
    #    is just kout_iso and the induced rate scales as (1 + Emax*Cc/(EC50+Cc)).
    # -----------------------------------------------------------------------
    R_2b6 <- kout_2b6 * (1 + emax_2b6 * Cc / (ec50_2b6 + Cc))
    R_2a6 <- kout_2a6 * (1 + emax_2a6 * Cc / (ec50_2a6 + Cc))

    # -----------------------------------------------------------------------
    # 4. EFV ODE system.
    #    CYP2B6 and CYP2A6 metabolic arms are multiplied by the
    #    corresponding relative enzyme amount (enzyme_2b6 / enzyme_2a6,
    #    both starting at 1). The UGT arm is constant in time. The 2-cpt
    #    distribution is parameterised by Q-EFV and Vp-EFV.
    # -----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot
                          - (cl_2b6 * enzyme_2b6
                             + cl_2a6 * enzyme_2a6
                             + cl_ugt) / vc * central
                          - q / vc * central
                          + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # Sequential zero+first-order absorption: dose enters depot uniformly
    # over [0, d1] (zero-order) then drains via ka (first-order).
    dur(depot) <- d1

    # -----------------------------------------------------------------------
    # 5. 8-OH EFV ODE system. Formation flux from EFV via CYP2B6 (scaled by
    #    enzyme_2b6); further elimination via CYP2B6 (scaled by enzyme_2b6)
    #    plus UGT (constant). Vc-8OH = Vc-EFV per the paper.
    # -----------------------------------------------------------------------
    Cc_8oh <- central_8oh / vc_8oh
    d/dt(central_8oh)     <-  cl_2b6 * enzyme_2b6 / vc * central
                              - (cl_2b6_8oh * enzyme_2b6 + cl_ugt_8oh) / vc_8oh * central_8oh
                              - q_8oh / vc_8oh * central_8oh
                              + q_8oh / vp_8oh * peripheral1_8oh
    d/dt(peripheral1_8oh) <-  q_8oh / vc_8oh * central_8oh
                              - q_8oh / vp_8oh * peripheral1_8oh

    # -----------------------------------------------------------------------
    # 6. 7-OH EFV ODE system. Formation flux from EFV via CYP2A6 (scaled by
    #    enzyme_2a6); single lumped CL-7OH (constant in time).
    # -----------------------------------------------------------------------
    Cc_7oh <- central_7oh / vc_7oh
    d/dt(central_7oh)     <-  cl_2a6 * enzyme_2a6 / vc * central
                              - cl_7oh / vc_7oh * central_7oh
                              - q_7oh / vc_7oh * central_7oh
                              + q_7oh / vp_7oh * peripheral1_7oh
    d/dt(peripheral1_7oh) <-  q_7oh / vc_7oh * central_7oh
                              - q_7oh / vp_7oh * peripheral1_7oh

    # -----------------------------------------------------------------------
    # 7. Enzyme turnover ODEs (Heathman 2024 Figure 1 ODE box).
    # -----------------------------------------------------------------------
    d/dt(enzyme_2b6) <- R_2b6 - kout_2b6 * enzyme_2b6
    d/dt(enzyme_2a6) <- R_2a6 - kout_2a6 * enzyme_2a6

    # Initial enzyme amounts normalised to baseline = 1.
    enzyme_2b6(0) <- 1
    enzyme_2a6(0) <- 1

    # -----------------------------------------------------------------------
    # 8. Observation models. Proportional residual error for each analyte.
    # -----------------------------------------------------------------------
    Cc     ~ prop(propSd)
    Cc_8oh ~ prop(propSd_8oh)
    Cc_7oh ~ prop(propSd_7oh)
  })
}
