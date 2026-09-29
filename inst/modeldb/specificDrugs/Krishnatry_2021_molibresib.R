Krishnatry_2021_molibresib <- function() {
  description <- "Semimechanistic liver-compartment population PK model for oral molibresib (GSK525762, a BET bromodomain inhibitor) and its active metabolite composite GSK3529246 in adults with advanced solid tumours. Molibresib is absorbed after a lag time by first-order absorption into a physiologic liver compartment (plasma flow 55 L/h, volume 1.5 L, both fixed), distributes in a two-compartment systemic model, and is eliminated only by hepatic extraction. The extraction ratio is proportional to a relative hepatic enzyme amount, whose production rate rises linearly with the molibresib concentration in the liver, so the model reproduces the autoinduction that lowers molibresib exposure on repeated dosing. Extracted molibresib is converted on a 1:1 molar basis into GSK3529246 through one transit compartment, and GSK3529246 has two-compartment disposition with first-order elimination. Body weight scales both central volumes with one shared exponent, and time-varying aspartate aminotransferase lowers GSK3529246 clearance. Amounts are in umol and concentrations in umol/L (molibresib 424 g/mol, GSK3529246 396 g/mol)."
  reference <- paste(
    "Krishnatry AS, Voelkner A, Dhar A, Prohn M, Ferron-Brady G (2021).",
    "Population pharmacokinetic modeling of molibresib and its active",
    "metabolites in patients with solid tumors: A semimechanistic",
    "autoinduction model. CPT Pharmacometrics Syst Pharmacol",
    "10(7):709-722. doi:10.1002/psp4.12639.",
    "Structure from the final NONMEM control stream (Run 32) in",
    "Supplementary Text S2; parameter values from Table 3.",
    sep = " "
  )
  vignette <- "Krishnatry_2021_molibresib"
  # Methods: 'Plasma concentrations were converted in molar concentrations
  # using the molecular weights of molibresib (424 g/mol) and GSK3529246
  # (396 g/mol) for the modeling.' The model therefore runs in molar units,
  # which keeps the 1:1 parent-to-metabolite conversion mass-balanced and fixes
  # the unit of the induction slope (per umol/L of molibresib in the liver).
  # Dose in umol = dose in mg / 424 * 1000.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling normalised at 70 kg (Supplementary Text S2:",
        "COVWT_V = (WT/70)**THETA(16)). One estimated exponent, 0.717, is",
        "shared by the molibresib central volume V1/F and the GSK3529246",
        "central volume mV1/F (Table 3 'WT on V1/F and mV1/F'); the control",
        "stream multiplies both volumes by the same COVWT_V. The cohort",
        "median weight was 69.6 kg (range 34-120 kg).",
        sep = " "
      ),
      source_name = "WT"
    ),
    AST = list(
      description = "Aspartate aminotransferase, time-varying",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power scaling normalised at 28 U/L (Supplementary Text S2:",
        "COVAST_MCL = (AST/28)**THETA(17)) on GSK3529246 clearance only;",
        "exponent -0.194, so higher AST lowers mCL/F. The source calls the",
        "covariate time-varying AST and reports it in IU/L (Table 1 mean",
        "33.7, SD 19 IU/L at baseline). The control stream reads the",
        "time-varying AST column, not the baseline BAST column.",
        sep = " "
      ),
      source_name = "AST"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Selected by the GAM screen for CL/F and mCL/F and entered in the full model (Run 17), then removed from both in backward elimination (Table S4 steps 1-2). Not in the final model.",
      source_name = "AGE"
    ),
    ALT = list(
      description = "Alanine aminotransferase, time-varying",
      units = "U/L",
      type = "continuous",
      notes = "Screened in the GAM analysis (time-varying and baseline ALT); not carried into the full model.",
      source_name = "ALT"
    ),
    ALB = list(
      description = "Serum albumin at baseline",
      units = "g/L",
      type = "continuous",
      notes = "Baseline albumin screened in the GAM analysis; not carried into the full model.",
      source_name = "BALB"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the GAM analysis; not carried into the full model. 53% of the analysis population was female.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "molibresib", units = "umol", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "molibresib", units = "umol", specimen = "tissue", verified = TRUE),
    transit1_gsk3529246 = list(
      analyte = "GSK3529246 (active metabolite composite of molibresib)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(analyte = "molibresib", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "molibresib", units = "umol", specimen = "plasma", verified = TRUE),
    enzyme = list(
      analyte = "relative amount of the hepatic enzyme (CYP3A4) that clears molibresib",
      units = "(fraction of baseline)",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_gsk3529246 = list(
      analyte = "GSK3529246 (active metabolite composite of molibresib)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_gsk3529246 = list(
      analyte = "GSK3529246 (active metabolite composite of molibresib)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 193L,
    n_studies = 1L,
    age_range = "16-86 years (median 58; mean 55.7, SD 14; Table 1)",
    weight_range = "34-120 kg (median 69.6; mean 71.7, SD 17; Table 1)",
    sex_female_pct = 53,
    race_ethnicity = c(White = 83, Asian = 6, Black = 6, Missing = 5),
    disease_state = "Advanced solid tumours: NUT (nuclear protein in testis) midline carcinoma, colorectal, breast (including triple-negative and ER-positive), prostate (castration-resistant), lung (including small-cell), gastrointestinal stromal tumour, neuroblastoma and multiple myeloma (Table 1).",
    dose_range = "Part 1: 2-100 mg once daily or 20-40 mg twice daily of the amorphous free-base formulation, and 80 mg once daily of the besylate salt in a 10-patient besylate substudy; Part 2: 75 mg once daily of the besylate salt (the recommended phase II dose). Table S2.",
    regions = "Not reported.",
    notes = paste(
      "First-time-in-human phase I/II study BET115521 (Part 1 dose escalation,",
      "n = 94; Part 2 expansion cohorts, n = 99). 2681 molibresib and 814",
      "GSK3529246 plasma concentrations above the lower limit of",
      "quantification (Table S2); 260 and 144 post-dose BLQ samples were",
      "excluded. GSK3529246 was assayed only in Part 2 and the Part 1 80 mg",
      "cohort (n = 131). The two major active metabolites (GSK3536835,",
      "ethyl-hydroxy, and GSK3529246, N-desethyl) were measured together after",
      "full conversion of one to the other and reported as the composite",
      "GSK3529246. Model development truncated the data at 1000 h. NONMEM 7.4,",
      "FOCE with interaction.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Molibresib (parent). Krishnatry 2021 Table 3 'Estimate' column
    # (final model, Run 32). Supplementary Text S2 estimates every THETA on
    # the log scale ('; LOG' tag on each $THETA line), so the values are
    # carried as logs. The $THETA values printed in Text S2 are the run's
    # initial values and differ slightly from Table 3; Table 3 is used.
    # ------------------------------------------------------------------
    lka <- log(4.08)
    label("Molibresib absorption rate constant ka (1/h)") # Table 3 'Molibresib absorption rate' ka 4.08 /h (95% CI 3.27-5.1)
    lcl <- log(9.02)
    label("Molibresib clearance CL/F at baseline enzyme amount (L/h)") # Table 3 'Molibresib clearance' CL/F 9.02 L/h (95% CI 8.34-9.76); Text S2 EH = CL*A(6)/QH
    lvc <- log(53.1)
    label("Molibresib central volume V1/F at 70 kg (L)") # Table 3 'Molibresib central volume of distribution' V1/F 53.1 L (95% CI 49.4-57)
    lqh <- fixed(log(55.0))
    label("Liver plasma flow Qh (L/h)") # Table 3 'Liver plasma flow' Qh 55.0 L/h Fixed; Methods: hepatic blood flow 100 L/h x (1 - haematocrit 0.45); Text S2 $THETA 4.007333 FIX = log(55.0)
    lvh <- fixed(log(1.50))
    label("Liver volume Vh (L)") # Table 3 'Liver volume' Vh/F 1.50 L Fixed; Text S2 $THETA 0.4054651 FIX = log(1.5)
    lq <- log(0.999)
    label("Molibresib intercompartmental clearance Q1 (L/h)") # Table 3 'Molibresib intercompartmental clearance' Q1 0.999 L/h (95% CI 0.774-1.29)
    lvp <- log(17.4)
    label("Molibresib peripheral volume V2/F (L)") # Table 3 'Molibresib peripheral volume of distribution' V2/F 17.4 L (95% CI 14.3-21.2)
    ltlag <- log(0.132)
    label("Absorption lag time ALAG (h)") # Table 3 'Lagtime' ALAG 0.132 h (95% CI 0.118-0.149)

    # ------------------------------------------------------------------
    # Autoinduction. Text S2 $DES:
    #   CH = A(2)/VH; IND = SLP*CH
    #   DADT(6) = KIN*(1 + IND) - KOUT*A(6), KOUT = KIN/BASE, BASE = 1
    # ------------------------------------------------------------------
    lslope <- log(0.922)
    label("Linear induction slope on enzyme production (L/umol, per umol/L of molibresib in the liver)") # Table 3 'Induction slope' Slope 0.922 (95% CI 0.687-1.24)
    lkin <- fixed(log(0.00550))
    label("Enzyme production rate kin; also the enzyme turnover rate kout (1/h)") # Table 3 'Enzyme production rate' kin 0.00550 /h Fixed; Text S1: fixed to 0.0055/h, an enzyme half-life of 126 h; Text S2 $THETA -5.203 FIX

    # ------------------------------------------------------------------
    # GSK3529246 (active metabolite composite). Table 3.
    # ------------------------------------------------------------------
    lktr_gsk3529246 <- log(11.5)
    label("GSK3529246 transit rate constant mka, transit compartment to GSK3529246 central (1/h)") # Table 3 'GSK3529246 transit rate' mka 11.5 /h (95% CI 7.76-17.1)
    lcl_gsk3529246 <- log(12.8)
    label("GSK3529246 clearance mCL/F at AST 28 U/L (L/h)") # Table 3 'GSK3529246 clearance' mCL/F 12.8 L/h (95% CI 11.5-14.4)
    lvc_gsk3529246 <- log(62.1)
    label("GSK3529246 central volume mV1/F at 70 kg (L)") # Table 3 'GSK3529246 central volume of distribution' mV1/F 62.1 L (95% CI 55.2-69.7)
    lq_gsk3529246 <- log(5.63)
    label("GSK3529246 intercompartmental clearance mQ1 (L/h)") # Table 3 'GSK3529246 intercompartmental clearance' mQ1 5.63 L/h (95% CI 4.04-7.84)
    lvp_gsk3529246 <- log(140)
    label("GSK3529246 peripheral volume mV2/F (L)") # Table 3 'GSK3529246 peripheral volume of distribution' mV2/F 140 L (95% CI 104-188)

    # ------------------------------------------------------------------
    # Covariate effects. Table 3; reference values from Text S2.
    # ------------------------------------------------------------------
    e_wt_vc <- 0.717
    label("Power exponent of body weight on V1/F and mV1/F, normalised at 70 kg (unitless)") # Table 3 'WT on V1/F and mV1/F' 0.717 (RSE 5.70%, 95% CI 0.637-0.796)
    e_ast_cl_gsk3529246 <- -0.194
    label("Power exponent of AST on GSK3529246 clearance, normalised at 28 U/L (unitless)") # Table 3 'AST on mCL/F' -0.194 (RSE 34.3%, 95% CI -0.324 to -0.0635)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 3 lists omega^2 (variances); the
    # footnote derives IIV% as sqrt(omega^2)*100 (e.g. sqrt(1.65) = 128%).
    # Text S2 fixes the omegas of Qh, Vh, Q1, V2, ALAG, slope, kin, mQ1 and
    # mV2 to 0; those etas are omitted here.
    # ------------------------------------------------------------------
    etalka ~ 1.65 # Table 3 'IIV on absorption rate' omega2 ka 1.65 (IIV 128%)
    etalcl ~ 0.205 # Table 3 'IIV on molibresib clearance' omega2 CL/F 0.205 (IIV 45.3%)
    etalvc ~ 0.0687 # Table 3 'IIV on molibresib central volume of distribution' omega2 V1/F 0.0687 (IIV 26.2%)
    etalktr_gsk3529246 ~ 0.488 # Table 3 'IIV on GSK3529246 transit rate' omega2 mka 0.488 (IIV 69.9%)
    # Table 3 prints the covariance as '174'; its 95% CI (0.107-0.24) and the
    # $OMEGA BLOCK(2) of Text S2 place it at 0.174 (correlation 0.80).
    etalcl_gsk3529246 + etalvc_gsk3529246 ~ c(0.227, 0.174, 0.206) # Table 3 omega2 mCL/F 0.227 (IIV 47.6%), covariance mCL-mV 0.174 (printed '174'; 95% CI 0.107-0.24), omega2 mV1/F 0.206 (IIV 45.4%)

    # ------------------------------------------------------------------
    # Residual error. Text S2: Y = IPRED + IPRED*EPS(n) with a $SIGMA
    # BLOCK(2), so the Table 3 values are variances and the proportional
    # SD is their square root. The covariance 0.0918 between the two
    # endpoints (correlation 0.56) cannot be expressed in nlmixr2 and is
    # dropped; see the vignette.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.207)
    label("Proportional residual error SD, molibresib (fraction)") # Table 3 'Proportional error molibresib' sigma2 0.207 (RSE 7.40%); SD = sqrt(0.207) = 0.455
    propSd_gsk3529246 <- sqrt(0.132)
    label("Proportional residual error SD, GSK3529246 (fraction)") # Table 3 'Proportional error GSK3529246' sigma2 0.132 (RSE 10.9%); SD = sqrt(0.132) = 0.363
  })

  model({
    # Covariate terms (Text S2 $PK)
    cov_wt_vc <- (WT / 70)^e_wt_vc
    cov_ast_cl_gsk3529246 <- (AST / 28)^e_ast_cl_gsk3529246

    # Molibresib individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * cov_wt_vc
    qh <- exp(lqh)
    vh <- exp(lvh)
    q <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag)
    slope <- exp(lslope)
    kin <- exp(lkin)
    # Text S2: BASE = 1, KOUT = KIN/BASE, A_0(6) = BASE
    kout <- kin

    # GSK3529246 individual parameters
    ktr_gsk3529246 <- exp(lktr_gsk3529246 + etalktr_gsk3529246)
    cl_gsk3529246 <- exp(lcl_gsk3529246 + etalcl_gsk3529246) * cov_ast_cl_gsk3529246
    vc_gsk3529246 <- exp(lvc_gsk3529246 + etalvc_gsk3529246) * cov_wt_vc
    q_gsk3529246 <- exp(lq_gsk3529246)
    vp_gsk3529246 <- exp(lvp_gsk3529246)

    # Micro-constants
    k12 <- q / vc
    k21 <- q / vp
    kel_gsk3529246 <- cl_gsk3529246 / vc_gsk3529246
    k12_gsk3529246 <- q_gsk3529246 / vc_gsk3529246
    k21_gsk3529246 <- q_gsk3529246 / vp_gsk3529246

    # Liver model (Text S2 $DES). The hepatic extraction ratio is
    # proportional to the relative enzyme amount, EH = CL * enzyme / Qh, and
    # the fraction escaping the liver is FH = 1 - EH. Qh * FH / Vh returns
    # molibresib to the systemic circulation, Qh * EH / Vh converts it to
    # GSK3529246, and systemic molibresib re-enters the liver at Qh / V1.
    # Text S2 adds 1E-9 to CH as a numerical guard; it shifts the enzyme
    # baseline by about 1e-9 and is omitted.
    ch <- liver / vh
    eh <- cl * enzyme / qh
    fh <- 1 - eh
    k_liver_central <- qh * fh / vh
    k_liver_transit1 <- qh * eh / vh
    k_central_liver <- qh / vc

    enzyme(0) <- 1

    d/dt(depot) <- -ka * depot
    d/dt(liver) <- ka * depot - k_liver_central * liver + k_central_liver * central - k_liver_transit1 * liver
    d/dt(transit1_gsk3529246) <- k_liver_transit1 * liver - ktr_gsk3529246 * transit1_gsk3529246
    d/dt(central) <- k_liver_central * liver - k_central_liver * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(enzyme) <- kin * (1 + slope * ch) - kout * enzyme
    d/dt(central_gsk3529246) <- ktr_gsk3529246 * transit1_gsk3529246 - k12_gsk3529246 * central_gsk3529246 + k21_gsk3529246 * peripheral1_gsk3529246 - kel_gsk3529246 * central_gsk3529246
    d/dt(peripheral1_gsk3529246) <- k12_gsk3529246 * central_gsk3529246 - k21_gsk3529246 * peripheral1_gsk3529246

    alag(depot) <- tlag

    # Observations (umol/L)
    Cc <- central / vc
    Cc_gsk3529246 <- central_gsk3529246 / vc_gsk3529246

    Cc ~ prop(propSd)
    Cc_gsk3529246 ~ prop(propSd_gsk3529246)
  })
}
