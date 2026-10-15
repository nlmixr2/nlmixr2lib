Gupta_2022_mobocertinib <- function() {
  description <- "Joint semimechanistic population PK model for oral mobocertinib (an irreversible EGFR exon 20 insertion tyrosine kinase inhibitor) and its two active metabolites AP32960 and AP32914 in healthy adult volunteers and adults with metastatic non-small cell lung cancer. Mobocertinib is absorbed from the depot through three transit compartments (transit rate equal to ka) into a two-compartment disposition model. Its clearance forms AP32960 (fixed fraction 0.62, two-compartment disposition) and AP32914 (fixed fraction 0.08, one-compartment disposition), and the remaining 0.30 is eliminated. Auto-induction is a turnover enzyme pool whose production rate is stimulated by an Emax function of the molar sum of the three plasma concentrations; the relative enzyme amount multiplies the clearances of all three moieties. Healthy-volunteer status raises the mobocertinib central volume and all three clearances. Amounts are in umol and concentrations in umol/L."
  reference <- paste(
    "Gupta N, Pierrillas PB, Hanley MJ, Zhang S, Diderichsen PM (2022).",
    "Population pharmacokinetics of mobocertinib in healthy volunteers and",
    "patients with non-small cell lung cancer. CPT Pharmacometrics Syst",
    "Pharmacol 11(6):731-744. doi:10.1002/psp4.12785.",
    "Structure from the final NONMEM control stream (Supporting Information,",
    "supplementary data file 'Final model control stream'); parameter values",
    "from Table 3.",
    sep = " "
  )
  vignette <- "Gupta_2022_mobocertinib"
  # Control stream $ERROR: CP = A(2)/V2, CM1 = A(4)/VM60, CM2 = A(5)/VM14, each
  # commented 'concentration in uM', and EC50 is reported in nM (Table 3). The
  # model therefore runs in molar units: amounts in umol, volumes in L,
  # concentrations in umol/L. Metabolite formation FM60 * K20 * A(2) is a molar
  # flux, so the fixed fractions 0.62 and 0.08 are molar fractions. Dose in
  # umol = dose in mg / 585.7 g/mol * 1000 (mobocertinib free base,
  # C32H39N7O4, PubChem CID 118607832); 160 mg = 273.2 umol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy volunteer, 0 = patient with mNSCLC)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with metastatic NSCLC)",
      notes = paste(
        "Control stream $PK: IF(HV.EQ.1) V2HV = (1 + THETA(18)), CLHV = (1 +",
        "THETA(15)), CLM60HV = (1 + THETA(17)), CLM14HV = (1 + THETA(16)),",
        "each equal to 1 when HV = 0. Linear relative effects 0.787 on the",
        "mobocertinib central volume and 0.900 / 0.738 / 0.909 on the",
        "mobocertinib / AP32960 / AP32914 clearances (Table 3). 110 of 427",
        "participants (25.8%) were healthy volunteers, all of whom received",
        "single doses only (Table 1).",
        sep = " "
      ),
      source_name = "HV"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Table 2 mean 52.9, SD 17.0; range 18-86 years); not significant in the stepwise covariate search. Simulated AUC24 at the 5th / 95th percentiles differed by at most 12% from the median (Figure 4b).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 2 mean 70.4, SD 15.9; range 37.3-132 kg); not retained. Continuous covariates were tested as power functions normalised to the median.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (57.1% female, Table 2); not retained.",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race screened (White 66.3%, Asian 26.9%, Black 4.4%, Table 2); not retained.",
      source_name = "RACE"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race screened; not retained.",
      source_name = "RACE"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Table 2 mean 41.2, SD 20.7 g/L); not retained.",
      source_name = "ALB"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Table 2 mean 23.7, SD 12.9 U/L); not retained.",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Table 2 mean 21.8, SD 15.9 U/L); not retained.",
      source_name = "ALT"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (Table 2 mean 9.50, SD 4.68 umol/L); not retained.",
      source_name = "BILI"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD equation)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened (Table 2 mean 92.7, SD 26.7 mL/min/1.73 m^2); not retained. The abstract reports no clinically meaningful effect of mild-to-moderate renal impairment (eGFR 30-89 mL/min/1.73 m^2).",
      source_name = "EGFR"
    ),
    SMOKE = list(
      description = "Smoking status",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened (never 74.9%, former 24.1%, current 0.7%, Table 2); not retained.",
      source_name = "SMOK"
    ),
    FORM_MOBOCERTINIB_CAPSULE_C = list(
      description = "Mobocertinib drug product capsule C (vs capsules A/B)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a covariate (capsule C 30.0%, Table 2); not retained. A separate model-based relative bioavailability refit (Supporting Information Table S4) estimated a capsule C effect of -0.0523 on ka and 0.0875 on F, both within the bioequivalence bounds; that refit is not the final model and is not encoded.",
      source_name = "capsule"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mobocertinib", units = "umol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "mobocertinib", units = "umol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "mobocertinib", units = "umol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "mobocertinib", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mobocertinib", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mobocertinib", units = "umol", specimen = "plasma", verified = TRUE),
    central_ap32960 = list(
      analyte = "AP32960 (active metabolite of mobocertinib)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_ap32960 = list(
      analyte = "AP32960 (active metabolite of mobocertinib)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    central_ap32914 = list(
      analyte = "AP32914 (active metabolite of mobocertinib)",
      units = "umol",
      specimen = "plasma",
      verified = TRUE
    ),
    enzyme = list(
      analyte = "relative amount of the inducible enzyme (CYP3A) clearing mobocertinib, AP32960 and AP32914",
      units = "(fraction of baseline)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 427L,
    n_studies = 4L,
    age_range = "18-86 years (mean 52.9, SD 17.0; Table 2)",
    weight_range = "37.3-132 kg (mean 70.4, SD 15.9; Table 2)",
    sex_female_pct = 57.1,
    race_ethnicity = c(
      White = 66.3,
      Asian = 26.9,
      Black = 4.4,
      `American Indian/Alaskan native` = 0.2,
      `Other/multiple` = 2.1
    ),
    disease_state = "317 adults with metastatic non-small cell lung cancer (EGFR or HER2 mutations; 289 with at least one prior treatment) and 110 healthy adult volunteers (Table 2).",
    dose_range = "Patients: 5-180 mg once daily or 40-60 mg twice daily (recommended dose 160 mg once daily); healthy volunteers: single doses of 20-160 mg (Table 1).",
    regions = "Global phase I/II study (NCT02716116), Japanese phase I/II study (NCT03807778), and two phase I healthy-volunteer studies (NCT03482453, NCT03928327).",
    notes = paste(
      "5880 mobocertinib, 5880 AP32960 and 5879 AP32914 post-dose",
      "observations (Supporting Information Table S1). BLQ samples handled",
      "by a modified M6 method. Only treatment periods without CYP3A",
      "modulators were used from the drug-drug interaction study. NONMEM",
      "7.3, log-transform-both-sides.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Mobocertinib. Table 3 'Untransformed parameter' column; the control
    # stream estimates every structural THETA on the log scale (TVKA =
    # EXP(THETA(1)), etc.). The untransformed column carries one more
    # significant figure than the log-scale 'Estimate' column and is used.
    # ------------------------------------------------------------------
    lka <- log(2.12)
    label("Absorption rate constant ka, also the transit rate constant ktr (1/h)") # Table 3 'Ka' 0.752 (log), untransformed 2.12 h-1; control stream KTR = KA
    lvc <- log(2340)
    label("Mobocertinib apparent central volume Vc/F in patients (L)") # Table 3 'Vc mobo/F' 7.76 (log), untransformed 2340 L
    lcl <- log(108)
    label("Mobocertinib apparent clearance CL/F in patients at baseline enzyme (L/h)") # Table 3 'CL mobo/F' 4.68 (log), untransformed 108 L/h
    lq <- log(11.7)
    label("Mobocertinib apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q mobo/F' 2.46 (log), untransformed 11.7 L/h
    lvp <- log(1110)
    label("Mobocertinib apparent peripheral volume Vp/F (L)") # Table 3 'Vp mobo/F' 7.01, untransformed 1110 L; control stream TVV3 = EXP(THETA(5))

    # ------------------------------------------------------------------
    # Metabolites. Control stream: FM60 = 0.62, FM14 = 0.08 hardcoded in
    # $PK ('Fraction of M60 and M14 is ~62% and ~8%'); Methods: 'set to 62%
    # and 8%, respectively' based on clinical data. M60 = AP32960 and M14 =
    # AP32914 (the M1 compartment has the peripheral compartment, matching
    # the two-compartment AP32960 model of Figure 1).
    # ------------------------------------------------------------------
    fm_ap32960 <- fixed(0.62)
    label("Fraction of mobocertinib clearance forming AP32960 (molar fraction)") # Methods 'Population PK modeling'; control stream FM60 = 0.62
    fm_ap32914 <- fixed(0.08)
    label("Fraction of mobocertinib clearance forming AP32914 (molar fraction)") # Methods 'Population PK modeling'; control stream FM14 = 0.08
    lcl_ap32960 <- log(117)
    label("AP32960 apparent clearance CL/F in patients at baseline enzyme (L/h)") # Table 3 'CL AP32960/F' 4.76 (log), untransformed 117 L/h
    lvc_ap32960 <- log(12.8)
    label("AP32960 apparent central volume Vc/F (L)") # Table 3 'Vc AP32960/F' 2.55 (log), untransformed 12.8 L
    lq_ap32960 <- log(26.5)
    label("AP32960 apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q AP32960/F' 3.28 (log), untransformed 26.5 L/h
    lvp_ap32960 <- log(1090)
    label("AP32960 apparent peripheral volume Vp/F (L)") # Table 3 'Vp AP32960/F' 6.99, untransformed 1090 L
    lcl_ap32914 <- log(124)
    label("AP32914 apparent clearance CL/F in patients at baseline enzyme (L/h)") # Table 3 'CL AP32914/F' 4.82 (log), untransformed 124 L/h
    lvc_ap32914 <- log(30.6)
    label("AP32914 apparent volume Vc/F (L)") # Table 3 'Vc AP32914/F' 3.42 (log), untransformed 30.6 L

    # ------------------------------------------------------------------
    # Healthy-volunteer effects, linear relative form (1 + theta * HV).
    # ------------------------------------------------------------------
    e_healthy_vc <- 0.787
    label("Linear relative effect of healthy-volunteer status on mobocertinib Vc/F (unitless)") # Table 3 'HV on Vc mobo/F' 0.787; control stream THETA(18)
    e_healthy_cl <- 0.900
    label("Linear relative effect of healthy-volunteer status on mobocertinib CL/F (unitless)") # Table 3 'HV on CL mobo/F' 0.900; control stream THETA(15)
    e_healthy_cl_ap32960 <- 0.738
    label("Linear relative effect of healthy-volunteer status on AP32960 CL/F (unitless)") # Table 3 'HV on CL AP32960/F' 0.738; control stream THETA(17)
    e_healthy_cl_ap32914 <- 0.909
    label("Linear relative effect of healthy-volunteer status on AP32914 CL/F (unitless)") # Table 3 'HV on CL AP32914/F' 0.909; control stream THETA(16)

    # ------------------------------------------------------------------
    # Enzyme auto-induction. Control stream $DES:
    #   CTOT = C2 + C4 + C5; EFF = CTOT*EMAX/(EC50 + CTOT)
    #   DADT(10) = KENZ*(1 + EFF) - KENZ*A(10), A_0(10) = 1
    # ------------------------------------------------------------------
    lkenz <- log(0.00392)
    label("Enzyme synthesis and degradation rate constant kenz (1/h)") # Table 3 'K enz' -5.54 (log), untransformed 0.00392 h-1
    lec50 <- log(0.213)
    label("Molar-sum plasma concentration giving half-maximal enzyme induction EC50 (umol/L)") # Table 3 'EC 50' -1.55 (log), untransformed 213 nM = 0.213 umol/L
    lemax <- log(0.781)
    label("Maximal fractional increase in enzyme production Emax (unitless)") # Table 3 'E max' 0.781 (estimated untransformed; control stream EMAX = THETA(21))

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 3 reports variances (omega^2):
    # the untransformed column is sqrt(omega^2), e.g. sqrt(0.209) = 45.7%,
    # and the covariance rows give correlations, e.g. 0.197 /
    # sqrt(0.237 * 0.246) = 0.814. ETA(1) KA; ETA(2) V2; ETA(3) CL; ETA(4)
    # CLM60; ETA(5) CLM14. Q, Vp, metabolite volumes and enzyme parameters
    # carry no IIV in the control stream.
    # ------------------------------------------------------------------
    etalka ~ 0.209 # Table 3 'IIV K a' 0.209 (45.7%)
    # Block order Vc mobo/F, CL mobo/F, CL AP32960/F, CL AP32914/F. Table 3:
    # variances 0.237, 0.246, 0.157, 0.295; covariances CL-Vc 0.197,
    # CL AP32960-Vc 0.153, CL AP32960-CL 0.182, CL AP32914-Vc 0.194,
    # CL AP32914-CL 0.221, CL AP32914-CL AP32960 0.171.
    etalvc + etalcl + etalcl_ap32960 + etalcl_ap32914 ~ c(
      0.237,
      0.197, 0.246,
      0.153, 0.182, 0.157,
      0.194, 0.221, 0.171, 0.295
    )

    # ------------------------------------------------------------------
    # Residual error. Control stream: Y = LOG(C) + RES*EPS(n) with $SIGMA 1
    # FIX, so THETA(6)-THETA(8) are SDs on the log scale (Table 3 'additive
    # error on log scale'), i.e. log-normal error in linear space.
    # ------------------------------------------------------------------
    expSd <- 0.414
    label("Mobocertinib additive residual SD on the log scale") # Table 3 residual 'Mobocertinib' 0.414; control stream RESCP = THETA(6)
    expSd_ap32960 <- 0.373
    label("AP32960 additive residual SD on the log scale") # Table 3 residual 'AP32960' 0.373; control stream RESM1 = THETA(7)
    expSd_ap32914 <- 0.405
    label("AP32914 additive residual SD on the log scale") # Table 3 residual 'AP32914' 0.405; control stream RESM2 = THETA(8)
  })

  model({
    # Healthy-volunteer covariate terms (control stream $PK)
    cov_healthy_vc <- 1 + e_healthy_vc * DIS_HEALTHY
    cov_healthy_cl <- 1 + e_healthy_cl * DIS_HEALTHY
    cov_healthy_cl_ap32960 <- 1 + e_healthy_cl_ap32960 * DIS_HEALTHY
    cov_healthy_cl_ap32914 <- 1 + e_healthy_cl_ap32914 * DIS_HEALTHY

    # Individual parameters
    ka <- exp(lka + etalka)
    ktr <- ka
    vc <- exp(lvc + etalvc) * cov_healthy_vc
    cl <- exp(lcl + etalcl) * cov_healthy_cl
    q <- exp(lq)
    vp <- exp(lvp)

    cl_ap32960 <- exp(lcl_ap32960 + etalcl_ap32960) * cov_healthy_cl_ap32960
    vc_ap32960 <- exp(lvc_ap32960)
    q_ap32960 <- exp(lq_ap32960)
    vp_ap32960 <- exp(lvp_ap32960)
    cl_ap32914 <- exp(lcl_ap32914 + etalcl_ap32914) * cov_healthy_cl_ap32914
    vc_ap32914 <- exp(lvc_ap32914)

    kenz <- exp(lkenz)
    ec50 <- exp(lec50)
    emax <- exp(lemax)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_ap32960 <- cl_ap32960 / vc_ap32960
    k12_ap32960 <- q_ap32960 / vc_ap32960
    k21_ap32960 <- q_ap32960 / vp_ap32960
    kel_ap32914 <- cl_ap32914 / vc_ap32914

    # Plasma concentrations (umol/L) and their molar sum, which drives
    # enzyme production with equal molar potency for the three moieties.
    Cc <- central / vc
    Cc_ap32960 <- central_ap32960 / vc_ap32960
    Cc_ap32914 <- central_ap32914 / vc_ap32914
    c_molar_sum <- Cc + Cc_ap32960 + Cc_ap32914
    eff_enzyme <- emax * c_molar_sum / (ec50 + c_molar_sum)

    enzyme(0) <- 1

    # Control stream $DES. The relative enzyme amount multiplies all three
    # elimination rates; 1 - fm_ap32960 - fm_ap32914 = 0.30 of mobocertinib
    # clearance leaves the system directly.
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ktr * transit3 - k12 * central + k21 * peripheral1 - kel * central * enzyme
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_ap32960) <- fm_ap32960 * kel * central * enzyme - kel_ap32960 * central_ap32960 * enzyme - k12_ap32960 * central_ap32960 + k21_ap32960 * peripheral1_ap32960
    d/dt(peripheral1_ap32960) <- k12_ap32960 * central_ap32960 - k21_ap32960 * peripheral1_ap32960
    d/dt(central_ap32914) <- fm_ap32914 * kel * central * enzyme - kel_ap32914 * central_ap32914 * enzyme
    d/dt(enzyme) <- kenz * (1 + eff_enzyme) - kenz * enzyme

    Cc ~ lnorm(expSd)
    Cc_ap32960 ~ lnorm(expSd_ap32960)
    Cc_ap32914 ~ lnorm(expSd_ap32914)
  })
}
