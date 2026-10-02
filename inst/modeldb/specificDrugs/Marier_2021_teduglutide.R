Marier_2021_teduglutide <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and an absorption lag time for the GLP-2 analog teduglutide in adult and pediatric patients with short bowel syndrome (SBS) and in non-SBS subjects (healthy volunteers and subjects with renal or hepatic impairment), pooled from 17 studies (Marier 2021). Power effects of body weight on Ka, CL/F and Vc/F; power effects of capped creatinine clearance on CL/F and of age on Vc/F; categorical effects of disease status (non-SBS) and sex on CL/F; injection site (non-abdomen) on Ka, lag time and relative bioavailability; formulation strength and supra-therapeutic dose on lag time."
  reference <- "Marier JF, Jomphe C, Peyret T, Wang Y. Population pharmacokinetics and exposure-response analyses of teduglutide in adult and pediatric patients with short bowel syndrome. Clin Transl Sci. 2021;14(6):2497-2509. doi:10.1111/cts.13117"
  vignette <- "Marier_2021_teduglutide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds. Verified against Marier 2021: teduglutide is
  # given by subcutaneous injection into the abdomen, arm or thigh (Methods,
  # 'Data sources and clinical studies'; Table S4 'Site of injection') and
  # measured in plasma (Results, 'Population PK analysis of teduglutide';
  # Figure 1 'plasma concentration').
  compartmentData <- list(
    depot = list(analyte = "teduglutide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "teduglutide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL/F (exponent 0.488), Vc/F (1.35) and Ka (-0.798), each normalized to 70 kg (Table S8; supplement control stream '(WT/70)**THETA'). The estimated (not fixed allometric) exponents were retained because the pooled cohort spans 5.15-127 kg (Methods, 'Population PK analysis of teduglutide'); Tables S11-S13 show that fixed 0.75/1 exponents with or without a renal maturation function did not improve the fit.",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated creatinine clearance, raw mL/min (NOT BSA-normalized): Cockcroft-Gault for subjects > 12 years and the modified Schwartz equation 'normalized to weight' for subjects < 12 years",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as the power function (min(CRCL, 150) / 99.35)^0.341. The source column is CRCLT, 'estimated creatinine clearance rate capped at 150 mL/min' (Table S8 footnote); the cap is applied inside model() so the user supplies the uncapped estimate. The reference 99.35 mL/min is the overall cohort median (Results, 'Baseline characteristics': 99.4 mL/min; Table S5). Estimating equations per Methods, 'Population PK analysis of teduglutide'.",
      source_name = "CRCLT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc/F, (AGE / 34)^-0.312 (Table S8; control stream 'VAGE = ((AGE/34)**THETA(18))'). The reference 34 years is close to the cohort median of 34.5 years (Table S5). Cohort range 0.380-80.0 years; the Discussion quotes the resulting Vc/F of 137 L at 0.38 years and 26.0 L at 80 years for an otherwise-typical subject (weight held at 70 kg).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "CL/F x 0.932 in females (Table S8). The control stream codes SEXN = 1 (male, 'Most common') as the reference and SEXN = 0 as '(1 + THETA(19))', so SEXF = 1 - SEXN and no sign inversion is needed.",
      source_name = "SEXN"
    ),
    DIS_SBS = list(
      description = "Short bowel syndrome patient indicator: 1 = patient with SBS, 0 = non-SBS subject (healthy volunteer or subject with renal or hepatic impairment)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (SBS patient)",
      notes = "CL/F x 0.668 in non-SBS subjects (Table S8 'x 0.668 if non-SBS'; Results: 'The typical CL/F in healthy subjects was ~33% lower than in patients with SBS'). The control stream codes POP1 = 1 as SBS (reference) and POP1 = 0 as '(1 + THETA(15))', so DIS_SBS = POP1. The typical values of the model are therefore those of an SBS patient.",
      source_name = "POP1"
    ),
    INJSITE_ARM = list(
      description = "SC injection-site indicator: 1 = arm, 0 = abdomen",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (abdomen)",
      notes = "Pooled with INJSITE_THIGH into a single non-abdomen indicator (control stream EXLOC1 = 0), which multiplies Ka by 0.766, the lag time by 1.458 and the relative bioavailability by 0.936 (Table S8). Per-dose-record covariate; the site of injection was tested as time-varying (Results). The 24 SBS patients whose injection site was not recorded (Table S4 'Missing') have no stated coding.",
      source_name = "EXLOC1"
    ),
    INJSITE_THIGH = list(
      description = "SC injection-site indicator: 1 = thigh, 0 = abdomen",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (abdomen)",
      notes = "Pooled with INJSITE_ARM into the single non-abdomen stratum; see the INJSITE_ARM notes. Both indicators 0 denotes the abdomen reference.",
      source_name = "EXLOC1"
    ),
    FORM_TEDUGLUTIDE_GE10MGVIAL = list(
      description = "Teduglutide formulation-strength indicator: 1 = vial strength of 10 mg/vial or higher, 0 = below 10 mg/vial",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (1.25, 2.5 or 5 mg/vial)",
      notes = "Absorption lag time x 0.476 for formulation strength >= 10 mg/vial (Table S8). The control stream carries a four-level STRENGTH2 variable whose levels 1, 2 and 3 all take the same '(1 + THETA(12))' factor (the run description reads 'Combined Strength cat'), so the four strength categories of the forward-inclusion step (Table S7 step 1) were collapsed to two. SBS patients received 1.25-10 mg/vials and non-SBS subjects 5-65 mg/vials (Results; Table S4).",
      source_name = "STRENGTH2"
    ),
    DOSE_HIGH = list(
      description = "Supra-therapeutic dose indicator: 1 = supra-therapeutic dose, 0 = therapeutic dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (therapeutic dose)",
      notes = "Absorption lag time x 1.783 for the supra-therapeutic dose level (Table S8). Methods list the extrinsic covariate as 'dose (therapeutic vs. supra-therapeutic [20 mg])', i.e. the 20 mg supra-therapeutic arm of the thorough-QT study C09-001 (Table S1). The paper does not say whether the other >= 20 mg single- and multiple-ascending-dose arms of the phase 1 studies were also flagged; code DOSE_HIGH = 1 only for a 20 mg supra-therapeutic dose unless better information is available.",
      source_name = "SUPRA"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese vs non-Japanese race",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F, Vc/F and Ka and not retained (Results: 'Race (Japanese vs. non-Japanese) did not have an impact on CL/F ... V/F ... K a')."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F and Vc/F and not retained (Results; Table S7)."
    ),
    ADA_POS = list(
      description = "Anti-drug antibody status",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a time-varying covariate and not retained (Results; Table S7)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 478L,
    n_studies = 17L,
    age_range = "0.380 - 80.0 years",
    age_median = "34.5 years",
    weight_range = "5.15 - 127 kg",
    weight_median = "65.5 kg",
    sex_female_pct = 36.8,
    race_ethnicity = c(White = 81.6, Black = 9.0, Japanese = 5.2, Asian = 1.3, Other = 2.9),
    disease_state = "219 patients with short bowel syndrome dependent on parenteral support (194 non-Japanese: 106 adults and 88 pediatric; 25 Japanese: 14 adults and 11 pediatric; 5 patients < 1 year, 86 aged 1-11 years, 8 aged 12-17 years, 120 adults) and 259 non-SBS subjects (healthy volunteers and subjects with renal or moderate hepatic impairment).",
    dose_range = "Subcutaneous teduglutide 0.0125-0.15 mg/kg once daily in SBS patients (all Japanese patients 0.05 mg/kg); single and multiple fixed doses up to 80 mg in phase 1 studies; formulation strengths 1.25-65 mg/vial.",
    renal_function = "Creatinine clearance 11.4-251 mL/min, median 99.4 mL/min. 346 (72.4%) normal, 78 (16.3%) mild, 40 (8.4%) moderate and 7 (1.5%) severe renal impairment, 6 (1.3%) ESRD (Table S4).",
    regions = "Multinational, including Japan.",
    notes = "Baseline characteristics from Tables S4 and S5. 6775 plasma samples, 768 (10.2%) below the limit of quantitation and handled with the M3 method. NONMEM 7.4.3, Laplacian estimation with interaction. The final model (supplement control stream, run7a) carries a full 3 x 3 OMEGA block on CL/F, Vc/F and Ka, but only the variances are reported (Table S8)."
  )

  ini({
    # ---- Structural PK parameters (Marier 2021 Table S8) ----
    # Reference subject: SBS patient (DIS_SBS = 1), male, 70 kg, age 34 years,
    # capped creatinine clearance 99.35 mL/min, abdominal injection, vial
    # strength < 10 mg/vial, therapeutic dose.
    lka <- log(0.330); label("First-order SC absorption rate constant (1/h)") # Table S8: Ka = 0.330 1/h (RSE 5.6%)
    lcl <- log(16.0); label("Apparent clearance CL/F (L/h)") # Table S8: CL/F = 16.0 L/h (RSE 18.1%)
    lvc <- log(33.9); label("Apparent volume of distribution Vc/F (L)") # Table S8: Vc/F = 33.9 L (RSE 6.2%)
    ltlag <- log(0.299); label("Absorption lag time (h)") # Table S8: ALAG = 0.299 h (RSE 40.5%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability for abdominal injection (unitless)") # Table S8: F1 = 1, Fixed

    # ---- Covariate effects (Marier 2021 Table S8) ----
    e_wt_cl <- 0.488; label("Power exponent of WT/70 on CL/F (unitless)") # Table S8: (body weight/70)^0.488 (RSE 18.1%)
    e_crcl_cl <- 0.341; label("Power exponent of capped CRCL/99.35 on CL/F (unitless)") # Table S8: (CRCLT/99.35)^0.341 (RSE 27.6%)
    e_nonsbs_cl <- 0.668; label("Multiplier on CL/F for non-SBS subjects (unitless)") # Table S8: x 0.668 if non-SBS (RSE 41.0%)
    e_sexf_cl <- 0.932; label("Multiplier on CL/F for female subjects (unitless)") # Table S8: x 0.932 if female (RSE 65.1%)
    e_wt_vc <- 1.35; label("Power exponent of WT/70 on Vc/F (unitless)") # Table S8: (body weight/70)^1.35 (RSE 8.1%)
    e_age_vc <- -0.312; label("Power exponent of AGE/34 on Vc/F (unitless)") # Table S8: (age/34.0)^-0.312 (RSE 15.3%)
    e_wt_ka <- -0.798; label("Power exponent of WT/70 on Ka (unitless)") # Table S8: (body weight/70)^-0.798 (RSE 17.7%)
    e_injsite_ka <- 0.766; label("Multiplier on Ka for injection outside the abdomen (unitless)") # Table S8: x 0.766 for SC administration other than abdomen (RSE 91.0%)
    e_injsite_tlag <- 1.458; label("Multiplier on lag time for injection outside the abdomen (unitless)") # Table S8: x 1.458 for SC administration other than abdomen (RSE 105%)
    e_form_tlag <- 0.476; label("Multiplier on lag time for vial strength >= 10 mg/vial (unitless)") # Table S8: x 0.476 for formulation strength >= 10 mg/vial (RSE 9.7%)
    e_dosehigh_tlag <- 1.783; label("Multiplier on lag time for the supra-therapeutic dose (unitless)") # Table S8: x 1.783 for supra-therapeutic dose level (RSE 89.1%)
    e_injsite_fdepot <- 0.936; label("Multiplier on relative bioavailability for injection outside the abdomen (unitless)") # Table S8: x 0.936 for SC administration other than abdomen (RSE 129%)

    # ---- IIV (Marier 2021 Table S8 'BSV') ----
    # Exponential random effects (Methods, 'Population PK analysis of
    # teduglutide'), BSV reported as %CV, so omega^2 = log(CV^2 + 1). The
    # control stream estimates a full OMEGA BLOCK(3) on CL, V and KA but the
    # covariances are not reported, so the block is encoded as diagonal. The
    # lag time carries no IIV (control stream '$OMEGA 0 FIX ; IIV ALAG1';
    # Table S8 '0, fixed').
    etalcl ~ 0.04769 # Table S8 BSV CL/F = 22.1% -> log(0.221^2 + 1)
    etalvc ~ 0.08895 # Table S8 BSV Vc/F = 30.5% -> log(0.305^2 + 1)
    etalka ~ 0.05111 # Table S8 BSV Ka = 22.9% -> log(0.229^2 + 1)

    # ---- Residual error (Marier 2021 Table S8 'Error model') ----
    # Combined proportional + additive on the linear scale (control stream
    # $ERROR: SD = SQRT((IPRED * THETA(8)/100)^2 + THETA(9)^2)).
    propSd <- 0.243; label("Proportional residual error (fraction)") # Table S8: proportional error = 24.3% (RSE 4.8%)
    addSd <- 6.51; label("Additive residual error (ng/mL)") # Table S8: additive error = 6.51 ng/mL (RSE 9.8%)
  })

  model({
    # ---- 1. Derived covariate terms ----
    # Non-abdominal injection site: arm and thigh share one stratum (control
    # stream EXLOC1 = 0). The two indicators are mutually exclusive.
    injsite_nonabd <- INJSITE_ARM + INJSITE_THIGH
    # Creatinine clearance capped at 150 mL/min (Table S8 footnote 'CRCLT').
    crcl_capped <- min(CRCL, 150)

    # ---- 2. Individual PK parameters ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (crcl_capped / 99.35)^e_crcl_cl *
      e_nonsbs_cl^(1 - DIS_SBS) * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (AGE / 34)^e_age_vc
    ka <- exp(lka + etalka) * (WT / 70)^e_wt_ka * e_injsite_ka^injsite_nonabd
    tlag <- exp(ltlag) * e_injsite_tlag^injsite_nonabd *
      e_form_tlag^FORM_TEDUGLUTIDE_GE10MGVIAL * e_dosehigh_tlag^DOSE_HIGH
    fdepot <- exp(lfdepot) * e_injsite_fdepot^injsite_nonabd

    # ---- 3. Micro-constants ----
    kel <- cl / vc

    # ---- 4. ODE system (control stream $DES) ----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 5. Bioavailability and lag time ----
    f(depot) <- fdepot
    alag(depot) <- tlag

    # ---- 6. Observation and error ----
    # central / vc is mg/L; x 1000 gives ng/mL, the units of addSd.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
