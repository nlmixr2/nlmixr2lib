Tanaka_2026_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in hospitalised Japanese children (<18 years; Tanaka 2026). Clearance scales allometrically with body weight (exponent 0.75, reference 70 kg), matures sigmoidally with postmenstrual age (TM50 49.1 weeks, Hill 2.17), follows a power function (exponent 0.423) of maturation-adjusted modified-Schwartz eGFR (eGFR divided by the Rhodin 2009 renal maturation fraction times 114.1 mL/min/1.73 m^2), and is 27% higher on days with a daily maximum body temperature >= 38 degC. Central volume, peripheral volume and intercompartmental clearance are fixed literature values scaled by weight (linear for volumes, 0.75 for Q). Proportional residual error."
  reference <- "Tanaka R, Sawaguchi A, Mizuno T, Irie K, Ono H, Amano E, Goto H, Tashiro S, Tatsuta R, Yoshikawa N, Hirano N, Ihara K, Itoh H. Impact of Febrile State on Vancomycin Clearance in Pediatric Patients: Insights From Population Pharmacokinetic Modeling. Pharm Res. 2026;43(8):2721-2732. doi:10.1007/s11095-026-04145-w"
  vignette <- "Tanaka_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Vancomycin is given by IV infusion into `central`. The Methods state that
  # SERUM vancomycin was assayed (particle-enhanced turbidimetric inhibition
  # immunoassay, Dimension system), so `central` is verified as serum. The
  # paper never names a matrix for the peripheral distribution compartment.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size term on CL and Q (exponent 0.75) and on Vc and Vp (exponent 1), all with reference weight 70 kg (Methods 'Covariate Analysis'; supplement 2 control stream column BBWT). Table I baseline median 16.9 kg (IQR 9.3-28.8).",
      source_name = "BBWT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Postmenstrual age in WEEKS (supplement 2 control stream column PMAW; Table III gives TM50 in weeks). Drives two sigmoid functions: the estimated clearance maturation PAGE^2.17 / (49.1^2.17 + PAGE^2.17) and the fixed Rhodin 2009 renal maturation PAGE^3.4 / (47.7^3.4 + PAGE^3.4) that adjusts eGFR. Computed from gestational age plus postnatal age; where gestational age was missing the source assigned 280 days (40 weeks). Table I baseline median 6.8 years (IQR 2.0-12.2).",
      source_name = "PMAW"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate from the modified (bedside) Schwartz equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = k * height (cm) / serum creatinine (mg/dL), with k = 0.244 for low-birth-weight (<2.5 kg) infants under 1 year, 0.277 for other infants, 0.296 for children aged 1-12 years and 0.413 above 12 years",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL as (CRCL / (Fmat_eGFR * 114.1))^0.423, where Fmat_eGFR = PAGE^3.4 / (47.7^3.4 + PAGE^3.4) is the Rhodin 2009 renal maturation fraction, so the eGFR reference is 114.1 mL/min/1.73 m^2 at full renal maturity and lower in young infants. Patients with eGFR >= 400 mL/min/1.73 m^2 were excluded. Table I baseline median 118.4 mL/min/1.73 m^2 (IQR 96.8-151.9). Laboratory values not measured on the sampling day were taken from the nearest available date.",
      source_name = "eGFR"
    ),
    BODYTEMP_FEBRILE = list(
      description = "Febrile-day indicator: 1 when the daily maximum body temperature is >= 38 degC, 0 when it is < 38 degC",
      units = "(binary)",
      type = "binary",
      reference_category = "daily maximum body temperature < 38 degC",
      notes = "Time-varying (one value per day). Threshold 38.0 degC on the DAILY MAXIMUM axillary temperature (Table II footnote; Table III footnote 'BT38 ... >= 38 degC, flag = 1; < 38 degC, flag = 0'). Source column BT38 in the supplement 2 control stream. Multiplies CL by 1.27 when 1. Chosen over continuous daily maximum temperature and over the febrile-neutropenia indicator, which were collinear with it, because it gave the largest OFV drop (Table II).",
      source_name = "BT38"
    )
  )

  # Screened in the forward-addition step on CL (Table II) but not retained in
  # the final model; documentation only.
  covariatesDataExcluded <- list(
    HGB = list(
      description = "Hemoglobin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL as a normalised power model; dOFV 0 against the base model (Table II). Not retained.",
      source_name = "HGB"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened on CL as a binary indicator (source column GENDER; coding direction not stated); dOFV -0.987 (Table II). Not retained.",
      source_name = "GENDER"
    ),
    BODYTEMP = list(
      description = "Daily maximum body temperature",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL as a normalised power model; dOFV -16.238 (Table II). Its relationship to the CL random effect was biphasic with an inflection at 38 degC (supplement Figure S3), so it was dichotomised into BODYTEMP_FEBRILE instead. Table I baseline median 38.6 degC (IQR 37.7-39.3).",
      source_name = "BT"
    ),
    FEBRILE_NEUTROPENIA = list(
      description = "Febrile neutropenia indicator: axillary temperature >= 37.5 degC with neutrophil count < 500/uL between vancomycin start and the first TDM sample",
      units = "(binary)",
      type = "binary",
      reference_category = "no febrile neutropenia",
      notes = "Screened on CL; dOFV -15.962 alone (Table II). Collinear with BODYTEMP_FEBRILE; added on top of it, dOFV -3.039 with theta_FN 1.14 (95% CI 0.987-1.293), below the retention criterion (supplement Table S2). Not retained. 30.2% of patients met the definition (Table I).",
      source_name = "FN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 129L,
    n_studies = 1L,
    n_sites = 1L,
    n_concentrations = 214L,
    age_range = "<18 years; Table I: <1 month 6.2%, 1-<3 months 5.4%, 3-<12 months 10.1%, 1-<6 years 27.9%, 6-<12 years 25.6%, 12-<18 years 24.8%",
    weight_median = "16.9 kg (IQR 9.3-28.8)",
    height_median = "104.0 cm (IQR 73.5-138.5)",
    pma_median = "6.8 years (IQR 2.0-12.2)",
    sex_female_pct = 47.3,
    race_ethnicity = "Japanese (single Japanese centre; not tabulated)",
    disease_state = "Paediatric inpatients receiving IV vancomycin with therapeutic drug monitoring; 30.2% met the febrile-neutropenia definition, MRSA was cultured in 14.7%. Median daily maximum body temperature 38.6 degC (IQR 37.7-39.3). Renal replacement therapy and eGFR >= 400 mL/min/1.73 m^2 excluded.",
    dose_range = "Clinician-chosen IV regimens (not tabulated)",
    regions = "Japan (Department of Pediatrics, Oita University Hospital, April 2011 to March 2021)",
    renal_function = "Modified-Schwartz eGFR median 118.4 mL/min/1.73 m^2 (IQR 96.8-151.9); <60: 11.6%, 60-<90: 9.3%, 90-<120: 29.5%, 120-<150: 23.3%, 150-<180: 13.2%, 180-<210: 4.7%, >=210: 8.5%",
    notes = "Retrospective single-centre TDM data. 214 concentrations from 129 patients, mostly troughs (supplement Figure S1); 70 patients (54.3%) contributed a single concentration. Data >= 14 days after vancomycin start, and data after a transition from febrile-neutropenia to non-febrile-neutropenia state, were excluded. Fit in NONMEM 7.3.0 with FOCE-I; 1,000-replicate bootstrap."
  )

  ini({
    # Clearance (Table III). Reference subject: 70 kg, fully matured
    # (Fmat -> 1), eGFR equal to 114.1 mL/min/1.73 m^2 at full renal maturity,
    # daily maximum temperature < 38 degC.
    lcl <- log(6.50); label("Clearance for a 70 kg mature afebrile subject at eGFR 114.1 mL/min/1.73 m^2 (L/h)") # Table III theta_CL = 6.50 (RSE 5.0%); control stream THETA(3)

    # Volumes and intercompartmental clearance were fixed to literature values
    # (Le 2013, Kloprogge 2019) because the trough-dominant data could not
    # inform them (Methods 'Population Pharmacokinetic Modeling').
    lvc <- fixed(log(39.9)); label("Central volume for a 70 kg subject (L)") # Table III Vc = 39.9 L/70 kg, no RSE; control stream TVV1 = 39.9*(BBWT/70)
    lq <- fixed(log(3.85)); label("Intercompartmental clearance for a 70 kg subject (L/h)") # Table III Q = 3.85 L/h/70 kg, no RSE; control stream TVQ = 3.85*(BBWT/70)**0.75
    lvp <- fixed(log(37.8)); label("Peripheral volume for a 70 kg subject (L)") # Table III Vp = 37.8 L/70 kg, no RSE; control stream TVV2 = 37.8*(BBWT/70)

    # Allometric exponents fixed by theory (Methods equations for CL, Q and V).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless)") # Methods: CL, Q = theta * (BW/70)^0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/70) for Vc and Vp (unitless)") # Methods: V = theta * (BW/70)

    # Sigmoid postmenstrual-age maturation of CL (estimated).
    ltmat50 <- log(49.1); label("Postmenstrual age at 50% CL maturation, TM50 (weeks)") # Table III theta_TM50 = 49.1 weeks (RSE 8.5%); control stream THETA(2)
    lhill_mat <- log(2.17); label("Hill coefficient of the CL maturation function (unitless)") # Table III theta_Hill = 2.17 (RSE 14.7%); control stream THETA(1)

    # Renal-function and fever effects on CL.
    e_crcl_cl <- 0.423; label("Power exponent on maturation-adjusted eGFR ratio for CL (unitless)") # Table III theta_eGFR = 0.423 (RSE 27.7%); control stream THETA(4)
    e_bodytemp_febrile_cl <- 1.27; label("Multiplicative factor on CL when daily maximum body temperature >= 38 degC (unitless)") # Table III theta_BT38 = 1.27 (RSE 3.8%); control stream THETA(5)

    # IIV on CL: Table III 31.0% = sqrt(omega^2); supplement Table S3 omega^2.
    etalcl ~ 0.0963 # Table S3 'IIV CL' = 0.0963; Table III 31.0% (RSE 19.5%)

    # Proportional residual error (control stream Y = F*(1+ERR(1))).
    propSd <- 0.303; label("Proportional residual SD (fraction)") # Table III 30.3% (RSE 20.0%); Table S3 sigma^2 = 0.0912
  })
  model({
    # Sigmoid maturation of CL on postmenstrual age in weeks.
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    # Rhodin 2009 renal maturation fraction with its published constants
    # (TM50 47.7 weeks, Hill 3.4), fixed inside the eGFR adjustment.
    fmat_egfr <- PAGE^3.4 / (47.7^3.4 + PAGE^3.4)

    # Table III: CL = theta_CL * (BW/70)^0.75 * Fmat *
    #   [eGFR / (Fmat_eGFR * 114.1)]^theta_eGFR * theta_BT38^flag
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * fmat *
      (CRCL / (fmat_egfr * 114.1))^e_crcl_cl *
      e_bodytemp_febrile_cl^BODYTEMP_FEBRILE
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L, so central/vc is mg/L == ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
