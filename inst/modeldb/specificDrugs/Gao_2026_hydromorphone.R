Gao_2026_hydromorphone <- function() {
  description <- paste(
    "Three-compartment population PK model for a 72-hour continuous",
    "intravenous infusion of hydromorphone in critically ill adults supported",
    "by veno-venous (V-V) or veno-arterial (V-A) extracorporeal membrane",
    "oxygenation (ECMO) (Gao 2026). First-order elimination from the central",
    "compartment. Every disposition parameter (CL, V1, Q2, V2, Q3, V3) has its",
    "own typical value in each ECMO mode, estimated in one joint NONMEM fit that",
    "shares a single set of inter-individual variability terms (a CL-V1",
    "covariance block plus a diagonal term on Q2) across both modes. No",
    "covariate was retained. Combined additive and proportional residual error.",
    "Parameters are in L/min, so model time is in minutes."
  )
  reference <- paste(
    "Gao X, Yao L, Zeng F, Zhang H, Wu Y, Huang Y, Hong Z, Shu H, Zhang J,",
    "Yang X, Zou X, Shang Y. Population Pharmacokinetics of Intravenous",
    "Infusion of Hydromorphone in Critically Ill Patients Supported by",
    "Extracorporeal Membrane Oxygenation. Clin Transl Sci. 2026;19(7):e70666.",
    "doi:10.1111/cts.70666. PMC13358379."
  )
  vignette <- "Gao_2026_hydromorphone"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ECMO_VA_STATUS = list(
      description = "Veno-arterial ECMO mode indicator (1 = V-A ECMO, 0 = V-V ECMO); selects the mode-specific typical value of every disposition parameter",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (V-V ECMO; every subject in the cohort was on ECMO)",
      notes = paste(
        "REPLACEMENT rule, not a multiplicative effect: the indicator switches",
        "each of the six typical values (CL, V1, Q2, V2, Q3, V3) between its",
        "V-V and V-A estimate (Gao 2026 Table 2 and its footnote, 'CLVVi =",
        "exp(CLTVV + eta_i) and CLVAi = exp(CLTVA + eta_i)'). Time-fixed per",
        "subject: each patient stayed on one ECMO mode throughout the 72-h",
        "infusion and 12-h washout. 12 V-V and 14 V-A patients (Table 1; per",
        "patient in Table S2). Because the whole cohort was on ECMO, 0 means",
        "V-V ECMO here, not 'no ECMO'; the model carries no information about",
        "patients off ECMO."
      ),
      source_name = "MODE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise forward inclusion / backward elimination on the PK parameters (Gao 2026 'Covariate Analysis Development'); not retained. V-V 55.75 +/- 10.87, V-A 49.21 +/- 12.72 years (Table 1).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as GEND; not retained. Male 11/12 V-V and 6/14 V-A (Table 1).",
      source_name = "GEND"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as BW, and additionally tested as an estimated-exponent power model (WT/65) and as fixed allometric scaling (0.75 on CL/Q, 1 on volumes, 70 kg); neither improved the fit (dOFV < 3.84) and weight was not retained. Dosing was weight-based (0.03 mg/kg/h).",
      source_name = "BW"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "BMI"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. V-V 26.93 +/- 5.22, V-A 27.90 +/- 3.69 g/L (Table 1).",
      source_name = "ALB"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained. The paper's exploratory post hoc exposure comparison (Table S9) found 9.6-19.8% higher geometric-mean exposure above the median ALT, which it reports as hypothesis-generating only.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained (post hoc exposure comparison in Table S10).",
      source_name = "AST"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "BUN"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as CR; not retained.",
      source_name = "CR"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened; not retained. Raw mL/min (Table 1), estimating equation not stated, not BSA-normalized. Patients on continuous renal replacement therapy during ECMO were excluded.",
      source_name = "CRCL"
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units = "degC",
      type = "continuous",
      notes = "Screened as Temp; not retained.",
      source_name = "Temp"
    ),
    TPRO = list(
      description = "Total plasma protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened as TP; not retained.",
      source_name = "TP"
    ),
    GGT = list(
      description = "Gamma-glutamyl transpeptidase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained.",
      source_name = "GGT"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as TBIL; not retained.",
      source_name = "TBIL"
    )
  )

  compartmentData <- list(
    central = list(analyte = "hydromorphone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "hydromorphone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "hydromorphone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    age_range = "18-80 years (eligibility); observed 25-77 years (Table S2)",
    age_median = "V-V 56 years, V-A 49.5 years (Table S6 footnote)",
    weight_range = "45-102 kg (Table S2)",
    weight_median = "V-V 72.5 kg, V-A 65 kg (Table S5 footnote)",
    sex_female_pct = 34.6,
    race_ethnicity = "Chinese (single centre, Wuhan, China); not otherwise reported",
    disease_state = "Critically ill adults on ECMO for at least 48 h: 12 on veno-venous ECMO (viral pneumonia / ARDS, Pneumocystis pneumonia, post-traumatic ARDS, fibrotic interstitial lung disease) and 14 on veno-arterial ECMO (acute myocardial infarction or viral myocarditis with cardiogenic shock) (Table S2).",
    dose_range = "Hydromorphone 0.03 mg/kg/h continuous IV infusion for 72 h.",
    regions = "China",
    ecmo_mode = "V-V 12 (46.2%), V-A 14 (53.8%)",
    hepatic_function = "Severe hepatic dysfunction (Child-Pugh C) excluded. Median AST V-V 51, V-A 195.5 U/L; median ALT V-V 39.5, V-A 97 U/L (Table 1).",
    renal_function = "Creatinine clearance median V-V 51.90, V-A 57.07 mL/min (Table 1); patients on CRRT during ECMO excluded.",
    notes = paste(
      "Union Hospital, Tongji Medical College, Wuhan; enrolled September 2022",
      "to August 2024 (NCT05565495). 484 arterial plasma concentrations (227",
      "V-V, 257 V-A) sampled during the infusion and up to 12 h after it",
      "stopped; LC-MS/MS assay. Patient sex: 9 of 26 female (Table S2).",
      "Centrifugal-pump circuit (Maquet) with 400-500 mL priming volume;",
      "post- to pre-oxygenator concentration ratios were close to 1",
      "throughout, so no circuit sequestration term was modelled. NONMEM",
      "7.5.0, FOCE-I; 1000-sample nonparametric bootstrap and pcVPC. A pooled",
      "sensitivity model with ECMO mode as a covariate on V2 only (Table S4)",
      "is not the final model and is not reproduced here."
    )
  )

  ini({
    # Gao 2026 Table 2 'Final model estimates'. Paper subscripts: VV = veno-venous
    # ECMO, VA = veno-arterial ECMO. Paper Q2/V2 = peripheral compartment 1
    # (canonical q / vp); paper Q3/V3 = peripheral compartment 2 (q2 / vp2).
    lcl_vv <- log(0.955); label("Clearance, V-V ECMO (L/min)") # Table 2 CLVV 0.955 (RSE 9.13%)
    lcl_va <- log(0.657); label("Clearance, V-A ECMO (L/min)") # Table 2 CLVA 0.657 (RSE 6.76%)
    lvc_vv <- log(23.6); label("Central volume V1, V-V ECMO (L)") # Table 2 V1VV 23.6 (RSE 17.3%)
    lvc_va <- log(16.3); label("Central volume V1, V-A ECMO (L)") # Table 2 V1VA 16.3 (RSE 18.3%)
    lq_vv <- log(0.112); label("Intercompartmental clearance Q2 to peripheral 1, V-V ECMO (L/min)") # Table 2 Q2VV 0.112 (RSE 45.6%)
    lq_va <- log(0.075); label("Intercompartmental clearance Q2 to peripheral 1, V-A ECMO (L/min)") # Table 2 Q2VA 0.075 (RSE 57%)
    lvp_vv <- log(112); label("Peripheral volume V2, V-V ECMO (L)") # Table 2 V2VV 112 (RSE 26.2%)
    lvp_va <- log(72.8); label("Peripheral volume V2, V-A ECMO (L)") # Table 2 V2VA 72.8 (RSE 13.9%)
    lq2_vv <- log(1.48); label("Intercompartmental clearance Q3 to peripheral 2, V-V ECMO (L/min)") # Table 2 Q3VV 1.48 (RSE 15.1%)
    lq2_va <- log(0.538); label("Intercompartmental clearance Q3 to peripheral 2, V-A ECMO (L/min)") # Table 2 Q3VA 0.538 (RSE 10.6%)
    lvp2_vv <- log(87); label("Peripheral volume V3, V-V ECMO (L)") # Table 2 V3VV 87 (RSE 8.93%)
    lvp2_va <- log(43); label("Peripheral volume V3, V-A ECMO (L)") # Table 2 V3VA 43 (RSE 11.8%)

    # IIV shared across both modes. Table 2 footnote: 'reported as approximate
    # CV%', i.e. CV = 100 * sqrt(omega^2), so omega^2 = (CV / 100)^2. The
    # CL-V1 covariance 0.0545 is printed on the variance scale.
    etalcl + etalvc ~ c(0.276^2, 0.0545, 0.572^2) # Table 2 'IIV of CL' 27.6, 'Covariance between CL and V1' 0.0545, 'IIV of V1' 57.2
    etalq ~ 1.40^2 # Table 2 'IIV of Q2' 140 (approximate CV%)

    propSd <- 0.185; label("Proportional residual error (fraction)") # Table 2 proportional error 18.5%
    addSd <- 1.7; label("Additive residual error (ng/mL)") # Table 2 additive error 1.7 ng/mL
  })
  model({
    # Mode-specific typical values (Table 2 footnote): ECMO_VA_STATUS = 1
    # selects the V-A estimate, 0 the V-V estimate; the eta is shared.
    cl <- exp(lcl_vv * (1 - ECMO_VA_STATUS) + lcl_va * ECMO_VA_STATUS + etalcl)
    vc <- exp(lvc_vv * (1 - ECMO_VA_STATUS) + lvc_va * ECMO_VA_STATUS + etalvc)
    q <- exp(lq_vv * (1 - ECMO_VA_STATUS) + lq_va * ECMO_VA_STATUS + etalq)
    vp <- exp(lvp_vv * (1 - ECMO_VA_STATUS) + lvp_va * ECMO_VA_STATUS)
    q2 <- exp(lq2_vv * (1 - ECMO_VA_STATUS) + lq2_va * ECMO_VA_STATUS)
    vp2 <- exp(lvp2_vv * (1 - ECMO_VA_STATUS) + lvp2_va * ECMO_VA_STATUS)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg, volume in L: central / vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
