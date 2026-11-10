Zhang_2024_vancomycin <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous vancomycin in",
    "overweight and obese hospitalised adults (ward and ICU) with varying",
    "renal function, pooled with morbidly obese otherwise healthy subjects",
    "undergoing bariatric surgery (Zhang 2024). Clearance scales with",
    "indexed CKD-EPI eGFR by a power function, linearly with total body",
    "weight, and is 15.5% lower in ICU patients. Residual error differs",
    "between the hospital-patient data (additive) and the AMIGO bariatric",
    "surgery data (combined), selected by STUDY_AMIGO."
  )
  reference <- paste(
    "Zhang T, Krekels EHJ, Smit C, van Dongen EPA, Bruggemann RJM, Knibbe CAJ.",
    "How to Dose Vancomycin in Overweight and Obese Patients with Varying",
    "Renal (Dys)function in the Novel Era of AUC 400-600 mg*h/L-Targeted",
    "Dosing. Clin Pharmacokinet. 2024;63(1):79-91. doi:10.1007/s40262-023-01324-5."
  )
  vignette <- "Zhang_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Indexed CKD-EPI estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Creatinine-based CKD-EPI (2009) eGFR, BSA-indexed (Supplementary",
        "Table S1). Power covariate on CL, (CRCL / 75)^0.658 (Table 2).",
        "Time-varying: serum creatinine was measured repeatedly during",
        "treatment (median 11 per patient, Table 1). Selected over",
        "Cockcroft-Gault (TBW or LBW), MDRD and the de-indexed CKD-EPI and",
        "MDRD variants (Results)."
      ),
      source_name = "CKD-EPI"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear covariate on CL, (1 + 0.0106 * (WT - 90)) (Table 2",
        "equation). Selected over LBW, ABW, IBW and BMI (Results)."
      ),
      source_name = "TBW"
    ),
    DIS_CRITILL = list(
      description = "ICU admission status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ward patient or bariatric-surgery subject)",
      notes = paste(
        "1 = admitted to the intensive care unit, 0 = not admitted (Table 2",
        "footnote a: FICU = 1 for patients not admitted to ICU). Encoded as",
        "0.845^DIS_CRITILL on CL, i.e. CL 15.5% lower in ICU patients."
      ),
      source_name = "ICU admission status (FICU)"
    ),
    STUDY_AMIGO = list(
      description = "AMIGO bariatric-surgery dataset indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (hospital-patient dataset)",
      notes = paste(
        "1 = record from the AMIGO trial (20 morbidly obese otherwise",
        "healthy subjects undergoing bariatric surgery, rich sampling after",
        "a single 12.5 mg/kg dose; Smit 2020), Table 2 'HS data'; 0 = record",
        "from the St. Antonius Hospital overweight and obese patient TDM",
        "dataset, Table 2 'PA data'. Selects the residual-error model only."
      ),
      source_name = "PA / HS dataset (Table 2 residual-error rows)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 230,
    n_studies = 2,
    age_range = "23-92 years",
    age_median = "70 years",
    weight_range = "58.6-234.6 kg",
    weight_median = "91.8 kg",
    sex_female_pct = 39.1,
    race_ethnicity = NULL,
    disease_state = paste(
      "210 overweight or obese (BMI >= 25 kg/m^2) hospitalised patients",
      "receiving vancomycin for infection (157 ward, 53 ICU; renal",
      "replacement therapy excluded) and 20 morbidly obese (BMI >= 35",
      "kg/m^2) otherwise healthy subjects with eGFR >= 60 mL/min/1.73 m^2",
      "undergoing bariatric surgery (AMIGO trial)."
    ),
    dose_range = paste(
      "Patients: continuous infusion or intermittent dosing per routine",
      "TDM-guided practice, median dose 1250 mg (250-3500 mg). AMIGO:",
      "single IV infusion of 12.5 mg/kg (maximum 2500 mg)."
    ),
    regions = "The Netherlands (St. Antonius Hospital, Nieuwegein)",
    renal_function = "Indexed CKD-EPI median 74.5 (6.0-154.9) mL/min/1.73 m^2",
    bmi_range = "25.0-65.7 kg/m^2",
    notes = paste(
      "Demographics from Zhang 2024 Table 1 (140 male / 90 female). 1395",
      "concentrations (1188 TDM samples: 11 peaks, 416 troughs, 761 random",
      "samples under continuous infusion; 207 rich samples from AMIGO).",
      "NONMEM 7.4, ADVAN11 TRANS4, FOCE-I."
    )
  )

  ini({
    # Structural parameters: Zhang 2024 Table 2, final model column.
    lcl <- log(3.36); label("Clearance at CKD-EPI 75 mL/min/1.73 m^2 and 90 kg, non-ICU (L/h)") # Table 2: TVCL = 3.36 L/h (RSE 2.60%)
    lvc <- log(16.7); label("Central volume of distribution V1 (L)") # Table 2: V1 = 16.7 L (RSE 14.0%)
    lq <- log(16.1); label("Intercompartmental clearance V1-V2 (L/h)") # Table 2: Q2 = 16.1 L/h (RSE 11.5%)
    lvp <- log(23.1); label("Peripheral volume V2 (L)") # Table 2: V2 = 23.1 L (RSE 3.60%)
    lq2 <- log(2.44); label("Intercompartmental clearance V1-V3 (L/h)") # Table 2: Q3 = 2.44 L/h (RSE 8.00%)
    lvp2 <- log(83.0); label("Second peripheral volume V3 (L)") # Table 2: V3 = 83.0 L (RSE 4.50%)

    # Covariate effects on CL (Table 2 equation):
    # CL = TVCL * (CKDEPI/75)^theta1 * [1 + theta2 * (TBW - 90)] * FICU
    e_crcl_cl <- 0.658; label("Power exponent of (CRCL/75) on CL (unitless)") # Table 2: theta1 = 0.658 (RSE 3.10%)
    e_wt_cl <- 0.0106; label("Linear slope of (WT - 90) on CL (1/kg)") # Table 2: theta2 = 0.0106 (RSE 13.8%)
    e_dis_critill_cl <- 0.845; label("Multiplicative factor on CL for ICU patients (unitless)") # Table 2: FICU = 0.845 (RSE 2.30%)

    # IIV: Table 2 reports variances on the log scale.
    etalcl ~ 0.0730 # Table 2: IIV CL variance = 0.0730 (RSE 10.5%)
    etalvc ~ 0.978 # Table 2: IIV V1 variance = 0.978 (RSE 29.1%)

    # Residual error: Table 2 reports variances; SD = sqrt(variance).
    addSd <- 3.0529; label("Additive residual error, patient data (mg/L)") # Table 2: PA additive error variance = 9.32; sqrt(9.32) = 3.0529
    propSd_amigo <- 0.046476; label("Proportional residual error, AMIGO data (fraction)") # Table 2: HS proportional error variance = 0.00216; sqrt(0.00216) = 0.046476
    addSd_amigo <- 1.0488; label("Additive residual error, AMIGO data (mg/L)") # Table 2: HS additive error variance = 1.10; sqrt(1.10) = 1.0488
  })

  model({
    # Individual parameters (Table 2 equation); DIS_CRITILL = 0 gives FICU = 1.
    cl <- exp(lcl + etalcl) * (CRCL / 75)^e_crcl_cl * (1 + e_wt_cl * (WT - 90)) *
      e_dis_critill_cl^DIS_CRITILL
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc

    # Dataset-specific residual error: additive for the patient (PA) data,
    # combined additive + proportional for the AMIGO (HS) data.
    addSd_study <- addSd * (1 - STUDY_AMIGO) + addSd_amigo * STUDY_AMIGO
    propSd_study <- propSd_amigo * STUDY_AMIGO
    Cc ~ add(addSd_study) + prop(propSd_study)
  })
}
