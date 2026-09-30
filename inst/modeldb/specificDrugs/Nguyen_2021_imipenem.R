Nguyen_2021_imipenem <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in 44 Vietnamese",
    "adults hospitalised for acute exacerbations of chronic obstructive",
    "pulmonary disease (Nguyen 2021). Clearance scales as a power of",
    "Cockcroft-Gault creatinine clearance (reference 75.54 mL/min, the",
    "observation-weighted cohort mean), the only covariate retained; the",
    "volume of distribution carries no covariate. Inter-individual",
    "variability is exponential on clearance and volume, and residual error",
    "is proportional. Fitted in Monolix alongside a separate ceftazidime",
    "model from the same study (Nguyen_2021_ceftazidime).",
    sep = " "
  )
  reference <- paste(
    "Nguyen TM, Ngo TH, Truong AQ, Vu DH, Le DC, Vu NB, Can TN, Nguyen HA,",
    "Phan TP, Van Bambeke F, Vidaillac C, Ngo QC.",
    "Population pharmacokinetics and dose optimization of ceftazidime and",
    "imipenem in patients with acute exacerbations of chronic obstructive",
    "pulmonary disease.",
    "Pharmaceutics. 2021;13(4):456. doi:10.3390/pharmaceutics13040456.",
    "This model is also catalogued (as study 10) by Zhang P, Zhao Y, Zhu J,",
    "Yang Y, Liang G, Wang X, Yu Z. Population pharmacokinetics of imipenem",
    "in different populations for individualized dosing: a systematic",
    "review. Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055.",
    sep = " "
  )
  vignette <- "Nguyen_2021_ceftazidime_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated by the Cockcroft-Gault equation",
        "(Nguyen 2021 Methods 2.3, 'CLCRCG'), raw mL/min and NOT",
        "BSA-normalised: Table 1 reports it in mL/min and the paper applies",
        "no body-surface-area correction."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL as the power term (CRCL/75.54)^0.532 (Nguyen 2021 Results",
        "3.3). The paper centres every continuous covariate on its mean",
        "weighted by the number of observations per individual (Methods 2.3,",
        "Equation 2), so 75.54 mL/min is that weighted mean for the imipenem",
        "cohort, not its median (Table 1 median 76.6 mL/min, IQR",
        "57.5-96.6). Measured at the time of blood sampling. Stored under the",
        "canonical CRCL column, which accepts raw Cockcroft-Gault mL/min when",
        "the source does not BSA-normalise (precedent: Delattre 2010",
        "amikacin, Bai 2024 imipenem, Wang 2024 imipenem)."
      ),
      source_name = "CLCRCG"
    )
  )

  # Screened in the stepwise covariate search and not retained (Nguyen 2021
  # Supplementary Table S2, models PK_I_09 to PK_I_33). AGE passed forward
  # selection on CL on its own (dOFV 8.23) but added only 3.35 once CRCL was
  # in the model (PK_I_31).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on V (dOFV -0.12) and CL (8.23 alone; 3.35 on top of CRCL), not retained (Table S2). Cohort median 65 years, IQR 60-72 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 0.71) and CL (-0.09), not retained (Table S2). Cohort is 41 male / 3 female (Table 1), so the term is near-unidentifiable in this study."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on V (dOFV 4.11) and CL (-0.1), not retained (Table S2). Cohort median 50 kg, IQR 47-55 (Table 1)."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = "Screened on V (dOFV 5.5) and CL (-0.38), not retained (Table S2). Cohort median 43 kg, IQR 40-46 (Table 1); the estimating equation is not stated."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed as a tested covariate in Methods 2.3 but absent from the Table S2 step list; not retained. Cohort median 19.51 kg/m^2 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "not reported",
      type = "continuous",
      notes = "Screened on V (dOFV -0.47) and CL (6.37), not retained (Table S2)."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 1.42) and CL (3.73), not retained (Table S2). 8 of 44 patients (Table 1)."
    ),
    MECH_VENT = list(
      description = "Invasive ventilation",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 2.7) and CL (-0.11) as 'VENTILATOR', not retained (Table S2). 13 of 44 patients (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    n_concentrations = 84L,
    age_median = "65 years (IQR 60-72)",
    weight_median = "50 kg (IQR 47-55)",
    height_median = "160 cm (IQR 159-165)",
    bmi_median = "19.51 kg/m^2 (IQR as printed 18.22-19.51)",
    crcl_median = "76.6 mL/min (Cockcroft-Gault; IQR 57.5-96.6)",
    sex_female_pct = 6.8,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Adults hospitalised on a respiratory ward for acute exacerbations of",
      "chronic obstructive pulmonary disease (GOLD 2020 definition) and",
      "treated with imipenem-cilastatin for at least three consecutive",
      "days. 66% had respiratory distress and 30% invasive ventilation; none",
      "was admitted to intensive care."
    ),
    dose_range = paste(
      "Imipenem 0.5 g q12h (n = 1), 0.5 g q8h (n = 4), 0.5 g q6h (n = 7),",
      "1 g q12h (n = 8) or 1 g q8h (n = 24) by IV infusion (Table 1);",
      "infusion durations varied between patients (under 30 min to over",
      "120 min, Supplementary Figure S1)."
    ),
    regions = "Vietnam (Bach Mai Hospital, Hanoi)",
    notes = paste(
      "Prospective study, August 2018 - March 2019. Two sparse samples per",
      "patient: at least 30 min after the end of the third infusion and 1-2",
      "h before the fourth dose (steady state assumed); 7 of the 94 patients",
      "across both drug cohorts gave one sample only. Plasma was stabilised",
      "1:1 in MOPS buffer and assayed by HPLC-UV (LLOQ 0.5 mg/L). Fitted in",
      "Monolix 2019R1 by SAEM; one- and two-compartment structures with",
      "constant, proportional and combined error were compared by BIC",
      "(Table S1). Covariates were added by forward inclusion (dOFV >",
      "6.635) and backward elimination (dOFV > 10.828). In addition to the",
      "entries in covariatesDataExcluded, the Anthonisen score (dOFV 7.28",
      "alone on CL, 3.35 on top of CRCL), respiratory distress ('ARDS' in",
      "Table S2) and the MDRD-4 estimated GFR (7.84 alone, 0.02 on top of",
      "CRCL) were screened and not retained; they are noted here rather than",
      "registered because no coefficient is reported."
    )
  )

  ini({
    # ===== Structural PK -- Nguyen 2021 Table 2, IMIPENEM block. Typical
    # values at the reference CRCL of 75.54 mL/min. =====
    lcl <- log(7.88); label("Clearance at CRCL = 75.54 mL/min (L/h)") # Nguyen 2021 Table 2: CL = 7.88 L/h (RSE 5.35%)
    lvc <- log(15.1); label("Volume of distribution (L)") # Nguyen 2021 Table 2: Vd = 15.1 L (RSE 6.07%)

    # ===== Covariate effect -- Nguyen 2021 Results 3.3:
    #   CL_i = 7.88 * (CLCR_i/75.54)^0.532 * exp(eta_CL)
    e_crcl_cl <- 0.532; label("Power exponent on (CRCL/75.54) for CL (unitless)") # Nguyen 2021 Table 2: beta CLCRCG on CL = 0.532 (RSE 27.2%)

    # ===== Inter-individual variability =====
    # Nguyen 2021 Methods 2.3 defines omega_P as 'the standard deviation of
    # eta_Pi' (Monolix convention), so the printed percentage is the SD of
    # eta and the variance is omega^2. Results 3.3 prose gives 12.9% (V) and
    # 30% (CL); Table 2 values are used (see the vignette Errata).
    etalcl ~ 0.086436 # Nguyen 2021 Table 2: omega CL = 29.4% -> 0.294^2
    etalvc ~ 0.011449 # Nguyen 2021 Table 2: omega V = 10.7% -> 0.107^2

    # ===== Residual error =====
    propSd <- 0.233; label("Proportional residual error (fraction)") # Nguyen 2021 Table 2: b = 23.3%
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (CRCL / 75.54)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central) <- -kel * central

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
