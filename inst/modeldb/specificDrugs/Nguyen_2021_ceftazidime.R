Nguyen_2021_ceftazidime <- function() {
  description <- paste(
    "One-compartment IV population PK model for ceftazidime in 50 Vietnamese",
    "adults hospitalised for acute exacerbations of chronic obstructive",
    "pulmonary disease (Nguyen 2021). Clearance scales as a power of",
    "Cockcroft-Gault creatinine clearance (reference 69.02 mL/min, the",
    "observation-weighted cohort mean), the only covariate retained; the",
    "volume of distribution carries no covariate. Inter-individual",
    "variability is exponential on clearance and volume, and residual error",
    "is proportional. Fitted in Monolix alongside a separate imipenem model",
    "from the same study (Nguyen_2021_imipenem).",
    sep = " "
  )
  reference <- paste(
    "Nguyen TM, Ngo TH, Truong AQ, Vu DH, Le DC, Vu NB, Can TN, Nguyen HA,",
    "Phan TP, Van Bambeke F, Vidaillac C, Ngo QC.",
    "Population pharmacokinetics and dose optimization of ceftazidime and",
    "imipenem in patients with acute exacerbations of chronic obstructive",
    "pulmonary disease.",
    "Pharmaceutics. 2021;13(4):456. doi:10.3390/pharmaceutics13040456.",
    sep = " "
  )
  vignette <- "Nguyen_2021_ceftazidime_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
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
        "Enters CL as the power term (CRCL/69.02)^0.485 (Nguyen 2021 Results",
        "3.3). The paper centres every continuous covariate on its mean",
        "weighted by the number of observations per individual (Methods 2.3,",
        "Equation 2), so 69.02 mL/min is that weighted mean for the",
        "ceftazidime cohort, not its median (Table 1 median 62.9 mL/min, IQR",
        "49.0-76.8). Measured at the time of blood sampling. Stored under the",
        "canonical CRCL column, which accepts raw Cockcroft-Gault mL/min when",
        "the source does not BSA-normalise (precedent: Delattre 2010",
        "amikacin, Nguyen 2021 imipenem)."
      ),
      source_name = "CLCRCG"
    )
  )

  # Screened in the stepwise covariate search and not retained (Nguyen 2021
  # Supplementary Table S2, models PK_C_09 to PK_C_33). AGE and CREAT passed
  # forward selection on CL on their own (dOFV 8.2 and 12.98) but added
  # nothing once CRCL was in the model (PK_C_33: 0.31; PK_C_32: -1.16).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on V (dOFV 0.02) and CL (8.2 alone; 0.31 on top of CRCL), not retained (Table S2). Cohort median 69 years, IQR 63-77 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 2.78) and CL (2.19), not retained (Table S2). Cohort is 47 male / 3 female (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on V (dOFV 0.22) and CL (1.77), not retained (Table S2). Cohort median 51 kg, IQR 47-57 (Table 1)."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      notes = "Screened on V (dOFV 0.68) and CL (4.86), not retained (Table S2). Cohort median 45 kg, IQR 41-47 (Table 1); the estimating equation is not stated."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed as a tested covariate in Methods 2.3 but absent from the Table S2 step list; not retained. Cohort median 19.49 kg/m^2 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "not reported",
      type = "continuous",
      notes = "Screened on V (dOFV 0.21) and CL (12.98 alone; -1.16 on top of CRCL), not retained (Table S2)."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 1.4) and CL (1.28), not retained (Table S2). 9 of 50 patients (Table 1)."
    ),
    MECH_VENT = list(
      description = "Invasive ventilation",
      units = "(binary)",
      type = "binary",
      notes = "Screened on V (dOFV 0.53) and CL (1.37) as 'VENTILATOR', not retained (Table S2). 7 of 50 patients (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    n_concentrations = 97L,
    age_median = "69 years (IQR 63-77)",
    weight_median = "51 kg (IQR 47-57)",
    height_median = "162.5 cm (IQR 160-167)",
    bmi_median = "19.49 kg/m^2 (IQR 17.55-21.44)",
    crcl_median = "62.9 mL/min (Cockcroft-Gault; IQR 49.0-76.8)",
    sex_female_pct = 6,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Adults hospitalised on a respiratory ward for acute exacerbations of",
      "chronic obstructive pulmonary disease (GOLD 2020 definition) and",
      "treated with ceftazidime for at least three consecutive days. 34%",
      "had respiratory distress and 14% invasive ventilation; none was",
      "admitted to intensive care."
    ),
    dose_range = paste(
      "1 g q12h (n = 1), 1 g q8h (n = 39), 1 g q6h (n = 1), 2 g q12h",
      "(n = 3) or 2 g q8h (n = 6) by IV infusion (Table 1); infusion",
      "durations varied between patients (under 30 min to over 120 min,",
      "Supplementary Figure S1)."
    ),
    regions = "Vietnam (Bach Mai Hospital, Hanoi)",
    notes = paste(
      "Prospective study, August 2018 - March 2019. Two sparse samples per",
      "patient: at least 30 min after the end of the third infusion and 1-2",
      "h before the fourth dose (steady state assumed); 7 of the 94 patients",
      "across both drug cohorts gave one sample only. Assayed by HPLC-UV",
      "(LLOQ 2 mg/L). Fitted in Monolix 2019R1 by SAEM; one- and",
      "two-compartment structures with constant, proportional and combined",
      "error were compared by BIC (Table S1). Covariates were added by",
      "forward inclusion (dOFV > 6.635) and backward elimination (dOFV >",
      "10.828). In addition to the entries in covariatesDataExcluded, the",
      "Anthonisen score, respiratory distress ('ARDS' in Table S2) and the",
      "MDRD-4 estimated GFR were screened and not retained; they are noted",
      "here rather than registered because no coefficient is reported."
    )
  )

  ini({
    # ===== Structural PK -- Nguyen 2021 Table 2, CEFTAZIDIME block. Typical
    # values at the reference CRCL of 69.02 mL/min. =====
    lcl <- log(8.74); label("Clearance at CRCL = 69.02 mL/min (L/h)") # Nguyen 2021 Table 2: CL = 8.74 L/h (RSE 3.18%)
    lvc <- log(23.7); label("Volume of distribution (L)") # Nguyen 2021 Table 2: Vd = 23.7 L (RSE 2.96%)

    # ===== Covariate effect -- Nguyen 2021 Results 3.3:
    #   CL_i = 8.74 * (CLCR_i/69.02)^0.485 * exp(eta_CL)
    e_crcl_cl <- 0.485; label("Power exponent on (CRCL/69.02) for CL (unitless)") # Nguyen 2021 Table 2: beta CLCRCG on CL = 0.485 (RSE 17.2%)

    # ===== Inter-individual variability =====
    # Nguyen 2021 Methods 2.3 defines omega_P as 'the standard deviation of
    # eta_Pi' (Monolix convention), so the printed percentage is the SD of
    # eta and the variance is omega^2.
    etalcl ~ 0.043264 # Nguyen 2021 Table 2: omega CL = 20.8% -> 0.208^2
    etalvc ~ 0.0169 # Nguyen 2021 Table 2: omega V = 13% -> 0.13^2

    # ===== Residual error =====
    propSd <- 0.121; label("Proportional residual error (fraction)") # Nguyen 2021 Table 2: b = 12.1%
  })

  model({
    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (CRCL / 69.02)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Intravenous infusion into the central compartment; the infusion
    # duration comes from the event table's rate / dur column.
    d/dt(central) <- -kel * central

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
