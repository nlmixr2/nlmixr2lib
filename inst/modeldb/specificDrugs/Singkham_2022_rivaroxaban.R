Singkham_2022_rivaroxaban <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for oral rivaroxaban in Thai adults with",
    "non-valvular atrial fibrillation (Singkham 2022). Apparent clearance",
    "carries a power effect of Cockcroft-Gault creatinine clearance",
    "normalised to 57.5 mL/min and apparent volume a power effect of body",
    "weight normalised to 63 kg. Between-subject variability on CL/F and ka",
    "only (none on V/F); additive residual error. Concentrations are",
    "rivaroxaban-calibrated anti-factor-Xa activity.",
    sep = " "
  )
  reference <- paste(
    "Singkham N, Phrommintikul A, Pacharasupa P, Norasetthada L, Gunaparn S,",
    "Prasertwitayakij N, Wongcharoen W, Punyawudho B. Population",
    "Pharmacokinetics and Dose Optimization Based on Renal Function of",
    "Rivaroxaban in Thai Patients with Non-Valvular Atrial Fibrillation.",
    "Pharmaceutics. 2022;14(8):1744. doi:10.3390/pharmaceutics14081744.",
    "PMCID: PMC9414338.",
    sep = " "
  )
  vignette <- "Singkham_2022_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault, raw mL/min (NOT BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Singkham 2022 Methods section 2.3: 'creatinine clearance (CrCl;",
        "mL/min, calculated according to Cockcroft and Gault equation)', so",
        "the values are raw mL/min and are NOT BSA-normalised. Enters CL/F",
        "as the power function of the Table 2 footnote c,",
        "CL/F = 4.19 x (CrCl/57.5)^0.277. Cohort mean (SD) 59.0 (22.8)",
        "mL/min (Table 1); patients with CrCl < 15 mL/min were excluded",
        "(Methods section 2.1). The paper does not say which cohort statistic",
        "57.5 mL/min is; it is close to, and presumably, the cohort median.",
        sep = " "
      ),
      source_name = "CrCl"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V/F as the power function of the Table 2 footnote d,",
        "V/F = 37.5 x (WT/63)^0.412. Cohort mean (SD) 64.0 (14.1) kg",
        "(Table 1). The paper does not say which cohort statistic 63 kg is.",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Methods section 2.3) and not retained. Cohort mean (SD) 69.4 (9.2) years (Table 1). Age is a Cockcroft-Gault input, so it enters the model indirectly through CRCL."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not significant on either CL/F or V/F (p > 0.05; Discussion). Cohort median (Q1, Q3) 24.2 (21.5, 26.9) kg/m^2 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened (Methods section 2.3) and not retained; creatinine clearance was the renal descriptor kept. Cohort mean (SD) 1.1 (0.3) mg/dL (Table 1)."
    ),
    SEXF = list(
      description = "Female-sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods section 2.3) and not retained. 38 of 60 patients (63.3%) were male (Table 1). Sex is a Cockcroft-Gault input, so it enters indirectly through CRCL."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    n_observations = 240L,
    age_range = "mean (SD) 69.4 (9.2) years",
    weight_range = "mean (SD) 64.0 (14.1) kg",
    sex_female_pct = 36.7,
    race_ethnicity = "Thai (100%)",
    disease_state = "Adults with non-valvular atrial fibrillation eligible for a direct oral anticoagulant",
    renal_function = "Cockcroft-Gault creatinine clearance mean (SD) 59.0 (22.8) mL/min; CrCl < 15 mL/min excluded",
    dose_range = paste(
      "Rivaroxaban orally once daily: standard dose (20 mg for CrCl >= 50",
      "mL/min, 15 mg for CrCl 15-49 mL/min) for at least one week, then the",
      "Japan-specific dose (15 mg for CrCl >= 50 mL/min, 10 mg for CrCl 15-49",
      "mL/min) for at least one week",
      sep = " "
    ),
    regions = "Thailand (tertiary hospital, Chiang Mai)",
    notes = paste(
      "Table 1 lists baseline characteristics. Each patient contributed",
      "steady-state peak (2-4 h post-dose) and trough (22-24 h post-dose)",
      "samples on each of the two dosing occasions, i.e. four samples per",
      "patient. Concentrations are anti-factor-Xa activity (BIOPHEN Heparin",
      "LRT) calibrated with rivaroxaban-specific calibrators. Concomitant",
      "dronedarone in 3 (5.0%) and amiodarone in 1 (1.7%) patient. Estimation",
      "in NONMEM 7.4 with FOCE-I; ka and its IIV were stabilised with a",
      "frequentist $PRIOR from an earlier Japanese population PK study.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters -- Singkham 2022 Table 2, 'Final Model (NONMEM)'
    # column. The structure is a one-compartment model with first-order
    # absorption and elimination (Results section 3.1).
    # ------------------------------------------------------------------------
    lka <- log(0.697)
    label("First-order absorption rate constant ka (1/h)") # Table 2 'ka (h-1)' = 0.697 [RSE 10.7%]; estimated under a frequentist $PRIOR (Methods section 2.3)
    lcl <- log(4.19)
    label("Apparent clearance CL/F at CrCl 57.5 mL/min (L/h)") # Table 2 'CL/F (L/h)' = 4.19 [RSE 3.8%]
    lvc <- log(37.5)
    label("Apparent volume of distribution V/F at 63 kg (L)") # Table 2 'V/F (L)' = 37.5 [RSE 4.7%]

    # ------------------------------------------------------------------------
    # Covariate effects -- Table 2 footnotes c and d give the power forms
    # CL/F = 4.19 x (CrCl/57.5)^0.277 and V/F = 37.5 x (WT/63)^0.412.
    # ------------------------------------------------------------------------
    e_crcl_cl <- 0.277
    label("Power exponent on (CRCL / 57.5) for CL/F (unitless)") # Table 2 'CrCl on CL/F' = 0.277 [RSE 29%] and footnote c; Results text quotes 0.278
    e_wt_vc <- 0.412
    label("Power exponent on (WT / 63) for V/F (unitless)") # Table 2 'WT on V/F' = 0.412 [RSE 35.7%] and footnote d

    # ------------------------------------------------------------------------
    # Between-subject variability. Table 2 prints %CV; the variances below
    # are omega^2 = log(1 + CV^2), the exact log-normal relation. The Table 2
    # Wald intervals fix that reading: back-transforming
    # omega^2 x (1 +/- 1.96 x RSE) through CV = sqrt(exp(omega^2) - 1) gives
    # 16.66-26.25% for CL/F (printed 16.67-26.24) and 66.37-85.14% for ka
    # (printed 66.39-85.10), whereas the CV = sqrt(omega^2) reading gives
    # 16.75-26.12% and 67.98-83.08%. No IIV on V/F could be estimated
    # (Discussion).
    # ------------------------------------------------------------------------
    etalcl ~ 0.0470137 # Table 2 'IIV of CL/F (%CV)' = 21.94 [RSE 21.3%]; log(1 + 0.2194^2)
    etalka ~ 0.4550377 # Table 2 'IIV of ka (%CV)' = 75.91 [RSE 10.1%]; log(1 + 0.7591^2)

    # ------------------------------------------------------------------------
    # Residual error. Methods section 2.3: 'The residual unexplained
    # variability (RUV) was modeled using an additive function.' The value
    # carries concentration units in Table 2, so it is a standard deviation.
    # ------------------------------------------------------------------------
    addSd <- 0.092
    label("Additive residual error (mg/L)") # Table 2 'RUV, additive (mg/L)' = 0.092 [RSE 11.7%]
  })

  model({
    # Individual PK parameters. Normalising constants 57.5 mL/min and 63 kg
    # are those printed in Table 2 footnotes c and d.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 57.5)^e_crcl_cl
    vc <- exp(lvc) * (WT / 63)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and volume in L give mg/L, the unit of the Table 2 residual
    # error (1 mg/L = 1000 ng/mL).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
