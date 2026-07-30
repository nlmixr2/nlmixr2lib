Krens_2020_ganciclovir <- function() {
  description <- paste(
    "One-compartment population PK model for IV ganciclovir in critically ill",
    "adults (Krens 2020), with a CKD-EPI eGFR power effect on clearance and a",
    "fixed central volume.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Krens 2020 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Krens SD, Hodiamont CJ, Juffermans NP, Mathot RAA, van Hest RM. Population",
    "Pharmacokinetics of Ganciclovir in Critically Ill Patients.",
    "Ther Drug Monit. 2020;42(2):295-301. doi:10.1097/FTD.0000000000000689.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "CKD-EPI-estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Renal function estimated with the Chronic Kidney Disease Epidemiology",
        "Collaboration (CKD-EPI) equation, which returns mL/min/1.73 m^2. Power",
        "effect on CL referenced to 65 mL/min/1.73 m^2. CKD-EPI was the only",
        "covariate retained; weight, sex, age, ideal body weight, ethnicity,",
        "comorbidity, continuous venovenous hemofiltration (CVVH), serum creatinine",
        "and serum sodium were tested and not retained (Yang 2023 Table 4).",
        "Body weight is therefore NOT a covariate in this model -- Vc is a fixed",
        "42 L for every subject regardless of size.",
        sep = " "
      ),
      source_name        = "CKD-EPI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    n_observations = 128L,
    age_median     = "56 years (range 30-82)",
    weight_median  = "70 kg (range 44-140)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported; ethnicity was tested as a covariate and not retained.",
    disease_state  = "Critically ill (intensive-care) patients receiving ganciclovir.",
    dose_range     = "IV ganciclovir 2.8 mg/kg/day (range 0.7-20 mg/kg/day).",
    regions        = "Netherlands (retrospective study).",
    bioassay       = "HPLC, LLOQ 0.5 mg/L; LC-MS/MS, LLOQ 0.1 mg/L.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Sampling strategy not",
      "reported. Yang 2023 notes that this model's simulated concentration-time",
      "profiles showed high variability at 5 mg/kg q12h and differed from the",
      "other adult models, attributed to the critically ill population in whom PK",
      "is complicated by hemodynamic instability, varied volumes of distribution",
      "and fluctuating renal function. The simulation target in the primary study",
      "was a trough concentration above 1.5 mg/L.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Krens et al. (2020) row. Reference
    # subject: CKD-EPI eGFR = 65 mL/min/1.73 m^2. Clearance in L/h, volume in L.
    # One-compartment IV model; no absorption parameters and no weight scaling.
    lcl <- log(2.3); label("Clearance at CKD-EPI = 65 mL/min/1.73 m^2 (CL, L/h)") # Yang 2023 Table 3 (Krens 2020): CL = 2.3 * (CKD-EPI/65)^0.71
    lvc <- log(42) ; label("Central volume of distribution (Vc, L; no covariates)") # Yang 2023 Table 3 (Krens 2020): Vc = 42

    # Covariate effect. The exponent 0.71 is a non-canonical estimated value.
    e_crcl_cl <- 0.71; label("Power exponent of CKD-EPI eGFR on CL (unitless; reference 65 mL/min/1.73 m^2)") # Yang 2023 Table 3 (Krens 2020): (CKD-EPI/65)^0.71

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2.
    etalcl ~ 0.2209  # Yang 2023 Table 3 (Krens 2020): BSV CL = 47.0% -> 0.470^2
    etalvc ~ 0.64    # Yang 2023 Table 3 (Krens 2020): BSV Vc = 80.0% -> 0.800^2

    # Residual unexplained variability: proportional.
    propSd <- 0.43; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Krens 2020): 43% proportional error
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 65)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
