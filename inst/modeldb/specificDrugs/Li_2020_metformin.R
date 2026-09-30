Li_2020_metformin <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order oral absorption and an absorption-lag time for immediate-release metformin hydrochloride at steady state in Chinese adults with type 2 diabetes mellitus (Li 2020). Apparent oral clearance scales as a power function of body weight (reference 75 kg) and of CKD-EPI estimated glomerular filtration rate (reference 102.5 mL/min/1.73 m^2); between-subject variability is estimated on CL/F only. OCT1, OCT2 and MATE1 polymorphisms were screened and not retained."
  reference <- "Li L, Guan Z, Li R, Zhao W, Hao G, Yan Y, Xu Y, Liao L, Wang H, Gao L, Wu K, Gao Y, Li Y. Population pharmacokinetics and dosing optimization of metformin in Chinese patients with type 2 diabetes mellitus. Medicine (Baltimore). 2020;99(46):e23212. doi:10.1097/MD.0000000000023212"
  vignette <- "Li_2020_metformin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "metformin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "metformin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the pharmacokinetic sampling visit.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL/F: (WT/75)^0.688. The reference 75 kg is the cohort median (Li 2020 Table 1, 75 kg, range 51-113 kg). Time-fixed per subject (single sampling visit at steady state).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI serum-creatinine equation (Li 2020 Methods section 2.5), BSA-normalized, in mL/min/1.73 m^2.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL/F: (CRCL/102.5)^0.914. The reference 102.5 mL/min/1.73 m^2 is the cohort median (Li 2020 Table 1, range 46.9-137.7). CKD-EPI eGFR (not a Cockcroft-Gault creatinine clearance); patients had stable renal function, so time-fixed per subject.",
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F in forward selection and not retained (delta OFV < 3.84; Li 2020 Results section 3.4).",
      source_name = "Age"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F in forward selection and not retained (delta OFV < 3.84; Li 2020 Results section 3.4).",
      source_name = "BMI"
    ),
    SNP_SLC22A1_RS622342 = list(
      description = "SLC22A1 (OCT1) rs622342 A>C carrier indicator (1 = heterozygous or homozygous variant).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild type)",
      notes = "Screened on CL/F (variant carriers vs wild type) and not retained (delta OFV < 3.84; Li 2020 Results sections 3.2 and 3.4, Table 3).",
      source_name = "OCT1 rs622342"
    ),
    SNP_SLC22A2_RS316019 = list(
      description = "SLC22A2 (OCT2) rs316019 G>T carrier indicator (1 = heterozygous or homozygous variant).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild type)",
      notes = "Screened on CL/F (variant carriers vs wild type) and not retained (delta OFV < 3.84; Li 2020 Results sections 3.2 and 3.4, Table 3).",
      source_name = "OCT2 rs316019"
    ),
    SNP_SLC47A1_RS2289669 = list(
      description = "SLC47A1 (MATE1) rs2289669 G>A carrier indicator (1 = heterozygous or homozygous variant).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild type)",
      notes = "Screened on CL/F (variant carriers vs wild type) and not retained (delta OFV < 3.84; Li 2020 Results sections 3.2 and 3.4, Table 3).",
      source_name = "MATE1 rs2289669"
    ),
    SNP_SLC47A1_RS2252281 = list(
      description = "SLC47A1 (MATE1) rs2252281 T>C carrier indicator (1 = heterozygous or homozygous variant).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild type)",
      notes = "Screened on CL/F (variant carriers vs wild type) and not retained (delta OFV < 3.84; Li 2020 Results sections 3.2 and 3.4, Table 3).",
      source_name = "MATE1 rs2252281"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 125L,
    n_studies = 1L,
    n_observations = 160L,
    age_range = "27-83 years",
    age_median = "56 years",
    weight_range = "51-113 kg",
    weight_median = "75 kg",
    bmi_median = "26.4 kg/m^2 (range 18.1-35.3)",
    sex_female_pct = 32,
    race_ethnicity = c(Asian = 100),
    disease_state = "Type 2 diabetes mellitus (1999 WHO criteria), hospitalized, on metformin for at least 7 days (pharmacokinetic steady state)",
    renal_function = "CKD-EPI eGFR median 102.5 mL/min/1.73 m^2 (range 46.9-137.7); 5 / 26 / 73 / 10 subjects in the 45-59 / 60-89 / 90-120 / >=120 mL/min/1.73 m^2 strata (Table 2)",
    dose_range = "Metformin hydrochloride immediate-release film-coated tablets (Glucophage), 850-2000 mg/day: 1000 mg b.i.d. (n = 55), 500 mg t.i.d. (n = 29), 500 mg b.i.d. (n = 28), 500 mg q.i.d. (n = 12), 850 mg q.d. (n = 1)",
    regions = "China (Shandong Provincial Qianfoshan Hospital, Jinan)",
    notes = "Prospective open-label sparse-sampling study (1-2 samples per patient, 2-4 h and 10-12 h post-dose), February-September 2017; ChiCTR1800014273. 130 enrolled, 5 excluded for irregular regimens or non-adherence. Plasma metformin by HPLC-UV (LLOQ 0.2 ug/mL; no samples below LLOQ). Demographics from Li 2020 Table 1, Table 2 and Results section 3.1."
  )

  ini({
    lka <- log(1.4)
    label("Absorption rate constant (1/h)") # Li 2020 Table 4, ka = 1.4 1/h (RSE 51.5%)
    lcl <- log(53.0)
    label("Apparent clearance CL/F at WT = 75 kg and eGFR = 102.5 mL/min/1.73 m^2 (L/h)") # Li 2020 Table 4, theta1 = 53.0 L/h (RSE 4.6%)
    lvc <- log(438)
    label("Apparent volume of distribution V/F (L)") # Li 2020 Table 4, V/F = 438 L (RSE 15.0%)
    ltlag <- log(0.914)
    label("Absorption lag time (h)") # Li 2020 Table 4, tlag = 0.914 h (RSE 30.3%)

    e_wt_cl <- 0.688
    label("Power exponent of body weight on CL/F (unitless)") # Li 2020 Table 4, theta2 = 0.688 (RSE 24.6%)
    e_crcl_cl <- 0.914
    label("Power exponent of CKD-EPI eGFR on CL/F (unitless)") # Li 2020 Table 4, theta3 = 0.914 (RSE 19.9%)

    # IIV on CL/F only (Results section 3.3). Table 4 prints 18.0% and the
    # final CL/F equation prints the eta term as EXP(0.1797); 0.1797 is the
    # omega SD, so omega^2 = 0.1797^2 = 0.0323.
    etalcl ~ 0.0323

    # Residual error: exponential model (Results section 3.3), 35.07% in
    # Table 4; encoded as proportional (first-order equivalent of NONMEM
    # Y = F * EXP(EPS) under FOCE-I).
    propSd <- 0.3507
    label("Proportional residual error (fraction)") # Li 2020 Table 4, residual variability 35.07% (RSE 19.3%)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl * (CRCL / 102.5)^e_crcl_cl
    vc <- exp(lvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
