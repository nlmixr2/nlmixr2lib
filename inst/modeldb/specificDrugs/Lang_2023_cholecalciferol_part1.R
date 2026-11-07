Lang_2023_cholecalciferol_part1 <- function() {
  description <- paste(
    "Interim two-compartment population PK model of serum",
    "25-hydroxyvitamin D (25(OH)D) after oral vitamin D3 (cholecalciferol)",
    "in children aged 6-17 years with asthma and overweight or obesity,",
    "fitted to VDORA Part 1 only (Lang 2023 Table 2) and used to choose the",
    "Part 2 dosing regimen. First-order absorption of the vitamin D3 dose",
    "into the 25(OH)D central compartment, linear elimination, and a",
    "constant endogenous input that holds the pre-dose concentration at an",
    "estimated baseline (BASE). No covariates. Doses are in nmol of",
    "vitamin D3 (1 IU = 0.025 ug = 0.0650 nmol) and concentrations in",
    "nmol/L of 25(OH)D (1 ng/mL = 2.496 nmol/L). The final model fitted to",
    "Parts 1 and 2 is Lang_2023_cholecalciferol."
  )
  reference <- paste(
    "Lang JE, Gonzalez Ramirez R, Balevic S, Bickel S, Hornik CP, Majure JM,",
    "Venkatachalam S, Snowden J, O'Sullivan B, James L. Pharmacokinetics",
    "of Oral Vitamin D in Children with Obesity and Asthma. Clin",
    "Pharmacokinet. 2023;62(11):1567. doi:10.1007/s40262-023-01285-9."
  )
  vignette <- "Lang_2023_cholecalciferol"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  # Covariates screened in the Part 1 model (ESM Table 1) but not retained:
  # 'The final popPK model did not identify any significant covariates.'
  # Season (SUMMER on CL/F, delta OFV -6.1; FALL on BASE, -4.8), non-White
  # race, normal fat mass and BMI percentile were also screened; they are
  # described here rather than listed because no register name fits them.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (WT/70) and WT on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "WT"
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (FFM/56.1) and FFM on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "FFM"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (LBM/53) and LBM on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "LBM"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "BMI/31.1 on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "BMI"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "AGE/12.8 on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "MALE on CL/F and BASE screened (ESM Table 1; SEXF = 1 - MALE); not retained.",
      source_name = "MALE"
    ),
    DIS_OBESE = list(
      description = "Obese indicator (BMI at or above the 95th percentile for age and sex)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BMI below the 95th percentile)",
      notes = "OBESE on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "OBESE"
    ),
    DIS_OBESE_MORBID = list(
      description = "Morbidly obese indicator (BMI at or above the 99th percentile for age and sex)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BMI below the 99th percentile)",
      notes = "Morbidly OBESE on CL/F and BASE screened (ESM Table 1); not retained.",
      source_name = "MorbidlyOBESE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cholecalciferol", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "25-hydroxyvitamin D", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "25-hydroxyvitamin D", units = "nmol", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 44,
    n_studies = 1,
    age_range = "6.5-17.9 years",
    age_median = "11.2-13.4 years across cohorts A-D",
    weight_range = "BMI 20.8-53.1 kg/m^2 (at or above the 85th percentile for age and sex)",
    sex_female_pct = 50,
    race_ethnicity = c(
      White = 36.4,
      Black = 45.5,
      Asian = 4.5,
      Multiple = 4.5,
      Unknown = 9.1,
      Hispanic_Latino = 20.5
    ),
    disease_state = paste(
      "Physician-diagnosed asthma with overweight or obesity and vitamin D",
      "insufficiency (baseline serum 25(OH)D 10 to < 30 ng/mL). Weight",
      "groups: overweight 9, obese 11, morbidly obese 24."
    ),
    dose_range = paste(
      "Oral vitamin D3 for 16 weeks: cohort A 50,000 IU loading + 6000 IU",
      "daily; B 50,000 IU loading + 10,000 IU daily; C 6000 IU daily; D",
      "600 IU daily."
    ),
    regions = "United States (multi-centre; ClinicalTrials.gov NCT03686150)",
    notes = paste(
      "VDORA Part 1 interim analysis: 227 serum 25(OH)D concentrations from",
      "44 children (median 5, range 1-9 per child). Demographics from Table 1",
      "cohorts A-D. Fitted in NONMEM 7.4 with FOCE-I."
    )
  )

  ini({
    # Structural parameters -- Lang 2023 Table 2 (Part 1 population PK model).
    # Every parameter is apparent (relative to the unknown fraction of the
    # vitamin D3 dose that reaches serum as 25(OH)D).
    lka <- log(0.00178); label("Absorption rate constant (1/h)") # Table 2, Ka = 0.00178 1/h (RSE 28%)
    lcl <- log(0.0858); label("Apparent 25(OH)D clearance CL/F (L/h)") # Table 2, CL/FM = 0.0858 L/h (RSE 55%)
    lvc <- log(35.2); label("Apparent central volume Vc/F (L)") # Table 2, Vc/FM = 35.2 L (RSE 38%)
    lq <- log(0.189); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2, Q/FM = 0.189 L/h (RSE 16%)
    lvp <- log(724); label("Apparent peripheral volume Vp/F (L)") # Table 2, Vp/FM = 724 L (RSE 36%)
    lrbase <- log(45.8); label("Baseline (endogenous steady-state) serum 25(OH)D (nmol/L)") # Table 2, BASE = 45.8 nmol/L (RSE 5%)

    # IIV -- Table 2 %CV (exponential model). The printed RSEs match the
    # variance-scale RSE of CV = sqrt(omega^2), so omega^2 = CV^2.
    etalcl ~ 1.13^2 # Table 2, IIV CL/FM = 113 %CV (RSE 62%)
    etalq ~ 0.28^2 # Table 2, IIV Q/FM = 28 %CV (RSE 45%)
    etalrbase ~ 0.27^2 # Table 2, IIV BASE = 27 %CV (RSE 22%)

    # Residual error -- Table 2
    propSd <- 0.13; label("Proportional residual error (fraction)") # Table 2, proportional error = 13% (RSE 20%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Constant endogenous input that holds serum 25(OH)D at BASE before
    # dosing; at steady state the input equals CL/F * BASE.
    kin <- cl * rbase

    # Both disposition compartments start at the endogenous steady state.
    central(0) <- rbase * vc
    peripheral1(0) <- rbase * vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + kin - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
