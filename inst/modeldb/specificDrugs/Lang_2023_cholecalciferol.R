Lang_2023_cholecalciferol <- function() {
  description <- paste(
    "Two-compartment population PK model of serum 25-hydroxyvitamin D",
    "(25(OH)D) after oral vitamin D3 (cholecalciferol) in children aged",
    "6-17 years with asthma and overweight or obesity (VDORA Parts 1 and 2",
    "combined, Lang 2023). First-order absorption of the vitamin D3 dose",
    "into the 25(OH)D central compartment, linear elimination, and a",
    "constant endogenous input that holds the pre-dose concentration at an",
    "estimated baseline (BASE). Clearance, central volume and",
    "intercompartmental clearance scale with the BMI-for-age z-score",
    "(fixed power exponents 0.75, 1 and 1; reference z-score 2.36).",
    "Doses are in nmol of vitamin D3 (1 IU = 0.025 ug = 0.0650 nmol) and",
    "concentrations in nmol/L of 25(OH)D (1 ng/mL = 2.496 nmol/L)."
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

  covariateData <- list(
    BMIZ = list(
      description = "BMI-for-age z-score (age- and sex-standardised body mass index)",
      units = "unitless (z-score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent 0.75), Vc/F (exponent 1) and Q/F",
        "(exponent 1), normalised to the cohort median z-score of 2.36",
        "(Lang 2023 Eqs 7-9). Because the effect is a power of BMIZ / 2.36,",
        "the model is only defined for BMIZ > 0; every enrolled child had",
        "BMI at or above the 85th percentile for age and sex (z-score of",
        "about 1.04 or more). The paper does not name the growth reference",
        "used for the z-score; the 85th/95th/99th-percentile weight groups",
        "in Table 1 follow the CDC definitions. Measured at the time of first",
        "dose; covariate values were carried forward after dosing stopped",
        "(Methods)."
      ),
      source_name = "ZBMI"
    )
  )

  # Covariates screened in the combined Part 1 + Part 2 model (ESM Table 3)
  # but not retained. Weight-based allometric scaling (WT/70 and WT/79.5)
  # improved the fit but was superseded by the BMI z-score, which gave
  # better parameter precision (Results).
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL/F, Vc/F, Q/F and Vp/F (reference 70 kg,",
        "delta OFV -13.73; reference 79.5 kg, delta OFV -14.75) and WT on",
        "Ka and BASE screened (ESM Table 3); replaced by BMIZ in the final",
        "model."
      ),
      source_name = "WT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (BMI/30) and BMI on Ka and BASE screened (ESM Table 3); not retained.",
      source_name = "BMI"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (LBM/49.7) and LBM on Ka and BASE screened (ESM Table 3); not retained.",
      source_name = "LBM"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "AGE on CL/F (AGE/12) screened (ESM Table 3); not retained.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "SEX on CL/F, Ka, Q/F, Vc/F, Vp/F and BASE screened (ESM Table 3); not retained.",
      source_name = "SEX"
    ),
    DIS_OBESE = list(
      description = "Obese indicator (BMI at or above the 95th percentile for age and sex)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BMI below the 95th percentile)",
      notes = "OBESE on CL/F, Vc/F, Ka and BASE screened (ESM Table 3); not retained.",
      source_name = "OBESE"
    ),
    DIS_OBESE_MORBID = list(
      description = "Morbidly obese indicator (BMI at or above the 99th percentile for age and sex)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BMI below the 99th percentile)",
      notes = paste(
        "MORBID_OBESE on CL/F, Vc/F, Ka and BASE screened (ESM Table 3);",
        "on Vc/F and Ka it entered the weight-based multivariable step but",
        "was dropped in backward elimination."
      ),
      source_name = "MORBID_OBESE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cholecalciferol", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "25-hydroxyvitamin D", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "25-hydroxyvitamin D", units = "nmol", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 109,
    n_studies = 1,
    age_range = "6-17 years (Table 1 ranges 6.5-17.9 years)",
    age_median = "11.2-13.4 years across the six dosing cohorts",
    weight_range = "BMI at or above the 85th percentile for age and sex (BMI 17.5-53.1 kg/m^2)",
    sex_female_pct = 48.6,
    race_ethnicity = c(
      White = 41.1,
      Black = 41.1,
      AmericanIndian_AlaskaNative = 3.7,
      NativeHawaiian_PacificIslander = 0.9,
      Asian = 1.9,
      Multiple = 4.7,
      Unknown = 6.5,
      Hispanic_Latino = 17.8
    ),
    disease_state = paste(
      "Physician-diagnosed asthma with overweight or obesity (BMI at or",
      "above the 85th percentile) and vitamin D insufficiency (baseline",
      "serum 25(OH)D 10 to < 30 ng/mL). Weight groups (Table 1, n = 107",
      "tabulated): overweight 23, obese 35, morbidly obese 49."
    ),
    dose_range = paste(
      "Oral vitamin D3 for 16 weeks: 50,000 IU loading dose + 6000, 8000",
      "or 10,000 IU daily; 6000 IU daily without loading; or 600 IU daily",
      "(standard of care)."
    ),
    regions = "United States (multi-centre; ClinicalTrials.gov NCT03686150)",
    notes = paste(
      "VDORA study. 848 serum 25(OH)D concentrations from 109 children (46",
      "Part 1, 63 Part 2) measured by immunoassay at baseline, weeks 1, 1-3",
      "(optional), 4, 8, 12, 16 and follow-up weeks 20, 24 and 28.",
      "Demographics from Table 1 (sums over cohorts A-F, n = 107). Fitted",
      "in NONMEM 7.4 with FOCE-I."
    )
  )

  ini({
    # Structural parameters -- Lang 2023 Table 4, column 'VDORA Part 1 and 2'
    # (final combined model, the typical values at the median BMI z-score
    # of 2.36). Every parameter is apparent (relative to the unknown
    # fraction of the vitamin D3 dose that reaches serum as 25(OH)D).
    lka <- log(0.0022); label("Absorption rate constant (1/h)") # Table 4, Ka = 0.0022 1/h (RSE 35%)
    lcl <- log(0.179); label("Apparent 25(OH)D clearance CL/F at BMIZ 2.36 (L/h)") # Table 4, CL/FM = 0.1790 L/h (RSE 8%)
    lvc <- log(30.2); label("Apparent central volume Vc/F at BMIZ 2.36 (L)") # Table 4, Vc/FM = 30.2 L (RSE 55%)
    lq <- log(0.116); label("Apparent intercompartmental clearance Q/F at BMIZ 2.36 (L/h)") # Table 4, Q/FM = 0.116 L/h (RSE 33%)
    lvp <- log(170); label("Apparent peripheral volume Vp/F (L)") # Table 4, Vp/FM = 170 L (RSE 14%)
    lrbase <- log(42.6); label("Baseline (endogenous steady-state) serum 25(OH)D (nmol/L)") # Table 4, BASE = 42.6 nmol/L (RSE 3%)

    # BMI z-score power exponents, fixed (Eqs 7-9; Discussion 'we used fixed
    # exponents to characterize the effect of the BMI z-score')
    e_bmiz_cl <- fixed(0.75); label("Power exponent of BMIZ/2.36 on CL/F (unitless)") # Eq 7, exponent 0.75 fixed
    e_bmiz_vc <- fixed(1); label("Power exponent of BMIZ/2.36 on Vc/F (unitless)") # Eq 8, exponent 1 fixed
    e_bmiz_q <- fixed(1); label("Power exponent of BMIZ/2.36 on Q/F (unitless)") # Eq 9, exponent 1 fixed

    # IIV -- Table 4 %CV (exponential model, Eq 3). The printed RSEs match
    # the variance-scale RSE of CV = sqrt(omega^2), so omega^2 = CV^2.
    etalcl ~ 0.63^2 # Table 4, IIV CL/FM = 63 %CV (RSE 20%)
    etalq ~ 1.20^2 # Table 4, IIV Q/FM = 120 %CV (RSE 39%)
    etalrbase ~ 0.30^2 # Table 4, IIV BASE = 30 %CV (RSE 18%)

    # Residual error -- Table 4
    propSd <- 0.15; label("Proportional residual error (fraction)") # Table 4, proportional error = 15% (RSE 10%)
  })

  model({
    # BMI z-score scaling normalised to the cohort median 2.36 (Eqs 7-9)
    bmiz_ratio <- BMIZ / 2.36

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * bmiz_ratio^e_bmiz_cl
    vc <- exp(lvc) * bmiz_ratio^e_bmiz_vc
    q <- exp(lq + etalq) * bmiz_ratio^e_bmiz_q
    vp <- exp(lvp)
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Constant endogenous input that holds serum 25(OH)D at BASE before
    # dosing (Results, 'constant endogenous input rate to the central
    # compartment'); at steady state the input equals CL/F * BASE.
    kin <- cl * rbase

    # Both disposition compartments start at the endogenous steady state
    # (equal concentrations, so there is no distributional transient).
    central(0) <- rbase * vc
    peripheral1(0) <- rbase * vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + kin - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
