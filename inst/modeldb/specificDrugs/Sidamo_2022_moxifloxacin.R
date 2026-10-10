Sidamo_2022_moxifloxacin <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and linear elimination for oral moxifloxacin in Ethiopian adults with multidrug- or rifampicin-resistant tuberculosis (Sidamo 2022). Body mass index is retained as a power covariate on V/F, but its exponent is not reported and is carried as fixed(0); the additive residual-error magnitude is not reported and is carried as fixed(0)."
  reference <- "Sidamo T, Rao PS, Aklillu E, Shibeshi W, Park Y, Cho YS, Shin JG, Heysell SK, Mpagama SG, Engidawork E. Population Pharmacokinetics of Levofloxacin and Moxifloxacin, and the Probability of Target Attainment in Ethiopian Patients with Multidrug-Resistant Tuberculosis. Infect Drug Resist. 2022;15:6839-6852. doi:10.2147/IDR.S389442"
  vignette <- "Sidamo_2022_fluoroquinolones_mdr_tb"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on V/F (Results: 'BMI as the apparent volume distribution (V) covariate'; Methods Equation 1: TVP = theta1 * (COV / COVmean)^theta2). Centred on 17.8 kg/m^2, the value printed in the Table 2 'BMI (Kg.m2)' row, which matches the moxifloxacin cohort's mean BMI (Table 1: median 17.3, IQR 16.2-19.2 kg/m^2). The exponent theta2 is not reported and is carried as fixed(0); see the vignette for the alternative reading of that row.",
      source_name = "BMI"
    )
  )

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods 'Covariates') but not retained for moxifloxacin.",
      source_name = "Cr"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Screened ('gender') but not retained for moxifloxacin.",
      source_name = "gender"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained for moxifloxacin.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained for moxifloxacin.",
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened ('bilirubin') but not retained for moxifloxacin. Nutritional status, adverse drug effects and comorbidities were also screened and not retained.",
      source_name = "bilirubin"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "adults 18 years and older (median 26 years, IQR 20-37)",
    age_median = "26 years",
    weight_range = "IQR 46.6-65.0 kg",
    weight_median = "51 kg",
    bmi_median = "17.3 kg/m^2 (IQR 16.2-19.2)",
    sex_female_pct = 55.6,
    race_ethnicity = "Ethiopian",
    disease_state = "Multidrug- or rifampicin-resistant pulmonary tuberculosis (MDR/RR-TB) outpatients, with or without HIV co-infection, on a moxifloxacin-based standardized regimen without rifampicin; 55.6% malnourished (BMI < 16 kg/m^2), 27.8% with a co-existing illness, 83.3% with prior TB treatment.",
    dose_range = "Moxifloxacin 600 mg orally once daily for at least 8 days before PK sampling (steady state).",
    regions = "Southern Ethiopia (Butajira, Yirgalem, Arbaminch and Nigist Eleni Mohammed Memorial hospitals).",
    renal_function = "Serum creatinine median 1 mg/dL (IQR 0.8-1.4).",
    notes = "39 of 55 sampled patients had acceptable profiles (21 levofloxacin, 18 moxifloxacin). Plasma sampled pre-dose and at 2, 4, 6, 9, 12 and 24 h; LC-MS/MS assay. Phoenix NLME, first-order conditional estimation. Demographics from Table 1, moxifloxacin column."
  )

  ini({
    # Sidamo 2022 Table 2, MXF rows, 'GM (%RSE)' column. Apparent oral
    # parameters (no IV reference; F not estimated). No absorption lag
    # (Results: 'without lag time').
    lka <- log(0.5); label("Absorption rate constant ka (1/h)") # Table 2: Ka = 0.5 1/h (18.8 %RSE)
    lvc <- log(102.2); label("Apparent volume of distribution V/F at BMI 17.8 kg/m^2 (L)") # Table 2: V = 102.2 L (10.0 %RSE)
    lcl <- log(20); label("Apparent clearance CL/F (L/h)") # Table 2: Cl = 20 L/h (3.5 %RSE)

    # Power exponent of BMI on V/F (Methods Equation 1). Not reported:
    # Table 2 prints one value for the 'BMI (Kg.m2)' row, 17.8, read here as
    # the centring value COVmean; see the vignette Assumptions.
    e_bmi_vc <- fixed(0); label("Power exponent of (BMI/17.8) on V/F (unitless; unreported in source)")

    # IIV: Table 2 '%CV (Shrinkage, %)' column, converted to log-normal
    # variances omega^2 = log(1 + CV^2) (Methods Equation 3, exponential IIV).
    etalka ~ 0.2987 # Table 2: Ka 59 %CV (shrinkage 50%)
    etalvc ~ 0.0099503 # Table 2: V 10 %CV (shrinkage 40%)
    etalcl ~ 0.00072873 # Table 2: Cl 2.7 %CV (shrinkage 70%)

    # Additive residual error (Results: 'one-compartment model with additive
    # error models'); its magnitude is not reported and the paper has no
    # supplement.
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- magnitude not reported in the source)")
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (BMI / 17.8)^e_bmi_vc
    cl <- exp(lcl + etalcl)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
