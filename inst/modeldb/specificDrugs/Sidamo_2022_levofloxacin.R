Sidamo_2022_levofloxacin <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption, an absorption lag time and linear elimination for oral levofloxacin in Ethiopian adults with multidrug- or rifampicin-resistant tuberculosis (Sidamo 2022). Serum creatinine is retained as a power covariate on CL/F, but its exponent is not reported and is carried as fixed(0); the combined additive and proportional residual-error magnitudes are not reported and are carried as fixed(0)."
  reference <- "Sidamo T, Rao PS, Aklillu E, Shibeshi W, Park Y, Cho YS, Shin JG, Heysell SK, Mpagama SG, Engidawork E. Population Pharmacokinetics of Levofloxacin and Moxifloxacin, and the Probability of Target Attainment in Ethiopian Patients with Multidrug-Resistant Tuberculosis. Infect Drug Resist. 2022;15:6839-6852. doi:10.2147/IDR.S389442"
  vignette <- "Sidamo_2022_fluoroquinolones_mdr_tb"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on CL/F (Results: 'Cr was used as a covariate for clearance (CL) of LFX'; Methods Equation 1: TVP = theta1 * (COV / COVmean)^theta2). Centred on 1.1 mg/dL, the value printed in the Table 2 'Cr (mg/dL)' row, which matches the levofloxacin cohort's mean creatinine (Table 1: median 1, IQR 0.8-1.3 mg/dL). The exponent theta2 is not reported and is carried as fixed(0); see the vignette for the alternative reading of that row.",
      source_name = "Cr"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods 'Covariates') but not retained for levofloxacin.",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Screened ('gender') but not retained for levofloxacin.",
      source_name = "gender"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained for levofloxacin.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained for levofloxacin.",
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "not reported",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened ('bilirubin') but not retained for levofloxacin. Nutritional status, adverse drug effects and comorbidities were also screened and not retained.",
      source_name = "bilirubin"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    age_range = "adults 18 years and older (median 26 years, IQR 20-30)",
    age_median = "26 years",
    weight_range = "IQR 42.5-52.3 kg",
    weight_median = "45 kg",
    bmi_median = "16.7 kg/m^2 (IQR 15.6-18)",
    sex_female_pct = 42.9,
    race_ethnicity = "Ethiopian",
    disease_state = "Multidrug- or rifampicin-resistant pulmonary tuberculosis (MDR/RR-TB) outpatients, with or without HIV co-infection, on a levofloxacin-based standardized regimen; 42.9% malnourished (BMI < 16 kg/m^2), 38.1% with a co-existing illness, 71.4% with prior TB treatment.",
    dose_range = "Levofloxacin 750 or 1000 mg orally once daily for at least 8 days before PK sampling (steady state).",
    regions = "Southern Ethiopia (Butajira, Yirgalem, Arbaminch and Nigist Eleni Mohammed Memorial hospitals).",
    renal_function = "Serum creatinine median 1 mg/dL (IQR 0.8-1.3).",
    notes = "39 of 55 sampled patients had acceptable profiles (21 levofloxacin, 18 moxifloxacin). Plasma sampled pre-dose and at 2, 4, 6, 9, 12 and 24 h; LC-MS/MS assay. Phoenix NLME, first-order conditional estimation. Demographics from Table 1, levofloxacin column."
  )

  ini({
    # Sidamo 2022 Table 2, LFX rows, 'GM (%RSE)' column. Apparent oral
    # parameters (no IV reference; F not estimated).
    lka <- log(0.2); label("Absorption rate constant ka (1/h)") # Table 2: Ka = 0.2 1/h (4.0 %RSE)
    lvc <- log(122); label("Apparent volume of distribution V/F (L)") # Table 2: V = 122 L (21.5 %RSE)
    lcl <- log(10.4); label("Apparent clearance CL/F at serum creatinine 1.1 mg/dL (L/h)") # Table 2: Cl = 10.4 L/h (7.7 %RSE)
    ltlag <- log(0.2); label("Absorption lag time (h)") # Table 2: Tlag = 0.2 h (6.1 %RSE)

    # Power exponent of serum creatinine on CL/F (Methods Equation 1). Not
    # reported: Table 2 prints one value for the 'Cr (mg/dL)' row, 1.1, read
    # here as the centring value COVmean; see the vignette Assumptions.
    e_creat_cl <- fixed(0); label("Power exponent of (CREAT/1.1) on CL/F (unitless; unreported in source)")

    # IIV: Table 2 '%CV (Shrinkage, %)' column, converted to log-normal
    # variances omega^2 = log(1 + CV^2) (Methods Equation 3, exponential IIV).
    etalka ~ 0.01941 # Table 2: Ka 14 %CV (shrinkage 60%)
    etalvc ~ 0.065413 # Table 2: V 26 %CV (shrinkage 30%)
    etalcl ~ 0.36158 # Table 2: Cl 66 %CV (shrinkage 60%)
    etaltlag ~ 0.0011553 # Table 2: Tlag 3.4 %CV (shrinkage 60%)

    # Combined additive and proportional residual error (Results: 'both
    # additive and multiplicative error models'); neither magnitude is
    # reported and the paper has no supplement.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- magnitude not reported in the source)")
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- magnitude not reported in the source)")
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (CREAT / 1.1)^e_creat_cl
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
