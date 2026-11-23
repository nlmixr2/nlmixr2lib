Konofal_2014_mazindol <- function() {
  description <- paste(
    "One-compartment population PK model for oral mazindol in children aged",
    "9-12.8 years with attention deficit/hyperactivity disorder (ADHD), from",
    "an open-label pilot phase II study (NCT00508677; 21 children, 1 mg once",
    "daily for 7 days). First-order absorption with a fixed 0.3 h lag time",
    "and first-order elimination. Apparent clearance (CL/F) scales with",
    "fat-free mass by an estimated power exponent (reference 28 kg) and",
    "apparent volume (V/F) with age by an estimated power exponent",
    "(reference 10 years). Exponential between-subject variability on CL/F",
    "and V/F; additive (constant SD) residual error on the ng/mL scale."
  )
  reference <- paste(
    "Konofal E, Zhao W, Laouenan C, Lecendreux M, Kaguelidou F, Benadjaoud L,",
    "Mentre F, Jacqz-Aigrain E. Pilot phase II study of mazindol in children",
    "with attention deficit/hyperactivity disorder. Drug Des Devel Ther.",
    "2014;8:2321-2332. doi:10.2147/DDDT.S65495.",
    "Parameter values from Table 3 (final model); covariate equations and",
    "reference values from Results, 'Pharmacokinetics'.",
    sep = " "
  )
  vignette <- "Konofal_2014_mazindol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass, the size descriptor on apparent clearance",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 28 kg, the median of the study population (Results,",
        "'Pharmacokinetics': 'FFMref and ageref were the median values of the",
        "present population, which were 28 kg and 10 years'). Konofal 2014",
        "Methods compute FFM with the Janmahasatian et al. formula",
        "FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT), HT in metres, with",
        "WHSmax = 42.92 kg/m^2 and WHS50 = 30.93 kg/m^2 for boys and",
        "WHSmax = 37.99 kg/m^2 and WHS50 = 35.98 kg/m^2 for girls. The adult",
        "formula is applied unchanged to these 9-12.8-year-old children (no",
        "paediatric age correction). Cohort FFM range 22.9-47.2 kg",
        "(Discussion; printed there with the unit 'kg/m2', a typographical",
        "error for kg). FFM outperformed total, lean and normal-fat body",
        "weight and fixed allometric 3/4 scaling as the CL/F size descriptor.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    AGE = list(
      description = "Age, the covariate on apparent volume of distribution",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 10 years, the median of the study population",
        "(Results, 'Pharmacokinetics'). Range 9.0-12.8 years (Table 1). The",
        "authors caution that the age effect was estimated over a narrow age",
        "range and that age and FFM were correlated (Spearman 0.706), so the",
        "steep exponent should not be extrapolated outside 9-13 years.",
        sep = " "
      ),
      source_name = "age"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mazindol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mazindol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 21,
    n_studies = 1,
    age_range = "9.0-12.8 years",
    age_mean = "10.7 (SD 1.0) years",
    weight_range = "27.5-68.0 kg",
    weight_mean = "34.9 (SD 9.1) kg",
    height_range = "123.1-153.9 cm",
    height_mean = "140.3 (SD 7.6) cm",
    sex_female_pct = 9.5,
    race_ethnicity = "not reported",
    disease_state = paste(
      "ADHD (DSM-IV-TR; combined 76.2%, inattentive 14.3%,",
      "hyperactive-impulsive 9.5%), low responders to methylphenidate for at",
      "least 6 months, ADHD RS-IV total score 43.0 (SD 6.0) at baseline"
    ),
    dose_range = "mazindol 1 mg orally once daily at 7 am for 7 days",
    regions = "France (Robert Debre Hospital, Paris)",
    notes = paste(
      "Demographics in Table 1. 24 children enrolled, 3 excluded before the",
      "first dose; 19 of the 21 dosed children had a full PK profile. Samples",
      "at 0.5, 2, 4, 6 and 12 h after the first dose and 12, 24, 36 and 48 h",
      "after the seventh dose; 183 concentrations, 73 below the 0.5 ng/mL",
      "LLOQ, handled with the M3 method. Dosing times recorded by MEMS.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Table 3, 'Final model' column. CL/F at FFM 28 kg,
    # V/F at age 10 years.
    ltlag <- fixed(log(0.3)); label("Absorption lag time (h)")             # Konofal 2014 Table 3 'Lag-time (h) 0.3' (no RSE); Results: 'the value of 0.3 hour was fixed'
    lka <- log(1.02);         label("First-order absorption rate constant (1/h)") # Konofal 2014 Table 3 'Ka (h-1) 1.02, RSE 6.6%'
    lcl <- log(27.9);         label("Apparent oral clearance CL/F at FFM 28 kg (L/h)") # Konofal 2014 Table 3 'CL/F (L h-1) 27.9, RSE 6.5%'
    lvc <- log(234);          label("Apparent volume of distribution V/F at age 10 years (L)") # Konofal 2014 Table 3 'V/F (L) 234, RSE 4.0%'

    e_ffm_cl <- 1.24; label("Power exponent of fat-free mass on CL/F (unitless)") # Konofal 2014 Table 3 'theta1 1.24, RSE 32.7%' in 'CL/F = CLref * (FFMi/28)^theta1'
    e_age_vc <- 2.28; label("Power exponent of age on V/F (unitless)")            # Konofal 2014 Table 3 'theta2 2.28, RSE 15.4%' in 'V/F = V/Fref * (Agei/10)^theta2'

    # Exponential IIV reported as CV% (Table 3, 'Interindividual variability
    # (%)'); converted with omega^2 = log(1 + CV^2).
    etalcl ~ 0.058293 # Konofal 2014 Table 3 IIV 'CL/F 24.5%' -> log(1 + 0.245^2)
    etalvc ~ 0.012912 # Konofal 2014 Table 3 IIV 'V/F 11.4%' -> log(1 + 0.114^2)

    addSd <- 0.31; label("Additive residual error (ng/mL)") # Konofal 2014 Table 3 'Residual variability (constant SD) ng/mL 0.31, RSE 14.8%'
  })

  model({
    # Covariate equations: Results, 'Pharmacokinetics' and Table 3.
    cl <- exp(lcl + etalcl) * (FFM / 28)^e_ffm_cl
    vc <- exp(lvc + etalvc) * (AGE / 10)^e_age_vc
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg and V/F in L give mg/L; x 1000 converts to ng/mL, the unit
    # of the assay, the LLOQ (0.5 ng/mL) and the additive residual SD.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
