Vos_2025_nortriptyline <- function() {
  description <- "One-compartment population PK model for oral nortriptyline in adults with major depressive disorder (Vos 2025): first-order absorption with the rate constant fixed from Kvist 2001, and a well-stirred liver in which the intrinsic clearance CLint (with interindividual variability) sets both the systemic clearance Q * CLint / (Q + CLint) and the hepatic availability Q / (Q + CLint), with hepatic plasma flow Q fixed to 36 L/h."
  reference <- paste(
    "Vos CF, Coenen MJH, ter Hark SE, Schellekens AFA, Aarnoutse RE,",
    "Janzing JGE, ter Heine R. Optimizing Nortriptyline Dosing: A Comparison",
    "between Pharmacogenetics-Based, Phenotype-Based, and Standard Dosing.",
    "Clin Pharmacokinet. 2025;64(7):1071-1078. doi:10.1007/s40262-025-01528-x",
    sep = " "
  )
  vignette <- "Vos_2025_nortriptyline"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "nortriptyline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nortriptyline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 85L,
    n_studies = 1L,
    age_range = "18-65 years (inclusion criterion); mean 42.0 (SD 13.8) years in the 48-participant dosing-comparison subset (Table 1)",
    weight_range = "not reported",
    sex_female_pct = 48,
    race_ethnicity = "not reported",
    disease_state = "Major depressive disorder, free of other psychotropic medication (benzodiazepines up to 4 mg lorazepam-equivalent per day allowed) and of interacting co-medication",
    dose_range = "Oral nortriptyline; standard-dosing arm 100-150 mg/day, genotype-guided arm 50-200 mg/day by CYP2D6 phenotype (ESM 4 mean doses 102.4 and 90.9 mg/day)",
    regions = "Netherlands",
    cyp2d6_phenotype = "48-participant subset: 5 PM, 20 IM, 21 NM, 2 UM (Table 1)",
    notes = "Secondary analysis of the PITA randomized controlled trial (Vos 2023, JAMA Netw Open). The PK model was fitted to the 85 nortriptyline-treated participants with at least one measured plasma concentration (Methods; ESM 1); the demographics in Table 1 describe only the 48 of them with both a single-dose (8-24 h post dose) and a steady-state (10-16 h post dose) sample. No absorption-phase samples were collected, which is why ka was fixed. CYP2D6 genotype was not a covariate of the PK model. NONMEM 7.5, FOCE-I."
  )

  ini({
    lka <- fixed(log(0.533)); label("Absorption rate constant ka (1/h)") # ESM 1 Table A row 'Ka (h-1)' = 0.533 (Fixed); text: fixed to the Kvist 2001 estimate because no absorption-phase samples were collected
    lclint <- log(42.3); label("Hepatic intrinsic clearance CLint (L/h)") # ESM 1 Table A row 'Cl int (L/h)' = 42.3 (RSE 5%); equation 2
    lvc <- log(1080); label("Volume of distribution V (L)") # ESM 1 Table A row 'V (L)' = 1080 (RSE 4%); equation 5
    lqh <- fixed(log(36)); label("Hepatic plasma flow Q (L/h)") # ESM 1 text below equation 5: 'Q is the hepatic plasma flow of 36 L/h'

    etalclint ~ 0.17293 # ESM 1 Table A row 'Interindividual variability in CL (%)' = 43.5% (RSE 16%, shrinkage 4.6%), on CLint per equation 2; log(1 + 0.435^2) = 0.17293

    propSd <- 0.156; label("Proportional residual error (fraction)") # ESM 1 Table A row 'Residual proportional error (%)' = 15.6% (RSE 19%, shrinkage 27%)
  })
  model({
    # Individual parameters (ESM 1 equations 1, 2, 5)
    ka <- exp(lka)
    clint <- exp(lclint + etalclint)
    vc <- exp(lvc)
    qh <- exp(lqh)

    # Well-stirred liver. The systemic clearance (equation 3) and the hepatic
    # availability of the oral dose (equation 4) both follow from CLint, so the
    # apparent oral clearance CL / F reduces to CLint.
    cl <- qh * clint / (qh + clint)
    fh <- qh / (qh + clint)

    # Micro-constant (equation 7)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fh

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L as reported.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
