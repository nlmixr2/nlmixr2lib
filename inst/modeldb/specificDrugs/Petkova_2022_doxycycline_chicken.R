Petkova_2022_doxycycline_chicken <- function() {
  description <- "Preclinical (chicken). Two-compartment population PK model with first-order absorption, absorption lag and first-order elimination, parameterised in micro-constants, for doxycycline hyclate given via drinking water at a nominal 20 mg/kg/24 h for 5 days to healthy and Mycoplasma gallisepticum-infected Ross broiler chickens with or without feed N-acetylcysteine; no covariates retained (Petkova 2022). The printed V/F and ke (CL/F = 0.728 L/kg/h) give about half of the paper's own NCA AUC under once-daily dose records; see the vignette before using the absolute exposure scale"
  reference <- "Petkova T, Yordanova A, Milanova A. Population Pharmacokinetics of Doxycycline, Administered Alone or with N-Acetylcysteine, in Chickens with Experimental Mycoplasma gallisepticum Infection. Pharmaceutics. 2022;14(11):2440. doi:10.3390/pharmaceutics14112440"
  vignette <- "Petkova_2022_doxycycline_chicken"
  # Every structural parameter is normalised to body weight exactly as
  # published: the dose is in mg/kg, compartment amounts are mg/kg and V/F is
  # L/kg, so central / vc is mg/L = ug/mL and no unit conversion is applied.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "doxycycline", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "doxycycline", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doxycycline", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    DIS_INFECT_ACTIVE = list(
      description = "Experimental Mycoplasma gallisepticum infection (1 = infected, 0 = healthy)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Health status (healthy vs M. gallisepticum Rlow-challenged at 1 day of age, with clinical airsacculitis at the 40-day PK study) was tested as a covariate and not retained in the final model (Petkova 2022 Section 2.6 and Results: 'The tested co-variables did not improve the model and were excluded').",
      source_name = "health status"
    ),
    CONMED_ACETYLCYSTEINE = list(
      description = "N-acetylcysteine co-administration via feed at 100 mg/kg/24 h (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Drug administration (doxycycline alone or with feed N-acetylcysteine) was tested as a covariate and not retained in the final model (Petkova 2022 Section 2.6 and Results).",
      source_name = "drug administration"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not tested as a covariate. Group means 2.79, 2.84, 2.47 and 2.48 kg (Section 2.2). Because the dose is in mg/kg and V/F is in L/kg, whole-bird V/F and CL/F are proportional to body weight.",
      source_name = "body weight"
    )
  )

  population <- list(
    species = "chicken (Ross broiler)",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "40 days at the start of treatment",
    weight_range = "group means 2.47-2.84 kg (Section 2.2: 2.79 +/- 0.17, 2.84 +/- 0.29, 2.47 +/- 0.27, 2.48 +/- 0.43 kg)",
    sex_female_pct = NA_real_,
    disease_state = "12 healthy and 20 Mycoplasma gallisepticum Rlow-infected (air-sac challenge at 1 day of age) broilers",
    dose_range = "Doxycycline hyclate via drinking water at a nominal 20 mg/kg/24 h for 5 days (water freshly prepared every 12 h); actually received 24.41, 23.70, 22.18 and 23.78 mg/kg/24 h in the four groups. Two groups also received N-acetylcysteine 100 mg/kg/24 h via feed.",
    regions = "Bulgaria (Trakia University, Stara Zagora)",
    sampling = "0, 0.5, 1, 2, 3, 4, 6, 9, 12, 14, 24, 120, 122, 124, 126, 132, 144, 152, 168 and 174 h after the start of treatment; 19 samples per healthy bird and 9-10 per infected bird. Plasma doxycycline by LC-MS/MS, LOQ 0.08 ug/mL; 4.6% of samples were below the LOQ and excluded.",
    notes = "Four groups: healthy doxycycline (n = 6), healthy doxycycline + NAC (n = 6), infected doxycycline (n = 10), infected doxycycline + NAC (n = 10). Sex mix not reported ('from both sexes'). Fit in Phoenix 8.3.4 with FOCE-ELS. Individual concentrations are in Supplementary Tables S2 and S3."
  )

  ini({
    # Table 1 typical values; Table 2 omega variances. All values are as
    # printed. The printed V/F and ke give CL/F = 0.728 L/kg/h, which under the
    # stated 20 mg/kg/24 h regimen reproduces only about half of the paper's own
    # Table S4 NCA AUC0-inf; see the vignette Errata.
    lka <- log(0.093); label("Absorption rate constant ka (1/h)") # Table 1, tvka = 0.093 1/h (printed 95% CI 0.096-0.096)
    lvc <- log(4.73); label("Apparent central volume V/F per kg body weight (L/kg)") # Table 1, tvV/F = 4.73 L/kg
    lkel <- log(0.154); label("Elimination rate constant ke from the central compartment (1/h)") # Table 1, tvke = 0.154 1/h
    lk12 <- log(4.59); label("Distribution rate constant central to peripheral k12 (1/h)") # Table 1, tvk12 = 4.59 1/h
    lk21 <- log(7.06); label("Distribution rate constant peripheral to central k21 (1/h)") # Table 1, tvk21 = 7.06 1/h
    ltlag <- log(0.80); label("Absorption lag time (h)") # Table 1, tvtlag = 0.80 h

    etalka ~ 0.002 # Table 2, eta ka variance 0.002 (printed BSV 0.20 %CV)
    etalvc ~ 0.06 # Table 2, eta V variance 0.06 (printed BSV 6.40 %CV)
    etalkel ~ 1.02e-7 # Table 2, eta ke variance 1.02e-7 (printed BSV 0.0001 %CV)
    etalk12 ~ 0.32 # Table 2, eta k12 variance 0.32 (printed BSV 32.82 %CV)
    etalk21 ~ 0.14 # Table 2, eta k21 variance 0.14 (printed BSV 14.27 %CV)
    etaltlag ~ 0.79 # Table 2, eta tlag variance 0.79 (printed BSV 93.17 %CV)

    propSd <- 0.388; label("Proportional residual error (fraction)") # Table 1, tvCMultStdev = 0.388
    addSd <- 0.121; label("Additive residual error (ug/mL)") # Table 1, stdev0 = 0.121 ug/mL
  })
  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    tlag <- exp(ltlag + etaltlag)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    Cc <- central / vc
    # Phoenix additive + multiplicative error: Var = stdev0^2 + (CMultStdev * C)^2
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
