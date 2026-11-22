Mag_2026_enrofloxacin_plasma_dog <- function() {
  description <- paste0(
    "Preclinical (beagle dog). One-compartment population PK model with first-order absorption, ",
    "an absorption lag time and first-order elimination for total enrofloxacin in PLASMA after a ",
    "single fixed 225 mg oral dose of enrofloxacin tablets (target 20 mg/kg) in healthy male Beagle ",
    "dogs (Mag 2026). Exploratory empirical model fitted separately from the interstitial-fluid data; ",
    "CL and V are apparent (CL/F, V/F). No covariates. Residual error magnitudes are not reported ",
    "by the source and are fixed to zero. ",
    "The printed Table 3 CL/F and V/F give about half the exposure observed and simulated by the ",
    "authors at 225 mg (see the vignette Errata). ",
    "Sibling models from the same paper: ",
    "Mag_2026_enrofloxacin_isf_dog, Mag_2026_ciprofloxacin_plasma_dog and Mag_2026_ciprofloxacin_isf_dog."
  )
  reference <- paste(
    "Mag P, Somogyi Z, Kertenics A, Kovacs M, Szabo P, Jerzsele A. (2026).",
    "Exploratory Pharmacokinetic Characterization of Enrofloxacin and Ciprofloxacin in Plasma",
    "and Interstitial Fluid of Dogs: Nonlinear Mixed-Effects Modeling and PK/PD Target",
    "Attainment Analysis. Antibiotics (Basel) 15(8):716. doi:10.3390/antibiotics15080716.",
    "Parameter estimates: Table 3 (structural) and Table 4 (IIV).",
    sep = " "
  )
  vignette <- "Mag_2026_enrofloxacin_ciprofloxacin_dog"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a median-centred power covariate on CL/F and V/F (Mag 2026 Section 4.9) but not retained: the improvement was not consistent across analytes and matrices. No estimate is reported.",
      source_name = "body weight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "enrofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "adult (ages not reported)",
    weight_range = "9.3-14.7 kg (mean 11.7, SD 1.5) across the 9 enrolled dogs",
    sex_female_pct = 0,
    disease_state = "healthy, purpose-bred, intact male Beagle dogs",
    dose_range = "single oral dose of 225 mg enrofloxacin (1.5 x 150 mg Baytril Flavour tablets) after a 12 h fast; 15.3-21.0 mg/kg (mean 19.0 mg/kg) in the analysed dogs",
    regions = "Hungary (Aurigon Labs, Dunakeszi)",
    notes = "Nine dogs enrolled; one vomited 20 min after dosing and was excluded, leaving 8 dogs with plasma data (Mag 2026 Sections 2.1 and 4.1). Plasma sampled at 0, 10, 20, 40 min and 1, 2, 4, 6, 8, 10, 24 h; LC-MS/MS total plasma concentrations. Phoenix NLME 8.5.2.4, FOCE-ELS."
  )

  ini({
    # Mag 2026 Table 3, rows 'Enrofloxacin / Plasma'.
    lka <- log(0.78); label("First-order absorption rate constant (1/h)") # Table 3 Enrofloxacin Plasma: Ka = 0.78 1/h (RSE 11.2%)
    lcl <- log(27.06); label("Apparent clearance CL/F (L/h)") # Table 3 Enrofloxacin Plasma: CL/F = 27.06 L/h (RSE 32.5%)
    lvc <- log(131.35); label("Apparent volume of distribution V/F (L)") # Table 3 Enrofloxacin Plasma: V/F = 131.35 L (RSE 17.1%)
    ltlag <- log(0.40); label("Absorption lag time (h)") # Table 3 Enrofloxacin Plasma: Tlag = 0.40 h (RSE 31.9%)

    # Mag 2026 Table 4 reports IIV as %CV of an exponential random effect,
    # with %CV = sqrt(exp(omega^2) - 1) * 100 (Section 4.9), so
    # omega^2 = log(1 + CV^2). No covariances are reported.
    etalvc ~ 0.140235 # Table 4 Enrofloxacin Plasma: V/F IIV = 38.8 %CV (shrinkage 6.4%)
    etalcl ~ 0.824800 # Table 4 Enrofloxacin Plasma: CL/F IIV = 113.2 %CV (shrinkage 1.2%)
    etaltlag ~ 0.804932 # Table 4 Enrofloxacin Plasma: Tlag IIV = 111.2 %CV (shrinkage 2.7%)

    # Section 4.9: combined additive and multiplicative residual error. The
    # magnitudes are not reported anywhere in the paper or its supplement,
    # so both are fixed to zero (typical-value and IIV simulation only).
    addSd <- fixed(0); label("Additive residual error (ug/mL); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Total enrofloxacin plasma concentration (mg/L = ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
