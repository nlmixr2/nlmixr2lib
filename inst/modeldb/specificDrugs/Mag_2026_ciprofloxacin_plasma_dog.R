Mag_2026_ciprofloxacin_plasma_dog <- function() {
  description <- paste0(
    "Preclinical (beagle dog). One-compartment population PK model with first-order absorption, ",
    "an absorption lag time and first-order elimination for CIPROFLOXACIN, the active metabolite, ",
    "in PLASMA after a single fixed 225 mg oral dose of ENROFLOXACIN tablets (target 20 mg/kg) in ",
    "healthy male Beagle dogs (Mag 2026). Exploratory empirical model: the enrofloxacin dose (mg) is ",
    "the input and the formed ciprofloxacin is described as if absorbed from it, so CL/F and V/F ",
    "absorb the unknown fraction of enrofloxacin converted to ciprofloxacin and are not intrinsic ",
    "ciprofloxacin disposition parameters (Mag 2026 Section 3.3). Dose enrofloxacin into depot; Cc is ",
    "the ciprofloxacin plasma concentration. No covariates. The additive residual error magnitude ",
    "is not reported and is fixed to zero. Sibling models: Mag_2026_enrofloxacin_plasma_dog, ",
    "Mag_2026_enrofloxacin_isf_dog and Mag_2026_ciprofloxacin_isf_dog."
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

  # The depot holds the administered enrofloxacin dose; the central state
  # holds formed ciprofloxacin expressed on the enrofloxacin-dose scale.
  compartmentData <- list(
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ciprofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
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
    # Mag 2026 Table 3, rows 'Ciprofloxacin / Plasma'.
    lka <- log(0.46); label("First-order apparent formation/absorption rate constant (1/h)") # Table 3 Ciprofloxacin Plasma: Ka = 0.46 1/h (RSE 59.1%)
    lcl <- log(88.62); label("Apparent clearance CL/F of ciprofloxacin per enrofloxacin dose (L/h)") # Table 3 Ciprofloxacin Plasma: CL/F = 88.62 L/h (RSE 29.7%)
    lvc <- log(804.76); label("Apparent volume of distribution V/F of ciprofloxacin per enrofloxacin dose (L)") # Table 3 Ciprofloxacin Plasma: V/F = 804.76 L (RSE 40.4%)
    ltlag <- log(0.40); label("Lag time (h)") # Table 3 Ciprofloxacin Plasma: Tlag = 0.40 h (RSE 41.0%)

    # Table 4 %CV of an exponential random effect; omega^2 = log(1 + CV^2)
    # (Section 4.9). No covariances are reported.
    etalvc ~ 0.371730 # Table 4 Ciprofloxacin Plasma: V/F IIV = 67.1 %CV (shrinkage 4.2%)
    etalcl ~ 0.488847 # Table 4 Ciprofloxacin Plasma: CL/F IIV = 79.4 %CV (shrinkage 10.7%)

    # Section 4.9: additive residual error for ciprofloxacin in plasma;
    # magnitude unreported, fixed to zero.
    addSd <- fixed(0); label("Additive residual error (ug/mL); not reported") # Section 4.9: additive residual error model; value unreported
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Total ciprofloxacin plasma concentration (ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
