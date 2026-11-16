Mag_2026_ciprofloxacin_isf_dog <- function() {
  description <- paste0(
    "Preclinical (beagle dog). One-compartment population PK model with first-order absorption, ",
    "an absorption lag time and first-order elimination for CIPROFLOXACIN, the active metabolite, ",
    "in subcutaneous INTERSTITIAL FLUID (protein-free in vivo ultrafiltrate) after a single fixed ",
    "225 mg oral dose of ENROFLOXACIN tablets (target 20 mg/kg) in healthy male Beagle dogs ",
    "(Mag 2026). Exploratory empirical model fitted to the ISF data alone: the enrofloxacin dose ",
    "(mg) is the input, and CL/F and V/F absorb both the unknown fraction converted to ciprofloxacin ",
    "and the plasma-to-ISF distribution, so they are descriptive only. Dose enrofloxacin into depot; ",
    "Cc is the ciprofloxacin ISF concentration, not a plasma concentration. No covariates. Residual ",
    "error magnitudes are not reported by the source and are fixed to zero. Sibling models: ",
    "Mag_2026_enrofloxacin_plasma_dog, Mag_2026_enrofloxacin_isf_dog and Mag_2026_ciprofloxacin_plasma_dog."
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
  # holds formed ciprofloxacin in ISF on the enrofloxacin-dose scale.
  # specimen 'tissue' is the library's matrix for interstitial fluid outside
  # the brain; the closed specimen vocabulary has 'brain ISF' only.
  compartmentData <- list(
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "adult (ages not reported)",
    weight_range = "9.3-14.7 kg (mean 11.7, SD 1.5) across the 9 enrolled dogs",
    sex_female_pct = 0,
    disease_state = "healthy, purpose-bred, intact male Beagle dogs",
    dose_range = "single oral dose of 225 mg enrofloxacin (1.5 x 150 mg Baytril Flavour tablets) after a 12 h fast; 15.3-21.0 mg/kg (mean 19.0 mg/kg) in the analysed dogs",
    regions = "Hungary (Aurigon Labs, Dunakeszi)",
    notes = "Nine dogs enrolled; one excluded after vomiting and two lost to ultrafiltration-probe failure, leaving 6 dogs with ISF data (Mag 2026 Section 2.1). ISF collected by a subcutaneous interscapular UF-3-12 ultrafiltration probe over intervals ending at 1, 2, 4, 6, 8, 10, 22.5 and 24 h; each concentration was assigned to the interval midpoint minus the probe dead-volume lag tUF = 70 uL / mean flow = 0.84 h (Section 4.8), and the model was fitted on that corrected clock. No probe-recovery correction. Phoenix NLME 8.5.2.4, FOCE-ELS."
  )

  ini({
    # Mag 2026 Table 3, rows 'Ciprofloxacin / ISF'.
    lka <- log(0.24); label("First-order apparent formation/absorption rate constant (1/h)") # Table 3 Ciprofloxacin ISF: Ka = 0.24 1/h (RSE 50.8%)
    lcl <- log(34.41); label("Apparent clearance CL/F of ciprofloxacin per enrofloxacin dose (L/h)") # Table 3 Ciprofloxacin ISF: CL/F = 34.41 L/h (RSE 28.0%)
    lvc <- log(334.12); label("Apparent volume of distribution V/F of ciprofloxacin per enrofloxacin dose (L)") # Table 3 Ciprofloxacin ISF: V/F = 334.12 L (RSE 25.4%)
    ltlag <- log(1.60); label("Lag time (h)") # Table 3 Ciprofloxacin ISF: Tlag = 1.60 h (RSE 8.7%)

    # Table 4 %CV of an exponential random effect; omega^2 = log(1 + CV^2)
    # (Section 4.9). No covariances are reported.
    etalcl ~ 0.418606 # Table 4 Ciprofloxacin ISF: CL/F IIV = 72.1 %CV (shrinkage 5.3%)
    etalka ~ 0.732138 # Table 4 Ciprofloxacin ISF: Ka IIV = 103.9 %CV (shrinkage 4.3%)

    # Section 4.9: combined additive and multiplicative residual error;
    # magnitudes unreported, fixed to zero.
    addSd <- fixed(0); label("Additive residual error (ug/mL); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Ciprofloxacin concentration in interstitial-fluid ultrafiltrate (ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
