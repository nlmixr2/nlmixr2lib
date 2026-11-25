Mag_2026_enrofloxacin_isf_dog <- function() {
  description <- paste0(
    "Preclinical (beagle dog). One-compartment population PK model with first-order absorption, ",
    "an absorption lag time and first-order elimination for enrofloxacin in subcutaneous ",
    "INTERSTITIAL FLUID (protein-free in vivo ultrafiltrate) after a single fixed 225 mg oral dose ",
    "of enrofloxacin tablets (target 20 mg/kg) in healthy male Beagle dogs (Mag 2026). Exploratory ",
    "empirical model fitted to the ISF data alone, not linked to plasma: the oral dose is absorbed ",
    "directly into a compartment whose concentration is the ISF concentration, so CL/F and V/F are ",
    "descriptive and the slow Ka (flip-flop) and long lag carry the delayed distribution into ISF. ",
    "The output Cc is the ISF concentration, not a plasma concentration. No covariates. Residual ",
    "error magnitudes are not reported by the source and are fixed to zero. ",
    "The printed Table 3 CL/F and V/F give about half the exposure observed and simulated by the ",
    "authors at 225 mg (see the vignette Errata). ",
    "Sibling models: ",
    "Mag_2026_enrofloxacin_plasma_dog, Mag_2026_ciprofloxacin_plasma_dog and Mag_2026_ciprofloxacin_isf_dog."
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

  # specimen 'tissue' is the library's matrix for interstitial fluid outside
  # the brain; the closed specimen vocabulary has 'brain ISF' only.
  compartmentData <- list(
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE)
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
    # Mag 2026 Table 3, rows 'Enrofloxacin / ISF'.
    lka <- log(0.13); label("First-order absorption rate constant (1/h)") # Table 3 Enrofloxacin ISF: Ka = 0.13 1/h (RSE 16.3%)
    lcl <- log(12.10); label("Apparent clearance CL/F (L/h)") # Table 3 Enrofloxacin ISF: CL/F = 12.10 L/h (RSE 38.8%)
    lvc <- log(21.70); label("Apparent volume of distribution V/F (L)") # Table 3 Enrofloxacin ISF: V/F = 21.70 L (RSE 33.8%)
    ltlag <- log(1.47); label("Absorption lag time (h)") # Table 3 Enrofloxacin ISF: Tlag = 1.47 h (RSE 11.9%)

    # Table 4 %CV of an exponential random effect; omega^2 = log(1 + CV^2)
    # (Section 4.9). No covariances are reported.
    etalvc ~ 0.074960 # Table 4 Enrofloxacin ISF: V/F IIV = 27.9 %CV (shrinkage 27.8%)
    etalcl ~ 0.884125 # Table 4 Enrofloxacin ISF: CL/F IIV = 119.2 %CV (shrinkage 1.3%)

    # Section 4.9: combined additive and multiplicative residual error;
    # magnitudes unreported, fixed to zero.
    addSd <- fixed(0); label("Additive residual error (ug/mL); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported") # Section 4.9: combined additive + multiplicative error; values unreported
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

    # Enrofloxacin concentration in interstitial-fluid ultrafiltrate (ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
