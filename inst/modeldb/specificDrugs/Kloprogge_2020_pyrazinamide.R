Kloprogge_2020_pyrazinamide <- function() {
  description <- "One-compartment population PK model for oral pyrazinamide at steady state in Malawian adults with drug-sensitive pulmonary tuberculosis on standard fixed-dose-combination therapy (Kloprogge 2020). First-order absorption with ka fixed to Alsultan 2017; CL/F and V/F estimated with allometric body-weight scaling referenced to 70 kg; IIV on CL/F, V/F and ka; combined proportional and additive residual error."
  reference <- paste(
    "Kloprogge F, Mwandumba HC, Banda G, Kamdolozi M, Shani D, Corbett EL,",
    "Kontogianni N, Ward S, Khoo SH, Davies GR, Sloan DJ. (2020).",
    "Longitudinal pharmacokinetic-pharmacodynamic biomarkers correlate with",
    "treatment outcome in drug-sensitive pulmonary tuberculosis: a population",
    "pharmacokinetic-pharmacodynamic analysis.",
    "Open Forum Infect Dis 7(7):ofaa218. doi:10.1093/ofid/ofaa218.",
    sep = " "
  )
  vignette <- "Kloprogge_2020_tuberculosis"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1) referenced to 70 kg (Supplementary Materials section 2: 'centralised around 70 kg patient for rifampicin and pyrazinamide'). Exponents not printed; the standard values of the Alsultan 2017 source model (same 70 kg reference) are used. Cohort median 52 kg, range 34-74 kg (Table 1).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrazinamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 154L,
    n_studies = 1L,
    age_range = "17-61 years",
    age_median = "30 years",
    weight_range = "34-74 kg",
    weight_median = "52 kg",
    sex_female_pct = 31,
    disease_state = "Smear-positive, drug-sensitive pulmonary tuberculosis; 58% HIV co-infected (antiretroviral therapy per national protocol).",
    dose_range = "Daily fixed-dose-combination tablets (pyrazinamide 400 mg per tablet, intensive phase) by weight band: 2 tablets 30-37 kg, 3 tablets 37-54 kg, 4 tablets 54-74 kg, 5 tablets > 74 kg (800-2000 mg/day).",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes = "Steady-state sampling on day 14 or 21 of treatment, predose and 2 and 6 h after an observed fasted dose (Methods 'Antibiotic Plasma Concentration Measurement'); pyrazinamide measured by HPLC-UV. Demographics from Table 1 'Pharmacokinetic Data' column."
  )

  ini({
    # Structural parameters (Supplementary S1 Table, pyrazinamide column) at
    # the 70 kg reference weight; ka fixed to Alsultan 2017 (Supplementary
    # Materials section 2).
    lcl <- log(3.86); label("Apparent clearance CL/F at 70 kg (L/h)")               # S1 Table 'Cl (l/hr)' pyrazinamide = 3.86 (CI 3.66-4.04)
    lvc <- log(45.2); label("Apparent volume of distribution V/F at 70 kg (L)")     # S1 Table 'Vc (l)' pyrazinamide = 45.2 (CI 43.47-47.5)
    lka <- fixed(log(3.94)); label("First-order absorption rate constant ka (1/h)") # S1 Table 'ka (hr-1)' pyrazinamide = 3.94 fixed

    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)") # Supplementary section 2 'using allometry'; exponent value from Alsultan 2017
    e_wt_vc <- fixed(1); label("Allometric exponent on V/F (unitless)")     # Supplementary section 2 'using allometry'; exponent value from Alsultan 2017

    # IIV: S1 Table footnote 'IIV was calculated as 100 x sqrt(e^eta - 1)',
    # so omega^2 = log(1 + CV^2). The ka IIV was estimated despite ka fixed.
    etalcl ~ 0.07863  # S1 Table 'IIV on Cl (%CV)' = 28.60; log(1 + 0.2860^2) = 0.07863
    etalvc ~ 0.001559 # S1 Table 'IIV on Vc (%CV)' = 3.95; log(1 + 0.0395^2) = 0.001559
    etalka ~ 2.820    # S1 Table 'IIV on ka (%CV)' = 397.2 (CI 174.38-1005.27); log(1 + 3.972^2) = 2.820

    # Combined residual error. The additive entry (74.5, no unit printed) is
    # read as the NONMEM variance on the umol/L data scale (Figure 1 y-axis),
    # so SD = sqrt(74.5) = 8.631 umol/L = 8.631 * 123.11 / 1000 = 1.063 mg/L.
    # Read as an SD, 74.5 umol/L (9.2 mg/L) puts the simulated 2.5th
    # percentile of the predose concentration below zero, against the Figure 1
    # VPC.
    propSd <- 0.0306; label("Proportional residual error (fraction)") # S1 Table 'Proportional residual variability (%)' pyrazinamide = 3.06
    addSd <- 1.063; label("Additive residual error (mg/L)")           # S1 Table 'Additive residual variability' pyrazinamide = 74.5 (CI 12.65-462.32); variance in (umol/L)^2 -> sqrt(74.5) umol/L x 123.11 g/mol
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
