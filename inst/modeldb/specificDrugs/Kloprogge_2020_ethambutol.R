Kloprogge_2020_ethambutol <- function() {
  description <- "Two-compartment population PK model for oral ethambutol at steady state in Malawian adults with drug-sensitive pulmonary tuberculosis on standard fixed-dose-combination therapy (Kloprogge 2020). One transit compartment followed by first-order absorption, with MTT, ka, Q/F and Vp/F fixed to Jonsson 2011; CL/F and Vc/F estimated with allometric body-weight scaling referenced to 50 kg; IIV on CL/F and (fixed) ka and MTT; proportional residual error."
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
      notes = "Allometric scaling of all clearance (CL/F, Q/F; exponent 0.75) and volume (Vc/F, Vp/F; exponent 1) terms referenced to 50 kg (Supplementary Materials section 2: 'centralised around ... 50 kg for ... ethambutol'). Exponents not printed; the standard values of the Jonsson 2011 source model (same 50 kg reference, scaling on all clearance and volume terms) are used. Cohort median 52 kg, range 34-74 kg (Table 1).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "not applicable", verified = TRUE)
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
    dose_range = "Daily fixed-dose-combination tablets (ethambutol 275 mg per tablet, intensive phase) by weight band: 2 tablets 30-37 kg, 3 tablets 37-54 kg, 4 tablets 54-74 kg, 5 tablets > 74 kg (550-1375 mg/day).",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes = "Steady-state sampling on day 14 or 21 of treatment, predose and 2 and 6 h after an observed fasted dose (Methods 'Antibiotic Plasma Concentration Measurement'). Demographics from Table 1 'Pharmacokinetic Data' column."
  )

  ini({
    # Structural parameters (Supplementary S1 Table, ethambutol column) at the
    # 50 kg reference weight; absorption and peripheral-distribution terms
    # fixed to Jonsson 2011 (Supplementary Materials section 2).
    lcl <- log(45.50); label("Apparent clearance CL/F at 50 kg (L/h)")                         # S1 Table 'Cl (l/hr)' ethambutol = 45.50 (CI 43.47-47.6)
    lvc <- log(124.0); label("Apparent central volume Vc/F at 50 kg (L)")                      # S1 Table 'Vc (l)' ethambutol = 124.0 (CI 109.26-142.45)
    lq <- fixed(log(34.3)); label("Apparent inter-compartmental clearance Q/F at 50 kg (L/h)") # S1 Table 'Q (l/hr)' ethambutol = 34.3 fixed
    lvp <- fixed(log(623)); label("Apparent peripheral volume Vp/F at 50 kg (L)")             # S1 Table 'Vp (l)' ethambutol = 623 fixed
    lka <- fixed(log(0.474)); label("First-order absorption rate constant ka (1/h)")           # S1 Table 'ka (hr-1)' ethambutol = 0.474 fixed
    lmtt <- fixed(log(0.789)); label("Mean transit time MTT through the transit compartment (h)") # S1 Table 'MTT (hr)' ethambutol = 0.789 fix

    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)") # Supplementary section 2 'using allometry'; exponent value from Jonsson 2011
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc/F and Vp/F (unitless)")   # Supplementary section 2 'using allometry'; exponent value from Jonsson 2011

    # IIV: S1 Table footnote 'IIV was calculated as 100 x sqrt(e^eta - 1)',
    # so omega^2 = log(1 + CV^2).
    etalcl ~ 0.04312         # S1 Table 'IIV on Cl (%CV)' = 20.99; log(1 + 0.2099^2) = 0.04312
    etalka ~ fixed(0.3760)   # S1 Table 'IIV on ka (%CV)' = 67.56 (held constant); log(1 + 0.6756^2) = 0.3760
    etalmtt ~ fixed(0.7890)  # S1 Table 'IIV on MTT (%CV)' = 109.6 (held constant); log(1 + 1.096^2) = 0.7890

    propSd <- 0.128; label("Proportional residual error (fraction)") # S1 Table 'Proportional residual variability (%)' ethambutol = 12.80
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 50)^e_wt_cl
    vc <- exp(lvc) * (WT / 50)^e_wt_vc
    q <- exp(lq) * (WT / 50)^e_wt_cl
    vp <- exp(lvp) * (WT / 50)^e_wt_vc
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt + etalmtt)
    ktr <- 1 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The transit count is not printed for ethambutol; the one transit
    # compartment ahead of the absorption compartment (ktr = 1 / MTT) of the
    # Jonsson 2011 source model is used, since MTT and ka were fixed to it.
    d/dt(transit1) <- -ktr * transit1
    d/dt(depot) <- ktr * transit1 - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
