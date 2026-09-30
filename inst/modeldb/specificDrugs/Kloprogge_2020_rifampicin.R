Kloprogge_2020_rifampicin <- function() {
  description <- "One-compartment population PK model for oral rifampicin at steady state in Malawian adults with drug-sensitive pulmonary tuberculosis on standard fixed-dose-combination therapy (Kloprogge 2020). Savic transit-compartment absorption (NN = 1.5, MTT and ka fixed to Sloan 2017) followed by first-order absorption; CL/F and V/F estimated with allometric body-weight scaling referenced to 70 kg and a male-sex effect on CL/F; IIV on CL/F, V/F and (fixed) MTT; proportional residual error."
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
      notes = "Allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1) referenced to 70 kg (Supplementary Materials section 2: 'centralised around 70 kg patient for rifampicin'). Exponents not printed; the standard values of the Sloan 2017 source model (same Malawian cohort, same 70 kg reference) are used. Cohort median 52 kg, range 34-74 kg (Table 1).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) -- the typical CL/F = 16.9 L/h applies to females; males take exp(0.183) = 1.20-fold higher CL/F.",
      notes = "S1 Table row 'Cl~Male' = 0.183 (CI -0.036 to 0.385). Supplementary section 2 says 'Rifampicin clearance estimates were centralised around male patients', which the maintainers read as describing the male-sex covariate rather than a male reference category: the parameter is named for males, 0.183 = log(1.20) reproduces the male/female ratio of 1.2 that Sloan 2017 estimated in the same cohort with a female reference, and the female-reference reading reproduces the published median AUC0-24 more closely than the male-reference reading. Entered on the log scale (MU-referenced, Supplementary section 2 'Pi = e^(log(theta) + eta)'); the linear (1 + 0.183) alternative differs by 1.5%. Applied as (1 - SEXF).",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE)
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
    dose_range = "Daily fixed-dose-combination tablets (rifampicin 150 mg per tablet) by weight band: 2 tablets 30-37 kg, 3 tablets 37-54 kg, 4 tablets 54-74 kg, 5 tablets > 74 kg (300-750 mg/day).",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes = "Steady-state sampling on day 14 or 21 of treatment, predose and 2 and 6 h after an observed fasted dose (Methods 'Antibiotic Plasma Concentration Measurement'). Demographics from Table 1 'Pharmacokinetic Data' column. The same cohort supplied the bacillary-load PKPD model Kloprogge_2020_hrze_bacterialload."
  )

  ini({
    # Structural parameters (Supplementary S1 Table, rifampicin column) at the
    # 70 kg female reference. Absorption parameters were fixed to Sloan 2017
    # (Supplementary Materials section 2).
    lcl <- log(16.90); label("Apparent clearance CL/F at 70 kg, female reference (L/h)") # S1 Table 'Cl (l/hr)' rifampicin = 16.90 (CI 14.73-20.16)
    lvc <- log(31.3); label("Apparent volume of distribution V/F at 70 kg (L)")         # S1 Table 'Vc (l)' rifampicin = 31.3 (CI 23.19-39.06)
    lka <- fixed(log(0.277)); label("First-order absorption rate constant ka (1/h)")     # S1 Table 'ka (hr-1)' rifampicin = 0.277 fixed
    lmtt <- fixed(log(0.326)); label("Mean transit time MTT (h)")                        # S1 Table 'MTT (hr)' rifampicin = 0.326 fix
    lnn <- fixed(log(1.5)); label("Number of absorption transit compartments NN (dimensionless)") # S1 Table 'Transit compartments (n)' rifampicin = 1.5 fix
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (CL/F and V/F are F-relative)") # not estimated; apparent parameters

    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F (unitless)") # Supplementary section 2 'using allometry'; exponent value from Sloan 2017
    e_wt_vc <- fixed(1); label("Allometric exponent on V/F (unitless)")     # Supplementary section 2 'using allometry'; exponent value from Sloan 2017

    e_sex_cl <- 0.183; label("Log-scale effect of male sex on CL/F (unitless; applied as (1 - SEXF))") # S1 Table 'Cl~Male' = 0.183 (CI -0.036 to 0.385)

    # IIV: S1 Table footnote 'IIV was calculated as 100 x sqrt(e^eta - 1)',
    # so omega^2 = log(1 + CV^2). The fixed MTT value back-transforms to the
    # Sloan 2017 variance 0.0706, confirming the convention.
    etalcl ~ 0.1490          # S1 Table 'IIV on Cl (%CV)' = 40.08; log(1 + 0.4008^2) = 0.1490
    etalvc ~ 0.5389          # S1 Table 'IIV on Vc (%CV)' = 84.52; log(1 + 0.8452^2) = 0.5389
    etalmtt ~ fixed(0.0706)  # S1 Table 'IIV on MTT (%CV)' = 27.05 (held constant); log(1 + 0.2705^2) = 0.0706

    propSd <- 0.194; label("Proportional residual error (fraction)") # S1 Table 'Proportional residual variability (%)' rifampicin = 19.40
  })
  model({
    cl <- exp(lcl + etalcl + e_sex_cl * (1 - SEXF)) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka)
    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot)

    kel <- cl / vc

    # Savic 2007 transit absorption (analytical gamma input for non-integer
    # NN) into a virtual depot, then first-order absorption, as in the Sloan
    # 2017 source model. f(depot) <- 0 suppresses the bolus; transit() reads
    # the dose from podo(depot).
    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
