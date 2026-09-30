Kloprogge_2020_isoniazid <- function() {
  description <- "Two-compartment population PK model for oral isoniazid at steady state in Malawian adults with drug-sensitive pulmonary tuberculosis on standard fixed-dose-combination therapy (Kloprogge 2020). First-order absorption; CL/F and Vc/F estimated, with ka, Q/F and Vp/F fixed to healthy-volunteer literature values (Seng 2015). Allometric body-weight scaling referenced to 63 kg; IIV on CL/F and (fixed) ka; proportional residual error. NAT2 acetylator status was not accounted for."
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
      notes = "Allometric scaling of all clearance (CL/F, Q/F; exponent 0.75) and volume (Vc/F, Vp/F; exponent 1) terms referenced to 63 kg (Supplementary Materials section 2: 'Bodyweight was incorporated as a covariate on clearance and volume estimates using allometry, centralised around ... 63 kg ... for isoniazid'). The exponents are not printed; the standard 0.75 / 1 values of the Seng 2015 source model, whose 63 kg reference weight this model adopts, are used. Cohort median 52 kg, range 34-74 kg (Table 1).",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    NAT2 = list(
      description = "NAT2 acetylator phenotype / genotype.",
      units = "(categorical)",
      type = "categorical",
      notes = "Not accounted for (Supplementary Materials section 2: 'isoniazid elimination clearance parameters were not accounted for NAT2 acetylator status'). The large CL/F IIV (46.7% CV) absorbs the acetylator polymorphism."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "not applicable", verified = TRUE)
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
    dose_range = "Daily fixed-dose-combination tablets (isoniazid 75 mg per tablet) by weight band: 2 tablets 30-37 kg, 3 tablets 37-54 kg, 4 tablets 54-74 kg, 5 tablets > 74 kg (150-375 mg/day).",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes = "Steady-state sampling on day 14 or 21 of treatment, predose and 2 and 6 h after an observed fasted dose (Methods 'Antibiotic Plasma Concentration Measurement'). Demographics from Table 1 'Pharmacokinetic Data' column. The same cohort supplied the bacillary-load PKPD model Kloprogge_2020_hrze_bacterialload."
  )

  ini({
    # Structural parameters (Supplementary S1 Table, isoniazid column) at the
    # 63 kg reference weight. ka, Q and Vp were fixed to the healthy-volunteer
    # Seng 2015 values (Supplementary Materials section 2).
    lcl <- log(13.70); label("Apparent clearance CL/F at 63 kg (L/h)")                        # S1 Table 'Cl (l/hr)' isoniazid = 13.70 (CI 12.8-14.88)
    lvc <- log(39.7); label("Apparent central volume Vc/F at 63 kg (L)")                      # S1 Table 'Vc (l)' isoniazid = 39.7 (CI 34.68-45.43)
    lq <- fixed(log(2.9)); label("Apparent inter-compartmental clearance Q/F at 63 kg (L/h)") # S1 Table 'Q (l/hr)' isoniazid = 2.9 fixed
    lvp <- fixed(log(16.5)); label("Apparent peripheral volume Vp/F at 63 kg (L)")            # S1 Table 'Vp (l)' isoniazid = 16.5 fixed
    lka <- fixed(log(0.6)); label("First-order absorption rate constant ka (1/h)")            # S1 Table 'ka (hr-1)' isoniazid = 0.6 fixed

    # Allometric exponents: not printed; the standard values of the Seng 2015
    # source model (which also uses the 63 kg reference weight).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)") # Supplementary section 2 'using allometry'; exponent value from Seng 2015
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc/F and Vp/F (unitless)")   # Supplementary section 2 'using allometry'; exponent value from Seng 2015

    # IIV: S1 Table footnote 'IIV was calculated as 100 x sqrt(e^eta - 1)',
    # so omega^2 = log(1 + CV^2).
    etalcl ~ 0.1970         # S1 Table 'IIV on Cl (%CV)' = 46.66; log(1 + 0.4666^2) = 0.1970
    etalka ~ fixed(0.1260)  # S1 Table 'IIV on ka (%CV)' = 36.64 (held constant); log(1 + 0.3664^2) = 0.1260

    propSd <- 0.211; label("Proportional residual error (fraction)") # S1 Table 'Proportional residual variability (%)' isoniazid = 21.10
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 63)^e_wt_cl
    vc <- exp(lvc) * (WT / 63)^e_wt_vc
    q <- exp(lq) * (WT / 63)^e_wt_cl
    vp <- exp(lvp) * (WT / 63)^e_wt_vc
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L. The source dataset was in umol/L
    # (Figure 1 y-axis); CL and V are unit-independent.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
