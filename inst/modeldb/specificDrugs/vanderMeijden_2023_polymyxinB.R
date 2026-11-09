vanderMeijden_2023_polymyxinB <- function() {
  description <- paste(
    "Preclinical (mouse, neutropenic female CD-1). One-compartment population",
    "PK model with first-order absorption for subcutaneous polymyxin B in",
    "neutropenic mice infected in the thigh or lung (van der Meijden 2023).",
    "Parameters are for a virtual 1 kg mouse, so doses are in mg/kg and CL/F",
    "and V/F are per kg. The PK is non-linear: apparent clearance and the",
    "absorption rate constant both fall with the dose given at each",
    "administration, as power functions of DOSE / 10 mg/kg. Between-subject",
    "variability on ka and CL/F; log-additive residual error. Infection site",
    "(thigh vs lung) was tested and had no effect on PK."
  )
  reference <- paste(
    "van der Meijden A, Aranzana-Climent V, van der Spek H, de Winter BCM,",
    "Couet W, Meletiadis J, Muller AE, van den Berg S. Pharmacokinetic and",
    "pharmacodynamic properties of polymyxin B in Escherichia coli and",
    "Klebsiella pneumoniae murine infection models. J Antimicrob Chemother.",
    "2023;78(3):832-839. doi:10.1093/jac/dkad022"
  )
  vignette <- "vanderMeijden_2023_polymyxinB"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list(
    DOSE = list(
      description = "Polymyxin B dose given at each subcutaneous administration",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE register entry, in mg/kg rather than mg",
        "because the model is parameterised for a virtual 1 kg mouse and",
        "the dosing amount is itself in mg/kg. Enters CL/F and ka as power",
        "functions of DOSE / 10 mg/kg (Table 1 footnote). It is the dose per",
        "ADMINISTRATION, not the total daily dose: the paper notes that",
        "fractionating the same total daily dose changes total exposure",
        "for this reason. Supply it on every record, equal to the",
        "amount of the most recent (or only) dose; place the event columns",
        "before DOSE in the event table."
      ),
      source_name = "DOSE (polymyxin B dose in mg/kg; Table 1 footnote)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "polymyxin B", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "polymyxin B", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, outbred, female; neutropenic)",
    n_subjects = 160L,
    n_studies = 1L,
    age_range = "7-8 weeks",
    weight_range = "20-25 g",
    sex_female_pct = 100,
    disease_state = paste(
      "Neutropenic (cyclophosphamide 150 and 100 mg/kg intraperitoneally 4",
      "and 1 days before infection) mice infected 2 h before dosing with",
      "about 2.5 x 10^6 cfu, either intramuscularly in each thigh or",
      "intranasally (lung)"
    ),
    dose_range = paste(
      "Single subcutaneous dose of 0.5, 2, 8 or 32 mg/kg (thigh model) or",
      "1, 4, 16 or 64 mg/kg (lung model)"
    ),
    regions = "The Netherlands (Erasmus MC, Rotterdam)",
    notes = paste(
      "Plasma sampled at 10 time points per dose level by orbital sinus",
      "bleeding followed by cervical dislocation, two mice per time point,",
      "so each mouse contributed one sample (n_subjects is 8 dose levels x",
      "10 time points x 2 mice, derived from the Methods). Polymyxin B",
      "measured by LC-MS/MS (LLOQ 0.2 mg/L in plasma); samples below the",
      "LLOQ were omitted. Fitted in NONMEM 7.4.2 with FOCE-I on",
      "log-transformed concentrations. Plasma unbound fraction 20% (4.0%",
      "SD) by ultrafiltration, independent of concentration."
    )
  )

  ini({
    # Structural parameters for a virtual 1 kg mouse at the 10 mg/kg
    # reference dose. Bioavailability was not estimable, so CL and V are
    # apparent (CL/F, V/F).
    lka <- log(1.63)
    label("Absorption rate constant ka at a 10 mg/kg dose (1/h)") # Table 1 'ka (h-1)' = 1.63
    lcl <- log(0.384)
    label("Apparent clearance CL/F at a 10 mg/kg dose (L/h/kg)") # Table 1 'CL/F (L/h/kg)' = 0.384
    lvc <- log(0.774)
    label("Apparent central volume of distribution Vc/F (L/kg)") # Table 1 'Vc/F (L/kg)' = 0.774

    # Power effects of the dose per administration (Table 1 footnote:
    # CL = 0.384 x (DOSE/10)^-0.26 and ka = 1.63 x (DOSE/10)^-0.407).
    e_dose_ka <- -0.407
    label("Power exponent of dose per administration (DOSE/10 mg/kg) on ka (unitless)") # Table 1 'Dose on ka' = -0.407
    e_dose_cl <- -0.26
    label("Power exponent of dose per administration (DOSE/10 mg/kg) on CL/F (unitless)") # Table 1 'Dose on CL' = -0.26

    # Between-subject variability, exponential model. Table 1 reports BSV as
    # percent; converted with omega^2 = log(CV^2 + 1).
    etalka ~ 0.48496 # Table 1 'BSV ka (%)' = 79; log(0.79^2 + 1)
    etalcl ~ 0.18440 # Table 1 'BSV CL (%)' = 45; log(0.45^2 + 1)

    # Residual error: additive on log-transformed concentrations, reported
    # as a variance.
    expSd <- 0.30496
    label("Log-additive residual SD (log mg/L)") # Table 1 'Residual error' variance = 0.093; sqrt(0.093)
  })

  model({
    ka <- exp(lka + etalka) * (DOSE / 10)^e_dose_ka
    cl <- exp(lcl + etalcl) * (DOSE / 10)^e_dose_cl
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
