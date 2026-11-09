Yang_2024_ganciclovir <- function() {
  description <- paste(
    "One-compartment population PK model for ganciclovir after IV ganciclovir and",
    "enteral valganciclovir in immunocompromised children (Yang 2024), with",
    "first-order absorption, bioavailability fixed from Acosta 2007, fixed",
    "allometric body-weight scaling to 70 kg on clearance (0.75) and volume (1),",
    "a power effect of time-varying serum creatinine on clearance, and a combined",
    "proportional + additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Yang W, Irwin A, Weerdenburg H, McWhinney B, Cole T, Lei A, Han B, Zhu X,",
    "Gwee A. Serum ganciclovir drug exposure in children receiving standard",
    "ganciclovir dosing. Antimicrob Agents Chemother. 2024;68(10):e00525-24.",
    "doi:10.1128/aac.00525-24.",
    "Bioavailability fixed from Acosta et al. (2007); see",
    "modellib('Acosta_2007_ganciclovir').",
    sep = " "
  )
  vignette <- "Yang_2024_ganciclovir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "valganciclovir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ganciclovir", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling to a 70 kg reference: exponent 0.75 on CL and",
        "1 on V (Yang 2024 Table 3 and its footnote equations). Cohort median",
        "28.2 kg, range 5.6-73.3 kg (Table 1).",
        sep = " "
      ),
      source_name = "weight"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying (Methods: 'SCR was tested a priori as a time-varying",
        "covariate through a power function'). Enters CL as (30 / CREAT)^0.997,",
        "so CL falls almost in proportion to a rise in serum creatinine; 30 umol/L",
        "is the cohort median baseline value (Table 1, range 19-167 umol/L).",
        sep = " "
      ),
      source_name = "SCR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    n_observations = 161L,
    age_median = "9.0 years (range 0.4-17.0)",
    weight_median = "28.2 kg (range 5.6-73.3)",
    sex_female_pct = 47.8,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Immunocompromised children (haematological or solid-organ malignancy,",
      "haematopoietic stem cell or solid-organ transplant, primary",
      "immunodeficiency, or systemic immunosuppression) receiving ganciclovir or",
      "valganciclovir for CMV prophylaxis or treatment; neonates with congenital",
      "CMV excluded. Four children received renal replacement therapy.",
      sep = " "
    ),
    dose_range = paste(
      "IV ganciclovir, standard 5 mg/kg twice daily as a 1-h infusion (median",
      "5, range 0.6-10.9 mg/kg/dose); enteral valganciclovir median 17.8 (range",
      "6.0-47.4) mg/kg/dose.",
      sep = " "
    ),
    regions = paste(
      "Australia (Royal Children's Hospital Melbourne and Queensland Children's",
      "Hospital Brisbane; retrospective audit, April 2018 to January 2023).",
      sep = " "
    ),
    bioassay = "UPLC with QDa mass detection; LLOQ 0.2 mg/L.",
    notes = paste(
      "Demographics from Table 1. 24 children contributed 176 serum samples;",
      "15 below the LLOQ (9%) were excluded, leaving 161 samples from 23 children",
      "(113 after IV ganciclovir, 48 after valganciclovir). Routine TDM sampling:",
      "a trough 0.5-1 h before the dose and a peak 1-4 h after it. Median",
      "baseline serum creatinine 30 umol/L (range 19-167). Renal replacement",
      "therapy was tested on CL and not retained (Discussion).",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2024 Table 3 (final model). Reference subject:
    # WT = 70 kg, serum creatinine = 30 umol/L.
    lcl <- log(14.1); label("Clearance at WT = 70 kg and CREAT = 30 umol/L (L/h)") # Table 3: CL = 14.1 L/h (RSE 9%); footnote CL = 14.10*(30/SCR)^0.997*(weight/70)^0.75
    lvc <- log(62.5); label("Volume of distribution at WT = 70 kg (L)") # Table 3: V = 62.5 L (RSE 4%); footnote V = 62.5*weight/70
    lka <- log(0.66); label("First-order absorption rate constant of enteral valganciclovir (1/h)") # Table 3: Ka = 0.66 1/h (RSE 25%)

    # Bioavailability of ganciclovir per mg of valganciclovir dosed. It absorbs
    # the GCV/VGCV molecular-weight ratio 0.72, so the absolute bioavailability
    # is 0.536 / 0.72 = 0.744 (Results).
    lfdepot <- fixed(log(0.536)); label("Bioavailability of ganciclovir from a valganciclovir mg dose (fraction)") # Table 3: F = 0.536 FIX, taken from Acosta et al. (Results)

    # Covariate effects.
    e_creat_cl <- 0.997; label("Power exponent on (30 umol/L / CREAT) for CL (unitless)") # Table 3: effect of SCR on CL = 0.997 (RSE 7%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL (unitless)") # Table 3: effect of weight on CL = 0.75 FIX
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on V (unitless)") # Table 3: effect of weight on V = 1 FIX

    # Between-subject variability on CL only; IIV on V, F and Ka were not
    # estimable (Supplement, population PK results). Table 3 footnote:
    # CV% = sqrt(omega^2) * 100%, so omega^2 = 0.366^2.
    etalcl ~ 0.133956 # Table 3: BSV CL = 36.6% CV (RSE 15%, shrinkage 10%)

    # Residual unexplained variability: combined proportional + additive.
    propSd <- 0.392; label("Proportional residual error (fraction)") # Table 3: prop.error.sd = 39.2%
    addSd <- 0.144; label("Additive residual error (mg/L)") # Table 3: add.error.sd = 0.144 mg/L
  })

  model({
    cl <- exp(lcl + etalcl) * (30 / CREAT)^e_creat_cl * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka)

    kel <- cl / vc

    # IV ganciclovir doses (mg) go to central; enteral valganciclovir doses
    # (mg of valganciclovir) go to depot.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Dose in mg, volume in L -> serum concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
