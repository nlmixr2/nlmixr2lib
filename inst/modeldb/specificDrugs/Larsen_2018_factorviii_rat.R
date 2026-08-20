Larsen_2018_factorviii_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "One-compartment population PK model for B-domain truncated",
    "recombinant factor VIII (rFVIII) in male and female",
    "Sprague-Dawley rats with concentration-dependent (nonlinear)",
    "clearance, from the Larsen 2018 preclinical scaling",
    "programme (n=60; 30 males, 30 females). Clearance rises",
    "exponentially with total FVIII activity (Larsen 2018 Eq. 1:",
    "CL(C) = CL_0 * exp(beta_cl * C)) to capture saturation of the",
    "von Willebrand factor buffer at supraphysiological doses. Sex",
    "enters as a fractional multiplier on CL (females clear",
    "rFVIII 39% faster than males), and body-weight-normalised",
    "allometric scaling (exponents 1 for V, 0.75 for CL; scaling",
    "principle I) is applied within-species around the median",
    "weight of 0.32 kg. Endogenous FVIII (1.12 IU/mL) is carried as",
    "an additive baseline on the total observed activity.",
    "Parameter values from Larsen 2018 Table 3 (Rat column)."
  )
  reference <- paste(
    "Larsen MS, Juul RV, Groth AV, Simonsson USH, Kristensen AT,",
    "Knudsen T, Agerso H, Kreilgaard M. (2018).",
    "Prediction of human pharmacokinetics of activated recombinant",
    "factor VII and B-domain truncated factor VIII from animal",
    "population pharmacokinetic models of haemophilia.",
    "European Journal of Pharmaceutical Sciences 119, 265-273.",
    "doi:10.1016/j.ejps.2018.01.035.",
    sep = " "
  )
  vignette <- "Larsen_2018_haemophilia_animal_popPK"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "factorviii", units = "IU", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (rat)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the rat median weight of",
        "0.32 kg (Larsen 2018 Table 1). Range 0.22-0.41 kg."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Gates the fractional female-",
        "multiplier on CL (Larsen 2018 Table 3 Rat Cov",
        "CL,Female = 1.39)."
      ),
      source_name        = "SEXF"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 60L,
    n_studies      = NA_integer_,
    weight_range   = "0.22-0.41 kg",
    weight_median  = "0.32 kg",
    sex_female_pct = 50,
    disease_state  = "healthy (non-haemophilic; endogenous FVIII present)",
    dose_range     = "single IV bolus, vehicle, 50, 250, or 1250 IU/kg rFVIII",
    regions        = "Denmark (Novo Nordisk in-house)",
    notes          = paste(
      "60 Sprague-Dawley rats (30 M, 30 F). Sampling up to 24",
      "h post-dose. Endogenous FVIII was measurable in rats",
      "and is carried as a fitted additive baseline of 1.12",
      "IU/mL. Nonlinear elimination was observed at",
      "supraphysiological doses. Data sources: Novo Nordisk",
      "historical unpublished in-house studies. See Larsen",
      "2018 Table 1 (Rat rFVIII column) for the study",
      "inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 0.32 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 3, Rat column, scaling principle I.
    # NB: `lcl` here is CL_0 (clearance at C = 0) in the exponential nonlinear
    # elimination form CL(C) = CL_0 * exp(beta_cl * C) -- Larsen 2018 Eq. 1.
    lvc      <- log(9.52)        ; label("Central volume of distribution (mL)")            # Table 3 Rat V     = 9.52 mL [RSE 14%]
    lcl      <- log(1.89)        ; label("Clearance at C=0 (mL/h)")                        # Table 3 Rat CL    = 1.89 mL/h [RSE 20%]
    lrbase   <- log(1.12)        ; label("Endogenous baseline FVIII activity (IU/mL)")     # Table 3 Rat Base FVIII = 1.12 IU/mL [RSE 8%]
    lbeta_cl <- log(0.162)       ; label("Nonlinear CL slope on FVIII activity (mL/IU)")   # Table 3 Rat beta = 0.162 [RSE 13%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL_0 (unitless)")         # Larsen 2018 Eq. 2 + Sec 2.5 note

    # Sex covariate: CL_female = CL * 1.39 (Larsen 2018 Table 3 Rat footnote).
    e_sexf_cl <- 1.39            ; label("CL multiplier for female vs male rats (fraction)")  # Table 3 Rat Cov CL,Female = 1.39 [RSE 17%]

    # IIV: only Base FVIII was reported (Table 3 Rat column). No IIV on CL was
    # reported despite the female/male difference being captured structurally.
    # 24.2% CV -> variance on log scale = log(0.242^2 + 1) = 0.056913
    etalrbase ~ 0.056913         # Table 3 Rat IIV on Base FVIII = 24.2% [RSE 33%, shrinkage 30%]

    propSd <- 0.305              ; label("Proportional residual error (fraction)")         # Table 3 Rat Residual error = 0.305 [RSE 9%]
  })

  model({
    sex_cl <- 1 + (e_sexf_cl - 1) * SEXF

    vc      <- exp(lvc)                * (WT / 0.32)^e_wt_vc
    cl_0    <- exp(lcl)                * (WT / 0.32)^e_wt_cl * sex_cl
    rbase   <- exp(lrbase + etalrbase)
    beta_cl <- exp(lbeta_cl)

    # Total observed FVIII activity: dosed + endogenous baseline
    Cc <- central / vc + rbase

    # Concentration-dependent CL (Larsen 2018 Eq. 1) evaluated at total activity
    cl <- cl_0 * exp(beta_cl * Cc)

    d/dt(central) <- -cl * central / vc

    Cc ~ prop(propSd)
  })
}
