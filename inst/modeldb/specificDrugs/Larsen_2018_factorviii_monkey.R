Larsen_2018_factorviii_monkey <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey).",
    "One-compartment population PK model for B-domain truncated",
    "recombinant factor VIII (rFVIII) in male cynomolgus monkeys",
    "with concentration-dependent (nonlinear) clearance, from the",
    "Larsen 2018 preclinical scaling programme (n=35). Clearance",
    "rises exponentially with total FVIII activity (Larsen 2018",
    "Eq. 1: CL(C) = CL_0 * exp(beta_cl * C)) to capture saturation",
    "of the von Willebrand factor buffer at supraphysiological",
    "doses. Body-weight-normalised allometric scaling (exponents 1",
    "for V, 0.75 for CL; scaling principle I) is applied within-",
    "species around the median weight of 2.96 kg. IIV is estimated",
    "on CL and the endogenous baseline. Endogenous FVIII (1.75",
    "IU/mL) is carried as an additive baseline on the total observed",
    "activity. Parameter values from Larsen 2018 Table 3 (Monkey",
    "column)."
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
      description        = "Body weight (cynomolgus monkey)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the monkey median weight of",
        "2.96 kg (Larsen 2018 Table 1). Range 2.11-3.97 kg."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "monkey (cynomolgus)",
    n_subjects     = 35L,
    n_studies      = NA_integer_,
    weight_range   = "2.11-3.97 kg",
    weight_median  = "2.96 kg",
    sex_female_pct = 0,
    disease_state  = "healthy (non-haemophilic; endogenous FVIII present)",
    dose_range     = "single IV bolus, vehicle, 50, 250, 500, 1000, 1250, or 2500 IU/kg rFVIII",
    regions        = "Denmark (Novo Nordisk in-house)",
    notes          = paste(
      "35 male cynomolgus monkeys. Sampling up to 480 h",
      "post-dose (Table 1 lists 480 h). Endogenous FVIII was",
      "measurable in monkeys and is carried as a fitted",
      "additive baseline of 1.75 IU/mL. Nonlinear elimination",
      "was observed at supraphysiological doses. Data sources:",
      "Novo Nordisk historical unpublished in-house studies.",
      "See Larsen 2018 Table 1 (Monkey rFVIII column) for the",
      "study inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 2.96 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 3, Monkey column, scaling principle I.
    # NB: `lcl` here is CL_0 (clearance at C = 0) in the exponential nonlinear
    # elimination form CL(C) = CL_0 * exp(beta_cl * C) -- Larsen 2018 Eq. 1.
    lvc      <- log(145)         ; label("Central volume of distribution (mL)")            # Table 3 Monkey V     = 145 mL [RSE 5%]
    lcl      <- log(13.1)        ; label("Clearance at C=0 (mL/h)")                        # Table 3 Monkey CL    = 13.1 mL/h [RSE 13%]
    lrbase   <- log(1.75)        ; label("Endogenous baseline FVIII activity (IU/mL)")     # Table 3 Monkey Base FVIII = 1.75 IU/mL [RSE 7%]
    lbeta_cl <- log(0.0355)      ; label("Nonlinear CL slope on FVIII activity (mL/IU)")   # Table 3 Monkey beta = 0.0355 [RSE 8%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL_0 (unitless)")         # Larsen 2018 Eq. 2 + Sec 2.5 note

    # IIV: CL (55% CV) and Base FVIII (39.5% CV).
    # Variance on log scale = log(CV^2 + 1).
    etalcl    ~ 0.264285         # Table 3 Monkey IIV on CL = 55%  [RSE 18%, shrinkage 21%]
    etalrbase ~ 0.144987         # Table 3 Monkey IIV on Base FVIII = 39.5% [RSE 13%, shrinkage 6%]

    propSd <- 0.229              ; label("Proportional residual error (fraction)")         # Table 3 Monkey Residual error = 0.229 [RSE 5%]
  })

  model({
    vc      <- exp(lvc)                * (WT / 2.96)^e_wt_vc
    cl_0    <- exp(lcl + etalcl)       * (WT / 2.96)^e_wt_cl
    rbase   <- exp(lrbase + etalrbase)
    beta_cl <- exp(lbeta_cl)

    Cc <- central / vc + rbase

    cl <- cl_0 * exp(beta_cl * Cc)

    d/dt(central) <- -cl * central / vc

    Cc ~ prop(propSd)
  })
}
