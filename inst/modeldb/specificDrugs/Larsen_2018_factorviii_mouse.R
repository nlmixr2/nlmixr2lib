Larsen_2018_factorviii_mouse <- function() {
  description <- paste(
    "Preclinical (mouse).",
    "One-compartment population PK model for B-domain truncated",
    "recombinant factor VIII (rFVIII) in FVIII-knockout (FVIII KO)",
    "mice following single IV bolus administration, from the",
    "Larsen 2018 preclinical scaling programme (n=132; 71 males,",
    "61 females). Sex enters as a fractional multiplier on CL",
    "(females clear rFVIII 24% faster than males), and body-",
    "weight-normalised allometric scaling (exponents 1 for V, 0.75",
    "for CL; scaling principle I) is applied within-species around",
    "the median weight of 0.021 kg. IIV is estimated on CL.",
    "Endogenous FVIII is not present (KO). Parameter values from",
    "Larsen 2018 Table 3 (Mouse column)."
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
      description        = "Body weight (mouse)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the mouse median weight of",
        "0.021 kg (Larsen 2018 Table 1). Range 0.019-0.021 kg."
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
        "multiplier on CL (Larsen 2018 Table 3 Mouse Cov",
        "CL,Female = 1.24)."
      ),
      source_name        = "SEXF"
    )
  )

  population <- list(
    species        = "mouse (FVIII knockout)",
    n_subjects     = 132L,
    n_studies      = NA_integer_,
    weight_range   = "0.019-0.021 kg",
    weight_median  = "0.021 kg",
    sex_female_pct = 46.2,
    disease_state  = "FVIII knockout (haemophilia A model)",
    dose_range     = "single IV bolus, 17.5, 35, 70, 140, 280, or 1000 IU/kg rFVIII",
    regions        = "Denmark (Novo Nordisk in-house and Elm 2012 / Stennicke 2013)",
    notes          = paste(
      "132 FVIII-knockout mice (71 M, 61 F). Sampling up to",
      "48 h post-dose. Endogenous FVIII is absent by design.",
      "Data sources: Elm 2012, Stennicke 2013 (N8-GP), and",
      "Novo Nordisk historical unpublished in-house studies.",
      "See Larsen 2018 Table 1 (Mouse rFVIII column) for the",
      "study inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 0.021 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 3, Mouse column, scaling principle I.
    lvc <- log(0.716)            ; label("Central volume of distribution (mL)")            # Table 3 Mouse V  = 0.716 mL [RSE 2%]
    lcl <- log(0.0639)           ; label("Clearance (mL/h)")                               # Table 3 Mouse CL = 0.0639 mL/h [RSE 3%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2

    # Sex covariate: CL_female = CL * 1.24 (Larsen 2018 Table 3 Mouse footnote).
    e_sexf_cl <- 1.24            ; label("CL multiplier for female vs male mice (fraction)")  # Table 3 Mouse Cov CL,Female = 1.24 [RSE 4%]

    # IIV: only CL was reported (Table 3 Mouse column).
    # 19% CV -> variance on log scale = log(0.19^2 + 1) = 0.035464
    etalcl ~ 0.035464            # Table 3 Mouse IIV on CL = 19% [RSE 8%, shrinkage 16%]

    propSd <- 0.231              ; label("Proportional residual error (fraction)")         # Table 3 Mouse Residual error = 0.231 [RSE 5%]
  })

  model({
    sex_cl <- 1 + (e_sexf_cl - 1) * SEXF

    vc <- exp(lvc)          * (WT / 0.021)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 0.021)^e_wt_cl * sex_cl

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
