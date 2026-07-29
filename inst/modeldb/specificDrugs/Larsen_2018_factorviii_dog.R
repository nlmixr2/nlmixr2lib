Larsen_2018_factorviii_dog <- function() {
  description <- paste(
    "Preclinical (haemophilia A dog).",
    "Two-compartment population PK model for B-domain truncated",
    "recombinant factor VIII (rFVIII) in haemophilia A dogs (1",
    "male, 2 female) following single IV bolus administration, from",
    "the Larsen 2018 preclinical scaling programme (n=3). Body-",
    "weight-normalised allometric scaling (exponents 1 for V/V2,",
    "0.75 for CL/Q; scaling principle I) is applied within-species",
    "around the median weight of 19.55 kg. No IIV was estimated",
    "(only 3 subjects). Endogenous FVIII is not present (haemophilia",
    "A background). Parameter values from Larsen 2018 Table 3 (Dog",
    "column). Note: Q is estimated with high RSE (83%), reflecting",
    "the small sample size."
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

  covariateData <- list(
    WT = list(
      description        = "Body weight (haemophilia A dog)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the dog median weight of",
        "19.55 kg (Larsen 2018 Table 1). Range 19.55-19.8 kg",
        "(very narrow, n=3)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "dog (haemophilia A)",
    n_subjects     = 3L,
    n_studies      = NA_integer_,
    weight_range   = "19.55-19.8 kg",
    weight_median  = "19.55 kg",
    sex_female_pct = 66.7,
    disease_state  = "haemophilia A (endogenous FVIII absent)",
    dose_range     = "single IV bolus, 100 IU/kg rFVIII",
    regions        = "Denmark / United States (Novo Nordisk in-house and Agerso 2012)",
    notes          = paste(
      "3 haemophilia A dogs (1 M, 2 F). Sampling up to 80 h",
      "post-dose. Endogenous FVIII is absent by disease.",
      "Data sources: Agerso 2012 (rFVIII in dogs). See Larsen",
      "2018 Table 1 (Dog rFVIII column) for the study",
      "inventory. Only 3 subjects; parameters (particularly",
      "Q, RSE 83%) reflect the small sample size and should",
      "be used accordingly."
    )
  )

  ini({
    # Structural PK parameters, reference weight 19.55 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 3, Dog column, scaling principle I.
    lvc <- log(765)              ; label("Central volume of distribution (mL)")            # Table 3 Dog V  = 765 mL [RSE 5%]
    lcl <- log(67.7)             ; label("Clearance (mL/h)")                               # Table 3 Dog CL = 67.7 mL/h [RSE 4%]
    lvp <- log(92.4)             ; label("Peripheral volume of distribution (mL)")         # Table 3 Dog V2 = 92.4 mL [RSE 32%]
    lq  <- log(4.92)             ; label("Inter-compartmental clearance (mL/h)")           # Table 3 Dog Q  = 4.92 mL/h [RSE 83%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2
    e_wt_vp <- fixed(1)          ; label("Allometric exponent on Vp (unitless)")           # Larsen 2018 Eq. 2
    e_wt_q  <- fixed(0.75)       ; label("Allometric exponent on Q (unitless)")            # Larsen 2018 Eq. 2

    # No IIV parameters reported for the dog rFVIII model (Table 3 Dog column
    # shows blank IIV rows -- consistent with n = 3). Encoded as fixed near-zero
    # IIVs so downstream users can flip them on if desired.
    etalvc ~ fixed(0)            # Table 3 Dog IIV on V  not estimated (n=3)
    etalcl ~ fixed(0)            # Table 3 Dog IIV on CL not estimated (n=3)

    propSd <- 0.209              ; label("Proportional residual error (fraction)")         # Table 3 Dog Residual error = 0.209 [RSE 11%]
  })

  model({
    vc <- exp(lvc + etalvc) * (WT / 19.55)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 19.55)^e_wt_cl
    vp <- exp(lvp)          * (WT / 19.55)^e_wt_vp
    q  <- exp(lq)           * (WT / 19.55)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
