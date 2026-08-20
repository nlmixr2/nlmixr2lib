Larsen_2018_factorviia_monkey <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey).",
    "Two-compartment population PK model for activated recombinant",
    "factor VII (rFVIIa) in cynomolgus monkeys following single IV",
    "bolus administration, from the Larsen 2018 preclinical scaling",
    "programme (n=27; 16 males and 11 females). Endogenous rFVIIa",
    "was detectable in monkeys and is carried as an additive",
    "baseline (0.287 IU/mL) on the observation. Body-weight-",
    "normalised allometric scaling (exponents 1 for V/V2, 0.75 for",
    "CL/Q; scaling principle I) is applied within-species around",
    "the median weight of 2.78 kg. IIV is estimated on V, CL, and",
    "the endogenous baseline. Parameter values from Larsen 2018",
    "Table 2 (Monkey column)."
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
    central     = list(analyte = "factorviia", units = "IU", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "factorviia", units = "IU", specimen = "plasma", verified = FALSE)
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
        "2.78 kg (Larsen 2018 Table 1). Range 2.2-5.2 kg."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "monkey (cynomolgus)",
    n_subjects     = 27L,
    n_studies      = NA_integer_,
    weight_range   = "2.2-5.2 kg",
    weight_median  = "2.78 kg",
    sex_female_pct = 40.7,
    disease_state  = "healthy (non-haemophilic)",
    dose_range     = "single IV bolus, vehicle, 90, 270, or 3000 ug/kg rFVIIa (mass dose; user supplies the corresponding IU dose)",
    regions        = "Denmark (Novo Nordisk in-house)",
    notes          = paste(
      "27 cynomolgus monkeys (16 M, 11 F). Sampling up to 13",
      "h post-dose. Endogenous rFVIIa was detectable in the",
      "monkey and is modelled as an additive baseline of 0.287",
      "IU/mL on the total observed activity (Larsen 2018",
      "Table 2 Base FVIIa). Data sources: Karpf 2011 and Novo",
      "Nordisk historical unpublished in-house studies. See",
      "Larsen 2018 Table 1 (Monkey column) for the study",
      "inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 2.78 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 2, Monkey column, scaling principle I.
    lvc    <- log(229)           ; label("Central volume of distribution (mL)")            # Table 2 Monkey V  = 229 mL [RSE 12%]
    lcl    <- log(148)           ; label("Clearance (mL/h)")                               # Table 2 Monkey CL = 148 mL/h [RSE 8%]
    lvp    <- log(62.1)          ; label("Peripheral volume of distribution (mL)")         # Table 2 Monkey V2 = 62.1 mL [RSE 11%]
    lq     <- log(42.3)          ; label("Inter-compartmental clearance (mL/h)")           # Table 2 Monkey Q  = 42.3 mL/h [RSE 19%]
    lrbase <- log(0.287)         ; label("Endogenous baseline FVIIa activity (IU/mL)")     # Table 2 Monkey Base FVIIa = 0.287 IU/mL [RSE 12%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2
    e_wt_vp <- fixed(1)          ; label("Allometric exponent on Vp (unitless)")           # Larsen 2018 Eq. 2
    e_wt_q  <- fixed(0.75)       ; label("Allometric exponent on Q (unitless)")            # Larsen 2018 Eq. 2

    # IIV: V (54.7% CV), CL (35.1% CV), Base FVIIa (42.8% CV).
    # Variance on log scale = log(CV^2 + 1).
    etalvc    ~ 0.261756         # Table 2 Monkey IIV on V  = 54.7% [RSE 17%, shrinkage 13%]
    etalcl    ~ 0.116183         # Table 2 Monkey IIV on CL = 35.1% [RSE 16%, shrinkage 12%]
    etalrbase ~ 0.168209         # Table 2 Monkey IIV on Base FVIIa = 42.8% [RSE 22%, shrinkage 32%]

    propSd <- 0.0343             ; label("Proportional residual error (fraction)")         # Table 2 Monkey Residual error = 0.0343 [RSE 14%]
  })

  model({
    vc    <- exp(lvc + etalvc) * (WT / 2.78)^e_wt_vc
    cl    <- exp(lcl + etalcl) * (WT / 2.78)^e_wt_cl
    vp    <- exp(lvp)          * (WT / 2.78)^e_wt_vp
    q     <- exp(lq)           * (WT / 2.78)^e_wt_q
    rbase <- exp(lrbase + etalrbase)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    Cc <- central / vc + rbase
    Cc ~ prop(propSd)
  })
}
