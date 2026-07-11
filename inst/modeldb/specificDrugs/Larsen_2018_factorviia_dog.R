Larsen_2018_factorviia_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog).",
    "Two-compartment population PK model for activated recombinant",
    "factor VII (rFVIIa) in male beagle dogs following single IV",
    "bolus administration, from the Larsen 2018 preclinical scaling",
    "programme (n=10). Body-weight-normalised allometric scaling",
    "(exponents 1 for V/V2, 0.75 for CL/Q; scaling principle I) is",
    "applied within-species around the median weight of 12.8 kg.",
    "IIV is estimated on V and CL. Endogenous rFVIIa was not",
    "detectable in dogs. Parameter values from Larsen 2018 Table 2",
    "(Dog column)."
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
      description        = "Body weight (beagle dog)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the dog median weight of",
        "12.8 kg (Larsen 2018 Table 1). Range 10-21.5 kg."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "beagle dog",
    n_subjects     = 10L,
    n_studies      = NA_integer_,
    weight_range   = "10-21.5 kg",
    weight_median  = "12.8 kg",
    sex_female_pct = 0,
    disease_state  = "healthy (non-haemophilic)",
    dose_range     = "single IV bolus, 270 ug/kg rFVIIa (mass dose; user supplies the corresponding IU dose)",
    regions        = "Denmark (Novo Nordisk in-house)",
    notes          = paste(
      "10 male beagle dogs. Sampling up to 12 h post-dose.",
      "Endogenous rFVIIa was not detected in dogs. Data",
      "sources: Agerso 2011b (rFVIIa in dogs) and Novo Nordisk",
      "historical unpublished in-house studies. See Larsen",
      "2018 Table 1 (Dog column) for the study inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 12.8 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 2, Dog column, scaling principle I.
    lvc <- log(1300)             ; label("Central volume of distribution (mL)")            # Table 2 Dog V  = 1300 mL [RSE 9%]
    lcl <- log(761)              ; label("Clearance (mL/h)")                               # Table 2 Dog CL = 761 mL/h [RSE 11%]
    lvp <- log(200)              ; label("Peripheral volume of distribution (mL)")         # Table 2 Dog V2 = 200 mL [RSE 11%]
    lq  <- log(55.6)             ; label("Inter-compartmental clearance (mL/h)")           # Table 2 Dog Q  = 55.6 mL/h [RSE 31%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2
    e_wt_vp <- fixed(1)          ; label("Allometric exponent on Vp (unitless)")           # Larsen 2018 Eq. 2
    e_wt_q  <- fixed(0.75)       ; label("Allometric exponent on Q (unitless)")            # Larsen 2018 Eq. 2

    # IIV: V (25.8% CV) and CL (34.7% CV). Independent (no block correlation reported).
    etalvc ~ 0.064442            # Table 2 Dog IIV on V  = 25.8% [RSE 23%, shrinkage 0.1%]
    etalcl ~ 0.113694            # Table 2 Dog IIV on CL = 34.7% [RSE 25%, shrinkage 0.1%]

    propSd <- 0.0236             ; label("Proportional residual error (fraction)")         # Table 2 Dog Residual error = 0.0236 [RSE 11%]
  })

  model({
    vc <- exp(lvc + etalvc) * (WT / 12.8)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 12.8)^e_wt_cl
    vp <- exp(lvp)          * (WT / 12.8)^e_wt_vp
    q  <- exp(lq)           * (WT / 12.8)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
