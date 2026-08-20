Larsen_2018_factorviia_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model for activated recombinant",
    "factor VII (rFVIIa) in male Sprague-Dawley rats following",
    "single IV bolus administration, from the Larsen 2018 preclinical",
    "scaling programme (n=37 rats). Body-weight-normalised allometric",
    "scaling (exponents 1 for V/V2, 0.75 for CL/Q; scaling principle",
    "I) is applied within-species around the median weight of 0.24",
    "kg. IIV is estimated on V and CL. Endogenous rFVIIa was not",
    "detectable in rats. Parameter values from Larsen 2018 Table 2",
    "(Rat column)."
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
      description        = "Body weight (rat)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the rat median weight of",
        "0.24 kg (Larsen 2018 Table 1). Range 0.18-0.34 kg."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 37L,
    n_studies      = NA_integer_,
    weight_range   = "0.18-0.34 kg",
    weight_median  = "0.24 kg",
    sex_female_pct = 0,
    disease_state  = "healthy (non-haemophilic)",
    dose_range     = "single IV bolus, 1000 or 3000 ug/kg rFVIIa (mass dose; user supplies the corresponding IU dose)",
    regions        = "Denmark (Novo Nordisk in-house and published Karpf 2011)",
    notes          = paste(
      "37 male Sprague-Dawley rats. Sampling up to 14 h",
      "post-dose. Endogenous rFVIIa was not detected in rats.",
      "Data sources: Karpf 2011 and Novo Nordisk historical",
      "unpublished in-house studies. See Larsen 2018 Table 1",
      "(Rat column) for the study inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 0.24 kg (Larsen 2018 Table 1).
    # Values from Larsen 2018 Table 2, Rat column, scaling principle I.
    lvc <- log(18.6)             ; label("Central volume of distribution (mL)")            # Table 2 Rat V  = 18.6 mL [RSE 3%]
    lcl <- log(16.5)             ; label("Clearance (mL/h)")                               # Table 2 Rat CL = 16.5 mL/h [RSE 3%]
    lvp <- log(2.68)             ; label("Peripheral volume of distribution (mL)")         # Table 2 Rat V2 = 2.68 mL [RSE 9%]
    lq  <- log(1.63)             ; label("Inter-compartmental clearance (mL/h)")           # Table 2 Rat Q  = 1.63 mL/h [RSE 17%]

    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2
    e_wt_vp <- fixed(1)          ; label("Allometric exponent on Vp (unitless)")           # Larsen 2018 Eq. 2
    e_wt_q  <- fixed(0.75)       ; label("Allometric exponent on Q (unitless)")            # Larsen 2018 Eq. 2

    # IIV: V (16.3% CV) and CL (14.6% CV). Independent (no block correlation reported).
    # Variance on log scale = log(CV^2 + 1).
    etalvc ~ 0.026222            # Table 2 Rat IIV on V  = 16.3% [RSE 17%, shrinkage 15%]
    etalcl ~ 0.021092            # Table 2 Rat IIV on CL = 14.6% [RSE 14%, shrinkage 5%]

    propSd <- 0.161              ; label("Proportional residual error (fraction)")         # Table 2 Rat Residual error = 0.161 [RSE 6%]
  })

  model({
    vc <- exp(lvc + etalvc) * (WT / 0.24)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 0.24)^e_wt_cl
    vp <- exp(lvp)          * (WT / 0.24)^e_wt_vp
    q  <- exp(lq)           * (WT / 0.24)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
