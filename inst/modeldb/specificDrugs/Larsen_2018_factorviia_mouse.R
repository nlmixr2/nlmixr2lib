Larsen_2018_factorviia_mouse <- function() {
  description <- paste(
    "Preclinical (mouse).",
    "Two-compartment population PK model for activated recombinant",
    "factor VII (rFVIIa) in male C57BI/6 and NMRI mice following",
    "single IV bolus administration, developed as part of the",
    "Larsen 2018 preclinical scaling programme (n=51 mice pooled",
    "across two strains). Strain enters as a fractional multiplier",
    "on both V and CL (C57BI/6 = 0.588 x V, 0.87 x CL relative to",
    "the NMRI reference), and body-weight-normalised allometric",
    "scaling (exponents 1 for V/V2, 0.75 for CL/Q; scaling",
    "principle I) is applied within-species around the median",
    "weight of 0.026 kg. IIV is estimated only on CL. Endogenous",
    "rFVIIa was not detectable in mice. Parameter values from",
    "Larsen 2018 Table 2 (Mouse column)."
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
      description        = "Body weight (mouse)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal body weight. Used for within-species",
        "allometric scaling around the mouse median weight of",
        "0.026 kg (Larsen 2018 Table 1). Range 0.021-0.036 kg."
      ),
      source_name        = "WT"
    ),
    STRAIN_C57BI6 = list(
      description        = "C57BI/6 mouse-strain indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = NMRI mouse (reference strain)",
      notes              = paste(
        "1 = C57BI/6 mouse; 0 = NMRI mouse. Gates the",
        "fractional strain multipliers on V and CL",
        "(Larsen 2018 Table 2 footnote and Cov V,C57BI/6 /",
        "Cov CL,C57BI/6 rows)."
      ),
      source_name        = "STRAIN_C57BI6"
    )
  )

  population <- list(
    species        = "mouse (C57BI/6 and NMRI, pooled)",
    n_subjects     = 51L,
    n_studies      = NA_integer_,
    weight_range   = "0.021-0.036 kg",
    weight_median  = "0.026 kg",
    sex_female_pct = 0,
    disease_state  = "healthy (non-haemophilic)",
    dose_range     = "single IV bolus, 100 or 1000 ug/kg rFVIIa (mass dose; user supplies the corresponding IU dose)",
    regions        = "Denmark (Novo Nordisk in-house and published Karpf 2011)",
    notes          = paste(
      "51 male mice (C57BI/6 + NMRI, pooled). Sampling up to 8",
      "h post-dose. Endogenous rFVIIa was not detected in mice",
      "and is not carried in the model. Data sources: Karpf",
      "2011 (glycopegylated rFVIIa variants) and Novo Nordisk",
      "historical unpublished in-house studies. See Larsen 2018",
      "Table 1 (Mouse column) for the study inventory."
    )
  )

  ini({
    # Structural PK parameters, reference weight 0.026 kg (Larsen 2018 Table 1),
    # NMRI reference strain (STRAIN_C57BI6 = 0). Values from Larsen 2018
    # Table 2, Mouse column, scaling principle I. RSE% shown in trailing
    # comments.
    lvc <- log(2.72)             ; label("Central volume of distribution (mL)")            # Table 2 Mouse V = 2.72 mL [RSE 4%]
    lcl <- log(3.2)              ; label("Clearance (mL/h)")                               # Table 2 Mouse CL = 3.2 mL/h [RSE 3%]
    lvp <- log(0.384)            ; label("Peripheral volume of distribution (mL)")         # Table 2 Mouse V2 = 0.384 mL [RSE 10%]
    lq  <- log(0.244)            ; label("Inter-compartmental clearance (mL/h)")           # Table 2 Mouse Q = 0.244 mL/h [RSE 21%]

    # Allometric exponents (fixed at theory-based values; Larsen 2018 Sec 2.5 Eq. 2).
    e_wt_vc <- fixed(1)          ; label("Allometric exponent on Vc (unitless)")           # Larsen 2018 Eq. 2
    e_wt_cl <- fixed(0.75)       ; label("Allometric exponent on CL (unitless)")           # Larsen 2018 Eq. 2
    e_wt_vp <- fixed(1)          ; label("Allometric exponent on Vp (unitless)")           # Larsen 2018 Eq. 2
    e_wt_q  <- fixed(0.75)       ; label("Allometric exponent on Q (unitless)")            # Larsen 2018 Eq. 2

    # Strain covariate effects: V_C57BI6 = V * 0.588; CL_C57BI6 = CL * 0.87
    # (Larsen 2018 Table 2 Mouse footnote). Encoded as fractional multipliers
    # applied when STRAIN_C57BI6 == 1.
    e_strain_c57bi6_vc <- 0.588  ; label("V multiplier for C57BI/6 vs NMRI (fraction)")    # Table 2 Mouse Cov V,C57BI/6 = 0.588 [RSE 4%]
    e_strain_c57bi6_cl <- 0.87   ; label("CL multiplier for C57BI/6 vs NMRI (fraction)")   # Table 2 Mouse Cov CL,C57BI/6 = 0.87 [RSE 4%]

    # IIV: only CL was reported (Table 2 Mouse column).
    # 7.4 %CV -> variance on log scale = log(0.074^2 + 1) = 0.005461
    etalcl ~ 0.005461            # Table 2 Mouse IIV on CL = 7.4% [RSE 15%, shrinkage 26%]

    # Residual proportional error (Table 2 Mouse Residual error 0.131 fraction)
    propSd <- 0.131              ; label("Proportional residual error (fraction)")         # Table 2 Mouse Residual error = 0.131 [RSE 9%]
  })

  model({
    # Strain-fractional multipliers (STRAIN_C57BI6 == 0 -> multiplier = 1 = NMRI reference).
    strain_vc <- 1 + (e_strain_c57bi6_vc - 1) * STRAIN_C57BI6
    strain_cl <- 1 + (e_strain_c57bi6_cl - 1) * STRAIN_C57BI6

    # Individual PK parameters (within-species allometric scaling around 0.026 kg).
    vc <- exp(lvc)          * (WT / 0.026)^e_wt_vc * strain_vc
    cl <- exp(lcl + etalcl) * (WT / 0.026)^e_wt_cl * strain_cl
    vp <- exp(lvp)          * (WT / 0.026)^e_wt_vp
    q  <- exp(lq)           * (WT / 0.026)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                    k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
