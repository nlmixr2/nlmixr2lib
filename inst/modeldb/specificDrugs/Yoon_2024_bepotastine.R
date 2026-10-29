Yoon_2024_bepotastine <- function() {
  description <- paste(
    "One-compartment population PK model for oral bepotastine in healthy",
    "Korean adults and Korean children aged 2-6 years with allergic rhinitis",
    "or urticaria (Yoon 2024). First-order absorption with an absorption lag",
    "time feeds a one-compartment disposition model. Body weight is the only",
    "covariate, entering apparent clearance and apparent volume as fixed",
    "allometric power functions (exponents 0.75 and 1) normalised to 70 kg.",
    "IIV on ka, CL/F, V/F and lag time with CL-V and CL-lag correlations;",
    "combined proportional and additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Yoon S, Jin BH, Kim CO, Park K, Park MS, Chae D. Pharmacokinetic",
    "Modeling of Bepotastine for Determination of Optimal Dosage Regimen in",
    "Pediatric Patients with Allergic Rhinitis or Urticaria. Pharmaceutics.",
    "2024;16(3):334. doi:10.3390/pharmaceutics16030334. PMCID: PMC10975960.",
    sep = " "
  )
  vignette <- "Yoon_2024_bepotastine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "bepotastine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bepotastine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Yoon 2024 Methods section 2.2: 'Individual CL = Typical value of CL",
        "x (Body weight/70)^0.75' and 'Individual V = Typical value of V x",
        "(Body weight/70)'. Table 2 rows 'Exponent of CL' = 0.75 fixed and",
        "'Exponent of V' = 1.0 fixed. The 70 kg reference is the printed",
        "normalising constant (close to the adult median of 70.8 kg, Table 1).",
        "Weight was the only covariate retained (Results 3.1).",
        sep = " "
      ),
      source_name = "Body weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested on CL and V (Methods 2.2) and not retained. Table 1: adults",
        "median 24 (19-37) years, children median 5 (2-6) years.",
        sep = " "
      ),
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Selected in forward selection but removed in backward elimination",
        "(Results 3.1). Table 1: adults median 0.85 (0.69-1.10) mg/dL,",
        "children median 0.36 (0.28-0.52) mg/dL.",
        sep = " "
      ),
      source_name = "Creatinine"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Tested on CL and V (Methods 2.2) and not retained. Table 1: adults",
        "median 112.0 (80-133), children median 121.5 (90-152)",
        "mL/min/1.73 m^2. The paper does not name the eGFR equation.",
        sep = " "
      ),
      source_name = "eGFR"
    ),
    SEXF = list(
      description = "Female-sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Gender was tested on CL and V (Methods 2.2) and not retained. All 32",
        "adults were male; 11 of 30 children (37%) were female (Table 1).",
        sep = " "
      ),
      source_name = "Gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 62,
    n_studies = 2,
    age_range = "2-37 years (adults 19-37, median 24; children 2-6, median 5)",
    weight_range = "13-89.6 kg (adults 56.8-89.6, median 70.8; children 13-26, median 17.8)",
    sex_female_pct = 17.7,
    race_ethnicity = "Korean (62/62)",
    disease_state = paste(
      "32 healthy adult men plus 30 children aged 2-6 years with allergic",
      "rhinitis or urticaria",
      sep = " "
    ),
    dose_range = paste(
      "Adults: single 10 mg oral tablet. Children: 3 mg oral dry syrup",
      "(1.5 mL of 2 mg/mL) twice daily for at least four days.",
      sep = " "
    ),
    regions = "Republic of Korea (Severance Hospital, Seoul)",
    notes = paste(
      "Table 1 gives the baseline demographics by study. Adults were sampled",
      "intensively (13 samples over 24 h; the 0 h and 24 h samples were BLQ",
      "and excluded); children were sampled sparsely (two samples, one at",
      "0.5-2 h and one at 6-12 h post-dose). Plasma LC-MS/MS assay, LLOQ",
      "1 ng/mL. The paper also reports a PK-Sim PBPK model, which is not",
      "part of this file.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Yoon 2024 Table 2 (final popPK model, FOCE),
    # also quoted in Results 3.1: 'The population estimates of KA, CL, V, and
    # lag time (ALAG) were 4.21 h-1, 28.0 L/h, 103.0 L, and 0.27 h'. All
    # values are apparent (oral) quantities for a 70 kg subject; the paper
    # estimates no bioavailability.
    lka <- log(4.21)
    label("First-order absorption rate constant (1/h)")
    # Table 2 row 'KA (h-1)' = 4.21 (RSE 5%).
    lcl <- log(28.0)
    label("Apparent clearance CL/F at 70 kg (L/h)")
    # Table 2 row 'CL (L/h)' = 28.0 (RSE 4%).
    lvc <- log(103.0)
    label("Apparent volume of distribution V/F at 70 kg (L)")
    # Table 2 row 'V (L)' = 103.0 (RSE 5%).
    ltlag <- log(0.27)
    label("Absorption lag time (h)")
    # Table 2 row 'ALAG (h)' = 0.27 (RSE 5%).

    # Allometric body-weight exponents, held constant in the fit.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT/70) for CL/F (unitless)")
    # Table 2 row 'Exponent of CL' = 0.75; Methods 2.2 equation for CL.
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on (WT/70) for V/F (unitless)")
    # Table 2 row 'Exponent of V' = 1.0; Methods 2.2 equation for V.

    # IIV, log-normal (Methods 2.2). Table 2 reports CV%; converted with
    # omega^2 = log(1 + CV^2). Covariances were estimated for the CL-V and
    # CL-ALAG pairs only (Results 3.1), so the V-lag element is zero.
    # cov = r * sqrt(omega_i^2 * omega_j^2).
    etalka ~ 0.819837
    # Table 2 row 'IIV of KA (CV%)' = 112.7 (RSE 13%); log(1 + 1.127^2).
    etalcl + etalvc + etaltlag ~ c(
      0.0493852,
      0.0314038, 0.0472652,
      0.0329999, 0, 0.0998236
    )
    # Table 2 rows 'IIV of CL (CV%)' = 22.5, 'IIV of V (CV%)' = 22.0, 'IIV of
    # ALAG' = 32.4, 'CORR of CL-V' = 0.65 and 'CORR of CL-ALAG' = 0.47.

    # Residual error, combined (Results 3.1).
    propSd <- 0.103
    label("Proportional residual error (fraction)")
    # Table 2 row 'Proportional error' = 10.3 (RSE 9%), read as a CV% (SD 0.103).
    addSd <- fixed(0.316228)
    label("Additive residual error (ng/mL)")
    # Table 2 row 'Additive error' = 0.1 (not estimated). Results 3.1: 'the
    # additive error variance fixed at 0.1 ng/mL', so SD = sqrt(0.1).
  })

  model({
    # Individual PK parameters with allometric weight scaling (Methods 2.2).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    # One compartment with first-order absorption and lag time (Results 3.1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg, volume in L: mg/L x 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
