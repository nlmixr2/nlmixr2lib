Munn_2024_carprofen_moxidectinCorrected <- function() {
  description <- paste(
    "Preclinical (sheep).",
    "Three-compartment population PK model for carprofen after a single",
    "4 mg/kg intravenous bolus in merino wethers, fitted to plasma and",
    "subcutaneous tissue-cage fluid simultaneously. Structurally identical",
    "to Munn_2024_carprofen_raw, but refitted to carprofen concentrations",
    "that were CORRECTED using moxidectin as an in vivo internal standard:",
    "each carprofen result was multiplied by the mean moxidectin",
    "concentration for that tissue and animal divided by the concurrent",
    "moxidectin result. Comparing this fit against the raw fit is the",
    "central question of the source paper; the authors conclude the",
    "correction did not alter the estimated parameters or their",
    "uncertainty.",
    sep = " ")
  reference <- paste(
    "Munn R, Whittem T.",
    "Moxidectin is a candidate for use as an in vivo internal standard in",
    "pharmacokinetic studies, as demonstrated with use in simultaneous",
    "tissue cage and ultrafiltration fluid collection.",
    "Front Vet Sci. 2024;11:1332974.",
    "doi:10.3389/fvets.2024.1332974.",
    "Parameter values from Supplementary Table S1,",
    "'Moxidectin Corrected' panel.",
    "The tissue-cage model structure is inherited from",
    "Munn R, Whittem T, Woodward AP.",
    "The surface area to volume ratio changes the pharmacokinetic and",
    "pharmacodynamic parameters in the subcutaneous tissue cage model:",
    "as illustrated by carprofen in sheep.",
    "Front Vet Sci. 2022;9:905797. doi:10.3389/fvets.2022.905797.",
    sep = " ")
  vignette <- "Munn_2024_carprofen_moxidectin_internal_standard"

  # See Munn_2024_carprofen_raw for why the tissue cage is a
  # paper-mechanistic state rather than the canonical `isf`.
  paper_specific_compartments <- c("cage6", "cage10")

  units <- list(
    time          = "h",
    dosing        = "mg/kg",
    concentration = "ug/mL"
  )

  # Per-kg normalisation, as in the raw fit: Vc is L/kg, CL is L/h/kg and
  # the dose is 4 mg/kg, so central/vc is mg/L == ug/mL.
  compartmentData <- list(
    central     = list(analyte = "carprofen", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "carprofen", units = "mg/kg", specimen = "plasma", verified = TRUE),
    cage6       = list(analyte = "carprofen", units = "ug/mL", specimen = "tissue", verified = TRUE),
    cage10      = list(analyte = "carprofen", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "sheep (merino wether)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "approximately 18 months",
    weight_range   = "42-51.5 kg",
    sex_female_pct = 0,
    disease_state  = "Healthy (veterinary clinical examination plus routine haematology and biochemistry before enrolment)",
    dose_range     = paste(
      "Single 4 mg/kg carprofen intravenous bolus into a cephalic vein at",
      "time zero, with 0.2 mg/kg moxidectin given subcutaneously 14 days",
      "earlier as the in vivo internal standard."
    ),
    regions        = "Australia (University of Melbourne, Werribee, Victoria)",
    notes          = paste(
      "Same eight animals and same samples as the raw fit; only the",
      "carprofen concentrations differ. Correction was applied only where",
      "a valid concurrent moxidectin result existed (Results 3). Plasma",
      "moxidectin was stable across the 72 h window at a mean of",
      "8.55-8.57 ng/mL with a within-animal CV of 0.05-0.15, and",
      "within-cage means were 8.25-8.58 ng/mL with CV 0.03-14%, which is",
      "the paper's evidence that moxidectin holds a pseudo-steady state",
      "and is therefore usable as an internal standard."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Munn 2024 Supplementary Table S1,
    # "Moxidectin Corrected" panel, Fixed Effects. Monolix 2023R1 SAEM.
    #
    # Vc and CL are identical to the raw fit to the reported precision, so
    # the Discussion's half-life identity ln(2)/(CL/Vc) = 0.693/(0.0016 /
    # 0.045) = 19.5 h holds for this fit too.
    # ------------------------------------------------------------------
    lvc  <- log(0.045)   ; label("Central volume of distribution Vc (log L/kg)")            # Suppl Table S1 corrected: Central Volume = 0.045 L/kg (SE 0.0021, RSE 4.80%)
    lcl  <- log(0.0016)  ; label("Clearance from central CL (log L/h/kg)")                  # Suppl Table S1 corrected: Clearance = 0.0016 L/h.kg (SE 0.000087, RSE 5.52%)
    lk12 <- log(0.12)    ; label("Rate constant central -> peripheral1 k12 (log 1/h)")      # Suppl Table S1 corrected: k12 = 0.12 /h (SE 0.048, RSE 40.9%)
    lk21 <- log(0.32)    ; label("Rate constant peripheral1 -> central k21 (log 1/h)")      # Suppl Table S1 corrected: k21 = 0.32 /h (SE 0.1, RSE 31.8%)

    # ------------------------------------------------------------------
    # Tissue-cage rate constants; the 6 cm cage is the reference level
    # (printed as 0* in Supplementary Table S1).
    # ------------------------------------------------------------------
    lk13 <- log(0.053)   ; label("Rate constant central -> tissue cage k13, 6 cm cage (log 1/h)")   # Suppl Table S1 corrected: k13 = 0.053 /h (SE 0.013, RSE 25.3%)
    lk31 <- log(0.21)    ; label("Rate constant out of tissue cage k31, 6 cm cage (log 1/h)")       # Suppl Table S1 corrected: k31 = 0.21 /h (SE 0.076, RSE 35.9%)

    # ------------------------------------------------------------------
    # Cage-length effects, 10 cm vs the 6 cm reference, applied on the LOG
    # scale (k_10cm = k_6cm * exp(beta)). See Munn_2024_carprofen_raw for
    # the full justification: the additive reading drives k13 negative
    # (0.053 - 0.057 = -0.004 /h) and contradicts Munn 2022 Table 5.
    # ------------------------------------------------------------------
    e_cage10_k13 <- -0.057 ; label("Log-scale effect of a 10 cm (vs 6 cm) tissue cage on k13 (unitless)")   # Suppl Table S1 corrected: Covariate for k13 for Cage Size 10 cm = -0.057 (SE 0.03, RSE 53.1%)
    e_cage10_k31 <- -0.17  ; label("Log-scale effect of a 10 cm (vs 6 cm) tissue cage on k31 (unitless)")   # Suppl Table S1 corrected: Covariate for k31 for Cage Size 10 cm = -0.17 (SE 0.045, RSE 27.0%)

    # ------------------------------------------------------------------
    # IIV -- Supplementary Table S1, "Standard Deviation of the Random
    # Effects" block of the Moxidectin Corrected panel. Log-scale SDs,
    # confirmed by the table's own C.V.(%) column via
    # CV = sqrt(exp(omega^2) - 1); nlmixr2 takes variances, so each entry
    # is SD^2. As in the raw fit, Vc and k21 are numerically zero and not
    # identifiable (RSE 2.25e+6% and 8.82e+6%); they are reproduced as
    # reported.
    # ------------------------------------------------------------------
    etalvc  ~ 6.76e-08   # Suppl Table S1 corrected: Volume SD = 0.00026 (CV 0.026%, RSE 2.25e+6%) -> 0.00026^2
    etalk12 ~ 0.1681     # Suppl Table S1 corrected: k12 SD = 0.41 (CV 43.05%, RSE 42.8%) -> 0.41^2
    etalk21 ~ 7.921e-07  # Suppl Table S1 corrected: k21 SD = 0.00089 (CV 0.089%, RSE 8.82e+6%) -> 0.00089^2
    etalcl  ~ 0.0169     # Suppl Table S1 corrected: Clearance SD = 0.13 (CV 13.46%, RSE 31.7%) -> 0.13^2

    # ------------------------------------------------------------------
    # Residual error -- NOT REPORTED, exactly as in the raw fit.
    # Supplementary Table S1 carries no error-model rows for either panel.
    # Fixed to zero rather than borrowing the upstream Munn 2022 values,
    # which come from a different dataset. See the vignette "Assumptions
    # and deviations" section.
    # ------------------------------------------------------------------
    propSd         <- fixed(0) ; label("Proportional residual error, plasma (fraction; not reported by the source)")        # Not reported; see note above
    propSd_Ccage6  <- fixed(0) ; label("Proportional residual error, 6 cm tissue cage (fraction; not reported)")            # Not reported; see note above
    propSd_Ccage10 <- fixed(0) ; label("Proportional residual error, 10 cm tissue cage (fraction; not reported)")           # Not reported; see note above
  })

  model({
    # 1. Individual plasma-disposition parameters (log-normal IIV).
    vc  <- exp(lvc  + etalvc)
    cl  <- exp(lcl  + etalcl)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. Cage rate constants, one pair per implanted cage length.
    k13_cage6  <- exp(lk13)
    k31_cage6  <- exp(lk31)
    k13_cage10 <- exp(lk13 + e_cage10_k13)
    k31_cage10 <- exp(lk31 + e_cage10_k31)

    kel <- cl / vc

    # 3. Plasma disposition -- self-contained two-compartment model.
    d/dt(central)     <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central     - k21 * peripheral1

    # 4. Plasma concentration: mg/kg over L/kg gives mg/L == ug/mL.
    Cc <- central / vc

    # 5. Tissue cages, driven by Cc with no back-transfer into central
    # (Munn 2022: the cage kinetics are "driven by the central compartment
    # concentrations without altering the central compartment
    # concentrations"). Both cages solve simultaneously because each sheep
    # carried both a 6 cm and a 10 cm cage.
    d/dt(cage6)  <- k13_cage6  * Cc - k31_cage6  * cage6
    d/dt(cage10) <- k13_cage10 * Cc - k31_cage10 * cage10

    Ccage6  <- cage6
    Ccage10 <- cage10

    Cc      ~ prop(propSd)
    Ccage6  ~ prop(propSd_Ccage6)
    Ccage10 ~ prop(propSd_Ccage10)
  })
}
