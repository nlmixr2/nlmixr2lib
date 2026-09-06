DominguezMore_2024_rutin_rabbit <- function() {
  description <- paste(
    "Preclinical (rabbit). Two-compartment intravenous population PK model for",
    "the flavonoid rutin (quercetin-3-O-rutinoside) in male New Zealand White",
    "rabbits, given either as the pure compound or within a hydroethanolic",
    "extract of Physalis peruviana calyces (Dominguez More 2024). The model is",
    "parameterised in micro-constants (V, k, k12, k21) exactly as published.",
    "The extract matrix acts as a categorical covariate on V, k and k12: it",
    "raises the central volume and the elimination rate constant and lowers the",
    "central-to-peripheral distribution rate constant, which the authors",
    "attribute to displacement of rutin from plasma protein binding by other",
    "extract constituents. V is carried on the LOGIT scale (Eq. 7) and the",
    "remaining parameters on the log scale (Eqs. 8-10), following the",
    "Monolix parameter transformations the authors selected."
  )
  reference <- paste(
    "Dominguez More GP, Rey DP, Valderrama IH, Ospina LF, Aragon DM.",
    "Rutin and Physalis peruviana extract: population pharmacokinetics in",
    "New Zealand rabbits. Pharmaceutics. 2024;16(10):1241.",
    "doi:10.3390/pharmaceutics16101241",
    sep = " "
  )
  vignette <- "DominguezMore_2024_rutin_physalis"

  units <- list(time = "h", dosing = "ug/kg", concentration = "ng/mL")

  compartmentData <- list(
    central     = list(analyte = "rutin", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rutin", units = "ug/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_RUTIN_EXTRACT = list(
      description        = "Source of the administered rutin: within the Physalis peruviana calyx extract matrix versus the isolated pure compound",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pure rutin)",
      notes              = paste(
        "Time-fixed per animal; each rabbit received only one treatment.",
        "Dominguez More 2024 Section 2.3.2 describes it as 'the source of rutin",
        "as categorical covariate (pure compound or extract)'. The intravenous",
        "extract dose (100 mg/kg of extract, containing 14.80 ug rutin per mg)",
        "delivers 1.48 mg/kg of rutin, four times the 0.37 mg/kg pure-rutin",
        "dose, so the covariate is confounded with dose level in this study",
        "design; the authors did not fit a separate dose effect."
      ),
      source_name        = "source of rutin"
    )
  )

  population <- list(
    species        = "rabbit (New Zealand White, male)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "9-10 weeks",
    weight_range   = "1.8-2.2 kg",
    sex_female_pct = 0,
    race_ethnicity = "not applicable",
    disease_state  = "healthy",
    dose_range     = paste(
      "Single intravenous marginal-ear-vein dose (0.5 mL/kg): pure rutin",
      "0.37 mg/kg, or P. peruviana calyx extract 100 mg/kg equivalent to",
      "1.48 mg/kg of rutin"
    ),
    regions        = "Colombia (Universidad Nacional de Colombia, Bogota)",
    notes          = paste(
      "Twenty rabbits in total were split across four arms (n = 5 each):",
      "intravenous pure rutin, intravenous extract, oral pure rutin and oral",
      "extract. Only the two intravenous arms (n = 10) contributed to this",
      "model; the two oral arms fed the companion metabolite model",
      "DominguezMore_2024_quercetin_rabbit. Plasma was sampled at 0, 0.083,",
      "0.166, 0.333, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12 and 24 h and assayed",
      "by a validated UHPLC-UV method (260 nm) over 100-10,000 ng/mL.",
      "Free quercetin was not detectable in any intravenous sample. Fitted in",
      "Monolix 2024R1 (MonolixSuite); AIC 1159; all RSEs below 20%; validated",
      "by VPC (1000 simulations) and by a 1000-run bootstrap (177 runs did not",
      "converge). See Dominguez More 2024 Sections 2.2.2-2.2.3 and Table 3."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural disposition -- Table 3, "Model Estimations / Population
    # Value" column. The model is parameterised in micro-constants, as
    # stated in Section 2.3.2 ("Models were parametrized in microconstants,
    # i.e., elimination rate constant (k), distribution rate constant from
    # compartment 1 to compartment 2 (k12), and distribution rate constant
    # from compartment 2 to compartment 1 (k21)"). Clearance is a derived
    # quantity, Cl = k * V (Eq. 5).
    # ----------------------------------------------------------------------
    lvc  <- log(0.048) ; label("Central volume of distribution V (L/kg)")                                # Table 3: V = 0.048 L/kg (RSE 5.2%)
    lkel <- log(1.924) ; label("Elimination rate constant k (1/h)")                                      # Table 3: k = 1.924 1/h (RSE 9.1%)
    lk12 <- log(3.666) ; label("Central-to-peripheral distribution rate constant k12 (1/h)")             # Table 3: k12 = 3.666 1/h (RSE 19.3%)
    lk21 <- log(3.777) ; label("Peripheral-to-central distribution rate constant k21 (1/h)")             # Table 3: k21 = 3.777 1/h (RSE 9.7%)

    # ----------------------------------------------------------------------
    # Extract-matrix covariate effects (beta in Eq. 4). Section 3.3: "the
    # matrix of the extract introduces significant variability in the subject
    # population for V, k, and k12. No correlations among these parameters
    # were detected."
    #
    # SCALE MATTERS AND IS NOT THE SAME FOR ALL THREE. Section 3.3 states the
    # final model used "lognormal transformation of the data, except for V,
    # where a logit-normal transformation was more appropriated", and the
    # printed final equations make the scale explicit:
    #   Eq. 7   logit(V)   = logit(Vpop)   + beta_V                (no eta on V)
    #   Eq. 8   log(k)     = log(kpop)     + beta_k    + eta_k
    #   Eq. 9   log(k12)   = log(k12pop)   + beta_k12  + eta_k12
    #   Eq. 10  log(k21)   = log(k21pop)               + eta_k21
    # so beta_V is an additive shift on the LOGIT of V (bounded in (0, 1) by
    # Monolix's default logit-normal support, which V in L/kg satisfies),
    # while beta_k and beta_k12 are additive on the log scale.
    # ----------------------------------------------------------------------
    e_form_rutin_extract_vc  <-  0.678 ; label("Extract-matrix shift on logit(V) (logit units)")         # Table 3: beta_V = 0.678 (RSE 9.6%)
    e_form_rutin_extract_kel <-  0.625 ; label("Extract-matrix shift on log(k) (log units)")             # Table 3: beta_k = 0.625 (RSE 19.7%)
    e_form_rutin_extract_k12 <- -0.634 ; label("Extract-matrix shift on log(k12) (log units)")           # Table 3: beta_k12 = -0.634 (RSE 43.6%)

    # ----------------------------------------------------------------------
    # Inter-individual variability. Table 3's sub-heading is "Standard
    # deviation of the Random Effects", so the tabulated omega values are SDs
    # and are squared here because nlmixr2's `~` takes a VARIANCE. No random
    # effect was estimated on V (Eq. 7 carries no eta term and Table 3 lists
    # no Omega_V). No correlations were retained (Section 3.3).
    # ----------------------------------------------------------------------
    etalkel ~ 0.029929  # Table 3: Omega_k   = 0.173 (SD, RSE 23.0%) -> 0.173^2
    etalk12 ~ 0.131044  # Table 3: Omega_k12 = 0.362 (SD, RSE 26.4%) -> 0.362^2
    etalk21 ~ 0.070756  # Table 3: Omega_k21 = 0.266 (SD, RSE 27.3%) -> 0.266^2

    # ----------------------------------------------------------------------
    # Residual error -- proportional only (Section 3.3: "The final model used
    # a proportional error approach"). Monolix's proportional error model is
    # y = f + b * f * e, so the tabulated b is directly the proportional SD.
    # ----------------------------------------------------------------------
    propSd <- 0.076 ; label("Proportional residual error (fraction)")                                    # Table 3: b = 0.076 (RSE 9.8%)
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # V follows Eq. 7 on the logit scale: the population value is shifted by
    # beta_V on logit(V) and back-transformed. Writing it any other way (for
    # instance exp(lvc + beta * COV)) would misstate the extract-arm volume.
    # ----------------------------------------------------------------------
    vcpop   <- exp(lvc)
    logitvc <- log(vcpop / (1 - vcpop)) +
      e_form_rutin_extract_vc * FORM_RUTIN_EXTRACT
    vc <- 1 / (1 + exp(-logitvc))

    kel <- exp(lkel + e_form_rutin_extract_kel * FORM_RUTIN_EXTRACT + etalkel)
    k12 <- exp(lk12 + e_form_rutin_extract_k12 * FORM_RUTIN_EXTRACT + etalk12)
    k21 <- exp(lk21 + etalk21)

    # ----------------------------------------------------------------------
    # 2. Two-compartment disposition with first-order elimination from the
    #    central compartment (Figure 3a). X1 is the amount in the central
    #    compartment and X2 the amount in the peripheral compartment.
    # ----------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------------
    # 3. Observation. Doses are in ug/kg and vc is in L/kg, so central / vc is
    #    in ug/L = ng/mL, the unit the assay reports.
    # ----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
