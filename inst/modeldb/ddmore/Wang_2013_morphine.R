Wang_2013_morphine <- function() {
  description <- "Two-compartment population PK model for morphine across the entire paediatric age range and adults using a bodyweight-dependent allometric exponent (BDE) on clearance, with adolescent-specific intercompartmental clearance and central volume and an adult-stratum oral-bioavailability adjustment, as packaged in DDMORE Foundation Model Repository entry DDMODEL00000269 (Wang 2013 Model I)."
  reference <- paste(
    "Wang C, Sadhasivam S, Krekels EHJ, Dahan A, Tibboel D, Danhof M, Vinks AA, Knibbe CAJ (2013).",
    "Developmental changes in morphine clearance across the entire paediatric age range are best described by a bodyweight-dependent exponent model.",
    "Clinical Drug Investigation 33(7):523-534.",
    "doi:10.1007/s40261-013-0097-6.",
    "PMID:23754691.",
    "DDMORE Foundation Model Repository: DDMODEL00000269 (Model I, morphine alone).",
    sep = " "
  )
  vignette <- "Wang_2013_morphine"
  units <- list(time = "minute", dosing = "ug", concentration = "ug/L")

  ddmore_id    <- "DDMODEL00000269"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline (one row per subject in the bundled simulated dataset; the .mod $INPUT comment labels BW as 'bodyweight in kg').",
        "Drives the bodyweight-dependent allometric exponent (BDE) on morphine clearance:",
        "KBDE(WT) = (KDEC + (KMAX - KDEC)) - KDEC * WT^GAMMA / (KHAL^GAMMA + WT^GAMMA),",
        "with KDEC = 0.594, KMAX - KDEC = 0.872, KHAL = 4.01 kg, GAMMA = 4.62.",
        "Reference weight for allometric scaling is 70 kg (TVCL = THETA(5) * (BW/70)^KBDE in the .mod $PK)."
      ),
      source_name        = "BW"
    ),
    CHILD = list(
      description        = "Age-stratum indicator: 1 if the subject is in the 0-3 year (newborn / infant / young child) stratum, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult baseline together with ADOLESCENT = 0)",
      notes              = paste(
        "Derived from the source POP column documented in the .mod $INPUT comments:",
        "POP = 1 (0-3 years) -> CHILD = 1; POP != 1 -> CHILD = 0.",
        "CHILD = 0 AND ADOLESCENT = 0 identifies the adult stratum (POP = 3, 18-36 years)",
        "that triggers the F1 = 0.88 oral-bioavailability adjustment in the .mod $PK block."
      ),
      source_name        = "POP"
    ),
    ADOLESCENT = list(
      description        = "Age-stratum indicator: 1 if the subject is in the 6-15 year (children-adolescents) stratum, 0 otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-adolescent: 0-3 year paediatric stratum or adult stratum)",
      notes              = paste(
        "Derived from the source POP column documented in the .mod $INPUT comments:",
        "POP = 2 (6-15 years) -> ADOLESCENT = 1; POP != 2 -> ADOLESCENT = 0.",
        "Adolescents take a different intercompartmental clearance Q (lq_adolescent, no BW scaling; THETA(10))",
        "and a different central volume V1 (lvc_adolescent, still BW-scaled; THETA(11))",
        "via the IF (POP.EQ.2) overrides in the .mod $PK block."
      ),
      source_name        = "POP"
    )
  )

  population <- list(
    n_subjects     = 475L,
    n_studies      = "Pooled across multiple paediatric and adult studies (full study list reported in the publication; original Wang 2013 PDF is not on disk under literature/).",
    age_range      = "0-3 years (newborns / infants / young children stratum, POP = 1), 6-15 years (children-adolescents stratum, POP = 2), 18-36 years (adults stratum, POP = 3) per the .mod $INPUT comments.",
    weight_range   = "0.6-85 kg in the bundle's simulated dataset (across the three POP strata).",
    sex_female_pct = "Not extractable from DDMORE bundle (Wang 2013 PDF not on disk under literature/).",
    race_ethnicity = "Not extractable from DDMORE bundle.",
    disease_state  = "Mixed paediatric and adult cohorts receiving morphine. Specific clinical indication not extractable from DDMORE bundle alone; the publication abstract describes 358 neonates / infants / children / adults plus 117 adolescents.",
    dose_range     = "Variable by study; doses in the bundle's simulated dataset are in micrograms with infusion rates ranging from short bolus-equivalent infusions (~10 s) to 60-min infusions (e.g., 1850 ug at 30.83 ug/min in adults; 53 ug at 265 ug/min in 0.6-kg neonates).",
    regions        = "Not extractable from DDMORE bundle.",
    notes          = paste(
      "Demographic counts are summarised from the publication abstract (358 neonates / infants / children / adults + 117 adolescents = 475)",
      "and the .mod $INPUT comments (POP age-stratum definitions).",
      "Original Wang 2013 PDF is not on disk under the literature tree at extraction time;",
      "see the vignette's 'Assumptions and deviations' section for the items that could not be cross-checked against the publication."
    )
  )

  ini({
    # ----- Bodyweight-dependent exponent (BDE) on CL -----
    # KBDE = (KDEC + (KMAX - KDEC)) - KDEC * WT^GAMMA / (KHAL^GAMMA + WT^GAMMA).
    # Mechanistic constants; no IIV (the .mod has no ETA on these).
    # All four are positive, so log-transform for positivity.
    lkdec  <- log(0.594); label("Decrease in BDE clearance exponent across the paediatric BW range (KDEC, unitless)")           # Output_real_ModelI_Morphine.lst FINAL PARAMETER ESTIMATE TH 1
    lkmin  <- log(0.872); label("Adult-asymptote BDE clearance exponent (KMAX - KDEC, unitless)")                                 # Output_real_ModelI_Morphine.lst TH 2
    lkhal  <- log(4.01) ; label("BDE half-maximal-effect bodyweight (KHAL, kg)")                                                  # Output_real_ModelI_Morphine.lst TH 3
    lhill <- log(4.62) ; label("BDE Hill coefficient (GAMMA, unitless)")                                                          # Output_real_ModelI_Morphine.lst TH 4

    # ----- Structural PK typical values (reference 70 kg) -----
    lcl  <- log(1.62); label("Population clearance for a 70-kg adult (CL, L/min)")                                                # Output_real_ModelI_Morphine.lst TH 5
    lq   <- log(1.90); label("Population intercompartmental clearance for a 70-kg subject in non-adolescent strata (Q, L/min)")    # Output_real_ModelI_Morphine.lst TH 6
    lvc  <- log(81.2); label("Population central volume for a 70-kg subject in non-adolescent strata (V1, L)")                     # Output_real_ModelI_Morphine.lst TH 7
    lvp  <- log(128) ; label("Population peripheral volume for a 70-kg subject (V2, L)")                                           # Output_real_ModelI_Morphine.lst TH 8

    # ----- Adolescent-specific (POP = 2) overrides -----
    lq_adolescent  <- log(0.500); label("Adolescent (POP = 2, 6-15 yrs) intercompartmental clearance, NOT BW-scaled (Q, L/min)") # Output_real_ModelI_Morphine.lst TH10
    lvc_adolescent <- log(46.0) ; label("Adolescent (POP = 2, 6-15 yrs) central volume for a 70-kg subject, BW-scaled (V1, L)")   # Output_real_ModelI_Morphine.lst TH11

    # ----- Adult-stratum oral bioavailability -----
    # The .mod sets F1 = 0.88 if POP == 3. Encoded multiplicatively as an adult-stratum
    # ratio applied to subjects with CHILD = 0 AND ADOLESCENT = 0, paralleling the
    # CarlssonPetri_2021_liraglutide age-stratum-ratio encoding.
    e_age_adult_f <- 0.88; label("Adult-stratum (POP = 3, oral-route arm) F1 ratio applied as e_age_adult_f^(1 - CHILD - ADOLESCENT) (unitless)") # .mod $PK: IF (POP.EQ.3) F1 = 0.88

    # ----- Inter-individual variability -----
    # Final OMEGA estimates from the .lst; the .mod has $OMEGA 0 FIX on Q and V2,
    # so only CL and V1 carry IIV.
    etalcl ~ 0.159  # Output_real_ModelI_Morphine.lst OMEGA ETA1, IIV CL
    etalvc ~ 0.253  # Output_real_ModelI_Morphine.lst OMEGA ETA3, IIV V1

    # ----- Residual error -----
    # The .mod $ERROR uses NONMEM log-transform-both-sides:
    #   IPRED = LOG(F); W = THETA(9); Y = IPRED + ERR(1)*W with $SIGMA EPS1 FIXED at 1.
    # NONMEM "additive on log-scale" with SIGMA fixed at 1 maps to proportional residual error
    # in linear nlmixr2 space with propSd = THETA(9). See naming-conventions.md NONMEM-syntax-translation table.
    propSd <- 0.432; label("Proportional residual error (fraction)")  # Output_real_ModelI_Morphine.lst TH 9 (THETA-encoded log-transform-both-sides residual SD)
  })

  model({
    # 1. Bodyweight-dependent allometric exponent (BDE) on CL.
    kdec  <- exp(lkdec)
    kmin  <- exp(lkmin)
    khal  <- exp(lkhal)
    hill <- exp(lhill)
    kbde  <- (kdec + kmin) - kdec * (WT ^ hill) / (khal ^ hill + WT ^ hill)

    # 2. Individual PK parameters with adolescent-specific overrides.
    cl <- exp(lcl + etalcl) * (WT / 70) ^ kbde
    q  <- exp(lq)            * (WT / 70) * (1 - ADOLESCENT) +
          exp(lq_adolescent)               * ADOLESCENT
    vc <- (exp(lvc            + etalvc) * (1 - ADOLESCENT) +
           exp(lvc_adolescent + etalvc) * ADOLESCENT) * (WT / 70)
    vp <- exp(lvp) * (WT / 70)

    # 3. Micro-constants for the 2-compartment ODE system (NONMEM ADVAN5 K10/K12/K21).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system. Dose enters the central compartment ($MODEL COMP=CENTRAL DEFDOSE).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability: F1 = 1 for paediatric strata (CHILD = 1 or ADOLESCENT = 1);
    #    F1 = e_age_adult_f for adults (CHILD = 0 AND ADOLESCENT = 0).
    f(central) <- e_age_adult_f ^ (1 - CHILD - ADOLESCENT)

    # 6. Observation and residual error model (proportional in linear space; see ini() comment).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
