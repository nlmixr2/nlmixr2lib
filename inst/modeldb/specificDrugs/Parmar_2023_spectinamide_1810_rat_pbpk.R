Parmar_2023_spectinamide_1810_rat_pbpk <- function() {
  description <- paste(
    "Minimal PBPK (mPBPK; author-coded ODEs, fitted in Monolix 2021R1).",
    "Preclinical (Sprague-Dawley rat, 225 g). Spectinamide 1810 after a",
    "single intravenous dose (Parmar 2023, Pharmaceutics). Rat-physiology",
    "counterpart of Parmar_2023_spectinamide_1810_mouse_pbpk: the same",
    "structure and the same distribution rate constants, with the",
    "Sprague-Dawley physiology and rat-specific blood-to-plasma ratio and",
    "unbound fraction substituted. Venous and arterial blood plus five",
    "tissues (lung, spleen, liver, kidney, and a lumped 'other'), each split",
    "into a rapid-equilibrium extracellular pool (vascular + interstitial)",
    "and a slow cellular pool coupled by a first-order influx K(I->C) on the",
    "unbound fraction and a back flux K(C->I). Lung in series between venous",
    "and arterial blood, spleen draining portally into the liver, and",
    "elimination by glomerular filtration (GFR x fu) on the kidney",
    "extracellular pool. Intravenous only, and plasma is the only endpoint",
    "with a reported residual error: the rat study sampled plasma alone. Two",
    "mass-balance defects in the published equations plus a flow-rounding",
    "inconsistency are corrected here; see the vignette Errata."
  )
  reference <- paste(
    "Parmar KR, Lukka PB, Wagh S, Temrikar ZH, Liu J, Lee RE, Braunstein M,",
    "Hickey AJ, Robertson GT, Gonzalez-Juarrero M, Edginton A, Meibohm B.",
    "(2023). Development of a Minimalistic Physiologically Based",
    "Pharmacokinetic (mPBPK) Model for the Preclinical Development of",
    "Spectinamide Antibiotics.",
    "Pharmaceutics 15(6):1759. doi:10.3390/pharmaceutics15061759.",
    "PMCID PMC10305115.",
    sep = " "
  )
  vignette <- "Parmar_2023_spectinamides_pbpk"
  units <- list(
    time          = "h",
    dosing        = "mg (absolute amount; the paper's mg/kg dose x 0.225 kg rat body weight, e.g. 10 mg/kg = 2.25 mg)",
    concentration = "mg/L (equivalently ug/mL for plasma and ug/g for tissue under the paper's unit-tissue-density assumption)"
  )

  # No covariates. The rat physiology (Parmar 2023 Table 2, Rat (225 g)
  # column) is held fixed and no demographic or disease covariate effect is
  # estimated. Species is a separate model file rather than a covariate
  # switch, following the An_2012_mitoxantrone_{mouse,human}_pbpk precedent.
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    venous               = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "blood cell", verified = FALSE),
    arterial             = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "blood cell", verified = FALSE),
    lung_extracellular   = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    lung_cellular        = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen_extracellular = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen_cellular      = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    liver_extracellular  = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    liver_cellular       = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    kidney_extracellular = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    kidney_cellular      = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    other_extracellular  = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE),
    other_cellular       = list(analyte = "Spectinamide 1810", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley, male)",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = "225 g reference body weight (Parmar 2023 Table 2)",
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy (uninfected) double-catheterized male Sprague-Dawley rats.",
    dose_range     = "Spectinamide 1810 10 mg/kg single intravenous dose (Parmar 2023 Table 1)",
    regions        = "USA (University of Tennessee Health Science Center)",
    notes          = paste(
      "18 male Sprague-Dawley rats across three separate studies, dosed at",
      "10 mg/kg intravenously via the femoral-vein catheter (10 mg/mL in",
      "PlasmaLyte / water 9:1) with serial jugular sampling at 0.08, 0.25,",
      "0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 10, 24 and 48 h -- 13 plasma data",
      "points per animal, 18 animals per sampling time point (Parmar 2023",
      "Table 1, 'Generated in this study as described'). Urine was collected",
      "over 0-6, 6-10, 10-24 and 24-48 h for fraction-excreted calculations.",
      "No rat tissue samples were taken for spectinamide 1810, so plasma is",
      "the only endpoint carrying a residual-error model; the tissue",
      "concentrations remain available as derived outputs. Parmar 2023",
      "Section 3.8 uses this dataset to check the species scaling of the",
      "mouse model (Figure 13, Table 8 last row)."
    )
  )

  ini({
    # =========================================================================
    # Tissue distribution rate constants (Parmar 2023 Table 7, Intravenous
    # column). Shared with Parmar_2023_spectinamide_1810_mouse_pbpk; only the
    # physiology and the drug-specific k(b/p) / fu differ between the two
    # species files.
    # =========================================================================
    lkic_lung <- log(0.13)
    label("Lung uptake rate constant K_Lung(I->C) (1/h)")
    # Parmar 2023 Table 7 IV: 0.13 1/h (RSE 16.0%).

    lkci_lung <- log(0.076)
    label("Lung cellular back-flux rate constant K_Lung(C->I) (1/h)")
    # Parmar 2023 Table 7 IV: 0.076 1/h (RSE 17.6%).

    lkic_spleen <- log(0.059)
    label("Spleen uptake rate constant K_Spleen(I->C) (1/h)")
    # Parmar 2023 Table 7 IV: 0.059 1/h (RSE 23.1%).

    lkci_spleen <- log(0.013)
    label("Spleen cellular back-flux rate constant K_Spleen(C->I) (1/h)")
    # Parmar 2023 Table 7 IV: 0.013 1/h (RSE 120%).

    lkic_liver <- log(1.19)
    label("Liver uptake rate constant K_Liver(I->C) (1/h)")
    # Parmar 2023 Table 7 IV: 1.19 1/h (RSE 12.0%).

    lkci_liver <- log(0.051)
    label("Liver cellular back-flux rate constant K_Liver(C->I) (1/h)")
    # Parmar 2023 Table 7 IV: 0.051 1/h (RSE 17.2%).

    lkic_kidney <- log(3.94)
    label("Kidney uptake rate constant K_Kidney(I->C) (1/h)")
    # Parmar 2023 Table 7 IV: 3.94 1/h (RSE 43.9%).

    lkci_kidney <- log(0.097)
    label("Kidney cellular back-flux rate constant K_Kidney(C->I) (1/h)")
    # Parmar 2023 Table 7 IV: 0.097 1/h (RSE 45.1%).

    lkic_other <- log(4.77)
    label("Other-tissue uptake rate constant K_Others(I->C) (1/h)")
    # Parmar 2023 Table 7 IV: 4.77 1/h (RSE 7.88%).

    lkci_other <- log(7.0e-5)
    label("Other-tissue cellular back-flux rate constant K_Others(C->I) (1/h)")
    # Parmar 2023 Table 7 IV: 7.0 x 10^-5 1/h (RSE 0.00278%).

    # =========================================================================
    # Inter-animal variability. Table 7 reports a single omega for
    # spectinamide 1810, on the other-tissue uptake constant, estimated by
    # simultaneously fitting the mouse and rat data. As for 1599 the
    # tabulated value is a log-scale STANDARD DEVIATION rather than a
    # variance (Parmar 2023 Section 3.7 maps 0.48-0.87 to 50.9-106% CV, which
    # only holds under the SD reading); nlmixr2 ini() takes variances, so the
    # entry below is the published SD squared.
    # =========================================================================
    etalkic_other ~ 0.0961
    # Parmar 2023 Table 7 IV: omega K_Others(I->C) = 0.31 (RSE 14.7%) -> 0.31^2.

    # =========================================================================
    # Residual unexplained variability. Only plasma was sampled in the rat,
    # and the plasma intravenous residual in Table 7 is clean (unlike the
    # three tissue rows, which are transcription duplicates of the
    # rate-constant rows above them -- see
    # Parmar_2023_spectinamide_1810_mouse_pbpk and the vignette Errata).
    # =========================================================================
    propSd <- 0.43
    label("Proportional residual SD for venous plasma Cc (fraction)")
    # Parmar 2023 Table 7: eps_Plasma IV 0.43 (RSE 14.7%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Sprague-Dawley rat physiology, held fixed (Parmar 2023 Table 2,
    #    Rat (225 g) column; flows and GFR from ref [16], volumes from
    #    ref [17]). Flows in L/h, volumes in L.
    #
    #    q_lung is DERIVED as the sum of the systemic organ flows
    #    (0.901 + 0.0412 + 0.601 + 3.29 = 4.8332 L/h) rather than taken as
    #    the printed 4.83 L/h: pulmonary flow is the cardiac output and must
    #    equal that sum exactly, otherwise the arterial node creates drug.
    #    See the vignette Errata and mass-balance check.
    # -----------------------------------------------------------------------
    q_spleen <- 0.0412   # Parmar 2023 Table 2, rat Q_Spleen (L/h)
    q_liver  <- 0.901    # Parmar 2023 Table 2, rat Q_Liver (L/h)
    q_kidney <- 0.601    # Parmar 2023 Table 2, rat Q_Kidney (L/h)
    q_other  <- 3.29     # Parmar 2023 Table 2, rat Q_Other (L/h)
    q_lung   <- q_liver + q_spleen + q_kidney + q_other  # = 4.8332 L/h;
    # Parmar 2023 Table 2 prints rat Q_Lung = 4.83 L/h (0.066% lower).
    gfr      <- 0.088    # Parmar 2023 Table 2, rat GFR (L/h)

    v_venous   <- 0.0115    # Parmar 2023 Table 2, rat V_Venous blood (L)
    v_arterial <- 0.00494   # Parmar 2023 Table 2, rat V_Arterial blood (L)
    v_lung     <- 0.00140   # Parmar 2023 Table 2, rat V_Lung (L)
    v_spleen   <- 0.00277   # Parmar 2023 Table 2, rat V_Spleen (L)
    v_liver    <- 0.0157    # Parmar 2023 Table 2, rat V_Liver (L)
    v_kidney   <- 0.00241   # Parmar 2023 Table 2, rat V_Kidney (L)
    v_other    <- 0.245     # Parmar 2023 Table 2, rat V_Other (L)

    # Drug-specific parameters for spectinamide 1810 in the rat (Parmar 2023
    # Table 2, "Generated in this study as described under Methods").
    kbp <- 0.785  # Parmar 2023 Table 2, rat spectinamide 1810 k(b/p)
    fu  <- 0.607  # Parmar 2023 Table 2, rat spectinamide 1810 fu Plasma

    # -----------------------------------------------------------------------
    # 2. Tissue sub-compartment volume fractions (Parmar 2023 Table 3;
    #    identical in mouse and rat and across the two compounds). The
    #    rapid-equilibrium ("extracellular") pool lumps the vascular and
    #    interstitial spaces, assumed to be in instantaneous equilibrium with
    #    blood, so a single blood-scale concentration Cb describes both with
    #    effective volume V_vascular + V_interstitial / k(b/p).
    # -----------------------------------------------------------------------
    fvas_lung   <- 0.26   # Parmar 2023 Table 3, lung fraction vascular
    fint_lung   <- 0.19   # Parmar 2023 Table 3, lung fraction interstitial
    fcell_lung  <- 0.55   # Parmar 2023 Table 3, lung fraction cellular
    fvas_spleen  <- 0.22  # Parmar 2023 Table 3, spleen fraction vascular
    fint_spleen  <- 0.20  # Parmar 2023 Table 3, spleen fraction interstitial
    fcell_spleen <- 0.58  # Parmar 2023 Table 3, spleen fraction cellular
    fvas_liver   <- 0.15  # Parmar 2023 Table 3, liver fraction vascular
    fint_liver   <- 0.20  # Parmar 2023 Table 3, liver fraction interstitial
    fcell_liver  <- 0.64  # Parmar 2023 Table 3, liver fraction cellular
    fvas_kidney  <- 0.10  # Parmar 2023 Table 3, kidney fraction vascular
    fint_kidney  <- 0.15  # Parmar 2023 Table 3, kidney fraction interstitial
    fcell_kidney <- 0.75  # Parmar 2023 Table 3, kidney fraction cellular
    fvas_other   <- 0.040 # Parmar 2023 Table 3, other fraction vascular
    fint_other   <- 0.19  # Parmar 2023 Table 3, other fraction interstitial
    fcell_other  <- 0.77  # Parmar 2023 Table 3, other fraction cellular

    veff_lung   <- v_lung   * (fvas_lung   + fint_lung   / kbp)
    veff_spleen <- v_spleen * (fvas_spleen + fint_spleen / kbp)
    veff_liver  <- v_liver  * (fvas_liver  + fint_liver  / kbp)
    veff_kidney <- v_kidney * (fvas_kidney + fint_kidney / kbp)
    veff_other  <- v_other  * (fvas_other  + fint_other  / kbp)

    # -----------------------------------------------------------------------
    # 3. Individual parameters. Inter-animal variability sits on the
    #    other-tissue uptake rate constant, the only one for which Table 7
    #    reports an omega.
    # -----------------------------------------------------------------------
    kic_lung   <- exp(lkic_lung)
    kci_lung   <- exp(lkci_lung)
    kic_spleen <- exp(lkic_spleen)
    kci_spleen <- exp(lkci_spleen)
    kic_liver  <- exp(lkic_liver)
    kci_liver  <- exp(lkci_liver)
    kic_kidney <- exp(lkic_kidney)
    kci_kidney <- exp(lkci_kidney)
    kic_other  <- exp(lkic_other + etalkic_other)
    kci_other  <- exp(lkci_other)

    # -----------------------------------------------------------------------
    # 4. Concentrations. States are amounts (mg); Cb_<organ> is the
    #    blood-scale concentration of the rapid-equilibrium pool.
    # -----------------------------------------------------------------------
    Cven <- venous   / v_venous
    Cart <- arterial / v_arterial
    Cb_lung   <- lung_extracellular   / veff_lung
    Cb_spleen <- spleen_extracellular / veff_spleen
    Cb_liver  <- liver_extracellular  / veff_liver
    Cb_kidney <- kidney_extracellular / veff_kidney
    Cb_other  <- other_extracellular  / veff_other

    # -----------------------------------------------------------------------
    # 5. Blood compartments, with the two supplement-S1 / S10 mass-balance
    #    corrections: the venous loss to the lung is Q_Lung x C_VenousBlood
    #    (S3 has the lung gaining exactly that), and the hepatic venous
    #    return carries the portal spleen flow, (Q_Liver + Q_Spleen) x
    #    Cb_Liver. See the vignette Errata.
    # -----------------------------------------------------------------------
    d/dt(venous) <- (q_liver + q_spleen) * Cb_liver +
      q_kidney * Cb_kidney +
      q_other * Cb_other -
      q_lung * Cven
    d/dt(arterial) <- q_lung * (Cb_lung - Cart)  # supplement S6

    # -----------------------------------------------------------------------
    # 6. Lung, in series between venous and arterial blood (supplement
    #    S3 / S4). K(I->C) x fu x Cb x V_effective reduces to
    #    K(I->C) x fu x <organ>_extracellular in amount form, and
    #    K(C->I) x Ccell x V_cellular to K(C->I) x <organ>_cellular.
    # -----------------------------------------------------------------------
    d/dt(lung_extracellular) <- q_lung * (Cven - Cb_lung) -
      kic_lung * fu * lung_extracellular +
      kci_lung * lung_cellular
    d/dt(lung_cellular) <- kic_lung * fu * lung_extracellular -
      kci_lung * lung_cellular

    # -----------------------------------------------------------------------
    # 7. Spleen (supplement S7 / S8), draining portally into the liver.
    # -----------------------------------------------------------------------
    d/dt(spleen_extracellular) <- q_spleen * (Cart - Cb_spleen) -
      kic_spleen * fu * spleen_extracellular +
      kci_spleen * spleen_cellular
    d/dt(spleen_cellular) <- kic_spleen * fu * spleen_extracellular -
      kci_spleen * spleen_cellular

    # -----------------------------------------------------------------------
    # 8. Liver (supplement S10 / S11).
    # -----------------------------------------------------------------------
    d/dt(liver_extracellular) <- q_liver * Cart +
      q_spleen * Cb_spleen -
      (q_liver + q_spleen) * Cb_liver -
      kic_liver * fu * liver_extracellular +
      kci_liver * liver_cellular
    d/dt(liver_cellular) <- kic_liver * fu * liver_extracellular -
      kci_liver * liver_cellular

    # -----------------------------------------------------------------------
    # 9. Kidney (supplement S13 / S14). The only elimination pathway:
    #    glomerular filtration of unbound drug from the rapid-equilibrium
    #    pool (Parmar 2023 Section 2.8).
    # -----------------------------------------------------------------------
    d/dt(kidney_extracellular) <- q_kidney * (Cart - Cb_kidney) -
      kic_kidney * fu * kidney_extracellular +
      kci_kidney * kidney_cellular -
      gfr * fu * Cb_kidney
    d/dt(kidney_cellular) <- kic_kidney * fu * kidney_extracellular -
      kci_kidney * kidney_cellular

    # -----------------------------------------------------------------------
    # 10. Lumped remainder tissue (supplement S16 / S17).
    # -----------------------------------------------------------------------
    d/dt(other_extracellular) <- q_other * (Cart - Cb_other) -
      kic_other * fu * other_extracellular +
      kci_other * other_cellular
    d/dt(other_cellular) <- kic_other * fu * other_extracellular -
      kci_other * other_cellular

    # -----------------------------------------------------------------------
    # 11. Observations. Cc is the venous PLASMA concentration
    #     (supplement S2: C_VenousPlasma = C_VenousBlood / k(b/p)) and is the
    #     only endpoint with a residual-error model, matching the rat study's
    #     plasma-only sampling. The tissue concentrations (supplement S5, S9,
    #     S12, S15) are derived outputs.
    # -----------------------------------------------------------------------
    # Renal elimination rate (mg/h), exposed as a derived output so a
    # simulation can verify mass balance without re-hardcoding GFR and fu:
    # the amount remaining in all states plus the time integral of elimRate
    # must equal the absorbed dose. See the vignette mass-balance check.
    elimRate <- gfr * fu * Cb_kidney

    Cc      <- Cven / kbp
    Clung   <- (lung_extracellular + lung_cellular) / v_lung
    Cspleen <- (spleen_extracellular + spleen_cellular) / v_spleen
    Cliver  <- (liver_extracellular + liver_cellular) / v_liver
    Ckidney <- (kidney_extracellular + kidney_cellular) / v_kidney

    Cc ~ prop(propSd)
  })
}
