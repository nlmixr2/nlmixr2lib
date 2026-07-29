Parmar_2023_spectinamide_1599_rat_pbpk <- function() {
  description <- paste(
    "Minimal PBPK (mPBPK; author-coded ODEs, fitted in Monolix 2021R1).",
    "Preclinical (Sprague-Dawley rat, 225 g). Spectinamide 1599, an",
    "anti-tuberculosis spectinomycin analogue, after a single intravenous",
    "dose (Parmar 2023, Pharmaceutics). Rat-physiology counterpart of",
    "Parmar_2023_spectinamide_1599_mouse_pbpk: the same structure and the",
    "same jointly estimated distribution rate constants, with the",
    "Sprague-Dawley physiology and rat-specific blood-to-plasma ratio and",
    "unbound fraction substituted. Venous and arterial blood plus five",
    "tissues (lung, spleen, liver, kidney, and a lumped 'other'), each split",
    "into a rapid-equilibrium extracellular pool (vascular + interstitial)",
    "and a slow cellular pool coupled by a first-order influx K(I->C) on the",
    "unbound fraction and a back flux K(C->I). Lung in series between venous",
    "and arterial blood, spleen draining portally into the liver, and",
    "elimination by glomerular filtration (GFR x fu) on the kidney",
    "extracellular pool. Intravenous only: the paper reports no",
    "subcutaneous or aerosol rat data and estimated Ka / F from mouse data",
    "alone, so no absorption depot is included. Two mass-balance defects in",
    "the published equations plus a flow-rounding inconsistency are",
    "corrected here; see the vignette Errata."
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
    time          = "hour",
    dosing        = "mg (absolute amount; the paper's mg/kg dose x 0.225 kg rat body weight, e.g. 10 mg/kg = 2.25 mg)",
    concentration = "mg/L (equivalently ug/mL for plasma and ug/g for tissue under the paper's unit-tissue-density assumption)"
  )

  # No covariates. The rat physiology (Parmar 2023 Table 2, Rat (225 g)
  # column) is held fixed and no demographic or disease covariate effect is
  # estimated. Species is a separate model file rather than a covariate
  # switch, following the An_2012_mitoxantrone_{mouse,human}_pbpk precedent.
  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 19L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = "225 g reference body weight (Parmar 2023 Table 2)",
    sex_female_pct = 52.6,
    race_ethnicity = NA_character_,
    disease_state  = "Healthy (uninfected) double-catheterized Sprague-Dawley rats.",
    dose_range     = "Spectinamide 1599 10 mg/kg single intravenous dose (Parmar 2023 Table 1)",
    regions        = "USA (University of Tennessee Health Science Center)",
    notes          = paste(
      "Two rat datasets (Parmar 2023 Table 1). (i) Serial-sampling plasma",
      "study: 6 female + 5 male rats, 13 plasma samples per animal at 0.08,",
      "0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6, 8, 10, 24 and 48 h, with urine",
      "collected over 0-6, 6-10, 10-24 and 24-48 h for fraction-excreted",
      "calculations. (ii) Destructive-sampling tissue study: 4 male + 4",
      "female rats, groups of four (two male, two female) euthanized at",
      "15 min or 4 h with blood, lung, liver, spleen and kidney collected",
      "(1 data point per animal, 4 animals per sampling time point).",
      "Spectinamide 1599 was formulated at 10 mg/mL in PlasmaLyte / water",
      "(9:1) and injected via the femoral-vein catheter. The distribution",
      "rate constants and inter-animal variances in Table 4 were estimated",
      "by simultaneously fitting the mouse and rat plasma and tissue data,",
      "so they are shared with",
      "Parmar_2023_spectinamide_1599_mouse_pbpk."
    )
  )

  ini({
    # =========================================================================
    # Tissue distribution rate constants (Parmar 2023 Table 4, Intravenous
    # column; verified against the typeset PDF pages 9-10). Estimated by
    # simultaneously fitting the mouse and rat single-dose intravenous plasma
    # and tissue data, so the same values apply to both species files; only
    # the physiology and the drug-specific k(b/p) / fu differ.
    # =========================================================================
    lkic_lung <- log(0.068)
    label("Lung uptake rate constant K_Lung(I->C) (1/h)")
    # Parmar 2023 Table 4 IV: 0.068 1/h (RSE 15.3%).

    lkci_lung <- log(0.028)
    label("Lung cellular back-flux rate constant K_Lung(C->I) (1/h)")
    # Parmar 2023 Table 4 IV: 0.028 1/h (RSE 41.5%).

    lkic_spleen <- log(0.048)
    label("Spleen uptake rate constant K_Spleen(I->C) (1/h)")
    # Parmar 2023 Table 4 IV: 0.048 1/h (RSE 16.5%).

    lkci_spleen <- log(0.01)
    label("Spleen cellular back-flux rate constant K_Spleen(C->I) (1/h)")
    # Parmar 2023 Table 4 IV: 0.01 1/h (RSE 106%).

    lkic_liver <- log(0.87)
    label("Liver uptake rate constant K_Liver(I->C) (1/h)")
    # Parmar 2023 Table 4 IV: 0.87 1/h (RSE 10.1%).

    lkci_liver <- log(0.061)
    label("Liver cellular back-flux rate constant K_Liver(C->I) (1/h)")
    # Parmar 2023 Table 4 IV: 0.061 1/h (RSE 13.7%).

    lkic_kidney <- log(12.1)
    label("Kidney uptake rate constant K_Kidney(I->C) (1/h)")
    # Parmar 2023 Table 4 IV: 12.1 1/h (RSE 19.7%).

    lkci_kidney <- log(0.15)
    label("Kidney cellular back-flux rate constant K_Kidney(C->I) (1/h)")
    # Parmar 2023 Table 4 IV: 0.15 1/h (RSE 31.0%).

    lkic_other <- log(5.4)
    label("Other-tissue uptake rate constant K_Others(I->C) (1/h)")
    # Parmar 2023 Table 4 IV: 5.4 1/h (RSE 4.92%).

    lkci_other <- log(7.0e-5)
    label("Other-tissue cellular back-flux rate constant K_Others(C->I) (1/h)")
    # Parmar 2023 Table 4 IV: 7.0 x 10^-5 1/h (RSE 142%).

    # =========================================================================
    # Inter-animal variability on the four tissue uptake rate constants
    # (Parmar 2023 Table 4, Intravenous column). The tabulated omega values
    # are log-scale STANDARD DEVIATIONS, not variances: the paper's own
    # Section 3.7 states that 0.48-0.87 "correspond to 50.9 to 106% CV", and
    # sqrt(exp(0.48^2) - 1) = 50.9% while sqrt(exp(0.87^2) - 1) = 106% (the
    # variance reading would give 78.5-106% instead). nlmixr2 ini() takes
    # variances, so each entry is the published SD squared.
    # =========================================================================
    etalkic_lung ~ 0.7569
    # Parmar 2023 Table 4 IV: omega K_Lung(I->C) = 0.87 (RSE 30.6%) -> 0.87^2.

    etalkic_spleen ~ 0.4356
    # Parmar 2023 Table 4 IV: omega K_Spleen(I->C) = 0.66 (RSE 41.0%) -> 0.66^2.

    etalkic_liver ~ 0.3721
    # Parmar 2023 Table 4 IV: omega K_Liver(I->C) = 0.61 (RSE 29.8%) -> 0.61^2.

    etalkic_other ~ 0.2304
    # Parmar 2023 Table 4 IV: omega K_Others(I->C) = 0.48 (RSE 21.6%) -> 0.48^2.

    # =========================================================================
    # Residual unexplained variability (Parmar 2023 Table 4, Intravenous
    # column -- the only route with rat data). Proportional per matrix. No
    # residual error is reported for kidney, so Ckidney is exposed as a
    # derived concentration without an error model.
    # =========================================================================
    propSd <- 0.32
    label("Proportional residual SD for venous plasma Cc (fraction)")
    # Parmar 2023 Table 4: eps_Plasma IV 0.32 (RSE 15.1%).

    propSd_Clung <- 0.35
    label("Proportional residual SD for lung concentration (fraction)")
    # Parmar 2023 Table 4: eps_Lung IV 0.35 (RSE 14.4%).

    propSd_Cliver <- 0.28
    label("Proportional residual SD for liver concentration (fraction)")
    # Parmar 2023 Table 4: eps_Liver IV 0.28 (RSE 14.4%).

    propSd_Cspleen <- 0.53
    label("Proportional residual SD for spleen concentration (fraction)")
    # Parmar 2023 Table 4: eps_Spleen IV 0.53 (RSE 14.4%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Sprague-Dawley rat physiology, held fixed (Parmar 2023 Table 2,
    #    Rat (225 g) column; flows and GFR from ref [16], volumes from
    #    ref [17]). Flows in L/h, volumes in L.
    #
    #    q_lung is DERIVED as the sum of the systemic organ flows
    #    (0.901 + 0.0412 + 0.601 + 3.29 = 4.8332 L/h) rather than taken as
    #    the printed 4.83 L/h. Pulmonary flow is the cardiac output and must
    #    equal the sum exactly, otherwise the arterial node creates or
    #    destroys drug: with the printed 4.83 the sum exceeds it by 0.066%,
    #    which is a small fraction of cardiac output but a material fraction
    #    of total drug removal because glomerular filtration is itself a
    #    small clearance. See the vignette Errata and mass-balance check.
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

    # Drug-specific parameters for spectinamide 1599 in the rat (Parmar 2023
    # Table 2, "Generated in this study as described under Methods").
    kbp <- 0.812  # Parmar 2023 Table 2, rat spectinamide 1599 k(b/p)
    fu  <- 0.563  # Parmar 2023 Table 2, rat spectinamide 1599 fu Plasma

    # -----------------------------------------------------------------------
    # 2. Tissue sub-compartment volume fractions (Parmar 2023 Table 3, from
    #    ref [17]; identical in mouse and rat). The rapid-equilibrium
    #    ("extracellular") pool lumps the vascular and interstitial spaces,
    #    assumed to be in instantaneous equilibrium with blood, so a single
    #    blood-scale concentration Cb describes both with effective volume
    #    V_vascular + V_interstitial / k(b/p) (supplement S3-S18 LHS).
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
    # 3. Individual parameters. Inter-animal variability sits on the four
    #    tissue uptake rate constants for which Table 4 reports an omega.
    # -----------------------------------------------------------------------
    kic_lung   <- exp(lkic_lung + etalkic_lung)
    kci_lung   <- exp(lkci_lung)
    kic_spleen <- exp(lkic_spleen + etalkic_spleen)
    kci_spleen <- exp(lkci_spleen)
    kic_liver  <- exp(lkic_liver + etalkic_liver)
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
    # 5. Blood compartments. Supplement S1 is printed as
    #      V_ven dC_ven/dt = Q_Liver Cb_Liver + Q_Kidney Cb_Kidney
    #                        + Q_Other Cb_Other - Q_Lung Cb_Lung
    #    Two corrections are applied (see the vignette Errata):
    #      (a) the venous LOSS to the lung must be Q_Lung x C_VenousBlood,
    #          not Q_Lung x Cb_Lung, because S3 has the lung GAINING
    #          Q_Lung x C_VenousBlood;
    #      (b) the hepatic venous return must carry the portal spleen flow,
    #          (Q_Liver + Q_Spleen) x Cb_Liver, because S10 adds
    #          Q_Spleen x Cb_Spleen to the liver while removing only
    #          Q_Liver x Cb_Liver.
    #    Both follow uniquely from mass balance and from the tabulated-flow
    #    identity in block 1.
    # -----------------------------------------------------------------------
    d/dt(venous) <- (q_liver + q_spleen) * Cb_liver +
      q_kidney * Cb_kidney +
      q_other * Cb_other -
      q_lung * Cven
    d/dt(arterial) <- q_lung * (Cb_lung - Cart)  # supplement S6

    # -----------------------------------------------------------------------
    # 6. Lung, in series between venous and arterial blood (supplement
    #    S3 / S4). Because K(I->C) multiplies Cb x V_effective it reduces to
    #    K(I->C) x fu x <organ>_extracellular in amount form, and
    #    K(C->I) x Ccell x V_cellular to K(C->I) x <organ>_cellular.
    # -----------------------------------------------------------------------
    d/dt(lung_extracellular) <- q_lung * (Cven - Cb_lung) -
      kic_lung * fu * lung_extracellular +
      kci_lung * lung_cellular
    d/dt(lung_cellular) <- kic_lung * fu * lung_extracellular -
      kci_lung * lung_cellular

    # -----------------------------------------------------------------------
    # 7. Spleen (supplement S7 / S8). Perfused from arterial blood; its
    #    venous outflow drains portally into the liver.
    # -----------------------------------------------------------------------
    d/dt(spleen_extracellular) <- q_spleen * (Cart - Cb_spleen) -
      kic_spleen * fu * spleen_extracellular +
      kci_spleen * spleen_cellular
    d/dt(spleen_cellular) <- kic_spleen * fu * spleen_extracellular -
      kci_spleen * spleen_cellular

    # -----------------------------------------------------------------------
    # 8. Liver (supplement S10 / S11). Hepatic-arterial inflow q_liver plus
    #    the portal spleen outflow q_spleen x Cb_spleen; returns
    #    (q_liver + q_spleen) x Cb_liver to venous blood.
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
    #     (supplement S2: C_VenousPlasma = C_VenousBlood / k(b/p)); the
    #     tissue concentrations are supplement S5, S9, S12 and S15.
    #     Proportional residual error per matrix; kidney has no reported
    #     residual so Ckidney carries no error model.
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

    Cc      ~ prop(propSd)
    Clung   ~ prop(propSd_Clung)
    Cliver  ~ prop(propSd_Cliver)
    Cspleen ~ prop(propSd_Cspleen)
  })
}
