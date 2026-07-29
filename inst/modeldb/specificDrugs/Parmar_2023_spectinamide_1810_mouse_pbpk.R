Parmar_2023_spectinamide_1810_mouse_pbpk <- function() {
  description <- paste(
    "Minimal PBPK (mPBPK; author-coded ODEs, fitted in Monolix 2021R1).",
    "Preclinical (BALB/c mouse, 20 g). Spectinamide 1810, a structurally",
    "close analogue of spectinamide 1599, after intravenous and subcutaneous",
    "dosing (Parmar 2023, Pharmaceutics). Same structure as",
    "Parmar_2023_spectinamide_1599_mouse_pbpk with the drug-specific",
    "blood-to-plasma ratio and unbound fraction updated and the distribution",
    "rate constants, Ka, F and residual errors re-estimated for 1810",
    "(Table 7). Venous and arterial blood plus five tissues (lung, spleen,",
    "liver, kidney, and a lumped 'other'), each split into a",
    "rapid-equilibrium extracellular pool (vascular + interstitial) and a",
    "slow cellular pool coupled by a first-order influx K(I->C) on the",
    "unbound fraction and a back flux K(C->I). Lung in series between venous",
    "and arterial blood, spleen draining portally into the liver, and",
    "elimination by glomerular filtration (GFR x fu) on the kidney",
    "extracellular pool. Subcutaneous doses reach venous blood through depot",
    "with ka / F. No intrapulmonary-aerosol route: the paper studied 1810 by",
    "the intravenous and subcutaneous routes only. Two mass-balance defects",
    "in the published equations plus a flow-rounding inconsistency are",
    "corrected here, and Table 7's three intravenous tissue residual errors",
    "are unrecoverable transcription duplicates so they are fixed at zero;",
    "see the vignette Errata."
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
    dosing        = "mg (absolute amount; the paper's mg/kg dose x 0.02 kg mouse body weight, e.g. 10 mg/kg = 0.2 mg)",
    concentration = "mg/L (equivalently ug/mL for plasma and ug/g for tissue under the paper's unit-tissue-density assumption)"
  )

  # No covariates: physiology is tabulated per species and held fixed, and no
  # demographic, weight or disease-status effect is estimated (Parmar 2023
  # Section 3.8 predicts the infected-mouse plasma profiles with the
  # healthy-animal parameters unchanged). The rat counterpart is a separate
  # model file (Parmar_2023_spectinamide_1810_rat_pbpk.R).
  covariateData <- list()

  population <- list(
    species        = "mouse (BALB/c, 20 g)",
    n_subjects     = 303L,
    n_studies      = 3L,
    age_range      = NA_character_,
    weight_range   = "20 g reference body weight (Parmar 2023 Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Healthy BALB/c mice for the model-parameterisation and",
      "model-qualification data (plasma for all routes; lung, liver and",
      "spleen for the intravenous studies).",
      "Mycobacterium tuberculosis-infected BALB/c mice (studies 2A and 2B,",
      "plasma only) were used as an external predictive check of the",
      "subcutaneous model; Parmar 2023 found the healthy-animal parameters",
      "predicted the infected-animal profiles, so no disease-status effect",
      "is encoded."
    ),
    dose_range     = paste(
      "Intravenous 10 mg/kg single dose and QD5; subcutaneous 46, 50 and",
      "200 mg/kg single dose and 50 / 200 mg/kg QD5 (Parmar 2023 Table 1).",
      "Infected-mouse subcutaneous regimens spanned 10-500 mg/kg per dose",
      "at BID / QD / TIW / BIW / QW for 4 weeks."
    ),
    regions        = "USA (University of Tennessee Health Science Center; Colorado State University BSL-3 for infected animals)",
    notes          = paste(
      "n_subjects counts the BALB/c mice in Parmar 2023 Table 1 that",
      "contributed spectinamide 1810 data: 153 healthy mice (IV single dose",
      "24, IV QD5 24, SC 46 mg/kg 21, SC 50 / 200 mg/kg 48, SC QD5 36) plus",
      "150 Mtb-infected mice across studies 2A and 2B for the external",
      "subcutaneous check. Destructive sampling with 1 data point per animal",
      "and 3 animals per sampling time point for the healthy studies; 2",
      "plasma points per animal and 5-6 animals per time point for the",
      "infected studies. The intravenous studies collected plasma, lung,",
      "liver and spleen; the subcutaneous studies collected plasma only.",
      "Structural parameters are the mouse + rat joint fit reported in",
      "Table 7."
    )
  )

  ini({
    # =========================================================================
    # Tissue distribution rate constants (Parmar 2023 Table 7, Intravenous
    # column). Estimated from the healthy-mouse and rat single-dose
    # intravenous plasma and tissue data for spectinamide 1810, then held
    # fixed for the subcutaneous fit ("Fixed" in the Subcutaneous column),
    # so they are the same parameters here rather than route-specific ones.
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
    # Subcutaneous absorption (Parmar 2023 Table 7, Subcutaneous column;
    # Section 2.8: "the dose was administered as an external compartment,
    # absorbed into the venous compartment via a first-order absorption rate
    # constant and a bioavailability component").
    # =========================================================================
    lka <- log(8.26)
    label("Subcutaneous first-order absorption rate constant Ka (1/h)")
    # Parmar 2023 Table 7 Subcutaneous: Ka = 8.26 (RSE 14.5%).

    lfdepot <- log(1.00)
    label("Subcutaneous bioavailability F (fraction)")
    # Parmar 2023 Table 7 Subcutaneous: F = 1.00 (RSE 0.463%). Reported at
    # the upper boundary of the F parameterisation; kept as estimated (not
    # fixed) because the table gives it an RSE.

    # =========================================================================
    # Inter-animal variability. For spectinamide 1810 Table 7 reports a
    # single omega, on the other-tissue uptake constant. As for 1599 the
    # tabulated omega is a log-scale STANDARD DEVIATION rather than a
    # variance (Parmar 2023 Section 3.7 maps 0.48-0.87 to 50.9-106% CV for
    # the 1599 fit, which only holds under the SD reading); nlmixr2 ini()
    # takes variances, so the entry below is the published SD squared.
    # =========================================================================
    etalkic_other ~ 0.0961
    # Parmar 2023 Table 7 IV: omega K_Others(I->C) = 0.31 (RSE 14.7%) -> 0.31^2.

    # =========================================================================
    # Residual unexplained variability. The bare parameters carry the
    # INTRAVENOUS column of Parmar 2023 Table 7, which is the dataset that
    # parameterised the structural model.
    #
    # PUBLISHER TRANSCRIPTION ERROR: the three intravenous TISSUE residuals
    # in Table 7 are printed as eps_Lung 0.13 (16.0), eps_Liver 0.076 (17.6)
    # and eps_Spleen 1.19 (12.0) -- verbatim duplicates, value AND %RSE, of
    # the K_Lung(I->C) 0.13 (16.0), K_Lung(C->I) 0.076 (17.6) and
    # K_Liver(I->C) 1.19 (12.0) rows above them. Confirmed against the
    # typeset PDF (page 20), and a 119% proportional error on spleen
    # alongside a 43% plasma error is not credible. The three values are
    # therefore treated as UNREPORTED and fixed at zero, per the standing
    # "unreported residual variability -> fixed(0) + erratum; never invent
    # variances" rule. A user wanting stochastic tissue simulations can
    # override them with the clean SUBCUTANEOUS-column values noted on each
    # line. The intravenous plasma residual (0.43) is clean and is used as
    # printed.
    # =========================================================================
    propSd <- 0.43
    label("Proportional residual SD for venous plasma Cc (fraction)")
    # Parmar 2023 Table 7: eps_Plasma IV 0.43 (RSE 14.7%);
    # subcutaneous 0.36 (15.6%).

    propSd_Clung <- fixed(0)
    label("Proportional residual SD for lung concentration (FIXED AT ZERO - published value is a transcription duplicate)")
    # Parmar 2023 Table 7 eps_Lung IV is printed as 0.13 (16.0), identical to
    # the K_Lung(I->C) row; unrecoverable. Subcutaneous column: 0.30 (14.7%).

    propSd_Cliver <- fixed(0)
    label("Proportional residual SD for liver concentration (FIXED AT ZERO - published value is a transcription duplicate)")
    # Parmar 2023 Table 7 eps_Liver IV is printed as 0.076 (17.6), identical
    # to the K_Lung(C->I) row; unrecoverable. Subcutaneous column: 0.31 (14.7%).

    propSd_Cspleen <- fixed(0)
    label("Proportional residual SD for spleen concentration (FIXED AT ZERO - published value is a transcription duplicate)")
    # Parmar 2023 Table 7 eps_Spleen IV is printed as 1.19 (12.0), identical
    # to the K_Liver(I->C) row; unrecoverable. Subcutaneous column: 0.43 (14.7%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Mouse physiology, held fixed (Parmar 2023 Table 2, Mouse (20 g)
    #    column; flows and GFR from ref [16], volumes from ref [17]).
    #    Flows in L/h, volumes in L.
    #
    #    q_lung is DERIVED as the sum of the systemic organ flows
    #    (0.139 + 0.00695 + 0.100 + 0.371 = 0.61695 L/h) rather than taken as
    #    the printed 0.618 L/h. Pulmonary flow is the cardiac output and must
    #    equal that sum exactly, otherwise the arterial node destroys the
    #    difference; with the printed value 2.8% of the dose goes missing by
    #    72 h. See the vignette Errata and mass-balance check.
    # -----------------------------------------------------------------------
    q_spleen <- 0.00695  # Parmar 2023 Table 2, mouse Q_Spleen (L/h)
    q_liver  <- 0.139    # Parmar 2023 Table 2, mouse Q_Liver (L/h)
    q_kidney <- 0.100    # Parmar 2023 Table 2, mouse Q_Kidney (L/h)
    q_other  <- 0.371    # Parmar 2023 Table 2, mouse Q_Other (L/h)
    q_lung   <- q_liver + q_spleen + q_kidney + q_other  # = 0.61695 L/h;
    # Parmar 2023 Table 2 prints mouse Q_Lung = 0.618 L/h (0.17% higher).
    gfr      <- 0.0168   # Parmar 2023 Table 2, mouse GFR (L/h)

    v_venous   <- 0.00120    # Parmar 2023 Table 2, mouse V_Venous blood (L)
    v_arterial <- 0.000515   # Parmar 2023 Table 2, mouse V_Arterial blood (L)
    v_lung     <- 0.000194   # Parmar 2023 Table 2, mouse V_Lung (L)
    v_spleen   <- 0.000127   # Parmar 2023 Table 2, mouse V_Spleen (L)
    v_liver    <- 0.00193    # Parmar 2023 Table 2, mouse V_Liver (L)
    v_kidney   <- 0.000525   # Parmar 2023 Table 2, mouse V_Kidney (L)
    v_other    <- 0.0235     # Parmar 2023 Table 2, mouse V_Other (L)

    # Drug-specific parameters for spectinamide 1810 in the mouse (Parmar
    # 2023 Table 2, "Generated in this study as described under Methods").
    kbp <- 0.604  # Parmar 2023 Table 2, mouse spectinamide 1810 k(b/p)
    fu  <- 0.693  # Parmar 2023 Table 2, mouse spectinamide 1810 fu Plasma

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
    ka         <- exp(lka)

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
    #    Subcutaneous absorption enters venous blood (supplement S19-S21).
    # -----------------------------------------------------------------------
    d/dt(venous) <- (q_liver + q_spleen) * Cb_liver +
      q_kidney * Cb_kidney +
      q_other * Cb_other -
      q_lung * Cven +
      ka * depot
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
    #    the portal spleen outflow; returns (q_liver + q_spleen) x Cb_liver
    #    to venous blood.
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
    # 11. Subcutaneous absorption depot (supplement S19-S21). The paper's
    #     50 uL depot volume only converts the depot amount to a depot
    #     concentration in the authors' MlxTran code and cancels out of the
    #     first-order transfer, so it is not needed in the amount-based form.
    # -----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    f(depot) <- exp(lfdepot)

    # -----------------------------------------------------------------------
    # 12. Observations. Cc is the venous PLASMA concentration
    #     (supplement S2: C_VenousPlasma = C_VenousBlood / k(b/p)); the
    #     tissue concentrations are supplement S5, S9, S12 and S15. Kidney
    #     has no reported residual so Ckidney carries no error model.
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
