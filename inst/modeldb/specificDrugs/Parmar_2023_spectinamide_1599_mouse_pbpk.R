Parmar_2023_spectinamide_1599_mouse_pbpk <- function() {
  description <- paste(
    "Minimal PBPK (mPBPK; author-coded ODEs, fitted in Monolix 2021R1).",
    "Preclinical (BALB/c mouse, 20 g). Spectinamide 1599, an anti-tuberculosis",
    "spectinomycin analogue, after intravenous, subcutaneous, and",
    "intrapulmonary-aerosol dosing (Parmar 2023, Pharmaceutics). Venous and",
    "arterial blood plus five tissues (lung, spleen, liver, kidney, and a",
    "lumped 'other'), each tissue split into a rapid-equilibrium",
    "extracellular pool (vascular + interstitial, in instantaneous",
    "equilibrium with blood through the blood-to-plasma ratio) and a slow",
    "cellular pool coupled by a first-order influx K(I->C) acting on the",
    "unbound fraction and a first-order back flux K(C->I). The lung sits in",
    "series between venous and arterial blood; the spleen drains portally",
    "into the liver; elimination is glomerular filtration (GFR x fu) acting",
    "on the kidney extracellular pool. Subcutaneous doses reach venous blood",
    "through depot with ka / F; intrapulmonary-aerosol doses reach an",
    "epithelial-lining-fluid compartment (elf) through depot2 with its own",
    "ka / F, and the ELF exchanges with the lung cellular and extracellular",
    "pools. Tissue concentrations (lung, spleen, liver, kidney) are",
    "simultaneous outputs alongside venous plasma Cc. Two mass-balance",
    "defects in the published equations are corrected here (venous outflow",
    "to the lung, and hepatic outflow of the portal spleen flow); see the",
    "vignette Errata."
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
    dosing        = "mg (absolute amount; the paper's mg/kg dose x 0.02 kg mouse body weight, e.g. 10 mg/kg = 0.2 mg)",
    concentration = "mg/L (equivalently ug/mL for plasma and ug/g for tissue under the paper's unit-tissue-density assumption)"
  )

  # No covariates. Parmar 2023 does not carry any subject-level covariate on
  # the model parameters: the physiology is tabulated per species (Table 2 /
  # Table 3) and held fixed, and no demographic, weight, or disease-status
  # effect is estimated (the paper explicitly concludes there is no relevant
  # PK difference between healthy and Mtb-infected mice, Section 3.6). The
  # Sprague-Dawley rat physiology is a separate model file
  # (Parmar_2023_spectinamide_1599_rat_pbpk.R) rather than a covariate switch,
  # following the An_2012_mitoxantrone_{mouse,human}_pbpk precedent.
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    venous               = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "blood cell", verified = FALSE),
    arterial             = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "blood cell", verified = FALSE),
    lung_extracellular   = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    lung_cellular        = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    elf                  = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "blood cell", verified = FALSE),
    spleen_extracellular = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    spleen_cellular      = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    liver_extracellular  = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    liver_cellular       = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    kidney_extracellular = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    kidney_cellular      = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    other_extracellular  = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    other_cellular       = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "tissue", verified = FALSE),
    depot                = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "administration site", verified = FALSE),
    depot2               = list(analyte = "Spectinamide 1599", units = NA_character_, specimen = "administration site", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (BALB/c, 20 g)",
    n_subjects     = 552L,
    n_studies      = 4L,
    age_range      = NA_character_,
    weight_range   = "20 g reference body weight (Parmar 2023 Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Healthy BALB/c mice for all model-parameterisation and",
      "model-qualification data (plasma + lung + liver + spleen + kidney +",
      "ELF). Mycobacterium tuberculosis-infected BALB/c mice (studies 1A /",
      "1B / 1C, plasma only) were used as an external predictive check of",
      "the subcutaneous model; Parmar 2023 Section 3.6 found no relevant PK",
      "difference between healthy and infected animals, so no",
      "disease-status effect is encoded."
    ),
    dose_range     = paste(
      "Intravenous 10 mg/kg single dose and QD5; subcutaneous 50 and",
      "200 mg/kg single dose, QD5 / BIW / TIW; intrapulmonary aerosol 10,",
      "50 and 150 mg/kg single dose, QD5 / BIW / TIW (Parmar 2023 Table 1).",
      "Infected-mouse subcutaneous regimens spanned 1-200 mg/kg per dose at",
      "BID / QD / TIW / BIW / QW for 4 weeks."
    ),
    regions        = "USA (University of Tennessee Health Science Center; Colorado State University BSL-3 for infected animals)",
    notes          = paste(
      "n_subjects counts the BALB/c mice in Parmar 2023 Table 1 that",
      "contributed spectinamide 1599 data: 168 healthy mice for the IV and",
      "SC datasets, 234 healthy mice for the intrapulmonary-aerosol",
      "datasets, and 150 Mtb-infected mice for the external subcutaneous",
      "check. Destructive sampling: 1 data point per animal (plasma + lung +",
      "liver + spleen, plus ELF in the aerosol studies), 3 animals per",
      "sampling time point, except the infected studies with 2 plasma points",
      "per animal and 5-6 animals per time point. n_studies counts the four",
      "healthy-mouse route groups plus infected studies 1A-1C as one",
      "external-validation group. Parameter estimates are the mouse + rat",
      "joint fit reported in Table 4 (the inter-animal variances were",
      "'estimated by simultaneously fitting the mouse and rat plasma and",
      "tissue data obtained after intravenous administration')."
    )
  )

  ini({
    # =========================================================================
    # Tissue distribution rate constants (Parmar 2023 Table 4, Intravenous
    # column; verified against the typeset PDF pages 9-10 because the
    # machine-extracted markdown scrambles Table 4). All ten constants are
    # estimated from the healthy-mouse + rat single-dose intravenous plasma
    # and tissue data and are then held fixed for the intratracheal and
    # subcutaneous fits ("Fixed" in the Intratracheal / Subcutaneous columns),
    # so they are the same parameters here rather than route-specific ones.
    #
    # K(I->C) is the first-order uptake from the rapid-equilibrium
    # (vascular + interstitial) sub-compartment into the cellular
    # sub-compartment and always multiplies the unbound fraction fu;
    # K(C->I) is the first-order back flux out of the cellular pool and acts
    # on the total cellular amount.
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
    # Extravascular absorption. Two independent routes with their own
    # first-order absorption rate constant and bioavailability component:
    # subcutaneous (depot -> venous blood) and intrapulmonary aerosol
    # (depot2 -> epithelial lining fluid). Parmar 2023 Section 2.8: "For
    # subcutaneous administration, the dose was administered as an external
    # compartment, absorbed into the venous compartment via a first-order
    # absorption rate constant and a bioavailability component. ... The dose
    # was administered as an external compartment, absorbed into the ELF
    # compartment via a first-order absorption rate constant and a
    # bioavailability component."
    # =========================================================================
    lka <- log(4.36)
    label("Subcutaneous first-order absorption rate constant Ka (1/h)")
    # Parmar 2023 Table 4 Subcutaneous: Ka = 4.36 (RSE 7.86%).

    lfdepot <- log(0.86)
    label("Subcutaneous bioavailability F (fraction)")
    # Parmar 2023 Table 4 Subcutaneous: F = 0.86 (RSE 6.15%).

    lka2 <- log(5.03)
    label("Intrapulmonary-aerosol first-order absorption rate constant Ka (1/h)")
    # Parmar 2023 Table 4 Intratracheal: Ka = 5.03 (RSE 4.53%).

    lfdepot2 <- log(0.33)
    label("Intrapulmonary-aerosol bioavailability F (fraction)")
    # Parmar 2023 Table 4 Intratracheal: F = 0.33 (RSE 2.03%).

    # =========================================================================
    # Inter-animal variability on the four tissue uptake rate constants
    # (Parmar 2023 Table 4, Intravenous column). The table's omega symbols
    # are STANDARD DEVIATIONS on the log scale, not variances, despite the
    # footnote calling omega "the inter-animal variance": the paper's own
    # Section 3.7 states the values 0.48-0.87 "correspond to 50.9 to 106% CV",
    # and sqrt(exp(0.48^2) - 1) = 50.9% while sqrt(exp(0.87^2) - 1) = 106%
    # (the variance reading would give 78.5-106% instead). nlmixr2 ini()
    # takes variances, so each entry below is the published SD squared.
    # No inter-animal variability was estimated on the kidney uptake
    # constant, on any back-flux constant, or on Ka / F.
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
    # Residual unexplained variability. Parmar 2023 Table 4 reports a
    # proportional residual per matrix and per route. The bare parameters
    # below carry the INTRAVENOUS column, which is the dataset that
    # parameterised the structural model; the subcutaneous and intratracheal
    # columns are recorded in the comments and in the vignette Errata so a
    # user can override them when simulating those routes. No residual error
    # is reported for kidney, so Ckidney is exposed as a derived
    # concentration without a residual-error endpoint.
    # =========================================================================
    propSd <- 0.32
    label("Proportional residual SD for venous plasma Cc (fraction)")
    # Parmar 2023 Table 4: eps_Plasma IV 0.32 (RSE 15.1%);
    # subcutaneous 0.54 (17.6%); intratracheal 6.57 (8.62%).

    propSd_Clung <- 0.35
    label("Proportional residual SD for lung concentration (fraction)")
    # Parmar 2023 Table 4: eps_Lung IV 0.35 (RSE 14.4%);
    # subcutaneous 0.50 (15.3%); intratracheal 0.95 (9.92%).

    propSd_Cliver <- 0.28
    label("Proportional residual SD for liver concentration (fraction)")
    # Parmar 2023 Table 4: eps_Liver IV 0.28 (RSE 14.4%);
    # subcutaneous 0.29 (19.1%); intratracheal 2.55 (10.3%).

    propSd_Cspleen <- 0.53
    label("Proportional residual SD for spleen concentration (fraction)")
    # Parmar 2023 Table 4: eps_Spleen IV 0.53 (RSE 14.4%);
    # subcutaneous 0.97 (16.5%); intratracheal 8.48 (9.86%).
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Mouse physiology, held fixed (Parmar 2023 Table 2, Mouse (20 g)
    #    column; blood flows and GFR from ref [16], volumes from ref [17],
    #    V_ELF from ref [18]). Flows in L/h, volumes in L.
    #
    #    The tabulated flows satisfy Q_Liver + Q_Spleen + Q_Kidney + Q_Other
    #    = 0.139 + 0.00695 + 0.100 + 0.371 = 0.61695 L/h against a tabulated
    #    Q_Lung = 0.618 L/h, i.e. Q_Liver is the hepatic-ARTERIAL flow and
    #    the spleen is a parallel arterial branch that drains portally into
    #    the liver. That identity is what fixes the hepatic-outflow
    #    correction in block 5 below.
    #
    #    q_lung is DERIVED as the sum rather than taken as the printed
    #    0.618 L/h. The pulmonary flow is the cardiac output and must equal
    #    the sum of the systemic organ flows exactly, otherwise the arterial
    #    node removes q_lung x Cart while the tissues receive only
    #    (q_liver + q_spleen + q_kidney + q_other) x Cart and the difference
    #    is destroyed. With the printed 0.618 that leak is 0.17% of cardiac
    #    output per pass, which is ~10% of total drug removal here because
    #    glomerular filtration is itself a small clearance -- 2.8% of the
    #    dose goes missing by 72 h. Using the sum (a 0.17% change, i.e.
    #    rounding-level) closes the loop to solver tolerance. See the
    #    vignette Errata and mass-balance check.
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
    v_elf      <- 0.0000100  # Parmar 2023 Table 2, mouse V_ELF (L)

    # Drug-specific parameters for spectinamide 1599 in the mouse
    # (Parmar 2023 Table 2, "Generated in this study as described under
    # Methods"). k(b/p) converts blood to plasma concentrations (S2);
    # fu is the plasma unbound fraction that gates cellular uptake and
    # glomerular filtration; fu_elf is the unbound fraction in epithelial
    # lining fluid, computed by the authors from fu and the plasma/ELF
    # albumin ratio via Poulin and Theil.
    kbp    <- 0.552  # Parmar 2023 Table 2, mouse spectinamide 1599 k(b/p)
    fu     <- 0.602  # Parmar 2023 Table 2, mouse spectinamide 1599 fu Plasma
    fu_elf <- 0.948  # Parmar 2023 Table 2, mouse spectinamide 1599 fu ELF

    # -----------------------------------------------------------------------
    # 2. Tissue sub-compartment volume fractions (Parmar 2023 Table 3, from
    #    ref [17]; identical in mouse and rat).
    #
    #    The rapid-equilibrium ("extracellular") pool lumps the vascular and
    #    interstitial spaces, which the paper assumes to be in instantaneous
    #    equilibrium with blood. Its effective volume is
    #    V_vascular + V_interstitial / k(b/p) so that a single blood-scale
    #    concentration Cb describes both spaces (supplement S3-S18 LHS).
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
    ka         <- exp(lka)
    ka2        <- exp(lka2)

    # ELF transfer constants for the intrapulmonary-aerosol route
    # (supplement S23-S25). Their numeric values are never tabulated, but
    # Figure 2 labels every ELF arrow with an already-reported constant:
    #   * the ELF <-> cellular arrows carry "K_B2C x fup_ELF" and "K3",
    #     the same K_B2C symbol used for the interstitium -> cellular arrow,
    #     so K(ELF->C) is the lung uptake constant K_Lung(I->C);
    #   * the ELF <-> rapid-pool arrows carry "Ka x fup" and
    #     "Ka x fup_ELF", matching the fu pattern of the two K(ELF->B) terms
    #     in S23 / S25, so K(ELF->B) is the intrapulmonary Ka;
    #   * S23 writes the cellular -> ELF back flux with the lung's own
    #     K_Lung(C->I), so reusing the lung constants in the ELF equations is
    #     the authors' own pattern rather than an inference.
    # The identification is checked empirically in the vignette against the
    # authors' published intrapulmonary predicted AUCinf (Table 5).
    kelfc      <- kic_lung  # K(ELF->C) = K_Lung(I->C), Figure 2 "K_B2C"
    kelfb      <- ka2       # K(ELF->B) = intrapulmonary Ka, Figure 2 "Ka"
    kcelf_lung <- kci_lung  # cellular -> ELF back flux, supplement S23

    # -----------------------------------------------------------------------
    # 4. Concentrations. States are amounts (mg); Cb_<organ> is the
    #    blood-scale concentration of the rapid-equilibrium pool and
    #    Ccell_<organ> the cellular concentration.
    # -----------------------------------------------------------------------
    Cven <- venous   / v_venous
    Cart <- arterial / v_arterial
    Cb_lung   <- lung_extracellular   / veff_lung
    Cb_spleen <- spleen_extracellular / veff_spleen
    Cb_liver  <- liver_extracellular  / veff_liver
    Cb_kidney <- kidney_extracellular / veff_kidney
    Cb_other  <- other_extracellular  / veff_other

    # -----------------------------------------------------------------------
    # 5. Blood compartments.
    #
    #    Supplement S1 is printed as
    #      V_ven dC_ven/dt = Q_Liver Cb_Liver + Q_Kidney Cb_Kidney
    #                        + Q_Other Cb_Other - Q_Lung Cb_Lung
    #    Two corrections are applied (both documented in the vignette
    #    Errata):
    #      (a) the venous LOSS to the lung must be Q_Lung x C_VenousBlood,
    #          not Q_Lung x Cb_Lung, because S3 has the lung GAINING
    #          Q_Lung x C_VenousBlood; as printed the pair does not
    #          conserve drug instant by instant.
    #      (b) the hepatic venous return must carry the portal spleen flow,
    #          (Q_Liver + Q_Spleen) x Cb_Liver, because S10 adds
    #          Q_Spleen x Cb_Spleen to the liver while removing only
    #          Q_Liver x Cb_Liver.
    #    Both follow uniquely from mass balance and from the tabulated-flow
    #    identity in block 1; with them the system conserves drug exactly.
    #    Both corrections leave the AUC-level balance unchanged and alter only
    #    the transient shape: integrating S3 to infinity gives
    #    AUC(Cb_Lung) = AUC(C_VenousBlood), so (a) is AUC-neutral, and the
    #    printed liver equation inflates AUC(Cb_Liver) by exactly
    #    (Q_Liver + Q_Spleen) / Q_Liver, which cancels against the smaller
    #    outflow coefficient in (b). The rounding fix in block 1, by contrast,
    #    DOES change AUC -- see the vignette mass-balance check.
    #
    #    Subcutaneous absorption enters venous blood (Parmar 2023
    #    Section 2.8); supplement S19-S21.
    # -----------------------------------------------------------------------
    d/dt(venous) <- (q_liver + q_spleen) * Cb_liver +
      q_kidney * Cb_kidney +
      q_other * Cb_other -
      q_lung * Cven +
      ka * depot
    d/dt(arterial) <- q_lung * (Cb_lung - Cart)  # supplement S6

    # -----------------------------------------------------------------------
    # 6. Lung, in series between venous and arterial blood (supplement S3 /
    #    S4, extended to S25 / S24 for the intrapulmonary-aerosol route).
    #    Because K(I->C) multiplies Cb x V_effective it reduces to
    #    K(I->C) x fu x <organ>_extracellular in amount form, and
    #    K(C->I) x Ccell x V_cellular reduces to
    #    K(C->I) x <organ>_cellular.
    #
    #    The ELF terms are inert when no aerosol dose has been given
    #    (elf = 0), so the same equations serve the IV / SC and the aerosol
    #    routes.
    # -----------------------------------------------------------------------
    d/dt(lung_extracellular) <- q_lung * (Cven - Cb_lung) -
      kic_lung * fu * lung_extracellular +
      kci_lung * lung_cellular +
      kelfb * fu_elf * elf -
      kelfb * fu * lung_extracellular
    d/dt(lung_cellular) <- kic_lung * fu * lung_extracellular -
      kci_lung * lung_cellular -
      kcelf_lung * lung_cellular +
      kelfc * fu_elf * elf

    # Epithelial lining fluid (supplement S23). The printed S23 / S24 / S25
    # triple adds the cellular back flux K_Lung(C->I) x A_cellular to BOTH
    # the ELF pool (S23) and the lung rapid pool (S25) while removing it from
    # the cellular pool only once (S24): as printed, the aerosol sub-model
    # creates drug. Figure 2 draws TWO distinct arrows out of the lung
    # cellular pool ("K2" to the interstitium and "K3" to the ELF), so the
    # missing S24 loss term is restored here as kcelf_lung x lung_cellular,
    # with kcelf_lung identified with K_Lung(C->I) exactly as S23 writes it.
    # This restores mass balance without inventing a number; the alternative
    # reading (split the single printed efflux between the two destinations
    # with an invented fraction) is recorded in the vignette Errata.
    d/dt(elf) <- ka2 * depot2 +
      kcelf_lung * lung_cellular -
      kelfc * fu_elf * elf -
      kelfb * fu_elf * elf +
      kelfb * fu * lung_extracellular

    # -----------------------------------------------------------------------
    # 7. Spleen (supplement S7 / S8). Perfused from arterial blood; its
    #    venous outflow drains portally into the liver rather than into
    #    venous blood.
    # -----------------------------------------------------------------------
    d/dt(spleen_extracellular) <- q_spleen * (Cart - Cb_spleen) -
      kic_spleen * fu * spleen_extracellular +
      kci_spleen * spleen_cellular
    d/dt(spleen_cellular) <- kic_spleen * fu * spleen_extracellular -
      kci_spleen * spleen_cellular

    # -----------------------------------------------------------------------
    # 8. Liver (supplement S10 / S11). Receives hepatic-arterial flow
    #    q_liver from arterial blood plus the portal spleen outflow
    #    q_spleen x Cb_spleen, and returns (q_liver + q_spleen) x Cb_liver to
    #    venous blood (correction (b) of block 5).
    # -----------------------------------------------------------------------
    d/dt(liver_extracellular) <- q_liver * Cart +
      q_spleen * Cb_spleen -
      (q_liver + q_spleen) * Cb_liver -
      kic_liver * fu * liver_extracellular +
      kci_liver * liver_cellular
    d/dt(liver_cellular) <- kic_liver * fu * liver_extracellular -
      kci_liver * liver_cellular

    # -----------------------------------------------------------------------
    # 9. Kidney (supplement S13 / S14). The only elimination pathway in the
    #    model: glomerular filtration of the unbound drug out of the
    #    rapid-equilibrium pool. Parmar 2023 Section 2.8: "the unbound
    #    spectinamide 1599 and 1810 concentrations in plasma were assumed to
    #    be eliminated entirely by the kidney via glomerular filtration".
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
    # 11. Absorption depots. depot is the subcutaneous site (supplement
    #     S19-S21, absorbed into venous blood); depot2 is the intratracheal
    #     aerosol dose (supplement S22, absorbed into the ELF). The paper's
    #     50 uL subcutaneous depot volume (S19) only converts the depot
    #     amount to a depot concentration in the authors' MlxTran code and
    #     cancels out of the first-order transfer, so it is not needed in the
    #     amount-based form used here.
    # -----------------------------------------------------------------------
    d/dt(depot)  <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2

    f(depot)  <- exp(lfdepot)
    f(depot2) <- exp(lfdepot2)

    # -----------------------------------------------------------------------
    # 12. Observations.
    #     Cc      -- venous PLASMA concentration, supplement S2:
    #                C_VenousPlasma = C_VenousBlood / k(b/p).
    #     Clung   -- whole-lung concentration including the ELF pool,
    #                supplement S26 (reduces to S5 when elf = 0).
    #     Cspleen, Cliver, Ckidney -- supplement S9, S12, S15.
    #     Proportional residual error per matrix; kidney has no reported
    #     residual so Ckidney carries no error model.
    # -----------------------------------------------------------------------
    # Renal elimination rate (mg/h), exposed as a derived output so a
    # simulation can verify mass balance without re-hardcoding GFR and fu:
    # the amount remaining in all states plus the time integral of elimRate
    # must equal the absorbed dose. See the vignette mass-balance check.
    elimRate <- gfr * fu * Cb_kidney

    Cc      <- Cven / kbp
    Clung   <- (lung_extracellular + lung_cellular + elf) / v_lung
    Celf    <- elf / v_elf
    Cspleen <- (spleen_extracellular + spleen_cellular) / v_spleen
    Cliver  <- (liver_extracellular + liver_cellular) / v_liver
    Ckidney <- (kidney_extracellular + kidney_cellular) / v_kidney

    Cc      ~ prop(propSd)
    Clung   ~ prop(propSd_Clung)
    Cliver  ~ prop(propSd_Cliver)
    Cspleen ~ prop(propSd_Cspleen)
  })
}
