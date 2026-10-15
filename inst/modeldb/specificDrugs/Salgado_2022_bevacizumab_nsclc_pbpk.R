Salgado_2022_bevacizumab_nsclc_pbpk <- function() {
  description <- paste(
    "PBPK (minimal, tumor-to-lymph-node axis). Proof-of-concept",
    "physiologically based pharmacokinetic model of the anti-VEGF-A antibody bevacizumab distribution",
    "from plasma into an organ-specific primary tumor (non-small-cell lung cancer (lung)) and its",
    "tumor-draining lymph node (TDLN), calibrated to digitized human",
    "immuno-PET (89Zr-bevacizumab) standardized-uptake-value data.",
    "Plasma is a single well-stirred pool with linear systemic clearance;",
    "antibody enters the tumor interstitium by vascular convection, drains",
    "by afferent lymph to the TDLN with an empirical, distance-dependent",
    "FcRn-salvage fraction, and returns to plasma by efferent lymph.",
    "Target binding in the tumor and the TDLN follows quasi-equilibrium",
    "TMDD with target turnover and complex internalization.",
    "One of eleven organ- and antibody-specific calibrations in the source",
    "paper; deterministic (no between-subject or residual variability).",
    sep = " "
  )
  reference <- paste(
    "Salgado E, Cao Y. A Physiologically Based Pharmacokinetic Framework",
    "for Quantifying Antibody Distribution Gradients from Tumors to",
    "Tumor-Draining Lymph Nodes. Antibodies (Basel). 2022;11(2):28.",
    "doi:10.3390/antib11020028.",
    "Calibration data digitized by the authors from Bahce I et al., EJNMMI Res 2014;4:35.",
    sep = " "
  )
  vignette <- "Salgado_2022_antibody_tumor_lymph_node_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  covariateData <- list()

  # State amounts are nmol for the three antibody pools (Appendix A,
  # Eqs. A1, A2 and A4 are written as mass balances, volume x dC/dt);
  # the two target states are total-target concentrations in nM
  # (Eqs. A3 and A5).
  compartmentData <- list(
    central = list(analyte = "bevacizumab", units = "nmol", specimen = "plasma", verified = TRUE),
    is_tumor = list(
      analyte = "bevacizumab (total: free + target-bound)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    total_target_tumor = list(
      analyte = "VEGF-A (total: free + antibody-bound)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    lnode = list(
      analyte = "bevacizumab (total: free + target-bound)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    total_target_lnode = list(
      analyte = "VEGF-A (total: free + antibody-bound)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # Total-target concentration states in the tumor and the TDLN (Salgado
  # 2022 Eqs. A3 and A5, R_T,PT and R_T,TDLN); the registered location
  # suffixes for target states do not include a tumor or a lymph node.
  paper_specific_compartments <- c("total_target_tumor", "total_target_lnode")

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Advanced non-small-cell lung cancer",
    dose_range = paste(
      "Single intravenous tracer dose of 89Zr-bevacizumab in the source",
      "immuno-PET study; the dose is not restated in Salgado 2022."
    ),
    regions = NA_character_,
    notes = paste(
      "Calibrated to mean (+/- SD) standardized uptake values (SUV) in",
      "plasma, primary tumor and TDLN digitized with WebPlotDigitizer from",
      "Bahce I et al., EJNMMI Res 2014;4:35 (Salgado 2022 Section 2.2 and Figure 2b).",
      "Individual-level data were not used, so the number of subjects is",
      "not stated in Salgado 2022. One of eleven organ-specific",
      "calibrations drawn from eight immuno-PET studies."
    )
  )

  ini({
    # Calibrated organ- and antibody-specific parameters. Salgado 2022
    # Results 3.1: 'Tables 1 and 2 denote the optimized range of parameters
    # for calibrating each organ-specific tumor model ... (sigma_V, R01,
    # R02, Kd, dln, CLp)'. No standard errors are reported.
    lcl <- log(0.07); label("Systemic antibody clearance from plasma CLp (L/h)") # Table 2 row 'NSCLC [17]' CLp = 0.07 L/h
    sigma_v <- 0.85; label("Tumor vascular reflection coefficient sigma_V (unitless)") # Table 1 row 'NSCLC [17]' sigma_V = 0.85
    r0_tumor <- 10; label("Baseline target concentration in the primary-tumor interstitium R01 (nM)") # Table 2 row 'NSCLC [17]' R01 = 10 nM
    r0_lnode <- 10; label("Baseline target concentration in the TDLN interstitium R02 (nM)") # Table 2 row 'NSCLC [17]' R02 = 10 nM
    dln <- 21; label("Shape factor for the tumor-to-TDLN lymphatic distance in the FcRn-salvage fraction (unitless)") # Table 2 row 'NSCLC [17]' dln = 21; derivation in Appendix B

    # Antibody-specific literature parameters (Table 2 footnote references)
    kd <- fixed(0.058); label("Equilibrium dissociation constant of the antibody for its target Kd (nM)") # Table 2 row 'NSCLC [17]' Kd = 0.058 nM (literature, refs 30-32)
    kd_fcrn <- fixed(2400); label("Equilibrium dissociation constant of the antibody for FcRn Kd,FcRn (nM)") # Table 2 row 'NSCLC [17]' Kd,FcRn = 2400 nM (literature, ref 33)

    # Organ-specific physiological parameters (Table 1)
    lvc <- fixed(log(5.0)); label("Plasma volume Vp (L)") # Table 1 row 'NSCLC [17]' Vp = 5.0 L
    l_organ <- fixed(0.012); label("Organ-specific lymph flow Lorgan, 0.2% of organ blood flow (L/h)") # Table 1 row 'NSCLC [17]' Lorgan = 0.012 L/h (footnote a)
    v_isf_tumor <- fixed(0.175); label("Primary-tumor interstitial fluid volume VISF,PT, 20% of organ volume (L)") # Table 1 row 'NSCLC [17]' VISF,PT = 0.175 L (footnote b)

    # Parameters shared by every organ system (Table 1 and Table 2 footnotes)
    sigma_l <- fixed(0.2); label("Lymphatic reflection coefficient sigma_L (unitless)") # Table 1 footnote: sigma_L = 0.2 (ref 27)
    l_aff <- fixed(0.004); label("Afferent lymph flow from the tumor to the TDLN Laff (L/h)") # Table 1 footnote: Laff = 0.004 L/h (ref 28)
    l_eff <- fixed(0.004); label("Efferent lymph flow from the TDLN to plasma Leff (L/h)") # Table 1 footnote: Leff = 0.004 L/h (ref 28)
    lv_lnode <- fixed(log(0.0000584)); label("TDLN interstitial fluid volume VISF,TDLN, 20% of an average lymph node (L)") # Table 1 footnote: VISF,TDLN = 0.0000584 L (ref 28)
    fcrn <- fixed(40000); label("FcRn concentration in lymphatic endothelial cells (nM; 40 uM)") # Table 1 footnote: [FcRn] = 40 uM (ref 29)
    f_fcrn_base <- fixed(0.418); label("Fraction of antibody trafficked through the lymphatics independently of FcRn (unitless)") # Appendix B: baseline 0.418 added to Eq. A15 (ref 47, FcRn-knockout mice)
    kdeg <- fixed(0.01); label("Target degradation rate constant kdeg (1/h)") # Table 2 footnote: kdeg = 0.01 1/h
    kint <- fixed(0.01); label("Antibody-target complex internalization rate constant kint (1/h)") # Table 2 footnote: kint = 0.01 1/h

    # Residual error: SUV calibration only, no residual-error model reported
    propSd <- fixed(0); label("Proportional residual error, plasma (fraction; not reported)") # Not reported in Salgado 2022; deterministic model
    propSd_Ctumor <- fixed(0); label("Proportional residual error, primary tumor (fraction; not reported)") # Not reported in Salgado 2022; deterministic model
    propSd_Clnode <- fixed(0); label("Proportional residual error, TDLN (fraction; not reported)") # Not reported in Salgado 2022; deterministic model
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    v_lnode <- exp(lv_lnode)

    # Target synthesis rate. Appendix A Eqs. A3 and A5 print R01 (and R02)
    # as the zero-order synthesis term with initial condition R01; with
    # kdeg = 0.01 1/h that literal form is not at steady state and drives
    # total target up about 100-fold over 200 h. The target baseline is
    # held at steady state (ksyn = kdeg * R0, the synthesis rate 'ksyn'
    # that the Figure 1 caption names), which reproduces the Figure 2
    # calibration curves; see the vignette.
    ksyn_tumor <- kdeg * r0_tumor
    ksyn_lnode <- kdeg * r0_lnode

    # Total antibody concentrations (nM)
    Cc <- central / vc
    Ctumor <- is_tumor / v_isf_tumor
    Clnode <- lnode / v_lnode

    # Free antibody concentrations: quasi-equilibrium root, Eqs. A6 (tumor)
    # and A7 (TDLN): Cf = 0.5 * ((C - RT - Kd) + sqrt((C - RT - Kd)^2 +
    # 4 * C * Kd)). When C - RT - Kd < 0 the printed form subtracts two
    # nearly equal numbers; the algebraically identical rationalized form
    # 2 * C * Kd / (sqrt(...) - (C - RT - Kd)) is used there instead.
    b_tumor <- Ctumor - total_target_tumor - kd
    disc_tumor <- sqrt(b_tumor^2 + 4 * Ctumor * kd)
    if (b_tumor >= 0) {
      cf_tumor <- 0.5 * (b_tumor + disc_tumor)
    } else {
      cf_tumor <- 2 * Ctumor * kd / (disc_tumor - b_tumor)
    }
    b_lnode <- Clnode - total_target_lnode - kd
    disc_lnode <- sqrt(b_lnode^2 + 4 * Clnode * kd)
    if (b_lnode >= 0) {
      cf_lnode <- 0.5 * (b_lnode + disc_lnode)
    } else {
      cf_lnode <- 2 * Clnode * kd / (disc_lnode - b_lnode)
    }

    # Antibody-target complex concentrations (nM), Eqs. A8 and A9
    complex_tumor <- total_target_tumor * cf_tumor / (kd + cf_tumor)
    complex_lnode <- total_target_lnode * cf_lnode / (kd + cf_lnode)

    # Fraction of antibody surviving lymphatic transit through FcRn salvage,
    # Eq. A10 (FcRnb1 leaving the tumor, FcRnb2 leaving the TDLN)
    fcrnb_tumor <- f_fcrn_base + (fcrn / (cf_tumor + kd_fcrn + fcrn))^dln
    fcrnb_lnode <- f_fcrn_base + (fcrn / (cf_lnode + kd_fcrn + fcrn))^dln

    # Eq. A1: plasma (single tumor / TDLN channel; the sum over channels
    # applies only when metastases are added)
    d/dt(central) <- -(1 - sigma_v) * l_organ * Cc - cl * Cc +
      l_eff * cf_lnode * fcrnb_lnode
    # Eq. A2: total antibody in the primary-tumor interstitium
    d/dt(is_tumor) <- (1 - sigma_v) * l_organ * Cc -
      (1 - sigma_l) * l_aff * cf_tumor * fcrnb_tumor -
      kint * complex_tumor * v_isf_tumor
    # Eq. A3: total target in the primary tumor (IC = R01)
    d/dt(total_target_tumor) <- ksyn_tumor -
      kdeg * (total_target_tumor - complex_tumor) - kint * complex_tumor
    total_target_tumor(0) <- r0_tumor
    # Eq. A4: total antibody in the TDLN interstitium
    d/dt(lnode) <- (1 - sigma_l) * l_aff * cf_tumor * fcrnb_tumor -
      l_eff * cf_lnode * fcrnb_lnode - kint * complex_lnode * v_lnode
    # Eq. A5: total target in the TDLN (IC = R02)
    d/dt(total_target_lnode) <- ksyn_lnode -
      kdeg * (total_target_lnode - complex_lnode) - kint * complex_lnode
    total_target_lnode(0) <- r0_lnode

    Cc ~ prop(propSd)
    Ctumor ~ prop(propSd_Ctumor)
    Clnode ~ prop(propSd_Clnode)
  })
}
