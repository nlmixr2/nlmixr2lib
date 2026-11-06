Pressly_2023_trastuzumab_mpbpk <- function() {
  description <- paste(
    "PBPK (minimal, TMDD). Model 2 of a four-model tutorial hierarchy for",
    "trastuzumab: the basic TMDD model with tissue distribution resolved into",
    "a vascular space, an endothelial space, the interstitial space and the",
    "lymph, after Li et al. (2014). Antibody reaches the vascular space with",
    "the cardiac output, is taken up into the endothelium (where it is",
    "cleared) from the vascular and interstitial spaces, crosses to the",
    "interstitium by convection with the lymph flow (vascular reflection",
    "coefficient), and returns to plasma through the lymph (interstitial",
    "reflection coefficient). Target binding, target turnover and complex",
    "internalisation are confined to plasma. No FcRn recycling. Parameter",
    "values are the authors' trastuzumab set, mostly from Li et al. (2014),",
    "not fitted. Deterministic: no between-subject or residual variability.",
    "Siblings: modellib('Pressly_2023_trastuzumab_tmdd'),",
    "modellib('Pressly_2023_trastuzumab_mpbpk_fcrn') and",
    "modellib('Pressly_2023_trastuzumab_mpbpk_shedding').",
    sep = " "
  )
  reference <- paste(
    "Pressly MA, Peletier LA, Zheng S, Sharma VD, Lien YTK, Wang W, Zhou H,",
    "Schmidt S. The quest for balance between capturing data and model",
    "complexity: A quantitative clinical pharmacology approach applied to",
    "monoclonal antibodies. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(5):639-652. doi:10.1002/psp4.12927.",
    "Equations 7-12 (main text); parameter values Supplementary Material S2",
    "Table S2; deposited MATLAB code Data S1 (tmdd2.m, tmdd2_data_impulse.m).",
    sep = " "
  )
  vignette <- "Pressly_2023_trastuzumab_tmdd_tutorial"

  # The ODEs are written in concentrations (nM) in the source. Here the five
  # antibody states hold amounts (nmol) so that a dose is given as an amount;
  # the target and complex states keep the source's plasma concentrations (nM).
  units <- list(time = "h", dosing = "nmol", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    tissue_vasc = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    endosome = list(analyte = "trastuzumab", units = "nmol", specimen = "endosome", verified = TRUE),
    tissue_isf = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "lymph", verified = TRUE),
    target = list(analyte = "target (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex = list(analyte = "trastuzumab-target complex", units = "nM", specimen = "plasma", verified = TRUE)
  )

  dosing <- "central"

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "HER2-positive metastatic breast cancer (comparison data); the model itself is a parameterised tutorial, not a fit",
    dose_range = "1, 2, 4 and 8 mg/kg single 90-min i.v. infusion (Tokuda 1999 comparison data); 650 nM i.v. bolus in the didactic simulations",
    weight_median = "70 kg (reference body weight used by the authors to convert mg/kg doses)",
    notes = paste(
      "TUTORIAL model. Parameter values are the authors' trastuzumab set",
      "(Supplementary Table S2, mostly from Li et al. 2014 CPT PSP), shared",
      "by all four models; nothing was estimated. The authors overlaid the",
      "simulations on mean concentrations digitised from Tokuda et al.",
      "1999 (18 Japanese patients with HER2-overexpressing metastatic",
      "breast cancer, single 1-8 mg/kg 90-min i.v. infusions). Figure 9e",
      "shows an adjusted variant with kup = 0.00536 1/h. There is no",
      "inter-individual variability and no residual error; propSd is fixed(0)."
    )
  )

  ini({
    lvc <- fixed(log(3)); label("Plasma volume Vp (L)") # Table S2 'Vp' = 3 L
    lv_tissue_vasc <- fixed(log(1.61)); label("Vascular-space volume Vv (L)") # Table S2 'Vv' = 1.61 L
    lv_endosome <- fixed(log(0.08)); label("Endothelial-space volume Ve (L)") # Table S2 'Ve' = 0.08 L
    lv_tissue_isf <- fixed(log(11.2)); label("Interstitial-space volume Vi (L)") # Table S2 'Vi' = 11.2 L
    lv_lymph <- fixed(log(0.2702)); label("Lymph volume Vs (L)") # Table S2 'Vs' = 0.2702 L
    lqc <- fixed(log(300)); label("Cardiac output Q (L/h)") # Table S2 'Q' = 300 L/h
    lq_lymph <- fixed(log(0.119)); label("Lymph flow S (L/h)") # Table S2 'S' = 0.119 L/h
    sigma_v <- fixed(0.9); label("Vascular reflection coefficient sigma_v (unitless)") # Table S2 'sigma_v' = 0.9
    sigma_i <- fixed(0.2); label("Interstitial (lymphatic) reflection coefficient sigma_i (unitless)") # Table S2 'sigma_i' = 0.2
    lkup <- fixed(log(0.0298)); label("Endothelial uptake rate constant from the vascular space kup (1/h)") # Table S2 'kup' = 0.0298 1/h
    delta_kup <- fixed(1); label("Multiplier on kup for uptake from the interstitial space delta (unitless)") # Table S2 'delta' = 1
    lclint <- fixed(log(0.01751)); label("Clearance of antibody from the endothelial space CLint (L/h)") # Table S2 'CLint' = 0.01751 L/h
    kel <- fixed(0); label("First-order elimination rate constant of free antibody from plasma kel,p (1/h)") # Table S2 'kel,p' = 0; tmdd2.m 'k10=0'

    lkon <- fixed(log(2.52)); label("Association rate constant of antibody and target kon,p (1/(nM*h))") # Table S2 'kon,p' = 2.52
    lkoff <- fixed(log(1.26)); label("Dissociation rate constant of the antibody-target complex koff,p (1/h)") # Table S2 'koff,p' = 1.26
    lksyn <- fixed(log(0.83)); label("Zero-order synthesis rate of plasma target ksyn,p (nM/h)") # Table S2 'ksyn,p' = 0.83
    lkdeg <- fixed(log(0.1)); label("First-order degradation rate constant of free target kdeg,p (1/h)") # Table S2 'kdeg,p' = 0.1
    lkint <- fixed(log(0.01)); label("Internalisation rate constant of the antibody-target complex kint,p (1/h)") # Table S2 'kint,p' = 0.01

    mw <- fixed(145532); label("Molecular weight of trastuzumab (g/mol)") # Table S2 'MW' = 145,532 Da

    # The source is a deterministic simulation with no residual error;
    # fixed(0) keeps the model simulable without inventing a value.
    propSd <- fixed(0); label("Proportional residual error (fraction; zero - none in the source)")
  })

  model({
    vc <- exp(lvc)
    v_tissue_vasc <- exp(lv_tissue_vasc)
    v_endosome <- exp(lv_endosome)
    v_tissue_isf <- exp(lv_tissue_isf)
    v_lymph <- exp(lv_lymph)
    qc <- exp(lqc)
    q_lymph <- exp(lq_lymph)
    kup <- exp(lkup)
    clint <- exp(lclint)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)

    # Equations 5 and 12: target at its antibody-free baseline at time zero.
    target(0) <- ksyn / kdeg

    lig_p <- central / vc # Lp (nM)
    lig_v <- tissue_vasc / v_tissue_vasc # Lv (nM)
    lig_e <- endosome / v_endosome # Le (nM)
    lig_i <- tissue_isf / v_tissue_isf # Li (nM)
    lig_s <- lymph / v_lymph # Ls (nM)

    # Equations 7-11, each antibody equation multiplied through by its volume.
    d/dt(central) <- (qc - q_lymph) * lig_v + q_lymph * lig_s - qc * lig_p +
      vc * (-kon * lig_p * target + koff * complex - kel * lig_p)
    d/dt(target) <- ksyn - kdeg * target - kon * lig_p * target + koff * complex
    d/dt(complex) <- kon * lig_p * target - (koff + kint) * complex
    d/dt(tissue_vasc) <- qc * lig_p - (qc - q_lymph) * lig_v - (1 - sigma_v) * q_lymph * lig_v -
      kup * tissue_vasc
    d/dt(endosome) <- kup * tissue_vasc + delta_kup * kup * tissue_isf - clint * lig_e
    d/dt(tissue_isf) <- -delta_kup * kup * tissue_isf + (1 - sigma_v) * q_lymph * lig_v -
      (1 - sigma_i) * q_lymph * lig_i
    d/dt(lymph) <- (1 - sigma_i) * q_lymph * lig_i - q_lymph * lig_s

    # Free trastuzumab in plasma: nM * g/mol = ng/L, / 1e6 = ug/mL.
    Cc <- lig_p * mw / 1e6
    Cc ~ prop(propSd)
  })
}
