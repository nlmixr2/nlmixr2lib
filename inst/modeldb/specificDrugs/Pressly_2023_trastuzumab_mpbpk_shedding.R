Pressly_2023_trastuzumab_mpbpk_shedding <- function() {
  description <- paste(
    "PBPK (minimal, TMDD). Model 4 of a four-model tutorial hierarchy for",
    "trastuzumab: Model 3 (plasma, vascular, endothelial, interstitial and",
    "lymph spaces with quasi-steady-state FcRn recycling in the endothelium,",
    "after Li et al. 2014) with the target moved to the interstitium and",
    "target shedding added. Membrane-bound HER2 is synthesised and degraded",
    "in the interstitial space and sheds its ectodomain at a first-order",
    "rate; the shed (soluble) target is degraded in the interstitium or",
    "transferred to plasma, where it is degraded. Antibody binds membrane",
    "and shed target in the interstitium and shed target in plasma; each",
    "complex is internalised. Parameter values are the authors' trastuzumab",
    "set, mostly from Li et al. (2014), not fitted. Deterministic: no",
    "between-subject or residual variability.",
    "Siblings: modellib('Pressly_2023_trastuzumab_tmdd'),",
    "modellib('Pressly_2023_trastuzumab_mpbpk') and",
    "modellib('Pressly_2023_trastuzumab_mpbpk_fcrn').",
    sep = " "
  )
  reference <- paste(
    "Pressly MA, Peletier LA, Zheng S, Sharma VD, Lien YTK, Wang W, Zhou H,",
    "Schmidt S. The quest for balance between capturing data and model",
    "complexity: A quantitative clinical pharmacology approach applied to",
    "monoclonal antibodies. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(5):639-652. doi:10.1002/psp4.12927.",
    "Equations 23-36 (main text), Equation 18 and Supplementary Material S4;",
    "parameter values Supplementary Material S2 Table S2; deposited MATLAB code",
    "Data S1 (tmdd4.m, tmdd4_data_impulse.m).",
    sep = " "
  )
  vignette <- "Pressly_2023_trastuzumab_tmdd_tutorial"

  # The ODEs are written in concentrations (nM) in the source. Here the five
  # antibody states hold amounts (nmol) so that a dose is given as an amount;
  # the target and complex states keep the source's concentrations (nM) in
  # plasma or in the interstitial space.
  units <- list(time = "h", dosing = "nmol", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    tissue_vasc = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    endosome = list(
      analyte = "trastuzumab (total: free + FcRn-bound)",
      units = "nmol",
      specimen = "endosome",
      verified = TRUE
    ),
    tissue_isf = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "lymph", verified = TRUE),
    target = list(
      analyte = "shed HER2 extracellular domain (free)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    complex = list(analyte = "trastuzumab-shed HER2 complex", units = "nM", specimen = "plasma", verified = TRUE),
    target_isf = list(
      analyte = "shed HER2 extracellular domain (free)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_isf = list(analyte = "trastuzumab-shed HER2 complex", units = "nM", specimen = "tissue", verified = TRUE),
    target_mem_isf = list(analyte = "membrane-bound HER2 (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_mem_isf = list(
      analyte = "trastuzumab-membrane HER2 complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # Membrane-bound target and its antibody complex in the interstitial space
  # (Equation 28, R_m,i and LR_m,i). The registered location-suffixed target
  # family (`target_isf`, `complex_isf`) holds the shed, soluble target here;
  # the membrane-bound species has no canonical name.
  paper_specific_compartments <- c("target_mem_isf", "complex_mem_isf")

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
      "breast cancer, single 1-8 mg/kg 90-min i.v. infusions). Figure 9g",
      "shows an adjusted variant with krc = 3.13 1/h, lambda = 0.5 1/h and",
      "kup = 0.0507 1/h. There is no",
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
    lclint <- fixed(log(0.01751)); label("Clearance of FcRn-unbound antibody from the endothelial space CLint (L/h)") # Table S2 'CLint' = 0.01751 L/h
    kel <- fixed(0); label("First-order elimination rate constant of free antibody from plasma kel,p (1/h)") # Table S2 'kel,p' = 0; tmdd4.m 'k10=0'

    lc_fcrn_t <- fixed(log(40000)); label("Total FcRn concentration in the endothelial space R_FcRn (nM)") # Table S2 'FcRn' = 40000 nM
    lkon_fcrn <- fixed(log(1)); label("Association rate constant of antibody and FcRn kon,F (1/(nM*h))") # Table S2 'kon,F' = 1
    lkoff_fcrn <- fixed(log(774)); label("Dissociation rate constant of the antibody-FcRn complex koff,F (1/h)") # Table S2 'koff,F' = 774
    fr <- fixed(0.715); label("Fraction of recycled FcRn-bound antibody returned to the vascular space FR (unitless)") # Table S2 'FR' = 0.715 (Simcyp V21)
    lkrec <- fixed(log(16.7)); label("Recycling rate constant of FcRn-bound antibody krc (1/h)") # Table S2 'krc' = 16.7 1/h (Simcyp V21)

    # One association / dissociation pair serves all three binding reactions:
    # Table S2 gives KD,sh = KD,m = koff,p / kon,p and the text sets
    # koff,sh = koff,m and kon,sh = kon,m.
    lkon <- fixed(log(2.52)); label("Association rate constant of antibody and target, kon,p = kon,sh = kon,m (1/(nM*h))") # Table S2 'kon,p' = 2.52; 'KD,sh' = 'KD,m' = 1.26/2.52
    lkoff <- fixed(log(1.26)); label("Dissociation rate constant of the antibody-target complexes, koff,p = koff,sh = koff,m (1/h)") # Table S2 'koff,p' = 1.26
    lkdeg <- fixed(log(0.1)); label("First-order degradation rate constant of free shed target in plasma kdeg,p (1/h)") # Table S2 'kdeg,p' = 0.1
    lkint <- fixed(log(0.01)); label("Internalisation rate constant of the antibody-shed target complex in plasma kint,p (1/h)") # Table S2 'kint,p' = 0.01

    lksyn_mem <- fixed(log(0.83)); label("Zero-order synthesis rate of membrane-bound target in the interstitium ksyn,m (nM/h)") # Table S2 'ksyn,m' = 0.83
    lkdeg_mem <- fixed(log(0.1)); label("First-order degradation rate constant of membrane-bound target kdeg,m (1/h)") # Table S2 'kdeg,m' = 0.1
    lkint_mem <- fixed(log(0.1)); label("Internalisation rate constant of the antibody-membrane target complex kint,m (1/h)") # Table S2 'kint,m' = 0.1
    lksh <- fixed(log(0.138)); label("First-order shedding rate constant of membrane-bound target ksh (1/h)") # Table S2 'ksh' = 0.138
    lkdeg_shed <- fixed(log(0.1)); label("First-order degradation rate constant of free shed target in the interstitium kdeg,sh (1/h)") # Table S2 'kdeg,sh' = 0.1
    lkint_shed <- fixed(log(0.01)); label("Internalisation rate constant of the antibody-shed target complex in the interstitium kint,sh (1/h)") # Table S2 'kint,sh' = 0.01
    llambda_shed <- fixed(log(0.01)); label("Transfer rate constant of free shed target from the interstitium to plasma lambda (1/h)") # Table S2 'lambda' = 0.01 1/h
    alpha_shed <- fixed(0); label("Fraction of membrane complex returned as free shed target by shedding alpha (unitless)") # Table S2 'alpha' = 0

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
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    ksyn_mem <- exp(lksyn_mem)
    kdeg_mem <- exp(lkdeg_mem)
    kint_mem <- exp(lkint_mem)
    ksh <- exp(lksh)
    kdeg_shed <- exp(lkdeg_shed)
    kint_shed <- exp(lkint_shed)
    lambda_shed <- exp(llambda_shed)
    c_fcrn_t <- exp(lc_fcrn_t)
    kd_fcrn <- exp(lkoff_fcrn) / exp(lkon_fcrn) # KD,F = koff,F / kon,F
    krec <- exp(lkrec)

    # Equations 29 and 30: antibody-free steady state of the three free
    # target species at time zero.
    target_mem_isf(0) <- ksyn_mem / (kdeg_mem + ksh)
    target_isf(0) <- ksh / (kdeg_shed + lambda_shed) * ksyn_mem / (kdeg_mem + ksh)
    target(0) <- lambda_shed / kdeg * ksh / (kdeg_shed + lambda_shed) * ksyn_mem / (kdeg_mem + ksh) *
      v_tissue_isf / vc

    lig_p <- central / vc # Lp (nM)
    lig_v <- tissue_vasc / v_tissue_vasc # Lv (nM)
    lig_e <- endosome / v_endosome # Le (nM)
    lig_i <- tissue_isf / v_tissue_isf # Li (nM)
    lig_s <- lymph / v_lymph # Ls (nM)

    # Equation 18 (QSS FcRn binding): FcRn-bound antibody in the endothelium,
    # (1 - fu) * Le, as the smaller root of the binding quadratic.
    sum_fcrn <- lig_e + c_fcrn_t + kd_fcrn
    lig_e_bound <- 0.5 * (sum_fcrn - sqrt(sum_fcrn^2 - 4 * lig_e * c_fcrn_t))
    lig_e_free <- lig_e - lig_e_bound # fu * Le

    # Equations 23, 24, 25, 28 and 36, each antibody equation multiplied
    # through by its volume. The plasma target is now the shed target, which
    # arrives only from the interstitium (Equation 23 has no ksyn,p).
    d/dt(central) <- (qc - q_lymph) * lig_v + q_lymph * lig_s - qc * lig_p +
      vc * (-kon * lig_p * target + koff * complex - kel * lig_p)
    d/dt(target) <- lambda_shed * target_isf * v_tissue_isf / vc - kdeg * target -
      kon * lig_p * target + koff * complex
    d/dt(complex) <- kon * lig_p * target - (koff + kint) * complex
    d/dt(tissue_vasc) <- qc * lig_p - (qc - q_lymph) * lig_v - (1 - sigma_v) * q_lymph * lig_v -
      kup * tissue_vasc + fr * krec * lig_e_bound * v_endosome
    d/dt(endosome) <- kup * tissue_vasc + delta_kup * kup * tissue_isf - clint * lig_e_free -
      krec * lig_e_bound * v_endosome
    d/dt(tissue_isf) <- -delta_kup * kup * tissue_isf + (1 - sigma_v) * q_lymph * lig_v -
      (1 - sigma_i) * q_lymph * lig_i + (1 - fr) * krec * lig_e_bound * v_endosome +
      v_tissue_isf * (-kon * lig_i * target_isf + koff * complex_isf -
        kon * lig_i * target_mem_isf + koff * complex_mem_isf)
    d/dt(target_isf) <- ksh * target_mem_isf - (kdeg_shed + lambda_shed) * target_isf -
      kon * lig_i * target_isf + koff * complex_isf + alpha_shed * ksh * complex_mem_isf
    d/dt(complex_isf) <- kon * lig_i * target_isf - (koff + kint_shed) * complex_isf
    d/dt(target_mem_isf) <- ksyn_mem - (kdeg_mem + ksh) * target_mem_isf -
      kon * lig_i * target_mem_isf + koff * complex_mem_isf
    d/dt(complex_mem_isf) <- kon * lig_i * target_mem_isf - (koff + kint_mem) * complex_mem_isf -
      alpha_shed * ksh * complex_mem_isf
    d/dt(lymph) <- (1 - sigma_i) * q_lymph * lig_i - q_lymph * lig_s

    # Free trastuzumab in plasma: nM * g/mol = ng/L, / 1e6 = ug/mL.
    Cc <- lig_p * mw / 1e6
    Cc ~ prop(propSd)
  })
}
