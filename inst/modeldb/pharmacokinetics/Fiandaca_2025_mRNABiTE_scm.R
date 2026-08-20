Fiandaca_2025_mRNABiTE_scm <- function() {
  description <- paste0(
    "PBPK / QSP (whole-body two-pore, 15 tissues + tumour, 53 ODE states). ",
    "Preclinical (mouse, 8-10 week old immunodeficient tumour-bearing). ",
    "Multiscale model for an mRNA-encoded B7H3 x CD3 bispecific T-cell engager (BiTE) ",
    "delivered in a liver-tropic ionizable lipid nanoparticle. ",
    "A protein-transport module (Li and Shah 2019 two-pore whole-body PBPK, no FcRn) is ",
    "driven by an mRNA-transport module in which intravenously dosed mRNA migrates to the ",
    "liver and passes through 2 sequential intracellular states (short-chain model, SCM), ",
    "each translating BiTE into the hepatic interstitium at its own rate. ",
    "Dose plasma for the recombinant BiTE and mrna_plasma for the mRNA-encoded product; ",
    "concentrations are ng/mL and doses ug."
  )
  reference <- paste0(
    "Fiandaca G, Campanile E, Leonardelli L, Pettina E, Giampiccolo S, Carstens EJ, ",
    "Dasti L, Zangani N, Marchetti L. A multiscale physiologically based pharmacokinetic ",
    "model to support mRNA-encoded BiTE therapy in cancer treatment. Mol Ther Nucleic Acids. ",
    "2025;36(3):102606. doi:10.1016/j.omtn.2025.102606. ",
    "MATLAB source deposited by the authors at https://github.com/cosbi-research/PBPKmRNABiTE ",
    "(Data and code availability); every physiological constant and fitted estimate in this ",
    "file is traced to that deposit and/or paper Table 1."
  )
  vignette <- "Fiandaca_2025_mRNABiTE"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  # The mRNA staging chain is phenomenological: the authors introduce the
  # minimum number of liver mRNA states needed to reproduce the observed delay
  # in BiTE appearance, and explicitly decline to interpret the individual
  # rates ("we refrain from discussing the estimated values of rates").
  # These states therefore do not map onto a canonical compartment role.
  paper_specific_compartments <- c("mrna_plasma", "mrna_liver1", "mrna_liver2")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. All states are amounts in ug; concentrations are formed in
  # model() by dividing by the corresponding sub-compartment volume.
  compartmentData <- list(
    vp_heart            = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_lung             = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_muscle           = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_skin             = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_adipose          = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_bone             = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_brain            = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_kidney           = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_liver            = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_small_intestine  = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_large_intestine  = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_pancreas         = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_thymus           = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_spleen           = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    vp_other            = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    eu_heart            = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_lung             = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_muscle           = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_skin             = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_adipose          = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_bone             = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_brain            = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_kidney           = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_liver            = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_small_intestine  = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_large_intestine  = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_pancreas         = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_thymus           = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_spleen           = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    eu_other            = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    is_heart            = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_lung             = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_muscle           = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_skin             = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_adipose          = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_bone             = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_brain            = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_kidney           = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver            = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_small_intestine  = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_large_intestine  = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_pancreas         = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_thymus           = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_spleen           = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    is_other            = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    plasma              = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    lnode               = list(analyte = "BiTE", units = "ug", specimen = "lymph", verified = TRUE),
    vp_tumor            = list(analyte = "BiTE", units = "ug", specimen = "plasma", verified = TRUE),
    eu_tumor            = list(analyte = "BiTE", units = "ug", specimen = "endosome", verified = TRUE),
    is_tumor            = list(analyte = "BiTE", units = "ug", specimen = "tissue", verified = TRUE),
    mrna_plasma         = list(analyte = "mRNA", units = "ug", specimen = "plasma", verified = TRUE),
    mrna_liver1         = list(analyte = "mRNA", units = "ug", specimen = "tissue", verified = TRUE),
    mrna_liver2         = list(analyte = "mRNA", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (8-10 weeks old, immunodeficient, tumour-bearing)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "8-10 weeks",
    weight_range   = paste0(
      "0.028 kg used for the physiological parameter set (deposit pars.BW); the dose ",
      "calculation in the deposit simulation drivers uses 0.02 kg (6 mg/kg = 120 ug BiTE, ",
      "1.5 mg/kg = 30 ug mRNA-BiTE), matching the absolute doses reported by Huang 2023."
    ),
    disease_state  = "Subcutaneous tumour-bearing immunodeficient mice; a 0.472 mL tumour is carried as a 16th perfused tissue.",
    dose_range     = paste0(
      "Single IV bolus. Recombinant BiTE 6 mg/kg (120 ug). mRNA-BiTE 1.5 mg/kg (30 ug) for ",
      "training; 0.5, 1, 1.5 and 2 mg/kg for the dose-escalation validation."
    ),
    notes          = paste0(
      "Training and validation data are the published mouse serum BiTE profiles of ",
      "Huang et al. 2023 (the B7H3 x CD3 BiTE and its liver-tropic ionizable LNP-mRNA ",
      "formulation), re-digitised in the authors deposit. Serum was treated as equivalent ",
      "to plasma. Sampling at 1, 4, 6, 12, 24, 48, 72, 144 and 168 h for the time course; ",
      "a single 24 h sample for each dose-escalation arm."
    ),
    scope_note     = paste0(
      "Deterministic mechanistic platform: no between-subject random effects and no ",
      "residual-error magnitude are reported. The published uncertainty bands (Figures 4 ",
      "and 6B) come from a 1,000-subject virtual population in which the six (SCM) or ",
      "eight (LCM) fitted parameters were drawn from normal distributions centred on the ",
      "estimate with a relative standard deviation of 0.2; that is a parameter-uncertainty ",
      "sweep, not an estimated random-effects model, so it is reproduced in the vignette ",
      "rather than encoded as etas here."
    )
  )

  ini({
    # --- mRNA-transport module: fitted to the mRNA-BiTE profile ------------
    lkbl2lv          <- fixed(log(0.68035623)); label("Migration rate of mRNA from blood to liver (1/h)")            # Table 1 SCM 6.804e-1; deposit SCM_parameters_estimated.mat[2]
    lks12            <- fixed(log(0.48074765)); label("Transition rate, first to second liver mRNA state (1/h)")     # Table 1 SCM 4.807e-1; deposit .mat[6]
    lkrna            <- fixed(log(0.01083587)); label("mRNA clearance (degradation) rate (1/h)")                     # Table 1 SCM 1.08e-2; deposit .mat[3]
    lktr_mrna1       <- fixed(log(0.96823094)); label("BiTE translation rate from the first liver mRNA state (1/h)")  # Table 1 SCM 9.682e-1; deposit .mat[4]
    lktr_mrna2       <- fixed(log(0.31818579)); label("BiTE translation rate from the second liver mRNA state (1/h)") # Table 1 SCM 3.182e-1; deposit .mat[5]

    # --- Protein-transport module -----------------------------------------
    lclup            <- fixed(log(0.09999872)); label("Fluid-phase pinocytosis rate per unit endosomal volume (1/h)") # Table 1 SCM 1e-1; deposit .mat[1]

    # Fraction of the dosed mRNA delivered to organs other than the liver.
    # Fixed, not estimated: 'The parameter for the proportion of mRNA released
    # in organs other than the liver was set to 0.011 based on Huang et al.
    # (2023), indicating that approximately 99% of administered mRNA
    # accumulates in the liver after 6 h.'
    lfrac_mrna_other <- fixed(log(0.011));      label("Fraction of dosed mRNA delivered outside the liver (unitless)")   # Results, 'Model training and validation'; deposit pars.ratio_in_other

    # --- System parameters inherited from the upstream platform ------------
    # Supplied by the authors deposit rather than re-derived; the deposit cites
    # Li & Shah 2019 (itself built on Shah & Betts 2012) for each.
    lkdeg            <- fixed(log(32.2));       label("Endosomal (lysosomal) degradation rate of the BiTE (1/h)")        # deposit SCM_simulation.m K_deg, cited to Li & Shah 2019
    lclnlf           <- fixed(log(9.1));        label("Lymph-node-to-plasma flow scaling C_LNLF (unitless)")             # deposit SCM_simulation.m C_LNLF, cited to Shah & Betts 2012

    # No residual-error magnitude is reported anywhere in the paper or the
    # deposit: the model is fitted by least squares on relative deviations and
    # simulated deterministically.
    propSd           <- fixed(0);               label("Proportional residual error (fraction; magnitude not reported)")
  })

  model({
    # === Back-transform the fitted / inherited parameters ==================
    kbl2lv    <- exp(lkbl2lv)
    ks12      <- exp(lks12)
    krna      <- exp(lkrna)
    ktr_mrna1 <- exp(lktr_mrna1)
    ktr_mrna2 <- exp(lktr_mrna2)
    clup      <- exp(lclup)
    frac_mrna_other <- exp(lfrac_mrna_other)
    kdeg      <- exp(lkdeg)
    c_lnlf    <- exp(lclnlf)

    # === Provenance note on three deposit-vs-print discrepancies ============
    # 1. Two-pore convection. The convective coefficients below are written
    #    (1 - alpha_l) and (1 - alpha_s) -- the fractional hydraulic
    #    conductances -- rather than (1 - sigma_l) / (1 - sigma_s), the
    #    osmotic reflection coefficients that two-pore theory calls for. That
    #    is the authors' own deposit form: sigma_S = 0.906 and sigma_L =
    #    0.0877 are defined in Parameters_collection_BiTEs_55kDa.m and never
    #    referenced. Table 1 was fitted with this form, so it is reproduced
    #    here; substituting the reflection coefficients changes the plasma
    #    profile by at most 0.2%.
    # 2. Large intestine. The deposit drives the large-intestine two-pore
    #    terms off C_V_HEART instead of C_V_LG_INT (a copy-paste slip present
    #    in all three of its ODE files). The paper's Equations 3 and 5 are
    #    per-organ, so cv_lr is used below; the difference on plasma is
    #    < 0.03%.
    # 3. Printed Equations 7 and 8 (LCM) each carry a subscript typo
    #    (-k_s23 * mRNA_LV1 and -k_d * mRNA_LV2 respectively). Mass balance
    #    and the deposit both give -k_s23 * mRNA_LV2 and -k_d * mRNA_LV3;
    #    the corrected forms are used.
    # See the vignette's "Assumptions and deviations" section.

    # === Mouse physiology (Fiandaca 2025 deposit,
    # SCM_parameters_collection_mouse.m; inherited from Li & Shah 2019 /
    # Shah & Betts 2012). Volumes L, plasma flows L/h. ======================
    v_plasma <- 0.00094435          # V_TOT_PLASMA
    v_lnode  <- 0.000113            # V_TOT_LYMPH

    # heart
    v_ht_vp <- 5.852e-06; v_ht_is <- 2.1736e-05; v_ht_e <- 7.6e-07
    q_ht <- 3.6498e-02
    # lung
    v_lu_vp <- 2.945509193e-05; v_lu_is <- 3.84285724e-05; v_lu_e <- 1.0220365e-06
    # muscle
    v_mu_vp <- 2.49018e-04; v_mu_is <- 1.47147e-03; v_mu_e <- 5.6595e-05
    q_mu <- 8.613e-02
    # skin
    v_sk_vp <- 1.877106e-04; v_sk_is <- 1.65627e-03; v_sk_e <- 2.5095e-05
    q_sk <- 2.7819e-02
    # adipose
    v_ad_vp <- 2.1802e-05; v_ad_is <- 3.3694e-04; v_ad_e <- 9.91e-06
    q_ad <- 1.3431e-02
    # bone
    v_bo_vp <- 6.2128e-05; v_bo_is <- 5.25264e-04; v_bo_e <- 1.412e-05
    q_bo <- 1.518e-02
    # brain
    v_br_vp <- 1.067e-05; v_br_is <- 8.73e-05; v_br_e <- 2.425e-06
    q_br <- 1.1781e-02
    # kidney
    v_ki_vp <- 2.8875e-05; v_ki_is <- 7.875e-05; v_ki_e <- 2.625e-06
    q_ki <- 6.8508e-02
    # liver
    v_li_vp <- 1.6410625e-04; v_li_is <- 3.84595163e-04; v_li_e <- 9.625e-06
    q_li <- 1.0263e-02
    # small_intestine
    v_si_vp <- 1.16116e-05; v_si_is <- 1.26672e-04; v_si_e <- 3.64e-06
    q_si <- 5.808e-02
    # large_intestine
    v_lr_vp <- 5.000325e-06; v_lr_is <- 5.4549e-05; v_lr_e <- 1.5675e-06
    q_lr <- 1.7292e-02
    # pancreas
    v_pa_vp <- 5.335e-06; v_pa_is <- 1.6878e-05; v_pa_e <- 4.85e-07
    q_pa <- 6.237e-03
    # thymus
    v_th_vp <- 4.95e-07; v_th_is <- 1.53e-06; v_th_e <- 4.5e-08
    q_th <- 1.188e-03
    # spleen
    v_sp_vp <- 1.5367e-05; v_sp_is <- 2.54e-05; v_sp_e <- 6.35e-07
    q_sp <- 8.184e-03
    # other
    v_ot_vp <- 1.9535609097e-05; v_ot_is <- 7.97245183e-05; v_ot_e <- 2.326591e-06
    q_ot <- 1.254e-02

    # Tumour: V_TUMOUR with fixed vascular / endosomal / interstitial fractions
    v_tumour <- 0.000472
    v_tu_vp  <- 0.07  * v_tumour; v_tu_e <- 0.005 * v_tumour; v_tu_is <- 0.55 * v_tumour
    q_tu     <- 0.0127

    # Lung plasma flow closes the arterial mass balance: the sum of every
    # tissue's arterial supply, tumour included (deposit: pars.PLQ_LUNG).
    q_lu <- q_ht + q_mu + q_sk + q_ad + q_bo + q_br + q_ki + q_li + q_si + q_lr + q_pa + q_th + q_sp + q_ot + q_tu

    # Lymph flow is 0.2% of the tissue's plasma flow (deposit: L_ORG =
    # 0.002 * PLQ_ORG). The liver takes 0.2% of its TOTAL (arterial +
    # portal) inflow so its vascular outflow closes like every other organ.
    l_ht <- 0.002 * q_ht; l_mu <- 0.002 * q_mu; l_sk <- 0.002 * q_sk; l_ad <- 0.002 * q_ad; l_bo <- 0.002 * q_bo; l_br <- 0.002 * q_br; l_ki <- 0.002 * q_ki; l_si <- 0.002 * q_si; l_lr <- 0.002 * q_lr; l_pa <- 0.002 * q_pa; l_th <- 0.002 * q_th; l_sp <- 0.002 * q_sp; l_ot <- 0.002 * q_ot; l_lu <- 0.002 * q_lu; l_tu <- 0.002 * q_tu
    q_li_upstream <- (q_si - l_si) + (q_lr - l_lr) + (q_pa - l_pa) + (q_sp - l_sp)
    l_li <- 0.002 * (q_li + q_li_upstream)

    # Lymph-node return flow (deposit: L_LYMPH = PLQ_LUNG * C_LNLF)
    l_lnode <- c_lnlf * q_lu

    # Endosomal pinocytosis clearance per organ: CL_UP_ORG = CL_up * V_E_ORG
    clu_ht <- clup * v_ht_e; clu_lu <- clup * v_lu_e; clu_mu <- clup * v_mu_e; clu_sk <- clup * v_sk_e; clu_ad <- clup * v_ad_e; clu_bo <- clup * v_bo_e; clu_br <- clup * v_br_e; clu_ki <- clup * v_ki_e; clu_li <- clup * v_li_e; clu_si <- clup * v_si_e; clu_lr <- clup * v_lr_e; clu_pa <- clup * v_pa_e; clu_th <- clup * v_th_e; clu_sp <- clup * v_sp_e; clu_ot <- clup * v_ot_e; clu_tu <- clup * v_tu_e

    # === 55 kDa BiTE two-pore transport constants ==========================
    # Parameters_collection_BiTEs_55kDa.m (Li & Shah 2019 Table 4, 55 kDa).
    a_e     <- 3.26        # Stokes-Einstein radius (nm)
    theta   <- 0.0986      # glomerular sieving coefficient
    aona0_s <- 2.45e-03    # fractional accessible pore area, small pores
    aona0_l <- 0.522       # fractional accessible pore area, large pores
    rr_s    <- a_e / 0.735 # protein-to-small-pore size ratio
    rr_l    <- a_e / 0.143 # protein-to-large-pore size ratio
    pe_s    <- 0.113               # small-pore Peclet number
    pe_l    <- 2.254425255333927   # large-pore Peclet number
    alpha_l <- 0.042       # fractional hydraulic conductance, large pores
    alpha_s <- 0.958       # fractional hydraulic conductance, small pores
    x_j     <- 0.38        # relative pore-abundance constant
    xp      <- 13197       # pore-size / conductance constant (nm^3)
    gfr     <- 0.0167      # mouse glomerular filtration rate (L/h)
    sig_is  <- 0.2         # interstitial (lymphatic) reflection coefficient

    cl_r <- gfr * theta    # renal clearance of the BiTE (L/h)
    xp_s <- xp * (1 / a_e) * aona0_s * (alpha_s / rr_s^2)
    xp_l <- xp * (1 / a_e) * aona0_l * (alpha_l / rr_l^2)
    pe_s_ratio <- pe_s / (exp(pe_s) - 1)
    pe_l_ratio <- pe_l / (exp(pe_l) - 1)

    # === Concentration aliases (amount / volume, ug/L = ng/mL) ==============
    cv_p  <- plasma / v_plasma
    cv_ln <- lnode  / v_lnode
    cv_ht <- vp_heart / v_ht_vp; cis_ht <- is_heart / v_ht_is
    cv_lu <- vp_lung / v_lu_vp; cis_lu <- is_lung / v_lu_is
    cv_mu <- vp_muscle / v_mu_vp; cis_mu <- is_muscle / v_mu_is
    cv_sk <- vp_skin / v_sk_vp; cis_sk <- is_skin / v_sk_is
    cv_ad <- vp_adipose / v_ad_vp; cis_ad <- is_adipose / v_ad_is
    cv_bo <- vp_bone / v_bo_vp; cis_bo <- is_bone / v_bo_is
    cv_br <- vp_brain / v_br_vp; cis_br <- is_brain / v_br_is
    cv_ki <- vp_kidney / v_ki_vp; cis_ki <- is_kidney / v_ki_is
    cv_li <- vp_liver / v_li_vp; cis_li <- is_liver / v_li_is
    cv_si <- vp_small_intestine / v_si_vp; cis_si <- is_small_intestine / v_si_is
    cv_lr <- vp_large_intestine / v_lr_vp; cis_lr <- is_large_intestine / v_lr_is
    cv_pa <- vp_pancreas / v_pa_vp; cis_pa <- is_pancreas / v_pa_is
    cv_th <- vp_thymus / v_th_vp; cis_th <- is_thymus / v_th_is
    cv_sp <- vp_spleen / v_sp_vp; cis_sp <- is_spleen / v_sp_is
    cv_ot <- vp_other / v_ot_vp; cis_ot <- is_other / v_ot_is
    cv_tu <- vp_tumor / v_tu_vp; cis_tu <- is_tumor / v_tu_is

    # === Two-pore transcapillary clearance (Li & Shah 2019 Eq 9 + 10) ======
    # Per organ: PS_L = xp_l * L_ORG, PS_S = xp_s * L_ORG, and the
    # isogravimetric-flow split J_L = x_j*L + alpha_l*L, J_S = -x_j*L + alpha_s*L.
    # The convective coefficients are the deposit's own form (note 1 above).
    cltp_ht <- xp_l * l_ht * (cv_ht - cis_ht) * pe_l_ratio + (x_j * l_ht + alpha_l * l_ht) * (1 - alpha_l) * cv_ht +
      xp_s * l_ht * (cv_ht - cis_ht) * pe_s_ratio + (-x_j * l_ht + alpha_s * l_ht) * (1 - alpha_s) * cv_ht
    cltp_lu <- xp_l * l_lu * (cv_lu - cis_lu) * pe_l_ratio + (x_j * l_lu + alpha_l * l_lu) * (1 - alpha_l) * cv_lu +
      xp_s * l_lu * (cv_lu - cis_lu) * pe_s_ratio + (-x_j * l_lu + alpha_s * l_lu) * (1 - alpha_s) * cv_lu
    cltp_mu <- xp_l * l_mu * (cv_mu - cis_mu) * pe_l_ratio + (x_j * l_mu + alpha_l * l_mu) * (1 - alpha_l) * cv_mu +
      xp_s * l_mu * (cv_mu - cis_mu) * pe_s_ratio + (-x_j * l_mu + alpha_s * l_mu) * (1 - alpha_s) * cv_mu
    cltp_sk <- xp_l * l_sk * (cv_sk - cis_sk) * pe_l_ratio + (x_j * l_sk + alpha_l * l_sk) * (1 - alpha_l) * cv_sk +
      xp_s * l_sk * (cv_sk - cis_sk) * pe_s_ratio + (-x_j * l_sk + alpha_s * l_sk) * (1 - alpha_s) * cv_sk
    cltp_ad <- xp_l * l_ad * (cv_ad - cis_ad) * pe_l_ratio + (x_j * l_ad + alpha_l * l_ad) * (1 - alpha_l) * cv_ad +
      xp_s * l_ad * (cv_ad - cis_ad) * pe_s_ratio + (-x_j * l_ad + alpha_s * l_ad) * (1 - alpha_s) * cv_ad
    cltp_bo <- xp_l * l_bo * (cv_bo - cis_bo) * pe_l_ratio + (x_j * l_bo + alpha_l * l_bo) * (1 - alpha_l) * cv_bo +
      xp_s * l_bo * (cv_bo - cis_bo) * pe_s_ratio + (-x_j * l_bo + alpha_s * l_bo) * (1 - alpha_s) * cv_bo
    cltp_br <- xp_l * l_br * (cv_br - cis_br) * pe_l_ratio + (x_j * l_br + alpha_l * l_br) * (1 - alpha_l) * cv_br +
      xp_s * l_br * (cv_br - cis_br) * pe_s_ratio + (-x_j * l_br + alpha_s * l_br) * (1 - alpha_s) * cv_br
    cltp_ki <- xp_l * l_ki * (cv_ki - cis_ki) * pe_l_ratio + (x_j * l_ki + alpha_l * l_ki) * (1 - alpha_l) * cv_ki +
      xp_s * l_ki * (cv_ki - cis_ki) * pe_s_ratio + (-x_j * l_ki + alpha_s * l_ki) * (1 - alpha_s) * cv_ki
    cltp_li <- xp_l * l_li * (cv_li - cis_li) * pe_l_ratio + (x_j * l_li + alpha_l * l_li) * (1 - alpha_l) * cv_li +
      xp_s * l_li * (cv_li - cis_li) * pe_s_ratio + (-x_j * l_li + alpha_s * l_li) * (1 - alpha_s) * cv_li
    cltp_si <- xp_l * l_si * (cv_si - cis_si) * pe_l_ratio + (x_j * l_si + alpha_l * l_si) * (1 - alpha_l) * cv_si +
      xp_s * l_si * (cv_si - cis_si) * pe_s_ratio + (-x_j * l_si + alpha_s * l_si) * (1 - alpha_s) * cv_si
    cltp_lr <- xp_l * l_lr * (cv_lr - cis_lr) * pe_l_ratio + (x_j * l_lr + alpha_l * l_lr) * (1 - alpha_l) * cv_lr +
      xp_s * l_lr * (cv_lr - cis_lr) * pe_s_ratio + (-x_j * l_lr + alpha_s * l_lr) * (1 - alpha_s) * cv_lr
    cltp_pa <- xp_l * l_pa * (cv_pa - cis_pa) * pe_l_ratio + (x_j * l_pa + alpha_l * l_pa) * (1 - alpha_l) * cv_pa +
      xp_s * l_pa * (cv_pa - cis_pa) * pe_s_ratio + (-x_j * l_pa + alpha_s * l_pa) * (1 - alpha_s) * cv_pa
    cltp_th <- xp_l * l_th * (cv_th - cis_th) * pe_l_ratio + (x_j * l_th + alpha_l * l_th) * (1 - alpha_l) * cv_th +
      xp_s * l_th * (cv_th - cis_th) * pe_s_ratio + (-x_j * l_th + alpha_s * l_th) * (1 - alpha_s) * cv_th
    cltp_sp <- xp_l * l_sp * (cv_sp - cis_sp) * pe_l_ratio + (x_j * l_sp + alpha_l * l_sp) * (1 - alpha_l) * cv_sp +
      xp_s * l_sp * (cv_sp - cis_sp) * pe_s_ratio + (-x_j * l_sp + alpha_s * l_sp) * (1 - alpha_s) * cv_sp
    cltp_ot <- xp_l * l_ot * (cv_ot - cis_ot) * pe_l_ratio + (x_j * l_ot + alpha_l * l_ot) * (1 - alpha_l) * cv_ot +
      xp_s * l_ot * (cv_ot - cis_ot) * pe_s_ratio + (-x_j * l_ot + alpha_s * l_ot) * (1 - alpha_s) * cv_ot
    cltp_tu <- xp_l * l_tu * (cv_tu - cis_tu) * pe_l_ratio + (x_j * l_tu + alpha_l * l_tu) * (1 - alpha_l) * cv_tu +
      xp_s * l_tu * (cv_tu - cis_tu) * pe_s_ratio + (-x_j * l_tu + alpha_s * l_tu) * (1 - alpha_s) * cv_tu

    # === Tissue vascular spaces (deposit dydt(1:15), Eq 3) =================
    # Arterial supply is the lung vascular concentration for every organ
    # except the lung (supplied from central plasma) and the liver (which
    # also receives the portal return from gut, pancreas and spleen).
    d/dt(vp_heart) <- q_ht * cv_lu - (q_ht - l_ht) * cv_ht - cltp_ht - clu_ht * cv_ht
    d/dt(vp_lung) <- (q_lu + l_lu) * cv_p - q_lu * cv_lu - cltp_lu - clu_lu * cv_lu
    d/dt(vp_muscle) <- q_mu * cv_lu - (q_mu - l_mu) * cv_mu - cltp_mu - clu_mu * cv_mu
    d/dt(vp_skin) <- q_sk * cv_lu - (q_sk - l_sk) * cv_sk - cltp_sk - clu_sk * cv_sk
    d/dt(vp_adipose) <- q_ad * cv_lu - (q_ad - l_ad) * cv_ad - cltp_ad - clu_ad * cv_ad
    d/dt(vp_bone) <- q_bo * cv_lu - (q_bo - l_bo) * cv_bo - cltp_bo - clu_bo * cv_bo
    d/dt(vp_brain) <- q_br * cv_lu - (q_br - l_br) * cv_br - cltp_br - clu_br * cv_br
    # Kidney additionally clears by glomerular filtration (Eq 7)
    d/dt(vp_kidney) <- q_ki * cv_lu - (q_ki - l_ki) * cv_ki - cltp_ki -
      (clu_ki + cl_r) * cv_ki
    qvout_li <- (q_li - l_li) + q_li_upstream
    d/dt(vp_liver) <- q_li * cv_lu + (q_si - l_si) * cv_si + (q_lr - l_lr) * cv_lr + (q_pa - l_pa) * cv_pa + (q_sp - l_sp) * cv_sp -
      qvout_li * cv_li - cltp_li - clu_li * cv_li
    d/dt(vp_small_intestine) <- q_si * cv_lu - (q_si - l_si) * cv_si - cltp_si - clu_si * cv_si
    d/dt(vp_large_intestine) <- q_lr * cv_lu - (q_lr - l_lr) * cv_lr - cltp_lr - clu_lr * cv_lr
    d/dt(vp_pancreas) <- q_pa * cv_lu - (q_pa - l_pa) * cv_pa - cltp_pa - clu_pa * cv_pa
    d/dt(vp_thymus) <- q_th * cv_lu - (q_th - l_th) * cv_th - cltp_th - clu_th * cv_th
    d/dt(vp_spleen) <- q_sp * cv_lu - (q_sp - l_sp) * cv_sp - cltp_sp - clu_sp * cv_sp
    d/dt(vp_other) <- q_ot * cv_lu - (q_ot - l_ot) * cv_ot - cltp_ot - clu_ot * cv_ot

    # === Tissue endosomal spaces (deposit dydt(16:30), Eq 4) ===============
    # Fluid-phase pinocytosis from both the vascular and interstitial side,
    # then lysosomal degradation. No FcRn: the BiTE carries no Fc sequence,
    # so there is no recycling arm and the endosome is a pure sink.
    d/dt(eu_heart) <- clu_ht * (cv_ht + cis_ht) - kdeg * eu_heart
    d/dt(eu_lung) <- clu_lu * (cv_lu + cis_lu) - kdeg * eu_lung
    d/dt(eu_muscle) <- clu_mu * (cv_mu + cis_mu) - kdeg * eu_muscle
    d/dt(eu_skin) <- clu_sk * (cv_sk + cis_sk) - kdeg * eu_skin
    d/dt(eu_adipose) <- clu_ad * (cv_ad + cis_ad) - kdeg * eu_adipose
    d/dt(eu_bone) <- clu_bo * (cv_bo + cis_bo) - kdeg * eu_bone
    d/dt(eu_brain) <- clu_br * (cv_br + cis_br) - kdeg * eu_brain
    d/dt(eu_kidney) <- clu_ki * (cv_ki + cis_ki) - kdeg * eu_kidney
    d/dt(eu_liver) <- clu_li * (cv_li + cis_li) - kdeg * eu_liver
    d/dt(eu_small_intestine) <- clu_si * (cv_si + cis_si) - kdeg * eu_small_intestine
    d/dt(eu_large_intestine) <- clu_lr * (cv_lr + cis_lr) - kdeg * eu_large_intestine
    d/dt(eu_pancreas) <- clu_pa * (cv_pa + cis_pa) - kdeg * eu_pancreas
    d/dt(eu_thymus) <- clu_th * (cv_th + cis_th) - kdeg * eu_thymus
    d/dt(eu_spleen) <- clu_sp * (cv_sp + cis_sp) - kdeg * eu_spleen
    d/dt(eu_other) <- clu_ot * (cv_ot + cis_ot) - kdeg * eu_other

    # === Tissue interstitial spaces (deposit dydt(31:45), Eq 5) ============
    d/dt(is_heart) <- cltp_ht - (1 - sig_is) * l_ht * cis_ht - clu_ht * cis_ht
    d/dt(is_lung) <- cltp_lu - (1 - sig_is) * l_lu * cis_lu - clu_lu * cis_lu
    d/dt(is_muscle) <- cltp_mu - (1 - sig_is) * l_mu * cis_mu - clu_mu * cis_mu
    d/dt(is_skin) <- cltp_sk - (1 - sig_is) * l_sk * cis_sk - clu_sk * cis_sk
    d/dt(is_adipose) <- cltp_ad - (1 - sig_is) * l_ad * cis_ad - clu_ad * cis_ad
    d/dt(is_bone) <- cltp_bo - (1 - sig_is) * l_bo * cis_bo - clu_bo * cis_bo
    d/dt(is_brain) <- cltp_br - (1 - sig_is) * l_br * cis_br - clu_br * cis_br
    d/dt(is_kidney) <- cltp_ki - (1 - sig_is) * l_ki * cis_ki - clu_ki * cis_ki
    d/dt(is_liver) <- cltp_li - (1 - sig_is) * l_li * cis_li - clu_li * cis_li +
      ktr_mrna1 * mrna_liver1 + ktr_mrna2 * mrna_liver2
    d/dt(is_small_intestine) <- cltp_si - (1 - sig_is) * l_si * cis_si - clu_si * cis_si
    d/dt(is_large_intestine) <- cltp_lr - (1 - sig_is) * l_lr * cis_lr - clu_lr * cis_lr
    d/dt(is_pancreas) <- cltp_pa - (1 - sig_is) * l_pa * cis_pa - clu_pa * cis_pa
    d/dt(is_thymus) <- cltp_th - (1 - sig_is) * l_th * cis_th - clu_th * cis_th
    d/dt(is_spleen) <- cltp_sp - (1 - sig_is) * l_sp * cis_sp - clu_sp * cis_sp
    d/dt(is_other) <- cltp_ot - (1 - sig_is) * l_ot * cis_ot - clu_ot * cis_ot

    # === Central plasma (deposit dydt(46), Eq 1) ===========================
    # Venous return from gut, pancreas and spleen goes to the liver, not here.
    d/dt(plasma) <- (q_ht - l_ht) * cv_ht +
      (q_mu - l_mu) * cv_mu +
      (q_sk - l_sk) * cv_sk +
      (q_ad - l_ad) * cv_ad +
      (q_bo - l_bo) * cv_bo +
      (q_br - l_br) * cv_br +
      (q_ki - l_ki) * cv_ki +
      (q_th - l_th) * cv_th +
      (q_ot - l_ot) * cv_ot +
      qvout_li * cv_li + (q_tu - l_tu) * cv_tu + l_lnode * cv_ln -
      (q_lu + l_lu) * cv_p

    # === Lymph node (deposit dydt(47), Eq 2) ===============================
    # The tumour's interstitial lymph is NOT returned here -- see model note.
    d/dt(lnode) <- (1 - sig_is) * (
      l_ht * cis_ht + l_lu * cis_lu + l_mu * cis_mu + l_sk * cis_sk + l_ad * cis_ad + l_bo * cis_bo + l_br * cis_br + l_ki * cis_ki + l_li * cis_li + l_si * cis_si + l_lr * cis_lr + l_pa * cis_pa + l_th * cis_th + l_sp * cis_sp + l_ot * cis_ot
    ) - l_lnode * cv_ln

    # === Tumour (deposit dydt(48:50)) ======================================
    d/dt(vp_tumor) <- q_tu * cv_lu - (q_tu - l_tu) * cv_tu - cltp_tu - clu_tu * cv_tu
    d/dt(eu_tumor) <- clu_tu * (cv_tu + cis_tu) - kdeg * eu_tumor
    d/dt(is_tumor) <- cltp_tu - (1 - sig_is) * l_tu * cis_tu - clu_tu * cis_tu

    # === mRNA-transport module (Equations 1-3, deposit dydt(51:53)) ================
    # Amounts (ug). Writing the chain in amounts rather than the deposit's
    # concentrations makes the V_BL^TOT / V_LV^IS ratio of Equation 2 / 6
    # cancel exactly, so no volume term appears here.
    d/dt(mrna_plasma) <- -(krna + kbl2lv) * mrna_plasma
    d/dt(mrna_liver1) <- kbl2lv * (1 - frac_mrna_other) * mrna_plasma -
      (ks12 + krna) * mrna_liver1
    d/dt(mrna_liver2) <- ks12 * mrna_liver1 - krna * mrna_liver2

    # === Observation =======================================================
    # Central plasma BiTE concentration in ng/mL (ug amount / L volume).
    Cc <- cv_p
    Cc ~ prop(propSd)
  })
}
