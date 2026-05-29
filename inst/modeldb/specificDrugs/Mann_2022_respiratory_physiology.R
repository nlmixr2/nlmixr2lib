Mann_2022_respiratory_physiology <- function() {
  description <- paste(
    "QSP. Magosso / Ursino respiratory and cerebrovascular physiology",
    "with Mann 2022 extensions for opioid-induced ventilatory",
    "depression and cardiovascular-collapse / cardiac-arrest dynamics.",
    "Encodes the 11-state physiological submodel from the FDA",
    "delaymymod.c implementation, plus the cardiac-arrest event rule",
    "(PaO2 below 15 mm Hg sustained 220 s -> cardiac output decays",
    "toward 0.01 L/min). The CAR (fraction of opioid receptors bound",
    "by an agonist) input drives reductions in wakefulness drive",
    "(W - Wmax * CAR^P3) and chemoreflex drives (factor 1 - CAR^P1).",
    "The Spencer dissociation algebra for blood gas exchange is",
    "carried inline. The original FDA implementation uses delay-",
    "differential equations for peripheral and central chemoreflex",
    "filtering with delays of roughly K_Dp/(Qb+Qt) ~ 7 s and",
    "K_Dc/(Qb+Qt) ~ 11 s; this model deploys the limit of zero delay",
    "(Plag_X = X, Clag_X = X), which preserves the steady-state",
    "structure and longer-time-scale overdose dynamics but does not",
    "reproduce the second-scale delay artefacts of the original.",
    "Composes downstream of Mann_2022_mu_receptor_binding (CAR_OPIOID",
    "input)."
  )
  reference <- paste(
    "Mann J, Samieegohar M, Chaturbedi A, Zirkle J, Han X, Ahmadi SF,",
    "Eshleman A, Janowsky A, Wolfrum K, Swanson T, Bloom S, Dahan A,",
    "Olofsen E, Florian J, Strauss DG, Li Z.",
    "Development of a Translational Model to Assess the Impact of Opioid",
    "Overdose and Naloxone Dosing on Respiratory Depression and",
    "Cardiac Arrest. Clin Pharmacol Ther. 2022;112(5):1020-1032.",
    "doi:10.1002/cpt.2696. PMID 35766413.",
    "Upstream Magosso / Ursino respiratory + cerebrovascular base model:",
    "Magosso E, Ursino M, van Oostrom JH. Opioid-induced respiratory",
    "depression: a mathematical model for fentanyl.",
    "IEEE Trans Biomed Eng. 2004;51(7):1115-1128.",
    "doi:10.1109/TBME.2004.827330. Plus Ursino M, Magosso E,",
    "Avanzolini G. An integrated model of the human ventilatory control",
    "system: the response to hypercapnia / hypoxia. Clin Physiol.",
    "2001;21(4):447-477.",
    "Upstream cerebral-blood-flow control: Duffin J, Hare GMT, Fisher",
    "JA. A mathematical model of cerebral blood flow control in anaemia",
    "and hypoxia. J Physiol. 2020;598(4):717-730.",
    "Authoritative implementation: FDA delaymymod.c at",
    "https://github.com/FDA/Mechanistic-PK-PD-Model-to-Rescue-Opioid-Overdose,",
    "with delaypars.R (parameter table) and fentanyl_pars.txt (P1, P2",
    "overrides for naive vs chronic patient types)."
  )
  units <- list(
    time          = "min",
    dosing        = "(not applicable; CAR is a time-varying covariate input from the binding layer)",
    concentration = "L/min (primary output Venti is minute ventilation; the model has no plasma-concentration analyte and the canonical Cc observation is not used)"
  )

  covariateData <- list(
    CAR_OPIOID = list(
      description        = "Time-varying fraction of mu-opioid receptors bound by an opioid agonist (RL_op output of Mann_2022_mu_receptor_binding); drives the ventilatory-depression reductions of wakefulness and chemoreflex drives.",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Range 0..1. In the composed Mann 2022 overdose chain, this is",
        "the output of the binding layer Mann_2022_mu_receptor_binding",
        "for the OPIOID slot (RL_op). At CAR_OPIOID = 0 the model",
        "reduces to baseline Magosso / Ursino respiratory physiology",
        "(eupnoea); at CAR_OPIOID = 1 the wakefulness drive is",
        "completely abolished and the chemoreflex drives are also",
        "fully suppressed (1 - 1^P1 = 0), driving terminal apnoea."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    ),
    OPIOID_PATIENT_TYPE = list(
      description        = "Categorical 0 = healthy opioid-naive volunteer, 1 = chronic opioid user; selects the P1 (chemoreflex-drive reduction exponent) and P3 (wakefulness-drive reduction exponent) sensitivity parameters.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (opioid-naive)",
      notes              = paste(
        "Mann 2022 fits two pharmacodynamic parameter sets:",
        "naive (P1 = 2.875, P3 = 0.9) and chronic (P1 = 4.226,",
        "P3 = 1.323). Both share the same metabolism exponent",
        "P2 = 0.06319. Sources for the naive vs chronic split: Algera",
        "MH et al. Clin Pharmacol Ther 2021;109(3):637-645",
        "(chronic-user respiratory-depression tolerance data).",
        "Implementation source for the numeric values: FDA",
        "simulateToGetOD_IM.R lines 185-192 (allpatients hand-",
        "adjustment after fentanyl_pars.txt load)."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = "70 kg adult (Magosso / Ursino reference)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy adult ventilatory physiology calibrated to clinical",
      "ventilation responses to hypercapnia, hypoxia, anesthesia, and",
      "extended cardiovascular collapse from severe acute hypoxemia."
    ),
    dose_range     = NA_character_,
    regions        = NA_character_,
    notes          = paste(
      "Eleven physiological states are integrated:",
      "palv_co2, palv_o2 (alveolar / arterial gas partial pressures);",
      "cb_co2, cb_o2 (brain blood-gas concentrations);",
      "ct_co2, ct_o2 (peripheral-tissue blood-gas concentrations);",
      "yco2, yo2 (filtered peripheral-chemoreflex inputs);",
      "dp_state, dc_state (filtered peripheral and central chemoreflex drives);",
      "alpha_h (central ventilatory depression factor).",
      "Two cardiac-arrest event states are added:",
      "t_pao2_below (accumulated time arterial O2 partial pressure has",
      "spent below the critical 15 mm Hg threshold, in minutes), and",
      "im_arrest (multiplier on cardiac output Qb + Qt during",
      "post-trigger decay).",
      "Spencer / Spence dissociation algebra is computed inline as",
      "intermediate quantities in model() (not as states).",
      "DEVIATION from the FDA reference implementation:",
      "the original delaymymod.c uses delay-differential equations",
      "for the chemoreflex drives via deSolve::dede + lagvalue() with",
      "peripheralDelay = K_Dp/(Qb+Qt) ~ 7 s and centralDelay =",
      "K_Dc/(Qb+Qt) ~ 11 s. This nlmixr2lib encoding instead applies",
      "the zero-delay limit Plag_X = X, Clag_X = X. Effect on overdose",
      "dynamics: short-time-scale chemoreflex perturbations are",
      "smoothed; minute-to-hour-scale outcomes (time to cardiac arrest,",
      "PaO2 trajectory, ventilation suppression) are well-preserved.",
      "P1 (chemoreflex sensitivity exponent) and P3 (wakefulness",
      "sensitivity exponent) labelling note: the manuscript text says",
      "'P1 for wakefulness, P3 for chemoreflex' but the FDA code",
      "implements the opposite mapping (P3 acts on wakefulness via Kf",
      "and P1 acts on chemoreflex drives via AV1, AV2). The code is",
      "the authoritative implementation; the model() block here matches",
      "the code mapping."
    )
  )

  ini({
    # ===== Pharmacodynamic sensitivity exponents (FDA fentanyl_pars.txt
    # + simulateToGetOD_IM.R hand-adjustments lines 185-192) =====
    #
    # Two patient subtypes: opioid-naive (PATIENT_TYPE = 0) and chronic
    # opioid user (PATIENT_TYPE = 1). The naive / chronic split is
    # implemented as paired parameters selected by OPIOID_PATIENT_TYPE
    # at runtime.

    p1_naive    <- fixed(2.875)
    label("Chemoreflex-drive sensitivity exponent P1 in opioid-naive subjects (unitless, FIXED)")  # FDA simulateToGetOD_IM.R line 188: P1 = 2.5 * 1.15 = 2.875 (naive)
    p1_chronic  <- fixed(4.22625)
    label("Chemoreflex-drive sensitivity exponent P1 in chronic opioid users (unitless, FIXED)")    # FDA simulateToGetOD_IM.R line 191: P1 = 2.5 * 1.47 * 1.15 = 4.22625 (chronic)
    p3_naive    <- fixed(0.9)
    label("Wakefulness-drive sensitivity exponent P3 in opioid-naive subjects (unitless, FIXED)")   # FDA simulateToGetOD_IM.R line 189: P3 = 0.9 (naive)
    p3_chronic  <- fixed(1.323)
    label("Wakefulness-drive sensitivity exponent P3 in chronic opioid users (unitless, FIXED)")    # FDA simulateToGetOD_IM.R line 191: P3 = 0.9 * 1.47 = 1.323 (chronic)
    p2          <- fixed(0.06319)
    label("Metabolism scaling exponent P2 (unitless, FIXED, both patient types)")                   # FDA simulateToGetOD_IM.R line 186: P2 = 0.06319

    # ===== Wakefulness drive (Magosso / Ursino) =====
    wmax_drive  <- fixed(6.62)
    label("Maximum wakefulness drive Wmax (L/min, FIXED)")                                          # FDA delaypars.R: Wmax = 6.62
    w_baseline  <- fixed(6.62)
    label("Baseline wakefulness drive W (L/min, FIXED)")                                            # FDA delaypars.R: W = 6.62
    # NOTE: FDA delaypars.R also defines Bmax = 0.66 used internally as
    # Venti = totalVentilation * Bmax / 60 -- this scaling is undone at
    # the output stage (yout[0] = Venti * 60 / 0.66 = totalVentilation),
    # so on the L/min time scale used here the net Bmax factor is 1 and
    # is not encoded as an ini() parameter.

    # ===== Alveolar / gas-exchange parameters =====
    v_a         <- fixed(3.28)
    label("Alveolar volume V_A (L, FIXED)")                                                         # FDA delaypars.R: V_A = 3.28
    p_i_co2     <- fixed(0.0)
    label("Inspired CO2 partial pressure P_I_co2 (mm Hg, FIXED room-air = 0)")                     # FDA delaypars.R: P_I_co2 = 0
    p_i_o2_init <- fixed(149.0)
    label("Inspired O2 partial pressure P_I_o2 (mm Hg, FIXED room-air sea-level)")                 # FDA delaystates.R initial state P_I_o2 = 149
    lumbda      <- fixed(863.0)
    label("Blood-to-gas partition coefficient lambda (unitless, FIXED, Henry's law family)")        # FDA delaypars.R: lumbda = 863
    s1_shunt    <- fixed(0.024)
    label("Pulmonary venous-admixture shunt fraction s1 (unitless, FIXED)")                        # FDA delaypars.R: s1 = 0.024

    # ===== Blood-gas (Spencer) dissociation curve constants =====
    alpha_co2   <- fixed(6.67e-4)
    label("CO2 plasma solubility coefficient alpha_co2 (FIXED, Spencer)")                          # FDA delaypars.R: alpha_co2 = 6.67e-4
    alpha_o2    <- fixed(3.17e-5)
    label("O2 plasma solubility coefficient alpha_o2 (FIXED, Spencer)")                            # FDA delaypars.R: alpha_o2 = 3.17e-5
    c_b_hco3    <- fixed(0.5824)
    label("Brain bicarbonate buffer concentration C_B_hco3 (FIXED, Spencer)")                      # FDA delaypars.R: C_B_hco3 = 0.5824
    c_t_hco3    <- fixed(0.5824)
    label("Tissue bicarbonate buffer concentration C_T_hco3 (FIXED, Spencer)")                     # FDA delaypars.R: C_T_hco3 = 0.5824
    k1_spencer  <- fixed(14.99)
    label("Spencer constant K1 for O2 dissociation (FIXED)")                                       # FDA delaypars.R: K1 = 14.99
    k2_spencer  <- fixed(194.4)
    label("Spencer constant K2 for CO2 dissociation (FIXED)")                                      # FDA delaypars.R: K2 = 194.4
    a1_spencer  <- fixed(0.3836)
    label("Spencer exponent a1 for O2 dissociation (FIXED)")                                       # FDA delaypars.R: a1 = 0.3836
    a2_spencer  <- fixed(1.819)
    label("Spencer exponent a2 for CO2 dissociation (FIXED)")                                      # FDA delaypars.R: a2 = 1.819
    alpha1_sp   <- fixed(0.03198)
    label("Spencer cross-term alpha1 (CO2 effect on O2 binding, FIXED)")                           # FDA delaypars.R: alpha1 = 0.03198
    beta1_sp    <- fixed(0.008275)
    label("Spencer cross-term beta1 (CO2 effect on O2 binding, FIXED)")                            # FDA delaypars.R: beta1 = 0.008275
    alpha2_sp   <- fixed(0.05591)
    label("Spencer cross-term alpha2 (O2 effect on CO2 binding, FIXED)")                           # FDA delaypars.R: alpha2 = 0.05591
    beta2_sp    <- fixed(0.03255)
    label("Spencer cross-term beta2 (O2 effect on CO2 binding, FIXED)")                            # FDA delaypars.R: beta2 = 0.03255
    z_sp        <- fixed(0.02272)
    label("Spencer scale Z (FIXED)")                                                                # FDA delaypars.R: Z = 0.02272
    c1_spencer  <- fixed(9.0)
    label("Spencer O2-asymptote C1Spencer (FIXED)")                                                # FDA delaypars.R: C1Spencer = 9
    c2_spencer  <- fixed(86.11)
    label("Spencer CO2-asymptote C2Spencer (FIXED)")                                               # FDA delaypars.R: C2Spencer = 86.11

    # ===== Brain and tissue physiology =====
    qb0_lmin    <- fixed(0.75)
    label("Baseline cerebral blood flow Qb0 (L/min, FIXED)")                                       # FDA delaypars.R: Qb0 = 0.75/60 L/s = 0.75 L/min
    qt0_lmin    <- fixed(4.25)
    label("Baseline non-brain (peripheral tissue) blood flow Qt0 (L/min, FIXED)")                  # FDA delaypars.R: Qt0 = 4.25/60 L/s = 4.25 L/min
    rou_factor  <- fixed(0.32)
    label("Tissue blood-flow O2-sensitivity coefficient rou (unitless, FIXED)")                    # FDA delaypars.R: rou = 0.32
    v_b         <- fixed(1.32)
    label("Brain blood volume V_B (L, FIXED)")                                                     # FDA delaypars.R: V_B = 1.32
    v_t         <- fixed(38.68)
    label("Peripheral-tissue blood volume V_T (L, FIXED)")                                         # FDA delaypars.R: V_T = 38.68
    m_b_co2_lmin <- fixed(0.04)
    label("Brain CO2 metabolic production rate M_B_co2 at baseline (L/min, FIXED)")                # FDA delaypars.R: M_B_co2 = 0.04/60 L/s = 0.04 L/min
    m_b_o2_0_lmin <- fixed(-0.05)
    label("Brain O2 metabolic consumption rate M_B_o2_0 at baseline (L/min, FIXED, negative = consumption)") # FDA delaypars.R: M_B_o2_0 = -0.05/60 L/s = -0.05 L/min
    m_t_co2_lmin <- fixed(0.16)
    label("Tissue CO2 metabolic production rate M_T_co2 at baseline (L/min, FIXED)")               # FDA delaypars.R: M_T_co2 = 0.16/60 L/s = 0.16 L/min
    m_t_o2_0_lmin <- fixed(-0.2)
    label("Tissue O2 metabolic consumption rate M_T_o2_0 at baseline (L/min, FIXED, negative)")    # FDA delaypars.R: M_T_o2_0 = -0.2/60 L/s = -0.2 L/min
    fl_minfrac  <- fixed(0.87)
    label("Minimum metabolic fraction fL during chronic suppression (FIXED)")                      # FDA delaypars.R: fL = 0.87
    fn_maxfrac  <- fixed(0.93)
    label("Maximum metabolic fraction fN during compensation (FIXED)")                             # FDA delaypars.R: fN = 0.93

    # ===== Chemoreflex drives =====
    aco2_drive  <- fixed(-30.0)
    label("CO2 chemoreflex sigmoid asymptote A_co2 (FIXED)")                                       # FDA delaypars.R: Aco2 = -30
    bco2_drive  <- fixed(60.0)
    label("CO2 chemoreflex sigmoid amplitude B_co2 (FIXED)")                                       # FDA delaypars.R: Bco2 = 60
    cco2_drive  <- fixed(40.0)
    label("CO2 chemoreflex sigmoid midpoint C_co2 (mm Hg, FIXED)")                                 # FDA delaypars.R: Cco2 = 40
    dco2_drive  <- fixed(4.5)
    label("CO2 chemoreflex sigmoid steepness D_co2 (mm Hg, FIXED)")                                # FDA delaypars.R: Dco2 = 4.5
    # NOTE: FDA delaypars.R also defines P_a_co2_0 = 40 mm Hg used as the
    # baseline arterial CO2 partial pressure. The central drive
    # ODE uses brain CO2 baseline P_B_co2_0 instead (delaymymod.c line
    # 518), so the arterial baseline is not separately required in this
    # encoding.
    tau_co2_min <- fixed(0.3333)
    label("CO2 chemoreflex filtering time-constant tau_co2 (min, FIXED, = 20 s)")                 # FDA delaypars.R: tau_co2 = 20 s = 0.333 min
    c1o2_drive  <- fixed(17.0)
    label("O2 chemoreflex exponential prefactor c1o2 (FIXED)")                                     # FDA delaypars.R: c1o2 = 17
    c2o2_drive  <- fixed(11.0)
    label("O2 chemoreflex exponential length scale c2o2 (mm Hg, FIXED)")                           # FDA delaypars.R: c2o2 = 11
    p_a_o2_0    <- fixed(95.0)
    label("Baseline arterial O2 partial pressure P_a_o2_0 (mm Hg, FIXED)")                         # FDA delaypars.R: P_a_o2_0 = 95
    tau_o2_min  <- fixed(0.1667)
    label("O2 chemoreflex filtering time-constant tau_o2 (min, FIXED, = 10 s)")                   # FDA delaypars.R: tau_o2 = 10 s = 0.1667 min

    # ===== Peripheral and central chemoreflex drive parameters =====
    bp_log      <- fixed(18.0)
    label("Peripheral-drive CO2 reference Bp (mm Hg, FIXED)")                                      # FDA delaypars.R: Bp = 18
    k_fpc       <- fixed(1.738)
    label("Peripheral-drive scale K_fpc (FIXED)")                                                  # FDA delaypars.R: K_fpc = 1.738
    f_pc_max    <- fixed(12.3)
    label("Peripheral-drive O2-modulation upper f_pc_max (FIXED)")                                 # FDA delaypars.R: f_pc_max = 12.3
    f_pc_min    <- fixed(0.8352)
    label("Peripheral-drive O2-modulation lower f_pc_min (FIXED)")                                 # FDA delaypars.R: f_pc_min = 0.8352
    f_pc_0      <- fixed(3.67)
    label("Peripheral-drive baseline f_pc_0 (FIXED)")                                              # FDA delaypars.R: f_pc_0 = 3.67
    p_a_o2_c    <- fixed(45.0)
    label("Peripheral-drive O2-modulation midpoint P_a_o2_c (mm Hg, FIXED)")                       # FDA delaypars.R: P_a_o2_c = 45
    k_pc_steep  <- fixed(29.27)
    label("Peripheral-drive O2-modulation steepness K_pc (mm Hg, FIXED)")                          # FDA delaypars.R: K_pc = 29.27
    tau_dp_min  <- fixed(0.1167)
    label("Peripheral-drive filtering time-constant tau_Dp (min, FIXED, = 7 s)")                  # FDA delaypars.R: tau_Dp = 7 s = 0.1167 min
    tau_dc_min  <- fixed(1.0)
    label("Central-drive filtering time-constant tau_Dc (min, FIXED, = 60 s)")                    # FDA delaypars.R: tau_Dc = 60 s = 1 min
    g_dp        <- fixed(2.5)
    label("Peripheral-drive gain G_Dp (FIXED)")                                                    # FDA delaypars.R: G_Dp = 2.5
    g_dc        <- fixed(2.0)
    label("Central-drive gain G_Dc (FIXED)")                                                       # FDA delaypars.R: G_Dc = 2

    # ===== Brain CO2 / O2 baselines for central-drive feedback =====
    p_b_co2_0   <- fixed(45.27)
    label("Baseline brain CO2 partial pressure P_B_co2_0 (mm Hg, FIXED)")                          # FDA delaypars.R: P_B_co2_0 = 45.27
    p_b_o2_0    <- fixed(32.0)
    label("Baseline brain O2 partial pressure P_B_o2_0 (mm Hg, FIXED)")                            # FDA delaypars.R: P_B_o2_0 = 32

    # ===== Central ventilatory depression factor alphaH (Mann 2022 eqs 6-9) =====
    gh          <- fixed(10.0)
    label("Hypoxic ventilatory depression slope G_H (unitless, FIXED)")                            # FDA delaypars.R: Gh = 10; Mann 2022 Supp eq 6-8
    theta_hmin  <- fixed(29.8)
    label("alphaH lower-saturation brain O2 threshold theta_Hmin (mm Hg, FIXED)")                  # FDA delaypars.R: theta_Hmin = 29.8; Mann 2022 Supp eq 6-8
    theta_hmax  <- fixed(37.0)
    label("alphaH upper-saturation brain O2 threshold theta_Hmax (mm Hg, FIXED, Mann updated from Ursino 35)") # FDA delaypars.R: theta_Hmax = 37; Mann 2022 Supp text (paragraph after eq 8) -- updated from Ursino 35 to 37 mm Hg
    tau_h_min   <- fixed(5.0)
    label("alphaH filtering time-constant tau_h (min, FIXED, = 300 s)")                            # FDA delaypars.R: tau_h = 300 s = 5 min

    # ===== Cardiac arrest event parameters (Mann 2022 cardiac-arrest section) =====
    ca_pao2_threshold <- fixed(15.0)
    label("PaO2 threshold for cardiac-arrest trigger (mm Hg, FIXED, Mann 2022 + Hobler 1973)")     # Mann 2022 Supplement 1 'Cardiovascular Collapse' section: PaO2 < 15 mm Hg triggers collapse; threshold inherited from Hobler 1973 [17]
    ca_sustain_min    <- fixed(3.6667)
    label("Sustained-below-threshold time before collapse-decay activates (min, FIXED, = 220 s)")   # Mann 2022 Supplement 1 'Cardiovascular Collapse' section: 220 s sustained delay
    ca_decay_rate     <- fixed(0.024)
    label("Cardiovascular-collapse decay rate Cim25 (1/min after trigger, FIXED)")                  # FDA delaypars.R: Cim25 = 0.004 * 6 = 0.024
  })

  model({
    # ===== Per-patient-type sensitivity exponents =====
    # OPIOID_PATIENT_TYPE = 0 -> opioid-naive; = 1 -> chronic opioid user.
    p1 <- (1 - OPIOID_PATIENT_TYPE) * p1_naive + OPIOID_PATIENT_TYPE * p1_chronic
    p3 <- (1 - OPIOID_PATIENT_TYPE) * p3_naive + OPIOID_PATIENT_TYPE * p3_chronic

    # Opioid-driven reductions on the three ventilatory drives.
    # P3 acts on the wakefulness drive (FDA delaymymod.c line 274: Kf = Wmax*pow(CAR, P3)).
    # P1 acts on both the peripheral and central chemoreflex drives
    # (FDA delaymymod.c lines 511, 517: AV1 = AV2 = 1 - CAR^P1).
    kf  <- wmax_drive * (CAR_OPIOID + 1e-12)^p3
    av1 <- 1 - (CAR_OPIOID + 1e-12)^p1
    av2 <- av1
    # The (CAR_OPIOID + 1e-12)^p3 / ^p1 guards against 0^p3 evaluating
    # to NaN when CAR_OPIOID is identically zero during initial conditions.

    # Residual fraction of wakefulness drive after opioid suppression
    fracw <- (w_baseline - kf) / 6.62
    # Clamp non-negative (matches FDA delaymymod.c lines 276-278)
    fracw_pos <- (fracw > 0) * fracw

    # ===== Metabolic rate scaling by residual wakefulness =====
    # FDA delaymymod.c lines 279-294: M_real = M_baseline * fracw^P2, with
    # lower clamp at M_baseline * fL and upper clamp at M_baseline * fN.
    metab_scale <- (fracw_pos)^p2
    metab_scale_co2 <- (metab_scale < fl_minfrac) * fl_minfrac +
                       (metab_scale >= fl_minfrac) * metab_scale
    metab_scale_o2  <- (metab_scale > fn_maxfrac) * fn_maxfrac +
                       (metab_scale <= fn_maxfrac) * metab_scale
    m_b_co2_real <- m_b_co2_lmin * metab_scale_co2
    m_t_co2_real <- m_t_co2_lmin * metab_scale_co2
    m_b_o2_real  <- m_b_o2_0_lmin * metab_scale_o2
    m_t_o2_real  <- m_t_o2_0_lmin * metab_scale_o2

    # ===== Spencer gas dissociation (forward direction: P -> C) =====
    # Brain partial pressures from blood gas concentrations.
    c_b_co2_diss <- cb_co2 - c_b_hco3
    c_b_co2_pos  <- (c_b_co2_diss > 0) * c_b_co2_diss
    p_b_co2      <- c_b_co2_pos / alpha_co2

    c_b_o2_pos   <- (cb_o2 > 0) * cb_o2
    p_b_o2       <- c_b_o2_pos / alpha_o2

    c_t_co2_diss <- ct_co2 - c_t_hco3
    c_t_co2_pos  <- (c_t_co2_diss > 0) * c_t_co2_diss
    p_t_co2      <- c_t_co2_pos / alpha_co2

    c_t_o2_pos   <- (ct_o2 > 0) * ct_o2
    p_t_o2       <- c_t_o2_pos / alpha_o2

    # Spencer-curve mapping from brain venous P_co2, P_o2 -> C_Vb_co2, C_Vb_o2
    f2_b   <- p_b_co2 * (1 + beta2_sp * p_b_o2) /
              (k2_spencer * (1 + alpha2_sp * p_b_o2))
    f2_b_p <- (f2_b > 0) * f2_b
    c_vb_co2 <- z_sp * c2_spencer * (f2_b_p^(1 / a2_spencer)) /
                (1 + f2_b_p^(1 / a2_spencer))
    f1_b   <- p_b_o2 * (1 + beta1_sp * p_b_co2) /
              (k1_spencer * (1 + alpha1_sp * p_b_co2))
    f1_b_p <- (f1_b > 0) * f1_b
    c_vb_o2  <- z_sp * c1_spencer * (f1_b_p^(1 / a1_spencer)) /
                (1 + f1_b_p^(1 / a1_spencer))

    # Tissue venous concentrations from P_t_co2 and P_t_o2
    f2_t   <- p_t_co2 * (1 + beta2_sp * p_t_o2) /
              (k2_spencer * (1 + alpha2_sp * p_t_o2))
    f2_t_p <- (f2_t > 0) * f2_t
    c_vt_co2 <- z_sp * c2_spencer * (f2_t_p^(1 / a2_spencer)) /
                (1 + f2_t_p^(1 / a2_spencer))
    f1_t   <- p_t_o2 * (1 + beta1_sp * p_t_co2) /
              (k1_spencer * (1 + alpha1_sp * p_t_co2))
    f1_t_p <- (f1_t > 0) * f1_t
    c_vt_o2  <- z_sp * c1_spencer * (f1_t_p^(1 / a1_spencer)) /
                (1 + f1_t_p^(1 / a1_spencer))

    # End-tidal (alveolar) blood-gas concentrations from palv_co2 / palv_o2
    f2_e   <- palv_co2 * (1 + beta2_sp * palv_o2) /
              (k2_spencer * (1 + alpha2_sp * palv_o2))
    f2_e_p <- (f2_e > 0) * f2_e
    c_e_co2 <- z_sp * c2_spencer * (f2_e_p^(1 / a2_spencer)) /
               (1 + f2_e_p^(1 / a2_spencer))
    f1_e   <- palv_o2 * (1 + beta1_sp * palv_co2) /
              (k1_spencer * (1 + alpha1_sp * palv_co2))
    f1_e_p <- (f1_e > 0) * f1_e
    c_e_o2  <- z_sp * c1_spencer * (f1_e_p^(1 / a1_spencer)) /
               (1 + f1_e_p^(1 / a1_spencer))

    # ===== Blood flow control (Magosso / Ursino with Mann cardiac-arrest decay) =====
    # Baseline blood flows are modulated by chemoreflex states yco2, yo2.
    # The Mann cardiac-arrest decay enters as im_arrest (multiplier 0..1).
    qb <- qb0_lmin * (1 + (yco2 + yo2)) * im_arrest
    qt <- qt0_lmin * (1 + rou_factor * yo2) * im_arrest
    q_total <- qb + qt

    # Mixed venous concentrations (weighted by per-territory blood flow)
    qb_safe <- qb + 1e-9
    qt_safe <- qt + 1e-9
    c_v_co2 <- (qb * c_vb_co2 + qt * c_vt_co2) / (qb_safe + qt_safe)
    c_v_o2  <- (qb * c_vb_o2  + qt * c_vt_o2)  / (qb_safe + qt_safe)

    # Arterial concentrations = pulmonary-shunted blend of end-tidal +
    # mixed venous (FDA delaymymod.c lines 181-182)
    c_a_co2 <- (1 - s1_shunt) * c_e_co2 + s1_shunt * c_v_co2
    c_a_o2  <- (1 - s1_shunt) * c_e_o2  + s1_shunt * c_v_o2

    # ===== Spencer reverse mapping (C -> P) for arterial partial pressures =====
    # FDA delaymymod.c lines 183-191 give the closed-form solution.
    asat_co2 <- z_sp * c2_spencer - c_a_co2
    asat_co2_pos <- (asat_co2 > 1e-9) * asat_co2 + (asat_co2 <= 1e-9) * 1e-9
    asat_o2  <- z_sp * c1_spencer - c_a_o2
    asat_o2_pos  <- (asat_o2  > 1e-9) * asat_o2  + (asat_o2  <= 1e-9) * 1e-9
    c_a_co2_pos <- (c_a_co2 > 1e-9) * c_a_co2 + (c_a_co2 <= 1e-9) * 1e-9
    c_a_o2_pos  <- (c_a_o2  > 1e-9) * c_a_o2  + (c_a_o2  <= 1e-9) * 1e-9

    d2 <- k2_spencer * (c_a_co2_pos / asat_co2_pos)^a2_spencer
    d1 <- k1_spencer * (c_a_o2_pos  / asat_o2_pos )^a1_spencer
    s2_rev <- -(d2 + alpha2_sp * d2 * d1) / (beta1_sp + alpha1_sp * beta2_sp * d1)
    s1_rev <- -(d1 + alpha1_sp * d1 * d2) / (beta2_sp + alpha2_sp * beta1_sp * d2)
    r1_rev <- -(1 + beta1_sp * d2 - beta2_sp * d1 - alpha1_sp * alpha2_sp * d1 * d2) /
              (2 * (beta2_sp + alpha2_sp * beta1_sp * d2))
    r2_rev <- -(1 + beta2_sp * d1 - beta1_sp * d2 - alpha2_sp * alpha1_sp * d2 * d1) /
              (2 * (beta1_sp + alpha1_sp * beta2_sp * d1))

    p_a_co2 <- r2_rev + sqrt(r2_rev * r2_rev - s2_rev + 1e-9)
    p_a_o2  <- r1_rev + sqrt(r1_rev * r1_rev - s1_rev + 1e-9)

    # ===== Chemoreflex input states (Magosso / Ursino) =====
    # psai_co2: CO2 sigmoid input to peripheral chemoreflex.
    psai_co2 <- ((aco2_drive + bco2_drive / (1 + exp(-(p_a_co2 - cco2_drive) / dco2_drive))) *
                 1500 / 100 / 1000 + qb0_lmin) / qb0_lmin - 1
    # FDA delaymymod.c line 502: 1500/100/1000/60 = 0.00025; in min-time we
    # drop the /60 because our blood-flow units are L/min.

    # psai_o2: O2 exponential input to peripheral chemoreflex.
    psai_o2_raw <- c1o2_drive * (exp(-p_a_o2 / c2o2_drive) -
                                 exp(-p_a_o2_0 / c2o2_drive))
    psai_o2     <- (psai_o2_raw > 2) * 2 + (psai_o2_raw <= 2) * psai_o2_raw

    # ===== Peripheral chemoreflex drive (Plag_X = X under zero-delay limit) =====
    f_pc_o2_mod_num <- (f_pc_max + f_pc_min * exp((p_a_o2 - p_a_o2_c) / k_pc_steep))
    f_pc_o2_mod_den <- 1 + exp((p_a_o2 - p_a_o2_c) / k_pc_steep)
    p_a_co2_safe <- (p_a_co2 > 1e-6) * p_a_co2 + (p_a_co2 <= 1e-6) * 1e-6
    log_arg <- p_a_co2_safe / bp_log
    log_arg_safe <- (log_arg > 1e-6) * log_arg + (log_arg <= 1e-6) * 1e-6
    plag_f_pc <- k_fpc * log(log_arg_safe) * f_pc_o2_mod_num / f_pc_o2_mod_den

    # ===== alphaH dynamics (Mann 2022 eqs 6-9) =====
    hstat_low  <- 1 + gh * (theta_hmin - p_b_o2_0) / p_b_o2_0
    hstat_mid  <- 1 + gh * (p_b_o2 - p_b_o2_0) / p_b_o2_0
    hstat_high <- 1 + gh * (theta_hmax - p_b_o2_0) / p_b_o2_0
    hstat <- (p_b_o2 <= theta_hmin) * hstat_low +
             (p_b_o2 >  theta_hmin) * (p_b_o2 <= theta_hmax) * hstat_mid +
             (p_b_o2 >  theta_hmax) * hstat_high

    # ===== Total ventilation drive synthesis =====
    # Mann 2022 Supp eq 10-11: V_E = W + max(0, alphaH*D_P + D_C) with
    # alphaH = 1 if D_p < 0 (eq 11 clamp).
    alpha_h_effective <- (dp_state < 0) * 1 + (dp_state >= 0) * alpha_h
    chemoreflex_drive_raw <- alpha_h_effective * dp_state + dc_state
    chemoreflex_drive <- (chemoreflex_drive_raw > 0) * chemoreflex_drive_raw
    # Total ventilation (L/min) = residual wakefulness + chemoreflex
    total_drive <- (w_baseline - kf) + chemoreflex_drive
    venti <- (total_drive > 0) * total_drive

    # ===== Cardiac arrest tripwire dynamics =====
    # t_pao2_below accumulates time PaO2 has spent below the critical
    # threshold (in minutes). When t_pao2_below crosses ca_sustain_min
    # (= 220 s / 60), the collapse decay activates on the im_arrest
    # multiplier on cardiac output.
    pao2_below <- (p_a_o2 < ca_pao2_threshold) * 1.0
    ca_trigger <- (t_pao2_below > ca_sustain_min) * 1.0

    # ===== Gas-exchange ODEs (Magosso / Ursino) =====
    # Alveolar / arterial CO2 and O2 (FDA delaymymod.c lines 597-598)
    d/dt(palv_co2) <- (venti * (p_i_co2 - palv_co2) +
                      lumbda * q_total * (1 - s1_shunt) * (c_v_co2 - c_e_co2)) / v_a
    d/dt(palv_o2)  <- (venti * (p_i_o2_init - palv_o2) +
                      lumbda * q_total * (1 - s1_shunt) * (c_v_o2 - c_e_o2)) / v_a

    # Brain blood-gas concentrations (FDA lines 599-600)
    d/dt(cb_co2) <- (qb * (c_a_co2 - c_vb_co2) + m_b_co2_real) / v_b
    d/dt(cb_o2)  <- (qb * (c_a_o2  - c_vb_o2 ) + m_b_o2_real ) / v_b

    # Tissue blood-gas concentrations (FDA lines 601-602)
    d/dt(ct_co2) <- (qt * (c_a_co2 - c_vt_co2) + m_t_co2_real) / v_t
    d/dt(ct_o2)  <- (qt * (c_a_o2  - c_vt_o2 ) + m_t_o2_real ) / v_t

    # Chemoreflex input filters (FDA lines 604-605, time-constants in min)
    d/dt(yco2) <- (psai_co2 - yco2) / tau_co2_min
    d/dt(yo2)  <- (psai_o2  - yo2)  / tau_o2_min

    # Peripheral and central chemoreflex drives, opioid-attenuated
    # (FDA lines 514-518 + 606-607, zero-delay limit)
    d/dt(dp_state) <- (-dp_state + av1 * g_dp * (plag_f_pc - f_pc_0)) / tau_dp_min
    d/dt(dc_state) <- (-dc_state + av2 * g_dc * (p_b_co2 - p_b_co2_0)) / tau_dc_min

    # alphaH dynamics (Mann 2022 eq 9, in min)
    d/dt(alpha_h) <- (hstat - alpha_h) / tau_h_min

    # Time-below-threshold accumulator for cardiac arrest trigger.
    # pao2_below evaluates 1 while P_a_o2 < threshold, 0 otherwise.
    # The state grows at 1 per min while below threshold and resets
    # slowly while above (so brief excursions above threshold do not
    # reset the trigger; if a sustained recovery occurs, the slow
    # negative drift eventually re-arms).
    d/dt(t_pao2_below) <- pao2_below - (1 - pao2_below) * t_pao2_below * 0.1

    # Cardiac-arrest decay multiplier on cardiac output. ca_trigger is 1
    # once t_pao2_below has accumulated past ca_sustain_min; im_arrest
    # then decays toward 0 at rate ca_decay_rate (1/min), driving
    # cardiac output (qb + qt) toward the Mann 2022 cardiac-arrest
    # floor of 0.01 L/min.
    d/dt(im_arrest) <- -ca_trigger * ca_decay_rate * im_arrest

    # ===== Initial conditions =====
    palv_co2(0)      <- 40.28
    palv_o2(0)       <- 100.2
    cb_co2(0)      <- 0.645
    cb_o2(0)       <- 9.78e-4
    ct_co2(0)      <- 0.605
    ct_o2(0)       <- 13e-4
    yco2(0)         <- 0
    yo2(0)          <- 0
    dp_state(0)     <- 0
    dc_state(0)     <- 0
    alpha_h(0)      <- 1
    t_pao2_below(0) <- 0
    im_arrest(0)    <- 1

    # ===== Outputs =====
    Venti         <- venti                 # minute ventilation (L/min)
    P_a_O2        <- p_a_o2                # arterial O2 partial pressure (mm Hg)
    P_a_CO2       <- p_a_co2               # arterial CO2 partial pressure (mm Hg)
    Q_brain       <- qb                    # cerebral blood flow (L/min)
    Q_tissue      <- qt                    # peripheral-tissue blood flow (L/min)
    Q_total       <- q_total               # total cardiac output (L/min)
    P_B_O2        <- p_b_o2                # brain O2 partial pressure (mm Hg)

    # Arterial oxygen saturation (Severinghaus-style approximation as
    # used in FDA delaymymod.c lines 555-559)
    nSat <- 2.6
    k3_sat <- 26.6
    p_a_o2_safe <- (p_a_o2 > 1e-6) * p_a_o2 + (p_a_o2 <= 1e-6) * 1e-6
    po2_virt <- p_a_o2_safe * (40 / p_a_o2_safe)^0.3
    SAT_O2 <- 100 * po2_virt^nSat / (k3_sat^nSat + po2_virt^nSat)

    CA_active     <- ca_trigger            # 1 = cardiac-arrest decay has been activated
    T_pao2_below  <- t_pao2_below          # accumulated time PaO2 spent below critical threshold (min)
  })
}
