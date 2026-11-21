Chakravartula_2024_fentanyl_pediatric <- function() {
  description <- paste(
    "QSP (pediatric). Age-adjusted extension of the Mann 2022 respiratory /",
    "cerebrovascular physiology model (Mann_2022_respiratory_physiology) to",
    "children aged 2-12 years, from Chakravartula 2024. The receptor-binding",
    "and ventilatory-depression pharmacodynamics are unchanged from the adult",
    "model; the physiological component is re-parameterized by age for body,",
    "brain and tissue volumes, baseline cerebral and tissue blood flows,",
    "oxygen / carbon-dioxide metabolic rates, wakefulness and chemoreflex",
    "set-points (paper Tables S1-S4), and the cerebrovascular reactivity",
    "(CBF response to hypercapnia and hypoxia) is reduced in children.",
    "Structurally this is the Mann model plus one extra state: a brain-",
    "specific O2 cerebral-blood-flow response (y_child_o2) with child",
    "parameters, while the peripheral-tissue O2 response keeps adult",
    "parameters. With the age parameters set to the paper's adult values the",
    "model reduces exactly to Mann_2022_respiratory_physiology. Parameter",
    "defaults here are the 10-year-old child (the paper's best-validated",
    "pediatric age; Figures 5-6). Driven by CAR_OPIOID (fraction of mu-opioid",
    "receptors bound), it reproduces the paper's finding that children",
    "desaturate (arterial and brain PO2) faster than adults at the same",
    "plasma fentanyl concentration. Composes downstream of",
    "Mann_2022_mu_receptor_binding."
  )
  reference <- paste(
    "Chakravartula S, Thrasher B, Mann J, Chaturbedi A, Han X, Dahan A,",
    "Florian J, Strauss D, Li Z. Physiologically based modeling reveals",
    "different risk of respiratory depression after fentanyl overdose",
    "between adults and children. Clin Transl Sci. 2024;17(4):e13780.",
    "doi:10.1111/cts.13780. PMID 38618722.",
    "Adult structural respiratory / cerebrovascular model and FDA reference",
    "implementation inherited from Mann J et al. Clin Pharmacol Ther.",
    "2022;112(5):1020-1032 (doi:10.1002/cpt.2696; Supplement 1 'Physiological",
    "Component'), itself building on Magosso E, Ursino M, van Oostrom JH.",
    "IEEE Trans Biomed Eng. 2004;51(7):1115-1128. Pediatric parameter",
    "sources cited by Chakravartula 2024 Supplement 1: brain / cardiac blood",
    "flow Wu C et al. J Am Heart Assoc. 2016;5(1); brain metabolism Kuzawa",
    "CW et al. PNAS. 2014;111(36):13010-5; cerebrovascular reactivity Tallon",
    "CM et al. Exp Physiol. 2020;105(7):1090-1101."
  )
  vignette <- "Chakravartula_2024_fentanyl_pediatric_respiratory_depression"
  units <- list(
    time = "min",
    dosing = "(not applicable; CAR_OPIOID is a time-varying covariate input from the binding layer)",
    concentration = "L/min (primary output Venti is minute ventilation; the model has no plasma-concentration analyte and the canonical Cc observation is not used)"
  )

  compartmentData <- list(
    palv_co2 = list(
      analyte = "CO2 partial pressure in alveoli",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    palv_o2 = list(
      analyte = "O2 partial pressure in alveoli",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    cb_co2 = list(
      analyte = "CO2 partial pressure in brain capillaries",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    cb_o2 = list(
      analyte = "O2 partial pressure in brain capillaries",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    ct_co2 = list(
      analyte = "CO2 partial pressure in tissues",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    ct_o2 = list(
      analyte = "O2 partial pressure in tissues",
      units = NA_character_,
      specimen = "blood cell",
      verified = FALSE
    ),
    yco2 = list(
      analyte = "CBF response to CO2 (y)",
      units = NA_character_,
      specimen = "administration site",
      verified = FALSE
    ),
    yo2 = list(
      analyte = "tissue CBF response to O2 (y)",
      units = NA_character_,
      specimen = "administration site",
      verified = FALSE
    ),
    y_child_o2 = list(
      analyte = "brain CBF response to O2 (child, y)",
      units = NA_character_,
      specimen = "administration site",
      verified = FALSE
    ),
    dp_state = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    dc_state = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    alpha_h = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    crossed_latch = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    first_cross_time = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    check_complete = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    ca_arrest_latch = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE),
    im_arrest = list(analyte = "none", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    CAR_OPIOID = list(
      description = "Time-varying fraction of mu-opioid receptors bound by an opioid agonist (RL_op output of Mann_2022_mu_receptor_binding); drives the ventilatory-depression reductions of wakefulness and chemoreflex drives. The binding / pharmacodynamic relationship is unchanged from the adult Mann 2022 model (Chakravartula 2024 assumes the same CAR-to-drive mapping in children).",
      units = "fraction",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Range 0..1. At CAR_OPIOID = 0 the model is at eupnoeic steady",
        "state for the configured age; at CAR_OPIOID = 1 wakefulness drive",
        "is abolished and the chemoreflex drives are fully suppressed."
      ),
      source_name = "(inherited from Mann_2022_respiratory_physiology)"
    ),
    OPIOID_PATIENT_TYPE = list(
      description = "Categorical 0 = opioid-naive, 1 = chronic opioid user; selects the P1 (chemoreflex-drive reduction exponent) and P3 (wakefulness-drive reduction exponent) sensitivity parameters. Inherited unchanged from Mann 2022; pediatric overdose simulations use the naive setting.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (opioid-naive)",
      notes = paste(
        "Mann 2022 naive (P1 = 2.875, P3 = 0.9) and chronic (P1 = 4.226,",
        "P3 = 1.323) parameter sets; shared metabolism exponent P2 = 0.06319."
      ),
      source_name = "(inherited from Mann_2022_respiratory_physiology)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "2-12 years (defaults are the 10-year-old child); reduces to the adult Mann model at adult parameter values",
    weight_range = "Table S1 typical pediatric WT 12.5-42.0 kg (ages 2-12); 34.0 kg at age 10",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state = paste(
      "Age-adjusted ventilatory and cerebrovascular physiology calibrated to",
      "pediatric clinical steady-state minute ventilation, arterial blood-gas",
      "tensions, and cerebral-blood-flow reactivity to hypercapnia; used to",
      "study opioid-induced respiratory depression and hypoxia in children",
      "versus adults."
    ),
    dose_range = NA_character_,
    regions = NA_character_,
    notes = paste(
      "Seventeen physiological + event states are integrated. The eleven",
      "Mann 2022 physiology states (palv_co2, palv_o2, cb_co2, cb_o2,",
      "ct_co2, ct_o2, yco2, yo2, dp_state, dc_state, alpha_h) are carried",
      "unchanged in STRUCTURE; one pediatric state is added, y_child_o2 -",
      "the brain-specific cerebral-blood-flow response to O2 with child",
      "parameters (paper Supplement eqs 3-4). Brain blood flow is modulated",
      "by the CO2 response (yco2) plus the brain O2 response (y_child_o2);",
      "peripheral-tissue blood flow keeps the adult O2 response (yo2, eqs",
      "5-6). Five cardiac-arrest event states (crossed_latch,",
      "first_cross_time, check_complete, ca_arrest_latch, im_arrest) are",
      "inherited verbatim from Mann_2022_respiratory_physiology.",
      "AGE PARAMETERIZATION. The parameter defaults below are the",
      "10-year-old values from paper Tables S1-S4. To simulate another age",
      "(2-12, or adult), override the age-varying parameters via the params",
      "argument of rxode2::rxSolve. The age-varying parameters are: qb0_lmin,",
      "qt0_lmin (Qb0, Qt0; Table S1); v_a, v_b, v_t (lung / brain / other-",
      "tissue volumes; Table S1); brain_weight_g (Table S1); m_b_co2_lmin,",
      "m_b_o2_0_lmin, m_t_co2_lmin, m_t_o2_0_lmin (Table S2); wmax_drive,",
      "w_baseline (W; Table S4), p_b_co2_0 (Table S4), f_pc_0 (Table S4);",
      "aco2_drive, bco2_drive, tau_co2_min (A_CO2, B_CO2, tau_CO2; Table S3,",
      "shared across ages 2-12); c1o2_child, tau_child_o2_min (c1,child,",
      "tau,child,O2; Table S3). Setting all of these to their Table adult-row",
      "values recovers Mann_2022_respiratory_physiology exactly (then",
      "y_child_o2 collapses onto yo2 and brain_weight_g = 1500).",
      "The paper assumes ages 2-9 and 11-12 share the 10-year-old",
      "cerebrovascular-reactivity parameters (Table S3) because dynamic",
      "hypercapnia CBF data exist only for ~10-year-olds (Tallon et al.).",
      "PRE-EQUILIBRATION. As for the adult Mann model, the ini() initial",
      "conditions are the FDA delaystates.R nominal values, NOT the",
      "age-specific pre-dose equilibrium; downstream overdose simulations",
      "MUST first run a no-drug pre-equilibration (helper Mann2022Equilibrate)",
      "and pass the result as the inits argument to rxode2::rxSolve, or the",
      "chemoreflex transient biases the overdose trajectory.",
      "DEVIATIONS inherited from Mann_2022_respiratory_physiology (zero-delay",
      "chemoreflex via lagState, single-im cardiac-arrest decay, one-shot",
      "mtime arrest trigger, P1 / P3 code-vs-text mapping) are documented in",
      "that model's population notes and apply unchanged here."
    )
  )

  ini({
    # ===== Pharmacodynamic sensitivity exponents (unchanged from Mann 2022,
    # FDA fentanyl_pars.txt + simulateToGetOD_IM.R; Chakravartula 2024 keeps
    # the adult CAR-to-drive mapping in children) =====
    p1_naive    <- fixed(2.875)
    label("Chemoreflex-drive sensitivity exponent P1 in opioid-naive subjects (unitless)")  # Mann 2022 / FDA simulateToGetOD_IM.R: 2.5 * 1.15 = 2.875
    p1_chronic  <- fixed(4.22625)
    label("Chemoreflex-drive sensitivity exponent P1 in chronic opioid users (unitless)")    # Mann 2022 / FDA: 2.5 * 1.47 * 1.15 = 4.22625
    p3_naive    <- fixed(0.9)
    label("Wakefulness-drive sensitivity exponent P3 in opioid-naive subjects (unitless)")   # Mann 2022 / FDA: 0.9
    p3_chronic  <- fixed(1.323)
    label("Wakefulness-drive sensitivity exponent P3 in chronic opioid users (unitless)")    # Mann 2022 / FDA: 0.9 * 1.47 = 1.323
    p2          <- fixed(0.06319)
    label("Metabolism scaling exponent P2 (unitless,, both patient types)")                   # Mann 2022 / FDA: 0.06319

    # ===== Wakefulness drive (age-varying W; paper Table S4) =====
    wmax_drive  <- fixed(6.142)
    label("Maximum wakefulness drive Wmax (L/min)")                                          # Table S4 age 10: W = 6.142 (adult 6.620)
    w_baseline  <- fixed(6.142)
    label("Baseline wakefulness drive W (L/min)")                                            # Table S4 age 10: W = 6.142 (adult 6.620)

    # ===== Alveolar / gas-exchange parameters (v_a age-varying; Table S1) =====
    v_a         <- fixed(1.691)
    label("Alveolar / lung volume V_A (L)")                                                  # Table S1 age 10: lung volume = 1.691 L (adult 3.28)
    p_i_co2     <- fixed(0.0)
    label("Inspired CO2 partial pressure P_I_co2 (mm Hg, room-air = 0)")                     # Mann 2022 / FDA delaypars.R: P_I_co2 = 0
    p_i_o2_init <- fixed(149.0)
    label("Inspired O2 partial pressure P_I_o2 (mm Hg, room-air sea-level)")                 # Mann 2022 / FDA delaystates.R: P_I_o2 = 149
    lumbda      <- fixed(863.0)
    label("Blood-to-gas partition coefficient lambda (unitless,, Henry's law family)")        # Mann 2022 / FDA delaypars.R: lumbda = 863
    s1_shunt    <- fixed(0.024)
    label("Pulmonary venous-admixture shunt fraction s1 (unitless)")                        # Mann 2022 / FDA delaypars.R: s1 = 0.024
    b_max       <- fixed(0.66)
    label("Alveolar-fraction (dead-space-corrected) Bmax (unitless)")                       # Mann 2022 / FDA delaypars.R: Bmax = 0.66

    # ===== Blood-gas (Spencer) dissociation curve constants (unchanged) =====
    alpha_co2   <- fixed(6.67e-4)
    label("CO2 plasma solubility coefficient alpha_co2 (Spencer)")                          # Mann 2022 / FDA delaypars.R: alpha_co2 = 6.67e-4
    alpha_o2    <- fixed(3.17e-5)
    label("O2 plasma solubility coefficient alpha_o2 (Spencer)")                            # Mann 2022 / FDA delaypars.R: alpha_o2 = 3.17e-5
    c_b_hco3    <- fixed(0.5824)
    label("Brain bicarbonate buffer concentration C_B_hco3 (Spencer)")                      # Mann 2022 / FDA delaypars.R: C_B_hco3 = 0.5824
    c_t_hco3    <- fixed(0.5824)
    label("Tissue bicarbonate buffer concentration C_T_hco3 (Spencer)")                     # Mann 2022 / FDA delaypars.R: C_T_hco3 = 0.5824
    k1_spencer  <- fixed(14.99)
    label("Spencer constant K1 for O2 dissociation")                                       # Mann 2022 / FDA delaypars.R: K1 = 14.99
    k2_spencer  <- fixed(194.4)
    label("Spencer constant K2 for CO2 dissociation")                                      # Mann 2022 / FDA delaypars.R: K2 = 194.4
    a1_spencer  <- fixed(0.3836)
    label("Spencer exponent a1 for O2 dissociation")                                       # Mann 2022 / FDA delaypars.R: a1 = 0.3836
    a2_spencer  <- fixed(1.819)
    label("Spencer exponent a2 for CO2 dissociation")                                      # Mann 2022 / FDA delaypars.R: a2 = 1.819
    alpha1_sp   <- fixed(0.03198)
    label("Spencer cross-term alpha1 (CO2 effect on O2 binding)")                           # Mann 2022 / FDA delaypars.R: alpha1 = 0.03198
    beta1_sp    <- fixed(0.008275)
    label("Spencer cross-term beta1 (CO2 effect on O2 binding)")                            # Mann 2022 / FDA delaypars.R: beta1 = 0.008275
    alpha2_sp   <- fixed(0.05591)
    label("Spencer cross-term alpha2 (O2 effect on CO2 binding)")                           # Mann 2022 / FDA delaypars.R: alpha2 = 0.05591
    beta2_sp    <- fixed(0.03255)
    label("Spencer cross-term beta2 (O2 effect on CO2 binding)")                            # Mann 2022 / FDA delaypars.R: beta2 = 0.03255
    z_sp        <- fixed(0.02272)
    label("Spencer scale Z")                                                                # Mann 2022 / FDA delaypars.R: Z = 0.02272
    c1_spencer  <- fixed(9.0)
    label("Spencer O2-asymptote C1Spencer")                                                # Mann 2022 / FDA delaypars.R: C1Spencer = 9
    c2_spencer  <- fixed(86.11)
    label("Spencer CO2-asymptote C2Spencer")                                               # Mann 2022 / FDA delaypars.R: C2Spencer = 86.11

    # ===== Brain and tissue physiology (age-varying; Tables S1, S2) =====
    qb0_lmin    <- fixed(0.900)
    label("Baseline cerebral blood flow Qb0 (L/min)")                                       # Table S1 age 10: Qb0 = 0.900 (adult 0.750)
    qt0_lmin    <- fixed(3.700)
    label("Baseline non-brain (peripheral tissue) blood flow Qt0 (L/min)")                  # Table S1 age 10: Qt0 = 3.700 (adult 4.250)
    rou_factor  <- fixed(0.32)
    label("Tissue blood-flow O2-sensitivity coefficient rou (unitless)")                    # Mann 2022 / FDA delaypars.R: rou = 0.32
    v_b         <- fixed(1.205)
    label("Brain volume V_B (L)")                                                           # Table S1 age 10: brain volume = 1.205 L (adult 1.320)
    v_t         <- fixed(19.165)
    label("Peripheral-tissue (other tissue) volume V_T (L)")                                # Table S1 age 10: other tissue volume = 19.165 L (adult 38.680)
    m_b_co2_lmin <- fixed(0.065)
    label("Brain CO2 metabolic production rate M_B_co2 at baseline (L/min)")                # Table S2 age 10: M_B_co2 = 0.065 (adult 0.040)
    m_b_o2_0_lmin <- fixed(-0.081)
    label("Brain O2 metabolic consumption rate M_B_o2_0 at baseline (L/min,, negative = consumption)") # Table S2 age 10: M_B_o2_0 = 0.081 consumption (adult 0.050)
    m_t_co2_lmin <- fixed(0.107)
    label("Tissue CO2 metabolic production rate M_T_co2 at baseline (L/min)")               # Table S2 age 10: M_T_co2 = 0.107 (adult 0.160)
    m_t_o2_0_lmin <- fixed(-0.134)
    label("Tissue O2 metabolic consumption rate M_T_o2_0 at baseline (L/min,, negative)")    # Table S2 age 10: M_T_o2_0 = 0.134 consumption (adult 0.200)
    fl_minfrac  <- fixed(0.87)
    label("Minimum metabolic fraction fL during chronic suppression")                      # Mann 2022 / FDA delaypars.R: fL = 0.87
    fn_maxfrac  <- fixed(0.93)
    label("Maximum metabolic fraction fN during compensation")                             # Mann 2022 / FDA delaypars.R: fN = 0.93
    brain_weight_g <- fixed(1318)
    label("Brain weight scaling the CBF-CO2 response amplitude (g)")                        # Table S1 age 10: brain weight = 1318 g (adult 1500); Mann hardcodes 1500

    # ===== Chemoreflex drives (CO2 CBF-reactivity age-varying; Table S3) =====
    aco2_drive  <- fixed(-46.3)
    label("CO2 chemoreflex / CBF sigmoid asymptote A_co2 (ml/min/100ml)")                   # Table S3 ages 2-12: A_CO2 = -46.3 (adult -30.0)
    bco2_drive  <- fixed(92.6)
    label("CO2 chemoreflex / CBF sigmoid amplitude B_co2 (ml/min/100ml)")                   # Table S3 ages 2-12: B_CO2 = 92.6 (adult 60.0)
    cco2_drive  <- fixed(40.0)
    label("CO2 chemoreflex sigmoid midpoint C_co2 (mm Hg)")                                 # Mann 2022 / FDA delaypars.R: Cco2 = 40 (paper C_CO2 unchanged by age)
    dco2_drive  <- fixed(4.5)
    label("CO2 chemoreflex sigmoid steepness D_co2 (mm Hg)")                                # Mann 2022 / FDA delaypars.R: Dco2 = 4.5 (paper D_CO2 unchanged by age)
    tau_co2_min <- fixed(0.75833)
    label("CO2 CBF-response time constant tau_co2 (min,, = 45.5 s)")                         # Table S3 ages 2-12: tau_CO2 = 45.5 s (adult 20 s)
    c1o2_drive  <- fixed(17.0)
    label("Tissue O2 CBF exponential prefactor c1o2 (adult, shared)")                       # Mann 2022 / FDA delaypars.R: c1o2 = 17 (tissue uses adult params at all ages)
    c2o2_drive  <- fixed(11.0)
    label("O2 CBF exponential length scale c2o2 (mm Hg, shared)")                            # Mann 2022 / FDA delaypars.R: c2o2 = 11
    p_a_o2_0    <- fixed(95.0)
    label("Baseline arterial O2 partial pressure P_a_o2_0 (mm Hg)")                         # Mann 2022 / FDA delaypars.R: P_a_o2_0 = 95
    tau_o2_min  <- fixed(0.1667)
    label("Tissue O2 CBF-response time constant tau_o2 (min,, = 10 s, adult)")               # Mann 2022 / FDA delaypars.R: tau_o2 = 10 s (tissue uses adult params)

    # ===== Brain O2 CBF response, child-specific (paper Supplement eqs 3-4) =====
    c1o2_child  <- fixed(15.32)
    label("Brain O2 CBF exponential prefactor c1,child (ml/min/100ml family)")              # Table S3 ages 2-12: c1,child = 15.32 (adult brain uses tissue c1 = 17)
    tau_child_o2_min <- fixed(2.3983)
    label("Brain O2 CBF-response time constant tau,child,O2 (min,, = 143.9 s)")              # Table S3 ages 2-12: tau,child,O2 = 143.9 s (adult brain uses tissue tau = 10 s)

    # ===== Peripheral and central chemoreflex drive parameters =====
    bp_log      <- fixed(18.0)
    label("Peripheral-drive CO2 reference Bp (mm Hg)")                                      # Mann 2022 / FDA delaypars.R: Bp = 18
    k_fpc       <- fixed(1.738)
    label("Peripheral-drive scale K_fpc")                                                  # Mann 2022 / FDA delaypars.R: K_fpc = 1.738
    f_pc_max    <- fixed(12.3)
    label("Peripheral-drive O2-modulation upper f_pc_max")                                 # Mann 2022 / FDA delaypars.R: f_pc_max = 12.3
    f_pc_min    <- fixed(0.8352)
    label("Peripheral-drive O2-modulation lower f_pc_min")                                 # Mann 2022 / FDA delaypars.R: f_pc_min = 0.8352
    f_pc_0      <- fixed(3.246)
    label("Peripheral-drive baseline f_pc_0 (age-varying)")                                 # Table S4 age 10: f_pc_0 = 3.246 (adult 3.670)
    p_a_o2_c    <- fixed(45.0)
    label("Peripheral-drive O2-modulation midpoint P_a_o2_c (mm Hg)")                       # Mann 2022 / FDA delaypars.R: P_a_o2_c = 45
    k_pc_steep  <- fixed(29.27)
    label("Peripheral-drive O2-modulation steepness K_pc (mm Hg)")                          # Mann 2022 / FDA delaypars.R: K_pc = 29.27
    tau_dp_min  <- fixed(0.1167)
    label("Peripheral-drive filtering time-constant tau_Dp (min,, = 7 s)")                  # Mann 2022 / FDA delaypars.R: tau_Dp = 7 s
    tau_dc_min  <- fixed(1.0)
    label("Central-drive filtering time-constant tau_Dc (min,, = 60 s)")                    # Mann 2022 / FDA delaypars.R: tau_Dc = 60 s
    g_dp        <- fixed(2.5)
    label("Peripheral-drive gain G_Dp")                                                    # Mann 2022 / FDA delaypars.R: G_Dp = 2.5
    g_dc        <- fixed(2.0)
    label("Central-drive gain G_Dc")                                                       # Mann 2022 / FDA delaypars.R: G_Dc = 2

    # ===== Chemoreflex transport-delay durations (unchanged from Mann 2022) =====
    peripheral_delay_min <- fixed(7 / 60)
    label("Peripheral chemoreflex transport delay (min,, = 7 s / 60)")                       # Mann 2022 / FDA: K_Dp = 0.588, Q_baseline => 7 s
    central_delay_min    <- fixed(11 / 60)
    label("Central chemoreflex transport delay (min,, = 11 s / 60)")                         # Mann 2022 / FDA: K_Dc = 0.9239, Q_baseline => 11 s

    # ===== Brain CO2 / O2 baselines for central-drive feedback (p_b_co2_0 age-varying) =====
    p_b_co2_0   <- fixed(44.37)
    label("Baseline brain CO2 partial pressure P_B_co2_0 (mm Hg, age-varying)")             # Table S4 age 10: P_B_co2_0 = 44.37 (adult 45.27)
    p_b_o2_0    <- fixed(32.0)
    label("Baseline brain O2 partial pressure P_B_o2_0 (mm Hg)")                            # Mann 2022 / FDA delaypars.R: P_B_o2_0 = 32 (reference set-point, not age-tabulated)

    # ===== Central ventilatory depression factor alphaH (unchanged from Mann 2022) =====
    gh          <- fixed(10.0)
    label("Hypoxic ventilatory depression slope G_H (unitless)")                            # Mann 2022 / FDA delaypars.R: Gh = 10
    theta_hmin  <- fixed(29.8)
    label("alphaH lower-saturation brain O2 threshold theta_Hmin (mm Hg)")                  # Mann 2022 / FDA delaypars.R: theta_Hmin = 29.8
    theta_hmax  <- fixed(37.0)
    label("alphaH upper-saturation brain O2 threshold theta_Hmax (mm Hg)")                  # Mann 2022 / FDA delaypars.R: theta_Hmax = 37
    tau_h_min   <- fixed(5.0)
    label("alphaH filtering time-constant tau_h (min,, = 300 s)")                            # Mann 2022 / FDA delaypars.R: tau_h = 300 s

    # ===== Cardiac arrest event parameters (unchanged from Mann 2022) =====
    ca_pao2_threshold <- fixed(15.0)
    label("PaO2 threshold for cardiac-arrest trigger (mm Hg)")                              # Mann 2022 Supplement 1 'Cardiovascular Collapse': PaO2 < 15 mm Hg
    ca_sustain_min    <- fixed(3.6667)
    label("Sustained-below-threshold time before collapse-decay activates (min,, = 220 s)")   # Mann 2022 Supplement 1: 220 s sustained delay
    ca_decay_rate     <- fixed(0.004 * 6 * 60)
    label("Cardiovascular-collapse decay rate Cim25 (1/min after trigger)")                  # Mann 2022 / FDA delaypars.R: Cim25 = 0.024 /s * 60 = 1.44 /min
  })

  model({
    # ===== Per-patient-type sensitivity exponents =====
    p1 <- (1 - OPIOID_PATIENT_TYPE) * p1_naive + OPIOID_PATIENT_TYPE * p1_chronic
    p3 <- (1 - OPIOID_PATIENT_TYPE) * p3_naive + OPIOID_PATIENT_TYPE * p3_chronic

    # Opioid-driven reductions on the ventilatory drives (CAR-to-drive
    # mapping unchanged from Mann 2022; Chakravartula 2024 assumes the same
    # relationship in children).
    kf  <- wmax_drive * (CAR_OPIOID + 1e-12)^p3
    av1 <- 1 - (CAR_OPIOID + 1e-12)^p1
    av2 <- av1

    # Residual fraction of wakefulness drive after opioid suppression.
    # Normalized by the age-specific Wmax (= wmax_drive) so the fraction is
    # unit-consistent across ages; at adult Wmax = 6.62 this matches Mann.
    fracw <- (w_baseline - kf) / wmax_drive
    fracw_pos <- (fracw > 0) * fracw

    # ===== Metabolic rate scaling by residual wakefulness (unchanged) =====
    metab_scale <- (fracw_pos)^p2
    metab_scale_co2 <- (metab_scale < fl_minfrac) * fl_minfrac +
                       (metab_scale >= fl_minfrac) * metab_scale
    metab_scale_o2  <- (metab_scale < fn_maxfrac) * fn_maxfrac +
                       (metab_scale >= fn_maxfrac) * metab_scale
    m_b_co2_real <- m_b_co2_lmin * metab_scale_co2
    m_t_co2_real <- m_t_co2_lmin * metab_scale_co2
    m_b_o2_real  <- m_b_o2_0_lmin * metab_scale_o2
    m_t_o2_real  <- m_t_o2_0_lmin * metab_scale_o2

    # ===== Spencer gas dissociation (forward direction: P -> C) =====
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

    # ===== Low-PaO2 metabolism safety clamps (unchanged) =====
    m_o2_safety_brain  <- (p_b_o2 < 6) * (p_b_o2 / 6) + (p_b_o2 >= 6) * 1.0
    m_o2_safety_tissue <- (p_t_o2 < 6) * (p_t_o2 / 6) + (p_t_o2 >= 6) * 1.0

    # Spencer-curve mapping brain venous P -> C
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

    # Tissue venous concentrations
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

    # End-tidal (alveolar) concentrations
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

    # ===== Blood flow control (PEDIATRIC: brain uses yco2 + y_child_o2) =====
    # Chakravartula 2024 Supplement: brain CBF responds to CO2 (yco2) and to
    # O2 through a child-specific response (y_child_o2, eqs 3-4); peripheral
    # tissue CBF keeps the adult O2 response (yo2, eqs 5-6). At adult
    # parameter values y_child_o2 collapses onto yo2 and this reduces to the
    # Mann adult qb = qb0 * (1 + yco2 + yo2).
    qb <- qb0_lmin * (1 + (yco2 + y_child_o2)) * im_arrest
    qt <- qt0_lmin * (1 + rou_factor * yo2) * im_arrest
    q_total <- qb + qt

    # Mixed venous concentrations
    qb_safe <- qb + 1e-9
    qt_safe <- qt + 1e-9
    c_v_co2 <- (qb * c_vb_co2 + qt * c_vt_co2) / (qb_safe + qt_safe)
    c_v_o2  <- (qb * c_vb_o2  + qt * c_vt_o2)  / (qb_safe + qt_safe)

    # Arterial concentrations
    c_a_co2 <- (1 - s1_shunt) * c_e_co2 + s1_shunt * c_v_co2
    c_a_o2  <- (1 - s1_shunt) * c_e_o2  + s1_shunt * c_v_o2

    # ===== Spencer reverse mapping (C -> P) for arterial partial pressures =====
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

    # ===== CBF response inputs (Chakravartula 2024 Supplement eqs 2, 4, 6) =====
    # CO2 CBF response amplitude scales with brain weight (Mann hardcodes the
    # adult 1500 g; here brain_weight_g is age-specific). psai_co2 is the
    # fractional CBF change at the current arterial CO2.
    psai_co2 <- ((aco2_drive + bco2_drive / (1 + exp(-(p_a_co2 - cco2_drive) / dco2_drive))) *
                 brain_weight_g / 100 / 1000 + qb0_lmin) / qb0_lmin - 1

    # Tissue O2 CBF response (adult parameters; eq 6), capped as in Mann.
    psai_o2_raw <- c1o2_drive * (exp(-p_a_o2 / c2o2_drive) -
                                 exp(-p_a_o2_0 / c2o2_drive))
    psai_o2     <- (psai_o2_raw > 2) * 2 + (psai_o2_raw <= 2) * psai_o2_raw

    # Brain O2 CBF response (child parameters; eq 4), capped likewise. At
    # adult parameters c1o2_child = c1o2_drive so this equals psai_o2.
    psai_o2_child_raw <- c1o2_child * (exp(-p_a_o2 / c2o2_drive) -
                                       exp(-p_a_o2_0 / c2o2_drive))
    psai_o2_child <- (psai_o2_child_raw > 2) * 2 + (psai_o2_child_raw <= 2) * psai_o2_child_raw

    # ===== Chemoreflex transport-delay lookups via lagState (true delay) =====
    plag_p_a_o2  <- lagState(t, p_a_o2,  peripheral_delay_min, 0, p_a_o2_0)
    plag_p_a_co2 <- lagState(t, p_a_co2, peripheral_delay_min, 1, 40)
    clag_p_b_co2 <- lagState(t, p_b_co2, central_delay_min,    2, p_b_co2_0)

    # ===== Peripheral chemoreflex drive (uses lagged arterial gases) =====
    f_pc_o2_mod_num <- (f_pc_max + f_pc_min * exp((plag_p_a_o2 - p_a_o2_c) / k_pc_steep))
    f_pc_o2_mod_den <- 1 + exp((plag_p_a_o2 - p_a_o2_c) / k_pc_steep)
    plag_p_a_co2_safe <- (plag_p_a_co2 > 1e-6) * plag_p_a_co2 + (plag_p_a_co2 <= 1e-6) * 1e-6
    log_arg <- plag_p_a_co2_safe / bp_log
    log_arg_safe <- (log_arg > 1e-6) * log_arg + (log_arg <= 1e-6) * 1e-6
    plag_f_pc <- k_fpc * log(log_arg_safe) * f_pc_o2_mod_num / f_pc_o2_mod_den

    # ===== alphaH dynamics (unchanged from Mann 2022) =====
    hstat_low  <- 1 + gh * (theta_hmin - p_b_o2_0) / p_b_o2_0
    hstat_mid  <- 1 + gh * (p_b_o2 - p_b_o2_0) / p_b_o2_0
    hstat_high <- 1 + gh * (theta_hmax - p_b_o2_0) / p_b_o2_0
    hstat <- (p_b_o2 <= theta_hmin) * hstat_low +
             (p_b_o2 >  theta_hmin) * (p_b_o2 <= theta_hmax) * hstat_mid +
             (p_b_o2 >  theta_hmax) * hstat_high

    # ===== Total ventilation drive synthesis =====
    alpha_h_effective <- (dp_state < 0) * 1 + (dp_state >= 0) * alpha_h
    chemoreflex_drive_raw <- alpha_h_effective * dp_state + dc_state
    chemoreflex_drive <- (chemoreflex_drive_raw > 0) * chemoreflex_drive_raw
    total_drive <- (w_baseline - kf) + chemoreflex_drive
    venti <- (total_drive > 0) * total_drive

    # ===== Cardiac arrest tripwire dynamics (unchanged from Mann 2022) =====
    pao2_below <- (p_a_o2 < ca_pao2_threshold) * 1.0
    mtime(t_arrest_check) <- first_cross_time + ca_sustain_min
    ca_arrest_decision <- (t >= t_arrest_check) * (crossed_latch > 0.5) *
                          (check_complete < 0.5) * pao2_below

    # ===== Gas-exchange ODEs (Magosso / Ursino) =====
    venti_alveolar <- venti * b_max
    d/dt(palv_co2) <- (venti_alveolar * (p_i_co2 - palv_co2) +
                      lumbda * q_total * (1 - s1_shunt) * (c_v_co2 - c_e_co2)) / v_a
    d/dt(palv_o2)  <- (venti_alveolar * (p_i_o2_init - palv_o2) +
                      lumbda * q_total * (1 - s1_shunt) * (c_v_o2 - c_e_o2)) / v_a

    d/dt(cb_co2) <- (qb * (c_a_co2 - c_vb_co2) + m_b_co2_real * m_o2_safety_brain) / v_b
    d/dt(cb_o2)  <- (qb * (c_a_o2  - c_vb_o2 ) + m_b_o2_real  * m_o2_safety_brain) / v_b

    d/dt(ct_co2) <- (qt * (c_a_co2 - c_vt_co2) + m_t_co2_real * m_o2_safety_tissue) / v_t
    d/dt(ct_o2)  <- (qt * (c_a_o2  - c_vt_o2 ) + m_t_o2_real  * m_o2_safety_tissue) / v_t

    # CBF response filters: yco2 (CO2, age tau), yo2 (tissue O2, adult tau),
    # y_child_o2 (brain O2, child tau). (Chakravartula 2024 eqs 1, 5, 3.)
    d/dt(yco2) <- (psai_co2 - yco2) / tau_co2_min
    d/dt(yo2)  <- (psai_o2  - yo2)  / tau_o2_min
    d/dt(y_child_o2) <- (psai_o2_child - y_child_o2) / tau_child_o2_min

    d/dt(dp_state) <- (-dp_state + av1 * g_dp * (plag_f_pc - f_pc_0)) / tau_dp_min
    d/dt(dc_state) <- (-dc_state + av2 * g_dc * (clag_p_b_co2 - p_b_co2_0)) / tau_dc_min

    d/dt(alpha_h) <- (hstat - alpha_h) / tau_h_min

    d/dt(crossed_latch) <- (1 - crossed_latch) * pao2_below * 100
    d/dt(first_cross_time) <- (crossed_latch < 0.5) * 1.0
    d/dt(check_complete) <- (1 - check_complete) * (t >= t_arrest_check) *
                            (crossed_latch > 0.5) * 10
    d/dt(ca_arrest_latch) <- (1 - ca_arrest_latch) * ca_arrest_decision * 1000
    d/dt(im_arrest) <- -ca_arrest_latch * ca_decay_rate * im_arrest

    # ===== Initial conditions (FDA delaystates.R nominal; see population
    # notes - pre-equilibrate for overdose simulations). y_child_o2 starts
    # at 0 like the other CBF-response filters. =====
    palv_co2(0)      <- 40.28
    palv_o2(0)       <- 100.2
    cb_co2(0)        <- 0.645
    cb_o2(0)         <- 9.78e-4
    ct_co2(0)        <- 0.605
    ct_o2(0)         <- 13e-4
    yco2(0)          <- 0
    yo2(0)           <- 0
    y_child_o2(0)    <- 0
    dp_state(0)      <- 0
    dc_state(0)      <- 0
    alpha_h(0)       <- 1
    crossed_latch(0)   <- 0
    first_cross_time(0) <- 0
    check_complete(0)  <- 0
    ca_arrest_latch(0) <- 0
    im_arrest(0)       <- 1

    # ===== Outputs =====
    Venti         <- venti                 # minute ventilation (L/min)
    P_a_O2        <- p_a_o2                # arterial O2 partial pressure (mm Hg)
    P_a_CO2       <- p_a_co2               # arterial CO2 partial pressure (mm Hg)
    Q_brain       <- qb                    # cerebral blood flow (L/min)
    Q_tissue      <- qt                    # peripheral-tissue blood flow (L/min)
    Q_total       <- q_total               # total cardiac output (L/min)
    P_B_O2        <- p_b_o2                # brain O2 partial pressure (mm Hg)

    # Arterial oxygen saturation (Severinghaus-style approximation)
    nSat <- 2.6
    k3_sat <- 26.6
    p_a_o2_safe <- (p_a_o2 > 1e-6) * p_a_o2 + (p_a_o2 <= 1e-6) * 1e-6
    po2_virt <- p_a_o2_safe * (40 / p_a_o2_safe)^0.3
    SAT_O2 <- 100 * po2_virt^nSat / (k3_sat^nSat + po2_virt^nSat)

    CA_active        <- ca_arrest_latch
    First_cross_time <- first_cross_time
    T_arrest_check   <- t_arrest_check
  })
}
