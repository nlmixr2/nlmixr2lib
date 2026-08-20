Luo_2024_perindopril_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, custom WinNonlin 8.1 implementation). Joint",
    "perindopril + perindoprilat disposition in healthy adults and in liver",
    "cirrhosis (Child-Pugh A/B/C). Perindopril is an inactive ester prodrug",
    "hydrolysed by hepatic CES1 to the active diacid perindoprilat, which is",
    "eliminated renally. Perindopril has a second, parallel hepatic pathway",
    "(UGT glucuronidation) that does NOT form perindoprilat; because the",
    "responsible UGT isoenzyme is unidentified, its cirrhosis scaling uses",
    "the metoprolol AUC0-inf change rate of 0.62 (Section 3.1.4). The",
    "semi-PBPK circuit is stomach, three small-intestinal lumen segments",
    "(duodenum / jejunum / ileum), the matching three gut-wall segments,",
    "portal vein, liver, kidney and a two-compartment systemic compartment.",
    "Cirrhosis is applied by switching the Child-Pugh-specific physiology of",
    "Table 1 (organ blood flows, functional liver volume, GI transit rates,",
    "GFR, plasma-binding protein concentrations, hepatic CES1 content) and",
    "rescaling CLint, CLint,K, Peff and Vsys through Eq 1-6. Deterministic:",
    "the paper's virtual populations are uniform 80-120% draws on the drug",
    "parameters, not lognormal random effects, so no IIV is encoded."
  )
  reference <- paste(
    "Luo X, Zhang Z, Mu R, Hu G, Liu L, Liu X (2024). Simultaneously",
    "Predicting the Pharmacokinetics of CES1-Metabolized Drugs and Their",
    "Metabolites Using Physiologically Based Pharmacokinetic Model in",
    "Cirrhosis Subjects. Pharmaceutics 16(2):234.",
    "doi:10.3390/pharmaceutics16020234.",
    "ODE system from Supplementary Material Eq S1-S20; drug parameters",
    "from Table 2; system physiology and Child-Pugh scaling from Table 1",
    "and Eq 1-6; validation targets from Table 7.",
    sep = " "
  )
  vignette <- "Luo_2024_CES1_cirrhosis"
  units    <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  # Segment-resolved gut-WALL states. The gut LUMEN segments use the
  # canonical stomach / duodenum / jejunum / ileum names; the tissue (wall)
  # states between lumen and portal vein have no canonical yet, so they are
  # declared paper-specific pending a second segment-resolved absorption
  # PBPK to ratify a canonical trio. portal_vein follows
  # vandenBerg_2021_uprifosbuvir_pbpk.R.
  paper_specific_compartments <- c(
    "wall_duodenum", "wall_jejunum", "wall_ileum", "portal_vein", "wall_duodenum_perat", "wall_jejunum_perat", "wall_ileum_perat", "portal_vein_perat", "liver_perat", "kidney_perat"
  )

  covariateData <- list(
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (1 = Child-Pugh class A).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy hepatic function when HEPIMP_MOD and HEPIMP_SEV are also 0)",
      notes              = paste(
        "Classification scheme is Child-Pugh, NOT NCI ODWG: HEPIMP_MILD = 1",
        "selects the Child-Pugh A column of Luo 2024 Table 1. The three",
        "HEPIMP_* indicators are mutually exclusive; all three 0 selects the",
        "Normal column. In Child-Pugh A the hepatic CES1 content is unchanged",
        "from healthy (fCES1 = 1) while functional liver volume is already",
        "reduced to 81% of normal, so the CP-A effect is driven by liver",
        "volume, blood-flow redistribution, GFR, albumin and GI transit."
      ),
      source_name        = "Child-Pugh A"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = Child-Pugh class B).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy hepatic function when HEPIMP_MILD and HEPIMP_SEV are also 0)",
      notes              = paste(
        "Classification scheme is Child-Pugh, NOT NCI ODWG: HEPIMP_MOD = 1",
        "selects the Child-Pugh B column of Luo 2024 Table 1 (hepatic CES1",
        "content 1.715 mg/g liver = 70% of healthy; functional liver volume",
        "65% of normal)."
      ),
      source_name        = "Child-Pugh B"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment indicator (1 = Child-Pugh class C).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy hepatic function when HEPIMP_MILD and HEPIMP_MOD are also 0)",
      notes              = paste(
        "Classification scheme is Child-Pugh, NOT NCI ODWG: HEPIMP_SEV = 1",
        "selects the Child-Pugh C column of Luo 2024 Table 1 (hepatic CES1",
        "content 0.735 mg/g liver = 30% of healthy; functional liver volume",
        "53% of normal; hepatic arterial flow raised to 1020 mL/min)."
      ),
      source_name        = "Child-Pugh C"
    )
  )

  compartmentData <- list(
    stomach                  = list(analyte = "perindopril", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum                 = list(analyte = "perindopril", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum                  = list(analyte = "perindopril", units = "mg", specimen = "administration site", verified = TRUE),
    ileum                    = list(analyte = "perindopril", units = "mg", specimen = "administration site", verified = TRUE),
    wall_duodenum            = list(analyte = "perindopril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum             = list(analyte = "perindopril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum               = list(analyte = "perindopril", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein              = list(analyte = "perindopril", units = "mg", specimen = "plasma", verified = TRUE),
    liver                    = list(analyte = "perindopril", units = "mg", specimen = "tissue", verified = TRUE),
    kidney                   = list(analyte = "perindopril", units = "mg", specimen = "tissue", verified = TRUE),
    central                  = list(analyte = "perindopril", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1              = list(analyte = "perindopril", units = "mg", specimen = "plasma", verified = TRUE),
    wall_duodenum_perat      = list(analyte = "perindoprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum_perat       = list(analyte = "perindoprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum_perat         = list(analyte = "perindoprilat", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein_perat        = list(analyte = "perindoprilat", units = "mg", specimen = "plasma", verified = TRUE),
    liver_perat              = list(analyte = "perindoprilat", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_perat             = list(analyte = "perindoprilat", units = "mg", specimen = "tissue", verified = TRUE),
    central_perat            = list(analyte = "perindoprilat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_perat        = list(analyte = "perindoprilat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 4L,
    age_range      = "adults",
    disease_state  = paste(
      "Pooled healthy volunteers and liver cirrhosis patients. Healthy: Lecocq",
      "1990 (n = 12), Lees 1988 (n = 8). Cirrhosis: Tsai 1989 Child-Pugh A (n",
      "= 8), Thiollet 1992 Child-Pugh B (n = 10). Section 3.1.4 notes the",
      "cirrhosis reports supply pharmacokinetic parameters only, with no",
      "concentration-time profiles, so only a parameter comparison was",
      "possible."
    ),
    dose_range     = "Perindopril tert-butylamine 4 and 8 mg single oral doses",
    notes          = paste(
      "Luo 2024 Table 3. Literature-digitised clinical data; the authors",
      "simulated 1000 virtual individuals per population by drawing CLint,",
      "CLint,K, fu,b, Vsys, Peff, ka, KL:P, KG:P and KK:P uniformly over",
      "80-120% of the Table 2 point values (Section 2.3). Validation targets",
      "are in Table 7."
    )
  )

  ini({
    # ================= PERINDOPRIL (parent prodrug) =================
    # Luo 2024 Table 2. Clearances converted mL/min -> L/min.
    lcl_int_h            <- fixed(log(1.01115));      label("perindopril hepatic CES1 intrinsic clearance (log L/min)")  # Table 2: CLint 1011.15 mL/min [ref 68] (CES1-mediated)
    lcl_int_k            <- fixed(log(0.1302));       label("perindopril renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 130.2 mL/min [ref 71]
    lvc                  <- fixed(log(13.119));       label("perindopril systemic apparent volume of distribution (log L)")  # Table 2: Vsys 13.119 L [ref 70] (footnote g)
    lk12                 <- fixed(log(0.0028));       label("perindopril central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.0028 1/min [ref 70] (footnote g)
    lk21                 <- fixed(log(0.0024));       label("perindopril peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.0024 1/min [ref 70] (footnote g)
    lpeff                <- fixed(log(0.000134));     label("perindopril effective permeability Peff,A-B (log cm/s)")  # Table 2: Peff,A-B 1.34e-4 cm/s [ref 43]
    kp_liver             <- fixed(0.665);             label("perindopril liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut               <- fixed(0.633);             label("perindopril gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney            <- fixed(0.742);             label("perindopril kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr                  <- fixed(1);                 label("perindopril blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 1 (footnote e, assumed)
    fu_p                 <- fixed(0.4);               label("perindopril fraction unbound in plasma (unitless)")  # Section 3.1.4 percent binding 60% -> fu,p 0.40 [ref 72]; Table 2 fu,b 0.4 with Rb 1 via Eq S8
    lcl_int_ugt          <- fixed(log(0.15647));      label("perindopril hepatic UGT intrinsic clearance (log L/min)")  # Table 2: 156.47 mL/min (footnotes c and j, UGT intrinsic clearance of perindopril) [ref 69]

    # ================= PERINDOPRILAT (active metabolite) =================
    # Luo 2024 Table 2, perindoprilat row.
    lcl_int_k_perat      <- fixed(log(0.23178));      label("perindoprilat renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 231.78 mL/min [ref 74]
    lvc_perat            <- fixed(log(53.44));        label("perindoprilat systemic apparent volume of distribution (log L)")  # Table 2: Vsys 53.44 L [ref 73] (footnote f)
    lk12_perat           <- fixed(log(0.271));        label("perindoprilat central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.271 1/min [ref 73] (footnote f)
    lk21_perat           <- fixed(log(0.0996));       label("perindoprilat peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.0996 1/min [ref 73] (footnote f)
    kp_liver_perat       <- fixed(1.45);              label("perindoprilat liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut_perat         <- fixed(1.38);              label("perindoprilat gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney_perat      <- fixed(1.61);              label("perindoprilat kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr_perat            <- fixed(1);                 label("perindoprilat blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 1 (footnote e, assumed)
    fu_p_perat           <- fixed(0.85);              label("perindoprilat fraction unbound in plasma (unitless)")  # Section 3.1.4 mean percent binding 15% -> fu,p 0.85 [ref 72]; Table 2 fu,b 0.85 with Rb 1 via Eq S8

    # ================= Residual error =================
    # Luo 2024 report no residual-error model: the semi-PBPK model is a
    # forward predictor assessed by 0.5-2-fold AUC/Cmax agreement and
    # 5th-95th percentile coverage (Section 2.4), not by a fitted sigma.
    # Non-fitted placeholders per the in-repo PBPK convention
    # (An_2012_mitoxantrone_human_pbpk.R, Gaohua_2012_pregnancy_pbpk_midazolam.R).
    propSd               <- fixed(0.10);  label("Proportional residual error placeholder, perindopril (fraction)")  # not reported by Luo 2024
    propSd_perat         <- fixed(0.10);  label("Proportional residual error placeholder, perindoprilat (fraction)")  # not reported by Luo 2024
  })

  model({
    # ============================================================
    # 1. Child-Pugh decode. The three HEPIMP_* indicators are mutually
    #    exclusive; cp_0 = 1 selects the Normal column of Table 1.
    # ============================================================
    cp_a <- HEPIMP_MILD
    cp_b <- HEPIMP_MOD
    cp_c <- HEPIMP_SEV
    cp_0 <- 1 - cp_a - cp_b - cp_c
    cp_lc <- cp_a + cp_b + cp_c

    # ============================================================
    # 2. System physiology, Luo 2024 Table 1 (Normal | CP-A | CP-B | CP-C).
    #    Flows mL/min -> L/min, volumes mL -> L. QL = QLA + QPV holds
    #    exactly in every column (Table 1 footnote a).
    # ============================================================
    q_liver <- cp_0 * 1.4500 + cp_a * 1.43650 + cp_b * 1.17690 + cp_c * 1.65630  # liver flow 1450 / 1436.5 / 1176.9 / 1656.3 mL/min
    q_ha    <- cp_0 * 0.3000 + cp_a * 0.39000 + cp_b * 0.48690 + cp_c * 1.02000  # hepatic arterial 300 / 390 / 486.9 / 1020 mL/min
    q_pv    <- cp_0 * 1.1500 + cp_a * 1.04650 + cp_b * 0.69000 + cp_c * 0.63630  # portal vein 1150 / 1046.5 / 690 / 636.3 mL/min
    q_kid   <- cp_0 * 1.2400 + cp_a * 1.09120 + cp_b * 0.80600 + cp_c * 0.59520  # kidney 1240 / 1091.2 / 806 / 595.2 mL/min
    q_duo   <- 0.045                                                             # duodenum 45 mL/min (unchanged in cirrhosis, footnote b)
    q_jej   <- 0.173                                                             # jejunum 173 mL/min (footnote b)
    q_ile   <- 0.102                                                             # ileum 102 mL/min (footnote b)
    q_gut   <- q_duo + q_jej + q_ile

    v_liver <- cp_0 * 1.6900 + cp_a * 1.36890 + cp_b * 1.09850 + cp_c * 0.89570  # liver volume 1690 / 1368.9 / 1098.5 / 895.7 mL
    v_pv    <- 0.070                                                             # portal vein 70 mL (footnote b)
    v_kid   <- 0.280                                                             # kidney 280 mL (footnote b)
    v_duo   <- 0.021                                                             # duodenum wall 21 mL (footnote b)
    v_jej   <- 0.063                                                             # jejunum wall 63 mL (footnote b)
    v_ile   <- 0.042                                                             # ileum wall 42 mL (footnote b)

    # ============================================================
    # 3. Cirrhosis scaling ratios (Table 1 -> Eq 1-6). Every ratio is
    #    exactly 1 in the Normal column, so Eq 1-6 apply unconditionally
    #    with no branching.
    # ============================================================
    f_ces1  <- cp_0 * 1 + cp_a * 1.00 + cp_b * 0.70 + cp_c * 0.30               # CES1 2.45 / 2.45 / 1.715 / 0.735 mg/g liver
    f_liver <- cp_0 * 1 + cp_a * 0.81 + cp_b * 0.65 + cp_c * 0.53               # functional liver volume ratio vs 1690 mL
    r_gfr   <- cp_0 * 1 + cp_lc * (82 / 105)                                    # GFR 105 -> 82 mL/min
    r_alb   <- cp_0 * 1 + cp_a * (36.2 / 44.7) + cp_b * (30.4 / 44.7) + cp_c * (26.3 / 44.7)  # albumin 44.7 -> 36.2 / 30.4 / 26.3 g/L
    r_lr    <- cp_0 * 1 + cp_a * (0.046 / 0.037) + cp_b * (0.052 / 0.037) + cp_c * (0.057 / 0.037)  # lactulose/rhamnose ratio (intestinal permeability)
    # Section 3.1.4: the UGT isoenzyme metabolising perindopril is unidentified,
    # so the metoprolol AUC0-inf change rate (0.62) is used as the cirrhosis
    # scaling factor `fother` of Eq 2 for all three Child-Pugh classes.
    f_other <- cp_0 * 1 + cp_lc * 0.62

    r_duo   <- 2.00                                                              # gut radius r1 2 cm (footnote b)
    r_jej   <- 1.63                                                              # gut radius r2 1.63 cm (footnote b)
    r_ile   <- 1.45                                                              # gut radius r3 1.45 cm (footnote b)
    kt_sto  <- cp_0 * 0.04 + cp_lc * 0.0504                                      # stomach transit 0.04 -> 0.0504 1/min (footnote c)
    kt_duo  <- cp_0 * 0.07 + cp_lc * 0.0889                                      # duodenum transit 0.07 -> 0.0889 1/min (footnote c)
    kt_jej  <- cp_0 * 0.03 + cp_lc * 0.0381                                      # jejunum transit 0.03 -> 0.0381 1/min (footnote c)
    kt_ile  <- cp_0 * 0.04 + cp_lc * 0.0508                                      # ileum transit 0.04 -> 0.0508 1/min (footnote c)

    # ============================================================
    # 4. Drug parameters after cirrhosis scaling.
    #    Eq 3: fu,p,CI = 1 / (1 + (1 - fu,p,HT) * Pratio / fu,p,HT)
    #    Eq S8: fu,b = fu,p / Rb          Eq 4: Vsys,CI = (fu,p,CI/fu,p,HT) * Vsys
    #    Eq 1: CLint,CES1 * fCES1 * fliver   Eq 5: CLint,K * GFR ratio
    #    Eq 6: Peff,CI = Peff,HT * LR ratio
    #    perindopril binds mainly to albumin (Section 2.3).
    # ============================================================
    fu_p_i   <- 1 / (1 + (1 - fu_p) * r_alb / fu_p)
    fu_b_i   <- fu_p_i / bpr
    v_sys    <- exp(lvc) * (fu_p_i / fu_p)
    cl_int_h <- exp(lcl_int_h) * f_ces1 * f_liver
    cl_int_k <- exp(lcl_int_k) * r_gfr
    cl_int_ugt <- exp(lcl_int_ugt) * f_other * f_liver          # Eq 2 (non-CES1 hepatic route)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Eq S3: ka,i = 2 * Peff / ri. Peff is cm/s and ri cm, so 60 converts 1/s -> 1/min.
    peff   <- exp(lpeff) * r_lr
    ka_duo <- 2 * peff * 60 / r_duo
    ka_jej <- 2 * peff * 60 / r_jej
    ka_ile <- 2 * peff * 60 / r_ile

    fu_p_i_perat  <- 1 / (1 + (1 - fu_p_perat) * r_alb / fu_p_perat)
    fu_b_i_perat  <- fu_p_i_perat / bpr_perat
    v_sys_perat   <- exp(lvc_perat) * (fu_p_i_perat / fu_p_perat)
    cl_int_k_perat <- exp(lcl_int_k_perat) * r_gfr
    k12_perat <- exp(lk12_perat)
    k21_perat <- exp(lk21_perat)

    # ============================================================
    # 5. Concentrations from state amounts (every state is an amount in mg)
    #    and the blood-side effluent concentration Rb / Kp leaving each tissue.
    # ============================================================
    c_wall_duo <- wall_duodenum / v_duo
    c_wall_jej <- wall_jejunum  / v_jej
    c_wall_ile <- wall_ileum    / v_ile
    c_pv       <- portal_vein   / v_pv
    c_liver    <- liver         / v_liver
    c_kid      <- kidney        / v_kid
    c_sys      <- central       / v_sys
    eff_wall_duo <- c_wall_duo * bpr / kp_gut
    eff_wall_jej <- c_wall_jej * bpr / kp_gut
    eff_wall_ile <- c_wall_ile * bpr / kp_gut
    eff_liver    <- c_liver    * bpr / kp_liver
    eff_kid      <- c_kid      * bpr / kp_kidney

    c_wall_duo_perat <- wall_duodenum_perat / v_duo
    c_wall_jej_perat <- wall_jejunum_perat  / v_jej
    c_wall_ile_perat <- wall_ileum_perat    / v_ile
    c_pv_perat       <- portal_vein_perat   / v_pv
    c_liver_perat    <- liver_perat         / v_liver
    c_kid_perat      <- kidney_perat        / v_kid
    c_sys_perat      <- central_perat       / v_sys_perat
    eff_wall_duo_perat <- c_wall_duo_perat * bpr_perat / kp_gut_perat
    eff_wall_jej_perat <- c_wall_jej_perat * bpr_perat / kp_gut_perat
    eff_wall_ile_perat <- c_wall_ile_perat * bpr_perat / kp_gut_perat
    eff_liver_perat    <- c_liver_perat    * bpr_perat / kp_liver_perat
    eff_kid_perat      <- c_kid_perat      * bpr_perat / kp_kidney_perat

    # Hepatic CES1 hydrolysis flux (mg/min): Eq S7 loss term for the parent
    # and Eq S13 formation term for the metabolite.
    hydrolysis <- cl_int_h * fu_b_i * eff_liver
    # Parallel UGT glucuronidation removes parent WITHOUT forming perindoprilat.
    ugt_loss   <- cl_int_ugt * fu_b_i * eff_liver

    # ============================================================
    # 6. ODE system (Luo 2024 Supplementary Material).
    #    Parent: Eq S1 stomach, S2 lumen, S5 gut wall, S6 portal vein,
    #    S7 liver, S14 kidney, S16-S17 (2-cmt) systemic.
    #    perindoprilat: S5/S6 without absorption, S13 liver, S14 kidney, S16-S17 (2-cmt).
    # ============================================================
    d/dt(stomach)  <- -kt_sto * stomach
    d/dt(duodenum) <-  kt_sto * stomach  - kt_duo * duodenum - ka_duo * duodenum
    d/dt(jejunum)  <-  kt_duo * duodenum - kt_jej * jejunum  - ka_jej * jejunum
    d/dt(ileum)    <-  kt_jej * jejunum  - kt_ile * ileum    - ka_ile * ileum

    d/dt(wall_duodenum) <- ka_duo * duodenum + q_duo * c_sys - q_duo * eff_wall_duo
    d/dt(wall_jejunum)  <- ka_jej * jejunum  + q_jej * c_sys - q_jej * eff_wall_jej
    d/dt(wall_ileum)    <- ka_ile * ileum    + q_ile * c_sys - q_ile * eff_wall_ile

    d/dt(portal_vein) <- q_duo * eff_wall_duo + q_jej * eff_wall_jej +
                         q_ile * eff_wall_ile - q_pv * c_pv

    d/dt(liver) <- q_pv * c_pv + q_ha * c_sys - q_liver * eff_liver - hydrolysis - ugt_loss

    d/dt(kidney) <- q_kid * c_sys - (q_kid + fu_b_i * cl_int_k) * eff_kid

    d/dt(central) <- q_liver * eff_liver + q_kid * eff_kid -
                     (q_ha + q_kid) * c_sys - q_gut * c_sys
                     + k21 * peripheral1 - k12 * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- perindoprilat: no lumen states (it is not absorbed from the gut lumen) ----
    d/dt(wall_duodenum_perat) <- q_duo * c_sys_perat - q_duo * eff_wall_duo_perat
    d/dt(wall_jejunum_perat)  <- q_jej * c_sys_perat - q_jej * eff_wall_jej_perat
    d/dt(wall_ileum_perat)    <- q_ile * c_sys_perat - q_ile * eff_wall_ile_perat

    d/dt(portal_vein_perat) <- q_duo * eff_wall_duo_perat + q_jej * eff_wall_jej_perat +
                            q_ile * eff_wall_ile_perat - q_pv * c_pv_perat

    d/dt(liver_perat) <- q_pv * c_pv_perat + q_ha * c_sys_perat + hydrolysis -
                      q_liver * eff_liver_perat

    d/dt(kidney_perat) <- q_kid * c_sys_perat -
                       (q_kid + fu_b_i_perat * cl_int_k_perat) * eff_kid_perat

    d/dt(central_perat) <- q_liver * eff_liver_perat + q_kid * eff_kid_perat -
                        (q_ha + q_kid) * c_sys_perat - q_gut * c_sys_perat
                        + k21_perat * peripheral1_perat - k12_perat * central_perat
    d/dt(peripheral1_perat) <- k12_perat * central_perat - k21_perat * peripheral1_perat

    # Table 2 also reports an oral bioavailability F = 0.66 for this drug
    # (Table 2: F 0.66 [ref 64]; Section 3.1.4 bioavailability 66%).
    # It is deliberately NOT applied as a dose multiplier: F is a reported
    # literature property that this mechanistic model PREDICTS from the
    # segment absorption competition (Eq S2-S3) and the well-stirred hepatic
    # first pass (Eq S7), not an input it consumes. Applying it again
    # double-counts -- for benazepril the model's own predicted
    # Fa * Fh = 0.441 * 0.878 = 0.387, reproducing the tabulated 0.35, and
    # applying 0.35 on top drove AUC and Cmax to ~0.31-0.36 of the Table 5
    # predictions. Dropping it brings benazepril, temocapril and perindopril
    # all into line with the five drugs that have no F entry. See the
    # vignette's Assumptions and deviations section.

    # ============================================================
    # 7. Observations. Amounts are mg and volumes L, so concentrations are
    #    mg/L = ug/mL, matching the AUC unit (ug*h/mL) of Table 7.
    # ============================================================
    Cc <- c_sys
    Cc_perat <- c_sys_perat

    Cc ~ prop(propSd)
    Cc_perat ~ prop(propSd_perat)
  })
}

