Luo_2024_enalapril_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, custom WinNonlin 8.1 implementation). Joint",
    "enalapril + enalaprilat disposition in healthy adults and in liver",
    "cirrhosis (Child-Pugh A/B/C). Enalapril is an inactive ester prodrug",
    "hydrolysed by hepatic carboxylesterase 1 (CES1) to the active diacid",
    "enalaprilat, which is cleared renally. The semi-PBPK circuit is",
    "stomach, three small-intestinal lumen segments (duodenum / jejunum /",
    "ileum), the matching three gut-wall segments, portal vein, liver,",
    "kidney and a one-compartment systemic compartment. Cirrhosis is applied",
    "by switching the Child-Pugh-specific physiology of Table 1 (organ blood",
    "flows, functional liver volume, GI transit rates, GFR, plasma-binding",
    "protein concentrations, hepatic CES1 content) and rescaling CLint,",
    "CLint,K, Peff and Vsys through Eq 1-6. Deterministic: the paper's",
    "virtual populations are uniform 80-120% draws on the drug parameters,",
    "not lognormal random effects, so no IIV is encoded."
  )
  reference <- paste(
    "Luo X, Zhang Z, Mu R, Hu G, Liu L, Liu X (2024). Simultaneously",
    "Predicting the Pharmacokinetics of CES1-Metabolized Drugs and Their",
    "Metabolites Using Physiologically Based Pharmacokinetic Model in",
    "Cirrhosis Subjects. Pharmaceutics 16(2):234.",
    "doi:10.3390/pharmaceutics16020234.",
    "ODE system from Supplementary Material Eq S1-S20; drug parameters",
    "from Table 2; system physiology and Child-Pugh scaling from Table 1",
    "and Eq 1-6; validation targets from Table 4.",
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
    "wall_duodenum", "wall_jejunum", "wall_ileum", "portal_vein", "wall_duodenum_enaat", "wall_jejunum_enaat", "wall_ileum_enaat", "portal_vein_enaat", "liver_enaat", "kidney_enaat"
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
    stomach                  = list(analyte = "enalapril", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum                 = list(analyte = "enalapril", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum                  = list(analyte = "enalapril", units = "mg", specimen = "administration site", verified = TRUE),
    ileum                    = list(analyte = "enalapril", units = "mg", specimen = "administration site", verified = TRUE),
    wall_duodenum            = list(analyte = "enalapril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum             = list(analyte = "enalapril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum               = list(analyte = "enalapril", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein              = list(analyte = "enalapril", units = "mg", specimen = "plasma", verified = TRUE),
    liver                    = list(analyte = "enalapril", units = "mg", specimen = "tissue", verified = TRUE),
    kidney                   = list(analyte = "enalapril", units = "mg", specimen = "tissue", verified = TRUE),
    central                  = list(analyte = "enalapril", units = "mg", specimen = "plasma", verified = TRUE),
    wall_duodenum_enaat      = list(analyte = "enalaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum_enaat       = list(analyte = "enalaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum_enaat         = list(analyte = "enalaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein_enaat        = list(analyte = "enalaprilat", units = "mg", specimen = "plasma", verified = TRUE),
    liver_enaat              = list(analyte = "enalaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_enaat             = list(analyte = "enalaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    central_enaat            = list(analyte = "enalaprilat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_enaat        = list(analyte = "enalaprilat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 5L,
    age_range      = "adults",
    disease_state  = paste(
      "Pooled healthy volunteers and liver cirrhosis patients. Healthy:",
      "Ohnishi 1989 (n = 7), Todd 1986 (n = 12), Weisser 1991 (n = 8),",
      "Dickstein 1987 (n = 10). Cirrhosis: Baba 1990 Child-Pugh B (n = 7),",
      "Ohnishi 1989 Child-Pugh C (n = 7)."
    ),
    dose_range     = "Enalapril maleate 10 mg single oral dose",
    notes          = paste(
      "Luo 2024 Table 3. Literature-digitised clinical data; the authors",
      "simulated 1000 virtual individuals per population by drawing CLint,",
      "CLint,K, fu,b, Vsys, Peff, ka, KL:P, KG:P and KK:P uniformly over",
      "80-120% of the Table 2 point values (Section 2.3). Validation targets",
      "are in Table 4."
    )
  )

  ini({
    # ================= ENALAPRIL (parent prodrug) =================
    # Luo 2024 Table 2. Clearances converted mL/min -> L/min.
    lcl_int_h            <- fixed(log(0.784));        label("enalapril hepatic CES1 intrinsic clearance (log L/min)")  # Table 2: CLint 784 mL/min [ref 28]
    lcl_int_k            <- fixed(log(0.6246));       label("enalapril renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 624.6 mL/min [ref 31]
    lvc                  <- fixed(log(40));           label("enalapril systemic apparent volume of distribution (log L)")  # Table 2: Vsys 40 L [ref 29]
    lpeff                <- fixed(log(0.00016));      label("enalapril effective permeability Peff,A-B (log cm/s)")  # Table 2: Peff,A-B 1.60e-4 cm/s [ref 30]
    kp_liver             <- fixed(1.66);              label("enalapril liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut               <- fixed(2.29);              label("enalapril gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney            <- fixed(1.79);              label("enalapril kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr                  <- fixed(0.74);              label("enalapril blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 0.74 [ref 32]
    fu_p                 <- fixed(0.55);              label("enalapril fraction unbound in plasma (unitless)")  # Section 3.1.1 free fraction in plasma 0.55 [ref 33]; Table 2 fu,b 0.74 = 0.55/0.74 via Eq S8

    # ================= ENALAPRILAT (active metabolite) =================
    # Luo 2024 Table 2, enalaprilat row.
    lcl_int_k_enaat      <- fixed(log(0.1864));       label("enalaprilat renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 186.4 mL/min [ref 35]
    lvc_enaat            <- fixed(log(46.1));         label("enalaprilat systemic apparent volume of distribution (log L)")  # Table 2: Vsys 46.1 L [ref 34]
    lk12_enaat           <- fixed(log(0.001));        label("enalaprilat central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.001 1/min [ref 34]
    lk21_enaat           <- fixed(log(0.0009));       label("enalaprilat peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.0009 1/min [ref 34]
    kp_liver_enaat       <- fixed(1.12);              label("enalaprilat liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut_enaat         <- fixed(1.04);              label("enalaprilat gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney_enaat      <- fixed(1.25);              label("enalaprilat kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr_enaat            <- fixed(0.73);              label("enalaprilat blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 0.73 [ref 32]
    fu_p_enaat           <- fixed(0.5);               label("enalaprilat fraction unbound in plasma (unitless)")  # Section 3.1.1 free fraction in plasma 0.50 [ref 33]; Table 2 fu,b 0.68 = 0.50/0.73 via Eq S8

    # ================= Residual error =================
    # Luo 2024 report no residual-error model: the semi-PBPK model is a
    # forward predictor assessed by 0.5-2-fold AUC/Cmax agreement and
    # 5th-95th percentile coverage (Section 2.4), not by a fitted sigma.
    # Non-fitted placeholders per the in-repo PBPK convention
    # (An_2012_mitoxantrone_human_pbpk.R, Gaohua_2012_pregnancy_pbpk_midazolam.R).
    propSd               <- fixed(0.10);  label("Proportional residual error placeholder, enalapril (fraction)")  # not reported by Luo 2024
    propSd_enaat         <- fixed(0.10);  label("Proportional residual error placeholder, enalaprilat (fraction)")  # not reported by Luo 2024
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
    #    enalapril binds mainly to albumin (Section 2.3).
    # ============================================================
    fu_p_i   <- 1 / (1 + (1 - fu_p) * r_alb / fu_p)
    fu_b_i   <- fu_p_i / bpr
    v_sys    <- exp(lvc) * (fu_p_i / fu_p)
    cl_int_h <- exp(lcl_int_h) * f_ces1 * f_liver
    cl_int_k <- exp(lcl_int_k) * r_gfr

    # Eq S3: ka,i = 2 * Peff / ri. Peff is cm/s and ri cm, so 60 converts 1/s -> 1/min.
    peff   <- exp(lpeff) * r_lr
    ka_duo <- 2 * peff * 60 / r_duo
    ka_jej <- 2 * peff * 60 / r_jej
    ka_ile <- 2 * peff * 60 / r_ile

    fu_p_i_enaat  <- 1 / (1 + (1 - fu_p_enaat) * r_alb / fu_p_enaat)
    fu_b_i_enaat  <- fu_p_i_enaat / bpr_enaat
    v_sys_enaat   <- exp(lvc_enaat) * (fu_p_i_enaat / fu_p_enaat)
    cl_int_k_enaat <- exp(lcl_int_k_enaat) * r_gfr
    k12_enaat <- exp(lk12_enaat)
    k21_enaat <- exp(lk21_enaat)

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

    c_wall_duo_enaat <- wall_duodenum_enaat / v_duo
    c_wall_jej_enaat <- wall_jejunum_enaat  / v_jej
    c_wall_ile_enaat <- wall_ileum_enaat    / v_ile
    c_pv_enaat       <- portal_vein_enaat   / v_pv
    c_liver_enaat    <- liver_enaat         / v_liver
    c_kid_enaat      <- kidney_enaat        / v_kid
    c_sys_enaat      <- central_enaat       / v_sys_enaat
    eff_wall_duo_enaat <- c_wall_duo_enaat * bpr_enaat / kp_gut_enaat
    eff_wall_jej_enaat <- c_wall_jej_enaat * bpr_enaat / kp_gut_enaat
    eff_wall_ile_enaat <- c_wall_ile_enaat * bpr_enaat / kp_gut_enaat
    eff_liver_enaat    <- c_liver_enaat    * bpr_enaat / kp_liver_enaat
    eff_kid_enaat      <- c_kid_enaat      * bpr_enaat / kp_kidney_enaat

    # Hepatic CES1 hydrolysis flux (mg/min): Eq S7 loss term for the parent
    # and Eq S13 formation term for the metabolite.
    hydrolysis <- cl_int_h * fu_b_i * eff_liver

    # ============================================================
    # 6. ODE system (Luo 2024 Supplementary Material).
    #    Parent: Eq S1 stomach, S2 lumen, S5 gut wall, S6 portal vein,
    #    S7 liver, S14 kidney, S15 (1-cmt) systemic.
    #    enalaprilat: S5/S6 without absorption, S13 liver, S14 kidney, S16-S17 (2-cmt).
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

    d/dt(liver) <- q_pv * c_pv + q_ha * c_sys - q_liver * eff_liver - hydrolysis

    d/dt(kidney) <- q_kid * c_sys - (q_kid + fu_b_i * cl_int_k) * eff_kid

    d/dt(central) <- q_liver * eff_liver + q_kid * eff_kid -
                     (q_ha + q_kid) * c_sys - q_gut * c_sys

    # ---- enalaprilat: no lumen states (it is not absorbed from the gut lumen) ----
    d/dt(wall_duodenum_enaat) <- q_duo * c_sys_enaat - q_duo * eff_wall_duo_enaat
    d/dt(wall_jejunum_enaat)  <- q_jej * c_sys_enaat - q_jej * eff_wall_jej_enaat
    d/dt(wall_ileum_enaat)    <- q_ile * c_sys_enaat - q_ile * eff_wall_ile_enaat

    d/dt(portal_vein_enaat) <- q_duo * eff_wall_duo_enaat + q_jej * eff_wall_jej_enaat +
                            q_ile * eff_wall_ile_enaat - q_pv * c_pv_enaat

    d/dt(liver_enaat) <- q_pv * c_pv_enaat + q_ha * c_sys_enaat + hydrolysis -
                      q_liver * eff_liver_enaat

    d/dt(kidney_enaat) <- q_kid * c_sys_enaat -
                       (q_kid + fu_b_i_enaat * cl_int_k_enaat) * eff_kid_enaat

    d/dt(central_enaat) <- q_liver * eff_liver_enaat + q_kid * eff_kid_enaat -
                        (q_ha + q_kid) * c_sys_enaat - q_gut * c_sys_enaat
                        + k21_enaat * peripheral1_enaat - k12_enaat * central_enaat
    d/dt(peripheral1_enaat) <- k12_enaat * central_enaat - k21_enaat * peripheral1_enaat

    # ============================================================
    # 7. Observations. Amounts are mg and volumes L, so concentrations are
    #    mg/L = ug/mL, matching the AUC unit (ug*h/mL) of Table 4.
    # ============================================================
    Cc <- c_sys
    Cc_enaat <- c_sys_enaat

    Cc ~ prop(propSd)
    Cc_enaat ~ prop(propSd_enaat)
  })
}

