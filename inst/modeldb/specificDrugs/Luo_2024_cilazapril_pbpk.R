Luo_2024_cilazapril_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, custom WinNonlin 8.1 implementation). Joint",
    "cilazapril + cilazaprilat disposition in healthy adults and in liver",
    "cirrhosis (Child-Pugh A/B/C). Cilazapril is an inactive ester prodrug",
    "hydrolysed by hepatic CES1 to the active diacid cilazaprilat, which is",
    "eliminated renally. Absorption is described by a directly reported",
    "first-order ka rather than by a permeability. The semi-PBPK circuit is",
    "stomach, three small-intestinal lumen segments (duodenum / jejunum /",
    "ileum), the matching three gut-wall segments, portal vein, liver,",
    "kidney and a two-compartment systemic compartment. Cirrhosis is applied",
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
    "and Eq 1-6; validation targets from Table 6.",
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
    "wall_duodenum", "wall_jejunum", "wall_ileum", "portal_vein", "wall_duodenum_cilat", "wall_jejunum_cilat", "wall_ileum_cilat", "portal_vein_cilat", "liver_cilat", "kidney_cilat"
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
    stomach                  = list(analyte = "cilazapril", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum                 = list(analyte = "cilazapril", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum                  = list(analyte = "cilazapril", units = "mg", specimen = "administration site", verified = TRUE),
    ileum                    = list(analyte = "cilazapril", units = "mg", specimen = "administration site", verified = TRUE),
    wall_duodenum            = list(analyte = "cilazapril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum             = list(analyte = "cilazapril", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum               = list(analyte = "cilazapril", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein              = list(analyte = "cilazapril", units = "mg", specimen = "plasma", verified = TRUE),
    liver                    = list(analyte = "cilazapril", units = "mg", specimen = "tissue", verified = TRUE),
    kidney                   = list(analyte = "cilazapril", units = "mg", specimen = "tissue", verified = TRUE),
    central                  = list(analyte = "cilazapril", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1              = list(analyte = "cilazapril", units = "mg", specimen = "plasma", verified = TRUE),
    wall_duodenum_cilat      = list(analyte = "cilazaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum_cilat       = list(analyte = "cilazaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum_cilat         = list(analyte = "cilazaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein_cilat        = list(analyte = "cilazaprilat", units = "mg", specimen = "plasma", verified = TRUE),
    liver_cilat              = list(analyte = "cilazaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_cilat             = list(analyte = "cilazaprilat", units = "mg", specimen = "tissue", verified = TRUE),
    central_cilat            = list(analyte = "cilazaprilat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cilat        = list(analyte = "cilazaprilat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 6L,
    age_range      = "adults",
    disease_state  = paste(
      "Pooled healthy volunteers and liver cirrhosis patients. Healthy:",
      "Massarella 1989 (n = 24), Williams 1990 (n = 13), Williams 1989 (n =",
      "12), Massarella 1989 food study (n = 16), Francis 1987 (n = 12), Gross",
      "1993 (n = 10). Cirrhosis: Gross 1993 Child-Pugh B (n = 9)."
    ),
    dose_range     = "Cilazapril 1, 1.25, 2.5, 5 and 10 mg single oral doses",
    notes          = paste(
      "Luo 2024 Table 3. Literature-digitised clinical data; the authors",
      "simulated 1000 virtual individuals per population by drawing CLint,",
      "CLint,K, fu,b, Vsys, Peff, ka, KL:P, KG:P and KK:P uniformly over",
      "80-120% of the Table 2 point values (Section 2.3). Validation targets",
      "are in Table 6."
    )
  )

  ini({
    # ================= CILAZAPRIL (parent prodrug) =================
    # Luo 2024 Table 2. Clearances converted mL/min -> L/min.
    lcl_int_h            <- fixed(log(0.1997));       label("cilazapril hepatic CES1 intrinsic clearance (log L/min)")  # Table 2: CLint 199.7 mL/min (footnote c, recalculated from CL_L,b); CLb 205 mL/min
    lcl_int_k            <- fixed(log(0.118095));     label("cilazapril renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 118.095 mL/min [ref 52]
    lvc                  <- fixed(log(18.23));        label("cilazapril systemic apparent volume of distribution (log L)")  # Table 2: Vsys 18.23 L [ref 51] (footnote f, simulated by WinNonlin)
    lk12                 <- fixed(log(0.00325));      label("cilazapril central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.00325 1/min [ref 51] (footnote f)
    lk21                 <- fixed(log(0.00155));      label("cilazapril peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.00155 1/min [ref 51] (footnote f)
    lka                  <- fixed(log(0.099));        label("cilazapril gut-lumen absorption rate constant (log 1/min)")  # Table 2: ka 0.099 1/min [ref 53] (footnote g, calculated by WinNonlin)
    kp_liver             <- fixed(1.32);              label("cilazapril liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut               <- fixed(1.31);              label("cilazapril gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney            <- fixed(1.43);              label("cilazapril kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr                  <- fixed(1);                 label("cilazapril blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 1 (footnote e, assumed)
    fu_p                 <- fixed(0.7);               label("cilazapril fraction unbound in plasma (unitless)")  # Section 3.1.3 free fraction in plasma 0.70 [ref 50]; Table 2 fu,b 0.7 with Rb 1 via Eq S8

    # ================= CILAZAPRILAT (active metabolite) =================
    # Luo 2024 Table 2, cilazaprilat row.
    lcl_int_k_cilat      <- fixed(log(0.07548));      label("cilazaprilat renal intrinsic clearance (log L/min)")  # Table 2: CLint,K 75.48 mL/min [ref 52]
    lvc_cilat            <- fixed(log(10.3517));      label("cilazaprilat systemic apparent volume of distribution (log L)")  # Table 2: Vsys 10.3517 L [ref 51] (footnote f)
    lk12_cilat           <- fixed(log(0.00084));      label("cilazaprilat central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.00084 1/min [ref 51] (footnote f)
    lk21_cilat           <- fixed(log(0.008));        label("cilazaprilat peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.008 1/min [ref 51] (footnote f)
    kp_liver_cilat       <- fixed(1.28);              label("cilazaprilat liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut_cilat         <- fixed(1.22);              label("cilazaprilat gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney_cilat      <- fixed(1.42);              label("cilazaprilat kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr_cilat            <- fixed(1);                 label("cilazaprilat blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 1 (footnote e, assumed)
    fu_p_cilat           <- fixed(0.76);              label("cilazaprilat fraction unbound in plasma (unitless)")  # Section 3.1.3 free fraction in plasma 0.76 [ref 50]; Table 2 fu,b 0.76 with Rb 1 via Eq S8

    # ================= Residual error =================
    # Luo 2024 report no residual-error model: the semi-PBPK model is a
    # forward predictor assessed by 0.5-2-fold AUC/Cmax agreement and
    # 5th-95th percentile coverage (Section 2.4), not by a fitted sigma.
    # Non-fitted placeholders per the in-repo PBPK convention
    # (An_2012_mitoxantrone_human_pbpk.R, Gaohua_2012_pregnancy_pbpk_midazolam.R).
    propSd               <- fixed(0.10);  label("Proportional residual error placeholder, cilazapril (fraction)")  # not reported by Luo 2024
    propSd_cilat         <- fixed(0.10);  label("Proportional residual error placeholder, cilazaprilat (fraction)")  # not reported by Luo 2024
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
    #    cilazapril binds mainly to albumin (Section 2.3).
    # ============================================================
    fu_p_i   <- 1 / (1 + (1 - fu_p) * r_alb / fu_p)
    fu_b_i   <- fu_p_i / bpr
    v_sys    <- exp(lvc) * (fu_p_i / fu_p)
    cl_int_h <- exp(lcl_int_h) * f_ces1 * f_liver
    cl_int_k <- exp(lcl_int_k) * r_gfr
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Table 2 reports ka directly for this drug rather than Peff, so Eq S3 is
    # not used; the single ka applies to all three segments and carries the
    # Eq 6 intestinal-permeability correction in cirrhosis.
    ka_duo <- exp(lka) * r_lr
    ka_jej <- ka_duo
    ka_ile <- ka_duo

    fu_p_i_cilat  <- 1 / (1 + (1 - fu_p_cilat) * r_alb / fu_p_cilat)
    fu_b_i_cilat  <- fu_p_i_cilat / bpr_cilat
    v_sys_cilat   <- exp(lvc_cilat) * (fu_p_i_cilat / fu_p_cilat)
    cl_int_k_cilat <- exp(lcl_int_k_cilat) * r_gfr
    k12_cilat <- exp(lk12_cilat)
    k21_cilat <- exp(lk21_cilat)

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

    c_wall_duo_cilat <- wall_duodenum_cilat / v_duo
    c_wall_jej_cilat <- wall_jejunum_cilat  / v_jej
    c_wall_ile_cilat <- wall_ileum_cilat    / v_ile
    c_pv_cilat       <- portal_vein_cilat   / v_pv
    c_liver_cilat    <- liver_cilat         / v_liver
    c_kid_cilat      <- kidney_cilat        / v_kid
    c_sys_cilat      <- central_cilat       / v_sys_cilat
    eff_wall_duo_cilat <- c_wall_duo_cilat * bpr_cilat / kp_gut_cilat
    eff_wall_jej_cilat <- c_wall_jej_cilat * bpr_cilat / kp_gut_cilat
    eff_wall_ile_cilat <- c_wall_ile_cilat * bpr_cilat / kp_gut_cilat
    eff_liver_cilat    <- c_liver_cilat    * bpr_cilat / kp_liver_cilat
    eff_kid_cilat      <- c_kid_cilat      * bpr_cilat / kp_kidney_cilat

    # Hepatic CES1 hydrolysis flux (mg/min): Eq S7 loss term for the parent
    # and Eq S13 formation term for the metabolite.
    hydrolysis <- cl_int_h * fu_b_i * eff_liver

    # ============================================================
    # 6. ODE system (Luo 2024 Supplementary Material).
    #    Parent: Eq S1 stomach, S2 lumen, S5 gut wall, S6 portal vein,
    #    S7 liver, S14 kidney, S16-S17 (2-cmt) systemic.
    #    cilazaprilat: S5/S6 without absorption, S13 liver, S14 kidney, S16-S17 (2-cmt).
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
                     + k21 * peripheral1 - k12 * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- cilazaprilat: no lumen states (it is not absorbed from the gut lumen) ----
    d/dt(wall_duodenum_cilat) <- q_duo * c_sys_cilat - q_duo * eff_wall_duo_cilat
    d/dt(wall_jejunum_cilat)  <- q_jej * c_sys_cilat - q_jej * eff_wall_jej_cilat
    d/dt(wall_ileum_cilat)    <- q_ile * c_sys_cilat - q_ile * eff_wall_ile_cilat

    d/dt(portal_vein_cilat) <- q_duo * eff_wall_duo_cilat + q_jej * eff_wall_jej_cilat +
                            q_ile * eff_wall_ile_cilat - q_pv * c_pv_cilat

    d/dt(liver_cilat) <- q_pv * c_pv_cilat + q_ha * c_sys_cilat + hydrolysis -
                      q_liver * eff_liver_cilat

    d/dt(kidney_cilat) <- q_kid * c_sys_cilat -
                       (q_kid + fu_b_i_cilat * cl_int_k_cilat) * eff_kid_cilat

    d/dt(central_cilat) <- q_liver * eff_liver_cilat + q_kid * eff_kid_cilat -
                        (q_ha + q_kid) * c_sys_cilat - q_gut * c_sys_cilat
                        + k21_cilat * peripheral1_cilat - k12_cilat * central_cilat
    d/dt(peripheral1_cilat) <- k12_cilat * central_cilat - k21_cilat * peripheral1_cilat

    # ============================================================
    # 7. Observations. Amounts are mg and volumes L, so concentrations are
    #    mg/L = ug/mL, matching the AUC unit (ug*h/mL) of Table 6.
    # ============================================================
    Cc <- c_sys
    Cc_cilat <- c_sys_cilat

    Cc ~ prop(propSd)
    Cc_cilat ~ prop(propSd_cilat)
  })
}

