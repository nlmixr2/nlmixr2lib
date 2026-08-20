Luo_2024_remimazolam_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, custom WinNonlin 8.1 implementation).",
    "remimazolam disposition in healthy adults and in liver cirrhosis",
    "(Child-Pugh A/B/C). Remimazolam is a direct-acting ultrashort-acting",
    "sedative that CES1 INACTIVATES to a carboxylic-acid metabolite, so no",
    "active metabolite is tracked. It is administered intravenously only,",
    "and Table 2 reports no CLint,K, so the kidney compartment carries",
    "distribution but no elimination. The semi-PBPK circuit is three",
    "gut-wall segments, portal vein, liver, kidney and a three-compartment",
    "systemic compartment. Cirrhosis is applied by switching the",
    "Child-Pugh-specific physiology of Table 1 (organ blood flows,",
    "functional liver volume, GI transit rates, GFR, plasma-binding protein",
    "concentrations, hepatic CES1 content) and rescaling CLint, CLint,K,",
    "Peff and Vsys through Eq 1-6. Deterministic: the paper's virtual",
    "populations are uniform 80-120% draws on the drug parameters, not",
    "lognormal random effects, so no IIV is encoded."
  )
  reference <- paste(
    "Luo X, Zhang Z, Mu R, Hu G, Liu L, Liu X (2024). Simultaneously",
    "Predicting the Pharmacokinetics of CES1-Metabolized Drugs and Their",
    "Metabolites Using Physiologically Based Pharmacokinetic Model in",
    "Cirrhosis Subjects. Pharmaceutics 16(2):234.",
    "doi:10.3390/pharmaceutics16020234.",
    "ODE system from Supplementary Material Eq S1-S20; drug parameters",
    "from Table 2; system physiology and Child-Pugh scaling from Table 1",
    "and Eq 1-6; validation targets from Table 12.",
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
    "wall_duodenum", "wall_jejunum", "wall_ileum", "portal_vein"
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
    wall_duodenum            = list(analyte = "remimazolam", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum             = list(analyte = "remimazolam", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum               = list(analyte = "remimazolam", units = "mg", specimen = "tissue", verified = TRUE),
    portal_vein              = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    liver                    = list(analyte = "remimazolam", units = "mg", specimen = "tissue", verified = TRUE),
    kidney                   = list(analyte = "remimazolam", units = "mg", specimen = "tissue", verified = TRUE),
    central                  = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1              = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2              = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 2L,
    age_range      = "adults",
    disease_state  = paste(
      "Two clinical reports. Healthy: Sheng 2020 single-ascending-dose and",
      "continuous-infusion cohorts in healthy Chinese volunteers (n = 3 to 10",
      "per dose level). Cirrhosis: Stohr 2021 Child-Pugh B (n = 8) and",
      "Child-Pugh C (n = 3)."
    ),
    dose_range     = "Remimazolam besylate 0.05-0.4 mg/kg and 3.315-24.6 mg intravenous doses",
    notes          = paste(
      "Luo 2024 Table 3. Literature-digitised clinical data; the authors",
      "simulated 1000 virtual individuals per population by drawing CLint,",
      "CLint,K, fu,b, Vsys, Peff, ka, KL:P, KG:P and KK:P uniformly over",
      "80-120% of the Table 2 point values (Section 2.3). Validation targets",
      "are in Table 12."
    )
  )

  ini({
    # ================= REMIMAZOLAM =================
    # Luo 2024 Table 2. Clearances converted mL/min -> L/min.
    lcl_int_h            <- fixed(log(79.21296));     label("remimazolam hepatic CES1 intrinsic clearance (log L/min)")  # Table 2: CLint 79,212.96 mL/min (footnote c, recalculated from CL_L,b); CLb 1180 mL/min [ref 75]
    lvc                  <- fixed(log(15.0768));      label("remimazolam systemic apparent volume of distribution (log L)")  # Table 2: Vsys 15.0768 L [ref 76] (footnote f)
    lk12                 <- fixed(log(0.01638));      label("remimazolam central -> peripheral1 rate K12 (log 1/min)")  # Table 2: K12 0.01638 1/min [ref 76] (footnote f)
    lk21                 <- fixed(log(0.000476));     label("remimazolam peripheral1 -> central rate K21 (log 1/min)")  # Table 2: K21 0.000476 1/min [ref 76] (footnote f)
    lk13                 <- fixed(log(0.3117));       label("remimazolam central -> peripheral2 rate K13 (log 1/min)")  # Table 2: K13 0.3117 1/min [ref 76] (footnote f)
    lk31                 <- fixed(log(0.5057));       label("remimazolam peripheral2 -> central rate K31 (log 1/min)")  # Table 2: K31 0.5057 1/min [ref 76] (footnote f)
    kp_liver             <- fixed(36.34);             label("remimazolam liver:plasma partition coefficient KL:P (unitless)")  # Table 2 KL:P (Rodgers-Rowland, footnote d)
    kp_gut               <- fixed(63.19);             label("remimazolam gut:plasma partition coefficient KG:P (unitless)")  # Table 2 KG:P (Rodgers-Rowland, footnote d)
    kp_kidney            <- fixed(31.2);              label("remimazolam kidney:plasma partition coefficient KK:P (unitless)")  # Table 2 KK:P (Rodgers-Rowland, footnote d)
    bpr                  <- fixed(1);                 label("remimazolam blood:plasma concentration ratio Rb (unitless)")  # Table 2: Rb 1 (footnote e, assumed)
    fu_p                 <- fixed(0.08);              label("remimazolam fraction unbound in plasma (unitless)")  # Section 3.1.9 plasma protein binding ~92% -> fu,p 0.08 [ref 77]; Table 2 fu,b 0.08 with Rb 1 via Eq S8

    # ================= Residual error =================
    # Luo 2024 report no residual-error model: the semi-PBPK model is a
    # forward predictor assessed by 0.5-2-fold AUC/Cmax agreement and
    # 5th-95th percentile coverage (Section 2.4), not by a fitted sigma.
    # Non-fitted placeholders per the in-repo PBPK convention
    # (An_2012_mitoxantrone_human_pbpk.R, Gaohua_2012_pregnancy_pbpk_midazolam.R).
    propSd               <- fixed(0.10);  label("Proportional residual error placeholder, remimazolam (fraction)")  # not reported by Luo 2024
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

    # ============================================================
    # 4. Drug parameters after cirrhosis scaling.
    #    Eq 3: fu,p,CI = 1 / (1 + (1 - fu,p,HT) * Pratio / fu,p,HT)
    #    Eq S8: fu,b = fu,p / Rb          Eq 4: Vsys,CI = (fu,p,CI/fu,p,HT) * Vsys
    #    Eq 1: CLint,CES1 * fCES1 * fliver   Eq 5: CLint,K * GFR ratio
    #    Eq 6: Peff,CI = Peff,HT * LR ratio
    #    remimazolam binds mainly to albumin (Section 2.3).
    # ============================================================
    fu_p_i   <- 1 / (1 + (1 - fu_p) * r_alb / fu_p)
    fu_b_i   <- fu_p_i / bpr
    v_sys    <- exp(lvc) * (fu_p_i / fu_p)
    cl_int_h <- exp(lcl_int_h) * f_ces1 * f_liver
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)

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

    # Hepatic CES1 hydrolysis flux (mg/min): Eq S7 loss term for the parent
    # only; CES1 inactivates this drug so no active metabolite is formed.
    hydrolysis <- cl_int_h * fu_b_i * eff_liver

    # ============================================================
    # 6. ODE system (Luo 2024 Supplementary Material).
    #    Parent: Eq S5 gut wall, S6 portal vein,
    #    S7 liver, S14 kidney, S18-S20 (3-cmt) systemic.
    # ============================================================
    d/dt(wall_duodenum) <- q_duo * c_sys - q_duo * eff_wall_duo
    d/dt(wall_jejunum)  <- q_jej * c_sys - q_jej * eff_wall_jej
    d/dt(wall_ileum)    <- q_ile * c_sys - q_ile * eff_wall_ile

    d/dt(portal_vein) <- q_duo * eff_wall_duo + q_jej * eff_wall_jej +
                         q_ile * eff_wall_ile - q_pv * c_pv

    d/dt(liver) <- q_pv * c_pv + q_ha * c_sys - q_liver * eff_liver - hydrolysis

    # Table 2 reports no CLint,K for this drug: the kidney distributes but does not eliminate.
    d/dt(kidney) <- q_kid * c_sys - q_kid * eff_kid

    d/dt(central) <- q_liver * eff_liver + q_kid * eff_kid -
                     (q_ha + q_kid) * c_sys - q_gut * c_sys
                     + k21 * peripheral1 - k12 * central
                     + k31 * peripheral2 - k13 * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ============================================================
    # 7. Observations. Amounts are mg and volumes L, so concentrations are
    #    mg/L = ug/mL, matching the AUC unit (ug*h/mL) of Table 12.
    # ============================================================
    Cc <- c_sys

    Cc ~ prop(propSd)
  })
}

