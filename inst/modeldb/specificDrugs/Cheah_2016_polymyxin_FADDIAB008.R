Cheah_2016_polymyxin_FADDIAB008 <- function() {
  description <- "In vitro (Acinetobacter baumannii FADDI-AB008; clinical heteroresistant isolate; loss-of-LPS resistance mechanism). Mechanism-based PK/PD model for polymyxin B and colistin against A. baumannii in a dynamic one-compartment in vitro infection model (IVM). Identical structural model to Cheah_2016_polymyxin_ATCC19606 (three bacterial subpopulations bact_s/bact_r/bact_d with logistic carrying capacity, Bulitta 2010 lipid-A receptor-occupancy submodel, Hill-function killing, single-compartment adaptive-resistance turnover, IVM one-compartment PK); strain-specific parameter values from Cheah 2016 Table 1 column FADDI-AB008. For this strain the experimental data supported inclusion of a fitness cost f_cost on susceptible replication (Eq 10): G_inhib_max = 0.994 (Table 1). SC50 is inherited as a FIXED proxy from Bulitta 2015 (Cheah 2016 does not report it); see vignette Errata."
  reference <- paste(
    "Cheah S-E, Li J, Tsuji BT, Forrest A, Bulitta JB, Nation RL. (2016).",
    "Colistin and polymyxin B dosage regimens against Acinetobacter baumannii:",
    "differences in activity and the emergence of resistance.",
    "Antimicrobial Agents and Chemotherapy 60(7):3921-3933.",
    "doi:10.1128/AAC.02927-15.",
    sep = " "
  )
  vignette <- "Cheah_2016_polymyxin_Abaumannii_dynamics"
  units <- list(
    time          = "hour",
    dosing        = "mg (polymyxin B or colistin base, IV bolus or 1-h infusion into IVM central reservoir)",
    concentration = "log10 CFU/mL (Cc, observed viable count on drug-free agar)"
  )

  paper_specific_compartments <- c("bact_s", "bact_r", "bact_d", "r_adapt")

  covariateData <- list()

  population <- list(
    species         = "in vitro (Acinetobacter baumannii FADDI-AB008)",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    organism        = "A. baumannii FADDI-AB008 (clinical heteroresistant isolate; described in reference 19 as isolate 8; polymyxin B and colistin MIC 0.5 mg/L; polymyxin resistance via loss of lipopolysaccharide from the outer membrane)",
    system          = "Dynamic one-compartment in vitro infection model (IVM), 80 mL central reservoir, 37 C, CAMHB circulated at 4.8 mL/h to simulate elimination half-life 11.6 h and average steady-state polymyxin concentration 3 mg/L",
    duration        = "96 h with viable counting at 0, 0.5, 1, 2, 4, 8, 11, 13, 23, 25, 26, 28, 47, 49, 50, 52, 71, 73, 74, 76, and 96 h on drug-free and drug-containing (6.6 mg/L polymyxin B base) plates; PAPs at 0, 23, 47, 71, and 96 h",
    inoculum_target = "approximately 10^6 CFU/mL (Table 1: log10 CFU_total,0 = 5.93)",
    regimens        = "R1 (gradual rise of colistin), R2 (polymyxin B 1-h infusion every 12 h, no loading dose), R3 (R2 + conventional loading dose), R4 (R2 + augmented loading dose attaining initial peak 6 mg/L)",
    notes           = "In-vitro pharmacodynamic study; no human or animal subjects. Random effects (etas) are not included (typical-value fit per strain). For FADDI-AB008 the experimental data supported a fitness cost G_inhib_max = 0.994 (Cheah 2016 Table 1, Results paragraph 5); the F_cost term in Eq 1 slows the replication rate of the susceptible population as R_adaptive rises. Cheah 2016 notes that an augmented polymyxin B loading dose (R4) delayed bacterial regrowth in this strain (Discussion paragraph 5)."
  )

  ini({
    # ===============================================================
    # Bacterial growth (Cheah 2016 Table 1, FADDI-AB008)
    # ===============================================================
    mgt_s <- 35.2
    label("Mean generation time of polymyxin-susceptible bacteria (min; MGT_S)")  # Cheah 2016 Table 1, FADDI-AB008
    mgt_r <- 48.1
    label("Mean generation time of constitutively polymyxin-resistant bacteria (min; MGT_R)")  # Cheah 2016 Table 1, FADDI-AB008

    log10_cfumax <- 7.66
    label("log10 maximal bacterial population (log10 CFU/mL; CFU_max)")  # Cheah 2016 Table 1, FADDI-AB008

    log10_cfu0_total <- 5.93
    label("log10 initial inoculum total viable count (log10 CFU/mL; CFU_total,0)")  # Cheah 2016 Table 1, FADDI-AB008
    log10_cfu0_r <- -0.44
    label("log10 initial constitutively-resistant inoculum (log10 CFU/mL; CFU_R,0)")  # Cheah 2016 Table 1, FADDI-AB008

    k_sd <- 0.00246
    label("Transition rate constant from susceptible into dormant state (1/h; k_SD)")  # Cheah 2016 Table 1, FADDI-AB008
    k_ds <- 5.45
    label("Transition rate constant from dormant back into susceptible state (1/h; k_DS)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # IVM pharmacokinetics (Cheah 2016 Table 1, FADDI-AB008)
    # ===============================================================
    cl_ivm <- 0.00585
    label("IVM flow rate / clearance from central reservoir (L/h; CL)")  # Cheah 2016 Table 1, FADDI-AB008
    v_ivm <- 0.119
    label("IVM central reservoir volume (L; V)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # Polymyxin target-site binding (Eqs 4-5)
    # ===============================================================
    hill_binding <- 3.25
    label("Hill coefficient for polymyxin binding to lipid-A LPS receptor (Hill_binding)")  # Cheah 2016 Table 1, FADDI-AB008
    ec50 <- 0.0237
    label("Fraction of receptors not occupied by Mg2+/Ca2+ at 50% F_polymyxin_eff (unitless; EC50)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # Bacterial killing (Eq 7)
    # ===============================================================
    kill_max <- fixed(100)
    label("Maximal polymyxin-induced killing rate constant (1/h; Kill_max; FIXED per Cheah 2016)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    hill_killing <- 1.75
    label("Hill coefficient for polymyxin-induced killing (Hill_killing)")  # Cheah 2016 Table 1, FADDI-AB008
    killc50 <- 0.693
    label("EC50 of effective polymyxin concentration for killing (mg/L; KillC50)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # Adaptive resistance (Eqs 8-9)
    # ===============================================================
    s_max <- fixed(300)
    label("Maximal fold reduction in effective polymyxin via adaptive resistance (S_max; FIXED per Cheah 2016)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    sc50 <- fixed(36.5)
    label("Adaptive-resistance half-saturation polymyxin concentration (mg/L; SC50; FIXED, inherited from Bulitta 2015)")  # Bulitta JB et al. 2015 AAC 59:2315-2327, Table 1 PAO1-RH (operator-approved fixed-from-class proxy; not in Cheah 2016)
    k_adapt <- 8.09
    label("Adaptation rate constant for R_adaptive turnover (1/h; k_adapt)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # Fitness cost (Eq 10) -- ESTIMATED for FADDI-AB008
    # ===============================================================
    g_inhib_max <- 0.994
    label("Maximal fitness cost associated with adaptive resistance (G_inhib,max; unitless)")  # Cheah 2016 Table 1, FADDI-AB008

    # ===============================================================
    # Medium / physical constants (inherited; see ATCC 19606 file)
    # ===============================================================
    kd_cations <- fixed(200)
    label("Receptor dissociation constant for Mg2+/Ca2+ (umol/L; Kd_cations; FIXED per Bulitta 2010)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Table 1 footnote (g)
    kd_polymyxin <- fixed(0.3)
    label("Receptor dissociation constant for polymyxin (umol/L; Kd_polymyxin; FIXED per Bulitta 2010)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Table 1 footnote (g)
    mw_polymyxin <- fixed(1.163)
    label("Mean polymyxin molar mass (mg/umol; MW_polymyxin; FIXED at 1163 g/mol per Bulitta 2010 colistin reference)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Methods + Table 1 footnote
    c_cations <- fixed(1138)
    label("Sum of Mg2+ and Ca2+ in CAMHB (umol/L; C_cations; FIXED per Bulitta 2010 Table 1 footnote f)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Table 1 footnote (f)

    # ===============================================================
    # Residual error -- not reported in Cheah 2016
    # ===============================================================
    addSd <- fixed(0.01)
    label("Additive residual SD on log10 CFU/mL (FIXED tiny placeholder -- not reported in Cheah 2016; vignette uses zeroRe())")  # Cheah 2016 does NOT report residual error SD
  })

  model({
    kg_s       <- log(2) * 60 / mgt_s
    kg_r       <- log(2) * 60 / mgt_r
    cfumax     <- 10 ^ log10_cfumax
    cfu0_total <- 10 ^ log10_cfu0_total
    cfu0_r     <- 10 ^ log10_cfu0_r
    cfu0_s     <- cfu0_total - cfu0_r
    kel_ivm    <- cl_ivm / v_ivm

    c_polymyxin     <- central / v_ivm
    c_polymyxin_um  <- c_polymyxin / mw_polymyxin
    # Denominator split with parens to work around an rxode2 parser
    # false-positive mu-reference detection that mis-flags
    # 'pop_param + pop_param' as 'THETA + ETA' syntax.
    kd_ratio        <- kd_cations / kd_polymyxin
    denom_baseline  <- (kd_cations) + (c_cations)
    denom_drug      <- kd_ratio * c_polymyxin_um
    denom_total     <- (denom_baseline) + (denom_drug)
    f_bound_cations <- c_cations / denom_total

    not_occ          <- 1 - f_bound_cations
    not_occ_h        <- not_occ ^ hill_binding
    ec50_h           <- ec50 ^ hill_binding
    f_polymyxin_eff  <- not_occ_h / ((ec50_h) + (not_occ_h))
    c_polymyxin_eff  <- f_polymyxin_eff * c_polymyxin / (1 + r_adapt)

    c_eff_hk   <- c_polymyxin_eff ^ hill_killing
    killc50_hk <- killc50 ^ hill_killing
    kill       <- kill_max * c_eff_hk / ((killc50_hk) + (c_eff_hk))

    stim   <- s_max * c_polymyxin / ((sc50) + (c_polymyxin))
    f_cost <- g_inhib_max * r_adapt / s_max

    cfu_total <- bact_s + bact_r

    d/dt(bact_s) <- bact_s * (
        kg_s * (1 - cfu_total / cfumax) * (1 - f_cost) -
        kill -
        k_sd -
        kel_ivm
      ) + k_ds * bact_d
    d/dt(bact_r) <- bact_r * (
        kg_r * (1 - cfu_total / cfumax) -
        kel_ivm
      )
    d/dt(bact_d) <- k_sd * bact_s - ((k_ds) + (kel_ivm)) * bact_d
    d/dt(r_adapt) <- k_adapt * (stim - r_adapt)
    d/dt(central) <- -kel_ivm * central

    bact_s(0)  <- cfu0_s
    bact_r(0)  <- cfu0_r
    bact_d(0)  <- 0
    r_adapt(0) <- 0
    central(0) <- 0

    cfu_obs <- bact_s + bact_r + 1
    Cc      <- log10(cfu_obs)
    Cc      ~ add(addSd)
  })
}
