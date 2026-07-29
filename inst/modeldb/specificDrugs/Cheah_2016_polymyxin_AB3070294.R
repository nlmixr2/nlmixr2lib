Cheah_2016_polymyxin_AB3070294 <- function() {
  description <- "In vitro (Acinetobacter baumannii AB307-0294; clinical heteroresistant isolate). Mechanism-based PK/PD model for polymyxin B and colistin against A. baumannii in a dynamic one-compartment in vitro infection model (IVM). Identical structural model to Cheah_2016_polymyxin_ATCC19606 (three bacterial subpopulations bact_s/bact_r/bact_d with logistic carrying capacity, Bulitta 2010 lipid-A receptor-occupancy submodel, Hill-function killing, single-compartment adaptive-resistance turnover, IVM one-compartment PK); strain-specific parameter values from Cheah 2016 Table 1 column AB307-0294. Fitness cost G_inhib_max was not estimated for this strain (Table 1 'NE') so f_cost is held at 0. SC50 is inherited as a FIXED proxy from Bulitta 2015 (Cheah 2016 does not report it); see vignette Errata for the full inheritance list."
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
    species         = "in vitro (Acinetobacter baumannii AB307-0294)",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    organism        = "A. baumannii AB307-0294 (clinical heteroresistant isolate; polymyxin B and colistin MIC 1 mg/L; population-analysis-profile evidence of heteroresistance)",
    system          = "Dynamic one-compartment in vitro infection model (IVM), 80 mL central reservoir, 37 C, cation-adjusted Mueller-Hinton broth (CAMHB) circulated at 4.8 mL/h to simulate elimination half-life 11.6 h and average steady-state polymyxin concentration 3 mg/L",
    duration        = "96 h with viable counting at 0, 0.5, 1, 2, 4, 8, 11, 13, 23, 25, 26, 28, 47, 49, 50, 52, 71, 73, 74, 76, and 96 h on drug-free and drug-containing (6.6 mg/L polymyxin B base) plates; PAPs at 0, 23, 47, 71, and 96 h",
    inoculum_target = "approximately 10^6 CFU/mL (Table 1: log10 CFU_total,0 = 6.14)",
    regimens        = "R1 (gradual rise of colistin), R2 (polymyxin B 1-h infusion every 12 h, no loading dose), R3 (R2 + conventional loading dose), R4 (R2 + augmented loading dose attaining initial peak 6 mg/L)",
    notes           = "In-vitro pharmacodynamic study; no human or animal subjects. Random effects (etas) are not included (typical-value fit per strain). Polymyxin B and colistin share structural binding/killing parameters per strain; only the simulated PK profile (R1 vs R2-R4) differs. Note: Table 1 reports a relatively slow adaptive-resistance rise of bacterial regrowth for this strain (consistent with the high k_DS / low k_DS combination and the moderate k_adapt 14.2 h^-1)."
  )

  ini({
    # ===============================================================
    # Bacterial growth (Cheah 2016 Table 1, AB307-0294)
    # ===============================================================
    mgt_s <- 32
    label("Mean generation time of polymyxin-susceptible bacteria (min; MGT_S)")  # Cheah 2016 Table 1, AB307-0294
    mgt_r <- 81.0
    label("Mean generation time of constitutively polymyxin-resistant bacteria (min; MGT_R)")  # Cheah 2016 Table 1, AB307-0294

    log10_cfumax <- 8.05
    label("log10 maximal bacterial population (log10 CFU/mL; CFU_max)")  # Cheah 2016 Table 1, AB307-0294

    log10_cfu0_total <- 6.14
    label("log10 initial inoculum total viable count (log10 CFU/mL; CFU_total,0)")  # Cheah 2016 Table 1, AB307-0294
    log10_cfu0_r <- -0.503
    label("log10 initial constitutively-resistant inoculum (log10 CFU/mL; CFU_R,0)")  # Cheah 2016 Table 1, AB307-0294

    k_sd <- 0.0367
    label("Transition rate constant from susceptible into dormant state (1/h; k_SD)")  # Cheah 2016 Table 1, AB307-0294
    k_ds <- 0.0332
    label("Transition rate constant from dormant back into susceptible state (1/h; k_DS)")  # Cheah 2016 Table 1, AB307-0294

    # ===============================================================
    # IVM pharmacokinetics (Cheah 2016 Table 1, AB307-0294)
    # ===============================================================
    cl_ivm <- 0.00506
    label("IVM flow rate / clearance from central reservoir (L/h; CL)")  # Cheah 2016 Table 1, AB307-0294
    v_ivm <- 0.101
    label("IVM central reservoir volume (L; V)")  # Cheah 2016 Table 1, AB307-0294

    # ===============================================================
    # Polymyxin target-site binding (Eqs 4-5)
    # ===============================================================
    hill_binding <- 4.13
    label("Hill coefficient for polymyxin binding to lipid-A LPS receptor (Hill_binding)")  # Cheah 2016 Table 1, AB307-0294
    ec50 <- 0.0989
    label("Fraction of receptors not occupied by Mg2+/Ca2+ at 50% F_polymyxin_eff (unitless; EC50)")  # Cheah 2016 Table 1, AB307-0294

    # ===============================================================
    # Bacterial killing (Eq 7)
    # ===============================================================
    kill_max <- fixed(100)
    label("Maximal polymyxin-induced killing rate constant (1/h; Kill_max; FIXED per Cheah 2016)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    hill_killing <- 2.79
    label("Hill coefficient for polymyxin-induced killing (Hill_killing)")  # Cheah 2016 Table 1, AB307-0294
    killc50 <- 0.561
    label("EC50 of effective polymyxin concentration for killing (mg/L; KillC50)")  # Cheah 2016 Table 1, AB307-0294

    # ===============================================================
    # Adaptive resistance (Eqs 8-9)
    # ===============================================================
    s_max <- fixed(300)
    label("Maximal fold reduction in effective polymyxin via adaptive resistance (S_max; FIXED per Cheah 2016)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    sc50 <- fixed(36.5)
    label("Adaptive-resistance half-saturation polymyxin concentration (mg/L; SC50; FIXED, inherited from Bulitta 2015)")  # Bulitta JB et al. 2015 AAC 59:2315-2327, Table 1 PAO1-RH (operator-approved fixed-from-class proxy; not in Cheah 2016)
    k_adapt <- 14.2
    label("Adaptation rate constant for R_adaptive turnover (1/h; k_adapt)")  # Cheah 2016 Table 1, AB307-0294

    # ===============================================================
    # Fitness cost -- NOT estimated for AB307-0294 (Table 1 NE)
    # ===============================================================
    g_inhib_max <- fixed(0)
    label("Maximal fitness cost (G_inhib,max; FIXED at 0 -- Cheah 2016 Table 1 reports NE for AB307-0294)")  # Cheah 2016 Table 1, AB307-0294 (NE, not estimated)

    # ===============================================================
    # Medium / physical constants (inherited; see ATCC 19606 file for
    # the same provenance)
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
