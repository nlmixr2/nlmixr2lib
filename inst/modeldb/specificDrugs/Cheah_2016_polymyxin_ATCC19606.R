Cheah_2016_polymyxin_ATCC19606 <- function() {
  description <- "In vitro (Acinetobacter baumannii ATCC 19606; heteroresistant reference strain). Mechanism-based PK/PD model for polymyxin B and colistin against A. baumannii in a dynamic one-compartment in vitro infection model (IVM). The bacterial system is partitioned into three subpopulations -- polymyxin-susceptible (bact_s, CFU_S in the paper), constitutively polymyxin-resistant (bact_r, CFU_R; killing rate fixed at zero), and dormant or extremely slowly replicating cells (bact_d, Pop_D; nonobservable on viable-count plates) -- with a logistic carrying capacity CFU_max constraining the total observable population and a first-order bidirectional transition between susceptible and dormant states (k_SD, k_DS). Polymyxin in the IVM reservoir (central compartment) follows one-compartment first-order kinetics (CL_IVM, V_IVM) with simulated elimination half-life 11.6 h. Polymyxin target-site binding follows the Bulitta 2010 lipid-A LPS receptor-occupancy model: competitive displacement of bound divalent cations (Ca2+ and Mg2+) by polymyxin gives F_bound_cations (Eq 4), and a Hill function of the unoccupied fraction (Hill_binding, EC50) gives F_polymyxin_eff (Eq 5). The effective polymyxin concentration C_polymyxin_eff is the F_polymyxin_eff-weighted reservoir concentration divided by (1 + R_adaptive) to encode adaptive resistance attenuation (Eq 6). Bacterial killing is a Hill function of C_polymyxin_eff (Eq 7; Kill_max fixed at 100/h, Hill_killing, KillC50). Adaptive resistance R_adaptive is a single-compartment turnover whose driver Stim is a Hill-1 of raw reservoir polymyxin concentration (Eq 8; S_max fixed at 300, SC50 inherited from Bulitta 2015 since Cheah 2016 does not report it), with rate constant k_adapt (Eq 9). For ATCC 19606 the fitness cost G_inhib_max was not estimated by the authors (Table 1 'NE') so f_cost is held at 0. Observation is the log10 of the drug-free agar viable count CFU_S + CFU_R (Eq 11). The model has no inter-experiment IIV (typical-value fit per strain) and the residual error is set to a tiny fixed value because Cheah 2016 does not report it -- see vignette Errata for the full inheritance / approximation list."
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
    species         = "in vitro (Acinetobacter baumannii ATCC 19606)",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    organism        = "A. baumannii ATCC 19606 (heteroresistant reference strain; polymyxin B and colistin MIC 0.5 mg/L; resistance occurs via lipid-A phosphoethanolamine modification or loss of LPS from the outer membrane)",
    system          = "Dynamic one-compartment in vitro infection model (IVM), 80 mL central reservoir, 37 C, cation-adjusted Mueller-Hinton broth (CAMHB) circulated at 4.8 mL/h to simulate elimination half-life 11.6 h and average steady-state polymyxin concentration 3 mg/L",
    duration        = "96 h with viable counting at 0, 0.5, 1, 2, 4, 8, 11, 13, 23, 25, 26, 28, 47, 49, 50, 52, 71, 73, 74, 76, and 96 h on drug-free and drug-containing (6.6 mg/L polymyxin B base) plates; population analysis profiles (PAPs) at 0, 23, 47, 71, and 96 h on 1.7, 3.3, and 6.6 mg/L polymyxin B agar plates",
    inoculum_target = "approximately 10^6 CFU/mL (Table 1: log10 CFU_total,0 = 6.34)",
    regimens        = "R1 (gradual rise of colistin, mimicking the Plachouras 2009 ref 8 predicted profile for a patient with no loading CMS dose; continuous infusion approximating accumulation to Css 3 mg/L), R2 (polymyxin B, 1-h infusion every 12 h, no loading dose), R3 (R2 plus a conventional loading dose to attain Css 3 mg/L rapidly), R4 (R2 plus an augmented loading dose attaining initial peak 6 mg/L)",
    notes           = "In-vitro pharmacodynamic study; no human or animal subjects. Random effects (etas) are not included: the paper reports per-strain typical-value fits with %SE values that reflect parameter precision, not between-replicate IIV. Polymyxin B and colistin were modelled with identical structural binding/killing parameters per strain -- only the simulated PK profile (R1 vs R2-R4) differs. Heteroresistance is encoded via the constitutively-resistant subpopulation CFU_R, whose initial fraction log10 CFU_R,0 = -0.249 of the total inoculum (~0.56 CFU/mL) is reported in Table 1."
  )

  ini({
    # ===============================================================
    # Bacterial growth (per-strain, Cheah 2016 Table 1, ATCC 19606)
    # ===============================================================
    # Mean generation times (MIN per Table 1); the model() block
    # converts MGT (min) into growth rate kg (1/h) via kg = ln(2) * 60
    # / MGT_min.  MGT is treated as a typical-value structural fit.
    mgt_s <- 75.2
    label("Mean generation time of polymyxin-susceptible bacteria (min; MGT_S)")  # Cheah 2016 Table 1, ATCC 19606
    mgt_r <- 71.4
    label("Mean generation time of constitutively polymyxin-resistant bacteria (min; MGT_R)")  # Cheah 2016 Table 1, ATCC 19606

    # log10 of the logistic carrying capacity (CFU_max in the paper);
    # the model() block exponentiates this to plain CFU/mL.
    log10_cfumax <- 7.35
    label("log10 maximal bacterial population (log10 CFU/mL; CFU_max)")  # Cheah 2016 Table 1, ATCC 19606

    # Initial inoculum -- log10 of total CFU/mL and log10 of CFU_R/mL.
    # The model() block back-computes log10 CFU_S,0 = log10(CFU_total,0
    # - CFU_R,0) (CFU_R,0 << CFU_total,0 so the floor never bites).
    log10_cfu0_total <- 6.34
    label("log10 initial inoculum total viable count (log10 CFU/mL; CFU_total,0)")  # Cheah 2016 Table 1, ATCC 19606
    log10_cfu0_r <- -0.249
    label("log10 initial constitutively-resistant inoculum (log10 CFU/mL; CFU_R,0)")  # Cheah 2016 Table 1, ATCC 19606

    # Susceptible<->dormant transition rate constants (1/h).
    k_sd <- 0.00635
    label("Transition rate constant from susceptible into dormant state (1/h; k_SD)")  # Cheah 2016 Table 1, ATCC 19606
    k_ds <- 3.19
    label("Transition rate constant from dormant back into susceptible state (1/h; k_DS)")  # Cheah 2016 Table 1, ATCC 19606

    # ===============================================================
    # IVM pharmacokinetics (per-strain Table 1; one-compartment,
    # first-order elimination CL_IVM / V_IVM giving t1/2 = 11.6 h)
    # ===============================================================
    cl_ivm <- 0.00434
    label("IVM flow rate / clearance from central reservoir (L/h; CL)")  # Cheah 2016 Table 1, ATCC 19606
    v_ivm <- 0.111
    label("IVM central reservoir volume (L; V)")  # Cheah 2016 Table 1, ATCC 19606

    # ===============================================================
    # Polymyxin target-site (lipid-A LPS) binding (Eqs 4-5)
    # ===============================================================
    hill_binding <- 3.06
    label("Hill coefficient for polymyxin binding to lipid-A LPS receptor (Hill_binding)")  # Cheah 2016 Table 1, ATCC 19606
    # EC50 here is unitless: fraction of receptors NOT occupied by
    # Mg2+/Ca2+ at which F_polymyxin_eff = 0.5 (per Bulitta 2010 ref 30
    # model). Linear-space numbers in Table 1 are 0.0173-0.0989, all <<
    # 1, confirming the unitless / fractional interpretation.
    ec50 <- 0.0483
    label("Fraction of receptors not occupied by Mg2+/Ca2+ at 50% F_polymyxin_eff (unitless; EC50)")  # Cheah 2016 Table 1, ATCC 19606

    # ===============================================================
    # Bacterial killing by polymyxin (Eq 7)
    # ===============================================================
    # Kill_max FIXED at 100 /h per Cheah 2016 Results paragraph 3 ("the
    # maximal bacterial killing rate (Kill_max) was fixed to 100 h-1
    # (equivalent to 11 log10 killing in 15 min)") and Table 1.
    kill_max <- fixed(100)
    label("Maximal polymyxin-induced killing rate constant (1/h; Kill_max; FIXED per Cheah 2016 Results / Table 1)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    hill_killing <- 3.85
    label("Hill coefficient for polymyxin-induced killing (Hill_killing)")  # Cheah 2016 Table 1, ATCC 19606
    killc50 <- 0.521
    label("EC50 of effective polymyxin concentration for killing (mg/L; KillC50)")  # Cheah 2016 Table 1, ATCC 19606

    # ===============================================================
    # Adaptive resistance turnover (Eqs 8-9)
    # ===============================================================
    # S_max FIXED at 300 per Cheah 2016 Results paragraph 3 ("the
    # maximal attenuation in killing (i.e., extent of adaptive
    # resistance; S_max) [was fixed] to 300") and Table 1.
    s_max <- fixed(300)
    label("Maximal fold reduction in effective polymyxin concentration via adaptive resistance (S_max; FIXED per Cheah 2016 Results / Table 1)")  # Cheah 2016 Table 1 (fixed) + Results paragraph 3
    # SC50 is the half-saturation polymyxin concentration in the Stim
    # equation. Cheah 2016 does NOT report SC50 anywhere; Table 1 lists
    # only k_adapt and the fixed S_max. The closest source for SC50 is
    # the same lab's prior Bulitta 2015 mechanism-based adaptive-
    # resistance fit (Bulitta JB et al. 2015. Two mechanisms of killing
    # of Pseudomonas aeruginosa by tobramycin assessed at multiple
    # inocula via mechanism-based modeling. Antimicrob Agents Chemother
    # 59:2315-2327. doi:10.1128/AAC.04099-14), which reports SC50,Adapt
    # = 36.5 mg/L for tobramycin against PAO1-RH. SC50 is inherited as
    # FIXED here; operator-approved fixed-from-class proxy (see vignette
    # Errata).
    sc50 <- fixed(36.5)
    label("Adaptive-resistance half-saturation polymyxin concentration (mg/L; SC50; FIXED, inherited from Bulitta 2015 ref 31 framework -- not reported in Cheah 2016)")  # Bulitta JB et al. 2015 AAC 59:2315-2327, Table 1 PAO1-RH (operator-approved fixed-from-class proxy)
    k_adapt <- 7.20
    label("Adaptation rate constant for R_adaptive turnover (1/h; k_adapt)")  # Cheah 2016 Table 1, ATCC 19606

    # ===============================================================
    # Fitness cost (Eq 10) -- NOT estimated for ATCC 19606 (Table 1 NE)
    # ===============================================================
    # Holds f_cost identically at 0 in the susceptible-growth term of
    # Eq 1. The operator-approved policy for unreported / not-estimated
    # structural shifts is fixed(0) with a clear inline note.
    g_inhib_max <- fixed(0)
    label("Maximal fitness cost associated with adaptive resistance (G_inhib,max; FIXED at 0 -- Cheah 2016 Table 1 reports NE for ATCC 19606)")  # Cheah 2016 Table 1, ATCC 19606 (NE, not estimated)

    # ===============================================================
    # Medium / physical constants (NOT in Cheah 2016; operator-approved
    # inheritance from Bulitta 2010 ref 30 + CAMHB CLSI specification)
    # ===============================================================
    # Kd_cations and Kd_polymyxin are receptor dissociation constants
    # for Mg2+/Ca2+ and polymyxin on lipid-A LPS, fixed per Bulitta 2010
    # Methods Receptor-occupancy paragraph (citing Schindler & Osborn
    # 1979) and the same paper's Table 1 footnote (g).
    kd_cations <- fixed(200)
    label("Receptor dissociation constant for Mg2+/Ca2+ on lipid-A LPS (umol/L; Kd_cations; FIXED per Bulitta 2010 Table 1 footnote g)")  # Bulitta JB et al. 2010 AAC 54:2051-2062 (ref 30 in Cheah 2016), Table 1 footnote (g)
    kd_polymyxin <- fixed(0.3)
    label("Receptor dissociation constant for polymyxin on lipid-A LPS (umol/L; Kd_polymyxin; FIXED per Bulitta 2010 Table 1 footnote g)")  # Bulitta JB et al. 2010 AAC 54:2051-2062 (ref 30 in Cheah 2016), Table 1 footnote (g)
    # Polymyxin molecular mass: polymyxin B (B1+B2 mix, ~1203 g/mol)
    # and colistin (E1+E2 mix, ~1163 g/mol) are within 3%. Inherited as
    # 1163 g/mol = 1.163 mg/umol from Bulitta 2010 (colistin reference)
    # because Cheah 2016 uses identical structural binding parameters
    # for both polymyxins per strain.
    mw_polymyxin <- fixed(1.163)
    label("Mean polymyxin molar mass (mg/umol = g/mmol; MW_polymyxin; FIXED at 1163 g/mol per Bulitta 2010 colistin reference; polymyxin B and colistin within 3%)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Methods paragraph after Eq. 1 + Table 1 footnote
    # Sum of Mg2+ and Ca2+ molar concentration in CAMHB: 0.514 mmol/L
    # Mg2+ + 0.624 mmol/L Ca2+ = 1.138 mmol/L = 1138 umol/L, per
    # Bulitta 2010 Table 1 footnote (f), which describes the same
    # CAMHB used by Cheah 2016 (Methods: "cation-adjusted Mueller-
    # Hinton broth").
    c_cations <- fixed(1138)
    label("Sum of Mg2+ and Ca2+ molar concentration in CAMHB (umol/L; C_cations; FIXED per Bulitta 2010 Table 1 footnote f)")  # Bulitta JB et al. 2010 AAC 54:2051-2062, Table 1 footnote (f)

    # ===============================================================
    # Residual error (NOT reported in Cheah 2016)
    # ===============================================================
    # Cheah 2016 does not report a residual error SD. Set to a tiny
    # fixed value so the model is well-formed for typical-value
    # rxode2 simulation (the vignette uses zeroRe() to suppress
    # residual noise entirely). See vignette Errata.
    addSd <- fixed(0.01)
    label("Additive residual SD on log10 CFU/mL (FIXED tiny placeholder -- not reported in Cheah 2016; vignette uses zeroRe())")  # Cheah 2016 does NOT report residual error SD (typical-value fit per strain)
  })

  model({
    # ---------------------------------------------------------------
    # Derived rate constants and starting populations
    # ---------------------------------------------------------------
    # Growth rate kg (1/h) from mean generation time MGT (min):
    #   kg = ln(2) / MGT_h = ln(2) * 60 / MGT_min.
    kg_s <- log(2) * 60 / mgt_s
    kg_r <- log(2) * 60 / mgt_r

    # Logistic carrying capacity and inoculum components (CFU/mL).
    cfumax     <- 10 ^ log10_cfumax
    cfu0_total <- 10 ^ log10_cfu0_total
    cfu0_r     <- 10 ^ log10_cfu0_r
    cfu0_s     <- cfu0_total - cfu0_r  # CFU_R,0 << CFU_total,0

    # IVM elimination rate constant (1/h).
    kel_ivm <- cl_ivm / v_ivm

    # ---------------------------------------------------------------
    # Polymyxin concentration in the IVM central reservoir (mg/L).
    # ---------------------------------------------------------------
    c_polymyxin <- central / v_ivm

    # ---------------------------------------------------------------
    # Cation displacement at the lipid-A receptor (Eq 4):
    #   F_bound,cations = C_cations / (Kd_cations + C_cations
    #                       + (Kd_cations / Kd_polymyxin) * C_polymyxin / MW_polymyxin)
    # C_polymyxin is mg/L; dividing by MW_polymyxin (mg/umol) gives
    # umol/L so the term is unit-consistent with C_cations (umol/L).
    # The denominator is split onto its own lines, with parentheses
    # around each pair of named parameters in a sum, to work around an
    # rxode2 parser false-positive mu-reference detection that
    # mis-flags 'pop_param + pop_param' as 'THETA + ETA' syntax.
    # ---------------------------------------------------------------
    c_polymyxin_um  <- c_polymyxin / mw_polymyxin
    kd_ratio        <- kd_cations / kd_polymyxin
    denom_baseline  <- (kd_cations) + (c_cations)
    denom_drug      <- kd_ratio * c_polymyxin_um
    denom_total     <- (denom_baseline) + (denom_drug)
    f_bound_cations <- c_cations / denom_total

    # ---------------------------------------------------------------
    # Effective polymyxin concentration (Eqs 5, 6):
    #   F_polymyxin,eff = (1 - F_bound,cations)^H / (EC50^H + (1 - F_bound,cations)^H)
    #   C_polymyxin,eff = F_polymyxin,eff * C_polymyxin / (1 + R_adaptive)
    # ---------------------------------------------------------------
    not_occ          <- 1 - f_bound_cations
    not_occ_h        <- not_occ ^ hill_binding
    ec50_h           <- ec50 ^ hill_binding
    f_polymyxin_eff  <- not_occ_h / ((ec50_h) + (not_occ_h))
    c_polymyxin_eff  <- f_polymyxin_eff * c_polymyxin / (1 + r_adapt)

    # ---------------------------------------------------------------
    # Bacterial killing (Eq 7):
    #   Kill = Kill_max * C_eff^Hk / (KillC50^Hk + C_eff^Hk)
    # ---------------------------------------------------------------
    c_eff_hk  <- c_polymyxin_eff ^ hill_killing
    killc50_hk <- killc50 ^ hill_killing
    kill      <- kill_max * c_eff_hk / ((killc50_hk) + (c_eff_hk))

    # ---------------------------------------------------------------
    # Adaptive-resistance turnover (Eqs 8, 9):
    #   Stim       = S_max * C_polymyxin / (SC50 + C_polymyxin)
    #   d/dt(r_adapt) = k_adapt * (Stim - r_adapt)
    # Stim uses RAW reservoir C_polymyxin per the printed equation
    # (not effective C); the paper's prose calls it "accounting for
    # effective polymyxin concentration (Stim)" but the printed Eq 8
    # uses C_polymyxin -- trust the equation per the operator policy.
    # ---------------------------------------------------------------
    stim <- s_max * c_polymyxin / ((sc50) + (c_polymyxin))

    # ---------------------------------------------------------------
    # Fitness cost (Eq 10): only active when G_inhib_max > 0.
    # For ATCC 19606, G_inhib_max is fixed at 0 (Table 1 NE).
    # ---------------------------------------------------------------
    f_cost <- g_inhib_max * r_adapt / s_max

    # ---------------------------------------------------------------
    # Total OBSERVABLE viable count CFU_total (Eq 11 = CFU_S + CFU_R)
    # used for the logistic growth constraint per Eqs 1-2. The dormant
    # population is non-replicating and non-observable and does NOT
    # contribute to the logistic constraint.
    # ---------------------------------------------------------------
    cfu_total <- bact_s + bact_r

    # ---------------------------------------------------------------
    # ODE system (Eqs 1-3 + Eq 9 + one-compartment IVM PK).
    # ---------------------------------------------------------------
    # Eq 1: dCFU_S/dt = CFU_S * [k_growth,S * (1 - CFU_total/CFU_max)
    #                              * (1 - F_cost) - Kill - k_SD - CL/V]
    #                  + k_DS * Pop_D
    d/dt(bact_s) <- bact_s * (
        kg_s * (1 - cfu_total / cfumax) * (1 - f_cost) -
        kill -
        k_sd -
        kel_ivm
      ) + k_ds * bact_d

    # Eq 2: dCFU_R/dt = CFU_R * [k_growth,R * (1 - CFU_total/CFU_max) - CL/V]
    d/dt(bact_r) <- bact_r * (
        kg_r * (1 - cfu_total / cfumax) -
        kel_ivm
      )

    # Eq 3: dPop_D/dt = k_SD * CFU_S - (k_DS + CL/V) * Pop_D
    d/dt(bact_d) <- k_sd * bact_s - ((k_ds) + (kel_ivm)) * bact_d

    # Eq 9: dR_adaptive/dt = k_adapt * (Stim - R_adaptive)
    d/dt(r_adapt) <- k_adapt * (stim - r_adapt)

    # IVM PK: one-compartment first-order elimination from central
    # reservoir; doses (R1-R4) are supplied via the event table.
    d/dt(central) <- -kel_ivm * central

    # ---------------------------------------------------------------
    # Initial conditions
    # ---------------------------------------------------------------
    bact_s(0)  <- cfu0_s
    bact_r(0)  <- cfu0_r
    bact_d(0)  <- 0
    r_adapt(0) <- 0
    central(0) <- 0

    # ---------------------------------------------------------------
    # Observation (Eq 11: drug-free agar viable count = CFU_S + CFU_R).
    # Add a 1 CFU/mL floor to keep log10 finite when killing drives
    # CFU very low (Bulitta 2010 / Wicha 2017 / Landersdorfer 2018
    # convention for in-vitro PD).
    # ---------------------------------------------------------------
    cfu_obs <- bact_s + bact_r + 1
    Cc      <- log10(cfu_obs)
    Cc      ~ add(addSd)
  })
}
