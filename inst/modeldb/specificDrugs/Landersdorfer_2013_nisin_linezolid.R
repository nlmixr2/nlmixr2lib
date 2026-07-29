Landersdorfer_2013_nisin_linezolid <- function() {
  description <- "In vitro (methicillin-resistant Staphylococcus aureus USA300). Mechanism-based pharmacodynamic model for nisin plus linezolid in a 48-h static-concentration time-kill assay (S-ADAPT and NONMEM analyses; subpopulation synergy concept). Three pre-existing bacterial populations (Nis-S/Lin-S, Nis-I/Lin-S, Nis-R/Lin-I), each following a Bulitta two-state life-cycle growth model. Nisin kills with a second-order function (k2*Cnis); linezolid inhibits protein synthesis (turnover of a protein pool P), which (i) raises Inh_Rep = 1 - P and therefore reduces successful replication for the Lin-S populations and (ii) inhibits the slow state 1 -> state 2 growth-rate transition via a steep Hill function (Inh_k12) in all three populations. The Nis-R/Lin-I population has no Inh_Rep effect (only Inh_k12). A saturating carrying-capacity replication factor (REP = 2*CFUmax/(CFUmax + CFUall)) caps the population. Nisin and linezolid concentrations are external time-varying inputs (covariates Cnis and Clin); the model contains no human PK component."
  reference <- paste(
    "Landersdorfer CB, Ly NS, Xu H, Tsuji BT, Bulitta JB. (2013).",
    "Quantifying subpopulation synergy for antibiotic combinations via mechanism-based modeling and a sequential dosing design.",
    "Antimicrobial Agents and Chemotherapy 57(5):2343-2351.",
    "doi:10.1128/AAC.00092-13.",
    "Parameter values are the S-ADAPT estimates of Table 1 (nisin + linezolid column);",
    "the NONMEM estimates in the adjacent column are not used here.",
    "Structural equations are Eqs 1, 3, and 5-9 of the main text.",
    sep = " "
  )
  vignette <- "Landersdorfer_2013_nisin_amikacin_linezolid"
  units <- list(time = "hour", dosing = "mg/L (drug input concentration)", concentration = "log10 CFU/mL (observation); mg/L (drug covariates)")

  # Cnis / Clin are the experimentally-controlled nisin and linezolid
  # broth concentrations in the static time-kill (and sequential
  # combination) experiments. They are supplied as time-varying
  # covariates from the event data, not estimated PK profiles.
  depends <- c("Cnis", "Clin")

  # prot_pool is a paper-mechanistic compartment (the protein
  # constituent pool P; Eq 6) that linezolid depletes; depletion of
  # the pool then raises Inh_Rep (Eq 7) and reduces successful
  # replication for the Lin-S populations.
  paper_specific_compartments <- c("prot_pool")

  covariateData <- list(
    Cnis = list(
      description        = "Nisin concentration in the broth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally. Static time-kill: held at 4, 8, 32, or 128 mg/L for monotherapy and at 8, 16, or 32 mg/L for simultaneous combinations with linezolid. Sequential combinations: held at 8 or 32 mg/L for the 1.5-h pretreatment then set to 0 at ~1.75 h. In-vitro experimental input -- not in inst/references/covariate-columns.md (the canonical register is for human pop-PK covariates and does not apply to this in-vitro PD model).",
      source_name        = "Nisin concentration (paper Methods, Combinations and Time-kill experiments)"
    ),
    Clin = list(
      description        = "Linezolid concentration in the broth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally. Static time-kill: held at 2, 8, or 32 mg/L for monotherapy. Simultaneous combinations: 2, 8, or 32 mg/L throughout. Sequential combinations: 0 during the 1.5-h nisin pretreatment then 8 or 32 mg/L after resuspension at ~1.75 h. In-vitro experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Linezolid concentration (paper Methods, Combinations and Time-kill experiments)"
    )
  )

  population <- list(
    species             = "in vitro (methicillin-resistant Staphylococcus aureus, USA300 strain from the Network on Antimicrobial Resistance)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "MRSA USA300 (nisin MIC 16 mg/L; linezolid MIC 2 mg/L)",
    system              = "Static-concentration time-kill assay in 20-mL conical tubes; serial viable counts over 48 h",
    medium              = "Luria-Bertani broth supplemented with 12.5 mg/L Mg2+ and 25 mg/L Ca2+; viable counts on Luria-Bertani agar",
    temperature         = "37 C",
    duration            = "48 h",
    inoculum            = "~10^7.7 CFU/mL diluted from a ~10^9 CFU/mL spectrophotometric saline suspension, then grown for 60 min to ~10^8.0 CFU/mL before antibiotic dosing",
    mic_values          = c(nisin = "16 mg/L", linezolid = "2 mg/L"),
    regimens            = "Nisin monotherapy at 4, 8, 32, and 128 mg/L; linezolid monotherapy at 2, 8, and 32 mg/L; nisin pretreatment (8 or 32 mg/L for 1.5 h) followed by removal and resuspension in 0/8/32 mg/L linezolid; simultaneous combinations of nisin (16 or 32 mg/L) with linezolid (2, 8, or 32 mg/L); 20 viable-count profiles total",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. Fitted simultaneously in S-ADAPT (importance-sampling Monte Carlo expectation-maximization, pmethod=4) and NONMEM (FOCE+I, ADVAN9); r = 0.99 and slope = 1.00 for observed vs individually fitted log10 viable counts (S-ADAPT). The S-ADAPT analysis included biological between-curve variability fixed to CV = 15% (linear-scale parameters), CV = 10% (Hill coefficient), and variance = 0.25 (log10-scale parameters) per Table 1 footnotes b-c; the packaged model omits between-curve variability and is intended for typical-value simulation only. See Landersdorfer 2013 Methods (page 2344) and Table 1."
  )

  ini({
    # =============================================================
    # Bacterial life-cycle parameters (Bulitta two-state LCGM)
    # =============================================================
    # Replication rate constant k21 (state 2 -> 2*state 1 doubling step),
    # shared across subpopulations and fixed per Table 1.
    lk21 <- fixed(log(50))
    label("Replication rate constant k21 (1/h; state 2 -> 2*state 1 doubling; FIXED)")  # Landersdorfer 2013 Table 1, "k21 = 50 (fixed)"

    # Mean generation time for the slow state 1 -> state 2 transition.
    # Table 1 row "Nis_s/Lin_s and Nis_i/Lin_s = 83.9 min" (the two
    # Lin-S subpopulations share the base growth rate). The
    # Nis-R/Lin-I population carries reduced biofitness via fk12_r.
    mgt_base <- 83.9
    label("Mean generation time, Nis-S/Lin-S and Nis-I/Lin-S subpopulations (min)")  # Landersdorfer 2013 Table 1

    # Initial total inoculum and carrying capacity (log10 CFU/mL).
    log10cfu0   <- 7.88
    label("Initial total inoculum (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (Log10 initial inoculum)
    log10cfumax <- 9.38
    label("Maximum population size (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (Log10 CFUmax)

    # =============================================================
    # Initial subpopulation sizes (log10 mutation frequencies)
    # =============================================================
    # Per Methods (Initial conditions): the product of total inoculum
    # and mutation frequency seeds each less-susceptible population at
    # t = 0; the Nis-S/Lin-S population is the remainder.
    log10mf_is <- -2.88
    label("log10 mutation frequency for Nis-I/Lin-S subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1
    log10mf_ri <- -4.15
    label("log10 mutation frequency for Nis-R/Lin-I subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1

    # =============================================================
    # Biofitness reduction factor (multiply k12_base) for Nis-R/Lin-I
    # =============================================================
    fk12_r <- 0.144
    label("Growth-rate factor for Nis-R/Lin-I subpopulation (unitless)")  # Landersdorfer 2013 Table 1 (S-ADAPT)

    # =============================================================
    # Nisin killing -- second-order in nisin concentration
    # =============================================================
    # Note: the nisin killing-rate constants are RE-ESTIMATED in the
    # linezolid combination model and differ from the nisin + amikacin
    # column of Table 1 (k2S 4.49 vs 5.67, k2I 0.0209 vs 0.0664, etc.);
    # this is intrinsic to the two independent fits described in the
    # paper Methods.
    lk2s <- log(4.49)
    label("Second-order nisin killing rate constant, Nis-S populations (L/(mg*h))")  # Landersdorfer 2013 Table 1
    lk2i <- log(0.0209)
    label("Second-order nisin killing rate constant, Nis-I populations (L/(mg*h))")  # Landersdorfer 2013 Table 1
    lk2r <- log(0.00318)
    label("Second-order nisin killing rate constant, Nis-R populations (L/(mg*h))")  # Landersdorfer 2013 Table 1

    # =============================================================
    # Linezolid inhibition of successful replication (Inh_Rep)
    # via depletion of the protein pool P (Eqs 6, 7)
    # =============================================================
    # P obeys dP/dt = k_Prot * [(1 - C_Lin/(C_Lin + IC_50_Prot)) - P]
    # with P(0) = 1 (baseline 100%). Linezolid depletes P; Inh_Rep =
    # Imax_Rep * (1 - P) raises the probability of unsuccessful
    # replication. Imax_Rep is fixed to 1.0 (Table 1, footnote).
    imax_rep    <- fixed(1.0)
    label("Maximum inhibition of successful replication by linezolid (unitless; FIXED)")  # Landersdorfer 2013 Table 1
    ic50_prot   <- 3.92
    label("Linezolid concn for half-maximum inhibition of protein synthesis (mg/L)")  # Landersdorfer 2013 Table 1 (IC_50,Prot)
    lk_prot     <- log(0.72)
    label("Turnover rate constant for protein pool P (1/h)")  # Landersdorfer 2013 Table 1 (k_Prot)

    # =============================================================
    # Linezolid inhibition of the slow growth-rate transition (Inh_k12)
    # =============================================================
    # Steep Hill function (Hill_k12 = 10 fixed) acting on Cnis via
    # Eq 8: Inh_k12 = Imax_k12 * Clin^Hill / (Clin^Hill + IC50_k12^Hill).
    # Applied as the multiplicative growth-rate-reduction (1 - Inh_k12)
    # on k12 -- see the Assumptions section of the validation vignette
    # for discussion of the parallel-structure reading of Eq 5.
    imax_k12    <- 0.858
    label("Maximum inhibition of state 1 -> state 2 growth rate by linezolid (unitless)")  # Landersdorfer 2013 Table 1
    ic50_k12    <- 4.25
    label("Linezolid concn for 50% of Imax_k12 effect (mg/L)")  # Landersdorfer 2013 Table 1 (IC_50,k12)
    lhill_k12   <- fixed(log(10))
    label("Hill coefficient for the Inh_k12 function (unitless; FIXED)")  # Landersdorfer 2013 Table 1, "Hill_k12 = 10 (fixed)"

    # =============================================================
    # Residual error (additive on log10 viable count)
    # =============================================================
    # Per Methods (Observation model): all log10 viable counts fitted
    # with an additive residual error on the log10 scale. For counts
    # below 100 CFU/mL the paper used a per-colony-count residual
    # model (ref 19); only the additive log10 component is packaged
    # here, matching the published Table 1 value.
    addSd <- 0.307
    label("Additive residual SD on log10 total viable count (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (SD_CFU)
  })

  model({
    # 0. Back-transform log-scale parameters.
    k21       <- exp(lk21)
    k2s       <- exp(lk2s)
    k2i       <- exp(lk2i)
    k2r       <- exp(lk2r)
    k_prot    <- exp(lk_prot)
    hill_k12  <- exp(lhill_k12)

    # 1. Base growth rate k12 (1/h) from mean generation time (min).
    k12_base <- 60 / mgt_base

    # 2. Linezolid Inh_k12 inhibition factor (Eq 8): steep Hill on
    #    linezolid concentration. Effective growth rate is reduced by
    #    multiplying k12 by (1 - Inh_k12); see vignette Assumptions
    #    for the reading of Eq 5.
    clin_h    <- Clin^hill_k12
    inh_k12   <- imax_k12 * clin_h / (clin_h + ic50_k12^hill_k12)
    growth_mult <- 1 - inh_k12
    k12_ss    <- k12_base * growth_mult
    k12_is    <- k12_base * growth_mult
    k12_ri    <- k12_base * fk12_r * growth_mult

    # 3. Linezolid Inh_Rep effect on successful replication (Eq 7):
    #    depletion of the protein pool P raises Inh_Rep. With
    #    Imax_Rep = 1 (fixed), Inh_Rep = 1 - P. (1 - Inh_Rep) = P
    #    therefore enters the REP-modulated state-2 -> state-1 term
    #    for the two Lin-S populations.
    inh_rep <- imax_rep * (1 - prot_pool)
    rep_mult_lin_s <- 1 - inh_rep
    # Nis-R/Lin-I population has no Inh_Rep effect (paper text after
    # Eq 9: "lacked the effect on the probability of successful
    # replication (Inh_Rep)").
    rep_mult_ri    <- 1.0

    # 4. Total population (CFUall = sum of all 6 bacterial states; Eq 1).
    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 +
              bact_intermediate_susceptible1 + bact_intermediate_susceptible2 +
              bact_resistant_intermediate1 + bact_resistant_intermediate2

    # 5. Replication factor REP (Eq 3, saturating form):
    #    REP = 2 * CFUmax / (CFUmax + CFUall)
    cfumax <- 10^log10cfumax
    REP    <- 2 * cfumax / (cfumax + CFUall)

    # 6. Per-population per-state killing rates from nisin (1/h).
    #    Linezolid does not contribute a direct kill term; killing of
    #    the Lin-S populations is driven by Inh_Rep (replication
    #    failure) and the slowing of the growth rate via Inh_k12.
    kill_ss <- k2s * Cnis
    kill_is <- k2i * Cnis
    kill_ri <- k2r * Cnis

    # 7. Protein pool turnover (Eq 6). Baseline P = 1; linezolid
    #    depletes P by reducing the saturation factor on synthesis.
    d/dt(prot_pool) <- k_prot * ((1 - Clin / (Clin + ic50_prot)) - prot_pool)

    # 8. Two-state life-cycle ODEs per population (Eqs 5, 9 for the
    #    Lin-S populations; the Nis-R/Lin-I population drops the
    #    (1 - Inh_Rep) factor on the replication term per the text
    #    following Eq 9).
    d/dt(bact_susceptible_susceptible1)  <- rep_mult_lin_s * REP * k21 * bact_susceptible_susceptible2 -
                                              (k12_ss + kill_ss) * bact_susceptible_susceptible1
    d/dt(bact_susceptible_susceptible2)  <- k12_ss * bact_susceptible_susceptible1 -
                                              (k21 + kill_ss) * bact_susceptible_susceptible2

    d/dt(bact_intermediate_susceptible1) <- rep_mult_lin_s * REP * k21 * bact_intermediate_susceptible2 -
                                              (k12_is + kill_is) * bact_intermediate_susceptible1
    d/dt(bact_intermediate_susceptible2) <- k12_is * bact_intermediate_susceptible1 -
                                              (k21 + kill_is) * bact_intermediate_susceptible2

    d/dt(bact_resistant_intermediate1)   <- rep_mult_ri * REP * k21 * bact_resistant_intermediate2 -
                                              (k12_ri + kill_ri) * bact_resistant_intermediate1
    d/dt(bact_resistant_intermediate2)   <- k12_ri * bact_resistant_intermediate1 -
                                              (k21 + kill_ri) * bact_resistant_intermediate2

    # 9. Initial conditions (Methods: Initial conditions).
    #    Less-susceptible subpopulations are seeded by mutation
    #    frequencies; the Nis-S/Lin-S population takes the remainder.
    #    All bacteria are placed in state 1 (state 2 = 0) since k21
    #    is rapid. Protein pool starts at its hypothetical baseline
    #    of 1.0 (100%).
    total0  <- 10^log10cfu0
    frac_is <- 10^log10mf_is
    frac_ri <- 10^log10mf_ri
    frac_ss <- 1 - frac_is - frac_ri

    bact_susceptible_susceptible1(0)  <- total0 * frac_ss
    bact_intermediate_susceptible1(0) <- total0 * frac_is
    bact_resistant_intermediate1(0)   <- total0 * frac_ri
    prot_pool(0)                      <- 1.0

    # 10. Observation: log10 total viable count (CFU/mL).
    #     A 1 CFU/mL floor mirrors the experimental limit of counting
    #     (paper Figures 5-6: counts of zero colonies plotted as
    #     0 log10 CFU/mL). Cc carries log10 CFU/mL in this MBM, not a
    #     drug concentration (matches Rees 2018 / Landersdorfer 2018
    #     HFIM-MBM convention).
    cfu_obs <- CFUall + 1
    Cc <- log10(cfu_obs)
    Cc ~ add(addSd)
  })
}
