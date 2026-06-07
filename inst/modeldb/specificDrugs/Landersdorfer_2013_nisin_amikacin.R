Landersdorfer_2013_nisin_amikacin <- function() {
  description <- "In vitro (methicillin-resistant Staphylococcus aureus USA300). Mechanism-based pharmacodynamic model for nisin plus amikacin in a 48-h static-concentration time-kill assay (S-ADAPT and NONMEM analyses; subpopulation synergy concept). Six pre-existing bacterial populations crossing nisin (susceptible Nis-S, intermediate Nis-I, resistant Nis-R) with amikacin (susceptible Ami-S, resistant Ami-R) susceptibility, each following a Bulitta two-state life-cycle growth model (state 1 -> state 2 -> 2*state 1 with replication rate k21 fixed). Nisin kills with a second-order function (k2*Cnis) and amikacin with a saturating Hill function; a saturating carrying-capacity replication factor (REP = 2*CFUmax/(CFUmax + CFUall)) caps the population. Nisin- and amikacin-cross-resistant populations carry reduced biofitness (multiplicative growth-rate factor fk12). Nisin and amikacin concentrations are external time-varying inputs (covariates Cnis and Cami); the model contains no human PK component."
  reference <- paste(
    "Landersdorfer CB, Ly NS, Xu H, Tsuji BT, Bulitta JB. (2013).",
    "Quantifying subpopulation synergy for antibiotic combinations via mechanism-based modeling and a sequential dosing design.",
    "Antimicrobial Agents and Chemotherapy 57(5):2343-2351.",
    "doi:10.1128/AAC.00092-13.",
    "Parameter values are the S-ADAPT estimates of Table 1 (nisin + amikacin column);",
    "the NONMEM estimates in the adjacent column are not used here.",
    "Structural equations are Eqs 1-4 of the main text.",
    sep = " "
  )
  vignette <- "Landersdorfer_2013_nisin_amikacin_linezolid"
  units <- list(time = "hour", dosing = "mg/L (drug input concentration)", concentration = "log10 CFU/mL (observation); mg/L (drug covariates)")

  # Cnis / Cami are the experimentally-controlled nisin and amikacin
  # broth concentrations in the static time-kill (and sequential
  # combination) experiments. They are supplied as time-varying
  # covariates from the event data, not estimated PK profiles.
  depends <- c("Cnis", "Cami")

  covariateData <- list(
    Cnis = list(
      description        = "Nisin concentration in the broth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally. Static time-kill: held at a constant assay concentration after addition (4, 8, 32, or 128 mg/L for monotherapy; 8, 16, or 32 mg/L for simultaneous combinations). Sequential combinations: held at 8 or 32 mg/L for 1.5 h pretreatment, then set to 0 at ~1.75 h when bacteria were resuspended in nisin-free broth. In-vitro experimental input -- not in inst/references/covariate-columns.md (the canonical register is for human pop-PK covariates and does not apply to this in-vitro PD model).",
      source_name        = "Nisin concentration (paper Methods, Combinations and Time-kill experiments)"
    ),
    Cami = list(
      description        = "Amikacin concentration in the broth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally. Static time-kill: held at 1, 4, 16, or 64 mg/L for monotherapy. Simultaneous combinations: 4, 8, or 16 mg/L throughout. Sequential combinations: 0 during the 1.5-h nisin pretreatment then 8 or 16 mg/L after resuspension at ~1.75 h. In-vitro experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Amikacin concentration (paper Methods, Combinations and Time-kill experiments)"
    )
  )

  population <- list(
    species             = "in vitro (methicillin-resistant Staphylococcus aureus, USA300 strain from the Network on Antimicrobial Resistance)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "MRSA USA300 (nisin MIC 16 mg/L; amikacin MIC 8 mg/L)",
    system              = "Static-concentration time-kill assay in 20-mL conical tubes; serial viable counts over 48 h",
    medium              = "Luria-Bertani broth supplemented with 12.5 mg/L Mg2+ and 25 mg/L Ca2+; viable counts on Luria-Bertani agar",
    temperature         = "37 C",
    duration            = "48 h",
    inoculum            = "~10^7.7 CFU/mL diluted from a ~10^9 CFU/mL spectrophotometric saline suspension, then grown for 60 min to ~10^8.0 CFU/mL before antibiotic dosing",
    mic_values          = c(nisin = "16 mg/L", amikacin = "8 mg/L"),
    regimens            = "Nisin monotherapy at 4, 8, 32, and 128 mg/L; amikacin monotherapy at 1, 4, 16, and 64 mg/L; nisin pretreatment (8 or 32 mg/L for 1.5 h) followed by removal and resuspension in 0/8/16 mg/L amikacin; simultaneous combinations of nisin (8, 16, 32 mg/L) with amikacin (4, 8, 16 mg/L); 20 viable-count profiles total",
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
    # Table 1 row "Nis_s/Ami_s and Nis_i/Ami_s = 57.3 min" (the two
    # Ami-S subpopulations carrying the susceptible amikacin phenotype
    # share the base growth rate). The Nis-R/Ami-S, Nis-S/Ami-R,
    # Nis-I/Ami-R, and Nis-R/Ami-R subpopulations carry reduced
    # biofitness through the fk12 multipliers below.
    mgt_base <- 57.3
    label("Mean generation time, Nis-S/Ami-S and Nis-I/Ami-S subpopulations (min)")  # Landersdorfer 2013 Table 1

    # Initial total inoculum and carrying capacity (log10 CFU/mL).
    log10cfu0   <- 7.89
    label("Initial total inoculum (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (Log10 initial inoculum)
    log10cfumax <- 9.23
    label("Maximum population size (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (Log10 CFUmax)

    # =============================================================
    # Initial subpopulation sizes (log10 mutation frequencies)
    # =============================================================
    # Per Methods (Initial conditions): the product of total inoculum
    # and mutation frequency seeds each less-susceptible population at
    # t = 0; the Nis-S/Ami-S population is the remainder.
    log10mf_is <- -4.25
    label("log10 mutation frequency for Nis-I/Ami-S subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1
    log10mf_rs <- -3.25
    label("log10 mutation frequency for Nis-R/Ami-S subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1
    log10mf_sr <- -2.57
    label("log10 mutation frequency for Nis-S/Ami-R subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1
    log10mf_ir <- -4.37
    label("log10 mutation frequency for Nis-I/Ami-R subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1
    log10mf_rr <- -7.42
    label("log10 mutation frequency for Nis-R/Ami-R subpopulation (log10 fraction; unitless)")  # Landersdorfer 2013 Table 1

    # =============================================================
    # Biofitness reduction factors (multiply k12_base)
    # =============================================================
    # Footnote a: replication of the Nis-R/Ami-S and Nis-S/Ami-R
    # populations was estimated to be very slow and eventually fixed
    # to 0. fk12_IR and fk12_RR were estimated via a logistic
    # transformation in S-ADAPT (footnote d).
    fk12_rs_sr <- fixed(0)
    label("Growth-rate factor for Nis-R/Ami-S and Nis-S/Ami-R subpopulations (unitless; FIXED)")  # Landersdorfer 2013 Table 1 footnote a
    fk12_ir    <- 0.532
    label("Growth-rate factor for Nis-I/Ami-R subpopulation (unitless)")  # Landersdorfer 2013 Table 1 (S-ADAPT)
    fk12_rr    <- 0.413
    label("Growth-rate factor for Nis-R/Ami-R subpopulation (unitless)")  # Landersdorfer 2013 Table 1 (S-ADAPT)

    # =============================================================
    # Nisin killing -- second-order in nisin concentration
    # =============================================================
    lk2s <- log(5.67)
    label("Second-order nisin killing rate constant, Nis-S populations (L/(mg*h))")  # Landersdorfer 2013 Table 1
    lk2i <- log(0.0664)
    label("Second-order nisin killing rate constant, Nis-I populations (L/(mg*h))")  # Landersdorfer 2013 Table 1
    lk2r <- log(0.00691)
    label("Second-order nisin killing rate constant, Nis-R populations (L/(mg*h))")  # Landersdorfer 2013 Table 1

    # =============================================================
    # Amikacin killing -- saturating Hill function (shared KC50, Hill)
    # =============================================================
    lkmax_s <- log(10.1)
    label("Maximum amikacin killing rate constant, Ami-S populations (1/h)")  # Landersdorfer 2013 Table 1 (Kmax,S)
    lkmax_r <- log(0.771)
    label("Maximum amikacin killing rate constant, Ami-R populations (1/h)")  # Landersdorfer 2013 Table 1 (Kmax,R)
    lkc50   <- log(14.7)
    label("Amikacin concentration causing 50% of Kmax (mg/L)")  # Landersdorfer 2013 Table 1 (KC50)
    lhill   <- log(2.45)
    label("Amikacin Hill coefficient (unitless)")  # Landersdorfer 2013 Table 1 (Hill)

    # =============================================================
    # Residual error (additive on log10 viable count)
    # =============================================================
    # Per Methods (Observation model): all log10 viable counts fitted
    # with an additive residual error on the log10 scale. For counts
    # below 100 CFU/mL the paper used a per-colony-count residual
    # model (ref 19); only the additive log10 component is packaged
    # here, matching the published Table 1 value.
    addSd <- 0.395
    label("Additive residual SD on log10 total viable count (log10 CFU/mL)")  # Landersdorfer 2013 Table 1 (SD_CFU)
  })

  model({
    # 0. Back-transform log-scale parameters.
    k21    <- exp(lk21)
    k2s    <- exp(lk2s)
    k2i    <- exp(lk2i)
    k2r    <- exp(lk2r)
    kmax_s <- exp(lkmax_s)
    kmax_r <- exp(lkmax_r)
    kc50   <- exp(lkc50)
    hill   <- exp(lhill)

    # 1. Base growth rate k12 (1/h) from mean generation time (min).
    k12_base <- 60 / mgt_base

    # 2. Per-subpopulation k12 with biofitness multipliers (Table 1).
    #    The two Ami-S / Nis-{S, I} populations share k12_base; the
    #    Nis-R/Ami-S and Nis-S/Ami-R populations carry fk12_rs_sr = 0
    #    (no replication); the Nis-I/Ami-R and Nis-R/Ami-R populations
    #    use fk12_ir and fk12_rr respectively.
    k12_ss <- k12_base
    k12_is <- k12_base
    k12_rs <- k12_base * fk12_rs_sr
    k12_sr <- k12_base * fk12_rs_sr
    k12_ir <- k12_base * fk12_ir
    k12_rr <- k12_base * fk12_rr

    # 3. Total population (CFUall = sum of all 12 bacterial states; Eq 1).
    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 +
              bact_intermediate_susceptible1 + bact_intermediate_susceptible2 +
              bact_resistant_susceptible1 + bact_resistant_susceptible2 +
              bact_susceptible_resistant1 + bact_susceptible_resistant2 +
              bact_intermediate_resistant1 + bact_intermediate_resistant2 +
              bact_resistant_resistant1 + bact_resistant_resistant2

    # 4. Replication factor REP (Eq 3, saturating form):
    #    REP = 2 * (1 - CFUall/(CFUmax + CFUall))
    #        = 2 * CFUmax / (CFUmax + CFUall)
    #    REP -> 2 as CFUall -> 0 (doubling); REP -> 1 as
    #    CFUall -> CFUmax (50% successful replication; net stasis).
    cfumax <- 10^log10cfumax
    REP    <- 2 * cfumax / (cfumax + CFUall)

    # 5. Amikacin saturating killing factor (shared across populations).
    cami_h <- Cami^hill
    Fami   <- cami_h / (cami_h + kc50^hill)

    # 6. Per-population per-state killing rates (1/h).
    #    Eq 2: nisin contributes a second-order term k2*Cnis and
    #    amikacin contributes the Kmax * Fami Hill term; the same
    #    killing rate applies to state 1 and state 2 of each
    #    population (Eq 4 for state 2 includes the same killing term
    #    as Eq 2 for state 1).
    kill_ss <- k2s * Cnis + kmax_s * Fami
    kill_is <- k2i * Cnis + kmax_s * Fami
    kill_rs <- k2r * Cnis + kmax_s * Fami
    kill_sr <- k2s * Cnis + kmax_r * Fami
    kill_ir <- k2i * Cnis + kmax_r * Fami
    kill_rr <- k2r * Cnis + kmax_r * Fami

    # 7. Two-state life-cycle ODEs per population (Eqs 2 and 4 for SS;
    #    the analogous equations for the other 5 populations follow
    #    the same structure with the appropriate k12, k2, and Kmax).
    d/dt(bact_susceptible_susceptible1)  <- REP * k21 * bact_susceptible_susceptible2  - (k12_ss + kill_ss) * bact_susceptible_susceptible1
    d/dt(bact_susceptible_susceptible2)  <- k12_ss   * bact_susceptible_susceptible1   - (k21    + kill_ss) * bact_susceptible_susceptible2

    d/dt(bact_intermediate_susceptible1) <- REP * k21 * bact_intermediate_susceptible2 - (k12_is + kill_is) * bact_intermediate_susceptible1
    d/dt(bact_intermediate_susceptible2) <- k12_is   * bact_intermediate_susceptible1  - (k21    + kill_is) * bact_intermediate_susceptible2

    d/dt(bact_resistant_susceptible1)    <- REP * k21 * bact_resistant_susceptible2    - (k12_rs + kill_rs) * bact_resistant_susceptible1
    d/dt(bact_resistant_susceptible2)    <- k12_rs   * bact_resistant_susceptible1     - (k21    + kill_rs) * bact_resistant_susceptible2

    d/dt(bact_susceptible_resistant1)    <- REP * k21 * bact_susceptible_resistant2    - (k12_sr + kill_sr) * bact_susceptible_resistant1
    d/dt(bact_susceptible_resistant2)    <- k12_sr   * bact_susceptible_resistant1     - (k21    + kill_sr) * bact_susceptible_resistant2

    d/dt(bact_intermediate_resistant1)   <- REP * k21 * bact_intermediate_resistant2   - (k12_ir + kill_ir) * bact_intermediate_resistant1
    d/dt(bact_intermediate_resistant2)   <- k12_ir   * bact_intermediate_resistant1    - (k21    + kill_ir) * bact_intermediate_resistant2

    d/dt(bact_resistant_resistant1)      <- REP * k21 * bact_resistant_resistant2      - (k12_rr + kill_rr) * bact_resistant_resistant1
    d/dt(bact_resistant_resistant2)      <- k12_rr   * bact_resistant_resistant1       - (k21    + kill_rr) * bact_resistant_resistant2

    # 8. Initial conditions (Methods: Initial conditions).
    #    Less-susceptible subpopulations are seeded by mutation
    #    frequencies; the Nis-S/Ami-S population takes the remainder.
    #    All bacteria are placed in state 1 (state 2 = 0) since k21
    #    is rapid.
    total0  <- 10^log10cfu0
    frac_is <- 10^log10mf_is
    frac_rs <- 10^log10mf_rs
    frac_sr <- 10^log10mf_sr
    frac_ir <- 10^log10mf_ir
    frac_rr <- 10^log10mf_rr
    frac_ss <- 1 - frac_is - frac_rs - frac_sr - frac_ir - frac_rr

    bact_susceptible_susceptible1(0)  <- total0 * frac_ss
    bact_intermediate_susceptible1(0) <- total0 * frac_is
    bact_resistant_susceptible1(0)    <- total0 * frac_rs
    bact_susceptible_resistant1(0)    <- total0 * frac_sr
    bact_intermediate_resistant1(0)   <- total0 * frac_ir
    bact_resistant_resistant1(0)      <- total0 * frac_rr

    # 9. Observation: log10 total viable count (CFU/mL).
    #    A 1 CFU/mL floor mirrors the experimental limit of counting
    #    (paper Figures 5-6: counts of zero colonies plotted as
    #    0 log10 CFU/mL). Cc carries log10 CFU/mL in this MBM, not a
    #    drug concentration (matches Rees 2018 / Landersdorfer 2018
    #    HFIM-MBM convention).
    cfu_obs <- CFUall + 1
    Cc <- log10(cfu_obs)
    Cc ~ add(addSd)
  })
}
