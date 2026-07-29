PiresdeMello_2018_zika_FAV_HFIM <- function() {
  description <- paste(
    "In vitro (HUH-7 human hepatoma cells, hollow-fiber infection model).",
    "Refined translational mechanism-based pharmacodynamic (MBM) model of",
    "Zika virus replication and inhibition by favipiravir (FAV) under",
    "dynamic, human-like FAV concentration-time profiles. Twelve-state",
    "model: uninfected host cells (uninfected) with logistic-growth",
    "replication limited by carrying capacity HOSTmax; five sequential",
    "infected host cell stages (infected1..infected5) representing the",
    "delay from infection to virus release; five intracellular virus",
    "transit compartments (vi1..vi5) for viral maturation; and",
    "extracellular virus (vextra) as the observation output (log10",
    "PFU/mL). FAV inhibits viral RNA release between vi4 and vi5 via a",
    "simple Imax/IC50 inhibition function (Eq 8). FAV concentration is a",
    "time-varying covariate (CONC_FAV_UM) driven externally by the",
    "user-supplied clinical PK profile. Parameters are the HFIM column of",
    "Table 1; drug-effect parameters (Imax_FAV, IC50_FAV) and the",
    "additive residual SD are shared estimates with the parallel plate",
    "assay co-fit reported in the same paper.",
    sep = " "
  )
  reference <- paste(
    "Pires de Mello CP, Tao X, Kim TH, Vicchiarelli M, Bulitta JB,",
    "Kaushik A, Brown AN. (2018). Clinical regimens of favipiravir",
    "inhibit Zika virus replication in the hollow-fiber infection model.",
    "Antimicrob Agents Chemother 62(9):e00967-18.",
    "doi:10.1128/AAC.00967-18."
  )
  vignette <- "PiresdeMello_2018_zika_FAV_HFIM"

  paper_specific_compartments <- c(
    "uninfected",
    "infected1", "infected2", "infected3", "infected4", "infected5",
    "vi1", "vi2", "vi3", "vi4", "vi5",
    "vextra"
  )

  units <- list(
    time          = "hour",
    dosing        = "no dosing events -- FAV concentration is supplied via the CONC_FAV_UM time-varying covariate",
    concentration = "log10(PFU/mL) for the model observation Cc"
  )

  covariateData <- list(
    CONC_FAV_UM = list(
      description        = "Favipiravir extracellular (free-drug) concentration in the HFIM medium (uM)",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying covariate driving Eq 8 (INH_FAV). The HFIM system",
        "simulates the free-drug concentration-time profiles associated",
        "with clinically utilized FAV regimens (Pires de Mello 2018 Fig 4),",
        "delivered via syringe pumps to the central reservoir. Two",
        "regimens were prospectively validated: (a) low-dose influenza",
        "(1,800 mg at 0 and 12 h on day 1, then 800 mg q12h from day 2)",
        "and (b) high-dose Ebola (2,400 mg at 0 and 8 h, 1,800 mg at 16 h",
        "on day 1, then 1,200 mg q12h from day 2). Static dose-ranging",
        "studies used continuous infusion at 31.25-500 uM. Set to 0 for",
        "the no-treatment growth control. Paper-specific covariate not in",
        "inst/references/covariate-columns.md because the canonical",
        "concentration concept in nlmixr2lib is a state-derived plasma",
        "concentration (Cc), not an exogenous-drug-concentration covariate",
        "used to drive an in vitro PD model."
      ),
      source_name        = "C_FAV (Pires de Mello 2018 Eq 8)"
    )
  )

  population <- list(
    species          = "in vitro (HUH-7 human hepatoma cell line)",
    n_subjects       = NA_integer_,
    n_studies        = 1L,
    cell_line        = "HUH-7 (human hepatoma)",
    virus_strain     = "Zika virus, 2015 human Puerto Rican strain PRVABC59 (BEI Resources)",
    inoculum         = paste(
      "10^8 HUH-7 cells mixed with 10^5 PFU ZIKV (MOI ~ 0.001 PFU/cell)",
      "inoculated into the extracapillary space of a cellulosic hollow-",
      "fiber cartridge. Initial uninfected host cell density U(0) =",
      "10^6.82 ~ 6.61e6 cells/mL (Log_U fixed, Table 1). No pre-existing",
      "infected cells: I1(0) = 0 (Log_I fixed = 0 for HFIM, Table 1).",
      "Initial extracellular virus V_extra(0) = 6,670 PFU/mL",
      "(Virus_Load fixed, Table 1)."
    ),
    disease_state    = paste(
      "Zika virus in vitro pharmacodynamic experiments in HUH-7 cells",
      "cultured in Dulbecco's modified Eagle's medium (high glucose)",
      "with 5% fetal bovine serum and 1% penicillin-streptomycin in the",
      "hollow-fiber infection model (HFIM) system at 37 C in 5% CO2.",
      "1% DMSO maintained in all medium. ECS sampled daily for 7 days.",
      "Plaque assay limit of detection 100 PFU/mL (log10 = 2). The Beal",
      "M3 method was used in S-ADAPT to handle BLQ samples (e.g. at",
      "time zero for the plate assay)."
    ),
    dose_range       = paste(
      "Dose-ranging continuous-infusion experiments: FAV 31.25-500 uM",
      "(constant). Prospective clinical-regimen validation: low-dose",
      "influenza (1,800/800 mg q12h) and high-dose Ebola",
      "(2,400/1,800/1,200 mg) -- supplied as time-varying free-drug",
      "concentration profiles."
    ),
    notes            = paste(
      "Two replicate experiments were performed for the dose-ranging",
      "studies; the prospective regimen validation used three cartridges",
      "(low-dose, high-dose, growth control). Between-curve variability",
      "was small; the variance was eventually fixed at 0.01 in S-ADAPT",
      "(Table 1 footnote a). The Table 1 RSE percentages reflect this",
      "small between-curve variability and are NOT encoded as etas in",
      "this model file -- the registry uses the typical-value mechanism.",
      "Drug-effect parameters Imax_FAV (0.9998) and IC50_FAV (61.6 uM)",
      "and the additive residual SD (0.286) are shared estimates from",
      "the simultaneous comodel fit of the plate assay and HFIM data;",
      "the plate-assay-specific host-cell parameters (Log_U = 5.82,",
      "Log_I = 3.61, k_infect log10 = -7.36, k_syn = 39.9 1/h,",
      "T_delay = 32.4 h, MST_virus = 15.7 h, no PLAT term) are reported",
      "in Table 1 but not encoded here -- this registry entry models the",
      "HFIM system used for clinical-regimen validation. Model fit and",
      "parameter estimation performed in importance-sampling parallelized",
      "S-ADAPT (version 1.57) with the SADAPT-TRAN facilitator tool;",
      "deterministic simulations in Berkeley Madonna (version 8.23.3.0)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Host-cell dynamics (Pires de Mello 2018 Table 1 HFIM column, Eq 1).
    # log10kinfect is kept on the log10 scale (matching the paper's
    # reported form) and back-transformed in model() via
    # kinfect <- 10^log10kinfect.
    # ---------------------------------------------------------------------
    log10kinfect  <- fixed(-7.03)
    label("log10 of 2nd-order virus-host infection rate constant (mL/PFU/h)")
    # Table 1 row 1 (HFIM): log10(k_infect) = -7.03 (RSE 1.92 percent).

    lksyn         <- fixed(log(39.7))
    label("ksyn -- virus synthesis rate constant in infected cells (1/h)")
    # Table 1 row 2 (HFIM): k_syn = 39.7 1/h (RSE 31.4 percent).

    ltdelay       <- fixed(log(94.0))
    label("T_delay = 5/ktr -- mean delay to virus release / infected-cell survival (h)")
    # Table 1 row 3 (HFIM): T_delay = 5/k_tr = 94.0 h (RSE 18.8 percent);
    # five-transit chain so k_tr = 5/T_delay = 5/94 ~ 0.0532 1/h.

    lmstvirus     <- fixed(log(10.6))
    label("MST_virus = 1/klossvirus -- mean survival time for extracellular virus (h)")
    # Table 1 row 4 (HFIM): MST_virus = 1/k_loss,virus = 10.6 h
    # (RSE 27.9 percent); k_loss = 1/10.6 ~ 0.0943 1/h.

    log10U0       <- fixed(6.82)
    label("log10 initial number of uninfected host cells (cells/mL)")
    # Table 1 row 5 (HFIM): Log_U = 6.82 (fixed); U(0) = 10^6.82
    # ~ 6.61e6 cells/mL.
    # Table 1 row 6 (HFIM): Log_I = 0 (fixed). Footnote c: "no
    # infected cells, i.e. 0 cells on a linear scale." Encoded as
    # infected1(0) = 0 by relying on rxode2's default zero initial
    # condition (no ini() parameter needed; the plate-assay row of
    # Table 1 reports the estimated Log_I = 3.61, applicable only to
    # the plate co-fit and not used in this HFIM extraction).

    log10hostmax  <- fixed(7.39)
    label("log10 plateau (carrying capacity) of total host cell density HOSTmax (cells/mL)")
    # Table 1 row 7 (HFIM): log_max = 7.39 (RSE 0.849 percent);
    # HOSTmax = 10^7.39 ~ 2.45e7 cells/mL. HFIM-only -- the plate
    # assay has no PLAT term (Eq 1: PLAT = 0 in plate assay).

    ltrepl        <- fixed(log(24))
    label("T_repl -- mean doubling time for host cell replication (h)")
    # Table 1 row 8: T_repl = 24 h (fixed); k21 = 1/T_repl = 1/24
    # ~ 0.0417 1/h (Eq 1 cell-replication rate constant).

    lvirusload    <- fixed(log(6670))
    label("Virus_Load -- initial extracellular virus loading V_extra(0) (PFU/mL)")
    # Table 1 row 9 (HFIM): Virus_Load = 6,670 PFU/mL (fixed) -- the
    # HFIM viral inoculum (= 10^3.82 per the prose; Table 1 rounds to
    # 6670). The plate assay row reports 0; the prose of Materials
    # and methods notes ~ 10^2 PFU/mL residual cell-free virus after
    # washing -- a discrepancy noted in the vignette Errata. This
    # HFIM extraction uses the table value.

    # ---------------------------------------------------------------------
    # Drug effect: FAV (Pires de Mello 2018 Table 1, Eq 8). Shared
    # estimates between the plate and HFIM assays in the comodel fit.
    # Imax_FAV is a fraction (0-1); the paper used a logistic transform
    # during estimation (Table 1 footnote b) but reported the normal-
    # scale point estimate.
    # ---------------------------------------------------------------------
    imax_fav      <- fixed(0.9998)
    label("Imax_FAV -- maximum extent of inhibition by FAV (fraction)")
    # Table 1 (shared plate+HFIM): Imax_FAV = 0.9998
    # (95 percent CI on the logistically transformed scale ~0.983-1.00).

    lic50_fav     <- fixed(log(61.6))
    label("IC50_FAV -- FAV concentration causing 50 percent of Imax_FAV (uM)")
    # Table 1 (shared plate+HFIM): IC50_FAV = 61.6 uM (RSE 18.1 percent).

    # ---------------------------------------------------------------------
    # Residual error (Pires de Mello 2018 Table 1 last row; System
    # outputs and residual-error model section: additive error on
    # log10 PFU/mL scale).
    # ---------------------------------------------------------------------
    addSd         <- fixed(0.286)
    label("addSd -- additive residual SD on log10 PFU/mL scale (Beal M3 method for BLQ)")
    # Table 1 (shared plate+HFIM): SDin = 0.286 (RSE 5.09 percent).
  })

  model({
    # ================================================================
    # 1. Back-transform structural parameters from log / log10 scale.
    # ================================================================
    kinfect    <- 10^log10kinfect
    ksyn       <- exp(lksyn)
    tdelay     <- exp(ltdelay)
    ktr        <- 5 / tdelay                # five-transit chain: T_delay = 5/k_tr
    mstvirus   <- exp(lmstvirus)
    klossvirus <- 1 / mstvirus
    trepl      <- exp(ltrepl)
    k21        <- 1 / trepl                 # cell-replication rate (1/h)
    hostmax    <- 10^log10hostmax
    virusload  <- exp(lvirusload)
    ic50_fav   <- exp(lic50_fav)

    # ================================================================
    # 2. FAV inhibition (Eq 8). INH_FAV is the transmission fraction
    #    of virus from vi4 -> vi5: INH_FAV = 1 - Imax_FAV when FAV is
    #    saturating, INH_FAV = 1 when FAV = 0. The +1e-30 floor avoids
    #    a 0/0 limit when CONC_FAV_UM is exactly zero (no drug).
    # ================================================================
    cfav    <- CONC_FAV_UM + 1e-30
    inh_fav <- 1 - imax_fav * cfav / (cfav + ic50_fav)

    # ================================================================
    # 3. PLAT logistic-growth plateau (HFIM only; Eq 1). Total host
    #    cells = U + sum of infected stages. The plateau caps growth
    #    when total cell density approaches HOSTmax.
    # ================================================================
    cells_total <- uninfected + infected1 + infected2 + infected3 +
                   infected4 + infected5
    plat        <- 1 - cells_total / hostmax

    # ================================================================
    # 4. Initial conditions (Pires de Mello 2018 Table 1 HFIM column;
    #    Materials and methods).
    #    uninfected(0)  = 10^Log_U = 10^6.82
    #    infected1..5(0) = 0  (Log_I = 0 -- no infected cells at t=0)
    #    vi1..vi5(0)     = 0
    #    vextra(0)       = Virus_Load = 6,670 PFU/mL
    # ================================================================
    uninfected(0) <- 10^log10U0
    vextra(0)     <- virusload

    # ================================================================
    # 5. ODE system (Pires de Mello 2018 Eqs 1-7, 9).
    #    Note: Eq 6 simplifies algebraically. As written:
    #      dVi4/dt = ktr * (Vi3 - Vi4 * INH_FAV) - ktr * Vi4 * (1 - INH_FAV)
    #             = ktr * Vi3 - ktr * Vi4 * INH_FAV - ktr * Vi4 + ktr * Vi4 * INH_FAV
    #             = ktr * Vi3 - ktr * Vi4
    #    so Vi4 has standard transit dynamics; only Vi5 receives the
    #    fraction INH_FAV of what leaves Vi4 (Eq 7).
    # ================================================================
    # Eq 1: dU/dt = -kinfect * Vextra * U + PLAT * k21 * U (HFIM has PLAT)
    d/dt(uninfected) <- -kinfect * vextra * uninfected + plat * k21 * uninfected

    # Eq 2: infected-cell life-cycle (5 transit stages)
    d/dt(infected1) <- kinfect * vextra * uninfected - ktr * infected1
    d/dt(infected2) <- ktr * (infected1 - infected2)
    d/dt(infected3) <- ktr * (infected2 - infected3)
    d/dt(infected4) <- ktr * (infected3 - infected4)
    d/dt(infected5) <- ktr * (infected4 - infected5)

    # Eqs 3-7: intracellular-virus maturation chain. Vi1 is produced
    # from the total infected-cell population at rate ksyn (consistent
    # with the companion Pires de Mello 2018 Vero MBM, Eq 3, where the
    # symbol "I" in dVi1/dt = ksyn * I - ktr * Vi1 denotes total
    # infected cells -- a single I compartment in the companion paper,
    # extended here to the sum across the 5 infected stages).
    infected_total <- infected1 + infected2 + infected3 + infected4 + infected5
    d/dt(vi1) <- ksyn * infected_total - ktr * vi1
    d/dt(vi2) <- ktr * (vi1 - vi2)
    d/dt(vi3) <- ktr * (vi2 - vi3)
    d/dt(vi4) <- ktr * (vi3 - vi4)
    d/dt(vi5) <- ktr * (vi4 * inh_fav - vi5)

    # Eq 9: extracellular virus
    d/dt(vextra) <- ktr * vi5 - klossvirus * vextra - kinfect * vextra * uninfected

    # ================================================================
    # 6. Observation: log10 PFU/mL of extracellular virus.
    #    Plaque assay LoD is 100 PFU/mL (log10 = 2); the Beal M3
    #    method was used in S-ADAPT to handle BLQ samples. A small
    #    floor (1e-30) avoids -Inf when vextra is exactly 0.
    # ================================================================
    Cc <- log10(vextra + 1e-30)
    Cc ~ add(addSd)
  })
}
