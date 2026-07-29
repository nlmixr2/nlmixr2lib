Bulitta_2009_ceftazidime <- function() {
  description <- "In vitro (Pseudomonas aeruginosa PAO1). Mechanism-based PD model (model D) for the pronounced inoculum effect of ceftazidime, fit to static time-kill experiments at initial inocula 10^6 / 10^7 / 10^8 CFU/mL with ceftazidime 0-128 mg/L. Two co-existing bacterial subpopulations (genotypically susceptible and resistant) each follow a two-state life-cycle (S1 -> S2 -> 2 S1) with a shared slow S1 -> S2 generation step (k12) and a shared fast doubling step (k21 fixed at 50 /h). A logistic-saturation replication-efficiency factor Rep gates the doubling step toward the carrying capacity CFUmax. Ceftazidime stimulates an autolysin effect (alys_s / alys_r) that decreases the per-doubling success probability via (1 - alys); the resistant population has a smaller maximum autolysin effect (Smax_r < Smax_s). Freely diffusible signal molecules (csig1 / csig2) synthesised by all viable bacteria mediate quorum-sensing-like phenotypic tolerance at high inocula: high csig1 suppresses the autolysin loss (the inoculum effect) and prolongs the generation time (drug-conditional Inhk12 factor). Ceftazidime broth concentration cb degrades first-order at fixed half-life 45.9 h (Viaene 1973). The model has no human PK component; the experimental drug exposure is supplied as a dose into cb at time zero, and the bacterial / signal initial conditions are derived from the inoculum log10_cfuo parameter. Random effects (eta) are NOT estimated in the NONMEM analysis -- the paper reports the additive log10-scale residual SD only (between-run variability was negligible). The packaged model is therefore intended for typical-value simulation."
  reference <- paste(
    "Bulitta JB, Ly NS, Yang JC, Forrest A, Jusko WJ, Tsuji BT (2009).",
    "Development and qualification of a pharmacodynamic model for the pronounced inoculum effect",
    "of ceftazidime against Pseudomonas aeruginosa.",
    "Antimicrobial Agents and Chemotherapy 53(1):46-56.",
    "doi:10.1128/AAC.00489-08.",
    sep = " "
  )
  vignette <- "Bulitta_2009_ceftazidime"
  units <- list(
    time = "hour",
    dosing = "mg/L (initial ceftazidime concentration in broth)",
    concentration = "log10 CFU/mL (observation); mg/L (drug state cb)"
  )

  paper_specific_compartments <- c(
    "bact_s1", "bact_s2",
    "bact_r1", "bact_r2",
    "alys_s", "alys_r",
    "csig1", "csig2",
    "cb"
  )

  covariateData <- list()

  population <- list(
    species             = "in vitro (Pseudomonas aeruginosa, PAO1 clinical isolate from the R.E.W. Hancock laboratory)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Pseudomonas aeruginosa PAO1 (ceftazidime MIC 2 mg/L); the same model was also fit to P. aeruginosa ATCC 27853 in an external-qualification arm with the ATCC parameter set reported in Table 2 (not packaged here; only the primary PAO1 NONMEM fit is reproduced)",
    system              = "Static time-kill experiments at five inocula (10^5 / 10^6 / 10^7 / 10^8 / 10^9 CFU/mL); samples at 0, 0.5, 1, 2, 4, 8, 24 h (10^7 also at 37 and 48 h); duplicate per concentration",
    medium              = "Luria-Bertani broth supplemented with calcium 25 mg/L and magnesium 12.5 mg/L",
    temperature         = "37 C",
    duration            = "24 h (48 h for 10^7 CFU/mL)",
    mic_values          = c(ceftazidime = "2 mg/L"),
    concentration_range = c(ceftazidime = "0-128 mg/L (0-64 x MIC)"),
    regimens            = "Ceftazidime monotherapy at static concentrations 0, 2, 4, 8, 16, 32, 64, 128 mg/L applied at t=0; broth concentration degrades first-order at fixed half-life 45.9 h (Viaene 1973)",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. Final NONMEM fit reports the additive residual SD on log10(CFU/mL) only (sigma = 0.224); IIV / between-replicate etas were not estimated because the duplicate-run variability was small. An S-ADAPT confirmation run with BSV on every parameter is reported in Table 2 but those BSV variances are not packaged. Parameters here reproduce the PAO1 NONMEM column of Table 2."
  )

  ini({
    # ---- Initial inoculum and carrying capacity ------------------------------
    # Default log10 CFUo set to the middle of the studied range (10^7.37); the
    # paper's Table 2 lists fit values 6.01 / 7.37 / 8.10 for the three studied
    # inocula. Override this parameter to simulate other CFUo.
    log10_cfuo <- 7.37
    label("Initial total inoculum CFUo (log10 CFU/mL)")  # Bulitta 2009 Table 2, "Log 10 IC for 10^7 CFU/ml"
    log10_cfumax <- 9.78
    label("Maximum population size CFUmax (log10 CFU/mL)")  # Bulitta 2009 Table 2, "Log 10 CFU max"
    log10_frr <- -3.54
    label("log10 fraction of resistant subpopulation at t = 0 (Log10 FrR)")  # Bulitta 2009 Table 2, "Log 10 fraction of resistant population at time zero"

    # ---- Bacterial life cycle ------------------------------------------------
    # k12 = slow generation step (1/MTT12); k21 = fast doubling step (fixed at
    # 50/h per Table 2 footnote: assumed not rate-limiting). MTT12 is parameterised
    # in minutes per Table 2 column "Unit"; conversion to /h happens in model().
    lmtt12 <- log(28.3)
    label("Mean generation time MTT12 (min)")  # Bulitta 2009 Table 2, "Generation time at low CFUo"
    lk21 <- fixed(log(50))
    label("Fast doubling rate constant k21 (1/h; FIXED -- assumed not rate-limiting)")  # Bulitta 2009 Table 2 footnote a

    # ---- Signal-molecule kinetics -------------------------------------------
    # csig1 = central signal pool; csig2 = peripheral signal pool. MTT_S10 is
    # the production-and-loss MTT in csig1 (synthesis term tracks SRALL with
    # rate kS10 = 60/MTT_S10 [/h]). MTT_S12 governs central -> peripheral, fixed
    # at 1 min per Table 2 footnote b (signal-molecule onset is very fast based
    # on starvation-experiment evidence). MTT_S21 governs peripheral -> central,
    # fixed at 24 h per Table 2 footnote c (low identifiability, fixed at a value
    # >= 48 h estimate).
    lmtt_s10 <- log(2.33)
    label("Signal-molecule synthesis/loss MTT in csig1 (min)")  # Bulitta 2009 Table 2, "MTT for elimination and synthesis"
    lmtt_s12 <- fixed(log(1))
    label("Signal-molecule MTT csig1 -> csig2 (min; FIXED)")  # Bulitta 2009 Table 2 footnote b
    lmtt_s21 <- fixed(log(24))
    label("Signal-molecule MTT csig2 -> csig1 (h; FIXED)")  # Bulitta 2009 Table 2 footnote c

    # ---- Autolysin / drug effect --------------------------------------------
    # Smax,S = 1 fixed (full inhibition of successful replication achievable
    # in susceptible population at high drug); Smax,R = 0.560 estimated (the
    # resistant population is partially protected). SC50 shared across subpops.
    # The kout governing alys turnover is parameterised as a ratio to k12.
    smax_s <- fixed(1)
    label("Maximum autolysin effect Smax,S in susceptible subpop (FIXED at 1)")  # Bulitta 2009 Table 2 footnote e
    smax_r <- 0.560
    label("Maximum autolysin effect Smax,R in resistant subpop (unitless)")  # Bulitta 2009 Table 2
    lsc50 <- log(0.294)
    label("Drug concentration at 50% stimulation of autolysin SC50 (mg/L)")  # Bulitta 2009 Table 2
    ratio_kout_k12 <- 0.438
    label("Ratio kout/k12 of autolysin loss rate to slow growth-step rate (unitless)")  # Bulitta 2009 Table 2

    # ---- Inoculum-effect (signal-mediated) parameters ------------------------
    log10_c50sig <- 7.24
    label("Signal-molecule concentration at 50% of maximum inoculum effect (Log10 C50,Sig in CFU/mL)")  # Bulitta 2009 Table 2
    smax_loss <- 1.18
    label("Maximum stimulation of autolysin loss by signal molecules Smax,loss (unitless)")  # Bulitta 2009 Table 2

    # ---- Drug-dependent prolongation of generation time ---------------------
    # EC50,drug and Smax,k12 govern the Inhk12 multiplier that couples high
    # drug + high signal to slower replication. Smax,k12 fixed at 10 per Table 2
    # footnote d (low identifiability in this data set).
    lec50_drug <- log(35.3)
    label("Drug concentration at 50% stimulation of generation-time prolongation EC50,drug (mg/L)")  # Bulitta 2009 Table 2
    smax_k12 <- fixed(10)
    label("Maximum stimulation of generation-time prolongation Smax,k12 (FIXED)")  # Bulitta 2009 Table 2 footnote d

    # ---- Ceftazidime broth degradation --------------------------------------
    # First-order degradation in growth medium; half-life fixed at 45.9 h from
    # Viaene et al. 1973 (Methods, "the ceftazidime concentration in broth (CB)
    # was assumed to degrade with a fixed half-life of 45.9 h").
    lkdeg <- fixed(log(log(2) / 45.9))
    label("Ceftazidime first-order degradation rate kdeg (1/h; FIXED, t1/2 = 45.9 h)")  # Bulitta 2009 Methods (cites Viaene 1973)

    # ---- Residual error -----------------------------------------------------
    # Paper: "additive error model on the log scale was applied" (Methods,
    # log-both-sides). Observation is log10(SRALL); additive residual SD on
    # log10 CFU/mL per Table 2 (sigma = 0.224).
    addSd <- 0.224
    label("Additive residual SD on log10 total viable count (log10 CFU/mL)")  # Bulitta 2009 Table 2, "SD of residual error on log 10 scale"
  })

  model({
    # ---- Effective parameters (delog and rate-constant conversions) ---------
    mtt12   <- exp(lmtt12)
    k12     <- 60 / mtt12              # min -> /h
    k21     <- exp(lk21)               # /h
    mtt_s10 <- exp(lmtt_s10)
    kS10    <- 60 / mtt_s10            # min -> /h
    mtt_s12 <- exp(lmtt_s12)
    kS12    <- 60 / mtt_s12            # min -> /h
    mtt_s21 <- exp(lmtt_s21)
    kS21    <- 1  / mtt_s21            # h -> /h
    sc50    <- exp(lsc50)
    kout    <- ratio_kout_k12 * k12    # /h (Bulitta 2009 parameterises kout via its ratio to k12)
    c50sig  <- 10 ^ log10_c50sig
    ec50dr  <- exp(lec50_drug)
    kdeg    <- exp(lkdeg)
    cfumax  <- 10 ^ log10_cfumax
    cfuo    <- 10 ^ log10_cfuo
    frr     <- 10 ^ log10_frr

    # ---- Total viable count and replication-efficiency factor (Eqs 7-8) -----
    # Rep approaches 2 at zero density (perfect doubling) and decreases as
    # SRALL approaches CFUmax, encoding logistic-saturation toward carrying
    # capacity. The doubling step in the bacterial ODEs uses Rep * k21.
    sr_all <- bact_s1 + bact_s2 + bact_r1 + bact_r2
    rep_eff <- 2 * (1 - sr_all / (cfumax + sr_all))

    # ---- Signal-molecule effect on generation time (Eqs 13-14) --------------
    # csig_k12 is the effective signal-molecule level after drug-induced
    # prolongation; Inhk12 in [0, 1] scales the slow S1 -> S2 transition.
    csig_k12 <- csig1 * (1 + smax_k12 * cb / (ec50dr + cb))
    inhk12   <- 1 - csig_k12 / (c50sig + csig_k12)

    # ---- Drug stimulation of autolysin (Eq 10) ------------------------------
    # StimDrug shared form; per Table 2 the susceptible population has Smax_s
    # and the resistant population Smax_r, with SC50 shared.
    stim_s <- smax_s * cb / (sc50 + cb)
    stim_r <- smax_r * cb / (sc50 + cb)

    # ---- Signal-mediated loss of autolysin effect (Eq 9 inner term) ---------
    # The loss-amplifier (1 + Smax,loss * csig1 / (C50,Sig + csig1)) is shared
    # between subpops; it depends on the central signal pool only.
    loss_amp <- 1 + smax_loss * csig1 / (c50sig + csig1)

    # ---- Bacterial life-cycle ODEs ------------------------------------------
    # Susceptible (Eqs 5-6): drug-stimulated autolysin alys_s reduces the
    # successful-doubling probability via (1 - alys_s); both the slow S1 -> S2
    # transition and the fast doubling step are common to both subpops, with
    # the slow transition further modulated by Inhk12.
    d/dt(bact_s1) <-  rep_eff * (1 - alys_s) * k21 * bact_s2 - inhk12 * k12 * bact_s1
    d/dt(bact_s2) <- -k21 * bact_s2 + inhk12 * k12 * bact_s1

    # Resistant (analogous to susceptible per paper text "the differential
    # equations for the genotypically resistant population were the same as
    # that for the susceptible population"; only Smax differs -- captured
    # downstream via alys_r).
    d/dt(bact_r1) <-  rep_eff * (1 - alys_r) * k21 * bact_r2 - inhk12 * k12 * bact_r1
    d/dt(bact_r2) <- -k21 * bact_r2 + inhk12 * k12 * bact_r1

    # ---- Autolysin-activity ODEs (Eq 9) -------------------------------------
    d/dt(alys_s) <- (stim_s - loss_amp * alys_s) * kout
    d/dt(alys_r) <- (stim_r - loss_amp * alys_r) * kout

    # ---- Signal-molecule ODEs (Eqs 11-12) -----------------------------------
    # Central pool csig1 is synthesised by all viable bacteria with a tracking
    # rate kS10 and exchanges with peripheral pool csig2 at rates kS12 / kS21.
    d/dt(csig1) <- (sr_all - csig1) * kS10 - kS12 * csig1 + kS21 * csig2
    d/dt(csig2) <-  kS12 * csig1 - kS21 * csig2

    # ---- Ceftazidime broth concentration ------------------------------------
    # First-order chemical degradation only (no PK -- this is in-vitro).
    d/dt(cb) <- -kdeg * cb

    # ---- Initial conditions (Eqs 5-6, 11-12) --------------------------------
    # Per Table 2 the susceptible pool starts at CFUo*(1 - FrR) entirely in
    # state S1; the resistant pool starts at CFUo*FrR entirely in state R1.
    # The signal-molecule peripheral pool starts at its quasi-steady-state
    # ratio CFUo*kS12/kS21 (paper Eq 12 IC) so that the signal system begins
    # at balance. Autolysin pools start at zero (no drug pre-exposure).
    bact_s1(0) <- cfuo * (1 - frr)
    bact_s2(0) <- 0
    bact_r1(0) <- cfuo * frr
    bact_r2(0) <- 0
    alys_s(0)  <- 0
    alys_r(0)  <- 0
    csig1(0)   <- cfuo
    csig2(0)   <- cfuo * kS12 / kS21
    # cb is set by the dosing event at time 0 (user supplies the broth
    # ceftazidime concentration in mg/L as the dose amount into cb).

    # ---- Observation: log10 total viable count (CFU/mL) ---------------------
    # The single observation Cc carries log10(SRALL) per the nlmixr2lib
    # single-output convention (see units$concentration). A 1e-6 floor avoids
    # log10(0) when drug kill drives all states toward zero. Matches the
    # Wicha 2017 / Landersdorfer 2018 pharmacodynamics/ pattern. Additive
    # residual SD on log10 CFU/mL per Bulitta 2009 Table 2.
    Cc <- log10(sr_all + 1e-6)
    Cc ~ add(addSd)
  })
}
