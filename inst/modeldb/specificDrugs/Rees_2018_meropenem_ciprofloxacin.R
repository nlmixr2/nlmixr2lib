Rees_2018_meropenem_ciprofloxacin <- function() {
  description <- "In vitro (hollow-fiber infection model). Mechanism-based PK/PD (life-cycle growth) model of bacterial killing and resistance for meropenem plus ciprofloxacin against hypermutable Pseudomonas aeruginosa CW44, with three pre-existing subpopulations and subpopulation plus mechanistic synergy"
  reference <- "Rees VE, Yadav R, Rogers KE, Bulitta JB, Wirth V, Oliver A, Boyce JD, Peleg AY, Nation RL, Landersdorfer CB. Meropenem combined with ciprofloxacin combats hypermutable Pseudomonas aeruginosa from respiratory infections of cystic fibrosis patients. Antimicrob Agents Chemother. 2018 Oct 24;62(11):e01150-18. doi:10.1128/AAC.01150-18. Model differential equations (Eqs 1-4) and static-time-kill parameters (Table S1) are in the supplemental material."
  vignette <- "Rees_2018_meropenem_ciprofloxacin"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  # No patient covariates: this is an in vitro mechanism-based model. The drug
  # exposures (meropenem and ciprofloxacin concentrations) are state variables
  # (cmem, ccip) dosed by the user, not covariate columns.
  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa CW44, hypermutable cystic fibrosis isolate)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Hypermutable P. aeruginosa respiratory infection in cystic fibrosis",
    model_system     = "Dynamic hollow-fiber infection model (HFIM) over 8 days, simulating cystic-fibrosis epithelial-lining-fluid (ELF) antibiotic exposure",
    initial_inoculum = "~10^7.4 CFU/mL",
    dose_range       = "Meropenem 1-2 g q8h (3-h infusions) or 3 g/day continuous infusion; ciprofloxacin 400 mg q8h (1-h infusions); alone and combined",
    notes            = paste(
      "Mechanism-based model (S-ADAPT) fit to total and resistant viable counts for the double-susceptible hypermutable isolate CW44.",
      "Parameters are the HFIM estimates (Rees 2018 Table 2); the static-concentration time-kill (SCTK) counterpart is Table S1 in the supplement.",
      "Antibiotic ELF concentrations were simulated in Berkeley Madonna with elimination half-lives of 0.8 h (meropenem) and 2.9 h (ciprofloxacin); ELF penetration 30%/60% (meropenem) and 85% (ciprofloxacin).",
      "Strains PAO1, PAOmutS, CW8, and CW35 were studied only in static time-kill assays; the HFIM (and this model) describe CW44."
    )
  )

  ini({
    # --- Bacterial growth and subpopulations (Table 2, HFIM) ---
    log10cfu0   <- 7.37;  label("Initial inoculum (log10 CFU/mL)")                                  # Table 2: Log10CFU0 = 7.37 (SE 2.53%)
    log10cfumax <- 8.80;  label("Maximum population size (log10 CFU/mL)")                           # Table 2: Log10CFUmax = 8.80 (SE 0.796%)
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")        # Table 2: k21 = 50.0 (fixed; fast replication, ref 69)

    # Mean generation time per subpopulation (minutes); growth rate k12 = 60/MGT (1/h)
    mgt_ss <- 181;  label("Mean generation time, double-susceptible MEMs/CIPs (min)")              # Table 2: k12,ss row = 181 (SE 6.32%)
    mgt_ri <- 86.3; label("Mean generation time, MEM-resistant/CIP-intermediate (min)")           # Table 2: k12,ri row = 86.3 (SE 4.40%)
    mgt_ir <- 121;  label("Mean generation time, MEM-intermediate/CIP-resistant (min)")           # Table 2: k12,ir row = 121 (SE 5.82%)

    # Log10 mutation frequencies seeding the two pre-existing resistant subpopulations
    log10mf_mem <- -7.70; label("Log10 meropenem mutation frequency (seeds MEM-resistant subpop)")  # Table 2: Log10MFMEM = -7.70 (SE 1.39%)
    log10mf_cip <- -4.31; label("Log10 ciprofloxacin mutation frequency (seeds CIP-resistant subpop)") # Table 2: Log10MFCIP = -4.31 (SE 3.78%)

    # --- Killing by meropenem (Table 2) ---
    kmax_mem    <- 0.975; label("Maximum meropenem killing rate constant (1/h)")                   # Table 2: Kmax,MEM = 0.975 (SE 12.1%)
    kc50_ss_mem <- 2.08;  label("Meropenem KC50, double-susceptible (mg/L)")                       # Table 2: KC50,ss,MEM = 2.08 (SE 26.6%)
    kc50_ri_mem <- 24.5;  label("Meropenem KC50, MEM-resistant/CIP-intermediate (mg/L)")           # Table 2: KC50,ri,MEM = 24.5 (SE 13.0%)
    kc50_ir_mem <- 5.48;  label("Meropenem KC50, MEM-intermediate/CIP-resistant (mg/L)")           # Table 2: KC50,ir,MEM = 5.48 (SE 8.76%)
    hill_mem    <- 1.60;  label("Hill coefficient for meropenem killing (unitless)")               # Table 2: HILLMEM = 1.60 (SE 11.1%)

    # --- Killing by ciprofloxacin (Table 2) ---
    kmax_cip    <- 5.28;  label("Maximum ciprofloxacin killing rate constant (1/h)")               # Table 2: Kmax,CIP = 5.28 (SE 6.58%)
    kc50_ss_cip <- 1.67;  label("Ciprofloxacin KC50, double-susceptible (mg/L)")                   # Table 2: KC50,ss,CIP = 1.67 (SE 13.9%)
    kc50_ri_cip <- 14.5;  label("Ciprofloxacin KC50, MEM-resistant/CIP-intermediate (mg/L)")       # Table 2: KC50,ri,CIP = 14.5 (SE 11.8%)
    kc50_ir_cip <- 45.8;  label("Ciprofloxacin KC50, MEM-intermediate/CIP-resistant (mg/L)")       # Table 2: KC50,ir,CIP = 45.8 (SE 10.3%)

    # --- Mechanistic synergy: ciprofloxacin lowers the effective meropenem KC50 (Eq 4) ---
    imax_syn <- fixed(1); label("Max fractional decrease of KC50,MEM via mechanistic synergy (unitless)") # Table 2: Imax,SYN = 1 (fixed)
    ic50_syn <- 0.914;    label("Ciprofloxacin concn for 50% of Imax,SYN (mg/L)")                  # Table 2: IC50,SYN = 0.914 (SE 27.8%)

    # --- Residual error (additive on log10 scale) ---
    addSd <- 0.370; label("Additive residual SD on log10 scale (log10 CFU/mL)")                    # Table 2: SDCFU = 0.370 (SE 6.97%)

    # --- Simulated antibiotic disposition in the HFIM (fixed, from clinical PK literature) ---
    thalf_mem <- fixed(0.8); label("Simulated meropenem half-life in HFIM (h)")                    # Methods (refs 48, 49, 56); not an MBM estimate
    thalf_cip <- fixed(2.9); label("Simulated ciprofloxacin half-life in HFIM (h)")                # Methods (refs 48, 49, 56); not an MBM estimate
  })

  model({
    # 0. Back-transform the log replication rate constant.
    k21 <- exp(lk21)

    # 1. Carrying capacity and total inoculum (log10 -> linear, CFU/mL)
    cfumax <- 10^log10cfumax
    cfu0   <- 10^log10cfu0

    # 2. First-order growth rate constants from mean generation time
    #    k12 (1/h) = 60 / MGT(min)  (supplement: "k12 = 60/MGT")
    k12ss <- 60 / mgt_ss
    k12ri <- 60 / mgt_ri
    k12ir <- 60 / mgt_ir

    # 3. Antibiotic elimination rate constants from the simulated HFIM half-lives
    kel_mem <- log(2) / thalf_mem    # 1/h
    kel_cip <- log(2) / thalf_cip    # 1/h

    # 4. Mechanistic synergy factor (Eq 4): ciprofloxacin scales down KC50,MEM.
    #    syn ranges from 1 (no CIP) to 1 - imax_syn (saturating CIP); the
    #    effective meropenem KC50 in each killing term is syn * KC50,MEM, so
    #    1/syn is the fold-decrease in KC50,MEM (~4.3-fold at 3 mg/L CIP).
    syn <- 1 - imax_syn * ccip / (ccip + ic50_syn)

    # 5. Subpopulation killing rates (Eq 2 killing term): Hill-type meropenem
    #    killing (with synergy on KC50,MEM) plus Emax ciprofloxacin killing.
    kill_ss <- kmax_mem * cmem^hill_mem / (cmem^hill_mem + (syn * kc50_ss_mem)^hill_mem) +
      kmax_cip * ccip / (ccip + kc50_ss_cip)
    kill_ri <- kmax_mem * cmem^hill_mem / (cmem^hill_mem + (syn * kc50_ri_mem)^hill_mem) +
      kmax_cip * ccip / (ccip + kc50_ri_cip)
    kill_ir <- kmax_mem * cmem^hill_mem / (cmem^hill_mem + (syn * kc50_ir_mem)^hill_mem) +
      kmax_cip * ccip / (ccip + kc50_ir_cip)

    # 6. Plateau factor (logistic carrying-capacity limit on successful replication)
    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 + bact_resistant_intermediate1 + bact_resistant_intermediate2 + bact_intermediate_resistant1 + bact_intermediate_resistant2
    plat <- 1 - CFUall / cfumax

    # 7. Life-cycle growth model: two states per subpopulation (Eqs 2-3).
    #    State 1 = preparing for replication; state 2 = immediately before
    #    replication. State-2 cells replicate at k21 producing 2 daughter
    #    state-1 cells with success probability plat.
    d/dt(bact_susceptible_susceptible1) <- 2 * plat * k21 * bact_susceptible_susceptible2 - k12ss * bact_susceptible_susceptible1 - kill_ss * bact_susceptible_susceptible1
    d/dt(bact_susceptible_susceptible2) <- -k21 * bact_susceptible_susceptible2 + k12ss * bact_susceptible_susceptible1 - kill_ss * bact_susceptible_susceptible2
    d/dt(bact_resistant_intermediate1) <- 2 * plat * k21 * bact_resistant_intermediate2 - k12ri * bact_resistant_intermediate1 - kill_ri * bact_resistant_intermediate1
    d/dt(bact_resistant_intermediate2) <- -k21 * bact_resistant_intermediate2 + k12ri * bact_resistant_intermediate1 - kill_ri * bact_resistant_intermediate2
    d/dt(bact_intermediate_resistant1) <- 2 * plat * k21 * bact_intermediate_resistant2 - k12ir * bact_intermediate_resistant1 - kill_ir * bact_intermediate_resistant1
    d/dt(bact_intermediate_resistant2) <- -k21 * bact_intermediate_resistant2 + k12ir * bact_intermediate_resistant1 - kill_ir * bact_intermediate_resistant2

    # 8. Antibiotic concentrations (mg/L): dosed by the user, decline per half-life.
    d/dt(cmem) <- -kel_mem * cmem
    d/dt(ccip) <- -kel_cip * ccip

    # 9. Initial conditions: pre-existing subpopulations seeded by mutation
    #    frequency; all inoculum placed in state 1 (state-2 equilibrates in
    #    < ~1% of the population within ~1/k12). cmem/ccip start at 0.
    bact_susceptible_susceptible1(0) <- cfu0
    bact_resistant_intermediate1(0) <- cfu0 * 10^log10mf_mem
    bact_intermediate_resistant1(0) <- cfu0 * 10^log10mf_cip

    # 10. Outputs: total and resistant viable counts; fit on the log10 scale.
    CFUrmem <- bact_resistant_intermediate1 + bact_resistant_intermediate2    # meropenem-resistant subpopulation (CFU/mL)
    CFUrcip <- bact_intermediate_resistant1 + bact_intermediate_resistant2    # ciprofloxacin-resistant subpopulation (CFU/mL)
    Cc <- log10(CFUall)
    Cc ~ add(addSd)
  })
}
