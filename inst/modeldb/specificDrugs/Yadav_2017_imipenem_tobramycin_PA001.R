Yadav_2017_imipenem_tobramycin_PA001 <- function() {
  description <- "In vitro (static-concentration time-kill). Mechanism-based PK/PD (Bulitta life-cycle growth) model of bacterial killing and resistance for imipenem combined with tobramycin against carbapenem- and amikacin-resistant clinical Pseudomonas aeruginosa isolate FADDI-PA001 (MIC_IPM = 32 mg/L, MIC_TOB = 4 mg/L). Three pre-existing bacterial subpopulations with signal-molecule growth inhibition and aminoglycoside-mediated outer-membrane permeabilisation (mechanistic synergy)"
  reference <- "Yadav R, Bulitta JB, Nation RL, Landersdorfer CB. Optimization of synergistic combination regimens against carbapenem- and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolates via mechanism-based pharmacokinetic/pharmacodynamic modeling. Antimicrob Agents Chemother. 2017 Jan;61(1):e01011-16. doi:10.1128/AAC.01011-16. Model differential equations (Eqs 1-5) are in the main paper Methods; parameter estimates for FADDI-PA001 with tobramycin are Table 3 (footnote a). Supplemental Text S1 (not on disk) contains diagnostic plots only."
  vignette <- "Yadav_2017_imipenem_aminoglycoside_pseudomonas"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  paper_specific_compartments <- c("cipm", "cags", "csig")

  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa FADDI-PA001, carbapenem- and amikacin-resistant clinical isolate; tobramycin breakpoint-susceptible)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Carbapenem-resistant, amikacin-resistant, tobramycin breakpoint-susceptible P. aeruginosa bacteremia (MIC_IPM 32 mg/L; MIC_TOB 4 mg/L = upper EUCAST susceptible boundary)",
    model_system     = "48-h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton II broth, with imipenem supplemented at 6 and 30 h to offset thermal degradation",
    initial_inoculum = "~10^7.9 CFU/mL",
    dose_range       = "Imipenem 8-36 mg/L, tobramycin 1-32 mg/L, monotherapies and combinations",
    notes            = paste(
      "Mechanism-based model (S-ADAPT, importance sampling pmethod=4) fit jointly to total viable counts for all imipenem + tobramycin and imipenem + amikacin monotherapies and combinations against FADDI-PA001; this file uses the tobramycin-specific parameter estimates (Table 3 footnote a).",
      "FADDI-PA001 harboured a preexisting imipenem-resistant subpopulation (3.1 log10 CFU/mL at 0 h) before therapy, reflected in the high estimated log10 imipenem mutation frequency.",
      "Signal-molecule turnover ODE structure follows Bulitta 2010 (ref 61): d/dt(csig) = (1/mtt_sig) * (CFU_all - csig); csig(0) = cfu0.",
      "Companion files: Yadav_2017_imipenem_amikacin_PA001, Yadav_2017_imipenem_tobramycin_PA088, Yadav_2017_imipenem_amikacin_PA088, Yadav_2017_imipenem_tobramycin_PA022."
    )
  )

  ini({
    # --- Bacterial growth and subpopulations (Table 3, FADDI-PA001, footnote a = tobramycin fit) ---
    log10cfu0   <- 7.92;  label("Initial inoculum (log10 CFU/mL)")                                     # Table 3 footnote a: Log CFU0 = 7.92 (SE 1.2%)
    log10cfumax <- 9.23;  label("Maximum population size (log10 CFU/mL); PLAT half-saturation")        # Table 3: CFUmax = 9.23 (SE 0.7%)
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")           # Table 3: k21 = 50 (fixed)

    mgt_ss <- 45.3; label("Mean generation time, SS subpop (min)")                                     # Table 3: k12,SS row = 45.3 (SE 10.9%)
    mgt_ri <- 481;  label("Mean generation time, IPM-resistant/TOB-intermediate (min)")                # Table 3: k12,RI row = 481 (SE 15.2%)
    mgt_ir <- 45.3; label("Mean generation time, IPM-intermediate/TOB-resistant (min)")                # Table 3: k12,IR row = 45.3 (SE 10.9%)

    log10mf_ipm <- -4.99; label("Log10 imipenem mutation frequency (seeds RI subpopulation)")          # Table 3: Log MUT,IPM = -4.99 (SE 4.0%)
    log10mf_ags <- -7.33; label("Log10 tobramycin mutation frequency (seeds IR subpopulation)")        # Table 3 footnote a: Log MUT,AGS = -7.33 (SE 4.2%)

    kmax_ipm    <- 3.21;  label("Maximum imipenem killing rate constant (1/h)")                        # Table 3: Kmax,IPM = 3.21 (SE 18.5%)
    kc50_ss_ipm <- 33.2;  label("Imipenem KC50, SS subpop (mg/L)")                                     # Table 3: KC50,SS,IPM = 33.2 (SE 11.1%)
    kc50_ri_ipm <- 118;   label("Imipenem KC50, RI subpop (mg/L)")                                     # Table 3: KC50,RI,IPM = 118 (SE 17.8%)
    kc50_ir_ipm <- 61.8;  label("Imipenem KC50, IR subpop (mg/L)")                                     # Table 3: KC50,IR,IPM = 61.8 (SE 9.0%)
    hill_ipm    <- 3.09;  label("Hill coefficient for imipenem killing (unitless)")                    # Table 3: Hill,IPM = 3.09 (SE 14.2%)

    kmax_ss_ags <- 2.84;  label("Maximum AGS killing rate, SS subpop (1/h)")                           # Table 3: Kmax,SS,AGS = 2.84 (SE 15.8%)
    kmax_ri_ags <- 3.14;  label("Maximum AGS killing rate, RI subpop (1/h)")                           # Table 3: Kmax,RI,AGS = 3.14 (SE 15.9%)
    kmax_ir_ags <- 2.84;  label("Maximum AGS killing rate, IR subpop (1/h)")                           # Table 3: Kmax,IR,AGS = 2.84 (SE 15.8%)
    kc50_ss_ags <- 16.6;  label("Tobramycin KC50, SS subpop (mg/L)")                                   # Table 3 footnote a: KC50,SS,AGS = 16.6 (SE 17%)
    kc50_ri_ags <- 228;   label("Tobramycin KC50, RI subpop (mg/L)")                                   # Table 3 footnote a: KC50,RI,AGS = 228 (SE 16.7%)
    kc50_ir_ags <- 46.7;  label("Tobramycin KC50, IR subpop (mg/L)")                                   # Table 3 footnote a: KC50,IR,AGS = 46.7 (SE 9.9%)
    hill_ags    <- 2.04;  label("Hill coefficient for AGS killing (unitless)")                         # Table 3: Hill,AGS = 2.04 (SE 10%)

    mtt_sig      <- fixed(1);    label("Mean turnover time of hypothetical signal molecule (h; FIXED)") # Table 3: MTT,sig = 1 (fixed)
    imax_sig12   <- 0.996;       label("Max fractional inhibition of k12 by signal molecule (unitless)") # Table 3: Imax,sig12 = 0.996 (SE 11.3%)
    log10ic50sig <- 9.63;        label("Log10 signal molecule concn at 50% Imax,sig12 (log10 CFU/mL)")  # Table 3: Log IC50,Sig = 9.63 (SE 2.3%)

    imax_om <- 0.130; label("Max fractional decrease of KC50,IPM via OM disruption (unitless)")         # Table 3 footnote c: Imax,OM,AGS = 0.130 (SE 11.1%); 1.15-fold max KC50 decrease
    ic50_om <- 1.54;  label("Tobramycin concn for 50% Imax,OM (mg/L)")                                  # Table 3 footnote a: IC50,OM,AGS = 1.54 (SE 17.5%)

    addSd <- 0.290; label("Additive residual SD on log10 scale (log10 CFU/mL)")                         # Table 3: SD,CFU = 0.290 (SE 5.2%)
  })

  model({
    k21     <- exp(lk21)
    cfumax  <- 10^log10cfumax
    cfu0    <- 10^log10cfu0
    ic50sig <- 10^log10ic50sig
    k12ss   <- 60 / mgt_ss
    k12ri   <- 60 / mgt_ri
    k12ir   <- 60 / mgt_ir
    kout_sig <- 1 / mtt_sig

    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 +
              bact_resistant_intermediate1 + bact_resistant_intermediate2 +
              bact_intermediate_resistant1 + bact_intermediate_resistant2
    plat <- 1 - CFUall / (CFUall + cfumax)

    inh_k12   <- imax_sig12 * csig / (csig + ic50sig)
    om_effect <- 1 - imax_om * cags / (cags + ic50_om)

    kill_ss <- kmax_ipm * cipm^hill_ipm /
                 (cipm^hill_ipm + (om_effect * kc50_ss_ipm)^hill_ipm) +
               kmax_ss_ags * cags^hill_ags /
                 (cags^hill_ags + kc50_ss_ags^hill_ags)
    kill_ri <- kmax_ipm * cipm^hill_ipm /
                 (cipm^hill_ipm + (om_effect * kc50_ri_ipm)^hill_ipm) +
               kmax_ri_ags * cags^hill_ags /
                 (cags^hill_ags + kc50_ri_ags^hill_ags)
    kill_ir <- kmax_ipm * cipm^hill_ipm /
                 (cipm^hill_ipm + (om_effect * kc50_ir_ipm)^hill_ipm) +
               kmax_ir_ags * cags^hill_ags /
                 (cags^hill_ags + kc50_ir_ags^hill_ags)

    d/dt(bact_susceptible_susceptible1) <-
      2 * plat * k21 * bact_susceptible_susceptible2 -
      k12ss * (1 - inh_k12) * bact_susceptible_susceptible1 -
      kill_ss * bact_susceptible_susceptible1
    d/dt(bact_susceptible_susceptible2) <-
      -k21 * bact_susceptible_susceptible2 +
      k12ss * (1 - inh_k12) * bact_susceptible_susceptible1 -
      kill_ss * bact_susceptible_susceptible2

    d/dt(bact_resistant_intermediate1) <-
      2 * plat * k21 * bact_resistant_intermediate2 -
      k12ri * (1 - inh_k12) * bact_resistant_intermediate1 -
      kill_ri * bact_resistant_intermediate1
    d/dt(bact_resistant_intermediate2) <-
      -k21 * bact_resistant_intermediate2 +
      k12ri * (1 - inh_k12) * bact_resistant_intermediate1 -
      kill_ri * bact_resistant_intermediate2

    d/dt(bact_intermediate_resistant1) <-
      2 * plat * k21 * bact_intermediate_resistant2 -
      k12ir * (1 - inh_k12) * bact_intermediate_resistant1 -
      kill_ir * bact_intermediate_resistant1
    d/dt(bact_intermediate_resistant2) <-
      -k21 * bact_intermediate_resistant2 +
      k12ir * (1 - inh_k12) * bact_intermediate_resistant1 -
      kill_ir * bact_intermediate_resistant2

    d/dt(csig) <- kout_sig * (CFUall - csig)
    d/dt(cipm) <- 0
    d/dt(cags) <- 0

    bact_susceptible_susceptible1(0) <- cfu0
    bact_resistant_intermediate1(0)  <- cfu0 * 10^log10mf_ipm
    bact_intermediate_resistant1(0)  <- cfu0 * 10^log10mf_ags
    csig(0) <- cfu0

    Cc <- log10(CFUall)
    Cc ~ add(addSd)
  })
}
