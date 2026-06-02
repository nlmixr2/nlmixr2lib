Yadav_2017_imipenem_tobramycin_PA022 <- function() {
  description <- "In vitro (static-concentration time-kill). Mechanism-based PK/PD (Bulitta life-cycle growth) model of bacterial killing and resistance for imipenem combined with tobramycin against carbapenem-resistant and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolate FADDI-PA022 (MIC_IPM = 16 mg/L, MIC_TOB = 8 mg/L). Three pre-existing bacterial subpopulations with signal-molecule growth inhibition and aminoglycoside-mediated outer-membrane permeabilisation (mechanistic synergy)"
  reference <- "Yadav R, Bulitta JB, Nation RL, Landersdorfer CB. Optimization of synergistic combination regimens against carbapenem- and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolates via mechanism-based pharmacokinetic/pharmacodynamic modeling. Antimicrob Agents Chemother. 2017 Jan;61(1):e01011-16. doi:10.1128/AAC.01011-16. Model differential equations (Eqs 1-5) are in the main paper Methods; parameter estimates for FADDI-PA022 with tobramycin are Table 3 (FADDI-PA022 column). Amikacin combinations against FADDI-PA022 were not studied (Table 1 footnote 'NS'). Supplemental Text S1 (not on disk) contains diagnostic plots only."
  vignette <- "Yadav_2017_imipenem_aminoglycoside_pseudomonas"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  paper_specific_compartments <- c("cipm", "cags", "csig")

  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa FADDI-PA022, carbapenem-resistant, tobramycin-resistant, amikacin-resistant clinical isolate)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Carbapenem-resistant, tobramycin-resistant P. aeruginosa bacteremia (MIC_IPM 16 mg/L; MIC_TOB 8 mg/L; MIC_AMK >32 mg/L per EUCAST)",
    model_system     = "48-h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton II broth, with imipenem supplemented at 6 and 30 h to offset thermal degradation",
    initial_inoculum = "~10^7.0 CFU/mL",
    dose_range       = "Imipenem 8-36 mg/L, tobramycin 4-32 mg/L, monotherapies and combinations",
    notes            = paste(
      "Mechanism-based model fit jointly to total viable counts for imipenem + tobramycin monotherapies and combinations against FADDI-PA022. Amikacin was not studied against FADDI-PA022 (Table 1 'NS') so no amikacin companion file exists.",
      "Note that for FADDI-PA022 the signal-molecule turnover time is estimated (MTT_sig = 0.087 h), not fixed at 1 h as for FADDI-PA088 and FADDI-PA001. The fast turnover gives a much stronger growth-inhibition effect via the signal-molecule mechanism for this isolate.",
      "Signal-molecule turnover ODE structure follows Bulitta 2010 (ref 61): d/dt(csig) = (1/mtt_sig) * (CFU_all - csig); csig(0) = cfu0.",
      "Companion files: Yadav_2017_imipenem_tobramycin_PA088, Yadav_2017_imipenem_amikacin_PA088, Yadav_2017_imipenem_tobramycin_PA001, Yadav_2017_imipenem_amikacin_PA001."
    )
  )

  ini({
    log10cfu0   <- 6.99;  label("Initial inoculum (log10 CFU/mL)")                                     # Table 3: Log CFU0 = 6.99 (SE 1.3%)
    log10cfumax <- 9.54;  label("Maximum population size (log10 CFU/mL); PLAT half-saturation")        # Table 3: CFUmax = 9.54 (SE 2.9%)
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")           # Table 3: k21 = 50 (fixed)

    mgt_ss <- 63.8; label("Mean generation time, SS subpop (min)")                                     # Table 3: k12,SS row = 63.8 (SE 10.8%)
    mgt_ri <- 565;  label("Mean generation time, IPM-resistant/TOB-intermediate (min)")                # Table 3: k12,RI row = 565 (SE 19.9%)
    mgt_ir <- 63.8; label("Mean generation time, IPM-intermediate/TOB-resistant (min)")                # Table 3: k12,IR row = 63.8 (SE 10.8%)

    log10mf_ipm <- -3.90; label("Log10 imipenem mutation frequency (seeds RI subpopulation)")          # Table 3: Log MUT,IPM = -3.90 (SE 9.4%)
    log10mf_ags <- -7.55; label("Log10 tobramycin mutation frequency (seeds IR subpopulation)")        # Table 3: Log MUT,AGS = -7.55 (SE 4.0%)

    kmax_ipm    <- 2.23;  label("Maximum imipenem killing rate constant (1/h)")                        # Table 3: Kmax,IPM = 2.23 (SE 10%)
    kc50_ss_ipm <- 15.2;  label("Imipenem KC50, SS subpop (mg/L)")                                     # Table 3: KC50,SS,IPM = 15.2 (SE 25.7%)
    kc50_ri_ipm <- 67.3;  label("Imipenem KC50, RI subpop (mg/L)")                                     # Table 3: KC50,RI,IPM = 67.3 (SE 24.2%)
    kc50_ir_ipm <- 52.7;  label("Imipenem KC50, IR subpop (mg/L)")                                     # Table 3: KC50,IR,IPM = 52.7 (SE 33.8%)
    hill_ipm    <- 3.08;  label("Hill coefficient for imipenem killing (unitless)")                    # Table 3: Hill,IPM = 3.08 (SE 27.4%)

    kmax_ss_ags <- 2.93;  label("Maximum tobramycin killing rate, SS subpop (1/h)")                    # Table 3: Kmax,SS,AGS = 2.93 (SE 11.8%)
    kmax_ri_ags <- 2.80;  label("Maximum tobramycin killing rate, RI subpop (1/h)")                    # Table 3: Kmax,RI,AGS = 2.80 (SE 21.2%)
    kmax_ir_ags <- 2.93;  label("Maximum tobramycin killing rate, IR subpop (1/h)")                    # Table 3: Kmax,IR,AGS = 2.93 (SE 11.8%)
    kc50_ss_ags <- 54.9;  label("Tobramycin KC50, SS subpop (mg/L)")                                   # Table 3: KC50,SS,AGS = 54.9 (SE 8.1%)
    kc50_ri_ags <- 568;   label("Tobramycin KC50, RI subpop (mg/L)")                                   # Table 3: KC50,RI,AGS = 568 (SE 10.9%)
    kc50_ir_ags <- 159;   label("Tobramycin KC50, IR subpop (mg/L)")                                   # Table 3: KC50,IR,AGS = 159 (SE 34.7%)
    hill_ags    <- 0.998; label("Hill coefficient for tobramycin killing (unitless)")                  # Table 3: Hill,AGS = 0.998 (SE 11.9%)

    # MTT_sig is ESTIMATED (not fixed) for FADDI-PA022; the 0.087 h turnover is much
    # faster than the 1-h fixed value used for FADDI-PA088 and FADDI-PA001.
    mtt_sig      <- 0.087;  label("Mean turnover time of hypothetical signal molecule (h)")            # Table 3: MTT,sig = 0.087 (SE 31.3%)
    imax_sig12   <- 0.881;  label("Max fractional inhibition of k12 by signal molecule (unitless)")    # Table 3: Imax,sig12 = 0.881 (SE 56.2%)
    log10ic50sig <- 8.95;   label("Log10 signal molecule concn at 50% Imax,sig12 (log10 CFU/mL)")      # Table 3: Log IC50,Sig = 8.95 (SE 26%)

    imax_om <- 0.604; label("Max fractional decrease of KC50,IPM via OM disruption (unitless)")         # Table 3 footnote c: Imax,OM,AGS = 0.604 (SE 126%); 2.53-fold max KC50 decrease
    ic50_om <- 4.29;  label("Tobramycin concn for 50% Imax,OM (mg/L)")                                  # Table 3: IC50,OM,AGS = 4.29 (SE 21.6%)

    addSd <- 0.475; label("Additive residual SD on log10 scale (log10 CFU/mL)")                         # Table 3: SD,CFU = 0.475 (SE 10.4%)
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
