Yadav_2017_imipenem_amikacin_PA088 <- function() {
  description <- "In vitro (static-concentration time-kill). Mechanism-based PK/PD (Bulitta life-cycle growth) model of bacterial killing and resistance for imipenem combined with amikacin against carbapenem- and tobramycin-resistant clinical Pseudomonas aeruginosa isolate FADDI-PA088 (MIC_IPM = 16 mg/L, MIC_AMK = 4 mg/L). Three pre-existing bacterial subpopulations with signal-molecule growth inhibition and aminoglycoside-mediated outer-membrane permeabilisation (mechanistic synergy)"
  reference <- "Yadav R, Bulitta JB, Nation RL, Landersdorfer CB. Optimization of synergistic combination regimens against carbapenem- and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolates via mechanism-based pharmacokinetic/pharmacodynamic modeling. Antimicrob Agents Chemother. 2017 Jan;61(1):e01011-16. doi:10.1128/AAC.01011-16. Model differential equations (Eqs 1-5) are in the main paper Methods; parameter estimates for FADDI-PA088 with amikacin are Table 3 (footnote b). Supplemental Text S1 (not on disk) contains diagnostic plots only."
  vignette <- "Yadav_2017_imipenem_aminoglycoside_pseudomonas"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  paper_specific_compartments <- c("cipm", "cags", "csig")

  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa FADDI-PA088, carbapenem- and tobramycin-resistant clinical isolate; amikacin-susceptible)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Carbapenem-resistant, tobramycin-resistant, amikacin-susceptible P. aeruginosa bacteremia (MIC_IPM 16 mg/L; MIC_AMK 4 mg/L)",
    model_system     = "48-h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton II broth, with imipenem supplemented at 6 and 30 h to offset thermal degradation",
    initial_inoculum = "~10^7.5 CFU/mL",
    dose_range       = "Imipenem 8-36 mg/L, amikacin 0.5-64 mg/L, monotherapies and combinations",
    notes            = paste(
      "Mechanism-based model (S-ADAPT, importance sampling pmethod=4) fit jointly to total viable counts for all imipenem + tobramycin and imipenem + amikacin monotherapies and combinations against FADDI-PA088; this file uses the amikacin-specific parameter estimates (Table 3 footnote b).",
      "Three pre-existing bacterial subpopulations: SS (IPM-susceptible, AGS-susceptible), RI (IPM-resistant, AGS-intermediate), IR (IPM-intermediate, AGS-resistant), each described by a two-state Bulitta life-cycle growth model.",
      "Outer-membrane synergy: amikacin permeabilises the outer membrane and reduces the effective imipenem KC50 to OM_effect * KC50_IPM. For FADDI-PA088, KC50_IPM decreases up to 3.11-fold at high amikacin (Table 3 footnote c).",
      "Signal-molecule turnover ODE structure follows Bulitta 2010 (ref 61): d/dt(csig) = (1/mtt_sig) * (CFU_all - csig); csig(0) = cfu0.",
      "Companion files: Yadav_2017_imipenem_tobramycin_PA088, Yadav_2017_imipenem_tobramycin_PA001, Yadav_2017_imipenem_amikacin_PA001, Yadav_2017_imipenem_tobramycin_PA022."
    )
  )

  ini({
    # --- Bacterial growth and subpopulations (Table 3, FADDI-PA088, footnote b = amikacin fit) ---
    log10cfu0   <- 7.46;  label("Initial inoculum (log10 CFU/mL)")                                     # Table 3 footnote b: Log CFU0 = 7.46 (SE 1.8%)
    log10cfumax <- 9.56;  label("Maximum population size (log10 CFU/mL); PLAT half-saturation")        # Table 3: CFUmax = 9.56 (SE 1.8%) [strain-shared]
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")           # Table 3: k21 = 50 (fixed)

    mgt_ss <- 26.9; label("Mean generation time, SS subpop (min)")                                     # Table 3: k12,SS row = 26.9 (SE 6.2%) [strain-shared]
    mgt_ri <- 873;  label("Mean generation time, IPM-resistant/AGS-intermediate (min)")                # Table 3: k12,RI row = 873 (SE 8.5%) [strain-shared]
    mgt_ir <- 26.9; label("Mean generation time, IPM-intermediate/AGS-resistant (min)")                # Table 3: k12,IR row = 26.9 (SE 6.2%) [strain-shared]

    log10mf_ipm <- -4.73; label("Log10 imipenem mutation frequency (seeds RI subpopulation)")          # Table 3: Log MUT,IPM = -4.73 (SE 5.8%) [strain-shared]
    log10mf_ags <- -6.68; label("Log10 amikacin mutation frequency (seeds IR subpopulation)")          # Table 3 footnote b: Log MUT,AGS = -6.68 (SE 3.4%)

    kmax_ipm    <- 3.34;  label("Maximum imipenem killing rate constant (1/h)")                        # Table 3: Kmax,IPM = 3.34 (SE 5.5%) [strain-shared]
    kc50_ss_ipm <- 0.992; label("Imipenem KC50, SS subpop (mg/L)")                                     # Table 3: KC50,SS,IPM = 0.992 (SE 33%) [strain-shared]
    kc50_ri_ipm <- 264;   label("Imipenem KC50, RI subpop (mg/L)")                                     # Table 3: KC50,RI,IPM = 264 (SE 9.7%) [strain-shared]
    kc50_ir_ipm <- 23.5;  label("Imipenem KC50, IR subpop (mg/L)")                                     # Table 3: KC50,IR,IPM = 23.5 (SE 11%) [strain-shared]
    hill_ipm    <- 3.00;  label("Hill coefficient for imipenem killing (unitless)")                    # Table 3: Hill,IPM = 3.00 (SE 13.8%) [strain-shared]

    kmax_ss_ags <- 11.8;  label("Maximum AGS killing rate, SS subpop (1/h)")                           # Table 3: Kmax,SS,AGS = 11.8 (SE 7.9%) [shared across TOB/AMK]
    kmax_ri_ags <- 3.28;  label("Maximum AGS killing rate, RI subpop (1/h)")                           # Table 3: Kmax,RI,AGS = 3.28 (SE 26.3%) [shared]
    kmax_ir_ags <- 11.8;  label("Maximum AGS killing rate, IR subpop (1/h)")                           # Table 3: Kmax,IR,AGS = 11.8 (SE 7.9%) [shared]
    kc50_ss_ags <- 4.88;  label("Amikacin KC50, SS subpop (mg/L)")                                     # Table 3 footnote b: KC50,SS,AGS = 4.88 (SE 25.3%)
    kc50_ri_ags <- 280;   label("Amikacin KC50, RI subpop (mg/L)")                                     # Table 3 footnote b: KC50,RI,AGS = 280 (SE 15.7%)
    kc50_ir_ags <- 329;   label("Amikacin KC50, IR subpop (mg/L)")                                     # Table 3 footnote b: KC50,IR,AGS = 329 (SE 10.8%)
    hill_ags    <- 1.11;  label("Hill coefficient for AGS killing (unitless)")                         # Table 3: Hill,AGS = 1.11 (SE 5.5%) [shared]

    mtt_sig      <- fixed(1);    label("Mean turnover time of hypothetical signal molecule (h; FIXED)") # Table 3: MTT,sig = 1 (fixed) [strain-shared]
    imax_sig12   <- 0.997;       label("Max fractional inhibition of k12 by signal molecule (unitless)") # Table 3: Imax,sig12 = 0.997 (SE 13.7%) [strain-shared]
    log10ic50sig <- 10.2;        label("Log10 signal molecule concn at 50% Imax,sig12 (log10 CFU/mL)")  # Table 3: Log IC50,Sig = 10.2 (SE 4.3%) [strain-shared]

    imax_om <- 0.678; label("Max fractional decrease of KC50,IPM via OM disruption (unitless)")         # Table 3 footnote c: Imax,OM,AGS = 0.678 (SE 31.8%) [strain-shared]
    ic50_om <- 1.13;  label("Amikacin concn for 50% Imax,OM (mg/L)")                                   # Table 3 footnote b: IC50,OM,AGS = 1.13 (SE 49.0%)

    addSd <- 0.378; label("Additive residual SD on log10 scale (log10 CFU/mL)")                         # Table 3: SD,CFU = 0.378 (SE 6.0%) [strain-shared]
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
