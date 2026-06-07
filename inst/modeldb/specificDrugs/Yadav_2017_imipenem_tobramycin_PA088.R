Yadav_2017_imipenem_tobramycin_PA088 <- function() {
  description <- "In vitro (static-concentration time-kill). Mechanism-based PK/PD (Bulitta life-cycle growth) model of bacterial killing and resistance for imipenem combined with tobramycin against carbapenem- and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolate FADDI-PA088 (MIC_IPM = 16 mg/L, MIC_TOB = 32 mg/L). Three pre-existing bacterial subpopulations with signal-molecule growth inhibition and aminoglycoside-mediated outer-membrane permeabilisation (mechanistic synergy)"
  reference <- "Yadav R, Bulitta JB, Nation RL, Landersdorfer CB. Optimization of synergistic combination regimens against carbapenem- and aminoglycoside-resistant clinical Pseudomonas aeruginosa isolates via mechanism-based pharmacokinetic/pharmacodynamic modeling. Antimicrob Agents Chemother. 2017 Jan;61(1):e01011-16. doi:10.1128/AAC.01011-16. Model differential equations (Eqs 1-5) are in the main paper Methods; parameter estimates for FADDI-PA088 with tobramycin are Table 3 (footnote a). Supplemental Text S1 (not on disk) contains diagnostic plots only."
  vignette <- "Yadav_2017_imipenem_aminoglycoside_pseudomonas"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  paper_specific_compartments <- c("cipm", "cags", "csig")

  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa FADDI-PA088, carbapenem- and tobramycin-resistant clinical isolate)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Double-resistant P. aeruginosa bacteremia (MIC_IPM 16 mg/L = 98th-percentile imipenem; MIC_TOB 32 mg/L = 98th-percentile tobramycin per EUCAST)",
    model_system     = "48-h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton II broth, with imipenem supplemented at 6 and 30 h to offset thermal degradation",
    initial_inoculum = "~10^7.5 CFU/mL",
    dose_range       = "Imipenem 8-36 mg/L, tobramycin 1-32 mg/L, monotherapies and combinations",
    notes            = paste(
      "Mechanism-based model (S-ADAPT, importance sampling pmethod=4) fit jointly to total viable counts for all imipenem + tobramycin and imipenem + amikacin monotherapies and combinations against FADDI-PA088; this file uses the tobramycin-specific parameter estimates (Table 3 footnote a).",
      "Three pre-existing bacterial subpopulations: SS (IPM-susceptible, AGS-susceptible), RI (IPM-resistant, AGS-intermediate), IR (IPM-intermediate, AGS-resistant), each described by a two-state Bulitta life-cycle growth model (state 1 preparing for replication, state 2 immediately before replication, replication rate k21 fixed at 50/h).",
      "Outer-membrane synergy (Eq 5): aminoglycoside permeabilises the outer membrane and reduces the effective imipenem KC50 to OM_effect * KC50_IPM, where OM_effect ranges from 1 (no AGS) to 1 - imax_om at saturating AGS. For FADDI-PA088, KC50_IPM decreases up to 3.11-fold at high tobramycin (Table 3 footnote c).",
      "Signal-molecule growth inhibition (Eq 2): hypothetical signal molecules track CFU_all with mean turnover time MTT_sig and inhibit the k12 state-1->state-2 transition. The Csig turnover ODE structure is not explicitly given in the main paper; it cites Bulitta 2010 (ref 61). The standard form d/dt(csig) = (1/mtt_sig) * (CFU_all - csig), giving csig = CFU_all at steady state, is assumed here. Csig(0) = cfu0 (steady state with inoculum).",
      "Drug concentrations cipm and cags (here tobramycin) are state variables dosed by the user; in the original SCTK experiment they are essentially constant (with thermal-degradation supplementation for imipenem). For Monte Carlo simulation under clinical PK, the authors coupled the MBM with previously reported critically-ill popPK models (ref 36 tobramycin, ref 37 imipenem; not extracted here).",
      "Companion files: Yadav_2017_imipenem_amikacin_PA088 (same isolate, amikacin parameters; Table 3 footnote b), Yadav_2017_imipenem_tobramycin_PA001, Yadav_2017_imipenem_amikacin_PA001, Yadav_2017_imipenem_tobramycin_PA022."
    )
  )

  ini({
    # --- Bacterial growth and subpopulations (Table 3, FADDI-PA088, footnote a = tobramycin fit) ---
    log10cfu0   <- 7.54;  label("Initial inoculum (log10 CFU/mL)")                                     # Table 3: Log CFU0 = 7.54 (SE 3.8%)
    log10cfumax <- 9.56;  label("Maximum population size (log10 CFU/mL); PLAT half-saturation")        # Table 3: CFUmax = 9.56 (SE 1.8%)
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")           # Table 3: k21 = 50 (fixed; rapid replication, ref 64)

    # Mean generation time per subpopulation (minutes); growth rate k12 = 60/MGT (1/h).
    # Per Table 3, MGT for SS and IR are equal in this paper (same first-state turnover);
    # RI is much slower (resistant-mutant fitness cost).
    mgt_ss <- 26.9; label("Mean generation time, SS subpop (min)")                                     # Table 3: k12,SS row = 26.9 (SE 6.2%)
    mgt_ri <- 873;  label("Mean generation time, IPM-resistant/TOB-intermediate (min)")                # Table 3: k12,RI row = 873 (SE 8.5%)
    mgt_ir <- 26.9; label("Mean generation time, IPM-intermediate/TOB-resistant (min)")                # Table 3: k12,IR row = 26.9 (SE 6.2%)

    # Log10 mutation frequencies seeding the two pre-existing resistant subpopulations
    log10mf_ipm <- -4.73; label("Log10 imipenem mutation frequency (seeds RI subpopulation)")          # Table 3: Log MUT,IPM = -4.73 (SE 5.8%)
    log10mf_ags <- -7.00; label("Log10 tobramycin mutation frequency (seeds IR subpopulation)")        # Table 3 footnote a: Log MUT,AGS = -7.00 (SE 5.9%)

    # --- Killing by imipenem (Table 3) ---
    kmax_ipm    <- 3.34;  label("Maximum imipenem killing rate constant (1/h)")                        # Table 3: Kmax,IPM = 3.34 (SE 5.5%)
    kc50_ss_ipm <- 0.992; label("Imipenem KC50, SS subpop (mg/L)")                                     # Table 3: KC50,SS,IPM = 0.992 (SE 33%)
    kc50_ri_ipm <- 264;   label("Imipenem KC50, RI subpop (mg/L)")                                     # Table 3: KC50,RI,IPM = 264 (SE 9.7%)
    kc50_ir_ipm <- 23.5;  label("Imipenem KC50, IR subpop (mg/L)")                                     # Table 3: KC50,IR,IPM = 23.5 (SE 11%)
    hill_ipm    <- 3.00;  label("Hill coefficient for imipenem killing (unitless)")                    # Table 3: Hill,IPM = 3.00 (SE 13.8%)

    # --- Killing by aminoglycoside (here tobramycin; Table 3 footnote a) ---
    # Kmax,AGS and Hill,AGS are shared across tobramycin and amikacin within a strain
    # (single value reported in Table 3); KC50,AGS values are AGS-specific.
    kmax_ss_ags <- 11.8;  label("Maximum AGS killing rate, SS subpop (1/h)")                           # Table 3: Kmax,SS,AGS = 11.8 (SE 7.9%)
    kmax_ri_ags <- 3.28;  label("Maximum AGS killing rate, RI subpop (1/h)")                           # Table 3: Kmax,RI,AGS = 3.28 (SE 26.3%)
    kmax_ir_ags <- 11.8;  label("Maximum AGS killing rate, IR subpop (1/h)")                           # Table 3: Kmax,IR,AGS = 11.8 (SE 7.9%)
    kc50_ss_ags <- 18.6;  label("Tobramycin KC50, SS subpop (mg/L)")                                   # Table 3 footnote a: KC50,SS,AGS = 18.6 (SE 39.6%)
    kc50_ri_ags <- 849;   label("Tobramycin KC50, RI subpop (mg/L)")                                   # Table 3 footnote a: KC50,RI,AGS = 849 (SE 5.7%)
    kc50_ir_ags <- 615;   label("Tobramycin KC50, IR subpop (mg/L)")                                   # Table 3 footnote a: KC50,IR,AGS = 615 (SE 4.1%)
    hill_ags    <- 1.11;  label("Hill coefficient for AGS killing (unitless)")                         # Table 3: Hill,AGS = 1.11 (SE 5.5%)

    # --- Signal-molecule growth inhibition (Eq 2) ---
    mtt_sig      <- fixed(1);    label("Mean turnover time of hypothetical signal molecule (h; FIXED)") # Table 3: MTT,sig = 1 (fixed)
    imax_sig12   <- 0.997;       label("Max fractional inhibition of k12 by signal molecule (unitless)") # Table 3: Imax,sig12 = 0.997 (SE 13.7%)
    log10ic50sig <- 10.2;        label("Log10 signal molecule concn at 50% Imax,sig12 (log10 CFU/mL)")  # Table 3: Log IC50,Sig = 10.2 (SE 4.3%)

    # --- Outer-membrane permeabilisation by AGS (Eq 5) ---
    # imax_om is the maximum fractional decrease of KC50_IPM at saturating AGS:
    # at C_AGS >> ic50_om, OM_effect -> 1 - imax_om, so effective KC50,IPM -> (1-imax_om)*KC50,IPM.
    # Footnote c: imax_om = 0.678 maps to up to 3.11-fold KC50_IPM decrease (1/(1-0.678) = 3.11).
    imax_om <- 0.678; label("Max fractional decrease of KC50,IPM via OM disruption (unitless)")         # Table 3 footnote c: Imax,OM,AGS = 0.678 (SE 31.8%)
    ic50_om <- 5.55;  label("Tobramycin concn for 50% Imax,OM (mg/L)")                                  # Table 3 footnote a: IC50,OM,AGS = 5.55 (SE 29.3%)

    # --- Residual error (additive on log10 scale) ---
    addSd <- 0.378; label("Additive residual SD on log10 scale (log10 CFU/mL)")                         # Table 3: SD,CFU = 0.378 (SE 6.0%)
  })

  model({
    # 0. Back-transform the log replication rate constant.
    k21 <- exp(lk21)

    # 1. Carrying capacity and total inoculum (log10 -> linear, CFU/mL)
    cfumax <- 10^log10cfumax
    cfu0   <- 10^log10cfu0
    ic50sig <- 10^log10ic50sig

    # 2. First-order growth rate constants from mean generation time (k12 = 60/MGT(min))
    k12ss <- 60 / mgt_ss
    k12ri <- 60 / mgt_ri
    k12ir <- 60 / mgt_ir

    # 3. Signal-molecule turnover rate constant
    kout_sig <- 1 / mtt_sig

    # 4. Plateau factor (Yadav 2017 definition; Eq 3-4 footnote):
    #    PLAT = 1 - CFUall/(CFUall + CFUmax). PLAT = 1 when CFUall = 0;
    #    PLAT = 0.5 when CFUall = CFUmax; PLAT -> 0 as CFUall -> infinity.
    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 +
              bact_resistant_intermediate1 + bact_resistant_intermediate2 +
              bact_intermediate_resistant1 + bact_intermediate_resistant2
    plat <- 1 - CFUall / (CFUall + cfumax)

    # 5. Signal-molecule inhibition of k12 (Eq 2):
    #    Inh_k12 = imax_sig12 * csig / (csig + IC50_sig); applied as (1-Inh) on k12.
    inh_k12 <- imax_sig12 * csig / (csig + ic50sig)

    # 6. Outer-membrane synergy (Eq 5): AGS permeabilises the OM and reduces the
    #    effective imipenem KC50 to OM_effect * KC50_IPM. Without AGS, OM_effect=1.
    om_effect <- 1 - imax_om * cags / (cags + ic50_om)

    # 7. Per-subpopulation killing terms (Eq 3 right-hand-side bracket):
    #    Hill-type imipenem killing with synergy on KC50,IPM plus Hill-type AGS killing.
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

    # 8. Bacterial life cycle (Eq 3-4 per subpopulation): state 2 cells replicate
    #    via k21 producing 2 daughters in state 1 with success probability plat;
    #    state 1 -> state 2 at rate k12*(1-inh_k12).
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

    # 9. Signal-molecule turnover (CFU-tracking with delay MTT_sig).
    #    Standard Bulitta life-cycle form; explicit ODE not in main paper text
    #    (Methods cites ref 61 for the formulation; supplement Text S1 not on disk).
    d/dt(csig) <- kout_sig * (CFUall - csig)

    # 10. Antibiotic bath concentrations (mg/L): dosed by the user. Default to
    #     constant (SCTK assumption); user can attach a clinical PK model by
    #     adding elimination terms or by feeding the states via dosing events.
    d/dt(cipm) <- 0
    d/dt(cags) <- 0

    # 11. Initial conditions: total inoculum CFU0 placed in subpopulation 1
    #     state-1; resistant subpopulations seeded by mutation frequency.
    #     csig starts at steady state with the inoculum.
    bact_susceptible_susceptible1(0) <- cfu0
    bact_resistant_intermediate1(0)  <- cfu0 * 10^log10mf_ipm
    bact_intermediate_resistant1(0)  <- cfu0 * 10^log10mf_ags
    csig(0) <- cfu0

    # 12. Outputs: total viable count on the log10 scale (fitted endpoint).
    Cc <- log10(CFUall)
    Cc ~ add(addSd)
  })
}
