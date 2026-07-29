Landersdorfer_2018_meropenem_tobramycin_PAOmutS <- function() {
  description <- "In vitro (Pseudomonas aeruginosa PAOdelta-mutS hypermutable strain). Mechanism-based PK/PD (life-cycle growth) model of bacterial killing and resistance for meropenem plus tobramycin, with three pre-existing subpopulations (susceptible, MEM-resistant/TOB-intermediate, MEM-intermediate/TOB-resistant) and mechanistic synergy via tobramycin-induced outer-membrane permeabilisation of meropenem"
  reference <- "Landersdorfer CB, Rees VE, Yadav R, Rogers KE, Kim TH, Bergen PJ, Cheah SE, Boyce JD, Peleg AY, Oliver A, Shin BS, Nation RL, Bulitta JB. Optimization of a meropenem-tobramycin combination dosage regimen against hypermutable and nonhypermutable Pseudomonas aeruginosa via mechanism-based modeling and the hollow-fiber infection model. Antimicrob Agents Chemother. 2018 Mar 27;62(4):e02055-17. doi:10.1128/AAC.02055-17. Model differential equations (Eqs 2-5) and Table S3 population parameter estimates are in the supplemental material."
  vignette <- "Landersdorfer_2018_meropenem_tobramycin"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  paper_specific_compartment_pattern <- "^bact_"

  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa PAOdelta-mutS hypermutable strain; isogenic mutS deletion of PAO1 by Mena et al.)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Hypermutable P. aeruginosa relevant to cystic fibrosis chronic lung infection; baseline MIC meropenem 1 mg/L, MIC tobramycin 0.5 mg/L (identical pre-treatment MICs to PAO1)",
    model_system     = "96-h static-concentration time-kill (SCTK) experiments; meropenem 2, 8, 16 mg/L and tobramycin 1, 4, 8 mg/L alone and in combination",
    initial_inoculum = "~10^7.8 CFU/mL targeted; 10^7.72 CFU/mL estimated",
    dose_range       = "Meropenem 2-16 mg/L static concentrations; tobramycin 1-8 mg/L static concentrations; alone and combined",
    notes            = paste(
      "Mechanism-based model (S-ADAPT, importance sampling) co-modeled SCTK data for PAOdelta-mutS and its isogenic wild-type sibling PAO1.",
      "Strain-specific parameter estimates for PAOdelta-mutS are reported in Table S3 column 'b'.",
      "PAOdelta-mutS has a defective DNA-mismatch repair system (~1000-fold higher mutation rate than PAO1); resistance emerges via rapid ascent of less-susceptible mutants.",
      "Predicted in silico for hollow-fiber infection model (HFIM) regimens using clinically relevant CF-patient meropenem and tobramycin concentration-time profiles (simulated CL_MEM = 15.9 L/h, t1/2 = 0.8 h; CL_TOB = 4.9 L/h, t1/2 = 2.5 h).",
      "Sibling model: Landersdorfer_2018_meropenem_tobramycin_PAO1."
    )
  )

  ini({
    # --- Bacterial growth and subpopulations (Table S3) ---
    log10cfu0    <- 7.72;  label("Initial inoculum (log10 CFU/mL)")                                # Table S3: Log10CFU0 = 7.72 (SE 1.0%); shared between strains
    log10cfumax  <- 9.57;  label("Maximum population size CFUmax (log10 CFU/mL)")                  # Table S3: Log10CFUmax = 9.57 (SE 1.4%); shared between strains
    lk21         <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)")      # Same-group precedent (Rees 2018, AAC); supplement states "k21 was assumed to be fast"

    # Mean generation time per subpopulation (minutes); growth rate k12 = 60/MGT (1/h)
    mgt_ss <- 57.3;  label("Mean generation time, double-susceptible MEMs/TOBs (min)")             # Table S3: MGT_SS = 57.3 (SE 11.2%); shared between strains
    mgt_ri <- 2341;  label("Mean generation time, MEM-resistant/TOB-intermediate (min)")           # Table S3: MGT_RI = 2341 (SE 10.5%); shared between strains
    mgt_ir <- 78.4;  label("Mean generation time, MEM-intermediate/TOB-resistant (min)")           # Table S3: MGT_IR = 78.4 (SE 4.7%); shared between strains

    # Log10 pre-existing fraction of less-susceptible subpopulations (seeds the initial conditions)
    log10prb_ri  <- -7.63; label("Log10 proportion of MEM-resistant/TOB-intermediate bacteria (seeds RI subpop)") # Table S3: Log10PRB_RI = -7.63 (SE 2.2%) for PAOdelta-mutS
    log10prb_ir  <- -2.99; label("Log10 proportion of MEM-intermediate/TOB-resistant bacteria (seeds IR subpop)") # Table S3: Log10PRB_IR = -2.99 (SE 4.3%) for PAOdelta-mutS

    # --- Killing by meropenem (Table S3) ---
    kmax_mem_ss  <- 2.81;  label("Maximum meropenem killing rate constant, double-susceptible (1/h)")    # Table S3: Kmax_MEM_SS = 2.81 (SE 17.9%) for PAOdelta-mutS
    kmax_mem_ri  <- 0.458; label("Maximum meropenem killing rate constant, MEM-resistant/TOB-int (1/h)") # Table S3: Kmax_MEM_RI = 0.458 (SE 9.5%) for PAOdelta-mutS
    kmax_mem_ir  <- 1.69;  label("Maximum meropenem killing rate constant, MEM-int/TOB-resistant (1/h)") # Table S3: Kmax_MEM_IR = 1.69 (SE 12.4%) for PAOdelta-mutS
    kc50_ss_mem  <- 2.07;  label("Meropenem KC50, double-susceptible (mg/L)")                            # Table S3: KC50_MEM_SS = 2.07 (SE 16%); shared between strains
    kc50_ri_mem  <- 81.2;  label("Meropenem KC50, MEM-resistant/TOB-intermediate (mg/L)")                # Table S3: KC50_MEM_RI = 81.2 (SE 17.2%); shared between strains
    kc50_ir_mem  <- 40.0;  label("Meropenem KC50, MEM-intermediate/TOB-resistant (mg/L)")                # Table S3: KC50_MEM_IR = 40.0 (SE 13.9%) for PAOdelta-mutS
    hill_mem     <- 0.563; label("Hill coefficient for meropenem killing (unitless)")                    # Table S3: HillMEM = 0.563 (SE 18.9%); shared between strains

    # --- Killing by tobramycin (Table S3) ---
    kmax_tob_ss  <- 5.73;  label("Maximum tobramycin killing rate constant, double-susceptible (1/h)")     # Table S3: Kmax_TOB_SS = 5.73 (SE 12.3%) for PAOdelta-mutS
    kmax_tob_ri  <- 0.182; label("Maximum tobramycin killing rate constant, MEM-resistant/TOB-int (1/h)")  # Table S3: Kmax_TOB_RI = 0.182 (SE 23.4%) for PAOdelta-mutS
    kmax_tob_ir  <- 0.156; label("Maximum tobramycin killing rate constant, MEM-int/TOB-resistant (1/h)")  # Table S3: Kmax_TOB_IR = 0.156 (SE 21.8%); shared between strains
    kc50_ss_tob  <- 1.62;  label("Tobramycin KC50, double-susceptible (mg/L)")                             # Table S3: KC50_TOB_SS = 1.62 (SE 20.2%); shared between strains
    kc50_ri_tob  <- 31.0;  label("Tobramycin KC50, MEM-resistant/TOB-intermediate (mg/L)")                 # Table S3: KC50_TOB_RI = 31.0 (SE 14.8%); shared between strains
    kc50_ir_tob  <- 92.8;  label("Tobramycin KC50, MEM-intermediate/TOB-resistant (mg/L)")                 # Table S3: KC50_TOB_IR = 92.8 (SE 7.4%); shared between strains
    hill_tob     <- 3.03;  label("Hill coefficient for tobramycin killing (unitless)")                     # Table S3: HillTOB = 3.03 (SE 12.9%); shared between strains

    # --- Mechanistic synergy: tobramycin disrupts the bacterial outer membrane, lowering effective KC50,MEM (Eq 5) ---
    imax_om <- fixed(1.0); label("Max fractional decrease of KC50,MEM via outer-membrane disruption (unitless; FIXED)") # Table S3: Imax_OM = 1.0 (footnote e: fixed)
    ic50_om <- 0.587;      label("Tobramycin concentration for 50% of Imax_OM (mg/L)")                                  # Table S3: IC50_OM = 0.587 (SE 19.5%); shared between strains

    # --- Residual error (additive on log10 scale) ---
    addSd <- 0.41; label("Additive residual SD on log10 scale (log10 CFU/mL)") # Table S3: SDCFU = 0.41 (SE 5.6%); shared between strains

    # --- Simulated antibiotic disposition (fixed, from clinical CF-patient PK literature) ---
    thalf_mem <- fixed(0.8); label("Simulated meropenem half-life (h; FIXED)") # Main text Methods: t1/2,MEM = 0.8 h (refs 55, 56); not an MBM estimate
    thalf_tob <- fixed(2.5); label("Simulated tobramycin half-life (h; FIXED)") # Main text Methods: t1/2,TOB = 2.5 h (refs 55, 56); not an MBM estimate
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

    # 3. Antibiotic elimination rate constants from the simulated half-lives
    kel_mem <- log(2) / thalf_mem    # 1/h
    kel_tob <- log(2) / thalf_tob    # 1/h

    # 4. Mechanistic synergy factor (Eq 5): tobramycin lowers effective KC50,MEM.
    #    OM_effect ranges from 1 (no TOB) to 1 - imax_om (saturating TOB); the
    #    effective KC50,MEM in each killing term is OM_effect * KC50,MEM, so
    #    1/OM_effect is the fold-decrease in KC50,MEM. HillOM is taken as 1
    #    (Emax form); the supplement renders Eq 5 with a HillOM exponent but
    #    Table S3 reports no estimate (sibling Rees 2018 also uses Hill=1).
    om_effect <- 1 - imax_om * ctob / (ctob + ic50_om)

    # 5. Subpopulation killing rates (Eq 2 killing term): Hill-type killing by
    #    both meropenem (with synergy on KC50,MEM) and tobramycin.
    kill_ss <- kmax_mem_ss * cmem^hill_mem / (cmem^hill_mem + (om_effect * kc50_ss_mem)^hill_mem) +
      kmax_tob_ss * ctob^hill_tob / (ctob^hill_tob + kc50_ss_tob^hill_tob)
    kill_ri <- kmax_mem_ri * cmem^hill_mem / (cmem^hill_mem + (om_effect * kc50_ri_mem)^hill_mem) +
      kmax_tob_ri * ctob^hill_tob / (ctob^hill_tob + kc50_ri_tob^hill_tob)
    kill_ir <- kmax_mem_ir * cmem^hill_mem / (cmem^hill_mem + (om_effect * kc50_ir_mem)^hill_mem) +
      kmax_tob_ir * ctob^hill_tob / (ctob^hill_tob + kc50_ir_tob^hill_tob)

    # 6. Replication factor (Eq 3): REP = 2 at low CFU, REP = 1 when CFU = CFUmax.
    #    Equivalent to 2 * (1 - CFUall / (2 * CFUmax)); stationary plateau at CFUmax.
    CFUall <- bact_susceptible_susceptible1 + bact_susceptible_susceptible2 +
      bact_resistant_intermediate1 + bact_resistant_intermediate2 +
      bact_intermediate_resistant1 + bact_intermediate_resistant2
    rep_factor <- 2 - CFUall / cfumax

    # 7. Life-cycle growth model: two states per subpopulation (Eqs 2, 4).
    #    State 1 = preparing for replication; state 2 = immediately before
    #    replication. State-2 cells replicate at k21 producing rep_factor daughter
    #    state-1 cells (rep_factor in [1, 2] over the population trajectory).
    d/dt(bact_susceptible_susceptible1) <- rep_factor * k21 * bact_susceptible_susceptible2 - k12ss * bact_susceptible_susceptible1 - kill_ss * bact_susceptible_susceptible1
    d/dt(bact_susceptible_susceptible2) <- -k21 * bact_susceptible_susceptible2 + k12ss * bact_susceptible_susceptible1 - kill_ss * bact_susceptible_susceptible2
    d/dt(bact_resistant_intermediate1)  <- rep_factor * k21 * bact_resistant_intermediate2  - k12ri * bact_resistant_intermediate1  - kill_ri * bact_resistant_intermediate1
    d/dt(bact_resistant_intermediate2)  <- -k21 * bact_resistant_intermediate2 + k12ri * bact_resistant_intermediate1 - kill_ri * bact_resistant_intermediate2
    d/dt(bact_intermediate_resistant1)  <- rep_factor * k21 * bact_intermediate_resistant2  - k12ir * bact_intermediate_resistant1  - kill_ir * bact_intermediate_resistant1
    d/dt(bact_intermediate_resistant2)  <- -k21 * bact_intermediate_resistant2 + k12ir * bact_intermediate_resistant1 - kill_ir * bact_intermediate_resistant2

    # 8. Antibiotic concentrations (mg/L): dosed by the user, decline per half-life.
    d/dt(cmem) <- -kel_mem * cmem
    d/dt(ctob) <- -kel_tob * ctob

    # 9. Initial conditions: pre-existing subpopulations seeded by Log10PRB;
    #    all inoculum placed in state 1 (state-2 equilibrates within ~1/k12).
    #    cmem and ctob start at 0.
    bact_susceptible_susceptible1(0) <- cfu0
    bact_resistant_intermediate1(0)  <- cfu0 * 10^log10prb_ri
    bact_intermediate_resistant1(0)  <- cfu0 * 10^log10prb_ir

    # 10. Outputs: total viable counts on log10 scale; fit on the log10 scale.
    CFUri <- bact_resistant_intermediate1 + bact_resistant_intermediate2    # MEM-resistant subpopulation (CFU/mL)
    CFUir <- bact_intermediate_resistant1 + bact_intermediate_resistant2    # TOB-resistant subpopulation (CFU/mL)
    Cc <- log10(CFUall)
    Cc ~ add(addSd)
  })
}
