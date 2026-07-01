Nielsen_2007_semimechanistic_antibiotic_pd <- function() {
  description <- "In vitro (Streptococcus pyogenes M12 NCTC P1800). Semimechanistic time-kill pharmacokinetic/pharmacodynamic model of five antibiotics (benzylpenicillin, cefuroxime, erythromycin, moxifloxacin, vancomycin) against S. pyogenes. The bacterial system has two states: a growing drug-susceptible population (bact_susceptible) and a resting drug-insusceptible population (bact_resting) that captures phenotypic persister-cell switching. Bacteria grow in the susceptible state at first-order rate kgrowth, die in both states at first-order rate kdeath, and transfer from susceptible to resting at rate kSR = (kgrowth - kdeath) * (bact_susceptible + bact_resting) / Bmax so the total population reaches Bmax at stationary phase (reverse transfer kRS is fixed to 0). Each drug is dosed into its own PK compartment (pen, cxm, ery, mxf, van; the compartment state IS the bath concentration in mg/L per the in-vitro convention). Drug decays first-order via degradation (kdeg fixed from stability experiments; nonzero only for benzylpenicillin and cefuroxime). A biophase (effect) compartment (pen_e, cxm_e, ery_e, mxf_e, van_e) equilibrates with the PK compartment at first-order rate ke and drives the killing effect through a sigmoidal Emax function DRUG = Emax * Ce^gamma / (Ce^gamma + EC50^gamma). DRUG adds to the natural death rate on susceptible bacteria only (paper equation 6). Multi-drug DRUG contributions sum; in monotherapy (as fitted) only one drug is active. The starting inoculum is treated as mix1 (all bacteria in the growing state) per the paper's mixture-model dominant mode; parameters fmix1 (0.747) and fpers (0.0529) are retained as fixed documentation for the mix2 alternative starting condition."
  reference <- "Nielsen EI, Viberg A, Lowdin E, Cars O, Karlsson MO, Sandstrom M. Semimechanistic pharmacokinetic/pharmacodynamic model for assessment of activity of antibacterial agents from time-kill curve experiments. Antimicrob Agents Chemother. 2007 Jan;51(1):128-136. doi:10.1128/AAC.00604-06. PMID: 17060527."
  vignette <- "Nielsen_2007_semimechanistic_antibiotic_pd"
  units <- list(time = "hour", dosing = "mg/L (initial concentration)", concentration = "log10 CFU/mL (observation); mg/L (drug states)")

  # No patient covariates: this is an in-vitro semi-mechanistic PD model with
  # static drug exposures. Drug initial concentrations are applied via dosing
  # events into the pen / cxm / ery / mxf / van compartments at time 0 with
  # amt interpreted as mg/L (consistent with the in-vitro convention where
  # the compartment state IS the bath concentration). Paper-studied ranges:
  # benzylpenicillin 0.00075-0.77 mg/L (0.0625-64x MIC 0.012 mg/L),
  # cefuroxime 0.00196-2.00 mg/L (0.0625-64x MIC 0.0313 mg/L),
  # erythromycin 0.0078-8.00 mg/L (0.0625-64x MIC 0.125 mg/L),
  # moxifloxacin 0.0313-8.00 mg/L (0.25-64x MIC 0.125 mg/L),
  # vancomycin 0.0625-16.0 mg/L (0.25-64x MIC 0.25 mg/L).
  covariateData <- list()

  # Paper-specific compartments: five drug PK compartments (one per antibiotic
  # studied), five biophase / effect compartments, and two bacterial-state
  # compartments (susceptible vs. resting). The drug PK / biophase compartments
  # are not classical central / peripheral compartments because each represents
  # the in-vitro medium concentration of a single antibiotic; there is no
  # subject-level PK.
  paper_specific_compartments <- c(
    "pen", "cxm", "ery", "mxf", "van",
    "pen_e", "cxm_e", "ery_e", "mxf_e", "van_e",
    "bact_susceptible", "bact_resting"
  )

  population <- list(
    species             = "in vitro (Streptococcus pyogenes group A M12 strain NCTC P1800)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Streptococcus pyogenes group A M12 NCTC P1800 (single reference strain)",
    system              = "Static time-kill curve experiments in 10 mL glass tubes with 4 mL Todd-Hewitt broth; 24-hour incubation with dense sampling at 0, 1, 2, 4, 6, 9, 12, 15, 18, and 24 h",
    medium              = "Todd-Hewitt broth (bacteria); MIC determinations on Iso-Sensitest agar",
    temperature         = "35 C incubation; 5% CO2 for colony counting",
    duration            = "24 h",
    starting_inoculum   = "10^6 CFU/mL (standard); additional lower inocula used only for antibiotic-free growth controls",
    limit_of_detection  = "10 CFU/mL",
    mic_values          = c(
      benzylpenicillin = "0.012 mg/L",
      cefuroxime       = "0.0313 mg/L",
      erythromycin     = "0.125 mg/L",
      moxifloxacin     = "0.125 mg/L",
      vancomycin       = "0.25 mg/L"
    ),
    concentration_range = c(
      benzylpenicillin = "0.0625x to 64x MIC (0.00075-0.77 mg/L); 455 observations",
      cefuroxime       = "0.0625x to 64x MIC (0.00196-2.00 mg/L); 427 observations",
      erythromycin     = "0.0625x to 64x MIC (0.0078-8.00 mg/L); 455 observations",
      moxifloxacin     = "0.25x to 64x MIC (0.0313-8.00 mg/L); 376 observations",
      vancomycin       = "0.25x to 64x MIC (0.0625-16.0 mg/L); 409 observations"
    ),
    regimens            = "Monotherapy only; each experiment exposed the bacterial inoculum to a single antibiotic at a single fixed concentration for 24 h. Each experiment was run in duplicate or triplicate on separate days; at least one drug-free growth control was included per day.",
    notes               = "Total of 135 time-kill experiments were fitted simultaneously in a single NONMEM ADVAN9 / FOCE analysis. Population-analysis structure used a 2-component residual model (Karlsson 1995) with replicate-specific (repl) and consistent-across-replicates (eps) components. See Nielsen 2007 Materials and Methods (p 129) and Tables 1-3."
  )

  ini({
    # --- Shared bacterial-system parameters (Table 2, typical values with RSE %) ---
    lkgrowth <- log(1.35)
    label("Log growth rate constant of susceptible bacteria (kgrowth, 1/h)")  # Table 2: kgrowth = 1.35 (RSE 5.4%)
    lkdeath <- log(0.179)
    label("Log natural death rate constant, both bacterial states (kdeath, 1/h)")  # Table 2: kdeath = 0.179 (RSE 6.5%)
    lbmax <- log(4.15e8)
    label("Log maximum bacterial concentration at stationary phase (Bmax, CFU/mL)")  # Table 2: Bmax = 4.15e8 (RSE 9.2%)

    # Mixture-model parameters describing starting-inoculum heterogeneity.
    # These are fixed documentation of the paper's mixture module: mix1 (log
    # phase, all susceptible) contained fmix1 = 74.7% of the 23 starting
    # inocula, and mix2 (approaching stationary, some resting cells) had fpers
    # = 5.29% of cells already in the resting state at t = 0. The default
    # simulation initial condition below assumes mix1 (dominant mode); a user
    # can override via the compartment initial-state ini fields if a mix2
    # start is required.
    fmix1 <- fixed(0.747)
    label("Fraction of starting inocula in mix1 -- log-phase, no resting cells (unitless)")  # Table 2: fmix1 = 0.747 (RSE 16%)
    fpers <- fixed(0.0529)
    label("Fraction of resting cells in the mix2 starting inoculum (unitless)")  # Table 2: fpers = 0.0529 (RSE 48%)

    # --- Drug-specific parameters (Table 3, typical values with RSE %) ---

    # Benzylpenicillin (Emax = 2.44 1/h, EC50 = 0.00438 mg/L, gamma = 1.29, ke = 1.00 1/h)
    lemax_pen <- log(2.44)
    label("Log benzylpenicillin maximum killing rate constant (Emax_pen, 1/h)")  # Table 3: Emax = 2.44 (RSE 8.6%)
    lec50_pen <- log(0.00438)
    label("Log benzylpenicillin half-maximum-effect biophase concentration (EC50_pen, mg/L)")  # Table 3: EC50 = 0.00438 (RSE 7.7%)
    lhill_pen <- log(1.29)
    label("Log benzylpenicillin Hill sigmoidicity exponent (gamma_pen, unitless)")  # Table 3: gamma = 1.29 (RSE 10%)
    lke0_pen <- log(1.00)
    label("Log benzylpenicillin biophase equilibration rate constant (ke_pen, 1/h)")  # Table 3: ke = 1.00 (RSE 9.6%)
    lkdeg_pen <- fixed(log(0.020))
    label("Log benzylpenicillin first-order degradation rate constant in broth (kdeg_pen, 1/h; FIXED from stability study)")  # Results p 132: kdeg = 0.020 h-1 (FIXED)

    # Cefuroxime (Emax = 3.30 1/h, EC50 = 0.00829 mg/L, gamma = 1.69, ke = 0.861 1/h)
    lemax_cxm <- log(3.30)
    label("Log cefuroxime maximum killing rate constant (Emax_cxm, 1/h)")  # Table 3: Emax = 3.30 (RSE 6.1%)
    lec50_cxm <- log(0.00829)
    label("Log cefuroxime half-maximum-effect biophase concentration (EC50_cxm, mg/L)")  # Table 3: EC50 = 0.00829 (RSE 6.6%)
    lhill_cxm <- log(1.69)
    label("Log cefuroxime Hill sigmoidicity exponent (gamma_cxm, unitless)")  # Table 3: gamma = 1.69 (RSE 8.5%)
    lke0_cxm <- log(0.861)
    label("Log cefuroxime biophase equilibration rate constant (ke_cxm, 1/h)")  # Table 3: ke = 0.861 (RSE 17%)
    lkdeg_cxm <- fixed(log(0.026))
    label("Log cefuroxime first-order degradation rate constant in broth (kdeg_cxm, 1/h; FIXED from stability study)")  # Results p 132: kdeg = 0.026 h-1 (FIXED)

    # Erythromycin (Emax = 2.03 1/h, EC50 = 0.0276 mg/L, gamma = 0.769, ke fixed at 100 1/h,
    # kdeg = 0 because no significant degradation over 24 h per the stability study)
    lemax_ery <- log(2.03)
    label("Log erythromycin maximum killing rate constant (Emax_ery, 1/h)")  # Table 3: Emax = 2.03 (RSE 6.4%)
    lec50_ery <- log(0.0276)
    label("Log erythromycin half-maximum-effect biophase concentration (EC50_ery, mg/L)")  # Table 3: EC50 = 0.0276 (RSE 15%)
    lhill_ery <- log(0.769)
    label("Log erythromycin Hill sigmoidicity exponent (gamma_ery, unitless)")  # Table 3: gamma = 0.769 (RSE 19%)
    lke0_ery <- fixed(log(100))
    label("Log erythromycin biophase equilibration rate constant (ke_ery, 1/h; FIXED -- data did not support estimation)")  # Table 3: ke = 100 (FIXED); Results p 132
    lkdeg_ery <- fixed(log(1e-6))
    label("Log erythromycin first-order degradation rate constant in broth (kdeg_ery, 1/h; FIXED negligible)")  # Results p 132: no significant degradation over 24 h

    # Moxifloxacin (Emax = 3.20 1/h, EC50 = 0.0747 mg/L, gamma = 1.59, ke = 0.644 1/h)
    lemax_mxf <- log(3.20)
    label("Log moxifloxacin maximum killing rate constant (Emax_mxf, 1/h)")  # Table 3: Emax = 3.20 (RSE 4.6%)
    lec50_mxf <- log(0.0747)
    label("Log moxifloxacin half-maximum-effect biophase concentration (EC50_mxf, mg/L)")  # Table 3: EC50 = 0.0747 (RSE 3.0%)
    lhill_mxf <- log(1.59)
    label("Log moxifloxacin Hill sigmoidicity exponent (gamma_mxf, unitless)")  # Table 3: gamma = 1.59 (RSE 7.2%)
    lke0_mxf <- log(0.644)
    label("Log moxifloxacin biophase equilibration rate constant (ke_mxf, 1/h)")  # Table 3: ke = 0.644 (RSE 20%)
    lkdeg_mxf <- fixed(log(1e-6))
    label("Log moxifloxacin first-order degradation rate constant in broth (kdeg_mxf, 1/h; FIXED negligible)")  # Results p 132: no significant degradation over 24 h

    # Vancomycin (Emax = 1.36 1/h, EC50 = 0.384 mg/L, gamma = 20 FIXED, ke fixed at 100 1/h,
    # kdeg = 0 because no significant degradation over 24 h)
    lemax_van <- log(1.36)
    label("Log vancomycin maximum killing rate constant (Emax_van, 1/h)")  # Table 3: Emax = 1.36 (RSE 5.5%)
    lec50_van <- log(0.384)
    label("Log vancomycin half-maximum-effect biophase concentration (EC50_van, mg/L)")  # Table 3: EC50 = 0.384 (RSE 0.9%)
    lhill_van <- fixed(log(20))
    label("Log vancomycin Hill sigmoidicity exponent (gamma_van, unitless; FIXED -- very steep all-or-nothing effect)")  # Table 3: gamma = 20 (FIXED at lowest value that did not harm the fit)
    lke0_van <- fixed(log(100))
    label("Log vancomycin biophase equilibration rate constant (ke_van, 1/h; FIXED -- data did not support estimation)")  # Table 3: ke = 100 (FIXED); Results p 132
    lkdeg_van <- fixed(log(1e-6))
    label("Log vancomycin first-order degradation rate constant in broth (kdeg_van, 1/h; FIXED negligible)")  # Results p 132: no significant degradation over 24 h

    # --- Residual error ---
    # Nielsen 2007 used the Karlsson 1995 two-component residual model with
    # replicate-specific (eps_repl = 47%) and consistent-across-replicates
    # (eps = 98%) components estimated on natural-log-transformed CFU/mL. For
    # single-observation simulation the combined SD on the natural-log scale is
    # sqrt(0.98^2 + 0.47^2) = 1.087; converted to the log10-CFU/mL observation
    # scale used here (Cc = log10(bact_susceptible + bact_resting)) this is
    # 1.087 / log(10) = 0.472. Held FIXED as a simulation-oriented summary of
    # the two-component fit rather than re-estimated.
    addSd <- fixed(0.472)
    label("Additive residual SD on log10(CFU/mL); combined eps + eps_repl on log10 scale")  # Table 2: eps = 98%, eps_repl = 47% (both natural-log SDs); combined SD_log10 = sqrt(0.98^2 + 0.47^2) / log(10)
  })

  model({
    # 1. Back-transform log-parameterised primary parameters.
    kgrowth <- exp(lkgrowth)
    kdeath  <- exp(lkdeath)
    bmax    <- exp(lbmax)

    emax_pen <- exp(lemax_pen); ec50_pen <- exp(lec50_pen); hill_pen <- exp(lhill_pen)
    ke0_pen  <- exp(lke0_pen);  kdeg_pen <- exp(lkdeg_pen)

    emax_cxm <- exp(lemax_cxm); ec50_cxm <- exp(lec50_cxm); hill_cxm <- exp(lhill_cxm)
    ke0_cxm  <- exp(lke0_cxm);  kdeg_cxm <- exp(lkdeg_cxm)

    emax_ery <- exp(lemax_ery); ec50_ery <- exp(lec50_ery); hill_ery <- exp(lhill_ery)
    ke0_ery  <- exp(lke0_ery);  kdeg_ery <- exp(lkdeg_ery)

    emax_mxf <- exp(lemax_mxf); ec50_mxf <- exp(lec50_mxf); hill_mxf <- exp(lhill_mxf)
    ke0_mxf  <- exp(lke0_mxf);  kdeg_mxf <- exp(lkdeg_mxf)

    emax_van <- exp(lemax_van); ec50_van <- exp(lec50_van); hill_van <- exp(lhill_van)
    ke0_van  <- exp(lke0_van);  kdeg_van <- exp(lkdeg_van)

    # 2. Per-drug sigmoidal Emax DRUG term (paper equation 4).
    #    DRUG = Emax * Ce^gamma / (Ce^gamma + EC50^gamma).
    #    Each drug's contribution is zero when its biophase concentration is
    #    zero, so in monotherapy (as the paper fitted) only one term is
    #    non-zero. In hypothetical multi-drug simulation the contributions
    #    sum additively into kdeath -- not validated against paper data.
    drug_pen <- emax_pen * pen_e^hill_pen / (pen_e^hill_pen + ec50_pen^hill_pen)
    drug_cxm <- emax_cxm * cxm_e^hill_cxm / (cxm_e^hill_cxm + ec50_cxm^hill_cxm)
    drug_ery <- emax_ery * ery_e^hill_ery / (ery_e^hill_ery + ec50_ery^hill_ery)
    drug_mxf <- emax_mxf * mxf_e^hill_mxf / (mxf_e^hill_mxf + ec50_mxf^hill_mxf)
    drug_van <- emax_van * van_e^hill_van / (van_e^hill_van + ec50_van^hill_van)
    drug_total <- drug_pen + drug_cxm + drug_ery + drug_mxf + drug_van

    # 3. Transfer rate from susceptible to resting (paper Materials and Methods,
    #    p 130). kSR is linear in the total bacterial concentration; the
    #    proportionality constant is (kgrowth - kdeath) / Bmax so that at
    #    stationary phase (bact_susceptible + bact_resting = Bmax) the
    #    susceptible-population net growth rate is zero. kRS is fixed to 0
    #    per the paper's assumption of negligible back-transfer during a
    #    constant-antibiotic-pressure experiment.
    cfu_total <- bact_susceptible + bact_resting
    kSR <- (kgrowth - kdeath) * cfu_total / bmax
    kRS <- 0

    # 4. Drug PK compartments (in-vitro bath concentration; first-order
    #    degradation with drug-specific kdeg; nonzero only for benzylpenicillin
    #    and cefuroxime).
    d/dt(pen) <- -kdeg_pen * pen
    d/dt(cxm) <- -kdeg_cxm * cxm
    d/dt(ery) <- -kdeg_ery * ery
    d/dt(mxf) <- -kdeg_mxf * mxf
    d/dt(van) <- -kdeg_van * van

    # 5. Biophase / effect compartments (Sheiner 1979 link model, cited as
    #    ref 22). Introduced without affecting the mass balance of the drug
    #    PK compartment. The biophase concentration Ce chases the PK
    #    concentration with first-order rate ke.
    d/dt(pen_e) <- ke0_pen * (pen - pen_e)
    d/dt(cxm_e) <- ke0_cxm * (cxm - cxm_e)
    d/dt(ery_e) <- ke0_ery * (ery - ery_e)
    d/dt(mxf_e) <- ke0_mxf * (mxf - mxf_e)
    d/dt(van_e) <- ke0_van * (van - van_e)

    # 6. Bacterial system (paper equations 5-7; the final model uses
    #    equation 6 -- DRUG additive to kdeath on susceptible bacteria only).
    d/dt(bact_susceptible) <- kgrowth * bact_susceptible -
                              (kdeath + drug_total) * bact_susceptible -
                              kSR * bact_susceptible + kRS * bact_resting
    d/dt(bact_resting)     <- kSR * bact_susceptible -
                              kRS * bact_resting -
                              kdeath * bact_resting

    # 7. Initial conditions: expected-value marginalization over the paper's
    #    mixture model. Mix1 (probability fmix1) contributes all bacteria to
    #    the susceptible state. Mix2 (probability 1 - fmix1) contributes a
    #    fraction fpers to the resting state. The expected resting fraction
    #    across all starting inocula is therefore (1 - fmix1) * fpers ~ 1.3%.
    #    Biophase compartments start at zero (no drug delay before dosing);
    #    drug compartments are dosed by the user event table at t = 0.
    inoc0        <- 1e6
    fpers_marg   <- (1 - fmix1) * fpers
    bact_susceptible(0) <- inoc0 * (1 - fpers_marg)
    bact_resting(0)     <- inoc0 * fpers_marg

    # 8. Observation. Cc is the log10 of the total bacterial concentration
    #    (susceptible + resting), matching the y-axis of Nielsen 2007 Fig 2/3.
    #    A 1e-6 floor prevents log10(0) if the ODE integrator drives the
    #    total to (numerically) zero. Additive residual on log10-CFU scale
    #    is the combined eps + eps_repl summarised above (see ini()).
    Cc <- log10(bact_susceptible + bact_resting + 1e-6)
    Cc ~ add(addSd)
  })
}
