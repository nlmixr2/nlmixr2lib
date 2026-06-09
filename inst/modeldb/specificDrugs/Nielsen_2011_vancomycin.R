Nielsen_2011_vancomycin <- function() {
  description <- "In vitro (Streptococcus pyogenes M12 NCTC P1800). Semimechanistic PKPD model of vancomycin time-kill kinetics; two-stage bacterial life-cycle (proliferating drug-sensitive S and non-growing drug-insensitive R) with sigmoidal Emax killing of S via an effect compartment; first-order drug elimination (ke set per in vitro kinetic-system flow rate); drug-specific degradation kdeg fixed at zero. Parameter values are from the combined static and dynamic estimation in Table 3."
  reference <- paste(
    "Nielsen EI, Cars O, Friberg LE. (2011).",
    "Predicting in vitro antibacterial efficacy across experimental designs with a semimechanistic pharmacokinetic-pharmacodynamic model.",
    "Antimicrobial Agents and Chemotherapy 55(4):1571-1579.",
    "doi:10.1128/AAC.01286-10.",
    "Model structure originally developed in Nielsen EI, Viberg A, Lowdin E, Cars O, Karlsson MO, Sandstrom M. (2007). Semimechanistic pharmacokinetic/pharmacodynamic model for assessment of activity of antibacterial agents from time-kill curve experiments. Antimicrob Agents Chemother 51(1):128-136 (reference 31 in Nielsen 2011); the 2011 paper extends the same structure to dynamic concentration-time profiles and reports re-estimated parameters.",
    sep = " "
  )
  vignette <- "Nielsen_2011_antibacterial_efficacy"
  units <- list(time = "hour", dosing = "mg/L (drug input concentration)", concentration = "natural log CFU/mL (observation); mg/L (drug compartment)")

  paper_specific_compartments <- c("bact_sensitive", "bact_resting")

  covariateData <- list()

  population <- list(
    species       = "in vitro (Streptococcus pyogenes group A strain M12 NCTC P1800)",
    n_subjects    = NA_integer_,
    n_studies     = 1L,
    organism      = "Streptococcus pyogenes group A M12 NCTC P1800 (National Culture Type Collection); vancomycin MIC 0.25 mg/L",
    medium        = "Todd-Hewitt broth (35 C, 5% CO2 for plating)",
    inoculum      = "Target 10^6 CFU/mL at t = 0 (logarithmic-growth-phase culture)",
    system        = "Static 10-mL test tubes (4 mL broth) and a dynamic in vitro kinetic system (110-mL open-bottom spinner flask with pump-driven dilution producing first-order antibiotic elimination at flow-rate / volume)",
    duration      = "24 h (some experiments 48 h)",
    drug          = "Vancomycin (Vancomycin Abbott; ElectraBox Pharma)",
    dose_range    = "Static initial concentrations 0.25, 0.5, 1, 1.5, 2, 4, 16, 64 x MIC; dynamic initial concentrations 2 and 16 x MIC with simulated half-lives 0 (constant infusion), 5.1 (human), and 15.2 (one-third human) h",
    n_observations = 580L,
    notes         = "Parameters are typical values from the simultaneous static + dynamic estimation (Table 3 'Static and dynamic' column). Random effects (eta) are NOT included: the source NONMEM run estimated only inter-experiment variability in the resting fraction at t=0 via a two-component mixture model, plus replicate and common additive residual error on the natural-log viable count. The packaged file folds the mixture into a typical starting fraction (E[f_resting] = (1 - fmix1) * fpers) and combines the replicate and common residual SDs into a single additive SD. Note: the paper reports that the model under-predicts the sustained suppression of bacterial growth observed for vancomycin at 2 x MIC in dynamic experiments (see Discussion); users simulating that regimen should expect a discrepancy with the paper's measured time-kill curves at low concentrations."
  )

  ini({
    # =============================================================
    # Bacterial life-cycle parameters (shared across the five drugs)
    # =============================================================
    lkgrowth <- log(1.46)
    label("Log bacterial multiplication rate kgrowth (1/h)")          # Table 3 "Static and dynamic" column: kgrowth = 1.46 (RSE 4.1%)
    lkdeath  <- log(0.187)
    label("Log natural-death rate kdeath (1/h)")                      # Table 3 "Static and dynamic" column: kdeath = 0.187 (RSE 6.2%)
    lbmax    <- log(5.00e8)
    label("Log carrying-capacity Bmax (CFU/mL)")                      # Table 3 "Static and dynamic" column: Bmax = 5.00e8 (RSE 8.7%)
    fpers <- 0.0652
    label("Fraction of bacteria in the resting stage at t = 0 (subpopulation 2; unitless)")  # Table 3 "Static and dynamic" column: fpers = 0.0652 (RSE 49%)
    fmix1 <- 0.880
    label("Proportion of experiments in subpopulation 1 (no resting bacteria at t = 0; unitless)")  # Table 3 "Static and dynamic" column: fMix1 = 0.880 (RSE 7.5%)

    # =============================================================
    # Drug-specific PK parameters (vancomycin)
    # =============================================================
    kdeg <- fixed(0)
    label("Drug degradation rate kdeg in broth (1/h; FIXED at 0 per Methods 'Semimechanistic PKPD model')")  # Methods: "For the other drugs, the kdeg was set to zero" (vancomycin = 0)

    # ke is the in vitro kinetic-system elimination rate. The
    # default below corresponds to a 5.1-h simulated human
    # half-life (Table 1 dynamic-vancomycin column); override via
    # params = list(lke = ...) for a different half-life or to ~ 0
    # (e.g., log(1e-9)) for a static-concentration simulation.
    lke <- fixed(log(log(2) / 5.1))
    label("Log in vitro kinetic-system elimination rate ke (1/h; default ln(2)/5.1 = human half-life; override per experiment)")  # Table 1: simulated human half-life for vancomycin = 5.1 h

    # =============================================================
    # Drug-specific PD parameters (vancomycin; Table 3)
    # =============================================================
    lemax <- log(1.52)
    label("Log maximum killing rate Emax (1/h)")                                          # Table 3 "Static and dynamic" column: Emax VAN = 1.52 (RSE 5.3%)
    lec50 <- log(0.304)
    label("Log effect-compartment concentration for 50% of Emax (mg/L)")                  # Table 3 "Static and dynamic" column: EC50 VAN = 0.304 (RSE 13%)
    lhill <- log(4.99)
    label("Log Hill (sigmoidicity) exponent gamma on the killing function (unitless)")    # Table 3 "Static and dynamic" column: gamma VAN = 4.99 (RSE 24%)
    lke0  <- fixed(log(100))
    label("Log effect-compartment delay rate ke0 (1/h; FIXED at 100 per Table 3 'Static and dynamic' column)")  # Table 3 "Static and dynamic" column: ke0 VAN = 100 FIX (effectively instantaneous equilibrium)

    # =============================================================
    # Starting inoculum
    # =============================================================
    lcfu0 <- fixed(log(1e6))
    label("Log starting bacterial inoculum (CFU/mL; FIXED at 10^6 per Methods)")  # Methods: "starting inoculum of 10^6 CFU/ml"

    # =============================================================
    # Residual error
    # =============================================================
    addSd <- 1.476
    label("Additive residual SD on natural-log CFU/mL (combined common + replicate error)")  # Derived from Table 3 "Static and dynamic" column: sqrt(1.40^2 + 0.468^2)
  })

  model({
    kgrowth <- exp(lkgrowth)
    kdeath  <- exp(lkdeath)
    bmax    <- exp(lbmax)
    ke      <- exp(lke)
    emax    <- exp(lemax)
    ec50    <- exp(lec50)
    hill    <- exp(lhill)
    ke0     <- exp(lke0)
    cfu0    <- exp(lcfu0)

    btot <- bact_sensitive + bact_resting
    ksr <- (kgrowth - kdeath) / bmax * btot
    drug_kill <- emax * effect^hill / (ec50^hill + effect^hill)

    d/dt(central) <- -ke * central - kdeg * central
    d/dt(effect) <- ke0 * central - ke0 * effect
    d/dt(bact_sensitive) <- kgrowth * bact_sensitive - (kdeath + drug_kill) * bact_sensitive - ksr * bact_sensitive
    d/dt(bact_resting) <- ksr * bact_sensitive - kdeath * bact_resting

    f_resting_t0 <- (1 - fmix1) * fpers
    bact_sensitive(0) <- cfu0 * (1 - f_resting_t0)
    bact_resting(0)   <- cfu0 * f_resting_t0

    cfu_obs_floor <- btot + 1
    Cc <- log(cfu_obs_floor)
    Cc ~ add(addSd)
  })
}
