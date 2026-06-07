Nielsen_2011_cefuroxime <- function() {
  description <- "In vitro (Streptococcus pyogenes M12 NCTC P1800). Semimechanistic PKPD model of cefuroxime time-kill kinetics; two-stage bacterial life-cycle (proliferating drug-sensitive S and non-growing drug-insensitive R) with sigmoidal Emax killing of S via an effect compartment; first-order drug elimination (ke set per in vitro kinetic-system flow rate) plus drug-specific degradation kdeg. Parameter values are from the combined static and dynamic estimation in Table 3."
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
    organism      = "Streptococcus pyogenes group A M12 NCTC P1800 (National Culture Type Collection); cefuroxime MIC 0.0313 mg/L",
    medium        = "Todd-Hewitt broth (35 C, 5% CO2 for plating)",
    inoculum      = "Target 10^6 CFU/mL at t = 0 (logarithmic-growth-phase culture)",
    system        = "Static 10-mL test tubes (4 mL broth) and a dynamic in vitro kinetic system (110-mL open-bottom spinner flask with pump-driven dilution producing first-order antibiotic elimination at flow-rate / volume)",
    duration      = "24 h (some experiments 48 h)",
    drug          = "Cefuroxime (Zinacef; GlaxoSmithKline)",
    dose_range    = "Static initial concentrations 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 16, 64 x MIC; dynamic initial concentrations 2 and 16 x MIC with simulated half-lives 0 (constant infusion), 1.7 (human), and 5.1 (one-third human) h",
    n_observations = 588L,
    notes         = "Parameters are typical values from the simultaneous static + dynamic estimation (Table 3 'Static and dynamic' column). Random effects (eta) are NOT included: the source NONMEM run estimated only inter-experiment variability in the resting fraction at t=0 via a two-component mixture model, plus replicate and common additive residual error on the natural-log viable count. The packaged file folds the mixture into a typical starting fraction (E[f_resting] = (1 - fmix1) * fpers) and combines the replicate and common residual SDs into a single additive SD."
  )

  ini({
    # =============================================================
    # Bacterial life-cycle parameters (shared across the five drugs)
    # =============================================================
    # Values are from Table 3 column "Static and dynamic" (the
    # combined-data joint estimation across all five antibiotics).
    lkgrowth <- log(1.46)
    label("Log bacterial multiplication rate kgrowth (1/h)")          # Table 3 "Static and dynamic" column: kgrowth = 1.46 (RSE 4.1%)
    lkdeath  <- log(0.187)
    label("Log natural-death rate kdeath (1/h)")                      # Table 3 "Static and dynamic" column: kdeath = 0.187 (RSE 6.2%)
    lbmax    <- log(5.00e8)
    label("Log carrying-capacity Bmax (CFU/mL)")                      # Table 3 "Static and dynamic" column: Bmax = 5.00e8 (RSE 8.7%)

    # Mixture-model parameters: NONMEM allowed two subpopulations of
    # experiments differing in the resting-fraction at t = 0. Subpop
    # 1 has no resting bacteria; subpop 2 has fraction fpers in the
    # resting stage. fmix1 is the proportion of experiments in
    # subpop 1. The packaged file uses the mixture-mean fraction as
    # a single deterministic typical-value initial condition.
    fpers <- 0.0652
    label("Fraction of bacteria in the resting stage at t = 0 (subpopulation 2; unitless)")  # Table 3 "Static and dynamic" column: fpers = 0.0652 (RSE 49%)
    fmix1 <- 0.880
    label("Proportion of experiments in subpopulation 1 (no resting bacteria at t = 0; unitless)")  # Table 3 "Static and dynamic" column: fMix1 = 0.880 (RSE 7.5%)

    # =============================================================
    # Drug-specific PK parameters (cefuroxime)
    # =============================================================
    # kdeg is the in-broth degradation rate of the drug; Nielsen
    # 2011 fixed it for benzylpenicillin and cefuroxime based on
    # 24-h incubation experiments (reference 31). For the other
    # three drugs (erythromycin, moxifloxacin, vancomycin) kdeg
    # was set to zero. Encoded on the linear scale so the same
    # ini() shape works across the five Nielsen 2011 drug files.
    kdeg <- fixed(0.026)
    label("Drug degradation rate kdeg in broth (1/h; FIXED per Methods 'Semimechanistic PKPD model')")  # Methods: "for these drugs the kdeg values were fixed to 0.020 and 0.026 h-1, respectively" (cefuroxime = 0.026)

    # ke is the in vitro kinetic-system elimination rate (flow rate
    # / volume). For dynamic experiments it is set to ln(2) /
    # simulated-half-life; for static experiments it is 0. The
    # default value below corresponds to a 1.7-h simulated human
    # half-life (Table 1 dynamic-cefuroxime column); override via
    # params = list(lke = ...) for a different half-life or to ~ 0
    # (e.g., log(1e-9)) for a static-concentration simulation.
    lke <- fixed(log(log(2) / 1.7))
    label("Log in vitro kinetic-system elimination rate ke (1/h; default ln(2)/1.7 = human half-life; override per experiment)")  # Table 1: simulated human half-life for cefuroxime = 1.7 h

    # =============================================================
    # Drug-specific PD parameters (cefuroxime; Table 3)
    # =============================================================
    lemax <- log(2.72)
    label("Log maximum killing rate Emax (1/h)")                                          # Table 3 "Static and dynamic" column: Emax CXM = 2.72 (RSE 5.1%)
    lec50 <- log(0.00787)
    label("Log effect-compartment concentration for 50% of Emax (mg/L)")                  # Table 3 "Static and dynamic" column: EC50 CXM = 0.00787 (RSE 8.3%)
    lhill <- log(1.35)
    label("Log Hill (sigmoidicity) exponent gamma on the killing function (unitless)")    # Table 3 "Static and dynamic" column: gamma CXM = 1.35 (RSE 14%)
    lke0  <- fixed(log(100))
    label("Log effect-compartment delay rate ke0 (1/h; FIXED at 100 per Table 3 'Static and dynamic' column)")  # Table 3 "Static and dynamic" column: ke0 CXM = 100 FIX (effectively instantaneous equilibrium)

    # =============================================================
    # Starting inoculum
    # =============================================================
    lcfu0 <- fixed(log(1e6))
    label("Log starting bacterial inoculum (CFU/mL; FIXED at 10^6 per Methods)")  # Methods: "starting inoculum of 10^6 CFU/ml"

    # =============================================================
    # Residual error
    # =============================================================
    # See Nielsen_2011_benzylpenicillin.R for the derivation: the
    # paper reports common (eps = 1.40) and replicate (epsrepl =
    # 0.468) additive SDs on natural log; combined in quadrature.
    addSd <- 1.476
    label("Additive residual SD on natural-log CFU/mL (combined common + replicate error)")  # Derived from Table 3 "Static and dynamic" column: sqrt(1.40^2 + 0.468^2)
  })

  model({
    # =============================================================
    # Effective (linear-scale) parameters
    # =============================================================
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
