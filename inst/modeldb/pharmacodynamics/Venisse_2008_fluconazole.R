Venisse_2008_fluconazole <- function() {
  description <- "In vitro (Candida albicans, ATCC 3153). Mechanism-based PK-PD model of fluconazole fungistatic activity in a dynamic broth-renewal flask. Candida population dynamics follow Mouton-type saturable exponential growth with a maximum carrying capacity Nmax minus a first-order natural-death rate Kd minus a fixed broth-renewal rate Ke (set by the experimenters to 0.231 1/h, equal to the drug 1-compartment kel and yielding a 3 h elimination half-life). Fluconazole multiplicatively inhibits growth via a sigmoid-Emax (Hill=1) Imax * C / (IC50 + C) term where Imax < 1 represents the maximum fractional reduction of the growth rate. Drug concentration is the central/Vc concentration from a one-compartment PK with V and CL fixed to the apparatus geometry (V = 0.4 L bulk volume of the glass flask, CL = 1.54 mL/min consistent with t1/2 = 3 h). Candida growth/death/Nmax parameters were estimated from all six experiments (3 fluconazole + 3 caspofungin) in a single Nonmem run, so the Kg / Kd / Nmax values reproduced here are the joint typical-value estimates shared with the companion model Venisse_2008_caspofungin. Interindividual variability was reported on Nmax (225% CV exponential); the residual error was reported as exponential at 197% CV on the linear CFU/mL scale, encoded here as an additive SD on natural-log CFU/mL via sigma_ln = sqrt(log(CV^2 + 1))."
  reference <- paste(
    "Venisse N, Gregoire N, Marliat M, Couet W. (2008).",
    "Mechanism-based pharmacokinetic-pharmacodynamic models of in vitro",
    "fungistatic and fungicidal effects against Candida albicans.",
    "Antimicrobial Agents and Chemotherapy 52(3):937-943.",
    "doi:10.1128/AAC.01030-07."
  )
  vignette <- "Venisse_2008_candida_albicans"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L (drug central; numerically equal to ug/mL used in the paper); CFU/mL (Candida count); log CFU/mL (Cc observation)"
  )

  paper_specific_compartments <- c("candida")

  covariateData <- list()

  population <- list(
    species                = "in vitro (Candida albicans, ATCC 3153 reference strain)",
    n_subjects             = NA_integer_,
    n_studies              = 1L,
    organism               = "Candida albicans ATCC 3153 (NCCLS-susceptible reference strain)",
    system                 = "Dynamic in vitro infection model: one-compartment glass flask (400 mL working volume) with a peristaltic pump (Masterflex L/S) continuously supplying and removing RPMI 1640 + 2% sorbitol at a flow rate adjusted to simulate a drug elimination half-life of 3 h",
    medium                 = "RPMI 1640 + 2% sorbitol",
    temperature            = "37 C",
    duration               = "31 h post-inoculation (4 h pre-drug control growth + 27 h post-drug observation)",
    mic                    = "fluconazole MIC = 1.56 ug/mL on this isolate (NCCLS-susceptible)",
    concentration_range    = "fluconazole initial concentrations 0.01-100 ug/mL across three experiments (Table 1)",
    regimens               = "Three fluconazole-monotherapy experiments + three caspofungin-monotherapy experiments (six experiments total); each experiment includes a drug-free growth control and three drug-treated arms differing by initial concentration",
    notes                  = "Initial Candida inoculum 5e3 CFU/mL introduced at t = 0 h; fluconazole introduced 4 h later as a bolus to the bath. Sampling at 4 h (before drug) and approximately 4, 6, 8, 15, 18, 22, 27 h after drug introduction. The joint Nonmem run pooled all 6 experiments to estimate Candida-specific parameters (Kg, Kd, Nmax) shared across drugs and drug-specific parameters (IC50, Imax for fluconazole) per experiment set."
  )

  ini({
    # ----- Drug PK (one-compartment, both V and CL fixed by apparatus geometry) -----
    # Methods (page 938): V = 0.4 L (bulk volume of the glass compartment),
    # CL = 1.54 mL/min = 0.0924 L/h, half-life = 3 h. Both fixed.
    lcl <- fixed(log(0.0924))
    label("Drug clearance (L/h; CL = 1.54 mL/min, fixed by broth-pump flow rate)")
    # Venisse 2008 Methods (PK-PD analysis paragraph): CL = 1.54 mL/min fixed.
    lvc <- fixed(log(0.4))
    label("Drug central volume (L; V = 0.4 L bulk bath volume, fixed)")
    # Venisse 2008 Methods (PK-PD analysis paragraph): V = 0.4 liters fixed.

    # ----- Candida population dynamics (shared across all 6 experiments) -----
    lkg <- log(0.864)
    label("Candida natural growth rate constant Kg (1/h)")
    # Venisse 2008 Table 2: Kg = 0.864 1/h, SE = 0.0321.
    lkd <- log(0.0414)
    label("Candida natural-death rate constant Kd (1/h)")
    # Venisse 2008 Table 2: Kd = 0.0414 1/h, SE = 0.0245.
    lnmax <- log(1.48e6)
    label("Maximum Candida carrying capacity Nmax (CFU/mL)")
    # Venisse 2008 Table 2: Nmax = 1.48e6 CFU/mL, SE = 5.54e5.
    lke <- fixed(log(0.231))
    label("Broth-renewal elimination rate constant Ke (1/h; fixed at 0.231 = ln(2)/3 h)")
    # Venisse 2008 Methods (page 938): Ke fixed at 0.231 1/h, equal to the drug
    # 1-compartment kel = CL/V = 0.0924/0.4. Set by experimenters (pump flow rate).

    # ----- Fluconazole growth-inhibition PD parameters -----
    lic50 <- log(0.0929)
    label("Fluconazole IC50 (ug/mL); concentration for 50% of maximal growth inhibition")
    # Venisse 2008 Table 2: IC50 = 0.0929 ug/mL = 9.29e-2, SE = 5.67e-2.
    imax <- 0.675
    label("Fluconazole Imax (unitless, bounded 0-1); maximum fractional growth inhibition")
    # Venisse 2008 Table 2: Imax = 0.675, SE = 0.0342. Confined in [0, 1] by the
    # paper's parameterisation; Imax < 1 reflects the observation that even at
    # very high fluconazole concentrations Candida growth is not fully inhibited.

    # ----- Interindividual variability (between-experiment) -----
    # Paper Table 2 reports eta(Nmax) = 225% CV exponential model.
    # omega^2 = log(CV^2 + 1) = log(2.25^2 + 1) = log(6.0625) = 1.802.
    etalnmax ~ 1.802
    # Venisse 2008 Table 2: eta(Nmax) = 225% CV; converted via log(CV^2 + 1).

    # ----- Residual variability -----
    # Paper Table 2: epsilon_1 fluconazole = 197% CV (exponential residual on
    # CFU/mL). Mapped to additive SD on natural-log CFU/mL:
    # sigma_ln = sqrt(log(CV^2 + 1)) = sqrt(log(1.97^2 + 1)) = sqrt(1.5854) = 1.259.
    addSd <- 1.259
    label("Additive residual SD on natural-log CFU/mL (197% CV exponential -> log-scale SD)")
    # Venisse 2008 Table 2: epsilon_1 fluconazole = 197% CV exponential.
  })

  model({
    # ----- Back-transform structural parameters -----
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    kg   <- exp(lkg)
    kd   <- exp(lkd)
    nmax <- exp(lnmax + etalnmax)
    ke   <- exp(lke)
    ic50 <- exp(lic50)

    # ----- Drug PK micro-constant -----
    kel <- cl / vc

    # ----- Bath drug concentration (drives the fluconazole inhibition term) -----
    cc <- central / vc

    # ----- Fluconazole growth-inhibition fractional effect (Imax * C / (IC50 + C)) -----
    e_flu <- imax * cc / (ic50 + cc)

    # ----- ODE system -----
    # Eq 4 (Venisse 2008):
    #   dN/dt = Kg * (1 - N/Nmax) * (1 - Imax*C/(IC50+C)) * N - Kd*N - Ke*N
    # with drug PK in a separate 1-compartment central state (Methods PK-PD
    # analysis paragraph and Eq 3 / Eq 4).
    d/dt(central) <- -kel * central
    d/dt(candida) <-
      kg * (1 - candida / nmax) * (1 - e_flu) * candida -
      kd * candida -
      ke * candida

    # Initial Candida inoculum (Methods page 938: "starting concentration of
    # 5 x 10^3 organisms per ml ... starting inoculum ... was 5 x 10^3 cells
    # per ml"). Drug enters via a dosing event at the experimental
    # drug-introduction time (4 h post inoculation).
    candida(0) <- 5000

    # ----- Observation: log(Candida) with additive residual on the log scale -----
    # Equivalent to the paper's exponential residual on the linear CFU/mL scale;
    # the +1 floor prevents log(0) when drug effects drive the count toward 0
    # (does not arise for fluconazole, which is fungistatic).
    Cc <- log(candida + 1)
    Cc ~ add(addSd)
  })
}
