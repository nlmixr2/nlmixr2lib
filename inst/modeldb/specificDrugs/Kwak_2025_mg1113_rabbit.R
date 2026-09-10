Kwak_2025_mg1113_rabbit <- function() {
  description <- "Preclinical (rabbit). QSP (full TMDD). Allometrically scaled rabbit projection of the refined two-target target-mediated drug disposition model for MG1113, a humanized anti-tissue-factor-pathway-inhibitor (anti-TFPI) IgG4 antibody. Structure is identical to the cynomolgus-monkey fit (explicit bimolecular binding to both soluble TFPI-alpha and membrane-bound TFPI, target turnover, complex elimination, and a transit compartment for delayed subcutaneous absorption); every parameter is a PREDICTION obtained by allometric scaling of the monkey estimates to a 2.5 kg rabbit, not a rabbit fit. Externally validated against observed rabbit MG1113 profiles (Kwak 2025 Figure 3, AAFE 1.5-1.9 in the high-dose groups)."
  reference   <- "Kwak H, Jeong YS, Kim J, Lee M, Byoun S, Aoki Y, Chung SJ, Lee W. Refined target-mediated drug disposition modeling of the anti-tissue factor pathway inhibitor antibody MG1113 in cynomolgus monkeys and rabbits. Front Pharmacol. 2025;16:1745702. doi:10.3389/fphar.2025.1745702. PMCID PMC12819659. Model equations from the Supplementary Material section 4 (Supplementary Methods, Model equations); rabbit parameter values from Table 3 (Rabbit column). Allometric exponents (0.75 for clearance, -0.25 for rate constants, 1.0 for volumes) from Germovsek E et al. MAbs. 2021;13(1):1964935. Rabbit sTFPI-alpha baseline and KD from Kwak H et al. Res Pract Thromb Haemost. 2020;4(8):1301-1312 (doi:10.1002/rth2.12438). Monkey parameters this projection scales from: see modellib('Kwak_2025_mg1113_monkey')."
  vignette    <- "Kwak_2025_mg1113"

  paper_specific_compartments <- c("stfpi", "mtfpi", "astfpi", "amtfpi")

  units       <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482: the four TFPI-species states are carried by the paper as
  # CONCENTRATIONS (nM) in the central compartment, not as amounts.
  compartmentData <- list(
    depot       = list(analyte = "MG1113",                                units = "nmol", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "MG1113",                                units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "free MG1113",                           units = "nmol", specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "free MG1113",                           units = "nmol", specimen = "plasma",              verified = TRUE),
    stfpi       = list(analyte = "free soluble TFPI-alpha",               units = "nM",   specimen = "plasma",              verified = TRUE),
    mtfpi       = list(analyte = "free membrane-bound TFPI",              units = "nM",   specimen = "plasma",              verified = TRUE),
    astfpi      = list(analyte = "MG1113 / soluble TFPI-alpha complex",   units = "nM",   specimen = "plasma",              verified = TRUE),
    amtfpi      = list(analyte = "MG1113 / membrane-bound TFPI complex",  units = "nM",   specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "rabbit",
    n_subjects     = NA,
    n_studies      = 1,
    age_range      = "2-3 months",
    weight_range   = "2.5 kg representative body weight (Kwak 2025 section 2.6, citing Zhao 2015)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "Healthy male rabbits (Orient Bio); the 90-minute sTFPI-alpha comparison points come from a separate MG1113 efficacy study in a rabbit model of hemophilia A (Kwak 2020)",
    dose_range     = "Single i.v. at 2.75, 17.2 or 34.4 nmol/kg; single s.c. at 17.2, 34.4, 68.8 or 137.6 nmol/kg; n = 2-10 per group",
    regions        = "Republic of Korea",
    notes          = paste0(
      "PREDICTION, NOT A FIT. No parameter in this file was estimated from ",
      "rabbit data: every value is the monkey rank 1 estimate scaled ",
      "allometrically to a 2.5 kg rabbit (exponents 0.75 for CLD, -0.25 for ",
      "rate constants, 1.0 for volumes), except KD and kon, which were ",
      "assumed conserved across species on the strength of K2-domain sequence ",
      "similarity (rabbit vs human 92% identity), and the sTFPI-alpha ",
      "baseline, which is a literature value (1.114 nM, Kwak 2020). Every ",
      "parameter is therefore wrapped in fixed(). Predictive performance was ",
      "good in the high-dose groups (i.v. 34.4 nmol/kg; s.c. 68.8 and 137.6 ",
      "nmol/kg, AAFE 1.5-1.9) but the model OVER-PREDICTS the low-dose groups ",
      "(i.v. 2.75 and 17.2 nmol/kg; s.c. 17.2 and 34.4 nmol/kg, AAFE ",
      "2.8-4.2), not capturing the rapid decline of MG1113 there (Kwak 2025 ",
      "section 3.2). Typical-value simulator: no IIV and no residual-error ",
      "model are reported anywhere in the paper."
    ),
    model_class    = "QSP / full TMDD (allometrically scaled interspecies projection)",
    n_states       = 8
  )

  ini({
    # Every value below is an allometric projection or a conserved / literature
    # constant, NOT a rabbit estimate; hence fixed() throughout.

    # -- Binding parameters (assumed conserved across species; Kwak 2025 section 2.6) --
    kd     <- fixed(0.04665); label("Equilibrium dissociation constant of the MG1113/TFPI complex (nM)")                 # Kwak 2025 Table 3, Rabbit column; conserved from Kwak 2020 in-vitro measurement
    kon    <- fixed(28.50);   label("Association rate constant of MG1113 and TFPI (1/(nM*day))")                         # Kwak 2025 Table 3, Rabbit column; conserved across species, not scaled

    # -- Distribution (Kwak 2025 Table 3, Rabbit column) --
    lq     <- fixed(log(0.1130));   label("Inter-compartmental clearance of MG1113 (L/day)")                        # Kwak 2025 Table 3, Rabbit column; monkey 0.1454 scaled by (2.5/3.5)^0.75; paper symbol CLD
    lvc    <- fixed(log(0.07944));  label("Central compartment volume (L)")                                          # Kwak 2025 Table 3, Rabbit column; monkey 0.1112 scaled by (2.5/3.5)^1.0; paper symbol V2
    lvp    <- fixed(log(0.09615));  label("Peripheral compartment volume (L)")                                       # Kwak 2025 Table 3, Rabbit column; monkey 0.1346 scaled by (2.5/3.5)^1.0; paper symbol V3

    # -- Subcutaneous absorption (Kwak 2025 Table 3, Rabbit column) --
    k_depot_transit1   <- fixed(22.10);      label("Depot-to-transit transfer rate constant (1/day)")               # Kwak 2025 Table 3, Rabbit column; monkey 20.32 scaled by (2.5/3.5)^-0.25; paper symbol k01
    k_depot_central    <- fixed(9.637e-8);   label("Depot-to-central direct absorption rate constant (1/day)")      # Kwak 2025 Table 3, Rabbit column; scaled from a monkey value judged non-identifiable; paper symbol k02
    k_transit1_central <- fixed(0.4387);     label("Transit-to-central transfer rate constant (1/day)")             # Kwak 2025 Table 3, Rabbit column; monkey 0.4033 scaled by (2.5/3.5)^-0.25; paper symbol k12
    kloss              <- fixed(3.425e-7);   label("Loss rate constant of MG1113 at the depot compartment (1/day)")      # Kwak 2025 Table 3, Rabbit column; scaled from a monkey value judged non-identifiable

    # -- MG1113 elimination (Kwak 2025 Table 3, Rabbit column) --
    lkel   <- fixed(log(0.4942)); label("First-order elimination rate constant of MG1113 from the central compartment (1/day)")  # Kwak 2025 Table 3, Rabbit column; monkey 0.4543 scaled by (2.5/3.5)^-0.25

    # -- Target turnover and baselines (Kwak 2025 Table 3, Rabbit column) --
    kdegs   <- fixed(82.13);  label("Degradation rate constant of soluble TFPI-alpha (1/day)")                           # Kwak 2025 Table 3, Rabbit column; monkey 75.50 scaled by (2.5/3.5)^-0.25
    kdegm   <- fixed(1.251);  label("Degradation rate constant of membrane-bound TFPI (1/day)")                          # Kwak 2025 Table 3, Rabbit column; monkey 1.150 scaled by (2.5/3.5)^-0.25
    stfpi_b <- fixed(1.114);  label("Baseline concentration of soluble TFPI-alpha (nM)")                                 # Kwak 2025 Table 3, Rabbit column; literature value taken from Kwak 2020, not scaled
    mtfpi_b <- fixed(14.19);  label("Baseline concentration of membrane-bound TFPI (nM)")                                # Kwak 2025 Table 3, Rabbit column; from the monkey mTFPI:sTFPI-alpha ratio applied to the rabbit sTFPI-alpha baseline (section 2.6)

    # -- Drug-target complex elimination / internalisation (Kwak 2025 Table 3, Rabbit column) --
    kints  <- fixed(0.3366);   label("Elimination rate constant of the MG1113/soluble TFPI-alpha complex (1/day)")   # Kwak 2025 Table 3, Rabbit column; monkey 0.3094 scaled by (2.5/3.5)^-0.25; paper symbol kel,CS
    kintm  <- fixed(0.007452); label("Elimination rate constant of the MG1113/membrane-bound TFPI complex (1/day)")  # Kwak 2025 Table 3, Rabbit column; monkey 0.006851 scaled by (2.5/3.5)^-0.25; paper symbol kel,CM
  })

  model({
    # 1. Back-transform log-scale parameters
    q   <- exp(lq)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    kel <- exp(lkel)

    # 2. Derived rate constants (Kwak 2025 Supplementary Material section 4,
    #    'Related equations'). koff is pinned by KD; the two synthesis rates
    #    are pinned by the steady-state balance at the reported baselines.
    koff  <- kd * kon
    ksyns <- kdegs * stfpi_b
    ksynm <- kdegm * mtfpi_b

    # 3. Free MG1113 concentrations (nM). See the monkey model file for why
    #    the two drug states are carried as amounts rather than as the
    #    paper's concentrations; the ODEs below are the paper's equations
    #    multiplied through by their compartment volumes.
    cfree <- central / vc
    cper  <- peripheral1 / vp

    # 4. ODE system
    d/dt(depot)       <- -(k_depot_transit1 + k_depot_central + kloss) * depot
    d/dt(transit1)    <- k_depot_transit1 * depot - k_transit1_central * transit1

    d/dt(central)     <- k_depot_central * depot +
                         k_transit1_central * transit1 +
                         q * (cper - cfree) +
                         vc * koff * (amtfpi + astfpi) -
                         vc * kon * cfree * (mtfpi + stfpi) -
                         kel * central

    d/dt(peripheral1) <- q * (cfree - cper)

    d/dt(mtfpi)       <- ksynm - kdegm * mtfpi + koff * amtfpi - kon * mtfpi * cfree
    d/dt(stfpi)       <- ksyns - kdegs * stfpi + koff * astfpi - kon * stfpi * cfree

    d/dt(amtfpi)      <- kon * mtfpi * cfree - (koff + kintm) * amtfpi
    d/dt(astfpi)      <- kon * stfpi * cfree - (koff + kints) * astfpi

    # 5. Endogenous initial conditions (Supplementary Material section 4)
    mtfpi(0) <- mtfpi_b
    stfpi(0) <- stfpi_b

    # 6. Outputs (nM). Both assays detect the FREE form (Kwak 2025 section
    #    2.1). The paper reports no residual-error model, so none is applied.
    Cc    <- cfree
    sTFPI <- stfpi
  })
}
