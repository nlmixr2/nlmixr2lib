Kwak_2025_mg1113_human <- function() {
  description <- "QSP (full TMDD). Allometrically scaled human projection of the refined two-target target-mediated drug disposition model for MG1113, a humanized anti-tissue-factor-pathway-inhibitor (anti-TFPI) IgG4 antibody in development for hemophilia. Structure is identical to the cynomolgus-monkey fit (explicit bimolecular binding to both soluble TFPI-alpha and membrane-bound TFPI, target turnover, complex elimination, and a transit compartment for delayed subcutaneous absorption); every parameter is a PREDICTION obtained by allometric scaling of the monkey estimates to a 70 kg adult, not a fit to human data. Used to predict that weekly s.c. dosing at 3.3 mg/kg holds sTFPI-alpha below 25% of baseline for most of the dosing interval (Kwak 2025 Figure 4)."
  reference   <- "Kwak H, Jeong YS, Kim J, Lee M, Byoun S, Aoki Y, Chung SJ, Lee W. Refined target-mediated drug disposition modeling of the anti-tissue factor pathway inhibitor antibody MG1113 in cynomolgus monkeys and rabbits. Front Pharmacol. 2025;16:1745702. doi:10.3389/fphar.2025.1745702. PMCID PMC12819659. Model equations from the Supplementary Material section 4 (Supplementary Methods, Model equations); human parameter values from Table 3 (Human column). Allometric exponents (0.75 for clearance, -0.25 for rate constants, 1.0 for volumes) from Germovsek E et al. MAbs. 2021;13(1):1964935. Human sTFPI-alpha baseline from Kwak EY et al. J Thromb Haemost. 2021;19(6):1425-1435 (doi:10.1111/jth.15244); KD from Kwak H et al. Res Pract Thromb Haemost. 2020;4(8):1301-1312. Monkey parameters this projection scales from: see modellib('Kwak_2025_mg1113_monkey')."
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

  # Screened / discussed but NOT retained as a covariate in the model. The
  # paper varies the sTFPI-alpha baseline over its reported 95% prediction
  # interval in a LOCAL SENSITIVITY ANALYSIS (section 3.3, Figure 5) rather
  # than fitting a covariate relationship, so there is no coefficient to
  # encode. To reproduce that analysis, set stfpi_b directly and rescale
  # mtfpi_b as mtfpi_b = 12.7 * stfpi_b (Supplementary Material item 7).
  covariatesDataExcluded <- list(
    STFPI_BASE = list(
      description = "Individual baseline plasma soluble TFPI-alpha concentration",
      units       = "nM",
      type        = "continuous",
      notes       = paste0(
        "95% prediction interval 1.3-2.9 nM in human plasma (Dahm 2003), ",
        "around the 2.3 nM reference used here. Explored by local sensitivity ",
        "analysis in Kwak 2025 section 3.3 / Figure 5, not by an estimated ",
        "covariate effect. Higher baselines reduce MG1113 AUC0-30days (more so ",
        "after s.c. than i.v. dosing) and shorten the time sTFPI-alpha stays ",
        "below 25% of baseline; at 1.3 nM, 3.3 mg/kg suppresses sTFPI-alpha ",
        "below 25% for about 6-7 days by either route. Baseline sTFPI-alpha is ",
        "reported not to differ between healthy adults and patients with ",
        "hemophilia (Gu 2015)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA,
    n_studies      = 0,
    age_range      = "adult",
    weight_range   = "70 kg reference body weight (Kwak 2025 section 2.6, citing Zhao 2015)",
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state  = "Healthy adults and patients with hemophilia (simulated phase I dosing regimens; no human MG1113 data were available to the authors)",
    dose_range     = "Simulated single and weekly (Q7d) s.c. doses of 0.5, 1.7 and 3.3 mg/kg and i.v. doses of 3.3 mg/kg (3.44, 11.7 and 22.7 nmol/kg), matching the phase I trial NCT03855696",
    regions        = NA,
    notes          = paste0(
      "PREDICTION, NOT A FIT, AND NOT VALIDATED. No human MG1113 data were ",
      "available, so nothing in this file was estimated from or checked ",
      "against human observations (Kwak 2025 section 4: 'Due to the ",
      "unavailability of MG1113 clinical trial data, we could not validate the ",
      "simulated profiles of MG1113 and sTFPI-alpha in humans'). Every value ",
      "is the monkey rank 1 estimate scaled allometrically to 70 kg ",
      "(exponents 0.75 for CLD, -0.25 for rate constants, 1.0 for volumes), ",
      "except KD and kon, which were assumed conserved across species on the ",
      "strength of K2-domain sequence similarity (human vs monkey 96% ",
      "identity), and the sTFPI-alpha baseline, which is a literature value ",
      "(2.3 nM). Every parameter is therefore wrapped in fixed(). Typical-value ",
      "simulator: no IIV and no residual-error model are reported anywhere in ",
      "the paper. Clinical context: sTFPI-alpha below 25% of baseline is the ",
      "efficacy-associated target (Chowdary 2015, Eichler 2018)."
    ),
    model_class    = "QSP / full TMDD (allometrically scaled interspecies projection)",
    n_states       = 8
  )

  ini({
    # Every value below is an allometric projection or a conserved / literature
    # constant, NOT a human estimate; hence fixed() throughout.

    # -- Binding parameters (assumed conserved across species; Kwak 2025 section 2.6) --
    kd     <- fixed(0.04665); label("Equilibrium dissociation constant of the MG1113/TFPI complex (nM)")                 # Kwak 2025 Table 3, Human column; conserved from Kwak 2020 in-vitro measurement
    kon    <- fixed(28.50);   label("Association rate constant of MG1113 and TFPI (1/(nM*day))")                         # Kwak 2025 Table 3, Human column; conserved across species, not scaled

    # -- Distribution (Kwak 2025 Table 3, Human column) --
    lq     <- fixed(log(1.375));   label("Inter-compartmental clearance of MG1113 (L/day)")                         # Kwak 2025 Table 3, Human column; monkey 0.1454 scaled by (70/3.5)^0.75; paper symbol CLD
    lvc    <- fixed(log(2.224));   label("Central compartment volume (L)")                                           # Kwak 2025 Table 3, Human column; monkey 0.1112 scaled by (70/3.5)^1.0; paper symbol V2
    lvp    <- fixed(log(2.692));   label("Peripheral compartment volume (L)")                                        # Kwak 2025 Table 3, Human column; monkey 0.1346 scaled by (70/3.5)^1.0; paper symbol V3

    # -- Subcutaneous absorption (Kwak 2025 Table 3, Human column) --
    k_depot_transit1   <- fixed(9.609);      label("Depot-to-transit transfer rate constant (1/day)")               # Kwak 2025 Table 3, Human column; monkey 20.32 scaled by (70/3.5)^-0.25; paper symbol k01
    k_depot_central    <- fixed(4.189e-8);   label("Depot-to-central direct absorption rate constant (1/day)")      # Kwak 2025 Table 3, Human column; scaled from a monkey value judged non-identifiable; paper symbol k02
    k_transit1_central <- fixed(0.1907);     label("Transit-to-central transfer rate constant (1/day)")             # Kwak 2025 Table 3, Human column; monkey 0.4033 scaled by (70/3.5)^-0.25; paper symbol k12
    kloss              <- fixed(1.489e-7);   label("Loss rate constant of MG1113 at the depot compartment (1/day)")      # Kwak 2025 Table 3, Human column; scaled from a monkey value judged non-identifiable

    # -- MG1113 elimination (Kwak 2025 Table 3, Human column) --
    lkel   <- fixed(log(0.2148)); label("First-order elimination rate constant of MG1113 from the central compartment (1/day)")  # Kwak 2025 Table 3, Human column; monkey 0.4543 scaled by (70/3.5)^-0.25

    # -- Target turnover and baselines (Kwak 2025 Table 3, Human column) --
    kdegs   <- fixed(35.70);  label("Degradation rate constant of soluble TFPI-alpha (1/day)")                           # Kwak 2025 Table 3, Human column; monkey 75.50 scaled by (70/3.5)^-0.25
    kdegm   <- fixed(0.5440); label("Degradation rate constant of membrane-bound TFPI (1/day)")                          # Kwak 2025 Table 3, Human column; monkey 1.150 scaled by (70/3.5)^-0.25
    stfpi_b <- fixed(2.3);    label("Baseline concentration of soluble TFPI-alpha (nM)")                                 # Kwak 2025 Table 3, Human column; literature value taken from Kwak 2021, not scaled
    mtfpi_b <- fixed(29.11);  label("Baseline concentration of membrane-bound TFPI (nM)")                                # Kwak 2025 Table 3, Human column; from the monkey mTFPI:sTFPI-alpha ratio applied to the human sTFPI-alpha baseline (section 2.6)

    # -- Drug-target complex elimination / internalisation (Kwak 2025 Table 3, Human column) --
    kints  <- fixed(0.1463);   label("Elimination rate constant of the MG1113/soluble TFPI-alpha complex (1/day)")   # Kwak 2025 Table 3, Human column; monkey 0.3094 scaled by (70/3.5)^-0.25; paper symbol kel,CS
    kintm  <- fixed(0.003240); label("Elimination rate constant of the MG1113/membrane-bound TFPI complex (1/day)")  # Kwak 2025 Table 3, Human column; monkey 0.006851 scaled by (70/3.5)^-0.25; paper symbol kel,CM
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
