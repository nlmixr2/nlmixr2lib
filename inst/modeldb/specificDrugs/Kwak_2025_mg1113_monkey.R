Kwak_2025_mg1113_monkey <- function() {
  description <- "Preclinical (cynomolgus monkey). QSP (full TMDD). Refined two-target target-mediated drug disposition model for MG1113, a humanized anti-tissue-factor-pathway-inhibitor (anti-TFPI) IgG4 antibody, describing explicit bimolecular binding to BOTH soluble TFPI-alpha (sTFPI-alpha) and membrane-bound TFPI (mTFPI) in a two-compartment PK framework, with zero-order synthesis / first-order degradation turnover of each target, first-order elimination of each drug-target complex, and a single transit compartment for delayed subcutaneous absorption. Parameters fitted to monkey MG1113 and sTFPI-alpha plasma profiles by the Cluster Gauss-Newton Method (rank 1 accepted parameter set, Kwak 2025 Table 1)."
  reference   <- "Kwak H, Jeong YS, Kim J, Lee M, Byoun S, Aoki Y, Chung SJ, Lee W. Refined target-mediated drug disposition modeling of the anti-tissue factor pathway inhibitor antibody MG1113 in cynomolgus monkeys and rabbits. Front Pharmacol. 2025;16:1745702. doi:10.3389/fphar.2025.1745702. PMCID PMC12819659. Model equations from the Supplementary Material section 4 (Supplementary Methods, Model equations); parameter values from Table 1 (rank 1 column). KD fixed from Kwak H et al. Res Pract Thromb Haemost. 2020;4(8):1301-1312 (doi:10.1002/rth2.12438). Monkey MG1113 and sTFPI-alpha data reanalysed from Kwak EY et al. J Thromb Haemost. 2021;19(6):1425-1435 (doi:10.1111/jth.15244)."
  vignette    <- "Kwak_2025_mg1113"

  paper_specific_compartments <- c("stfpi", "mtfpi", "astfpi", "amtfpi")

  units       <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The four TFPI-species states are carried by the paper
  # as CONCENTRATIONS (nM) in the central compartment, not as amounts; the
  # `units` field below records that, matching Kwak 2025 Supplementary
  # Material section 4, which labels each of them "The concentration (nM) of
  # ... in the central compartment".
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
    species        = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects     = 18,
    n_studies      = 1,
    age_range      = NA,
    weight_range   = "3.5 kg representative body weight; individual weights were not available (Kwak 2025 section 2.1, citing Zhao 2015)",
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state  = "Healthy cynomolgus monkeys (non-disease preclinical PK/PD study)",
    dose_range     = "Single i.v. or s.c. MG1113 at 2.5, 5.0 or 10.0 mg/kg (17.2, 34.4 or 68.8 nmol/kg), i.e. 60.2, 120.4 or 240.8 nmol per 3.5 kg monkey; n = 1-3 per dose/route group",
    regions        = NA,
    notes          = paste0(
      "Typical-value mechanistic (QSP / full TMDD) simulator. The Cluster ",
      "Gauss-Newton Method (CGNM) was fitted to the MEAN observed profiles, so ",
      "the paper reports NO inter-individual variability and NO residual-error ",
      "model; none is encoded here. Of 1,000 CGNM parameter sets, 136 met the ",
      "acceptance criterion SSR <= 6.34; the rank 1 set is encoded. Five ",
      "parameters (k02, kel_CM, kel_CS, kloss and kon) were judged ",
      "NON-IDENTIFIABLE by approximate profile likelihood (Kwak 2025 ",
      "Supplementary Table S1); see the vignette for the profile-likelihood ",
      "intervals. Plasma MG1113 and sTFPI-alpha were measured by ELISAs that ",
      "detect the FREE forms, so the model outputs are free species. Data ",
      "reanalysed from Kwak EY 2021 (J Thromb Haemost 19:1425-1435)."
    ),
    model_class    = "QSP / full TMDD (explicit bimolecular binding to two target pools, soluble and membrane-bound, with target turnover and complex internalisation)",
    n_states       = 8
  )

  ini({
    # -- Binding parameters (conserved across species; Kwak 2025 section 2.6) --
    kd     <- fixed(0.04665); label("Equilibrium dissociation constant of the MG1113/TFPI complex (nM)")                 # Kwak 2025 Table 1, 'Fixed parameter' block; adopted from Kwak 2020 in-vitro measurement
    kon    <- 28.50;          label("Association rate constant of MG1113 and TFPI (1/(nM*day))")                         # Kwak 2025 Table 1, rank 1; koff is derived as kd * kon

    # -- Distribution (Kwak 2025 Table 1, rank 1) --
    lq     <- log(0.1454);    label("Inter-compartmental clearance of MG1113 (L/day)")                              # Kwak 2025 Table 1, rank 1; paper symbol CLD
    lvc    <- log(0.1112);    label("Central compartment volume (L)")                                                # Kwak 2025 Table 1, rank 1; paper symbol V2
    lvp    <- log(0.1346);    label("Peripheral compartment volume (L)")                                             # Kwak 2025 Table 1, rank 1; paper symbol V3

    # -- Subcutaneous absorption (Kwak 2025 Table 1, rank 1) --
    k_depot_transit1   <- 20.32;      label("Depot-to-transit transfer rate constant (1/day)")                      # Kwak 2025 Table 1, rank 1; paper symbol k01
    k_depot_central    <- 8.859e-8;   label("Depot-to-central direct absorption rate constant (1/day)")             # Kwak 2025 Table 1, rank 1; non-identifiable, minute value retained per section 3.1; paper symbol k02
    k_transit1_central <- 0.4033;     label("Transit-to-central transfer rate constant (1/day)")                    # Kwak 2025 Table 1, rank 1; paper symbol k12
    kloss              <- 3.149e-7;   label("Loss rate constant of MG1113 at the depot compartment (1/day)")             # Kwak 2025 Table 1, rank 1; non-identifiable, minute value retained per section 3.1

    # -- MG1113 elimination (Kwak 2025 Table 1, rank 1) --
    lkel   <- log(0.4543);    label("First-order elimination rate constant of MG1113 from the central compartment (1/day)")  # Kwak 2025 Table 1, rank 1

    # -- Target turnover and baselines (Kwak 2025 Table 1, rank 1) --
    kdegs   <- 75.50;         label("Degradation rate constant of soluble TFPI-alpha (1/day)")                           # Kwak 2025 Table 1, rank 1
    kdegm   <- 1.150;         label("Degradation rate constant of membrane-bound TFPI (1/day)")                          # Kwak 2025 Table 1, rank 1
    stfpi_b <- 0.9456;        label("Baseline concentration of soluble TFPI-alpha (nM)")                                 # Kwak 2025 Table 1, rank 1; experimentally measured value was 0.9768 nM
    mtfpi_b <- 12.04;         label("Baseline concentration of membrane-bound TFPI (nM)")                                # Kwak 2025 Table 1, rank 1

    # -- Drug-target complex elimination / internalisation (Kwak 2025 Table 1, rank 1) --
    kints  <- 0.3094;         label("Elimination rate constant of the MG1113/soluble TFPI-alpha complex (1/day)")   # Kwak 2025 Table 1, rank 1; paper symbol kel,CS
    kintm  <- 0.006851;       label("Elimination rate constant of the MG1113/membrane-bound TFPI complex (1/day)")  # Kwak 2025 Table 1, rank 1; paper symbol kel,CM
  })

  model({
    # 1. Back-transform log-scale parameters
    q   <- exp(lq)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    kel <- exp(lkel)

    # 2. Derived rate constants (Kwak 2025 Supplementary Material section 4,
    #    'Related equations'). koff is not estimated: it is pinned by the
    #    in-vitro KD so that kd = koff / kon. The two synthesis rates are
    #    pinned by the steady-state balance at the reported baselines.
    koff  <- kd * kon
    ksyns <- kdegs * stfpi_b
    ksynm <- kdegm * mtfpi_b

    # 3. Free MG1113 concentrations (nM) from the central and peripheral
    #    amounts. The paper writes its central / peripheral ODEs directly on
    #    concentration (C2MG1113, C3MG1113); carrying the two drug states as
    #    AMOUNTS instead lets a dose land on either `central` (i.v.) or
    #    `depot` (s.c.) through the ordinary rxode2 event table. The two forms
    #    are identical: each ODE below is the paper's equation multiplied
    #    through by its compartment volume.
    cfree <- central / vc
    cper  <- peripheral1 / vp

    # 4. ODE system
    #    -- Subcutaneous absorption (Supplementary Material section 4,
    #       'Subcutaneous administration'). Amounts in nmol.
    d/dt(depot)       <- -(k_depot_transit1 + k_depot_central + kloss) * depot
    d/dt(transit1)    <- k_depot_transit1 * depot - k_transit1_central * transit1

    #    -- Free MG1113 in the central compartment. Paper form (s.c.):
    #       dC2/dt = (k02*A0 + k12*A1 + CLD*(C3 - C2))/V2
    #                + koff*(C2CM + C2CS) - kon*C2*(mTFPI + sTFPI_alpha)
    #                - kel_MG1113*C2
    #       multiplied through by V2. Under i.v. dosing depot and transit1 are
    #       empty, so the first two terms vanish and the equation reduces to
    #       the paper's 'Intravenous administration' form.
    d/dt(central)     <- k_depot_central * depot +
                         k_transit1_central * transit1 +
                         q * (cper - cfree) +
                         vc * koff * (amtfpi + astfpi) -
                         vc * kon * cfree * (mtfpi + stfpi) -
                         kel * central

    #    -- Free MG1113 in the peripheral compartment.
    #       Paper form: dC3/dt = CLD*(C2 - C3)/V3, multiplied through by V3.
    d/dt(peripheral1) <- q * (cfree - cper)

    #    -- Free membrane-bound TFPI (nM) in the central compartment.
    d/dt(mtfpi)       <- ksynm - kdegm * mtfpi + koff * amtfpi - kon * mtfpi * cfree

    #    -- Free soluble TFPI-alpha (nM) in the central compartment.
    d/dt(stfpi)       <- ksyns - kdegs * stfpi + koff * astfpi - kon * stfpi * cfree

    #    -- MG1113 / membrane-bound TFPI complex (C2CM, nM).
    d/dt(amtfpi)      <- kon * mtfpi * cfree - (koff + kintm) * amtfpi

    #    -- MG1113 / soluble TFPI-alpha complex (C2CS, nM).
    d/dt(astfpi)      <- kon * stfpi * cfree - (koff + kints) * astfpi

    # 5. Endogenous initial conditions (Supplementary Material section 4).
    #    The two complexes and all drug states start empty (rxode2 default).
    mtfpi(0) <- mtfpi_b
    stfpi(0) <- stfpi_b

    # 6. Outputs (nM). Both assays detect the FREE form (Kwak 2025 section
    #    2.1), so neither output adds the bound species back in. The paper
    #    reports no residual-error model, so none is applied.
    Cc    <- cfree
    sTFPI <- stfpi
  })
}
