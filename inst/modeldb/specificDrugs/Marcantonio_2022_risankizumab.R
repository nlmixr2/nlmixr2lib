Marcantonio_2022_risankizumab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "risankizumab-p19 binding in adults with plaque psoriasis (Marcantonio",
    "2022 Early Feasibility Assessment, Case Study 3). Risankizumab binds",
    "the p19 subunit of IL-23 (specific to IL-23; distinguished from",
    "ustekinumab which binds the shared p40). Bivalent binding (assumed",
    "valency = 2 as a monoclonal antibody). Structure identical to the",
    "Marcantonio 2022 anti-ligand family; parameters FIXED from paper",
    "Table S4.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 3 (risankizumab IL-23 p19,",
    "plaque psoriasis). Drug-specific parameters from paper Supplementary",
    "Table S4.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Risankizumab dose amount into depot (SC) in nmol; MW = 145610 Da so 150 mg = 1030 nmol.",
    concentration = "Free risankizumab plasma concentration Cc = Ab_00 / V in nM; V = 5 L."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "risankizumab", units = NA_character_, specimen = "administration site", verified = FALSE),
    Ab_00 = list(analyte = "free risankizumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_0L = list(analyte = "bound risankizumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_L0 = list(analyte = "risankizumab-p19 complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_LL = list(analyte = "dimerized risankizumab-p19 complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    L1    = list(analyte = "p19 subunit of IL-23", units = NA_character_, specimen = "plasma", verified = FALSE),
    R1    = list(analyte = "risankizumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    L1R1  = list(analyte = "risankizumab-p19 complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    S1    = list(analyte = "IL-23", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with moderate-to-severe plaque psoriasis.",
    dose_range     = "150 mg SC at weeks 0 and 4, then every 12 weeks (Skyrizi USPI; Ph3 dose per Gordon 2018). Paper Table 5 reports both a Q4W and Q12W ID90 prediction.",
    regions        = NA_character_,
    notes          = "See sibling Marcantonio 2022 anti-ligand models for shared methodology."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption rate constant (1/day; from t1/2 = 2.5 days)")             # Marcantonio 2022 Assess default
    lkclearAb <- fixed(log(log(2) / 27));          label("First-order risankizumab elimination rate constant (1/day; from t1/2 = 27 days)")   # Marcantonio 2022 Table S4 Drug Half Life (Suleiman 2019)
    kd_drug   <- fixed(0.001);                     label("Risankizumab-p19 equilibrium dissociation constant (nM)")                            # Marcantonio 2022 Table S4 KD = 1 pM (Singh 2015; assumed < 10 pM)
    valency   <- fixed(2);                         label("Effective drug valency for p19 binding (dimensionless)")                             # Marcantonio 2022 Table S4 Drug Valency = 2 (assumed monoclonal antibody)

    L1_conc      <- fixed(0.174);                  label("Baseline plasma IL-23 p19 concentration (nM)")                                         # Marcantonio 2022 Assess run file lig_css_1_central = 0.174 nM; paper Table S4 reports 872 pM but the Assess JSON that produced Table 5's 273 mg / 37.1 mg predictions uses 0.174 nM = 174 pM. Following JSON as authoritative.
    R1_conc      <- fixed(0.88);                   label("Baseline p19-receptor concentration (nM; 4.4 nmols total-body per 5 L)")             # Marcantonio 2022 Table S4 Receptor Concentration = 4.4 nmols; 4.4/5 = 0.88 nM
    S1_conc      <- fixed(0);                      label("Baseline soluble-shed receptor concentration (nM)")                                     # Marcantonio 2022 Assess run file default
    lkclearL1    <- fixed(log(log(2) / (90 / 1440)));  label("First-order p19 elimination rate constant (1/day; from t1/2 = 1.5 hr)")             # Marcantonio 2022 Table S4 p-19 Half Life (Lotze 1985; assumed IL-2-like)
    lkclearR1    <- fixed(log(log(2) / (120 / 1440))); label("First-order p19-receptor elimination rate constant (1/day; from t1/2 = 2 hr)")     # Marcantonio 2022 Table S4 Receptor Internalization = 2 h (standard assumption)
    lkclearS1    <- fixed(log(log(2) / (30 / 1440)));  label("First-order shed-receptor elimination rate constant (1/day; unused)")               # Assess default; unused when S1_conc = 0
    kd_lr        <- fixed(1.6);                    label("p19:Receptor equilibrium dissociation constant (nM)")                                   # Marcantonio 2022 Table S4 p-19:Receptor KD = 1.6 nM (Parham 2002, midpoint 0.3-3 nM)

    V           <- fixed(5);                       label("Central compartment volume of distribution (L)")                                        # Marcantonio 2022 Table S4 Volume
    kon         <- fixed(0.001 * 86400);           label("Bimolecular association rate constant (L/nmol/day)")                                    # Assess default

    propSd    <- fixed(0.01);                      label("Placeholder proportional residual error on Cc (not paper-derived)")                    # Placeholder
  })

  model({
    ka         <- exp(lka)
    kclearAb   <- exp(lkclearAb)
    kclearL1   <- exp(lkclearL1)
    kclearR1   <- exp(lkclearR1)
    kclearS1   <- exp(lkclearS1)
    kclearL1R1 <- kclearR1
    kon1Ab     <- kon
    kon2Ab     <- floor(valency / 2) * kon
    koffAb     <- kon * kd_drug
    koffL1R1   <- kon * kd_lr

    total_R1 <- R1_conc * V
    L1_0     <- L1_conc * V
    S1_0     <- S1_conc * V
    L1R1_0   <- (kon * L1_0 * total_R1 / V) /
                  (kon * L1_0 / V + koffL1R1 + kclearR1)
    R1_0     <- total_R1 - L1R1_0
    ksynth_L1 <- kon * L1_0 * R1_0 / V - koffL1R1 * L1R1_0 + kclearL1 * L1_0
    kshed_R1  <- kclearS1 * S1_0 / R1_0
    ksynth_R1 <- kon * L1_0 * R1_0 / V - koffL1R1 * L1R1_0 +
                   kclearR1 * R1_0 + kshed_R1 * R1_0

    d/dt(depot) <- -ka * depot

    d/dt(Ab_00) <-  ka * depot -
                    kclearAb * Ab_00 -
                    kon2Ab * L1 * Ab_00 / V + koffAb * Ab_0L -
                    kon1Ab * L1 * Ab_00 / V + koffAb * Ab_L0

    d/dt(Ab_0L) <-  kon2Ab * L1 * Ab_00 / V - koffAb * Ab_0L -
                    kclearAb * Ab_0L -
                    kon1Ab * L1 * Ab_0L / V + koffAb * Ab_LL

    d/dt(Ab_L0) <-  kon1Ab * L1 * Ab_00 / V - koffAb * Ab_L0 -
                    kclearAb * Ab_L0 -
                    kon2Ab * L1 * Ab_L0 / V + koffAb * Ab_LL

    d/dt(Ab_LL) <-  kon1Ab * L1 * Ab_0L / V - koffAb * Ab_LL +
                    kon2Ab * L1 * Ab_L0 / V - koffAb * Ab_LL -
                    kclearAb * Ab_LL

    d/dt(L1)   <-  ksynth_L1 - kclearL1 * L1 -
                    kon * L1 * R1 / V + koffL1R1 * L1R1 -
                    kon2Ab * L1 * Ab_00 / V + koffAb * Ab_0L -
                    kon1Ab * L1 * Ab_00 / V + koffAb * Ab_L0 -
                    kon1Ab * L1 * Ab_0L / V + koffAb * Ab_LL -
                    kon2Ab * L1 * Ab_L0 / V + koffAb * Ab_LL

    d/dt(R1)   <-  ksynth_R1 - kclearR1 * R1 - kshed_R1 * R1 -
                    kon * L1 * R1 / V + koffL1R1 * L1R1

    d/dt(L1R1) <-  kon * L1 * R1 / V - koffL1R1 * L1R1 -
                    kclearL1R1 * L1R1

    d/dt(S1)   <-  kshed_R1 * R1 - kclearS1 * S1

    L1(0)   <- L1_0
    R1(0)   <- R1_0
    L1R1(0) <- L1R1_0
    S1(0)   <- S1_0

    Cc <- Ab_00 / V
    Cc ~ prop(propSd)
  })
}
