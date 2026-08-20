Marcantonio_2022_belimumab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "belimumab-BAFF binding in adults with systemic lupus erythematosus (SLE)",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 3). Belimumab",
    "binds soluble BAFF (BLyS); model treats BAFF as the ligand and BAFF-R",
    "as the cognate receptor (BAFF also binds BCMA and TACI in vivo but this",
    "model uses only the BAFF-R affinity 15 nM per Day 2005 / Cachero 2006",
    "/ Hymowitz 2005 as reported in Table S5). Belimumab valency = 1",
    "(BAFF trimer binds 3:3 per Shin 2018). Parameters FIXED from paper",
    "Table S5.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 3 (belimumab BAFF, SLE).",
    "Drug-specific parameters from paper Supplementary Table S5.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Belimumab dose amount into depot (SC) or Ab_00 (IV) in nmol; MW = 147000 Da so 200 mg = 1361 nmol.",
    concentration = "Free belimumab plasma concentration Cc = Ab_00 / V in nM; V = 5 L."
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "Belimumab", units = NA_character_, specimen = "administration site", verified = FALSE),
    Ab_00 = list(analyte = "Belimumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_0L = list(analyte = "Belimumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_L0 = list(analyte = "Belimumab-BAFF complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    Ab_LL = list(analyte = "Belimumab-BAFF complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    L1    = list(analyte = "soluble BAFF (BLyS)", units = NA_character_, specimen = "plasma", verified = FALSE),
    R1    = list(analyte = "BAFF-R", units = NA_character_, specimen = "not applicable", verified = FALSE),
    L1R1  = list(analyte = "Belimumab-BAFF-BAFF-R complex", units = NA_character_, specimen = "plasma", verified = FALSE),
    S1    = list(analyte = "Free soluble BAFF (BLyS)", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with systemic lupus erythematosus (SLE).",
    dose_range     = "200 mg SC every week or 10 mg/kg IV Q4W (Benlysta USPI). Model prediction targets the weekly SC schedule; a Q4W IV prediction is also reported.",
    regions        = NA_character_,
    notes          = "See sibling Marcantonio 2022 anti-ligand models for shared methodology."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption rate constant (1/day; from t1/2 = 2.5 days)")            # Assess default
    lkclearAb <- fixed(log(log(2) / 19.4));        label("First-order belimumab elimination rate constant (1/day; from t1/2 = 19.4 days)")   # Marcantonio 2022 Table S5 Drug Half Life (SumR BLA; linear PK)
    kd_drug   <- fixed(0.274);                     label("Belimumab-BAFF equilibrium dissociation constant (nM)")                             # Marcantonio 2022 Table S5 KD = 274 pM (PharmR BLA)
    valency   <- fixed(1);                         label("Effective drug valency for BAFF binding (dimensionless)")                            # Marcantonio 2022 Table S5 Drug Valency = 1 (BAFF trimer 3:3 stoichiometry, Shin 2018)

    L1_conc      <- fixed(0.294);                  label("Baseline plasma BAFF concentration (nM)")                                             # Marcantonio 2022 Table S5 BAFF Concentration = 1.47 nmols; 1.47/5 = 0.294 nM
    R1_conc      <- fixed(0.058);                  label("Baseline BAFF-R receptor concentration (nM)")                                          # Marcantonio 2022 Table S5 Receptor Concentration = 290 pmols; 290e-3/5 = 0.058 nM
    S1_conc      <- fixed(0);                      label("Baseline soluble-shed BAFF-R concentration (nM)")                                       # Assess run file default
    lkclearL1    <- fixed(log(log(2) / (60 / 1440)));  label("First-order BAFF elimination rate constant (1/day; from t1/2 = 60 min)")           # Marcantonio 2022 Table S5 BAFF Half Life (Moritz 1989, TNF-alpha analog)
    lkclearR1    <- fixed(log(log(2) / (120 / 1440))); label("First-order BAFF-R elimination rate constant (1/day; from t1/2 = 2 hr)")           # Marcantonio 2022 Table S5 Receptor Internalization = 2 h (standard)
    lkclearS1    <- fixed(log(log(2) / (30 / 1440)));  label("First-order shed-receptor elimination rate constant (1/day; unused)")               # Assess default
    kd_lr        <- fixed(15);                     label("BAFF:BAFF-R equilibrium dissociation constant (nM)")                                   # Marcantonio 2022 Table S5 BAFF:Receptor KD = 15 nM for BAFF-R (Day 2005); BCMA 1550 nM and TACI 1.3 nM not modelled

    V           <- fixed(5);                       label("Central compartment volume of distribution (L)")                                        # Table S5 Volume
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
