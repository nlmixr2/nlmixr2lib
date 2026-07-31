Marcantonio_2022_omalizumab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "omalizumab-IgE binding in adults with moderate-to-severe allergic asthma",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 3). Omalizumab",
    "binds soluble IgE; model treats IgE as the ligand and FcepsilonRI",
    "(high-affinity IgE receptor on basophils and mast cells) as the cognate",
    "receptor. Bivalent binding (crystallography suggests 2:1 stoichiometry",
    "with some 2:2 and 3:3 populations per Davies 2017). Parameters FIXED",
    "from paper Table S6; the receptor half-life used is the Assess default",
    "of 15 min rather than the Table S6 alternative of 2 hr (paper explicitly",
    "tests both).",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 3 (omalizumab IgE, asthma).",
    "Drug-specific parameters from paper Supplementary Table S6.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Omalizumab dose amount into depot (SC) in nmol; MW = 149000 Da so 225 mg = 1510 nmol.",
    concentration = "Free omalizumab plasma concentration Cc = Ab_00 / V in nM; V = 5 L."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with moderate-to-severe persistent allergic asthma; the IgE concentration used here (40 nmols total = 8 nM at V = 5 L) corresponds to a patient with serum IgE 222 IU/mL = 533 ng/mL per Table S6.",
    dose_range     = "150 mg or 300 mg SC every 2 or 4 weeks per body weight and serum IgE (Xolair USPI). Model prediction targets 225 mg Q2W SC.",
    regions        = NA_character_,
    notes          = "See sibling Marcantonio 2022 anti-ligand models for shared methodology."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));           label("First-order SC absorption rate constant (1/day; from t1/2 = 2.5 days)")            # Assess default
    lkclearAb <- fixed(log(log(2) / 26));            label("First-order omalizumab elimination rate constant (1/day; from t1/2 = 26 days)")    # Marcantonio 2022 Table S6 Drug Half Life (Xolair USPI)
    kd_drug   <- fixed(0.02);                        label("Omalizumab-IgE equilibrium dissociation constant (nM)")                            # Marcantonio 2022 Table S6 KD = 20 pM (Chang 2000; assumed 20 pM given avidity effects)
    valency   <- fixed(2);                           label("Effective drug valency for IgE binding (dimensionless)")                             # Marcantonio 2022 Table S6 Drug Valency = 2 (Davies 2017 crystallography)

    L1_conc      <- fixed(8);                        label("Baseline plasma free IgE concentration (nM)")                                         # Marcantonio 2022 Table S6 IgE Concentration = 40 nmols total-body; 40/5 = 8 nM (corresponds to 533 ng/mL)
    R1_conc      <- fixed(0.78);                     label("Baseline FcepsilonRI receptor concentration (nM)")                                     # Marcantonio 2022 Table S6 Receptor Concentration = 3.9 nmols; 3.9/5 = 0.78 nM
    S1_conc      <- fixed(0);                        label("Baseline soluble-shed FcepsilonRI concentration (nM)")                                 # Assess default
    lkclearL1    <- fixed(log(log(2) / (2880 / 1440))); label("First-order IgE elimination rate constant (1/day; from t1/2 = 2 days)")             # Marcantonio 2022 Table S6 IgE Half Life = 2 days (Corne 1997); 2880 min
    lkclearR1    <- fixed(log(log(2) / (15 / 1440)));   label("First-order FcepsilonRI elimination rate constant (1/day; from t1/2 = 15 min)")     # Marcantonio 2022 Assess run file rec_half_1 = 15 min; Table S6 lists a 2 hr alternative that the paper explicitly says was 'tried both' (occupancy of soluble IgE vs FcepsilonRI inhibition); the Assess-JSON value of 15 min drives the paper's Table 5 dose prediction
    lkclearS1    <- fixed(log(log(2) / (30 / 1440)));   label("First-order shed-receptor elimination rate constant (1/day; unused)")               # Assess default
    kd_lr        <- fixed(0.12);                     label("IgE:FcepsilonRI equilibrium dissociation constant (nM)")                              # Marcantonio 2022 Table S6 IgE:Receptor KD = 0.12 nM (Miller 1989)

    V           <- fixed(5);                         label("Central compartment volume of distribution (L)")                                        # Table S6 Volume
    kon         <- fixed(0.001 * 86400);             label("Bimolecular association rate constant (L/nmol/day)")                                    # Assess default

    propSd    <- fixed(0.01);                        label("Placeholder proportional residual error on Cc (not paper-derived)")                    # Placeholder
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
