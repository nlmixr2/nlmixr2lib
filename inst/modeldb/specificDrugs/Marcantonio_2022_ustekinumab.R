Marcantonio_2022_ustekinumab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "ustekinumab-p40 binding in adults with plaque psoriasis (Marcantonio",
    "2022 Early Feasibility Assessment, Case Study 3). Ustekinumab binds",
    "the shared p40 subunit of IL-12 and IL-23 with bivalent binding",
    "(valency = 2 per Benson 2011); model treats p40 as the soluble target",
    "and its cognate receptor (IL-12R-beta1) as the membrane binding",
    "partner. Structure identical to the Marcantonio 2022 anti-ligand",
    "family; parameters FIXED from paper Table S3.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 3 (ustekinumab IL-12/IL-23",
    "p40, plaque psoriasis). Drug-specific parameters from paper",
    "Supplementary Table S3.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Ustekinumab dose amount into depot (SC) in nmol; MW = 148600 Da so 45 mg = 303.5 nmol.",
    concentration = "Free ustekinumab plasma concentration Cc = Ab_00 / V in nM; V = 5 L."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with moderate-to-severe plaque psoriasis.",
    dose_range     = "45 mg SC (60-100 kg subjects) or 90 mg SC (>100 kg) at weeks 0 and 4, then every 12 weeks (Stelara USPI). Model prediction targets the 12-week maintenance schedule.",
    regions        = NA_character_,
    notes          = "See sibling Marcantonio 2022 anti-ligand models for shared methodology."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption rate constant (1/day; from t1/2 = 2.5 days)")             # Marcantonio 2022 Assess default; typical mAb SC t1/2 (Kagan 2014)
    lkclearAb <- fixed(log(log(2) / 21.6));        label("First-order ustekinumab elimination rate constant (1/day; from t1/2 = 21.6 days)") # Marcantonio 2022 Table S3 Drug Half Life (Zhu 2009)
    kd_drug   <- fixed(0.001);                     label("Ustekinumab-p40 equilibrium dissociation constant (nM)")                            # Marcantonio 2022 Table S3 KD = 1 pM (Luo 2010; assumed < 10 pM)
    valency   <- fixed(2);                         label("Effective drug valency for p40 binding (dimensionless)")                             # Marcantonio 2022 Table S3 Drug Valency = 2 (Benson 2011)

    L1_conc      <- fixed(1.9e-3 / 5);             label("Baseline p40 concentration (nM; 1.9 pmols total-body divided by 5 L)")               # Marcantonio 2022 Table S3 p-40 Concentration = 1.9 pmols (bottom-up); Assess JSON lig_css_1_central = 0.00038 nM
    R1_conc      <- fixed(61e-3 / 5);              label("Baseline p40-receptor concentration (nM; 61 pmols total-body divided by 5 L)")       # Marcantonio 2022 Table S3 Receptor Concentration = 61 pmols (bottom-up); Assess JSON rec_css_1_central = 0.0122 nM
    S1_conc      <- fixed(0);                      label("Baseline soluble-shed receptor concentration (nM)")                                     # Marcantonio 2022 Assess run file default
    lkclearL1    <- fixed(log(log(2) / (90 / 1440)));  label("First-order p40 elimination rate constant (1/day; from t1/2 = 1.5 hr)")             # Marcantonio 2022 Table S3 p-40 Half Life (Lotze 1985; assumed IL-2-like)
    lkclearR1    <- fixed(log(log(2) / (120 / 1440))); label("First-order p40-receptor elimination rate constant (1/day; from t1/2 = 2 hr)")     # Marcantonio 2022 Table S3 Receptor Internalization = 2 h (standard assumption)
    lkclearS1    <- fixed(log(log(2) / (30 / 1440)));  label("First-order shed-receptor elimination rate constant (1/day; unused)")               # Marcantonio 2022 Assess default; unused when S1_conc = 0
    kd_LR        <- fixed(0.01);                   label("p40:Receptor equilibrium dissociation constant (nM)")                                   # Marcantonio 2022 Table S3 p-40:Receptor KD = 10 pM (Ma 2001)

    V           <- fixed(5);                       label("Central compartment volume of distribution (L)")                                        # Marcantonio 2022 Table S3 Volume
    kon         <- fixed(0.001 * 86400);           label("Bimolecular association rate constant (L/nmol/day)")                                    # Marcantonio 2022 Assess default (0.001 nM-1 s-1)

    propSd    <- fixed(0.01);                      label("Placeholder proportional residual error on Cc (not paper-derived)")                    # Not paper-derived; placeholder
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
    koffL1R1   <- kon * kd_LR

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
