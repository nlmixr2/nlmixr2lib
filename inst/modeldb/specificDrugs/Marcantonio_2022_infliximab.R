Marcantonio_2022_infliximab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "infliximab-TNF-alpha binding in adults with rheumatoid arthritis",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 1).",
    "Structure identical to the Marcantonio 2022 adalimumab model; differs",
    "only in the drug-specific parameters (KD, half-life, molecular weight,",
    "route). Infliximab is administered IV (dose directly into Ab_00), binds",
    "soluble TNF-alpha monovalently (effective valency = 1 per Tran 2017 /",
    "Lim 2018), and blocks the reversible TNF:TNFR interaction. All",
    "parameters are FIXED from paper Table 2 (no fitting, no IIV, no",
    "residual error reported). Run in the Applied BioMath Assess",
    "Monospecific Anti-Ligand model.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH, Matteson A, Presler M, Burke JM, Hagen DR, Hua F,",
    "Apgar JF (2022). Early Feasibility Assessment: A Method for",
    "Accurately Predicting Biotherapeutic Dosing to Inform Early Drug",
    "Discovery Decisions. Frontiers in Pharmacology 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 1 (infliximab TNF-alpha,",
    "rheumatoid arthritis, IV Q8W maintenance dosing). Drug-specific",
    "parameters from paper Table 2.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = paste(
      "Infliximab dose amount into Ab_00 (IV bolus) must be in nmol",
      "(convert from mg via amt_nmol = amt_mg * 1e6 / MW_Da; infliximab",
      "MW = 149100 Da, so 210 mg (=3 mg/kg x 70 kg) = 1408 nmol,",
      "441 mg (=6.3 mg/kg x 70 kg, the model-predicted ID90) = 2957 nmol.",
      sep = " "
    ),
    concentration = "Free infliximab plasma concentration Cc = Ab_00 / V in nM; typical plasma volume V = 5 L."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with active rheumatoid arthritis (RA); clinically approved indication for infliximab.",
    dose_range     = paste(
      "Clinically approved RA regimen: 3 mg/kg IV at weeks 0, 2, and 6,",
      "then maintenance every 8 weeks; may escalate to 10 mg/kg every 4",
      "weeks (Infliximab USPI 2013). Model prediction here targets the",
      "maintenance every-8-weeks dose.",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "See the sibling model Marcantonio_2022_adalimumab (Case Study 1) for",
      "shared target-side parameters (TNF-alpha soluble ligand and TNFR1",
      "membrane receptor turnover), the joint parameter-scan sensitivity",
      "analysis, and the paper's rationale for a mechanistic rather than",
      "exposure-only dose prediction.",
      sep = " "
    )
  )

  ini({
    # -------------------------------------------------------------------------
    # Drug-specific parameters (infliximab). All FIXED (no fitting; the paper
    # reports no uncertainty and no IIV). Values from Marcantonio 2022 Table 2.
    # -------------------------------------------------------------------------
    lka       <- fixed(log(log(2) / 2.5));    label("First-order SC absorption rate constant placeholder (1/day; unused for IV dosing)") # Not paper-derived for IV infliximab; kept identical to sibling SC model for structural consistency
    lkclearAb <- fixed(log(log(2) / 14));     label("First-order infliximab elimination rate constant (1/day; from t1/2 = 14 days)")     # Marcantonio 2022 Table 2 Infliximab Half-Life (Hemperly 2018)
    kd_drug   <- fixed(0.0042);               label("Infliximab-TNF equilibrium dissociation constant (nM)")                              # Marcantonio 2022 Table 2 Infliximab KD = 4.2 pM (Kaymakcalan 2009)
    valency   <- fixed(1);                    label("Effective drug valency for TNF-alpha binding (dimensionless)")                       # Marcantonio 2022 Table 2 Drug Valency = 1 (Tran 2017, Lim 2018)

    # -------------------------------------------------------------------------
    # Target-specific parameters (TNF-alpha and TNFR1; identical to the
    # Marcantonio 2022 adalimumab model since both drugs share the target).
    # -------------------------------------------------------------------------
    L1_conc      <- fixed(5.73e-5);             label("Baseline plasma soluble TNF-alpha concentration in RA patients (nM)")              # Marcantonio 2022 Table 2 TNF Concentration (Takeuchi 2011)
    R1_conc      <- fixed(0.23);                label("Baseline total TNFR1 concentration (nM; bottom-up estimate)")                      # Marcantonio 2022 Table 2 TNFR Concentration (bottom-up per methods)
    S1_conc      <- fixed(0);                   label("Baseline soluble-shed TNFR concentration (nM)")                                      # Marcantonio 2022 Assess run file default (shed_css_1_central = 0)
    lkclearL1    <- fixed(log(log(2) / (30 / 1440)));  label("First-order TNF-alpha elimination rate constant (1/day; from t1/2 = 30 min)") # Marcantonio 2022 Table 2 TNF Half-Life (Moritz 1989)
    lkclearR1    <- fixed(log(log(2) / (540 / 1440))); label("First-order TNFR1 elimination rate constant (1/day; from t1/2 = 9 hr)")       # Marcantonio 2022 Table 2 TNFR receptor half-life (Higuchi 1994)
    lkclearS1    <- fixed(log(log(2) / (30 / 1440)));  label("First-order soluble-shed-TNFR elimination rate constant (1/day)")             # Marcantonio 2022 Assess run file default (shed_half_1 = 30 min); unused when S1_conc = 0
    kd_lr        <- fixed(0.019);               label("TNF:TNFR1 equilibrium dissociation constant (nM)")                                   # Marcantonio 2022 Table 2 TNF:TNFR KD = 19 pM (Grell 1998)

    # -------------------------------------------------------------------------
    # System parameters.
    # -------------------------------------------------------------------------
    V           <- fixed(5);                    label("Central compartment volume of distribution (L)")                                     # Marcantonio 2022 Table 2 Volume; typical mAb Vd (Pearson 1995, Ovacik 2018)
    kon         <- fixed(0.001 * 86400);        label("Bimolecular association rate constant for all reversible binding (L/nmol/day)")     # Marcantonio 2022 Assess default (kon = 0.001 nM-1 s-1 converted to 1/nM/day)

    # Placeholder proportional residual error.
    propSd    <- fixed(0.01);                 label("Placeholder proportional residual error on Cc (not paper-derived)")                # Not paper-derived; placeholder for nlmixr2 error-model contract
  })

  model({
    # -------------------------------------------------------------------------
    # Rate constants and equilibrium constants (all derived from ini() values).
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Steady-state initial conditions (identical algebra as adalimumab).
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # ODE system (identical structure as adalimumab; monovalent so Ab_0L and
    # Ab_LL remain at 0 throughout).
    # -------------------------------------------------------------------------
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
