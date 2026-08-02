Marcantonio_2022_adalimumab <- function() {
  description <- paste(
    "QSP. One-compartment monospecific anti-ligand mechanistic PKPD model of",
    "adalimumab-TNF-alpha binding in adults with rheumatoid arthritis",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 1). Drug",
    "administered SC (first-order absorption) or IV, binds soluble TNF-alpha",
    "(monovalent, effective valency = 1 per Tran 2017 / Lim 2018), and blocks",
    "the reversible TNF:TNFR interaction. All species are eliminated by",
    "first-order processes; the drug-target-ligand complex eliminates at the",
    "same rate as the free drug; the ligand-receptor complex eliminates at",
    "the same rate as the free receptor (paper Model Assumptions). All",
    "parameters are FIXED from paper Table 2 (no fitting, no IIV, no",
    "residual error reported). Initial conditions for L1 (free TNF), R1",
    "(free TNFR), L1R1 (TNF-TNFR complex) and S1 (soluble shed receptor)",
    "are derived analytically from the reported steady-state concentrations",
    "and the quasi-equilibrium ligand-receptor binding balance. Simulation",
    "output activity (paper active_1_central) is the L1R1 complex amount;",
    "target inhibition is computed as 1 - L1R1(t)/L1R1(0). Run in the",
    "Applied BioMath Assess Monospecific Anti-Ligand model; the ODE system",
    "encoded here reproduces the Assess reaction network verbatim (see",
    "Supplementary Material Data Sheet 2, one_compartment_anti_ligand.pdf).",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH, Matteson A, Presler M, Burke JM, Hagen DR, Hua F,",
    "Apgar JF (2022). Early Feasibility Assessment: A Method for",
    "Accurately Predicting Biotherapeutic Dosing to Inform Early Drug",
    "Discovery Decisions. Frontiers in Pharmacology 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 1 (adalimumab TNF-alpha,",
    "rheumatoid arthritis). Drug-specific parameters from paper Table 2;",
    "absorption half-life from paper Case Study 1 text.",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = paste(
      "Dose amount into the depot compartment (SC) or Ab_00 compartment (IV)",
      "must be in nmol (convert from mg via amt_nmol = amt_mg * 1e6 / MW_Da;",
      "for adalimumab MW = 148000 Da, 40 mg = 270.3 nmol).",
      sep = " "
    ),
    concentration = paste(
      "Free adalimumab plasma concentration Cc = Ab_00 / V in nM.",
      "Soluble TNF and target-receptor amounts are in nmol; target complex",
      "L1R1 is in nmol; typical plasma volume V = 5 L.",
      sep = " "
    )
  )

  covariateData <- list()
  covariatesDataExcluded <- list(
    BW = list(
      description = "Body weight; the paper assumes a typical 70 kg adult for mg/kg-per-kg dose conversion but body weight does not appear in the ODE system.",
      units       = "kg",
      type        = "continuous",
      notes       = "Paper Table 2 lists BW = 70 kg only for converting mg-per-kg dosing to mg. Not a model covariate."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with active rheumatoid arthritis (RA); clinically approved indication for adalimumab.",
    dose_range     = "40 mg SC every other week (approved RA dose); some patients not receiving methotrexate benefit from 40 mg every week (Adalimumab USPI 2021).",
    regions        = NA_character_,
    notes          = paste(
      "Marcantonio 2022 is a methodological Early Feasibility Assessment",
      "(EFA) paper that uses mechanistic PKPD models parameterised entirely",
      "from literature (no fitting, no clinical dataset). The 'population'",
      "is therefore the target indication rather than a fitted cohort: an",
      "adult RA patient with typical body weight 70 kg and typical",
      "monoclonal-antibody volume of distribution 5 L. Baseline soluble",
      "TNF-alpha concentration is the median RA-patient plasma level",
      "(Takeuchi 2011), TNFR1 receptor expression is a bottom-up estimate",
      "from Scatchard analysis (Imamura 1987, Michishita 1990) across a high",
      "percentage of nucleated cells (Sender 2016), and TNFR1 turnover is",
      "from Higuchi 1994.",
      sep = " "
    )
  )

  ini({
    # -------------------------------------------------------------------------
    # Drug-specific parameters (adalimumab). All FIXED (no fitting; the paper
    # reports no uncertainty and no IIV). Values from Marcantonio 2022 Table 2.
    # -------------------------------------------------------------------------
    lka       <- fixed(log(log(2) / 2.5));    label("First-order SC absorption rate constant (1/day; from t1/2 = 2.5 days)")           # Marcantonio 2022 Case Study 1 text; assumed typical mAb SC t1/2 (Kagan 2014)
    lkclearAb <- fixed(log(log(2) / 20));     label("First-order adalimumab elimination rate constant (1/day; from t1/2 = 20 days)")   # Marcantonio 2022 Table 2 Adalimumab Half-Life (Weisman 2003, Ternant 2015)
    kd_drug   <- fixed(0.0086);               label("Adalimumab-TNF equilibrium dissociation constant (nM)")                            # Marcantonio 2022 Table 2 Adalimumab KD = 8.6 pM (Kaymakcalan 2009); Assess JSON stores 0.008 nM (~5% rounding)
    valency   <- fixed(1);                    label("Effective drug valency for TNF-alpha binding (dimensionless)")                     # Marcantonio 2022 Table 2 (Tran 2017, Lim 2018); floor(valency/2) gates second-arm kinetics

    # -------------------------------------------------------------------------
    # Target-specific parameters (TNF-alpha ligand and TNFR1 receptor).
    # -------------------------------------------------------------------------
    L1_conc      <- fixed(5.73e-5);             label("Baseline plasma soluble TNF-alpha concentration in RA patients (nM)")            # Marcantonio 2022 Table 2 TNF Concentration (Takeuchi 2011)
    R1_conc      <- fixed(0.23);                label("Baseline total TNFR1 concentration (nM; bottom-up estimate)")                    # Marcantonio 2022 Table 2 TNFR Concentration (bottom-up per methods)
    S1_conc      <- fixed(0);                   label("Baseline soluble-shed TNFR concentration (nM)")                                    # Marcantonio 2022 Assess run file default (shed_css_1_central = 0)
    lkclearL1    <- fixed(log(log(2) / (30 / 1440))); label("First-order TNF-alpha elimination rate constant (1/day; from t1/2 = 30 min)") # Marcantonio 2022 Table 2 TNF Half-Life (Moritz 1989)
    lkclearR1    <- fixed(log(log(2) / (540 / 1440))); label("First-order TNFR1 elimination rate constant (1/day; from t1/2 = 9 hr)")     # Marcantonio 2022 Table 2 TNFR receptor half-life (Higuchi 1994); 540 min = 9 hr
    lkclearS1    <- fixed(log(log(2) / (30 / 1440))); label("First-order soluble-shed-TNFR elimination rate constant (1/day)")           # Marcantonio 2022 Assess run file default (shed_half_1 = 30 min); unused when S1_conc = 0
    kd_lr        <- fixed(0.019);               label("TNF:TNFR1 equilibrium dissociation constant (nM)")                                 # Marcantonio 2022 Table 2 TNF:TNFR KD = 19 pM (Grell 1998)

    # -------------------------------------------------------------------------
    # System parameters.
    # -------------------------------------------------------------------------
    V           <- fixed(5);                    label("Central compartment volume of distribution (L)")                                   # Marcantonio 2022 Table 2 Volume; typical mAb Vd (Pearson 1995, Ovacik 2018)
    kon         <- fixed(0.001 * 86400);        label("Bimolecular association rate constant for all reversible binding (L/nmol/day)")   # Marcantonio 2022 Assess default (kon = 0.001 nM-1 s-1 converted to 1/nM/day)

    # -------------------------------------------------------------------------
    # Placeholder residual error. The paper reports no residual variability
    # (deterministic simulation of typical trajectories). A small proportional
    # SD is retained here for the nlmixr2 UI observation contract; simulations
    # that use rxode2::zeroRe() suppress it.
    # -------------------------------------------------------------------------
    propSd    <- fixed(0.01);                 label("Placeholder proportional residual error on Cc (not paper-derived)")               # Not paper-derived; placeholder for nlmixr2 error-model contract
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
    # Steady-state initial conditions for the target-side species (paper's
    # Assess % relationships block; the quadratic-form-free quasi-equilibrium
    # formulas apply because at t=0 no drug is bound).
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
    # ODE system (Marcantonio 2022 % reactions block for the 1-compartment
    # monospecific anti-ligand model; time in days, amounts in nmol).
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

    # -------------------------------------------------------------------------
    # Steady-state initial values.
    # -------------------------------------------------------------------------
    L1(0)   <- L1_0
    R1(0)   <- R1_0
    L1R1(0) <- L1R1_0
    S1(0)   <- S1_0

    # -------------------------------------------------------------------------
    # Observations. Cc = free adalimumab concentration in central (nM).
    # -------------------------------------------------------------------------
    Cc <- Ab_00 / V
    Cc ~ prop(propSd)
  })
}
