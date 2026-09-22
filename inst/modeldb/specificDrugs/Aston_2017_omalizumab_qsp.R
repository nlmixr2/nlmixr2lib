Aston_2017_omalizumab_qsp <- function() {
  description <- paste(
    "QSP. Target-mediated drug disposition (TMDD) with HOMEOSTATIC FEEDBACK on",
    "the free-receptor synthesis rate, parameterised for the anti-IgE mAb",
    "omalizumab (Example 1, Sect. 7.1 of Aston et al. 2017). The classical",
    "Levy / Mager-Jusko one-compartment TMDD system (free ligand L, free receptor",
    "R, ligand-receptor complex P) is extended with a fourth state, a feedback",
    "moderator F that scales the zero-order receptor synthesis rate and relaxes",
    "first-order toward a feedback function H(R) of the free-receptor level.",
    "The paper's result is that the moderator's response SPEED decides whether",
    "free receptor rebounds above baseline after a dose: with these omalizumab",
    "parameters ke(L) < ke(P), so the no-feedback and fast-feedback (direct,",
    "quasi-equilibrium) limits can never rebound, yet a slow moderator produces",
    "a rebound peaking near 117% of baseline. ktol is the paper's swept quantity",
    "(alpha); its default is the analytic bound of Theorem 5.12, which is also",
    "close to the maximum-rebound alpha of Fig. 10. Deterministic illustration:",
    "no subjects were fitted, so there is no inter-individual variability and no",
    "residual-error model.",
    sep = " "
  )
  reference <- paste(
    "Aston PJ, Derks G, Agoram BM, van der Graaf PH.",
    "A mathematical analysis of rebound in a target-mediated drug disposition",
    "model: II. With feedback.",
    "J Math Biol. 2017;75(1):39-73. doi:10.1007/s00285-016-1073-6.",
    "Companion model from the same paper: modellib('Aston_2017_efalizumab_qsp').",
    "The omalizumab parameter values are attributed by Aston et al. to",
    "Sun T (2001), poster, Advanced Methods of PKPD Systems Analysis, Los Angeles,",
    "and Agoram BM, Martin SW, van der Graaf PH, Drug Discov Today.",
    "2007;12(23-24):1018-1024, doi:10.1016/j.drudis.2007.10.002;",
    "they are reproduced in full in Sect. 7.1 of the present paper.",
    sep = " "
  )
  vignette <- "Aston_2017_receptor_rebound"

  # Aston et al. work in nM and days throughout (Sect. 2 and Sect. 7.1): every
  # state is a CONCENTRATION in nM, and the "initial drug dose L0 = 14.8148 nM"
  # of Sect. 7.1 is added directly to the free-ligand state as the initial
  # condition L(0) = L0 of Eq. (4). There is no volume anywhere in the model,
  # so a dose is supplied in nM straight into `central`.
  units <- list(time = "day", dosing = "nM", concentration = "nM")

  # `moderator1` is the canonical Gabrielsson-Hjorth moderator state: a
  # first-order delay driven by a system state (here the free receptor) with NO
  # mass transfer, whose value scales the production rate it modulates (here
  # ksyn). It is Aston's F, which is dimensionless with baseline F = 1.
  compartmentData <- list(
    central = list(analyte = "omalizumab (free)", units = "nM", specimen = "plasma", verified = TRUE),
    target = list(analyte = "IgE receptor / antigen (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex = list(analyte = "omalizumab-receptor complex", units = "nM", specimen = "plasma", verified = TRUE),
    moderator1 = list(
      analyte = "receptor-synthesis feedback moderator F (dimensionless, baseline 1)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = paste(
      "not applicable; the parameter values are literature estimates for the",
      "anti-IgE monoclonal antibody omalizumab, whose clinical indication is",
      "moderate-to-severe allergic asthma"
    ),
    dose_range = paste(
      "a single bolus of L0 = 14.8148 nM added to the free-ligand state",
      "(Sect. 7.1); the analysis is qualitative in the dose"
    ),
    notes = paste(
      "Aston et al. (2017) is a MATHEMATICAL ANALYSIS paper. Sect. 7.1 is a worked",
      "example, not a fit: the omalizumab parameters are taken from Sun (2001) and",
      "Agoram et al. (2007) and are reproduced verbatim in the text of Sect. 7.1.",
      "No subject-level data are fitted anywhere in the paper, so there is no",
      "demographic table, no inter-individual variability and no residual-error",
      "model. Every parameter is therefore encoded with fixed(). The paper's own",
      "numerical anchors for this example are the non-dimensional quantities",
      "lambda1 = -0.084, k3 = 0.517 (Sect. 7.1), the analytic rebound bound",
      "alpha < 0.135, the loss of rebound near alpha = 0.98, and the Fig. 10",
      "sweep of Rmax/R0 and tmax against alpha; the vignette reproduces all of them."
    )
  )

  ini({
    # ---- Elimination rates (Eqs. 1-3; values in Sect. 7.1) -------------------
    # Aston's ke(L): first-order elimination of the FREE ligand (drug).
    lkel <- fixed(log(0.024)); label("Free-drug (ligand) elimination rate ke(L) (1/day)") # Sect. 7.1: ke(L) = 0.024 /day
    # Aston's ke(P): first-order elimination of the ligand-receptor COMPLEX.
    # Named kint after the canonical TMDD complex-internalisation rate, which is
    # the role this rate plays in Eq. (3).
    kint <- fixed(0.201); label("Drug-receptor complex elimination rate ke(P) (1/day)") # Sect. 7.1: ke(P) = 0.201 /day
    # Aston's kout: first-order elimination of the FREE receptor.
    kdeg <- fixed(0.823); label("Free-receptor elimination rate kout (1/day)") # Sect. 7.1: kout = 0.823 /day

    # ---- Receptor baseline ---------------------------------------------------
    # R0 = kin / kout is the drug-free steady state of Eq. (2); Aston prints both
    # R0 and the implied kin = kout * R0 = 2.212224 nM/day, which model() derives.
    bl_target <- fixed(2.688); label("Baseline free-receptor concentration R0 (nM)") # Sect. 7.1: R0 = 2.688 nM

    # ---- Binding (Eqs. 1-3) ---------------------------------------------------
    kon <- fixed(0.592); label("Ligand-receptor association rate constant kon (1/(nM*day))") # Sect. 7.1: kon = 0.592 (nM day)^-1
    koff <- fixed(0.900); label("Ligand-receptor dissociation rate constant koff (1/day)") # Sect. 7.1: koff = 0.900 /day

    # ---- Feedback moderator (Eqs. 6 and 14) -----------------------------------
    # Aston's alpha: the response SPEED of the feedback moderator, and the swept
    # quantity of Sect. 7.1 / Fig. 10. alpha = 0 recovers the no-feedback TMDD
    # model and alpha -> Inf the direct (quasi-equilibrium) feedback limit;
    # neither can rebound with these elimination rates. The default below is the
    # analytic bound of Theorem 5.12 quoted in Sect. 7.1, 'rebound will happen
    # for alpha < 0.135'; Fig. 10 shows the numerical maximum of Rmax/R0 sits
    # close to it and that rebound in fact persists up to alpha = 0.98. Named
    # ktol after the canonical moderator-chain turnover rate constant, whose
    # d/dt(moderator) = ktol * (driver - moderator) shape this is exactly.
    ktol <- fixed(0.135); label("Feedback-moderator response rate alpha (1/day)") # Sect. 7.1: rebound for alpha < 0.135 (Theorem 5.12 bound)

    # Aston's H0: the slope of the 'mainly linear' feedback function of Eq. (34),
    #   H(R) = 1 + H0 * (R0 - R),
    # used in Sect. 7.1 with H0 = 1. The canonical sstim_<driver>_<target> family
    # is written with the driver's own deviation from baseline,
    #   1 + sstim * (driver - driver_baseline),
    # so sstim = -H0: the SIGN is negative because the feedback is negative
    # (synthesis rises when the receptor falls below R0). Sect. 7.1 notes that
    # the rebound is small enough here that the nonlinear branch H1(R) of
    # Eq. (34), which applies only for R > R0 + beta/H0, is never reached.
    sstim_target_moderator1 <- fixed(-1); label("Feedback slope on the moderator setpoint, -H0 (1/nM)") # Sect. 7.1 / Eq. (34): H0 = 1 (per nM), negative feedback

    # The source is a deterministic mathematical analysis and reports no
    # residual-error model and no inter-individual variability.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    kel <- exp(lkel)

    # Zero-order receptor synthesis holding the free receptor at R0 when the
    # moderator sits at its baseline of 1 (Eq. 4: R0 = kin / kout).
    # Sect. 7.1 prints the same product: kin = kout * R0 = 2.212224 nM/day.
    ksyn <- kdeg * bl_target

    # Drug-free steady state of Eqs. (11)-(14), with the extra condition
    # F(0) = 1 stated immediately after Eq. (14).
    target(0) <- bl_target
    complex(0) <- 0
    moderator1(0) <- 1

    # The 'mainly linear' feedback function H(R), Eq. (34) with H0 = 1; the
    # nonlinear branch is not required at this rebound magnitude (Sect. 7.1).
    hfb <- 1 + sstim_target_moderator1 * (target - bl_target)

    # Eq. (11): free ligand. Eliminated, bound, released.
    d/dt(central) <- -kel * central - kon * central * target + koff * complex

    # Eq. (13): free receptor. Synthesised at kin * F, eliminated, bound, released.
    d/dt(target) <- ksyn * moderator1 - kdeg * target -
      kon * central * target + koff * complex

    # Eq. (12): ligand-receptor complex. Formed by binding, lost by dissociation
    # and by elimination of the complex.
    d/dt(complex) <- kon * central * target - koff * complex - kint * complex

    # Eq. (14): the feedback moderator relaxes first-order toward H(R) at rate
    # alpha. No mass is exchanged with any other state.
    d/dt(moderator1) <- ktol * (hfb - moderator1)

    # Free drug is the paper's L; the free receptor R and its ratio to baseline
    # are the quantities Figs. 10 and 11 are drawn in.
    Cc <- central
    freeTarget <- target
    totalTarget <- target + complex
    targetRatio <- target / bl_target

    Cc ~ prop(propSd)
  })
}
