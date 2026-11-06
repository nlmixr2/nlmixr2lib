Pressly_2023_trastuzumab_tmdd <- function() {
  description <- paste(
    "Mechanistic TMDD model (Model 1 of a four-model tutorial hierarchy) for",
    "trastuzumab: the basic two-compartment target-mediated drug disposition",
    "model of Mager and Jusko. Free antibody in plasma binds a plasma target",
    "that is synthesised at zero order and degraded at first order; the",
    "antibody-target complex is internalised; free antibody distributes to a",
    "tissue compartment. The linear elimination rate of free antibody is zero",
    "as tabulated (Table S2), so elimination is target-mediated only; the",
    "authors' Figure 9a simulation used kel = 0.004 1/h. Parameter values are the",
    "authors' trastuzumab set taken mostly from Li et al. (2014), not fitted;",
    "the authors overlaid the simulation on digitised Tokuda (1999) phase-I",
    "data (1-8 mg/kg). Deterministic: no between-subject or residual",
    "variability. Siblings: modellib('Pressly_2023_trastuzumab_mpbpk'),",
    "modellib('Pressly_2023_trastuzumab_mpbpk_fcrn') and",
    "modellib('Pressly_2023_trastuzumab_mpbpk_shedding').",
    sep = " "
  )
  reference <- paste(
    "Pressly MA, Peletier LA, Zheng S, Sharma VD, Lien YTK, Wang W, Zhou H,",
    "Schmidt S. The quest for balance between capturing data and model",
    "complexity: A quantitative clinical pharmacology approach applied to",
    "monoclonal antibodies. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(5):639-652. doi:10.1002/psp4.12927.",
    "Equations 1-6 (main text); parameter values Supplementary Material S2",
    "Table S2; deposited MATLAB code Data S1 (tmdd1.m, tmdd1_data_impulse.m).",
    sep = " "
  )
  vignette <- "Pressly_2023_trastuzumab_tmdd_tutorial"

  # The ODEs are written in concentrations (nM) in the source. Here the two
  # antibody states hold amounts (nmol) so that a dose is given as an amount;
  # the target and complex states keep the source's plasma concentrations (nM).
  units <- list(time = "h", dosing = "nmol", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "trastuzumab (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    target = list(analyte = "target (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex = list(analyte = "trastuzumab-target complex", units = "nM", specimen = "plasma", verified = TRUE)
  )

  dosing <- "central"

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "HER2-positive metastatic breast cancer (comparison data); the model itself is a parameterised tutorial, not a fit",
    dose_range = "1, 2, 4 and 8 mg/kg single 90-min i.v. infusion (Tokuda 1999 comparison data); 650 nM i.v. bolus in the didactic simulations",
    weight_median = "70 kg (reference body weight used by the authors to convert mg/kg doses)",
    notes = paste(
      "TUTORIAL model. Parameter values are the authors' trastuzumab set",
      "(Supplementary Table S2, mostly from Li et al. 2014 CPT PSP), shared",
      "by all four models; nothing was estimated. The authors overlaid the",
      "simulations on mean concentrations digitised from Tokuda et al.",
      "1999 (18 Japanese patients with HER2-overexpressing metastatic",
      "breast cancer, single 1-8 mg/kg 90-min i.v. infusions). There is no",
      "inter-individual variability and no residual error; propSd is fixed(0)."
    )
  )

  ini({
    lvc <- fixed(log(3)); label("Plasma volume Vp (L)") # Table S2 'Vc' and 'Vp' = 3 L (plasma)
    lvp <- fixed(log(3)); label("Tissue volume Vt (L)") # Table S2 'Vp' = 3 L; tmdd1.m 'Vp=3' (kpt = Q/Vc, ktp = Q/Vp)
    lq <- fixed(log(300)); label("Plasma-tissue intercompartmental clearance Cldist (L/h)") # Table S2 'kpt = Q/Vc', 'ktp = Q/Vp' with Q = 300 L/h
    # Table S2 lists kel,p = 0 for the basic TMDD model, so free antibody
    # leaves plasma only through the target. The deposited Model 1 code
    # (tmdd1.m, tmdd1_data_impulse.m) uses kel = 0.004 1/h instead; that is
    # the value behind Figure 9a.
    kel <- fixed(0); label("First-order elimination rate constant of free antibody from plasma kel,p (1/h)") # Table S2 'kel,p' = 0

    lkon <- fixed(log(2.52)); label("Association rate constant of antibody and target kon,p (1/(nM*h))") # Table S2 'kon,p' = 2.52
    lkoff <- fixed(log(1.26)); label("Dissociation rate constant of the antibody-target complex koff,p (1/h)") # Table S2 'koff,p' = 1.26
    lksyn <- fixed(log(0.83)); label("Zero-order synthesis rate of plasma target ksyn,p (nM/h)") # Table S2 'ksyn,p' = 0.83
    lkdeg <- fixed(log(0.1)); label("First-order degradation rate constant of free target kdeg,p (1/h)") # Table S2 'kdeg,p' = 0.1
    lkint <- fixed(log(0.01)); label("Internalisation rate constant of the antibody-target complex kint,p (1/h)") # Table S2 'kint,p' = 0.01

    mw <- fixed(145532); label("Molecular weight of trastuzumab (g/mol)") # Table S2 'MW' = 145,532 Da

    # The source is a deterministic simulation with no residual error;
    # fixed(0) keeps the model simulable without inventing a value.
    propSd <- fixed(0); label("Proportional residual error (fraction; zero - none in the source)")
  })

  model({
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)

    # Equation 5: baseline target in the absence of antibody.
    target(0) <- ksyn / kdeg

    lig_p <- central / vc # free antibody in plasma, Lp (nM)
    lig_t <- peripheral1 / vp # free antibody in tissue, Lt (nM)

    # Equations 1-4 multiplied through by the plasma or tissue volume for the
    # two antibody states; kpt * Vp = ktp * Vt = Cldist.
    d/dt(central) <- vc * (-kon * lig_p * target + koff * complex - kel * lig_p) - q * lig_p + q * lig_t
    d/dt(peripheral1) <- q * lig_p - q * lig_t
    d/dt(target) <- ksyn - kon * lig_p * target + koff * complex - kdeg * target
    d/dt(complex) <- kon * lig_p * target - (kint + koff) * complex

    # Free trastuzumab in plasma: nM * g/mol = ng/L, / 1e6 = ug/mL.
    Cc <- lig_p * mw / 1e6
    Cc ~ prop(propSd)
  })
}
