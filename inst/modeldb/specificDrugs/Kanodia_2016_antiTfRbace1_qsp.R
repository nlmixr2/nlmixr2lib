Kanodia_2016_antiTfRbace1_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (cynomolgus monkey). Kanodia 2016 mechanistic",
    "pharmacokinetic-pharmacodynamic model for anti-transferrin-receptor",
    "(anti-TfR)/BACE1 bispecific antibodies delivered into the brain. Ten",
    "ODEs: a two-compartment serum PK with target-mediated disposition on",
    "soluble/membrane TfR in the central compartment; a brain-uptake module",
    "that reversibly binds antibody to TfR on the blood-brain-barrier",
    "endothelium, internalizes the antibody-TfR complex, dissociates and",
    "transcytoses free antibody into the brain cortex, plus a parallel",
    "non-specific uptake route; and an amyloid-beta (Abeta) turnover module",
    "whose synthesis is inhibited by brain-cortex antibody through a Hill",
    "term, giving the CSF-Abeta pharmacodynamic readout. Antibody-TfR",
    "association rate (kon) and affinity (KD) are fixed from in vitro",
    "measurements; the remaining parameters were fit by particle-swarm",
    "optimization. Deterministic: the source reports no IIV and no",
    "residual-error magnitudes. Three brain-module constants are not cleanly",
    "tabulated and are resolved from the deposited SimBiology model code; see",
    "the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Kanodia JS, Gadkar K, Bumbaca D, Zhang Y, Tong RK, Luk W, Hoyte K,",
    "Lu Y, Wildsmith KR, Couch JA, Watts RJ, Dennis MS, Ernst JA,",
    "Scearce-Levie K, Atwal JK, Ramanujan S, Joseph S (2016). Prospective",
    "Design of Anti-Transferrin Receptor Bispecific Antibodies for Optimal",
    "Delivery into the Human Brain. CPT Pharmacometrics Syst Pharmacol",
    "5(5):283-291. doi:10.1002/psp4.12081. PMCID: PMC4879477.",
    "Parameter values from Table 2 (final estimates) and Table 3 (anti-TfR",
    "affinities); model equations from the Supplementary Material; reaction",
    "topology and the brain-module constants not tabulated in Table 2 from",
    "the deposited SimBiology project (Supporting Information,",
    "TfRBACE_model_ver3p5.sbproj).",
    sep = " "
  )
  vignette <- "Kanodia_2016_antiTfRbace1_qsp"

  # Concentration-based QSP model, exactly as the published ODEs and the
  # deposited SimBiology code are written: every state is a concentration
  # (nM), volumes enter the rate laws explicitly as CL/Vc, CLd/Vc, CLd/Vp,
  # and the intravenous antibody dose is supplied as the central-compartment
  # concentration jump (nM). The vignette converts a mg/kg dose to nM via the
  # antibody molecular weight and Vc.
  units <- list(
    time = "day",
    dosing = "nM",
    concentration = "nM"
  )

  # No subject-level covariates. The only molecule-specific input is the
  # anti-TfR affinity KD (an in vitro property of the antibody, carried as an
  # ini parameter), not a population covariate.
  covariateData <- list()

  # Every state holds a CONCENTRATION (nM), not an amount; the model works in
  # a unit-volume frame with the physiologic volumes appearing explicitly in
  # the rate laws (see units above). analyte/specimen verified against the
  # Supplementary equations, Table 1 (list of model species) and the Methods.
  compartmentData <- list(
    central = list(analyte = "anti-TfR/BACE1 bispecific antibody", units = "nM", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "anti-TfR/BACE1 bispecific antibody",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    tfr_central = list(analyte = "transferrin receptor (free)", units = "nM", specimen = "serum", verified = TRUE),
    complex_central = list(analyte = "antibody-TfR complex", units = "nM", specimen = "serum", verified = TRUE),
    tfr_brain = list(analyte = "transferrin receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_brain = list(analyte = "antibody-TfR complex", units = "nM", specimen = "tissue", verified = TRUE),
    complex_brain_int = list(
      analyte = "antibody-TfR complex (internalized)",
      units = "nM",
      specimen = "endosome",
      verified = TRUE
    ),
    ab_brain_int = list(
      analyte = "anti-TfR/BACE1 bispecific antibody (internalized)",
      units = "nM",
      specimen = "endosome",
      verified = TRUE
    ),
    ab_cortex = list(
      analyte = "anti-TfR/BACE1 bispecific antibody",
      units = "nM",
      specimen = "brain ISF",
      verified = TRUE
    ),
    abeta = list(analyte = "amyloid-beta (Abeta)", units = "nM", specimen = "CSF", verified = TRUE)
  )

  # All ten states are paper-mechanistic TfR-trafficking / antibody-uptake /
  # Abeta species. central and peripheral1 use the canonical serum-PK names;
  # the remaining eight carry role-based lower-case snake_case names derived
  # from the Supplementary equation symbols (R_c -> tfr_central,
  # [Ab:R]_bi -> complex_brain_int, Ab_bc -> ab_cortex, and so on).
  paper_specific_compartments <- c(
    "tfr_central",
    "complex_central",
    "tfr_brain",
    "complex_brain",
    "complex_brain_int",
    "ab_brain_int",
    "ab_cortex",
    "abeta"
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 5,
    n_studies = 1,
    disease_state = "healthy (brain-delivery / amyloid-beta biomarker study)",
    dose_range = "intravenous bispecific antibody (calibration doses from Yu et al. 2014; human projections used 30 mg/kg i.v.)",
    notes = paste(
      "Model calibrated to published cynomolgus-monkey PK-PD data (Yu et al.",
      "2014, Sci Transl Med 6:261ra154) for two anti-TfR/BACE1 bispecifics",
      "(anti-TfR1, anti-TfR2) and a non-binding anti-gD antibody, with five",
      "animals per arm; validated against three further bispecifics",
      "(anti-TfR1/nBACE1, anti-TfR52/nBACE1, anti-TfR53/nBACE1). Parameters",
      "were fit by particle-swarm optimization, which yields no parameter",
      "precision, so no IIV or residual error is reported."
    )
  )

  ini({
    # =====================================================================
    # All parameters FIXED: this is a deterministic QSP model fit by
    # particle-swarm optimization, which gives no parameter precision; the
    # source reports no IIV and no residual-error magnitudes. Values are the
    # Table 2 final estimates (and Table 3 affinity) unless noted. The
    # deposited SimBiology project (ver3p5) carries DRAFT placeholder values
    # that disagree with Table 2 (e.g. kint = 1 vs 5.72, ktrans = 1 vs 200,
    # kon = 120 vs 17.9) plus an unused oral-absorption module, so it is used
    # only to confirm the reaction topology and to resolve the three
    # brain-module constants Table 2 does not cleanly report (see below).
    #
    # Units: time in days; concentrations in nM; volumes in L/kg; clearances
    # in L/day/kg. Only volume-normalized ratios (CL/Vc etc.) appear in the
    # rate laws, so the L/kg bases cancel.
    # =====================================================================

    # --- Serum two-compartment PK (log-transformed per library convention,
    #     though fixed here: the particle-swarm fit gives no precision) ---
    lcl <- fixed(log(0.0039)) ; label("Non-specific elimination clearance (L/day/kg)")  # Table 2 CL = 0.0039 L/day/kg
    lq <- fixed(log(0.0188)) ; label("Distribution clearance (L/day/kg)")               # Table 2 CLd = 0.0188 L/day/kg
    lvc <- fixed(log(0.0322)) ; label("Central compartment volume (L/kg)")              # Table 2 Vi central = 0.0322 L/kg
    lvp <- fixed(log(0.0289)) ; label("Peripheral compartment volume (L/kg)")           # Table 2 Vi peripheral = 0.0289 L/kg

    # --- Antibody-TfR binding (kon and KD fixed from in vitro; koff = KD*kon) ---
    kon <- fixed(17.9) ; label("Antibody-TfR association rate (1/(nM*day))")      # Table 2 kon = 17.9 1/nM/day, fixed from in vitro
    kd <- fixed(37) ; label("Antibody-TfR dissociation constant KD (nM)")         # Table 3 anti-TfR1 KD = 37 nM (default arm); other arms: TfR2 1650, TfR52 343, TfR53 143

    # --- Transferrin receptor turnover (central and brain) ---
    syn_tfr_c <- fixed(886) ; label("TfR synthesis rate, central (nM/day)")       # Table 2 SynC/bR central = 886 nM/day
    syn_tfr_b <- fixed(341) ; label("TfR synthesis rate, brain (nM/day)")         # Table 2 SynC/bR brain = 341 nM/day
    deg_tfr_c <- fixed(2.1) ; label("TfR degradation rate, central (1/day)")      # Table 2 DegC/bR central; printed '2-1 / 7' in the published PDF too; 2.1, 21 or 2^-1 = 0.5 all possible -- provisional 2.1 pending author query (Errata)
    deg_tfr_b <- fixed(7) ; label("TfR degradation rate, brain (1/day)")          # Table 2 DegC/bR brain = 7 1/day

    # --- Complex degradation ---
    deg_complex_c <- fixed(0.305) ; label("Antibody-TfR complex degradation, central (1/day)")  # Table 2 DegCAbR = 0.305 1/day
    deg_complex_b <- fixed(7) ; label("Antibody-TfR complex degradation, brain internalized (1/day)")  # not in Table 2; deposited code = tfrdeg*fdeg with fdeg = 1, so = deg_tfr_b = 7 (see Errata)

    # --- Brain uptake and transcytosis ---
    kint <- fixed(5.72) ; label("Antibody-TfR complex internalization rate (1/day)")  # Table 2 kint = 5.72 1/day
    kbuns <- fixed(0.016) ; label("Non-specific antibody uptake into brain (1/day)")  # Table 2 kbuns = 0.016 1/day
    ktrans <- fixed(200) ; label("Antibody transcytosis rate into brain cortex (1/day)")  # Table 2 ktrans = 200 1/day
    kbel <- fixed(291) ; label("Antibody clearance from brain cortex (1/day)")     # Table 2 DegCortexC = 291 1/day

    # --- Amyloid-beta (Abeta) turnover and BACE inhibition ---
    syn_abeta <- fixed(100) ; label("Abeta synthesis rate (nM/day)")              # Table 2 SynAbeta = 100 nM/day
    deg_abeta <- fixed(2.5) ; label("Abeta degradation rate (1/day)")             # Table 2 DegCortexAbeta = 2.5 1/day
    vm_bace <- fixed(0.75) ; label("Maximum fractional inhibition of Abeta synthesis by BACE engagement")  # Table 2 VmBACE = 0.75, fixed
    km_bace <- fixed(6.67) ; label("Brain-cortex antibody at half-maximal BACE inhibition (nM)")  # Table 2 KmBACE = 6.67 nM (4-fold lower for the higher-potency nBACE1 arm)
  })

  model({
    # Back-transform the serum-PK parameters (fixed; log scale per convention).
    cl <- exp(lcl)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)

    # Derived off rate (SimBiology initialAssignment koff = KD_tfr*kon).
    koff <- kd * kon

    # Steady-state (baseline) initial conditions. Table 1 lists the TfR
    # species initial condition as Syn/Deg and the Abeta initial condition as
    # Syn/Deg; every antibody / complex state starts at zero.
    tfr_c_base <- syn_tfr_c / deg_tfr_c
    tfr_b_base <- syn_tfr_b / deg_tfr_b
    abeta_base <- syn_abeta / deg_abeta

    tfr_central(0) <- tfr_c_base
    tfr_brain(0) <- tfr_b_base
    abeta(0) <- abeta_base

    # =====================================================================
    # ODE system (Supplementary Material equations). States are
    # concentrations (nM). The published equations carry volume ratios
    # Vb/Vc and Vc/VB on the brain<->central flux terms; the deposited
    # SimBiology model works in a unit-capacity compartment binding brain TfR
    # directly to the central concentration, i.e. these ratios are 1. No
    # numeric brain volume is reported anywhere, so the code's unit scaling
    # is used (Vb/Vc = Vc/VB = 1) -- see Errata.
    # =====================================================================

    # Systemic PK module -------------------------------------------------
    # Central antibody (Supp Eq d(Ab_c)/dt)
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * (central - peripheral1) -
      kon * central * tfr_central + koff * complex_central -
      kon * central * tfr_brain + koff * complex_brain -
      kbuns * central + kbel * ab_cortex
    # Peripheral antibody (Supp Eq d(Ab_p)/dt)
    d/dt(peripheral1) <- (q / vp) * (central - peripheral1)
    # Central free TfR (Supp Eq d(R_c)/dt)
    d/dt(tfr_central) <- syn_tfr_c - deg_tfr_c * tfr_central -
      kon * central * tfr_central + koff * complex_central
    # Central antibody-TfR complex (Supp Eq d[Ab:R]_c/dt)
    d/dt(complex_central) <- kon * central * tfr_central -
      koff * complex_central - deg_complex_c * complex_central

    # Brain-uptake module ------------------------------------------------
    # Brain (BBB) free TfR (Supp Eq d(R_b)/dt)
    d/dt(tfr_brain) <- syn_tfr_b - deg_tfr_b * tfr_brain -
      kon * ab_brain_int * tfr_brain + koff * complex_brain_int -
      kon * central * tfr_brain + koff * complex_brain
    # Brain luminal antibody-TfR complex (Supp Eq d[Ab:R]_b/dt)
    d/dt(complex_brain) <- kon * central * tfr_brain -
      koff * complex_brain - kint * complex_brain
    # Internalized antibody-TfR complex (Supp Eq d[Ab:R]_bi/dt)
    d/dt(complex_brain_int) <- kint * complex_brain -
      koff * complex_brain_int + kon * ab_brain_int * tfr_brain -
      deg_complex_b * complex_brain_int
    # Internalized free antibody in endothelium (Supp Eq d(Ab_bi)/dt)
    d/dt(ab_brain_int) <- kbuns * central + koff * complex_brain_int -
      kon * ab_brain_int * tfr_brain - ktrans * ab_brain_int
    # Transcytosed antibody in brain cortex (Supp Eq d(Ab_bc)/dt)
    d/dt(ab_cortex) <- ktrans * ab_brain_int - kbel * ab_cortex

    # Pharmacodynamic module ---------------------------------------------
    # Brain-cortex Abeta; synthesis inhibited by BACE engagement (Supp Eq dAbeta/dt)
    d/dt(abeta) <- syn_abeta *
      (1 - vm_bace * ab_cortex / (km_bace + ab_cortex)) -
      deg_abeta * abeta

    # Observations -------------------------------------------------------
    # Serum antibody concentration (the state is already a concentration).
    Cc <- central
    # CSF Abeta as a percentage of its drug-free baseline (SimBiology rule
    # abetaN = 100*abeta/Abeta_nominal); inhibition = 100 - AbetaPct.
    AbetaPct <- 100 * abeta / abeta_base
    # Deterministic simulation: no residual-error or IIV model, as the source
    # reports none (particle-swarm fit; no parameter precision).
  })
}
