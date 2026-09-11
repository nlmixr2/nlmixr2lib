Penney_2025_tce_tmdd_qsp <- function() {
  description <- "QSP. Mechanistic target-mediated drug disposition (TMDD) model predicting the pharmacokinetics of T-cell engagers (TCEs) as a function of target engagement. A two-compartment linear disposition backbone carries mass-action binding of free drug to two independent targets - the CD3 activating receptor on circulating T cells and a tumor-associated antigen (TAA) - each with its own synthesis, internalisation and drug-complex internalisation. Total elimination is the sum of an intrinsic (non-target) clearance plus the two target-mediated routes, which is what makes TCE half-lives much shorter than their IgG-like format would suggest and breaks the usual monkey-to-human allometric ranking. Trimer (CD3 + TAA simultaneous) formation is deliberately omitted by the authors. Defaults are the paper's generic HUMAN parameterisation in the CD3-binding-only configuration (bl_taa = 0), which is the configuration the authors themselves use for TCEs with no published TAA turnover data (DLL3, GPRC5D) and for Figure 2; set bl_taa / kint_taa / kd_taa to switch the TAA arm on, and thalf_intrinsic to 9 days for the cynomolgus monkey."
  reference <- paste(
    "Penney M, Ippolito A, Fevola E, Rata S, Brown L, Morentin Gutierrez P,",
    "Jones RDO. Predicting the Pharmacokinetics of T-Cell Engagers as a",
    "Function of Target-Mediated Drug Disposition.",
    "Clin Transl Sci. 2025;18(11):e70384. doi:10.1111/cts.70384.",
    "Structural siblings in this library: modellib('Betts_2019_pf_06671008_qsp')",
    "and modellib('Poels_2025_elranatamab_qsp').",
    sep = " "
  )
  vignette <- "Penney_2025_tce_tmdd"

  # Time is days and every concentration is nM, because the binding terms of
  # Equations 1-6 mix drug with receptor concentrations and the published
  # affinities (K_D) are quoted in nM. Every drug state holds a CONCENTRATION
  # (nM), exactly as printed in Equations 1-6 -- `central` is the paper's
  # D_Free and `peripheral1` is its D_P. The whole model is BODY-WEIGHT
  # NORMALISED, as the paper is (V_C = 40 mL/kg, V_P = 60 mL/kg, doses in
  # mg/kg): volumes are L/kg and a dose is given in nmol/kg, converted to the
  # central concentration by `f(central) <- 1 / vc`.
  #
  # The paper reports doses in mg/kg but never states a molecular weight for
  # the generic TCE, so mg/kg cannot be converted here. This does not affect
  # the paper's own quantitative anchors: at the 0.01 mg/kg dose of Figure 2
  # the system is in its dose-independent (linear, "fourth phase") regime, so
  # the predicted terminal half-life is identical for any TCE molecular weight
  # between 50 and 150 kDa. See the vignette for the demonstration.
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # Free circulating CD3 is carried as a dynamic turnover state here (Eq 2),
  # unlike its siblings Betts 2019 / Poels 2025 which derive CD3 algebraically
  # from a T-cell density. `target_cd3_central` is registered as a new canonical
  # in inst/references/compartment-names.md with this extraction as the founding
  # example. It joins the existing free-target family beside `target_bonemarrow`,
  # which keeps the free receptor pool typographically distinct from the
  # drug-bound `drug_cd3_central` dimer (Poels 2025) it binds to. The TAA arm
  # reuses the canonical generic TMDD pair `target` / `complex` because the
  # paper's TAA is deliberately unspecified (it is fitted per antigen).
  compartmentData <- list(
    central            = list(analyte = "T-cell engager (free)",            units = "nM", specimen = "serum",       verified = TRUE),
    peripheral1        = list(analyte = "T-cell engager (free)",            units = "nM", specimen = "serum",       verified = TRUE),
    target_cd3_central = list(analyte = "CD3 receptor (free)",              units = "nM", specimen = "whole blood", verified = TRUE),
    drug_cd3_central   = list(analyte = "T-cell engager-CD3 dimer",         units = "nM", specimen = "whole blood", verified = TRUE),
    target             = list(analyte = "tumor-associated antigen (free)",  units = "nM", specimen = "tumor",       verified = TRUE),
    complex            = list(analyte = "T-cell engager-TAA dimer",         units = "nM", specimen = "tumor",       verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "advanced solid and haematological malignancies (the indications of the surveyed clinical-stage TCEs)",
    dose_range = "0.01-10 mg/kg IV (the simulated range of Figure 2b); the surveyed TCEs span first-in-human microgram/kg doses to RP2Ds",
    notes = paste(
      "This is a PREDICTION method, not a model fitted to a subject-level dataset, so there is no",
      "subject count, no demographic table and no inter-individual variability. The generic",
      "disposition parameters (V_C = 40 mL/kg, V_P = 60 mL/kg, Q = 10 x CL, and CL set from an",
      "assumed intrinsic terminal half-life of 15 days in human / 9 days in the cynomolgus monkey)",
      "are assumptions stated in Sections 2.3, 2.6 and 2.7. The CD3 system parameters are generic",
      "literature values (Section 2.3). The method was evaluated against 29 published clinical-stage",
      "TCEs (Table 1); 18/22 cynomolgus-monkey and 16/18 human half-life predictions fell within",
      "two-fold (Table 2). Per-TCE K_D values to CD3 and to the TAA are tabulated in Supporting",
      "Information Tables S1 (cynomolgus monkey) and S2 (human)."
    )
  )

  ini({
    # ---- Generic linear (intrinsic) disposition -------------------------------
    # Body-weight-normalised throughout; the paper quotes mL/kg, encoded as L/kg.
    lvc <- fixed(log(0.040)); label("Central volume of distribution (L/kg)")                          # Section 2.3 / 2.7 (V_c = 40 mL/kg)
    lvp <- fixed(log(0.060)); label("Peripheral volume of distribution (L/kg)")                       # Section 2.3 / 2.7 (V_p = 60 mL/kg)
    # CL is not quoted directly; the paper sets it from the assumed intrinsic
    # terminal half-life as CL = ln(2) * (V_C + V_P) / half-life (Section 2.2).
    # 15 days is the human default (Section 2.7); use 9 days for the
    # cynomolgus monkey (Section 2.6).
    thalf_intrinsic <- fixed(15); label("Assumed intrinsic (non-target-mediated) terminal half-life (day)")  # Section 2.7
    qclratio <- fixed(10); label("Ratio of intercompartmental to intrinsic clearance, Q / CL (unitless)")    # Section 2.3 / 2.7 (Q = 10 x CL)

    # ---- Binding ---------------------------------------------------------------
    # A single association rate is assumed for BOTH targets, and each
    # dissociation rate is recovered as k_off = k_on * K_D (Section 2.3).
    # 1e6 /M/s * 1e-9 M/nM * 86400 s/day = 86.4 /nM/day.
    kon <- fixed(86.4); label("Association rate constant, shared by both targets (1/nM/day)")         # Section 2.3 (k_on = 1e6 /M/s)

    # ---- CD3 arm ---------------------------------------------------------------
    cd3_receptors <- fixed(50000); label("CD3 receptors per T cell (receptors/cell)")                 # Section 2.3
    tcell_blood <- fixed(1500); label("Circulating T-cell count (cells/uL)")                          # Section 2.3
    # Internalisation of CD3 and of the drug-CD3 dimer. The paper prints two
    # mutually inconsistent figures for this rate - "1.4%/minute, equating to a
    # half-life of about 36 min" - which cannot both hold (1.4%/min is a 49.5
    # min exponential half-life). The paper's OWN reported model output
    # adjudicates: only the 36-minute half-life reading reproduces the quoted
    # 4.4-day half-life at K_D,CD3 = 10 nM (Section 3.2). ln(2) / 36 min *
    # 1440 min/day = 27.726 /day. See the vignette Errata for the numerical
    # demonstration (the 1.4%/min reading gives 5.4 days).
    kint_cd3 <- fixed(27.726); label("CD3 and drug-CD3 dimer internalisation rate (1/day)")           # Section 2.3, resolved against the Section 3.2 output
    kd_cd3 <- fixed(10); label("Equilibrium dissociation constant of the TCE for CD3 (nM)")           # Section 3.2 (representative 'first generation' TCE affinity)

    # ---- TAA arm ---------------------------------------------------------------
    # OFF by default: bl_taa = 0 makes the whole TAA arm inert, which is the
    # paper's own configuration whenever no TAA turnover data exist (Section
    # 2.4: DLL3 and GPRC5D "predictions are made for CD3-binding only") and for
    # Figure 2. The authors state explicitly in the Table S2 caption that the
    # calibrated per-antigen TAA expression and turnover values are NOT
    # reported, so no non-zero default can be sourced; see the vignette Errata.
    bl_taa <- fixed(0); label("Baseline free tumor-associated-antigen concentration (nM); zero disables the TAA arm")  # Section 2.4 (values not reported by the authors)
    kint_taa <- fixed(1); label("TAA and drug-TAA dimer internalisation rate (1/day); placeholder, inert while bl_taa = 0")  # Section 2.4 (value not reported by the authors)
    kd_taa <- fixed(1); label("Equilibrium dissociation constant of the TCE for the TAA (nM); placeholder, inert while bl_taa = 0")  # per-TCE values in Supporting Information Tables S1 and S2

    # The source is a deterministic prediction method and reports no residual
    # error model and no inter-individual variability.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # Avogadro's number (1/mol), for the receptor-density to molar conversion.
    n_avogadro <- 6.02214076e23

    # ---- Derived system quantities ---------------------------------------------
    vc <- exp(lvc)
    vp <- exp(lvp)

    # Baseline free CD3: receptors/uL -> receptors/L -> mol/L -> nM.
    cd30 <- cd3_receptors * tcell_blood * 1e6 / n_avogadro * 1e9

    # Intrinsic clearance from the assigned terminal half-life (Section 2.2).
    cl <- log(2) * (vc + vp) / thalf_intrinsic
    q <- qclratio * cl

    # k_off = k_on * K_D for each target (Section 2.3).
    koff_cd3 <- kon * kd_cd3
    koff_taa <- kon * kd_taa

    # Zero-order receptor synthesis holding each free-target pool at its
    # baseline in the absence of drug (r_syn = k_int * baseline).
    ksyn_cd3 <- kint_cd3 * cd30
    ksyn_taa <- kint_taa * bl_taa

    target_cd3_central(0) <- cd30
    target(0) <- bl_taa

    # ---- Equations 1-6 ----------------------------------------------------------
    # Eq 1: free drug in the central compartment. Intrinsic clearance,
    # CD3 binding, distribution to the periphery, and TAA binding.
    d/dt(central) <-
      -(cl / vc) * central -
      (kon * target_cd3_central * central - koff_cd3 * drug_cd3_central) -
      (q / vc) * (central - peripheral1) -
      (kon * target * central - koff_taa * complex)

    # Eq 2: free CD3 - consumed by binding, replaced by synthesis, internalised.
    d/dt(target_cd3_central) <-
      -(kon * target_cd3_central * central - koff_cd3 * drug_cd3_central) +
      ksyn_cd3 - kint_cd3 * target_cd3_central

    # Eq 3: drug-CD3 dimer - formed by binding, lost by internalisation.
    d/dt(drug_cd3_central) <-
      (kon * target_cd3_central * central - koff_cd3 * drug_cd3_central) -
      kint_cd3 * drug_cd3_central

    # Eq 4: free TAA - the CD3 arm's structural mirror.
    d/dt(target) <-
      -(kon * target * central - koff_taa * complex) +
      ksyn_taa - kint_taa * target

    # Eq 5: drug-TAA dimer.
    d/dt(complex) <-
      (kon * target * central - koff_taa * complex) -
      kint_taa * complex

    # Eq 6: peripheral free drug.
    d/dt(peripheral1) <- q * (central - peripheral1) / vp

    # An IV dose is supplied in nmol/kg and enters as a central CONCENTRATION.
    f(central) <- 1 / vc

    # Cc is the free drug the paper's half-life is computed on (its D_Free);
    # CcTotal is the total-analyte assay equivalent in the central compartment.
    Cc <- central
    CcTotal <- central + drug_cd3_central + complex

    Cc ~ prop(propSd)
  })
}
