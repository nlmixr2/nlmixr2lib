Brenner_2016_pamrevlumab_rat <- function() {
  description <- "QSP. Preclinical (rat). Target-mediated drug disposition (TMDD) model for FG-3019 (pamrevlumab), a human anti-connective-tissue-growth-factor (CTGF) IgG1 monoclonal antibody, fit simultaneously to FG-3019, recombinant human CTGF and CTGF N-fragment kinetics in male Sprague-Dawley rats. Explicit binding of antibody to two constitutively produced target species -- intact CTGF (W) and its N-terminal half CTGF-N (N) -- each with a plasma and a tissue compartment, plus antibody-target complexes. Target-mediated elimination proceeds via tissue uptake of the antibody-CTGF complex, which reproduces the dose-dependent clearance of FG-3019 over 0.03-100 mg/kg. Typical-value mechanistic simulator: no IIV and no residual error are reported."
  reference <- "Brenner MC, Krzyzanski W, Chou JZ, Signore PE, Fung CK, Guzman D, Li D, Zhang W, Olsen DR, Nguyen VL, Koo CW, Sternlicht MD, Lipson KE. FG-3019, a Human Monoclonal Antibody Recognizing Connective Tissue Growth Factor, is Subject to Target-Mediated Drug Disposition. Pharm Res. 2016 Aug;33(8):1833-1849. doi:10.1007/s11095-016-1918-0. PMID 27059922. PMCID PMC4942499. Structural equations from the Electronic Supplementary Material (Kinetic Model development, Eqs. 1-20); parameter values from Table I. FG-3019 is the development code for the antibody later assigned the INN pamrevlumab."
  vignette <- "Brenner_2016_pamrevlumab_rat"

  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Two target species (intact CTGF and its N-terminal fragment), each free and
  # antibody-bound, each with a plasma and a tissue state. The canonical
  # target_/complex_ location suffixes cover a location only, so the species
  # token makes these paper-specific. elim_target / elim_nontarget are the
  # bookkeeping integrators of Supplement Eqs. (19)-(20) used to reproduce Fig. 7.
  paper_specific_compartments <- c(
    "target_ctgf",
    "target_ctgf_peripheral1",
    "target_ctgfn",
    "target_ctgfn_peripheral1",
    "complex_ctgf",
    "complex_ctgf_peripheral1",
    "complex_ctgfn",
    "complex_ctgfn_peripheral1",
    "elim_target",
    "elim_nontarget"
  )

  compartmentData <- list(
    central = list(analyte = "pamrevlumab (FG-3019), free", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pamrevlumab (FG-3019), free", units = "nmol", specimen = "tissue", verified = TRUE),
    target_ctgf = list(analyte = "CTGF (intact, W), free", units = "nmol", specimen = "plasma", verified = TRUE),
    target_ctgf_peripheral1 = list(analyte = "CTGF (intact, W)", units = "nmol", specimen = "tissue", verified = TRUE),
    target_ctgfn = list(
      analyte = "CTGF N-terminal fragment (N), free",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    target_ctgfn_peripheral1 = list(
      analyte = "CTGF N-terminal fragment (N)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_ctgf = list(
      analyte = "pamrevlumab-CTGF complex (AbW)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    complex_ctgf_peripheral1 = list(
      analyte = "pamrevlumab-CTGF complex (AbWT)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_ctgfn = list(
      analyte = "pamrevlumab-CTGF-N complex (AbN)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    complex_ctgfn_peripheral1 = list(
      analyte = "pamrevlumab-CTGF-N complex (AbNT)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    elim_target = list(
      analyte = "pamrevlumab (FG-3019) cleared by the target-mediated pathway, cumulative",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    elim_nontarget = list(
      analyte = "pamrevlumab (FG-3019) cleared by non-target-mediated pathways, cumulative",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 33,
    n_studies = 3,
    age_range = NA,
    weight_range = "average 336 g at dosing; model parameters expressed for a 0.3 kg animal",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state = "Healthy male Sprague-Dawley rats (no fibrotic disease model)",
    dose_range = "FG-3019 0.03, 0.3, 3, 10, 30 and 100 mg/kg IV bolus (n = 3-9 per dose); recombinant human CTGF 20 and 40 nmol/kg IV bolus (n = 3 per group); recombinant human CTGF N-fragment 20 and 40 nmol/kg IV bolus (n = 3 per group)",
    regions = NA,
    notes = paste0(
      "Cohort counts are the animals contributing to the three data sets fit ",
      "simultaneously (Brenner 2016 'Compartmental Pharmacokinetic Modeling' and ",
      "Results): FG-3019 single-dose PK, n = 3 + 3 + 9 + 6 + 6 + 6 = 33 across six ",
      "dose levels (Fig. 1 legend); rhCTGF and rhCTGF-N PK, 3 rats per group in four ",
      "groups (Fig. 3); endogenous CTGF-N response to FG-3019, 3 rats per group at ",
      "10, 30 and 100 mg/kg (Fig. 2, animals drawn from the same FG-3019 PK study). ",
      "Mean (not individual) data were fit by maximum likelihood in ADAPT 5, so the ",
      "CV% values in Table I are estimation precision, not between-animal variability. ",
      "The co-administration experiments of Figs. 4-5 (CTGF or RAP dosed after FG-3019) ",
      "were deliberately EXCLUDED from the fit -- the model treats FG-3019 as a ",
      "monovalent 75 kDa species, whereas co-dosed rhCTGF forms 2:1 complexes with the ",
      "bivalent 150 kDa antibody, and no rat RAP kinetics were available."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Plasma compartment shared by Ab, W, N, AbW and AbN (Supplement, text
    # following Eq. 10). All values are Brenner 2016 Table I; the CV% shown in
    # each comment is the estimation precision reported in that table.
    # ---------------------------------------------------------------------
    lvc <- log(0.01512); label("Plasma (central) volume shared by all plasma species (L)") # Table I, 'V, L' = 0.01512, CV 6.3%; 50.4 mL/kg for a 0.3 kg rat

    # --- Intact CTGF (W) disposition -------------------------------------
    clw <- 0.09243; label("Plasma clearance of intact CTGF (L/h)") # Table I, 'CLW, L/h' = 0.09243, CV 3.7%
    vwt <- 0.1160; label("Tissue volume for intact CTGF (L)") # Table I, 'VWT, L' = 0.1160, CV 13.0%
    clwt <- 1.380; label("Tissue clearance (elimination) of intact CTGF (L/h)") # Table I, 'CLWT, L/h' = 1.380, CV 9.5%
    cldw <- fixed(10); label("Distributional clearance of intact CTGF to tissue (L/h)") # Table I, 'CLdW, L/h' = 10, FIXED; set fast per Methods to force the very rapid initial phase of rhCTGF

    # --- CTGF N-fragment (N) disposition ---------------------------------
    vnt <- 0.04458; label("Tissue volume for the CTGF N-fragment (L)") # Table I, 'VNT, L' = 0.04458, CV 11.5%
    cldn <- 0.07998; label("Distributional clearance of the CTGF N-fragment to tissue (L/h)") # Table I, 'CLdN, L/h' = 0.07998, CV 9.4%

    # --- Antibody-target binding ------------------------------------------
    kdw <- 23.97; label("Equilibrium dissociation constant, FG-3019 to intact CTGF (nM)") # Table I, 'KDW, nM' = 23.97, CV 18.0%
    kdn <- 51.21; label("Equilibrium dissociation constant, FG-3019 to the CTGF N-fragment (nM)") # Table I, 'KDN, nM' = 51.21, CV 15.0%
    koffw <- fixed(100); label("Dissociation rate constant, FG-3019-CTGF complex (1/h)") # Table I, 'koffW, h-1' = 100, FIXED; enforces rapid-binding without the rapid-binding approximation
    koffn <- fixed(100); label("Dissociation rate constant, FG-3019-CTGF-N complex (1/h)") # Table I, 'koffN, h-1' = 100, FIXED; enforces rapid-binding without the rapid-binding approximation

    # --- Endogenous target baselines --------------------------------------
    bl_ctgf <- 0.02591; label("Baseline plasma concentration of intact CTGF (nM)") # Table I, 'CW0, nM' = 0.02591, CV 15.6%; 0.98 ng/mL, below the 5.6 ng/mL assay LLOQ
    bl_ctgfn <- 0.6572; label("Baseline plasma concentration of the CTGF N-fragment (nM)") # Table I, 'CN0, nM' = 0.6572, CV 10.6%; 12.5 ng/mL

    # --- FG-3019 (antibody) disposition -----------------------------------
    clab <- 0.0001321; label("Non-target-mediated plasma clearance of FG-3019 (L/h)") # Table I, 'CLAb, L/h' = 0.0001321, CV 6.4%; 0.45 mL/h/kg for a 0.3 kg rat
    vabt <- 0.01363; label("Tissue volume for FG-3019 (L)") # Table I, 'VAbT, L' = 0.01363, CV 12.5%
    cldab <- 0.0005692; label("Distributional clearance of FG-3019 to tissue (L/h)") # Table I, 'CLdAb, L/h' = 0.0005692, CV 27.7%
    cldabw <- fixed(10); label("Distributional clearance of the FG-3019-CTGF complex to tissue (L/h)") # Table I, 'CLdAbW, L/h' = 10, FIXED; set equal to CLdW per Methods
  })

  model({
    # -------------------------------------------------------------------
    # Constrained (secondary) parameters. Brenner 2016 Table I footnotes a-d
    # record the equality constraints imposed during the fit to keep the model
    # identifiable; they are derived here rather than duplicated as free
    # parameters so the constraints cannot drift apart.
    # -------------------------------------------------------------------
    vc <- exp(lvc)
    cln <- clw # Table I footnote a: CLN = CLW
    clabw <- clab # Table I footnote b: CLAb = CLAbW = CLAbN
    clabn <- clab # Table I footnote b: CLAb = CLAbW = CLAbN
    cldabn <- cldab # Table I footnote c: CLdAbN = CLdAb
    clabwt <- clwt # Table I footnote d: CLAbWT = CLWT
    vabwt <- vwt # Methods: 'CLdW = CLdAbW = 10 L/h with VAbWT = VWT'
    vabnt <- vabt # Methods: 'CLdAbN = CLdAb and VAbNT = VAbT'

    # Association rate constants, Supplement Eq. (15): kon = koff / KD (1/(nM*h))
    konw <- koffw / kdw
    konn <- koffn / kdn

    # Zero-order target production rates, Supplement Eqs. (11)-(12), written in
    # terms of the baseline plasma concentrations. Substituting Table I values
    # reproduces the tabulated secondary estimates (Table I footnote e):
    #   kw = 0.02591 * (0.09243 + 10 * 1.380 / 11.380) = 0.03382 nmol/h
    #   kn = 0.09243 * 0.6572                          = 0.06075 nmol/h
    kw <- bl_ctgf * (clw + cldw * clwt / (cldw + clwt)) # Supplement Eq. (11)
    kn <- cln * bl_ctgfn # Supplement Eq. (12)

    # -------------------------------------------------------------------
    # Free antibody plasma concentration drives every binding term. States are
    # AMOUNTS (nmol); dividing by the relevant volume gives nM.
    # -------------------------------------------------------------------
    cab <- central / vc

    # --- Intact CTGF, Supplement Eqs. (1)-(2) ---------------------------
    d/dt(target_ctgf) <- kw -
      (clw + cldw) / vc * target_ctgf -
      konw * cab * target_ctgf +
      koffw * complex_ctgf +
      cldw / vwt * target_ctgf_peripheral1
    d/dt(target_ctgf_peripheral1) <- cldw / vc * target_ctgf -
      (cldw + clwt) / vwt * target_ctgf_peripheral1

    # --- CTGF N-fragment, Supplement Eqs. (3)-(4) -----------------------
    # The N-fragment tissue compartment has no elimination term, so all
    # N-fragment loss is from plasma (renal filtration; Discussion).
    d/dt(target_ctgfn) <- kn -
      (cln + cldn) / vc * target_ctgfn -
      konn * cab * target_ctgfn +
      koffn * complex_ctgfn +
      cldn / vnt * target_ctgfn_peripheral1
    d/dt(target_ctgfn_peripheral1) <- cldn / vc * target_ctgfn -
      cldn / vnt * target_ctgfn_peripheral1

    # --- Free FG-3019, Supplement Eqs. (5)-(6) --------------------------
    d/dt(central) <- -(clab + cldab) / vc * central +
      cldab / vabt * peripheral1 -
      konw * cab * target_ctgf + koffw * complex_ctgf -
      konn * cab * target_ctgfn + koffn * complex_ctgfn
    d/dt(peripheral1) <- cldab / vc * central -
      cldab / vabt * peripheral1

    # --- FG-3019-CTGF complex, Supplement Eqs. (7)-(8) ------------------
    # This is the target-mediated elimination route: the complex distributes
    # rapidly to tissue (cldabw) and is eliminated there (clabwt).
    d/dt(complex_ctgf) <- konw * cab * target_ctgf -
      koffw * complex_ctgf -
      (clabw + cldabw) / vc * complex_ctgf +
      cldabw / vabwt * complex_ctgf_peripheral1
    d/dt(complex_ctgf_peripheral1) <- cldabw / vc * complex_ctgf -
      (cldabw + clabwt) / vabwt * complex_ctgf_peripheral1

    # --- FG-3019-CTGF-N complex, Supplement Eqs. (9)-(10) ---------------
    # The N-fragment complex follows free-antibody distribution kinetics and
    # has no tissue elimination, so binding FG-3019 rescues CTGF-N from renal
    # clearance and lets it accumulate in plasma.
    d/dt(complex_ctgfn) <- konn * cab * target_ctgfn -
      koffn * complex_ctgfn -
      (clabn + cldabn) / vc * complex_ctgfn +
      cldabn / vabnt * complex_ctgfn_peripheral1
    d/dt(complex_ctgfn_peripheral1) <- cldabn / vc * complex_ctgfn -
      cldabn / vabnt * complex_ctgfn_peripheral1

    # --- Elimination bookkeeping, Supplement Eqs. (17)-(20) -------------
    # Cumulative amount of FG-3019 (free and bound) removed by each pathway;
    # their ratio reproduces Fig. 7.
    d/dt(elim_target) <- clabwt / vabwt * complex_ctgf_peripheral1
    d/dt(elim_nontarget) <- clab / vc * central +
      clabw / vc * complex_ctgf +
      clabn / vc * complex_ctgfn

    # -------------------------------------------------------------------
    # Endogenous baselines, Supplement Eqs. (13), (14) and (16). Setting the
    # tissue states to their steady-state partners holds the untreated system
    # at baseline.
    # -------------------------------------------------------------------
    target_ctgf(0) <- bl_ctgf * vc # Supplement Eq. (16): W0 = CW0 * V
    target_ctgf_peripheral1(0) <- cldw * vwt * bl_ctgf / (cldw + clwt) # Supplement Eq. (13)
    target_ctgfn(0) <- bl_ctgfn * vc # Supplement Eq. (16): N0 = CN0 * V
    target_ctgfn_peripheral1(0) <- vnt * bl_ctgfn # Supplement Eq. (14)

    # -------------------------------------------------------------------
    # Observations, expressed as the assays measure them (nM).
    # -------------------------------------------------------------------
    # FG-3019 immunoassay reports total antibody: free plus both complexes,
    # in binding-site nM (MW 75 kDa = half of the 150 kDa bivalent IgG).
    Cc <- (central + complex_ctgf + complex_ctgfn) / vc
    # W-CTGF assay: forms with an intact hinge, i.e. intact CTGF free or bound.
    Cctgf_w <- (target_ctgf + complex_ctgf) / vc
    # N+W-CTGF assay: every Domain-2-containing form, i.e. intact CTGF and the
    # N-fragment, each free or antibody-bound.
    Cctgf_nw <- (target_ctgf + complex_ctgf + target_ctgfn + complex_ctgfn) / vc
    # Percent of cumulative FG-3019 elimination attributable to the
    # target-mediated pathway (Fig. 7).
    pctTargetMediated <- 100 * elim_target / (elim_target + elim_nontarget + 1e-30)
  })
}
