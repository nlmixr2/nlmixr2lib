MethvenJeknic_2026_attrcm_fpe_invitro <- function() {
  description <- "In vitro (human serum fluorescent probe exclusion assay, 25 C). Mass-action model of the fluorescent probe exclusion (FPE) assay used to estimate acoramidis and tafamidis binding to transthyretin (TTR) and albumin: the stabiliser binds the two tetramer thyroxine pockets with negative cooperativity and four equivalent albumin sites, while a covalent stilbene probe modifies only stabiliser-free tetramer (one bound molecule excludes the probe from the second pocket); fluorescence is a linear function of the fraction of probe that has reacted. Tetramer dissociation is neglected. Deterministic model without residual error."
  reference <- "Methven Jeknic S, Ahn-Horst TA, Ji AX, Ju JH, Paull M, Ji SG, Sinha U. Evaluation of the Current ATTR-CM Treatment Landscape via a Mathematical Model of TTR Dissociation and Amyloid Formation. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70277. doi:10.1002/psp4.70277. Equations from Supplemental Equations 1-7, initial conditions from Table S1, dissociation constants from Table 1, probe rate and readout slope from Table S5; probe depletion and the off-rate convention from the authors' code deposit https://github.com/BridgeBioAnalytics/ttr-mechanistic-model (commit 6686acd, ttr/fpe.py)."
  vignette <- "MethvenJeknic_2026_attrcm"
  # Assay species of Table S1: TTR tetramer free, singly and doubly
  # stabiliser-bound and covalently probe-modified, the unreacted probe, and
  # free / stabiliser-occupied albumin sites. Free stabiliser is held in
  # `central` so a bolus there sets the nominal drug concentration.
  paper_specific_compartments <- c(
    "ttr_tetramer",
    "ttr_tetramer_drug1",
    "ttr_tetramer_drug2",
    "ttr_tetramer_probe",
    "probe",
    "albumin",
    "complex_albumin"
  )
  units <- list(time = "h", dosing = "nM", concentration = "nM")

  covariateData <- list(
    TRT_ACORAMIDIS = list(
      description = "Acoramidis indicator; selects the acoramidis TTR- and albumin-binding constants for the stabiliser in `central`",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no acoramidis)",
      notes = "Mutually exclusive with TRT_TAFAMIDIS. The DMSO control is any drug indicator with no drug added.",
      source_name = "drug = 'acoramidis'"
    ),
    TRT_TAFAMIDIS = list(
      description = "Tafamidis indicator; selects the tafamidis TTR- and albumin-binding constants for the stabiliser in `central`",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no tafamidis)",
      notes = "Mutually exclusive with TRT_ACORAMIDIS.",
      source_name = "drug = 'tafamidis'"
    )
  )

  compartmentData <- list(
    ttr_tetramer = list(
      analyte = "transthyretin tetramer, unbound by stabiliser or probe",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    ttr_tetramer_drug1 = list(
      analyte = "transthyretin tetramer bound by one stabiliser molecule",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    ttr_tetramer_drug2 = list(
      analyte = "transthyretin tetramer bound by two stabiliser molecules",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    ttr_tetramer_probe = list(
      analyte = "transthyretin tetramer covalently modified by the fluorescent probe",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    probe = list(analyte = "unreacted covalent stilbene probe", units = "nM", specimen = "serum", verified = TRUE),
    central = list(
      analyte = "free stabiliser (acoramidis or tafamidis)",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    albumin = list(
      analyte = "free albumin drug-binding sites (four per albumin molecule)",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    complex_albumin = list(
      analyte = "stabiliser bound to an albumin binding site",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (pooled human serum)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Not applicable (in vitro assay)",
    dose_range = "Acoramidis 0, 1.25, 2.5, 5 and 10 uM; tafamidis 0, 1.25, 2.5, 5, 10, 20 and 40 uM, added with 3.6 uM probe to serum (TTR about 4.9 uM after dilution)",
    regions = "Not applicable",
    notes = "Fitted (LMFIT, Powell) to the human-serum FPE time courses of Miller et al. 2018 (J Med Chem 61:7862), plus the free fraction of tafamidis with albumin (Tess 2023), the 96.4% plasma-protein binding of acoramidis (prescribing information), and soft ITC targets on the TTR dissociation constants (Bulawa 2012, Penchala 2013) (Table S3). The albumin concentration (650 uM, four equivalent sites) is assumed; TTR (5 uM) and the 0.98 dilution are from the assay report."
  )

  ini({
    # ---- Association rates (Table S5, Assumed) ----
    kon_ttr <- fixed(3.6e10); label("Stabiliser association with either tetramer pocket, kaD1 = kaD2, both drugs (1/M/h)") # Table S5 kaD,1 and kaD,2 = 3.6e10 (Assumed)
    kon_alb <- fixed(3.6e7); label("Stabiliser association with an albumin site, kaDA, both drugs (1/M/h)") # Table S5 kaD,A = 3.6e7 (Assumed)

    # ---- Dissociation constants at 25 C (Table 1) ----
    lkd1_acoramidis <- log(6.00e-9); label("Acoramidis dissociation constant, first TTR binding event, Kd1 (M)") # Table 1 acoramidis Kd1 = 6.00 nM (4.18-8.51)
    lkd2_acoramidis <- log(190e-9); label("Acoramidis dissociation constant, second TTR binding event, Kd2 (M)") # Table 1 acoramidis Kd2 = 190 nM (139-322)
    lkd_alb_acoramidis <- log(250e-6); label("Acoramidis dissociation constant per albumin site, Kd,albumin (M)") # Table 1 acoramidis Kd,albumin = 250 uM (158-410)
    lkd1_tafamidis <- log(4.33e-9); label("Tafamidis dissociation constant, first TTR binding event, Kd1 (M)") # Table 1 tafamidis Kd1 = 4.33 nM (3.94-4.78)
    lkd2_tafamidis <- log(248e-9); label("Tafamidis dissociation constant, second TTR binding event, Kd2 (M)") # Table 1 tafamidis Kd2 = 248 nM (184-291)
    lkd_alb_tafamidis <- log(3.05e-6); label("Tafamidis dissociation constant per albumin site, Kd,albumin (M)") # Table 1 tafamidis Kd,albumin = 3.05 uM (2.94-3.16)

    # ---- Probe and optical readout (Table S5) ----
    lki <- log(1.926e5); label("Covalent probe modification of free tetramer, ki (1/M/h)") # Table S5 ki = 1.926e5 (Fit)
    s_fl <- 1.0589; label("Linear slope from fraction of probe reacted to fraction of total fluorescence, s (unitless)") # Table S5 s = 1.0589 (Fit)
    b_fl <- fixed(0); label("Linear intercept of the fluorescence readout, b (fraction)") # Table S5 b = 0 (Assumed)

    # ---- Assay initial conditions (Table S1) ----
    tet0 <- fixed(4.9); label("Initial TTR tetramer (uM)") # Table S1 T = 4.9 uM (5 uM x 0.98 dilution)
    alb0 <- fixed(2548); label("Albumin drug-binding sites, four per albumin molecule (uM)") # Table S1 A = 2548 uM (650 uM x 0.98 x 4)
    probe0 <- fixed(3.6); label("Initial covalent probe (uM)") # Table S1 P = 3.6 uM
  })
  model({
    # The indicator picks the drug's dissociation constants; with neither
    # set the drug in `central` does not bind.
    on_stab <- TRT_ACORAMIDIS + TRT_TAFAMIDIS
    kd1 <- TRT_ACORAMIDIS * exp(lkd1_acoramidis) + TRT_TAFAMIDIS * exp(lkd1_tafamidis)
    kd2 <- TRT_ACORAMIDIS * exp(lkd2_acoramidis) + TRT_TAFAMIDIS * exp(lkd2_tafamidis)
    kda <- TRT_ACORAMIDIS * exp(lkd_alb_acoramidis) + TRT_TAFAMIDIS * exp(lkd_alb_tafamidis)

    # Off-rates from the macroscopic Kds (code ttr/fpe.py): two equivalent
    # empty pockets give kdD1 = 2 * kon * Kd1 and the doubly bound tetramer
    # releases from either pocket, kdD2 = kon * Kd2 / 2. Rates in molar
    # units are converted to nM (second order 1/M/h x 1e-9).
    kon <- kon_ttr * 1e-9 * on_stab
    kona <- kon_alb * 1e-9 * on_stab
    koff1 <- 2 * kon_ttr * kd1
    koff2 <- 0.5 * kon_ttr * kd2
    koffa <- kon_alb * kda
    ki <- exp(lki) * 1e-9

    # Supplemental Equations 1-7 (nM). Eq 2 as printed omits the factor 2 on
    # kdD2 [T:D2] that Eqs 3 and 7 carry; the code (and the physiological
    # Eq 9) includes it, which conserves drug. The probe is consumed by the
    # covalent reaction (code: dP/dt = -ki [T][P]); the paper lists P in
    # Table S1 without printing its equation.
    d/dt(ttr_tetramer) <- koff1 * ttr_tetramer_drug1 - (2 * kon * central + ki * probe) * ttr_tetramer
    d/dt(ttr_tetramer_drug1) <- 2 * kon * ttr_tetramer * central + 2 * koff2 * ttr_tetramer_drug2 -
      (koff1 + kon * central) * ttr_tetramer_drug1
    d/dt(ttr_tetramer_drug2) <- kon * ttr_tetramer_drug1 * central - 2 * koff2 * ttr_tetramer_drug2
    d/dt(ttr_tetramer_probe) <- ki * ttr_tetramer * probe
    d/dt(probe) <- -ki * ttr_tetramer * probe
    d/dt(albumin) <- koffa * complex_albumin - kona * albumin * central
    d/dt(complex_albumin) <- kona * albumin * central - koffa * complex_albumin
    d/dt(central) <- koff1 * ttr_tetramer_drug1 + 2 * koff2 * ttr_tetramer_drug2 + koffa * complex_albumin -
      (2 * kon * ttr_tetramer + kon * ttr_tetramer_drug1 + kona * albumin) * central

    # Table S1 (uM -> nM); drug and probe are added together at time 0, so
    # no binding has equilibrated beforehand.
    ttr_tetramer(0) <- 1000 * tet0
    albumin(0) <- 1000 * alb0
    probe(0) <- 1000 * probe0

    # Fluorescence as a fraction of total, F(t) = s * [T:P] / ([T:P] + [P]) + b.
    fluor <- s_fl * ttr_tetramer_probe / (ttr_tetramer_probe + probe) + b_fl
  })
}
