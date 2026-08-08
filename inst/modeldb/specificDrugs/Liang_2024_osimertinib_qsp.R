Liang_2024_osimertinib_qsp <- function() {
  description <- paste(
    "QSP. EGFRm+ target-engagement (receptor-occupancy) model for osimertinib",
    "(OSI, AZD9291) in NSCLC, describing irreversible covalent binding of OSI",
    "to the two EGFR mutants studied by Liang 2024 -- T790M and L858R -- with",
    "target turnover. Two independent free-target / drug-target-complex pairs",
    "share a single osimertinib concentration driver, reproducing the pair of",
    "mutant inhibition curves in Liang 2024 Figure 1. The osimertinib",
    "concentration is NOT fitted here: Liang 2024 generated it with a",
    "whole-body PK-Sim 10.0 PBPK model that is a platform port (no ODEs, no",
    "organ volumes and no blood flows are published, and no .pksim5 project",
    "was deposited), so that layer is not reproducible from the on-disk",
    "sources and is deliberately NOT extracted. The total osimertinib",
    "concentration is instead supplied per record as the canonical",
    "time-varying covariate CEFFECT and multiplied by the fraction unbound fu",
    "inside model(), exactly as in the authors' own supplementary code",
    "listing. Deterministic mechanism model: Liang 2024 simulated ten virtual",
    "subjects per scenario from PK-Sim population physiology rather than",
    "fitting IIV or residual error, so no etas and no error model are encoded.",
    "IMPORTANT DEVIATION -- the complex ODE carries an inferred",
    "'- kdeg * complex' elimination term that is absent from the paper's",
    "printed Eq 5 and from the supplementary listing. As printed, the complex",
    "has no elimination at all and (with koff = 0) occupancy rises",
    "monotonically to 100% under any sustained exposure, which provably",
    "contradicts the flat sawtooth plateau of the paper's own Figure 1.",
    "Restoring the term conserves total target at rbase and yields the",
    "steady-state occupancy kon*Cfree / (kon*Cfree + kdeg); it uses only",
    "published values (kdeg = 0.025 /h is declared in the supplement and used",
    "on the free-target line), but it is an inferred correction rather than a",
    "transcription. Extraction performed under operator sidecar decision",
    "oare_PMC10946252 request-001 = option B (answered 2026-08-05). See the",
    "vignette 'Assumptions and deviations' section for this and for the",
    "residual ~400-fold kon / concentration scale discrepancy against Figure",
    "1's plasma band.",
    sep = " "
  )
  reference <- paste(
    "Liang F, Zhang Y, Xue Q, Yao N. (2024). Exploring inter-ethnic and",
    "inter-patient variability and optimal dosing of osimertinib: a",
    "physiologically based pharmacokinetic modeling approach.",
    "Front Pharmacol 15:1363259. doi:10.3389/fphar.2024.1363259. PMC10946252.",
    "Target-engagement equations from Section 2.3 (Eqs 5-7) and Supplementary",
    "Data Sheet 1 (the authors' Berkeley-Madonna-style listing).",
    sep = " "
  )
  vignette <- "Liang_2024_osimertinib"

  units <- list(
    time          = "h",
    dosing        = "(not applicable; the osimertinib concentration is supplied as the time-varying covariate CEFFECT)",
    concentration = "umol/L"
  )

  # T790M and L858R are EGFR point mutants, not anatomical locations or
  # registered metabolites, so the canonical `target` / `complex` TMDD names
  # carry a mutant suffix and are declared here rather than silently warned on.
  paper_specific_compartments <- c(
    "target_t790m", "complex_t790m",
    "target_l858r", "complex_l858r"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Liang 2024 works entirely in molar concentrations
  # (Em0 = 0.299 umol/L), so the states are concentrations rather than
  # amounts. Specimen is "tissue": the paper's headline readout is pulmonary
  # EGFRm+ inhibition in the target tissue (Section 3.3), though the same
  # system is also driven with plasma concentrations in Figure 1.
  compartmentData <- list(
    target_t790m = list(
      analyte  = "EGFR T790M mutant (free, unbound)",
      units    = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_t790m = list(
      analyte  = "osimertinib-EGFR T790M covalent complex",
      units    = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    target_l858r = list(
      analyte  = "EGFR L858R mutant (free, unbound)",
      units    = "umol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_l858r = list(
      analyte  = "osimertinib-EGFR L858R covalent complex",
      units    = "umol/L",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CEFFECT = list(
      description        = "Time-varying TOTAL osimertinib concentration (umol/L) driving EGFRm+ target engagement",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "For this model the canonical effect-site PD driver CEFFECT carries",
        "the TOTAL (not free) osimertinib concentration in umol/L; the free",
        "concentration that actually drives binding is computed inside",
        "model() as CEFFECT * fu. This follows the authors' own supplementary",
        "listing, whose binding term is 'kon * Cc * fu * Tfree' with Cc the",
        "total concentration. Liang 2024 Section 2.3 defines the driver of",
        "printed Eq 5 as Clung, 'the free OSI concentration in the lung';",
        "Clung therefore corresponds to CEFFECT * fu here. Supplying a total",
        "concentration is the more useful convention because that is what a",
        "PK model returns.",
        "In the source study CEFFECT is produced by a whole-body PK-Sim 10.0",
        "PBPK model that is not reproducible from the published sources and",
        "is not extracted (see description). Liang 2024 drove the system with",
        "both plasma and pulmonary concentrations; the paper's Table 1 gives a",
        "lung-to-plasma partition coefficient Klu,p of 28.5, so a pulmonary",
        "run uses a correspondingly higher CEFFECT than a plasma run.",
        "Published exposure anchors for an 80 mg once-daily regimen in NSCLC",
        "patients (Liang 2024 Tables 2 and 4) are Css,max approximately 0.577",
        "umol/L and Ctrough approximately 0.342 umol/L in Caucasians, with an",
        "efficacy/safety Ctrough window of 0.328-0.677 umol/L.",
        "For downstream simulations users may supply CEFFECT from any",
        "osimertinib PK model; note that Liang 2024 did NOT couple this",
        "target layer to any published popPK model, so pairing it with",
        "modellib('Brown_2017_osimertinib') is a user choice, not the",
        "authors' design."
      ),
      source_name        = "(none; Clung is computed by the PK-Sim PBPK model and is not a named NONMEM data column)"
    )
  )

  population <- list(
    species       = "human (in silico; PK-Sim virtual populations)",
    n_subjects    = 10L,
    n_studies     = 1L,
    age_range     = "44-83 years for the Caucasian NSCLC scenario; 53-78 years Japanese; 35-76 years Chinese (Liang 2024 Table 2)",
    sex_female_pct = 71,
    race_ethnicity = c(Caucasian = NA_real_, Japanese = NA_real_, Chinese = NA_real_),
    disease_state = "non-small cell lung cancer (NSCLC) with EGFR T790M and/or L858R activating mutations",
    dose_range    = "osimertinib 80 mg once daily for 14 consecutive days in the Figure 1 target-engagement simulation; 20-240 mg once daily across the wider validation set (Liang 2024 Table 2)",
    regions       = "simulated Caucasian, Japanese and Chinese populations",
    notes         = paste(
      "Liang 2024 Section 2.3: 'Each simulation consisted of ten virtual",
      "subjects.' The virtual populations' demographic characteristics and",
      "dosing regimens were taken from the clinical studies tabulated in",
      "Liang 2024 Table 2 (Planchard 2016, Vishwanathan 2018a/2019, Harvey",
      "2018, Grande 2019, Fujiwara 2023, Zhao 2018). The three ethnic groups",
      "differ in the PBPK layer only (liver volume and CYP abundances,",
      "Liang 2024 Table 1); the target-engagement parameters extracted here",
      "are shared across all three, so ethnicity enters this model solely",
      "through the CEFFECT concentration the user supplies.",
      "The target-side physiology is that of the diseased (NSCLC) population:",
      "fu = 0.019 rather than the healthy-volunteer 0.013 (Liang 2024",
      "Table 1)."
    )
  )

  ini({
    # ---- Association rate constants, one per EGFR mutant ------------------
    # Liang 2024 Section 2.3: "The kon values of 0.91 (binding to T790M) and
    # 0.44 (binding to L858R) uM^-1 s^-1 were obtained from the paper (Zhai
    # et al., 2020), equivalent to the ratio of Kinact/Ki."
    # Unit conversion to the model time base of hours (the paper's own
    # supplement, kturnover and Figure 1 all use hours):
    #   0.91 /uM/s * 3600 s/h = 3276 /uM/h
    #   0.44 /uM/s * 3600 s/h = 1584 /uM/h
    # Fixed: literature-inherited from Zhai 2020, not estimated by Liang 2024.
    lkon_t790m <- fixed(log(3276))
    label("Log association rate constant of osimertinib to EGFR T790M (1/(umol/L)/h)")
    # Liang 2024 Section 2.3: kon = 0.91 /uM/s (Zhai 2020) = 3276 /uM/h.

    lkon_l858r <- fixed(log(1584))
    label("Log association rate constant of osimertinib to EGFR L858R (1/(umol/L)/h)")
    # Liang 2024 Section 2.3: kon = 0.44 /uM/s (Zhai 2020) = 1584 /uM/h.

    # ---- Dissociation rate constant --------------------------------------
    # Liang 2024 Section 2.3: "koff was assumed to be 0 due to irreversible
    # covalent binding to EGFR for OSI." Supplementary Data Sheet 1 likewise
    # sets koff = 0.00. Kept as an explicit fixed structural zero (rather
    # than dropped from the equations) so the reversible-binding form of
    # Eqs 5-6 stays visible and a user can relax it.
    koff <- fixed(0)
    label("Dissociation rate constant of the osimertinib-EGFRm+ complex (1/h; zero for irreversible covalent binding)")
    # Liang 2024 Section 2.3 and Supplementary Data Sheet 1 (koff = 0.00).

    # ---- Target turnover --------------------------------------------------
    # Liang 2024 Section 2.3: "kturnover was obtained to be 0.025 h-1 from
    # the paper (Greig et al., 2015)." Supplementary Data Sheet 1 declares
    # Kturnover = 0.025 AND, separately, kdeg = 0.025, using kdeg only on the
    # free-target line. Both are fixed literature values, not estimated.
    lkturnover <- fixed(log(0.025))
    label("Log re-synthesis (turnover) rate constant of EGFRm+ (1/h)")
    # Liang 2024 Section 2.3 (Greig 2015); Supplementary Data Sheet 1 Kturnover = 0.025.

    lkdeg <- fixed(log(0.025))
    label("Log degradation rate constant of EGFRm+ and of the osimertinib-EGFRm+ complex (1/h)")
    # Supplementary Data Sheet 1 kdeg = 0.025 (free-target line). Its use on
    # the COMPLEX line is the inferred correction described in `description`
    # and in the vignette 'Assumptions and deviations'.

    # ---- Baseline target concentration ------------------------------------
    # Liang 2024 Section 2.3: "Em0 is starting concentration of EGFRm+, and
    # was set 0.299 uM based on the paper (Bartelink et al., 2022)."
    # Supplementary Data Sheet 1 repeated-dose block: T0 = 0.299, INIT
    # Tfree = 0.299. Fixed literature value.
    lrbase <- fixed(log(0.299))
    label("Log baseline (total) EGFRm+ concentration Em0 (umol/L)")
    # Liang 2024 Section 2.3 (Bartelink 2022); Supplementary Data Sheet 1 T0 = 0.299.

    # ---- Fraction unbound --------------------------------------------------
    # Liang 2024 Table 1: fup = 0.013 in the healthy population and 0.019 in
    # the diseased (NSCLC) population, the latter calculated from the patient
    # albumin level using the paper's Eq 1. The target-engagement simulations
    # are NSCLC-patient simulations, so 0.019 is the applicable value; it is
    # also the value hard-coded in Supplementary Data Sheet 1 (fu = 0.019).
    fu <- fixed(0.019)
    label("Fraction of osimertinib unbound in plasma in NSCLC patients (unitless)")
    # Liang 2024 Table 1 (diseased column); Supplementary Data Sheet 1 fu = 0.019.
  })

  model({
    # 1. Back-transform the fixed structural parameters.
    kon_t790m <- exp(lkon_t790m)
    kon_l858r <- exp(lkon_l858r)
    kturnover <- exp(lkturnover)
    kdeg      <- exp(lkdeg)
    rbase     <- exp(lrbase)

    # 2. Free osimertinib concentration driving covalent binding. CEFFECT is
    #    the TOTAL concentration supplied per record; Liang 2024's printed
    #    Eq 5 driver Clung ("the free OSI concentration in the lung") is this
    #    product. Matches the authors' supplementary listing term
    #    'kon * Cc * fu * Tfree'.
    cfree <- CEFFECT * fu

    # 3. Target-engagement ODEs, one independent pair per EGFR mutant.
    #
    #    Liang 2024 printed Eqs 5-6 (recovered from the PDF; the trimmed
    #    markdown dropped every display equation as `formula-not-decoded`):
    #       dOEm/dt    = kon * Clung * Emfree - koff * OEm
    #       dEmfree/dt = (Em0 - Emfree) * kturnover - kon * Clung * Emfree
    #                    + koff * OEm
    #
    #    The free-target line is transcribed verbatim. The complex line adds
    #    '- kdeg * complex', which is NOT in the paper. Rationale, in full,
    #    in `description` and the vignette; in brief: as printed the complex
    #    has no elimination, so with koff = 0 the complex is monotone
    #    non-decreasing while free target is bounded above by rbase, forcing
    #    occupancy to 100% under any sustained exposure. That contradicts
    #    Figure 1's flat sawtooth plateau. With the term restored, total
    #    target is conserved exactly at rbase (because kturnover = kdeg) and
    #    steady-state occupancy is kon*cfree / (kon*cfree + kdeg).
    #    Operator sidecar oare_PMC10946252 request-001, option B.
    d/dt(target_t790m) <-
      (rbase - target_t790m) * kturnover -
      kon_t790m * cfree * target_t790m +
      koff * complex_t790m
    d/dt(complex_t790m) <-
      kon_t790m * cfree * target_t790m -
      koff * complex_t790m -
      kdeg * complex_t790m

    d/dt(target_l858r) <-
      (rbase - target_l858r) * kturnover -
      kon_l858r * cfree * target_l858r +
      koff * complex_l858r
    d/dt(complex_l858r) <-
      kon_l858r * cfree * target_l858r -
      koff * complex_l858r -
      kdeg * complex_l858r

    # 4. Initial conditions. Supplementary Data Sheet 1 repeated-dose block:
    #    INIT Tfree = 0.299 (= T0), INIT TC = 0. The system starts drug-free
    #    with all EGFRm+ unbound.
    target_t790m(0)  <- rbase
    complex_t790m(0) <- 0
    target_l858r(0)  <- rbase
    complex_l858r(0) <- 0

    # 5. Outputs: per-mutant EGFRm+ inhibition (%), i.e. fractional receptor
    #    occupancy. Liang 2024 Eq 7:
    #       Inhibition (%) = OEm / (Emfree + OEm) * 100
    #    Deterministic mechanism model -- Liang 2024 reports no residual error
    #    and no IIV for this layer, so none is encoded (see description).
    inhT790M <- 100 * complex_t790m / (target_t790m + complex_t790m)
    inhL858R <- 100 * complex_l858r / (target_l858r + complex_l858r)

    # Total target concentration per mutant; a conserved quantity (= rbase)
    # that the vignette uses as a mass-balance check.
    totalT790M <- target_t790m + complex_t790m
    totalL858R <- target_l858r + complex_l858r
  })
}
