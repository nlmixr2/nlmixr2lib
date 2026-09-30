Abbiati_2021_avadomide_gbm <- function() {
  description <- paste(
    "QSP. Avadomide (CC-122)-induced neutropenia in adults with glioblastoma",
    "(GBM); GBM-cohort parameter set, the fit that also estimated the PD",
    "parameters EC50 and n.",
    "A two-compartment oral PK model with first-order absorption and a lag",
    "time drives a sigmoid Emax partial block of neutrophil maturation",
    "(cereblon-mediated Ikaros degradation), which acts on a Michaelis-Menten",
    "transfer from the second to the third bone-marrow maturation stage of a",
    "neutrophil life-cycle model: proliferating precursors -> three transit",
    "stages -> bone-marrow reservoir of mature neutrophils -> circulating",
    "neutrophils (ANC). Two power-law feedbacks regulate proliferation (on",
    "transit-2 level) and marrow egress (on ANC). Every rate constant is",
    "back-calculated from homeostasis. Deterministic mechanism model: the",
    "authors built virtual patients by resampling the empirical distributions",
    "of baseline ANC, gamma, the reservoir-to-circulation ratio and the KM",
    "fraction rather than by estimating IIV or residual error, so no etas and",
    "no error model are encoded. Sibling files carry the diffuse large",
    "B-cell lymphoma (Abbiati_2021_avadomide_dlbcl) and multiple myeloma",
    "(Abbiati_2021_avadomide_mm) parameter sets of the same structure.",
    sep = " "
  )
  reference <- paste(
    "Abbiati RA, Pourdehnad M, Carrancio S, Pierce DW, Kasibhatla S,",
    "McConnell M, Trotter MWB, Loos R, Santini CC, Ratushny AV. Quantitative",
    "Systems Pharmacology Modeling of Avadomide-Induced Neutropenia Enables",
    "Virtual Clinical Dose and Schedule Finding Studies. AAPS J.",
    "2021;23(5):103. doi:10.1208/s12248-021-00623-8.",
    "Correction: AAPS J. 2022;24(1):29. doi:10.1208/s12248-021-00673-y",
    "(replaces the supplementary material; no value in the main text changed).",
    sep = " "
  )
  vignette <- "Abbiati_2021_avadomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  # The bone-marrow reservoir of mature neutrophils awaiting egress (the
  # paper's 'Reserv') has no canonical compartment name. The proliferating
  # pool, the transit stages and the circulating pool use the registered
  # Friberg-chain names `prol`, `transit1`..`transit3` and `circ`.
  paper_specific_compartments <- c("reservoir")

  compartmentData <- list(
    depot = list(analyte = "avadomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "avadomide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avadomide", units = "mg", specimen = "tissue", verified = TRUE),
    prol = list(
      analyte = "committed proliferative neutrophil precursors",
      units = "10^9 cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    transit1 = list(analyte = "maturing neutrophils", units = "10^9 cells/L", specimen = "tissue", verified = TRUE),
    transit2 = list(analyte = "maturing neutrophils", units = "10^9 cells/L", specimen = "tissue", verified = TRUE),
    transit3 = list(analyte = "maturing neutrophils", units = "10^9 cells/L", specimen = "tissue", verified = TRUE),
    reservoir = list(
      analyte = "mature neutrophils (bone-marrow reservoir)",
      units = "10^9 cells/L",
      specimen = "tissue",
      verified = TRUE
    ),
    circ = list(analyte = "neutrophils", units = "10^9 cells/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    disease_state = "glioblastoma (GBM)",
    dose_range = paste(
      "Avadomide 3-6 mg orally once daily, continuously (QD) or on a",
      "5-days-on / 2-days-off (5/7) schedule (Figure 2: 4 mg QD, n = 4;",
      "5 mg 5/7, n = 5; 5 mg QD, n = 3; 6 mg 5/7, n = 3; 3 mg QD, n = 22)",
      sep = " "
    ),
    regions = "Multinational (trial NCT01421524, CC-122-ST-001)",
    notes = paste(
      "Absolute neutrophil counts from the first treatment cycle of the GBM",
      "cohorts of the first-in-human avadomide trial NCT01421524; ANC after",
      "the first G-CSF administration were removed, and individual ANC were",
      "resampled into 4-day windows before fitting (Supplementary Materials,",
      "Figure 3 data processing). GBM patients had not received prior",
      "marrow-depleting therapy and were taken as the closest match to a",
      "healthy marrow, so this cohort was fitted first: five parameters (EC50,",
      "n, gamma, the reservoir ratio and the KM fraction) were regressed",
      "simultaneously on all GBM dose groups (Results; Figure 3a). No",
      "demographic summary (age, sex, weight) is reported.",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Avadomide PK -- Supplementary Materials 'Avadomide PK' parameter table
    # (two-compartment, first-order absorption with lag time, first-order
    # elimination). The paper uses these values as inputs; none was
    # estimated in this analysis.
    # -----------------------------------------------------------------------
    lka <- fixed(log(5.5)); label("First-order absorption rate constant (1/h)") # Supplement, PK table: k_abs = 5.5 1/h
    ltlag <- fixed(log(0.406)); label("Absorption lag time (h)") # Supplement, PK table: Abs_lag = 0.406 h
    lkel <- fixed(log(0.0708)); label("Elimination rate constant from the central compartment (1/h)") # Supplement, PK table: k_el = 0.0708 1/h
    lk12 <- fixed(log(0.0265)); label("Central-to-peripheral rate constant (1/h)") # Supplement, PK table: k12 = 0.0265 1/h
    lk21 <- fixed(log(0.0195)); label("Peripheral-to-central rate constant (1/h)") # Supplement, PK table: k21 = 0.0195 1/h
    # The supplement prints only micro-constants, so the apparent central
    # volume is not stated. It is back-solved from the authors' deposited PK
    # simulation (Supplementary file MOESM2, avadomide 5 mg on the 5/7
    # schedule, concentrations in ng/mL): with the four rate constants above,
    # V/F = 48.70 L reproduces 97% of the deposited profile points to within
    # 0.5% and gives Cmax = 119 ng/mL and AUC(0-28 d) = 1181 ng/mL*day, the
    # 5 mg 5/7 row of Table II / Table S.III.
    lvc <- fixed(log(48.7)); label("Apparent central volume of distribution (L)") # Back-solved from Supplementary file MOESM2 (see comment)

    # -----------------------------------------------------------------------
    # Avadomide PD (Eq. 9) -- Table I
    # -----------------------------------------------------------------------
    # EC50 and n were regressed in this GBM fit, simultaneously with the three
    # disease-specific parameters and across all GBM dose groups (Results),
    # and then carried unchanged into the DLBCL and MM fits.
    lec50 <- log(15); label("Avadomide concentration giving half-maximal maturation block (ng/mL)") # Table I, EC50,PD = 15 ng/mL (regressed on GBM ANC)
    lhill <- log(2); label("Hill coefficient of the maturation-block function (unitless)") # Table I, n_PD = 2 (regressed on GBM ANC)
    emax <- fixed(0.9); label("Maximum avadomide-induced reduction of the transit-2 to transit-3 maturation rate (fraction)") # Table I, Emax,PD = 0.9 (type A)

    # -----------------------------------------------------------------------
    # Neutrophil life-cycle model (Eqs. 1-8) -- Table I, GBM column
    # -----------------------------------------------------------------------
    lgamma <- log(0.02); label("Exponent of the proliferation feedback on transit-2 level (unitless)") # Table I, gamma = 0.02 (GBM, regressed)
    lbeta <- fixed(log(20)); label("Exponent of the marrow-egress feedback on circulating ANC (unitless)") # Table I, beta = 20 (type A)
    # Baseline ANC is patient-specific: Table I prints 4.5E9 cell/L as an
    # example typical value and the virtual patients draw it from the
    # clinical distribution (Figure 4b). Carried here in 10^9 cells/L.
    lcirc0 <- fixed(log(4.5)); label("Baseline circulating neutrophil count (10^9 cells/L)") # Table I, Circ0 = 4.5E9 cell/L (example typical value)
    lthalf_circ <- fixed(log(30)); label("Half-life of circulating neutrophils (h)") # Table I, t1/2,Neutrophils = 30 h (literature, doubled for neutropenia)
    lkd_mat <- fixed(log(0.001)); label("Apoptosis rate constant of maturing bone-marrow neutrophils (1/h)") # Table I, kd = 0.001 1/h (type A)
    lratio_reserv <- log(3); label("Baseline ratio of bone-marrow reservoir to circulating neutrophils (unitless)") # Table I, Ratio Reserv0/Circ0 = 3 (GBM, regressed)
    lkm_frac <- log(0.6); label("Michaelis constant of the transit-2 to transit-3 transfer, as a fraction of baseline transit-2 level (fraction)") # Table I, KM,fraction = 0.6 (GBM, regressed)
  })

  model({
    # -----------------------------------------------------------------------
    # Individual parameters
    # -----------------------------------------------------------------------
    ka <- exp(lka)
    tlag <- exp(ltlag)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc <- exp(lvc)

    ec50 <- exp(lec50)
    hill <- exp(lhill)
    gamma <- exp(lgamma)
    beta <- exp(lbeta)
    circ0 <- exp(lcirc0)
    thalf_circ <- exp(lthalf_circ)
    kd_mat <- exp(lkd_mat)
    ratio_reserv <- exp(lratio_reserv)
    km_frac <- exp(lkm_frac)

    # -----------------------------------------------------------------------
    # Homeostatic back-calculation of the rate constants (Table I, type C
    # rows; identical to the supplementary R script). Every compartment but
    # the reservoir starts at the baseline ANC (Tran0 = Prol0 = Circ0).
    # -----------------------------------------------------------------------
    kelim <- log(2) / thalf_circ
    reserv0 <- ratio_reserv * circ0
    tran0 <- circ0
    kout <- kelim * circ0 / reserv0
    ktr4 <- (kd_mat + kout) * reserv0 / tran0
    ktr3 <- kd_mat + ktr4
    ktr2 <- kd_mat + ktr3
    ktr1 <- kd_mat + ktr2
    kprol <- ktr1
    # Michaelis-Menten form of k_tr,3 = Vmax / (KM + Transit2); Vmax is set
    # so that the homeostatic transfer rate equals ktr3.
    km <- km_frac * tran0
    vmax <- ktr3 * (km + tran0)

    # -----------------------------------------------------------------------
    # Avadomide PK
    # -----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc

    # -----------------------------------------------------------------------
    # PD: fractional maturation rate relative to homeostasis (Eq. 9)
    # -----------------------------------------------------------------------
    eff_mat <- 1 - emax * Cc^hill / (ec50^hill + Cc^hill)

    # -----------------------------------------------------------------------
    # Neutrophil life cycle (Eqs. 1-8)
    # -----------------------------------------------------------------------
    fb_prol <- (tran0 / transit2)^gamma # Eq. 7, feedback proliferation
    # With beta = 20 the egress feedback makes the system extremely stiff once
    # the ANC falls well below baseline: (circ0 / circ)^beta can reach ~1e14
    # and the reservoir is emptied to ~1e-14. The default liblsoda solver can
    # then fail for patients with a high baseline ANC; solve with a stiff BDF
    # solver, e.g. rxSolve(..., method = "cvode", atol = 1e-10, rtol = 1e-10).
    fb_egress <- (circ0 / circ)^beta # Eq. 8, feedback egress
    flux_mat <- vmax * eff_mat * transit2 / (km + transit2)

    d/dt(prol) <- kprol * fb_prol * prol - ktr1 * prol
    d/dt(transit1) <- ktr1 * prol - (ktr2 + kd_mat) * transit1
    d/dt(transit2) <- ktr2 * transit1 - flux_mat - kd_mat * transit2
    d/dt(transit3) <- flux_mat - (ktr4 + kd_mat) * transit3
    d/dt(reservoir) <- ktr4 * transit3 - (kd_mat + kout * fb_egress) * reservoir
    d/dt(circ) <- kout * fb_egress * reservoir - kelim * circ

    prol(0) <- tran0
    transit1(0) <- tran0
    transit2(0) <- tran0
    transit3(0) <- tran0
    reservoir(0) <- reserv0
    circ(0) <- circ0

    # Absolute neutrophil count (10^9 cells/L), the model's observable.
    ANC <- circ
  })
}
