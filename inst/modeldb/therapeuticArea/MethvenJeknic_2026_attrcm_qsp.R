MethvenJeknic_2026_attrcm_qsp <- function() {
  description <- "QSP. Mechanistic mass-action model of transthyretin (TTR) dissociation and amyloid formation in transthyretin amyloid cardiomyopathy (ATTR-CM), with tetramer stabilisation by acoramidis or tafamidis and TTR knockdown by vutrisiran. TTR tetramers are produced at a constant rate, cleared intact, or dissociate reversibly into monomers that unfold reversibly; unfolded monomers form amyloid irreversibly (monomer efflux, the model's disease-progression proxy). Stabilisers bind the two tetramer thyroxine pockets with negative cooperativity and four equivalent albumin sites at a constant total drug concentration; vutrisiran scales tetramer production by a time-varying silencing fraction refitted from the HELIOS-B (or HELIOS-A) TTR reductions. Genotype-specific kinetics for wild-type (ATTRwt), variant (ATTRv, parameterised on p.Val142Ile) and the protective p.Thr139Met TTR. Deterministic typical-cohort model without between-subject variability or residual error."
  reference <- "Methven Jeknic S, Ahn-Horst TA, Ji AX, Ju JH, Paull M, Ji SG, Sinha U. Evaluation of the Current ATTR-CM Treatment Landscape via a Mathematical Model of TTR Dissociation and Amyloid Formation. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70277. doi:10.1002/psp4.70277. Equations from Supplemental Equations 8-16, initial conditions from Table S2, rate constants from Table S6; implementation details (initial state, albumin-site count, drug-release on complex clearance, vutrisiran silencing fit) from the authors' code deposit https://github.com/BridgeBioAnalytics/ttr-mechanistic-model (commit 6686acd)."
  vignette <- "MethvenJeknic_2026_attrcm"
  # The TTR states (tetramer, singly / doubly stabiliser-bound tetramer,
  # folded and unfolded monomer), the accumulated amyloid, and the free and
  # stabiliser-occupied albumin binding sites are the paper's own species
  # (Table S2); none is a reusable canonical compartment. Free stabiliser is
  # held in `central` so a bolus there sets the constant total drug pool.
  paper_specific_compartments <- c(
    "ttr_tetramer",
    "ttr_tetramer_drug1",
    "ttr_tetramer_drug2",
    "ttr_monomer",
    "ttr_unfolded",
    "amyloid",
    "albumin",
    "complex_albumin"
  )
  units <- list(time = "h", dosing = "nM", concentration = "nM")

  covariateData <- list(
    SNP_TTR_RS76992529 = list(
      description = "TTR p.Val142Ile (V122I, rs76992529) variant indicator; selects the ATTRv parameter set",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type TTR, ATTRwt parameter set)",
      notes = "The paper parameterises its whole ATTRv cohort on p.Val142Ile thermodynamics and kinetics (Methods, 'Parameter Optimization'; Supplemental Methods 1.3), so 1 selects the ATTRv tetramer association / dissociation, monomer folding and baseline tetramer values of Table S6 / Table S2. Mutually exclusive with SNP_TTR_RS28933981.",
      source_name = "genotype = 'V122I' (code) / 'ATTRv' (paper)"
    ),
    SNP_TTR_RS28933981 = list(
      description = "TTR p.Thr139Met (T119M, rs28933981) protective variant indicator; selects the p.Thr139Met parameter set",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type TTR, ATTRwt parameter set)",
      notes = "1 selects the p.Thr139Met (kinetically stabilised) tetramer association / dissociation, monomer folding and baseline tetramer values of Table S6 / Table S2. Mutually exclusive with SNP_TTR_RS76992529.",
      source_name = "genotype = 'T119M' (code) / 'p.Thr139Met' (paper)"
    ),
    TRT_ACORAMIDIS = list(
      description = "Acoramidis treatment indicator; selects the acoramidis TTR- and albumin-binding rate constants for the drug in `central`",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no acoramidis)",
      notes = "The total acoramidis concentration is set by a bolus into `central` at time 0 (10000 nM = the 10 uM trough the paper holds constant). Mutually exclusive with TRT_TAFAMIDIS. With neither stabiliser indicator set the drug in `central` does not bind.",
      source_name = "drug = 'acoramidis'"
    ),
    TRT_TAFAMIDIS = list(
      description = "Tafamidis treatment indicator; selects the tafamidis TTR- and albumin-binding rate constants for the drug in `central`",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no tafamidis)",
      notes = "The total tafamidis concentration is set by a bolus into `central` at time 0 (26000 nM = the 26 uM peak the paper holds constant). Mutually exclusive with TRT_ACORAMIDIS.",
      source_name = "drug = 'tafamidis'"
    ),
    TRT_VUTRISIRAN = list(
      description = "Vutrisiran treatment indicator; scales tetramer production by (1 - silencing fraction) from time 0",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no knockdown)",
      notes = "Treatment starts at model time 0. The silencing fraction follows the single-exponential fit to the HELIOS-B (default) or HELIOS-A (STUDY_HELIOS_A = 1) TTR reductions, held at its first-measurement value before the first post-dose visit (Supplemental Methods 1.3.3).",
      source_name = "drug = 'vutrisiran_heliosb' / 'vutrisiran_heliosa'"
    ),
    STUDY_HELIOS_A = list(
      description = "Selects the HELIOS-A (hATTR polyneuropathy) instead of the HELIOS-B (ATTR-CM) vutrisiran silencing trajectory",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HELIOS-B trajectory, used for the main-text Figure 5)",
      notes = "Only acts when TRT_VUTRISIRAN = 1. HELIOS-A is used in the paper's Figures S11 and S12.",
      source_name = "drug = 'vutrisiran_heliosa'"
    )
  )

  compartmentData <- list(
    ttr_tetramer = list(
      analyte = "transthyretin tetramer, unbound by stabiliser",
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
    ttr_monomer = list(analyte = "folded transthyretin monomer", units = "nM", specimen = "serum", verified = TRUE),
    ttr_unfolded = list(
      analyte = "unfolded transthyretin monomer (amyloid precursor)",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    amyloid = list(
      analyte = "cumulative transthyretin amyloid formed (monomer units)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
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
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "Transthyretin amyloid cardiomyopathy, wild-type (ATTRwt) and variant (ATTRv); p.Thr139Met carriers as a protective-genotype reference",
    dose_range = "Acoramidis at a constant 10 uM total serum concentration (trough; 8 uM sensitivity case), tafamidis at a constant 26 uM (peak; 35 uM sensitivity case), vutrisiran as HELIOS-B or HELIOS-A TTR knockdown",
    regions = "Not applicable (calibrated to published cohort means)",
    notes = "Not a fit to individual patient data. Physiological parameters were calibrated by least squares (LMFIT, Powell, 50 random starts) to cohort-mean serum TTR in ATTRwt / ATTRv / p.Thr139Met cohorts (Alexander 2025, Tess 2023, Hanson 2018, Arvanitis 2017, Greve 2021), TTR turnover and half-life, tetramer dissociation and folding thermodynamics (Hurshman Babbes 2008, Yee 2019), and treatment-induced changes in serum TTR in ATTRibute-CM, ATTR-ACT and HELIOS-B (Table S4). Stabiliser binding constants come from a separate fit of the in vitro FPE model (MethvenJeknic_2026_attrcm_fpe_invitro) adjusted from 25 to 37 C with the van 't Hoff equation."
  )

  ini({
    # ---- TTR turnover (Table S6) ----
    lkin_tet <- log(5.849e-8); label("Tetramer production rate, kP (M/h)") # Table S6 kP = 5.849e-8 M/h (Fit)
    lkout_tet <- log(1.075e-2); label("First-order clearance of intact tetramer, free or stabiliser-bound, kLT (1/h)") # Table S6 kLT = 1.075e-2 1/h (Fit)
    lkout_mon <- log(1.3005); label("First-order clearance of folded and unfolded monomer, kLM (1/h)") # Table S6 kLM = 1.3005 1/h (Fit)

    # ---- Genotype-specific tetramer association (Table S6) ----
    lkat_attrwt <- log(2.791e19); label("Monomer association into tetramer, ATTRwt, kaT (1/M^3/h)") # Table S6 kaT ATTRwt = 2.791e19 (Fit)
    lkat_attrv <- log(8.912e15); label("Monomer association into tetramer, ATTRv, kaT (1/M^3/h)") # Table S6 kaT ATTRv = 8.912e15 (Fit)
    lkat_t139m <- log(1.101e26); label("Monomer association into tetramer, p.Thr139Met, kaT (1/M^3/h)") # Table S6 kaT p.Thr139Met = 1.101e26 (Fit)

    # ---- Genotype-specific tetramer dissociation (Table S6) ----
    lkdt_attrwt <- log(2.080e-3); label("Tetramer dissociation into monomers, ATTRwt, kdT (1/h)") # Table S6 kdT ATTRwt = 2.080e-3 (Fit)
    lkdt_attrv <- log(8.196e-3); label("Tetramer dissociation into monomers, ATTRv, kdT (1/h)") # Table S6 kdT ATTRv = 8.196e-3 (Fit)
    lkdt_t139m <- log(4.520e-4); label("Tetramer dissociation into monomers, p.Thr139Met, kdT (1/h)") # Table S6 kdT p.Thr139Met = 4.520e-4 (Fit)

    # ---- Monomer unfolding / folding (Table S6) ----
    lku <- log(8.567e-2); label("Monomer unfolding, ku, shared by all genotypes (1/h)") # Table S6 ku = 8.567e-2 for ATTRwt, ATTRv and p.Thr139Met (Fit)
    lkf_attrwt <- log(50.971); label("Unfolded-monomer refolding, ATTRwt, kf (1/h)") # Table S6 kf ATTRwt = 50.971 (Fit)
    lkf_attrv <- log(303.746); label("Unfolded-monomer refolding, ATTRv, kf (1/h)") # Table S6 kf ATTRv = 303.746 (Fit)
    lkf_t139m <- log(945.83); label("Unfolded-monomer refolding, p.Thr139Met, kf (1/h)") # Table S6 kf p.Thr139Met = 945.83 (Fit)
    ky <- fixed(3.6e3); label("Amyloid formation from unfolded monomer, kY (1/h)") # Table S6 kY = 3.6e3 (Assumed)

    # ---- Stabiliser binding at 37 C (Table S6) ----
    kon_ttr <- fixed(3.6e10); label("Stabiliser association with either tetramer pocket, kaD1 = kaD2, both drugs (1/M/h)") # Table S6 kaD,1 and kaD,2 = 3.6e10 for acoramidis and tafamidis (Assumed)
    kon_alb <- fixed(3.6e7); label("Stabiliser association with an albumin site, kaDA, both drugs (1/M/h)") # Table S6 kaD,A = 3.6e7 for acoramidis and tafamidis (Assumed)
    lkoff1_acoramidis <- log(1.050e3); label("Acoramidis unbinding from singly bound tetramer, kdD1 (1/h)") # Table S6 kdD,1 acoramidis = 1.050e3 (Fit)
    lkoff1_tafamidis <- log(4.3305e2); label("Tafamidis unbinding from singly bound tetramer, kdD1 (1/h)") # Table S6 kdD,1 tafamidis = 4.3305e2 (Fit)
    lkoff2_acoramidis <- log(8.319e3); label("Acoramidis unbinding from doubly bound tetramer, per site, kdD2 (1/h)") # Table S6 kdD,2 acoramidis = 8.319e3 (Fit)
    lkoff2_tafamidis <- log(6.205e3); label("Tafamidis unbinding from doubly bound tetramer, per site, kdD2 (1/h)") # Table S6 kdD,2 tafamidis = 6.205e3 (Fit)
    lkoff_alb_acoramidis <- log(9.004e3); label("Acoramidis unbinding from an albumin site, kdDA (1/h)") # Table S6 kdD,A acoramidis = 9.004e3 (Fit)
    lkoff_alb_tafamidis <- log(1.098e2); label("Tafamidis unbinding from an albumin site, kdDA (1/h)") # Table S6 kdD,A tafamidis = 1.098e2 (Fit)

    # ---- Dissociation of stabiliser-bound tetramer (Table S6) ----
    lkdtd1 <- log(5.671e-9); label("Dissociation of singly bound tetramer into monomers, kdT:D (1/h)") # Table S6 kdT:D = 5.671e-9 (Fit; table unit '1/M h' is a misprint, the term is first order)
    lkdtd2 <- log(2.817e-10); label("Dissociation of doubly bound tetramer into monomers, kdT:D2 (1/h)") # Table S6 kdT:D2 = 2.817e-10 (Fit; table unit '1/M h' is a misprint)

    # ---- Initial conditions (Table S2) ----
    tet0_attrwt <- fixed(4.21); label("Initial serum tetramer, ATTRwt (uM)") # Table S2 footnote d, ATTRwt 4.21 uM
    tet0_attrv <- fixed(3.09); label("Initial serum tetramer, ATTRv (uM)") # Table S2 footnote d, ATTRv 3.09 uM
    tet0_t139m <- fixed(5.04); label("Initial serum tetramer, p.Thr139Met (uM)") # Table S2 footnote d, p.Thr139Met 5.04 uM
    mon0 <- fixed(0.0411); label("Initial serum monomer (uM)") # Table S2, M = 0.0411 uM
    alb0 <- fixed(2600); label("Albumin drug-binding sites, four per albumin molecule (uM)") # Table S2, A = 2600 uM (code: 650 uM albumin x 4 sites)

    # ---- Vutrisiran silencing, s(t) = smax - (smax - s0) * exp(-ksil * t) ----
    # The paper fits but does not print these; values reproduced by refitting
    # the authors' deposited HELIOS data (data/external/*.csv) with the code's
    # model form (ttr/physiological.py::_silencing_fit, curve_fit LM).
    smax_heliosb <- 0.79688; label("Plateau TTR silencing fraction, HELIOS-B (fraction)") # refit of fontana_et_al_2025_nejm.csv (code deposit), A = 0.79688
    s0_heliosb <- 0.57114; label("Intercept TTR silencing fraction at time 0 of the fit, HELIOS-B (fraction)") # refit of fontana_et_al_2025_nejm.csv, B = 0.57114
    lksil_heliosb <- log(0.17898); label("Rate of approach to the silencing plateau, HELIOS-B (1/month)") # refit of fontana_et_al_2025_nejm.csv, k = 0.17898 1/month
    tsil_heliosb <- fixed(1.5); label("First post-dose HELIOS-B TTR visit; silencing held at its value here before (month)") # Results, 'first time point measured was 6 weeks'; fontana_et_al_2025_nejm.csv first month = 1.5
    smax_heliosa <- 0.85088; label("Plateau TTR silencing fraction, HELIOS-A (fraction)") # refit of adams_et_al_2023_amyloid.csv (code deposit), A = 0.85088
    s0_heliosa <- 0.63285; label("Intercept TTR silencing fraction at time 0 of the fit, HELIOS-A (fraction)") # refit of adams_et_al_2023_amyloid.csv, B = 0.63285
    lksil_heliosa <- log(0.18291); label("Rate of approach to the silencing plateau, HELIOS-A (1/month)") # refit of adams_et_al_2023_amyloid.csv, k = 0.18291 1/month
    tsil_heliosa <- fixed(0.75); label("First post-dose HELIOS-A TTR visit; silencing held at its value here before (month)") # adams_et_al_2023_amyloid.csv first month = 0.75

    # ---- Unit conversion for the reported serum TTR ----
    mw_ttr <- fixed(55000); label("TTR tetramer molecular weight (g/mol)") # code deposit plot.py TTR_UM_TO_MGDL, '55 kDa'
  })
  model({
    # Genotype selection; the reference (both indicators 0) is wild-type TTR.
    g_attrv <- SNP_TTR_RS76992529
    g_t139m <- SNP_TTR_RS28933981
    g_attrwt <- 1 - g_attrv - g_t139m

    # Rate constants converted from the paper's molar units to nM:
    # zero order M/h x 1e9, second order 1/M/h x 1e-9, fifth order
    # (tetramer from four monomers) 1/M^3/h x 1e-27.
    kin <- exp(lkin_tet) * 1e9
    kout_tet <- exp(lkout_tet)
    kout_mon <- exp(lkout_mon)
    kat <- (g_attrwt * exp(lkat_attrwt) + g_attrv * exp(lkat_attrv) + g_t139m * exp(lkat_t139m)) * 1e-27
    kdt <- g_attrwt * exp(lkdt_attrwt) + g_attrv * exp(lkdt_attrv) + g_t139m * exp(lkdt_t139m)
    ku <- exp(lku)
    kf <- g_attrwt * exp(lkf_attrwt) + g_attrv * exp(lkf_attrv) + g_t139m * exp(lkf_t139m)
    kdtd1 <- exp(lkdtd1)
    kdtd2 <- exp(lkdtd2)

    # Stabiliser binding: the indicator picks the drug's off-rates; with no
    # stabiliser indicator the association rates are zero and the drug in
    # `central` is inert.
    on_stab <- TRT_ACORAMIDIS + TRT_TAFAMIDIS
    kon <- kon_ttr * 1e-9 * on_stab
    kona <- kon_alb * 1e-9 * on_stab
    koff1 <- TRT_ACORAMIDIS * exp(lkoff1_acoramidis) + TRT_TAFAMIDIS * exp(lkoff1_tafamidis)
    koff2 <- TRT_ACORAMIDIS * exp(lkoff2_acoramidis) + TRT_TAFAMIDIS * exp(lkoff2_tafamidis)
    koffa <- TRT_ACORAMIDIS * exp(lkoff_alb_acoramidis) + TRT_TAFAMIDIS * exp(lkoff_alb_tafamidis)

    # Vutrisiran knockdown, alpha in Eq 8 = 1 - silencing. Model time is time
    # since the first dose; before the first post-dose visit the silencing is
    # evaluated at that visit (code: silencing_at clamps t to t_first).
    # A month is 30 days (code: HOURS_PER_MONTH = 720).
    smax <- (1 - STUDY_HELIOS_A) * smax_heliosb + STUDY_HELIOS_A * smax_heliosa
    s0 <- (1 - STUDY_HELIOS_A) * s0_heliosb + STUDY_HELIOS_A * s0_heliosa
    ksil <- (1 - STUDY_HELIOS_A) * exp(lksil_heliosb) + STUDY_HELIOS_A * exp(lksil_heliosa)
    tsil <- (1 - STUDY_HELIOS_A) * tsil_heliosb + STUDY_HELIOS_A * tsil_heliosa
    t_month <- t / 720
    if (t_month < tsil) {
      t_month <- tsil
    }
    silencing <- TRT_VUTRISIRAN * (smax - (smax - s0) * exp(-ksil * t_month))
    alpha <- 1 - silencing

    # Supplemental Equations 8-16 (concentrations in nM). Clearance of a
    # stabiliser-bound tetramer (kout_tet) and its dissociation (kdtd1,
    # kdtd2) release the bound drug back to `central`, so total drug is
    # conserved.
    d/dt(ttr_tetramer) <- alpha * kin + kat * ttr_monomer^4 + koff1 * ttr_tetramer_drug1 -
      (kdt + kout_tet + 2 * kon * central) * ttr_tetramer
    d/dt(ttr_tetramer_drug1) <- 2 * (kon * ttr_tetramer * central + koff2 * ttr_tetramer_drug2) -
      (koff1 + kon * central + kout_tet + kdtd1) * ttr_tetramer_drug1
    d/dt(ttr_tetramer_drug2) <- kon * ttr_tetramer_drug1 * central -
      (2 * koff2 + kout_tet + kdtd2) * ttr_tetramer_drug2
    d/dt(ttr_monomer) <- kf * ttr_unfolded +
      4 * (kdt * ttr_tetramer + kdtd1 * ttr_tetramer_drug1 + kdtd2 * ttr_tetramer_drug2) -
      (ku + kout_mon + 4 * kat * ttr_monomer^3) * ttr_monomer
    d/dt(ttr_unfolded) <- ku * ttr_monomer - (kout_mon + ky + kf) * ttr_unfolded
    d/dt(amyloid) <- ky * ttr_unfolded
    d/dt(albumin) <- koffa * complex_albumin - kona * albumin * central
    d/dt(complex_albumin) <- kona * albumin * central - koffa * complex_albumin
    d/dt(central) <- (koff1 + kdtd1 + kout_tet - kon * central) * ttr_tetramer_drug1 +
      2 * (koff2 + kdtd2 + kout_tet) * ttr_tetramer_drug2 + koffa * complex_albumin -
      2 * kon * ttr_tetramer * central - kona * albumin * central

    # Table S2 initial state (uM -> nM); every other species starts at 0 and
    # the stabiliser is added as a bolus into `central` at time 0.
    ttr_tetramer(0) <- 1000 * (g_attrwt * tet0_attrwt + g_attrv * tet0_attrv + g_t139m * tet0_t139m)
    ttr_monomer(0) <- 1000 * mon0
    albumin(0) <- 1000 * alb0

    # Outputs. Serum TTR is total tetramer (free + stabiliser-bound) in mg/dL
    # (nM x 1e-9 mol/nmol x MW g/mol x 100 = mg/dL). Monomer efflux is the
    # amyloid-formation rate, the paper's disease-progression proxy; the
    # `amyloid` state is its running integral (the efflux AUC).
    ttr_tet_total <- ttr_tetramer + ttr_tetramer_drug1 + ttr_tetramer_drug2
    ttr <- ttr_tet_total * 1e-7 * mw_ttr
    efflux <- ky * ttr_unfolded
    frac_tet_bound <- (ttr_tetramer_drug1 + ttr_tetramer_drug2) / ttr_tet_total
    drug_total <- central + ttr_tetramer_drug1 + 2 * ttr_tetramer_drug2 + complex_albumin
    frac_drug_tet <- 0
    frac_drug_alb <- 0
    if (drug_total > 0) {
      frac_drug_tet <- (ttr_tetramer_drug1 + 2 * ttr_tetramer_drug2) / drug_total
      frac_drug_alb <- complex_albumin / drug_total
    }
    # Share of TTR (tetramer equivalents) leaving as monomer through
    # clearance rather than intact or as amyloid (code plot.py
    # _steady_state_metrics, Figure 3b / Figure 4a).
    clear_tet <- kout_tet * ttr_tet_total
    clear_mon <- 0.25 * kout_mon * (ttr_monomer + ttr_unfolded)
    clear_amy <- 0.25 * efflux
    frac_mon_clear <- clear_mon / (clear_tet + clear_mon + clear_amy)
  })
}
