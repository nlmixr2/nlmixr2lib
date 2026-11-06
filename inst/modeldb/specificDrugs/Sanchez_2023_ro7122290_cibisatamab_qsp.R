Sanchez_2023_ro7122290_cibisatamab_qsp <- function() {
  description <- "QSP. Clinical translation of the in vitro FAP-4-1BBL (RO7122290) + cibisatamab co-culture models to a colorectal tumor: two-compartment FAP-4-1BBL PK with parallel linear and Michaelis-Menten clearance and correlated IIV, two-compartment linear cibisatamab PK, first-order plasma-to-tumor distribution of both drugs, cibisatamab-driven 4-1BB expression per CD8+ T cell (in vitro Model 1), FAP-4-1BBL / FAP / 4-1BB trimeric-complex formation in the tumor, and the in vitro Models 3 and 2 for the fold increase in tumor cell killing and the FAP-dependent loss of cibisatamab monotherapy killing (Sanchez 2023 Simulations 2 and 3; valid over the first week because the 4-1BB induction signal and synthesis stop after the in vitro stop time)."
  reference <- "Sanchez J, Claus C, Albrecht R, Gaillard BC, Marinho J, McIntyre C, Tanos T, Boehnke A, Friberg LE, Jonsson S, Frances N. A model-based approach leveraging in vitro data to support dose selection from the outset: A framework for bispecific antibodies in immuno-oncology. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1804-1818. doi:10.1002/psp4.13065"
  vignette <- "Sanchez_2023_ro7122290_cibisatamab"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    tumor = "nM",
    trimer_per_tcell = "trimeric complexes per T cell",
    fold_kill = "fold increase over cibisatamab monotherapy"
  )

  # Cibisatamab has its own plasma and tumor states beside the FAP-4-1BBL
  # ones; the receptor pools and their complexes are tumor concentrations;
  # expr_41bb counts 4-1BB per CD8+ T cell; auc_trimer integrates
  # trimer_per_tcell so Model 3 can read its time average.
  paper_specific_compartments <- c(
    "central_cibisatamab",
    "peripheral1_cibisatamab",
    "tumor_cibisatamab",
    "expr_41bb",
    "target_41bb",
    "target_fap",
    "complex_fap",
    "complex_41bb",
    "auc_trimer"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "FAP-4-1BBL (RO7122290)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "FAP-4-1BBL (RO7122290)", units = "mg", specimen = "plasma", verified = TRUE),
    central_cibisatamab = list(analyte = "cibisatamab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cibisatamab = list(analyte = "cibisatamab", units = "mg", specimen = "plasma", verified = TRUE),
    tumor_cibisatamab = list(analyte = "cibisatamab", units = "mg/L", specimen = "tumor", verified = TRUE),
    tumor = list(analyte = "free FAP-4-1BBL (RO7122290)", units = "nM", specimen = "tumor", verified = TRUE),
    transit1 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    expr_41bb = list(analyte = "4-1BB receptor", units = "receptors/cell", specimen = "tumor", verified = TRUE),
    target_41bb = list(analyte = "free 4-1BB receptor", units = "nM", specimen = "tumor", verified = TRUE),
    target_fap = list(analyte = "free FAP receptor", units = "nM", specimen = "tumor", verified = TRUE),
    complex_fap = list(analyte = "FAP-4-1BBL-FAP dimer", units = "nM", specimen = "tumor", verified = TRUE),
    complex_41bb = list(analyte = "FAP-4-1BBL-4-1BB dimer", units = "nM", specimen = "tumor", verified = TRUE),
    trimer = list(analyte = "FAP-FAP-4-1BBL-4-1BB trimeric complex", units = "nM", specimen = "tumor", verified = TRUE),
    auc_trimer = list(
      analyte = "cumulative trimeric complexes per T cell",
      units = "trimeric complexes per T cell * h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human (virtual patients; pharmacodynamics translated from in vitro co-culture data)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "Colorectal carcinoma (virtual patients); FAP expression from published immunohistochemistry H-scores",
    dose_range = "Cibisatamab 100 mg q3w 2-h IV infusion followed by FAP-4-1BBL 1-180 mg 2-h IV infusion q1w, q2w or q3w",
    notes = paste(
      "No patients were fitted. FAP-4-1BBL population PK comes from Micallef et al. (ACoP 2020) and cibisatamab PK from internal data;",
      "both are tabulated in Table S1. The pharmacodynamic layer is the in vitro Models 1-3 of the same paper.",
      "Tumor geometry, cell densities and the clinical median FAP expression (12% FAP-positive fibroblasts, 12,000 FAP/fibroblast)",
      "come from the Supplementary Materials and Methods and the deposited RxODE code."
    )
  )

  ini({
    # FAP-4-1BBL PK (Table S1, Micallef 2020); values carried in, not estimated here.
    lvc <- fixed(log(3.77)); label("FAP-4-1BBL central volume (L)") # Table S1 V1 3770 (2.42); deposited code s05 'mL for vol'
    lvp <- fixed(log(0.83)); label("FAP-4-1BBL peripheral volume (L)") # Table S1 V2 830 (2.97); deposited code s05 'mL for vol'
    lcl <- fixed(log(0.0154)); label("FAP-4-1BBL linear clearance (L/h)") # Table S1 CL 15.4 mL/h (3.89)
    lq <- fixed(log(0.0219)); label("FAP-4-1BBL intercompartmental clearance (L/h)") # Table S1 Q 21.9 mL/h (6.3)
    lvmax <- fixed(log(0.455)); label("FAP-4-1BBL maximum Michaelis-Menten elimination rate (mg/h)") # Table S1 Vmax 0.455 (5.62), printed ug/h; mg/h per deposited code s05
    lkm <- fixed(log(5.18)); label("FAP-4-1BBL Michaelis-Menten constant (mg/L)") # Table S1 Km 5.18 ug/mL (7.08)
    mw <- fixed(177000); label("FAP-4-1BBL molecular weight (g/mol)") # deposited code s04/s05 'Cc_FAP41BBL_nM = C1/0.177'

    etalvc + etalcl ~ fixed(c(0.485, 0.04, 0.104)) # deposited code s05 Sigma for V1 and CL
    etalvp ~ fixed(1.52) # deposited code s05 'V2 = 830 * exp(rnorm(NN, 0, sqrt(1.52)))'
    etalvmax + etalkm ~ fixed(c(0.306, 0.365, 0.798)) # deposited code s05 Sigma for Vmax and Km

    # Cibisatamab PK (Table S1, internal data).
    lvc_cibisatamab <- fixed(log(3.95)); label("Cibisatamab central volume (L)") # Table S1 V1_Cibi 3.95
    lvp_cibisatamab <- fixed(log(4.38)); label("Cibisatamab peripheral volume (L)") # Table S1 V2_Cibi 4.38
    lcl_cibisatamab <- fixed(log(0.063)); label("Cibisatamab clearance (L/h)") # Table S1 CL_Cibi 0.063
    lq_cibisatamab <- fixed(log(0.039)); label("Cibisatamab intercompartmental clearance (L/h)") # Table S1 Q_Cibi 0.039
    mw_cibisatamab <- fixed(194000); label("Cibisatamab molecular weight (g/mol)") # deposited code s04 'CEA = Tumor_Cibi * 5.1546'

    # Plasma-to-tumor distribution (Table S1); 2.2-fold plasma:tumor scenario.
    ktum_in <- fixed(0.06); label("FAP-4-1BBL plasma-to-tumor distribution rate (1/h)") # Table S1 kain 0.06 (2.2:1 scenario; 10:1 scenario 0.03)
    ktum_out <- fixed(0.137); label("FAP-4-1BBL tumor-to-plasma distribution rate (1/h)") # Table S1 kout 0.137 (2.2:1 scenario; 10:1 scenario 0.3)
    ktum_in_cibisatamab <- fixed(0.06); label("Cibisatamab plasma-to-tumor distribution rate (1/h)") # Table S1 kain 0.06 (cibisatamab)
    ktum_out_cibisatamab <- fixed(0.137); label("Cibisatamab tumor-to-plasma distribution rate (1/h)") # Table S1 kout 0.137 (cibisatamab)

    # In vitro Model 1 (Table 2, Table S1).
    lemax <- fixed(log(537.682)); label("Maximum Signal for cibisatamab, SignalMAX (receptors/cell)") # Table S1 'Emax 537.682'
    lec50 <- fixed(log(16.004)); label("Cibisatamab concentration giving half-maximal Signal (nM)") # Table S1 '16.004'
    lhill <- fixed(log(0.210)); label("Hill coefficient of Signal versus cibisatamab concentration (unitless)") # Table S1 0.210
    lt_stop <- fixed(log(16.18)); label("Stop time after which Signal and kin are zero (h)") # Table S1 '16.18'
    lktr <- fixed(log(0.2801)); label("Transit rate constant of the 7-compartment chain (1/h)") # Table S1 0.2801
    lkin <- fixed(log(2.952)); label("Basal 4-1BB synthesis rate (receptors/cell/h)") # Table S1 2.952
    lkout <- fixed(log(0.2832)); label("4-1BB elimination rate constant (1/h)") # Table S1 0.2832
    added_value <- fixed(1); label("Multiplier on the cibisatamab Signal from FAP-4-1BBL costimulation (unitless)") # deposited code s04 'Added_value = 1'

    # Binding (Table S1).
    kd_fap <- fixed(0.7); label("Equilibrium dissociation constant of FAP-4-1BBL for FAP (nM)") # Table S1 KD,FAP 0.7
    koff_fap <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from FAP (1/h)") # Table S1 koff,FAP 0.1 (assumed)
    kd_41bb <- fixed(0.2); label("Equilibrium dissociation constant of FAP-4-1BBL for 4-1BB (nM)") # Table S1 KD,4-1BB 0.2
    koff_41bb <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from 4-1BB (1/h)") # Table S1 koff,4-1BB 0.1 (assumed)

    # Tumor composition (deposited code s04; Supplementary Materials and Methods, Simulation 1).
    tumor_ldiam <- fixed(50); label("Longest tumor diameter of the ellipsoid tumor (mm)") # deposited code s04 'longest diameter = 50 mm, and other diameters = Longest/3'
    tumor_cells_cm3 <- fixed(1e7); label("Tumor cell density (cells/cm^3)") # deposited code s04 'Tumor_cells_number = 1e7 * Tumor_volume_L * 1000'
    fib_frac <- fixed(0.15); label("Fibroblasts per tumor cell (unitless)") # deposited code s04 'Proportion_Fibroblasts_in_tumor = 0.15'
    et_ratio <- fixed(0.1); label("Effector to tumor cell ratio in the tumor (unitless)") # deposited code s04 'E_T_ratio = 0.1'
    tcell_pct <- fixed(40); label("Percentage of effector cells that are T cells (%)") # deposited code s04 'Tcell_percentage_in_tumor = 40'
    cd8_frac <- fixed(0.4); label("Fraction of T cells carrying the modelled 4-1BB (CD8+) (unitless)") # deposited code s04 factor 0.4 in the T-cell count
    fap_per_fib <- fixed(12000); label("FAP receptors per FAP-positive fibroblast (receptors/cell)") # Supplement Simulation 1: median 12,000 FAP/fibroblast in colorectal carcinoma
    fap_pos_pct <- fixed(12); label("Percentage of tumor fibroblasts expressing FAP (%)") # Supplement Simulation 1: median 12% FAP-positive fibroblasts

    # In vitro Model 3 (Table 3).
    kill_max <- fixed(4.37); label("Maximum fold increase in tumor cell killing over cibisatamab monotherapy (unitless)") # Table 3 4.37
    trimer50 <- fixed(3.91e-3); label("Trimeric complexes per T cell giving half-maximal fold increase (complexes/cell)") # Table 3 3.91e-3
    hill_kill <- fixed(1.16); label("Hill coefficient of fold increase versus trimeric complexes (unitless)") # Table 3 1.16
    kill_base <- fixed(1.04); label("Fold killing without trimeric complex formation (unitless)") # Table 3 1.04

    # In vitro Model 2 (Table 3).
    kill_nofap <- fixed(101.9); label("Cibisatamab monotherapy killing without FAP-positive fibroblasts (%)") # Table 3 101.9
    kill_maxdec <- fixed(67.4); label("Maximum decrease in monotherapy killing due to FAP-positive fibroblasts (%)") # Table 3 67.4
    hill_fap <- fixed(0.86); label("Sigmoid coefficient of the generalized logistic model (unitless)") # Table 3 0.86
    logit_fap50 <- fixed(-1.1); label("Logit FAP-positive proportion at half-maximal decrease (logit units)") # Table 3 -1.1
  })
  model({
    # FAP-4-1BBL PK.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Cibisatamab PK.
    vc_cibisatamab <- exp(lvc_cibisatamab)
    vp_cibisatamab <- exp(lvp_cibisatamab)
    cl_cibisatamab <- exp(lcl_cibisatamab)
    q_cibisatamab <- exp(lq_cibisatamab)

    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    t_stop <- exp(lt_stop)
    ktr <- exp(lktr)
    kin <- exp(lkin)
    kout <- exp(lkout)
    kon_fap <- koff_fap / kd_fap
    kon_41bb <- koff_41bb / kd_41bb

    # Ellipsoid tumor with the two short diameters one third of the longest.
    vol_tumor <- 4 / 3 * pi * (tumor_ldiam / 2) * (tumor_ldiam / 6)^2 / 1e6
    n_tumor <- tumor_cells_cm3 * vol_tumor * 1000
    n_tcell <- et_ratio * n_tumor * tcell_pct / 100 * cd8_frac
    conv <- n_tcell / 6.022e23 * 1e9 / vol_tumor
    # Equation 22 with the tumor volume.
    fap0 <- n_tumor * fib_frac * fap_pos_pct / 100 * fap_per_fib / 6.022e23 * 1e9 / vol_tumor

    Cc <- central / vc
    Cc_nM <- Cc / mw * 1e6
    Cc_cibisatamab <- central_cibisatamab / vc_cibisatamab
    Ctumor_cibisatamab_nM <- tumor_cibisatamab / mw_cibisatamab * 1e6

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - vmax * Cc / (Cc + km)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_cibisatamab) <- -(cl_cibisatamab + q_cibisatamab) / vc_cibisatamab * central_cibisatamab +
      q_cibisatamab / vp_cibisatamab * peripheral1_cibisatamab
    d/dt(peripheral1_cibisatamab) <- q_cibisatamab / vc_cibisatamab * central_cibisatamab -
      q_cibisatamab / vp_cibisatamab * peripheral1_cibisatamab
    # Tumor uptake does not deplete plasma, as in the deposited code.
    d/dt(tumor_cibisatamab) <- ktum_in_cibisatamab * Cc_cibisatamab - ktum_out_cibisatamab * tumor_cibisatamab

    # Equations 10-13 driven by the tumor cibisatamab concentration; t is time
    # since the first cibisatamab dose.
    stepf <- 0.5 * (tanh(1e15 * (t_stop - t)) + 1)
    sig_mono <- emax * Ctumor_cibisatamab_nM^hill / (Ctumor_cibisatamab_nM^hill + ec50^hill)
    signal <- min(sig_mono * added_value, emax) * stepf
    kin_t <- kin * stepf

    d/dt(transit1) <- ktr * (signal - transit1)
    d/dt(transit2) <- ktr * (transit1 - transit2)
    d/dt(transit3) <- ktr * (transit2 - transit3)
    d/dt(transit4) <- ktr * (transit3 - transit4)
    d/dt(transit5) <- ktr * (transit4 - transit5)
    d/dt(transit6) <- ktr * (transit5 - transit6)
    d/dt(transit7) <- ktr * (transit6 - transit7)
    d/dt(expr_41bb) <- ktr * transit7 - kout * expr_41bb + kin_t
    expr_41bb(0) <- kin / kout

    # Equations 16-21 in the deposited-code form.
    d/dt(target_41bb) <- (ktr * transit7 - kout * expr_41bb + kin_t) * conv +
      koff_41bb * (complex_41bb + trimer) - kon_41bb * target_41bb * (tumor + complex_fap)
    target_41bb(0) <- kin / kout * conv
    d/dt(tumor) <- ktum_in * Cc_nM - ktum_out * tumor + koff_fap * complex_fap + koff_41bb * complex_41bb -
      kon_fap * target_fap * tumor - kon_41bb * target_41bb * tumor
    d/dt(complex_fap) <- kon_fap * target_fap * tumor + koff_41bb * trimer -
      kon_41bb * target_41bb * complex_fap - koff_fap * complex_fap
    d/dt(complex_41bb) <- kon_41bb * target_41bb * tumor + koff_fap * trimer -
      kon_fap * target_fap * complex_41bb - koff_41bb * complex_41bb
    d/dt(trimer) <- kon_41bb * target_41bb * complex_fap + kon_fap * target_fap * complex_41bb -
      (koff_fap + koff_41bb) * trimer
    d/dt(target_fap) <- koff_fap * (trimer + complex_fap) - kon_fap * target_fap * (tumor + complex_41bb)
    target_fap(0) <- fap0

    # Equations 23-26.
    ro_fap <- (complex_fap + trimer) / (target_fap + complex_fap + trimer) * 100
    ro_41bb <- (complex_41bb + trimer) / (target_41bb + complex_41bb + trimer) * 100
    trimer_per_tcell <- trimer / conv

    # Model 3 (Equation 27) on the trimeric complexes per T cell averaged from 0 to t.
    d/dt(auc_trimer) <- trimer_per_tcell
    trimer_avg <- trimer_per_tcell
    if (t > 0) {
      trimer_avg <- auc_trimer / t
    }
    fold_kill <- (kill_max - kill_base) * trimer_avg^hill_kill /
      (trimer_avg^hill_kill + trimer50^hill_kill) + kill_base

    # Model 2 (Equation 9) on the tumor FAP-positive fibroblast proportion.
    p_fap <- fap_pos_pct / 100
    logit_fap <- -7
    if (p_fap >= 1) {
      logit_fap <- 7
    } else if (p_fap > 0) {
      logit_fap <- log(p_fap / (1 - p_fap))
    }
    kill_mono_pct <- kill_nofap - kill_maxdec / (1 + exp(-hill_fap * (logit_fap - logit_fap50)))
    kill_abs <- fold_kill * kill_mono_pct / 100
  })
}
