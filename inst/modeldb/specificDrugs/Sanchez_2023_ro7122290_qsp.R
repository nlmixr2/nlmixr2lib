Sanchez_2023_ro7122290_qsp <- function() {
  description <- "QSP. Clinical dose-selection model for FAP-4-1BBL (RO7122290) given with cibisatamab in colorectal carcinoma: two-compartment FAP-4-1BBL PK with parallel linear and Michaelis-Menten clearance and correlated IIV, first-order plasma-to-tumor distribution, and FAP-4-1BBL / FAP / 4-1BB trimeric-complex formation in the tumor with 4-1BB held at 80 receptors per T cell (the level expected at steady-state cibisatamab 100 mg q3w) and FAP at the clinical median (Sanchez 2023 Simulation 1)."
  reference <- "Sanchez J, Claus C, Albrecht R, Gaillard BC, Marinho J, McIntyre C, Tanos T, Boehnke A, Friberg LE, Jonsson S, Frances N. A model-based approach leveraging in vitro data to support dose selection from the outset: A framework for bispecific antibodies in immuno-oncology. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1804-1818. doi:10.1002/psp4.13065"
  vignette <- "Sanchez_2023_ro7122290_cibisatamab"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    tumor = "nM",
    trimer_per_tcell = "trimeric complexes per T cell"
  )

  # The receptor pools and their FAP-4-1BBL complexes are tumor concentrations.
  paper_specific_compartments <- c(
    "target_41bb",
    "target_fap",
    "complex_fap",
    "complex_41bb"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "FAP-4-1BBL (RO7122290)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "FAP-4-1BBL (RO7122290)", units = "mg", specimen = "plasma", verified = TRUE),
    tumor = list(analyte = "free FAP-4-1BBL (RO7122290)", units = "nM", specimen = "tumor", verified = TRUE),
    target_41bb = list(analyte = "free 4-1BB receptor", units = "nM", specimen = "tumor", verified = TRUE),
    target_fap = list(analyte = "free FAP receptor", units = "nM", specimen = "tumor", verified = TRUE),
    complex_fap = list(analyte = "FAP-4-1BBL-FAP dimer", units = "nM", specimen = "tumor", verified = TRUE),
    complex_41bb = list(analyte = "FAP-4-1BBL-4-1BB dimer", units = "nM", specimen = "tumor", verified = TRUE),
    trimer = list(analyte = "FAP-FAP-4-1BBL-4-1BB trimeric complex", units = "nM", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species = "human (virtual patients; binding layer from in vitro data)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = "Colorectal carcinoma (virtual patients) receiving cibisatamab 100 mg q3w",
    dose_range = "FAP-4-1BBL 1-180 mg 2-h IV infusion q1w, q2w or q3w for 30 weeks",
    notes = paste(
      "No patients were fitted. FAP-4-1BBL population PK comes from Micallef et al. (ACoP 2020), tabulated in Table S1;",
      "the random effects are those of the deposited dose-selection code. 4-1BB is held at 80 receptors per T cell and FAP at the",
      "median of colorectal carcinoma (12% FAP-positive fibroblasts, 12,000 FAP/fibroblast; Supplementary Materials and Methods, Simulation 1)."
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
    mw <- fixed(177000); label("FAP-4-1BBL molecular weight (g/mol)") # deposited code s05 'Cc_FAP41BBL_nM = Cc_FAP41BBL/0.177'

    etalvc + etalcl ~ fixed(c(0.485, 0.04, 0.104)) # deposited code s05 Sigma for V1 and CL
    etalvp ~ fixed(1.52) # deposited code s05 'V2 = 830 * exp(rnorm(NN, 0, sqrt(1.52)))'
    etalvmax + etalkm ~ fixed(c(0.306, 0.365, 0.798)) # deposited code s05 Sigma for Vmax and Km

    # Plasma-to-tumor distribution (Table S1); 2.2-fold plasma:tumor scenario.
    ktum_in <- fixed(0.06); label("FAP-4-1BBL plasma-to-tumor distribution rate (1/h)") # Table S1 kain 0.06 (2.2:1 scenario; 10:1 scenario 0.03)
    ktum_out <- fixed(0.137); label("FAP-4-1BBL tumor-to-plasma distribution rate (1/h)") # Table S1 kout 0.137 (2.2:1 scenario; 10:1 scenario 0.3)

    # Binding (Table S1).
    kd_fap <- fixed(0.7); label("Equilibrium dissociation constant of FAP-4-1BBL for FAP (nM)") # Table S1 KD,FAP 0.7
    koff_fap <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from FAP (1/h)") # Table S1 koff,FAP 0.1 (assumed)
    kd_41bb <- fixed(0.2); label("Equilibrium dissociation constant of FAP-4-1BBL for 4-1BB (nM)") # Table S1 KD,4-1BB 0.2
    koff_41bb <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from 4-1BB (1/h)") # Table S1 koff,4-1BB 0.1 (assumed)

    # Tumor composition (deposited code s05; Supplementary Materials and Methods, Simulation 1).
    r41bb_tcell <- fixed(80); label("4-1BB receptors per T cell, held constant (receptors/cell)") # Supplement Simulation 1: 80 4-1BB receptors per T cell
    tumor_ldiam <- fixed(50); label("Longest tumor diameter of the ellipsoid tumor (mm)") # deposited code s05 'longest diameter = 50 mm, and other diameters = Longest/3'
    tumor_cells_cm3 <- fixed(1e7); label("Tumor cell density (cells/cm^3)") # deposited code s05 'Tumor_cells_number = 1e7 * Tumor_volume_L * 1000'
    fib_frac <- fixed(0.15); label("Fibroblasts per tumor cell (unitless)") # deposited code s05 'Proportion_Fibroblasts_in_tumor = 0.15'
    et_ratio <- fixed(0.1); label("Effector to tumor cell ratio in the tumor (unitless)") # deposited code s05 'E_T_ratio = 0.1'
    tcell_pct <- fixed(40); label("Percentage of effector cells that are T cells (%)") # deposited code s05 'Tcell_percentage_in_tumor = 40'
    cd8_frac <- fixed(0.4); label("Fraction of T cells carrying the modelled 4-1BB (CD8+) (unitless)") # deposited code s05 factor 0.4 in the T-cell count
    fap_per_fib <- fixed(12000); label("FAP receptors per FAP-positive fibroblast (receptors/cell)") # Supplement Simulation 1: median 12,000 FAP/fibroblast in colorectal carcinoma
    fap_pos_pct <- fixed(12); label("Percentage of tumor fibroblasts expressing FAP (%)") # Supplement Simulation 1: median 12% FAP-positive fibroblasts
  })
  model({
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kon_fap <- koff_fap / kd_fap
    kon_41bb <- koff_41bb / kd_41bb

    # Ellipsoid tumor with the two short diameters one third of the longest.
    vol_tumor <- 4 / 3 * pi * (tumor_ldiam / 2) * (tumor_ldiam / 6)^2 / 1e6
    n_tumor <- tumor_cells_cm3 * vol_tumor * 1000
    n_tcell <- et_ratio * n_tumor * tcell_pct / 100 * cd8_frac
    fap0 <- n_tumor * fib_frac * fap_pos_pct / 100 * fap_per_fib / 6.022e23 * 1e9 / vol_tumor

    Cc <- central / vc
    Cc_nM <- Cc / mw * 1e6

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - vmax * Cc / (Cc + km)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equations 17-21 with 4-1BB conserved (no turnover).
    d/dt(target_41bb) <- koff_41bb * (complex_41bb + trimer) - kon_41bb * target_41bb * (tumor + complex_fap)
    target_41bb(0) <- n_tcell * r41bb_tcell / 6.022e23 * 1e9 / vol_tumor
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

    # Equations 23, 24 and 26.
    ro_fap <- (complex_fap + trimer) / (target_fap + complex_fap + trimer) * 100
    ro_41bb <- (complex_41bb + trimer) / (target_41bb + complex_41bb + trimer) * 100
    trimer_per_tcell <- trimer / (target_41bb + complex_41bb + trimer) * r41bb_tcell
  })
}
