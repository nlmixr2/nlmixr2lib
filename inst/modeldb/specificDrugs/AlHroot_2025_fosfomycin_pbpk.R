AlHroot_2025_fosfomycin_pbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based). Intravenous fosfomycin (ZTI-01) in healthy adults:",
    "plasma, a lumped residual-tissue compartment, and a lung resolved into interstitial fluid,",
    "lung tissue, epithelial lining fluid (ELF) and alveolar macrophages (AMs). The lung",
    "interstitial fluid receives the whole cardiac output, the residual tissues a fitted",
    "fraction of it; lung tissue and ELF exchange with the interstitial fluid by",
    "permeability-surface products and partition coefficients, and drug moves between ELF and",
    "AMs by first-order uptake and efflux. Linear elimination from plasma. Physiological",
    "volumes and cardiac output are fixed 70 kg values; plasma, ELF and AM concentrations",
    "were fitted simultaneously with a single proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Al Hroot J, Reeder J, Yuan X, Gu K, Walter E, Boole L, Que L, An G.",
    "Determining an appropriate fosfomycin (ZTI-01) dosing regimen in pneumonia patients by",
    "utilizing minimal PBPK modeling and target attainment analysis.",
    "Antimicrob Agents Chemother. 2025;69(6):e01869-24. doi:10.1128/aac.01869-24.",
    sep = " "
  )
  vignette <- "AlHroot_2025_fosfomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The AM state has no registered compartment name.
  paper_specific_compartments <- c("alveolar_macrophage")

  compartmentData <- list(
    plasma = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    is_lung = list(analyte = "fosfomycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "fosfomycin", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "fosfomycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    alveolar_macrophage = list(
      analyte = "fosfomycin",
      units = "mg",
      specimen = "bronchoalveolar lavage",
      verified = TRUE
    ),
    other = list(analyte = "fosfomycin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "18-45 years",
    disease_state = "Healthy adult volunteers with normal renal function.",
    dose_range = paste(
      "Three doses of 6 g intravenous fosfomycin (ZTI-01) infused over 1 h every 8 h;",
      "bronchoalveolar lavage at one randomised time (0.5, 1.25, 2, 5 or 8 h) after the",
      "third dose in 30 participants (6 per time point).",
      sep = " "
    ),
    regions = "United States (Duke University Early Phase Clinical Research Unit)",
    notes = paste(
      "NIAID Phase 1 intrapulmonary PK study NCT03910673 (Methods 'Data source' and Fig. 1).",
      "Plasma sampled pre-dose and at 0.5, 1, 1.25, 2 and 5 h after doses 1 and 3, plus 8 and",
      "12 h after dose 3. The paper reports no weight, sex or race summary. Pneumonia patients",
      "were simulated, not studied: the paper lowers the ELF partition coefficient to 0.29 and",
      "raises the ELF volume to 0.039 L (Table 3) to represent them.",
      sep = " "
    )
  )

  ini({
    # Fixed physiology for a 70 kg human (Table 1).
    qc <- fixed(336); label("Cardiac output, QC (L/h)") # Table 1 QC = 336 L/hr
    v_plasma <- fixed(3); label("Plasma volume, Vp (L)") # Table 1 Vp = 3 L
    v_other <- fixed(66); label("Residual (all other organs) tissue volume, Vresidual (L)") # Table 1 Vresidual = 66 L
    v_is_lung <- fixed(0.3); label("Lung interstitial fluid volume, VISF (L)") # Table 1 VISF = 0.3 L
    v_lung <- fixed(0.7); label("Lung tissue volume, VLung (L)") # Table 1 VLung = 0.7 L
    v_elf <- fixed(0.026); label("Epithelial lining fluid volume, VELF (L)") # Table 1 VELF = 0.026 L
    v_alveolar_macrophage <- fixed(0.0026); label("Alveolar macrophage volume, VAMs (L)") # Table 1 VAMs = 0.0026 L

    # Estimated drug-specific parameters (Table 2).
    lcl <- log(6.59); label("Systemic clearance from plasma, CL (L/h)") # Results 'Estimated parameters': 'systematic clearance was estimated to be 6.59 L/h'; Table 2 prints CL = 7.0 (1.9%) -- see vignette Errata
    lfq_other <- log(0.396); label("Fraction of cardiac output perfusing the residual tissues, Fd (unitless)") # Table 2 Fd = 0.396 (2.8%)
    lkp_lung <- log(11.6); label("Lung tissue partition coefficient, Kp lung tissues (unitless)") # Table 2 KpLung tissues = 11.6 (6.6%)
    lkp_other <- log(0.142); label("Residual tissue partition coefficient, Kpresidual (unitless)") # Table 2 Kpresidual = 0.142 (7.8%)
    lkp_elf <- log(0.393); label("ELF partition coefficient, KpELF (unitless)") # Table 2 KpELF = 0.393 (5.1%)
    lkp_alveolar_macrophage <- log(0.000023); label("Alveolar macrophage partition coefficient, KpAM (unitless)") # Table 2 KpAM = 0.000023 (9.1%)
    lps_lung <- log(4.55); label("Permeability-surface product between lung ISF and lung tissue, PS1 (L/h)") # Table 2 PS1 = 4.55 (4.3%)
    lps_elf <- log(0.0316); label("Permeability-surface product between lung ISF and ELF, PS2 (L/h)") # Table 2 PS2 = 0.0316 (52.5%)
    lk_elf_alveolar_macrophage <- log(0.00515); label("ELF-to-AM uptake rate constant, K1 (1/h)") # Table 2 K1 = 0.00515 (11.9%)
    lk_alveolar_macrophage_elf <- log(0.000006); label("AM-to-ELF efflux rate constant, K2 (1/h)") # Table 2 K2 = 0.000006 (32.0%)

    # IIV: Table 2 reports CV%; omega^2 = log(1 + CV^2).
    etalcl ~ 0.01453 # Table 2 'IIV- CL' = 12.1% (RSE 4.0%, shrinkage 3.0%); log(1 + 0.121^2)
    etalk_alveolar_macrophage_elf ~ 0.14842 # Table 2 'IIV- K2' = 40.0% (RSE 9.0%, shrinkage 20.0%); log(1 + 0.40^2)

    # Residual error: the paper estimated ONE proportional variance for plasma,
    # ELF and AMs together. nlmixr2 needs a distinct parameter per endpoint, so
    # the shared value is carried three times; re-fits should constrain them
    # to be equal.
    propSd <- 0.17407; label("Proportional residual error, plasma (fraction)") # Table 2 sigma^2 (proportional) = 0.0303 (3.0%); sqrt(0.0303)
    propSd_Celf <- 0.17407; label("Proportional residual error, ELF (fraction)") # Table 2 sigma^2 (proportional) = 0.0303, shared with plasma and AMs
    propSd_Calveolar_macrophage <- 0.17407; label("Proportional residual error, alveolar macrophages (fraction)") # Table 2 sigma^2 (proportional) = 0.0303, shared with plasma and ELF
  })

  model({
    cl <- exp(lcl + etalcl)
    fq_other <- exp(lfq_other)
    kp_lung <- exp(lkp_lung)
    kp_other <- exp(lkp_other)
    kp_elf <- exp(lkp_elf)
    kp_alveolar_macrophage <- exp(lkp_alveolar_macrophage)
    ps_lung <- exp(lps_lung)
    ps_elf <- exp(lps_elf)
    k_elf_alveolar_macrophage <- exp(lk_elf_alveolar_macrophage)
    k_alveolar_macrophage_elf <- exp(lk_alveolar_macrophage_elf + etalk_alveolar_macrophage_elf)

    # Concentrations driving the exchange terms. Tissue terms are divided by
    # their partition coefficient, as in Equations 1-6.
    c_plasma <- plasma / v_plasma
    c_is_lung <- is_lung / v_is_lung
    c_lung_free <- lung / (v_lung * kp_lung)
    c_elf_free <- elf / (v_elf * kp_elf)
    c_other_free <- other / (v_other * kp_other)

    # Equations 1-6. Dose is an IV infusion into plasma.
    d/dt(plasma) <- qc * (c_is_lung - c_plasma) +
      fq_other * qc * (c_other_free - c_plasma) - cl * c_plasma
    d/dt(is_lung) <- qc * (c_plasma - c_is_lung) +
      ps_lung * (c_lung_free - c_is_lung) - ps_elf * (c_is_lung - c_elf_free)
    d/dt(lung) <- ps_lung * (c_is_lung - c_lung_free)
    d/dt(elf) <- ps_elf * (c_is_lung - c_elf_free) -
      k_elf_alveolar_macrophage * elf / kp_elf +
      k_alveolar_macrophage_elf * alveolar_macrophage / kp_alveolar_macrophage
    d/dt(alveolar_macrophage) <- k_elf_alveolar_macrophage * elf / kp_elf -
      k_alveolar_macrophage_elf * alveolar_macrophage / kp_alveolar_macrophage
    d/dt(other) <- fq_other * qc * (c_plasma - c_other_free)

    # Measured concentrations are amount over the physiological volume.
    Cc <- plasma / v_plasma
    Celf <- elf / v_elf
    Calveolar_macrophage <- alveolar_macrophage / v_alveolar_macrophage

    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
    Calveolar_macrophage ~ prop(propSd_Calveolar_macrophage)
  })
}
