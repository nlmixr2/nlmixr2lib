Tosca_2023_men1611_trastuzumab_hcc1954_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; HCC1954 cell-line-derived HER2-positive,",
    "trastuzumab-resistant breast-cancer xenograft, TGI study C). PK-PD",
    "tumor-growth-inhibition model for MEN1611 (oral PI3K inhibitor) combined",
    "with trastuzumab: the mouse two-compartment PK models of both drugs",
    "(held fixed) drive a Simeoni-type model with exponential-only (reduced",
    "Simeoni) unperturbed growth and a four-by-four grid of damaged-cell",
    "transit compartments, one axis per drug (synergistic 2D-TGI model, in",
    "which the trastuzumab potency is increased by a term proportional to the",
    "MEN1611 concentration). Typical-value fit to mean tumor weight with no",
    "between-animal variability."
  )
  reference <- paste(
    "Tosca EM, Borella E, Piana C, Bouchene S, Merlino G, Fiascarelli A,",
    "Mazzei P, Magni P. Model-based prediction of",
    "effective target exposure for MEN1611 in combination with trastuzumab in",
    "HER2-positive advanced or metastatic breast cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1626-1639.",
    "doi:10.1002/psp4.12910. Tumor parameters from Table 2; PK parameters from",
    "Supplementary Tables S3 (trastuzumab) and S4 (MEN1611); equations from",
    "Equations 2-3 and Supplementary Material S3. Growth model from Simeoni M,",
    "et al. Cancer Res. 2004;64(3):1094-1101. doi:10.1158/0008-5472.CAN-03-2524.",
    sep = " "
  )
  vignette <- "Tosca_2023_men1611_trastuzumab"

  units <- list(
    time = "day",
    dosing = "MEN1611 in mg/kg (depot_men1611); trastuzumab in mg per mouse (depot_trastuzumab IP, central_trastuzumab IV)",
    concentration = "ng/mL"
  )

  covariateData <- list()

  # The 2D damage grid of Tosca 2023 Material S3: damaged_cells_m<i>_t<j> is
  # paper compartment X_(i+1)(j+1), a tumor cell that has progressed i stages
  # along the MEN1611 mortality chain and j stages along the trastuzumab chain.
  paper_specific_compartments <- c(
    "damaged_cells_m0_t1",
    "damaged_cells_m0_t2",
    "damaged_cells_m0_t3",
    "damaged_cells_m1_t0",
    "damaged_cells_m1_t1",
    "damaged_cells_m1_t2",
    "damaged_cells_m1_t3",
    "damaged_cells_m2_t0",
    "damaged_cells_m2_t1",
    "damaged_cells_m2_t2",
    "damaged_cells_m2_t3",
    "damaged_cells_m3_t0",
    "damaged_cells_m3_t1",
    "damaged_cells_m3_t2",
    "damaged_cells_m3_t3"
  )

  compartmentData <- list(
    depot_men1611 = list(analyte = "MEN1611", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_men1611 = list(analyte = "MEN1611", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_men1611 = list(analyte = "MEN1611", units = "mg/kg", specimen = "tissue", verified = TRUE),
    depot_trastuzumab = list(analyte = "trastuzumab", units = "mg", specimen = "administration site", verified = TRUE),
    central_trastuzumab = list(analyte = "trastuzumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_trastuzumab = list(analyte = "trastuzumab", units = "mg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m0_t1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m0_t2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m0_t3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m1_t0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m1_t1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m1_t2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m1_t3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m2_t0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m2_t1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m2_t2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m2_t3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m3_t0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m3_t1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m3_t2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_m3_t3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1 bearing HCC1954 xenografts)",
    n_subjects = 21L,
    n_studies = 1L,
    disease_state = "subcutaneous HCC1954 cell-line-derived human breast-cancer xenograft: HER2-amplified, trastuzumab-resistant, PIK3CA H1047R mutation",
    dose_range = "MEN1611 6.5 mg/kg oral once daily for 12 days; trastuzumab 30 mg/kg IP once weekly for 2 weeks; given alone and in combination (Tosca 2023 Table 1)",
    regions = "preclinical",
    notes = paste(
      "Tosca 2023 TGI study C (Table 1): vehicle, MEN1611, trastuzumab and",
      "combination arms (5 vehicle, 5 MEN1611, 5 trastuzumab, 6 combination",
      "mice). Time zero is the start of the tumor-growth record, before",
      "treatment began; tumor weight was computed from caliper length and",
      "width assuming an oblate ellipsoid of unit density (Equation 1). The",
      "tumors did not reach the linear growth phase, so the reduced",
      "(exponential-only) Simeoni growth model was used. Control and",
      "single-agent arms were fitted simultaneously. The synergy factor was",
      "then estimated on the combination arm with every other parameter held",
      "at its single-agent estimate. Tosca 2023 Table 2 prints the line as",
      "'HCC194'; Table 1 names it HCC1954."
    )
  )

  ini({
    # --- MEN1611 mouse PK: Tosca 2023 Supplementary Table S4 (held fixed) ---
    # Printed in hours; model() multiplies rate constants and clearances by 24
    # for the day time base. CL and Q are per-kg (L/h/kg; see the vignette).
    lka_men1611 <- fixed(log(1.84)); label("MEN1611 absorption rate constant (1/h)") # Table S4: Kabs = 1.84 1/h
    lcl_men1611 <- fixed(log(1)); label("MEN1611 clearance (L/h/kg)") # Table S4: CL = 1
    lvc_men1611 <- fixed(log(0.997)); label("MEN1611 central volume (L/kg)") # Table S4: V1 = 0.997 L/kg
    lq_men1611 <- fixed(log(1.49)); label("MEN1611 intercompartmental clearance (L/h/kg)") # Table S4: Q = 1.49
    lvp_men1611 <- fixed(log(7.6)); label("MEN1611 peripheral volume (L/kg)") # Table S4: V2 = 7.6 L/kg
    # The xenograft vehicle was DMSO/Cremophor EL; the solution bioavailability
    # (fixed to 1) reproduces the MEN1611 peak concentrations implied by
    # Table S9, while the suspension value 0.559 does not (see the vignette).
    lfdepot_men1611 <- fixed(log(1)); label("MEN1611 oral bioavailability of the solution (fraction)") # Table S4: F solution = 1 (fixed)

    # --- Trastuzumab mouse PK: Tosca 2023 Supplementary Table S3 (held fixed) ---
    lka_trastuzumab <- fixed(log(0.151)); label("Trastuzumab IP absorption rate constant (1/h)") # Table S3: Kabs = 0.151 1/h
    lcl_trastuzumab <- fixed(log(1.53e-5)); label("Trastuzumab clearance (L/h)") # Table S3: CL = 1.53E-5 L/h
    lvc_trastuzumab <- fixed(log(2.45e-3)); label("Trastuzumab central volume (L)") # Table S3: V1 = 2.45E-3 L
    lq_trastuzumab <- fixed(log(2.58e-4)); label("Trastuzumab intercompartmental clearance (L/h)") # Table S3: Q = 2.58e-4 L/h
    lvp_trastuzumab <- fixed(log(2.22e-3)); label("Trastuzumab peripheral volume (L)") # Table S3: V2 = 2.22e-3 L
    lfdepot_trastuzumab <- fixed(log(0.5557)); label("Trastuzumab IP bioavailability (fraction)") # Table S3: F = 0.5557

    # --- Tumor growth and drug effects: Tosca 2023 Table 2, Study C (HCC1954) ---
    ltumorExpGrowth <- log(0.038); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Table 2 Study C: lambda0 = 0.038 1/day (19%)
    lrbase_tumor <- log(0.137); label("Initial tumor weight W0 (g)") # Table 2 Study C: W0 = 0.137 g (19%)
    ldamageTransit_trastuzumab <- log(84.5); label("Trastuzumab damaged-cell transit rate constant k1,TZB (1/day)") # Table 2 Study C: k1,TZB = 84.5 1/day (>100%)
    ldrugSlope_trastuzumab <- log(1.2e-8); label("Trastuzumab potency k2,TZB (mL/(ng*day))") # Table 2 Study C: k2,TZB = 1.2 E-8 ml*day/ng (>100%)
    ldamageTransit_men1611 <- log(55.1); label("MEN1611 damaged-cell transit rate constant k1,MEN1611 (1/day)") # Table 2 Study C: k1,MEN1611 = 55.1 1/day (>100%)
    ldrugSlope_men1611 <- log(1.13e-4); label("MEN1611 potency k2,MEN1611 (mL/(ng*day))") # Table 2 Study C: k2,MEN1611 = 1.13 E-4 ml*day/ng (>100%)
    # Table 2 labels gamma '(-)', but Equation 3 multiplies it by the MEN1611
    # concentration, so its unit is the reciprocal of that concentration.
    lsynergy <- log(0.263); label("Synergy factor gamma: increase in trastuzumab potency per unit MEN1611 concentration (mL/ng)") # Table 2 Study C: gamma = 0.263 (23%)

    # Constant residual error y = f + a*eps (Tosca 2023 Methods); the value of a
    # is not reported, so it is held at zero.
    addSd <- fixed(0); label("Additive residual error on tumor weight (g; not reported in the source)")
  })

  model({
    # --- PK (Tosca 2023 Supplementary Material S5); x 24 converts 1/h and L/h to per day ---
    ka_men1611 <- exp(lka_men1611) * 24
    cl_men1611 <- exp(lcl_men1611) * 24
    vc_men1611 <- exp(lvc_men1611)
    q_men1611 <- exp(lq_men1611) * 24
    vp_men1611 <- exp(lvp_men1611)
    ka_trastuzumab <- exp(lka_trastuzumab) * 24
    cl_trastuzumab <- exp(lcl_trastuzumab) * 24
    vc_trastuzumab <- exp(lvc_trastuzumab)
    q_trastuzumab <- exp(lq_trastuzumab) * 24
    vp_trastuzumab <- exp(lvp_trastuzumab)

    d/dt(depot_men1611) <- -ka_men1611 * depot_men1611
    d/dt(central_men1611) <- ka_men1611 * depot_men1611 - (cl_men1611 + q_men1611) / vc_men1611 * central_men1611 +
      q_men1611 / vp_men1611 * peripheral1_men1611
    d/dt(peripheral1_men1611) <- q_men1611 / vc_men1611 * central_men1611 - q_men1611 / vp_men1611 * peripheral1_men1611
    f(depot_men1611) <- exp(lfdepot_men1611)

    # Trastuzumab IP doses go to depot_trastuzumab, IV doses to central_trastuzumab.
    d/dt(depot_trastuzumab) <- -ka_trastuzumab * depot_trastuzumab
    d/dt(central_trastuzumab) <- ka_trastuzumab * depot_trastuzumab -
      (cl_trastuzumab + q_trastuzumab) / vc_trastuzumab * central_trastuzumab +
      q_trastuzumab / vp_trastuzumab * peripheral1_trastuzumab
    d/dt(peripheral1_trastuzumab) <- q_trastuzumab / vc_trastuzumab * central_trastuzumab -
      q_trastuzumab / vp_trastuzumab * peripheral1_trastuzumab
    f(depot_trastuzumab) <- exp(lfdepot_trastuzumab)

    # mg/kg / (L/kg) and mg / L are both mg/L; x 1000 gives ng/mL.
    Cc_men1611 <- central_men1611 / vc_men1611 * 1000
    Cc_trastuzumab <- central_trastuzumab / vc_trastuzumab * 1000

    # --- Tumor growth inhibition ---
    tumorExpGrowth <- exp(ltumorExpGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit_trastuzumab <- exp(ldamageTransit_trastuzumab)
    drugSlope_trastuzumab <- exp(ldrugSlope_trastuzumab)
    damageTransit_men1611 <- exp(ldamageTransit_men1611)
    drugSlope_men1611 <- exp(ldrugSlope_men1611)
    synergy <- exp(lsynergy)

    # Equation 3: k2,TZB,combo = k2,TZB + k2,TZB * gamma * C_MEN1611. The
    # synergistic model is the additive grid with k2,TZB replaced by k2,TZB,combo.
    drugSlope_trastuzumab_combo <- drugSlope_trastuzumab * (1 + synergy * Cc_men1611)
    killRate_trastuzumab <- drugSlope_trastuzumab_combo * Cc_trastuzumab
    killRate_men1611 <- drugSlope_men1611 * Cc_men1611

    tumor_size <- cycling_cells + damaged_cells_m0_t1 + damaged_cells_m0_t2 +
      damaged_cells_m0_t3 + damaged_cells_m1_t0 + damaged_cells_m1_t1 +
      damaged_cells_m1_t2 + damaged_cells_m1_t3 + damaged_cells_m2_t0 +
      damaged_cells_m2_t1 + damaged_cells_m2_t2 + damaged_cells_m2_t3 +
      damaged_cells_m3_t0 + damaged_cells_m3_t1 + damaged_cells_m3_t2 +
      damaged_cells_m3_t3

    # Reduced Simeoni growth: exponential only (Tosca 2023 Methods).
    growth <- tumorExpGrowth * cycling_cells

    # Material S3 mortality grid: a cell can be hit by either drug at any stage
    # of the other drug's chain, so damage progresses independently along both
    # axes; cells leaving the last stage of either chain die.
    d/dt(cycling_cells) <- growth - (killRate_men1611 + killRate_trastuzumab) * cycling_cells
    d/dt(damaged_cells_m0_t1) <- killRate_trastuzumab * cycling_cells - (killRate_men1611 + damageTransit_trastuzumab) * damaged_cells_m0_t1
    d/dt(damaged_cells_m0_t2) <- damageTransit_trastuzumab * damaged_cells_m0_t1 - (killRate_men1611 + damageTransit_trastuzumab) * damaged_cells_m0_t2
    d/dt(damaged_cells_m0_t3) <- damageTransit_trastuzumab * damaged_cells_m0_t2 - (killRate_men1611 + damageTransit_trastuzumab) * damaged_cells_m0_t3
    d/dt(damaged_cells_m1_t0) <- killRate_men1611 * cycling_cells - (damageTransit_men1611 + killRate_trastuzumab) * damaged_cells_m1_t0
    d/dt(damaged_cells_m1_t1) <- killRate_men1611 * damaged_cells_m0_t1 + killRate_trastuzumab * damaged_cells_m1_t0 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m1_t1
    d/dt(damaged_cells_m1_t2) <- killRate_men1611 * damaged_cells_m0_t2 + damageTransit_trastuzumab * damaged_cells_m1_t1 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m1_t2
    d/dt(damaged_cells_m1_t3) <- killRate_men1611 * damaged_cells_m0_t3 + damageTransit_trastuzumab * damaged_cells_m1_t2 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m1_t3
    d/dt(damaged_cells_m2_t0) <- damageTransit_men1611 * damaged_cells_m1_t0 - (damageTransit_men1611 + killRate_trastuzumab) * damaged_cells_m2_t0
    d/dt(damaged_cells_m2_t1) <- damageTransit_men1611 * damaged_cells_m1_t1 + killRate_trastuzumab * damaged_cells_m2_t0 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m2_t1
    d/dt(damaged_cells_m2_t2) <- damageTransit_men1611 * damaged_cells_m1_t2 + damageTransit_trastuzumab * damaged_cells_m2_t1 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m2_t2
    d/dt(damaged_cells_m2_t3) <- damageTransit_men1611 * damaged_cells_m1_t3 + damageTransit_trastuzumab * damaged_cells_m2_t2 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m2_t3
    d/dt(damaged_cells_m3_t0) <- damageTransit_men1611 * damaged_cells_m2_t0 - (damageTransit_men1611 + killRate_trastuzumab) * damaged_cells_m3_t0
    d/dt(damaged_cells_m3_t1) <- damageTransit_men1611 * damaged_cells_m2_t1 + killRate_trastuzumab * damaged_cells_m3_t0 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m3_t1
    d/dt(damaged_cells_m3_t2) <- damageTransit_men1611 * damaged_cells_m2_t2 + damageTransit_trastuzumab * damaged_cells_m3_t1 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m3_t2
    d/dt(damaged_cells_m3_t3) <- damageTransit_men1611 * damaged_cells_m2_t3 + damageTransit_trastuzumab * damaged_cells_m3_t2 - (damageTransit_men1611 + damageTransit_trastuzumab) * damaged_cells_m3_t3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ add(addSd)
  })
}
