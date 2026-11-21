Vasalou_2024_trastuzumabDeruxtecan_mouse_qsp <- function() {
  description <- "Preclinical (mouse, nude and NOD-SCID tumor-bearing). QSP. Mechanistic plasma and tumor PK/PD model of trastuzumab deruxtecan (T-DXd) and its released payload DXd in four human xenografts of differing HER2 expression (NCI-N87, Capan-1, JIMT-1, MDA-MB-468). Total antibody follows a two-compartment model; intact T-DXd is total antibody scaled by a first-order-decaying drug-to-antibody ratio (DAR); released DXd is a one-compartment model fed by deconjugation and by ADC elimination. A Krogh-cylinder tumor disposition model (Shah et al.) carries antibody across the vessel wall and tumor surface, HER2 binding, internalization, non-specific uptake and lysosomal degradation, with DXd diffusing between cytosol and tumor interstitium and back to plasma. Tumor DXd (intracellular + extracellular) drives an indirect-response model of gamma-H2AX (fold of vehicle). Cell line enters through HER2 receptors per cell and the internalization rate. Deterministic: naive-pooled fit, no IIV and no residual error reported."
  reference <- "Vasalou C, Proia TA, Kazlauskas L, Przybyla A, Sung M, Mamidi S, Maratea K, Griffin M, Sargeant R, Urosevic J, Rosenbaum AI, Yuan J, Aluri KC, Ramsden D, Hariparsad N, Jones RDO, Mettetal JT. Quantitative evaluation of trastuzumab deruxtecan pharmacokinetics and pharmacodynamics in mouse models of varying degrees of HER2 expression. CPT Pharmacometrics Syst Pharmacol. 2024;13(6):994. doi:10.1002/psp4.13133. Equations 10-27 and Table S2 of Appendix S1; full-precision values from the NONMEM control stream Figure5D_TumorFits.txt in Appendix S3."
  vignette <- "Vasalou_2024_trastuzumabDeruxtecan"
  units <- list(
    time = "day",
    dosing = "nmol/kg",
    concentration = "nM"
  )
  # The IV bolus enters the total-antibody central state, which is not named
  # central or depot, so buildModelDb() cannot infer it.
  dosing <- "central_tab"

  # Drug-to-antibody-ratio state and the five Shah tumor-disposition states
  # (Appendix S1 Eqs 11, 20-24) have no canonical compartment.
  paper_specific_compartments <- c(
    "dar",
    "tumor_ext_tab",
    "tumor_bound_tab",
    "tumor_endo_tab",
    "tumor_intra_dxd",
    "tumor_ext_dxd",
    "gh2ax"
  )

  compartmentData <- list(
    central_tab = list(
      analyte = "trastuzumab deruxtecan (total antibody)",
      units = "nmol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_tab = list(
      analyte = "trastuzumab deruxtecan (total antibody)",
      units = "nmol/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_dxd = list(analyte = "DXd (released payload)", units = "nmol/kg", specimen = "plasma", verified = TRUE),
    dar = list(analyte = "drug-to-antibody ratio", units = "(ratio)", specimen = "not applicable", verified = TRUE),
    tumor_ext_tab = list(
      analyte = "trastuzumab deruxtecan (unbound, tumor interstitium)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_bound_tab = list(
      analyte = "trastuzumab deruxtecan-HER2 complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_endo_tab = list(
      analyte = "trastuzumab deruxtecan (endosomal/lysosomal)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_intra_dxd = list(
      analyte = "DXd (released payload, intracellular)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor_ext_dxd = list(
      analyte = "DXd (released payload, extracellular)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    gh2ax = list(
      analyte = "gamma-H2AX (fold of vehicle)",
      units = "(ratio)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CELLLINE_CAPAN1 = list(
      description = "Xenograft cell-line indicator: 1 = Capan-1 (human pancreatic adenocarcinoma) tumor, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Mutually exclusive with CELLLINE_JIMT1 and CELLLINE_MDAMB468; all three = 0 is the NCI-N87 reference. Selects 2e4 HER2 receptors per cell and an internalization rate of 2.77 /day (NONMEM TYPE = 2).",
      source_name = "TYPE"
    ),
    CELLLINE_JIMT1 = list(
      description = "Xenograft cell-line indicator: 1 = JIMT-1 (human breast carcinoma) tumor, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Mutually exclusive with CELLLINE_CAPAN1 and CELLLINE_MDAMB468. Selects 3.8e5 HER2 receptors per cell and the re-estimated internalization rate of 0.14 /day (NONMEM TYPE = 3).",
      source_name = "TYPE"
    ),
    CELLLINE_MDAMB468 = list(
      description = "Xenograft cell-line indicator: 1 = MDA-MB-468 (human triple-negative breast carcinoma) tumor, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Mutually exclusive with CELLLINE_CAPAN1 and CELLLINE_JIMT1. Selects 4830 HER2 receptors per cell and an internalization rate of 4.158 /day (NONMEM TYPE = 4).",
      source_name = "TYPE"
    )
  )

  population <- list(
    species = "mouse (male/female nude and NOD-SCID, human tumor xenografts)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5-8 weeks at study start",
    weight_range = "20 g assumed for tumor-volume scaling (Table S2)",
    disease_state = "Subcutaneous human tumor xenografts of differing HER2 expression, approximately 300 mm^3 at dosing: NCI-N87 gastric (male nude, IHC 3+, 3.5e6 receptors/cell), JIMT-1 breast (female nude, IHC 2+, 3.8e5), Capan-1 pancreatic (male NOD-SCID, IHC 2+, ~2e4), MDA-MB-468 breast (female NOD-SCID, IHC 0, 4800).",
    dose_range = "Single IV dose of 1 mg/kg (NCI-N87, MDA-MB-468) or 10 mg/kg (all four xenografts) trastuzumab deruxtecan",
    regions = "Preclinical (Axis Bioservices, UK)",
    notes = "Sparse terminal sampling: three mice per treatment group per time point up to 14 days. Plasma total antibody, intact T-DXd and released DXd, tumor released DXd, and tumor gamma-H2AX by western blot. Naive-pooled NONMEM 7.4 FOCE fit performed in sequence: DAR (kdec), plasma PK, tumor PK (KUP), then PD."
  )

  ini({
    # Plasma pharmacokinetics (Appendix S1 Eqs 10-16; Table S2). Full
    # precision from the fixed THETAs of the tumor-fit control stream
    # (Figure5D_TumorFits.txt), which carry forward the plasma-fit estimates.
    lcl_tab <- log(0.02378); label("Total antibody clearance CL_mAb (L/kg/day)") # Table S2: 0.023 (SEM 0.002); control stream THETA(1) 0.02378
    lvc_tab <- log(0.0818); label("Total antibody central volume Vc_mAb (L/kg)") # Table S2: 0.082 (SEM 0.01); control stream THETA(2) 0.0818
    lq_tab <- log(0.078); label("Total antibody intercompartmental clearance Q_mAb (L/kg/day)") # Table S2: 0.078 (SEM 0.04); control stream THETA(3) 0.078
    lvp_tab <- log(0.0816); label("Total antibody peripheral volume Vper_mAb (L/kg)") # Table S2: 0.082 (SEM 0.007); control stream THETA(4) 0.0816
    lkdec <- log(0.034219); label("DXd deconjugation rate kdec (1/day)") # Table S2: 0.034 (SEM 0.002); control stream THETA(5) 0.034219; Results 'kdec = 0.034 day-1'
    dar0 <- fixed(8); label("Nominal drug-to-antibody ratio at dosing DAR0 (unitless)") # Table S2: DAR0 = 8, Fixed
    lcl_dxd <- log(96.574); label("Released DXd clearance CL_Drug (L/kg/day)") # Table S2: 96.6 (SEM 15.5); control stream THETA(6) 96.574
    lvc_dxd <- log(43.522); label("Released DXd volume Vc_Drug (L/kg)") # Table S2: 43.5 (SEM 24.9); control stream THETA(7) 43.522

    # Tumor geometry and transport (Table S2; Shah et al. Krogh-cylinder model)
    tumor_vol_mm3 <- fixed(300); label("Tumor volume at dosing (mm^3)") # Table S2 R_tumor and TV notes: 'assuming a 300 mm3 tumor'; control stream TV_mm = 300
    bw <- fixed(0.02); label("Mouse body weight (kg)") # Table S2 TV note: 'a 20 g mouse'; control stream BW = 20/1000
    p_tab <- fixed(243.648); label("Antibody vascular permeability P_mAb (um/day)") # Table S2: 243.6, Fixed; control stream PADC = 243.6480
    d_tab <- fixed(0.011); label("Antibody diffusivity D_mAb (cm^2/day)") # Table S2: 0.011, Fixed
    eps_tab <- fixed(0.24); label("Tumor void fraction for antibody epsilon_mAb (fraction)") # Table S2: 0.24, Fixed
    rcap <- fixed(8); label("Capillary radius R_cap (um)") # Table S2: 8
    rkrogh <- fixed(75); label("Krogh cylinder radius R_Krogh (um)") # Table S2: 75
    p_dxd <- fixed(21000); label("DXd vascular permeability P_Drug (um/day)") # Table S2: 21000, Fixed
    d_dxd <- fixed(0.25); label("DXd diffusivity D_Drug (cm^2/day)") # Table S2: 0.25, Fixed
    eps_dxd <- fixed(0.44); label("Tumor void fraction for DXd epsilon_Drug (fraction)") # Table S2: 0.44, Fixed

    # HER2 binding and cellular processing (Table S2)
    lkoff <- fixed(log(0.336)); label("T-DXd/HER2 dissociation rate koff (1/day)") # Table S2: 0.336, Fixed (= 0.014 1/h)
    kd <- fixed(0.0477); label("T-DXd/HER2 equilibrium dissociation constant KD (nM)") # Table S2: 0.047, Fixed; control stream KD = 47.7/(10**9) in mmol/L = 0.0477 nM
    lkdeg <- fixed(log(7.2)); label("Endosomal/lysosomal ADC degradation rate kdeg (1/day)") # Table S2: 7.2 (= in-vitro 0.3 1/h)
    lkdiff <- fixed(log(16.8)); label("DXd influx = efflux rate across the tumor cell membrane kdiff (1/day)") # Table S2: 16.8 (= in-vitro kinD 0.7 1/h, Table S1)
    lkup <- log(0.51); label("Non-specific antibody uptake into tumor cells KUP (1/day)") # Table S2: 0.51 (SEM 0.05), model derived
    vcell <- fixed(3e-12); label("Volume of a single tumor cell (L)") # Table S2: 3*10^-12, Fixed

    # Cell-line-specific HER2 receptor density (Table S2 '# Receptors/cell')
    rpc_n87 <- fixed(3.5e6); label("HER2 receptors per cell, NCI-N87") # Table S2: 3.5*10^6
    rpc_capan1 <- fixed(2e4); label("HER2 receptors per cell, Capan-1") # Table S2: 2*10^4 (published studies)
    rpc_jimt1 <- fixed(3.8e5); label("HER2 receptors per cell, JIMT-1") # Table S2: 3.8*10^5
    rpc_mdamb468 <- fixed(4830); label("HER2 receptors per cell, MDA-MB-468") # control stream Number_Rec = 4830; Table S2 and Figure 1b: 4.8*10^3

    # Cell-line-specific internalization rate (Table S2 'Kint')
    lkint_n87 <- fixed(log(0.36)); label("HER2/ADC internalization rate, NCI-N87 (1/day)") # Table S2: 0.36 (t1/2 46 h from in-vitro modeling)
    lkint_capan1 <- fixed(log(2.77)); label("HER2/ADC internalization rate, Capan-1 (1/day)") # Table S2: 2.8 (t1/2 6 h); control stream THETA(13) 2.77
    lkint_jimt1 <- log(0.14); label("HER2/ADC internalization rate, JIMT-1, re-estimated (1/day)") # Table S2 footnote: 0.14 day-1 (t1/2 ~120 h); literature value 5.5 (control stream THETA(14) 5.54) was replaced
    lkint_mdamb468 <- fixed(log(4.158)); label("HER2/ADC internalization rate, MDA-MB-468 (1/day)") # Table S2: 4.2 (t1/2 4 h); control stream THETA(15) 4.158

    # gamma-H2AX indirect-response PD (Appendix S1 Eq 27; Table S2)
    lkin <- log(2); label("Zero-order gamma-H2AX production rate kinp (1/day); koutp = kinp") # Table S2: kinp 2 (SEM 1.2), estimated; koutp '= kinp'
    emax <- fixed(5); label("Maximal stimulation of gamma-H2AX production PDmax (unitless)") # Table S2: PDmax 5, 'Fixed based on maximal observable response'
    lic50 <- log(96.6); label("Tumor DXd concentration at half-maximal stimulation IC50_PD (nM)") # Table S2: 96.6 (SEM 24.2), estimated
  })

  model({
    cl_tab <- exp(lcl_tab)
    vc_tab <- exp(lvc_tab)
    q_tab <- exp(lq_tab)
    vp_tab <- exp(lvp_tab)
    kdec <- exp(lkdec)
    cl_dxd <- exp(lcl_dxd)
    vc_dxd <- exp(lvc_dxd)
    koff <- exp(lkoff)
    kon <- koff / kd
    kdeg <- exp(lkdeg)
    kdiff <- exp(lkdiff)
    kup <- exp(lkup)
    kin <- exp(lkin)
    kout <- kin
    ic50 <- exp(lic50)

    # Cell line selection (control stream TYPE 1-4); NCI-N87 is the reference.
    f_n87 <- 1 - CELLLINE_CAPAN1 - CELLLINE_JIMT1 - CELLLINE_MDAMB468
    rpc <- rpc_n87 * f_n87 + rpc_capan1 * CELLLINE_CAPAN1 +
      rpc_jimt1 * CELLLINE_JIMT1 + rpc_mdamb468 * CELLLINE_MDAMB468
    kint <- exp(lkint_n87) * f_n87 + exp(lkint_capan1) * CELLLINE_CAPAN1 +
      exp(lkint_jimt1) * CELLLINE_JIMT1 + exp(lkint_mdamb468) * CELLLINE_MDAMB468

    # Total HER2 in tumor (Eq 26): receptors/cell / (Vcell * Avogadro), in nM
    ag <- rpc / (6.023e23 * vcell) * 1e9

    # Tumor geometry: spherical radius (cm) and tumor volume per body weight (L/kg)
    rtumor <- (tumor_vol_mm3 / ((4 / 3) * pi))^(1 / 3) / 10
    tvbw <- tumor_vol_mm3 * 1e-6 / bw

    # Exchange rate constants (1/day): vessel-wall permeation (um/day, um) plus
    # surface diffusion (cm^2/day, cm)
    kexch_tab <- 2 * p_tab * rcap / rkrogh^2 + 6 * d_tab / rtumor^2
    kexch_dxd <- 2 * p_dxd * rcap / rkrogh^2 + 6 * d_dxd / rtumor^2

    Cc_tab <- central_tab / vc_tab
    Cc_dxd <- central_dxd / vc_dxd

    dar(0) <- dar0
    gh2ax(0) <- kin / kout

    # Plasma (Eqs 17-19, 11)
    d/dt(central_tab) <- -cl_tab / vc_tab * central_tab -
      q_tab * (Cc_tab - peripheral1_tab / vp_tab) -
      tvbw * kexch_tab * (Cc_tab - tumor_ext_tab / eps_tab)
    d/dt(peripheral1_tab) <- q_tab * (Cc_tab - peripheral1_tab / vp_tab)
    d/dt(central_dxd) <- -cl_dxd / vc_dxd * central_dxd +
      (kdec + cl_tab / vc_tab) * dar * central_tab -
      tvbw * kexch_dxd * (Cc_dxd - tumor_ext_dxd / eps_dxd)
    d/dt(dar) <- -kdec * dar

    # Tumor (Eqs 20-24); states are tumor concentrations (nM)
    d/dt(tumor_ext_tab) <- kexch_tab * (Cc_tab - tumor_ext_tab / eps_tab) -
      kon * tumor_ext_tab * (ag - tumor_bound_tab) / eps_tab +
      koff * tumor_bound_tab - kup * tumor_ext_tab
    d/dt(tumor_bound_tab) <- kon * tumor_ext_tab * (ag - tumor_bound_tab) / eps_tab -
      koff * tumor_bound_tab - kint * tumor_bound_tab
    d/dt(tumor_endo_tab) <- kint * tumor_bound_tab + kup * tumor_ext_tab -
      kdeg * tumor_endo_tab
    d/dt(tumor_intra_dxd) <- kdeg * tumor_endo_tab * dar -
      kdiff * (tumor_intra_dxd - tumor_ext_dxd)
    d/dt(tumor_ext_dxd) <- kdiff * (tumor_intra_dxd - tumor_ext_dxd) +
      kdec * dar * (tumor_bound_tab + tumor_ext_tab) +
      kexch_dxd * (Cc_dxd - tumor_ext_dxd / eps_dxd)

    # Released DXd in tumor (Eq 25)
    Ctumor_dxd <- tumor_intra_dxd + tumor_ext_dxd

    # gamma-H2AX indirect response (Eq 27)
    d/dt(gh2ax) <- kin * (1 + emax * Ctumor_dxd / (ic50 + Ctumor_dxd)) - kout * gh2ax
    gh2ax_pct <- 100 * gh2ax

    # Intact T-DXd (Eq 16): total antibody scaled by DAR / DAR0
    Cc <- Cc_tab * dar / dar0
  })
}
