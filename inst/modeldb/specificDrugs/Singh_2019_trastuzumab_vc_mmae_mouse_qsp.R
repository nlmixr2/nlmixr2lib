Singh_2019_trastuzumab_vc_mmae_mouse_qsp <- function() {
  description <- "QSP (preclinical, mouse). Cell-level systems PK-PD model for the HER2-targeting antibody-drug conjugate trastuzumab-vc-MMAE (T-vc-MMAE, payload MMAE) in NCI-N87 (HER2-high) and GFP-MCF7 (HER2-low) xenograft-bearing SCID mice, fit jointly to both tumor lines. Plasma: two-compartment ADC with linear catabolic clearance and first-order non-specific deconjugation (declining average DAR), plus two-compartment unconjugated MMAE formed by both pathways; plasma total trastuzumab adds the deconjugated antibody, cleared like the ADC. Tumor: Krogh-cylinder vascular + surface exchange of ADC and MMAE into the tumor extracellular space, then a single-cell model (HER2 binding, internalization, lysosomal degradation, MMAE influx/efflux and tubulin binding, dilution by cell growth). Intracellular tubulin occupancy drives a sigmoid kill function that moves viable tumor volume through a three-compartment death transit chain; dying cells release their intracellular ADC and MMAE back to the tumor extracellular space. Cell-line-specific HER2 count, doubling time and cell volume are switched by CELLLINE_NCIN87. Time in days, ADC dose in nmol."
  reference <- paste(
    "Singh AP, Guo L, Verma A, Wong GG-L, Shah DK. A Cell-Level Systems PK-PD",
    "Model to Characterize In Vivo Efficacy of ADCs. Pharmaceutics.",
    "2019;11(2):98. doi:10.3390/pharmaceutics11020098.",
    "Model structure: Figure 1 and Equations (1)-(23); parameters: Table 1.",
    "Cell volumes and the molecule-to-nanomole scaling factor SF, which the",
    "2019 paper uses but does not print, are from its single-cell source:",
    "Singh AP, Shah DK. Drug Metab Dispos. 2017;45(11):1120-1132",
    "(doi:10.1124/dmd.117.076414), Table 2.",
    sep = " "
  )
  vignette <- "Singh_2019_trastuzumab_vc_mmae"
  units <- list(
    time = "day",
    dosing = "nmol (T-vc-MMAE, antibody basis)",
    concentration = "nM"
  )

  # What each ODE state holds. Plasma states are whole-mouse amounts (nmol) or
  # plasma concentrations (nM, as the paper writes Eqs 5-6); tumor extracellular
  # ADC is a concentration over the whole tumor volume (nM), tumor extracellular
  # MMAE an amount (nmol); the four cell states are molecules per single tumor
  # cell (Eqs 10-13). verified = TRUE: matched to Figure 1 and Eqs (3)-(23).
  compartmentData <- list(
    central = list(
      analyte = "T-vc-MMAE (ADC), amount",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "T-vc-MMAE (ADC), amount",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    mab_central = list(
      analyte = "deconjugated (naked) trastuzumab formed in plasma, amount",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    mab_peripheral = list(
      analyte = "deconjugated (naked) trastuzumab, amount",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    mmae_central = list(
      analyte = "unconjugated MMAE, concentration",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    mmae_peripheral = list(
      analyte = "unconjugated MMAE, concentration",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    dar = list(
      analyte = "average drug-to-antibody ratio of circulating ADC",
      units = "(unitless)",
      specimen = "plasma",
      verified = TRUE
    ),
    adc_ex_tumor = list(
      analyte = "T-vc-MMAE in tumor extracellular space, concentration over total tumor volume",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    mmae_ex_tumor = list(
      analyte = "unconjugated MMAE in tumor extracellular space, amount",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_bound_cell = list(
      analyte = "HER2-bound T-vc-MMAE per tumor cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    adc_lyso_cell = list(
      analyte = "internalized (endosomal/lysosomal) T-vc-MMAE per tumor cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    mmae_free_cell = list(
      analyte = "unbound cytoplasmic MMAE per tumor cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    mmae_tub_cell = list(
      analyte = "tubulin-bound MMAE per tumor cell",
      units = "molecules/cell",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor1 = list(
      analyte = "proliferating tumor volume",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor2 = list(
      analyte = "non-growing (dying) tumor volume, transit 1",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor3 = list(
      analyte = "non-growing (dying) tumor volume, transit 2",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    tumor4 = list(
      analyte = "non-growing (dying) tumor volume, transit 3",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    CELLLINE_NCIN87 = list(
      description = "Tumor cell-line indicator: 1 = NCI-N87 (HER2-high gastric carcinoma) xenograft, 0 = GFP-MCF7 (HER2-low breast carcinoma) xenograft",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = GFP-MCF7",
      notes = "Selects the cell-line-specific HER2 receptor count (Table 1 AgN87ex / AgMCF7ex), tumor doubling time (DTN87 / DTMCF7) and single-cell volume (Singh 2017 Table 2). Kill parameters are shared by both lines.",
      source_name = "tumor model (Group A GFP-MCF7 / Group B N87, Section 2.6)"
    ),
    TUM_VOL = list(
      description = "Tumor volume at the time of the ADC dose; sets the proliferating tumor state tumor1(0)",
      units = "mm^3",
      type = "continuous",
      source_name = "TVmm3(0), initial condition of Eq (1); ~500 mm^3 in the PK study, ~350-400 mm^3 in the TGI study (Sections 2.5, 2.6)"
    ),
    WT = list(
      description = "Body weight; converts the per-kg plasma clearances and volumes (Table 1) to whole-mouse values",
      units = "kg",
      type = "continuous",
      source_name = "not reported; Table 1 PK parameters are per kg and Eq (3) exchanges whole-tumor amounts with the central compartment"
    )
  )

  covariatesDataExcluded <- list()

  population <- list(
    species = "mouse (male SCID, NOD.CB17-Prkdcscid/J; NCI-N87 or GFP-MCF7 subcutaneous xenograft)",
    n_subjects = 74L,
    n_studies = 2L,
    disease_state = "Subcutaneous flank xenografts of NCI-N87 human gastric carcinoma (~950,000 HER2/cell in vitro) or GFP-transfected MCF7 human breast carcinoma (~55,000 HER2/cell in vitro); GFP-MCF7 mice supplemented with 17beta-estradiol valerate",
    age_range = "8 weeks at tumor implantation (6 weeks on arrival + 2 weeks acclimation); dosed 4-6 weeks later",
    weight_range = "not reported",
    sex_female_pct = 0,
    dose_range = "Single IV dose. Tumor PK: 10 mg/kg (n = 9 per tumor line). TGI: GFP-MCF7 0 (control), 3, 5, 10 mg/kg; N87 0 (control), 1, 3, 10 mg/kg (n = 7 per arm)",
    regions = "Preclinical (University at Buffalo, USA)",
    notes = paste(
      "18 mice for tumor PK (sparse terminal sampling at 24, 72, 168 h plus a",
      "10-min retro-orbital plasma sample) and 56 for TGI. Average DAR ~4",
      "(random cysteine conjugation). Plasma PK fit in ADAPT-5 (ML, additive +",
      "proportional variance); MMAE plasma disposition fixed from Shah 2012;",
      "single-cell constants fixed from Singh & Shah 2017; only HER2 counts",
      "estimated from tumor PK; kill parameters (Kmax, KC50, tau, gamma) and",
      "doubling times estimated from TGI in Monolix 8 (SAEM) with log-normal",
      "IIV on Kmax and tau."
    )
  )

  # Every state except central / peripheral1 is a mechanistic, model-specific
  # pool; whitelist them so checkModelConventions() does not flag them.
  paper_specific_compartments <- c(
    "mab_central",
    "mab_peripheral",
    "mmae_central",
    "mmae_peripheral",
    "dar",
    "adc_ex_tumor",
    "mmae_ex_tumor",
    "adc_bound_cell",
    "adc_lyso_cell",
    "mmae_free_cell",
    "mmae_tub_cell",
    "tumor1",
    "tumor2",
    "tumor3",
    "tumor4"
  )

  ini({
    # ---- Plasma PK of T-vc-MMAE (Table 1; per kg body weight) ----
    lcl <- log(0.033); label("ADC catabolic clearance CLADC (L/day/kg)") # Table 1: 0.033 (4.8%)
    lq <- log(0.0585); label("ADC distributional clearance CLDADC (L/day/kg)") # Table 1: 0.0585 (12.6%)
    lvc <- log(0.084); label("ADC central volume V1ADC (L/kg)") # Table 1: 0.084 (7.3%)
    lvp <- log(0.051); label("ADC peripheral volume V2ADC (L/kg)") # Table 1: 0.051 (5.2%)
    lkdec <- log(0.323); label("Non-specific deconjugation rate KdecP (1/day)") # Table 1: 0.323 (8.8%)
    dar0 <- fixed(4); label("Initial average drug-to-antibody ratio DAR(0) (unitless)") # Section 2.2: 'average drug: antibody ratio (DAR) of ~4'

    # ---- Plasma PK of unconjugated MMAE (Table 1, fixed from ref [11]) ----
    lcl_mmae <- fixed(log(18.40)); label("MMAE systemic clearance CLDrug (L/day/kg)") # Table 1: 18.40 [11]
    lq_mmae <- fixed(log(1.84)); label("MMAE distributional clearance CLDDrug (L/day/kg)") # Table 1: 1.84 [11]
    lvc_mmae <- fixed(log(0.136)); label("MMAE central volume V1Drug (L/kg)") # Table 1: 0.136 [11]
    lvp_mmae <- fixed(log(0.523)); label("MMAE peripheral volume V2Drug (L/kg)") # Table 1: 0.523 [11]

    # ---- Tumor distribution (Table 1, literature [18-20]) ----
    rcap <- fixed(8.0); label("Tumor capillary radius RCap (um)") # Table 1: 8.0
    rkrogh <- fixed(75.0); label("Krogh radius, mean inter-capillary distance RKrogh (um)") # Table 1: 75.0
    perm <- fixed(334); label("Vascular permeability of T-vc-MMAE PADC (um/day)") # Table 1: 334
    perm_mmae <- fixed(21000); label("Vascular permeability of MMAE PDrug (um/day)") # Table 1: 21000
    dif <- fixed(0.022); label("Diffusivity of T-vc-MMAE DADC (cm^2/day)") # Table 1: 0.022
    dif_mmae <- fixed(0.25); label("Diffusivity of MMAE DDrug (cm^2/day)") # Table 1: 0.25
    eps <- fixed(0.24); label("Tumor void fraction for T-vc-MMAE epsilonADC (unitless)") # Table 1: 0.24
    eps_mmae <- fixed(0.44); label("Tumor void fraction for MMAE epsilonDrug (unitless)") # Table 1: 0.44

    # ---- Single-cell disposition (Table 1, fixed from refs [11,14]; per hour) ----
    kon <- fixed(0.03); label("ADC-HER2 association rate KonADC (1/nM/h)") # Table 1: 0.03 [14]
    koff <- fixed(0.014); label("ADC-HER2 dissociation rate KoffADC (1/h)") # Table 1: 0.014 [14]
    kint <- fixed(0.11); label("HER2-ADC complex internalization rate KintADC (1/h)") # Table 1: 0.11 [14]
    kdeg <- fixed(0.353); label("Lysosomal ADC degradation rate KdegADC (1/h)") # Table 1: 0.353 [14]
    kon_tub <- fixed(0.0183); label("MMAE-tubulin association rate KonTub (1/nM/h)") # Table 1: 0.0183 [14]
    koff_tub <- fixed(0.545); label("MMAE-tubulin dissociation rate KoffTub (1/h)") # Table 1: 0.545 [14]
    tub_tot <- fixed(65); label("Total intracellular tubulin concentration Tubtot (nM)") # Table 1: 65 [11,14]
    kin_mmae <- fixed(8.33); label("MMAE influx rate into cells KinDrug (1/h)") # Table 1: 8.33 [14]
    kout_mmae <- fixed(0.046); label("MMAE efflux rate out of cells KoutDrug (1/h)") # Table 1: 0.046 [11]
    vcell_n87 <- fixed(3.12); label("Volume of one NCI-N87 cell VCell (pL)") # Singh 2017 DMD Table 2: VCell N87 = 3.12 pL (not printed in Singh 2019)
    vcell_mcf7 <- fixed(8.14); label("Volume of one GFP-MCF7 cell VCell (pL)") # Singh 2017 DMD Table 2: VCell MCF7 = 8.14 pL (not printed in Singh 2019)

    # ---- In-vivo HER2 receptor counts (estimated from tumor PK) ----
    lag_n87 <- log(185000); label("HER2 receptors per N87 tumor cell AgN87ex (molecules/cell)") # Table 1: 185,000 (2.8%)
    lag_mcf7 <- log(22400); label("HER2 receptors per GFP-MCF7 tumor cell AgMCF7ex (molecules/cell)") # Table 1: 22,400 (3.2%)

    # ---- Tumor growth and tubulin-occupancy-driven kill (estimated from TGI) ----
    ldt_n87 <- log(13.5); label("Doubling time of N87 tumors DTN87 (day)") # Table 1: 13.5 (11.4%)
    ldt_mcf7 <- log(10.6); label("Doubling time of GFP-MCF7 tumors DTMCF7 (day)") # Table 1: 10.6 (18.7%)
    lkmax <- log(1.03); label("Maximum kill rate Kmax (1/day)") # Table 1: 1.03 (31.3%)
    lkc50 <- log(96.8); label("Tubulin occupancy giving half-maximal kill KC50 (%)") # Table 1: 96.8 (13.2%)
    ltau <- log(2.03); label("Transit time of the death chain tau (day)") # Table 1: 2.03
    lhill <- log(15.02); label("Sigmoidicity of the occupancy-kill relationship gamma (unitless)") # Table 1: 15.02 (38.6%)

    # Between-mouse variability (Monolix log-normal; Table 1 reports the omega
    # as a percentage with its RSE in parentheses).
    etalkmax ~ 0.01032 # Table 1: IIVKmax 10.16 (47%) -> omega 0.1016, omega^2 = 0.01032
    etaltau ~ 0.03764 # Table 1: IIVtau 19.4 (32%) -> omega 0.194, omega^2 = 0.03764

    # Residual error magnitudes are not reported (Eq 24 gives only the ADAPT
    # variance-model form; the Monolix TGI error model is not described).
    propSd <- fixed(0); label("Proportional residual SD on plasma total trastuzumab (fraction); not reported")
    propSd_tumor_volume <- fixed(0); label("Proportional residual SD on tumor volume (fraction); not reported")
  })

  model({
    # ---- whole-mouse plasma PK parameters (per-kg Table 1 values x WT) ----
    cl <- exp(lcl) * WT
    q <- exp(lq) * WT
    vc <- exp(lvc) * WT
    vp <- exp(lvp) * WT
    kdec <- exp(lkdec)
    cl_mmae <- exp(lcl_mmae) * WT
    q_mmae <- exp(lq_mmae) * WT
    vc_mmae <- exp(lvc_mmae) * WT
    vp_mmae <- exp(lvp_mmae) * WT

    # ---- single-cell rate constants, 1/h -> 1/day ----
    kon_d <- kon * 24
    koff_d <- koff * 24
    kint_d <- kint * 24
    kdeg_d <- kdeg * 24
    koff_tub_d <- koff_tub * 24
    kin_d <- kin_mmae * 24
    kout_d <- kout_mmae * 24

    # ---- cell-line switches ----
    ag_ex <- exp(lag_n87) * CELLLINE_NCIN87 + exp(lag_mcf7) * (1 - CELLLINE_NCIN87)
    dt_tumor <- exp(ldt_n87) * CELLLINE_NCIN87 + exp(ldt_mcf7) * (1 - CELLLINE_NCIN87)
    vcell <- (vcell_n87 * CELLLINE_NCIN87 + vcell_mcf7 * (1 - CELLLINE_NCIN87)) * 1e-12 # pL -> L

    # SF converts molecules to nmol (Singh 2017 Table 2: 10^9 / 6.023e23).
    # One nM inside one cell is vcell / SF molecules; tubulin (nM) and its
    # association rate (1/nM/h) are put on that per-cell molecule scale, as the
    # single-cell source describes.
    sf <- 1e9 / 6.023e23
    molec_per_nm <- vcell / sf
    tub_tot_cell <- tub_tot * molec_per_nm
    kon_tub_cell <- kon_tub * 24 / molec_per_nm

    # ---- tumor growth and kill parameters ----
    kg <- 0.693 / dt_tumor # Eq (15); also the intracellular dilution rate
    kmax <- exp(lkmax + etalkmax)
    kc50 <- exp(lkc50)
    tau <- exp(ltau + etaltau)
    hill <- exp(lhill)

    # ---- tumor geometry (Eqs 2, 21) ----
    tv_mm3 <- tumor1 + tumor2 + tumor3 + tumor4
    tv <- tv_mm3 * 1e-6 # L
    rtumor <- (3 * tv * 1000 / (4 * 3.141592653589793))^(1 / 3) # cm
    # ~10^8 cells per gram of tumor (Section 2.7.3), i.e. 10^5 cells per mm^3 as
    # in Eqs (22)-(23).
    nc_tumor <- tv_mm3 * 1e5

    # ---- vascular (2 P Rcap / Rkrogh^2, um) + surface (6 D / R^2, cm) exchange, 1/day ----
    kx_adc <- 2 * perm * rcap / rkrogh^2 + 6 * dif / rtumor^2
    kx_mmae <- 2 * perm_mmae * rcap / rkrogh^2 + 6 * dif_mmae / rtumor^2

    c_adc <- central / vc # nM
    ex_adc <- (c_adc - adc_ex_tumor / eps) * kx_adc # nM/day over tumor volume
    c_mmae_tumor_ex <- mmae_ex_tumor / (tv * eps_mmae) # nM in tumor interstitium
    ex_mmae <- (mmae_central - c_mmae_tumor_ex) * kx_mmae # nM/day over tumor volume

    # HER2 binding per cell (molecules/cell/day), Eqs (8)/(10)
    bind_adc <- kon_d * adc_ex_tumor / eps * (ag_ex - adc_bound_cell) -
      koff_d * adc_bound_cell
    # Tubulin binding per cell (molecules/cell/day), Eqs (12)/(13)
    bind_tub <- kon_tub_cell * mmae_free_cell * (tub_tot_cell - mmae_tub_cell) -
      koff_tub_d * mmae_tub_cell
    # MMAE influx into one cell (molecules/cell/day), Eq (12)
    influx_cell <- kin_d * (vcell / (tv * eps_mmae)) * (mmae_ex_tumor / sf)
    # Cells leaving the last death-transit compartment (cells/day), Eqs (22)-(23)
    death_cells <- tumor4 / tau * 1e5

    # ---- plasma ADC (Eqs 3-4, amounts) ----
    d/dt(central) <- -cl / vc * central - q / vc * central + q / vp * peripheral1 -
      kdec * central - ex_adc * tv
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # ---- plasma deconjugated (naked) trastuzumab, output only ----
    # Eq (3) removes deconjugated antibody from the ADC state, yet the paper's
    # plasma 'total trastuzumab' fit (Figure 4A) declines as ADC plus the
    # antibody left behind by deconjugation, cleared with the ADC disposition
    # parameters. This pair carries that antibody; it feeds nothing else.
    d/dt(mab_central) <- kdec * central - cl / vc * mab_central -
      q / vc * mab_central + q / vp * mab_peripheral
    d/dt(mab_peripheral) <- q / vc * mab_central - q / vp * mab_peripheral

    # ---- plasma unconjugated MMAE (Eqs 5-6, concentrations) ----
    # Formation by deconjugation and by ADC catabolism, each releasing DAR
    # payloads. The tumor exchange term is scaled by tv / vc_mmae so that the
    # amount leaving plasma equals the amount entering the tumor in Eq (23);
    # the printed Eq (5) omits this factor (see the vignette's deviations).
    d/dt(mmae_central) <- -cl_mmae / vc_mmae * mmae_central -
      q_mmae / vc_mmae * mmae_central + q_mmae / vc_mmae * mmae_peripheral +
      kdec * central * dar / vc_mmae + cl * dar * c_adc / vc_mmae -
      ex_mmae * tv / vc_mmae
    d/dt(mmae_peripheral) <- q_mmae / vp_mmae * mmae_central -
      q_mmae / vp_mmae * mmae_peripheral

    # ---- average DAR (Eq 7) ----
    d/dt(dar) <- -kdec * dar

    # ---- tumor extracellular ADC (nM) and MMAE (nmol), Eqs (22)-(23) ----
    d/dt(adc_ex_tumor) <- ex_adc - bind_adc * nc_tumor * sf / tv -
      kdec * adc_ex_tumor +
      death_cells * (adc_bound_cell + adc_lyso_cell) * sf / tv
    d/dt(mmae_ex_tumor) <- ex_mmae * tv + kdec * adc_ex_tumor * dar * tv +
      (kdec * adc_bound_cell * dar + kout_d * mmae_free_cell) * nc_tumor * sf -
      influx_cell * nc_tumor * sf +
      death_cells * (mmae_free_cell + mmae_tub_cell) * sf

    # ---- single tumor cell (molecules/cell), Eqs (10)-(13) ----
    d/dt(adc_bound_cell) <- bind_adc - (kdec + kint_d) * adc_bound_cell -
      kg * adc_bound_cell
    d/dt(adc_lyso_cell) <- kint_d * adc_bound_cell - kdeg_d * adc_lyso_cell -
      kg * adc_lyso_cell
    d/dt(mmae_free_cell) <- kdeg_d * adc_lyso_cell * dar - kout_d * mmae_free_cell -
      bind_tub + influx_cell - kg * mmae_free_cell
    d/dt(mmae_tub_cell) <- bind_tub - kg * mmae_tub_cell

    # ---- tubulin occupancy (Eq 14) and kill (Eq 16) ----
    occ_tub <- mmae_tub_cell / tub_tot_cell * 100
    kill <- kmax * occ_tub^hill / (kc50^hill + occ_tub^hill)

    # ---- tumor growth inhibition, Eqs (17)-(20) ----
    d/dt(tumor1) <- (kg - kill) * tumor1
    d/dt(tumor2) <- kill * tumor1 - tumor2 / tau
    d/dt(tumor3) <- (tumor2 - tumor3) / tau
    d/dt(tumor4) <- (tumor3 - tumor4) / tau

    # ---- initial conditions ----
    dar(0) <- dar0
    tumor1(0) <- TUM_VOL

    # ---- outputs (nM unless noted) ----
    Cc <- (central + mab_central) / vc # plasma total trastuzumab (ADC + naked)
    plasma_adc <- c_adc
    plasma_total_mmae <- c_adc * dar + mmae_central
    plasma_free_mmae <- mmae_central
    tumor_total_trastuzumab <- adc_ex_tumor +
      (adc_bound_cell + adc_lyso_cell) * nc_tumor * sf / tv
    tumor_free_mmae <- (mmae_ex_tumor +
      (mmae_free_cell + mmae_tub_cell) * nc_tumor * sf) / tv
    tumor_total_mmae <- tumor_free_mmae + tumor_total_trastuzumab * dar
    tubulin_occupancy <- occ_tub # percent
    tumor_volume <- tv_mm3 # mm^3

    Cc ~ prop(propSd)
    tumor_volume ~ prop(propSd_tumor_volume)
  })
}
