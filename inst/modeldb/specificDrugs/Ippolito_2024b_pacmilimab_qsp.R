Ippolito_2024b_pacmilimab_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of immuno-oncology (QSP-IO)",
    "describing the antitumour immune response to pacmilimab, a conditionally",
    "activated (PROBODY) anti-PD-L1 antibody, in non-small cell lung cancer",
    "(NSCLC). It is the NSCLC QSP-IO platform of Wang 2023 with the PROBODY",
    "kinetics and out-of-synapse PD-L1 of Ippolito 2024 (CPT:PSP) added and 17",
    "parameters recalibrated to NSCLC, including the tumour protease cleavage",
    "rate. Four physiological compartments (central, peripheral, tumour,",
    "tumour-draining lymph node), APC endosomal and surface compartments, three",
    "immunological-synapse compartments (T cell-cancer cell, T cell-APC,",
    "macrophage-cancer cell) and the out-of-synapse PD-L1 of cancer cells and",
    "APCs. Covers logistic tumour growth, naive and activated CD8+ T cells,",
    "regulatory and helper CD4+ T cells, self- and neo-antigen release, uptake,",
    "processing and MHC presentation, TCR-mediated activation, PD-1/PD-L1/PD-L2",
    "and CTLA-4/CD28/CD80/CD86 checkpoint binding, MDSC recruitment with",
    "arginase-I and nitric-oxide suppression, M1/M2 macrophages with CD47/SIRP-alpha",
    "phagocytosis control, and IL-2 / IL-10 / IL-12 / IFN-gamma / TGF-beta / CCL2.",
    "The antibody carries six free states per compartment (both sites masked,",
    "one or both open, one or both cleaved) with reversible unmasking,",
    "irreversible protease cleavage (1% of the tumour rate outside the tumour)",
    "and linear plus saturable central clearance, and nine mono- and bivalently",
    "bound states for each of the two PD-L1 pools (cancer cell or APC, and T",
    "cell) in each T-cell synapse. Dosing the fully masked state",
    "(q_V_C_aPDL1cc) gives pacmilimab; dosing the fully cleaved state",
    "(q_V_C_aPDL1) gives the unmasked anti-PD-L1 antibody, the paper's",
    "reference arm. The platform's anti-PD-1 and anti-CTLA-4 antibody modules",
    "(same mask kinetics; only the fully unmasked form binds PD-1 or CTLA-4)",
    "are carried as in the export but are never dosed in the paper.",
    "Deterministic mechanism model: virtual patients come from",
    "Latin hypercube sampling of parameter distributions and selection against",
    "iAtlas immune-subset ratios, not from estimated IIV, so no etas and no",
    "residual error are encoded.",
    sep = " "
  )
  reference <- paste(
    "Ippolito A, Wang H, Zhang Y, Vakil V, Popel AS. Virtual clinical trials via",
    "a QSP immuno-oncology model to simulate the response to a conditionally",
    "activated PD-L1 targeting antibody in NSCLC. J Pharmacokinet Pharmacodyn.",
    "2024;51(6):747-757. doi:10.1007/s10928-024-09928-5.",
    "Model structure and every value from the paper's Supplementary Material 2",
    "(the SimBiology export of the model).",
    sep = " "
  )
  vignette <- "Ippolito_2024b_pacmilimab_qsp"
  # States keep the SimBiology export's compartment and species names
  # (q_<compartment>_<species>, an amount) so each one can be audited against
  # the Species sheet of Supplementary Material 2; none is a reusable
  # library compartment.
  paper_specific_compartments <- c(
    "q_V_C_nT0",
    "q_V_C_T0",
    "q_V_C_nT1",
    "q_V_C_T1",
    "q_V_C_aPDL1cc",
    "q_V_C_aPDL1co",
    "q_V_C_aPDL1oo",
    "q_V_C_aPDL1c",
    "q_V_C_aPDL1o",
    "q_V_C_aPDL1",
    "q_V_C_aPD1cc",
    "q_V_C_aPD1co",
    "q_V_C_aPD1oo",
    "q_V_C_aPD1c",
    "q_V_C_aPD1o",
    "q_V_C_aPD1",
    "q_V_C_aCTLA4cc",
    "q_V_C_aCTLA4co",
    "q_V_C_aCTLA4oo",
    "q_V_C_aCTLA4c",
    "q_V_C_aCTLA4o",
    "q_V_C_aCTLA4",
    "q_V_C_Th",
    "q_V_P_nT0",
    "q_V_P_T0",
    "q_V_P_nT1",
    "q_V_P_T1",
    "q_V_P_aPDL1cc",
    "q_V_P_aPDL1co",
    "q_V_P_aPDL1oo",
    "q_V_P_aPDL1c",
    "q_V_P_aPDL1o",
    "q_V_P_aPDL1",
    "q_V_P_aPD1cc",
    "q_V_P_aPD1co",
    "q_V_P_aPD1oo",
    "q_V_P_aPD1c",
    "q_V_P_aPD1o",
    "q_V_P_aPD1",
    "q_V_P_aCTLA4cc",
    "q_V_P_aCTLA4co",
    "q_V_P_aCTLA4oo",
    "q_V_P_aCTLA4c",
    "q_V_P_aCTLA4o",
    "q_V_P_aCTLA4",
    "q_V_P_Th",
    "q_V_T_C_x",
    "q_V_T_T1_exh",
    "q_V_T_Th_exh",
    "q_V_T_C1",
    "q_V_T_T0",
    "q_V_T_T1",
    "q_V_T_IFNg",
    "q_V_T_APC",
    "q_V_T_mAPC",
    "q_V_T_P0",
    "q_V_T_P1",
    "q_V_T_aPDL1cc",
    "q_V_T_aPDL1co",
    "q_V_T_aPDL1oo",
    "q_V_T_aPDL1c",
    "q_V_T_aPDL1o",
    "q_V_T_aPDL1",
    "q_V_T_aPD1cc",
    "q_V_T_aPD1co",
    "q_V_T_aPD1oo",
    "q_V_T_aPD1c",
    "q_V_T_aPD1o",
    "q_V_T_aPD1",
    "q_V_T_aCTLA4cc",
    "q_V_T_aCTLA4co",
    "q_V_T_aCTLA4oo",
    "q_V_T_aCTLA4c",
    "q_V_T_aCTLA4o",
    "q_V_T_aCTLA4",
    "q_V_T_Th",
    "q_V_T_TGFb",
    "q_V_T_MDSC",
    "q_V_T_NO",
    "q_V_T_ArgI",
    "q_V_T_CCL2",
    "q_V_T_Mac_M1",
    "q_V_T_Mac_M2",
    "q_V_T_IL12",
    "q_V_T_IL10",
    "q_V_LN_nT0",
    "q_V_LN_aT0",
    "q_V_LN_T0",
    "q_V_LN_IL2",
    "q_V_LN_nT1",
    "q_V_LN_aT1",
    "q_V_LN_T1",
    "q_V_LN_APC",
    "q_V_LN_mAPC",
    "q_V_LN_aPDL1cc",
    "q_V_LN_aPDL1co",
    "q_V_LN_aPDL1oo",
    "q_V_LN_aPDL1c",
    "q_V_LN_aPDL1o",
    "q_V_LN_aPDL1",
    "q_V_LN_aPD1cc",
    "q_V_LN_aPD1co",
    "q_V_LN_aPD1oo",
    "q_V_LN_aPD1c",
    "q_V_LN_aPD1o",
    "q_V_LN_aPD1",
    "q_V_LN_aCTLA4cc",
    "q_V_LN_aCTLA4co",
    "q_V_LN_aCTLA4oo",
    "q_V_LN_aCTLA4c",
    "q_V_LN_aCTLA4o",
    "q_V_LN_aCTLA4",
    "q_V_LN_aTh",
    "q_V_LN_Th",
    "q_V_e_P0",
    "q_V_e_p0",
    "q_V_e_P1",
    "q_V_e_p1",
    "q_A_e_M1",
    "q_A_e_M1p0",
    "q_A_e_M1p1",
    "q_A_s_M1",
    "q_A_s_M1p0",
    "q_A_s_M1p1",
    "q_syn_T_C1_PD1_PDL1",
    "q_syn_T_C1_PD1_PDL2",
    "q_syn_T_C1_PD1",
    "q_syn_T_C1_PDL1",
    "q_syn_T_C1_PDL2",
    "q_syn_T_C1_PD1_aPD1",
    "q_syn_T_C1_PD1_aPD1_PD1",
    "q_syn_T_C1_TPDL1",
    "q_syn_T_C1_PDL1_aPDL1cRo",
    "q_syn_T_C1_PDL1_aPDL1oRo",
    "q_syn_T_C1_PDL1_aPDL1Rmc",
    "q_syn_T_C1_PDL1_aPDL1mRo",
    "q_syn_T_C1_PDL1_aPDL1Rmo",
    "q_syn_T_C1_PDL1_aPDL1Rmm",
    "q_syn_T_C1_PDL1_aPDL1RoRo_PDL1",
    "q_syn_T_C1_PDL1_aPDL1RmRo_PDL1",
    "q_syn_T_C1_PDL1_aPDL1RmRm_PDL1",
    "q_syn_T_C1_TPDL1_aPDL1cRo",
    "q_syn_T_C1_TPDL1_aPDL1oRo",
    "q_syn_T_C1_TPDL1_aPDL1Rmo",
    "q_syn_T_C1_TPDL1_aPDL1Rmc",
    "q_syn_T_C1_TPDL1_aPDL1mRo",
    "q_syn_T_C1_TPDL1_aPDL1Rmm",
    "q_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1",
    "q_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1",
    "q_syn_T_C1_TPDL1_aPDL1RmRm_TPDL1",
    "q_syn_T_C1_CD28_CD80",
    "q_syn_T_C1_CD28_CD80_CD28",
    "q_syn_T_C1_CD28_CD86",
    "q_syn_T_C1_CD80_CTLA4",
    "q_syn_T_C1_CD80_CTLA4_CD80",
    "q_syn_T_C1_CTLA4_CD80_CTLA4",
    "q_syn_T_C1_CD80_CTLA4_CD80_CTLA4",
    "q_syn_T_C1_CD86_CTLA4",
    "q_syn_T_C1_CD86_CTLA4_CD86",
    "q_syn_T_C1_PDL1_CD80",
    "q_syn_T_C1_PDL1_CD80_CD28",
    "q_syn_T_C1_PDL1_CD80_CTLA4",
    "q_syn_T_C1_CD28",
    "q_syn_T_C1_CTLA4",
    "q_syn_T_C1_CD80",
    "q_syn_T_C1_CD80m",
    "q_syn_T_C1_CD86",
    "q_syn_T_C1_CTLA4_aCTLA4",
    "q_syn_T_C1_CTLA4_aCTLA4_CTLA4",
    "q_syn_T_APC_PD1_PDL1",
    "q_syn_T_APC_PD1_PDL2",
    "q_syn_T_APC_PD1",
    "q_syn_T_APC_PDL1",
    "q_syn_T_APC_PDL2",
    "q_syn_T_APC_PD1_aPD1",
    "q_syn_T_APC_PD1_aPD1_PD1",
    "q_syn_T_APC_TPDL1",
    "q_syn_T_APC_PDL1_aPDL1cRo",
    "q_syn_T_APC_PDL1_aPDL1oRo",
    "q_syn_T_APC_PDL1_aPDL1Rmc",
    "q_syn_T_APC_PDL1_aPDL1mRo",
    "q_syn_T_APC_PDL1_aPDL1Rmo",
    "q_syn_T_APC_PDL1_aPDL1Rmm",
    "q_syn_T_APC_PDL1_aPDL1RoRo_PDL1",
    "q_syn_T_APC_PDL1_aPDL1RmRo_PDL1",
    "q_syn_T_APC_PDL1_aPDL1RmRm_PDL1",
    "q_syn_T_APC_TPDL1_aPDL1cRo",
    "q_syn_T_APC_TPDL1_aPDL1oRo",
    "q_syn_T_APC_TPDL1_aPDL1Rmo",
    "q_syn_T_APC_TPDL1_aPDL1Rmc",
    "q_syn_T_APC_TPDL1_aPDL1mRo",
    "q_syn_T_APC_TPDL1_aPDL1Rmm",
    "q_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1",
    "q_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1",
    "q_syn_T_APC_TPDL1_aPDL1RmRm_TPDL1",
    "q_syn_T_APC_CD28_CD80",
    "q_syn_T_APC_CD28_CD80_CD28",
    "q_syn_T_APC_CD28_CD86",
    "q_syn_T_APC_CD80_CTLA4",
    "q_syn_T_APC_CD80_CTLA4_CD80",
    "q_syn_T_APC_CTLA4_CD80_CTLA4",
    "q_syn_T_APC_CD80_CTLA4_CD80_CTLA4",
    "q_syn_T_APC_CD86_CTLA4",
    "q_syn_T_APC_CD86_CTLA4_CD86",
    "q_syn_T_APC_PDL1_CD80",
    "q_syn_T_APC_PDL1_CD80_CD28",
    "q_syn_T_APC_PDL1_CD80_CTLA4",
    "q_syn_T_APC_CD28",
    "q_syn_T_APC_CTLA4",
    "q_syn_T_APC_CD80",
    "q_syn_T_APC_CD80m",
    "q_syn_T_APC_CD86",
    "q_syn_T_APC_CTLA4_aCTLA4",
    "q_syn_T_APC_CTLA4_aCTLA4_CTLA4",
    "q_nonsyn_C1_T_PDL1",
    "q_nonsyn_APC_LN_PDL1",
    "q_syn_M_C_CD47",
    "q_syn_M_C_SIRPa",
    "q_syn_M_C_CD47_SIRPa",
    "q_syn_M_C_PD1_PDL1",
    "q_syn_M_C_PD1_PDL2",
    "q_syn_M_C_PD1",
    "q_syn_M_C_PDL1",
    "q_syn_M_C_PDL2",
    "q_syn_M_C_PD1_aPD1",
    "q_syn_M_C_PD1_aPD1_PD1",
    "q_syn_M_C_PDL1_aPDL1",
    "q_syn_M_C_PDL1_aPDL1_PDL1",
    "q_syn_M_C_PDL1_CD80",
    "q_syn_M_C_CD80",
    "q_syn_M_C_CD80m"
  )

  units <- list(
    time = "day",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  population <- list(
    species = "human (in silico virtual cohort)",
    n_subjects = 460L,
    disease_state = "non-small cell lung cancer",
    dose_range = paste(
      "Pacmilimab 10 mg/kg every 2 weeks for 400 days, the TNBC phase I regimen",
      "applied prospectively to NSCLC; the unmasked anti-PD-L1 antibody",
      "(atezolizumab-like) at the same dose and schedule is the reference arm.",
      sep = " "
    ),
    notes = paste(
      "20,000 proposed patients drawn by Latin hypercube sampling and grown from",
      "a single cancer cell; plausible patients whose tumours reach the 3.7 cm",
      "pre-treatment diameter were then selected against the iAtlas NSCLC",
      "distributions of three bounded immune-subset fractions, M1/(M1+M2),",
      "Treg/(Treg+CD8) and CD4/(CD4+CD8),",
      "leaving 460 virtual patients. A further ten 20,000-patient cohorts",
      "restrict the T cell clone number to successive deciles of its",
      "distribution. No clinical trial of pacmilimab in NSCLC exists; the",
      "comparison anchor is the 18% response rate of anti-PD-L1 therapy.",
      sep = " "
    )
  )
  ini({
    # ---- compartment capacities (Supplementary Material 2, Compartments sheet) ----
    vol_V_C <- fixed(6); label("Capacity of compartment V_C [L]")  # ESM 2 Compartments (V_C = 6 liter)
    vol_V_P <- fixed(60); label("Capacity of compartment V_P [L]")  # ESM 2 Compartments (V_P = 60 liter)
    vol_V_LN <- fixed(0.00137444678594553); label("Capacity of compartment V_LN [L]")  # ESM 2 Compartments (V_LN = 1374.44678594553 millimeter^3)
    vol_V_e <- fixed(4e-16); label("Capacity of compartment V_e [L]")  # ESM 2 Compartments (V_e = 4e-16 litre)
    vol_A_e <- fixed(1.5e-09); label("Capacity of compartment A_e [dm^2]")  # ESM 2 Compartments (A_e = 15 micrometer^2)
    vol_A_s <- fixed(9e-08); label("Capacity of compartment A_s [dm^2]")  # ESM 2 Compartments (A_s = 900 micrometer^2)
    vol_syn_T_C1 <- fixed(1.26e-09); label("Capacity of compartment syn_T_C1 [dm^2]")  # ESM 2 Compartments (syn_T_C1 = 12.6 micrometer^2)
    vol_syn_T_APC <- fixed(1.26e-09); label("Capacity of compartment syn_T_APC [dm^2]")  # ESM 2 Compartments (syn_T_APC = 12.6 micrometer^2)
    vol_nonsyn_C1_T <- fixed(897.270277790521); label("Capacity of compartment nonsyn_C1_T [dm^2]")  # ESM 2 Compartments (nonsyn_C1_T = 8972702777905.21 micrometer^2)
    vol_nonsyn_APC_LN <- fixed(899.99999999874); label("Capacity of compartment nonsyn_APC_LN [dm^2]")  # ESM 2 Compartments (nonsyn_APC_LN = 8999999999987.4 micrometer^2)
    vol_syn_M_C <- fixed(1.26e-09); label("Capacity of compartment syn_M_C [dm^2]")  # ESM 2 Compartments (syn_M_C = 12.6 micrometer^2)

    NAVG_nmol <- fixed(6.0221409e+14); label("Avogadro constant [molecule/nmol]")  # physical constant; bridges the nanomole-denominated binding constants to the molecule-denominated surface densities

    # ---- parameters (Supplementary Material 2, Parameters sheet), converted to day / L / dm^2 / nmol ----
    k_cell_clear <- fixed(0.1); label("Rate of dead cell clearance from tumour compartment (estimated)")  # ESM 2 Parameters (k_cell_clear = 0.1 1/day)
    cell <- fixed(1); label("unit parameter for calculation")  # ESM 2 Parameters (cell = 1 cell)
    day <- fixed(1); label("unit parameter for calculation")  # ESM 2 Parameters (day = 1 day)
    vol_cell <- fixed(2.52731128244685e-12); label("Average volume of cancer cell Volume of a cancer cell calculated based on can...")  # ESM 2 Parameters (vol_cell = 2527.31128244685 micrometer^3/cell)
    vol_Tcell <- fixed(1.75015709799539e-13); label("Average volume of T cells Volume of a T cell calculated based on the average ...")  # ESM 2 Parameters (vol_Tcell = 175.015709799539 micrometer^3/cell)
    Ve_T <- fixed(0.75); label("Void fraction of the tumor (Zhang 2010, PMID: 20713014)")  # ESM 2 Parameters (Ve_T = 0.75 dimensionless)
    k_C1_growth <- fixed(0.007); label("Cancer cell growth rate (Usuda 1994, PMID: 7922975)")  # ESM 2 Parameters (k_C1_growth = 0.007 1/day)
    C_max <- fixed(1.5e+11); label("Cancer cell capacity estimated")  # ESM 2 Parameters (C_max = 1.5e+11 cell)
    k_C1_death <- fixed(1e-04); label("Cancer cell death rate from innate immune cells (Palsson 2013, PMID: 24074340)")  # ESM 2 Parameters (k_C1_death = 1e-04 1/day)
    initial_tumour_diameter <- fixed(3.7); label("Pre-treatment tumor diameter")  # ESM 2 Parameters (initial_tumour_diameter = 3.7 centimeter)
    div_T0 <- fixed(1160000); label("T Cell Diversity (Robins 2009, PMID: 19706884)")  # ESM 2 Parameters (div_T0 = 1160000 dimensionless)
    n_T0_clones <- fixed(100); label("Number of T cell clones (Jafarnejad 2019 PMID: 31236847)")  # ESM 2 Parameters (n_T0_clones = 100 dimensionless)
    q_nT0_LN_in <- fixed(0.1); label("Rate of naive T cell transport into the LN calculated based on the steady-sta...")  # ESM 2 Parameters (q_nT0_LN_in = 0.1 1/day)
    q_T0_LN_out <- fixed(24); label("Rate of activated T cell transport out of the LN (De Boer and Perelson 1995, ...")  # ESM 2 Parameters (q_T0_LN_out = 24 1/day)
    q_nT0_LN_out <- fixed(2.88); label("Rate of naive T cell transport out of the LN (Mandl 2012, PMID: 23071319)")  # ESM 2 Parameters (q_nT0_LN_out = 2.88 1/day)
    k_T0_act <- fixed(5); label("T0 activation rate (De Boer and Perelson 1995, PMID: 7475092; Robertson-Tessi...")  # ESM 2 Parameters (k_T0_act = 5 1/day)
    k_T0_pro <- fixed(1); label("T0 proliferation rate (Marchingo 2014, PMID: 25430770)")  # ESM 2 Parameters (k_T0_pro = 1 1/day)
    k_T0_death <- fixed(0.01); label("T0 death rate (De Boer and Perelson 1995, PMID: 7475092)")  # ESM 2 Parameters (k_T0_death = 0.01 1/day)
    q_T0_P_in <- fixed(3507.84); label("rate of T0 transport into the peripheral compartment calculated based on T ce...")  # ESM 2 Parameters (q_T0_P_in = 2.436 1/minute)
    q_T0_P_out <- fixed(24); label("rate of T0 transport out of the peripheral compartment (De Boer and Perelson ...")  # ESM 2 Parameters (q_T0_P_out = 24 1/day)
    q_T0_T_in <- fixed(208.8); label("rate of T0 transport into the tumor compartment calculated based on T cell tr...")  # ESM 2 Parameters (q_T0_T_in = 0.000145 1/(centimeter^3*minute))
    q_nT0_P_in <- fixed(254.016); label("rate of nT0 transport into the peripheral compartment calculated based on T c...")  # ESM 2 Parameters (q_nT0_P_in = 0.1764 1/minute)
    q_nT0_P_out <- fixed(5.1); label("rate of nT0 transport out of the peripheral compartment calculated based on s...")  # ESM 2 Parameters (q_nT0_P_out = 5.1 1/day)
    Q_nT0_thym <- fixed(7e+07); label("Thymic output of nT0 (Bains 2009, PMID: 19179300; Yeh 2017, PMID: 28912973)")  # ESM 2 Parameters (Q_nT0_thym = 7e+07 cell/day)
    k_nT0_pro <- fixed(3.2e+08); label("rate of nT0 proliferation (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (k_nT0_pro = 3.2e+08 cell/day)
    K_nT0_pro <- fixed(1e+09); label("half-maximal peripheral proliferation of of nT0 (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (K_nT0_pro = 1e+09 cell)
    k_nT0_death <- fixed(0.002); label("rate of nT0 death (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (k_nT0_death = 0.002 1/day)
    k_IL2_deg <- fixed(288); label("rate of IL2 degradation (Lotze 1985, PMID: 3871099)")  # ESM 2 Parameters (k_IL2_deg = 0.2 1/minute)
    k_IL2_cons <- fixed(0.000144); label("rate of IL2 consumption by T cells (Lotze 1985, PMID: 3871099)")  # ESM 2 Parameters (k_IL2_cons = 6e-06 nanomole/cell/hour)
    k_IL2_sec <- fixed(0.00072); label("rate of IL2 secretion from T cells (Han 2012, PMID: 22160692; Thurley 2015, P...")  # ESM 2 Parameters (k_IL2_sec = 3e-05 nanomole/cell/hour)
    IL2_50 <- fixed(0.32); label("T cell activation half-maximal IL2 concentration (Marchingo 2014, PMID: 25430...")  # ESM 2 Parameters (IL2_50 = 0.32 nanomolarity)
    IL2_50_Treg <- fixed(0.32); label("Treg activation half-maximal IL2 concentration (Wang and Smith 1987, PMID: 31...")  # ESM 2 Parameters (IL2_50_Treg = 0.32 nanomolarity)
    N0 <- fixed(2); label("numer of activated T cell generation by TCR signaling only(Marchingo 2014, PM...")  # ESM 2 Parameters (N0 = 2 dimensionless)
    N_costim <- fixed(3); label("numer of activated T cell generation by co-stimulatory signaling only (Marchi...")  # ESM 2 Parameters (N_costim = 3 dimensionless)
    N_IL2_CD8 <- fixed(11); label("maximum number of activated T cell generations due to IL2 (Marchingo 2014, PM...")  # ESM 2 Parameters (N_IL2_CD8 = 11 dimensionless)
    N_IL2_CD4 <- fixed(8); label("maximum number of activated T cell generations due to IL2 (Marchingo 2014, PM...")  # ESM 2 Parameters (N_IL2_CD4 = 8 dimensionless)
    k_Treg <- fixed(0.1); label("Rate of T cell death by Tregs (estimated)")  # ESM 2 Parameters (k_Treg = 0.1 1/day)
    div_T1 <- fixed(1110000); label("T Cell Diversity (Robins 2009, PMID: 19706884)")  # ESM 2 Parameters (div_T1 = 1110000 dimensionless)
    n_T1_clones <- fixed(92); label("Number of T cell clones (Jafarnejad 2019 PMID: 31236847, Forde 2018, PMID: 30...")  # ESM 2 Parameters (n_T1_clones = 92 dimensionless)
    q_nT1_LN_in <- fixed(0.076); label("Rate of naive T cell transport into the LN calculated based on the steady-sta...")  # ESM 2 Parameters (q_nT1_LN_in = 0.076 1/day)
    q_T1_LN_out <- fixed(24); label("Rate of activated T cell transport out of the LN (De Boer and Perelson 1995, ...")  # ESM 2 Parameters (q_T1_LN_out = 24 1/day)
    q_nT1_LN_out <- fixed(1.8); label("Rate of naive T cell transport out of the LN (Mandl 2012, PMID: 23071319)")  # ESM 2 Parameters (q_nT1_LN_out = 1.8 1/day)
    k_T1_act <- fixed(23); label("T1 activation rate (De Boer and Perelson 1995, PMID: 7475092; Robertson-Tessi...")  # ESM 2 Parameters (k_T1_act = 23 1/day)
    k_T1_pro <- fixed(1); label("T1 proliferation rate (Marchingo 2014, PMID: 25430770)")  # ESM 2 Parameters (k_T1_pro = 1 1/day)
    k_T1_death <- fixed(0.01); label("T1 death rate (De Boer and Perelson 1995, PMID: 7475092)")  # ESM 2 Parameters (k_T1_death = 0.01 1/day)
    q_T1_P_in <- fixed(3507.84); label("rate of T1 transport into the peripheral compartment calculated based on T ce...")  # ESM 2 Parameters (q_T1_P_in = 2.436 1/minute)
    q_T1_P_out <- fixed(24); label("rate of T1 transport out of the peripheral compartment (De Boer and Perelson ...")  # ESM 2 Parameters (q_T1_P_out = 24 1/day)
    q_T1_T_in <- fixed(208.8); label("rate of T1 transport into the tumor compartment calculated based on T cell tr...")  # ESM 2 Parameters (q_T1_T_in = 0.000145 1/(centimeter^3*minute))
    q_nT1_P_in <- fixed(254.016); label("rate of nT1 transport into the peripheral compartment calculated based on T c...")  # ESM 2 Parameters (q_nT1_P_in = 0.1764 1/minute)
    q_nT1_P_out <- fixed(5.1); label("rate of nT1 transport out of the peripheral compartment calculated based on s...")  # ESM 2 Parameters (q_nT1_P_out = 5.1 1/day)
    Q_nT1_thym <- fixed(3.5e+07); label("Thymic output of nT1 (Bains 2009, PMID: 19179300; Ye 2002, PMID: 11994448)")  # ESM 2 Parameters (Q_nT1_thym = 3.5e+07 cell/day)
    k_nT1_pro <- fixed(3.2e+08); label("rate of nT1 proliferation (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (k_nT1_pro = 3.2e+08 cell/day)
    K_nT1_pro <- fixed(1e+09); label("half-maximal peripheral proliferation of of nT1 (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (K_nT1_pro = 1e+09 cell)
    k_nT1_death <- fixed(0.002); label("rate of nT1 death (Braber 2012, PMID: 22365666)")  # ESM 2 Parameters (k_nT1_death = 0.002 1/day)
    k_T1 <- fixed(0.1); label("Rate of T cell exhaustion by cancer cells (estimated)")  # ESM 2 Parameters (k_T1 = 0.1 1/day)
    k_C_T1 <- fixed(0.9); label("Rate of cancer cell death by T cells (Robertson-Tessi 2012, PMID: 22051568)")  # ESM 2 Parameters (k_C_T1 = 0.9 1/day)
    k_IFNg_Tsec <- fixed(4e-13); label("rate of IFNg secretion from CD8 T cells (Autenshlyus 2021, PMID: 34059727; Bo...")  # ESM 2 Parameters (k_IFNg_Tsec = 4e-13 nanomole/cell/day)
    K_T_C <- fixed(1.2); label("Dependence of Teff killing rate on Teff/C ratio (Robertson-Tessi 2012, PMID: ...")  # ESM 2 Parameters (K_T_C = 1.2 dimensionless)
    K_T_Treg <- fixed(11); label("Dependence of Teff killing rate on Teff/Treg ratio (Robertson-Tessi 2012, PMI...")  # ESM 2 Parameters (K_T_Treg = 11 dimensionless)
    k_APC_mat <- fixed(1.5); label("Maximum rate of APC maturation (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_APC_mat = 1.5 1/day)
    k_APC_mig <- fixed(4); label("Rate of APC migration (Russo Halin 2016, PMID: 26876174)")  # ESM 2 Parameters (k_APC_mig = 4 1/day)
    k_APC_death <- fixed(0.01); label("Rate of APC death (Marino and Kirschner 2004, PMID: 15038983)")  # ESM 2 Parameters (k_APC_death = 0.01 1/day)
    k_mAPC_death <- fixed(0.02); label("Rate of mAPC death(Marino and Kirschner 2004, PMID: 15038983)")  # ESM 2 Parameters (k_mAPC_death = 0.02 1/day)
    APC0_T <- fixed(4e+08); label("Steady state density of CD103+ APCs in the tumor, calculated to fit basal APC...")  # ESM 2 Parameters (APC0_T = 4e+05 cell/milliliter)
    APC0_LN <- fixed(1.2e+09); label("Steady state density of all APC subtypes in the LN (Catron 2004, PMID: 153579...")  # ESM 2 Parameters (APC0_LN = 1200000 cell/milliliter)
    n_sites_APC <- fixed(10); label("Maxium number of T Cells an APC can interact with (De Boer and Perelson 1995,...")  # ESM 2 Parameters (n_sites_APC = 10 dimensionless)
    kin <- fixed(14.4); label("Rate of MHC internalization (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (kin = 14.4 1/day)
    kout <- fixed(28.8); label("Rate of MHC externalization (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (kout = 28.8 1/day)
    k_P0_up <- fixed(14.4); label("Rate of antigen uptake by APCs (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_P0_up = 14.4 1/day/cell)
    k_xP0_deg <- fixed(2); label("Rate of extracellular antigen degradation (Palsson 2013, PMID: 24074340)")  # ESM 2 Parameters (k_xP0_deg = 2 1/day)
    k_P0_deg <- fixed(17.28); label("Rate of endosomal antigen degradation (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_P0_deg = 17.28 1/day)
    k_p0_deg <- fixed(144); label("Rate of endosomal epitope degradation (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_p0_deg = 144 1/day)
    k_P0_on <- fixed(0.000144); label("Rate of antigen-MHC binding (Agrawal and Linderman 1996, PMID: 8944895)")  # ESM 2 Parameters (k_P0_on = 144000 1/day/molarity)
    k_P0_d1 <- fixed(100); label("Antigen-MHC kd (Stone 2015, PMID: 25618219)")  # ESM 2 Parameters (k_P0_d1 = 1e-07 molarity)
    p0_50 <- fixed(793650.793650794); label("calculated based on the number of molecules for half-maximal T cell activatio...")  # ESM 2 Parameters (p0_50 = 7.93650793650794e-05 molecule/micrometer^2)
    P0_C1 <- fixed(0.00054); label("Concentration of P0 in C1 (estimated by 2.7e6 proteins/um^3 in a cancer cell ...")  # ESM 2 Parameters (P0_C1 = 5.4e-13 mole/cell)
    A_Tcell <- fixed(1.51310411930437e-08); label("Surface area of the T cell calculated based on the average T cell diameter (C...")  # ESM 2 Parameters (A_Tcell = 151.310411930437 micrometer^2)
    A_cell <- fixed(8.97270277791781e-08); label("Surface area of the Cancer cell calculated based on the average cancer cell d...")  # ESM 2 Parameters (A_cell = 897.270277791781 micrometer^2)
    A_APC <- fixed(9e-08); label("Surface area of the APC (Agrawal and Linderman 1996, PMID: 8944895)")  # ESM 2 Parameters (A_APC = 900 micrometer^2)
    k_M1p0_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p0_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p0_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p0_TCR_off = 1 1/second)
    k_M1p0_TCR_p <- fixed(86400); label("Rate of MHC-peptide-TCR complex modification (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p0_TCR_p = 1 1/second)
    phi_M1p0_TCR <- fixed(7776); label("Rate of MHC-peptide-TCR complex with maximal modification that leads to non-s...")  # ESM 2 Parameters (phi_M1p0_TCR = 0.09 1/second)
    N_M1p0_TCR <- fixed(10); label("Number of modifcation steps for MHC-peptide-TCR complex (Lever 2014, PMID: 25...")  # ESM 2 Parameters (N_M1p0_TCR = 10 dimensionless)
    TCR_p0_tot <- fixed(1038130806703.61); label("Total number of TCR molecules per naive T cell calculated based on the number...")  # ESM 2 Parameters (TCR_p0_tot = 103.813080670361 molecule/micrometer^2)
    k_P1_up <- fixed(14.4); label("Rate of antigen uptake by APCs (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_P1_up = 14.4 1/day/cell)
    k_xP1_deg <- fixed(2); label("Rate of extracellular antigen degradation (Palsson 2013, PMID: 24074340)")  # ESM 2 Parameters (k_xP1_deg = 2 1/day)
    k_P1_deg <- fixed(17.28); label("Rate of endosomal antigen degradation (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_P1_deg = 17.28 1/day)
    k_p1_deg <- fixed(144); label("Rate of endosomal epitope degradation (Chen 2014, PMID: 25184733)")  # ESM 2 Parameters (k_p1_deg = 144 1/day)
    k_P1_on <- fixed(0.000144); label("Rate of antigen-MHC binding (Agrawal and Linderman 1996, PMID: 8944895)")  # ESM 2 Parameters (k_P1_on = 144000 1/day/molarity)
    k_P1_d1 <- fixed(10); label("Antigen-MHC kd (Stone 2015, PMID: 25618219)")  # ESM 2 Parameters (k_P1_d1 = 1e-08 molarity)
    p1_50 <- fixed(793650.793650794); label("calculated based on the number of molecules for half-maximal T cell activatio...")  # ESM 2 Parameters (p1_50 = 7.93650793650794e-05 molecule/micrometer^2)
    P1_C1 <- fixed(0.00054); label("Concentration of P1 in C1 (estimated by 2.7e6 proteins/um^3 in a cancer cell ...")  # ESM 2 Parameters (P1_C1 = 5.4e-13 mole/cell)
    k_M1p1_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p1_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p1_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p1_TCR_off = 1 1/second)
    k_M1p1_TCR_p <- fixed(86400); label("Rate of MHC-peptide-TCR complex modification (Lever 2014, PMID: 25145757)")  # ESM 2 Parameters (k_M1p1_TCR_p = 1 1/second)
    phi_M1p1_TCR <- fixed(7776); label("Rate of MHC-peptide-TCR complex with maximal modification that leads to non-s...")  # ESM 2 Parameters (phi_M1p1_TCR = 0.09 1/second)
    N_M1p1_TCR <- fixed(10); label("Number of modifcation steps for MHC-peptide-TCR complex (Lever 2014, PMID: 25...")  # ESM 2 Parameters (N_M1p1_TCR = 10 dimensionless)
    TCR_p1_tot <- fixed(1038130806703.61); label("Total number of TCR molecules per naive T cell calculated based on the number...")  # ESM 2 Parameters (TCR_p1_tot = 103.813080670361 molecule/micrometer^2)
    q_P_aPDL1 <- fixed(0.75168); label("Capillary filtration rate of aPDL1 between central and peripheral compartment...")  # ESM 2 Parameters (q_P_aPDL1 = 8.7e-06 liter/second)
    q_T_aPDL1 <- fixed(0.00736128); label("Capillary filtration rate of aPDL1 between central and tumor compartment(Wang...")  # ESM 2 Parameters (q_T_aPDL1 = 8.52e-05 milliliter/second)
    q_LN_aPDL1 <- fixed(0.0002808); label("Capillary filtration rate of aPDL1 between central and TDLN compartment(Meije...")  # ESM 2 Parameters (q_LN_aPDL1 = 3.25e-06 milliliter/second)
    q_LD_aPDL1 <- fixed(2.16); label("Rate of lymphatic drainage of aPDL1 from TDLN to central compartment(Zhu 1996...")  # ESM 2 Parameters (q_LD_aPDL1 = 0.0015 1/minute)
    k_cl_aPDL1 <- fixed(0.324); label("Clearance rate of aPDL1 from central compartment(Stroh 2017, PMID: 27981577)")  # ESM 2 Parameters (k_cl_aPDL1 = 0.324 liter/day)
    gamma_C_aPDL1 <- fixed(0.61); label("Volume fraction of interstitial space available to aPDL1 in central compartme...")  # ESM 2 Parameters (gamma_C_aPDL1 = 0.61 dimensionless)
    gamma_P_aPDL1 <- fixed(0.068); label("Volume fraction of interstitial space available to aPDL1 in peripheral compar...")  # ESM 2 Parameters (gamma_P_aPDL1 = 0.068 dimensionless)
    gamma_T_aPDL1 <- fixed(0.522); label("Volume fraction of interstitial space available to aPDL1 in tumor compartment...")  # ESM 2 Parameters (gamma_T_aPDL1 = 0.522 dimensionless)
    gamma_LN_aPDL1 <- fixed(0.2); label("Volume fraction of interstitial space available to aPDL1 in TDLN compartment(...")  # ESM 2 Parameters (gamma_LN_aPDL1 = 0.2 dimensionless)
    k_cvg <- fixed(0.00864); label("tumor cleavage rate")  # ESM 2 Parameters (k_cvg = 0.00864 1/day)
    kFac <- fixed(100); label("Reduced cleavage outside tumor")  # ESM 2 Parameters (kFac = 100 dimensionless)
    km <- fixed(1002.24); label("Mask opening rate of aPDL1 (Stroh 2019, doi: 10.1002/psp4.1244B)")  # ESM 2 Parameters (km = 1002.24 1/day)
    K_m <- fixed(0.00571428571428571); label("Mask opening rate over closing rate aPDL1 (Stroh 2019, doi: 10.1002/psp4.1244B)")  # ESM 2 Parameters (K_m = 0.00571428571428571 dimensionless)
    k_cln_aPDL1 <- fixed(5.6); label("Non-linear clearance rate of aPDL1 from central compartment(fitted)")  # ESM 2 Parameters (k_cln_aPDL1 = 5.6 nanomole/day)
    Kc_aPDL1 <- fixed(1.69); label("Half-maximal concentration of aPDL1 in central compartment for nonlinear clea...")  # ESM 2 Parameters (Kc_aPDL1 = 1.69 nanomole/liter)
    q_P_aPD1 <- fixed(0.791424); label("Capillary filtration rate of aPD1 between central and peripheral compartment(...")  # ESM 2 Parameters (q_P_aPD1 = 9.16e-06 liter/second)
    q_T_aPD1 <- fixed(0.00736128); label("Capillary filtration rate of aPD1 between central and tumor compartment(Wang ...")  # ESM 2 Parameters (q_T_aPD1 = 8.52e-05 milliliter/second)
    q_LN_aPD1 <- fixed(0.0002808); label("Capillary filtration rate of aPD1 between central and TDLN compartment(Meijer...")  # ESM 2 Parameters (q_LN_aPD1 = 3.25e-06 milliliter/second)
    q_LD_aPD1 <- fixed(2.16); label("Rate of lymphatic drainage of aPD1 from TDLN to central compartment(Zhu 1996,...")  # ESM 2 Parameters (q_LD_aPD1 = 0.0015 1/minute)
    k_cl_aPD1 <- fixed(0.315); label("Clearance rate of aPD1 from central compartment(Ahamadi 2016, PMID: 27863186)")  # ESM 2 Parameters (k_cl_aPD1 = 0.315 liter/day)
    gamma_C_aPD1 <- fixed(0.698); label("Volume fraction of interstitial space available to aPD1 in central compartmen...")  # ESM 2 Parameters (gamma_C_aPD1 = 0.698 dimensionless)
    gamma_P_aPD1 <- fixed(0.0789); label("Volume fraction of interstitial space available to aPD1 in peripheral compart...")  # ESM 2 Parameters (gamma_P_aPD1 = 0.0789 dimensionless)
    gamma_T_aPD1 <- fixed(0.522); label("Volume fraction of interstitial space available to aPD1 in tumor compartment(...")  # ESM 2 Parameters (gamma_T_aPD1 = 0.522 dimensionless)
    gamma_LN_aPD1 <- fixed(0.2); label("Volume fraction of interstitial space available to aPD1 in TDLN compartment(W...")  # ESM 2 Parameters (gamma_LN_aPD1 = 0.2 dimensionless)
    q_P_aCTLA4 <- fixed(0.4605984); label("Capillary filtration rate of aCTLA4 between central and peripheral compartmen...")  # ESM 2 Parameters (q_P_aCTLA4 = 5.331e-06 liter/second)
    q_T_aCTLA4 <- fixed(0.00736128); label("Capillary filtration rate of aCTLA4 between central and tumor compartment(Wan...")  # ESM 2 Parameters (q_T_aCTLA4 = 8.52e-05 milliliter/second)
    q_LN_aCTLA4 <- fixed(0.0002808); label("Capillary filtration rate of aCTLA4 between central and TDLN compartment(Meij...")  # ESM 2 Parameters (q_LN_aCTLA4 = 3.25e-06 milliliter/second)
    q_LD_aCTLA4 <- fixed(2.16); label("Rate of lymphatic drainage of aCTLA4 from TDLN to central compartment(Zhu 199...")  # ESM 2 Parameters (q_LD_aCTLA4 = 0.0015 1/minute)
    k_cl_aCTLA4 <- fixed(0.3276); label("Clearance rate of aCTLA4 from central compartment(Wang 2014, PMID: 24737343)")  # ESM 2 Parameters (k_cl_aCTLA4 = 0.3276 liter/day)
    gamma_C_aCTLA4 <- fixed(0.794); label("Volume fraction of interstitial space available to aCTLA4 in central compartm...")  # ESM 2 Parameters (gamma_C_aCTLA4 = 0.794 dimensionless)
    gamma_P_aCTLA4 <- fixed(0.062); label("Volume fraction of interstitial space available to aCTLA4 in peripheral compa...")  # ESM 2 Parameters (gamma_P_aCTLA4 = 0.062 dimensionless)
    gamma_T_aCTLA4 <- fixed(0.522); label("Volume fraction of interstitial space available to aCTLA4 in tumor compartmen...")  # ESM 2 Parameters (gamma_T_aCTLA4 = 0.522 dimensionless)
    gamma_LN_aCTLA4 <- fixed(0.2); label("Volume fraction of interstitial space available to aCTLA4 in TDLN compartment...")  # ESM 2 Parameters (gamma_LN_aCTLA4 = 0.2 dimensionless)
    kon_PD1_PDL1 <- fixed(518400000); label("kon of PD1-PDL1 binding kon for PD1-PDL1 (Cheng 2013, PMID: 23417675)")  # ESM 2 Parameters (kon_PD1_PDL1 = 0.06 1/(micromolarity*nanometer*second))
    k_out_PDL1 <- fixed(50000); label("Expression rate of PDL1 on tumor cells (Mimura 2018, PMID: 29034543)")  # ESM 2 Parameters (k_out_PDL1 = 50000 molecule/day)
    k_in_PDL1 <- fixed(1); label("Degradation rate of PDL1 on tumor cells (Hsu 2018, PMID: 30442814)")  # ESM 2 Parameters (k_in_PDL1 = 1 1/day)
    r_PDL1_IFNg <- fixed(6); label("Number of folds increase of PDL1 expression by IFNg (Shin 2017, PMID: 27903500)")  # ESM 2 Parameters (r_PDL1_IFNg = 6 dimensionless)
    kon_PD1_PDL2 <- fixed(691200000); label("kon of PD1-PDL2 binding (Cheng 2013, PMID: 23417675)")  # ESM 2 Parameters (kon_PD1_PDL2 = 0.08 1/(micromolarity*nanometer*second))
    kon_PD1_aPD1 <- fixed(57.888); label("kon of PD1-aPD1 binding (Wang 2014, PMID: 24872026; Brown 2020, PMID: 32134960)")  # ESM 2 Parameters (kon_PD1_aPD1 = 670000 1/(molarity*second))
    kon_PDL1_aPDL1 <- fixed(37.152); label("kon of PDL1-aPDL1 binding (Wang 2014, PMID: 24872026)")  # ESM 2 Parameters (kon_PDL1_aPDL1 = 430000 1/(molarity*second))
    kon_CD28_CD80 <- fixed(1.152e+09); label("kon of CD28-CD80 binding (van der Merwe 1997, PMID: 9053440)")  # ESM 2 Parameters (kon_CD28_CD80 = 0.133333333333333 1/(micromolarity*nanometer*second))
    kon_CD28_CD86 <- fixed(4.032e+09); label("kon of CD28-CD86 binding (Collins 2002, PMID: 12196291)")  # ESM 2 Parameters (kon_CD28_CD86 = 0.466666666666667 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD80 <- fixed(2.88e+09); label("kon of CTLA4-CD80 binding (van der Merwe 1997, PMID: 9053440)")  # ESM 2 Parameters (kon_CTLA4_CD80 = 0.333333333333333 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD86 <- fixed(5.76e+09); label("kon of CTLA4-CD86 binding (Jansson 2005, PMID: 16034096; Collins 2002, PMID: ...")  # ESM 2 Parameters (kon_CTLA4_CD86 = 0.666666666666667 1/(micromolarity*nanometer*second))
    kon_CD80_PDL1 <- fixed(921600000.000003); label("kon of CD80-PDL1 binding (Cheng 2014, PMID: 23417675)")  # ESM 2 Parameters (kon_CD80_PDL1 = 0.106666666666667 1/(micromolarity*nanometer*second))
    kon_CTLA4_aCTLA4 <- fixed(33.0912); label("kon of CTLA4-aCTLA4 binding (Wang 2014, PMID: 24872026)")  # ESM 2 Parameters (kon_CTLA4_aCTLA4 = 383000 1/(molarity*second))
    kon_CD80_CD80 <- fixed(1699200); label("kon of CD80 self-association (Ikemizu 2000, PMID: 10661405)")  # ESM 2 Parameters (kon_CD80_CD80 = 0.000196666666666667 1/(micromolarity*nanometer*second))
    koff_PD1_PDL1 <- fixed(127526.4); label("koff of PD1-PDL1 binding calculated based on the measured kd and kon (Cheng 2...")  # ESM 2 Parameters (koff_PD1_PDL1 = 1.476 1/second)
    koff_PD1_PDL2 <- fixed(47692.8); label("koff of PD1-PDL2 binding calculated based on the measured kd and kon (Cheng 2...")  # ESM 2 Parameters (koff_PD1_PDL2 = 0.552 1/second)
    koff_PD1_aPD1 <- fixed(231.552); label("koff of PD1-aPD1 binding calculated based on the measured kd and kon (Wang 20...")  # ESM 2 Parameters (koff_PD1_aPD1 = 0.00268 1/second)
    koff_PDL1_aPDL1 <- fixed(14.8608); label("koff of PDL1-aPDL1 binding calculated based on the measured kd and kon (Wang ...")  # ESM 2 Parameters (koff_PDL1_aPDL1 = 0.000172 1/second)
    koff_CD28_CD80 <- fixed(138240); label("koff of CD28-CD80 binding calculated based on the measured kd and kon (van de...")  # ESM 2 Parameters (koff_CD28_CD80 = 1.6 1/second)
    koff_CD28_CD86 <- fixed(2419200); label("koff of CD28-CD86 binding calculated based on the measured kd and kon (Collin...")  # ESM 2 Parameters (koff_CD28_CD86 = 28 1/second)
    koff_CTLA4_CD80 <- fixed(36288); label("koff of CTLA4-CD80 binding calculated based on the measured kd and kon (van d...")  # ESM 2 Parameters (koff_CTLA4_CD80 = 0.42 1/second)
    koff_CTLA4_CD86 <- fixed(449280); label("koff of CTLA4-CD86 binding calculated based on the measured kd and kon (Janss...")  # ESM 2 Parameters (koff_CTLA4_CD86 = 5.2 1/second)
    koff_CD80_PDL1 <- fixed(519782.4); label("koff of CD80-PDL1 binding calculated based on the measured kd and kon (Cheng ...")  # ESM 2 Parameters (koff_CD80_PDL1 = 6.016 1/second)
    koff_CTLA4_aCTLA4 <- fixed(602.25984); label("koff of CTLA4-aCTLA4 binding calculated based on the measured kd and kon (Wan...")  # ESM 2 Parameters (koff_CTLA4_aCTLA4 = 0.0069706 1/second)
    koff_CD80_CD80 <- fixed(866.592); label("koff of CD80 self-association calculated based on the measured kd and kon (Ik...")  # ESM 2 Parameters (koff_CD80_CD80 = 0.01003 1/second)
    Chi_PD1_aPD1 <- fixed(3333333309.53227); label("Antibody cross-arm binding efficiency that also includes the conversion of ko...")  # ESM 2 Parameters (Chi_PD1_aPD1 = 33.3333330953227 1/nanometer)
    Chi_PDL1_aPDL1 <- fixed(3333333309.53227); label("Antibody cross-arm binding efficiency that also includes the conversion of ko...")  # ESM 2 Parameters (Chi_PDL1_aPDL1 = 33.3333330953227 1/nanometer)
    Chi_CTLA4_aCTLA4 <- fixed(3333333309.53227); label("Antibody cross-arm binding efficiency that also includes the conversion of ko...")  # ESM 2 Parameters (Chi_CTLA4_aCTLA4 = 33.3333330953227 1/nanometer)
    PD1_50 <- fixed(6e+10); label("PD1/PDL1 concentration for half-maximal T cell inactivation (Jafarnejad 2019,...")  # ESM 2 Parameters (PD1_50 = 6 molecule/micrometer^2)
    n_PD1 <- fixed(2); label("Hill coefficient for PD1/PDL1 half-maximal T cell inactivation (estimated)")  # ESM 2 Parameters (n_PD1 = 2 dimensionless)
    CD28_CD8X_50 <- fixed(2e+12); label("CD28-CD80/CD28-CD86 concentration for half-maximal T cell co-estimulation est...")  # ESM 2 Parameters (CD28_CD8X_50 = 200 molecule/micrometer^2)
    n_CD28_CD8X <- fixed(2); label("Hill coefficient for CD28-CD80/CD28-CD86 half-maximal T cell co-estimulation ...")  # ESM 2 Parameters (n_CD28_CD8X = 2 dimensionless)
    T_PD1_total <- fixed(26660); label("concentration of PD1 on T cells (Cheng 2014, PMID: 23417675; Mkrtichyan 2012,...")  # ESM 2 Parameters (T_PD1_total = 26660 molecule)
    T_CD28_total <- fixed(184000); label("concentration of CD28 on T cells (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (T_CD28_total = 184000 molecule)
    T_CTLA4_syn <- fixed(8000); label("concentration of CTLA4 on T cells (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (T_CTLA4_syn = 8000 molecule)
    T_PDL1_total <- fixed(79980); label("concentration of PDL1 on T cells (Cheng 2014, PMID: 23417675; Mkrtichyan 2012...")  # ESM 2 Parameters (T_PDL1_total = 79980 molecule)
    C1_PDL1_base <- fixed(224000); label("baseline number of PDL1 molecules per C1 cell (Cheng 2014, PMID: 23417675; Mk...")  # ESM 2 Parameters (C1_PDL1_base = 224000 molecule)
    r_PDL2C1 <- fixed(0.07); label("PDL2/PDL1 molecule ratio on C1cell (Cheng 2014, PMID: 23417675)")  # ESM 2 Parameters (r_PDL2C1 = 0.07 dimensionless)
    C1_CD80_total <- fixed(80000); label("number of CD80 molecules per C1 cell (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (C1_CD80_total = 80000 molecule)
    C1_CD86_total <- fixed(860000); label("number of CD86 molecules per C1 cell (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (C1_CD86_total = 860000 molecule)
    APC_PDL1_base <- fixed(224000); label("baseline number of PDL1 molecules per APC cell (Cheng 2014, PMID: 23417675; M...")  # ESM 2 Parameters (APC_PDL1_base = 224000 molecule)
    r_PDL2APC <- fixed(0.07); label("PDL2/PDL1 molecule ratio on APCcell (Cheng 2014, PMID: 23417675)")  # ESM 2 Parameters (r_PDL2APC = 0.07 dimensionless)
    APC_CD80_total <- fixed(80000); label("number of CD80 molecules per APC cell (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (APC_CD80_total = 80000 molecule)
    APC_CD86_total <- fixed(860000); label("number of CD86 molecules per APC cell (Jansson 2005, PMID: 16034096)")  # ESM 2 Parameters (APC_CD86_total = 860000 molecule)
    Lig_fac <- fixed(1); label("Multiplicative factor for ligand increased density")  # ESM 2 Parameters (Lig_fac = 1 dimensionless)
    k_Th_act <- fixed(10); label("Rate of T helper cell activation (Robertson-Tessi 2012, PMID: 22051568)")  # ESM 2 Parameters (k_Th_act = 10 1/day)
    k_Th_Treg <- fixed(0.022); label("Th differentiation rate to Treg (Robertson-Tessi 2012, PMID: 22051568)")  # ESM 2 Parameters (k_Th_Treg = 0.022 1/day)
    k_TGFb_Tsec <- fixed(1.2e-10); label("TGF secretion rate by Treg (Liyanage 2002, PMID: 12193750)")  # ESM 2 Parameters (k_TGFb_Tsec = 1.2e-10 nanomole/cell/day)
    k_TGFb_deg <- fixed(14.3); label("TGF degradtion rate (Robertson-Tessi 2012, PMID: 22051568)")  # ESM 2 Parameters (k_TGFb_deg = 14.3 1/day)
    TGFb_50 <- fixed(0.07); label("Half-Maximal TGFb level for Th-to-Treg differentiation / chemoresistance deve...")  # ESM 2 Parameters (TGFb_50 = 0.07 nanomolarity)
    TGFb_50_Teff <- fixed(0.14); label("Half-Maximal TGFb level for CD8 T cell inhibition (Robertson-Tessi 2012, PMID...")  # ESM 2 Parameters (TGFb_50_Teff = 0.14 nanomolarity)
    Kc_rec <- fixed(20200000); label("Half-Maximal cancer cell number for T cell recruitment (Phillis 2006, PMID: 1...")  # ESM 2 Parameters (Kc_rec = 20200000 cell^2)
    TGFbase <- fixed(0.016); label("Baseline TGFb level in breast tumor (Ivanovic 2003, PMID: 12751375; Manuyakor...")  # ESM 2 Parameters (TGFbase = 0.016 nanomolarity)
    k_IFNg_Thsec <- fixed(80.2952114266667); label("IFNg secretion rate by T helper cell estimated assuming IFNg secretion is 3 t...")  # ESM 2 Parameters (k_IFNg_Thsec = 80.2952114266667 1/day)
    k_IFNg_deg <- fixed(7.68); label("IFNg degradation rate (Hofstra 1998, PMID: 9806748)")  # ESM 2 Parameters (k_IFNg_deg = 7.68 1/day)
    IFNg_50_ind <- fixed(0.00296); label("Half-Maximal IFNg level for PD-L1 induction (Shin 2017, PMID: 27903500)")  # ESM 2 Parameters (IFNg_50_ind = 2.96 picomolarity)
    k_CCL2_sec <- fixed(1.7e-12); label("rate of CCL2 secretion (Huang 2007, PMID: 17257744; Dutta 2018, PMID: 29594759)")  # ESM 2 Parameters (k_CCL2_sec = 1.7e-12 nanomole/cell/day)
    k_CCL2_deg <- fixed(1.44); label("rate of CCL2 degradation (Tanimoto 2007, PMID: 18089573)")  # ESM 2 Parameters (k_CCL2_deg = 0.06 1/hour)
    CCL2_50 <- fixed(0.23); label("Half-maximal CCL2 level of MDSC recruitment (Opalek 2007, PMID: 17888174; Ern...")  # ESM 2 Parameters (CCL2_50 = 0.23 nanomolarity)
    k_MDSC_rec <- fixed(2e+06); label("Rate of MDSC recruitment into the tumor (Parra 2021, PMID: 33875760)")  # ESM 2 Parameters (k_MDSC_rec = 2000 cell/(milliliter*day))
    k_MDSC_death <- fixed(0.015); label("Rate of MDSC death (Lai 2018, PMID: 29735668)")  # ESM 2 Parameters (k_MDSC_death = 0.015 1/day)
    k_NO_deg <- fixed(135); label("rate of NO degradation (Hakim 1996, PMID: 8953625)")  # ESM 2 Parameters (k_NO_deg = 135 1/day)
    k_ArgI_deg <- fixed(0.173); label("rate of ArgI degradation (Schimke 1964, PMID: 14257612)")  # ESM 2 Parameters (k_ArgI_deg = 0.173 1/day)
    k_NO_sec <- fixed(4.8e-07); label("rate of NO secretion from MDSCs (Serafini 2008, PMID: 18593947)")  # ESM 2 Parameters (k_NO_sec = 4.8e-07 nanomole/cell/day)
    k_ArgI_sec <- fixed(1.4e-08); label("rate of ArgI secretion (Serafini 2008, PMID: 18593947)")  # ESM 2 Parameters (k_ArgI_sec = 0.014 (mU*microliter)/cell/day)
    ArgI_50_Teff <- fixed(61.7); label("rate of ArgI-induced T cell death (Serafini 2008, PMID: 18593947)")  # ESM 2 Parameters (ArgI_50_Teff = 61.7 mU)
    NO_50_Teff <- fixed(0.75); label("rate of NO-induced T cell death (Serafini 2008, PMID: 18593947)")  # ESM 2 Parameters (NO_50_Teff = 0.75 nanomolarity)
    ArgI_50_Treg <- fixed(22.1); label("Half-maximal ArgI level of Treg expansion (estimated based on Serafini 2008, ...")  # ESM 2 Parameters (ArgI_50_Treg = 22.1 mU)
    k_Mac_rec <- fixed(2e+08); label("Recruitment rate of macrophage into tumor compartment (Welsh 2005, PMID: 1621...")  # ESM 2 Parameters (k_Mac_rec = 2e+05 cell/(milliliter*day))
    k_Mac_death <- fixed(0.02); label("Death rate of macrophages into tumor compartment (Ginhoux 2016, PMID: 26982352)")  # ESM 2 Parameters (k_Mac_death = 0.02 1/day)
    k_TGFb_Msec <- fixed(2e-11); label("Secretion rate of TGFb by macrophage in tumor compartment (Matsuo 2000, PMID:...")  # ESM 2 Parameters (k_TGFb_Msec = 2e-11 nanomole/cell/day)
    k_IL12_sec <- fixed(8.5e-12); label("Secretion rate of IL-12 by mAPC in tumor compartment (Heufler 1996, PMID: 860...")  # ESM 2 Parameters (k_IL12_sec = 8.5e-12 nanomole/cell/day)
    k_IL12_Msec <- fixed(5e-13); label("Secretion rate of IL-12 by macrophage in tumor compartment (Clough 2007, PMID...")  # ESM 2 Parameters (k_IL12_Msec = 5e-13 nanomole/cell/day)
    k_IL12_deg <- fixed(0.5544); label("Degradation rate of IL-12 in tumor compartment (Carreno 2000, PMID: 10707873)")  # ESM 2 Parameters (k_IL12_deg = 0.0231 1/hour)
    k_IL10_sec <- fixed(3e-12); label("Secretion rate of IL-10 by macrophage in tumor compartment (Autenshlyus 2021,...")  # ESM 2 Parameters (k_IL10_sec = 3e-12 nanomole/cell/day)
    k_IL10_deg <- fixed(4); label("Degradation rate of IL-10 in tumor compartment (Saxena 2015, PMID: 25481648)")  # ESM 2 Parameters (k_IL10_deg = 4 1/day)
    k_M2_pol <- fixed(0.25); label("Rate of M1 to M2 macrophage polarization (estimated assuming that polarizatio...")  # ESM 2 Parameters (k_M2_pol = 0.25 1/day)
    k_M1_pol <- fixed(0.02); label("Rate of M2 to M1 macrophage polarization (estimated based on the M1/M2 ratio)")  # ESM 2 Parameters (k_M1_pol = 0.02 1/day)
    IL10_50 <- fixed(0.008); label("Half-maximal IL10 level for M1 to M2 polarization / maintaining Treg function...")  # ESM 2 Parameters (IL10_50 = 8 picomolarity)
    IL12_50 <- fixed(0.00014); label("Half-maximal IL-12 level for M2 to M1 macrophage polarization (Altan-Bonnet 2...")  # ESM 2 Parameters (IL12_50 = 0.14 picomolarity)
    IFNg_50 <- fixed(0.0029); label("Half-maximal IFNg level for M2 to M1 macrophage polarization (Altan-Bonnet 20...")  # ESM 2 Parameters (IFNg_50 = 2.9 picomolarity)
    k_M1_phago <- fixed(0.33); label("Rate of M1 macrophage-mediated killing of cancer cell (Pan 2020, PMID: 33365025)")  # ESM 2 Parameters (k_M1_phago = 0.33 1/day)
    vol_Mcell <- fixed(4.84904826081585e-12); label("Volume of a single macrophage cell Volume of a macrophage calculated based on...")  # ESM 2 Parameters (vol_Mcell = 4849.04826081585 micrometer^3/cell)
    kon_CD47_SIRPa <- fixed(9e+06); label("kon of CD47-SIRPa binding (Hayes 2020, PMID: 31964705)")  # ESM 2 Parameters (kon_CD47_SIRPa = 0.0625 1/(micromolarity*minute*nanometer))
    koff_CD47_SIRPa <- fixed(432); label("koff of CD47-SIRPa binding (Hayes 2020, PMID: 31964705)")  # ESM 2 Parameters (koff_CD47_SIRPa = 0.3 1/minute)
    SIRPa_50 <- fixed(3.7e+11); label("SIRPa occupancy for half-maximal CD47/SIRPa inhibition on phagocytosis of can...")  # ESM 2 Parameters (SIRPa_50 = 37 molecule/micrometer^2)
    n_SIRPa <- fixed(2); label("Hill coefficient for CD47/SIRPa inhibition on phagocytosis of cancer cells by...")  # ESM 2 Parameters (n_SIRPa = 2 dimensionless)
    C_CD47 <- fixed(1e+12); label("CD47 expression on cancer cells (Morrissey 2020, PMID: 32768386; Zhang 2019, ...")  # ESM 2 Parameters (C_CD47 = 100 molecule/micrometer^2)
    M_PD1_total <- fixed(12500); label("PD-1 expression on macrophages estimated based on PD-1 expression on T cell (...")  # ESM 2 Parameters (M_PD1_total = 12500 molecule)
    M_SIRPa <- fixed(1e+12); label("SIRPa expression on macrophages (Subramanian 2006, PMID: 16291597)")  # ESM 2 Parameters (M_SIRPa = 100 molecule/micrometer^2)
    A_Mcell <- fixed(1.3854423602331e-07); label("Surface area of the macrophage calculated based on the average macrophage dia...")  # ESM 2 Parameters (A_Mcell = 1385.4423602331 micrometer^2)
    IL10_50_phago <- fixed(0.27); label("Half-maximal IL-10 level for phagocytosis by macrophage (Bian 2016, PMID: 275...")  # ESM 2 Parameters (IL10_50_phago = 270 picomolarity)
    K_Mac_C <- fixed(2); label("Association coefficient for macrophage and cancer cell (estimated)")  # ESM 2 Parameters (K_Mac_C = 2 dimensionless)
  })

  model({
    x_V_C_nT0 <- q_V_C_nT0
    x_V_C_T0 <- q_V_C_T0
    x_V_C_nT1 <- q_V_C_nT1
    x_V_C_T1 <- q_V_C_T1
    x_V_C_aPDL1cc <- q_V_C_aPDL1cc/vol_V_C
    x_V_C_aPDL1co <- q_V_C_aPDL1co/vol_V_C
    x_V_C_aPDL1oo <- q_V_C_aPDL1oo/vol_V_C
    x_V_C_aPDL1c <- q_V_C_aPDL1c/vol_V_C
    x_V_C_aPDL1o <- q_V_C_aPDL1o/vol_V_C
    x_V_C_aPDL1 <- q_V_C_aPDL1/vol_V_C
    x_V_C_aPD1cc <- q_V_C_aPD1cc/vol_V_C
    x_V_C_aPD1co <- q_V_C_aPD1co/vol_V_C
    x_V_C_aPD1oo <- q_V_C_aPD1oo/vol_V_C
    x_V_C_aPD1c <- q_V_C_aPD1c/vol_V_C
    x_V_C_aPD1o <- q_V_C_aPD1o/vol_V_C
    x_V_C_aPD1 <- q_V_C_aPD1/vol_V_C
    x_V_C_aCTLA4cc <- q_V_C_aCTLA4cc/vol_V_C
    x_V_C_aCTLA4co <- q_V_C_aCTLA4co/vol_V_C
    x_V_C_aCTLA4oo <- q_V_C_aCTLA4oo/vol_V_C
    x_V_C_aCTLA4c <- q_V_C_aCTLA4c/vol_V_C
    x_V_C_aCTLA4o <- q_V_C_aCTLA4o/vol_V_C
    x_V_C_aCTLA4 <- q_V_C_aCTLA4/vol_V_C
    x_V_C_Th <- q_V_C_Th
    x_V_P_nT0 <- q_V_P_nT0
    x_V_P_T0 <- q_V_P_T0
    x_V_P_nT1 <- q_V_P_nT1
    x_V_P_T1 <- q_V_P_T1
    x_V_P_aPDL1cc <- q_V_P_aPDL1cc/vol_V_P
    x_V_P_aPDL1co <- q_V_P_aPDL1co/vol_V_P
    x_V_P_aPDL1oo <- q_V_P_aPDL1oo/vol_V_P
    x_V_P_aPDL1c <- q_V_P_aPDL1c/vol_V_P
    x_V_P_aPDL1o <- q_V_P_aPDL1o/vol_V_P
    x_V_P_aPDL1 <- q_V_P_aPDL1/vol_V_P
    x_V_P_aPD1cc <- q_V_P_aPD1cc/vol_V_P
    x_V_P_aPD1co <- q_V_P_aPD1co/vol_V_P
    x_V_P_aPD1oo <- q_V_P_aPD1oo/vol_V_P
    x_V_P_aPD1c <- q_V_P_aPD1c/vol_V_P
    x_V_P_aPD1o <- q_V_P_aPD1o/vol_V_P
    x_V_P_aPD1 <- q_V_P_aPD1/vol_V_P
    x_V_P_aCTLA4cc <- q_V_P_aCTLA4cc/vol_V_P
    x_V_P_aCTLA4co <- q_V_P_aCTLA4co/vol_V_P
    x_V_P_aCTLA4oo <- q_V_P_aCTLA4oo/vol_V_P
    x_V_P_aCTLA4c <- q_V_P_aCTLA4c/vol_V_P
    x_V_P_aCTLA4o <- q_V_P_aCTLA4o/vol_V_P
    x_V_P_aCTLA4 <- q_V_P_aCTLA4/vol_V_P
    x_V_P_Th <- q_V_P_Th
    x_V_T_C_x <- q_V_T_C_x
    x_V_T_T1_exh <- q_V_T_T1_exh
    x_V_T_Th_exh <- q_V_T_Th_exh
    x_V_T_C1 <- q_V_T_C1
    x_V_T_T0 <- q_V_T_T0
    x_V_T_T1 <- q_V_T_T1
    x_V_T_APC <- q_V_T_APC
    x_V_T_mAPC <- q_V_T_mAPC
    x_V_T_Th <- q_V_T_Th
    x_V_T_MDSC <- q_V_T_MDSC
    x_V_T_Mac_M1 <- q_V_T_Mac_M1
    x_V_T_Mac_M2 <- q_V_T_Mac_M2
    x_V_LN_nT0 <- q_V_LN_nT0
    x_V_LN_aT0 <- q_V_LN_aT0
    x_V_LN_T0 <- q_V_LN_T0
    x_V_LN_IL2 <- q_V_LN_IL2/vol_V_LN
    x_V_LN_nT1 <- q_V_LN_nT1
    x_V_LN_aT1 <- q_V_LN_aT1
    x_V_LN_T1 <- q_V_LN_T1
    x_V_LN_APC <- q_V_LN_APC
    x_V_LN_mAPC <- q_V_LN_mAPC
    x_V_LN_aPDL1cc <- q_V_LN_aPDL1cc/vol_V_LN
    x_V_LN_aPDL1co <- q_V_LN_aPDL1co/vol_V_LN
    x_V_LN_aPDL1oo <- q_V_LN_aPDL1oo/vol_V_LN
    x_V_LN_aPDL1c <- q_V_LN_aPDL1c/vol_V_LN
    x_V_LN_aPDL1o <- q_V_LN_aPDL1o/vol_V_LN
    x_V_LN_aPDL1 <- q_V_LN_aPDL1/vol_V_LN
    x_V_LN_aPD1cc <- q_V_LN_aPD1cc/vol_V_LN
    x_V_LN_aPD1co <- q_V_LN_aPD1co/vol_V_LN
    x_V_LN_aPD1oo <- q_V_LN_aPD1oo/vol_V_LN
    x_V_LN_aPD1c <- q_V_LN_aPD1c/vol_V_LN
    x_V_LN_aPD1o <- q_V_LN_aPD1o/vol_V_LN
    x_V_LN_aPD1 <- q_V_LN_aPD1/vol_V_LN
    x_V_LN_aCTLA4cc <- q_V_LN_aCTLA4cc/vol_V_LN
    x_V_LN_aCTLA4co <- q_V_LN_aCTLA4co/vol_V_LN
    x_V_LN_aCTLA4oo <- q_V_LN_aCTLA4oo/vol_V_LN
    x_V_LN_aCTLA4c <- q_V_LN_aCTLA4c/vol_V_LN
    x_V_LN_aCTLA4o <- q_V_LN_aCTLA4o/vol_V_LN
    x_V_LN_aCTLA4 <- q_V_LN_aCTLA4/vol_V_LN
    x_V_LN_aTh <- q_V_LN_aTh
    x_V_LN_Th <- q_V_LN_Th
    x_V_e_P0 <- q_V_e_P0/vol_V_e
    x_V_e_p0 <- q_V_e_p0/vol_V_e
    x_V_e_P1 <- q_V_e_P1/vol_V_e
    x_V_e_p1 <- q_V_e_p1/vol_V_e
    x_A_e_M1 <- q_A_e_M1/vol_A_e
    x_A_e_M1p0 <- q_A_e_M1p0/vol_A_e
    x_A_e_M1p1 <- q_A_e_M1p1/vol_A_e
    x_A_s_M1 <- q_A_s_M1/vol_A_s
    x_A_s_M1p0 <- q_A_s_M1p0/vol_A_s
    x_A_s_M1p1 <- q_A_s_M1p1/vol_A_s
    x_syn_T_C1_PD1_PDL1 <- q_syn_T_C1_PD1_PDL1/vol_syn_T_C1
    x_syn_T_C1_PD1_PDL2 <- q_syn_T_C1_PD1_PDL2/vol_syn_T_C1
    x_syn_T_C1_PD1 <- q_syn_T_C1_PD1/vol_syn_T_C1
    x_syn_T_C1_PDL1 <- q_syn_T_C1_PDL1/vol_syn_T_C1
    x_syn_T_C1_PDL2 <- q_syn_T_C1_PDL2/vol_syn_T_C1
    x_syn_T_C1_PD1_aPD1 <- q_syn_T_C1_PD1_aPD1/vol_syn_T_C1
    x_syn_T_C1_PD1_aPD1_PD1 <- q_syn_T_C1_PD1_aPD1_PD1/vol_syn_T_C1
    x_syn_T_C1_TPDL1 <- q_syn_T_C1_TPDL1/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1cRo <- q_syn_T_C1_PDL1_aPDL1cRo/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1oRo <- q_syn_T_C1_PDL1_aPDL1oRo/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1Rmc <- q_syn_T_C1_PDL1_aPDL1Rmc/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1mRo <- q_syn_T_C1_PDL1_aPDL1mRo/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1Rmo <- q_syn_T_C1_PDL1_aPDL1Rmo/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1Rmm <- q_syn_T_C1_PDL1_aPDL1Rmm/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1RoRo_PDL1 <- q_syn_T_C1_PDL1_aPDL1RoRo_PDL1/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1RmRo_PDL1 <- q_syn_T_C1_PDL1_aPDL1RmRo_PDL1/vol_syn_T_C1
    x_syn_T_C1_PDL1_aPDL1RmRm_PDL1 <- q_syn_T_C1_PDL1_aPDL1RmRm_PDL1/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1cRo <- q_syn_T_C1_TPDL1_aPDL1cRo/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1oRo <- q_syn_T_C1_TPDL1_aPDL1oRo/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1Rmo <- q_syn_T_C1_TPDL1_aPDL1Rmo/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1Rmc <- q_syn_T_C1_TPDL1_aPDL1Rmc/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1mRo <- q_syn_T_C1_TPDL1_aPDL1mRo/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1Rmm <- q_syn_T_C1_TPDL1_aPDL1Rmm/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1 <- q_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1 <- q_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1/vol_syn_T_C1
    x_syn_T_C1_TPDL1_aPDL1RmRm_TPDL1 <- q_syn_T_C1_TPDL1_aPDL1RmRm_TPDL1/vol_syn_T_C1
    x_syn_T_C1_CD28_CD80 <- q_syn_T_C1_CD28_CD80/vol_syn_T_C1
    x_syn_T_C1_CD28_CD80_CD28 <- q_syn_T_C1_CD28_CD80_CD28/vol_syn_T_C1
    x_syn_T_C1_CD28_CD86 <- q_syn_T_C1_CD28_CD86/vol_syn_T_C1
    x_syn_T_C1_CD80_CTLA4 <- q_syn_T_C1_CD80_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD80_CTLA4_CD80 <- q_syn_T_C1_CD80_CTLA4_CD80/vol_syn_T_C1
    x_syn_T_C1_CTLA4_CD80_CTLA4 <- q_syn_T_C1_CTLA4_CD80_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD80_CTLA4_CD80_CTLA4 <- q_syn_T_C1_CD80_CTLA4_CD80_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD86_CTLA4 <- q_syn_T_C1_CD86_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD86_CTLA4_CD86 <- q_syn_T_C1_CD86_CTLA4_CD86/vol_syn_T_C1
    x_syn_T_C1_PDL1_CD80 <- q_syn_T_C1_PDL1_CD80/vol_syn_T_C1
    x_syn_T_C1_PDL1_CD80_CD28 <- q_syn_T_C1_PDL1_CD80_CD28/vol_syn_T_C1
    x_syn_T_C1_PDL1_CD80_CTLA4 <- q_syn_T_C1_PDL1_CD80_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD28 <- q_syn_T_C1_CD28/vol_syn_T_C1
    x_syn_T_C1_CTLA4 <- q_syn_T_C1_CTLA4/vol_syn_T_C1
    x_syn_T_C1_CD80 <- q_syn_T_C1_CD80/vol_syn_T_C1
    x_syn_T_C1_CD80m <- q_syn_T_C1_CD80m/vol_syn_T_C1
    x_syn_T_C1_CD86 <- q_syn_T_C1_CD86/vol_syn_T_C1
    x_syn_T_C1_CTLA4_aCTLA4 <- q_syn_T_C1_CTLA4_aCTLA4/vol_syn_T_C1
    x_syn_T_C1_CTLA4_aCTLA4_CTLA4 <- q_syn_T_C1_CTLA4_aCTLA4_CTLA4/vol_syn_T_C1
    x_syn_T_APC_PD1_PDL1 <- q_syn_T_APC_PD1_PDL1/vol_syn_T_APC
    x_syn_T_APC_PD1_PDL2 <- q_syn_T_APC_PD1_PDL2/vol_syn_T_APC
    x_syn_T_APC_PD1 <- q_syn_T_APC_PD1/vol_syn_T_APC
    x_syn_T_APC_PDL1 <- q_syn_T_APC_PDL1/vol_syn_T_APC
    x_syn_T_APC_PDL2 <- q_syn_T_APC_PDL2/vol_syn_T_APC
    x_syn_T_APC_PD1_aPD1 <- q_syn_T_APC_PD1_aPD1/vol_syn_T_APC
    x_syn_T_APC_PD1_aPD1_PD1 <- q_syn_T_APC_PD1_aPD1_PD1/vol_syn_T_APC
    x_syn_T_APC_TPDL1 <- q_syn_T_APC_TPDL1/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1cRo <- q_syn_T_APC_PDL1_aPDL1cRo/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1oRo <- q_syn_T_APC_PDL1_aPDL1oRo/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1Rmc <- q_syn_T_APC_PDL1_aPDL1Rmc/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1mRo <- q_syn_T_APC_PDL1_aPDL1mRo/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1Rmo <- q_syn_T_APC_PDL1_aPDL1Rmo/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1Rmm <- q_syn_T_APC_PDL1_aPDL1Rmm/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1RoRo_PDL1 <- q_syn_T_APC_PDL1_aPDL1RoRo_PDL1/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1RmRo_PDL1 <- q_syn_T_APC_PDL1_aPDL1RmRo_PDL1/vol_syn_T_APC
    x_syn_T_APC_PDL1_aPDL1RmRm_PDL1 <- q_syn_T_APC_PDL1_aPDL1RmRm_PDL1/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1cRo <- q_syn_T_APC_TPDL1_aPDL1cRo/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1oRo <- q_syn_T_APC_TPDL1_aPDL1oRo/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1Rmo <- q_syn_T_APC_TPDL1_aPDL1Rmo/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1Rmc <- q_syn_T_APC_TPDL1_aPDL1Rmc/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1mRo <- q_syn_T_APC_TPDL1_aPDL1mRo/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1Rmm <- q_syn_T_APC_TPDL1_aPDL1Rmm/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1 <- q_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1 <- q_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1/vol_syn_T_APC
    x_syn_T_APC_TPDL1_aPDL1RmRm_TPDL1 <- q_syn_T_APC_TPDL1_aPDL1RmRm_TPDL1/vol_syn_T_APC
    x_syn_T_APC_CD28_CD80 <- q_syn_T_APC_CD28_CD80/vol_syn_T_APC
    x_syn_T_APC_CD28_CD80_CD28 <- q_syn_T_APC_CD28_CD80_CD28/vol_syn_T_APC
    x_syn_T_APC_CD28_CD86 <- q_syn_T_APC_CD28_CD86/vol_syn_T_APC
    x_syn_T_APC_CD80_CTLA4 <- q_syn_T_APC_CD80_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD80_CTLA4_CD80 <- q_syn_T_APC_CD80_CTLA4_CD80/vol_syn_T_APC
    x_syn_T_APC_CTLA4_CD80_CTLA4 <- q_syn_T_APC_CTLA4_CD80_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD80_CTLA4_CD80_CTLA4 <- q_syn_T_APC_CD80_CTLA4_CD80_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD86_CTLA4 <- q_syn_T_APC_CD86_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD86_CTLA4_CD86 <- q_syn_T_APC_CD86_CTLA4_CD86/vol_syn_T_APC
    x_syn_T_APC_PDL1_CD80 <- q_syn_T_APC_PDL1_CD80/vol_syn_T_APC
    x_syn_T_APC_PDL1_CD80_CD28 <- q_syn_T_APC_PDL1_CD80_CD28/vol_syn_T_APC
    x_syn_T_APC_PDL1_CD80_CTLA4 <- q_syn_T_APC_PDL1_CD80_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD28 <- q_syn_T_APC_CD28/vol_syn_T_APC
    x_syn_T_APC_CTLA4 <- q_syn_T_APC_CTLA4/vol_syn_T_APC
    x_syn_T_APC_CD80 <- q_syn_T_APC_CD80/vol_syn_T_APC
    x_syn_T_APC_CD80m <- q_syn_T_APC_CD80m/vol_syn_T_APC
    x_syn_T_APC_CD86 <- q_syn_T_APC_CD86/vol_syn_T_APC
    x_syn_T_APC_CTLA4_aCTLA4 <- q_syn_T_APC_CTLA4_aCTLA4/vol_syn_T_APC
    x_syn_T_APC_CTLA4_aCTLA4_CTLA4 <- q_syn_T_APC_CTLA4_aCTLA4_CTLA4/vol_syn_T_APC
    x_nonsyn_C1_T_PDL1 <- q_nonsyn_C1_T_PDL1/vol_nonsyn_C1_T
    x_nonsyn_APC_LN_PDL1 <- q_nonsyn_APC_LN_PDL1/vol_nonsyn_APC_LN
    x_syn_M_C_CD47 <- q_syn_M_C_CD47/vol_syn_M_C
    x_syn_M_C_SIRPa <- q_syn_M_C_SIRPa/vol_syn_M_C
    x_syn_M_C_CD47_SIRPa <- q_syn_M_C_CD47_SIRPa/vol_syn_M_C
    x_syn_M_C_PD1_PDL1 <- q_syn_M_C_PD1_PDL1/vol_syn_M_C
    x_syn_M_C_PD1_PDL2 <- q_syn_M_C_PD1_PDL2/vol_syn_M_C
    x_syn_M_C_PD1 <- q_syn_M_C_PD1/vol_syn_M_C
    x_syn_M_C_PDL1 <- q_syn_M_C_PDL1/vol_syn_M_C
    x_syn_M_C_PDL2 <- q_syn_M_C_PDL2/vol_syn_M_C
    x_syn_M_C_PD1_aPD1 <- q_syn_M_C_PD1_aPD1/vol_syn_M_C
    x_syn_M_C_PD1_aPD1_PD1 <- q_syn_M_C_PD1_aPD1_PD1/vol_syn_M_C
    x_syn_M_C_PDL1_aPDL1 <- q_syn_M_C_PDL1_aPDL1/vol_syn_M_C
    x_syn_M_C_PDL1_aPDL1_PDL1 <- q_syn_M_C_PDL1_aPDL1_PDL1/vol_syn_M_C
    x_syn_M_C_PDL1_CD80 <- q_syn_M_C_PDL1_CD80/vol_syn_M_C
    x_syn_M_C_CD80 <- q_syn_M_C_CD80/vol_syn_M_C
    x_syn_M_C_CD80m <- q_syn_M_C_CD80m/vol_syn_M_C
    x_syn_T_C1_PDL1_aPDL1 <- 0   # no reaction forms this species in the export
    x_syn_T_C1_PDL1_aPDL1_PDL1 <- 0   # no reaction forms this species in the export
    x_syn_T_C1_TPDL1_aPDL1 <- 0   # no reaction forms this species in the export
    x_syn_T_APC_PDL1_aPDL1 <- 0   # no reaction forms this species in the export
    x_syn_T_APC_PDL1_aPDL1_PDL1 <- 0   # no reaction forms this species in the export
    x_syn_T_APC_TPDL1_aPDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1cRo <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1oRo <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1mRo <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1Rmo <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1Rmm <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1RoRo_PDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1RmRo_PDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_C1_T_PDL1_aPDL1RmRm_PDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1cRo <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1oRo <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1mRo <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1Rmo <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1Rmm <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1RoRo_PDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1RmRo_PDL1 <- 0   # no reaction forms this species in the export
    x_nonsyn_APC_LN_PDL1_aPDL1RmRm_PDL1 <- 0   # no reaction forms this species in the export
    k_C1_therapy <- 0/day   # ESM 2 Rules (k_C1_therapy)
    v6 <- Q_nT0_thym/div_T0   # R6 Thymic output of naive T cell to blood
    v32 <- Q_nT1_thym/div_T1   # R32 Thymic output of naive T cell to blood
    C_total <- 0 * cell + x_V_T_C1   # ESM 2 Rules (C_total)
    T_total <- 0 * cell + x_V_T_T0 + x_V_T_T1 + x_V_T_Th   # ESM 2 Rules (T_total)
    T_total_LN <- 0 * cell + x_V_LN_T0 + x_V_LN_T1 + x_V_LN_Th   # ESM 2 Rules (T_total_LN)
    Tregs_ <- x_V_T_T0   # ESM 2 Rules (Tregs_)
    H_APC <- n_sites_APC * x_V_LN_mAPC/(n_sites_APC * x_V_LN_mAPC + x_V_LN_nT0 * n_T0_clones + cell)   # ESM 2 Rules (H_APC)
    H_mAPC <- n_sites_APC * x_V_LN_mAPC/(n_sites_APC * x_V_LN_mAPC + x_V_LN_nT1 * n_T1_clones + cell)   # ESM 2 Rules (H_mAPC)
    H_APCh <- n_sites_APC * x_V_LN_mAPC/(n_sites_APC * x_V_LN_mAPC + x_V_LN_nT0 * n_T1_clones + cell)   # ESM 2 Rules (H_APCh)
    pTCR_p0_MHC_tot <- k_M1p0_TCR_off/(k_M1p0_TCR_off + phi_M1p0_TCR) * (k_M1p0_TCR_p/(k_M1p0_TCR_off + k_M1p0_TCR_p))^N_M1p0_TCR * 0.5 * (x_A_s_M1p0/n_T0_clones + TCR_p0_tot + k_M1p0_TCR_off/k_M1p0_TCR_on - TCR_p0_tot * sqrt(((x_A_s_M1p0/n_T0_clones + TCR_p0_tot + k_M1p0_TCR_off/k_M1p0_TCR_on)/TCR_p0_tot)^2 - 4 * x_A_s_M1p0/n_T0_clones/TCR_p0_tot))   # ESM 2 Rules (pTCR_p0_MHC_tot)
    pTCR_p1_MHC_tot <- k_M1p1_TCR_off/(k_M1p1_TCR_off + phi_M1p1_TCR) * (k_M1p1_TCR_p/(k_M1p1_TCR_off + k_M1p1_TCR_p))^N_M1p1_TCR * 0.5 * (x_A_s_M1p1/n_T1_clones + TCR_p1_tot + k_M1p1_TCR_off/k_M1p1_TCR_on - TCR_p1_tot * sqrt(((x_A_s_M1p1/n_T1_clones + TCR_p1_tot + k_M1p1_TCR_off/k_M1p1_TCR_on)/TCR_p1_tot)^2 - 4 * x_A_s_M1p1/n_T1_clones/TCR_p1_tot))   # ESM 2 Rules (pTCR_p1_MHC_tot)
    x_syn_T_C1_PDL1_total <- x_syn_T_C1_PDL1 + x_syn_T_C1_PD1_PDL1 + x_syn_T_C1_PDL1_aPDL1 + 2 * x_syn_T_C1_PDL1_aPDL1_PDL1 + x_syn_T_C1_PDL1_CD80 + x_syn_T_C1_PDL1_CD80_CD28 + x_syn_T_C1_PDL1_CD80_CTLA4   # ESM 2 Rules (syn_T_C1.PDL1_total)
    x_syn_T_C1_PDL2_total <- x_syn_T_C1_PD1_PDL2 + x_syn_T_C1_PDL2   # ESM 2 Rules (syn_T_C1.PDL2_total)
    H_PD1_C1 <- ((x_syn_T_C1_PD1_PDL1 + x_syn_T_C1_PD1_PDL2)/PD1_50)^n_PD1/(((x_syn_T_C1_PD1_PDL1 + x_syn_T_C1_PD1_PDL2)/PD1_50)^n_PD1 + 1)   # ESM 2 Rules (H_PD1_C1)
    H_CD28_C1 <- ((x_syn_T_C1_CD28_CD80 + x_syn_T_C1_CD28_CD86 + 2 * x_syn_T_C1_CD28_CD80_CD28 + x_syn_T_C1_PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X/(((x_syn_T_C1_CD28_CD80 + x_syn_T_C1_CD28_CD86 + 2 * x_syn_T_C1_CD28_CD80_CD28 + x_syn_T_C1_PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X + 1)   # ESM 2 Rules (H_CD28_C1)
    x_syn_T_APC_PDL1_total <- x_syn_T_APC_PDL1 + x_syn_T_APC_PD1_PDL1 + x_syn_T_APC_PDL1_aPDL1 + 2 * x_syn_T_APC_PDL1_aPDL1_PDL1 + x_syn_T_APC_PDL1_CD80 + x_syn_T_APC_PDL1_CD80_CD28 + x_syn_T_APC_PDL1_CD80_CTLA4   # ESM 2 Rules (syn_T_APC.PDL1_total)
    x_syn_T_APC_PDL2_total <- x_syn_T_APC_PD1_PDL2 + x_syn_T_APC_PDL2   # ESM 2 Rules (syn_T_APC.PDL2_total)
    H_PD1_APC <- ((x_syn_T_APC_PD1_PDL1 + x_syn_T_APC_PD1_PDL2)/PD1_50)^n_PD1/(((x_syn_T_APC_PD1_PDL1 + x_syn_T_APC_PD1_PDL2)/PD1_50)^n_PD1 + 1)   # ESM 2 Rules (H_PD1_APC)
    H_CD28_APC <- ((x_syn_T_APC_CD28_CD80 + x_syn_T_APC_CD28_CD86 + 2 * x_syn_T_APC_CD28_CD80_CD28 + x_syn_T_APC_PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X/(((x_syn_T_APC_CD28_CD80 + x_syn_T_APC_CD28_CD86 + 2 * x_syn_T_APC_CD28_CD80_CD28 + x_syn_T_APC_PDL1_CD80_CD28)/CD28_CD8X_50)^n_CD28_CD8X + 1)   # ESM 2 Rules (H_CD28_APC)
    x_nonsyn_C1_T_PDL1_total <- x_nonsyn_C1_T_PDL1 + x_nonsyn_C1_T_PDL1_aPDL1cRo + x_nonsyn_C1_T_PDL1_aPDL1oRo + 2 * x_nonsyn_C1_T_PDL1_aPDL1RoRo_PDL1 + x_nonsyn_C1_T_PDL1_aPDL1Rmo + x_nonsyn_C1_T_PDL1_aPDL1mRo + 2 * x_nonsyn_C1_T_PDL1_aPDL1RmRo_PDL1 + x_nonsyn_C1_T_PDL1_aPDL1Rmm + 2 * x_nonsyn_C1_T_PDL1_aPDL1RmRm_PDL1   # ESM 2 Rules (nonsyn_C1_T.PDL1_total)
    x_nonsyn_APC_LN_PDL1_total <- x_nonsyn_APC_LN_PDL1 + x_nonsyn_APC_LN_PDL1_aPDL1cRo + x_nonsyn_APC_LN_PDL1_aPDL1oRo + 2 * x_nonsyn_APC_LN_PDL1_aPDL1RoRo_PDL1 + x_nonsyn_APC_LN_PDL1_aPDL1Rmo + x_nonsyn_APC_LN_PDL1_aPDL1mRo + 2 * x_nonsyn_APC_LN_PDL1_aPDL1RmRo_PDL1 + x_nonsyn_APC_LN_PDL1_aPDL1Rmm + 2 * x_nonsyn_APC_LN_PDL1_aPDL1RmRm_PDL1   # ESM 2 Rules (nonsyn_APC_LN.PDL1_total)
    x_syn_M_C_PDL1_total <- x_syn_M_C_PDL1 + x_syn_M_C_PD1_PDL1 + x_syn_M_C_PDL1_aPDL1 + 2 * x_syn_M_C_PDL1_aPDL1_PDL1 + x_syn_M_C_PDL1_CD80   # ESM 2 Rules (syn_M_C.PDL1_total)
    x_syn_M_C_PDL2_total <- x_syn_M_C_PD1_PDL2 + x_syn_M_C_PDL2   # ESM 2 Rules (syn_M_C.PDL2_total)
    H_SIRPa <- (x_syn_M_C_CD47_SIRPa/SIRPa_50)^n_SIRPa/((x_syn_M_C_CD47_SIRPa/SIRPa_50)^n_SIRPa + 1)   # ESM 2 Rules (H_SIRPa)
    H_PD1_M <- ((x_syn_M_C_PD1_PDL1 + x_syn_M_C_PD1_PDL2)/PD1_50)^n_PD1/(((x_syn_M_C_PD1_PDL1 + x_syn_M_C_PD1_PDL2)/PD1_50)^n_PD1 + 1)   # ESM 2 Rules (H_PD1_M)
    M_total <- x_V_T_Mac_M1 + x_V_T_Mac_M2   # ESM 2 Rules (M_total)
    v1 <- k_cell_clear * x_V_T_C_x   # R1 Clearance of dead cancer cells from tumor
    v2 <- k_cell_clear * x_V_T_T1_exh   # R2 Clearance of exhausted CD8 T cells from tumor
    v3 <- k_cell_clear * x_V_T_Th_exh   # R3 Clearance of exhausted CD4 T cells from tumor
    v5 <- k_C1_death * x_V_T_C1   # R5 Cancer cell death
    v7 <- k_nT0_pro/div_T0 * x_V_P_nT0/(K_nT0_pro/div_T0 + x_V_P_nT0)   # R7 Naive T cell proliferation in the peripheral compartment
    v8 <- k_nT0_pro/div_T0 * x_V_LN_nT0/(K_nT0_pro/div_T0 + x_V_LN_nT0)   # R8 Naive T cell proliferation in the TDLN compartment
    v9 <- k_nT0_death * x_V_P_nT0   # R9 Naive T cell death in the peripheral compartment
    v10 <- k_nT0_death * x_V_C_nT0   # R10 Naive T cell death in the central compartment
    v11 <- k_nT0_death * x_V_LN_nT0   # R11 Naive T cell death in the TDLN compartment
    v12 <- q_nT0_P_in * x_V_C_nT0   # R12 Naive T cell entry into the peripheral compartment
    v13 <- q_nT0_P_out * x_V_P_nT0   # R13 Naive T cell exit from the peripheral compartment
    v14 <- q_nT0_LN_in * x_V_C_nT0   # R14 Naive T cell entry into the lymph node
    v15 <- q_nT0_LN_out * x_V_LN_nT0   # R15 Naive T cell exit from the lymph node
    v20 <- k_T0_death * x_V_C_T0   # R20 T cell death in the central compartment
    v21 <- k_T0_death * x_V_P_T0   # R21 T cell death in the peripheral compartment
    v22 <- k_T0_death * x_V_LN_T0   # R22 T cell death in the lymph node compartment
    v23 <- k_T0_death * x_V_T_T0   # R23 T cell death in the tumor compartment
    v25 <- q_T0_P_in * x_V_C_T0   # R25 T cell transport into the peripheral compartment
    v26 <- q_T0_P_out * x_V_P_T0   # R26 T cell transport out of the peripheral compartment
    v28 <- q_T0_LN_out * x_V_LN_T0   # R28 T cell transport out of the lymph node compartment
    v29 <- k_IL2_deg * x_V_LN_IL2 * vol_V_LN   # R29 IL2 degradation
    v30 <- k_IL2_cons * x_V_LN_T1 * x_V_LN_IL2/(IL2_50 + x_V_LN_IL2)   # R30 IL2 consumption by T cells
    v31 <- k_IL2_cons * x_V_LN_T0 * x_V_LN_IL2/(IL2_50_Treg + x_V_LN_IL2)   # R31 IL2 consumption by Tregs
    v33 <- k_nT1_pro/div_T1 * x_V_P_nT1/(K_nT1_pro/div_T1 + x_V_P_nT1)   # R33 Naive T cell proliferation in the peripheral compartment
    v34 <- k_nT1_pro/div_T1 * x_V_LN_nT1/(K_nT1_pro/div_T1 + x_V_LN_nT1)   # R34 Naive T cell proliferation in the TDLN compartment
    v35 <- k_nT1_death * x_V_P_nT1   # R35 Naive T cell death in the peripheral compartment
    v36 <- k_nT1_death * x_V_C_nT1   # R36 Naive T cell death in the central compartment
    v37 <- k_nT1_death * x_V_LN_nT1   # R37 Naive T cell death in the TDLN compartment
    v38 <- q_nT1_P_in * x_V_C_nT1   # R38 Naive T cell entry into the peripheral compartment
    v39 <- q_nT1_P_out * x_V_P_nT1   # R39 Naive T cell exit from the peripheral compartment
    v40 <- q_nT1_LN_in * x_V_C_nT1   # R40 Naive T cell entry into the lymph node
    v41 <- q_nT1_LN_out * x_V_LN_nT1   # R41 Naive T cell exit from the lymph node
    v46 <- k_T1_death * x_V_C_T1   # R46 T cell death in the central compartment
    v47 <- k_T1_death * x_V_P_T1   # R47 T cell death in the peripheral compartment
    v48 <- k_T1_death * x_V_LN_T1   # R48 T cell death in the lymph node compartment
    v49 <- k_T1_death * x_V_T_T1   # R49 T cell death in the tumor compartment
    v53 <- q_T1_P_in * x_V_C_T1   # R53 T cell transport into the peripheral compartment
    v54 <- q_T1_P_out * x_V_P_T1   # R54 T cell transport out of the peripheral compartment
    v56 <- q_T1_LN_out * x_V_LN_T1   # R56 T cell transport out of the lymph node compartment
    v57 <- k_IL2_sec * x_V_LN_aT1   # R57 IL2 secretion from activated T cells
    v58 <- k_IFNg_Tsec * x_V_T_T1   # R58 IFNg secretion from CD8 T cells in tumor
    v61 <- k_APC_death * (APC0_LN * vol_V_LN - x_V_LN_APC)   # R61 APC recruitment/death in LN
    v63 <- k_APC_mig * x_V_T_mAPC   # R63 APC migration to the lymph node
    v64 <- k_mAPC_death * x_V_T_mAPC   # R64 mAPC death in the tumour
    v65 <- k_mAPC_death * x_V_LN_mAPC   # R65 mAPC death in the lymph node
    v68 <- kout * x_A_e_M1 * vol_A_e - kin * x_A_s_M1 * vol_A_s   # R68 MHC translocation
    v73 <- k_P0_deg * x_V_e_P0 * vol_V_e   # R73 Antigen degradation in APC endosomes
    v74 <- k_p0_deg * x_V_e_p0 * vol_V_e   # R74 Epitope degradation in APC endosomes
    v75 <- k_P0_on * x_V_e_p0 * x_A_e_M1 * vol_A_e   # R75 Antigen-MHC binding in endosome
    v76 <- k_P0_d1 * k_P0_on * x_A_e_M1p0 * vol_A_e   # R76 Antigen-MHC unbinding in endosome
    v77 <- k_P0_d1 * k_P0_on * x_A_s_M1p0 * vol_A_s   # R77 Antigen-MHC unbinding on APC surface
    v78 <- kout * x_A_e_M1p0 * vol_A_e   # R78 Antigen-MHC translocation
    v83 <- k_P1_deg * x_V_e_P1 * vol_V_e   # R83 Antigen degradation in APC endosomes
    v84 <- k_p1_deg * x_V_e_p1 * vol_V_e   # R84 Epitope degradation in APC endosomes
    v85 <- k_P1_on * x_V_e_p1 * x_A_e_M1 * vol_A_e   # R85 Antigen-MHC binding in endosome
    v86 <- k_P1_d1 * k_P1_on * x_A_e_M1p1 * vol_A_e   # R86 Antigen-MHC unbinding in endosome
    v87 <- k_P1_d1 * k_P1_on * x_A_s_M1p1 * vol_A_s   # R87 Antigen-MHC unbinding on APC surface
    v88 <- kout * x_A_e_M1p1 * vol_A_e   # R88 Antigen-MHC translocation
    v89 <- q_P_aPDL1 * (x_V_C_aPDL1cc/gamma_C_aPDL1 - x_V_P_aPDL1cc/gamma_P_aPDL1)   # R89 aPDL1_cc diffusive transport to peripheral compartment
    v90 <- q_P_aPDL1 * (x_V_C_aPDL1co/gamma_C_aPDL1 - x_V_P_aPDL1co/gamma_P_aPDL1)   # R90 aPDL1_co diffusive transport to peripheral compartment
    v91 <- q_P_aPDL1 * (x_V_C_aPDL1c/gamma_C_aPDL1 - x_V_P_aPDL1c/gamma_P_aPDL1)   # R91 aPDL1_c diffusive transport to peripheral compartment
    v92 <- q_P_aPDL1 * (x_V_C_aPDL1o/gamma_C_aPDL1 - x_V_P_aPDL1o/gamma_P_aPDL1)   # R92 aPDL1_o diffusive transport to peripheral compartment
    v93 <- q_P_aPDL1 * (x_V_C_aPDL1/gamma_C_aPDL1 - x_V_P_aPDL1/gamma_P_aPDL1)   # R93 aPDL1 diffusive transport to peripheral compartment
    v94 <- q_P_aPDL1 * (x_V_C_aPDL1oo/gamma_C_aPDL1 - x_V_P_aPDL1oo/gamma_P_aPDL1)   # R94 aPDL1_oo diffusive transport to peripheral compartment
    v101 <- q_LN_aPDL1 * (x_V_C_aPDL1cc/gamma_C_aPDL1 - x_V_LN_aPDL1cc/gamma_LN_aPDL1)   # R101 aPDL1_cc diffusive transport to lymph node compartment
    v102 <- q_LN_aPDL1 * (x_V_C_aPDL1co/gamma_C_aPDL1 - x_V_LN_aPDL1co/gamma_LN_aPDL1)   # R102 aPDL1_co diffusive transport to lymph node compartment
    v103 <- q_LN_aPDL1 * (x_V_C_aPDL1c/gamma_C_aPDL1 - x_V_LN_aPDL1c/gamma_LN_aPDL1)   # R103 aPDL1_c diffusive transport to lymph node compartment
    v104 <- q_LN_aPDL1 * (x_V_C_aPDL1o/gamma_C_aPDL1 - x_V_LN_aPDL1o/gamma_LN_aPDL1)   # R104 aPDL1_o diffusive transport to lymph node compartment
    v105 <- q_LN_aPDL1 * (x_V_C_aPDL1/gamma_C_aPDL1 - x_V_LN_aPDL1/gamma_LN_aPDL1)   # R105 aPDL1 diffusive transport to lymph node compartment
    v106 <- q_LN_aPDL1 * (x_V_C_aPDL1oo/gamma_LN_aPDL1 - x_V_LN_aPDL1oo/gamma_LN_aPDL1)   # R106 aPDL1_oo diffusive transport to lymph node compartment
    v119 <- k_cl_aPDL1 * x_V_C_aPDL1cc   # R119 Drug clearance from central compartment
    v120 <- k_cl_aPDL1 * x_V_C_aPDL1co   # R120 Drug clearance from central compartment
    v121 <- k_cl_aPDL1 * x_V_C_aPDL1oo   # R121 Drug clearance from central compartment
    v122 <- k_cl_aPDL1 * x_V_C_aPDL1c   # R122 Drug clearance from central compartment
    v123 <- k_cl_aPDL1 * x_V_C_aPDL1o   # R123 Drug clearance from central compartment
    v124 <- k_cl_aPDL1 * x_V_C_aPDL1   # R124 Drug clearance from central compartment
    v134 <- km * (2 * x_V_C_aPDL1cc - x_V_C_aPDL1co/K_m) * vol_V_C   # R134 cc mask opening to co in C
    v135 <- 2 * k_cvg/kFac * x_V_C_aPDL1cc * vol_V_C   # R135 cc mask cleavage to c in C
    v136 <- km * (x_V_C_aPDL1co - 2 * x_V_C_aPDL1oo/K_m) * vol_V_C   # R136 co mask opening to oo in C
    v137 <- k_cvg/kFac * x_V_C_aPDL1co * vol_V_C   # R137 co mask cleavage to c in C
    v138 <- k_cvg/kFac * x_V_C_aPDL1co * vol_V_C   # R138 co mask cleavage to o in C
    v139 <- 2 * k_cvg/kFac * x_V_C_aPDL1oo * vol_V_C   # R139 oo mask cleavage to o in C
    v140 <- km * (x_V_C_aPDL1c - x_V_C_aPDL1o/K_m) * vol_V_C   # R140 c mask opening to o in C
    v141 <- k_cvg/kFac * x_V_C_aPDL1c * vol_V_C   # R141 c mask cleavage to free in C
    v142 <- k_cvg/kFac * x_V_C_aPDL1o * vol_V_C   # R142 o mask cleavage to free in C
    v143 <- km * (2 * x_V_P_aPDL1cc - x_V_P_aPDL1co/K_m) * vol_V_P   # R143 cc mask opening to co in P
    v144 <- 2 * k_cvg/kFac * x_V_P_aPDL1cc * vol_V_P   # R144 cc mask cleavage to c in P
    v145 <- km * (x_V_P_aPDL1co - 2 * x_V_P_aPDL1oo/K_m) * vol_V_P   # R145 co mask opening to oo in P
    v146 <- k_cvg/kFac * x_V_P_aPDL1co * vol_V_P   # R146 co mask cleavage to c in P
    v147 <- k_cvg/kFac * x_V_P_aPDL1co * vol_V_P   # R147 co mask cleavage to o in P
    v148 <- 2 * k_cvg/kFac * x_V_P_aPDL1oo * vol_V_P   # R148 oo mask cleavage to o in P
    v149 <- km * (x_V_P_aPDL1c - x_V_P_aPDL1o/K_m) * vol_V_P   # R149 c mask opening to o in P
    v150 <- k_cvg/kFac * x_V_P_aPDL1c * vol_V_P   # R150 c mask cleavage to free in P
    v151 <- k_cvg/kFac * x_V_P_aPDL1o * vol_V_P   # R151 o mask cleavage to free in P
    v152 <- km * (2 * x_V_LN_aPDL1cc - x_V_LN_aPDL1co/K_m) * vol_V_LN   # R152 cc mask opening to co in LN
    v153 <- 2 * k_cvg/kFac * x_V_LN_aPDL1cc * vol_V_LN   # R153 cc mask cleavage to c in LN
    v154 <- km * (x_V_LN_aPDL1co - 2 * x_V_LN_aPDL1oo/K_m) * vol_V_LN   # R154 co mask opening to oo in LN
    v155 <- k_cvg/kFac * x_V_LN_aPDL1co * vol_V_LN   # R155 co mask cleavage to c in LN
    v156 <- k_cvg/kFac * x_V_LN_aPDL1co * vol_V_LN   # R156 co mask cleavage to o in LN
    v157 <- 2 * k_cvg/kFac * x_V_LN_aPDL1oo * vol_V_LN   # R157 oo mask cleavage to o in LN
    v158 <- km * (x_V_LN_aPDL1c - x_V_LN_aPDL1o/K_m) * vol_V_LN   # R158 c mask opening to o in LN
    v159 <- k_cvg/kFac * x_V_LN_aPDL1c * vol_V_LN   # R159 c mask cleavage to free in LN
    v160 <- k_cvg/kFac * x_V_LN_aPDL1o * vol_V_LN   # R160 o mask cleavage to free in LN
    v161 <- k_cln_aPDL1 * x_V_C_aPDL1cc/(x_V_C_aPDL1cc + Kc_aPDL1)   # R161 Nonlinear clearance of aPDL1 from central compartment
    v162 <- k_cln_aPDL1 * x_V_C_aPDL1co/(x_V_C_aPDL1co + Kc_aPDL1)   # R162 Nonlinear clearance of aPDL1 from central compartment
    v163 <- k_cln_aPDL1 * x_V_C_aPDL1oo/(x_V_C_aPDL1oo + Kc_aPDL1)   # R163 Nonlinear clearance of aPDL1 from central compartment
    v164 <- k_cln_aPDL1 * x_V_C_aPDL1c/(x_V_C_aPDL1c + Kc_aPDL1)   # R164 Nonlinear clearance of aPDL1 from central compartment
    v165 <- k_cln_aPDL1 * x_V_C_aPDL1o/(x_V_C_aPDL1o + Kc_aPDL1)   # R165 Nonlinear clearance of aPDL1 from central compartment
    v166 <- k_cln_aPDL1 * x_V_C_aPDL1/(x_V_C_aPDL1 + Kc_aPDL1)   # R166 Nonlinear clearance of aPDL1 from central compartment
    v167 <- q_P_aPD1 * (x_V_C_aPD1cc/gamma_C_aPD1 - x_V_P_aPD1cc/gamma_P_aPD1)   # R167 aPD1_cc diffusive transport to peripheral compartment
    v168 <- q_P_aPD1 * (x_V_C_aPD1co/gamma_C_aPD1 - x_V_P_aPD1co/gamma_P_aPD1)   # R168 aPD1_co diffusive transport to peripheral compartment
    v169 <- q_P_aPD1 * (x_V_C_aPD1c/gamma_C_aPD1 - x_V_P_aPD1c/gamma_P_aPD1)   # R169 aPD1_c diffusive transport to peripheral compartment
    v170 <- q_P_aPD1 * (x_V_C_aPD1o/gamma_C_aPD1 - x_V_P_aPD1o/gamma_P_aPD1)   # R170 aPD1_o diffusive transport to peripheral compartment
    v171 <- q_P_aPD1 * (x_V_C_aPD1/gamma_C_aPD1 - x_V_P_aPD1/gamma_P_aPD1)   # R171 aPD1 diffusive transport to peripheral compartment
    v172 <- q_P_aPD1 * (x_V_C_aPD1oo/gamma_C_aPD1 - x_V_P_aPD1oo/gamma_P_aPD1)   # R172 aPD1_oo diffusive transport to peripheral compartment
    v179 <- q_LN_aPD1 * (x_V_C_aPD1cc/gamma_C_aPD1 - x_V_LN_aPD1cc/gamma_LN_aPD1)   # R179 aPD1_cc diffusive transport to lymph node compartment
    v180 <- q_LN_aPD1 * (x_V_C_aPD1co/gamma_C_aPD1 - x_V_LN_aPD1co/gamma_LN_aPD1)   # R180 aPD1_co diffusive transport to lymph node compartment
    v181 <- q_LN_aPD1 * (x_V_C_aPD1c/gamma_C_aPD1 - x_V_LN_aPD1c/gamma_LN_aPD1)   # R181 aPD1_c diffusive transport to lymph node compartment
    v182 <- q_LN_aPD1 * (x_V_C_aPD1o/gamma_C_aPD1 - x_V_LN_aPD1o/gamma_LN_aPD1)   # R182 aPD1_o diffusive transport to lymph node compartment
    v183 <- q_LN_aPD1 * (x_V_C_aPD1/gamma_C_aPD1 - x_V_LN_aPD1/gamma_LN_aPD1)   # R183 aPD1 diffusive transport to lymph node compartment
    v184 <- q_LN_aPD1 * (x_V_C_aPD1oo/gamma_LN_aPD1 - x_V_LN_aPD1oo/gamma_LN_aPD1)   # R184 aPD1_oo diffusive transport to lymph node compartment
    v197 <- k_cl_aPD1 * x_V_C_aPD1cc   # R197 Drug clearance from central compartment
    v198 <- k_cl_aPD1 * x_V_C_aPD1co   # R198 Drug clearance from central compartment
    v199 <- k_cl_aPD1 * x_V_C_aPD1oo   # R199 Drug clearance from central compartment
    v200 <- k_cl_aPD1 * x_V_C_aPD1c   # R200 Drug clearance from central compartment
    v201 <- k_cl_aPD1 * x_V_C_aPD1o   # R201 Drug clearance from central compartment
    v202 <- k_cl_aPD1 * x_V_C_aPD1   # R202 Drug clearance from central compartment
    v212 <- km * (2 * x_V_C_aPD1cc - x_V_C_aPD1co/K_m) * vol_V_C   # R212 cc mask opening to co in C
    v213 <- 2 * k_cvg/kFac * x_V_C_aPD1cc * vol_V_C   # R213 cc mask cleavage to c in C
    v214 <- km * (x_V_C_aPD1co - 2 * x_V_C_aPD1oo/K_m) * vol_V_C   # R214 co mask opening to oo in C
    v215 <- k_cvg/kFac * x_V_C_aPD1co * vol_V_C   # R215 co mask cleavage to c in C
    v216 <- k_cvg/kFac * x_V_C_aPD1co * vol_V_C   # R216 co mask cleavage to o in C
    v217 <- 2 * k_cvg/kFac * x_V_C_aPD1oo * vol_V_C   # R217 oo mask cleavage to o in C
    v218 <- km * (x_V_C_aPD1c - x_V_C_aPD1o/K_m) * vol_V_C   # R218 c mask opening to o in C
    v219 <- k_cvg/kFac * x_V_C_aPD1c * vol_V_C   # R219 c mask cleavage to free in C
    v220 <- k_cvg/kFac * x_V_C_aPD1o * vol_V_C   # R220 o mask cleavage to free in C
    v221 <- km * (2 * x_V_P_aPD1cc - x_V_P_aPD1co/K_m) * vol_V_P   # R221 cc mask opening to co in P
    v222 <- 2 * k_cvg/kFac * x_V_P_aPD1cc * vol_V_P   # R222 cc mask cleavage to c in P
    v223 <- km * (x_V_P_aPD1co - 2 * x_V_P_aPD1oo/K_m) * vol_V_P   # R223 co mask opening to oo in P
    v224 <- k_cvg/kFac * x_V_P_aPD1co * vol_V_P   # R224 co mask cleavage to c in P
    v225 <- k_cvg/kFac * x_V_P_aPD1co * vol_V_P   # R225 co mask cleavage to o in P
    v226 <- 2 * k_cvg/kFac * x_V_P_aPD1oo * vol_V_P   # R226 oo mask cleavage to o in P
    v227 <- km * (x_V_P_aPD1c - x_V_P_aPD1o/K_m) * vol_V_P   # R227 c mask opening to o in P
    v228 <- k_cvg/kFac * x_V_P_aPD1c * vol_V_P   # R228 c mask cleavage to free in P
    v229 <- k_cvg/kFac * x_V_P_aPD1o * vol_V_P   # R229 o mask cleavage to free in P
    v230 <- km * (2 * x_V_LN_aPD1cc - x_V_LN_aPD1co/K_m) * vol_V_LN   # R230 cc mask opening to co in LN
    v231 <- 2 * k_cvg/kFac * x_V_LN_aPD1cc * vol_V_LN   # R231 cc mask cleavage to c in LN
    v232 <- km * (x_V_LN_aPD1co - 2 * x_V_LN_aPD1oo/K_m) * vol_V_LN   # R232 co mask opening to oo in LN
    v233 <- k_cvg/kFac * x_V_LN_aPD1co * vol_V_LN   # R233 co mask cleavage to c in LN
    v234 <- k_cvg/kFac * x_V_LN_aPD1co * vol_V_LN   # R234 co mask cleavage to o in LN
    v235 <- 2 * k_cvg/kFac * x_V_LN_aPD1oo * vol_V_LN   # R235 oo mask cleavage to o in LN
    v236 <- km * (x_V_LN_aPD1c - x_V_LN_aPD1o/K_m) * vol_V_LN   # R236 c mask opening to o in LN
    v237 <- k_cvg/kFac * x_V_LN_aPD1c * vol_V_LN   # R237 c mask cleavage to free in LN
    v238 <- k_cvg/kFac * x_V_LN_aPD1o * vol_V_LN   # R238 o mask cleavage to free in LN
    v239 <- q_P_aCTLA4 * (x_V_C_aCTLA4cc/gamma_C_aCTLA4 - x_V_P_aCTLA4cc/gamma_P_aCTLA4)   # R239 aCTLA4_cc diffusive transport to peripheral compartment
    v240 <- q_P_aCTLA4 * (x_V_C_aCTLA4co/gamma_C_aCTLA4 - x_V_P_aCTLA4co/gamma_P_aCTLA4)   # R240 aCTLA4_co diffusive transport to peripheral compartment
    v241 <- q_P_aCTLA4 * (x_V_C_aCTLA4c/gamma_C_aCTLA4 - x_V_P_aCTLA4c/gamma_P_aCTLA4)   # R241 aCTLA4_c diffusive transport to peripheral compartment
    v242 <- q_P_aCTLA4 * (x_V_C_aCTLA4o/gamma_C_aCTLA4 - x_V_P_aCTLA4o/gamma_P_aCTLA4)   # R242 aCTLA4_o diffusive transport to peripheral compartment
    v243 <- q_P_aCTLA4 * (x_V_C_aCTLA4/gamma_C_aCTLA4 - x_V_P_aCTLA4/gamma_P_aCTLA4)   # R243 aCTLA4 diffusive transport to peripheral compartment
    v244 <- q_P_aCTLA4 * (x_V_C_aCTLA4oo/gamma_C_aCTLA4 - x_V_P_aCTLA4oo/gamma_P_aCTLA4)   # R244 aCTLA4_oo diffusive transport to peripheral compartment
    v251 <- q_LN_aCTLA4 * (x_V_C_aCTLA4cc/gamma_C_aCTLA4 - x_V_LN_aCTLA4cc/gamma_LN_aCTLA4)   # R251 aCTLA4_cc diffusive transport to lymph node compartment
    v252 <- q_LN_aCTLA4 * (x_V_C_aCTLA4co/gamma_C_aCTLA4 - x_V_LN_aCTLA4co/gamma_LN_aCTLA4)   # R252 aCTLA4_co diffusive transport to lymph node compartment
    v253 <- q_LN_aCTLA4 * (x_V_C_aCTLA4c/gamma_C_aCTLA4 - x_V_LN_aCTLA4c/gamma_LN_aCTLA4)   # R253 aCTLA4_c diffusive transport to lymph node compartment
    v254 <- q_LN_aCTLA4 * (x_V_C_aCTLA4o/gamma_C_aCTLA4 - x_V_LN_aCTLA4o/gamma_LN_aCTLA4)   # R254 aCTLA4_o diffusive transport to lymph node compartment
    v255 <- q_LN_aCTLA4 * (x_V_C_aCTLA4/gamma_C_aCTLA4 - x_V_LN_aCTLA4/gamma_LN_aCTLA4)   # R255 aCTLA4 diffusive transport to lymph node compartment
    v256 <- q_LN_aCTLA4 * (x_V_C_aCTLA4oo/gamma_LN_aCTLA4 - x_V_LN_aCTLA4oo/gamma_LN_aCTLA4)   # R256 aCTLA4_oo diffusive transport to lymph node compartment
    v269 <- k_cl_aCTLA4 * x_V_C_aCTLA4cc   # R269 Drug clearance from central compartment
    v270 <- k_cl_aCTLA4 * x_V_C_aCTLA4co   # R270 Drug clearance from central compartment
    v271 <- k_cl_aCTLA4 * x_V_C_aCTLA4oo   # R271 Drug clearance from central compartment
    v272 <- k_cl_aCTLA4 * x_V_C_aCTLA4c   # R272 Drug clearance from central compartment
    v273 <- k_cl_aCTLA4 * x_V_C_aCTLA4o   # R273 Drug clearance from central compartment
    v274 <- k_cl_aCTLA4 * x_V_C_aCTLA4   # R274 Drug clearance from central compartment
    v284 <- km * (2 * x_V_C_aCTLA4cc - x_V_C_aCTLA4co/K_m) * vol_V_C   # R284 cc mask opening to co in C
    v285 <- 2 * k_cvg/kFac * x_V_C_aCTLA4cc * vol_V_C   # R285 cc mask cleavage to c in C
    v286 <- km * (x_V_C_aCTLA4co - 2 * x_V_C_aCTLA4oo/K_m) * vol_V_C   # R286 co mask opening to oo in C
    v287 <- k_cvg/kFac * x_V_C_aCTLA4co * vol_V_C   # R287 co mask cleavage to c in C
    v288 <- k_cvg/kFac * x_V_C_aCTLA4co * vol_V_C   # R288 co mask cleavage to o in C
    v289 <- 2 * k_cvg/kFac * x_V_C_aCTLA4oo * vol_V_C   # R289 oo mask cleavage to o in C
    v290 <- km * (x_V_C_aCTLA4c - x_V_C_aCTLA4o/K_m) * vol_V_C   # R290 c mask opening to o in C
    v291 <- k_cvg/kFac * x_V_C_aCTLA4c * vol_V_C   # R291 c mask cleavage to free in C
    v292 <- k_cvg/kFac * x_V_C_aCTLA4o * vol_V_C   # R292 o mask cleavage to free in C
    v293 <- km * (2 * x_V_P_aCTLA4cc - x_V_P_aCTLA4co/K_m) * vol_V_P   # R293 cc mask opening to co in P
    v294 <- 2 * k_cvg/kFac * x_V_P_aCTLA4cc * vol_V_P   # R294 cc mask cleavage to c in P
    v295 <- km * (x_V_P_aCTLA4co - 2 * x_V_P_aCTLA4oo/K_m) * vol_V_P   # R295 co mask opening to oo in P
    v296 <- k_cvg/kFac * x_V_P_aCTLA4co * vol_V_P   # R296 co mask cleavage to c in P
    v297 <- k_cvg/kFac * x_V_P_aCTLA4co * vol_V_P   # R297 co mask cleavage to o in P
    v298 <- 2 * k_cvg/kFac * x_V_P_aCTLA4oo * vol_V_P   # R298 oo mask cleavage to o in P
    v299 <- km * (x_V_P_aCTLA4c - x_V_P_aCTLA4o/K_m) * vol_V_P   # R299 c mask opening to o in P
    v300 <- k_cvg/kFac * x_V_P_aCTLA4c * vol_V_P   # R300 c mask cleavage to free in P
    v301 <- k_cvg/kFac * x_V_P_aCTLA4o * vol_V_P   # R301 o mask cleavage to free in P
    v302 <- km * (2 * x_V_LN_aCTLA4cc - x_V_LN_aCTLA4co/K_m) * vol_V_LN   # R302 cc mask opening to co in LN
    v303 <- 2 * k_cvg/kFac * x_V_LN_aCTLA4cc * vol_V_LN   # R303 cc mask cleavage to c in LN
    v304 <- km * (x_V_LN_aCTLA4co - 2 * x_V_LN_aCTLA4oo/K_m) * vol_V_LN   # R304 co mask opening to oo in LN
    v305 <- k_cvg/kFac * x_V_LN_aCTLA4co * vol_V_LN   # R305 co mask cleavage to c in LN
    v306 <- k_cvg/kFac * x_V_LN_aCTLA4co * vol_V_LN   # R306 co mask cleavage to o in LN
    v307 <- 2 * k_cvg/kFac * x_V_LN_aCTLA4oo * vol_V_LN   # R307 oo mask cleavage to o in LN
    v308 <- km * (x_V_LN_aCTLA4c - x_V_LN_aCTLA4o/K_m) * vol_V_LN   # R308 c mask opening to o in LN
    v309 <- k_cvg/kFac * x_V_LN_aCTLA4c * vol_V_LN   # R309 c mask cleavage to free in LN
    v310 <- k_cvg/kFac * x_V_LN_aCTLA4o * vol_V_LN   # R310 o mask cleavage to free in LN
    v314 <- km * x_syn_T_C1_PDL1_aPDL1cRo * vol_syn_T_C1 - km/K_m * x_syn_T_C1_PDL1_aPDL1oRo * vol_syn_T_C1   # R314 breathing of probody on C1 surface in synapse
    v316 <- km * x_syn_T_C1_TPDL1_aPDL1cRo * vol_syn_T_C1 - km/K_m * x_syn_T_C1_TPDL1_aPDL1oRo * vol_syn_T_C1   # R316 breathing of probody on C1 surface in synapse
    v317 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_PDL1 * x_syn_T_C1_PDL1_aPDL1oRo * vol_syn_T_C1)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1RoRo_PDL1 * vol_syn_T_C1   # R317 binding and unbinding of PDL1 to PDL1-aPDL1oRo on C1 surface in synapse
    v318 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_TPDL1 * x_syn_T_C1_TPDL1_aPDL1oRo * vol_syn_T_C1)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1 * vol_syn_T_C1   # R318 binding and unbinding of PDL1 to PDL1-aPDL1oRo on C1 surface in synapse
    v320 <- k_cvg * x_syn_T_C1_PDL1_aPDL1cRo * vol_syn_T_C1   # R320 cleaving of probody on C1 surface in synapse
    v322 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1cRo * vol_syn_T_C1   # R322 cleaving of probody on C1 surface in synapse
    v324 <- km * x_syn_T_C1_PDL1_aPDL1Rmc * vol_syn_T_C1 - km/K_m * x_syn_T_C1_PDL1_aPDL1Rmo * vol_syn_T_C1   # R324 breathing of probody on C1 surface in synapse
    v325 <- k_cvg * x_syn_T_C1_PDL1_aPDL1oRo * vol_syn_T_C1   # R325 cleaving of probody on C1 surface in synapse
    v327 <- km * x_syn_T_C1_TPDL1_aPDL1Rmc * vol_syn_T_C1 - km/K_m * x_syn_T_C1_TPDL1_aPDL1Rmo * vol_syn_T_C1   # R327 breathing of probody on C1 surface in synapse
    v328 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1oRo * vol_syn_T_C1   # R328 cleaving of probody on C1 surface in synapse
    v330 <- k_cvg * x_syn_T_C1_PDL1_aPDL1cRo * vol_syn_T_C1   # R330 cleaving of probody on C1 surface in synapse
    v332 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1cRo * vol_syn_T_C1   # R332 cleaving of probody on C1 surface in synapse
    v333 <- 2 * k_cvg * x_syn_T_C1_PDL1_aPDL1RoRo_PDL1 * vol_syn_T_C1   # R333 cleaving of probody on C1 surface in synapse
    v334 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_PDL1 * x_syn_T_C1_PDL1_aPDL1mRo * vol_syn_T_C1)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_C1   # R334 binding and unbinding of PDL1 to PDL1-aPDL1mRo on C1 surface in synapse
    v335 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_PDL1 * x_syn_T_C1_PDL1_aPDL1Rmo * vol_syn_T_C1)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_C1   # R335 binding and unbinding of PDL1 to PDL1-aPDL1Rmo on C1 surface in synapse
    v336 <- 2 * k_cvg * x_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1 * vol_syn_T_C1   # R336 cleaving of probody on C1 surface in synapse
    v337 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_TPDL1 * x_syn_T_C1_TPDL1_aPDL1mRo * vol_syn_T_C1)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_C1   # R337 binding and unbinding of PDL1 to PDL1-aPDL1mRo on C1 surface in synapse
    v338 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_TPDL1 * x_syn_T_C1_TPDL1_aPDL1Rmo * vol_syn_T_C1)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_C1   # R338 binding and unbinding of PDL1 to PDL1-aPDL1Rmo on C1 surface in synapse
    v340 <- k_cvg * x_syn_T_C1_PDL1_aPDL1Rmo * vol_syn_T_C1   # R340 cleaving of probody on C1 surface in synapse
    v341 <- k_cvg * x_syn_T_C1_PDL1_aPDL1mRo * vol_syn_T_C1   # R341 cleaving of probody on C1 surface in synapse
    v342 <- k_cvg * x_syn_T_C1_PDL1_aPDL1Rmc * vol_syn_T_C1   # R342 cleaving of probody on C1 surface in synapse
    v344 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1Rmo * vol_syn_T_C1   # R344 cleaving of probody on C1 surface in synapse
    v345 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1mRo * vol_syn_T_C1   # R345 cleaving of probody on C1 surface in synapse
    v346 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1Rmc * vol_syn_T_C1   # R346 cleaving of probody on C1 surface in synapse
    v347 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_PDL1 * x_syn_T_C1_PDL1_aPDL1Rmm * vol_syn_T_C1)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1RmRm_PDL1 * vol_syn_T_C1   # R347 binding and unbinding of PDL1 to PDL1-aPDL1 on C1 surface in synapse
    v348 <- k_cvg * x_syn_T_C1_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_C1   # R348 cleaving of probody on C1 surface in synapse
    v349 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_C1_TPDL1 * x_syn_T_C1_TPDL1_aPDL1Rmm * vol_syn_T_C1)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1RmRm_TPDL1 * vol_syn_T_C1   # R349 binding and unbinding of PDL1 to PDL1-aPDL1 on C1 surface in synapse
    v350 <- k_cvg * x_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_C1   # R350 cleaving of probody on C1 surface in synapse
    v355 <- (kon_PD1_PDL1 * x_syn_T_C1_PD1 * x_syn_T_C1_PDL1 * vol_syn_T_C1)/NAVG_nmol - koff_PD1_PDL1 * x_syn_T_C1_PD1_PDL1 * vol_syn_T_C1   # R355 binding and unbinding of PD1 PDL1 in synapse
    v356 <- (kon_PD1_PDL2 * x_syn_T_C1_PD1 * x_syn_T_C1_PDL2 * vol_syn_T_C1)/NAVG_nmol - koff_PD1_PDL2 * x_syn_T_C1_PD1_PDL2 * vol_syn_T_C1   # R356 binding and unbinding of PD1 PDL2 in synapse
    v358 <- (Chi_PD1_aPD1 * kon_PD1_aPD1 * x_syn_T_C1_PD1 * x_syn_T_C1_PD1_aPD1 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_PD1_aPD1 * x_syn_T_C1_PD1_aPD1_PD1 * vol_syn_T_C1   # R358 binding and unbinding of PD1 to PD1-aPD1 on T surface in synapse
    v359 <- (2 * kon_CD28_CD80 * x_syn_T_C1_CD28 * x_syn_T_C1_CD80 * vol_syn_T_C1)/NAVG_nmol - koff_CD28_CD80 * x_syn_T_C1_CD28_CD80 * vol_syn_T_C1   # R359 binding and unbinding of CD28 and CD80 in synapse
    v360 <- (kon_CD28_CD80 * x_syn_T_C1_CD28 * x_syn_T_C1_CD28_CD80 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_CD28_CD80 * x_syn_T_C1_CD28_CD80_CD28 * vol_syn_T_C1   # R360 binding and unbinding of CD28-CD80 and CD28 in synapse
    v361 <- (kon_CD28_CD86 * x_syn_T_C1_CD28 * x_syn_T_C1_CD86 * vol_syn_T_C1)/NAVG_nmol - koff_CD28_CD86 * x_syn_T_C1_CD28_CD86 * vol_syn_T_C1   # R361 binding and unbinding of CD28 and CD86 in synapse
    v362 <- (4 * kon_CTLA4_CD80 * x_syn_T_C1_CTLA4 * x_syn_T_C1_CD80 * vol_syn_T_C1)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_C1_CD80_CTLA4 * vol_syn_T_C1   # R362 binding and unbinding of CTLA4 and CD80 in synapse
    v363 <- (2 * kon_CTLA4_CD80 * x_syn_T_C1_CTLA4 * x_syn_T_C1_CD80_CTLA4 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_CTLA4_CD80 * x_syn_T_C1_CTLA4_CD80_CTLA4 * vol_syn_T_C1   # R363 binding and unbinding of CTLA4 and CD80-CTLA4 in synapse
    v364 <- (4 * kon_CTLA4_CD80 * x_syn_T_C1_CD80 * x_syn_T_C1_CTLA4_CD80_CTLA4 * vol_syn_T_C1)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_C1_CD80_CTLA4_CD80_CTLA4 * vol_syn_T_C1   # R364 binding and unbinding of CD80 and CTLA4-CD80-CTLA4 in synapse
    v365 <- (2 * kon_CTLA4_CD80 * x_syn_T_C1_CD80_CTLA4 * x_syn_T_C1_CD80 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_CTLA4_CD80 * x_syn_T_C1_CD80_CTLA4_CD80 * vol_syn_T_C1   # R365 binding and unbinding of CD80-CTLA4 and CD80 in synapse
    v366 <- (4 * kon_CTLA4_CD80 * x_syn_T_C1_CTLA4 * x_syn_T_C1_CD80_CTLA4_CD80 * vol_syn_T_C1)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_C1_CD80_CTLA4_CD80_CTLA4 * vol_syn_T_C1   # R366 binding and unbinding of CTLA4 and CD80-CTLA4-CD80 in synapse
    v367 <- (2 * kon_CTLA4_CD86 * x_syn_T_C1_CTLA4 * x_syn_T_C1_CD86 * vol_syn_T_C1)/NAVG_nmol - koff_CTLA4_CD86 * x_syn_T_C1_CD86_CTLA4 * vol_syn_T_C1   # R367 binding and unbinding of CTLA4 and CD86 in synapse
    v368 <- (kon_CTLA4_CD86 * x_syn_T_C1_CD86_CTLA4 * x_syn_T_C1_CD86 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_CTLA4_CD86 * x_syn_T_C1_CD86_CTLA4_CD86 * vol_syn_T_C1   # R368 binding and unbinding of CD86-CTLA4 and CD86 in synapse
    v370 <- (2 * Chi_CTLA4_aCTLA4 * kon_CTLA4_aCTLA4 * x_syn_T_C1_CTLA4 * x_syn_T_C1_CTLA4_aCTLA4 * vol_syn_T_C1)/NAVG_nmol - 2 * koff_CTLA4_aCTLA4 * x_syn_T_C1_CTLA4_aCTLA4_CTLA4 * vol_syn_T_C1   # R370 binding and unbinding of CTLA4 to aCTLA4 on T surface in synapse
    v371 <- (kon_CD80_CD80 * x_syn_T_C1_CD80m * x_syn_T_C1_CD80m * vol_syn_T_C1)/NAVG_nmol - koff_CD80_CD80 * x_syn_T_C1_CD80 * vol_syn_T_C1   # R371 self-association and dissociation of CD80 monomers in synapse
    v372 <- (kon_CD80_PDL1 * x_syn_T_C1_CD80m * x_syn_T_C1_PDL1 * vol_syn_T_C1)/NAVG_nmol - koff_CD80_PDL1 * x_syn_T_C1_PDL1_CD80 * vol_syn_T_C1   # R372 binding and unbinding of CD80 and PDL1 in synapse
    v373 <- (kon_CD28_CD80 * x_syn_T_C1_PDL1_CD80 * x_syn_T_C1_CD28 * vol_syn_T_C1)/NAVG_nmol - koff_CD28_CD80 * x_syn_T_C1_PDL1_CD80_CD28 * vol_syn_T_C1   # R373 binding and unbinding of PDL1-CD80 and CD28 in synapse
    v374 <- (2 * kon_CTLA4_CD80 * x_syn_T_C1_PDL1_CD80 * x_syn_T_C1_CTLA4 * vol_syn_T_C1)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_C1_PDL1_CD80_CTLA4 * vol_syn_T_C1   # R374 binding and unbinding of PDL1-CD80 and CTLA4 in synapse
    v375 <- kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1co/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1cRo * vol_syn_T_APC   # R375 binding and unbinding of PDL1 to aPDL1co on APC surface in synapse
    v376 <- kon_PDL1_aPDL1 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1co/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1cRo * vol_syn_T_APC   # R376 binding and unbinding of PDL1 to aPDL1co on APC surface in synapse
    v377 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1oo/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1oRo * vol_syn_T_APC   # R377 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v378 <- km * x_syn_T_APC_PDL1_aPDL1cRo * vol_syn_T_APC - km/K_m * x_syn_T_APC_PDL1_aPDL1oRo * vol_syn_T_APC   # R378 breathing of probody on APC surface in synapse
    v379 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1oo/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1oRo * vol_syn_T_APC   # R379 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v380 <- km * x_syn_T_APC_TPDL1_aPDL1cRo * vol_syn_T_APC - km/K_m * x_syn_T_APC_TPDL1_aPDL1oRo * vol_syn_T_APC   # R380 breathing of probody on APC surface in synapse
    v381 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_PDL1 * x_syn_T_APC_PDL1_aPDL1oRo * vol_syn_T_APC)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1RoRo_PDL1 * vol_syn_T_APC   # R381 binding and unbinding of PDL1 to PDL1-aPDL1oRo on APC surface in synapse
    v382 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_TPDL1 * x_syn_T_APC_TPDL1_aPDL1oRo * vol_syn_T_APC)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1 * vol_syn_T_APC   # R382 binding and unbinding of PDL1 to PDL1-aPDL1oRo on APC surface in synapse
    v383 <- kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1c/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1Rmc * vol_syn_T_APC   # R383 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v384 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1cRo * vol_syn_T_APC   # R384 cleaving of probody on APC surface in synapse
    v385 <- kon_PDL1_aPDL1 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1c/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1Rmc * vol_syn_T_APC   # R385 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v386 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1cRo * vol_syn_T_APC   # R386 cleaving of probody on APC surface in synapse
    v387 <- kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1o/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1Rmo * vol_syn_T_APC   # R387 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v388 <- km * x_syn_T_APC_PDL1_aPDL1Rmc * vol_syn_T_APC - km/K_m * x_syn_T_APC_PDL1_aPDL1Rmo * vol_syn_T_APC   # R388 breathing of probody on APC surface in synapse
    v389 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1oRo * vol_syn_T_APC   # R389 cleaving of probody on APC surface in synapse
    v390 <- kon_PDL1_aPDL1 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1o/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1Rmo * vol_syn_T_APC   # R390 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v391 <- km * x_syn_T_APC_TPDL1_aPDL1Rmc * vol_syn_T_APC - km/K_m * x_syn_T_APC_TPDL1_aPDL1Rmo * vol_syn_T_APC   # R391 breathing of probody on APC surface in synapse
    v392 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1oRo * vol_syn_T_APC   # R392 cleaving of probody on APC surface in synapse
    v393 <- kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1o/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1mRo * vol_syn_T_APC   # R393 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v394 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1cRo * vol_syn_T_APC   # R394 cleaving of probody on APC surface in synapse
    v395 <- kon_PDL1_aPDL1 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1o/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1mRo * vol_syn_T_APC   # R395 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v396 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1cRo * vol_syn_T_APC   # R396 cleaving of probody on APC surface in synapse
    v397 <- 2 * k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1RoRo_PDL1 * vol_syn_T_APC   # R397 cleaving of probody on APC surface in synapse
    v398 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_PDL1 * x_syn_T_APC_PDL1_aPDL1mRo * vol_syn_T_APC)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_APC   # R398 binding and unbinding of PDL1 to PDL1-aPDL1mRo on APC surface in synapse
    v399 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_PDL1 * x_syn_T_APC_PDL1_aPDL1Rmo * vol_syn_T_APC)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_APC   # R399 binding and unbinding of PDL1 to PDL1-aPDL1Rmo on APC surface in synapse
    v400 <- 2 * k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1 * vol_syn_T_APC   # R400 cleaving of probody on APC surface in synapse
    v401 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_TPDL1 * x_syn_T_APC_TPDL1_aPDL1mRo * vol_syn_T_APC)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_APC   # R401 binding and unbinding of PDL1 to PDL1-aPDL1mRo on APC surface in synapse
    v402 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_TPDL1 * x_syn_T_APC_TPDL1_aPDL1Rmo * vol_syn_T_APC)/NAVG_nmol - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_APC   # R402 binding and unbinding of PDL1 to PDL1-aPDL1Rmo on APC surface in synapse
    v403 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_APC_PDL1 * x_V_LN_aPDL1/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1Rmm * vol_syn_T_APC   # R403 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v404 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1Rmo * vol_syn_T_APC   # R404 cleaving of probody on APC surface in synapse
    v405 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1mRo * vol_syn_T_APC   # R405 cleaving of probody on APC surface in synapse
    v406 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1Rmc * vol_syn_T_APC   # R406 cleaving of probody on APC surface in synapse
    v407 <- 2 * kon_PDL1_aPDL1 * 0 * (x_syn_T_APC_TPDL1 * x_V_LN_aPDL1/gamma_LN_aPDL1) * vol_syn_T_APC - koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1Rmm * vol_syn_T_APC   # R407 binding and unbinding of PDL1 to aPDL1 on APC surface in synapse
    v408 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1Rmo * vol_syn_T_APC   # R408 cleaving of probody on APC surface in synapse
    v409 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1mRo * vol_syn_T_APC   # R409 cleaving of probody on APC surface in synapse
    v410 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1Rmc * vol_syn_T_APC   # R410 cleaving of probody on APC surface in synapse
    v411 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_PDL1 * x_syn_T_APC_PDL1_aPDL1Rmm * vol_syn_T_APC)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_APC_PDL1_aPDL1RmRm_PDL1 * vol_syn_T_APC   # R411 binding and unbinding of PDL1 to PDL1-aPDL1 on APC surface in synapse
    v412 <- k_cvg/kFac * x_syn_T_APC_PDL1_aPDL1RmRo_PDL1 * vol_syn_T_APC   # R412 cleaving of probody on APC surface in synapse
    v413 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_T_APC_TPDL1 * x_syn_T_APC_TPDL1_aPDL1Rmm * vol_syn_T_APC)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_T_APC_TPDL1_aPDL1RmRm_TPDL1 * vol_syn_T_APC   # R413 binding and unbinding of PDL1 to PDL1-aPDL1 on APC surface in synapse
    v414 <- k_cvg/kFac * x_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1 * vol_syn_T_APC   # R414 cleaving of probody on APC surface in synapse
    v419 <- (kon_PD1_PDL1 * x_syn_T_APC_PD1 * x_syn_T_APC_PDL1 * vol_syn_T_APC)/NAVG_nmol - koff_PD1_PDL1 * x_syn_T_APC_PD1_PDL1 * vol_syn_T_APC   # R419 binding and unbinding of PD1 PDL1 in synapse
    v420 <- (kon_PD1_PDL2 * x_syn_T_APC_PD1 * x_syn_T_APC_PDL2 * vol_syn_T_APC)/NAVG_nmol - koff_PD1_PDL2 * x_syn_T_APC_PD1_PDL2 * vol_syn_T_APC   # R420 binding and unbinding of PD1 PDL2 in synapse
    v421 <- 2 * kon_PD1_aPD1 * (x_syn_T_APC_PD1 * x_V_LN_aPD1/gamma_LN_aPD1) * vol_syn_T_APC - koff_PD1_aPD1 * x_syn_T_APC_PD1_aPD1 * vol_syn_T_APC   # R421 binding and unbinding of PD1 to aPD1 on T surface in synapse
    v422 <- (Chi_PD1_aPD1 * kon_PD1_aPD1 * x_syn_T_APC_PD1 * x_syn_T_APC_PD1_aPD1 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_PD1_aPD1 * x_syn_T_APC_PD1_aPD1_PD1 * vol_syn_T_APC   # R422 binding and unbinding of PD1 to PD1-aPD1 on T surface in synapse
    v423 <- (2 * kon_CD28_CD80 * x_syn_T_APC_CD28 * x_syn_T_APC_CD80 * vol_syn_T_APC)/NAVG_nmol - koff_CD28_CD80 * x_syn_T_APC_CD28_CD80 * vol_syn_T_APC   # R423 binding and unbinding of CD28 and CD80 in synapse
    v424 <- (kon_CD28_CD80 * x_syn_T_APC_CD28 * x_syn_T_APC_CD28_CD80 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_CD28_CD80 * x_syn_T_APC_CD28_CD80_CD28 * vol_syn_T_APC   # R424 binding and unbinding of CD28-CD80 and CD28 in synapse
    v425 <- (kon_CD28_CD86 * x_syn_T_APC_CD28 * x_syn_T_APC_CD86 * vol_syn_T_APC)/NAVG_nmol - koff_CD28_CD86 * x_syn_T_APC_CD28_CD86 * vol_syn_T_APC   # R425 binding and unbinding of CD28 and CD86 in synapse
    v426 <- (4 * kon_CTLA4_CD80 * x_syn_T_APC_CTLA4 * x_syn_T_APC_CD80 * vol_syn_T_APC)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_APC_CD80_CTLA4 * vol_syn_T_APC   # R426 binding and unbinding of CTLA4 and CD80 in synapse
    v427 <- (2 * kon_CTLA4_CD80 * x_syn_T_APC_CTLA4 * x_syn_T_APC_CD80_CTLA4 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_CTLA4_CD80 * x_syn_T_APC_CTLA4_CD80_CTLA4 * vol_syn_T_APC   # R427 binding and unbinding of CTLA4 and CD80-CTLA4 in synapse
    v428 <- (4 * kon_CTLA4_CD80 * x_syn_T_APC_CD80 * x_syn_T_APC_CTLA4_CD80_CTLA4 * vol_syn_T_APC)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_APC_CD80_CTLA4_CD80_CTLA4 * vol_syn_T_APC   # R428 binding and unbinding of CD80 and CTLA4-CD80-CTLA4 in synapse
    v429 <- (2 * kon_CTLA4_CD80 * x_syn_T_APC_CD80_CTLA4 * x_syn_T_APC_CD80 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_CTLA4_CD80 * x_syn_T_APC_CD80_CTLA4_CD80 * vol_syn_T_APC   # R429 binding and unbinding of CD80-CTLA4 and CD80 in synapse
    v430 <- (4 * kon_CTLA4_CD80 * x_syn_T_APC_CTLA4 * x_syn_T_APC_CD80_CTLA4_CD80 * vol_syn_T_APC)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_APC_CD80_CTLA4_CD80_CTLA4 * vol_syn_T_APC   # R430 binding and unbinding of CTLA4 and CD80-CTLA4-CD80 in synapse
    v431 <- (2 * kon_CTLA4_CD86 * x_syn_T_APC_CTLA4 * x_syn_T_APC_CD86 * vol_syn_T_APC)/NAVG_nmol - koff_CTLA4_CD86 * x_syn_T_APC_CD86_CTLA4 * vol_syn_T_APC   # R431 binding and unbinding of CTLA4 and CD86 in synapse
    v432 <- (kon_CTLA4_CD86 * x_syn_T_APC_CD86_CTLA4 * x_syn_T_APC_CD86 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_CTLA4_CD86 * x_syn_T_APC_CD86_CTLA4_CD86 * vol_syn_T_APC   # R432 binding and unbinding of CD86-CTLA4 and CD86 in synapse
    v433 <- 4 * kon_CTLA4_aCTLA4 * (x_syn_T_APC_CTLA4 * x_V_LN_aCTLA4/gamma_LN_aCTLA4) * vol_syn_T_APC - koff_CTLA4_aCTLA4 * x_syn_T_APC_CTLA4_aCTLA4 * vol_syn_T_APC   # R433 binding and unbinding of CTLA4 to aCTLA4 on T surface in synapse
    v434 <- (2 * Chi_CTLA4_aCTLA4 * kon_CTLA4_aCTLA4 * x_syn_T_APC_CTLA4 * x_syn_T_APC_CTLA4_aCTLA4 * vol_syn_T_APC)/NAVG_nmol - 2 * koff_CTLA4_aCTLA4 * x_syn_T_APC_CTLA4_aCTLA4_CTLA4 * vol_syn_T_APC   # R434 binding and unbinding of CTLA4 to aCTLA4 on T surface in synapse
    v435 <- (kon_CD80_CD80 * x_syn_T_APC_CD80m * x_syn_T_APC_CD80m * vol_syn_T_APC)/NAVG_nmol - koff_CD80_CD80 * x_syn_T_APC_CD80 * vol_syn_T_APC   # R435 self-association and dissociation of CD80 monomers in synapse
    v436 <- (kon_CD80_PDL1 * x_syn_T_APC_CD80m * x_syn_T_APC_PDL1 * vol_syn_T_APC)/NAVG_nmol - koff_CD80_PDL1 * x_syn_T_APC_PDL1_CD80 * vol_syn_T_APC   # R436 binding and unbinding of CD80 and PDL1 in synapse
    v437 <- (kon_CD28_CD80 * x_syn_T_APC_PDL1_CD80 * x_syn_T_APC_CD28 * vol_syn_T_APC)/NAVG_nmol - koff_CD28_CD80 * x_syn_T_APC_PDL1_CD80_CD28 * vol_syn_T_APC   # R437 binding and unbinding of PDL1-CD80 and CD28 in synapse
    v438 <- (2 * kon_CTLA4_CD80 * x_syn_T_APC_PDL1_CD80 * x_syn_T_APC_CTLA4 * vol_syn_T_APC)/NAVG_nmol - koff_CTLA4_CD80 * x_syn_T_APC_PDL1_CD80_CTLA4 * vol_syn_T_APC   # R438 binding and unbinding of PDL1-CD80 and CTLA4 in synapse
    v448 <- k_T0_death * x_V_C_Th   # R448 T cell death in the central compartment
    v449 <- k_T0_death * x_V_P_Th   # R449 T cell death in the peripheral compartment
    v450 <- k_T0_death * x_V_LN_Th   # R450 T cell death in the lymph node compartment
    v451 <- k_T0_death * x_V_T_Th   # R451 T cell death in the tumor compartment
    v453 <- q_T0_P_in * x_V_C_Th   # R453 T cell transport into the peripheral compartment
    v454 <- q_T0_P_out * x_V_P_Th   # R454 T cell transport out of the peripheral compartment
    v456 <- q_T0_LN_out * x_V_LN_Th   # R456 T cell transport out of the lymph node compartment
    v457 <- k_IL2_sec * x_V_LN_aTh   # R457 IL2 secretion from activated T cells
    v459 <- k_TGFb_Tsec * x_V_T_T0   # R459 TGFb secretion from activated T cells in tumor
    v460 <- (k_IFNg_Thsec * x_V_T_Th)/NAVG_nmol   # R460 IFNg secretion from CD4 T cells in tumor
    v465 <- k_MDSC_death * x_V_T_MDSC   # R465 MDSC death in the tumor compartment
    v469 <- k_NO_sec * x_V_T_MDSC   # R469 NO secretion from MDSC
    v470 <- k_ArgI_sec * x_V_T_MDSC   # R470 ArgI secretion from MDSC
    v474 <- k_Mac_death * x_V_T_Mac_M1   # R474 Death of macrophage in tumor
    v475 <- k_Mac_death * x_V_T_Mac_M2   # R475 Death of macrophage in tumor
    v476 <- k_IL12_sec * x_V_T_mAPC   # R476 Secretion of IL-12 by APC in the tumor
    v477 <- k_IL12_Msec * x_V_T_Mac_M1   # R477 Secretion of IL-12 by macrophage in tumor
    v479 <- k_TGFb_Msec * x_V_T_Mac_M2   # R479 Secretion of TGFb by macrophage in tumor
    v480 <- k_IL10_sec * x_V_T_Mac_M2   # R480 Secretion of IL-10 by macrophage in tumor
    v488 <- (kon_CD47_SIRPa * x_syn_M_C_CD47 * x_syn_M_C_SIRPa * vol_syn_M_C)/NAVG_nmol - koff_CD47_SIRPa * x_syn_M_C_CD47_SIRPa * vol_syn_M_C   # R488 binding and unbinding of CD47-SIRPa in synapse
    v489 <- (kon_PD1_PDL1 * x_syn_M_C_PD1 * x_syn_M_C_PDL1 * vol_syn_M_C)/NAVG_nmol - koff_PD1_PDL1 * x_syn_M_C_PD1_PDL1 * vol_syn_M_C   # R489 binding and unbinding of PD1 PDL1 in synapse
    v490 <- (kon_PD1_PDL2 * x_syn_M_C_PD1 * x_syn_M_C_PDL2 * vol_syn_M_C)/NAVG_nmol - koff_PD1_PDL2 * x_syn_M_C_PD1_PDL2 * vol_syn_M_C   # R490 binding and unbinding of PD1 PDL2 in synapse
    v492 <- (Chi_PD1_aPD1 * kon_PD1_aPD1 * x_syn_M_C_PD1 * x_syn_M_C_PD1_aPD1 * vol_syn_M_C)/NAVG_nmol - 2 * koff_PD1_aPD1 * x_syn_M_C_PD1_aPD1_PD1 * vol_syn_M_C   # R492 binding and unbinding of PD1 to PD1-aPD1 on macrophage surface in synapse
    v494 <- (Chi_PDL1_aPDL1 * kon_PDL1_aPDL1 * x_syn_M_C_PDL1 * x_syn_M_C_PDL1_aPDL1 * vol_syn_M_C)/NAVG_nmol - 2 * koff_PDL1_aPDL1 * x_syn_M_C_PDL1_aPDL1_PDL1 * vol_syn_M_C   # R494 binding and unbinding of PDL1 to PDL1-aPDL1 on cancer cell surface in synapse
    v495 <- (kon_CD80_CD80 * x_syn_M_C_CD80m * x_syn_M_C_CD80m * vol_syn_M_C)/NAVG_nmol - koff_CD80_CD80 * x_syn_M_C_CD80 * vol_syn_M_C   # R495 self-association and dissociation of CD80 monomers in synapse
    v496 <- (kon_CD80_PDL1 * x_syn_M_C_CD80m * x_syn_M_C_PDL1 * vol_syn_M_C)/NAVG_nmol - koff_CD80_PDL1 * x_syn_M_C_PDL1_CD80 * vol_syn_M_C   # R496 binding and unbinding of CD80 and PDL1 in synapse
    H_P0 <- pTCR_p0_MHC_tot/(pTCR_p0_MHC_tot + p0_50)   # ESM 2 Rules (H_P0)
    H_P1 <- pTCR_p1_MHC_tot/(pTCR_p1_MHC_tot + p1_50)   # ESM 2 Rules (H_P1)
    N_aTh <- N0 + N_costim * H_CD28_APC + N_IL2_CD4 * x_V_LN_IL2/(IL2_50 + x_V_LN_IL2)   # ESM 2 Rules (N_aTh)
    H_Mac_C <- 1 - (1 - H_SIRPa) * (1 - H_PD1_M)   # ESM 2 Rules (H_Mac_C)
    vol_V_T <- (((x_V_T_C_x + C_total) * vol_cell + (x_V_T_T1_exh + x_V_T_Th_exh + T_total) * vol_Tcell)/Ve_T + M_total * vol_Mcell/Ve_T)   # ESM 2 Rules (V_T)
    N_aT <- N0 + N_costim * H_CD28_APC + N_IL2_CD8 * x_V_LN_IL2/(IL2_50 + x_V_LN_IL2)   # ESM 2 Rules (N_aT)
    N_aT0 <- N0 + N_costim * H_CD28_APC + N_IL2_CD4 * x_V_LN_IL2/(IL2_50 + x_V_LN_IL2)   # ESM 2 Rules (N_aT0)
    v4 <- k_C1_growth * x_V_T_C1 * (1 - C_total/C_max)   # R4 Cancer cell growth
    v24 <- k_cell_clear * x_V_T_T0 * (Kc_rec/(C_total^2 + Kc_rec))   # R24 T cell clearance upon antigen clearance
    v50 <- k_cell_clear * x_V_T_T1 * (Kc_rec/(C_total^2 + Kc_rec))   # R50 T cell clearance upon antigen clearance
    v52 <- k_T1 * x_V_T_T1 * C_total/(C_total + x_V_T_T1 + cell) * H_PD1_C1   # R52 T cell death from cancer
    v353 <- k_in_PDL1 * (C1_PDL1_base/A_cell - x_syn_T_C1_PDL1_total) * vol_syn_T_C1   # R353 Translocation of PDL1 between cell surface and cytoplasm
    v354 <- k_in_PDL1 * (C1_PDL1_base/A_cell * r_PDL2C1 - x_syn_T_C1_PDL2_total) * vol_syn_T_C1   # R354 Translocation of PDL2 between cell surface and cytoplasm
    v417 <- k_in_PDL1 * (APC_PDL1_base/A_cell - x_syn_T_APC_PDL1_total) * vol_syn_T_APC   # R417 Translocation of PDL1 between cell surface and cytoplasm
    v418 <- k_in_PDL1 * (APC_PDL1_base/A_cell * r_PDL2APC - x_syn_T_APC_PDL2_total) * vol_syn_T_APC   # R418 Translocation of PDL2 between cell surface and cytoplasm
    v440 <- k_in_PDL1 * (C1_PDL1_base/A_cell - x_nonsyn_C1_T_PDL1_total) * vol_nonsyn_C1_T   # R440 Translocation of PDL1 between cell surface and cytoplasm
    v442 <- k_in_PDL1 * (APC_PDL1_base/A_cell - x_nonsyn_APC_LN_PDL1_total) * vol_nonsyn_APC_LN   # R442 Translocation of PDL1 between cell surface and cytoplasm
    v452 <- k_cell_clear * x_V_T_Th * (Kc_rec/(C_total^2 + Kc_rec))   # R452 T cell clearance upon antigen clearance
    v462 <- k_CCL2_sec * C_total   # R462 CCL2 secretion by MDSCs and cancer cells
    v466 <- k_cell_clear * x_V_T_MDSC * (Kc_rec/(C_total^2 + Kc_rec))   # R466 MDSC death upon tumor eradication
    v472 <- k_cell_clear * x_V_T_Mac_M1 * (Kc_rec/(C_total^2 + Kc_rec))   # R472 Macrophage death upon tumor eradication
    v473 <- k_cell_clear * x_V_T_Mac_M2 * (Kc_rec/(C_total^2 + Kc_rec))   # R473 Macrophage death upon tumor eradication
    v486 <- k_in_PDL1 * (C1_PDL1_base/A_cell - x_syn_M_C_PDL1_total) * vol_syn_M_C   # R486 Translocation of PDL1 between cell surface and cytoplasm
    v487 <- k_in_PDL1 * (C1_PDL1_base/A_cell * r_PDL2C1 - x_syn_M_C_PDL2_total) * vol_syn_M_C   # R487 Translocation of PDL2 between cell surface and cytoplasm
    x_V_T_IFNg <- q_V_T_IFNg/vol_V_T
    x_V_T_P0 <- q_V_T_P0/vol_V_T
    x_V_T_P1 <- q_V_T_P1/vol_V_T
    x_V_T_aPDL1cc <- q_V_T_aPDL1cc/vol_V_T
    x_V_T_aPDL1co <- q_V_T_aPDL1co/vol_V_T
    x_V_T_aPDL1oo <- q_V_T_aPDL1oo/vol_V_T
    x_V_T_aPDL1c <- q_V_T_aPDL1c/vol_V_T
    x_V_T_aPDL1o <- q_V_T_aPDL1o/vol_V_T
    x_V_T_aPDL1 <- q_V_T_aPDL1/vol_V_T
    x_V_T_aPD1cc <- q_V_T_aPD1cc/vol_V_T
    x_V_T_aPD1co <- q_V_T_aPD1co/vol_V_T
    x_V_T_aPD1oo <- q_V_T_aPD1oo/vol_V_T
    x_V_T_aPD1c <- q_V_T_aPD1c/vol_V_T
    x_V_T_aPD1o <- q_V_T_aPD1o/vol_V_T
    x_V_T_aPD1 <- q_V_T_aPD1/vol_V_T
    x_V_T_aCTLA4cc <- q_V_T_aCTLA4cc/vol_V_T
    x_V_T_aCTLA4co <- q_V_T_aCTLA4co/vol_V_T
    x_V_T_aCTLA4oo <- q_V_T_aCTLA4oo/vol_V_T
    x_V_T_aCTLA4c <- q_V_T_aCTLA4c/vol_V_T
    x_V_T_aCTLA4o <- q_V_T_aCTLA4o/vol_V_T
    x_V_T_aCTLA4 <- q_V_T_aCTLA4/vol_V_T
    x_V_T_TGFb <- q_V_T_TGFb/vol_V_T
    x_V_T_NO <- q_V_T_NO/vol_V_T
    x_V_T_ArgI <- q_V_T_ArgI/vol_V_T
    x_V_T_CCL2 <- q_V_T_CCL2/vol_V_T
    x_V_T_IL12 <- q_V_T_IL12/vol_V_T
    x_V_T_IL10 <- q_V_T_IL10/vol_V_T
    v16 <- k_T0_act * H_APC * H_P0 * x_V_LN_nT0   # R16 Naive T cell activation
    v17 <- k_T0_act * H_APC * H_P0 * x_V_LN_nT0 * n_T0_clones   # R17 Naive T cell activation
    v18 <- k_T0_pro/N_aT0 * x_V_LN_aT0   # R18 aT0 cell proliferation
    v19 <- k_T0_pro/N_aT0 * 2^N_aT0 * x_V_LN_aT0   # R19 aT0 cell proliferation
    v27 <- q_T0_T_in * vol_V_T * x_V_C_T0 * (C_total^2/(C_total^2 + Kc_rec))   # R27 T cell transport into the tumor compartment
    v42 <- k_T1_act * H_mAPC * H_P1 * x_V_LN_nT1   # R42 Naive T cell activation
    v43 <- k_T1_act * H_mAPC * H_P1 * x_V_LN_nT1 * n_T1_clones   # R43 Naive T cell activation
    v44 <- k_T1_pro/N_aT * x_V_LN_aT1   # R44 aT1 cell proliferation
    v45 <- k_T1_pro/N_aT * 2^N_aT * x_V_LN_aT1   # R45 aT1 cell proliferation
    v55 <- q_T1_T_in * vol_V_T * x_V_C_T1 * (C_total^2/(C_total^2 + Kc_rec))   # R55 T cell transport into the tumor compartment
    v60 <- k_APC_death * (APC0_T * vol_V_T - x_V_T_APC)   # R60 APC recruitment/death in the tumour
    v113 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1cc/gamma_LN_aPDL1   # R113 aPDL1_cc convective transport from lymph node to central
    v114 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1co/gamma_LN_aPDL1   # R114 aPDL1_co convective transport from lymph node to central
    v115 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1oo/gamma_LN_aPDL1   # R115 aPDL1_oo convective transport from lymph node to central
    v116 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1c/gamma_LN_aPDL1   # R116 aPDL1_c convective transport from lymph node to central
    v117 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1o/gamma_LN_aPDL1   # R117 aPDL1_o convective transport from lymph node to central
    v118 <- q_LD_aPDL1 * vol_V_T * x_V_LN_aPDL1/gamma_LN_aPDL1   # R118 aPDL1 convective transport from lymph node to central
    v191 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1cc/gamma_LN_aPD1   # R191 aPD1_cc convective transport from lymph node to central
    v192 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1co/gamma_LN_aPD1   # R192 aPD1_co convective transport from lymph node to central
    v193 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1oo/gamma_LN_aPD1   # R193 aPD1_oo convective transport from lymph node to central
    v194 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1c/gamma_LN_aPD1   # R194 aPD1_c convective transport from lymph node to central
    v195 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1o/gamma_LN_aPD1   # R195 aPD1_o convective transport from lymph node to central
    v196 <- q_LD_aPD1 * vol_V_T * x_V_LN_aPD1/gamma_LN_aPD1   # R196 aPD1 convective transport from lymph node to central
    v263 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4cc/gamma_LN_aCTLA4   # R263 aCTLA4_cc convective transport from lymph node to central
    v264 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4co/gamma_LN_aCTLA4   # R264 aCTLA4_co convective transport from lymph node to central
    v265 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4oo/gamma_LN_aCTLA4   # R265 aCTLA4_oo convective transport from lymph node to central
    v266 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4c/gamma_LN_aCTLA4   # R266 aCTLA4_c convective transport from lymph node to central
    v267 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4o/gamma_LN_aCTLA4   # R267 aCTLA4_o convective transport from lymph node to central
    v268 <- q_LD_aCTLA4 * vol_V_T * x_V_LN_aCTLA4/gamma_LN_aCTLA4   # R268 aCTLA4 convective transport from lymph node to central
    v443 <- k_Th_act * H_APCh * H_P1 * x_V_LN_nT0   # R443 Naive T cell activation
    v444 <- k_Th_act * H_APCh * H_P1 * x_V_LN_nT0 * n_T1_clones   # R444 Naive T cell activation
    v445 <- k_T0_pro/N_aTh * x_V_LN_aTh   # R445 aTh cell proliferation
    v446 <- k_T0_pro/N_aTh * 2^N_aTh * x_V_LN_aTh   # R446 aTh cell proliferation
    v455 <- q_T0_T_in * vol_V_T * x_V_C_Th * (C_total^2/(C_total^2 + Kc_rec))   # R455 T cell transport into the tumor compartment
    H_TGFb <- x_V_T_TGFb/(x_V_T_TGFb + TGFb_50)   # ESM 2 Rules (H_TGFb)
    H_TGFb_Teff <- x_V_T_TGFb/(x_V_T_TGFb + TGFb_50_Teff)   # ESM 2 Rules (H_TGFb_Teff)
    H_NO <- x_V_T_NO/(NO_50_Teff + x_V_T_NO)   # ESM 2 Rules (H_NO)
    H_ArgI_Teff <- x_V_T_ArgI/(ArgI_50_Teff + x_V_T_ArgI)   # ESM 2 Rules (H_ArgI_Teff)
    H_ArgI_Treg <- x_V_T_ArgI/(ArgI_50_Treg + x_V_T_ArgI)   # ESM 2 Rules (H_ArgI_Treg)
    H_IL10 <- x_V_T_IL10/(IL10_50 + x_V_T_IL10)   # ESM 2 Rules (H_IL10)
    H_IL10_phago <- x_V_T_IL10/(IL10_50_phago + x_V_T_IL10)   # ESM 2 Rules (H_IL10_phago)
    H_IL12 <- x_V_T_IL12/(IL12_50 + x_V_T_IL12)   # ESM 2 Rules (H_IL12)
    v70 <- k_xP0_deg * x_V_T_P0 * vol_V_T   # R70 Free antigen degradation
    v71 <- k_P0_up * x_V_T_APC * x_V_T_P0 * vol_V_T   # R71 Antigen uptake by mature antigen presenting cells
    v72 <- k_P0_up * cell * x_V_T_P0 * vol_V_e   # R72 Antigen uptake by mature antigen presenting cells
    v80 <- k_xP1_deg * x_V_T_P1 * vol_V_T   # R80 Free antigen degradation
    v81 <- k_P1_up * x_V_T_APC * x_V_T_P1 * vol_V_T   # R81 Antigen uptake by mature antigen presenting cells
    v82 <- k_P1_up * cell * x_V_T_P1 * vol_V_e   # R82 Antigen uptake by mature antigen presenting cells
    v95 <- q_T_aPDL1 * (x_V_C_aPDL1cc/gamma_C_aPDL1 - x_V_T_aPDL1cc/gamma_T_aPDL1)   # R95 aPDL1_cc diffusive transport to tumor compartment
    v96 <- q_T_aPDL1 * (x_V_C_aPDL1co/gamma_C_aPDL1 - x_V_T_aPDL1co/gamma_T_aPDL1)   # R96 aPDL1_co diffusive transport to tumor compartment
    v97 <- q_T_aPDL1 * (x_V_C_aPDL1c/gamma_C_aPDL1 - x_V_T_aPDL1c/gamma_T_aPDL1)   # R97 aPDL1_c diffusive transport to tumor compartment
    v98 <- q_T_aPDL1 * (x_V_C_aPDL1o/gamma_C_aPDL1 - x_V_T_aPDL1o/gamma_T_aPDL1)   # R98 aPDL1_o diffusive transport to tumor compartment
    v99 <- q_T_aPDL1 * (x_V_C_aPDL1/gamma_C_aPDL1 - x_V_T_aPDL1/gamma_T_aPDL1)   # R99 aPDL1 diffusive transport to tumor compartment
    v100 <- q_T_aPDL1 * (x_V_C_aPDL1oo/gamma_C_aPDL1 - x_V_T_aPDL1oo/gamma_T_aPDL1)   # R100 aPDL1_oo diffusive transport to tumor compartment
    v107 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1cc/gamma_T_aPDL1   # R107 aPDL1_cc convective transport from tumor to lymph node
    v108 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1co/gamma_T_aPDL1   # R108 aPDL1_co convective transport from tumor to lymph node
    v109 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1oo/gamma_T_aPDL1   # R109 aPDL1_oo convective transport from tumor to lymph node
    v110 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1c/gamma_T_aPDL1   # R110 aPDL1_c convective transport from tumor to lymph node
    v111 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1o/gamma_T_aPDL1   # R111 aPDL1_o convective transport from tumor to lymph node
    v112 <- q_LD_aPDL1 * vol_V_T * x_V_T_aPDL1/gamma_T_aPDL1   # R112 aPDL1 convective transport from tumor to lymph node
    v125 <- km * (2 * x_V_T_aPDL1cc - x_V_T_aPDL1co/K_m) * vol_V_T   # R125 cc mask opening to co in T
    v126 <- 2 * k_cvg * x_V_T_aPDL1cc * vol_V_T   # R126 cc mask cleavage to c in T
    v127 <- km * (x_V_T_aPDL1co - 2 * x_V_T_aPDL1oo/K_m) * vol_V_T   # R127 co mask opening to oo in T
    v128 <- k_cvg * x_V_T_aPDL1co * vol_V_T   # R128 co mask cleavage to c in T
    v129 <- k_cvg * x_V_T_aPDL1co * vol_V_T   # R129 co mask cleavage to o in T
    v130 <- 2 * k_cvg * x_V_T_aPDL1oo * vol_V_T   # R130 oo mask cleavage to o in T
    v131 <- km * (x_V_T_aPDL1c - x_V_T_aPDL1o/K_m) * vol_V_T   # R131 c mask opening to o in T
    v132 <- k_cvg * x_V_T_aPDL1c * vol_V_T   # R132 c mask cleavage to free in T
    v133 <- k_cvg * x_V_T_aPDL1o * vol_V_T   # R133 o mask cleavage to free in T
    v173 <- q_T_aPD1 * (x_V_C_aPD1cc/gamma_C_aPD1 - x_V_T_aPD1cc/gamma_T_aPD1)   # R173 aPD1_cc diffusive transport to tumor compartment
    v174 <- q_T_aPD1 * (x_V_C_aPD1co/gamma_C_aPD1 - x_V_T_aPD1co/gamma_T_aPD1)   # R174 aPD1_co diffusive transport to tumor compartment
    v175 <- q_T_aPD1 * (x_V_C_aPD1c/gamma_C_aPD1 - x_V_T_aPD1c/gamma_T_aPD1)   # R175 aPD1_c diffusive transport to tumor compartment
    v176 <- q_T_aPD1 * (x_V_C_aPD1o/gamma_C_aPD1 - x_V_T_aPD1o/gamma_T_aPD1)   # R176 aPD1_o diffusive transport to tumor compartment
    v177 <- q_T_aPD1 * (x_V_C_aPD1/gamma_C_aPD1 - x_V_T_aPD1/gamma_T_aPD1)   # R177 aPD1 diffusive transport to tumor compartment
    v178 <- q_T_aPD1 * (x_V_C_aPD1oo/gamma_C_aPD1 - x_V_T_aPD1oo/gamma_T_aPD1)   # R178 aPD1_oo diffusive transport to tumor compartment
    v185 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1cc/gamma_T_aPD1   # R185 aPD1_cc convective transport from tumor to lymph node
    v186 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1co/gamma_T_aPD1   # R186 aPD1_co convective transport from tumor to lymph node
    v187 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1oo/gamma_T_aPD1   # R187 aPD1_oo convective transport from tumor to lymph node
    v188 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1c/gamma_T_aPD1   # R188 aPD1_c convective transport from tumor to lymph node
    v189 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1o/gamma_T_aPD1   # R189 aPD1_o convective transport from tumor to lymph node
    v190 <- q_LD_aPD1 * vol_V_T * x_V_T_aPD1/gamma_T_aPD1   # R190 aPD1 convective transport from tumor to lymph node
    v203 <- km * (2 * x_V_T_aPD1cc - x_V_T_aPD1co/K_m) * vol_V_T   # R203 cc mask opening to co in T
    v204 <- 2 * k_cvg * x_V_T_aPD1cc * vol_V_T   # R204 cc mask cleavage to c in T
    v205 <- km * (x_V_T_aPD1co - 2 * x_V_T_aPD1oo/K_m) * vol_V_T   # R205 co mask opening to oo in T
    v206 <- k_cvg * x_V_T_aPD1co * vol_V_T   # R206 co mask cleavage to c in T
    v207 <- k_cvg * x_V_T_aPD1co * vol_V_T   # R207 co mask cleavage to o in T
    v208 <- 2 * k_cvg * x_V_T_aPD1oo * vol_V_T   # R208 oo mask cleavage to o in T
    v209 <- km * (x_V_T_aPD1c - x_V_T_aPD1o/K_m) * vol_V_T   # R209 c mask opening to o in T
    v210 <- k_cvg * x_V_T_aPD1c * vol_V_T   # R210 c mask cleavage to free in T
    v211 <- k_cvg * x_V_T_aPD1o * vol_V_T   # R211 o mask cleavage to free in T
    v245 <- q_T_aCTLA4 * (x_V_C_aCTLA4cc/gamma_C_aCTLA4 - x_V_T_aCTLA4cc/gamma_T_aCTLA4)   # R245 aCTLA4_cc diffusive transport to tumor compartment
    v246 <- q_T_aCTLA4 * (x_V_C_aCTLA4co/gamma_C_aCTLA4 - x_V_T_aCTLA4co/gamma_T_aCTLA4)   # R246 aCTLA4_co diffusive transport to tumor compartment
    v247 <- q_T_aCTLA4 * (x_V_C_aCTLA4c/gamma_C_aCTLA4 - x_V_T_aCTLA4c/gamma_T_aCTLA4)   # R247 aCTLA4_c diffusive transport to tumor compartment
    v248 <- q_T_aCTLA4 * (x_V_C_aCTLA4o/gamma_C_aCTLA4 - x_V_T_aCTLA4o/gamma_T_aCTLA4)   # R248 aCTLA4_o diffusive transport to tumor compartment
    v249 <- q_T_aCTLA4 * (x_V_C_aCTLA4/gamma_C_aCTLA4 - x_V_T_aCTLA4/gamma_T_aCTLA4)   # R249 aCTLA4 diffusive transport to tumor compartment
    v250 <- q_T_aCTLA4 * (x_V_C_aCTLA4oo/gamma_C_aCTLA4 - x_V_T_aCTLA4oo/gamma_T_aCTLA4)   # R250 aCTLA4_oo diffusive transport to tumor compartment
    v257 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4cc/gamma_T_aCTLA4   # R257 aCTLA4_cc convective transport from tumor to lymph node
    v258 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4co/gamma_T_aCTLA4   # R258 aCTLA4_co convective transport from tumor to lymph node
    v259 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4oo/gamma_T_aCTLA4   # R259 aCTLA4_oo convective transport from tumor to lymph node
    v260 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4c/gamma_T_aCTLA4   # R260 aCTLA4_c convective transport from tumor to lymph node
    v261 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4o/gamma_T_aCTLA4   # R261 aCTLA4_o convective transport from tumor to lymph node
    v262 <- q_LD_aCTLA4 * vol_V_T * x_V_T_aCTLA4/gamma_T_aCTLA4   # R262 aCTLA4 convective transport from tumor to lymph node
    v275 <- km * (2 * x_V_T_aCTLA4cc - x_V_T_aCTLA4co/K_m) * vol_V_T   # R275 cc mask opening to co in T
    v276 <- 2 * k_cvg * x_V_T_aCTLA4cc * vol_V_T   # R276 cc mask cleavage to c in T
    v277 <- km * (x_V_T_aCTLA4co - 2 * x_V_T_aCTLA4oo/K_m) * vol_V_T   # R277 co mask opening to oo in T
    v278 <- k_cvg * x_V_T_aCTLA4co * vol_V_T   # R278 co mask cleavage to c in T
    v279 <- k_cvg * x_V_T_aCTLA4co * vol_V_T   # R279 co mask cleavage to o in T
    v280 <- 2 * k_cvg * x_V_T_aCTLA4oo * vol_V_T   # R280 oo mask cleavage to o in T
    v281 <- km * (x_V_T_aCTLA4c - x_V_T_aCTLA4o/K_m) * vol_V_T   # R281 c mask opening to o in T
    v282 <- k_cvg * x_V_T_aCTLA4c * vol_V_T   # R282 c mask cleavage to free in T
    v283 <- k_cvg * x_V_T_aCTLA4o * vol_V_T   # R283 o mask cleavage to free in T
    v311 <- kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1co/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1cRo * vol_syn_T_C1   # R311 binding and unbinding of PDL1 to aPDL1co on C1 surface in synapse
    v312 <- kon_PDL1_aPDL1 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1co/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1cRo * vol_syn_T_C1   # R312 binding and unbinding of PDL1 to aPDL1co on C1 surface in synapse
    v313 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1oo/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1oRo * vol_syn_T_C1   # R313 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v315 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1oo/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1oRo * vol_syn_T_C1   # R315 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v319 <- kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1c/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1Rmc * vol_syn_T_C1   # R319 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v321 <- kon_PDL1_aPDL1 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1c/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1Rmc * vol_syn_T_C1   # R321 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v323 <- kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1o/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1Rmo * vol_syn_T_C1   # R323 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v326 <- kon_PDL1_aPDL1 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1o/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1Rmo * vol_syn_T_C1   # R326 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v329 <- kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1o/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1mRo * vol_syn_T_C1   # R329 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v331 <- kon_PDL1_aPDL1 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1o/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1mRo * vol_syn_T_C1   # R331 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v339 <- 2 * kon_PDL1_aPDL1 * (x_syn_T_C1_PDL1 * x_V_T_aPDL1/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_PDL1_aPDL1Rmm * vol_syn_T_C1   # R339 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v343 <- 2 * kon_PDL1_aPDL1 * 0 * (x_syn_T_C1_TPDL1 * x_V_T_aPDL1/gamma_T_aPDL1) * vol_syn_T_C1 - koff_PDL1_aPDL1 * x_syn_T_C1_TPDL1_aPDL1Rmm * vol_syn_T_C1   # R343 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    v351 <- k_out_PDL1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_T_C1_PDL1_total/(C1_PDL1_base * r_PDL1_IFNg/A_cell))   # R351 Translocation of PDL1 between cell surface and cytoplasm
    v352 <- k_out_PDL1 * r_PDL2C1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_T_C1_PDL2_total/(C1_PDL1_base * r_PDL1_IFNg/A_cell * r_PDL2C1))   # R352 Translocation of PDL2 between cell surface and cytoplasm
    v357 <- 2 * kon_PD1_aPD1 * (x_syn_T_C1_PD1 * x_V_T_aPD1/gamma_T_aPD1) * vol_syn_T_C1 - koff_PD1_aPD1 * x_syn_T_C1_PD1_aPD1 * vol_syn_T_C1   # R357 binding and unbinding of PD1 to aPD1 on T surface in synapse
    v369 <- 4 * kon_CTLA4_aCTLA4 * (x_syn_T_C1_CTLA4 * x_V_T_aCTLA4/gamma_T_aCTLA4) * vol_syn_T_C1 - koff_CTLA4_aCTLA4 * x_syn_T_C1_CTLA4_aCTLA4 * vol_syn_T_C1   # R369 binding and unbinding of CTLA4 to aCTLA4 on T surface in synapse
    v415 <- k_out_PDL1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_T_APC_PDL1_total/(APC_PDL1_base * r_PDL1_IFNg/A_cell))   # R415 Translocation of PDL1 between cell surface and cytoplasm
    v416 <- k_out_PDL1 * r_PDL2APC * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_T_APC_PDL2_total/(APC_PDL1_base * r_PDL1_IFNg/A_cell * r_PDL2APC))   # R416 Translocation of PDL2 between cell surface and cytoplasm
    v439 <- k_out_PDL1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_nonsyn_C1_T_PDL1_total/(C1_PDL1_base * r_PDL1_IFNg/A_cell))   # R439 Translocation of PDL1 between cell surface and cytoplasm
    v441 <- k_out_PDL1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_nonsyn_APC_LN_PDL1_total/(APC_PDL1_base * r_PDL1_IFNg/A_cell))   # R441 Translocation of PDL1 between cell surface and cytoplasm
    v458 <- k_TGFb_deg * (TGFbase - x_V_T_TGFb) * vol_V_T   # R458 TGFb secretion by triple-negative breast cancer cells
    v461 <- k_IFNg_deg * x_V_T_IFNg * vol_V_T   # R461 IFNg degradtion in tumor
    v463 <- k_CCL2_deg * x_V_T_CCL2 * vol_V_T   # R463 CCL2 degradation
    v464 <- k_MDSC_rec * vol_V_T * (x_V_T_CCL2/(x_V_T_CCL2 + CCL2_50))   # R464 Recruitment of macrophage to tumor
    v467 <- k_NO_deg * x_V_T_NO * vol_V_T   # R467 NO degradation
    v468 <- k_ArgI_deg * x_V_T_ArgI * vol_V_T   # R468 ArgI degradation
    v471 <- k_Mac_rec * vol_V_T * (x_V_T_CCL2/(x_V_T_CCL2 + CCL2_50))   # R471 Recruitment of pro-inflammatory M1 macrophage to tumor
    v478 <- k_IL12_deg * x_V_T_IL12 * vol_V_T   # R478 Degradation of IL-12 in tumor
    v481 <- k_IL10_deg * x_V_T_IL10 * vol_V_T   # R481 Degradation of IL-10 in tumor
    v482 <- k_M2_pol * x_V_T_Mac_M1 * (x_V_T_TGFb/(x_V_T_TGFb + TGFb_50) + x_V_T_IL10/(x_V_T_IL10 + IL10_50))   # R482 M1 to M2 macrophage polarization in tumor
    v483 <- k_M1_pol * x_V_T_Mac_M2 * (x_V_T_IL12/(x_V_T_IL12 + IL12_50) + x_V_T_IFNg/(x_V_T_IFNg + IFNg_50))   # R483 M2 to M1 macrophage polarization in tumor
    v484 <- k_out_PDL1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_M_C_PDL1_total/(C1_PDL1_base * r_PDL1_IFNg/A_cell))   # R484 Translocation of PDL1 between cell surface and cytoplasm
    v485 <- k_out_PDL1 * r_PDL2C1 * x_V_T_IFNg/(x_V_T_IFNg + IFNg_50_ind) * (1 - x_syn_M_C_PDL2_total/(C1_PDL1_base * r_PDL1_IFNg/A_cell * r_PDL2C1))   # R485 Translocation of PDL2 between cell surface and cytoplasm
    v491 <- 2 * kon_PD1_aPD1 * (x_syn_M_C_PD1 * x_V_T_aPD1/gamma_T_aPD1) * vol_syn_M_C - koff_PD1_aPD1 * x_syn_M_C_PD1_aPD1 * vol_syn_M_C   # R491 binding and unbinding of PD1 to aPD1 on macrophage surface in synapse
    v493 <- 2 * kon_PDL1_aPDL1 * (x_syn_M_C_PDL1 * x_V_T_aPDL1/gamma_T_aPDL1) * vol_syn_M_C - koff_PDL1_aPDL1 * x_syn_M_C_PDL1_aPDL1 * vol_syn_M_C   # R493 binding and unbinding of PDL1 to aPDL1 on C1 surface in synapse
    H_MDSC <- 1 - (1 - H_NO) * (1 - H_ArgI_Teff)   # ESM 2 Rules (H_MDSC)
    v51 <- k_Treg * x_V_T_T1 * Tregs_/(x_V_T_T1 + Tregs_ + cell) * H_IL10   # R51 T cell death from Tregs
    v62 <- k_APC_mat * x_V_T_APC * H_IL12 * (1 - H_IL10)   # R62 APC maturation in the tumour
    v447 <- k_Th_Treg * x_V_T_Th * H_TGFb * H_ArgI_Treg   # R447 Differentiation between Treg and Th Cells
    v497 <- k_M1_phago * x_V_T_C1 * x_V_T_Mac_M1/(x_V_T_Mac_M1 + K_Mac_C * x_V_T_C1 + cell) * (1 - H_Mac_C) * (1 - H_IL10_phago)   # R497 M1 macrophage-mediated killing of cancer cell
    R_Tcell <- 0 * cell/day + (k_C1_death + k_C1_therapy) * x_V_T_C1 + k_C_T1 * x_V_T_C1 * x_V_T_T1/(K_T_C * x_V_T_C1 + x_V_T_T1 + cell) * x_V_T_T1/(x_V_T_T1 + K_T_Treg * Tregs_ + cell) * (1 - H_TGFb_Teff) * (1 - H_PD1_C1) * (1 - H_MDSC) + k_M1_phago * x_V_T_C1 * x_V_T_Mac_M1/(x_V_T_Mac_M1 + K_Mac_C * x_V_T_C1 + cell) * (1 - H_Mac_C) * (1 - H_IL10_phago)   # ESM 2 Rules (R_Tcell)
    v59 <- k_C_T1 * x_V_T_C1 * x_V_T_T1/(K_T_C * x_V_T_C1 + x_V_T_T1 + cell) * x_V_T_T1/(x_V_T_T1 + K_T_Treg * Tregs_ + cell) * (1 - H_TGFb_Teff) * (1 - H_PD1_C1) * (1 - H_MDSC)   # R59 Cancer cell killing by T cells
    v69 <- n_T0_clones * (P0_C1 * (k_C1_death + k_C1_therapy) * x_V_T_C1 + P0_C1 * k_C_T1 * x_V_T_C1 * x_V_T_T1/(K_T_C * x_V_T_C1 + x_V_T_T1 + cell) * x_V_T_T1/(x_V_T_T1 + K_T_Treg * Tregs_ + cell) * (1 - H_TGFb_Teff) * (1 - H_PD1_C1) * (1 - H_MDSC))   # R69 Antigen deposition from dying cancer cells
    v79 <- n_T1_clones * (P1_C1 * (k_C1_death + k_C1_therapy) * x_V_T_C1 + P1_C1 * k_C_T1 * x_V_T_C1 * x_V_T_T1/(K_T_C * x_V_T_C1 + x_V_T_T1 + cell) * x_V_T_T1/(x_V_T_T1 + K_T_Treg * Tregs_ + cell) * (1 - H_TGFb_Teff) * (1 - H_PD1_C1) * (1 - H_MDSC))   # R79 Antigen deposition from dying cancer cells

    # ---------- ODEs (states are AMOUNTS: nmol, molecule = cell, pg, mU*L) ----------
    d/dt(q_V_C_nT0) <- +v6-v10-v12+v13-v14+v15
    d/dt(q_V_C_T0) <- -v20-v25+v26-v27+v28
    d/dt(q_V_C_nT1) <- +v32-v36-v38+v39-v40+v41
    d/dt(q_V_C_T1) <- -v46-v53+v54-v55+v56
    d/dt(q_V_C_aPDL1cc) <- -v89-v95-v101+v113-v119-v134-v135-v161
    d/dt(q_V_C_aPDL1co) <- -v90-v96-v102+v114-v120+v134-v136-v137-v138-v162
    d/dt(q_V_C_aPDL1oo) <- -v94-v100-v106+v115-v121+v136-v139-v163
    d/dt(q_V_C_aPDL1c) <- -v91-v97-v103+v116-v122+v135+v137-v140-v141-v164
    d/dt(q_V_C_aPDL1o) <- -v92-v98-v104+v117-v123+v138+v139+v140-v142-v165
    d/dt(q_V_C_aPDL1) <- -v93-v99-v105+v118-v124+v141+v142-v166
    d/dt(q_V_C_aPD1cc) <- -v167-v173-v179+v191-v197-v212-v213
    d/dt(q_V_C_aPD1co) <- -v168-v174-v180+v192-v198+v212-v214-v215-v216
    d/dt(q_V_C_aPD1oo) <- -v172-v178-v184+v193-v199+v214-v217
    d/dt(q_V_C_aPD1c) <- -v169-v175-v181+v194-v200+v213+v215-v218-v219
    d/dt(q_V_C_aPD1o) <- -v170-v176-v182+v195-v201+v216+v217+v218-v220
    d/dt(q_V_C_aPD1) <- -v171-v177-v183+v196-v202+v219+v220
    d/dt(q_V_C_aCTLA4cc) <- -v239-v245-v251+v263-v269-v284-v285
    d/dt(q_V_C_aCTLA4co) <- -v240-v246-v252+v264-v270+v284-v286-v287-v288
    d/dt(q_V_C_aCTLA4oo) <- -v244-v250-v256+v265-v271+v286-v289
    d/dt(q_V_C_aCTLA4c) <- -v241-v247-v253+v266-v272+v285+v287-v290-v291
    d/dt(q_V_C_aCTLA4o) <- -v242-v248-v254+v267-v273+v288+v289+v290-v292
    d/dt(q_V_C_aCTLA4) <- -v243-v249-v255+v268-v274+v291+v292
    d/dt(q_V_C_Th) <- -v448-v453+v454-v455+v456
    d/dt(q_V_P_nT0) <- +v7-v9+v12-v13
    d/dt(q_V_P_T0) <- -v21+v25-v26
    d/dt(q_V_P_nT1) <- +v33-v35+v38-v39
    d/dt(q_V_P_T1) <- -v47+v53-v54
    d/dt(q_V_P_aPDL1cc) <- +v89-v143-v144
    d/dt(q_V_P_aPDL1co) <- +v90+v143-v145-v146-v147
    d/dt(q_V_P_aPDL1oo) <- +v94+v145-v148
    d/dt(q_V_P_aPDL1c) <- +v91+v144+v146-v149-v150
    d/dt(q_V_P_aPDL1o) <- +v92+v147+v148+v149-v151
    d/dt(q_V_P_aPDL1) <- +v93+v150+v151
    d/dt(q_V_P_aPD1cc) <- +v167-v221-v222
    d/dt(q_V_P_aPD1co) <- +v168+v221-v223-v224-v225
    d/dt(q_V_P_aPD1oo) <- +v172+v223-v226
    d/dt(q_V_P_aPD1c) <- +v169+v222+v224-v227-v228
    d/dt(q_V_P_aPD1o) <- +v170+v225+v226+v227-v229
    d/dt(q_V_P_aPD1) <- +v171+v228+v229
    d/dt(q_V_P_aCTLA4cc) <- +v239-v293-v294
    d/dt(q_V_P_aCTLA4co) <- +v240+v293-v295-v296-v297
    d/dt(q_V_P_aCTLA4oo) <- +v244+v295-v298
    d/dt(q_V_P_aCTLA4c) <- +v241+v294+v296-v299-v300
    d/dt(q_V_P_aCTLA4o) <- +v242+v297+v298+v299-v301
    d/dt(q_V_P_aCTLA4) <- +v243+v300+v301
    d/dt(q_V_P_Th) <- -v449+v453-v454
    d/dt(q_V_T_C_x) <- -v1+v5+v59+v497
    d/dt(q_V_T_T1_exh) <- -v2+v49+v50+v51+v52
    d/dt(q_V_T_Th_exh) <- -v3+v451+v452
    d/dt(q_V_T_C1) <- +v4-v5-v59-v497
    d/dt(q_V_T_T0) <- -v23-v24+v27+v447
    d/dt(q_V_T_T1) <- -v49-v50-v51-v52+v55
    d/dt(q_V_T_IFNg) <- +v58+v460-v461
    d/dt(q_V_T_APC) <- +v60-v62
    d/dt(q_V_T_mAPC) <- +v62-v63-v64
    d/dt(q_V_T_P0) <- +v69-v70-v71
    d/dt(q_V_T_P1) <- +v79-v80-v81
    d/dt(q_V_T_aPDL1cc) <- +v95-v107-v125-v126
    d/dt(q_V_T_aPDL1co) <- +v96-v108+v125-v127-v128-v129
    d/dt(q_V_T_aPDL1oo) <- +v100-v109+v127-v130
    d/dt(q_V_T_aPDL1c) <- +v97-v110+v126+v128-v131-v132
    d/dt(q_V_T_aPDL1o) <- +v98-v111+v129+v130+v131-v133
    d/dt(q_V_T_aPDL1) <- +v99-v112+v132+v133
    d/dt(q_V_T_aPD1cc) <- +v173-v185-v203-v204
    d/dt(q_V_T_aPD1co) <- +v174-v186+v203-v205-v206-v207
    d/dt(q_V_T_aPD1oo) <- +v178-v187+v205-v208
    d/dt(q_V_T_aPD1c) <- +v175-v188+v204+v206-v209-v210
    d/dt(q_V_T_aPD1o) <- +v176-v189+v207+v208+v209-v211
    d/dt(q_V_T_aPD1) <- +v177-v190+v210+v211
    d/dt(q_V_T_aCTLA4cc) <- +v245-v257-v275-v276
    d/dt(q_V_T_aCTLA4co) <- +v246-v258+v275-v277-v278-v279
    d/dt(q_V_T_aCTLA4oo) <- +v250-v259+v277-v280
    d/dt(q_V_T_aCTLA4c) <- +v247-v260+v276+v278-v281-v282
    d/dt(q_V_T_aCTLA4o) <- +v248-v261+v279+v280+v281-v283
    d/dt(q_V_T_aCTLA4) <- +v249-v262+v282+v283
    d/dt(q_V_T_Th) <- -v447-v451-v452+v455
    d/dt(q_V_T_TGFb) <- +v458+v459+v479
    d/dt(q_V_T_MDSC) <- +v464-v465-v466
    d/dt(q_V_T_NO) <- -v467+v469
    d/dt(q_V_T_ArgI) <- -v468+v470
    d/dt(q_V_T_CCL2) <- +v462-v463
    d/dt(q_V_T_Mac_M1) <- +v471-v472-v474-v482+v483
    d/dt(q_V_T_Mac_M2) <- -v473-v475+v482-v483
    d/dt(q_V_T_IL12) <- +v476+v477-v478
    d/dt(q_V_T_IL10) <- +v480-v481
    d/dt(q_V_LN_nT0) <- +v8-v11+v14-v15-v16-v443
    d/dt(q_V_LN_aT0) <- +v17-v18
    d/dt(q_V_LN_T0) <- +v19-v22-v28
    d/dt(q_V_LN_IL2) <- -v29-v30-v31+v57+v457
    d/dt(q_V_LN_nT1) <- +v34-v37+v40-v41-v42
    d/dt(q_V_LN_aT1) <- +v43-v44
    d/dt(q_V_LN_T1) <- +v45-v48-v56
    d/dt(q_V_LN_APC) <- +v61
    d/dt(q_V_LN_mAPC) <- +v63-v65
    d/dt(q_V_LN_aPDL1cc) <- +v101+v107-v113-v152-v153
    d/dt(q_V_LN_aPDL1co) <- +v102+v108-v114+v152-v154-v155-v156
    d/dt(q_V_LN_aPDL1oo) <- +v106+v109-v115+v154-v157
    d/dt(q_V_LN_aPDL1c) <- +v103+v110-v116+v153+v155-v158-v159
    d/dt(q_V_LN_aPDL1o) <- +v104+v111-v117+v156+v157+v158-v160
    d/dt(q_V_LN_aPDL1) <- +v105+v112-v118+v159+v160
    d/dt(q_V_LN_aPD1cc) <- +v179+v185-v191-v230-v231
    d/dt(q_V_LN_aPD1co) <- +v180+v186-v192+v230-v232-v233-v234
    d/dt(q_V_LN_aPD1oo) <- +v184+v187-v193+v232-v235
    d/dt(q_V_LN_aPD1c) <- +v181+v188-v194+v231+v233-v236-v237
    d/dt(q_V_LN_aPD1o) <- +v182+v189-v195+v234+v235+v236-v238
    d/dt(q_V_LN_aPD1) <- +v183+v190-v196+v237+v238
    d/dt(q_V_LN_aCTLA4cc) <- +v251+v257-v263-v302-v303
    d/dt(q_V_LN_aCTLA4co) <- +v252+v258-v264+v302-v304-v305-v306
    d/dt(q_V_LN_aCTLA4oo) <- +v256+v259-v265+v304-v307
    d/dt(q_V_LN_aCTLA4c) <- +v253+v260-v266+v303+v305-v308-v309
    d/dt(q_V_LN_aCTLA4o) <- +v254+v261-v267+v306+v307+v308-v310
    d/dt(q_V_LN_aCTLA4) <- +v255+v262-v268+v309+v310
    d/dt(q_V_LN_aTh) <- +v444-v445
    d/dt(q_V_LN_Th) <- +v446-v450-v456
    d/dt(q_V_e_P0) <- +v72-v73
    d/dt(q_V_e_p0) <- +v73-v74-v75/NAVG_nmol+v76/NAVG_nmol
    d/dt(q_V_e_P1) <- +v82-v83
    d/dt(q_V_e_p1) <- +v83-v84-v85/NAVG_nmol+v86/NAVG_nmol
    d/dt(q_A_e_M1) <- -v68-v75+v76-v85+v86
    d/dt(q_A_e_M1p0) <- +v75-v76-v78
    d/dt(q_A_e_M1p1) <- +v85-v86-v88
    d/dt(q_A_s_M1) <- +v68+v77+v87
    d/dt(q_A_s_M1p0) <- -v77+v78
    d/dt(q_A_s_M1p1) <- -v87+v88
    d/dt(q_syn_T_C1_PD1_PDL1) <- +v355
    d/dt(q_syn_T_C1_PD1_PDL2) <- +v356
    d/dt(q_syn_T_C1_PD1) <- -v355-v356-v357-v358
    d/dt(q_syn_T_C1_PDL1) <- -v311-v313-v317-v319-v323-v329-v334-v335-v339-v347+v351+v353-v355-v372
    d/dt(q_syn_T_C1_PDL2) <- +v352+v354-v356
    d/dt(q_syn_T_C1_PD1_aPD1) <- +v357-v358
    d/dt(q_syn_T_C1_PD1_aPD1_PD1) <- +v358
    d/dt(q_syn_T_C1_TPDL1) <- -v312-v315-v318-v321-v326-v331-v337-v338-v343-v349
    d/dt(q_syn_T_C1_PDL1_aPDL1cRo) <- +v311-v314-v320-v330
    d/dt(q_syn_T_C1_PDL1_aPDL1oRo) <- +v313+v314-v317-v325
    d/dt(q_syn_T_C1_PDL1_aPDL1Rmc) <- +v319+v320-v324-v342
    d/dt(q_syn_T_C1_PDL1_aPDL1mRo) <- +v329+v330-v334-v341
    d/dt(q_syn_T_C1_PDL1_aPDL1Rmo) <- +v323+v324+v325-v335-v340
    d/dt(q_syn_T_C1_PDL1_aPDL1Rmm) <- +v339+v340+v341+v342-v347
    d/dt(q_syn_T_C1_PDL1_aPDL1RoRo_PDL1) <- +v317-v333
    d/dt(q_syn_T_C1_PDL1_aPDL1RmRo_PDL1) <- +v333+v334+v335-v348
    d/dt(q_syn_T_C1_PDL1_aPDL1RmRm_PDL1) <- +v347+v348
    d/dt(q_syn_T_C1_TPDL1_aPDL1cRo) <- +v312-v316-v322-v332
    d/dt(q_syn_T_C1_TPDL1_aPDL1oRo) <- +v315+v316-v318-v328
    d/dt(q_syn_T_C1_TPDL1_aPDL1Rmo) <- +v326+v327+v328-v338-v344
    d/dt(q_syn_T_C1_TPDL1_aPDL1Rmc) <- +v321+v322-v327-v346
    d/dt(q_syn_T_C1_TPDL1_aPDL1mRo) <- +v331+v332-v337-v345
    d/dt(q_syn_T_C1_TPDL1_aPDL1Rmm) <- +v343+v344+v345+v346-v349
    d/dt(q_syn_T_C1_TPDL1_aPDL1RoRo_TPDL1) <- +v318-v336
    d/dt(q_syn_T_C1_TPDL1_aPDL1RmRo_TPDL1) <- +v336+v337+v338-v350
    d/dt(q_syn_T_C1_TPDL1_aPDL1RmRm_TPDL1) <- +v349+v350
    d/dt(q_syn_T_C1_CD28_CD80) <- +v359-v360
    d/dt(q_syn_T_C1_CD28_CD80_CD28) <- +v360
    d/dt(q_syn_T_C1_CD28_CD86) <- +v361
    d/dt(q_syn_T_C1_CD80_CTLA4) <- +v362-v363-v365
    d/dt(q_syn_T_C1_CD80_CTLA4_CD80) <- +v365-v366
    d/dt(q_syn_T_C1_CTLA4_CD80_CTLA4) <- +v363-v364
    d/dt(q_syn_T_C1_CD80_CTLA4_CD80_CTLA4) <- +v364+v366
    d/dt(q_syn_T_C1_CD86_CTLA4) <- +v367-v368
    d/dt(q_syn_T_C1_CD86_CTLA4_CD86) <- +v368
    d/dt(q_syn_T_C1_PDL1_CD80) <- +v372-v373-v374
    d/dt(q_syn_T_C1_PDL1_CD80_CD28) <- +v373
    d/dt(q_syn_T_C1_PDL1_CD80_CTLA4) <- +v374
    d/dt(q_syn_T_C1_CD28) <- -v359-v360-v361-v373
    d/dt(q_syn_T_C1_CTLA4) <- -v362-v363-v366-v367-v369-v370-v374
    d/dt(q_syn_T_C1_CD80) <- -v359-v362-v364-v365+v371
    d/dt(q_syn_T_C1_CD80m) <- -2*v371-v372
    d/dt(q_syn_T_C1_CD86) <- -v361-v367-v368
    d/dt(q_syn_T_C1_CTLA4_aCTLA4) <- +v369-v370
    d/dt(q_syn_T_C1_CTLA4_aCTLA4_CTLA4) <- +v370
    d/dt(q_syn_T_APC_PD1_PDL1) <- +v419
    d/dt(q_syn_T_APC_PD1_PDL2) <- +v420
    d/dt(q_syn_T_APC_PD1) <- -v419-v420-v421-v422
    d/dt(q_syn_T_APC_PDL1) <- -v375-v377-v381-v383-v387-v393-v398-v399-v403-v411+v415+v417-v419-v436
    d/dt(q_syn_T_APC_PDL2) <- +v416+v418-v420
    d/dt(q_syn_T_APC_PD1_aPD1) <- +v421-v422
    d/dt(q_syn_T_APC_PD1_aPD1_PD1) <- +v422
    d/dt(q_syn_T_APC_TPDL1) <- -v376-v379-v382-v385-v390-v395-v401-v402-v407-v413
    d/dt(q_syn_T_APC_PDL1_aPDL1cRo) <- +v375-v378-v384-v394
    d/dt(q_syn_T_APC_PDL1_aPDL1oRo) <- +v377+v378-v381-v389
    d/dt(q_syn_T_APC_PDL1_aPDL1Rmc) <- +v383+v384-v388-v406
    d/dt(q_syn_T_APC_PDL1_aPDL1mRo) <- +v393+v394-v398-v405
    d/dt(q_syn_T_APC_PDL1_aPDL1Rmo) <- +v387+v388+v389-v399-v404
    d/dt(q_syn_T_APC_PDL1_aPDL1Rmm) <- +v403+v404+v405+v406-v411
    d/dt(q_syn_T_APC_PDL1_aPDL1RoRo_PDL1) <- +v381-v397
    d/dt(q_syn_T_APC_PDL1_aPDL1RmRo_PDL1) <- +v397+v398+v399-v412
    d/dt(q_syn_T_APC_PDL1_aPDL1RmRm_PDL1) <- +v411+v412
    d/dt(q_syn_T_APC_TPDL1_aPDL1cRo) <- +v376-v380-v386-v396
    d/dt(q_syn_T_APC_TPDL1_aPDL1oRo) <- +v379+v380-v382-v392
    d/dt(q_syn_T_APC_TPDL1_aPDL1Rmo) <- +v390+v391+v392-v402-v408
    d/dt(q_syn_T_APC_TPDL1_aPDL1Rmc) <- +v385+v386-v391-v410
    d/dt(q_syn_T_APC_TPDL1_aPDL1mRo) <- +v395+v396-v401-v409
    d/dt(q_syn_T_APC_TPDL1_aPDL1Rmm) <- +v407+v408+v409+v410-v413
    d/dt(q_syn_T_APC_TPDL1_aPDL1RoRo_TPDL1) <- +v382-v400
    d/dt(q_syn_T_APC_TPDL1_aPDL1RmRo_TPDL1) <- +v400+v401+v402-v414
    d/dt(q_syn_T_APC_TPDL1_aPDL1RmRm_TPDL1) <- +v413+v414
    d/dt(q_syn_T_APC_CD28_CD80) <- +v423-v424
    d/dt(q_syn_T_APC_CD28_CD80_CD28) <- +v424
    d/dt(q_syn_T_APC_CD28_CD86) <- +v425
    d/dt(q_syn_T_APC_CD80_CTLA4) <- +v426-v427-v429
    d/dt(q_syn_T_APC_CD80_CTLA4_CD80) <- +v429-v430
    d/dt(q_syn_T_APC_CTLA4_CD80_CTLA4) <- +v427-v428
    d/dt(q_syn_T_APC_CD80_CTLA4_CD80_CTLA4) <- +v428+v430
    d/dt(q_syn_T_APC_CD86_CTLA4) <- +v431-v432
    d/dt(q_syn_T_APC_CD86_CTLA4_CD86) <- +v432
    d/dt(q_syn_T_APC_PDL1_CD80) <- +v436-v437-v438
    d/dt(q_syn_T_APC_PDL1_CD80_CD28) <- +v437
    d/dt(q_syn_T_APC_PDL1_CD80_CTLA4) <- +v438
    d/dt(q_syn_T_APC_CD28) <- -v423-v424-v425-v437
    d/dt(q_syn_T_APC_CTLA4) <- -v426-v427-v430-v431-v433-v434-v438
    d/dt(q_syn_T_APC_CD80) <- -v423-v426-v428-v429+v435
    d/dt(q_syn_T_APC_CD80m) <- -2*v435-v436
    d/dt(q_syn_T_APC_CD86) <- -v425-v431-v432
    d/dt(q_syn_T_APC_CTLA4_aCTLA4) <- +v433-v434
    d/dt(q_syn_T_APC_CTLA4_aCTLA4_CTLA4) <- +v434
    d/dt(q_nonsyn_C1_T_PDL1) <- +v439+v440
    d/dt(q_nonsyn_APC_LN_PDL1) <- +v441+v442
    d/dt(q_syn_M_C_CD47) <- -v488
    d/dt(q_syn_M_C_SIRPa) <- -v488
    d/dt(q_syn_M_C_CD47_SIRPa) <- +v488
    d/dt(q_syn_M_C_PD1_PDL1) <- +v489
    d/dt(q_syn_M_C_PD1_PDL2) <- +v490
    d/dt(q_syn_M_C_PD1) <- -v489-v490-v491-v492
    d/dt(q_syn_M_C_PDL1) <- +v484+v486-v489-v493-v494-v496
    d/dt(q_syn_M_C_PDL2) <- +v485+v487-v490
    d/dt(q_syn_M_C_PD1_aPD1) <- +v491-v492
    d/dt(q_syn_M_C_PD1_aPD1_PD1) <- +v492
    d/dt(q_syn_M_C_PDL1_aPDL1) <- +v493-v494
    d/dt(q_syn_M_C_PDL1_aPDL1_PDL1) <- +v494
    d/dt(q_syn_M_C_PDL1_CD80) <- +v496
    d/dt(q_syn_M_C_CD80) <- +v495
    d/dt(q_syn_M_C_CD80m) <- -2*v495-v496

    # ---------- initial conditions (ESM 2 Species initial values and initialAssignment rules; all other states start at 0) ----------
    # V_T(0) = vol_cell/Ve_T: the V_T rule evaluated at the single initial cancer cell
    q_V_C_nT0(0) <- 4448.27586206897   # ESM 2 Species (V_C.nT0 = 4448.27586206897 cell)
    q_V_C_nT1(0) <- 2729.72972972973   # ESM 2 Species (V_C.nT1 = 2729.72972972973 cell)
    q_V_P_nT0(0) <- 222413.793103448   # ESM 2 Species (V_P.nT0 = 222413.793103448 cell)
    q_V_P_nT1(0) <- 136486.486486486   # ESM 2 Species (V_P.nT1 = 136486.486486486 cell)
    q_V_T_C1(0) <- 1   # ESM 2 Species (V_T.C1 = 1 cell)
    q_V_T_P0(0) <- 1e-09*(vol_cell/Ve_T)   # ESM 2 Species (V_T.P0 = 1e-18 molarity)
    q_V_T_P1(0) <- 1e-09*(vol_cell/Ve_T)   # ESM 2 Species (V_T.P1 = 1e-18 molarity)
    q_V_LN_nT0(0) <- 186.827586206897   # ESM 2 Species (V_LN.nT0 = 186.827586206897 cell)
    q_V_LN_IL2(0) <- 0.00019*vol_V_LN   # ESM 2 Species (V_LN.IL2 = 0.00019 nanomolarity)
    q_V_LN_nT1(0) <- 136.486486486486   # ESM 2 Species (V_LN.nT1 = 136.486486486486 cell)
    q_V_e_P0(0) <- 1e-09*vol_V_e   # ESM 2 Species (V_e.P0 = 1e-18 molarity)
    q_V_e_p0(0) <- 1e-09*vol_V_e   # ESM 2 Species (V_e.p0 = 1e-18 molarity)
    q_V_e_P1(0) <- 1e-09*vol_V_e   # ESM 2 Species (V_e.P1 = 1e-18 molarity)
    q_V_e_p1(0) <- 1e-09*vol_V_e   # ESM 2 Species (V_e.p1 = 1e-18 molarity)
    q_A_e_M1(0) <- 21857923497267.8*vol_A_e   # ESM 2 Species (A_e.M1 = 2185.79234972678 molecule/micrometer^2)
    q_A_e_M1p0(0) <- 10000*vol_A_e   # ESM 2 Species (A_e.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_e_M1p1(0) <- 10000*vol_A_e   # ESM 2 Species (A_e.M1p1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1(0) <- 10000*vol_A_s   # ESM 2 Species (A_s.M1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p0(0) <- 10000*vol_A_s   # ESM 2 Species (A_s.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p1(0) <- 10000*vol_A_s   # ESM 2 Species (A_s.M1p1 = 1e-06 molecule/micrometer^2)
    q_syn_T_C1_PD1(0) <- (T_PD1_total/A_Tcell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.PD1 = T_PD1_total/A_Tcell)
    q_syn_T_C1_PDL1(0) <- (C1_PDL1_base/A_cell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.PDL1 = C1_PDL1_base/A_cell)
    q_syn_T_C1_PDL2(0) <- (C1_PDL1_base * r_PDL2C1/A_cell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.PDL2 = C1_PDL1_base*r_PDL2C1/A_cell)
    q_syn_T_C1_TPDL1(0) <- (T_PDL1_total/A_Tcell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.TPDL1 = T_PDL1_total/A_Tcell)
    q_syn_T_C1_CD28(0) <- (T_CD28_total/A_Tcell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.CD28 = T_CD28_total/A_Tcell)
    q_syn_T_C1_CTLA4(0) <- (T_CTLA4_syn/A_Tcell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.CTLA4 = T_CTLA4_syn/A_Tcell)
    q_syn_T_C1_CD80(0) <- (C1_CD80_total/A_cell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.CD80 = C1_CD80_total/A_cell)
    q_syn_T_C1_CD86(0) <- (C1_CD86_total/A_cell)*vol_syn_T_C1   # ESM 2 initialAssignment (syn_T_C1.CD86 = C1_CD86_total/A_cell)
    q_syn_T_APC_PD1(0) <- (T_PD1_total/A_Tcell)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.PD1 = T_PD1_total/A_Tcell)
    q_syn_T_APC_PDL1(0) <- (APC_PDL1_base/A_APC)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.PDL1 = APC_PDL1_base/A_APC)
    q_syn_T_APC_PDL2(0) <- (APC_PDL1_base * r_PDL2APC/A_APC)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.PDL2 = APC_PDL1_base*r_PDL2APC/A_APC)
    q_syn_T_APC_TPDL1(0) <- (T_PDL1_total/A_Tcell)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.TPDL1 = T_PDL1_total/A_Tcell)
    q_syn_T_APC_CD28(0) <- (T_CD28_total/A_Tcell)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.CD28 = T_CD28_total/A_Tcell)
    q_syn_T_APC_CTLA4(0) <- (T_CTLA4_syn/A_Tcell)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.CTLA4 = T_CTLA4_syn/A_Tcell)
    q_syn_T_APC_CD80(0) <- (APC_CD80_total/A_APC)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.CD80 = APC_CD80_total/A_APC)
    q_syn_T_APC_CD86(0) <- (APC_CD86_total/A_APC)*vol_syn_T_APC   # ESM 2 initialAssignment (syn_T_APC.CD86 = APC_CD86_total/A_APC)
    q_nonsyn_C1_T_PDL1(0) <- (C1_PDL1_base * Lig_fac/A_cell)*vol_nonsyn_C1_T   # ESM 2 initialAssignment (nonsyn_C1_T.PDL1 = C1_PDL1_base*Lig_fac/A_cell)
    q_nonsyn_APC_LN_PDL1(0) <- (APC_PDL1_base * Lig_fac/A_APC)*vol_nonsyn_APC_LN   # ESM 2 initialAssignment (nonsyn_APC_LN.PDL1 = APC_PDL1_base*Lig_fac/A_APC)
    q_syn_M_C_CD47(0) <- (C_CD47)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.CD47 = C_CD47)
    q_syn_M_C_SIRPa(0) <- (M_SIRPa)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.SIRPa = M_SIRPa)
    q_syn_M_C_PD1(0) <- (M_PD1_total/A_Mcell)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.PD1 = M_PD1_total/A_Mcell)
    q_syn_M_C_PDL1(0) <- (C1_PDL1_base/A_cell)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.PDL1 = C1_PDL1_base/A_cell)
    q_syn_M_C_PDL2(0) <- (C1_PDL1_base * r_PDL2C1/A_cell)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.PDL2 = C1_PDL1_base*r_PDL2C1/A_cell)
    q_syn_M_C_CD80(0) <- (C1_CD80_total/A_cell)*vol_syn_M_C   # ESM 2 initialAssignment (syn_M_C.CD80 = C1_CD80_total/A_cell)

    # ---------- outputs ----------
    # tumour diameter (cm) of the spherical tumour of volume V_T (L = 1000 cm^3), the RECIST 1.1 size
    tumour_diameter <- 2*(3*1000*vol_V_T/(4*pi))^(1/3)
    tumour_diameter_ratio <- tumour_diameter/initial_tumour_diameter
    # total free antibody in the central compartment (nmol/L), all six mask/cleavage states
    Cc <- x_V_C_aPDL1cc + x_V_C_aPDL1co + x_V_C_aPDL1oo + x_V_C_aPDL1c + x_V_C_aPDL1o + x_V_C_aPDL1
  })
}
