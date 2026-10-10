RuizMartinez_2022_nivolumab_qsp <- function() {
  description <- paste(
    "QSP. Whole-patient quantitative systems pharmacology model of",
    "triple-negative breast cancer (TNBC) under anti-PD-1 immunotherapy",
    "(nivolumab), the non-spatial QSP layer of the Ruiz-Martinez 2022 spatial",
    "QSP (spQSP) platform. Four physiological compartments (central,",
    "peripheral, tumour, tumour-draining lymph node) plus APC endosomal and",
    "surface compartments and two immunological-synapse compartments. The",
    "equations are those of Wang 2020 (HER2-negative breast cancer; packaged as",
    "Wang_2020_entinostat_nivolumab_ipilimumab_qsp), unchanged, with the",
    "parameter set recalibrated by Ruiz-Martinez et al. against the number of",
    "anti-PD-1 responders and the tumour T cell densities of TNBC: logistic",
    "cancer-cell growth from a single cell, naive regulatory and cytotoxic T",
    "cell priming in the lymph node, T cell trafficking, exhaustion and Treg",
    "suppression, APC maturation and antigen presentation, PD-1/PD-L1/PD-L2 and",
    "CTLA-4/CD28/CD80/CD86 synapse binding with cross-arm antibody binding, an",
    "MDSC module, and the nivolumab, ipilimumab, durvalumab and entinostat PK",
    "modules carried by the deposited model (only nivolumab is dosed in the",
    "paper). 120 ODE states, 21 repeated-assignment rules and 188 parameters,",
    "with every value taken from the paper's S1/S2 Tables. The agent-based",
    "(ABM) half of the spQSP model is stochastic and spatial and is not part of",
    "this file. Deterministic mechanism model: no etas and no error model.",
    sep = " "
  )
  reference <- paste(
    "Ruiz-Martinez A, Gong C, Wang H, Sove RJ, Mi H, Kimko H, Popel AS.",
    "Simulations of tumor growth and response to immunotherapy by coupling a",
    "spatial agent-based model with a whole-patient quantitative systems",
    "pharmacology model. PLoS Comput Biol. 2022;18(7):e1010254.",
    "doi:10.1371/journal.pcbi.1010254.",
    sep = " "
  )
  vignette <- "RuizMartinez_2022_nivolumab_qsp"

  # The 120 ODE states keep the SimBiology compartment.species naming
  # (q_<compartment>_<species>) so each traces to its S1 Table species row;
  # they are QSP cell / molecule pools, not PK compartments.
  paper_specific_compartments <- c(
    "q_V_C_T0",
    "q_V_C_T1",
    "q_V_C_nivo",
    "q_V_C_durv",
    "q_V_C_ipi",
    "q_V_C_ENT",
    "q_V_C_ENT_Buccal",
    "q_V_C_ENT_GI",
    "q_V_C_Dose2",
    "q_V_P_T0",
    "q_V_P_T1",
    "q_V_P_nivo",
    "q_V_P_durv",
    "q_V_P_ipi",
    "q_V_P_Treg_CTLA4",
    "q_V_P_Treg_CTLA4_ipi",
    "q_V_P_Treg_CTLA4_ipi_CTLA4",
    "q_V_P_ENT",
    "q_V_T_C_x",
    "q_V_T_T_exh",
    "q_V_T_C1",
    "q_V_T_T0",
    "q_V_T_T1",
    "q_V_T_APC",
    "q_V_T_mAPC",
    "q_V_T_c",
    "q_V_T_nivo",
    "q_V_T_durv",
    "q_V_T_ipi",
    "q_V_T_Treg_CTLA4",
    "q_V_T_Treg_CTLA4_ipi",
    "q_V_T_Treg_CTLA4_ipi_CTLA4",
    "q_V_T_MDSC",
    "q_V_T_CCL2",
    "q_V_T_NO",
    "q_V_T_ArgI",
    "q_V_T_ENT",
    "q_V_LN_nT0",
    "q_V_LN_aT0",
    "q_V_LN_T0",
    "q_V_LN_IL2",
    "q_V_LN_nT1",
    "q_V_LN_aT1",
    "q_V_LN_T1",
    "q_V_LN_APC",
    "q_V_LN_mAPC",
    "q_V_LN_P0",
    "q_V_LN_P1",
    "q_V_LN_nivo",
    "q_V_LN_durv",
    "q_V_LN_ipi",
    "q_V_LN_ENT",
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
    "q_syn_T_C1_PD1_nivo",
    "q_syn_T_C1_PD1_nivo_PD1",
    "q_syn_T_C1_PDL1_durv",
    "q_syn_T_C1_PDL1_durv_PDL1",
    "q_syn_T_C1_TPDL1",
    "q_syn_T_C1_TPDL1_durv",
    "q_syn_T_C1_TPDL1_durv_TPDL1",
    "q_syn_T_C1_CD28_CD80",
    "q_syn_T_C1_CD28_CD80_CD28",
    "q_syn_T_C1_CD28_CD86",
    "q_syn_T_C1_CD80_CTLA4",
    "q_syn_T_C1_CD80_CTLA4_CD80",
    "q_syn_T_C1_CTLA4_CD80_CTLA4",
    "q_syn_T_C1_CD80_CTLA4_CD80_CTLA4",
    "q_syn_T_C1_CD86_CTLA4",
    "q_syn_T_C1_CD86_CTLA4_CD86",
    "q_syn_T_C1_TPDL1_CD80",
    "q_syn_T_C1_TPDL1_CD80_TPDL1",
    "q_syn_T_C1_CD28",
    "q_syn_T_C1_CTLA4",
    "q_syn_T_C1_CD80",
    "q_syn_T_C1_CD86",
    "q_syn_T_C1_CTLA4_ipi",
    "q_syn_T_C1_CTLA4_ipi_CTLA4",
    "q_syn_T_APC_PD1_PDL1",
    "q_syn_T_APC_PD1_PDL2",
    "q_syn_T_APC_PD1",
    "q_syn_T_APC_PDL1",
    "q_syn_T_APC_PDL2",
    "q_syn_T_APC_PD1_nivo",
    "q_syn_T_APC_PD1_nivo_PD1",
    "q_syn_T_APC_PDL1_durv",
    "q_syn_T_APC_PDL1_durv_PDL1",
    "q_syn_T_APC_TPDL1",
    "q_syn_T_APC_TPDL1_durv",
    "q_syn_T_APC_TPDL1_durv_TPDL1",
    "q_syn_T_APC_CD28_CD80",
    "q_syn_T_APC_CD28_CD80_CD28",
    "q_syn_T_APC_CD28_CD86",
    "q_syn_T_APC_CD80_CTLA4",
    "q_syn_T_APC_CD80_CTLA4_CD80",
    "q_syn_T_APC_CTLA4_CD80_CTLA4",
    "q_syn_T_APC_CD80_CTLA4_CD80_CTLA4",
    "q_syn_T_APC_CD86_CTLA4",
    "q_syn_T_APC_CD86_CTLA4_CD86",
    "q_syn_T_APC_TPDL1_CD80",
    "q_syn_T_APC_TPDL1_CD80_TPDL1",
    "q_syn_T_APC_CD28",
    "q_syn_T_APC_CTLA4",
    "q_syn_T_APC_CD80",
    "q_syn_T_APC_CD86",
    "q_syn_T_APC_CTLA4_ipi",
    "q_syn_T_APC_CTLA4_ipi_CTLA4"
  )

  units <- list(
    time = "day",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  compartmentData <- list(
    q_V_C_T0 = list(
      analyte = "Number of T0 cells in the central compartment",
      units = "cell",
      specimen = "whole blood",
      verified = FALSE
    ),
    q_V_C_T1 = list(
      analyte = "Number of T1 cells in the central compartment",
      units = "cell",
      specimen = "whole blood",
      verified = FALSE
    ),
    q_V_C_nivo = list(
      analyte = "Concentration of nivo in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_C_durv = list(
      analyte = "Concentration of durv in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_C_ipi = list(
      analyte = "Concentration of ipi in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_C_ENT = list(
      analyte = "Concentration of ENT in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_C_ENT_Buccal = list(
      analyte = "Concentration of ENT in central compartment",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    q_V_C_ENT_GI = list(
      analyte = "Concentration of ENT in central compartment",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    q_V_C_Dose2 = list(
      analyte = "Concentration of ENT before entering depot 2",
      units = "nmol",
      specimen = "administration site",
      verified = FALSE
    ),
    q_V_P_T0 = list(
      analyte = "Number of T0 cells in the peripheral compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_T1 = list(
      analyte = "Number of T1 cells in the peripheral compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_nivo = list(
      analyte = "Concentration of nivo in peripheral compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_durv = list(
      analyte = "Concentration of durv in peripheral compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_ipi = list(
      analyte = "Concentration of ipi in peripheral compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_Treg_CTLA4 = list(
      analyte = "Number of free CTLA4 molecules on Treg in Peripheral compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_Treg_CTLA4_ipi = list(
      analyte = "Number of CTLA4-ipilimumab complex on Treg in Peripheral compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_Treg_CTLA4_ipi_CTLA4 = list(
      analyte = "Number of CTLA4-ipilimumab-CTLA4 complex on Treg in Peripheral compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_ENT = list(
      analyte = "Concentration of ENT in peripheral compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_T_C_x = list(
      analyte = "Dead cancer cells in the tumour compartment",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_T_exh = list(analyte = "Exhausted T cells", units = "cell", specimen = "tumor", verified = FALSE),
    q_V_T_C1 = list(analyte = "Number of cancer cells in tumour", units = "cell", specimen = "tumor", verified = FALSE),
    q_V_T_T0 = list(
      analyte = "Number of T0 cells in the tumour compartment",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_T1 = list(
      analyte = "Number of T1 cells in the tumour compartment",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_APC = list(
      analyte = "Number of naive antigen presenting cells in the tumour",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_mAPC = list(
      analyte = "Number of mature antigen presenting cells in the tumour",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_c = list(
      analyte = "Concentration of maturation cytokines in the tumour",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_nivo = list(
      analyte = "Concentration of nivo in tumour compartment",
      units = "nmol",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_durv = list(
      analyte = "Concentration of durv in tumour compartment",
      units = "nmol",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_ipi = list(
      analyte = "Concentration of ipi in tumour compartment",
      units = "nmol",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_Treg_CTLA4 = list(
      analyte = "Number of free CTLA4 molecules on Treg in Tumour compartment",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_Treg_CTLA4_ipi = list(
      analyte = "Number of CTLA4-ipilimumab complex on Treg in Tumour compartment",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_Treg_CTLA4_ipi_CTLA4 = list(
      analyte = "Number of CTLA4-ipilimumab-CTLA4 complex on Treg in Tumour compartment",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_MDSC = list(
      analyte = "Number of MDSCs in the tumour compartment",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_CCL2 = list(
      analyte = "Concentration of CCL2 in the tumor compartment",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_NO = list(
      analyte = "Concentration of NO in the tumor compartment",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_ArgI = list(
      analyte = "Concentration of Arg I in the tumor compartment",
      units = "molecule (arginase-I activity placeholder, 1 mU = 1 mol/L in the export)",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_ENT = list(
      analyte = "Concentration of ENT in tumour compartment",
      units = "nmol",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_LN_nT0 = list(
      analyte = "Number of naive T0 cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_aT0 = list(
      analyte = "Number of activated T0 cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_T0 = list(
      analyte = "Number of T0 cells in the lymph node compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_IL2 = list(
      analyte = "Concentration of IL2 in the lymph node compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_nT1 = list(
      analyte = "Number of naive T1 cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_aT1 = list(
      analyte = "Number of activated T1 cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_T1 = list(
      analyte = "Number of T1 cells in the lymph node compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_APC = list(
      analyte = "Number of naive antigen presenting cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_mAPC = list(
      analyte = "Number of mature antigen presenting cells in the lymph node",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_P0 = list(
      analyte = "Concentration of free antigen (P0) in the LN compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_P1 = list(
      analyte = "Concentration of free antigen (P1) in the LN compartment",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_nivo = list(
      analyte = "Concentration of nivo in lymph node compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_durv = list(
      analyte = "Concentration of durv in lymph node compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_ipi = list(
      analyte = "Concentration of ipi in lymph node compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_ENT = list(
      analyte = "Concentration of ENT in lymph node compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_e_P0 = list(
      analyte = "Concentration of antigen P0 in the APC endosomes",
      units = "molecule",
      specimen = "endosome",
      verified = FALSE
    ),
    q_V_e_p0 = list(
      analyte = "Concentration of epitope P0 in the APC endosomes",
      units = "molecule",
      specimen = "endosome",
      verified = FALSE
    ),
    q_V_e_P1 = list(
      analyte = "Concentration of antigen P1 in the APC endosomes",
      units = "molecule",
      specimen = "endosome",
      verified = FALSE
    ),
    q_V_e_p1 = list(
      analyte = "Concentration of epitope P1 in the APC endosomes",
      units = "molecule",
      specimen = "endosome",
      verified = FALSE
    ),
    q_A_e_M1 = list(
      analyte = "Amount of MHC per area on the cell surface",
      units = "molecule",
      specimen = "endosome",
      verified = FALSE
    ),
    q_A_e_M1p0 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "endosome", verified = FALSE),
    q_A_e_M1p1 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "endosome", verified = FALSE),
    q_A_s_M1 = list(
      analyte = "Amount of MHC per area on the cell surface",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_A_s_M1p0 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "tissue", verified = FALSE),
    q_A_s_M1p1 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "tissue", verified = FALSE),
    q_syn_T_C1_PD1_PDL1 = list(
      analyte = "concentration of PD1-PDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PD1_PDL2 = list(
      analyte = "concentration of PD1-PDL2 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PD1 = list(
      analyte = "concentration of PD1 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PDL1 = list(
      analyte = "concentration of PDL1 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PDL2 = list(
      analyte = "concentration of PDL2 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PD1_nivo = list(
      analyte = "concentration of PD1-nivolumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PD1_nivo_PD1 = list(
      analyte = "concentration of PD1-nivolumab-PD1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PDL1_durv = list(
      analyte = "concentration of PDL1-durvalumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PDL1_durv_PDL1 = list(
      analyte = "concentration of PDL1-durvalumab-PDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_TPDL1 = list(
      analyte = "concentration of PDL1 in synapse of T cell",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_TPDL1_durv = list(
      analyte = "concentration of TPDL1-durvalumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_TPDL1_durv_TPDL1 = list(
      analyte = "concentration of TPDL1-durvalumab-TPDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD28_CD80 = list(
      analyte = "concentration of CD28-CD80 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD28_CD80_CD28 = list(
      analyte = "concentration of CD28-CD80-CD28 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD28_CD86 = list(
      analyte = "concentration of CD28-CD86 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD80_CTLA4 = list(
      analyte = "concentration of CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD80_CTLA4_CD80 = list(
      analyte = "concentration of CD80-CTLA4-CD80 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CTLA4_CD80_CTLA4 = list(
      analyte = "concentration of CTLA4-CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD80_CTLA4_CD80_CTLA4 = list(
      analyte = "concentration of CD80-CTLA4-CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD86_CTLA4 = list(
      analyte = "concentration of CD86-CTLA4 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD86_CTLA4_CD86 = list(
      analyte = "concentration of CD86-CTLA4-CD86 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_TPDL1_CD80 = list(
      analyte = "concentration of TPDL1-CD80 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_TPDL1_CD80_TPDL1 = list(
      analyte = "concentration of TPDL1-CD80-TPDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD28 = list(
      analyte = "concentration of CD28 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CTLA4 = list(
      analyte = "concentration of CTLA4 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD80 = list(
      analyte = "concentration of CD80 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CD86 = list(
      analyte = "concentration of CD86 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CTLA4_ipi = list(
      analyte = "concentration of CTLA4-ipilimumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_CTLA4_ipi_CTLA4 = list(
      analyte = "concentration of CTLA4-ipilimumab-CTLA4 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_APC_PD1_PDL1 = list(
      analyte = "concentration of PD1-PDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PD1_PDL2 = list(
      analyte = "concentration of PD1-PDL2 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PD1 = list(
      analyte = "concentration of PD1 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PDL1 = list(
      analyte = "concentration of PDL1 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PDL2 = list(
      analyte = "concentration of PDL2 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PD1_nivo = list(
      analyte = "concentration of PD1-nivolumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PD1_nivo_PD1 = list(
      analyte = "concentration of PD1-nivolumab-PD1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PDL1_durv = list(
      analyte = "concentration of PDL1-durvalumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PDL1_durv_PDL1 = list(
      analyte = "concentration of PDL1-durvalumab-PDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_TPDL1 = list(
      analyte = "concentration of PDL1 in synapse of T cell",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_TPDL1_durv = list(
      analyte = "concentration of TPDL1-durvalumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_TPDL1_durv_TPDL1 = list(
      analyte = "concentration of TPDL1-durvalumab-TPDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD28_CD80 = list(
      analyte = "concentration of CD28-CD80 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD28_CD80_CD28 = list(
      analyte = "concentration of CD28-CD80-CD28 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD28_CD86 = list(
      analyte = "concentration of CD28-CD86 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD80_CTLA4 = list(
      analyte = "concentration of CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD80_CTLA4_CD80 = list(
      analyte = "concentration of CD80-CTLA4-CD80 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CTLA4_CD80_CTLA4 = list(
      analyte = "concentration of CTLA4-CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD80_CTLA4_CD80_CTLA4 = list(
      analyte = "concentration of CD80-CTLA4-CD80-CTLA4 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD86_CTLA4 = list(
      analyte = "concentration of CD86-CTLA4 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD86_CTLA4_CD86 = list(
      analyte = "concentration of CD86-CTLA4-CD86 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_TPDL1_CD80 = list(
      analyte = "concentration of TPDL1-CD80 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_TPDL1_CD80_TPDL1 = list(
      analyte = "concentration of TPDL1-CD80-TPDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD28 = list(
      analyte = "concentration of CD28 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CTLA4 = list(
      analyte = "concentration of CTLA4 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD80 = list(
      analyte = "concentration of CD80 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CD86 = list(
      analyte = "concentration of CD86 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CTLA4_ipi = list(
      analyte = "concentration of CTLA4-ipilimumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_CTLA4_ipi_CTLA4 = list(
      analyte = "concentration of CTLA4-ipilimumab-CTLA4 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    )
  )

  population <- list(
    species = "human (in silico virtual patients)",
    n_subjects = 100L,
    disease_state = "Triple-negative breast cancer (TNBC)",
    dose_range = paste(
      "Nivolumab 3 mg/kg every 2 weeks (80 kg patient, 144.6 kDa), started",
      "once the tumour reaches its pre-treatment diameter (2.255 cm for the",
      "default parameter set).",
      sep = " "
    ),
    notes = paste(
      "Methodological in silico study; no clinical data were fitted. The QSP",
      "parameters were recalibrated so that 100 virtual patients sampled from",
      "the ranges of S1 Table agree with the anti-PD-1 responder numbers and",
      "tumour T cell densities of the entinostat + nivolumab TNBC trial",
      "ETCTN-9844 (Torres 2021); S1 and S2 Figs show those 100 cases without",
      "and with nivolumab. The packaged values are the default QSP parameter",
      "set of S2 Table, which the spQSP simulations of the paper use (Figures",
      "4-8 override k_C1_growth, k_C_T1, k_T1 and Treg_max per growth-rate",
      "case).",
      sep = " "
    )
  )

  ini({
    # ---- compartment capacities (S1 Table, sheet 'Table S1. Compartment') ----
    vol_V_C <- fixed(6); label("Capacity of compartment V_C [L]")  # S1 Table, sheet 'Compartment' (V_C = 6 litre)
    vol_V_P <- fixed(61.321); label("Capacity of compartment V_P [L]")  # S1 Table, sheet 'Compartment' (V_P = 61.321 litre)
    vol_V_LN <- fixed(0.00111264739815); label("Capacity of compartment V_LN [L]")  # S1 Table, sheet 'Compartment' (V_LN = 1112.647398146385 millimeter^3)
    vol_V_e <- fixed(4e-16); label("Capacity of compartment V_e [L]")  # S1 Table, sheet 'Compartment' (V_e = 4.0000000000000004e-16 litre)
    vol_A_e <- fixed(1.5e-09); label("Capacity of compartment A_e [dm^2]")  # S1 Table, sheet 'Compartment' (A_e = 15 micrometer^2)
    vol_A_s <- fixed(9e-08); label("Capacity of compartment A_s [dm^2]")  # S1 Table, sheet 'Compartment' (A_s = 900 micrometer^2)
    vol_syn_T_C1 <- fixed(3.78e-09); label("Capacity of compartment syn_T_C1 [dm^2]")  # S1 Table, sheet 'Compartment' (syn_T_C1 = 37.8 micrometer^2)
    vol_syn_T_APC <- fixed(3.78e-09); label("Capacity of compartment syn_T_APC [dm^2]")  # S1 Table, sheet 'Compartment' (syn_T_APC = 37.8 micrometer^2)

    # ---- model parameters (S2 Table, sheet 'QSP default parameters'; S1 Table,
    # sheet 'Table S2.Parameters' for sampling ranges), converted to the model
    # unit system ----
    k_cell_clear <- fixed(0.1); label("Rate of dead cell clearance from tumour compartment [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_cell_clear = 0.1 1/day)
    vol_cell <- fixed(2.57244078451e-12); label("Average volume of cancer cell calculated based on cancer cell diameter [L/count]")  # S2 Table, sheet 'QSP default parameters' (vol_cell = 2572.4407845144424 micrometer^3/cell)
    vol_Tcell <- fixed(1.750157098e-13); label("Average volume of T cells calculated based on the average T cell diameter [L/count]")  # S2 Table, sheet 'QSP default parameters' (vol_Tcell = 175.01570979953922 micrometer^3/cell)
    V_Tmin <- fixed(1e-09); label("Cancer-Free Tumour compartment volume [L]")  # S2 Table, sheet 'QSP default parameters' (V_Tmin = 1e-06 milliliter)
    k_C1_growth <- fixed(0.005); label("Cancer cell growth rate. Estimated by tumor doubling time: 241?+/-?166 days for ER-positive BC, and ... [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_C1_growth = 0.005 1/day); S1 Table sampling range 'log normal [0.0087, 0.5]'
    C_max <- fixed(862890782186); label("Cancer cell capacity calculated from tumour maximum diameter and cancer cell density [count]")  # S2 Table, sheet 'QSP default parameters' (C_max = 862890782185.9965 cell)
    k_C1_death <- fixed(3e-05); label("Cancer cell death rate from innate immune cells [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_C1_death = 3e-05 1/day); S1 Table sampling range 'log unif [1e-5, 1e-3]'
    k_C1_therapy <- fixed(0); label("Rate of C1 killing by therapy [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_C1_therapy = 0 1/day)
    n_T0_clones <- fixed(100); label("Number of T cell clones [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_T0_clones = 100 dimensionless)
    Q_T0_in <- fixed(346.455379492); label("Rate of naive T cell transport into the lLN calculated based on naive T cell entry rate estimated based on ... [count/day]")  # S2 Table, sheet 'QSP default parameters' (Q_T0_in = 346.45537949178816 cell/day)
    Q_T0_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (Q_T0_out = 1.13 1/day)
    k_T0_act <- fixed(20); label("T0 activation rate calculated based on the rate of T cell activation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T0_act = 20 1/day)
    k_T0_pro <- fixed(1); label("T0 proliferation rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T0_pro = 1 1/day)
    k_T0_death <- fixed(0.01); label("T0 death rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T0_death = 0.01 1/day)
    q_T0_P_in <- fixed(295.812504); label("rate of T0 tranport into the peripheral compartment calculated based on T cell transmigration rate and T cell ... [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T0_P_in = 0.20542535 1/minute)
    q_T0_P_out <- fixed(24); label("rate of T0 tranport out of the peripheral compartment [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T0_P_out = 24 1/day)
    q_T0_T_in <- fixed(14.4); label("rate of T0 tranport into the tumour compartment calculated based on T cell transmigration rate and T cell ... [1/(L*day)]")  # S2 Table, sheet 'QSP default parameters' (q_T0_T_in = 1e-05 1/(centimeter^3*minute))
    q_T0_LN_out <- fixed(24); label("rate of T0 tranport out of the LN compartment [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T0_LN_out = 24 1/day)
    k_IL2_deg <- fixed(288); label("rate of IL2 degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_IL2_deg = 0.2 1/minute)
    k_IL2_cons <- fixed(86718844656); label("rate of IL2 consumption by T cells [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_IL2_cons = 6e-06 nanomole/cell/hour)
    k_IL2_sec <- fixed(433594223280); label("rate of IL2 secretion from T cells [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_IL2_sec = 3e-05 nanomole/cell/hour)
    IL2_50 <- fixed(1.9270854368e+14); label("T cell activation half-maximal IL2 concentration [count/L]")  # S2 Table, sheet 'QSP default parameters' (IL2_50 = 0.32 nanomolarity)
    IL2_50_Treg <- fixed(1.9270854368e+12); label("Treg activation half-maximal IL2 concentration [count/L]")  # S2 Table, sheet 'QSP default parameters' (IL2_50_Treg = 0.0032 nanomolarity)
    N0 <- fixed(5); label("numer of activated T cell generation by TCR signaling only [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (N0 = 5 dimensionless)
    N_costim <- fixed(6); label("numer of activated T cell generation by co-stimulatory signaling only [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (N_costim = 6 dimensionless)
    N_IL2 <- fixed(14); label("maximum number of activated T cell generations due to IL2 [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (N_IL2 = 14 dimensionless)
    k_Treg <- fixed(0.244); label("Rate of T cell death by Tregs [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_Treg = 0.244 1/day); S1 Table sampling range 'log unif [0.1, 0.5]'
    n_T1_clones <- fixed(61.36); label("Number of T cell clones [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_T1_clones = 61.36 dimensionless); S1 Table sampling range 'log normal [63 - 0.5]'
    Q_T1_in <- fixed(212.6058677); label("Rate of naive T cell transport into the lLN calculated based on naive T cell entry rate estimated based on ... [count/day]")  # S2 Table, sheet 'QSP default parameters' (Q_T1_in = 212.605867699863 cell/day)
    Q_T1_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (Q_T1_out = 1.13 1/day)
    k_T1_act <- fixed(20); label("T1 activation rate calculated based on the rate of T cell activation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T1_act = 20 1/day)
    k_T1_pro <- fixed(1); label("T1 proliferation rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T1_pro = 1 1/day)
    k_T1_death <- fixed(0.01); label("T1 death rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T1_death = 0.01 1/day)
    q_T1_P_in <- fixed(295.812504); label("rate of T1 tranport into the peripheral compartment calculated based on T cell transmigration rate and T cell ... [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T1_P_in = 0.20542535 1/minute)
    q_T1_P_out <- fixed(24); label("rate of T1 tranport out of the peripheral compartment [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T1_P_out = 24 1/day)
    q_T1_T_in <- fixed(14.4); label("rate of T1 tranport into the tumour compartment calculated based on T cell transmigration rate and T cell ... [1/(L*day)]")  # S2 Table, sheet 'QSP default parameters' (q_T1_T_in = 1e-05 1/(centimeter^3*minute))
    q_T1_LN_out <- fixed(24); label("rate of T1 tranport out of the LN compartment [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_T1_LN_out = 24 1/day)
    k_T1 <- fixed(0.134); label("Rate of T cell exhaustion by cancer cells [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_T1 = 0.134 1/day); S1 Table sampling range 'log normal [0.5, 0.5]'
    k_C_T1 <- fixed(0.736); label("Rate of cancer cell death by T cells [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_C_T1 = 0.736 1/day); S1 Table sampling range 'log normal [0.5, 0.5]'
    k_APC_mat <- fixed(1.5); label("Maximum rate of APC maturation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_APC_mat = 1.5 1/day)
    k_APC_mig <- fixed(4); label("Rate of APC migration [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_APC_mig = 4 1/day)
    k_APC_death <- fixed(0.01); label("Rate of APC death [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_APC_death = 0.01 1/day)
    k_mAPC_death <- fixed(0.02); label("Rate of mAPC death [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_mAPC_death = 0.02 1/day)
    APC0_T <- fixed(400000000); label("Steady state density of CD103+ APCs in the tumor, calculated to fit basal APC equal to 26000/cm^3 from [count/L]")  # S2 Table, sheet 'QSP default parameters' (APC0_T = 400000 cell/milliliter)
    APC0_LN <- fixed(1200000000); label("Steady state density of all APC subtypes in the LN- divided by 10 to only take CD103+ based on [count/L]")  # S2 Table, sheet 'QSP default parameters' (APC0_LN = 1200000 cell/milliliter)
    k_c <- fixed(2); label("Cytokine rate constant [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_c = 2 1/day)
    c0 <- fixed(6.02214199e+14); label("Baseline cytokine concentration [count/L]")  # S2 Table, sheet 'QSP default parameters' (c0 = 1e-09 molarity)
    c50 <- fixed(6.02214199e+14); label("Cytokine concentration for half-maximal APC maturation [count/L]")  # S2 Table, sheet 'QSP default parameters' (c50 = 1e-09 molarity)
    DAMPs <- fixed(8069670266.6); label("Concentration of cytokines released per dying cancer cell [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (DAMPs = 1.34e-14 mole/cell)
    n_sites_APC <- fixed(10); label("Maxium number of T Cells an APC can interact with [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_sites_APC = 10 dimensionless)
    kin <- fixed(14.4); label("Rate of MHC internalization [1/day]")  # S2 Table, sheet 'QSP default parameters' (kin = 14.4 1/day)
    kout <- fixed(28.8); label("Rate of MHC externalization [1/day]")  # S2 Table, sheet 'QSP default parameters' (kout = 28.8 1/day)
    k_P0_up <- fixed(14.4); label("Rate of antigen uptake by APCs [1/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_P0_up = 14.4 1/day/cell)
    k_xP0_deg <- fixed(2); label("Rate of extracellular antigen degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_xP0_deg = 2 1/day)
    k_P0_deg <- fixed(17.28); label("Rate of endosomal antigen degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_P0_deg = 17.28 1/day)
    k_p0_deg <- fixed(144); label("Rate of endosomal epitope degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_p0_deg = 144 1/day)
    k_P0_on <- fixed(2.39117576834e-19); label("Rate of antigen-MHC binding [L/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_P0_on = 144000 1/day/molarity)
    k_P0_d1 <- fixed(6.02214199e+16); label("Antigen-MHC kd [count/L]")  # S2 Table, sheet 'QSP default parameters' (k_P0_d1 = 1e-07 molarity)
    p0_50 <- fixed(264550.26455); label("calculated based on the number of molecules for half-maximal T cell activation and synapse surface area [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (p0_50 = 2.6455026455026456e-05 molecule/micrometer^2)
    P0_C1 <- fixed(6.02214199e+15); label("Concentration of P0 in C1 [1/L]")  # S2 Table, sheet 'QSP default parameters' (P0_C1 = 1e-08 molarity/cell)
    A_syn <- fixed(3.78e-09); label("Surface area of the synapse [dm^2]")  # S2 Table, sheet 'QSP default parameters' (A_syn = 37.8 micrometer^2)
    A_Tcell <- fixed(1.5131041193e-08); label("Surface area of the T cell calculated based on the average T cell diameter [dm^2]")  # S2 Table, sheet 'QSP default parameters' (A_Tcell = 151.31041193043737 micrometer^2)
    A_cell <- fixed(9.07920276887e-08); label("Surface area of the Cancer cell calculated based on the average Cancer cell diameter [dm^2]")  # S2 Table, sheet 'QSP default parameters' (A_cell = 907.9202768874502 micrometer^2)
    A_APC <- fixed(9e-08); label("Surface area of the APC [dm^2]")  # S2 Table, sheet 'QSP default parameters' (A_APC = 900 micrometer^2)
    k_M1p0_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_M1p0_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p0_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_M1p0_TCR_off = 1 1/second)
    TCR_p0_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell calculated based on the number of molecules for TCR and T cell ... [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (TCR_p0_tot = 103.81308067036068 molecule/micrometer^2)
    k_P1_up <- fixed(14.4); label("Rate of antigen uptake by APCs [1/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_P1_up = 14.4 1/day/cell)
    k_xP1_deg <- fixed(2); label("Rate of extracellular antigen degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_xP1_deg = 2 1/day)
    k_P1_deg <- fixed(17.28); label("Rate of endosomal antigen degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_P1_deg = 17.28 1/day)
    k_p1_deg <- fixed(144); label("Rate of endosomal epitope degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_p1_deg = 144 1/day)
    k_P1_on <- fixed(2.39117576834e-19); label("Rate of antigen-MHC binding [L/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_P1_on = 144000 1/day/molarity)
    k_P1_d1 <- fixed(3.34228880445e+16); label("Antigen-MHC kd [count/L]")  # S2 Table, sheet 'QSP default parameters' (k_P1_d1 = 5.55e-08 molarity); S1 Table sampling range 'log normal [2.7e-9, 1]'
    p1_50 <- fixed(264550.26455); label("calculated based on the number of molecules for half-maximal T cell activation and synapse surface area [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (p1_50 = 2.6455026455026456e-05 molecule/micrometer^2)
    P1_C1 <- fixed(6.02214199e+15); label("Concentration of P1 in C1 [1/L]")  # S2 Table, sheet 'QSP default parameters' (P1_C1 = 1e-08 molarity/cell)
    k_M1p1_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (k_M1p1_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p1_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_M1p1_TCR_off = 1 1/second)
    k_M1p1_TCR_p <- fixed(86400); label("Rate of MHC-peptide-TCR complex modification [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_M1p1_TCR_p = 1 1/second)
    phi_M1p1_TCR <- fixed(7776); label("Rate of MHC-peptide-TCR complex with maximal modification that leads to non-signaling [1/day]")  # S2 Table, sheet 'QSP default parameters' (phi_M1p1_TCR = 0.09 1/second)
    N_M1p1_TCR <- fixed(10); label("Number of modifcation steps for MHC-peptide-TCR complex [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (N_M1p1_TCR = 10 dimensionless)
    TCR_p1_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell calculated based on the number of molecules for TCR and T cell ... [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (TCR_p1_tot = 103.81308067036068 molecule/micrometer^2)
    kon_PD1_PDL1 <- fixed(8.3691151892e-07); label("kon of PD1-PDL1 binding kon for PD1-PDL1 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_PD1_PDL1 = 0.058333333333333334 1/(micromolarity*nanometer*second))
    q_P_nivo <- fixed(0.76032); label("Volumetric flow rate of nivo into peripheral tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_P_nivo = 8.8e-06 liter/second)
    q_T_nivo <- fixed(0.00736128); label("Volumetric flow rate of nivo into tumor tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_T_nivo = 8.52e-05 milliliter/second)
    q_LN_nivo <- fixed(0.0002808); label("Volumetric flow rate of nivo into lymph nodes [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_LN_nivo = 3.25e-06 milliliter/second)
    q_LD_nivo <- fixed(2.16); label("Lymphatic flow rate of nivo from tumor to LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_LD_nivo = 0.0015 1/minute)
    k_cl_nivo <- fixed(0.37); label("Clearance rate of nivo [L/day]")  # S2 Table, sheet 'QSP default parameters' (k_cl_nivo = 0.37 liter/day)
    gamma_C_nivo <- fixed(0.61); label("Volume fraction available to nivo in central compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_C_nivo = 0.61 dimensionless)
    gamma_P_nivo <- fixed(0.0452); label("Volume fraction available to nivo in peripheral compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_P_nivo = 0.0452 dimensionless)
    gamma_T_nivo <- fixed(0.522); label("Volume fraction available to nivo in tumor compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_T_nivo = 0.522 dimensionless)
    gamma_LN_nivo <- fixed(0.2); label("Volume fraction available to nivo in lymph node compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_LN_nivo = 0.2 dimensionless)
    q_P_durv <- fixed(0.76032); label("Volumetric flow rate of durv into peripheral tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_P_durv = 8.8e-06 1/second); the printed 1/second is inconsistent with the rate law, so the SBML-deposit unit (liter/second) is used
    q_T_durv <- fixed(0.00736128); label("Volumetric flow rate of durv into tumor tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_T_durv = 8.52e-05 1/second); the printed 1/second is inconsistent with the rate law, so the SBML-deposit unit (milliliter/second) is used
    q_LN_durv <- fixed(0.0002808); label("Volumetric flow rate of durv into lymph nodes [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_LN_durv = 3.25e-06 1/second); the printed 1/second is inconsistent with the rate law, so the SBML-deposit unit (milliliter/second) is used
    q_LD_durv <- fixed(2.16); label("Lymphatic flow rate of durv from tumor to LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_LD_durv = 0.0015 1/minute)
    k_cl_durv <- fixed(0.4); label("Clearance rate of durv [L/day]")  # S2 Table, sheet 'QSP default parameters' (k_cl_durv = 0.4 1/second); the printed 1/second is inconsistent with the rate law, so the SBML-deposit unit (liter/day) is used
    gamma_C_durv <- fixed(0.58); label("Volume fraction available to durv in central compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_C_durv = 0.58 dimensionless)
    gamma_P_durv <- fixed(0.056); label("Volume fraction available to durv in peripheral compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_P_durv = 0.056 dimensionless)
    gamma_T_durv <- fixed(0.522); label("Volume fraction available to durv in tumor compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_T_durv = 0.522 dimensionless)
    gamma_LN_durv <- fixed(0.2); label("Volume fraction available to durv in lymph node compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_LN_durv = 0.2 dimensionless)
    q_P_ipi <- fixed(0.9936); label("Volumetric flow rate of ipi into peripheral tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_P_ipi = 1.15e-05 liter/second)
    q_T_ipi <- fixed(0.00736128); label("Volumetric flow rate of ipi into tumor tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_T_ipi = 8.52e-05 milliliter/second)
    q_LN_ipi <- fixed(0.0002808); label("Volumetric flow rate of ipi into lymph nodes [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_LN_ipi = 3.25e-06 milliliter/second)
    q_LD_ipi <- fixed(2.16); label("Lymphatic flow rate of ipi from tumor to LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_LD_ipi = 0.0015 1/minute)
    k_cl_ipi <- fixed(0.524); label("Clearance rate of ipi [L/day]")  # S2 Table, sheet 'QSP default parameters' (k_cl_ipi = 0.524 liter/day)
    gamma_C_ipi <- fixed(0.686); label("Volume fraction available to ipi in central compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_C_ipi = 0.686 dimensionless)
    gamma_P_ipi <- fixed(0.0496); label("Volume fraction available to ipi in peripheral compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_P_ipi = 0.0496 dimensionless)
    gamma_T_ipi <- fixed(0.522); label("Volume fraction available to ipi in tumor compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_T_ipi = 0.522 dimensionless)
    gamma_LN_ipi <- fixed(0.2); label("Volume fraction available to ipi in lymph node compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_LN_ipi = 0.2 dimensionless)
    kon_PD1_PDL2 <- fixed(1.09994085344e-06); label("kon of PD1-PDL2 binding kon for PD1-PDL2 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_PD1_PDL2 = 0.07666666666666666 1/(micromolarity*nanometer*second))
    kon_PD1_nivo <- fixed(1.86511709931e-13); label("kon of PD1-nivolumab binding kon for PD1-nivolumab [L/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_PD1_nivo = 1300000 1/(molarity*second))
    kon_PDL1_durv <- fixed(6.16923348232e-14); label("kon of PDL1-durvalumab binding kon for PDL1-durvalumab [L/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_PDL1_durv = 430000 1/(molarity*second))
    kon_CD28_CD80 <- fixed(1.91294061467e-06); label("kon of CD28-CD80 binding kon for CD28-CD80 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CD28_CD80 = 0.13333333333333333 1/(micromolarity*nanometer*second))
    kon_CD28_CD86 <- fixed(6.69529215136e-06); label("kon of CD28-CD86 binding kon for CD28-CD86 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CD28_CD86 = 0.4666666666666667 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD80 <- fixed(1.05211733807e-05); label("kon of CTLA4-CD80 binding kon for CTLA4-CD80 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CTLA4_CD80 = 0.7333333333333333 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD86 <- fixed(9.56470307337e-06); label("kon of CTLA4-CD86 binding kon for CTLA4-CD86 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CTLA4_CD86 = 0.6666666666666666 1/(micromolarity*nanometer*second))
    kon_CD80_PDL1 <- fixed(1.51122308559e-05); label("kon of CD80-PDL1 binding kon for CD80-PDL1 [dm^2/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CD80_PDL1 = 1.0533333333333332 1/(micromolarity*nanometer*second))
    kon_CTLA4_ipi <- fixed(5.49492191565e-14); label("kon of CTLA4-ipilimumab binding kon for CTLA4-ipilimumab [L/(count*day)]")  # S2 Table, sheet 'QSP default parameters' (kon_CTLA4_ipi = 383000 1/(molarity*second))
    koff_PD1_PDL1 <- fixed(123984); label("koff of PD1-PDL1 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_PD1_PDL1 = 1.435 1/second)
    koff_PD1_PDL2 <- fixed(45705.6); label("koff of PD1-PDL2 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_PD1_PDL2 = 0.529 1/second)
    koff_PD1_nivo <- fixed(292.032); label("koff of PD1-nivolumab binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_PD1_nivo = 0.00338 1/second)
    koff_PDL1_durv <- fixed(2470.608); label("koff of PDL1-durvalumab binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_PDL1_durv = 0.028595 1/second)
    koff_CD28_CD80 <- fixed(138240); label("koff of CD28-CD80 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CD28_CD80 = 1.6 1/second)
    koff_CD28_CD86 <- fixed(2419200); label("koff of CD28-CD86 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CD28_CD86 = 28 1/second)
    koff_CTLA4_CD80 <- fixed(38016); label("koff of CTLA4-CD80 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CTLA4_CD80 = 0.44 1/second)
    koff_CTLA4_CD86 <- fixed(449280); label("koff of CTLA4-CD86 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CTLA4_CD86 = 5.2 1/second)
    koff_CD80_PDL1 <- fixed(518745.6); label("koff of CD80-PDL1 binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CD80_PDL1 = 6.004 1/second)
    koff_CTLA4_ipi <- fixed(602.25984); label("koff of CTLA4-ipilimumab binding calculated based on the measured kd and kon [1/day]")  # S2 Table, sheet 'QSP default parameters' (koff_CTLA4_ipi = 0.0069706 1/second)
    Chi_PD1_nivo <- fixed(333333330.953); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # S2 Table, sheet 'QSP default parameters' (Chi_PD1_nivo = 3.333333309532274 1/nanometer)
    Chi_PDL1_durv <- fixed(3333333309.53); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # S2 Table, sheet 'QSP default parameters' (Chi_PDL1_durv = 33.33333309532274 1/nanometer)
    Chi_CTLA4_ipi <- fixed(333333330.953); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # S2 Table, sheet 'QSP default parameters' (Chi_CTLA4_ipi = 3.333333309532274 1/nanometer)
    PD1_50 <- fixed(60000000000); label("PD1/PDL1 concentration for half-maximal T cell inactivation [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (PD1_50 = 6 molecule/micrometer^2)
    n_PD1 <- fixed(2); label("Hill coefficient for PD1/PDL1 half-maximal T cell inactivation [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_PD1 = 2 dimensionless)
    CD28_CD8X_50 <- fixed(1.70045e+13); label("CD28-CD80/CD28-CD86 concentration for half-maximal T cell co-estimulation [count/dm^2]")  # S2 Table, sheet 'QSP default parameters' (CD28_CD8X_50 = 1700.45 molecule/micrometer^2); S1 Table sampling range 'log unif [100, 4000]'
    n_CD28_CD8X <- fixed(2); label("Hill coefficient for CD28-CD80/CD28-CD86 half-maximal T cell co-estimulation [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_CD28_CD8X = 2 dimensionless)
    T_PD1_total <- fixed(60000); label("concentration of PD1 on T cells [count]")  # S2 Table, sheet 'QSP default parameters' (T_PD1_total = 60000 molecule)
    T_CD28_total <- fixed(184000); label("concentration of CD28 on T cells [count]")  # S2 Table, sheet 'QSP default parameters' (T_CD28_total = 184000 molecule)
    T_CTLA4_syn <- fixed(400); label("concentration of CTLA4 on T cells [count]")  # S2 Table, sheet 'QSP default parameters' (T_CTLA4_syn = 400 molecule)
    T_PDL1_total <- fixed(1600000); label("concentration of PDL1 on T cells [count]")  # S2 Table, sheet 'QSP default parameters' (T_PDL1_total = 1600000 molecule)
    C1_PDL1_total <- fixed(387151); label("number of PDL1 molecules per C1 cell [count]")  # S2 Table, sheet 'QSP default parameters' (C1_PDL1_total = 387151 molecule); S1 Table sampling range 'log unif [1.8e5, 1.6e6]'
    C1_PDL2_total <- fixed(40000); label("number of PDL2 molecules per C1 cell [count]")  # S2 Table, sheet 'QSP default parameters' (C1_PDL2_total = 40000 molecule)
    C1_CD80_total <- fixed(40000); label("number of CD80 molecules per C1 cell [count]")  # S2 Table, sheet 'QSP default parameters' (C1_CD80_total = 40000 molecule)
    C1_CD86_total <- fixed(860000); label("number of CD86 molecules per C1 cell [count]")  # S2 Table, sheet 'QSP default parameters' (C1_CD86_total = 860000 molecule)
    APC_PDL1_total <- fixed(1600000); label("number of PDL1 molecules per APC cell [count]")  # S2 Table, sheet 'QSP default parameters' (APC_PDL1_total = 1600000 molecule)
    APC_PDL2_total <- fixed(2000); label("number of PDL2 molecules per APC cell [count]")  # S2 Table, sheet 'QSP default parameters' (APC_PDL2_total = 2000 molecule)
    APC_CD80_total <- fixed(40000); label("number of CD80 molecules per APC cell [count]")  # S2 Table, sheet 'QSP default parameters' (APC_CD80_total = 40000 molecule)
    APC_CD86_total <- fixed(860000); label("number of CD86 molecules per APC cell [count]")  # S2 Table, sheet 'QSP default parameters' (APC_CD86_total = 860000 molecule)
    Treg_CTLA4_tot <- fixed(5000); label("Total number of CTLA4 on Treg cells [count]")  # S2 Table, sheet 'QSP default parameters' (Treg_CTLA4_tot = 5000 molecule)
    Treg_CTLA4_50 <- fixed(617.104); label("CTLA4 occupancy for half-maximal Treg inactivation by macrophages [count]")  # S2 Table, sheet 'QSP default parameters' (Treg_CTLA4_50 = 617.104 molecule); S1 Table sampling range 'log unif [500, 1000]'
    n_Treg_CTLA4 <- fixed(2); label("CTLA4 occupancy Hill coefficient for Treg inactivation by macrophages [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (n_Treg_CTLA4 = 2 dimensionless)
    k_CTLA4_ADCC <- fixed(0.289); label("Anti-CTLA4 ADCC (antibody-dependent cellular cytotoxicity) rate of Treg [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_CTLA4_ADCC = 0.289 1/day); S1 Table sampling range 'log unif [0.1, 1]'
    k_rec_MDSC <- fixed(1.2); label("Rate of MDSC recruitment into the tumor [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_rec_MDSC = 1.2 1/day)
    kd_MDSC <- fixed(0.015); label("Rate of MDSC death [1/day]")  # S2 Table, sheet 'QSP default parameters' (kd_MDSC = 0.015 1/day)
    IC50_ENT_C <- fixed(4.34256658899e+16); label("ENT concentration for half-maximal tumor cell death [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_ENT_C = 7.211e-08 molarity); S1 Table sampling range 'log unif [2.5e-9, 3.74e-7]'
    k_deg_CCL2 <- fixed(1.44); label("rate of CCL2 degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_deg_CCL2 = 0.06 1/hour)
    k_deg_NO <- fixed(135); label("rate of NO degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_deg_NO = 135 1/day)
    k_deg_ArgI <- fixed(0.173); label("rate of ArgI degradation [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_deg_ArgI = 0.173 1/day)
    k_sec_CCL2 <- fixed(81437.4261308); label("rate of CCL2 secretion [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_sec_CCL2 = 1.3523e-10 nanomole/cell/day); S1 Table sampling range 'linear unif [8.2e-11, 2.02e-10]'
    k_sec_NO <- fixed(240517124.51); label("rate of NO secretion from MDSCs [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_sec_NO = 3.99388e-07 nanomole/cell/day); S1 Table sampling range 'linear unif [3.95e-7, 5.65e-7]'
    k_sec_ArgI <- fixed(7.90370003336e+15); label("rate of ArgI secretion [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_sec_ArgI = 0.0131244 (mU*microliter)/cell/day); S1 Table sampling range 'linear unif [0.012, 0.016]'
    IC50_ENT_NO <- fixed(1.13063307005e+15); label("half-maximal ENT concentration for NO inhibition [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_ENT_NO = 1.87746e-09 molarity); S1 Table sampling range 'log unif [5.6e-11, 5.6e-9]'
    ki_Treg <- fixed(1.40217); label("rate of ArgI-induced Treg expension [1/day]")  # S2 Table, sheet 'QSP default parameters' (ki_Treg = 1.40217 1/day); S1 Table sampling range 'log unif [0.27, 2.7]'
    IC50_ArgI_CTL <- fixed(2.05462235986e+26); label("rate of ArgI-induced T cell death [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_ArgI_CTL = 341.178 mU); S1 Table sampling range 'log unif [100, 617]'
    IC50_NO_CTL <- fixed(1.61862530193e+15); label("rate of NO-induced T cell death [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_NO_CTL = 2.68779e-09 molarity); S1 Table sampling range 'log unif [7.5e-11, 7.5e-9]'
    EC50_CCL2_rec <- fixed(3.011070995e+14); label("Half-maximal CCL2 level of MDSC recruitment [count/L]")  # S2 Table, sheet 'QSP default parameters' (EC50_CCL2_rec = 5e-10 molarity)
    EC50_ArgI_Treg <- fixed(2.85281512564e+25); label("Half-maximal ArgI level of Treg expansion [count/L]")  # S2 Table, sheet 'QSP default parameters' (EC50_ArgI_Treg = 47.3721 mU); S1 Table sampling range 'log unif [22.1, 221]'
    MDSC_max <- fixed(166432000); label("Maximal MDSC density in the tumour [count/L]")  # S2 Table, sheet 'QSP default parameters' (MDSC_max = 166432 cell/milliliter); S1 Table sampling range 'log normal [1.637e5, 1]'
    Treg_max <- fixed(900000000); label("Maximal FoxP3+ T cell density in the tumor [count/L]")  # S2 Table, sheet 'QSP default parameters' (Treg_max = 900000 cell/milliliter); S1 Table sampling range 'log normal [9.55e3, 1]'
    IC50_ENT_CCL2 <- fixed(9.81500745814e+14); label("half-maximal ENT concentration for CCL2 inhibition [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_ENT_CCL2 = 1.62982e-09 molarity); S1 Table sampling range 'log unif [1.2e-10, 1.2e-8]'
    k_brec_MDSC <- fixed(0.0021); label("Baseline rate of MDSC migration into the tumor [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_brec_MDSC = 0.0021 1/day)
    IC50_ENT_ArgI <- fixed(2.25228110426e+17); label("half-maximal ENT concentration for Arg I inhibition [count/L]")  # S2 Table, sheet 'QSP default parameters' (IC50_ENT_ArgI = 3.74e-07 molarity)
    k_a1_ENT <- fixed(50.33472); label("Buccal absorption rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_a1_ENT = 2.09728 1/hour); S1 Table sampling range 'linear unif [1.9, 2.24]'
    k_a2_ENT <- fixed(122.98296); label("GI absorption rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_a2_ENT = 5.12429 1/hour); S1 Table sampling range 'linear unif [2.47, 7.47]'
    k_cln_ENT <- fixed(1.33664380273e+21); label("Non-linear clearance rate of entinostat [count/(L*day)]")  # S2 Table, sheet 'QSP default parameters' (k_cln_ENT = 9.24812e-05 mole/(liter*hour)); S1 Table sampling range 'linear unif [6.11e-5, 1.22e-4]'
    Kc_ENT <- fixed(5.64763049312e+18); label("Kd for non-linear clearnance [count/L]")  # S2 Table, sheet 'QSP default parameters' (Kc_ENT = 9.378109155345853e-06 mole/liter)
    lagP <- fixed(0.151966666667); label("Tlag (entinostat) [day]")  # S2 Table, sheet 'QSP default parameters' (lagP = 3.6472 hour)
    durP <- fixed(0.0074875); label("Duration of zero-order absorption (entinostat) [day]")  # S2 Table, sheet 'QSP default parameters' (durP = 0.1797 hour)
    k_dose2 <- fixed(0.3312); label("First-order absorption rate [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_dose2 = 0.0138 1/hour)
    q_P_ENT <- fixed(27648); label("Volumetric flow rate of entinostat into peripheral tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_P_ENT = 320 milliliter/second)
    q_T_ENT <- fixed(286112.736); label("Volumetric flow rate of entinostat into tumor tissues [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_T_ENT = 3311.49 milliliter/second); S1 Table sampling range 'linear unif [3200, 3420]'
    q_LN_ENT <- fixed(27648); label("Volumetric flow rate of entinostat into lymph nodes [L/day]")  # S2 Table, sheet 'QSP default parameters' (q_LN_ENT = 320 milliliter/second)
    q_LD_ENT <- fixed(2.16); label("Lymphatic flow rate of entinostat from tumor to LN [1/day]")  # S2 Table, sheet 'QSP default parameters' (q_LD_ENT = 0.0015 1/minute)
    k_cl_ENT <- fixed(14.500176); label("Clearance rate of entinostat [1/day]")  # S2 Table, sheet 'QSP default parameters' (k_cl_ENT = 0.604174 1/hour); S1 Table sampling range 'log unif [0.152, 1.27]'
    gamma_C_ENT <- fixed(0.55); label("Volume fraction available to entinostat in central compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_C_ENT = 0.55 dimensionless)
    gamma_P_ENT <- fixed(0.062); label("Volume fraction available to entinostat in peripheral compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_P_ENT = 0.062 dimensionless)
    gamma_T_ENT <- fixed(0.611); label("Volume fraction available to entinostat in tumor compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_T_ENT = 0.611 dimensionless)
    gamma_LN_ENT <- fixed(0.2); label("Volume fraction available to entinostat in lymph node compartment [dimensionless]")  # S2 Table, sheet 'QSP default parameters' (gamma_LN_ENT = 0.2 dimensionless)

    # ---- entinostat dose split (Wang 2020 Figure 2A; value not printed in either paper) ----
    F_ENT_buccal <- fixed(0.324); label("Fraction of an entinostat dose absorbed through the zero-order buccal route [dimensionless]")  # not printed in Ruiz-Martinez 2022 (entinostat is not dosed there); carried unchanged from Wang_2020_entinostat_nivolumab_ipilimumab_qsp, where it was back-solved by the maintainers from the Wang 2020 entinostat Cmax
  })
  model({
    # Unit system: amounts in counts (cells, molecules; the export defines
    # 1 cell = 1 molecule = 1.66053872801495e-24 mol), volumes in L, areas in
    # dm^2, time in days. Drug states are held in nmol so dose records are in
    # nmol; nmol_to_molecule converts them to counts inside the rate laws.
    nmol_to_molecule <- 602214199000000
    # ---- parameter aliases (rxode2 reads a bare theta +/- theta as a
    # mu-reference block; these mechanistic rate laws legitimately add two
    # fixed constants, so each such parameter is copied to a local first) ----
    a_k_C1_death <- k_C1_death
    a_k_C1_therapy <- k_C1_therapy
    a_k_M1p1_TCR_off <- k_M1p1_TCR_off
    a_k_M1p1_TCR_p <- k_M1p1_TCR_p
    a_phi_M1p1_TCR <- phi_M1p1_TCR

    # ---- compartment capacities ----
    m_V_C <- vol_V_C
    m_V_P <- vol_V_P
    m_V_LN <- vol_V_LN
    m_V_e <- vol_V_e
    m_A_e <- vol_A_e
    m_A_s <- vol_A_s
    m_syn_T_C1 <- vol_syn_T_C1
    m_syn_T_APC <- vol_syn_T_APC

    # ---- species quantities: states are AMOUNTS (molecules/cells, or nmol
    # for the drugs); x_* are the SimBiology species values in model units ----
    x_V_C_T0 <- q_V_C_T0
    x_V_C_T1 <- q_V_C_T1
    x_V_P_T0 <- q_V_P_T0
    x_V_P_T1 <- q_V_P_T1
    x_V_P_Treg_CTLA4 <- q_V_P_Treg_CTLA4
    x_V_P_Treg_CTLA4_ipi <- q_V_P_Treg_CTLA4_ipi
    x_V_P_Treg_CTLA4_ipi_CTLA4 <- q_V_P_Treg_CTLA4_ipi_CTLA4
    x_V_T_C_x <- q_V_T_C_x
    x_V_T_T_exh <- q_V_T_T_exh
    x_V_T_C1 <- q_V_T_C1
    x_V_T_T0 <- q_V_T_T0
    x_V_T_T1 <- q_V_T_T1
    x_V_T_APC <- q_V_T_APC
    x_V_T_mAPC <- q_V_T_mAPC
    x_V_T_Treg_CTLA4 <- q_V_T_Treg_CTLA4
    x_V_T_Treg_CTLA4_ipi <- q_V_T_Treg_CTLA4_ipi
    x_V_T_Treg_CTLA4_ipi_CTLA4 <- q_V_T_Treg_CTLA4_ipi_CTLA4
    x_V_T_MDSC <- q_V_T_MDSC
    x_V_LN_nT0 <- q_V_LN_nT0
    x_V_LN_aT0 <- q_V_LN_aT0
    x_V_LN_T0 <- q_V_LN_T0
    x_V_LN_nT1 <- q_V_LN_nT1
    x_V_LN_aT1 <- q_V_LN_aT1
    x_V_LN_T1 <- q_V_LN_T1
    x_V_LN_APC <- q_V_LN_APC
    x_V_LN_mAPC <- q_V_LN_mAPC
    vol_V_T <- V_Tmin+(vol_cell*x_V_T_C_x)+(vol_Tcell*x_V_T_T_exh)+(vol_cell*x_V_T_C1)+(vol_Tcell*x_V_T_T0)+(vol_Tcell*x_V_T_T1)   # S1 Table S4 Rules, rule 1 (tumour volume)
    x_V_C_nivo <- q_V_C_nivo*nmol_to_molecule/m_V_C
    x_V_C_durv <- q_V_C_durv*nmol_to_molecule/m_V_C
    x_V_C_ipi <- q_V_C_ipi*nmol_to_molecule/m_V_C
    x_V_C_ENT <- q_V_C_ENT*nmol_to_molecule/m_V_C
    x_V_C_ENT_Buccal <- q_V_C_ENT_Buccal*nmol_to_molecule/m_V_C
    x_V_C_ENT_GI <- q_V_C_ENT_GI*nmol_to_molecule/m_V_C
    x_V_C_Dose2 <- q_V_C_Dose2*nmol_to_molecule/m_V_C
    x_V_P_nivo <- q_V_P_nivo*nmol_to_molecule/m_V_P
    x_V_P_durv <- q_V_P_durv*nmol_to_molecule/m_V_P
    x_V_P_ipi <- q_V_P_ipi*nmol_to_molecule/m_V_P
    x_V_P_ENT <- q_V_P_ENT*nmol_to_molecule/m_V_P
    x_V_T_c <- q_V_T_c/vol_V_T
    x_V_T_nivo <- q_V_T_nivo*nmol_to_molecule/vol_V_T
    x_V_T_durv <- q_V_T_durv*nmol_to_molecule/vol_V_T
    x_V_T_ipi <- q_V_T_ipi*nmol_to_molecule/vol_V_T
    x_V_T_CCL2 <- q_V_T_CCL2/vol_V_T
    x_V_T_NO <- q_V_T_NO/vol_V_T
    x_V_T_ArgI <- q_V_T_ArgI/vol_V_T
    x_V_T_ENT <- q_V_T_ENT*nmol_to_molecule/vol_V_T
    x_V_LN_IL2 <- q_V_LN_IL2/m_V_LN
    x_V_LN_P0 <- q_V_LN_P0/m_V_LN
    x_V_LN_P1 <- q_V_LN_P1/m_V_LN
    x_V_LN_nivo <- q_V_LN_nivo*nmol_to_molecule/m_V_LN
    x_V_LN_durv <- q_V_LN_durv*nmol_to_molecule/m_V_LN
    x_V_LN_ipi <- q_V_LN_ipi*nmol_to_molecule/m_V_LN
    x_V_LN_ENT <- q_V_LN_ENT*nmol_to_molecule/m_V_LN
    x_V_e_P0 <- q_V_e_P0/m_V_e
    x_V_e_p0 <- q_V_e_p0/m_V_e
    x_V_e_P1 <- q_V_e_P1/m_V_e
    x_V_e_p1 <- q_V_e_p1/m_V_e
    x_A_e_M1 <- q_A_e_M1/m_A_e
    x_A_e_M1p0 <- q_A_e_M1p0/m_A_e
    x_A_e_M1p1 <- q_A_e_M1p1/m_A_e
    x_A_s_M1 <- q_A_s_M1/m_A_s
    x_A_s_M1p0 <- q_A_s_M1p0/m_A_s
    x_A_s_M1p1 <- q_A_s_M1p1/m_A_s
    x_syn_T_C1_PD1_PDL1 <- q_syn_T_C1_PD1_PDL1/m_syn_T_C1
    x_syn_T_C1_PD1_PDL2 <- q_syn_T_C1_PD1_PDL2/m_syn_T_C1
    x_syn_T_C1_PD1 <- q_syn_T_C1_PD1/m_syn_T_C1
    x_syn_T_C1_PDL1 <- q_syn_T_C1_PDL1/m_syn_T_C1
    x_syn_T_C1_PDL2 <- q_syn_T_C1_PDL2/m_syn_T_C1
    x_syn_T_C1_PD1_nivo <- q_syn_T_C1_PD1_nivo/m_syn_T_C1
    x_syn_T_C1_PD1_nivo_PD1 <- q_syn_T_C1_PD1_nivo_PD1/m_syn_T_C1
    x_syn_T_C1_PDL1_durv <- q_syn_T_C1_PDL1_durv/m_syn_T_C1
    x_syn_T_C1_PDL1_durv_PDL1 <- q_syn_T_C1_PDL1_durv_PDL1/m_syn_T_C1
    x_syn_T_C1_TPDL1 <- q_syn_T_C1_TPDL1/m_syn_T_C1
    x_syn_T_C1_TPDL1_durv <- q_syn_T_C1_TPDL1_durv/m_syn_T_C1
    x_syn_T_C1_TPDL1_durv_TPDL1 <- q_syn_T_C1_TPDL1_durv_TPDL1/m_syn_T_C1
    x_syn_T_C1_CD28_CD80 <- q_syn_T_C1_CD28_CD80/m_syn_T_C1
    x_syn_T_C1_CD28_CD80_CD28 <- q_syn_T_C1_CD28_CD80_CD28/m_syn_T_C1
    x_syn_T_C1_CD28_CD86 <- q_syn_T_C1_CD28_CD86/m_syn_T_C1
    x_syn_T_C1_CD80_CTLA4 <- q_syn_T_C1_CD80_CTLA4/m_syn_T_C1
    x_syn_T_C1_CD80_CTLA4_CD80 <- q_syn_T_C1_CD80_CTLA4_CD80/m_syn_T_C1
    x_syn_T_C1_CTLA4_CD80_CTLA4 <- q_syn_T_C1_CTLA4_CD80_CTLA4/m_syn_T_C1
    x_syn_T_C1_CD80_CTLA4_CD80_CTLA4 <- q_syn_T_C1_CD80_CTLA4_CD80_CTLA4/m_syn_T_C1
    x_syn_T_C1_CD86_CTLA4 <- q_syn_T_C1_CD86_CTLA4/m_syn_T_C1
    x_syn_T_C1_CD86_CTLA4_CD86 <- q_syn_T_C1_CD86_CTLA4_CD86/m_syn_T_C1
    x_syn_T_C1_TPDL1_CD80 <- q_syn_T_C1_TPDL1_CD80/m_syn_T_C1
    x_syn_T_C1_TPDL1_CD80_TPDL1 <- q_syn_T_C1_TPDL1_CD80_TPDL1/m_syn_T_C1
    x_syn_T_C1_CD28 <- q_syn_T_C1_CD28/m_syn_T_C1
    x_syn_T_C1_CTLA4 <- q_syn_T_C1_CTLA4/m_syn_T_C1
    x_syn_T_C1_CD80 <- q_syn_T_C1_CD80/m_syn_T_C1
    x_syn_T_C1_CD86 <- q_syn_T_C1_CD86/m_syn_T_C1
    x_syn_T_C1_CTLA4_ipi <- q_syn_T_C1_CTLA4_ipi/m_syn_T_C1
    x_syn_T_C1_CTLA4_ipi_CTLA4 <- q_syn_T_C1_CTLA4_ipi_CTLA4/m_syn_T_C1
    x_syn_T_APC_PD1_PDL1 <- q_syn_T_APC_PD1_PDL1/m_syn_T_APC
    x_syn_T_APC_PD1_PDL2 <- q_syn_T_APC_PD1_PDL2/m_syn_T_APC
    x_syn_T_APC_PD1 <- q_syn_T_APC_PD1/m_syn_T_APC
    x_syn_T_APC_PDL1 <- q_syn_T_APC_PDL1/m_syn_T_APC
    x_syn_T_APC_PDL2 <- q_syn_T_APC_PDL2/m_syn_T_APC
    x_syn_T_APC_PD1_nivo <- q_syn_T_APC_PD1_nivo/m_syn_T_APC
    x_syn_T_APC_PD1_nivo_PD1 <- q_syn_T_APC_PD1_nivo_PD1/m_syn_T_APC
    x_syn_T_APC_PDL1_durv <- q_syn_T_APC_PDL1_durv/m_syn_T_APC
    x_syn_T_APC_PDL1_durv_PDL1 <- q_syn_T_APC_PDL1_durv_PDL1/m_syn_T_APC
    x_syn_T_APC_TPDL1 <- q_syn_T_APC_TPDL1/m_syn_T_APC
    x_syn_T_APC_TPDL1_durv <- q_syn_T_APC_TPDL1_durv/m_syn_T_APC
    x_syn_T_APC_TPDL1_durv_TPDL1 <- q_syn_T_APC_TPDL1_durv_TPDL1/m_syn_T_APC
    x_syn_T_APC_CD28_CD80 <- q_syn_T_APC_CD28_CD80/m_syn_T_APC
    x_syn_T_APC_CD28_CD80_CD28 <- q_syn_T_APC_CD28_CD80_CD28/m_syn_T_APC
    x_syn_T_APC_CD28_CD86 <- q_syn_T_APC_CD28_CD86/m_syn_T_APC
    x_syn_T_APC_CD80_CTLA4 <- q_syn_T_APC_CD80_CTLA4/m_syn_T_APC
    x_syn_T_APC_CD80_CTLA4_CD80 <- q_syn_T_APC_CD80_CTLA4_CD80/m_syn_T_APC
    x_syn_T_APC_CTLA4_CD80_CTLA4 <- q_syn_T_APC_CTLA4_CD80_CTLA4/m_syn_T_APC
    x_syn_T_APC_CD80_CTLA4_CD80_CTLA4 <- q_syn_T_APC_CD80_CTLA4_CD80_CTLA4/m_syn_T_APC
    x_syn_T_APC_CD86_CTLA4 <- q_syn_T_APC_CD86_CTLA4/m_syn_T_APC
    x_syn_T_APC_CD86_CTLA4_CD86 <- q_syn_T_APC_CD86_CTLA4_CD86/m_syn_T_APC
    x_syn_T_APC_TPDL1_CD80 <- q_syn_T_APC_TPDL1_CD80/m_syn_T_APC
    x_syn_T_APC_TPDL1_CD80_TPDL1 <- q_syn_T_APC_TPDL1_CD80_TPDL1/m_syn_T_APC
    x_syn_T_APC_CD28 <- q_syn_T_APC_CD28/m_syn_T_APC
    x_syn_T_APC_CTLA4 <- q_syn_T_APC_CTLA4/m_syn_T_APC
    x_syn_T_APC_CD80 <- q_syn_T_APC_CD80/m_syn_T_APC
    x_syn_T_APC_CD86 <- q_syn_T_APC_CD86/m_syn_T_APC
    x_syn_T_APC_CTLA4_ipi <- q_syn_T_APC_CTLA4_ipi/m_syn_T_APC
    x_syn_T_APC_CTLA4_ipi_CTLA4 <- q_syn_T_APC_CTLA4_ipi_CTLA4/m_syn_T_APC

    # ---- repeated-assignment rules (S1 Table, sheet 'Table S4. Rules') ----
    C_total <- x_V_T_C1   # S1 Table S4 Rules (C_total)
    T_total <- x_V_T_T0+x_V_T_T1   # S1 Table S4 Rules (T_total)
    T_total_LN <- x_V_LN_T1   # S1 Table S4 Rules (T_total_LN)
    Tregs_ <- x_V_T_T0   # S1 Table S4 Rules (Tregs_)
    H_APC <- (n_sites_APC*x_V_LN_APC)/((n_sites_APC*x_V_LN_APC)+T_total_LN+1)   # S1 Table S4 Rules (H_APC)
    H_mAPC <- (n_sites_APC*x_V_LN_mAPC)/((n_sites_APC*x_V_LN_mAPC)+T_total_LN+1)   # S1 Table S4 Rules (H_mAPC)
    pTCR_p0_MHC_tot <- 0.5*(((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))-(TCR_p0_tot*sqrt((((((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))/TCR_p0_tot))^2-(((4*x_A_s_M1p0)/n_T0_clones)/TCR_p0_tot)))))   # S1 Table S4 Rules (pTCR_p0_MHC_tot)
    H_P0 <- pTCR_p0_MHC_tot/(pTCR_p0_MHC_tot+p0_50)   # S1 Table S4 Rules (H_P0)
    pTCR_p1_MHC_tot <- (a_k_M1p1_TCR_off/(a_k_M1p1_TCR_off+a_phi_M1p1_TCR))*((a_k_M1p1_TCR_p/(a_k_M1p1_TCR_off+a_k_M1p1_TCR_p)))^N_M1p1_TCR*0.5*(((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))-(TCR_p1_tot*sqrt((((((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))/TCR_p1_tot))^2-(((4*x_A_s_M1p1)/n_T1_clones)/TCR_p1_tot)))))   # S1 Table S4 Rules (pTCR_p1_MHC_tot)
    H_P1 <- pTCR_p1_MHC_tot/(pTCR_p1_MHC_tot+p1_50)   # S1 Table S4 Rules (H_P1)
    H_PD1_C1 <- (((x_syn_T_C1_PD1_PDL1+x_syn_T_C1_PD1_PDL2)/PD1_50))^n_PD1/((((x_syn_T_C1_PD1_PDL1+x_syn_T_C1_PD1_PDL2)/PD1_50))^n_PD1+1)   # S1 Table S4 Rules (H_PD1_C1)
    H_CD28_C1 <- (((x_syn_T_C1_CD28_CD80+x_syn_T_C1_CD28_CD86+(2*x_syn_T_C1_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X/((((x_syn_T_C1_CD28_CD80+x_syn_T_C1_CD28_CD86+(2*x_syn_T_C1_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X+1)   # S1 Table S4 Rules (H_CD28_C1)
    H_PD1_APC <- (((x_syn_T_APC_PD1_PDL1+x_syn_T_APC_PD1_PDL2)/PD1_50))^n_PD1/((((x_syn_T_APC_PD1_PDL1+x_syn_T_APC_PD1_PDL2)/PD1_50))^n_PD1+1)   # S1 Table S4 Rules (H_PD1_APC)
    H_CD28_APC <- (((x_syn_T_APC_CD28_CD80+x_syn_T_APC_CD28_CD86+(2*x_syn_T_APC_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X/((((x_syn_T_APC_CD28_CD80+x_syn_T_APC_CD28_CD86+(2*x_syn_T_APC_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X+1)   # S1 Table S4 Rules (H_CD28_APC)
    H_Treg_T <- (((x_V_T_Treg_CTLA4_ipi+(2*x_V_T_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4/((((x_V_T_Treg_CTLA4_ipi+(2*x_V_T_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4+1)   # S1 Table S4 Rules (H_Treg_T)
    H_Treg_P <- (((x_V_P_Treg_CTLA4_ipi+(2*x_V_P_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4/((((x_V_P_Treg_CTLA4_ipi+(2*x_V_P_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4+1)   # S1 Table S4 Rules (H_Treg_P)
    H_MDSC_C1 <- 1-((1-(x_V_T_NO/(IC50_NO_CTL+x_V_T_NO)))*(1-(x_V_T_ArgI/(IC50_ArgI_CTL+x_V_T_ArgI)))*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_ArgI))))   # S1 Table S4 Rules (H_MDSC_C1)
    H_ENT_C1 <- x_V_T_ENT/(x_V_T_ENT+IC50_ENT_C)   # S1 Table S4 Rules (H_ENT_C1)
    R_Tcell <- (((k_T1*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*(1-H_PD1_C1)*(1-H_MDSC_C1))   # S1 Table S4 Rules (R_Tcell)
    N_aT <- N0+(N_costim*H_CD28_APC)+((N_IL2*x_V_LN_IL2)/(IL2_50+x_V_LN_IL2))   # S1 Table S4 Rules (N_aT)

    # S1 Table, sheet 'Table S6.Events': cancer cells are set to zero once fewer than 0.9 cells
    # remain. rxode2 has no state-reset events, so proliferation is switched
    # off below that threshold instead and the remaining fraction of a cell
    # decays through the death terms.
    c1_alive <- 1*(x_V_T_C1 >= 0.9)

    # ---- reaction rates (S1 Table, sheet 'Table S3. Reactions', rows numbered
    # identically), amount per day ----
    v1 <- k_cell_clear*x_V_T_C_x   # R1 Clearance of dead cancer cells from tumor
    v2 <- k_cell_clear*x_V_T_T_exh   # R2 Clearance of dead T cells from tumor
    v3 <- (k_C1_growth*x_V_T_C1*(1-(C_total/C_max))*(1-H_ENT_C1))*c1_alive   # R3 Cancer cell growth (Logistic model based on Komarova and reviewed in Marusic)
    v4 <- a_k_C1_death*x_V_T_C1   # R4 Cancer cell death
    v5 <- Q_T0_in*n_T0_clones   # R5 Naive T cell entry into the lymph node
    v6 <- Q_T0_out*x_V_LN_nT0   # R6 Naive T cell exit from the lymph node
    v7 <- k_T0_act*H_APC*H_P0*x_V_LN_nT0   # R7 Naive T cell activation
    v8 <- (k_T0_pro/N_aT)*x_V_LN_aT0   # R8 aT0 cell proliferation
    v9 <- (k_T0_pro/N_aT)*((2)^N_aT-1)*x_V_LN_aT0   # R9 aT0 cell proliferation
    v10 <- k_T0_death*x_V_C_T0   # R10 T cell death in the central compartment
    v11 <- k_T0_death*x_V_P_T0   # R11 T cell death in the peripheral compartment
    v12 <- k_T0_death*x_V_T_T0   # R12 T cell death in the tumour compartment
    v13 <- k_T0_death*x_V_LN_T0   # R13 T cell death in the lymph node compartment
    v14 <- q_T0_P_in*x_V_C_T0   # R14 T cell transport into the peripheral compartment
    v15 <- q_T0_P_out*x_V_P_T0   # R15 T cell transport out of the peripheral compartment
    v16 <- q_T0_T_in*vol_V_T*x_V_C_T0*(1-(x_V_T_T0/(Treg_max*vol_V_T)))   # R16 T cell transport into the tumour compartment
    v17 <- q_T0_LN_out*x_V_LN_T0   # R17 T cell transport out of the lymph node compartment
    v18 <- k_IL2_deg*x_V_LN_IL2*m_V_LN   # R18 IL2 degradation
    v19 <- (k_IL2_cons*T_total_LN*x_V_LN_IL2)/(IL2_50+x_V_LN_IL2)   # R19 IL2 consumption by T cells
    v20 <- (k_IL2_cons*x_V_LN_T0*x_V_LN_IL2)/(IL2_50_Treg+x_V_LN_IL2)   # R20 IL2 consumption by Tregs
    v21 <- Q_T1_in*n_T1_clones   # R21 Naive T cell entry into the lymph node
    v22 <- Q_T1_out*x_V_LN_nT1   # R22 Naive T cell exit from the lymph node
    v23 <- k_T1_act*H_mAPC*H_P1*x_V_LN_nT1   # R23 Naive T cell activation
    v24 <- (k_T1_pro/N_aT)*x_V_LN_aT1   # R24 aT1 cell proliferation
    v25 <- (k_T1_pro/N_aT)*((2)^N_aT-1)*x_V_LN_aT1   # R25 aT1 cell proliferation
    v26 <- k_T1_death*x_V_C_T1   # R26 T cell death in the central compartment
    v27 <- k_T1_death*x_V_P_T1   # R27 T cell death in the peripheral compartment
    v28 <- k_T1_death*x_V_T_T1   # R28 T cell death in the tumour compartment
    v29 <- k_T1_death*x_V_LN_T1   # R29 T cell death in the lymph node compartment
    v30 <- (k_Treg*x_V_T_T1*Tregs_)/(C_total+T_total+1)   # R30 T cell death from Tregs
    v31 <- ((k_T1*x_V_T_T1*C_total)/(C_total+T_total+1))*H_PD1_C1   # R31 T cell death from cancer
    v32 <- q_T1_P_in*x_V_C_T1   # R32 T cell transport into the peripheral compartment
    v33 <- q_T1_P_out*x_V_P_T1   # R33 T cell transport out of the peripheral compartment
    v34 <- q_T1_T_in*vol_V_T*x_V_C_T1   # R34 T cell transport into the tumour compartment
    v35 <- q_T1_LN_out*x_V_LN_T1   # R35 T cell transport out of the lymph node compartment
    v36 <- k_IL2_sec*x_V_LN_aT1   # R36 IL2 secretion from activated T cells
    v37 <- ((k_C_T1*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*(1-H_PD1_C1)*(1-H_MDSC_C1)   # R37 Cancer cell killing by T cells
    v38 <- k_APC_death*((APC0_T*vol_V_T)-x_V_T_APC)   # R38 APC recruitment/death in the tumour
    v39 <- k_APC_death*((APC0_LN*m_V_LN)-x_V_LN_APC)   # R39 APC recruitment/death in LN
    v40 <- ((k_APC_mat*x_V_T_c)/(x_V_T_c+c50))*x_V_T_APC   # R40 APC maturation in the tumour
    v41 <- k_APC_mig*x_V_T_mAPC   # R41 APC migration to the lymph node
    v42 <- k_mAPC_death*x_V_T_mAPC   # R42 mAPC death in the tumour
    v43 <- k_mAPC_death*x_V_LN_mAPC   # R43 mAPC death in the lymph node
    v44 <- k_c*(c0-x_V_T_c)*vol_V_T   # R44 Baseline cytokine secretion/degradation
    v45 <- R_Tcell*DAMPs   # R45 Cytokine release in response to the tumour
    v46 <- (kout*x_A_e_M1*m_A_e)-(kin*x_A_s_M1*m_A_s)   # R46 MHC translocation
    v47 <- n_T0_clones*P0_C1*(a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1)))*x_V_T_C1*vol_V_T   # R47 Antigen deposition from dying cancer cells
    v48 <- k_xP0_deg*x_V_LN_P0*m_V_LN   # R48 Free antigen degradation
    v49 <- k_P0_up*x_V_LN_mAPC*x_V_LN_P0*m_V_LN   # R49 Antigen uptake by mature antigen presenting cells
    v50 <- k_P0_up*1*x_V_LN_P0*m_V_e   # R50 Antigen uptake by mature antigen presenting cells
    v51 <- k_P0_deg*x_V_e_P0*m_V_e   # R51 Antigen degradation in APC endosomes
    v52 <- k_p0_deg*x_V_e_p0*m_V_e   # R52 Epitope degradation in APC endosomes
    v53 <- k_P0_on*x_V_e_p0*x_A_e_M1*m_A_e   # R53 Antigen-MHC binding in endosome
    v54 <- k_P0_d1*k_P0_on*x_A_e_M1p0*m_A_e   # R54 Antigen-MHC unbinding in endosome
    v55 <- k_P0_d1*k_P0_on*x_A_s_M1p0*m_A_s   # R55 Antigen-MHC unbinding on APC surface
    v56 <- kout*x_A_e_M1p0*m_A_e   # R56 Antigen-MHC translocation
    v57 <- n_T1_clones*P1_C1*(a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1)))*x_V_T_C1*vol_V_T   # R57 Antigen deposition from dying cancer cells
    v58 <- k_xP1_deg*x_V_LN_P1*m_V_LN   # R58 Free antigen degradation
    v59 <- k_P1_up*x_V_LN_mAPC*x_V_LN_P1*m_V_LN   # R59 Antigen uptake by mature antigen presenting cells
    v60 <- k_P1_up*1*x_V_LN_P1*m_V_e   # R60 Antigen uptake by mature antigen presenting cells
    v61 <- k_P1_deg*x_V_e_P1*m_V_e   # R61 Antigen degradation in APC endosomes
    v62 <- k_p1_deg*x_V_e_p1*m_V_e   # R62 Epitope degradation in APC endosomes
    v63 <- k_P1_on*x_V_e_p1*x_A_e_M1*m_A_e   # R63 Antigen-MHC binding in endosome
    v64 <- k_P1_d1*k_P1_on*x_A_e_M1p1*m_A_e   # R64 Antigen-MHC unbinding in endosome
    v65 <- k_P1_d1*k_P1_on*x_A_s_M1p1*m_A_s   # R65 Antigen-MHC unbinding on APC surface
    v66 <- kout*x_A_e_M1p1*m_A_e   # R66 Antigen-MHC translocation
    v67 <- q_P_nivo*((x_V_C_nivo/gamma_C_nivo)-(x_V_P_nivo/gamma_P_nivo))   # R67 nivo diffusive transport to peripheral compartment
    v68 <- q_T_nivo*((x_V_C_nivo/gamma_C_nivo)-(x_V_T_nivo/gamma_T_nivo))   # R68 nivo diffusive transport to tumor compartment
    v69 <- q_LN_nivo*((x_V_C_nivo/gamma_C_nivo)-(x_V_LN_nivo/gamma_LN_nivo))   # R69 nivo diffusive transport to lymph node compartment
    v70 <- ((q_LD_nivo*x_V_T_nivo)/gamma_T_nivo)*vol_V_T   # R70 nivo convective transport from tumor to lymph node
    v71 <- ((q_LD_nivo*x_V_LN_nivo)/gamma_LN_nivo)*m_V_LN   # R71 nivo convective transport from lymph node to central
    v72 <- k_cl_nivo*x_V_C_nivo   # R72 Drug clearance from central compartment
    v73 <- q_P_durv*((x_V_C_durv/gamma_C_durv)-(x_V_P_durv/gamma_P_durv))   # R73 durv diffusive transport to peripheral compartment
    v74 <- q_T_durv*((x_V_C_durv/gamma_C_durv)-(x_V_T_durv/gamma_T_durv))   # R74 durv diffusive transport to tumor compartment
    v75 <- q_LN_durv*((x_V_C_durv/gamma_C_durv)-(x_V_LN_durv/gamma_LN_durv))   # R75 durv diffusive transport to lymph node compartment
    v76 <- ((q_LD_durv*x_V_T_durv)/gamma_T_durv)*vol_V_T   # R76 durv convective transport from tumor to lymph node
    v77 <- ((q_LD_durv*x_V_LN_durv)/gamma_LN_durv)*m_V_LN   # R77 durv convective transport from lymph node to central
    v78 <- k_cl_durv*x_V_C_durv   # R78 Drug clearance from central compartment
    v79 <- q_P_ipi*((x_V_C_ipi/gamma_C_ipi)-(x_V_P_ipi/gamma_P_ipi))   # R79 ipi diffusive transport to peripheral compartment
    v80 <- q_T_ipi*((x_V_C_ipi/gamma_C_ipi)-(x_V_T_ipi/gamma_T_ipi))   # R80 ipi diffusive transport to tumor compartment
    v81 <- q_LN_ipi*((x_V_C_ipi/gamma_C_ipi)-(x_V_LN_ipi/gamma_LN_ipi))   # R81 ipi diffusive transport to lymph node compartment
    v82 <- ((q_LD_ipi*x_V_T_ipi)/gamma_T_ipi)*vol_V_T   # R82 ipi convective transport from tumor to lymph node
    v83 <- ((q_LD_ipi*x_V_LN_ipi)/gamma_LN_ipi)*m_V_LN   # R83 ipi convective transport from lymph node to central
    v84 <- k_cl_ipi*x_V_C_ipi   # R84 Drug clearance from central compartment
    v85 <- ((kon_PD1_PDL1*x_syn_T_C1_PD1*x_syn_T_C1_PDL1)-(koff_PD1_PDL1*x_syn_T_C1_PD1_PDL1))*m_syn_T_C1   # R85 binding and unbinding of PD1 PDL1 in synapse
    v86 <- ((kon_PD1_PDL2*x_syn_T_C1_PD1*x_syn_T_C1_PDL2)-(koff_PD1_PDL2*x_syn_T_C1_PD1_PDL2))*m_syn_T_C1   # R86 binding and unbinding of PD1 PDL2 in synapse
    v87 <- ((2*kon_PD1_nivo*((x_syn_T_C1_PD1*x_V_T_nivo)/gamma_T_nivo))-(koff_PD1_nivo*x_syn_T_C1_PD1_nivo))*m_syn_T_C1   # R87 binding and unbinding of PD1 to Nivo on T surface in synapse
    v88 <- ((Chi_PD1_nivo*kon_PD1_nivo*x_syn_T_C1_PD1*x_syn_T_C1_PD1_nivo)-(2*koff_PD1_nivo*x_syn_T_C1_PD1_nivo_PD1))*m_syn_T_C1   # R88 binding and unbinding of PD1 to PD1-Nivo on T surface in synapse
    v89 <- ((2*kon_PDL1_durv*((x_syn_T_C1_PDL1*x_V_T_durv)/gamma_T_durv))-(koff_PDL1_durv*x_syn_T_C1_PDL1_durv))*m_syn_T_C1   # R89 binding and unbinding of PDL1 to Durv on C1 surface in synapse
    v90 <- ((Chi_PDL1_durv*kon_PDL1_durv*x_syn_T_C1_PDL1*x_syn_T_C1_PDL1_durv)-(2*koff_PDL1_durv*x_syn_T_C1_PDL1_durv_PDL1))*m_syn_T_C1   # R90 binding and unbinding of PDL1 to PDL1-Durv on C1 surface in synapse
    v91 <- ((2*kon_CD28_CD80*x_syn_T_C1_CD28*x_syn_T_C1_CD80)-(koff_CD28_CD80*x_syn_T_C1_CD28_CD80))*m_syn_T_C1   # R91 binding and unbinding of CD28 and CD80 in synapse
    v92 <- ((kon_CD28_CD80*x_syn_T_C1_CD28*x_syn_T_C1_CD28_CD80)-(2*koff_CD28_CD80*x_syn_T_C1_CD28_CD80_CD28))*m_syn_T_C1   # R92 binding and unbinding of CD28-CD80 and CD28 in synapse
    v93 <- ((kon_CD28_CD86*x_syn_T_C1_CD28*x_syn_T_C1_CD86)-(koff_CD28_CD86*x_syn_T_C1_CD28_CD86))*m_syn_T_C1   # R93 binding and unbinding of CD28 and CD86 in synapse
    v94 <- ((4*kon_CTLA4_CD80*x_syn_T_C1_CTLA4*x_syn_T_C1_CD80)-(koff_CTLA4_CD80*x_syn_T_C1_CD80_CTLA4))*m_syn_T_C1   # R94 binding and unbinding of CTLA4 and CD80 in synapse
    v95 <- ((kon_CTLA4_CD80*x_syn_T_C1_CTLA4*x_syn_T_C1_CD80_CTLA4)-(2*koff_CTLA4_CD80*x_syn_T_C1_CTLA4_CD80_CTLA4))*m_syn_T_C1   # R95 binding and unbinding of CTLA4 and CD80-CTLA4 in synapse
    v96 <- ((kon_CTLA4_CD80*x_syn_T_C1_CD80*x_syn_T_C1_CTLA4_CD80_CTLA4)-(koff_CTLA4_CD80*x_syn_T_C1_CD80_CTLA4_CD80_CTLA4))*m_syn_T_C1   # R96 binding and unbinding of CD80 and CTLA4-CD80-CTLA4 in synapse
    v97 <- ((kon_CTLA4_CD80*x_syn_T_C1_CD80_CTLA4*x_syn_T_C1_CD80)-(2*koff_CTLA4_CD80*x_syn_T_C1_CD80_CTLA4_CD80))*m_syn_T_C1   # R97 binding and unbinding of CD80-CTLA4 and CD80 in synapse
    v98 <- ((kon_CTLA4_CD80*x_syn_T_C1_CTLA4*x_syn_T_C1_CD80_CTLA4_CD80)-(koff_CTLA4_CD80*x_syn_T_C1_CD80_CTLA4_CD80_CTLA4))*m_syn_T_C1   # R98 binding and unbinding of CTLA4 and CD80-CTLA4-CD80 in synapse
    v99 <- ((2*kon_CTLA4_CD86*x_syn_T_C1_CTLA4*x_syn_T_C1_CD86)-(koff_CTLA4_CD86*x_syn_T_C1_CD86_CTLA4))*m_syn_T_C1   # R99 binding and unbinding of CTLA4 and CD86 in synapse
    v100 <- ((kon_CTLA4_CD86*x_syn_T_C1_CD86_CTLA4*x_syn_T_C1_CD86)-(2*koff_CTLA4_CD86*x_syn_T_C1_CD86_CTLA4_CD86))*m_syn_T_C1   # R100 binding and unbinding of CD86-CTLA4 and CD86 in synapse
    v101 <- ((4*kon_CTLA4_ipi*((x_syn_T_C1_CTLA4*x_V_T_ipi)/gamma_T_ipi))-(koff_CTLA4_ipi*x_syn_T_C1_CTLA4_ipi))*m_syn_T_C1   # R101 binding and unbinding of CTLA4 to Ipi on T surface in synapse
    v102 <- ((Chi_CTLA4_ipi*kon_CTLA4_ipi*x_syn_T_C1_CTLA4*x_syn_T_C1_CTLA4_ipi)-(2*koff_CTLA4_ipi*x_syn_T_C1_CTLA4_ipi_CTLA4))*m_syn_T_C1   # R102 binding and unbinding of CTLA4 to Ipi on T surface in synapse
    v103 <- ((2*kon_CD80_PDL1*x_syn_T_C1_CD80*x_syn_T_C1_TPDL1)-(koff_CD80_PDL1*x_syn_T_C1_TPDL1_CD80))*m_syn_T_C1   # R103 binding and unbinding of CD80 and PDL1 in synapse
    v104 <- ((kon_CD80_PDL1*x_syn_T_C1_TPDL1_CD80*x_syn_T_C1_TPDL1)-(2*koff_CD80_PDL1*x_syn_T_C1_TPDL1_CD80_TPDL1))*m_syn_T_C1   # R104 binding and unbinding of PDL1-CD80 and PDL1 in synapse
    v105 <- ((2*kon_PDL1_durv*((x_syn_T_C1_TPDL1*x_V_T_durv)/gamma_T_durv))-(koff_PDL1_durv*x_syn_T_C1_TPDL1_durv))*m_syn_T_C1   # R105 binding and unbinding of PDL1 to Durv on C1 surface in synapse
    v106 <- ((Chi_PDL1_durv*kon_PDL1_durv*x_syn_T_C1_TPDL1*x_syn_T_C1_TPDL1_durv)-(2*koff_PDL1_durv*x_syn_T_C1_TPDL1_durv_TPDL1))*m_syn_T_C1   # R106 binding and unbinding of PDL1 to PDL1-Durv on C1 surface in synapse
    v107 <- ((kon_PD1_PDL1*x_syn_T_APC_PD1*x_syn_T_APC_PDL1)-(koff_PD1_PDL1*x_syn_T_APC_PD1_PDL1))*m_syn_T_APC   # R107 binding and unbinding of PD1 PDL1 in synapse
    v108 <- ((kon_PD1_PDL2*x_syn_T_APC_PD1*x_syn_T_APC_PDL2)-(koff_PD1_PDL2*x_syn_T_APC_PD1_PDL2))*m_syn_T_APC   # R108 binding and unbinding of PD1 PDL2 in synapse
    v109 <- ((2*kon_PD1_nivo*((x_syn_T_APC_PD1*x_V_LN_nivo)/gamma_LN_nivo))-(koff_PD1_nivo*x_syn_T_APC_PD1_nivo))*m_syn_T_APC   # R109 binding and unbinding of PD1 to Nivo on T surface in synapse
    v110 <- ((Chi_PD1_nivo*kon_PD1_nivo*x_syn_T_APC_PD1*x_syn_T_APC_PD1_nivo)-(2*koff_PD1_nivo*x_syn_T_APC_PD1_nivo_PD1))*m_syn_T_APC   # R110 binding and unbinding of PD1 to PD1-Nivo on T surface in synapse
    v111 <- ((2*kon_PDL1_durv*((x_syn_T_APC_PDL1*x_V_LN_durv)/gamma_LN_durv))-(koff_PDL1_durv*x_syn_T_APC_PDL1_durv))*m_syn_T_APC   # R111 binding and unbinding of PDL1 to Durv on APC surface in synapse
    v112 <- ((Chi_PDL1_durv*kon_PDL1_durv*x_syn_T_APC_PDL1*x_syn_T_APC_PDL1_durv)-(2*koff_PDL1_durv*x_syn_T_APC_PDL1_durv_PDL1))*m_syn_T_APC   # R112 binding and unbinding of PDL1 to PDL1-Durv on APC surface in synapse
    v113 <- ((2*kon_CD28_CD80*x_syn_T_APC_CD28*x_syn_T_APC_CD80)-(koff_CD28_CD80*x_syn_T_APC_CD28_CD80))*m_syn_T_APC   # R113 binding and unbinding of CD28 and CD80 in synapse
    v114 <- ((kon_CD28_CD80*x_syn_T_APC_CD28*x_syn_T_APC_CD28_CD80)-(2*koff_CD28_CD80*x_syn_T_APC_CD28_CD80_CD28))*m_syn_T_APC   # R114 binding and unbinding of CD28-CD80 and CD28 in synapse
    v115 <- ((kon_CD28_CD86*x_syn_T_APC_CD28*x_syn_T_APC_CD86)-(koff_CD28_CD86*x_syn_T_APC_CD28_CD86))*m_syn_T_APC   # R115 binding and unbinding of CD28 and CD86 in synapse
    v116 <- ((4*kon_CTLA4_CD80*x_syn_T_APC_CTLA4*x_syn_T_APC_CD80)-(koff_CTLA4_CD80*x_syn_T_APC_CD80_CTLA4))*m_syn_T_APC   # R116 binding and unbinding of CTLA4 and CD80 in synapse
    v117 <- ((kon_CTLA4_CD80*x_syn_T_APC_CTLA4*x_syn_T_APC_CD80_CTLA4)-(2*koff_CTLA4_CD80*x_syn_T_APC_CTLA4_CD80_CTLA4))*m_syn_T_APC   # R117 binding and unbinding of CTLA4 and CD80-CTLA4 in synapse
    v118 <- ((kon_CTLA4_CD80*x_syn_T_APC_CD80*x_syn_T_APC_CTLA4_CD80_CTLA4)-(koff_CTLA4_CD80*x_syn_T_APC_CD80_CTLA4_CD80_CTLA4))*m_syn_T_APC   # R118 binding and unbinding of CD80 and CTLA4-CD80-CTLA4 in synapse
    v119 <- ((kon_CTLA4_CD80*x_syn_T_APC_CD80_CTLA4*x_syn_T_APC_CD80)-(2*koff_CTLA4_CD80*x_syn_T_APC_CD80_CTLA4_CD80))*m_syn_T_APC   # R119 binding and unbinding of CD80-CTLA4 and CD80 in synapse
    v120 <- ((kon_CTLA4_CD80*x_syn_T_APC_CTLA4*x_syn_T_APC_CD80_CTLA4_CD80)-(koff_CTLA4_CD80*x_syn_T_APC_CD80_CTLA4_CD80_CTLA4))*m_syn_T_APC   # R120 binding and unbinding of CTLA4 and CD80-CTLA4-CD80 in synapse
    v121 <- ((2*kon_CTLA4_CD86*x_syn_T_APC_CTLA4*x_syn_T_APC_CD86)-(koff_CTLA4_CD86*x_syn_T_APC_CD86_CTLA4))*m_syn_T_APC   # R121 binding and unbinding of CTLA4 and CD86 in synapse
    v122 <- ((kon_CTLA4_CD86*x_syn_T_APC_CD86_CTLA4*x_syn_T_APC_CD86)-(2*koff_CTLA4_CD86*x_syn_T_APC_CD86_CTLA4_CD86))*m_syn_T_APC   # R122 binding and unbinding of CD86-CTLA4 and CD86 in synapse
    v123 <- ((4*kon_CTLA4_ipi*((x_syn_T_APC_CTLA4*x_V_LN_ipi)/gamma_LN_ipi))-(koff_CTLA4_ipi*x_syn_T_APC_CTLA4_ipi))*m_syn_T_APC   # R123 binding and unbinding of CTLA4 to Ipi on T surface in synapse
    v124 <- ((Chi_CTLA4_ipi*kon_CTLA4_ipi*x_syn_T_APC_CTLA4*x_syn_T_APC_CTLA4_ipi)-(2*koff_CTLA4_ipi*x_syn_T_APC_CTLA4_ipi_CTLA4))*m_syn_T_APC   # R124 binding and unbinding of CTLA4 to Ipi on T surface in synapse
    v125 <- ((2*kon_CD80_PDL1*x_syn_T_APC_CD80*x_syn_T_APC_TPDL1)-(koff_CD80_PDL1*x_syn_T_APC_TPDL1_CD80))*m_syn_T_APC   # R125 binding and unbinding of CD80 and PDL1 in synapse
    v126 <- ((kon_CD80_PDL1*x_syn_T_APC_TPDL1_CD80*x_syn_T_APC_TPDL1)-(2*koff_CD80_PDL1*x_syn_T_APC_TPDL1_CD80_TPDL1))*m_syn_T_APC   # R126 binding and unbinding of PDL1-CD80 and PDL1 in synapse
    v127 <- ((2*kon_PDL1_durv*((x_syn_T_APC_TPDL1*x_V_LN_durv)/gamma_LN_durv))-(koff_PDL1_durv*x_syn_T_APC_TPDL1_durv))*m_syn_T_APC   # R127 binding and unbinding of PDL1 to Durv on APC surface in synapse
    v128 <- ((Chi_PDL1_durv*kon_PDL1_durv*x_syn_T_APC_TPDL1*x_syn_T_APC_TPDL1_durv)-(2*koff_PDL1_durv*x_syn_T_APC_TPDL1_durv_TPDL1))*m_syn_T_APC   # R128 binding and unbinding of PDL1 to PDL1-Durv on APC surface in synapse
    v129 <- (kon_CTLA4_ipi*((x_V_T_Treg_CTLA4*x_V_T_ipi)/gamma_C_ipi))-(koff_CTLA4_ipi*x_V_T_Treg_CTLA4_ipi)   # R129 binding and unbinding of CTLA4 to Ipi on Treg surface in Tumor
    v130 <- ((Chi_CTLA4_ipi*kon_CTLA4_ipi*x_V_T_Treg_CTLA4*x_V_T_Treg_CTLA4_ipi)/A_Tcell)-(koff_CTLA4_ipi*x_V_T_Treg_CTLA4_ipi_CTLA4)   # R130 binding and unbinding of CTLA4 to CTLA4-Ipi on Treg surface in Tumor
    v131 <- (kon_CTLA4_ipi*((x_V_P_Treg_CTLA4*x_V_P_ipi)/gamma_P_ipi))-(koff_CTLA4_ipi*x_V_P_Treg_CTLA4_ipi)   # R131 binding and unbinding of CTLA4 to Ipi on Treg surface in Peripheral compartment
    v132 <- ((Chi_CTLA4_ipi*kon_CTLA4_ipi*x_V_P_Treg_CTLA4*x_V_P_Treg_CTLA4_ipi)/A_Tcell)-(koff_CTLA4_ipi*x_V_P_Treg_CTLA4_ipi_CTLA4)   # R132 binding and unbinding of CTLA4 to CTLA4-Ipi on Treg surface in Peripheral compartment
    v133 <- k_CTLA4_ADCC*x_V_T_T0*H_Treg_T   # R133 CTLA4 ADCC in the tumor compartment
    v134 <- k_CTLA4_ADCC*x_V_P_T0*H_Treg_P   # R134 CTLA4 ADCC in the peripheral compartment
    v135 <- k_rec_MDSC*((MDSC_max*vol_V_T)-x_V_T_MDSC)*(x_V_T_CCL2/(EC50_CCL2_rec+x_V_T_CCL2))   # R135 MDSC recruitment into the tumor compartment
    v136 <- k_brec_MDSC*((MDSC_max*vol_V_T)-x_V_T_MDSC)   # R136 Baseline MDSC Migration into the tumor compartment
    v137 <- kd_MDSC*x_V_T_MDSC   # R137 MDSC death in the tumor compartment
    v138 <- k_deg_CCL2*x_V_T_CCL2*vol_V_T   # R138 CCL2 degradation
    v139 <- k_deg_NO*x_V_T_NO*vol_V_T   # R139 NO degradation
    v140 <- k_deg_ArgI*x_V_T_ArgI*vol_V_T   # R140 ArgI degradation
    v141 <- k_sec_CCL2*x_V_T_C1*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_CCL2)))   # R141 CCL2 secretion by MDSCs and cancer cells
    v142 <- k_sec_NO*x_V_T_MDSC*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_NO)))   # R142 NO secretion from MDSC
    v143 <- k_sec_ArgI*x_V_T_MDSC*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_ArgI)))   # R143 ArgI secretion from MDSC
    v144 <- ((ki_Treg*x_V_T_T0*(1-(x_V_T_T0/(Treg_max*vol_V_T)))*x_V_T_ArgI)/(EC50_ArgI_Treg+x_V_T_ArgI))*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_ArgI)))   # R144 ArgI-induced Treg Expension
    v145 <- k_dose2*x_V_C_Dose2*m_V_C   # R145 Administeration of ENT into depot 2
    v146 <- k_a1_ENT*x_V_C_ENT_Buccal*m_V_C   # R146 ENT Buccal absorption
    v147 <- k_a2_ENT*x_V_C_ENT_GI*m_V_C   # R147 ENT Gastro-intestinal absorption
    v148 <- ((k_cln_ENT*x_V_C_ENT)/(x_V_C_ENT+Kc_ENT))*m_V_C   # R148 ENT Non-linear clearance from central compartment
    v149 <- q_P_ENT*((x_V_C_ENT/gamma_C_ENT)-(x_V_P_ENT/gamma_P_ENT))   # R149 ENT diffusive transport to peripheral compartment
    v150 <- q_T_ENT*((x_V_C_ENT/gamma_C_ENT)-(x_V_T_ENT/gamma_T_ENT))   # R150 ENT diffusive transport to tumor compartment
    v151 <- q_LN_ENT*((x_V_C_ENT/gamma_C_ENT)-(x_V_LN_ENT/gamma_LN_ENT))   # R151 ENT diffusive transport to lymph node compartment
    v152 <- ((q_LD_ENT*x_V_T_ENT)/gamma_T_ENT)*vol_V_T   # R152 ENT convective transport from tumor to lymph node
    v153 <- ((q_LD_ENT*x_V_LN_ENT)/gamma_LN_ENT)*m_V_LN   # R153 ENT convective transport from lymph node to central
    v154 <- k_cl_ENT*x_V_C_ENT*m_V_C   # R154 Drug clearance from central compartment

    # ---- ODEs ----
    d/dt(q_V_C_T0) <- -v10-v14+v15-v16+v17
    d/dt(q_V_C_T1) <- -v26-v32+v33-v34+v35
    d/dt(q_V_C_nivo) <- (-v67-v68-v69+v71-v72)/nmol_to_molecule
    d/dt(q_V_C_durv) <- (-v73-v74-v75+v77-v78)/nmol_to_molecule
    d/dt(q_V_C_ipi) <- (-v79-v80-v81+v83-v84)/nmol_to_molecule
    d/dt(q_V_C_ENT) <- (v146+v147-v148-v149-v150-v151+v153-v154)/nmol_to_molecule
    d/dt(q_V_C_ENT_Buccal) <- (-v146)/nmol_to_molecule
    d/dt(q_V_C_ENT_GI) <- (v145-v147)/nmol_to_molecule
    d/dt(q_V_C_Dose2) <- (-v145)/nmol_to_molecule
    d/dt(q_V_P_T0) <- -v11+v14-v15-v134
    d/dt(q_V_P_T1) <- -v27+v32-v33
    d/dt(q_V_P_nivo) <- (v67)/nmol_to_molecule
    d/dt(q_V_P_durv) <- (v73)/nmol_to_molecule
    d/dt(q_V_P_ipi) <- (v79)/nmol_to_molecule
    d/dt(q_V_P_Treg_CTLA4) <- -v131-v132
    d/dt(q_V_P_Treg_CTLA4_ipi) <- v131-v132
    d/dt(q_V_P_Treg_CTLA4_ipi_CTLA4) <- v132
    d/dt(q_V_P_ENT) <- (v149)/nmol_to_molecule
    d/dt(q_V_T_C_x) <- -v1+v4+v37
    d/dt(q_V_T_T_exh) <- -v2+v12+v28+v30+v31
    d/dt(q_V_T_C1) <- v3-v4-v37
    d/dt(q_V_T_T0) <- -v12+v16-v133+v144
    d/dt(q_V_T_T1) <- -v28-v30-v31+v34
    d/dt(q_V_T_APC) <- v38-v40
    d/dt(q_V_T_mAPC) <- v40-v41-v42
    d/dt(q_V_T_c) <- v44+v45
    d/dt(q_V_T_nivo) <- (v68-v70)/nmol_to_molecule
    d/dt(q_V_T_durv) <- (v74-v76)/nmol_to_molecule
    d/dt(q_V_T_ipi) <- (v80-v82)/nmol_to_molecule
    d/dt(q_V_T_Treg_CTLA4) <- -v129-v130
    d/dt(q_V_T_Treg_CTLA4_ipi) <- v129-v130
    d/dt(q_V_T_Treg_CTLA4_ipi_CTLA4) <- v130
    d/dt(q_V_T_MDSC) <- v135+v136-v137
    d/dt(q_V_T_CCL2) <- -v138+v141
    d/dt(q_V_T_NO) <- -v139+v142
    d/dt(q_V_T_ArgI) <- -v140+v143
    d/dt(q_V_T_ENT) <- (v150-v152)/nmol_to_molecule
    d/dt(q_V_LN_nT0) <- v5-v6-v7
    d/dt(q_V_LN_aT0) <- v7-v8
    d/dt(q_V_LN_T0) <- v9-v13-v17
    d/dt(q_V_LN_IL2) <- -v18-v19-v20+v36
    d/dt(q_V_LN_nT1) <- v21-v22-v23
    d/dt(q_V_LN_aT1) <- v23-v24
    d/dt(q_V_LN_T1) <- v25-v29-v35
    d/dt(q_V_LN_APC) <- v39
    d/dt(q_V_LN_mAPC) <- v41-v43
    d/dt(q_V_LN_P0) <- v47-v48-v49
    d/dt(q_V_LN_P1) <- v57-v58-v59
    d/dt(q_V_LN_nivo) <- (v69+v70-v71)/nmol_to_molecule
    d/dt(q_V_LN_durv) <- (v75+v76-v77)/nmol_to_molecule
    d/dt(q_V_LN_ipi) <- (v81+v82-v83)/nmol_to_molecule
    d/dt(q_V_LN_ENT) <- (v151+v152-v153)/nmol_to_molecule
    d/dt(q_V_e_P0) <- v50-v51
    d/dt(q_V_e_p0) <- v51-v52-v53+v54
    d/dt(q_V_e_P1) <- v60-v61
    d/dt(q_V_e_p1) <- v61-v62-v63+v64
    d/dt(q_A_e_M1) <- -v46-v53+v54-v63+v64
    d/dt(q_A_e_M1p0) <- v53-v54-v56
    d/dt(q_A_e_M1p1) <- v63-v64-v66
    d/dt(q_A_s_M1) <- v46+v55+v65
    d/dt(q_A_s_M1p0) <- -v55+v56
    d/dt(q_A_s_M1p1) <- -v65+v66
    d/dt(q_syn_T_C1_PD1_PDL1) <- v85
    d/dt(q_syn_T_C1_PD1_PDL2) <- v86
    d/dt(q_syn_T_C1_PD1) <- -v85-v86-v87-v88
    d/dt(q_syn_T_C1_PDL1) <- -v85-v89-v90
    d/dt(q_syn_T_C1_PDL2) <- -v86
    d/dt(q_syn_T_C1_PD1_nivo) <- v87-v88
    d/dt(q_syn_T_C1_PD1_nivo_PD1) <- v88
    d/dt(q_syn_T_C1_PDL1_durv) <- v89-v90
    d/dt(q_syn_T_C1_PDL1_durv_PDL1) <- v90
    d/dt(q_syn_T_C1_TPDL1) <- -v103-v104-v105-v106
    d/dt(q_syn_T_C1_TPDL1_durv) <- v105-v106
    d/dt(q_syn_T_C1_TPDL1_durv_TPDL1) <- v106
    d/dt(q_syn_T_C1_CD28_CD80) <- v91-v92
    d/dt(q_syn_T_C1_CD28_CD80_CD28) <- v92
    d/dt(q_syn_T_C1_CD28_CD86) <- v93
    d/dt(q_syn_T_C1_CD80_CTLA4) <- v94-v95-v97
    d/dt(q_syn_T_C1_CD80_CTLA4_CD80) <- v97-v98
    d/dt(q_syn_T_C1_CTLA4_CD80_CTLA4) <- v95-v96
    d/dt(q_syn_T_C1_CD80_CTLA4_CD80_CTLA4) <- v96+v98
    d/dt(q_syn_T_C1_CD86_CTLA4) <- v99-v100
    d/dt(q_syn_T_C1_CD86_CTLA4_CD86) <- v100
    d/dt(q_syn_T_C1_TPDL1_CD80) <- v103-v104
    d/dt(q_syn_T_C1_TPDL1_CD80_TPDL1) <- v104
    d/dt(q_syn_T_C1_CD28) <- -v91-v92-v93
    d/dt(q_syn_T_C1_CTLA4) <- -v94-v95-v98-v99-v101-v102
    d/dt(q_syn_T_C1_CD80) <- -v91-v94-v96-v97-v103
    d/dt(q_syn_T_C1_CD86) <- -v93-v99-v100
    d/dt(q_syn_T_C1_CTLA4_ipi) <- v101-v102
    d/dt(q_syn_T_C1_CTLA4_ipi_CTLA4) <- v102
    d/dt(q_syn_T_APC_PD1_PDL1) <- v107
    d/dt(q_syn_T_APC_PD1_PDL2) <- v108
    d/dt(q_syn_T_APC_PD1) <- -v107-v108-v109-v110
    d/dt(q_syn_T_APC_PDL1) <- -v107-v111-v112
    d/dt(q_syn_T_APC_PDL2) <- -v108
    d/dt(q_syn_T_APC_PD1_nivo) <- v109-v110
    d/dt(q_syn_T_APC_PD1_nivo_PD1) <- v110
    d/dt(q_syn_T_APC_PDL1_durv) <- v111-v112
    d/dt(q_syn_T_APC_PDL1_durv_PDL1) <- v112
    d/dt(q_syn_T_APC_TPDL1) <- -v125-v126-v127-v128
    d/dt(q_syn_T_APC_TPDL1_durv) <- v127-v128
    d/dt(q_syn_T_APC_TPDL1_durv_TPDL1) <- v128
    d/dt(q_syn_T_APC_CD28_CD80) <- v113-v114
    d/dt(q_syn_T_APC_CD28_CD80_CD28) <- v114
    d/dt(q_syn_T_APC_CD28_CD86) <- v115
    d/dt(q_syn_T_APC_CD80_CTLA4) <- v116-v117-v119
    d/dt(q_syn_T_APC_CD80_CTLA4_CD80) <- v119-v120
    d/dt(q_syn_T_APC_CTLA4_CD80_CTLA4) <- v117-v118
    d/dt(q_syn_T_APC_CD80_CTLA4_CD80_CTLA4) <- v118+v120
    d/dt(q_syn_T_APC_CD86_CTLA4) <- v121-v122
    d/dt(q_syn_T_APC_CD86_CTLA4_CD86) <- v122
    d/dt(q_syn_T_APC_TPDL1_CD80) <- v125-v126
    d/dt(q_syn_T_APC_TPDL1_CD80_TPDL1) <- v126
    d/dt(q_syn_T_APC_CD28) <- -v113-v114-v115
    d/dt(q_syn_T_APC_CTLA4) <- -v116-v117-v120-v121-v123-v124
    d/dt(q_syn_T_APC_CD80) <- -v113-v116-v118-v119-v125
    d/dt(q_syn_T_APC_CD86) <- -v115-v121-v122
    d/dt(q_syn_T_APC_CTLA4_ipi) <- v123-v124
    d/dt(q_syn_T_APC_CTLA4_ipi_CTLA4) <- v124

    # ---- initial conditions (S1 Table, sheet 'Table S5. Species', and the SBML
    # deposit's initial assignments) ----
    q_V_P_Treg_CTLA4(0) <- Treg_CTLA4_tot
    q_V_T_C1(0) <- 1   # S1 Table S5 Species (V_T.C1 = 1 cell)
    q_V_T_Treg_CTLA4(0) <- Treg_CTLA4_tot
    q_V_LN_IL2(0) <- 670.052061644   # S1 Table S5 Species prints V_LN.IL2 = 1e-24 'molarity'; the SBML deposit holds 1e-24 mol per mm^3 of V_LN, i.e. 1e-18 mol/L
    q_V_LN_P0(0) <- 670.052061644   # S1 Table S5 Species prints V_LN.P0 = 1e-24 'molarity'; the SBML deposit holds 1e-24 mol per mm^3 of V_LN, i.e. 1e-18 mol/L
    q_V_LN_P1(0) <- 670.052061644   # S1 Table S5 Species prints V_LN.P1 = 1e-24 'molarity'; the SBML deposit holds 1e-24 mol per mm^3 of V_LN, i.e. 1e-18 mol/L
    q_V_e_P0(0) <- 2.408856796e-10   # S1 Table S5 Species (V_e.P0 = 1e-18 molarity)
    q_V_e_p0(0) <- 2.408856796e-10   # S1 Table S5 Species (V_e.p0 = 1e-18 molarity)
    q_V_e_P1(0) <- 2.408856796e-10   # S1 Table S5 Species (V_e.P1 = 1e-18 molarity)
    q_V_e_p1(0) <- 2.408856796e-10   # S1 Table S5 Species (V_e.p1 = 1e-18 molarity)
    q_A_e_M1(0) <- 32786.8852459   # S1 Table S5 Species (A_e.M1 = 2185.79235 molecule/micrometer^2)
    q_A_e_M1p0(0) <- 1.5e-05   # S1 Table S5 Species (A_e.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_e_M1p1(0) <- 1.5e-05   # S1 Table S5 Species (A_e.M1p1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1(0) <- 0.0009   # S1 Table S5 Species (A_s.M1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p0(0) <- 0.0009   # S1 Table S5 Species (A_s.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p1(0) <- 0.0009   # S1 Table S5 Species (A_s.M1p1 = 1e-06 molecule/micrometer^2)
    q_syn_T_C1_PD1(0) <- (T_PD1_total/A_Tcell)*m_syn_T_C1
    q_syn_T_C1_PDL1(0) <- (C1_PDL1_total/A_cell)*m_syn_T_C1
    q_syn_T_C1_PDL2(0) <- (C1_PDL2_total/A_cell)*m_syn_T_C1
    q_syn_T_C1_TPDL1(0) <- (T_PDL1_total/A_Tcell)*m_syn_T_C1
    q_syn_T_C1_CD28(0) <- (T_CD28_total/A_Tcell)*m_syn_T_C1
    q_syn_T_C1_CTLA4(0) <- (T_CTLA4_syn/A_syn)*m_syn_T_C1
    q_syn_T_C1_CD80(0) <- (C1_CD80_total/A_cell)*m_syn_T_C1
    q_syn_T_C1_CD86(0) <- (C1_CD86_total/A_cell)*m_syn_T_C1
    q_syn_T_APC_PD1(0) <- (T_PD1_total/A_Tcell)*m_syn_T_APC
    q_syn_T_APC_PDL1(0) <- (APC_PDL1_total/A_APC)*m_syn_T_APC
    q_syn_T_APC_PDL2(0) <- (APC_PDL2_total/A_APC)*m_syn_T_APC
    q_syn_T_APC_TPDL1(0) <- (T_PDL1_total/A_Tcell)*m_syn_T_APC
    q_syn_T_APC_CD28(0) <- (T_CD28_total/A_Tcell)*m_syn_T_APC
    q_syn_T_APC_CTLA4(0) <- (T_CTLA4_syn/A_syn)*m_syn_T_APC
    q_syn_T_APC_CD80(0) <- (APC_CD80_total/A_APC)*m_syn_T_APC
    q_syn_T_APC_CD86(0) <- (APC_CD86_total/A_APC)*m_syn_T_APC

    # ---- entinostat dosing: dose the same amount into both depots. A fraction
    # F_ENT_buccal enters the buccal depot as a zero-order input over durP; the
    # remainder enters Dose2 as a bolus after the lag lagP (Wang 2020 Figure 2A). ----
    f(q_V_C_ENT_Buccal) <- F_ENT_buccal
    dur(q_V_C_ENT_Buccal) <- durP
    f(q_V_C_Dose2) <- 1 - F_ENT_buccal
    alag(q_V_C_Dose2) <- lagP

    # ---- outputs ----
    # tumour diameter [cm] from the total tumour volume assuming a sphere
    tumour_diameter <- (6*vol_V_T*1000/pi)^(1/3)
    # central (plasma) drug concentrations [nmol/L]
    Cc_nivolumab <- x_V_C_nivo/nmol_to_molecule
    Cc_ipilimumab <- x_V_C_ipi/nmol_to_molecule
    Cc_entinostat <- x_V_C_ENT/nmol_to_molecule
    # tumour-infiltrating cell densities [cells/mL of tumour]
    teff_density <- x_V_T_T1/(vol_V_T*1000)
    treg_density <- x_V_T_T0/(vol_V_T*1000)
    mdsc_density <- x_V_T_MDSC/(vol_V_T*1000)
  })
}
