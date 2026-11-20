Yamada_2025_oncolyticVirus_human_qsp <- function() {
  description <- paste(
    "QSP. Clinical quantitative systems pharmacology model of an intratumorally",
    "injected oncolytic vaccinia virus carrying human IL-7 and IL-12 genes",
    "(hIL-7/hIL-12-VV) in combination with pembrolizumab (anti-PD-1). Built by",
    "extending the HER2-negative-breast-cancer immune-oncology platform of Wang",
    "et al. 2020 (Front Bioeng Biotechnol 8:141) with an oncolytic-virus",
    "mechanism of action taken from the authors' preclinical mouse QSP model",
    "(see Yamada_2025_oncolyticVirus_mouse_qsp). Four physiological compartments",
    "(central, peripheral, tumour, tumour-draining lymph node) plus APC endosomal",
    "and surface compartments and two immunological-synapse compartments (T cell",
    ": cancer cell and T cell : APC). Covers logistic cancer-cell growth from a",
    "single cell, naive regulatory and cytotoxic T cell priming in the lymph node",
    "(kinetic-proofreading TCR signal, CD28 co-stimulation and IL-2 driven",
    "division number), T cell trafficking, exhaustion and Treg suppression, APC",
    "maturation and migration, self- and neo-antigen release, uptake, endosomal",
    "processing and MHC presentation, PD-1/PD-L1/PD-L2 and CTLA-4/CD28/CD80/CD86",
    "synapse binding with cross-arm antibody binding, an MDSC module (CCL2-driven",
    "recruitment, arginase-I and nitric-oxide suppression of Teff killing,",
    "arginase-I driven Treg expansion) and anti-CTLA-4 ADCC depletion of Tregs.",
    "The oncolytic-virus layer splits the tumour cancer-cell pool into uninfected",
    "(Cuit) and infected (Cit) cells, adds an intratumoral virus state (Virus_T)",
    "with production proportional to infected-cell lysis and pembrolizumab-enhanced",
    "clearance, infection of uninfected cells, and IL-7 / IL-12 immunostimulation",
    "of naive-T-cell influx and differentiation driven by infected-cell count.",
    "The anti-PD-L1 antibody durvalumab, anti-CTLA-4 ipilimumab and the",
    "entinostat PK/PD module are carried over from the Wang 2020 platform but are",
    "not dosed in this paper.",
    "Deterministic mechanism model: the authors generated virtual patients by",
    "Latin hypercube sampling of parameter ranges rather than by fitting IIV or",
    "residual error, so no etas and no error model are encoded.",
    sep = " "
  )
  reference <- paste(
    "Yamada A, Tanaka S, Nakamura S, Asano N, Nomura T, Takahashi M, Caunt MP,",
    "Boyd PL, Brennan FA, Ortega F, Markham H, Orrell D, Chakrabarti C, Kato K.",
    "A Multiple-Model-Informed Drug-Development Approach for Optimal Regimen",
    "Selection of an Oncolytic Virus in Combination With Pembrolizumab.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):572-583.",
    "doi:10.1002/psp4.13297.",
    "Immune-oncology platform and inherited parameters from Wang H, Sove RJ,",
    "Jafarnejad M, et al. Conducting a Virtual Clinical Trial in HER2-Negative",
    "Breast Cancer Using a Quantitative Systems Pharmacology Model With an",
    "Epigenetic Modulator and Immune Checkpoint Inhibitors. Front Bioeng",
    "Biotechnol. 2020;8:141. doi:10.3389/fbioe.2020.00141.",
    sep = " "
  )
  vignette <- "Yamada_2025_oncolytic_virus_pembrolizumab"

  units <- list(
    time = "day",
    dosing = "nmol (antibodies, entinostat); pfu (virus, held as molecule-equivalents with 1 pfu = 1 molecule)",
    concentration = "nmol/L"
  )
  dosing <- c("q_V_T_Virus_T", "q_V_C_ENT_Buccal", "q_V_C_Dose2", "q_V_C_pembro")

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
    q_V_C_pembro = list(
      analyte = "Concentration of pembro in central compartment",
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
    q_V_P_pembro = list(
      analyte = "Concentration of pembro in peripheral compartment",
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
    q_V_T_Cuit = list(
      analyte = "Number of uninfected cancer cells in tumour",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_Cit = list(
      analyte = "Number of virus-infected cancer cells in tumour",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_Virus_T = list(
      analyte = "Intratumoral oncolytic-virus amount",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_C2 = list(
      analyte = "Parallel no-virus cancer-cell tracker in tumour",
      units = "cell",
      specimen = "tumor",
      verified = FALSE
    ),
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
    q_V_T_pembro = list(
      analyte = "Concentration of pembro in tumour compartment",
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
    q_V_LN_pembro = list(
      analyte = "Concentration of pembro in lymph node compartment",
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
    q_syn_T_C1_PD1_pembro = list(
      analyte = "concentration of PD1-pembrolizumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T_C1_PD1_pembro_PD1 = list(
      analyte = "concentration of PD1-pembrolizumab-PD1 complex",
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
    q_syn_T_APC_PD1_pembro = list(
      analyte = "concentration of PD1-pembrolizumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T_APC_PD1_pembro_PD1 = list(
      analyte = "concentration of PD1-pembrolizumab-PD1 complex",
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
    species = "human (in silico virtual cohort)",
    n_subjects = 100L,
    disease_state = "solid tumour amenable to intratumoral injection",
    dose_range = paste(
      "hIL-7/hIL-12-VV 1e7, 1e8 or 5e8 pfu/mL intratumorally every 2 weeks;",
      "pembrolizumab 400 mg every 6 weeks (combination arms starting on day 1,",
      "8, 15 or 29).",
      sep = " "
    ),
    notes = paste(
      "Prospective virtual clinical trial. 125 virtual-patient parameter sets",
      "were generated by Latin hypercube sampling of the ranges in Supplementary",
      "Table S6 (inherited from Wang 2020); each simulation started from 1e6",
      "cancer cells and grew untreated until it reached a preselected",
      "baseline tumour diameter (uniform 1.0-4.5 cm) within an 8000-day run-in,",
      "yielding 100 evaluable patients. Node values at the baseline diameter were",
      "used as per-patient initial conditions and therapy was simulated for up to",
      "400 days. The viral infection rate constant was adjusted downward relative",
      "to the mouse QSP to keep the single-agent objective response rate below",
      "50%.",
      sep = " "
    )
  )

  # Mechanistic QSP states that are not canonical nlmixr2lib compartments.
  paper_specific_compartments <- c(
    "q_V_C_T0",
    "q_V_C_T1",
    "q_V_C_pembro",
    "q_V_C_durv",
    "q_V_C_ipi",
    "q_V_C_ENT",
    "q_V_C_ENT_Buccal",
    "q_V_C_ENT_GI",
    "q_V_C_Dose2",
    "q_V_P_T0",
    "q_V_P_T1",
    "q_V_P_pembro",
    "q_V_P_durv",
    "q_V_P_ipi",
    "q_V_P_Treg_CTLA4",
    "q_V_P_Treg_CTLA4_ipi",
    "q_V_P_Treg_CTLA4_ipi_CTLA4",
    "q_V_P_ENT",
    "q_V_T_C_x",
    "q_V_T_T_exh",
    "q_V_T_Cuit",
    "q_V_T_Cit",
    "q_V_T_Virus_T",
    "q_V_T_C2",
    "q_V_T_T0",
    "q_V_T_T1",
    "q_V_T_APC",
    "q_V_T_mAPC",
    "q_V_T_c",
    "q_V_T_pembro",
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
    "q_V_LN_pembro",
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
    "q_syn_T_C1_PD1_pembro",
    "q_syn_T_C1_PD1_pembro_PD1",
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
    "q_syn_T_APC_PD1_pembro",
    "q_syn_T_APC_PD1_pembro_PD1",
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

  ini({
    # ---- compartment capacities (Yamada Table S5; Wang 2020 Table S1) ----
    vol_V_C <- fixed(6); label("Capacity of compartment V_C [L]")  # Yamada S5 / Wang 2020 S1 (V_C = 6 litre)
    vol_V_P <- fixed(61.321); label("Capacity of compartment V_P [L]")  # Yamada S5 / Wang 2020 S1 (V_P = 61.321 litre)
    vol_V_LN <- fixed(0.00111264739815); label("Capacity of compartment V_LN [L]")  # Yamada S5 / Wang 2020 S1 (V_LN = 1112.647398 millimeter^3)
    vol_V_e <- fixed(4e-16); label("Capacity of compartment V_e [L]")  # Yamada S5 / Wang 2020 S1 (V_e = 4e-16 litre)
    vol_A_e <- fixed(1.5e-09); label("Capacity of compartment A_e [dm^2]")  # Yamada S5 / Wang 2020 S1 (A_e = 15 micrometer^2)
    vol_A_s <- fixed(9e-08); label("Capacity of compartment A_s [dm^2]")  # Yamada S5 / Wang 2020 S1 (A_s = 900 micrometer^2)
    vol_syn_T_C1 <- fixed(3.78e-09); label("Capacity of compartment syn_T_C1 [dm^2]")  # Yamada S5 / Wang 2020 S1 (syn_T_C1 = 37.8 micrometer^2)
    vol_syn_T_APC <- fixed(3.78e-09); label("Capacity of compartment syn_T_APC [dm^2]")  # Yamada S5 / Wang 2020 S1 (syn_T_APC = 37.8 micrometer^2)

    # ---- model parameters (Yamada Tables S5-S6; Wang 2020 Table S2), converted to the model unit system ----
    k_cell_clear <- fixed(0.1); label("Rate of dead cell clearance from tumour compartment [1/day]")  # Yamada S5 / Wang 2020 S2 (k_cell_clear = 0.1 1/day)
    vol_cell <- fixed(2.57244078451e-12); label("Average volume of cancer cell calculated based on cancer cell diameter [L/count]")  # Yamada S5 / Wang 2020 S2 (vol_cell = 2572.4407845144424 micrometer^3/cell)
    vol_Tcell <- fixed(1.750157098e-13); label("Average volume of T cells calculated based on the average T cell diameter [L/count]")  # Yamada S5 / Wang 2020 S2 (vol_Tcell = 175.01570979953922 micrometer^3/cell)
    V_Tmin <- fixed(1e-09); label("Cancer-Free Tumour compartment volume [L]")  # Yamada S5 / Wang 2020 S2 (V_Tmin = 1e-06 milliliter)
    k_C1_growth <- fixed(0.0075); label("Cancer cell growth rate. Estimated by tumor doubling time: 241?+/-?166 days for ER-positive BC, and ... [1/day]")  # Yamada S6 (k_C1_growth = 0.0075 [0.0075, 1] 1/day, lognormal mean/SD for the virtual population; Wang 2020 S2 used 0.00673)
    C_max <- fixed(862890782186); label("Cancer cell capacity calculated from tumour maximum diameter and cancer cell density [count]")  # Yamada S5 / Wang 2020 S2 (C_max = 862890782185.9965 cell)
    k_C1_death <- fixed(0.0001); label("Cancer cell death rate from innate immune cells [1/day]")  # Yamada S5 / Wang 2020 S2 (k_C1_death = 0.0001 [1e-5 - 1e-3] 1/day)
    k_C1_therapy <- fixed(0); label("Rate of C1 killing by therapy [1/day]")  # Wang 2020 S2, not listed in Yamada S5 (k_C1_therapy = 0 1/day)
    n_T0_clones <- fixed(100); label("Number of T cell clones [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_T0_clones = 100 dimensionless)
    Q_T0_in <- fixed(346.455379492); label("Rate of naive T cell transport into the lLN calculated based on naive T cell entry rate estimated based on ... [count/day]")  # Yamada S5 / Wang 2020 S2 (Q_T0_in = 346.45537949178816 cell/day)
    Q_T0_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN [1/day]")  # Yamada S5 / Wang 2020 S2 (Q_T0_out = 1.13 1/day)
    k_T0_act <- fixed(20); label("T0 activation rate calculated based on the rate of T cell activation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T0_act = 20 1/day)
    k_T0_pro <- fixed(1); label("T0 proliferation rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T0_pro = 1 1/day)
    k_T0_death <- fixed(0.01); label("T0 death rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T0_death = 0.01 1/day)
    q_T0_P_in <- fixed(295.812504); label("rate of T0 tranport into the peripheral compartment calculated based on T cell transmigration rate and T cell ... [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T0_P_in = 0.20542535 1/minute)
    q_T0_P_out <- fixed(0.015); label("rate of T0 tranport out of the peripheral compartment [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T0_P_out = 0.015 1/day)
    q_T0_T_in <- fixed(4.824); label("rate of T0 tranport into the tumour compartment calculated based on T cell transmigration rate and T cell ... [1/(L*day)]")  # Yamada S5 / Wang 2020 S2 (q_T0_T_in = 3.35e-06 1/(centimeter^3*minute))
    q_T0_LN_out <- fixed(24); label("rate of T0 tranport out of the LN compartment [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T0_LN_out = 24 1/day)
    k_IL2_deg <- fixed(288); label("rate of IL2 degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_IL2_deg = 0.2 1/minute)
    k_IL2_cons <- fixed(86718844656); label("rate of IL2 consumption by T cells [1/day]")  # Yamada S5 / Wang 2020 S2 (k_IL2_cons = 6e-06 nanomole/cell/hour)
    k_IL2_sec <- fixed(433594223280); label("rate of IL2 secretion from T cells [1/day]")  # Yamada S5 / Wang 2020 S2 (k_IL2_sec = 3e-05 nanomole/cell/hour)
    IL2_50 <- fixed(1.9270854368e+14); label("T cell activation half-maximal IL2 concentration [count/L]")  # Yamada S5 / Wang 2020 S2 (IL2_50 = 0.32 nanomolarity)
    IL2_50_Treg <- fixed(1.9270854368e+12); label("Treg activation half-maximal IL2 concentration [count/L]")  # Yamada S5 / Wang 2020 S2 (IL2_50_Treg = 0.0032 nanomolarity)
    N0 <- fixed(2); label("numer of activated T cell generation by TCR signaling only [dimensionless]")  # Yamada S5 / Wang 2020 S2 (N0 = 2 dimensionless)
    N_costim <- fixed(3); label("numer of activated T cell generation by co-stimulatory signaling only [dimensionless]")  # Yamada S5 / Wang 2020 S2 (N_costim = 3 [2 - 5] dimensionless)
    N_IL2 <- fixed(11); label("maximum number of activated T cell generations due to IL2 [dimensionless]")  # Yamada S5 / Wang 2020 S2 (N_IL2 = 11 dimensionless)
    k_Treg <- fixed(1); label("Rate of T cell death by Tregs [1/day]")  # Yamada S5 / Wang 2020 S2 (k_Treg = 1 [0.1 - 1] 1/day); SBML deposit carries 2
    n_T1_clones <- fixed(100); label("Number of T cell clones [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_T1_clones = 100 [4 - 1.84e3] dimensionless)
    Q_T1_in <- fixed(212.6058677); label("Rate of naive T cell transport into the lLN calculated based on naive T cell entry rate estimated based on ... [count/day]")  # Yamada S5 / Wang 2020 S2 (Q_T1_in = 212.60586769986332 cell/day)
    Q_T1_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN [1/day]")  # Yamada S5 / Wang 2020 S2 (Q_T1_out = 1.13 1/day)
    k_T1_act <- fixed(20); label("T1 activation rate calculated based on the rate of T cell activation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T1_act = 20 1/day)
    k_T1_pro <- fixed(1); label("T1 proliferation rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T1_pro = 1 1/day)
    k_T1_death <- fixed(0.01); label("T1 death rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T1_death = 0.01 1/day)
    q_T1_P_in <- fixed(295.812504); label("rate of T1 tranport into the peripheral compartment calculated based on T cell transmigration rate and T cell ... [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T1_P_in = 0.20542535 1/minute)
    q_T1_P_out <- fixed(1); label("rate of T1 tranport out of the peripheral compartment [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T1_P_out = 1 1/day)
    q_T1_T_in <- fixed(4.824); label("rate of T1 tranport into the tumour compartment calculated based on T cell transmigration rate and T cell ... [1/(L*day)]")  # Yamada S5 / Wang 2020 S2 (q_T1_T_in = 3.35e-06 1/(centimeter^3*minute))
    q_T1_LN_out <- fixed(24); label("rate of T1 tranport out of the LN compartment [1/day]")  # Yamada S5 / Wang 2020 S2 (q_T1_LN_out = 24 1/day)
    k_T1 <- fixed(0.1); label("Rate of T cell exhaustion by cancer cells [1/day]")  # Yamada S5 / Wang 2020 S2 (k_T1 = 0.1 [0.05 - 0.5] 1/day)
    k_C_T1 <- fixed(2); label("Rate of cancer cell death by T cells [1/day]")  # Yamada S5 / Wang 2020 S2 (k_C_T1 = 2 [2 - 10] 1/day)
    k_APC_mat <- fixed(1.5); label("Maximum rate of APC maturation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_APC_mat = 1.5 1/day)
    k_APC_mig <- fixed(4); label("Rate of APC migration [1/day]")  # Yamada S5 / Wang 2020 S2 (k_APC_mig = 4 1/day)
    k_APC_death <- fixed(0.01); label("Rate of APC death [1/day]")  # Yamada S5 / Wang 2020 S2 (k_APC_death = 0.01 1/day)
    k_mAPC_death <- fixed(0.02); label("Rate of mAPC death [1/day]")  # Yamada S5 / Wang 2020 S2 (k_mAPC_death = 0.02 1/day)
    APC0_T <- fixed(400000000); label("Steady state density of CD103+ APCs in the tumor, calculated to fit basal APC equal to 26000/cm^3 from [count/L]")  # Yamada S5 / Wang 2020 S2 (APC0_T = 400000 cell/milliliter)
    APC0_LN <- fixed(1200000000); label("Steady state density of all APC subtypes in the LN- divided by 10 to only take CD103+ based on [count/L]")  # Yamada S5 / Wang 2020 S2 (APC0_LN = 1200000 cell/milliliter)
    k_c <- fixed(2); label("Cytokine rate constant [1/day]")  # Yamada S5 / Wang 2020 S2 (k_c = 2 1/day)
    c0 <- fixed(6.02214199e+14); label("Baseline cytokine concentration [count/L]")  # Yamada S5 / Wang 2020 S2 (c0 = 1e-09 molarity)
    c50 <- fixed(6.02214199e+14); label("Cytokine concentration for half-maximal APC maturation [count/L]")  # Yamada S5 / Wang 2020 S2 (c50 = 1e-09 molarity)
    DAMPs <- fixed(8069670266.6); label("Concentration of cytokines released per dying cancer cell [dimensionless]")  # Yamada S5 / Wang 2020 S2 (DAMPs = 1.34e-14 mole/cell)
    n_sites_APC <- fixed(10); label("Maxium number of T Cells an APC can interact with [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_sites_APC = 10 dimensionless)
    kin <- fixed(14.4); label("Rate of MHC internalization [1/day]")  # Yamada S5 / Wang 2020 S2 (kin = 14.4 1/day)
    kout <- fixed(28.8); label("Rate of MHC externalization [1/day]")  # Yamada S5 / Wang 2020 S2 (kout = 28.8 1/day)
    k_P0_up <- fixed(14.4); label("Rate of antigen uptake by APCs [1/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_P0_up = 14.4 1/day/cell)
    k_xP0_deg <- fixed(2); label("Rate of extracellular antigen degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_xP0_deg = 2 1/day)
    k_P0_deg <- fixed(17.28); label("Rate of endosomal antigen degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_P0_deg = 17.28 1/day)
    k_p0_deg <- fixed(144); label("Rate of endosomal epitope degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_p0_deg = 144 1/day)
    k_P0_on <- fixed(2.39117576834e-19); label("Rate of antigen-MHC binding [L/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_P0_on = 144000 1/day/molarity)
    k_P0_d1 <- fixed(6.02214199e+16); label("Antigen-MHC kd [count/L]")  # Yamada S5 / Wang 2020 S2 (k_P0_d1 = 1e-07 molarity)
    p0_50 <- fixed(264550.26455); label("calculated based on the number of molecules for half-maximal T cell activation and synapse surface area [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (p0_50 = 2.6455026455026456e-05 molecule/micrometer^2)
    P0_C1 <- fixed(6.02214199e+15); label("Concentration of P0 in C1 [1/L]")  # Yamada S5 / Wang 2020 S2 (P0_C1 = 1e-08 molarity/cell)
    A_syn <- fixed(3.78e-09); label("Surface area of the synapse [dm^2]")  # Yamada S5 / Wang 2020 S2 (A_syn = 37.8 micrometer^2)
    A_Tcell <- fixed(1.5131041193e-08); label("Surface area of the T cell calculated based on the average T cell diameter [dm^2]")  # Yamada S5 / Wang 2020 S2 (A_Tcell = 151.31041193043737 micrometer^2)
    A_cell <- fixed(9.07920276887e-08); label("Surface area of the Cancer cell calculated based on the average Cancer cell diameter [dm^2]")  # Yamada S5 / Wang 2020 S2 (A_cell = 907.9202768874502 micrometer^2)
    A_APC <- fixed(9e-08); label("Surface area of the APC [dm^2]")  # Yamada S5 / Wang 2020 S2 (A_APC = 900 micrometer^2)
    k_M1p0_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_M1p0_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p0_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex [1/day]")  # Yamada S5 / Wang 2020 S2 (k_M1p0_TCR_off = 1 1/second)
    TCR_p0_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell calculated based on the number of molecules for TCR and T cell ... [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (TCR_p0_tot = 103.81308067036068 molecule/micrometer^2)
    k_P1_up <- fixed(14.4); label("Rate of antigen uptake by APCs [1/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_P1_up = 14.4 1/day/cell)
    k_xP1_deg <- fixed(2); label("Rate of extracellular antigen degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_xP1_deg = 2 1/day)
    k_P1_deg <- fixed(17.28); label("Rate of endosomal antigen degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_P1_deg = 17.28 1/day)
    k_p1_deg <- fixed(144); label("Rate of endosomal epitope degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_p1_deg = 144 1/day)
    k_P1_on <- fixed(2.39117576834e-19); label("Rate of antigen-MHC binding [L/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_P1_on = 144000 1/day/molarity)
    k_P1_d1 <- fixed(2.408856796e+16); label("Antigen-MHC kd [count/L]")  # Yamada S5 / Wang 2020 S2 (k_P1_d1 = 4e-8 [2.7e-9 - 5.4e-7] molarity)
    p1_50 <- fixed(264550.26455); label("calculated based on the number of molecules for half-maximal T cell activation and synapse surface area [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (p1_50 = 2.6455026455026456e-05 molecule/micrometer^2)
    P1_C1 <- fixed(6.02214199e+15); label("Concentration of P1 in C1 [1/L]")  # Yamada S5 / Wang 2020 S2 (P1_C1 = 1e-08 molarity/cell)
    k_M1p1_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (k_M1p1_TCR_on = 1 1/(second*molecule/micrometer^2))
    k_M1p1_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex [1/day]")  # Yamada S5 / Wang 2020 S2 (k_M1p1_TCR_off = 1 1/second)
    k_M1p1_TCR_p <- fixed(86400); label("Rate of MHC-peptide-TCR complex modification [1/day]")  # Yamada S5 / Wang 2020 S2 (k_M1p1_TCR_p = 1 1/second)
    phi_M1p1_TCR <- fixed(7776); label("Rate of MHC-peptide-TCR complex with maximal modification that leads to non-signaling [1/day]")  # Yamada S5 / Wang 2020 S2 (phi_M1p1_TCR = 0.09 1/second)
    N_M1p1_TCR <- fixed(10); label("Number of modifcation steps for MHC-peptide-TCR complex [dimensionless]")  # Yamada S5 / Wang 2020 S2 (N_M1p1_TCR = 10 dimensionless)
    TCR_p1_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell calculated based on the number of molecules for TCR and T cell ... [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (TCR_p1_tot = 103.81308067036068 molecule/micrometer^2)
    kon_PD1_PDL1 <- fixed(8.3691151892e-07); label("kon of PD1-PDL1 binding kon for PD1-PDL1 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_PD1_PDL1 = 0.058333333333333334 1/(micromolarity*nanometer*second))
    q_P_pembro <- fixed(0.79488); label("Volumetric flow rate of pembro into peripheral tissues [L/day]")  # Yamada S5 (q_P_pembro = 9.2e-06 liter/second)
    q_T_pembro <- fixed(0.00736128); label("Volumetric flow rate of pembro into tumor tissues [L/day]")  # Yamada S5 / Wang 2020 S2 (q_T_pembro = 8.52e-05 milliliter/second)
    q_LN_pembro <- fixed(0.0002808); label("Volumetric flow rate of pembro into lymph nodes [L/day]")  # Yamada S5 / Wang 2020 S2 (q_LN_pembro = 3.25e-06 milliliter/second)
    q_LD_pembro <- fixed(2.16); label("Lymphatic flow rate of pembro from tumor to LN [1/day]")  # Yamada S5 / Wang 2020 S2 (q_LD_pembro = 0.0015 1/minute)
    k_cl_pembro <- fixed(0.22); label("Clearance rate of pembro [L/day]")  # Yamada S5 (k_cl_pembro = 0.22 liter/day)
    gamma_C_pembro <- fixed(0.58); label("Volume fraction available to pembro in central compartment [dimensionless]")  # Yamada S5 (gamma_C_pembro = 0.58 dimensionless)
    gamma_P_pembro <- fixed(0.0662); label("Volume fraction available to pembro in peripheral compartment [dimensionless]")  # Yamada S5 (gamma_P_pembro = 0.0662 dimensionless)
    gamma_T_pembro <- fixed(0.522); label("Volume fraction available to pembro in tumor compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_T_pembro = 0.522 dimensionless)
    gamma_LN_pembro <- fixed(0.2); label("Volume fraction available to pembro in lymph node compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_LN_pembro = 0.2 dimensionless)
    q_P_durv <- fixed(0.76032); label("Volumetric flow rate of durv into peripheral tissues [L/day]")  # Wang 2020 SBML deposit, not listed in Yamada S5 (q_P_durv = 8.8e-06 liter/second); Wang Table S2 prints 1.73655e-07 1/second, dimensionally inconsistent with the rate law
    q_T_durv <- fixed(0.00736128); label("Volumetric flow rate of durv into tumor tissues [L/day]")  # Wang 2020 SBML deposit, not listed in Yamada S5 (q_T_durv = 8.52e-05 milliliter/second); Wang Table S2 prints 8.52e-06 1/second, dimensionally inconsistent with the rate law
    q_LN_durv <- fixed(0.0002808); label("Volumetric flow rate of durv into lymph nodes [L/day]")  # Wang 2020 SBML deposit, not listed in Yamada S5 (q_LN_durv = 3.25e-06 milliliter/second); Wang Table S2 prints 1.73655e-07 1/second, dimensionally inconsistent with the rate law
    q_LD_durv <- fixed(2.16); label("Lymphatic flow rate of durv from tumor to LN [1/day]")  # Wang 2020 S2, not listed in Yamada S5 (q_LD_durv = 0.0015 1/minute)
    k_cl_durv <- fixed(0.4); label("Clearance rate of durv [L/day]")  # Wang 2020 SBML deposit, not listed in Yamada S5 (k_cl_durv = 0.4 liter/day); Wang Table S2 prints 7.1813e-07 1/second, dimensionally inconsistent with the rate law
    gamma_C_durv <- fixed(0.58); label("Volume fraction available to durv in central compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_C_durv = 0.58 dimensionless)
    gamma_P_durv <- fixed(0.056); label("Volume fraction available to durv in peripheral compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_P_durv = 0.056 dimensionless)
    gamma_T_durv <- fixed(0.522); label("Volume fraction available to durv in tumor compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_T_durv = 0.522 dimensionless)
    gamma_LN_durv <- fixed(0.2); label("Volume fraction available to durv in lymph node compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_LN_durv = 0.2 dimensionless)
    q_P_ipi <- fixed(0.9936); label("Volumetric flow rate of ipi into peripheral tissues [L/day]")  # Wang 2020 S2, not listed in Yamada S5 (q_P_ipi = 1.15e-05 liter/second)
    q_T_ipi <- fixed(0.00736128); label("Volumetric flow rate of ipi into tumor tissues [L/day]")  # Wang 2020 S2, not listed in Yamada S5 (q_T_ipi = 8.52e-05 milliliter/second)
    q_LN_ipi <- fixed(0.0002808); label("Volumetric flow rate of ipi into lymph nodes [L/day]")  # Wang 2020 S2, not listed in Yamada S5 (q_LN_ipi = 3.25e-06 milliliter/second)
    q_LD_ipi <- fixed(2.16); label("Lymphatic flow rate of ipi from tumor to LN [1/day]")  # Wang 2020 S2, not listed in Yamada S5 (q_LD_ipi = 0.0015 1/minute)
    k_cl_ipi <- fixed(0.524); label("Clearance rate of ipi [L/day]")  # Wang 2020 S2, not listed in Yamada S5 (k_cl_ipi = 0.524 liter/day)
    gamma_C_ipi <- fixed(0.686); label("Volume fraction available to ipi in central compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_C_ipi = 0.686 dimensionless)
    gamma_P_ipi <- fixed(0.0496); label("Volume fraction available to ipi in peripheral compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_P_ipi = 0.0496 dimensionless)
    gamma_T_ipi <- fixed(0.522); label("Volume fraction available to ipi in tumor compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_T_ipi = 0.522 dimensionless)
    gamma_LN_ipi <- fixed(0.2); label("Volume fraction available to ipi in lymph node compartment [dimensionless]")  # Wang 2020 S2, not listed in Yamada S5 (gamma_LN_ipi = 0.2 dimensionless)
    kon_PD1_PDL2 <- fixed(1.09994085344e-06); label("kon of PD1-PDL2 binding kon for PD1-PDL2 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_PD1_PDL2 = 0.07666666666666666 1/(micromolarity*nanometer*second))
    kon_PD1_pembro <- fixed(1.492094e-13); label("kon of PD1-pembrolizumab binding [L/(count*day)]")  # Yamada S5 (kon_PD1_pembro = 1040000 1/(molarity*second))
    kon_PDL1_durv <- fixed(6.16923348232e-14); label("kon of PDL1-durvalumab binding kon for PDL1-durvalumab [L/(count*day)]")  # Wang 2020 S2, not listed in Yamada S5 (kon_PDL1_durv = 430000 1/(molarity*second))
    kon_CD28_CD80 <- fixed(1.91294061467e-06); label("kon of CD28-CD80 binding kon for CD28-CD80 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_CD28_CD80 = 0.13333333333333333 1/(micromolarity*nanometer*second))
    kon_CD28_CD86 <- fixed(6.69529215136e-06); label("kon of CD28-CD86 binding kon for CD28-CD86 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_CD28_CD86 = 0.4666666666666667 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD80 <- fixed(1.05211733807e-05); label("kon of CTLA4-CD80 binding kon for CTLA4-CD80 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_CTLA4_CD80 = 0.7333333333333333 1/(micromolarity*nanometer*second))
    kon_CTLA4_CD86 <- fixed(9.56470307337e-06); label("kon of CTLA4-CD86 binding kon for CTLA4-CD86 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_CTLA4_CD86 = 0.6666666666666666 1/(micromolarity*nanometer*second))
    kon_CD80_PDL1 <- fixed(1.51122308559e-05); label("kon of CD80-PDL1 binding kon for CD80-PDL1 [dm^2/(count*day)]")  # Yamada S5 / Wang 2020 S2 (kon_CD80_PDL1 = 1.0533333333333332 1/(micromolarity*nanometer*second))
    kon_CTLA4_ipi <- fixed(5.49492191565e-14); label("kon of CTLA4-ipilimumab binding kon for CTLA4-ipilimumab [L/(count*day)]")  # Wang 2020 S2, not listed in Yamada S5 (kon_CTLA4_ipi = 383000 1/(molarity*second))
    koff_PD1_PDL1 <- fixed(123984); label("koff of PD1-PDL1 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_PD1_PDL1 = 1.435 1/second)
    koff_PD1_PDL2 <- fixed(45705.6); label("koff of PD1-PDL2 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_PD1_PDL2 = 0.529 1/second)
    koff_PD1_pembro <- fixed(2.6352); label("koff of PD1-pembrolizumab binding calculated from the measured kd and kon [1/day]")  # Yamada S5 (koff_PD1_pembro = 3.05e-05 1/second)
    koff_PDL1_durv <- fixed(2470.608); label("koff of PDL1-durvalumab binding calculated based on the measured kd and kon [1/day]")  # Wang 2020 S2, not listed in Yamada S5 (koff_PDL1_durv = 0.028595 1/second)
    koff_CD28_CD80 <- fixed(138240); label("koff of CD28-CD80 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_CD28_CD80 = 1.6 1/second)
    koff_CD28_CD86 <- fixed(2419200); label("koff of CD28-CD86 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_CD28_CD86 = 28 1/second)
    koff_CTLA4_CD80 <- fixed(38016); label("koff of CTLA4-CD80 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_CTLA4_CD80 = 0.44 1/second)
    koff_CTLA4_CD86 <- fixed(449280); label("koff of CTLA4-CD86 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_CTLA4_CD86 = 5.2 1/second)
    koff_CD80_PDL1 <- fixed(518745.6); label("koff of CD80-PDL1 binding calculated based on the measured kd and kon [1/day]")  # Yamada S5 / Wang 2020 S2 (koff_CD80_PDL1 = 6.004 1/second)
    koff_CTLA4_ipi <- fixed(602.25984); label("koff of CTLA4-ipilimumab binding calculated based on the measured kd and kon [1/day]")  # Wang 2020 S2, not listed in Yamada S5 (koff_CTLA4_ipi = 0.0069706 1/second)
    Chi_PD1_pembro <- fixed(333333330.953); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # Yamada S5 / Wang 2020 S2 (Chi_PD1_pembro = 3.333333309532274 1/nanometer)
    Chi_PDL1_durv <- fixed(3333333309.53); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # Wang 2020 S2, not listed in Yamada S5 (Chi_PDL1_durv = 33.33333309532274 1/nanometer)
    Chi_CTLA4_ipi <- fixed(333333330.953); label("Antibody cross-arm binding efficiency that also includes the conversion of kon from 3D to 2D [1/dm]")  # Wang 2020 S2, not listed in Yamada S5 (Chi_CTLA4_ipi = 3.333333309532274 1/nanometer)
    PD1_50 <- fixed(60000000000); label("PD1/PDL1 concentration for half-maximal T cell inactivation [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (PD1_50 = 6 molecule/micrometer^2)
    n_PD1 <- fixed(2); label("Hill coefficient for PD1/PDL1 half-maximal T cell inactivation [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_PD1 = 2 dimensionless)
    CD28_CD8X_50 <- fixed(2e+12); label("CD28-CD80/CD28-CD86 concentration for half-maximal T cell co-estimulation [count/dm^2]")  # Yamada S5 / Wang 2020 S2 (CD28_CD8X_50 = 200 [100 - 1000] molecule/micrometer^2)
    n_CD28_CD8X <- fixed(2); label("Hill coefficient for CD28-CD80/CD28-CD86 half-maximal T cell co-estimulation [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_CD28_CD8X = 2 dimensionless)
    T_PD1_total <- fixed(60000); label("concentration of PD1 on T cells [count]")  # Yamada S5 / Wang 2020 S2 (T_PD1_total = 60000 molecule)
    T_CD28_total <- fixed(184000); label("concentration of CD28 on T cells [count]")  # Yamada S5 / Wang 2020 S2 (T_CD28_total = 184000 molecule)
    T_CTLA4_syn <- fixed(400); label("concentration of CTLA4 on T cells [count]")  # Yamada S5 / Wang 2020 S2 (T_CTLA4_syn = 400 molecule)
    T_PDL1_total <- fixed(1600000); label("concentration of PDL1 on T cells [count]")  # Yamada S5 / Wang 2020 S2 (T_PDL1_total = 1600000 molecule)
    C1_PDL1_total <- fixed(1600000); label("number of PDL1 molecules per C1 cell [count]")  # Yamada S5 / Wang 2020 S2 (C1_PDL1_total = 1600000 [4e4 - 1.6e6] molecule)
    C1_PDL2_total <- fixed(2000); label("number of PDL2 molecules per C1 cell [count]")  # Yamada S5 / Wang 2020 S2 (C1_PDL2_total = 2000 [1000 - 40000] molecule)
    C1_CD80_total <- fixed(40000); label("number of CD80 molecules per C1 cell [count]")  # Yamada S5 / Wang 2020 S2 (C1_CD80_total = 40000 molecule)
    C1_CD86_total <- fixed(860000); label("number of CD86 molecules per C1 cell [count]")  # Yamada S5 / Wang 2020 S2 (C1_CD86_total = 860000 molecule)
    APC_PDL1_total <- fixed(1600000); label("number of PDL1 molecules per APC cell [count]")  # Yamada S5 / Wang 2020 S2 (APC_PDL1_total = 1600000 molecule)
    APC_PDL2_total <- fixed(2000); label("number of PDL2 molecules per APC cell [count]")  # Yamada S5 / Wang 2020 S2 (APC_PDL2_total = 2000 molecule)
    APC_CD80_total <- fixed(40000); label("number of CD80 molecules per APC cell [count]")  # Yamada S5 / Wang 2020 S2 (APC_CD80_total = 40000 molecule)
    APC_CD86_total <- fixed(860000); label("number of CD86 molecules per APC cell [count]")  # Yamada S5 / Wang 2020 S2 (APC_CD86_total = 860000 molecule)
    Treg_CTLA4_tot <- fixed(5000); label("Total number of CTLA4 on Treg cells [count]")  # Yamada S5 / Wang 2020 S2 (Treg_CTLA4_tot = 5000 molecule)
    Treg_CTLA4_50 <- fixed(1000); label("CTLA4 occupancy for half-maximal Treg inactivation by macrophages [count]")  # Yamada S5 / Wang 2020 S2 (Treg_CTLA4_50 = 1000 [100 - 1000] molecule)
    n_Treg_CTLA4 <- fixed(2); label("CTLA4 occupancy Hill coefficient for Treg inactivation by macrophages [dimensionless]")  # Yamada S5 / Wang 2020 S2 (n_Treg_CTLA4 = 2 dimensionless)
    k_CTLA4_ADCC <- fixed(0.1); label("Anti-CTLA4 ADCC (antibody-dependent cellular cytotoxicity) rate of Treg [1/day]")  # Yamada S5 / Wang 2020 S2 (k_CTLA4_ADCC = 0.1 [0.1 - 1] 1/day)
    k_rec_MDSC <- fixed(1.2); label("Rate of MDSC recruitment into the tumor [1/day]")  # Yamada S5 / Wang 2020 S2 (k_rec_MDSC = 1.2 [0.5 - 1.74] 1/day)
    kd_MDSC <- fixed(0.015); label("Rate of MDSC death [1/day]")  # Yamada S5 / Wang 2020 S2 (kd_MDSC = 0.015 1/day)
    IC50_ENT_C <- fixed(2.25228110426e+17); label("ENT concentration for half-maximal tumor cell death [count/L]")  # Yamada S5 / Wang 2020 S2 (IC50_ENT_C = 3.74e-7 [2.5e-9 - 3.74e-7] molarity)
    k_deg_CCL2 <- fixed(1.44); label("rate of CCL2 degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_deg_CCL2 = 0.06 1/hour)
    k_deg_NO <- fixed(135); label("rate of NO degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_deg_NO = 135 1/day)
    k_deg_ArgI <- fixed(0.173); label("rate of ArgI degradation [1/day]")  # Yamada S5 / Wang 2020 S2 (k_deg_ArgI = 0.173 1/day)
    k_sec_CCL2 <- fixed(85514.416258); label("rate of CCL2 secretion [1/day]")  # Yamada S5 / Wang 2020 S2 (k_sec_CCL2 = 1.42e-10 [8.2e-11 - 2.02e-10] nanomole/cell/day)
    k_sec_NO <- fixed(289062815.52); label("rate of NO secretion from MDSCs [1/day]")  # Yamada S5 / Wang 2020 S2 (k_sec_NO = 4.8e-7 [3.95e-7 - 5.65e-7] nanomole/cell/day)
    k_sec_ArgI <- fixed(8.430998786e+15); label("rate of ArgI secretion [1/day]")  # Yamada S5 / Wang 2020 S2 (k_sec_ArgI = 0.014 [0.012 - 0.016] (mU*microliter)/cell/day)
    IC50_ENT_NO <- fixed(3.3723995144e+14); label("half-maximal ENT concentration for NO inhibition [count/L]")  # Yamada S5 / Wang 2020 S2 (IC50_ENT_NO = 5.6e-10 [5.6e-11 - 5.6e-9] molarity)
    ki_Treg <- fixed(2.7); label("rate of ArgI-induced Treg expension [1/day]")  # Yamada S5 / Wang 2020 S2 (ki_Treg = 2.7 [0.27 - 2.7] 1/day)
    IC50_ArgI_CTL <- fixed(3.71566160783e+25); label("rate of ArgI-induced T cell death [count/L]")  # Yamada S5 / Wang 2020 S2 (IC50_ArgI_CTL = 61.7 [6.17 - 617] mU)
    IC50_NO_CTL <- fixed(4.5166064925e+14); label("rate of NO-induced T cell death [count/L]")  # Yamada S5 / Wang 2020 S2 (IC50_NO_CTL = 7.5e-10 [7.5e-11 - 7.5e-9] molarity)
    EC50_CCL2_rec <- fixed(3.011070995e+14); label("Half-maximal CCL2 level of MDSC recruitment [count/L]")  # Yamada S5 / Wang 2020 S2 (EC50_CCL2_rec = 5e-10 molarity)
    EC50_ArgI_Treg <- fixed(1.33089337979e+25); label("Half-maximal ArgI level of Treg expansion [count/L]")  # Yamada S5 / Wang 2020 S2 (EC50_ArgI_Treg = 22.1 [2.21 - 221] mU)
    MDSC_max <- fixed(163700000); label("Maximal MDSC density in the tumour [count/L]")  # Yamada S5 / Wang 2020 S2 (MDSC_max = 163700 [1e2 - 7.8e7] cell/milliliter)
    Treg_max <- fixed(955000); label("Maximal FoxP3+ T cell density in the tumor [count/L]")  # Yamada S5 / Wang 2020 S2 (Treg_max = 955 [1 - 1.65e7] cell/milliliter)
    IC50_ENT_CCL2 <- fixed(7.226570388e+14); label("half-maximal ENT concentration for CCL2 inhibition [count/L]")  # Yamada S5 / Wang 2020 S2 (IC50_ENT_CCL2 = 1.2e-9 [1.2e-10 - 1.2e-8] molarity)
    k_brec_MDSC <- fixed(0.0021); label("Baseline rate of MDSC migration into the tumor [1/day]")  # Yamada S5 / Wang 2020 S2 (k_brec_MDSC = 0.0021 [0.001 - 0.0031] 1/day)
    IC50_ENT_ArgI <- fixed(2.252281104e+17); label("half-maximal ENT concentration for Arg I inhibition [count/L]")  # Yamada S5 (IC50_ENT_ArgI = 3.74e-07 molarity, matching the Wang 2020 SBML deposit; Wang Table S2 prints 5e-07)
    k_a1_ENT <- fixed(45.6); label("Buccal absorption rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_a1_ENT = 1.9 [1.9 - 2.24] 1/hour)
    k_a2_ENT <- fixed(59.28); label("GI absorption rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_a2_ENT = 2.47 [2.47 - 7.47] 1/hour); SBML deposit carries 2.5
    k_cln_ENT <- fixed(1.76628443751e+21); label("Non-linear clearance rate of entinostat [count/(L*day)]")  # Yamada S5 / Wang 2020 S2 (k_cln_ENT = 1.22e-4 [6e-5 - 1.22e-4] mole/(liter*hour)); full precision 0.000122207654715555 from the SBML deposit
    Kc_ENT <- fixed(5.64763049312e+18); label("Kd for non-linear clearnance [count/L]")  # Yamada S5 / Wang 2020 S2 (Kc_ENT = 9.378109155345853e-06 mole/liter)
    lagP <- fixed(0.151966666667); label("Tlag (entinostat) [day]")  # Yamada S5 / Wang 2020 S2 (lagP = 3.6472 hour)
    durP <- fixed(0.0074875); label("Duration of zero-order absorption (entinostat) [day]")  # Yamada S5 / Wang 2020 S2 (durP = 0.1797 hour)
    k_dose2 <- fixed(0.3312); label("First-order absorption rate [1/day]")  # Yamada S5 / Wang 2020 S2 (k_dose2 = 0.0138 1/hour)
    q_P_ENT <- fixed(27648); label("Volumetric flow rate of entinostat into peripheral tissues [L/day]")  # Yamada S5 / Wang 2020 S2 (q_P_ENT = 320 milliliter/second)
    q_T_ENT <- fixed(276480); label("Volumetric flow rate of entinostat into tumor tissues [L/day]")  # Yamada S5 / Wang 2020 S2 (q_T_ENT = 3200 [3200 - 3420] milliliter/second)
    q_LN_ENT <- fixed(27648); label("Volumetric flow rate of entinostat into lymph nodes [L/day]")  # Yamada S5 / Wang 2020 S2 (q_LN_ENT = 320 milliliter/second)
    q_LD_ENT <- fixed(2.16); label("Lymphatic flow rate of entinostat from tumor to LN [1/day]")  # Yamada S5 / Wang 2020 S2 (q_LD_ENT = 0.0015 1/minute)
    k_cl_ENT <- fixed(3.648); label("Clearance rate of entinostat [1/day]")  # Yamada S5 / Wang 2020 S2 (k_cl_ENT = 0.152 [0.152 - 1.27] 1/hour)
    gamma_C_ENT <- fixed(0.55); label("Volume fraction available to entinostat in central compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_C_ENT = 0.55 dimensionless)
    gamma_P_ENT <- fixed(0.062); label("Volume fraction available to entinostat in peripheral compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_P_ENT = 0.062 dimensionless)
    gamma_T_ENT <- fixed(0.611); label("Volume fraction available to entinostat in tumor compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_T_ENT = 0.611 dimensionless)
    gamma_LN_ENT <- fixed(0.2); label("Volume fraction available to entinostat in lymph node compartment [dimensionless]")  # Yamada S5 / Wang 2020 S2 (gamma_LN_ENT = 0.2 dimensionless)

    # ---- entinostat dose split (Figure 2A; value not printed, see vignette) ----
    F_ENT_buccal <- fixed(0.324); label("Fraction of an entinostat dose absorbed through the zero-order buccal route [dimensionless]")  # not printed; back-solved by the maintainers so the simulated plasma Cmax after 2 mg/m^2 (1.7 m^2) equals the 15.4 ng/mL in the Results text (linear in F; 4 and 6 mg/m^2 follow)

    # ---- oncolytic-virus module (Yamada Table S5, derived from the mouse QSP) ----
    Beta_infect <- fixed(0.2); label("Viral infection rate constant [cell/day/molecule]")  # Yamada S5 (Beta_infect = 0.200 cell/day/molecule, Adjusted to reduce single-agent ORR below 50%)
    Alpha <- fixed(8.6); label("Viral production size [molecule/cell]")  # Yamada S5 (Alpha = 8.6 molecule/cell, Derived from mouse QSP)
    Delta_V <- fixed(2.39); label("Virus elimination rate constant [1/day]")  # Yamada S5 (Delta_V = 2.39 1/day, Derived from mouse QSP)
    Delta_I <- fixed(57.4); label("Infected-cell lysis rate constant [1/day]")  # Yamada S5 (Delta_I = 57.4 1/day, Derived from mouse QSP)
    IL12_max <- fixed(1.56875); label("Maximum activation fold for T1 activation rate by IL-12 [dimensionless]")  # Yamada S5 (IL12_max = 1.56875, Derived from mouse QSP)
    IL7_max <- fixed(0.696667); label("Maximum activation fold for naive-T-cell influx by IL-7 [dimensionless]")  # Yamada S5 (IL7_max = 0.696667, Derived from mouse QSP)
    IL12_50 <- fixed(3887.359); label("Michaelis-Menten constant for IL-12 activation of T1 activation [cell]")  # Yamada S5 (IL12_50 = 3887.359 cell, Derived from mouse QSP)
    IL7_50 <- fixed(38873.59); label("Michaelis-Menten constant for IL-7 activation of naive-T-cell influx [cell]")  # Yamada S5 (IL7_50 = 38873.59 cell, Derived from mouse QSP)
    Imax_v <- fixed(0.139); label("Maximum fold enhancement of virus clearance by pembrolizumab [dimensionless]")  # Yamada S5 (Imax_v = 0.139, Derived from mouse QSP)
    IC50_v <- fixed(0.029); label("Half-maximal pembrolizumab serum concentration enhancing virus clearance [nmol/L]")  # Yamada S5 (IC50_v = 0.029 nanomole/liter, Derived from mouse QSP)
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
    # Oncolytic-virus layer: the tumour cancer-cell pool is split into
    # uninfected (Cuit) and infected (Cit) cells; V_T.C1 = Cuit + Cit is the
    # total cancer-cell count the immune network sees (Yamada Table S4 rule).
    # C2 is a parallel no-virus cancer tracker; Virus_T is the intratumoral
    # virus amount.
    x_V_T_Cuit <- q_V_T_Cuit
    x_V_T_Cit <- q_V_T_Cit
    x_V_T_C1 <- x_V_T_Cuit + x_V_T_Cit
    x_V_T_C2 <- q_V_T_C2
    x_V_T_Virus_T <- q_V_T_Virus_T
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
    vol_V_T <- V_Tmin+(vol_cell*x_V_T_C_x)+(vol_Tcell*x_V_T_T_exh)+(vol_cell*x_V_T_C1)+(vol_Tcell*x_V_T_T0)+(vol_Tcell*x_V_T_T1)   # Table S4 rule 1 (tumour volume)
    x_V_C_pembro <- q_V_C_pembro*nmol_to_molecule/m_V_C
    # central pembrolizumab concentration in nmol/L (nM) for the virus-clearance
    # enhancement term, whose IC50_v is reported in nmol/L
    pembro_conc_nM <- q_V_C_pembro/m_V_C
    x_V_C_durv <- q_V_C_durv*nmol_to_molecule/m_V_C
    x_V_C_ipi <- q_V_C_ipi*nmol_to_molecule/m_V_C
    x_V_C_ENT <- q_V_C_ENT*nmol_to_molecule/m_V_C
    x_V_C_ENT_Buccal <- q_V_C_ENT_Buccal*nmol_to_molecule/m_V_C
    x_V_C_ENT_GI <- q_V_C_ENT_GI*nmol_to_molecule/m_V_C
    x_V_C_Dose2 <- q_V_C_Dose2*nmol_to_molecule/m_V_C
    x_V_P_pembro <- q_V_P_pembro*nmol_to_molecule/m_V_P
    x_V_P_durv <- q_V_P_durv*nmol_to_molecule/m_V_P
    x_V_P_ipi <- q_V_P_ipi*nmol_to_molecule/m_V_P
    x_V_P_ENT <- q_V_P_ENT*nmol_to_molecule/m_V_P
    x_V_T_c <- q_V_T_c/vol_V_T
    x_V_T_pembro <- q_V_T_pembro*nmol_to_molecule/vol_V_T
    x_V_T_durv <- q_V_T_durv*nmol_to_molecule/vol_V_T
    x_V_T_ipi <- q_V_T_ipi*nmol_to_molecule/vol_V_T
    x_V_T_CCL2 <- q_V_T_CCL2/vol_V_T
    x_V_T_NO <- q_V_T_NO/vol_V_T
    x_V_T_ArgI <- q_V_T_ArgI/vol_V_T
    x_V_T_ENT <- q_V_T_ENT*nmol_to_molecule/vol_V_T
    x_V_LN_IL2 <- q_V_LN_IL2/m_V_LN
    x_V_LN_P0 <- q_V_LN_P0/m_V_LN
    x_V_LN_P1 <- q_V_LN_P1/m_V_LN
    x_V_LN_pembro <- q_V_LN_pembro*nmol_to_molecule/m_V_LN
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
    x_syn_T_C1_PD1_pembro <- q_syn_T_C1_PD1_pembro/m_syn_T_C1
    x_syn_T_C1_PD1_pembro_PD1 <- q_syn_T_C1_PD1_pembro_PD1/m_syn_T_C1
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
    x_syn_T_APC_PD1_pembro <- q_syn_T_APC_PD1_pembro/m_syn_T_APC
    x_syn_T_APC_PD1_pembro_PD1 <- q_syn_T_APC_PD1_pembro_PD1/m_syn_T_APC
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

    # ---- repeated-assignment rules (Supplementary Table S4) ----
    C_total <- x_V_T_C1   # Table S4 (C_total)
    # IL-7 / IL-12 immunostimulation driven by infected-cell count (Yamada
    # Table S4 rules IL7_act, IL12_act)
    IL7_act <- IL7_max*x_V_T_Cit/(IL7_50+x_V_T_Cit)
    IL12_act <- IL12_max*x_V_T_Cit/(IL12_50+x_V_T_Cit)
    T_total <- x_V_T_T0+x_V_T_T1   # Table S4 (T_total)
    T_total_LN <- x_V_LN_T1   # Table S4 (T_total_LN)
    Tregs_ <- x_V_T_T0   # Table S4 (Tregs_)
    H_APC <- (n_sites_APC*x_V_LN_APC)/((n_sites_APC*x_V_LN_APC)+T_total_LN+1)   # Table S4 (H_APC)
    H_mAPC <- (n_sites_APC*x_V_LN_mAPC)/((n_sites_APC*x_V_LN_mAPC)+T_total_LN+1)   # Table S4 (H_mAPC)
    pTCR_p0_MHC_tot <- 0.5*(((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))-(TCR_p0_tot*sqrt((((((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))/TCR_p0_tot))^2-(((4*x_A_s_M1p0)/n_T0_clones)/TCR_p0_tot)))))   # Table S4 (pTCR_p0_MHC_tot)
    H_P0 <- pTCR_p0_MHC_tot/(pTCR_p0_MHC_tot+p0_50)   # Table S4 (H_P0)
    pTCR_p1_MHC_tot <- (a_k_M1p1_TCR_off/(a_k_M1p1_TCR_off+a_phi_M1p1_TCR))*((a_k_M1p1_TCR_p/(a_k_M1p1_TCR_off+a_k_M1p1_TCR_p)))^N_M1p1_TCR*0.5*(((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))-(TCR_p1_tot*sqrt((((((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))/TCR_p1_tot))^2-(((4*x_A_s_M1p1)/n_T1_clones)/TCR_p1_tot)))))   # Table S4 (pTCR_p1_MHC_tot)
    H_P1 <- pTCR_p1_MHC_tot/(pTCR_p1_MHC_tot+p1_50)   # Table S4 (H_P1)
    H_PD1_C1 <- (((x_syn_T_C1_PD1_PDL1+x_syn_T_C1_PD1_PDL2)/PD1_50))^n_PD1/((((x_syn_T_C1_PD1_PDL1+x_syn_T_C1_PD1_PDL2)/PD1_50))^n_PD1+1)   # Table S4 (H_PD1_C1)
    H_CD28_C1 <- (((x_syn_T_C1_CD28_CD80+x_syn_T_C1_CD28_CD86+(2*x_syn_T_C1_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X/((((x_syn_T_C1_CD28_CD80+x_syn_T_C1_CD28_CD86+(2*x_syn_T_C1_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X+1)   # Table S4 (H_CD28_C1)
    H_PD1_APC <- (((x_syn_T_APC_PD1_PDL1+x_syn_T_APC_PD1_PDL2)/PD1_50))^n_PD1/((((x_syn_T_APC_PD1_PDL1+x_syn_T_APC_PD1_PDL2)/PD1_50))^n_PD1+1)   # Table S4 (H_PD1_APC)
    H_CD28_APC <- (((x_syn_T_APC_CD28_CD80+x_syn_T_APC_CD28_CD86+(2*x_syn_T_APC_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X/((((x_syn_T_APC_CD28_CD80+x_syn_T_APC_CD28_CD86+(2*x_syn_T_APC_CD28_CD80_CD28))/CD28_CD8X_50))^n_CD28_CD8X+1)   # Table S4 (H_CD28_APC)
    H_Treg_T <- (((x_V_T_Treg_CTLA4_ipi+(2*x_V_T_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4/((((x_V_T_Treg_CTLA4_ipi+(2*x_V_T_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4+1)   # Table S4 (H_Treg_T)
    H_Treg_P <- (((x_V_P_Treg_CTLA4_ipi+(2*x_V_P_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4/((((x_V_P_Treg_CTLA4_ipi+(2*x_V_P_Treg_CTLA4_ipi_CTLA4))/Treg_CTLA4_50))^n_Treg_CTLA4+1)   # Table S4 (H_Treg_P)
    H_MDSC_C1 <- 1-((1-(x_V_T_NO/(IC50_NO_CTL+x_V_T_NO)))*(1-(x_V_T_ArgI/(IC50_ArgI_CTL+x_V_T_ArgI)))*(1-(x_V_T_ENT/(x_V_T_ENT+IC50_ENT_ArgI))))   # Table S4 (H_MDSC_C1)
    H_ENT_C1 <- x_V_T_ENT/(x_V_T_ENT+IC50_ENT_C)   # Table S4 (H_ENT_C1)
    R_Tcell <- (((k_T1*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*(1-H_PD1_C1)*(1-H_MDSC_C1))   # Table S4 (R_Tcell)
    N_aT <- N0+(N_costim*H_CD28_APC)+((N_IL2*x_V_LN_IL2)/(IL2_50+x_V_LN_IL2))   # Table S4 (N_aT)

    # Wang 2020 Table S6 event (inherited platform; Yamada tabulates no events):
    # cancer cells are set to zero once fewer than 0.9 cells
    # remain. rxode2 has no state-reset events, so proliferation is switched
    # off below that threshold instead and the remaining fraction of a cell
    # decays through the death terms.
    c1_alive <- 1*(x_V_T_C1 >= 0.9)

    # ---- reaction rates (Supplementary Table S3), amount per day ----
    v1 <- k_cell_clear*x_V_T_C_x   # R1 Clearance of dead cancer cells from tumor
    v2 <- k_cell_clear*x_V_T_T_exh   # R2 Clearance of dead T cells from tumor
    v3_1 <- (k_C1_growth*x_V_T_Cuit*(1-(C_total/C_max))*(1-H_ENT_C1))*c1_alive   # R3_1 Uninfected cancer cell growth (logistic)
    v4_1 <- a_k_C1_death*x_V_T_Cuit   # R4_1 Uninfected cancer cell death
    v4_2 <- a_k_C1_death*x_V_T_Cit   # R4_2 Infected cancer cell death
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
    v21 <- Q_T1_in*n_T1_clones*(1+IL7_act)   # R21 Naive T cell entry into the lymph node (IL-7 enhanced)
    v22 <- Q_T1_out*x_V_LN_nT1   # R22 Naive T cell exit from the lymph node
    v23 <- k_T1_act*H_mAPC*H_P1*x_V_LN_nT1*(1+IL12_act)   # R23 Naive T cell activation (IL-12 enhanced)
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
    v37_1 <- ((k_C_T1*x_V_T_T1*x_V_T_Cuit)/(C_total+T_total+1))*(1-H_PD1_C1*0.65)*(1-H_MDSC_C1)   # R37_1 Uninfected cancer cell killing by T cells (Yamada uses 0.65*H_PD1_C1)
    v37_2 <- ((k_C_T1*x_V_T_T1*x_V_T_Cit)/(C_total+T_total+1))*(1-H_PD1_C1*0.65)*(1-H_MDSC_C1)   # R37_2 Infected cancer cell killing by T cells
    v38 <- k_APC_death*((APC0_T*vol_V_T)-x_V_T_APC)   # R38 APC recruitment/death in the tumour
    v39 <- k_APC_death*((APC0_LN*m_V_LN)-x_V_LN_APC)   # R39 APC recruitment/death in LN
    v40 <- ((k_APC_mat*x_V_T_c)/(x_V_T_c+c50))*x_V_T_APC   # R40 APC maturation in the tumour
    v41 <- k_APC_mig*x_V_T_mAPC   # R41 APC migration to the lymph node
    v42 <- k_mAPC_death*x_V_T_mAPC   # R42 mAPC death in the tumour
    v43 <- k_mAPC_death*x_V_LN_mAPC   # R43 mAPC death in the lymph node
    v44 <- k_c*(c0-x_V_T_c)*vol_V_T   # R44 Baseline cytokine secretion/degradation
    v45 <- R_Tcell*DAMPs   # R45 Cytokine release in response to the tumour
    v46 <- (kout*x_A_e_M1*m_A_e)-(kin*x_A_s_M1*m_A_s)   # R46 MHC translocation
    v47 <- n_T0_clones*P0_C1*((a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1*0.65)))*x_V_T_C1+Delta_I*x_V_T_Cit)*vol_V_T   # R47 Antigen deposition from dying cancer cells (incl. virus-lysed infected cells)
    v48 <- k_xP0_deg*x_V_LN_P0*m_V_LN   # R48 Free antigen degradation
    v49 <- k_P0_up*x_V_LN_mAPC*x_V_LN_P0*m_V_LN   # R49 Antigen uptake by mature antigen presenting cells
    v50 <- k_P0_up*1*x_V_LN_P0*m_V_e   # R50 Antigen uptake by mature antigen presenting cells
    v51 <- k_P0_deg*x_V_e_P0*m_V_e   # R51 Antigen degradation in APC endosomes
    v52 <- k_p0_deg*x_V_e_p0*m_V_e   # R52 Epitope degradation in APC endosomes
    v53 <- k_P0_on*x_V_e_p0*x_A_e_M1*m_A_e   # R53 Antigen-MHC binding in endosome
    v54 <- k_P0_d1*k_P0_on*x_A_e_M1p0*m_A_e   # R54 Antigen-MHC unbinding in endosome
    v55 <- k_P0_d1*k_P0_on*x_A_s_M1p0*m_A_s   # R55 Antigen-MHC unbinding on APC surface
    v56 <- kout*x_A_e_M1p0*m_A_e   # R56 Antigen-MHC translocation
    v57 <- n_T1_clones*P1_C1*((a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1*0.65)))*x_V_T_C1+Delta_I*x_V_T_Cit)*vol_V_T   # R57 Antigen deposition from dying cancer cells (incl. virus-lysed infected cells)
    v58 <- k_xP1_deg*x_V_LN_P1*m_V_LN   # R58 Free antigen degradation
    v59 <- k_P1_up*x_V_LN_mAPC*x_V_LN_P1*m_V_LN   # R59 Antigen uptake by mature antigen presenting cells
    v60 <- k_P1_up*1*x_V_LN_P1*m_V_e   # R60 Antigen uptake by mature antigen presenting cells
    v61 <- k_P1_deg*x_V_e_P1*m_V_e   # R61 Antigen degradation in APC endosomes
    v62 <- k_p1_deg*x_V_e_p1*m_V_e   # R62 Epitope degradation in APC endosomes
    v63 <- k_P1_on*x_V_e_p1*x_A_e_M1*m_A_e   # R63 Antigen-MHC binding in endosome
    v64 <- k_P1_d1*k_P1_on*x_A_e_M1p1*m_A_e   # R64 Antigen-MHC unbinding in endosome
    v65 <- k_P1_d1*k_P1_on*x_A_s_M1p1*m_A_s   # R65 Antigen-MHC unbinding on APC surface
    v66 <- kout*x_A_e_M1p1*m_A_e   # R66 Antigen-MHC translocation
    v67 <- q_P_pembro*((x_V_C_pembro/gamma_C_pembro)-(x_V_P_pembro/gamma_P_pembro))   # R67 pembro diffusive transport to peripheral compartment
    v68 <- q_T_pembro*((x_V_C_pembro/gamma_C_pembro)-(x_V_T_pembro/gamma_T_pembro))   # R68 pembro diffusive transport to tumor compartment
    v69 <- q_LN_pembro*((x_V_C_pembro/gamma_C_pembro)-(x_V_LN_pembro/gamma_LN_pembro))   # R69 pembro diffusive transport to lymph node compartment
    v70 <- ((q_LD_pembro*x_V_T_pembro)/gamma_T_pembro)*vol_V_T   # R70 pembro convective transport from tumor to lymph node
    v71 <- ((q_LD_pembro*x_V_LN_pembro)/gamma_LN_pembro)*m_V_LN   # R71 pembro convective transport from lymph node to central
    v72 <- k_cl_pembro*x_V_C_pembro   # R72 Drug clearance from central compartment
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
    v87 <- ((2*kon_PD1_pembro*((x_syn_T_C1_PD1*x_V_T_pembro)/gamma_T_pembro))-(koff_PD1_pembro*x_syn_T_C1_PD1_pembro))*m_syn_T_C1   # R87 binding and unbinding of PD1 to Nivo on T surface in synapse
    v88 <- ((Chi_PD1_pembro*kon_PD1_pembro*x_syn_T_C1_PD1*x_syn_T_C1_PD1_pembro)-(2*koff_PD1_pembro*x_syn_T_C1_PD1_pembro_PD1))*m_syn_T_C1   # R88 binding and unbinding of PD1 to PD1-Nivo on T surface in synapse
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
    v109 <- ((2*kon_PD1_pembro*((x_syn_T_APC_PD1*x_V_LN_pembro)/gamma_LN_pembro))-(koff_PD1_pembro*x_syn_T_APC_PD1_pembro))*m_syn_T_APC   # R109 binding and unbinding of PD1 to Nivo on T surface in synapse
    v110 <- ((Chi_PD1_pembro*kon_PD1_pembro*x_syn_T_APC_PD1*x_syn_T_APC_PD1_pembro)-(2*koff_PD1_pembro*x_syn_T_APC_PD1_pembro_PD1))*m_syn_T_APC   # R110 binding and unbinding of PD1 to PD1-Nivo on T surface in synapse
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

    # ---- oncolytic-virus reactions (Yamada Table S3 rows 155-158) ----
    # R155 parallel no-virus cancer-cell tracker (net rate: growth - T kill - death)
    v155 <- k_C1_growth*x_V_T_C2*(1-(x_V_T_C2/C_max))*(1-H_ENT_C1) -
      ((k_C_T1*x_V_T_T1*x_V_T_C2)/(x_V_T_C2+T_total+1))*(1-0.65*H_PD1_C1)*(1-H_MDSC_C1) -
      k_C1_death*x_V_T_C2
    # R156 intratumoral virus: production proportional to infected-cell lysis,
    # pembrolizumab-enhanced clearance (net rate)
    v156 <- Alpha*Delta_I*x_V_T_Cit -
      Delta_V*(1+Imax_v*pembro_conc_nM/(IC50_v+pembro_conc_nM))*x_V_T_Virus_T
    # R157 infection of uninfected cancer cells by virus
    v157 <- Beta_infect*x_V_T_Cuit*x_V_T_Virus_T /
      (1+x_V_T_C_x+x_V_T_T_exh+x_V_T_C1+x_V_T_T0+x_V_T_T1)
    # R158 lysis of infected cancer cells
    v158 <- Delta_I*x_V_T_Cit

    # ---- ODEs ----
    d/dt(q_V_C_T0) <- -v10-v14+v15-v16+v17
    d/dt(q_V_C_T1) <- -v26-v32+v33-v34+v35
    d/dt(q_V_C_pembro) <- (-v67-v68-v69+v71-v72)/nmol_to_molecule
    d/dt(q_V_C_durv) <- (-v73-v74-v75+v77-v78)/nmol_to_molecule
    d/dt(q_V_C_ipi) <- (-v79-v80-v81+v83-v84)/nmol_to_molecule
    d/dt(q_V_C_ENT) <- (v146+v147-v148-v149-v150-v151+v153-v154)/nmol_to_molecule
    d/dt(q_V_C_ENT_Buccal) <- (-v146)/nmol_to_molecule
    d/dt(q_V_C_ENT_GI) <- (v145-v147)/nmol_to_molecule
    d/dt(q_V_C_Dose2) <- (-v145)/nmol_to_molecule
    d/dt(q_V_P_T0) <- -v11+v14-v15-v134
    d/dt(q_V_P_T1) <- -v27+v32-v33
    d/dt(q_V_P_pembro) <- (v67)/nmol_to_molecule
    d/dt(q_V_P_durv) <- (v73)/nmol_to_molecule
    d/dt(q_V_P_ipi) <- (v79)/nmol_to_molecule
    d/dt(q_V_P_Treg_CTLA4) <- -v131-v132
    d/dt(q_V_P_Treg_CTLA4_ipi) <- v131-v132
    d/dt(q_V_P_Treg_CTLA4_ipi_CTLA4) <- v132
    d/dt(q_V_P_ENT) <- (v149)/nmol_to_molecule
    d/dt(q_V_T_C_x) <- -v1+v4_1+v4_2+v37_1+v37_2+v158
    d/dt(q_V_T_T_exh) <- -v2+v12+v28+v30+v31
    d/dt(q_V_T_Cuit) <- v3_1-v4_1-v37_1-v157
    d/dt(q_V_T_Cit) <- v157-v4_2-v37_2-v158
    d/dt(q_V_T_Virus_T) <- v156
    d/dt(q_V_T_C2) <- v155
    d/dt(q_V_T_T0) <- -v12+v16-v133+v144
    d/dt(q_V_T_T1) <- -v28-v30-v31+v34
    d/dt(q_V_T_APC) <- v38-v40
    d/dt(q_V_T_mAPC) <- v40-v41-v42
    d/dt(q_V_T_c) <- v44+v45
    d/dt(q_V_T_pembro) <- (v68-v70)/nmol_to_molecule
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
    d/dt(q_V_LN_pembro) <- (v69+v70-v71)/nmol_to_molecule
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
    d/dt(q_syn_T_C1_PD1_pembro) <- v87-v88
    d/dt(q_syn_T_C1_PD1_pembro_PD1) <- v88
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
    d/dt(q_syn_T_APC_PD1_pembro) <- v109-v110
    d/dt(q_syn_T_APC_PD1_pembro_PD1) <- v110
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

    # ---- initial conditions (Supplementary Table S5 and its initial assignments) ----
    q_V_P_Treg_CTLA4(0) <- Treg_CTLA4_tot
    # untreated run-in starts from 1e6 uninfected cancer cells (Yamada Methods,
    # Human QSP Model Simulations: 'The initial cancer cell number was set as
    # 10^6 cells'; the Wang 2020 platform started from 1 cell). The C2 tracker
    # is initialised to the same count.
    q_V_T_Cuit(0) <- 1e6
    q_V_T_C2(0) <- 1e6
    q_V_T_Treg_CTLA4(0) <- Treg_CTLA4_tot
    q_V_LN_IL2(0) <- 670.052061644   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_LN.IL2 = 1e-18 molarity)
    q_V_LN_P0(0) <- 670.052061644   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_LN.P0 = 1e-18 molarity)
    q_V_LN_P1(0) <- 670.052061644   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_LN.P1 = 1e-18 molarity)
    q_V_e_P0(0) <- 2.408856796e-10   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_e.P0 = 1e-18 molarity)
    q_V_e_p0(0) <- 2.408856796e-10   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_e.p0 = 1e-18 molarity)
    q_V_e_P1(0) <- 2.408856796e-10   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_e.P1 = 1e-18 molarity)
    q_V_e_p1(0) <- 2.408856796e-10   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (V_e.p1 = 1e-18 molarity)
    q_A_e_M1(0) <- 32786.8852459   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_e.M1 = 2185.79235 molecule/micrometer^2)
    q_A_e_M1p0(0) <- 1.5e-05   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_e.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_e_M1p1(0) <- 1.5e-05   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_e.M1p1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1(0) <- 0.0009   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_s.M1 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p0(0) <- 0.0009   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_s.M1p0 = 1e-06 molecule/micrometer^2)
    q_A_s_M1p1(0) <- 0.0009   # Wang 2020 Table S5 initial amount, not tabulated in Yamada (A_s.M1p1 = 1e-06 molecule/micrometer^2)
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
    # remainder enters Dose2 as a bolus after the lag lagP (Figure 2A). ----
    f(q_V_C_ENT_Buccal) <- F_ENT_buccal
    dur(q_V_C_ENT_Buccal) <- durP
    f(q_V_C_Dose2) <- 1 - F_ENT_buccal
    alag(q_V_C_Dose2) <- lagP

    # ---- outputs ----
    # tumour diameter [cm] from the total tumour volume assuming a sphere
    tumour_diameter <- (6*vol_V_T*1000/pi)^(1/3)
    # central (plasma) drug concentrations [nmol/L]
    Cc_pembrolizumab <- x_V_C_pembro/nmol_to_molecule
    Cc_ipilimumab <- x_V_C_ipi/nmol_to_molecule
    Cc_entinostat <- x_V_C_ENT/nmol_to_molecule
    # tumour-infiltrating cell densities [cells/mL of tumour]
    teff_density <- x_V_T_T1/(vol_V_T*1000)
    treg_density <- x_V_T_T0/(vol_V_T*1000)
    mdsc_density <- x_V_T_MDSC/(vol_V_T*1000)
    # oncolytic-virus layer outputs
    cancer_cells <- x_V_T_C1   # total tumour cancer-cell count (uninfected + infected)
    infected_fraction <- x_V_T_Cit/(x_V_T_C1+1)   # fraction of cancer cells infected
    virus_amount <- x_V_T_Virus_T   # intratumoral virus
  })
}
