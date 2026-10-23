Ma_2020_atezolizumab_cibisatamab_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of microsatellite-stable",
    "metastatic colorectal cancer treated with the anti-PD-L1 antibody",
    "atezolizumab, the CEA x CD3 bispecific T cell engager cibisatamab",
    "(CEA-TCB, RO6958688), or their combination. Four physiological",
    "compartments (central, peripheral, tumour, tumour-draining lymph node),",
    "APC endosomal and surface compartments, two PD-1 checkpoint synapses",
    "(Teff-cancer cell and Teff-APC) and two TCE synapses (Teff-cancer cell and",
    "Treg-cancer cell). Covers logistic cancer-cell growth from a single cell,",
    "naive CD8 T cell and Treg priming in the lymph node (kinetic-proofreading",
    "TCR signal and IL-2 driven division number), T cell trafficking, Treg",
    "suppression and exhaustion, APC maturation and migration, self- and",
    "neo-antigen release, uptake, endosomal processing and MHC presentation,",
    "PD-1/PD-L1/PD-L2 binding on cancer cells and on APCs with bivalent",
    "atezolizumab binding to PD-L1 (the APC synapse limits Teff proliferation in",
    "the lymph node), and CEA-cibisatamab-CD3 trimer formation driving",
    "TCE-mediated cancer-cell killing and TCE-activated Treg suppression. 73 ODE",
    "states, 22 repeated-assignment rules, 12 initial assignments and 105",
    "reactions, translated from the authors' SimBiology SBML export",
    "(Supplementary file 2) with the printed values of Supplementary Tables S2",
    "and S3. Deterministic mechanism model: the authors generated virtual",
    "patients by Latin hypercube sampling of the Table S3 ranges rather than by",
    "fitting IIV or residual error, so no etas and no error model are encoded.",
    sep = " "
  )
  reference <- paste(
    "Ma H, Wang H, Sove RJ, Wang J, Giragossian C, Popel AS. Combination",
    "therapy with T cell engager and PD-L1 blockade enhances the antitumor",
    "potency of T cells as predicted by a QSP model. J Immunother Cancer.",
    "2020;8(2):e001141. doi:10.1136/jitc-2020-001141.",
    sep = " "
  )
  vignette <- "Ma_2020_atezolizumab_cibisatamab_qsp"

  units <- list(
    time = "day",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  covariateData <- list()

  compartmentData <- list(
    q_V_C_T1 = list(
      analyte = "Number of T1 cells in the central compartment",
      units = "cell",
      specimen = "whole blood",
      verified = FALSE
    ),
    q_V_C_T0 = list(
      analyte = "Number of T0 cells in the central compartment",
      units = "cell",
      specimen = "whole blood",
      verified = FALSE
    ),
    q_V_C_atezolizumab = list(
      analyte = "Concentration of atezolizumab in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_C_cibis = list(
      analyte = "Concentration of cibis in central compartment",
      units = "nmol",
      specimen = "plasma",
      verified = FALSE
    ),
    q_V_P_T1 = list(
      analyte = "Number of T1 cells in the peripheral compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_T0 = list(
      analyte = "Number of T0 cells in the peripheral compartment",
      units = "cell",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_atezolizumab = list(
      analyte = "Concentration of atezolizumab in peripheral compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_P_cibis = list(
      analyte = "Concentration of cibis in peripheral compartment",
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
    q_V_T_T1 = list(
      analyte = "Number of T1 cells in the tumour compartment",
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
    q_V_T_atezolizumab = list(
      analyte = "Concentration of atezolizumab in tumour compartment",
      units = "nmol",
      specimen = "tumor",
      verified = FALSE
    ),
    q_V_T_cibis = list(
      analyte = "Concentration of cibis in tumour compartment",
      units = "nmol",
      specimen = "tumor",
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
    q_V_LN_IL2 = list(
      analyte = "Concentration of IL2 in the lymph node compartment",
      units = "molecule",
      specimen = "tissue",
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
    q_V_LN_atezolizumab = list(
      analyte = "Concentration of atezolizumab in lymph node compartment",
      units = "nmol",
      specimen = "tissue",
      verified = FALSE
    ),
    q_V_LN_cibis = list(
      analyte = "Concentration of cibis in lymph node compartment",
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
    q_A_s_M1 = list(analyte = "M1 in A_s", units = "molecule", specimen = "tissue", verified = FALSE),
    q_A_s_M1p0 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "tissue", verified = FALSE),
    q_A_s_M1p1 = list(analyte = "Antigen-MHC complex", units = "molecule", specimen = "tissue", verified = FALSE),
    q_syn_T1_C1_PD1_PDL1 = list(
      analyte = "concentration of PD1-PDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PD1_PDL2 = list(
      analyte = "concentration of PD1-PDL2 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PD1 = list(
      analyte = "concentration of PD1 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PDL1 = list(
      analyte = "concentration of PDL1 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PDL2 = list(
      analyte = "concentration of PDL2 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PD1_aPD1 = list(
      analyte = "concentration of PD1-nivolumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PD1_aPD1_PD1 = list(
      analyte = "concentration of PD1-nivolumab-PD1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PDL1_aPDL1 = list(
      analyte = "concentration of PDL1-durvalumab complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_C1_PDL1_aPDL1_PDL1 = list(
      analyte = "concentration of PDL1-durvalumab-PDL1 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_syn_T1_APC_PD1_PDL1 = list(
      analyte = "concentration of PD1-PDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PD1_PDL2 = list(
      analyte = "concentration of PD1-PDL2 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PD1 = list(
      analyte = "concentration of PD1 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PDL1 = list(
      analyte = "concentration of PDL1 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PDL2 = list(
      analyte = "concentration of PDL2 in synapse",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PD1_aPD1 = list(
      analyte = "concentration of PD1-nivolumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PD1_aPD1_PD1 = list(
      analyte = "concentration of PD1-nivolumab-PD1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PDL1_aPDL1 = list(
      analyte = "concentration of PDL1-durvalumab complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_syn_T1_APC_PDL1_aPDL1_PDL1 = list(
      analyte = "concentration of PDL1-durvalumab-PDL1 complex",
      units = "molecule",
      specimen = "tissue",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CEA_cibis_CEA_CD3 = list(
      analyte = "concentration of CEA-cibis-CEA_CD3 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CEA_cibis_CEA = list(
      analyte = "concentration of CEA-cibis-CEA complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CEA_cibis_CD3 = list(
      analyte = "concentration of CEA-cibis-CD3 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CD3_cibis = list(
      analyte = "concentration of CD3-cibis complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CEA_cibis = list(
      analyte = "concentration of CEA-cibis complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CEA = list(
      analyte = "concentration of CEA in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T1_C1_CD3 = list(
      analyte = "concentration of CD3 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CEA_cibis_CEA_CD3 = list(
      analyte = "concentration of CEA-cibis-CEA_CD3 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CEA_cibis_CEA = list(
      analyte = "concentration of CEA-cibis-CEA complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CEA_cibis_CD3 = list(
      analyte = "concentration of CEA-cibis-CD3 complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CD3_cibis = list(
      analyte = "concentration of CD3-cibis complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CEA_cibis = list(
      analyte = "concentration of CEA-cibis complex",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CEA = list(
      analyte = "concentration of CEA in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    ),
    q_TCEsyn_T0_C1_CD3 = list(
      analyte = "concentration of CD3 in synapse",
      units = "molecule",
      specimen = "tumor",
      verified = FALSE
    )
  )

  population <- list(
    species = "human (in silico virtual cohort)",
    n_subjects = 2000L,
    disease_state = "microsatellite-stable (MSS) metastatic colorectal cancer",
    dose_range = paste(
      "Atezolizumab 1200 mg IV Q3W; cibisatamab 60 mg IV QW; and the",
      "combination of both regimens, each simulated for 400 days.",
      sep = " "
    ),
    notes = paste(
      "Virtual cohort of 2000 patients generated by Latin hypercube sampling of",
      "the parameter ranges in Supplementary Table S3 (tumour growth rate,",
      "initial tumour diameter 0.5-5 cm, tumour mutational burden, PD-1/PD-L1/",
      "PD-L2/CEA/CD3 expression, binding constants and cross-arm efficiencies).",
      "Patients were screened on tumour diameter, blood T cell density,",
      "activated T cell density in the tumour and the Teff/Treg ratio, leaving",
      "1312 (atezolizumab), 1325 (cibisatamab) and 1299 (combination) plausible",
      "patients. Response was classified by RECIST 1.1 at day 400; the simulated",
      "ORRs (8.2, 5.2 and 11.2 percent) were compared with NCT02324257 and",
      "NCT02650713. The packaged model carries the Table S2/S3 baseline values.",
      sep = " "
    )
  )

  ini({
    # ---- compartment capacities (SBML listOfCompartments) ----
    vol_V_C <- fixed(5); label("Central compartment volume [L]")  # SBML compartment V_C = 5 liter
    vol_V_P <- fixed(60); label("Peripheral compartment volume [L]")  # SBML compartment V_P = 60 liter
    vol_V_LN <- fixed(0.00137444678595); label("Tumour-draining lymph node volume [L]")  # SBML compartment V_LN = 1374.446786 millimeter^3
    vol_V_e <- fixed(4e-16); label("APC endosomal volume [L]")  # SBML compartment V_e = 4e-16 litre
    vol_A_e <- fixed(1.5e-09); label("APC endosomal surface area [dm^2]")  # SBML compartment A_e = 15 micrometer^2
    vol_A_s <- fixed(9e-08); label("APC surface area [dm^2]")  # SBML compartment A_s = 900 micrometer^2
    vol_syn_T1_C1 <- fixed(3.78e-09); label("Area of the Teff-cancer cell synapse [dm^2]")  # SBML compartment syn_T1_C1 = 37.8 micrometer^2
    vol_syn_T1_APC <- fixed(3.78e-09); label("Area of the Teff-APC synapse [dm^2]")  # SBML compartment syn_T1_APC = 37.8 micrometer^2
    vol_TCEsyn_T1_C1 <- fixed(7.56e-09); label("Area of the TCE synapse between Teff and cancer cell [dm^2]")  # SBML compartment TCEsyn_T1_C1 = 75.6 micrometer^2
    vol_TCEsyn_T0_C1 <- fixed(7.56e-09); label("Area of the TCE synapse between Treg and cancer cell [dm^2]")  # SBML compartment TCEsyn_T0_C1 = 75.6 micrometer^2

    # ---- model parameters (SBML listOfParameters; Tables S2/S3 where printed), in the model unit system ----
    k_cell_clear <- fixed(0.1); label("Rate of dead cell clearance from tumour compartment Rate of clearence of dead cancer cells [1/day]")  # SBML k_cell_clear = 0.1 1/day
    vol_cell <- fixed(2.52731128245e-12); label("Average volume of cancer cell Volume of a single cancer cell [L/count]")  # SBML vol_cell = 2527.311282 micrometer^3/cell
    vol_Tcell <- fixed(1.76e-13); label("Average volume of T cells Volume of a single T cell [L/count]")  # SBML vol_Tcell = 176 micrometer^3/cell
    V_Tmin <- fixed(1e-09); label("Cancer-free tumour compartment volume [L]")  # SBML V_Tmin = 0.001 microliter
    k_C1_growth <- fixed(0.005); label("Cancer cell growth rate [1/day]")  # SBML k_C1_growth = 0.005 1/day
    C_max <- fixed(4.18879020479e+14); label("Cancer cell capacity [count]")  # SBML C_max = 4.188790205e+14 cell
    k_C1_death <- fixed(1e-05); label("Cancer cell death rate from innate immune cells [1/day]")  # SBML k_C1_death = 1e-05 1/day
    k_C1_therapy <- fixed(0); label("Rate of cancer cell killing by other therapy (none in this paper) [1/day]")  # SBML k_C1_therapy = 0 1/day
    n_T1_clones <- fixed(10); label("Number of tumour-antigen-specific CD8 T cell clones (tumour mutational burden) [dimensionless]")  # Table S3 'Number of Ag Clones (TMB)' = 10 (range 0-2.5E4); SBML deposit carries 80
    Q_T1_in <- fixed(218.858981456); label("Rate of naive CD8 T cell entry into the lymph node per clone [count/day]")  # SBML Q_T1_in = 218.8589815 cell/day
    Q_T1_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN (Germain 2012) [1/day]")  # SBML Q_T1_out = 1.13 1/day
    k_T1_act <- fixed(20); label("T1 activation rate [1/day]")  # SBML k_T1_act = 20 1/day
    k_T1_pro <- fixed(1); label("T1 proliferation rate [1/day]")  # SBML k_T1_pro = 1 1/day
    k_T1_death <- fixed(0.01); label("T1 death rate [1/day]")  # SBML k_T1_death = 0.01 1/day
    q_T1_P_in <- fixed(289.44); label("rate of T1 tranport into the peripheral compartment [1/day]")  # SBML q_T1_P_in = 0.201 1/minute
    q_T1_P_out <- fixed(1); label("rate of T1 tranport out of the peripheral compartment [1/day]")  # SBML q_T1_P_out = 1 1/day
    q_T1_T_in <- fixed(4.824); label("rate of T1 tranport into the tumour compartment [1/(L*day)]")  # SBML q_T1_T_in = 3.35e-06 1/(centimeter^3*minute)
    q_T1_LN_out <- fixed(24); label("rate of T1 tranport out of the LN compartment [1/day]")  # SBML q_T1_LN_out = 24 1/day
    k_C_BTcell <- fixed(1); label("Rate of cancer cell death by TCE-activated Teff [1/day]")  # SBML k_C_BTcell = 1 1/day
    k_Treg_TCE <- fixed(2); label("Rate of TCE-activated Treg inhibition of Teff [1/day]")  # SBML k_Treg_TCE = 2 1/day
    k_C_BTreg <- fixed(1); label("Rate of cancer cell death by Teff, driven by TCE trimers in the Treg-cancer cell synapse [1/day]")  # SBML k_C_BTreg = 1 1/day
    k_IL2_deg <- fixed(288); label("rate of IL2 degradation (Lotze 1985) [1/day]")  # SBML k_IL2_deg = 0.2 1/minute
    k_IL2_cons <- fixed(86718844656); label("rate of IL2 consumption by T cells (Lotze 1985) [1/day]")  # SBML k_IL2_cons = 6e-06 nanomole/cell/hour
    k_IL2_sec <- fixed(433594223280); label("rate of IL2 secretion from T cells (Han 2012, Thurley 2015, Liu 2001) [1/day]")  # SBML k_IL2_sec = 3e-05 nanomole/cell/hour
    IL2_50 <- fixed(1.9270854368e+14); label("T cell activation half-maximal IL2 concentration (Marchingo 2014) [count/L]")  # SBML IL2_50 = 0.32 nanomolarity
    IL2_50_Treg <- fixed(1.9270854368e+12); label("Treg activation half-maximal IL2 concentration (Wang Smith 1987) [count/L]")  # SBML IL2_50_Treg = 0.0032 nanomolarity
    N0 <- fixed(5); label("Number of activated T cell generations without IL-2 [dimensionless]")  # SBML N0 = 5 dimensionless
    N_IL2 <- fixed(11); label("maximum number of activated T cell generations due to IL2 (Marchingo 2014) [dimensionless]")  # SBML N_IL2 = 11 dimensionless
    k_T1 <- fixed(0.1); label("Rate of Teff exhaustion by cancer cells [1/day]")  # SBML k_T1 = 0.1 1/day
    k_C_T1 <- fixed(3); label("Rate of cancer cell death by Teff [1/day]")  # SBML k_C_T1 = 3 1/day
    k_Treg <- fixed(1); label("Rate of Teff inhibition by Treg [1/day]")  # Table S3 'Rate of Teff inhibition by Treg' = 1 1/day (range 0.1-1); SBML deposit carries 2
    a0 <- fixed(0.692556926728); label("Coefficient a0 of the polynomial m(N_aT) in the T cell activation fraction [dimensionless]")  # SBML a0 = 0.6925569267 dimensionless
    a1 <- fixed(-0.0961148357298); label("Coefficient a1 of the polynomial m(N_aT) in the T cell activation fraction [dimensionless]")  # SBML a1 = -0.09611483573 dimensionless
    a2 <- fixed(0.00863812431352); label("Coefficient a2 of the polynomial m(N_aT) in the T cell activation fraction [dimensionless]")  # SBML a2 = 0.008638124314 dimensionless
    a3 <- fixed(-0.000403990701606); label("Coefficient a3 of the polynomial m(N_aT) in the T cell activation fraction [dimensionless]")  # SBML a3 = -0.0004039907016 dimensionless
    a4 <- fixed(7.59161205e-06); label("Coefficient a4 of the polynomial m(N_aT) in the T cell activation fraction [dimensionless]")  # SBML a4 = 7.59161205e-06 dimensionless
    n_T0_clones <- fixed(100); label("Number of self-antigen-specific Treg clones [dimensionless]")  # SBML n_T0_clones = 100 dimensionless
    Q_T0_in <- fixed(356.645243594); label("Rate of naive Treg entry into the lymph node per clone [count/day]")  # SBML Q_T0_in = 356.6452436 cell/day
    Q_T0_out <- fixed(1.13); label("Rate of naive T cell transport out of the LN (Germain 2012) [1/day]")  # SBML Q_T0_out = 1.13 1/day
    k_T0_act <- fixed(20); label("T0 activation rate [1/day]")  # SBML k_T0_act = 20 1/day
    k_T0_pro <- fixed(1); label("T0 proliferation rate [1/day]")  # SBML k_T0_pro = 1 1/day
    k_T0_death <- fixed(0.01); label("T0 death rate [1/day]")  # SBML k_T0_death = 0.01 1/day
    q_T0_P_in <- fixed(289.44); label("rate of T0 tranport into the peripheral compartment [1/day]")  # SBML q_T0_P_in = 0.201 1/minute
    q_T0_P_out <- fixed(0.015); label("rate of T0 tranport out of the peripheral compartment [1/day]")  # SBML q_T0_P_out = 0.015 1/day
    q_T0_T_in <- fixed(4.824); label("rate of T0 tranport into the tumour compartment [1/(L*day)]")  # SBML q_T0_T_in = 3.35e-06 1/(centimeter^3*minute)
    q_T0_LN_out <- fixed(24); label("rate of T0 tranport out of the LN compartment [1/day]")  # SBML q_T0_LN_out = 24 1/day
    k_APC_mat <- fixed(1.5); label("Maximum rate of APC maturation [1/day]")  # SBML k_APC_mat = 1.5 1/day
    k_APC_mig <- fixed(4); label("Rate of APC migration [1/day]")  # SBML k_APC_mig = 4 1/day
    k_APC_death <- fixed(0.01); label("Rate of APC death [1/day]")  # SBML k_APC_death = 0.01 1/day
    k_mAPC_death <- fixed(0.02); label("Rate of mAPC death [1/day]")  # SBML k_mAPC_death = 0.02 1/day
    APC0_T <- fixed(400000000); label("APC density in the tumour [count/L]")  # SBML APC0_T = 400000 cell/milliliter
    APC0_LN <- fixed(1200000000); label("APC density in the LN [count/L]")  # SBML APC0_LN = 1200000 cell/milliliter
    k_c <- fixed(2); label("Cytokine rate constant [1/day]")  # SBML k_c = 2 1/day
    c0 <- fixed(6.02214199e+14); label("Baseline cytokine concentration [count/L]")  # SBML c0 = 1e-09 molarity
    c50 <- fixed(6.02214199e+14); label("Cytokine concentration for half-maximal APC maturation [count/L]")  # SBML c50 = 1e-09 molarity
    DAMPs <- fixed(8069670266.6); label("Cytokine released per dying cancer cell [dimensionless]")  # SBML DAMPs = 1.34e-14 mole/cell
    n_sites_APC <- fixed(10); label("Maxium number of T Cells an APC can interact with Maximum Number of T Cells an APC can Interact with [dimensionless]")  # SBML n_sites_APC = 10 dimensionless
    kin <- fixed(14.4); label("Rate of MHC internalization [1/day]")  # SBML kin = 14.4 1/day
    kout <- fixed(28.8); label("Rate of MHC externalization [1/day]")  # SBML kout = 28.8 1/day
    k_P0_up <- fixed(14.4); label("Rate of antigen uptake by APCs Rate of Antigen Uptake [1/(count*day)]")  # SBML k_P0_up = 14.4 1/day/cell
    k_xP0_deg <- fixed(2); label("Rate of extracellular antigen degradation Rate of Antigen Degradation in the Extracellular Space [1/day]")  # SBML k_xP0_deg = 2 1/day
    k_P0_deg <- fixed(17.28); label("Rate of endosomal antigen degradation Rate of Antigen Degradation in Endosomes [1/day]")  # SBML k_P0_deg = 17.28 1/day
    k_p0_deg <- fixed(144); label("Rate of endosomal epitope degradation Rate of Epitope Degradation [1/day]")  # SBML k_p0_deg = 144 1/day
    k_P0_on <- fixed(2.39117576834e-19); label("Rate of antigen- MHC binding Antigen Association Constant [L/(count*day)]")  # SBML k_P0_on = 144000 1/day/molarity
    k_P0_d1 <- fixed(6.02214199e+16); label("Antigen-MHC kd [count/L]")  # SBML k_P0_d1 = 1e-07 molarity
    p0_50 <- fixed(264550.26455); label("MHC-Epitope-TCR Concentration for Half-Maximal T Cell Activation [count/dm^2]")  # SBML p0_50 = 2.645502646e-05 molecule/micrometer^2
    P0_C1 <- fixed(6.02214199e+15); label("Concentration of P0 in C1 [1/L]")  # SBML P0_C1 = 1e-08 molarity/cell
    A_Tcell <- fixed(1.5131041193e-08); label("Surface area of the T cell T Cell Surface Area [dm^2]")  # SBML A_Tcell = 151.3104119 micrometer^2
    A_cell <- fixed(8.97270277792e-08); label("Surface area of the Cancer cell Cancer Cell Surface Area [dm^2]")  # SBML A_cell = 897.2702778 micrometer^2
    A_APC <- fixed(9e-08); label("Surface area of the APC APC Surface Area [dm^2]")  # SBML A_APC = 900 micrometer^2
    k_M1p0_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex Rate of TCR Binding [dm^2/(count*day)]")  # SBML k_M1p0_TCR_on = 1 1/(second*molecule/micrometer^2)
    k_M1p0_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex Rate of TCR Unbinding [1/day]")  # SBML k_M1p0_TCR_off = 1 1/second
    TCR_p0_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell Number of TCR Molecules on Naive T Cells per Unit ... [count/dm^2]")  # SBML TCR_p0_tot = 103.8130807 molecule/micrometer^2
    k_P1_up <- fixed(14.4); label("Rate of antigen uptake by APCs Rate of Antigen Uptake [1/(count*day)]")  # SBML k_P1_up = 14.4 1/day/cell
    k_xP1_deg <- fixed(2); label("Rate of extracellular antigen degradation Rate of Antigen Degradation in the Extracellular Space [1/day]")  # SBML k_xP1_deg = 2 1/day
    k_P1_deg <- fixed(17.28); label("Rate of endosomal antigen degradation Rate of Antigen Degradation in Endosomes [1/day]")  # SBML k_P1_deg = 17.28 1/day
    k_p1_deg <- fixed(144); label("Rate of endosomal epitope degradation Rate of Epitope Degradation [1/day]")  # SBML k_p1_deg = 144 1/day
    k_P1_on <- fixed(2.39117576834e-19); label("Rate of antigen- MHC binding Antigen Association Constant [L/(count*day)]")  # SBML k_P1_on = 144000 1/day/molarity
    k_P1_d1 <- fixed(2.408856796e+16); label("Antigen-MHC kd [count/L]")  # SBML k_P1_d1 = 4e-08 molarity
    p1_50 <- fixed(264550.26455); label("MHC-Epitope-TCR Concentration for Half-Maximal T Cell Activation [count/dm^2]")  # SBML p1_50 = 2.645502646e-05 molecule/micrometer^2
    P1_C1 <- fixed(6.02214199e+15); label("Concentration of P1 in C1 [1/L]")  # SBML P1_C1 = 1e-08 molarity/cell
    k_M1p1_TCR_on <- fixed(8.64e-06); label("Rate of TCR binding to MHC-peptide complex Rate of TCR Binding [dm^2/(count*day)]")  # SBML k_M1p1_TCR_on = 1 1/(second*molecule/micrometer^2)
    k_M1p1_TCR_off <- fixed(86400); label("Rate of TCR unbinding from MHC-peptide complex Rate of TCR Unbinding [1/day]")  # SBML k_M1p1_TCR_off = 1 1/second
    k_M1p1_TCR_p <- fixed(86400); label("Rate of MHC-peptide-TCR complex modification [1/day]")  # SBML k_M1p1_TCR_p = 1 1/second
    phi_M1p1_TCR <- fixed(7776); label("Rate at which a maximally modified MHC-peptide-TCR complex becomes non-signalling [1/day]")  # SBML phi_M1p1_TCR = 0.09 1/second
    N_M1p1_TCR <- fixed(10); label("Number of modifcation steps for MHC-peptide-TCR complex Number of Intermediate Steps [dimensionless]")  # SBML N_M1p1_TCR = 10 dimensionless
    TCR_p1_tot <- fixed(1.0381308067e+12); label("Total number of TCR molecules per naive T cell Number of TCR Molecules on Naive T Cells per Unit ... [count/dm^2]")  # SBML TCR_p1_tot = 103.8130807 molecule/micrometer^2
    kon_PD1_PDL1 <- fixed(8.3691151892e-07); label("2D kon of PD-1-PD-L1 binding [dm^2/(count*day)]")  # SBML kon_PD1_PDL1 = 0.05833333333 1/(micromolarity*nanometer*second)
    aPD1 <- fixed(0); label("Anti-PD-1 antibody concentration (not dosed) [count/L]")  # SBML aPD1 = 0 molarity
    q_P_atezolizumab <- fixed(6.0285045744e-06); label("Atezolizumab central-peripheral volumetric flow [L/day]")  # SBML q_P_atezolizumab = 6.97743585e-08 milliliter/second
    q_T_atezolizumab <- fixed(0.736128); label("Atezolizumab central-tumour transport rate [1/day]")  # SBML q_T_atezolizumab = 8.52e-06 1/second
    q_LN_atezolizumab <- fixed(6.0285045744e-06); label("Atezolizumab central-lymph node volumetric flow [L/day]")  # SBML q_LN_atezolizumab = 6.97743585e-08 milliliter/second
    q_LD_atezolizumab <- fixed(2.16); label("Atezolizumab lymphatic drainage rate [1/day]")  # SBML q_LD_atezolizumab = 0.0015 1/minute
    k_cl_atezolizumab <- fixed(0.072); label("Atezolizumab elimination rate from the central compartment [1/day]")  # SBML k_cl_atezolizumab = 0.003 1/hour
    gamma_C_atezolizumab <- fixed(0.3643); label("Fraction of the central volume available to atezolizumab [dimensionless]")  # SBML gamma_C_atezolizumab = 0.3643 dimensionless
    gamma_P_atezolizumab <- fixed(0.0405); label("Fraction of the peripheral volume available to atezolizumab [dimensionless]")  # SBML gamma_P_atezolizumab = 0.0405 dimensionless
    gamma_T_atezolizumab <- fixed(0.522); label("Fraction of the tumour volume available to atezolizumab [dimensionless]")  # SBML gamma_T_atezolizumab = 0.522 dimensionless
    gamma_LN_atezolizumab <- fixed(0.2); label("Fraction of the lymph node volume available to atezolizumab [dimensionless]")  # SBML gamma_LN_atezolizumab = 0.2 dimensionless
    kon_PD1_PDL2 <- fixed(1.09994085344e-06); label("2D kon of PD-1-PD-L2 binding [dm^2/(count*day)]")  # SBML kon_PD1_PDL2 = 0.07666666667 1/(micromolarity*nanometer*second)
    kon_PD1_aPD1 <- fixed(1.86511709931e-13); label("kon of PD-1-anti-PD-1 binding (anti-PD-1 not dosed) [L/(count*day)]")  # SBML kon_PD1_aPD1 = 1300000 1/(molarity*second)
    kon_PDL1_aPDL1 <- fixed(6.16923348232e-14); label("kon of PD-L1-atezolizumab binding [L/(count*day)]")  # SBML kon_PDL1_aPDL1 = 430000 1/(molarity*second)
    koff_PD1_PDL1 <- fixed(123984); label("koff of PD-1-PD-L1 binding [1/day]")  # SBML koff_PD1_PDL1 = 1.435 1/second
    koff_PD1_PDL2 <- fixed(45705.6); label("koff of PD-1-PD-L2 binding [1/day]")  # SBML koff_PD1_PDL2 = 0.529 1/second
    koff_PD1_aPD1 <- fixed(292.032); label("koff of PD-1-anti-PD-1 binding (anti-PD-1 not dosed) [1/day]")  # SBML koff_PD1_aPD1 = 0.00338 1/second
    koff_PDL1_aPDL1 <- fixed(14.8608); label("koff of PD-L1-atezolizumab binding [1/day]")  # SBML koff_PDL1_aPDL1 = 0.000172 1/second
    Chi_PD1 <- fixed(333333333.333); label("Anti-PD-1 cross-arm binding efficiency divided by the synapse gap (anti-PD-1 not dosed) [1/dm]")  # SBML Chi_PD1 = 3.333333333 1/nanometer
    Chi_PDL1 <- fixed(3333333333.33); label("Atezolizumab cross-arm binding efficiency chi divided by the synapse gap [1/dm]")  # SBML Chi_PDL1 = 33.33333333 1/nanometer
    n_PD1 <- fixed(2); label("Hill coefficient for PD-1 complex inhibition of T cells [dimensionless]")  # SBML n_PD1 = 2 dimensionless
    T_PD1_total <- fixed(60000); label("Number of PD-1 molecules per Teff cell [count]")  # SBML T_PD1_total = 60000 molecule
    C_PDL1_total <- fixed(160000); label("Number of PD-L1 molecules per cancer cell [count]")  # Table S2 C_PDL1_total = 160000 sites/cell, Table S3 'Total PD-L1 on tumor cell' = 160000 (range 1-160000); SBML deposit carries 400000
    APC_PDL1_total <- fixed(1600000); label("Number of PD-L1 molecules per APC [count]")  # SBML APC_PDL1_total = 1600000 molecule
    r_PDL2C <- fixed(0.065); label("Ratio of PD-L2 to PD-L1 molecules on cancer cells [dimensionless]")  # Table S2 C_PDL2_total / C_PDL1_total = 10400 / 160000 sites/cell; Table S3 prints a 0.1 baseline outside its own 0-0.07 range; SBML deposit carries 0.05
    r_PDL2APC <- fixed(0.065); label("Ratio of PD-L2 to PD-L1 molecules on APCs [dimensionless]")  # Table S2 APC_PDL2_total / APC_PDL1_total = 104000 / 1600000 sites/cell; Table S3 prints a 0.1 baseline outside its own 0-0.07 range; SBML deposit carries 0.05
    q_P_cibis <- fixed(2.38032e-05); label("Cibisatamab central-peripheral volumetric flow [L/day]")  # SBML q_P_cibis = 2.755e-07 milliliter/second
    q_T_cibis <- fixed(0.736128); label("Cibisatamab central-tumour transport rate [1/day]")  # SBML q_T_cibis = 8.52e-06 1/second
    q_LN_cibis <- fixed(2.38032e-05); label("Cibisatamab central-lymph node volumetric flow [L/day]")  # SBML q_LN_cibis = 2.755e-07 milliliter/second
    q_LD_cibis <- fixed(2.16); label("Cibisatamab lymphatic drainage rate [1/day]")  # SBML q_LD_cibis = 0.0015 1/minute
    k_cl_cibis <- fixed(0.6); label("Cibisatamab elimination rate from the central compartment [1/day]")  # SBML k_cl_cibis = 0.025 1/hour
    gamma_C_cibis <- fixed(0.46); label("Fraction of the central volume available to cibisatamab [dimensionless]")  # SBML gamma_C_cibis = 0.46 dimensionless
    gamma_P_cibis <- fixed(0.067); label("Fraction of the peripheral volume available to cibisatamab [dimensionless]")  # SBML gamma_P_cibis = 0.067 dimensionless
    gamma_T_cibis <- fixed(0.522); label("Fraction of the tumour volume available to cibisatamab [dimensionless]")  # SBML gamma_T_cibis = 0.522 dimensionless
    gamma_LN_cibis <- fixed(0.2); label("Fraction of the lymph node volume available to cibisatamab [dimensionless]")  # SBML gamma_LN_cibis = 0.2 dimensionless
    dT_syn <- fixed(3e-08); label("Gap between the two cell membranes in the TCE synapse [dm]")  # SBML dT_syn = 3 nanometer
    SA_Tcell <- fixed(1.5131041193e-08); label("T cell surface area [dm^2]")  # SBML SA_Tcell = 151.3104119 micrometer^2
    SA_cell <- fixed(1.0636176088e-07); label("Cancer cell surface area [dm^2]")  # SBML SA_cell = 1063.617609 micrometer^2
    kon_CD3_cibis <- fixed(1.434705461e-15); label("kon of CD3-cibisatamab binding [L/(count*day)]")  # SBML kon_CD3_cibis = 10000 1/(molarity*second)
    kon_CEA_cibis <- fixed(1.434705461e-16); label("kon of CEA-cibisatamab binding [L/(count*day)]")  # SBML kon_CEA_cibis = 1000 1/(molarity*second)
    koff_CD3_cibis <- fixed(64.8); label("koff of CD3-cibisatamab binding [1/day]")  # SBML koff_CD3_cibis = 0.00075 1/second
    koff_CEA_cibis <- fixed(11.232); label("koff of CEA-cibisatamab binding [1/day]")  # SBML koff_CEA_cibis = 0.00013 1/second
    Chi_CEA_cibis <- fixed(33333333333.3); label("Cibisatamab cross-arm binding efficiency lambda divided by the synapse gap [1/dm]")  # Table S3 lambda = 1000 (range 0.001-1e5), divided by the synapse gap dT_syn = 3 nm as in the SBML rate laws; SBML deposit carries 0.3333333 1/nm, i.e. lambda = 1
    CEA_CD3_50 <- fixed(80000000000); label("Density of CEA-cibisatamab-CD3 trimers for half-maximal TCE-driven killing [count/dm^2]")  # SBML CEA_CD3_50 = 8 molecule/micrometer^2
    n_CEA_CD3 <- fixed(3); label("Hill coefficient for TCE-driven killing [dimensionless]")  # SBML n_CEA_CD3 = 3 dimensionless
    T1_CD3_total <- fixed(61000); label("Number of CD3 molecules per Teff cell [count]")  # Table S3 CD3 expression = 61000 molecule (range 30000-90000; also Ma 2020a Table S2 Teff_CD3_total); SBML deposit carries 62000
    C1_CEA_total <- fixed(20000); label("Number of CEA molecules per cancer cell [count]")  # SBML C1_CEA_total = 20000 molecule
    T0_CD3_total <- fixed(61000); label("Number of CD3 molecules per Treg cell [count]")  # Table S3 CD3 expression = 61000 molecule (range 30000-90000; also Ma 2020a Table S2 Treg_CD3_total); SBML deposit carries 62000
    K_PD1_PDLX <- fixed(250); label("Number of PD-1/PD-L1 (PD-L2) complexes for half-maximal inhibition of T cell killing [count]")  # Table S2 K_PD1_PDLX_CC = K_PD1_PDLX_APC = 250 molecule; SBML PD1_50 = 6 molecule/micrometer^2 over the 37.8 micrometer^2 synapse = 226.8
    delta_PDL2 <- fixed(1); label("Weight of PD-1/PD-L2 complexes in the checkpoint Hill functions [dimensionless]")  # Supplementary Eqs 10-11 (range 0-1); the SBML Hill rules carry no weight, i.e. delta = 1
  })
  model({
    # molecules per nmol (1 molecule = 1.66053872801495e-24 mol in the export)
    nmol_to_molecule <- 602214199000000

    # ---- parameter aliases (rxode2 reads a bare theta +/- theta as a
    # mu-reference block; these mechanistic rate laws legitimately add two
    # fixed constants, so each such parameter is copied to a local first) ----
    a_k_C1_death <- k_C1_death
    a_k_C1_therapy <- k_C1_therapy
    a_k_M1p1_TCR_off <- k_M1p1_TCR_off
    a_k_M1p1_TCR_p <- k_M1p1_TCR_p
    a_phi_M1p1_TCR <- phi_M1p1_TCR

    # ---- initial-assignment parameters ----
    k_ratio_1 <- k_T1_pro/k_T1_act
    k_ratio_0 <- k_T0_pro/k_T0_act

    # ---- compartment capacities ----
    m_V_C <- vol_V_C
    m_V_P <- vol_V_P
    m_V_LN <- vol_V_LN
    m_V_e <- vol_V_e
    m_A_e <- vol_A_e
    m_A_s <- vol_A_s
    m_syn_T1_C1 <- vol_syn_T1_C1
    m_syn_T1_APC <- vol_syn_T1_APC
    m_TCEsyn_T1_C1 <- vol_TCEsyn_T1_C1
    m_TCEsyn_T0_C1 <- vol_TCEsyn_T0_C1

    # ---- species quantities: states are AMOUNTS (counts, or nmol for the
    # two antibodies); x_* are the SimBiology species values in model units ----
    x_V_C_T1 <- q_V_C_T1
    x_V_C_T0 <- q_V_C_T0
    x_V_P_T1 <- q_V_P_T1
    x_V_P_T0 <- q_V_P_T0
    x_V_T_C_x <- q_V_T_C_x
    x_V_T_T_exh <- q_V_T_T_exh
    x_V_T_C1 <- q_V_T_C1
    x_V_T_T1 <- q_V_T_T1
    x_V_T_T0 <- q_V_T_T0
    x_V_T_APC <- q_V_T_APC
    x_V_T_mAPC <- q_V_T_mAPC
    x_V_LN_nT1 <- q_V_LN_nT1
    x_V_LN_aT1 <- q_V_LN_aT1
    x_V_LN_T1 <- q_V_LN_T1
    x_V_LN_nT0 <- q_V_LN_nT0
    x_V_LN_aT0 <- q_V_LN_aT0
    x_V_LN_T0 <- q_V_LN_T0
    x_V_LN_APC <- q_V_LN_APC
    x_V_LN_mAPC <- q_V_LN_mAPC
    vol_V_T <- V_Tmin+(vol_cell*x_V_T_C_x)+(vol_Tcell*x_V_T_T_exh)+(vol_cell*x_V_T_C1)+(vol_Tcell*x_V_T_T1)+(vol_Tcell*x_V_T_T0)   # tumour volume rule
    x_V_C_atezolizumab <- q_V_C_atezolizumab*nmol_to_molecule/m_V_C
    x_V_C_cibis <- q_V_C_cibis*nmol_to_molecule/m_V_C
    x_V_P_atezolizumab <- q_V_P_atezolizumab*nmol_to_molecule/m_V_P
    x_V_P_cibis <- q_V_P_cibis*nmol_to_molecule/m_V_P
    x_V_T_c <- q_V_T_c/vol_V_T
    x_V_T_atezolizumab <- q_V_T_atezolizumab*nmol_to_molecule/vol_V_T
    x_V_T_cibis <- q_V_T_cibis*nmol_to_molecule/vol_V_T
    x_V_LN_IL2 <- q_V_LN_IL2/m_V_LN
    x_V_LN_P0 <- q_V_LN_P0/m_V_LN
    x_V_LN_P1 <- q_V_LN_P1/m_V_LN
    x_V_LN_atezolizumab <- q_V_LN_atezolizumab*nmol_to_molecule/m_V_LN
    x_V_LN_cibis <- q_V_LN_cibis*nmol_to_molecule/m_V_LN
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
    x_syn_T1_C1_PD1_PDL1 <- q_syn_T1_C1_PD1_PDL1/m_syn_T1_C1
    x_syn_T1_C1_PD1_PDL2 <- q_syn_T1_C1_PD1_PDL2/m_syn_T1_C1
    x_syn_T1_C1_PD1 <- q_syn_T1_C1_PD1/m_syn_T1_C1
    x_syn_T1_C1_PDL1 <- q_syn_T1_C1_PDL1/m_syn_T1_C1
    x_syn_T1_C1_PDL2 <- q_syn_T1_C1_PDL2/m_syn_T1_C1
    x_syn_T1_C1_PD1_aPD1 <- q_syn_T1_C1_PD1_aPD1/m_syn_T1_C1
    x_syn_T1_C1_PD1_aPD1_PD1 <- q_syn_T1_C1_PD1_aPD1_PD1/m_syn_T1_C1
    x_syn_T1_C1_PDL1_aPDL1 <- q_syn_T1_C1_PDL1_aPDL1/m_syn_T1_C1
    x_syn_T1_C1_PDL1_aPDL1_PDL1 <- q_syn_T1_C1_PDL1_aPDL1_PDL1/m_syn_T1_C1
    x_syn_T1_APC_PD1_PDL1 <- q_syn_T1_APC_PD1_PDL1/m_syn_T1_APC
    x_syn_T1_APC_PD1_PDL2 <- q_syn_T1_APC_PD1_PDL2/m_syn_T1_APC
    x_syn_T1_APC_PD1 <- q_syn_T1_APC_PD1/m_syn_T1_APC
    x_syn_T1_APC_PDL1 <- q_syn_T1_APC_PDL1/m_syn_T1_APC
    x_syn_T1_APC_PDL2 <- q_syn_T1_APC_PDL2/m_syn_T1_APC
    x_syn_T1_APC_PD1_aPD1 <- q_syn_T1_APC_PD1_aPD1/m_syn_T1_APC
    x_syn_T1_APC_PD1_aPD1_PD1 <- q_syn_T1_APC_PD1_aPD1_PD1/m_syn_T1_APC
    x_syn_T1_APC_PDL1_aPDL1 <- q_syn_T1_APC_PDL1_aPDL1/m_syn_T1_APC
    x_syn_T1_APC_PDL1_aPDL1_PDL1 <- q_syn_T1_APC_PDL1_aPDL1_PDL1/m_syn_T1_APC
    x_TCEsyn_T1_C1_CEA_cibis_CEA_CD3 <- q_TCEsyn_T1_C1_CEA_cibis_CEA_CD3/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CEA_cibis_CEA <- q_TCEsyn_T1_C1_CEA_cibis_CEA/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CEA_cibis_CD3 <- q_TCEsyn_T1_C1_CEA_cibis_CD3/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CD3_cibis <- q_TCEsyn_T1_C1_CD3_cibis/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CEA_cibis <- q_TCEsyn_T1_C1_CEA_cibis/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CEA <- q_TCEsyn_T1_C1_CEA/m_TCEsyn_T1_C1
    x_TCEsyn_T1_C1_CD3 <- q_TCEsyn_T1_C1_CD3/m_TCEsyn_T1_C1
    x_TCEsyn_T0_C1_CEA_cibis_CEA_CD3 <- q_TCEsyn_T0_C1_CEA_cibis_CEA_CD3/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CEA_cibis_CEA <- q_TCEsyn_T0_C1_CEA_cibis_CEA/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CEA_cibis_CD3 <- q_TCEsyn_T0_C1_CEA_cibis_CD3/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CD3_cibis <- q_TCEsyn_T0_C1_CD3_cibis/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CEA_cibis <- q_TCEsyn_T0_C1_CEA_cibis/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CEA <- q_TCEsyn_T0_C1_CEA/m_TCEsyn_T0_C1
    x_TCEsyn_T0_C1_CD3 <- q_TCEsyn_T0_C1_CD3/m_TCEsyn_T0_C1

    # ---- repeated-assignment rules ----
    C_total <- x_V_T_C1
    T_total <- x_V_T_T1+x_V_T_T0
    T_total_LN <- x_V_LN_T1
    N_aT <- N0+((N_IL2*x_V_LN_IL2)/(IL2_50+x_V_LN_IL2))
    Tregs_ <- x_V_T_T0
    H_APC <- (n_sites_APC*x_V_LN_APC)/((n_sites_APC*x_V_LN_APC)+T_total_LN+1)
    H_mAPC <- (n_sites_APC*x_V_LN_mAPC)/((n_sites_APC*x_V_LN_mAPC)+T_total_LN+1)
    pTCR_p0_MHC_tot <- 0.5*(((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))-(TCR_p0_tot*sqrt((((((x_A_s_M1p0/n_T0_clones)+TCR_p0_tot+(k_M1p0_TCR_off/k_M1p0_TCR_on))/TCR_p0_tot))^2-(((4*x_A_s_M1p0)/n_T0_clones)/TCR_p0_tot)))))
    H_P0 <- pTCR_p0_MHC_tot/(pTCR_p0_MHC_tot+p0_50)
    pTCR_p1_MHC_tot <- (a_k_M1p1_TCR_off/(a_k_M1p1_TCR_off+a_phi_M1p1_TCR))*((a_k_M1p1_TCR_p/(a_k_M1p1_TCR_off+a_k_M1p1_TCR_p)))^N_M1p1_TCR*0.5*(((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))-(TCR_p1_tot*sqrt((((((x_A_s_M1p1/n_T1_clones)+TCR_p1_tot+(a_k_M1p1_TCR_off/k_M1p1_TCR_on))/TCR_p1_tot))^2-(((4*x_A_s_M1p1)/n_T1_clones)/TCR_p1_tot)))))
    H_P1 <- pTCR_p1_MHC_tot/(pTCR_p1_MHC_tot+p1_50)
    aPDL1 <- x_V_T_atezolizumab
    PD1_PDLX_syn_T1_C1 <- q_syn_T1_C1_PD1_PDL1 + delta_PDL2*q_syn_T1_C1_PD1_PDL2   # Supplementary Eqs 8 and 10: complexes in the synapse
    H_PD1_C1 <- (PD1_PDLX_syn_T1_C1/K_PD1_PDLX)^n_PD1/((PD1_PDLX_syn_T1_C1/K_PD1_PDLX)^n_PD1 + 1)
    PD1_PDLX_syn_T1_APC <- q_syn_T1_APC_PD1_PDL1 + delta_PDL2*q_syn_T1_APC_PD1_PDL2   # Supplementary Eqs 9 and 11: complexes in the synapse
    H_PD1_APC <- (PD1_PDLX_syn_T1_APC/K_PD1_PDLX)^n_PD1/((PD1_PDLX_syn_T1_APC/K_PD1_PDLX)^n_PD1 + 1)
    H_CEA_C1_T1 <- ((x_TCEsyn_T1_C1_CEA_cibis_CEA_CD3/CEA_CD3_50))^n_CEA_CD3/(((x_TCEsyn_T1_C1_CEA_cibis_CEA_CD3/CEA_CD3_50))^n_CEA_CD3+1)
    H_CEA_C1_T0 <- ((x_TCEsyn_T0_C1_CEA_cibis_CEA_CD3/CEA_CD3_50))^n_CEA_CD3/(((x_TCEsyn_T0_C1_CEA_cibis_CEA_CD3/CEA_CD3_50))^n_CEA_CD3+1)
    R_Tcell <- (((k_C_T1*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*(1-H_PD1_C1))+(((k_C_BTcell*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*H_CEA_C1_T1)+(((k_C_BTcell*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*H_CEA_C1_T0)
    m <- a0+(a1*N_aT)+(a2*(N_aT)^2)+(a3*(N_aT)^3)+(a4*(N_aT)^4)
    k <- -log((((exp(N_aT)-1)/((2*exp(N_aT))-1))/(1-exp((-N_aT)))))
    F_T1_act <- x_V_LN_aT1*(1-((1-exp((-N_aT)))*exp(((-k)*((((k_ratio_1*x_V_LN_aT1)/((H_mAPC*x_V_LN_nT1)+1e-09))/((2*exp(N_aT))-1)))^m))))
    F_T0_act <- x_V_LN_aT0*(1-((1-exp((-N_aT)))*exp(((-k)*((((k_ratio_0*x_V_LN_aT0)/((H_APC*x_V_LN_nT0)+1e-09))/((2*exp(N_aT))-1)))^m))))

    # SBML event: cancer cells are set to zero once fewer than 0.9 cells
    # remain. rxode2 has no state-reset events, so proliferation is switched
    # off below that threshold instead and the remaining fraction of a cell
    # decays through the death terms.
    c1_alive <- 1*(x_V_T_C1 >= 0.9)

    # ---- reaction rates (SBML listOfReactions, numbered in order), amount per day ----
    v1 <- k_cell_clear*x_V_T_C_x   # R1: V_T.C_x -> null
    v2 <- k_cell_clear*x_V_T_T_exh   # R2: V_T.T_exh -> null
    v3 <- (k_C1_growth*x_V_T_C1*(1-(C_total/C_max)))*c1_alive   # R3: null -> V_T.C1
    v4 <- a_k_C1_death*x_V_T_C1   # R4: V_T.C1 -> V_T.C_x
    v5 <- Q_T1_in*n_T1_clones   # R5: null -> V_LN.nT1
    v6 <- Q_T1_out*x_V_LN_nT1   # R6: V_LN.nT1 -> null
    v7 <- k_T1_act*H_mAPC*H_P1*x_V_LN_nT1   # R7: V_LN.nT1 -> V_LN.aT1
    v8 <- (k_T1_pro/N_aT)*x_V_LN_aT1   # R8: V_LN.aT1 -> null
    v9 <- (k_T1_pro/N_aT)*2*exp((N_aT-1))*x_V_LN_aT1*(1-H_PD1_APC)   # R9: null -> V_LN.T1
    v10 <- k_T1_death*x_V_C_T1   # R10: V_C.T1 -> null
    v11 <- k_T1_death*x_V_P_T1   # R11: V_P.T1 -> null
    v12 <- k_T1_death*x_V_T_T1   # R12: V_T.T1 -> V_T.T_exh
    v13 <- k_T1_death*x_V_LN_T1   # R13: V_LN.T1 -> null
    v14 <- (k_Treg*x_V_T_T1*Tregs_)/(C_total+T_total+1)   # R14: V_T.T1 -> V_T.T_exh
    v15 <- ((k_Treg_TCE*x_V_T_T1*Tregs_)/(C_total+T_total+1))*H_CEA_C1_T0   # R15: V_T.T1 -> V_T.T_exh
    v16 <- ((k_T1*x_V_T_T1*C_total)/(C_total+T_total+1))*H_PD1_C1   # R16: V_T.T1 -> V_T.T_exh
    v17 <- q_T1_P_in*x_V_C_T1   # R17: V_C.T1 -> V_P.T1
    v18 <- q_T1_P_out*x_V_P_T1   # R18: V_P.T1 -> V_C.T1
    v19 <- q_T1_T_in*vol_V_T*x_V_C_T1   # R19: V_C.T1 -> V_T.T1
    v20 <- q_T1_LN_out*x_V_LN_T1   # R20: V_LN.T1 -> V_C.T1
    v21 <- k_IL2_deg*x_V_LN_IL2*m_V_LN   # R21: V_LN.IL2 -> null
    v22 <- (k_IL2_cons*T_total_LN*x_V_LN_IL2)/(IL2_50+x_V_LN_IL2)   # R22: V_LN.IL2 -> null
    v23 <- k_IL2_sec*x_V_LN_aT1   # R23: null -> V_LN.IL2
    v24 <- (((k_C_T1*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*(1-H_PD1_C1))+(((k_C_BTcell*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*H_CEA_C1_T1)+(((k_C_BTreg*x_V_T_T1*x_V_T_C1)/(C_total+T_total+1))*H_CEA_C1_T0)   # R24: V_T.C1 -> V_T.C_x
    v25 <- Q_T0_in*n_T0_clones   # R25: null -> V_LN.nT0
    v26 <- Q_T0_out*x_V_LN_nT0   # R26: V_LN.nT0 -> null
    v27 <- k_T0_act*H_APC*H_P0*x_V_LN_nT0   # R27: V_LN.nT0 -> V_LN.aT0
    v28 <- (k_T0_pro/N_aT)*x_V_LN_aT0   # R28: V_LN.aT0 -> null
    v29 <- (k_T0_pro/N_aT)*2*exp((N_aT-1))*x_V_LN_aT0*(1-H_PD1_APC)   # R29: null -> V_LN.T0
    v30 <- k_T0_death*x_V_C_T0   # R30: V_C.T0 -> null
    v31 <- k_T0_death*x_V_P_T0   # R31: V_P.T0 -> null
    v32 <- k_T0_death*x_V_T_T0   # R32: V_T.T0 -> V_T.T_exh
    v33 <- k_T0_death*x_V_LN_T0   # R33: V_LN.T0 -> null
    v34 <- q_T0_P_in*x_V_C_T0   # R34: V_C.T0 -> V_P.T0
    v35 <- q_T0_P_out*x_V_P_T0   # R35: V_P.T0 -> V_C.T0
    v36 <- q_T0_T_in*vol_V_T*x_V_C_T0   # R36: V_C.T0 -> V_T.T0
    v37 <- q_T0_LN_out*x_V_LN_T0   # R37: V_LN.T0 -> V_C.T0
    v38 <- (k_IL2_cons*x_V_LN_T0*x_V_LN_IL2)/(IL2_50_Treg+x_V_LN_IL2)   # R38: V_LN.IL2 -> null
    v39 <- k_APC_death*((APC0_T*vol_V_T)-x_V_T_APC)   # R39: null -> V_T.APC
    v40 <- k_APC_death*((APC0_LN*m_V_LN)-x_V_LN_APC)   # R40: null -> V_LN.APC
    v41 <- ((k_APC_mat*x_V_T_c)/(x_V_T_c+c50))*x_V_T_APC   # R41: V_T.APC -> V_T.mAPC
    v42 <- k_APC_mig*x_V_T_mAPC   # R42: V_T.mAPC -> V_LN.mAPC
    v43 <- k_mAPC_death*x_V_T_mAPC   # R43: V_T.mAPC -> null
    v44 <- k_mAPC_death*x_V_LN_mAPC   # R44: V_LN.mAPC -> null
    v45 <- k_c*(c0-x_V_T_c)*vol_V_T   # R45: null -> V_T.c
    v46 <- R_Tcell*DAMPs   # R46: null -> V_T.c
    v47 <- (kout*x_A_e_M1*m_A_e)-(kin*x_A_s_M1*m_A_s)   # R47: A_e.M1 -> A_s.M1
    v48 <- n_T0_clones*P0_C1*(a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1)))*x_V_T_C1*vol_V_T   # R48: null -> V_LN.P0
    v49 <- k_xP0_deg*x_V_LN_P0*m_V_LN   # R49: V_LN.P0 -> null
    v50 <- k_P0_up*x_V_LN_mAPC*x_V_LN_P0*m_V_LN   # R50: V_LN.P0 -> null
    v51 <- k_P0_up*1*x_V_LN_P0*m_V_e   # R51: null -> V_e.P0
    v52 <- k_P0_deg*x_V_e_P0*m_V_e   # R52: V_e.P0 -> V_e.p0
    v53 <- k_p0_deg*x_V_e_p0*m_V_e   # R53: V_e.p0 -> null
    v54 <- k_P0_on*x_V_e_p0*x_A_e_M1*m_A_e   # R54: V_e.p0 + A_e.M1 -> A_e.M1p0
    v55 <- k_P0_d1*k_P0_on*x_A_e_M1p0*m_A_e   # R55: A_e.M1p0 -> V_e.p0 + A_e.M1
    v56 <- k_P0_d1*k_P0_on*x_A_s_M1p0*m_A_s   # R56: A_s.M1p0 -> A_s.M1
    v57 <- kout*x_A_e_M1p0*m_A_e   # R57: A_e.M1p0 -> A_s.M1p0
    v58 <- n_T1_clones*P1_C1*(a_k_C1_death+a_k_C1_therapy+(((k_C_T1*x_V_T_T1)/(C_total+T_total+1))*(1-H_PD1_C1)))*x_V_T_C1*vol_V_T   # R58: null -> V_LN.P1
    v59 <- k_xP1_deg*x_V_LN_P1*m_V_LN   # R59: V_LN.P1 -> null
    v60 <- k_P1_up*x_V_LN_mAPC*x_V_LN_P1*m_V_LN   # R60: V_LN.P1 -> null
    v61 <- k_P1_up*1*x_V_LN_P1*m_V_e   # R61: null -> V_e.P1
    v62 <- k_P1_deg*x_V_e_P1*m_V_e   # R62: V_e.P1 -> V_e.p1
    v63 <- k_p1_deg*x_V_e_p1*m_V_e   # R63: V_e.p1 -> null
    v64 <- k_P1_on*x_V_e_p1*x_A_e_M1*m_A_e   # R64: V_e.p1 + A_e.M1 -> A_e.M1p1
    v65 <- k_P1_d1*k_P1_on*x_A_e_M1p1*m_A_e   # R65: A_e.M1p1 -> V_e.p1 + A_e.M1
    v66 <- k_P1_d1*k_P1_on*x_A_s_M1p1*m_A_s   # R66: A_s.M1p1 -> A_s.M1
    v67 <- kout*x_A_e_M1p1*m_A_e   # R67: A_e.M1p1 -> A_s.M1p1
    v68 <- q_P_atezolizumab*((x_V_C_atezolizumab/gamma_C_atezolizumab)-(x_V_P_atezolizumab/gamma_P_atezolizumab))   # R68: V_C.atezolizumab -> V_P.atezolizumab
    v69 <- q_T_atezolizumab*((x_V_C_atezolizumab/gamma_C_atezolizumab)-(x_V_T_atezolizumab/gamma_T_atezolizumab))*m_V_C   # R69: V_C.atezolizumab -> V_T.atezolizumab
    v70 <- q_LN_atezolizumab*((x_V_C_atezolizumab/gamma_C_atezolizumab)-(x_V_LN_atezolizumab/gamma_LN_atezolizumab))   # R70: V_C.atezolizumab -> V_LN.atezolizumab
    v71 <- ((q_LD_atezolizumab*x_V_T_atezolizumab)/gamma_T_atezolizumab)*vol_V_T   # R71: V_T.atezolizumab -> V_LN.atezolizumab
    v72 <- ((q_LD_atezolizumab*x_V_LN_atezolizumab)/gamma_LN_atezolizumab)*m_V_LN   # R72: V_LN.atezolizumab -> V_C.atezolizumab
    v73 <- k_cl_atezolizumab*x_V_C_atezolizumab*m_V_C   # R73: V_C.atezolizumab -> null
    v74 <- ((kon_PD1_PDL1*x_syn_T1_C1_PD1*x_syn_T1_C1_PDL1)-(koff_PD1_PDL1*x_syn_T1_C1_PD1_PDL1))*m_syn_T1_C1   # R74: syn_T1_C1.PD1 + syn_T1_C1.PDL1 -> syn_T1_C1.PD1_PDL1
    v75 <- ((kon_PD1_PDL2*x_syn_T1_C1_PD1*x_syn_T1_C1_PDL2)-(koff_PD1_PDL2*x_syn_T1_C1_PD1_PDL2))*m_syn_T1_C1   # R75: syn_T1_C1.PD1 + syn_T1_C1.PDL2 -> syn_T1_C1.PD1_PDL2
    v76 <- ((2*kon_PD1_aPD1*((x_syn_T1_C1_PD1*aPD1)/1))-(koff_PD1_aPD1*x_syn_T1_C1_PD1_aPD1))*m_syn_T1_C1   # R76: syn_T1_C1.PD1 -> syn_T1_C1.PD1_aPD1
    v77 <- ((Chi_PD1*kon_PD1_aPD1*x_syn_T1_C1_PD1*x_syn_T1_C1_PD1_aPD1)-(2*koff_PD1_aPD1*x_syn_T1_C1_PD1_aPD1_PD1))*m_syn_T1_C1   # R77: syn_T1_C1.PD1_aPD1 + syn_T1_C1.PD1 -> syn_T1_C1.PD1_aPD1_PD1
    v78 <- ((2*kon_PDL1_aPDL1*((x_syn_T1_C1_PDL1*aPDL1)/gamma_C_atezolizumab))-(koff_PDL1_aPDL1*x_syn_T1_C1_PDL1_aPDL1))*m_syn_T1_C1   # R78: syn_T1_C1.PDL1 -> syn_T1_C1.PDL1_aPDL1
    v79 <- ((Chi_PDL1*kon_PDL1_aPDL1*x_syn_T1_C1_PDL1*x_syn_T1_C1_PDL1_aPDL1)-(2*koff_PDL1_aPDL1*x_syn_T1_C1_PDL1_aPDL1_PDL1))*m_syn_T1_C1   # R79: syn_T1_C1.PDL1_aPDL1 + syn_T1_C1.PDL1 -> syn_T1_C1.PDL1_aPDL1_PDL1
    v80 <- ((kon_PD1_PDL1*x_syn_T1_APC_PD1*x_syn_T1_APC_PDL1)-(koff_PD1_PDL1*x_syn_T1_APC_PD1_PDL1))*m_syn_T1_APC   # R80: syn_T1_APC.PD1 + syn_T1_APC.PDL1 -> syn_T1_APC.PD1_PDL1
    v81 <- ((kon_PD1_PDL2*x_syn_T1_APC_PD1*x_syn_T1_APC_PDL2)-(koff_PD1_PDL2*x_syn_T1_APC_PD1_PDL2))*m_syn_T1_APC   # R81: syn_T1_APC.PD1 + syn_T1_APC.PDL2 -> syn_T1_APC.PD1_PDL2
    v82 <- ((2*kon_PD1_aPD1*((x_syn_T1_APC_PD1*aPD1)/gamma_LN_cibis))-(koff_PD1_aPD1*x_syn_T1_APC_PD1_aPD1))*m_syn_T1_APC   # R82: syn_T1_APC.PD1 -> syn_T1_APC.PD1_aPD1
    v83 <- ((Chi_PD1*kon_PD1_aPD1*x_syn_T1_APC_PD1*x_syn_T1_APC_PD1_aPD1)-(2*koff_PD1_aPD1*x_syn_T1_APC_PD1_aPD1_PD1))*m_syn_T1_APC   # R83: syn_T1_APC.PD1_aPD1 + syn_T1_APC.PD1 -> syn_T1_APC.PD1_aPD1_PD1
    v84 <- ((2*kon_PDL1_aPDL1*((x_syn_T1_APC_PDL1*aPDL1)/gamma_LN_cibis))-(koff_PDL1_aPDL1*x_syn_T1_APC_PDL1_aPDL1))*m_syn_T1_APC   # R84: syn_T1_APC.PDL1 -> syn_T1_APC.PDL1_aPDL1
    v85 <- ((Chi_PDL1*kon_PDL1_aPDL1*x_syn_T1_APC_PDL1*x_syn_T1_APC_PDL1_aPDL1)-(2*koff_PDL1_aPDL1*x_syn_T1_APC_PDL1_aPDL1_PDL1))*m_syn_T1_APC   # R85: syn_T1_APC.PDL1_aPDL1 + syn_T1_APC.PDL1 -> syn_T1_APC.PDL1_aPDL1_PDL1
    v86 <- q_P_cibis*((x_V_C_cibis/gamma_C_cibis)-(x_V_P_cibis/gamma_P_cibis))   # R86: V_C.cibis -> V_P.cibis
    v87 <- q_T_cibis*((x_V_C_cibis/gamma_C_cibis)-(x_V_T_cibis/gamma_T_cibis))*m_V_C   # R87: V_C.cibis -> V_T.cibis
    v88 <- q_LN_cibis*((x_V_C_cibis/gamma_C_cibis)-(x_V_LN_cibis/gamma_LN_cibis))   # R88: V_C.cibis -> V_LN.cibis
    v89 <- ((q_LD_cibis*x_V_T_cibis)/gamma_T_cibis)*vol_V_T   # R89: V_T.cibis -> V_LN.cibis
    v90 <- ((q_LD_cibis*x_V_LN_cibis)/gamma_LN_cibis)*m_V_LN   # R90: V_LN.cibis -> V_C.cibis
    v91 <- k_cl_cibis*x_V_C_cibis*m_V_C   # R91: V_C.cibis -> null
    v92 <- ((2*kon_CEA_cibis*((x_TCEsyn_T1_C1_CEA*x_V_T_cibis)/gamma_T_cibis))-(koff_CEA_cibis*x_TCEsyn_T1_C1_CEA_cibis))*m_TCEsyn_T1_C1   # R92: TCEsyn_T1_C1.CEA -> TCEsyn_T1_C1.CEA_cibis
    v93 <- ((Chi_CEA_cibis*kon_CEA_cibis*x_TCEsyn_T1_C1_CEA*x_TCEsyn_T1_C1_CEA_cibis)-(2*koff_CEA_cibis*x_TCEsyn_T1_C1_CEA_cibis_CEA))*m_TCEsyn_T1_C1   # R93: TCEsyn_T1_C1.CEA_cibis + TCEsyn_T1_C1.CEA -> TCEsyn_T1_C1.CEA_cibis_CEA
    v94 <- (((kon_CD3_cibis*x_TCEsyn_T1_C1_CD3*x_TCEsyn_T1_C1_CEA_cibis)/dT_syn)-(koff_CD3_cibis*x_TCEsyn_T1_C1_CEA_cibis_CD3))*m_TCEsyn_T1_C1   # R94: TCEsyn_T1_C1.CEA_cibis + TCEsyn_T1_C1.CD3 -> TCEsyn_T1_C1.CEA_cibis_CD3
    v95 <- (((kon_CD3_cibis*x_TCEsyn_T1_C1_CD3*x_TCEsyn_T1_C1_CEA_cibis_CEA)/dT_syn)-(koff_CD3_cibis*x_TCEsyn_T1_C1_CEA_cibis_CEA_CD3))*m_TCEsyn_T1_C1   # R95: TCEsyn_T1_C1.CEA_cibis_CEA + TCEsyn_T1_C1.CD3 -> TCEsyn_T1_C1.CEA_cibis_CEA_CD3
    v96 <- ((kon_CD3_cibis*((x_TCEsyn_T1_C1_CD3*x_V_T_cibis)/gamma_T_cibis))-(koff_CD3_cibis*x_TCEsyn_T1_C1_CD3_cibis))*m_TCEsyn_T1_C1   # R96: TCEsyn_T1_C1.CD3 -> TCEsyn_T1_C1.CD3_cibis
    v97 <- (((2*kon_CEA_cibis*x_TCEsyn_T1_C1_CD3_cibis*x_TCEsyn_T1_C1_CEA)/dT_syn)-(koff_CEA_cibis*x_TCEsyn_T1_C1_CEA_cibis_CD3))*m_TCEsyn_T1_C1   # R97: TCEsyn_T1_C1.CD3_cibis + TCEsyn_T1_C1.CEA -> TCEsyn_T1_C1.CEA_cibis_CD3
    v98 <- ((Chi_CEA_cibis*kon_CEA_cibis*x_TCEsyn_T1_C1_CEA*x_TCEsyn_T1_C1_CEA_cibis_CD3)-(2*koff_CEA_cibis*x_TCEsyn_T1_C1_CEA_cibis_CEA_CD3))*m_TCEsyn_T1_C1   # R98: TCEsyn_T1_C1.CEA_cibis_CD3 + TCEsyn_T1_C1.CEA -> TCEsyn_T1_C1.CEA_cibis_CEA_CD3
    v99 <- ((2*kon_CEA_cibis*((x_TCEsyn_T0_C1_CEA*x_V_T_cibis)/gamma_T_cibis))-(koff_CEA_cibis*x_TCEsyn_T0_C1_CEA_cibis))*m_TCEsyn_T0_C1   # R99: TCEsyn_T0_C1.CEA -> TCEsyn_T0_C1.CEA_cibis
    v100 <- ((Chi_CEA_cibis*kon_CEA_cibis*x_TCEsyn_T0_C1_CEA*x_TCEsyn_T0_C1_CEA_cibis)-(2*koff_CEA_cibis*x_TCEsyn_T0_C1_CEA_cibis_CEA))*m_TCEsyn_T0_C1   # R100: TCEsyn_T0_C1.CEA_cibis + TCEsyn_T0_C1.CEA -> TCEsyn_T0_C1.CEA_cibis_CEA
    v101 <- (((kon_CD3_cibis*x_TCEsyn_T0_C1_CD3*x_TCEsyn_T0_C1_CEA_cibis)/dT_syn)-(koff_CD3_cibis*x_TCEsyn_T0_C1_CEA_cibis_CD3))*m_TCEsyn_T0_C1   # R101: TCEsyn_T0_C1.CEA_cibis + TCEsyn_T0_C1.CD3 -> TCEsyn_T0_C1.CEA_cibis_CD3
    v102 <- (((kon_CD3_cibis*x_TCEsyn_T0_C1_CD3*x_TCEsyn_T0_C1_CEA_cibis_CEA)/dT_syn)-(koff_CD3_cibis*x_TCEsyn_T0_C1_CEA_cibis_CEA_CD3))*m_TCEsyn_T0_C1   # R102: TCEsyn_T0_C1.CEA_cibis_CEA + TCEsyn_T0_C1.CD3 -> TCEsyn_T0_C1.CEA_cibis_CEA_CD3
    v103 <- ((kon_CD3_cibis*((x_TCEsyn_T0_C1_CD3*x_V_T_cibis)/gamma_T_cibis))-(koff_CD3_cibis*x_TCEsyn_T0_C1_CD3_cibis))*m_TCEsyn_T0_C1   # R103: TCEsyn_T0_C1.CD3 -> TCEsyn_T0_C1.CD3_cibis
    v104 <- (((2*kon_CEA_cibis*x_TCEsyn_T0_C1_CD3_cibis*x_TCEsyn_T0_C1_CEA)/dT_syn)-(koff_CEA_cibis*x_TCEsyn_T0_C1_CEA_cibis_CD3))*m_TCEsyn_T0_C1   # R104: TCEsyn_T0_C1.CD3_cibis + TCEsyn_T0_C1.CEA -> TCEsyn_T0_C1.CEA_cibis_CD3
    v105 <- ((Chi_CEA_cibis*kon_CEA_cibis*x_TCEsyn_T0_C1_CEA*x_TCEsyn_T0_C1_CEA_cibis_CD3)-(2*koff_CEA_cibis*x_TCEsyn_T0_C1_CEA_cibis_CEA_CD3))*m_TCEsyn_T0_C1   # R105: TCEsyn_T0_C1.CEA_cibis_CD3 + TCEsyn_T0_C1.CEA -> TCEsyn_T0_C1.CEA_cibis_CEA_CD3

    # ---- ODEs ----
    d/dt(q_V_C_T1) <- -v10-v17+v18-v19+v20
    d/dt(q_V_C_T0) <- -v30-v34+v35-v36+v37
    d/dt(q_V_C_atezolizumab) <- (-v68-v69-v70+v72-v73)/nmol_to_molecule
    d/dt(q_V_C_cibis) <- (-v86-v87-v88+v90-v91)/nmol_to_molecule
    d/dt(q_V_P_T1) <- -v11+v17-v18
    d/dt(q_V_P_T0) <- -v31+v34-v35
    d/dt(q_V_P_atezolizumab) <- (v68)/nmol_to_molecule
    d/dt(q_V_P_cibis) <- (v86)/nmol_to_molecule
    d/dt(q_V_T_C_x) <- -v1+v4+v24
    d/dt(q_V_T_T_exh) <- -v2+v12+v14+v15+v16+v32
    d/dt(q_V_T_C1) <- v3-v4-v24
    d/dt(q_V_T_T1) <- -v12-v14-v15-v16+v19
    d/dt(q_V_T_T0) <- -v32+v36
    d/dt(q_V_T_APC) <- v39-v41
    d/dt(q_V_T_mAPC) <- v41-v42-v43
    d/dt(q_V_T_c) <- v45+v46
    d/dt(q_V_T_atezolizumab) <- (v69-v71)/nmol_to_molecule
    d/dt(q_V_T_cibis) <- (v87-v89)/nmol_to_molecule
    d/dt(q_V_LN_nT1) <- v5-v6-v7
    d/dt(q_V_LN_aT1) <- v7-v8
    d/dt(q_V_LN_T1) <- v9-v13-v20
    d/dt(q_V_LN_IL2) <- -v21-v22+v23-v38
    d/dt(q_V_LN_nT0) <- v25-v26-v27
    d/dt(q_V_LN_aT0) <- v27-v28
    d/dt(q_V_LN_T0) <- v29-v33-v37
    d/dt(q_V_LN_APC) <- v40
    d/dt(q_V_LN_mAPC) <- v42-v44
    d/dt(q_V_LN_P0) <- v48-v49-v50
    d/dt(q_V_LN_P1) <- v58-v59-v60
    d/dt(q_V_LN_atezolizumab) <- (v70+v71-v72)/nmol_to_molecule
    d/dt(q_V_LN_cibis) <- (v88+v89-v90)/nmol_to_molecule
    d/dt(q_V_e_P0) <- v51-v52
    d/dt(q_V_e_p0) <- v52-v53-v54+v55
    d/dt(q_V_e_P1) <- v61-v62
    d/dt(q_V_e_p1) <- v62-v63-v64+v65
    d/dt(q_A_e_M1) <- -v47-v54+v55-v64+v65
    d/dt(q_A_e_M1p0) <- v54-v55-v57
    d/dt(q_A_e_M1p1) <- v64-v65-v67
    d/dt(q_A_s_M1) <- v47+v56+v66
    d/dt(q_A_s_M1p0) <- -v56+v57
    d/dt(q_A_s_M1p1) <- -v66+v67
    d/dt(q_syn_T1_C1_PD1_PDL1) <- v74
    d/dt(q_syn_T1_C1_PD1_PDL2) <- v75
    d/dt(q_syn_T1_C1_PD1) <- -v74-v75-v76-v77
    d/dt(q_syn_T1_C1_PDL1) <- -v74-v78-v79
    d/dt(q_syn_T1_C1_PDL2) <- -v75
    d/dt(q_syn_T1_C1_PD1_aPD1) <- v76-v77
    d/dt(q_syn_T1_C1_PD1_aPD1_PD1) <- v77
    d/dt(q_syn_T1_C1_PDL1_aPDL1) <- v78-v79
    d/dt(q_syn_T1_C1_PDL1_aPDL1_PDL1) <- v79
    d/dt(q_syn_T1_APC_PD1_PDL1) <- v80
    d/dt(q_syn_T1_APC_PD1_PDL2) <- v81
    d/dt(q_syn_T1_APC_PD1) <- -v80-v81-v82-v83
    d/dt(q_syn_T1_APC_PDL1) <- -v80-v84-v85
    d/dt(q_syn_T1_APC_PDL2) <- -v81
    d/dt(q_syn_T1_APC_PD1_aPD1) <- v82-v83
    d/dt(q_syn_T1_APC_PD1_aPD1_PD1) <- v83
    d/dt(q_syn_T1_APC_PDL1_aPDL1) <- v84-v85
    d/dt(q_syn_T1_APC_PDL1_aPDL1_PDL1) <- v85
    d/dt(q_TCEsyn_T1_C1_CEA_cibis_CEA_CD3) <- v95+v98
    d/dt(q_TCEsyn_T1_C1_CEA_cibis_CEA) <- v93-v95
    d/dt(q_TCEsyn_T1_C1_CEA_cibis_CD3) <- v94+v97-v98
    d/dt(q_TCEsyn_T1_C1_CD3_cibis) <- v96-v97
    d/dt(q_TCEsyn_T1_C1_CEA_cibis) <- v92-v93-v94
    d/dt(q_TCEsyn_T1_C1_CEA) <- -v92-v93-v97-v98
    d/dt(q_TCEsyn_T1_C1_CD3) <- -v94-v95-v96
    d/dt(q_TCEsyn_T0_C1_CEA_cibis_CEA_CD3) <- v102+v105
    d/dt(q_TCEsyn_T0_C1_CEA_cibis_CEA) <- v100-v102
    d/dt(q_TCEsyn_T0_C1_CEA_cibis_CD3) <- v101+v104-v105
    d/dt(q_TCEsyn_T0_C1_CD3_cibis) <- v103-v104
    d/dt(q_TCEsyn_T0_C1_CEA_cibis) <- v99-v100-v101
    d/dt(q_TCEsyn_T0_C1_CEA) <- -v99-v100-v104-v105
    d/dt(q_TCEsyn_T0_C1_CD3) <- -v101-v102-v103

    # ---- initial conditions (SBML species initial values and initial assignments) ----
    q_V_T_C1(0) <- 1   # SBML V_T.C1 = 1 cell
    q_V_LN_IL2(0) <- 827.711370266   # SBML V_LN.IL2 = 1e-24 mole/millimeter^3
    q_V_LN_P0(0) <- 827.711370266   # SBML V_LN.P0 = 1e-24 mole/millimeter^3
    q_V_LN_P1(0) <- 827.711370266   # SBML V_LN.P1 = 1e-24 mole/millimeter^3
    q_V_e_P0(0) <- 2.408856796e-10   # SBML V_e.P0 = 1e-18 mole/litre
    q_V_e_p0(0) <- 2.408856796e-10   # SBML V_e.p0 = 1e-18 mole/litre
    q_V_e_P1(0) <- 2.408856796e-10   # SBML V_e.P1 = 1e-18 mole/litre
    q_V_e_p1(0) <- 2.408856796e-10   # SBML V_e.p1 = 1e-18 mole/litre
    q_A_e_M1(0) <- 32786.8852459   # SBML A_e.M1 = 2185.79235 molecule/micrometer^2
    q_A_e_M1p0(0) <- 1.5e-05   # SBML A_e.M1p0 = 1e-06 molecule/micrometer^2
    q_A_e_M1p1(0) <- 1.5e-05   # SBML A_e.M1p1 = 1e-06 molecule/micrometer^2
    q_A_s_M1(0) <- 0.0009   # SBML A_s.M1 = 1e-06 molecule/micrometer^2
    q_A_s_M1p0(0) <- 0.0009   # SBML A_s.M1p0 = 1e-06 molecule/micrometer^2
    q_A_s_M1p1(0) <- 0.0009   # SBML A_s.M1p1 = 1e-06 molecule/micrometer^2
    q_syn_T1_C1_PD1(0) <- (T_PD1_total/A_Tcell)*m_syn_T1_C1
    q_syn_T1_C1_PDL1(0) <- (C_PDL1_total/A_cell)*m_syn_T1_C1
    q_syn_T1_C1_PDL2(0) <- ((C_PDL1_total*r_PDL2C)/A_cell)*m_syn_T1_C1
    q_syn_T1_APC_PD1(0) <- (T_PD1_total/A_Tcell)*m_syn_T1_APC
    q_syn_T1_APC_PDL1(0) <- (APC_PDL1_total/A_APC)*m_syn_T1_APC
    q_syn_T1_APC_PDL2(0) <- ((APC_PDL1_total*r_PDL2APC)/A_APC)*m_syn_T1_APC
    q_TCEsyn_T1_C1_CEA(0) <- (C1_CEA_total/SA_cell)*m_TCEsyn_T1_C1
    q_TCEsyn_T1_C1_CD3(0) <- (T1_CD3_total/SA_Tcell)*m_TCEsyn_T1_C1
    q_TCEsyn_T0_C1_CEA(0) <- (C1_CEA_total/SA_cell)*m_TCEsyn_T0_C1
    q_TCEsyn_T0_C1_CD3(0) <- (T0_CD3_total/SA_Tcell)*m_TCEsyn_T0_C1

    # ---- outputs ----
    # tumour diameter [cm] from the total tumour volume assuming a sphere
    tumour_diameter <- (6*vol_V_T*1000/pi)^(1/3)
    # central (plasma) antibody concentrations [nmol/L]
    Cc_atezolizumab <- q_V_C_atezolizumab/vol_V_C
    Cc_cibisatamab <- q_V_C_cibis/vol_V_C
    # tumour-infiltrating T cell densities [cells/mL of tumour]
    teff_density <- x_V_T_T1/(vol_V_T*1000)
    treg_density <- x_V_T_T0/(vol_V_T*1000)
  })
}
