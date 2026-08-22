Demin_2025_zanubrutinib_acalabrutinib_ibrutinib_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of Bruton tyrosine kinase",
    "(BTK) occupancy by the three approved covalent BTK inhibitors",
    "zanubrutinib, acalabrutinib and ibrutinib (plus the active metabolites",
    "ACP-5862 and PCI-45227) in patients with B-cell malignancies. Couples a",
    "two-compartment population-PK layer per moiety (taken from the published",
    "popPK models cited in Table S1) to a perfusion-limited distribution step",
    "into lymph-node and bone-marrow interstitial fluid, an equilibrium",
    "partition into the intracellular space of CLL cells, and an",
    "ATP-competitive, covalent BTK-binding module with BTK turnover. BTK",
    "occupancy is resolved simultaneously in PBMCs, lymph nodes and bone",
    "marrow. 67 ODE states and 26 subject-level random effects; the three",
    "inhibitors compete for the same BTK pool in each tissue, so the model is",
    "a single coupled system rather than three independent models. No",
    "residual-error model: the source generated 1000 virtual patients by",
    "sampling the parameter distributions in Table S1 and compared predicted",
    "medians and 95% CI bands against observed occupancy.",
    sep = " "
  )
  reference <- paste(
    "Demin O Jr, Ou Y, Kolesova G, Shchelokov D, Stepanov A, Musatova V,",
    "Sahasranaman S, Zhao Y, Liu X, Tang Z, Hanley WD. Quantitative Systems",
    "Pharmacology Model to Predict Target Occupancy by Bruton Tyrosine Kinase",
    "Inhibitors in Patients With B-Cell Malignancies. CPT Pharmacometrics Syst",
    "Pharmacol. 2025;14(4):706-717. doi:10.1002/psp4.13307.",
    sep = " "
  )
  vignette <- "Demin_2025_zanubrutinib_acalabrutinib_ibrutinib_qsp"

  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  dosing <- c("gut_zan", "gut_aca", "gut_ibr")

  paper_specific_compartment_pattern <- "^(btk|isf)_"

  compartmentData <- list(
    gut_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    gut_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    gut_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit5_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    depot_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    depot_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    depot_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    central_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_acp = list(analyte = "ACP-5862", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_acp = list(analyte = "ACP-5862", units = "mg", specimen = "plasma", verified = TRUE),
    central_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_pci = list(analyte = "PCI-45227", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_pci = list(analyte = "PCI-45227", units = "mg", specimen = "plasma", verified = TRUE),
    isf_ln_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "lymph", verified = TRUE),
    isf_ln_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "lymph", verified = TRUE),
    isf_ln_acp = list(analyte = "ACP-5862", units = "mg", specimen = "lymph", verified = TRUE),
    isf_ln_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "lymph", verified = TRUE),
    isf_ln_pci = list(analyte = "PCI-45227", units = "mg", specimen = "lymph", verified = TRUE),
    isf_bm_zan = list(analyte = "zanubrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    isf_bm_aca = list(analyte = "acalabrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    isf_bm_acp = list(analyte = "ACP-5862", units = "mg", specimen = "tissue", verified = TRUE),
    isf_bm_ibr = list(analyte = "ibrutinib", units = "mg", specimen = "tissue", verified = TRUE),
    isf_bm_pci = list(analyte = "PCI-45227", units = "mg", specimen = "tissue", verified = TRUE),
    btk_free_pbmc = list(analyte = "Bruton tyrosine kinase (unoccupied)", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_atp_pbmc = list(analyte = "BTK-ATP complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_zan_pbmc = list(analyte = "BTK-zanubrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_zan_inact_pbmc = list(analyte = "irreversibly inactivated BTK-zanubrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_aca_pbmc = list(analyte = "BTK-acalabrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_aca_inact_pbmc = list(analyte = "irreversibly inactivated BTK-acalabrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_acp_pbmc = list(analyte = "BTK-ACP-5862 complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_acp_inact_pbmc = list(analyte = "irreversibly inactivated BTK-ACP-5862 complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_ibr_pbmc = list(analyte = "BTK-ibrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_ibr_inact_pbmc = list(analyte = "irreversibly inactivated BTK-ibrutinib complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_pci_pbmc = list(analyte = "BTK-PCI-45227 complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_pci_inact_pbmc = list(analyte = "irreversibly inactivated BTK-PCI-45227 complex", units = "nmol", specimen = "blood cell", verified = TRUE),
    btk_free_ln = list(analyte = "Bruton tyrosine kinase (unoccupied)", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_atp_ln = list(analyte = "BTK-ATP complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_zan_ln = list(analyte = "BTK-zanubrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_zan_inact_ln = list(analyte = "irreversibly inactivated BTK-zanubrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_aca_ln = list(analyte = "BTK-acalabrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_aca_inact_ln = list(analyte = "irreversibly inactivated BTK-acalabrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_acp_ln = list(analyte = "BTK-ACP-5862 complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_acp_inact_ln = list(analyte = "irreversibly inactivated BTK-ACP-5862 complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_ibr_ln = list(analyte = "BTK-ibrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_ibr_inact_ln = list(analyte = "irreversibly inactivated BTK-ibrutinib complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_pci_ln = list(analyte = "BTK-PCI-45227 complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_pci_inact_ln = list(analyte = "irreversibly inactivated BTK-PCI-45227 complex", units = "nmol", specimen = "lymph", verified = TRUE),
    btk_free_bm = list(analyte = "Bruton tyrosine kinase (unoccupied)", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_atp_bm = list(analyte = "BTK-ATP complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_zan_bm = list(analyte = "BTK-zanubrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_zan_inact_bm = list(analyte = "irreversibly inactivated BTK-zanubrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_aca_bm = list(analyte = "BTK-acalabrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_aca_inact_bm = list(analyte = "irreversibly inactivated BTK-acalabrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_acp_bm = list(analyte = "BTK-ACP-5862 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_acp_inact_bm = list(analyte = "irreversibly inactivated BTK-ACP-5862 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_ibr_bm = list(analyte = "BTK-ibrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_ibr_inact_bm = list(analyte = "irreversibly inactivated BTK-ibrutinib complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_pci_bm = list(analyte = "BTK-PCI-45227 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    btk_pci_inact_bm = list(analyte = "irreversibly inactivated BTK-PCI-45227 complex", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 1000L,
    disease_state = paste(
      "B-cell malignancies (chronic lymphocytic leukaemia / small lymphocytic",
      "lymphoma, mantle cell lymphoma and Waldenstrom macroglobulinaemia);",
      "the BTK-binding module is parameterised for primary CLL cells.",
      sep = " "
    ),
    dose_range    = paste(
      "Zanubrutinib 160 mg BID or 320 mg QD (reduced doses 80 mg BID, 80 mg",
      "QD); acalabrutinib 100 mg BID, 100 mg QD or 200 mg QD; ibrutinib",
      "420 mg QD (reduced doses 280 mg QD, 140 mg QD).",
      sep = " "
    ),
    notes         = paste(
      "In silico cohort: 1000 virtual patients generated by sampling the",
      "parameter distributions in Table S1 (Methods section 2.5). The",
      "underlying PK parameters come from published population-PK analyses",
      "reproduced in Table S1: zanubrutinib Ou 2021 (PMID 33306268),",
      "acalabrutinib and ACP-5862 Edlund 2022 (PMID 34265100) and Edlund 2019",
      "(PMID 30556110), ibrutinib and PCI-45227 Gallais 2020 (PMID 32328976).",
      "Observed BTK-occupancy validation data are from patients with CLL.",
      sep = " "
    )
  )

  ini({
    # ---- Systemic / physiological parameters (Table S1 "Systemic parameters") ----
    q_ln <- fixed(151.776); label("Blood plasma flow to lymph nodes (L/day)")  # Table S1 Q_ln = 152 L/day; model code Const Q_ln = 151.776 (cardiac output x LN fraction, ICRP 89)
    q_bm <- fixed(267.84); label("Blood plasma flow to bone marrow (L/day)")  # Table S1 Q_bm = 268 L/day; model code Const Q_bm = 267.84 (cardiac output x BM fraction, ICRP 89)
    v_isf_ln <- fixed(0.051); label("Interstitial-fluid volume of lymph nodes (L)")  # model code Compartment LN_isf = 0.051 L (0.186 L total LN volume x 0.273 ISF fraction)
    v_isf_bm <- fixed(0.2); label("Interstitial-fluid volume of bone marrow (L)")  # model code Compartment BM_isf = 0.2 L (1.0 L total BM volume x 0.2 ISF fraction)
    bl_btk <- fixed(561.3); label("Baseline total intracellular BTK concentration in CLL cells (nmol/L)")  # Table S1 BTK_base = 561 nM; model code Const BTK_cell_base = 561.3
    koff_atp <- fixed(12096); label("Dissociation rate constant of ATP from BTK (1/day)")  # Table S1 k_off^atp = 1.21e4 /day; model code Const k_off_atp_btk = 12096
    kd_atp <- fixed(26000); label("Equilibrium dissociation constant of ATP for BTK (nmol/L)")  # Table S1 K_d^atp = 2.6e4 nM
    c_adj_pbmc <- fixed(7.219); label("Adjustment coefficient for intracellular BTKi accumulation (unitless)")  # Table S1 c_adj^pbmc = 7.219, fitted to intracellular ibrutinib in WBC (Figure S2)
    hct_ref <- fixed(0.43); label("Reference haematocrit (unitless)")  # model code Const Hct_ref = 0.43 (ICRP 89 text page 138)
    coef_kdeg_ln <- fixed(0.772414); label("Ratio of BTK degradation rate in lymph node to blood (unitless)")  # Table S1 coef_ln^deg_btk = 0.772 (11.2%/14.5% BTK turnover, LN vs PBMC)
    coef_kdeg_bm <- fixed(0.427586); label("Ratio of BTK degradation rate in bone marrow to blood (unitless)")  # Table S1 coef_bm^deg_btk = 0.428 (6.2%/14.5% BTK turnover, BM vs PBMC)

    # ---- Systemic parameters carrying between-subject variability ----
    latp <- fixed(log(2.271e+06)); label("Intracellular ATP concentration in CLL cells (nmol/L)")  # Table S1 ATP_base = 2.27e6 nM, IIV 40% CV
    lkdeg_btk <- fixed(log(1.03972)); label("BTK degradation rate constant in blood/PBMC (1/day)")  # Table S1 k_deg^btk = 1.04 /day (half-life ~16 h), IIV 30% CV

    # ---- Drug-specific binding and physicochemical parameters (Table S1) ----
    # zanubrutinib
    koff_zan <- fixed(8.1216); label("Dissociation rate constant of zanubrutinib from BTK (1/day)")  # Table S1 k_off^zan = 8.12, K_d^zan = 0.264 (estimated as Kd_ibr x IC50_zan/IC50_ibr), k_inact^zan = 2880, BP_ratio^zan = 0.804
    kd_zan <- fixed(0.264); label("Equilibrium dissociation constant of zanubrutinib for BTK (nmol/L)")  # Table S1 k_off^zan = 8.12, K_d^zan = 0.264 (estimated as Kd_ibr x IC50_zan/IC50_ibr), k_inact^zan = 2880, BP_ratio^zan = 0.804
    kinact_zan <- fixed(2877.12); label("Covalent BTK inactivation rate constant by zanubrutinib (1/day)")  # Table S1 k_off^zan = 8.12, K_d^zan = 0.264 (estimated as Kd_ibr x IC50_zan/IC50_ibr), k_inact^zan = 2880, BP_ratio^zan = 0.804
    bp_zan <- fixed(0.804); label("Blood-to-plasma ratio of zanubrutinib (unitless)")  # Table S1 k_off^zan = 8.12, K_d^zan = 0.264 (estimated as Kd_ibr x IC50_zan/IC50_ibr), k_inact^zan = 2880, BP_ratio^zan = 0.804
    kp_isf_zan <- fixed(1); label("ISF-to-unbound-plasma partition coefficient of zanubrutinib (unitless)")  # Table S1 Kp_isf_pu = 1 for all moieties (Rodgers 2006: Cu,ew = Cu,p for weak bases, neutrals and zwitterions)
    fup_zan <- fixed(0.055); label("Fraction unbound in plasma of zanubrutinib (unitless)")  # Table S1 fup^zan = 0.055 (FDA BRUKINSA multidisciplinary review)
    mr_zan <- fixed(471.55); label("Molecular weight of zanubrutinib (g/mol)")  # model code Const Mr_zan = 471.55 g/mol
    # acalabrutinib
    koff_aca <- fixed(7.92); label("Dissociation rate constant of acalabrutinib from BTK (1/day)")  # Table S1 k_off^aca = 7.92, K_d^aca = 0.816 (estimated as Kd_ibr x IC50_aca/IC50_ibr), k_inact^aca = 483, BP_ratio^aca = 0.79
    kd_aca <- fixed(0.816); label("Equilibrium dissociation constant of acalabrutinib for BTK (nmol/L)")  # Table S1 k_off^aca = 7.92, K_d^aca = 0.816 (estimated as Kd_ibr x IC50_aca/IC50_ibr), k_inact^aca = 483, BP_ratio^aca = 0.79
    kinact_aca <- fixed(482.976); label("Covalent BTK inactivation rate constant by acalabrutinib (1/day)")  # Table S1 k_off^aca = 7.92, K_d^aca = 0.816 (estimated as Kd_ibr x IC50_aca/IC50_ibr), k_inact^aca = 483, BP_ratio^aca = 0.79
    bp_aca <- fixed(0.79); label("Blood-to-plasma ratio of acalabrutinib (unitless)")  # Table S1 k_off^aca = 7.92, K_d^aca = 0.816 (estimated as Kd_ibr x IC50_aca/IC50_ibr), k_inact^aca = 483, BP_ratio^aca = 0.79
    kp_isf_aca <- fixed(1); label("ISF-to-unbound-plasma partition coefficient of acalabrutinib (unitless)")  # Table S1 Kp_isf_pu = 1 for all moieties (Rodgers 2006: Cu,ew = Cu,p for weak bases, neutrals and zwitterions)
    fup_aca <- fixed(0.025); label("Fraction unbound in plasma of acalabrutinib (unitless)")  # Table S1 fup^aca = 0.025
    mr_aca <- fixed(465.52); label("Molecular weight of acalabrutinib (g/mol)")  # model code Const Mr_aca = 465.52 g/mol
    # ACP-5862
    koff_acp <- fixed(7.92); label("Dissociation rate constant of ACP-5862 from BTK (1/day)")  # Table S1 k_off^acp = 7.92 (assumed equal to acalabrutinib), K_d^acp = 0.8, k_inact^acp = 268, BP_ratio^acp = 0.66
    kd_acp <- fixed(0.8); label("Equilibrium dissociation constant of ACP-5862 for BTK (nmol/L)")  # Table S1 k_off^acp = 7.92 (assumed equal to acalabrutinib), K_d^acp = 0.8, k_inact^acp = 268, BP_ratio^acp = 0.66
    kinact_acp <- fixed(267.84); label("Covalent BTK inactivation rate constant by ACP-5862 (1/day)")  # Table S1 k_off^acp = 7.92 (assumed equal to acalabrutinib), K_d^acp = 0.8, k_inact^acp = 268, BP_ratio^acp = 0.66
    bp_acp <- fixed(0.66); label("Blood-to-plasma ratio of ACP-5862 (unitless)")  # Table S1 k_off^acp = 7.92 (assumed equal to acalabrutinib), K_d^acp = 0.8, k_inact^acp = 268, BP_ratio^acp = 0.66
    kp_isf_acp <- fixed(1); label("ISF-to-unbound-plasma partition coefficient of ACP-5862 (unitless)")  # Table S1 Kp_isf_pu = 1 for all moieties (Rodgers 2006: Cu,ew = Cu,p for weak bases, neutrals and zwitterions)
    fup_acp <- fixed(0.014); label("Fraction unbound in plasma of ACP-5862 (unitless)")  # Table S1 fup^acp = 0.014
    mr_acp <- fixed(481.52); label("Molecular weight of ACP-5862 (g/mol)")  # model code Const Mr_acp = 481.52 g/mol
    # ibrutinib
    koff_ibr <- fixed(11.1744); label("Dissociation rate constant of ibrutinib from BTK (1/day)")  # Table S1 k_off^ibr = 11.2, K_d^ibr = 0.24 (taken directly from the literature), k_inact^ibr = 2300, BP_ratio^ibr = 0.7
    kd_ibr <- fixed(0.24); label("Equilibrium dissociation constant of ibrutinib for BTK (nmol/L)")  # Table S1 k_off^ibr = 11.2, K_d^ibr = 0.24 (taken directly from the literature), k_inact^ibr = 2300, BP_ratio^ibr = 0.7
    kinact_ibr <- fixed(2298.24); label("Covalent BTK inactivation rate constant by ibrutinib (1/day)")  # Table S1 k_off^ibr = 11.2, K_d^ibr = 0.24 (taken directly from the literature), k_inact^ibr = 2300, BP_ratio^ibr = 0.7
    bp_ibr <- fixed(0.7); label("Blood-to-plasma ratio of ibrutinib (unitless)")  # Table S1 k_off^ibr = 11.2, K_d^ibr = 0.24 (taken directly from the literature), k_inact^ibr = 2300, BP_ratio^ibr = 0.7
    kp_isf_ibr <- fixed(1); label("ISF-to-unbound-plasma partition coefficient of ibrutinib (unitless)")  # Table S1 Kp_isf_pu = 1 for all moieties (Rodgers 2006: Cu,ew = Cu,p for weak bases, neutrals and zwitterions)
    fup_ibr <- fixed(0.027); label("Fraction unbound in plasma of ibrutinib (unitless)")  # Table S1 fup^ibr = 0.027
    mr_ibr <- fixed(440.51); label("Molecular weight of ibrutinib (g/mol)")  # model code Const Mr_ibr = 440.51 g/mol
    # PCI-45227
    koff_pci <- fixed(11.1744); label("Dissociation rate constant of PCI-45227 from BTK (1/day)")  # Table S1 k_off^pci = 11.2 and k_inact^pci = 2300 (both assumed equal to ibrutinib), K_d^pci = 3.6 (15x less potent than ibrutinib), BP_ratio^pci = 1.1
    kd_pci <- fixed(3.6); label("Equilibrium dissociation constant of PCI-45227 for BTK (nmol/L)")  # Table S1 k_off^pci = 11.2 and k_inact^pci = 2300 (both assumed equal to ibrutinib), K_d^pci = 3.6 (15x less potent than ibrutinib), BP_ratio^pci = 1.1
    kinact_pci <- fixed(2298.24); label("Covalent BTK inactivation rate constant by PCI-45227 (1/day)")  # Table S1 k_off^pci = 11.2 and k_inact^pci = 2300 (both assumed equal to ibrutinib), K_d^pci = 3.6 (15x less potent than ibrutinib), BP_ratio^pci = 1.1
    bp_pci <- fixed(1.1); label("Blood-to-plasma ratio of PCI-45227 (unitless)")  # Table S1 k_off^pci = 11.2 and k_inact^pci = 2300 (both assumed equal to ibrutinib), K_d^pci = 3.6 (15x less potent than ibrutinib), BP_ratio^pci = 1.1
    kp_isf_pci <- fixed(1); label("ISF-to-unbound-plasma partition coefficient of PCI-45227 (unitless)")  # Table S1 Kp_isf_pu = 1 for all moieties (Rodgers 2006: Cu,ew = Cu,p for weak bases, neutrals and zwitterions)
    fup_pci <- fixed(0.09); label("Fraction unbound in plasma of PCI-45227 (unitless)")  # Table S1 fup^pci = 0.09 (FDA IMBRUVICA clinical pharmacology review)
    mr_pci <- fixed(474.52); label("Molecular weight of PCI-45227 (g/mol)")  # model code Const Mr_pci = 474.52 g/mol

    # ---- Zanubrutinib PK (Table S1; popPK model of Ou 2021, PMID 33306268) ----
    kabs_zan <- fixed(11.3857); label("First-order absorption rate constant of zanubrutinib from the depot (1/day)")  # Table S1 k_abs^zan = 11.4 /day
    f1_zan <- fixed(1); label("Bioavailability factor of zanubrutinib (unitless)")  # Table S1 F1^zan = 1, fixed because CL, Vc, Vp and Q are apparent (divided by F)
    ld1_zan <- fixed(log(0.0513416)); label("Duration of zero-order transfer into the zanubrutinib depot (day)")  # model code Const D1_zan_mean = 0.0513416 day (1.23 h), IIV 51.9% CV; Table S1 prints 0.513 which is a factor-of-10 typo (see vignette Errata)
    lvc_zan <- fixed(log(73.9952)); label("Apparent central volume of zanubrutinib, Vc/F (L)")  # Table S1 Vc^zan = 74 L, IIV 65% CV, IOV 59.2% CV
    lvp_zan <- fixed(log(489.312)); label("Apparent peripheral volume of zanubrutinib, Vp/F (L)")  # Table S1 Vp^zan = 489 L, IIV 72.8% CV
    lcl_zan <- fixed(log(2143.2)); label("Apparent clearance of zanubrutinib, CL/F (L/day)")  # Table S1 CL^zan = 2140 L/day, IIV 37.6% CV, IOV 34.2% CV
    lq_zan <- fixed(log(318.319)); label("Apparent inter-compartmental clearance of zanubrutinib, Q/F (L/day)")  # Table S1 Q^zan = 318 L/day, IIV 118.7% CV

    # ---- Acalabrutinib and ACP-5862 PK (Table S1; Edlund 2022 PMID 34265100, Edlund 2019 PMID 30556110) ----
    kabs_aca <- fixed(60); label("First-order absorption rate constant of acalabrutinib from the depot (1/day)")  # Table S1 k_abs^aca = 60 /day
    ntr_aca <- fixed(5); label("Number of transit compartments for acalabrutinib absorption (unitless)")  # Table S1 N_tr^aca = 5
    lmtt_aca <- fixed(log(0.0104167)); label("Mean transit time of acalabrutinib absorption (day)")  # Table S1 MTT^aca = 0.0104 day (15 min), IOV 118% CV
    lf1_aca <- fixed(log(1)); label("Bioavailability factor of acalabrutinib (unitless)")  # Table S1 F1^aca = 1, IOV 56.1% CV
    lvc_aca <- fixed(log(31)); label("Apparent central volume of acalabrutinib, Vc/F (L)")  # Table S1 Vc^aca = 31 L, IIV 270% CV
    lvp_aca <- fixed(log(213)); label("Apparent peripheral volume of acalabrutinib, Vp/F (L)")  # Table S1 Vp^aca = 213 L; the model code applies the Vc eta to this parameter (see vignette Errata)
    lcl_aca <- fixed(log(2450)); label("Apparent clearance of acalabrutinib, CL/F (L/day)")  # Table S1 CL^aca = 2450 L/day, IIV 23.8% CV
    lq_aca <- fixed(log(350.4)); label("Apparent inter-compartmental clearance of acalabrutinib, Q/F (L/day)")  # Table S1 Q^aca = 350 L/day, no variability reported
    fm_aca <- fixed(0.4); label("Fraction of acalabrutinib clearance forming ACP-5862 (unitless)")  # Table S1 Fm^aca = 0.4
    lvc_acp <- fixed(log(38.5)); label("Apparent central volume of ACP-5862, Vc/F (L)")  # Table S1 Vc^acp = 38.5 L, IIV 77.1% CV
    lvp_acp <- fixed(log(38.4)); label("Apparent peripheral volume of ACP-5862, Vp/F (L)")  # Table S1 Vp^acp = 38.4 L, IIV 19.2% CV
    lcl_acp <- fixed(log(523.2)); label("Apparent clearance of ACP-5862, CL/F (L/day)")  # Table S1 CL^acp = 523 L/day, IIV 31.4% CV
    lq_acp <- fixed(log(146.4)); label("Apparent inter-compartmental clearance of ACP-5862, Q/F (L/day)")  # Table S1 Q^acp = 146 L/day, IIV 40.7% CV

    # ---- Ibrutinib and PCI-45227 PK (Table S1; Gallais 2020 PMID 32328976) ----
    kabs_ibr <- fixed(37.44); label("First-order absorption rate constant of ibrutinib from the depot (1/day)")  # Table S1 k_abs^ibr = 37.4 /day
    alag_ibr <- fixed(0.00991667); label("Absorption lag time of ibrutinib (day)")  # Table S1 ALAG^ibr = 0.00992 day (14.3 min); the model code holds it constant (see vignette Errata)
    f1_ibr <- fixed(1); label("Bioavailability factor of ibrutinib (unitless)")  # Table S1 F1^ibr = 1, fixed because CL, Vc, Vp and Q are apparent
    ld1_ibr <- fixed(log(0.0412083)); label("Duration of zero-order transfer into the ibrutinib depot (day)")  # Table S1 D1^ibr = 0.0412 day (0.99 h), IIV 115.2% CV
    lvc_ibr <- fixed(log(1010)); label("Apparent central volume of ibrutinib and PCI-45227, Vc/F (L)")  # Table S1 Vc^ibr = 1010 L, IIV 81.8% CV; PCI-45227 shares this volume
    lvp_ibr <- fixed(log(1480)); label("Apparent peripheral volume of ibrutinib and PCI-45227, Vp/F (L)")  # Table S1 Vp^ibr = 1480 L, IIV 76.9% CV; PCI-45227 shares this volume
    lcl_ibr <- fixed(log(5808)); label("Apparent clearance of ibrutinib, CL/F (L/day)")  # Table S1 CL^ibr = 5808 L/day, IIV 66.5% CV, IOV 46.7% CV
    lq_ibr <- fixed(log(4104)); label("Apparent inter-compartmental clearance of ibrutinib, Q/F (L/day)")  # Table S1 Q^ibr = 4104 L/day, no variability reported
    lkabs_pci <- fixed(log(29.04)); label("First-order absorption rate constant of PCI-45227 from the ibrutinib depot (1/day)")  # Table S1 k_abs^pci = 29 /day, IIV 64.2% CV
    lclmet_ibr <- fixed(log(3600)); label("Apparent formation clearance of ibrutinib to PCI-45227 (L/day)")  # Table S1 CL_met^ibr = 3600 L/day, IIV 64.4% CV
    lcl_pci <- fixed(log(4344)); label("Apparent clearance of PCI-45227, CL/F (L/day)")  # Table S1 CL^pci = 4344 L/day, IIV 50.7% CV, IOV 25.7% CV
    lq_pci <- fixed(log(1200)); label("Apparent inter-compartmental clearance of PCI-45227, Q/F (L/day)")  # Table S1 Q^pci = 1200 L/day, no variability reported

    # ---- Between-subject variability (Table S1 "IIV / IOV, CV%" column) ----
    etald1_zan ~ 0.238514  # Table S1 D1^zan IIV 51.9% CV; log(1 + 0.519^2)
    etalvc_zan ~ 0.352416  # Table S1 Vc^zan IIV 65% CV; log(1 + 0.65^2)
    etalvp_zan ~ 0.425257  # Table S1 Vp^zan IIV 72.8% CV; log(1 + 0.728^2)
    etalcl_zan ~ 0.132235  # Table S1 CL^zan IIV 37.6% CV; log(1 + 0.376^2)
    etalq_zan ~ 0.879199  # Table S1 Q^zan IIV 118.7% CV; log(1 + 1.187^2)
    etalvc_aca ~ 2.11505  # Table S1 Vc^aca IIV 270% CV; log(1 + 2.7^2)
    etalcl_aca ~ 0.0550978  # Table S1 CL^aca IIV 23.8% CV; log(1 + 0.238^2)
    etalvc_acp ~ 0.466523  # Table S1 Vc^acp IIV 77.1% CV; log(1 + 0.771^2)
    etalvp_acp ~ 0.0362008  # Table S1 Vp^acp IIV 19.2% CV; log(1 + 0.192^2)
    etalcl_acp ~ 0.094033  # Table S1 CL^acp IIV 31.4% CV; log(1 + 0.314^2)
    etalq_acp ~ 0.153278  # Table S1 Q^acp IIV 40.7% CV; log(1 + 0.407^2)
    etald1_ibr ~ 0.844625  # Table S1 D1^ibr IIV 115.2% CV; log(1 + 1.152^2)
    etalvc_ibr ~ 0.512299  # Table S1 Vc^ibr IIV 81.8% CV; log(1 + 0.818^2)
    etalvp_ibr ~ 0.46459  # Table S1 Vp^ibr IIV 76.9% CV; log(1 + 0.769^2)
    etalcl_ibr ~ 0.366187  # Table S1 CL^ibr IIV 66.5% CV; log(1 + 0.665^2)
    etalkabs_pci ~ 0.345123  # Table S1 k_abs^pci IIV 64.2% CV; log(1 + 0.642^2)
    etalclmet_ibr ~ 0.346943  # Table S1 CL_met^ibr IIV 64.4% CV; log(1 + 0.644^2)
    etalcl_pci ~ 0.228767  # Table S1 CL^pci IIV 50.7% CV; log(1 + 0.507^2)
    etalatp ~ 0.14842  # Table S1 ATP_base IIV 40% CV; log(1 + 0.4^2)
    etalkdeg_btk ~ 0.0861777  # Table S1 k_deg^btk IIV 30% CV; log(1 + 0.3^2)

    # ---- Inter-occasion variability, encoded as a second subject-level random
    #      effect because the source model code carries each IOV deviate as a
    #      per-virtual-patient constant (see vignette Errata) ----
    etaiovlvc_zan ~ 0.300448  # Table S1 Vc^zan IOV 59.2% CV; log(1 + 0.592^2)
    etaiovlcl_zan ~ 0.110614  # Table S1 CL^zan IOV 34.2% CV; log(1 + 0.342^2)
    etaiovlmtt_aca ~ 0.872297  # Table S1 MTT^aca IOV 118% CV; log(1 + 1.18^2)
    etaiovlf1_aca ~ 0.273624  # Table S1 F1^aca IOV 56.1% CV; log(1 + 0.561^2)
    etaiovlcl_ibr ~ 0.197283  # Table S1 CL^ibr IOV 46.7% CV; log(1 + 0.467^2)
    etaiovlcl_pci ~ 0.0639593  # Table S1 CL^pci IOV 25.7% CV; log(1 + 0.257^2)
  })

  model({
    # ---------------- Individual parameter values ----------------
    # The IOV deviates are routed through named intermediates because rxode2
    # mu-referencing does not accept theta + eta_iiv + eta_iov on one line.
    iov_vc_zan <- etaiovlvc_zan
    iov_cl_zan <- etaiovlcl_zan
    iov_cl_ibr <- etaiovlcl_ibr
    iov_cl_pci <- etaiovlcl_pci

    d1_zan  <- exp(ld1_zan + etald1_zan)
    vc_zan  <- exp(lvc_zan + etalvc_zan + iov_vc_zan)
    vp_zan  <- exp(lvp_zan + etalvp_zan)
    cl_zan  <- exp(lcl_zan + etalcl_zan + iov_cl_zan)
    q_zan   <- exp(lq_zan + etalq_zan)

    mtt_aca <- exp(lmtt_aca + etaiovlmtt_aca)
    ktr_aca <- (ntr_aca + 1) / mtt_aca
    f1_aca  <- exp(lf1_aca + etaiovlf1_aca)
    vc_aca  <- exp(lvc_aca + etalvc_aca)
    # The source model code applies the Vc eta (not a Vp eta) to Vp^aca.
    vp_aca  <- exp(lvp_aca + etalvc_aca)
    cl_aca  <- exp(lcl_aca + etalcl_aca)
    q_aca   <- exp(lq_aca)
    vc_acp  <- exp(lvc_acp + etalvc_acp)
    vp_acp  <- exp(lvp_acp + etalvp_acp)
    cl_acp  <- exp(lcl_acp + etalcl_acp)
    q_acp   <- exp(lq_acp + etalq_acp)

    d1_ibr    <- exp(ld1_ibr + etald1_ibr)
    vc_ibr    <- exp(lvc_ibr + etalvc_ibr)
    vp_ibr    <- exp(lvp_ibr + etalvp_ibr)
    cl_ibr    <- exp(lcl_ibr + etalcl_ibr + iov_cl_ibr)
    q_ibr     <- exp(lq_ibr)
    kabs_pci  <- exp(lkabs_pci + etalkabs_pci)
    clmet_ibr <- exp(lclmet_ibr + etalclmet_ibr)
    cl_pci    <- exp(lcl_pci + etalcl_pci + iov_cl_pci)
    q_pci     <- exp(lq_pci)

    atp      <- exp(latp + etalatp)
    kdeg_btk <- exp(lkdeg_btk + etalkdeg_btk)
    kdeg_pbmc <- kdeg_btk
    kdeg_ln   <- kdeg_btk * coef_kdeg_ln
    kdeg_bm   <- kdeg_btk * coef_kdeg_bm

    # ---------------- Partition coefficients (equations 8 and 9) ----------------
    kp_rbc_zan  <- (bp_zan - (1 - hct_ref)) / hct_ref
    kp_rbc_aca  <- (bp_aca - (1 - hct_ref)) / hct_ref
    kp_rbc_acp  <- (bp_acp - (1 - hct_ref)) / hct_ref
    kp_rbc_ibr  <- (bp_ibr - (1 - hct_ref)) / hct_ref
    kp_rbc_pci  <- (bp_pci - (1 - hct_ref)) / hct_ref
    kp_cell_zan <- (kp_rbc_zan / fup_zan) * c_adj_pbmc
    kp_cell_aca <- (kp_rbc_aca / fup_aca) * c_adj_pbmc
    kp_cell_acp <- (kp_rbc_acp / fup_acp) * c_adj_pbmc
    kp_cell_ibr <- (kp_rbc_ibr / fup_ibr) * c_adj_pbmc
    kp_cell_pci <- (kp_rbc_pci / fup_pci) * c_adj_pbmc

    # ---------------- Plasma concentrations (mg/L) ----------------
    cpls_zan <- central_zan / vc_zan
    cprf_zan <- peripheral1_zan / vp_zan
    cpls_aca <- central_aca / vc_aca
    cprf_aca <- peripheral1_aca / vp_aca
    cpls_acp <- central_acp / vc_acp
    cprf_acp <- peripheral1_acp / vp_acp
    cpls_ibr <- central_ibr / vc_ibr
    cprf_ibr <- peripheral1_ibr / vp_ibr
    cpls_pci <- central_pci / vc_ibr
    cprf_pci <- peripheral1_pci / vp_ibr

    # ---------------- PK ODEs ----------------
    # Zanubrutinib: gut -> depot -> central <-> peripheral
    d/dt(gut_zan)         <- -(1 / d1_zan) * f1_zan * gut_zan
    d/dt(depot_zan)       <-  (1 / d1_zan) * f1_zan * gut_zan - kabs_zan * depot_zan
    d/dt(central_zan)     <-  kabs_zan * depot_zan - cl_zan * cpls_zan - q_zan * (cpls_zan - cprf_zan)
    d/dt(peripheral1_zan) <-  q_zan * (cpls_zan - cprf_zan)

    # Acalabrutinib: gut -> 5 transit compartments -> depot -> central <-> peripheral
    d/dt(gut_aca)       <- -ktr_aca * f1_aca * gut_aca
    d/dt(transit1_aca)  <-  ktr_aca * f1_aca * gut_aca - ktr_aca * transit1_aca
    d/dt(transit2_aca)  <-  ktr_aca * transit1_aca - ktr_aca * transit2_aca
    d/dt(transit3_aca)  <-  ktr_aca * transit2_aca - ktr_aca * transit3_aca
    d/dt(transit4_aca)  <-  ktr_aca * transit3_aca - ktr_aca * transit4_aca
    d/dt(transit5_aca)  <-  ktr_aca * transit4_aca - ktr_aca * transit5_aca
    d/dt(depot_aca)     <-  ktr_aca * transit5_aca - kabs_aca * depot_aca
    d/dt(central_aca)     <-  kabs_aca * depot_aca - cl_aca * (1 - fm_aca) * cpls_aca -
      cl_aca * fm_aca * cpls_aca - q_aca * (cpls_aca - cprf_aca)
    d/dt(peripheral1_aca) <-  q_aca * (cpls_aca - cprf_aca)
    d/dt(central_acp)     <-  cl_aca * fm_aca * cpls_aca - cl_acp * cpls_acp -
      q_acp * (cpls_acp - cprf_acp)
    d/dt(peripheral1_acp) <-  q_acp * (cpls_acp - cprf_acp)

    # Ibrutinib: gut -> depot -> central <-> peripheral, with PCI-45227 formed
    # both pre-systemically (from the depot) and from systemic ibrutinib
    d/dt(gut_ibr)         <- -(1 / d1_ibr) * f1_ibr * gut_ibr
    d/dt(depot_ibr)       <-  (1 / d1_ibr) * f1_ibr * gut_ibr - kabs_ibr * depot_ibr -
      kabs_pci * depot_ibr
    d/dt(central_ibr)     <-  kabs_ibr * depot_ibr - cl_ibr * cpls_ibr -
      clmet_ibr * cpls_ibr - q_ibr * (cpls_ibr - cprf_ibr)
    d/dt(peripheral1_ibr) <-  q_ibr * (cpls_ibr - cprf_ibr)
    d/dt(central_pci)     <-  kabs_pci * depot_ibr + clmet_ibr * cpls_ibr -
      cl_pci * cpls_pci - q_pci * (cpls_pci - cprf_pci)
    d/dt(peripheral1_pci) <-  q_pci * (cpls_pci - cprf_pci)

    alag(gut_ibr) <- alag_ibr

    # ---------------- Perfusion-limited transport into LN and BM ISF (equation 1) ----------------
    cisf_ln_zan <- isf_ln_zan / v_isf_ln
    cisf_ln_aca <- isf_ln_aca / v_isf_ln
    cisf_ln_acp <- isf_ln_acp / v_isf_ln
    cisf_ln_ibr <- isf_ln_ibr / v_isf_ln
    cisf_ln_pci <- isf_ln_pci / v_isf_ln
    cisf_bm_zan <- isf_bm_zan / v_isf_bm
    cisf_bm_aca <- isf_bm_aca / v_isf_bm
    cisf_bm_acp <- isf_bm_acp / v_isf_bm
    cisf_bm_ibr <- isf_bm_ibr / v_isf_bm
    cisf_bm_pci <- isf_bm_pci / v_isf_bm
    d/dt(isf_ln_zan) <- q_ln * (cpls_zan * fup_zan - cisf_ln_zan / kp_isf_zan)
    d/dt(isf_ln_aca) <- q_ln * (cpls_aca * fup_aca - cisf_ln_aca / kp_isf_aca)
    d/dt(isf_ln_acp) <- q_ln * (cpls_acp * fup_acp - cisf_ln_acp / kp_isf_acp)
    d/dt(isf_ln_ibr) <- q_ln * (cpls_ibr * fup_ibr - cisf_ln_ibr / kp_isf_ibr)
    d/dt(isf_ln_pci) <- q_ln * (cpls_pci * fup_pci - cisf_ln_pci / kp_isf_pci)
    d/dt(isf_bm_zan) <- q_bm * (cpls_zan * fup_zan - cisf_bm_zan / kp_isf_zan)
    d/dt(isf_bm_aca) <- q_bm * (cpls_aca * fup_aca - cisf_bm_aca / kp_isf_aca)
    d/dt(isf_bm_acp) <- q_bm * (cpls_acp * fup_acp - cisf_bm_acp / kp_isf_acp)
    d/dt(isf_bm_ibr) <- q_bm * (cpls_ibr * fup_ibr - cisf_bm_ibr / kp_isf_ibr)
    d/dt(isf_bm_pci) <- q_bm * (cpls_pci * fup_pci - cisf_bm_pci / kp_isf_pci)

    # ---------------- Intracellular concentrations in CLL cells (nmol/L, equation 2) ----------------
    # PBMC: driven by total plasma concentration via the erythrocyte-to-plasma
    # partition coefficient; LN and BM: driven by unbound ISF concentration.
    ccell_pbmc_zan <- cpls_zan * kp_rbc_zan * c_adj_pbmc * 1e6 / mr_zan
    ccell_pbmc_aca <- cpls_aca * kp_rbc_aca * c_adj_pbmc * 1e6 / mr_aca
    ccell_pbmc_acp <- cpls_acp * kp_rbc_acp * c_adj_pbmc * 1e6 / mr_acp
    ccell_pbmc_ibr <- cpls_ibr * kp_rbc_ibr * c_adj_pbmc * 1e6 / mr_ibr
    ccell_pbmc_pci <- cpls_pci * kp_rbc_pci * c_adj_pbmc * 1e6 / mr_pci
    ccell_ln_zan <- cisf_ln_zan * kp_cell_zan * 1e6 / mr_zan
    ccell_ln_aca <- cisf_ln_aca * kp_cell_aca * 1e6 / mr_aca
    ccell_ln_acp <- cisf_ln_acp * kp_cell_acp * 1e6 / mr_acp
    ccell_ln_ibr <- cisf_ln_ibr * kp_cell_ibr * 1e6 / mr_ibr
    ccell_ln_pci <- cisf_ln_pci * kp_cell_pci * 1e6 / mr_pci
    ccell_bm_zan <- cisf_bm_zan * kp_cell_zan * 1e6 / mr_zan
    ccell_bm_aca <- cisf_bm_aca * kp_cell_aca * 1e6 / mr_aca
    ccell_bm_acp <- cisf_bm_acp * kp_cell_acp * 1e6 / mr_acp
    ccell_bm_ibr <- cisf_bm_ibr * kp_cell_ibr * 1e6 / mr_ibr
    ccell_bm_pci <- cisf_bm_pci * kp_cell_pci * 1e6 / mr_pci

    # ---------------- BTK turnover, ATP binding and covalent inhibition (equations 3-7) ----------------
    # PBMC
    vsyn_btk_pbmc <- bl_btk * kdeg_pbmc
    rbind_zan_pbmc <- koff_zan * (ccell_pbmc_zan * btk_free_pbmc / kd_zan - btk_zan_pbmc)
    rbind_aca_pbmc <- koff_aca * (ccell_pbmc_aca * btk_free_pbmc / kd_aca - btk_aca_pbmc)
    rbind_acp_pbmc <- koff_acp * (ccell_pbmc_acp * btk_free_pbmc / kd_acp - btk_acp_pbmc)
    rbind_ibr_pbmc <- koff_ibr * (ccell_pbmc_ibr * btk_free_pbmc / kd_ibr - btk_ibr_pbmc)
    rbind_pci_pbmc <- koff_pci * (ccell_pbmc_pci * btk_free_pbmc / kd_pci - btk_pci_pbmc)
    rbind_atp_pbmc <- koff_atp * (atp * btk_free_pbmc / kd_atp - btk_atp_pbmc)
    d/dt(btk_free_pbmc) <- vsyn_btk_pbmc - kdeg_pbmc * btk_free_pbmc - rbind_atp_pbmc -
      rbind_zan_pbmc -
      rbind_aca_pbmc -
      rbind_acp_pbmc -
      rbind_ibr_pbmc -
      rbind_pci_pbmc
    d/dt(btk_atp_pbmc) <- rbind_atp_pbmc - kdeg_pbmc * btk_atp_pbmc
    d/dt(btk_zan_pbmc) <- rbind_zan_pbmc - kdeg_pbmc * btk_zan_pbmc - kinact_zan * btk_zan_pbmc
    d/dt(btk_zan_inact_pbmc) <- kinact_zan * btk_zan_pbmc - kdeg_pbmc * btk_zan_inact_pbmc
    d/dt(btk_aca_pbmc) <- rbind_aca_pbmc - kdeg_pbmc * btk_aca_pbmc - kinact_aca * btk_aca_pbmc
    d/dt(btk_aca_inact_pbmc) <- kinact_aca * btk_aca_pbmc - kdeg_pbmc * btk_aca_inact_pbmc
    d/dt(btk_acp_pbmc) <- rbind_acp_pbmc - kdeg_pbmc * btk_acp_pbmc - kinact_acp * btk_acp_pbmc
    d/dt(btk_acp_inact_pbmc) <- kinact_acp * btk_acp_pbmc - kdeg_pbmc * btk_acp_inact_pbmc
    d/dt(btk_ibr_pbmc) <- rbind_ibr_pbmc - kdeg_pbmc * btk_ibr_pbmc - kinact_ibr * btk_ibr_pbmc
    d/dt(btk_ibr_inact_pbmc) <- kinact_ibr * btk_ibr_pbmc - kdeg_pbmc * btk_ibr_inact_pbmc
    d/dt(btk_pci_pbmc) <- rbind_pci_pbmc - kdeg_pbmc * btk_pci_pbmc - kinact_pci * btk_pci_pbmc
    d/dt(btk_pci_inact_pbmc) <- kinact_pci * btk_pci_pbmc - kdeg_pbmc * btk_pci_inact_pbmc

    # lymph node
    vsyn_btk_ln <- bl_btk * kdeg_ln
    rbind_zan_ln <- koff_zan * (ccell_ln_zan * btk_free_ln / kd_zan - btk_zan_ln)
    rbind_aca_ln <- koff_aca * (ccell_ln_aca * btk_free_ln / kd_aca - btk_aca_ln)
    rbind_acp_ln <- koff_acp * (ccell_ln_acp * btk_free_ln / kd_acp - btk_acp_ln)
    rbind_ibr_ln <- koff_ibr * (ccell_ln_ibr * btk_free_ln / kd_ibr - btk_ibr_ln)
    rbind_pci_ln <- koff_pci * (ccell_ln_pci * btk_free_ln / kd_pci - btk_pci_ln)
    rbind_atp_ln <- koff_atp * (atp * btk_free_ln / kd_atp - btk_atp_ln)
    d/dt(btk_free_ln) <- vsyn_btk_ln - kdeg_ln * btk_free_ln - rbind_atp_ln -
      rbind_zan_ln -
      rbind_aca_ln -
      rbind_acp_ln -
      rbind_ibr_ln -
      rbind_pci_ln
    d/dt(btk_atp_ln) <- rbind_atp_ln - kdeg_ln * btk_atp_ln
    d/dt(btk_zan_ln) <- rbind_zan_ln - kdeg_ln * btk_zan_ln - kinact_zan * btk_zan_ln
    d/dt(btk_zan_inact_ln) <- kinact_zan * btk_zan_ln - kdeg_ln * btk_zan_inact_ln
    d/dt(btk_aca_ln) <- rbind_aca_ln - kdeg_ln * btk_aca_ln - kinact_aca * btk_aca_ln
    d/dt(btk_aca_inact_ln) <- kinact_aca * btk_aca_ln - kdeg_ln * btk_aca_inact_ln
    d/dt(btk_acp_ln) <- rbind_acp_ln - kdeg_ln * btk_acp_ln - kinact_acp * btk_acp_ln
    d/dt(btk_acp_inact_ln) <- kinact_acp * btk_acp_ln - kdeg_ln * btk_acp_inact_ln
    d/dt(btk_ibr_ln) <- rbind_ibr_ln - kdeg_ln * btk_ibr_ln - kinact_ibr * btk_ibr_ln
    d/dt(btk_ibr_inact_ln) <- kinact_ibr * btk_ibr_ln - kdeg_ln * btk_ibr_inact_ln
    d/dt(btk_pci_ln) <- rbind_pci_ln - kdeg_ln * btk_pci_ln - kinact_pci * btk_pci_ln
    d/dt(btk_pci_inact_ln) <- kinact_pci * btk_pci_ln - kdeg_ln * btk_pci_inact_ln

    # bone marrow
    vsyn_btk_bm <- bl_btk * kdeg_bm
    rbind_zan_bm <- koff_zan * (ccell_bm_zan * btk_free_bm / kd_zan - btk_zan_bm)
    rbind_aca_bm <- koff_aca * (ccell_bm_aca * btk_free_bm / kd_aca - btk_aca_bm)
    rbind_acp_bm <- koff_acp * (ccell_bm_acp * btk_free_bm / kd_acp - btk_acp_bm)
    rbind_ibr_bm <- koff_ibr * (ccell_bm_ibr * btk_free_bm / kd_ibr - btk_ibr_bm)
    rbind_pci_bm <- koff_pci * (ccell_bm_pci * btk_free_bm / kd_pci - btk_pci_bm)
    rbind_atp_bm <- koff_atp * (atp * btk_free_bm / kd_atp - btk_atp_bm)
    d/dt(btk_free_bm) <- vsyn_btk_bm - kdeg_bm * btk_free_bm - rbind_atp_bm -
      rbind_zan_bm -
      rbind_aca_bm -
      rbind_acp_bm -
      rbind_ibr_bm -
      rbind_pci_bm
    d/dt(btk_atp_bm) <- rbind_atp_bm - kdeg_bm * btk_atp_bm
    d/dt(btk_zan_bm) <- rbind_zan_bm - kdeg_bm * btk_zan_bm - kinact_zan * btk_zan_bm
    d/dt(btk_zan_inact_bm) <- kinact_zan * btk_zan_bm - kdeg_bm * btk_zan_inact_bm
    d/dt(btk_aca_bm) <- rbind_aca_bm - kdeg_bm * btk_aca_bm - kinact_aca * btk_aca_bm
    d/dt(btk_aca_inact_bm) <- kinact_aca * btk_aca_bm - kdeg_bm * btk_aca_inact_bm
    d/dt(btk_acp_bm) <- rbind_acp_bm - kdeg_bm * btk_acp_bm - kinact_acp * btk_acp_bm
    d/dt(btk_acp_inact_bm) <- kinact_acp * btk_acp_bm - kdeg_bm * btk_acp_inact_bm
    d/dt(btk_ibr_bm) <- rbind_ibr_bm - kdeg_bm * btk_ibr_bm - kinact_ibr * btk_ibr_bm
    d/dt(btk_ibr_inact_bm) <- kinact_ibr * btk_ibr_bm - kdeg_bm * btk_ibr_inact_bm
    d/dt(btk_pci_bm) <- rbind_pci_bm - kdeg_bm * btk_pci_bm - kinact_pci * btk_pci_bm
    d/dt(btk_pci_inact_bm) <- kinact_pci * btk_pci_bm - kdeg_bm * btk_pci_inact_bm

    # ---------------- Initial conditions ----------------
    btk_free_pbmc(0) <- bl_btk
    btk_free_ln(0) <- bl_btk
    btk_free_bm(0) <- bl_btk

    # ---------------- Outputs ----------------
    # Plasma concentrations in ng/mL
    Cc_zan <- cpls_zan * 1000
    Cc_aca <- cpls_aca * 1000
    Cc_acp <- cpls_acp * 1000
    Cc_ibr <- cpls_ibr * 1000
    Cc_pci <- cpls_pci * 1000

    # BTK occupancy (%), i.e. the fraction of total BTK bound to any inhibitor
    # or its inactivated adduct. The 1e-10 term reproduces the source model's
    # antizero_denom_nM guard against a zero denominator.
    btk_bound_pbmc <- btk_zan_pbmc + btk_zan_inact_pbmc +
      btk_aca_pbmc + btk_aca_inact_pbmc +
      btk_acp_pbmc + btk_acp_inact_pbmc +
      btk_ibr_pbmc + btk_ibr_inact_pbmc +
      btk_pci_pbmc + btk_pci_inact_pbmc
    btk_total_pbmc <- btk_bound_pbmc + btk_free_pbmc + btk_atp_pbmc + 1e-10
    occupancy_pbmc <- 100 * btk_bound_pbmc / btk_total_pbmc
    btk_bound_ln <- btk_zan_ln + btk_zan_inact_ln +
      btk_aca_ln + btk_aca_inact_ln +
      btk_acp_ln + btk_acp_inact_ln +
      btk_ibr_ln + btk_ibr_inact_ln +
      btk_pci_ln + btk_pci_inact_ln
    btk_total_ln <- btk_bound_ln + btk_free_ln + btk_atp_ln + 1e-10
    occupancy_ln <- 100 * btk_bound_ln / btk_total_ln
    btk_bound_bm <- btk_zan_bm + btk_zan_inact_bm +
      btk_aca_bm + btk_aca_inact_bm +
      btk_acp_bm + btk_acp_inact_bm +
      btk_ibr_bm + btk_ibr_inact_bm +
      btk_pci_bm + btk_pci_inact_bm
    btk_total_bm <- btk_bound_bm + btk_free_bm + btk_atp_bm + 1e-10
    occupancy_bm <- 100 * btk_bound_bm / btk_total_bm
  })
}
