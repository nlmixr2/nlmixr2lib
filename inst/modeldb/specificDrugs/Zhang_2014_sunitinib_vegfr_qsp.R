Zhang_2014_sunitinib_vegfr_qsp <- function() {
  description <- "QSP (human, deterministic typical-value tumor-cell model). PK / enhanced-PD (ePD) model of VEGFR2 signaling in a tumor, linked to two-compartment PK of up to five anticancer drugs (Zhang 2014). The base VEGFR network (40 biochemical reactions, mass-action and Michaelis-Menten kinetics, literature-derived parameters) covers VEGF-VEGFR binding, receptor dimerisation and autophosphorylation; recruitment of PI3K, PLC-gamma and Grb2-Sos to the phosphorylated dimer; PIP2/PIP3 turnover and PDK/Akt activation; DAG/IP3, PKC, RKIP and RasGRP signaling; the Ras-Raf-MEK-ERK-RSK cascade; and three negative feedbacks (RSK-mediated Raf inactivation, Akt-mediated PI3K inactivation, and Akt-induced inhibition of receptor synthesis). Drugs: sunitinib (oral, binds free, ligand-bound and dimerised receptor competitively, preventing phosphorylation), bevacizumab (IV infusion, sequesters VEGF), buparlisib / BKM120 (oral, competitive PI3K inhibitor), selumetinib / AZD6244 (oral, non-competitive inhibitor of ERK phosphorylation by MEK) and a hypothetical PLC-gamma inhibitor (IV, competitive inhibitor of PIP2 hydrolysis). Outputs are active ERK (pERK, ErkP) and active Akt (pAkt, PIP3-Akt-P) in nM and the plasma concentration of each drug in nM. Time in minutes. The system must be equilibrated before dosing (the authors ran 20 days); see the vignette."
  reference <- "Zhang X-Y, Birtwistle MR, Gallo JM. A General Network Pharmacodynamic Model-Based Design Pipeline for Customized Cancer Therapy Applied to the VEGFR Pathway. CPT Pharmacometrics Syst Pharmacol. 2014;3(1):e92. doi:10.1038/psp.2013.65. Parameter values and equations follow the authors' published MATLAB code (Supplementary Data: VEGFRModelDrugs.m, GetInitialConditions.m, physiosetbasis.mat row 1); Supplementary Tables S1-S4."
  vignette <- "Zhang_2014_sunitinib_vegfr_qsp"

  paper_specific_compartments <- c(
    "vr",
    "vegf_vr",
    "vegf_vr2",
    "vegf_vr2_p",
    "plcg",
    "pi3k_p",
    "pi3k",
    "pip3",
    "pip2",
    "pdk",
    "akt",
    "pip3_pdk",
    "pip3_akt",
    "pip3_akt_p",
    "pkc",
    "dag_pkc",
    "raf",
    "rkip",
    "raf_rkip",
    "rkip_p",
    "raf_p",
    "raf_i",
    "rsk_p",
    "rsk",
    "erk_p",
    "erk",
    "grb2sos",
    "ras_gtp",
    "ras_gdp",
    "mek_p",
    "mek",
    "vegf_vr2_p_pi3k",
    "vegf_vr2_p_plcg",
    "vegf_vr2_p_grb2sos",
    "ip3",
    "dag",
    "rasgrp",
    "dag_rasgrp",
    "sunitinib_vr",
    "sunitinib_vegf_vr",
    "sunitinib_vegf_vr2",
    "bevacizumab_vegf"
  )

  units <- list(
    time = "min",
    dosing = "nmol",
    concentration = "nM"
  )

  # depot: oral sunitinib; central_bevacizumab: 30-min IV infusion;
  # depot_buparlisib, depot_selumetinib: oral; central_plcginhibitor: IV bolus
  # entered as the resulting plasma concentration (nM).
  dosing <- c("depot", "central_bevacizumab", "depot_buparlisib", "depot_selumetinib", "central_plcginhibitor")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib", units = "nmol", specimen = "tissue", verified = TRUE),
    central_bevacizumab = list(analyte = "bevacizumab", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_bevacizumab = list(analyte = "bevacizumab", units = "nmol", specimen = "tissue", verified = TRUE),
    depot_buparlisib = list(analyte = "buparlisib", units = "nmol", specimen = "administration site", verified = TRUE),
    central_buparlisib = list(analyte = "buparlisib", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_buparlisib = list(analyte = "buparlisib", units = "nmol", specimen = "tissue", verified = TRUE),
    depot_selumetinib = list(
      analyte = "selumetinib",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_selumetinib = list(analyte = "selumetinib", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_selumetinib = list(analyte = "selumetinib", units = "nmol", specimen = "tissue", verified = TRUE),
    central_plcginhibitor = list(
      analyte = "PLC-gamma inhibitor (hypothetical)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    vr = list(analyte = "free VEGFR2", units = "nM", specimen = "tumor", verified = TRUE),
    vegf = list(analyte = "free VEGF-A (tumor extracellular)", units = "nM", specimen = "tumor", verified = TRUE),
    vegf_vr = list(analyte = "VEGF-VEGFR2 complex", units = "nM", specimen = "tumor", verified = TRUE),
    vegf_vr2 = list(analyte = "VEGF-VEGFR2 dimer", units = "nM", specimen = "tumor", verified = TRUE),
    vegf_vr2_p = list(analyte = "phosphorylated VEGF-VEGFR2 dimer", units = "nM", specimen = "tumor", verified = TRUE),
    plcg = list(analyte = "PLC-gamma", units = "nM", specimen = "tumor", verified = TRUE),
    pi3k_p = list(analyte = "phosphorylated (inactivated) PI3K", units = "nM", specimen = "tumor", verified = TRUE),
    pi3k = list(analyte = "PI3K", units = "nM", specimen = "tumor", verified = TRUE),
    pip3 = list(analyte = "PIP3", units = "nM", specimen = "tumor", verified = TRUE),
    pip2 = list(analyte = "PIP2", units = "nM", specimen = "tumor", verified = TRUE),
    pdk = list(analyte = "PDK1", units = "nM", specimen = "tumor", verified = TRUE),
    akt = list(analyte = "Akt", units = "nM", specimen = "tumor", verified = TRUE),
    pip3_pdk = list(analyte = "PIP3-PDK1 complex", units = "nM", specimen = "tumor", verified = TRUE),
    pip3_akt = list(analyte = "PIP3-Akt complex", units = "nM", specimen = "tumor", verified = TRUE),
    pip3_akt_p = list(
      analyte = "phosphorylated PIP3-Akt (active Akt, pAkt)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    pkc = list(analyte = "PKC", units = "nM", specimen = "tumor", verified = TRUE),
    dag_pkc = list(analyte = "DAG-PKC complex", units = "nM", specimen = "tumor", verified = TRUE),
    raf = list(analyte = "Raf", units = "nM", specimen = "tumor", verified = TRUE),
    rkip = list(analyte = "RKIP", units = "nM", specimen = "tumor", verified = TRUE),
    raf_rkip = list(analyte = "Raf-RKIP complex", units = "nM", specimen = "tumor", verified = TRUE),
    rkip_p = list(analyte = "phosphorylated RKIP", units = "nM", specimen = "tumor", verified = TRUE),
    raf_p = list(analyte = "phosphorylated (active) Raf", units = "nM", specimen = "tumor", verified = TRUE),
    raf_i = list(analyte = "feedback-inactivated Raf", units = "nM", specimen = "tumor", verified = TRUE),
    rsk_p = list(analyte = "phosphorylated RSK", units = "nM", specimen = "tumor", verified = TRUE),
    rsk = list(analyte = "RSK", units = "nM", specimen = "tumor", verified = TRUE),
    erk_p = list(analyte = "phosphorylated ERK (pERK)", units = "nM", specimen = "tumor", verified = TRUE),
    erk = list(analyte = "ERK", units = "nM", specimen = "tumor", verified = TRUE),
    grb2sos = list(analyte = "Grb2-Sos", units = "nM", specimen = "tumor", verified = TRUE),
    ras_gtp = list(analyte = "Ras-GTP", units = "nM", specimen = "tumor", verified = TRUE),
    ras_gdp = list(analyte = "Ras-GDP", units = "nM", specimen = "tumor", verified = TRUE),
    mek_p = list(analyte = "phosphorylated MEK", units = "nM", specimen = "tumor", verified = TRUE),
    mek = list(analyte = "MEK", units = "nM", specimen = "tumor", verified = TRUE),
    vegf_vr2_p_pi3k = list(
      analyte = "phosphorylated receptor dimer-PI3K complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    vegf_vr2_p_plcg = list(
      analyte = "phosphorylated receptor dimer-PLC-gamma complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    vegf_vr2_p_grb2sos = list(
      analyte = "phosphorylated receptor dimer-Grb2-Sos complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    ip3 = list(analyte = "IP3", units = "nM", specimen = "tumor", verified = TRUE),
    dag = list(analyte = "DAG", units = "nM", specimen = "tumor", verified = TRUE),
    rasgrp = list(analyte = "RasGRP", units = "nM", specimen = "tumor", verified = TRUE),
    dag_rasgrp = list(analyte = "DAG-RasGRP complex", units = "nM", specimen = "tumor", verified = TRUE),
    sunitinib_vr = list(analyte = "sunitinib-VEGFR2 complex", units = "nM", specimen = "tumor", verified = TRUE),
    sunitinib_vegf_vr = list(
      analyte = "sunitinib-VEGF-VEGFR2 complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    sunitinib_vegf_vr2 = list(
      analyte = "sunitinib-VEGF-VEGFR2 dimer complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    bevacizumab_vegf = list(analyte = "bevacizumab-VEGF complex", units = "nM", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "Generic solid tumor driven by VEGFR2 signaling (a 20-mL tumor, half cells by volume); base network plus NF1-loss, activating-PI3K and PTEN-loss mutation scenarios",
    dose_range = "Sunitinib 10-100 mg PO QD; bevacizumab 350-1050 mg IV over 30 min, at most two doses per 28 days; buparlisib (BKM120) 50-100 mg PO QD; selumetinib (AZD6244) 50-100 mg PO Q12H; PLC-gamma inhibitor 2-10 umol/L IV QD (Table 1 optimization bounds)",
    notes = "Deterministic typical-value systems model; no subjects were fitted. Signaling parameters and protein abundances were taken from the literature (Supplementary Tables S1-S2). Sunitinib and bevacizumab PK are published human population PK models (Houk 2009; Li 2013); buparlisib and selumetinib PK were fitted by the authors to digitised mean plasma profiles (Bendell 2012; O'Neil 2011). The PLC-gamma inhibitor is hypothetical (Ki 1 umol/L, first-order loss 0.000963 1/min)."
  )

  ini({
    # ---- Base VEGFR network (Table S1; values from the authors' code, physiosetbasis.mat row 1, k(1)..k(77)) ----
    lk1f <- fixed(log(0.48)); label("VEGF-VEGFR2 association rate constant k1f (1/(nM*min))") # mat k(1) = 0.48; Table S1 rxn 1 k1f = 4.8x10-1
    lk1r <- fixed(log(0.06)); label("VEGF-VEGFR2 dissociation rate constant k1r (1/min)") # mat k(2) = 0.06; Table S1 rxn 1 k1r = 6.0x10-2
    lk2f <- fixed(log(0.036)); label("VEGF-VEGFR2 dimerisation rate constant k2f (1/(nM*min))") # mat k(3) = 0.036; Table S1 rxn 2 k2f = 3.6x10-2
    lk2r <- fixed(log(0.06)); label("Receptor dimer dissociation rate constant k2r (1/min)") # mat k(4) = 0.06; Table S1 rxn 2 k2r = 6.0x10-2
    lk3f <- fixed(log(60)); label("Receptor dimer autophosphorylation rate constant k3f (1/min)") # mat k(5) = 60; Table S1 rxn 3 k3f = 6.0x10^1
    lk3r <- fixed(log(0.6)); label("Receptor dimer dephosphorylation rate constant k3r (1/min)") # mat k(6) = 0.6; Table S1 rxn 3 k3r = 6.0x10-1
    lk4f <- fixed(log(0.0912)); label("PLC-gamma binding to phosphorylated dimer k4f (1/(nM*min))") # mat k(11) = 0.0912; Table S1 rxn 4 k4f = 9.1x10-2
    lk4r <- fixed(log(6)); label("PLC-gamma unbinding from phosphorylated dimer k4r (1/min)") # mat k(12) = 6; Table S1 rxn 4 k4r = 6.0
    lkpg <- fixed(log(6)); label("PLC-gamma catalysed PIP2 hydrolysis kcat kpg (1/min)") # mat k(13) = 6; Table S1 kpg = 6.0
    lkmpg <- fixed(log(500)); label("Km of PLC-gamma catalysed PIP2 hydrolysis Kmpg (nM)") # mat k(14) = 500; Table S1 prints Kmpg = 2.0x10^3 (code value used)
    lk6f <- fixed(log(0.84)); label("PI3K binding to phosphorylated dimer k6f (1/(nM*min))") # mat k(15) = 0.84; Table S1 rxn 6 k6f = 8.4x10-1
    lk6r <- fixed(log(0.2)); label("PI3K unbinding from phosphorylated dimer k6r (1/min)") # mat k(16) = 0.2; Table S1 prints k6r = 1.0x10^1 (code value used)
    lk7 <- fixed(log(35.7)); label("SHP1-catalysed PI3K-P dephosphorylation kcat k7 (1/min)") # mat k(17) = 35.7; Table S1 rxn 7 k7 = 3.6x10^1
    lkm7 <- fixed(log(632)); label("Km of PI3K-P dephosphorylation Km7 (nM)") # mat k(18) = 632; Table S1 Km7 = 6.3x10^2
    lk8 <- fixed(log(1500)); label("PI3K-catalysed PIP2 phosphorylation kcat k8 (1/min)") # mat k(19) = 1500; Table S1 rxn 8 k8 = 1.5x10^3
    lkm8 <- fixed(log(1140)); label("Km of PI3K-catalysed PIP2 phosphorylation Km8 (nM)") # mat k(20) = 1140; Table S1 Km8 = 1.1x10^3
    lk9 <- fixed(log(300)); label("PTEN-catalysed PIP3 dephosphorylation kcat k9 (1/min)") # mat k(21) = 300; Table S1 rxn 9 k9 = 3.0x10^2
    lkm9 <- fixed(log(1100)); label("Km of PIP3 dephosphorylation Km9 (nM)") # mat k(22) = 1100; Table S1 Km9 = 1.1x10^3
    lk10f <- fixed(log(0.179)); label("PIP3-PDK1 association k10f (1/(nM*min))") # mat k(23) = 0.179; Table S1 rxn 10 k10f = 1.8x10-1
    lk10r <- fixed(log(60)); label("PIP3-PDK1 dissociation k10r (1/min)") # mat k(24) = 60; Table S1 k10r = 6.0x10^1
    lk11f <- fixed(log(17.9)); label("PIP3-Akt association k11f (1/(nM*min))") # mat k(25) = 17.9; Table S1 prints k11f = 1.8x10-1 (code value used)
    lk11r <- fixed(log(60)); label("PIP3-Akt dissociation k11r (1/min)") # mat k(26) = 60; Table S1 k11r = 6.0x10^1
    lk12 <- fixed(log(180)); label("PDK1-catalysed Akt phosphorylation kcat k12 (1/min)") # mat k(27) = 180; Table S1 rxn 12 k12 = 1.8x10^2
    lkm12 <- fixed(log(100)); label("Km of Akt phosphorylation Km12 (nM)") # mat k(28) = 100; Table S1 prints Km12 = 0.99 (code value used)
    lk13 <- fixed(log(90)); label("PP2A-catalysed Akt-P dephosphorylation kcat k13 (1/min)") # mat k(29) = 90; Table S1 rxn 13 k13 = 9.0x10^1
    lkm13 <- fixed(log(100)); label("Km of Akt-P dephosphorylation Km13 (nM)") # mat k(30) = 100; Table S1 prints Km13 = 0.1 (code value used)
    lk14f <- fixed(log(0.107)); label("DAG-PKC association k14f (1/(nM*min))") # mat k(31) = 0.107; Table S1 rxn 14 k14f = 1.1x10-1
    lk14r <- fixed(log(6)); label("DAG-PKC dissociation k14r (1/min)") # mat k(32) = 6; Table S1 k14r = 6.0
    lk15f <- fixed(log(0.1)); label("DAG-RasGRP association k15f (1/(nM*min))") # mat k(33) = 0.1; Table S1 rxn 15 k15f = 1.0x10-1
    lk15r <- fixed(log(0.05)); label("DAG-RasGRP dissociation k15r (1/min)") # mat k(34) = 0.05; Table S1 k15r = 5.0x10-2
    lk16f <- fixed(log(17.9)); label("Raf-RKIP association k16f (1/(nM*min))") # mat k(35) = 17.9; Table S1 prints k16f = 1.8x10-1 (code value used)
    lk16r <- fixed(log(0.3)); label("Raf-RKIP dissociation k16r (1/min)") # mat k(36) = 0.3; Table S1 k16r = 3.0x10-1
    lk17 <- fixed(log(180)); label("DAG-PKC catalysed RKIP phosphorylation kcat k17 (1/min)") # mat k(37) = 180; Table S1 rxn 17 k17 = 1.8x10^2
    lkm17 <- fixed(log(301)); label("Km of RKIP phosphorylation Km17 (nM)") # mat k(38) = 301; Table S1 Km17 = 3.0x10^2
    lvm18 <- fixed(log(60.3)); label("Vmax of RKIP-P dephosphorylation Vm18 (nM/min)") # mat k(39) = 60.3; Table S1 prints Vm18 = 1.8x10^1 (code value used)
    lkm18 <- fixed(log(99.9)); label("Km of RKIP-P dephosphorylation Km18 (nM)") # mat k(40) = 99.9; Table S1 Km18 = 9.9x10^1
    lk19 <- fixed(log(27)); label("Ras-GTP catalysed Raf activation kcat k19 (1/min)") # mat k(41) = 27; Table S1 rxn 19 k19 = 2.7x10^1
    lkm19 <- fixed(log(41.5)); label("Km of Raf activation Km19 (nM)") # mat k(42) = 41.5; Table S1 Km19 = 4.2x10^1
    lvm20 <- fixed(log(161)); label("Vmax of Raf-P dephosphorylation Vm20 (nM/min)") # mat k(43) = 161; Table S1 Vm20 = 1.6x10^2
    lkm20 <- fixed(log(9.96)); label("Km of Raf-P dephosphorylation Km20 (nM)") # mat k(44) = 9.96; Table S1 Km20 = 9.96
    lk21 <- fixed(log(50)); label("Raf-P catalysed MEK phosphorylation kcat k21 (1/min)") # mat k(45) = 50; Table S1 rxn 21 k21 = 5.0x10^1
    lkm21 <- fixed(log(14.9)); label("Km of MEK phosphorylation Km21 (nM)") # mat k(46) = 14.9; Table S1 Km21 = 1.5x10^1
    lvm22 <- fixed(log(1530)); label("Vmax of MEK-P dephosphorylation Vm22 (nM/min)") # mat k(47) = 1530; Table S1 Vm22 = 1.5x10^3
    lkm22 <- fixed(log(996)); label("Km of MEK-P dephosphorylation Km22 (nM)") # mat k(48) = 996; Table S1 Km22 = 9.96x10^2
    lk23 <- fixed(log(8.3)); label("MEK-P catalysed ERK phosphorylation kcat k23 (1/min)") # mat k(49) = 8.3; Table S1 rxn 23 k23 = 8.3
    lkm23 <- fixed(log(14.9)); label("Km of ERK phosphorylation Km23 (nM)") # mat k(50) = 14.9; Table S1 Km23 = 1.5x10^1
    lvm24 <- fixed(log(332)); label("Vmax of ERK-P dephosphorylation Vm24 (nM/min)") # mat k(51) = 332; Table S1 Vm24 = 3.3x10^2
    lkm24 <- fixed(log(996)); label("Km of ERK-P dephosphorylation Km24 (nM)") # mat k(52) = 996; Table S1 Km24 = 9.96x10^2
    lk25 <- fixed(log(18)); label("ERK-P catalysed RSK phosphorylation kcat k25 (1/min)") # mat k(53) = 18; Table S1 rxn 25 k25 = 1.8x10^1
    lkm25 <- fixed(log(9.99)); label("Km of RSK phosphorylation Km25 (nM)") # mat k(54) = 9.99; Table S1 Km25 = 9.99
    lk26 <- fixed(log(18)); label("RSK-P dephosphorylation rate constant k26 (1/min)") # mat k(55) = 18; Table S1 rxn 26 k26 = 1.8x10^1
    lk27f <- fixed(log(0.06)); label("Grb2-Sos binding to phosphorylated dimer k27f (1/(nM*min))") # mat k(56) = 0.06; Table S1 rxn 27 k27f = 6.0x10-2
    lk27r <- fixed(log(2.7)); label("Grb2-Sos unbinding from phosphorylated dimer k27r (1/min)") # mat k(57) = 2.7; Table S1 k27r = 2.7
    lkm28 <- fixed(log(200)); label("Km of Ras-GDP activation by Sos and RasGRP Km28 (nM)") # mat k(58) = 200; Table S1 prints Km28 = 2.56 (code value used)
    lk28 <- fixed(log(15)); label("Sos / RasGRP catalysed Ras activation kcat k28 (1/min)") # mat k(59) = 15; Table S1 rxn 28 k28 = 1.5x10^1
    lk29 <- fixed(log(130)); label("NF1 (GAP) catalysed Ras-GTP inactivation kcat k29 (1/min)") # mat k(60) = 130; Table S1 prints k29 = 1.0x10-2 (code value used)
    lkm29 <- fixed(log(259)); label("Km of Ras-GTP inactivation Km29 (nM)") # mat k(61) = 259; Table S1 Km29 = 2.6x10^2
    lk32 <- fixed(log(2000)); label("Akt-P catalysed PI3K inactivation kcat k32 (1/min)") # mat k(62) = 2000; Table S1 rxn 32 k32 = 2.0x10^3
    lkm32 <- fixed(log(2000)); label("Km of PI3K inactivation Km32 (nM)") # mat k(63) = 2000; Table S1 prints Km32 = 2.4x10^3 (code value used)
    lk30 <- fixed(log(27)); label("RSK-P catalysed Raf inactivation kcat k30 (1/min)") # mat k(64) = 27; Table S1 rxn 30 k30 = 2.7x10^1
    lkm30 <- fixed(log(41.5)); label("Km of Raf inactivation Km30 (nM)") # mat k(65) = 41.5; Table S1 Km30 = 4.2x10^1
    lvm31 <- fixed(log(161)); label("Vmax of Raf reactivation Vm31 (nM/min)") # mat k(66) = 161; Table S1 Vm31 = 1.6x10^2
    lkm31 <- fixed(log(9.96)); label("Km of Raf reactivation Km31 (nM)") # mat k(67) = 9.96; Table S1 Km31 = 9.96
    lksyn <- fixed(log(1.66)); label("Maximal VEGFR2 synthesis rate ksyn (nM/min)") # mat k(68) = 1.66; Table S1 ksyn = 1.7
    lkdeg <- fixed(log(0.01)); label("Degradation rate constant of unphosphorylated receptor species kdeg (1/min)") # mat k(69) = 0.01; Table S1 kdeg = 0.01; Table S3 kdeg = 0.01 for the sunitinib complexes
    lkdeg2 <- fixed(log(0.1)); label("Degradation rate constant of the phosphorylated receptor dimer kdeg2 (1/min)") # mat k(70) = 0.1; Table S1 kdeg2 = 0.1
    lklsyn <- fixed(log(0.00225)); label("VEGF synthesis rate in the tumor interstitium klsyn (nM/min)") # mat k(71) = 0.00225; Table S1 klsyn = 2.3x10-3
    lkpip2syn <- fixed(log(2.88)); label("PIP2 synthesis rate kPIP2syn (nM/min)") # mat k(72) = 2.88; Table S1 ksynpip2 = 2.9
    lkpip2deg <- fixed(log(0.00288)); label("PIP2 degradation rate constant kPIP2deg (1/min)") # mat k(73) = 0.00288; Table S1 Kdegpip2 = 2.9x10-3
    lkdegdag <- fixed(log(0.065)); label("DAG degradation rate constant kdegdag (1/min)") # mat k(74) = 0.065; Table S1 kdegdag = 6.5x10-2
    lkdegip3 <- fixed(log(4.8)); label("IP3 degradation rate constant kdegip3 (1/min)") # mat k(75) = 4.8; Table S1 kdegip3 = 4.8
    lkh <- fixed(log(50)); label("pAkt concentration giving half-maximal inhibition of receptor synthesis kh (nM)") # mat k(76) = 50; Table S1 kh = 50
    n_fb <- fixed(3); label("Hill coefficient of pAkt feedback on receptor synthesis n (unitless)") # mat k(77) = 3; Table S1 lists n as unitless without a printed value

    # ---- Constant enzyme abundances and tumor geometry (VEGFRModelDrugs.m; Table S1) ----
    shp1 <- fixed(100); label("SHP1 phosphatase concentration (nM)") # VEGFRModelDrugs.m SHP1 = 100; Table S1 rxn 7 SHP1 = 100
    pten <- fixed(100); label("PTEN phosphatase concentration (nM)") # VEGFRModelDrugs.m PTEN = 100; Table S1 rxn 9 PTEN = 100
    pp2a <- fixed(20); label("PP2A phosphatase concentration (nM)") # VEGFRModelDrugs.m PP2A = 20; Table S1 rxn 12 PP2A = 20
    nf1 <- fixed(25); label("NF1 GTPase-activating protein concentration (nM)") # VEGFRModelDrugs.m NF1 = 25; Table S1 rxn 29 NF1 = 25
    f_cell <- fixed(0.5); label("Fraction of tumor volume occupied by cells (unitless)") # VEGFRModelDrugs.m f = 0.5 (Vc = f x V_tumor, Ve = (1 - f) x V_tumor); Table S1 note Vc/Ve = 1

    # ---- Sunitinib PK (Table S4; VEGFRModelDrugs.m; Houk 2009) ----
    lka <- fixed(log(0.0033)); label("Sunitinib absorption rate constant (1/min)") # VEGFRModelDrugs.m ka = 0.0033; Table S4 Ka = 0.0033
    lcl <- fixed(log(0.863)); label("Sunitinib clearance (L/min)") # VEGFRModelDrugs.m cls = 0.863; Table S4 CL = 0.86
    lvc <- fixed(log(2030)); label("Sunitinib central volume (L)") # VEGFRModelDrugs.m vcs = 2030; Table S4 Vc = 2030
    lq <- fixed(log(0.120)); label("Sunitinib intercompartmental clearance (L/min)") # VEGFRModelDrugs.m clds = 0.120; Table S4 CLd = 0.12
    lvp <- fixed(log(583)); label("Sunitinib peripheral volume (L)") # VEGFRModelDrugs.m vps = 583
    lfdepot <- fixed(log(0.95)); label("Sunitinib oral bioavailability (fraction)") # VEGFRModelDrugs.m F = 0.95; Table S4 F = 0.95
    lkon <- fixed(log(1)); label("Sunitinib-receptor association rate constant kton (1/(nM*min))") # VEGFRModelDrugs.m kton = 1; Table S3 sunitinib kton = 1
    lkoff <- fixed(log(10)); label("Sunitinib-receptor dissociation rate constant ktoff (1/min)") # VEGFRModelDrugs.m ktoff = 10; Table S3 sunitinib ktoff = 10

    # ---- Bevacizumab PK and VEGF binding (Table S4, Table S3; Li 2013) ----
    lcl_bevacizumab <- fixed(log(0.00016)); label("Bevacizumab clearance (L/min)") # VEGFRModelDrugs.m cl = 0.000160; Table S4 CL = 0.00016 (Bev2comp.m used 0.000150 for Figure S1 only)
    lvc_bevacizumab <- fixed(log(2.83)); label("Bevacizumab central volume (L)") # VEGFRModelDrugs.m vc = 2.83; Table S4 Vc = 2.83
    lq_bevacizumab <- fixed(log(0.0004416)); label("Bevacizumab intercompartmental clearance (L/min)") # VEGFRModelDrugs.m cld = 0.0004416; Table S4 CLd = 0.00044
    lvp_bevacizumab <- fixed(log(2.49)); label("Bevacizumab peripheral volume (L)") # VEGFRModelDrugs.m vp = 2.49; Table S4 Vp = 2.49
    lkon_bevacizumab <- fixed(log(0.0021)); label("Bevacizumab-VEGF association rate constant kton3 (1/(nM*min))") # VEGFRModelDrugs.m kton3 = 0.0021; Table S3 bevacizumab kton = 0.0021
    lkoff_bevacizumab <- fixed(log(0.00378)); label("Bevacizumab-VEGF dissociation rate constant ktoff3 (1/min)") # VEGFRModelDrugs.m ktoff3 = 0.00378; Table S3 bevacizumab ktoff = 0.00378

    # ---- Buparlisib (BKM120) PK and PI3K inhibition (Table S4, Table S3; fitted by the authors to Bendell 2012) ----
    lka_buparlisib <- fixed(log(0.0136)); label("Buparlisib absorption rate constant (1/min)") # VEGFRModelDrugs.m ka_b = 0.0136; Table S4 Ka = 0.0136
    lcl_buparlisib <- fixed(log(0.0144)); label("Buparlisib clearance (L/min)") # VEGFRModelDrugs.m cls_b = 0.0144; Table S4 CL = 0.014
    lvc_buparlisib <- fixed(log(43.1)); label("Buparlisib central volume (L)") # VEGFRModelDrugs.m vcs_b = 43.1; Table S4 Vc = 43.05
    lq_buparlisib <- fixed(log(0.711)); label("Buparlisib intercompartmental clearance (L/min)") # VEGFRModelDrugs.m clds_b = 0.711; Table S4 CLd = 0.71
    lvp_buparlisib <- fixed(log(899.95)); label("Buparlisib peripheral volume (L)") # VEGFRModelDrugs.m vps_b = 899.95; Table S4 Vp = 899.95
    lki_buparlisib <- fixed(log(30)); label("Buparlisib competitive inhibition constant on PI3K Ki (nM)") # VEGFRModelDrugs.m ki = 30; Table S3 BKM-120 ki = 30

    # ---- Selumetinib (AZD6244) PK and MEK inhibition (Table S4, Table S3; fitted by the authors to O'Neil 2011) ----
    lka_selumetinib <- fixed(log(0.022)); label("Selumetinib absorption rate constant (1/min)") # VEGFRModelDrugs.m ka_m = 0.022; Table S4 Ka = 0.022
    lcl_selumetinib <- fixed(log(0.4084)); label("Selumetinib clearance (L/min)") # VEGFRModelDrugs.m cls_m = 0.4084; Table S4 CL = 0.40
    lvc_selumetinib <- fixed(log(136.13)); label("Selumetinib central volume (L)") # VEGFRModelDrugs.m vcs_m = 136.13; Table S4 Vc = 136.13
    lq_selumetinib <- fixed(log(0.1134)); label("Selumetinib intercompartmental clearance (L/min)") # VEGFRModelDrugs.m clds_m = 0.1134; Table S4 CLd = 0.11
    lvp_selumetinib <- fixed(log(30.9)); label("Selumetinib peripheral volume (L)") # VEGFRModelDrugs.m vps_m = 30.9; Table S4 prints Vp = 32.09 (code value used)
    lki_selumetinib <- fixed(log(14)); label("Selumetinib non-competitive inhibition constant on ERK phosphorylation Ki (nM)") # VEGFRModelDrugs.m kim = 14; Table S3 AZD6244 kim = 14

    # ---- Hypothetical PLC-gamma inhibitor (one-compartment, concentration state) ----
    lkel_plcginhibitor <- fixed(log(0.000963)); label("PLC-gamma inhibitor first-order elimination rate constant (1/min)") # VEGFRModelDrugs.m kplc = 0.000963
    lki_plcginhibitor <- fixed(log(1000)); label("PLC-gamma inhibitor competitive inhibition constant Ki (nM)") # VEGFRModelDrugs.m ki2 = 1000; FigS1A_PLCgPK.m 'Ki (1000 nM)'
  })

  model({
    # Individual parameters (typical values; the model has no random effects)
    k1f <- exp(lk1f)
    k1r <- exp(lk1r)
    k2f <- exp(lk2f)
    k2r <- exp(lk2r)
    k3f <- exp(lk3f)
    k3r <- exp(lk3r)
    k4f <- exp(lk4f)
    k4r <- exp(lk4r)
    kpg <- exp(lkpg)
    kmpg <- exp(lkmpg)
    k6f <- exp(lk6f)
    k6r <- exp(lk6r)
    k7 <- exp(lk7)
    km7 <- exp(lkm7)
    k8 <- exp(lk8)
    km8 <- exp(lkm8)
    k9 <- exp(lk9)
    km9 <- exp(lkm9)
    k10f <- exp(lk10f)
    k10r <- exp(lk10r)
    k11f <- exp(lk11f)
    k11r <- exp(lk11r)
    k12 <- exp(lk12)
    km12 <- exp(lkm12)
    k13 <- exp(lk13)
    km13 <- exp(lkm13)
    k14f <- exp(lk14f)
    k14r <- exp(lk14r)
    k15f <- exp(lk15f)
    k15r <- exp(lk15r)
    k16f <- exp(lk16f)
    k16r <- exp(lk16r)
    k17 <- exp(lk17)
    km17 <- exp(lkm17)
    vm18 <- exp(lvm18)
    km18 <- exp(lkm18)
    k19 <- exp(lk19)
    km19 <- exp(lkm19)
    vm20 <- exp(lvm20)
    km20 <- exp(lkm20)
    k21 <- exp(lk21)
    km21 <- exp(lkm21)
    vm22 <- exp(lvm22)
    km22 <- exp(lkm22)
    k23 <- exp(lk23)
    km23 <- exp(lkm23)
    vm24 <- exp(lvm24)
    km24 <- exp(lkm24)
    k25 <- exp(lk25)
    km25 <- exp(lkm25)
    k26 <- exp(lk26)
    k27f <- exp(lk27f)
    k27r <- exp(lk27r)
    km28 <- exp(lkm28)
    k28 <- exp(lk28)
    k29 <- exp(lk29)
    km29 <- exp(lkm29)
    k32 <- exp(lk32)
    km32 <- exp(lkm32)
    k30 <- exp(lk30)
    km30 <- exp(lkm30)
    vm31 <- exp(lvm31)
    km31 <- exp(lkm31)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)
    kdeg2 <- exp(lkdeg2)
    klsyn <- exp(lklsyn)
    kpip2syn <- exp(lkpip2syn)
    kpip2deg <- exp(lkpip2deg)
    kdegdag <- exp(lkdegdag)
    kdegip3 <- exp(lkdegip3)
    kh <- exp(lkh)

    # Cell-to-extracellular volume ratio in the tumor (Vc / Ve)
    vc_ve <- f_cell / (1 - f_cell)

    # PK parameters
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    kon <- exp(lkon)
    koff <- exp(lkoff)

    cl_bevacizumab <- exp(lcl_bevacizumab)
    vc_bevacizumab <- exp(lvc_bevacizumab)
    q_bevacizumab <- exp(lq_bevacizumab)
    vp_bevacizumab <- exp(lvp_bevacizumab)
    kon_bevacizumab <- exp(lkon_bevacizumab)
    koff_bevacizumab <- exp(lkoff_bevacizumab)

    ka_buparlisib <- exp(lka_buparlisib)
    cl_buparlisib <- exp(lcl_buparlisib)
    vc_buparlisib <- exp(lvc_buparlisib)
    q_buparlisib <- exp(lq_buparlisib)
    vp_buparlisib <- exp(lvp_buparlisib)
    ki_buparlisib <- exp(lki_buparlisib)

    ka_selumetinib <- exp(lka_selumetinib)
    cl_selumetinib <- exp(lcl_selumetinib)
    vc_selumetinib <- exp(lvc_selumetinib)
    q_selumetinib <- exp(lq_selumetinib)
    vp_selumetinib <- exp(lvp_selumetinib)
    ki_selumetinib <- exp(lki_selumetinib)

    kel_plcginhibitor <- exp(lkel_plcginhibitor)
    ki_plcginhibitor <- exp(lki_plcginhibitor)

    # Plasma drug concentrations (nM = nmol/L). Plasma drug drives binding /
    # inhibition in the tumor without being depleted by it (VEGFRModelDrugs.m).
    Cc <- central / vc
    Cc_bevacizumab <- central_bevacizumab / vc_bevacizumab
    Cc_buparlisib <- central_buparlisib / vc_buparlisib
    Cc_selumetinib <- central_selumetinib / vc_selumetinib
    Cc_plcginhibitor <- central_plcginhibitor

    # Reaction rates (nM/min; Table S1 reaction numbers; VEGFRModelDrugs.m)
    v1f <- k1f * vr * vc_ve * vegf
    v1r <- k1r * vegf_vr
    v2f <- k2f * vegf_vr * vegf_vr
    v2r <- k2r * vegf_vr2
    v3f <- k3f * vegf_vr2
    v3r <- k3r * vegf_vr2_p
    v4f <- k4f * vegf_vr2_p * plcg
    v4r <- k4r * vegf_vr2_p_plcg
    # PLC-gamma inhibitor: competitive inhibition of PIP2 hydrolysis
    vpg <- kpg * vegf_vr2_p_plcg * pip2 / (pip2 + kmpg * (1 + Cc_plcginhibitor / ki_plcginhibitor))
    v6f <- k6f * vegf_vr2_p * pi3k
    v6r <- k6r * vegf_vr2_p_pi3k
    v7 <- k7 * shp1 * pi3k_p / (pi3k_p + km7)
    # Buparlisib: competitive inhibition of PI3K-catalysed PIP2 phosphorylation
    v8 <- k8 * vegf_vr2_p_pi3k * pip2 / (km8 * (1 + Cc_buparlisib / ki_buparlisib) + pip2)
    v9 <- k9 * pten * pip3 / (km9 + pip3)
    v10f <- k10f * pip3 * pdk
    v10r <- k10r * pip3_pdk
    v11f <- k11f * pip3 * akt
    v11r <- k11r * pip3_akt
    v12 <- k12 * pip3_pdk * pip3_akt / (km12 + pip3_akt)
    v13 <- k13 * pp2a * pip3_akt_p / (km13 + pip3_akt_p)
    v14f <- k14f * dag * pkc
    v14r <- k14r * dag_pkc
    v15f <- k15f * rasgrp * dag
    v15r <- k15r * dag_rasgrp
    v16f <- k16f * raf * rkip
    v16r <- k16r * raf_rkip
    v17 <- k17 * dag_pkc * rkip / (km17 + rkip)
    v18 <- vm18 * rkip_p / (km18 + rkip_p)
    v19 <- k19 * ras_gtp * raf / (km19 + raf)
    v20 <- vm20 * raf_p / (km20 + raf_p)
    v21 <- k21 * raf_p * mek / (km21 + mek)
    v22 <- vm22 * mek_p / (km22 + mek_p)
    # Selumetinib: non-competitive inhibition of ERK phosphorylation
    v23 <- k23 * mek_p * erk / (1 + Cc_selumetinib / ki_selumetinib) / (km23 + erk)
    v24 <- vm24 * erk_p / (km24 + erk_p)
    v25 <- k25 * erk_p * rsk / (km25 + rsk)
    v26 <- k26 * rsk_p
    v27f <- k27f * vegf_vr2_p * grb2sos
    v27r <- k27r * vegf_vr2_p_grb2sos
    v28 <- k28 * (vegf_vr2_p_grb2sos + dag_rasgrp) * ras_gdp / (km28 + ras_gdp)
    v29 <- k29 * nf1 * ras_gtp / (km29 + ras_gtp)
    v30 <- k30 * rsk_p * raf / (km30 + raf)
    v31 <- vm31 * raf_i / (km31 + raf_i)
    v32 <- k32 * pip3_akt_p * pi3k / (pi3k + km32)
    vdip3 <- kdegip3 * ip3
    vddag <- kdegdag * dag
    vdr <- kdeg * vr
    vdrl <- kdeg * vegf_vr
    vdrl2 <- kdeg * vegf_vr2
    vdrl2p <- kdeg2 * vegf_vr2_p
    # Receptor synthesis under pAkt negative feedback (Table S1 rxn 36)
    v36 <- ksyn * (1 - pip3_akt_p^n_fb / (kh^n_fb + pip3_akt_p^n_fb))
    vdegpip2 <- kpip2deg * pip2
    # Sunitinib binding to free, ligand-bound and dimerised receptor (Table S3)
    v69 <- kon * Cc * vr
    v70 <- koff * sunitinib_vr
    v73 <- kon * Cc * vegf_vr
    v74 <- koff * sunitinib_vegf_vr
    v75 <- kon * Cc * vegf_vr2
    v76 <- koff * sunitinib_vegf_vr2
    # Bevacizumab sequestration of VEGF (Table S3)
    v71 <- kon_bevacizumab * Cc_bevacizumab * vegf
    v72 <- koff_bevacizumab * bevacizumab_vegf

    # Signaling network ODEs (Table S2; VEGFRModelDrugs.m)
    d/dt(vr) <- v1r - v1f / vc_ve + v36 - vdr - v69 + v70
    d/dt(vegf) <- v1r * vc_ve - v1f + klsyn - v71 + v72
    d/dt(vegf_vr) <- v1f / vc_ve - v1r - 2 * v2f + 2 * v2r - vdrl - v73 + v74
    d/dt(vegf_vr2) <- v2f - v2r - v3f + v3r - vdrl2 - v75 + v76
    d/dt(vegf_vr2_p) <- v3f - v3r - vdrl2p - v6f + v6r - v4f + v4r - v27f + v27r
    d/dt(plcg) <- -v4f + v4r
    d/dt(pi3k_p) <- v32 - v7
    d/dt(pi3k) <- v7 - v32 + v6r - v6f
    d/dt(pip3) <- v8 - v9 - v10f + v10r - v11f + v11r
    d/dt(pip2) <- v9 - v8 + kpip2syn - vdegpip2 - vpg
    d/dt(pdk) <- v10r - v10f
    d/dt(akt) <- v11r - v11f
    d/dt(pip3_pdk) <- v10f - v10r
    d/dt(pip3_akt) <- v11f - v11r + v13 - v12
    d/dt(pip3_akt_p) <- v12 - v13
    d/dt(pkc) <- -v14f + v14r
    d/dt(dag_pkc) <- v14f - v14r
    d/dt(raf) <- -v16f + v16r - v19 + v20 - v30 + v31
    d/dt(rkip) <- -v16f + v16r - v17 + v18
    d/dt(raf_rkip) <- v16f - v16r
    d/dt(rkip_p) <- v17 - v18
    d/dt(raf_p) <- v19 - v20
    d/dt(raf_i) <- v30 - v31
    d/dt(rsk_p) <- v25 - v26
    d/dt(rsk) <- -v25 + v26
    d/dt(erk_p) <- v23 - v24
    d/dt(erk) <- -v23 + v24
    d/dt(grb2sos) <- -v27f + v27r
    d/dt(ras_gtp) <- v28 - v29
    d/dt(ras_gdp) <- v29 - v28
    d/dt(mek_p) <- v21 - v22
    d/dt(mek) <- -v21 + v22
    d/dt(vegf_vr2_p_pi3k) <- v6f - v6r
    d/dt(vegf_vr2_p_plcg) <- v4f - v4r
    d/dt(vegf_vr2_p_grb2sos) <- v27f - v27r
    d/dt(ip3) <- vpg - vdip3
    d/dt(dag) <- vpg - vddag - v15f + v15r - v14f + v14r
    d/dt(rasgrp) <- -v15f + v15r
    d/dt(dag_rasgrp) <- v15f - v15r
    d/dt(sunitinib_vr) <- v69 - v70 - kdeg * sunitinib_vr
    d/dt(sunitinib_vegf_vr) <- v73 - v74 - kdeg * sunitinib_vegf_vr
    d/dt(sunitinib_vegf_vr2) <- v75 - v76 - kdeg * sunitinib_vegf_vr2
    d/dt(bevacizumab_vegf) <- v71 - v72

    # Drug PK (amounts in nmol; Table S4)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(central_bevacizumab) <- -cl_bevacizumab / vc_bevacizumab * central_bevacizumab -
      q_bevacizumab / vc_bevacizumab * central_bevacizumab +
      q_bevacizumab / vp_bevacizumab * peripheral1_bevacizumab
    d/dt(peripheral1_bevacizumab) <- q_bevacizumab / vc_bevacizumab * central_bevacizumab -
      q_bevacizumab / vp_bevacizumab * peripheral1_bevacizumab
    d/dt(depot_buparlisib) <- -ka_buparlisib * depot_buparlisib
    d/dt(central_buparlisib) <- ka_buparlisib * depot_buparlisib -
      cl_buparlisib / vc_buparlisib * central_buparlisib -
      q_buparlisib / vc_buparlisib * central_buparlisib +
      q_buparlisib / vp_buparlisib * peripheral1_buparlisib
    d/dt(peripheral1_buparlisib) <- q_buparlisib / vc_buparlisib * central_buparlisib -
      q_buparlisib / vp_buparlisib * peripheral1_buparlisib
    d/dt(depot_selumetinib) <- -ka_selumetinib * depot_selumetinib
    d/dt(central_selumetinib) <- ka_selumetinib * depot_selumetinib -
      cl_selumetinib / vc_selumetinib * central_selumetinib -
      q_selumetinib / vc_selumetinib * central_selumetinib +
      q_selumetinib / vp_selumetinib * peripheral1_selumetinib
    d/dt(peripheral1_selumetinib) <- q_selumetinib / vc_selumetinib * central_selumetinib -
      q_selumetinib / vp_selumetinib * peripheral1_selumetinib
    d/dt(central_plcginhibitor) <- -kel_plcginhibitor * central_plcginhibitor

    f(depot) <- exp(lfdepot)

    # Pre-equilibration initial conditions (nM; GetInitialConditions.m, Table S2).
    # Every state not listed starts at 0.
    vr(0) <- 200
    vegf(0) <- 0.5
    plcg(0) <- 250
    pi3k(0) <- 200
    pip2(0) <- 1000
    pdk(0) <- 100.534
    akt(0) <- 100
    pkc(0) <- 200
    raf(0) <- 20
    rkip(0) <- 20
    rsk(0) <- 100
    erk(0) <- 400
    grb2sos(0) <- 82.3
    ras_gdp(0) <- 95.7
    mek(0) <- 200
    rasgrp(0) <- 82.3

    # Outputs: active ERK and active Akt (nM)
    pERK <- erk_p
    pAkt <- pip3_akt_p
  })
}
