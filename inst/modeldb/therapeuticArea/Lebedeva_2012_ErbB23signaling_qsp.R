Lebedeva_2012_ErbB23signaling_qsp <- function() {
  description <- paste(
    "QSP. In vitro (PE04 ovarian carcinoma cell line). Kinetic model of",
    "heregulin-beta (HRG) induced ErbB2/ErbB3 signalling with",
    "anti-ErbB2 antibody (pertuzumab, 2C4) binding, used for global and",
    "local sensitivity analysis of Akt phosphorylation (Lebedeva 2012;",
    "a modification of the Faratian 2009 model). 53 species from the",
    "authors' SBML file plus one bookkeeping integral: HRG binding to",
    "ErbB3 and two-step ErbB2/ErbB3 heterodimerisation and",
    "transphosphorylation; Shc / Grb2-Sos / Ras / Raf / MEK / ERK",
    "(MAPK) cascade; PI3K activation, PIP2 to PIP3 conversion, the PTEN",
    "catalytic and phosphorylation cycles and two-step Akt",
    "phosphorylation by PDK1; PP2A as the shared phosphatase of",
    "phospho-MEK and phospho-Akt (Akt/MAPK crosstalk), with ppAkt",
    "inactivating Raf; and two-step pertuzumab-ErbB2 binding. Mass-action",
    "and Michaelis-Menten rate laws, time in seconds, concentrations in nM",
    "referenced to the cytoplasmic volume (HRG 1 nM and pertuzumab 100 nM",
    "in the medium correspond to model concentrations of 3000 and",
    "300000). Nominal parameter values are the paper's best calibration",
    "fit to PE04 time courses of pErbB3, pAkt and pErk; the model is",
    "deterministic (no between-cell variability, no residual error).",
    "Outputs: pAkt (phosphorylated fraction of total Akt), its running",
    "integral auc_pakt (the paper's sensitivity characteristic S_pAkt),",
    "and total phosphorylated receptor, MEK and ERK.",
    sep = " "
  )
  reference <- paste(
    "Lebedeva G, Sorokin A, Faratian D, Mullen P, Goltsov A, Langdon SP,",
    "Harrison DJ, Goryanin I (2012). Model-based global sensitivity",
    "analysis as applied to identification of anti-cancer drug targets",
    "and biomarkers of drug resistance in the ErbB2/3 network.",
    "Eur J Pharm Sci 46(4):244-258. doi:10.1016/j.ejps.2011.10.026.",
    "Rate laws from the authors' SBML file (Additional File 4) and",
    "Additional File 1 (system of ODEs); nominal parameter values from",
    "Supplementary Table S2 (Additional File 2). The model-unit",
    "pertuzumab concentration follows the curated BioModels entry",
    "BIOMD0000000424 of the predecessor model (Faratian D et al. 2009,",
    "Cancer Res 69(16):6713-6720, doi:10.1158/0008-5472.CAN-09-0777).",
    sep = " "
  )
  vignette <- "Lebedeva_2012_ErbB23signaling_qsp"
  units <- list(
    time = "s",
    dosing = "nM",
    concentration = "nM"
  )

  # HRG and pertuzumab are added to the culture as concentration steps
  # (bolus doses into unit-volume states).
  dosing <- c("hrg", "per")

  # Every state is a named molecular species of the signalling network (or
  # the integral of pAkt); none has a canonical analogue. Names are the
  # lower-cased SBML species ids, except that the SBML pair E2_Per / E2Per
  # (loose and tight pertuzumab-ErbB2 complexes, 'E2-2C4' and 'E2-2C4c' in
  # Additional File 1) becomes e2_per / e2_per_c, matching e3h_c and e23h_c.
  paper_specific_compartments <- c(
    "hrg",
    "per",
    "e3",
    "e3h",
    "e3h_c",
    "e2",
    "e23h",
    "e23h_c",
    "e23hp",
    "e2_per",
    "e2_per_c",
    "shc",
    "e23hp_shc",
    "e23hp_shcp",
    "gs",
    "e23hp_shgs",
    "shgs",
    "shcp",
    "rasgdp",
    "rasgtp",
    "raf",
    "rafa",
    "mek",
    "mekp",
    "mekpp",
    "pp2a",
    "mekp_pp2a",
    "mek_pp2a",
    "mekpp_pp2a",
    "erk",
    "erkp",
    "erkpp",
    "pi3k",
    "e23hp_pi3k",
    "e23hp_pi3ka",
    "pi3ka",
    "pi2",
    "pi3ka_pi",
    "pi3ka_pip3",
    "pip3",
    "pten",
    "ptenp",
    "pten_pip3",
    "pten_pi",
    "pten_bpv",
    "pi3k_ly",
    "akt",
    "akt_pip3",
    "akt_pi_p",
    "akt_pi_pp",
    "akt_pi_p_pp2a",
    "akt_pi_pp_pp2a",
    "akt_pip3_pp2a",
    "auc_pakt"
  )

  compartmentData <- list(
    hrg = list(analyte = "heregulin-beta (free ligand)", units = "nM", specimen = "not applicable", verified = TRUE),
    per = list(analyte = "pertuzumab (free antibody)", units = "nM", specimen = "not applicable", verified = TRUE),
    e3 = list(analyte = "ErbB3 receptor (free)", units = "nM", specimen = "tumor", verified = TRUE),
    e3h = list(analyte = "ErbB3-HRG intermediate complex", units = "nM", specimen = "tumor", verified = TRUE),
    e3h_c = list(analyte = "ErbB3-HRG final complex", units = "nM", specimen = "tumor", verified = TRUE),
    e2 = list(analyte = "ErbB2 receptor (free)", units = "nM", specimen = "tumor", verified = TRUE),
    e23h = list(
      analyte = "ErbB2-ErbB3-HRG heterodimer, intermediate form",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    e23h_c = list(
      analyte = "ErbB2-ErbB3-HRG heterodimer, final form",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    e23hp = list(
      analyte = "phosphorylated ErbB2-ErbB3-HRG heterodimer (pE23H)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    e2_per = list(
      analyte = "ErbB2-pertuzumab primary (loose) complex",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    e2_per_c = list(analyte = "ErbB2-pertuzumab tight complex", units = "nM", specimen = "tumor", verified = TRUE),
    shc = list(analyte = "Shc", units = "nM", specimen = "tumor", verified = TRUE),
    e23hp_shc = list(analyte = "pE23H-Shc complex", units = "nM", specimen = "tumor", verified = TRUE),
    e23hp_shcp = list(analyte = "pE23H-pShc complex", units = "nM", specimen = "tumor", verified = TRUE),
    gs = list(analyte = "Grb2-Sos complex (GS)", units = "nM", specimen = "tumor", verified = TRUE),
    e23hp_shgs = list(analyte = "pE23H-pShc-GS complex", units = "nM", specimen = "tumor", verified = TRUE),
    shgs = list(analyte = "pShc-GS complex", units = "nM", specimen = "tumor", verified = TRUE),
    shcp = list(analyte = "phosphorylated Shc", units = "nM", specimen = "tumor", verified = TRUE),
    rasgdp = list(analyte = "Ras-GDP", units = "nM", specimen = "tumor", verified = TRUE),
    rasgtp = list(analyte = "Ras-GTP", units = "nM", specimen = "tumor", verified = TRUE),
    raf = list(analyte = "Raf (inactive)", units = "nM", specimen = "tumor", verified = TRUE),
    rafa = list(analyte = "activated Raf (Raf*)", units = "nM", specimen = "tumor", verified = TRUE),
    mek = list(analyte = "MEK", units = "nM", specimen = "tumor", verified = TRUE),
    mekp = list(analyte = "phosphorylated MEK (pMEK)", units = "nM", specimen = "tumor", verified = TRUE),
    mekpp = list(analyte = "doubly phosphorylated MEK (ppMEK)", units = "nM", specimen = "tumor", verified = TRUE),
    pp2a = list(analyte = "protein phosphatase 2A (free)", units = "nM", specimen = "tumor", verified = TRUE),
    mekp_pp2a = list(analyte = "pMEK-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    mek_pp2a = list(analyte = "MEK-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    mekpp_pp2a = list(analyte = "ppMEK-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    erk = list(analyte = "ERK", units = "nM", specimen = "tumor", verified = TRUE),
    erkp = list(analyte = "phosphorylated ERK (pERK)", units = "nM", specimen = "tumor", verified = TRUE),
    erkpp = list(analyte = "doubly phosphorylated ERK (ppERK)", units = "nM", specimen = "tumor", verified = TRUE),
    pi3k = list(analyte = "PI3K (inactive, free)", units = "nM", specimen = "tumor", verified = TRUE),
    e23hp_pi3k = list(analyte = "pE23H-PI3K complex", units = "nM", specimen = "tumor", verified = TRUE),
    e23hp_pi3ka = list(analyte = "pE23H-PI3K* complex", units = "nM", specimen = "tumor", verified = TRUE),
    pi3ka = list(analyte = "activated PI3K (PI3K*)", units = "nM", specimen = "tumor", verified = TRUE),
    pi2 = list(
      analyte = "phosphatidylinositol-4,5-bisphosphate (PIP2, 'PI')",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    pi3ka_pi = list(analyte = "PI3K*-PIP2 complex", units = "nM", specimen = "tumor", verified = TRUE),
    pi3ka_pip3 = list(analyte = "PI3K*-PIP3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    pip3 = list(
      analyte = "phosphatidylinositol-3,4,5-trisphosphate (PIP3)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    pten = list(analyte = "PTEN (active, unphosphorylated, free)", units = "nM", specimen = "tumor", verified = TRUE),
    ptenp = list(analyte = "phosphorylated PTEN (pPTEN)", units = "nM", specimen = "tumor", verified = TRUE),
    pten_pip3 = list(analyte = "PTEN-PIP3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    pten_pi = list(analyte = "PTEN-PIP2 complex", units = "nM", specimen = "tumor", verified = TRUE),
    pten_bpv = list(analyte = "PTEN-bpV(pic) inhibitor complex", units = "nM", specimen = "tumor", verified = TRUE),
    pi3k_ly = list(analyte = "PI3K-LY294002 inhibitor complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt = list(analyte = "Akt (free)", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pip3 = list(analyte = "Akt-PIP3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pi_p = list(analyte = "pAkt-PIP3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pi_pp = list(analyte = "ppAkt-PIP3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pi_p_pp2a = list(analyte = "pAkt-PIP3-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pi_pp_pp2a = list(analyte = "ppAkt-PIP3-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    akt_pip3_pp2a = list(analyte = "Akt-PIP3-PP2A complex", units = "nM", specimen = "tumor", verified = TRUE),
    auc_pakt = list(
      analyte = "running integral of pAkt (phosphorylated fraction of total Akt)",
      units = "s",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (PE04 ovarian carcinoma cell line)",
    disease_state = "Human ovarian carcinoma cell line PE04 (pertuzumab-responsive)",
    dose_range = "Heregulin-beta 1 nM, alone or 15 min after pertuzumab 100 nM; time courses over 60 min",
    notes = paste(
      "Calibrated (Faratian 2009 data, recalibrated in Lebedeva 2012) to",
      "western-blot time courses of pErbB3, pAkt and pErk in PE04 cells",
      "stimulated with heregulin-beta in the absence and presence of",
      "pertuzumab (Additional File 1 Fig. S6). The model was not fully",
      "identifiable; Table S2 holds the parameters of one of the best fits.",
      "Lebedeva 2012 also tested model predictions experimentally in the",
      "OVCAR4 cell line (pertuzumab-resistant) without refitting."
    )
  )

  ini({
    # All values: Supplementary Table S2 (Additional File 2), column
    # 'Nominal parameter value' - the parameters 'identified in one of the
    # best fittings'. Species concentrations, Kd, Ke and Km in nM; rate
    # constants in 1/s or 1/(nM*s). Names follow Table S2 (and the labels of
    # Figs. 3 and 6), except that each total protein concentration carries a
    # '_tot' suffix so it does not collide with the free-species state, and
    # that k13 and k31 enter on the log scale (lk13, lk31): those two names
    # coincide with the canonical two-compartment PK micro-constants, which
    # the package convention log-transforms. Their meaning here is the
    # paper's (Raf activation; PI3K*-PIP2 binding).

    # ---- Total protein concentrations (initial free-species amounts) ----
    E3_tot <- 80; label("Total ErbB3 (nM)") # Table S2 'E3' = 80
    E2_tot <- 100; label("Total ErbB2 (nM)") # Table S2 'E2' = 100
    Shc_tot <- 100; label("Total Shc (nM)") # Table S2 'Shc' = 100
    GS_tot <- 100; label("Total Grb2-Sos (nM)") # Table S2 'GS' = 100
    RasGDP_tot <- 120; label("Total Ras (initially Ras-GDP) (nM)") # Table S2 'RasGDP' = 120
    Raf_tot <- 50; label("Total Raf (nM)") # Table S2 'Raf' = 50
    MEK_tot <- 10; label("Total MEK (nM)") # Table S2 'MEK' = 10
    PP2A_tot <- 15; label("Total PP2A (nM)") # Table S2 'PP2A' = 15
    ERK_tot <- 10; label("Total ERK (nM)") # Table S2 'ERK' = 10
    PI3K_tot <- 200; label("Total PI3K (nM)") # Table S2 'PI3K' = 200
    PI_tot <- 300; label("Total phosphoinositol, PIP2 + PIP3 (nM)") # Table S2 'PI' = 300
    Akt_tot <- 50; label("Total Akt (nM)") # Table S2 'Akt' = 50
    PTEN_tot <- 30; label("Total PTEN, PTEN + pPTEN at rest (nM)") # Table S2 'PTEN' = 30

    # ---- Receptor module ----
    k1 <- 0.005; label("ErbB3-heregulin binding rate constant (1/(nM*s))") # Table S2 'k1' = 0.005
    Kd_1 <- 600; label("ErbB3-heregulin dissociation constant (nM)") # Table S2 'Kd_1' = 600
    k2 <- 10; label("ErbB2/3 receptor dimerisation rate constant (1/(nM*s))") # Table S2 'k2' = 10
    Kd_2 <- 10; label("ErbB2/3 receptor dimer dissociation constant (nM)") # Table S2 'Kd_2' = 10
    k3 <- 1; label("ErbB2/3 transphosphorylation rate constant (1/s)") # Table S2 'k3' = 1
    Ke_3 <- 0.1; label("ErbB2/3 transphosphorylation equilibrium constant (unitless)") # Table S2 'Ke_3' = 0.1 (SBML id Kd_3)
    V4 <- 10; label("Maximal rate of ErbB2/3 dimer dephosphorylation (nM/s)") # Table S2 'V4' = 10
    K4 <- 50; label("Michaelis constant of ErbB2/3 dimer dephosphorylation (nM)") # Table S2 'K4' = 50
    k48 <- 0.001; label("Degradation rate constant of the phosphorylated dimer pE23H (1/s)") # Table S2 'k48' = 0.001
    k51 <- 0.01; label("ErbB3-HRG intermediate to final complex conversion rate constant (1/s)") # Table S2 'k51' = 0.01
    k53 <- 0.01; label("ErbB2/3 dimer intermediate to final form conversion rate constant (1/s)") # Table S2 'k53' = 0.01

    # ---- Shc / Grb2-Sos / Ras / Raf ----
    k5 <- 0.06; label("pE23H-Shc binding rate constant (1/(nM*s))") # Table S2 'k5' = 0.06
    Kd_5 <- 1; label("pE23H-Shc dissociation constant (nM)") # Table S2 'Kd_5' = 1
    k6 <- 12; label("Shc phosphorylation rate constant on the receptor (1/s)") # Table S2 'k6' = 12
    k_6 <- 3; label("Reverse rate constant of Shc phosphorylation on the receptor (1/s)") # Table S2 'k_6' = 3
    k7 <- 36; label("pE23H-pShc binding with Grb2-Sos rate constant (1/(nM*s))") # Table S2 'k7' = 36
    Kd_7 <- 9; label("pE23H-pShc-GS dissociation constant (nM)") # Table S2 'Kd_7' = 9
    k8 <- 12; label("Release rate constant of pShc-GS from the receptor (1/s)") # Table S2 'k8' = 12
    Kd_8 <- 0.1; label("Dissociation constant of pShc-GS from the receptor (nM)") # Table S2 'Kd_8' = 0.1
    k9 <- 35; label("Dissociation rate constant of the pShc-GS complex (1/s)") # Table S2 'k9' = 35
    k_9 <- fixed(0); label("Reverse (re-association) constant of pShc-GS, as written in the SBML rate law (1/nM)") # SBML parameter 'k_9' = 0; absent from Table S2 (V9 is irreversible in the calibrated model)
    V10 <- 0.0154; label("Maximal rate of pShc dephosphorylation (nM/s)") # Table S2 'V10' = 0.0154
    K10 <- 340; label("Michaelis constant of pShc dephosphorylation (nM)") # Table S2 'K10' = 340
    k11 <- 6; label("Ras activation catalytic rate constant (1/s)") # Table S2 'k11' = 6
    K11 <- 0.18; label("Michaelis constant of Ras activation (nM)") # Table S2 'K11' = 0.18
    V12 <- 3; label("Maximal rate of Ras deactivation (nM/s)") # Table S2 'V12' = 3
    K12 <- 0.1; label("Michaelis constant of Ras deactivation (nM)") # Table S2 'K12' = 0.1
    lk13 <- log(5); label("Raf activation catalytic rate constant, k13 (log 1/s)") # Table S2 'k13' = 5
    K13 <- 11.7; label("Michaelis constant of Raf activation (nM)") # Table S2 'K13' = 11.7
    k14 <- 0.1; label("Raf* deactivation by ppAkt-PIP3, catalytic rate constant (1/s)") # Table S2 'k14' = 0.1
    K14 <- 50; label("Michaelis constant of Raf* deactivation by ppAkt-PIP3 (nM)") # Table S2 'K14' = 50
    k14_1 <- 0.2; label("Raf* deactivation by the Raf phosphatase E_raf, catalytic rate constant (1/s)") # Table S2 'k14_1' = 0.2
    K14_1 <- 50; label("Michaelis constant of Raf* deactivation by E_raf (nM)") # Table S2 'K14_1' = 50
    E_raf <- 7; label("Concentration of the Raf*-dephosphorylating enzyme E_raf (nM)") # Table S2 'E_raf' = 7

    # ---- MEK / ERK and PP2A ----
    k15 <- 2.1; label("MEK phosphorylation by Raf*, catalytic rate constant (1/s)") # Table S2 'k15' = 2.1
    K15 <- 1; label("Michaelis constant of MEK phosphorylation (nM)") # Table S2 'K15' = 1
    k16 <- 1.5; label("pMEK and ppMEK binding with PP2A rate constant (1/(nM*s))") # Table S2 'k16' = 1.5
    Kd_16 <- 0.1; label("Dissociation constant of the ppMEK-PP2A complex (nM)") # Table S2 'Kd_16' = 0.1
    k16_cat <- 5; label("Catalytic rate constant of pMEK and ppMEK dephosphorylation by PP2A (1/s)") # Table S2 'k16_cat' = 5 (SBML id k16_kat)
    k18 <- 5; label("Release rate constant of MEK from MEK-PP2A (1/s)") # Table S2 'k18' = 5
    k22 <- 5; label("Release rate constant of pMEK from pMEK-PP2A (1/s)") # Table S2 'k22' = 5
    k23 <- 1.2; label("ERK phosphorylation by ppMEK, catalytic rate constant (1/s)") # Table S2 'k23' = 1.2
    K23 <- 10; label("Michaelis constant of ERK phosphorylation (nM)") # Table S2 'K23' = 10
    V24 <- 10; label("Maximal rate of pERK and ppERK dephosphorylation (nM/s)") # Table S2 'V24' = 10
    K24 <- 10; label("Michaelis constant of pERK and ppERK dephosphorylation (nM)") # Table S2 'K24' = 10

    # ---- PI3K / PIP3 / PTEN ----
    k27 <- 3; label("PI3K binding to the ErbB2/3 complex rate constant (1/(nM*s))") # Table S2 'k27' = 3
    Kd_27 <- 1; label("pE23H-PI3K dissociation constant (nM)") # Table S2 'Kd_27' = 1
    k28 <- 300; label("PI3K activation on the receptor rate constant (1/s)") # Table S2 'k28' = 300
    k_28 <- 0.01; label("Reverse ratio of PI3K activation on the receptor (unitless, SBML rate-law form)") # Table S2 'k_28' = 0.01 (SBML file value 0)
    k29 <- 13520; label("PI3K*-ErbB2/3 dissociation rate constant (1/s)") # Table S2 'k29' = 13520
    k_29 <- 0.01; label("PI3K* re-binding to pE23H rate constant (1/(nM*s))") # Table S2 'k_29' = 0.01 (SBML file value 0)
    V30 <- 900; label("PI3K* inactivation rate constant (1/s)") # Table S2 'V30' = 900
    lk31 <- log(0.03); label("PI3K* and PIP2 binding rate constant, k31 (log 1/(nM*s))") # Table S2 'k31' = 0.03
    K_d31 <- 100; label("PI3K*-PIP2 dissociation constant (nM)") # Table S2 'K_d31' = 100
    k55 <- 30; label("PIP2 phosphorylation by PI3K* rate constant (1/s)") # Table S2 'k55' = 30
    k56 <- 30; label("Release rate constant of PIP3 from PI3K*-PIP3 (1/s)") # Table S2 'k56' = 30
    k32 <- 8000; label("PIP3 binding with PTEN rate constant (1/(nM*s))") # Table S2 'k32' = 8000
    Kd_32 <- 0.01; label("PTEN-PIP3 dissociation constant (nM)") # Table S2 'Kd_32' = 0.01
    k33 <- 15; label("PIP3 dephosphorylation by PTEN rate constant (1/s)") # Table S2 'k33' = 15
    k34 <- 3.6; label("PTEN-PIP2 dissociation rate constant (1/s)") # Table S2 'k34' = 3.6
    V35 <- 200; label("Maximal rate of PTEN phosphorylation (nM/s)") # Table S2 'V35' = 200
    K35 <- 2; label("Michaelis constant of PTEN phosphorylation (nM)") # Table S2 'K35' = 2
    V_35 <- 150; label("Maximal rate of pPTEN dephosphorylation (nM/s)") # Table S2 'V_35' = 150
    K_35 <- 2; label("Michaelis constant of pPTEN dephosphorylation (nM)") # Table S2 'K_35' = 2

    # ---- Akt / PDK1 / PP2A ----
    k39 <- 1000; label("Akt-PIP3 binding rate constant (1/(nM*s))") # Table S2 'k39' = 1000
    Kd_39 <- 20; label("Akt-PIP3 dissociation constant (nM)") # Table S2 'Kd_39' = 20
    V40 <- 100; label("Maximal rate of Akt phosphorylation by PDK1 (nM/s)") # Table S2 'V40' = 100
    K40 <- 0.1; label("Michaelis constant of Akt phosphorylation by PDK1 (nM)") # Table S2 'K40' = 0.1
    k41 <- 3; label("PP2A binding with pAkt-PIP3 and ppAkt-PIP3 rate constant (1/(nM*s))") # Table S2 'k41' = 3
    Kd_41 <- 0.1; label("Dissociation constant of the ppAkt-PIP3-PP2A complex (nM)") # Table S2 'Kd_41' = 0.1
    k42 <- 10; label("pAkt and ppAkt dephosphorylation by PP2A, catalytic rate constant (1/s)") # Table S2 'k42' = 10
    k43 <- 10; label("Release rate constant of Akt-PIP3 from Akt-PIP3-PP2A (1/s)") # Table S2 'k43' = 10
    k47 <- 10; label("Release rate constant of pAkt-PIP3 from pAkt-PIP3-PP2A (1/s)") # Table S2 'k47' = 10

    # ---- Pertuzumab (2C4) binding to ErbB2 ----
    k49 <- 0.003; label("ErbB2-pertuzumab primary complex formation rate constant (1/(nM*s))") # Table S2 'k49' = 0.003
    Kd_49 <- 40000; label("ErbB2-pertuzumab primary complex dissociation constant (nM)") # Table S2 'Kd_49' = 40000
    k50 <- 0.1; label("ErbB2-pertuzumab tight complex formation rate constant (1/s)") # Table S2 'k50' = 0.1
    k_50 <- 0.012; label("ErbB2-pertuzumab tight complex reversal rate constant (1/s)") # Table S2 'k_50' = 0.012

    # ---- Inhibitor reactions carried in the SBML file (Faratian 2009) ----
    # Not among the 91 parameters of Table S2; the SBML file sets both
    # inhibitor concentrations to 0, so these reactions are inactive in every
    # Lebedeva 2012 simulation.
    k57 <- fixed(100); label("PTEN binding with bpV(pic) rate constant (1/(nM*s))") # SBML parameter 'k57' = 100
    Kd_57 <- fixed(10); label("PTEN-bpV(pic) dissociation constant (nM)") # SBML parameter 'Kd_57' = 10
    k58 <- fixed(100); label("PI3K binding with LY294002 rate constant (1/(nM*s))") # SBML parameter 'k58' = 100
    Kd_58 <- fixed(80); label("PI3K-LY294002 dissociation constant (nM)") # SBML parameter 'Kd_58' = 80
    bpV <- fixed(0); label("Concentration of the PTEN inhibitor bpV(pic), model units (nM)") # SBML parameter 'bpV' = 0
    LY <- fixed(0); label("Concentration of the PI3K inhibitor LY294002, model units (nM)") # SBML parameter 'LY' = 0
  })

  model({
    # Additional File 1 numbers the reactions V1-V58 (V37 and V38 retired);
    # the SBML file names them R1-R44. Each rate below is the SBML kinetic
    # law (the compartment 'Default' has size 1 and cancels); the comment
    # names both ids. Where the SBML and Additional File 1 differ the SBML
    # form is used (see the vignette).

    k13 <- exp(lk13)
    k31 <- exp(lk31)

    # ---- Initial conditions ----
    # The PTEN phosphorylation cycle (V35/V36) is the only reaction active
    # before stimulation, so free PTEN starts at its Michaelis-Menten
    # steady state: V35*x/(K35 + x) = V_35*(T - x)/(K_35 + T - x), with
    # T = PTEN_tot. The quadratic a*x^2 + b*x + c = 0 has exactly one root
    # in (0, T); the form 2c/(-b - sqrt(b^2 - 4ac)) returns it for either
    # sign of a (and for a = 0). At Table S2 values this gives the SBML
    # initial concentrations PTEN = 4.5644 nM and PTENP = 25.4356 nM.
    a_pten <- -(V35 - V_35)
    b_pten <- V35 * K_35 + V_35 * K35 + (V35 - V_35) * PTEN_tot
    c_pten <- -V_35 * K35 * PTEN_tot
    pten0 <- 2 * c_pten / (-b_pten - sqrt(b_pten^2 - 4 * a_pten * c_pten))

    e3(0) <- E3_tot
    e2(0) <- E2_tot
    shc(0) <- Shc_tot
    gs(0) <- GS_tot
    rasgdp(0) <- RasGDP_tot
    raf(0) <- Raf_tot
    mek(0) <- MEK_tot
    pp2a(0) <- PP2A_tot
    erk(0) <- ERK_tot
    pi3k(0) <- PI3K_tot
    pi2(0) <- PI_tot
    akt(0) <- Akt_tot
    pten(0) <- pten0
    ptenp(0) <- PTEN_tot - pten0

    # ---- Reaction rates (nM/s) ----
    # Receptor module
    v1 <- k1 * (e3 * hrg - Kd_1 * e3h) # V1 / R1
    v2 <- k2 * (e3h * e2 - Kd_2 * e23h) # V2 / R2
    v3 <- k3 * (e23h - Ke_3 * e23hp) # V3 / R3
    v4 <- V4 * e23hp / (K4 + e23hp) # V4 / R4
    v48 <- k48 * e23hp # V48 / R34
    v51 <- k51 * e3h # V51 / R37
    v52 <- k2 * (e3h_c * e2 - Kd_2 * e23h) # V52 / R38
    v53 <- k53 * e23h # V53 / R39
    v54 <- k3 * (e23h_c - Ke_3 * e23hp) # V54 / R40
    # Pertuzumab
    v49 <- k49 * (per * e2 - Kd_49 * e2_per) # V49 / R35
    v50 <- k50 * e2_per - k_50 * e2_per_c # V50 / R36
    # Shc / Grb2-Sos / Ras / Raf
    v5 <- k5 * (e23hp * shc - Kd_5 * e23hp_shc) # V5 / R5
    v6 <- k6 * e23hp_shc - k_6 * e23hp_shcp # V6 / R6
    v7 <- k7 * (e23hp_shcp * gs - Kd_7 * e23hp_shgs) # V7 / R7
    v8 <- k8 * (e23hp_shgs - Kd_8 * e23hp * shgs) # V8 / R8
    v9 <- k9 * (shgs - k_9 * shcp * gs) # V9 / R9
    v10 <- V10 * shcp / (K10 + shcp) # V10 / R10
    v11 <- k11 * rasgdp * shgs / (K11 + rasgdp) # V11 / R11
    v12 <- V12 * rasgtp / (K12 + rasgtp) # V12 / R12
    v13 <- k13 * raf * rasgtp / (K13 + raf) # V13 / R13
    v14 <- rafa * (k14 * akt_pi_pp / (rafa + K14) + k14_1 * E_raf / (rafa + K14_1)) # V14 / R14
    # MEK / ERK and PP2A
    v15 <- k15 * mek * rafa / (K15 + mek) # V15 / R15
    v16 <- k16 * mekp * pp2a - Kd_16 * mekp_pp2a # V16 / R16_1 (SBML form; see vignette)
    v17 <- k16_cat * mekp_pp2a # V17 / R16_2
    v18 <- k18 * mek_pp2a # V18 / R16_3
    v19 <- k15 * mekp * rafa / (K15 + mekp) # V19 / R17_1
    v20 <- k16 * (pp2a * mekpp - Kd_16 * mekpp_pp2a) # V20 / R18_1
    v21 <- k16_cat * mekpp_pp2a # V21 / R18_2
    v22 <- k22 * mekp_pp2a # V22 / R18_3
    v23 <- k23 * erk * mekpp / (K23 + erk) # V23 / R19
    v24 <- V24 * erkp / (K24 + erkp) # V24 / R20
    v25 <- k23 * erkp * mekpp / (K23 + erkp) # V25 / R21
    v26 <- V24 * erkpp / (K24 + erkpp) # V26 / R22
    # PI3K / PIP3 / PTEN
    v27 <- k27 * (e23hp * pi3k - Kd_27 * e23hp_pi3k) # V27 / R23
    v28 <- k28 * (e23hp_pi3k - k_28 * e23hp_pi3ka) # V28 / R24 (SBML form; see vignette)
    v29 <- k29 * e23hp_pi3ka - k_29 * e23hp * pi3ka # V29 / R25
    v30 <- V30 * pi3ka # V30 / R26
    v31 <- k31 * (pi2 * pi3ka - K_d31 * pi3ka_pi) # V31 / R27_1
    v55 <- k55 * pi3ka_pi # V55 / R41
    v56 <- k56 * pi3ka_pip3 # V56 / R42
    v32 <- k32 * (pip3 * pten - Kd_32 * pten_pip3) # V32 / R28_1
    v33 <- k33 * pten_pip3 # V33 / R28_2
    v34 <- k34 * pten_pi # V34 / R28_3
    v35 <- V35 * pten / (K35 + pten) # V35 / R28_4
    v36 <- V_35 * ptenp / (K_35 + ptenp) # V36 / R28_5
    v57 <- k57 * (pten * bpV - Kd_57 * pten_bpv) # V57 / R43
    v58 <- k58 * (pi3k * LY - Kd_58 * pi3k_ly) # V58 / R44
    # Akt / PDK1 / PP2A
    v39 <- k39 * (pip3 * akt - Kd_39 * akt_pip3) # V39 / R29
    v40 <- V40 * akt_pip3 / (K40 + akt_pip3) # V40 / R30
    v41 <- k41 * akt_pi_p * pp2a - Kd_41 * akt_pi_p_pp2a # V41 / R31_1 (SBML form; see vignette)
    v42 <- k42 * akt_pi_p_pp2a # V42 / R31_2
    v43 <- k43 * akt_pip3_pp2a # V43 / R31_3
    v44 <- V40 * akt_pi_p / (K40 + akt_pi_p) # V44 / R32
    v45 <- k41 * (akt_pi_pp * pp2a - Kd_41 * akt_pi_pp_pp2a) # V45 / R33_1
    v46 <- k42 * akt_pi_pp_pp2a # V46 / R33_2
    v47 <- k47 * akt_pi_p_pp2a # V47 / R33_3

    # ---- ODEs (Additional File 1, Eqs. S1.1-S1.54) ----
    d/dt(hrg) <- -v1
    d/dt(per) <- -v49
    d/dt(e3) <- -v1
    d/dt(e3h) <- v1 - v2 - v51
    d/dt(e3h_c) <- v51 - v52
    d/dt(e2) <- -v2 - v49 - v52
    d/dt(e23h) <- v2 - v3 + v4 + v52 - v53
    d/dt(e23h_c) <- v53 - v54
    d/dt(e23hp) <- v3 - v4 - v5 + v8 - v27 + v29 - v48 + v54
    d/dt(e2_per) <- v49 - v50
    d/dt(e2_per_c) <- v50
    d/dt(shc) <- -v5 + v10
    d/dt(e23hp_shc) <- v5 - v6
    d/dt(e23hp_shcp) <- v6 - v7
    d/dt(gs) <- -v7 + v9
    d/dt(e23hp_shgs) <- v7 - v8
    d/dt(shgs) <- v8 - v9
    d/dt(shcp) <- v9 - v10
    d/dt(rasgdp) <- -v11 + v12
    d/dt(rasgtp) <- v11 - v12
    d/dt(raf) <- -v13 + v14
    d/dt(rafa) <- v13 - v14
    d/dt(mek) <- -v15 + v18
    d/dt(mekp) <- v15 - v16 - v19 + v22
    d/dt(mekpp) <- v19 - v20
    d/dt(pp2a) <- -v16 + v18 - v20 + v22 - v41 + v43 - v45 + v47
    d/dt(mekp_pp2a) <- v16 - v17 + v21 - v22
    d/dt(mek_pp2a) <- v17 - v18
    d/dt(mekpp_pp2a) <- v20 - v21
    d/dt(erk) <- -v23 + v24
    d/dt(erkp) <- v23 - v24 - v25 + v26
    d/dt(erkpp) <- v25 - v26
    d/dt(pi3k) <- -v27 + v30 - v58
    d/dt(e23hp_pi3k) <- v27 - v28
    d/dt(e23hp_pi3ka) <- v28 - v29
    d/dt(pi3ka) <- v29 - v30 - v31 + v56
    d/dt(pi2) <- -v31 + v34
    d/dt(pi3ka_pi) <- v31 - v55
    d/dt(pi3ka_pip3) <- v55 - v56
    d/dt(pip3) <- -v32 - v39 + v56
    d/dt(pten) <- -v32 + v34 - v35 + v36 - v57
    d/dt(ptenp) <- v35 - v36
    d/dt(pten_pip3) <- v32 - v33
    d/dt(pten_pi) <- v33 - v34
    d/dt(pten_bpv) <- v57
    d/dt(pi3k_ly) <- v58
    d/dt(akt) <- -v39
    d/dt(akt_pip3) <- v39 - v40 + v43
    d/dt(akt_pi_p) <- v40 - v41 - v44 + v47
    d/dt(akt_pi_pp) <- v44 - v45
    d/dt(akt_pi_p_pp2a) <- v41 - v42 + v46 - v47
    d/dt(akt_pi_pp_pp2a) <- v45 - v46
    d/dt(akt_pip3_pp2a) <- v42 - v43

    # ---- Outputs ----
    # pAkt as defined in Methods S.1.4: the phosphorylated forms of Akt
    # normalised by total Akt (a fraction). Its integral over 0-T is the
    # sensitivity characteristic S_pAkt; S_pAkt,n = S_pAkt / T (Methods
    # S.3.1, T = 3600 s).
    pAkt <- (akt_pi_p + akt_pi_pp + akt_pi_p_pp2a + akt_pi_pp_pp2a) / Akt_tot
    d/dt(auc_pakt) <- pAkt
    # Total phosphorylated receptor, MEK and ERK (nM), the SBML assignment
    # rules tE3P, tMEKP and tERKP.
    tE3P <- e23hp + e23hp_pi3k + e23hp_pi3ka + e23hp_shc + e23hp_shcp + e23hp_shgs
    tMEKP <- mekp + mekpp + mekp_pp2a + mekpp_pp2a
    tERKP <- erkp + erkpp
  })
}
