Volkova_2025_SLE_ifn_qsp <- function() {
  description <- paste(
    "QSP. Volkova 2025 quantitative systems pharmacology model of type I",
    "interferon (IFN1)-mediated inflammation in adults with systemic lupus",
    "erythematosus (SLE), coupled to the PK of four IFN1-pathway-targeting",
    "monoclonal antibodies. Twenty ODEs. The biological core (7 ODEs) carries",
    "IFN1-producing immune cells (IC, proportional to plasmacytoid dendritic",
    "cells), free IFN-alpha and IFN-beta, free IFNAR1 receptor, the",
    "IFNAR1-IFN-alpha and IFNAR1-IFN-beta complexes and the 21-gene IFN gene",
    "signature (IFNGS), which is stimulated log-linearly by the summed",
    "IFN1-IFNAR1 complexes above a residual IFN1-independent floor. Two",
    "feedbacks: IFNGS relative to baseline amplifies IFN1 production per IC",
    "and accelerates IC loss from blood. Drugs: anifrolumab (anti-IFNAR1;",
    "two-compartment PK with linear plus saturable elimination, IFNAR1",
    "binding with complex internalisation, IV or SC), sifalimumab",
    "(anti-IFN-alpha; two-compartment PK with IFN-alpha binding, IV or SC),",
    "daxdilimab (anti-ILT7, pDC-depleting; one-compartment SC PK, Emax",
    "stimulation of IC loss) and litifilimab (anti-BDCA2; two-compartment SC",
    "PK with lag and Michaelis-Menten elimination, Hill inhibition of IFN1",
    "production). The system starts at its steady-state baseline. IIV on the",
    "IFNGS synthesis rate constant and on anifrolumab bioavailability, linear",
    "elimination and Vmax; residual error magnitudes are not reported, so",
    "none is encoded. Doses are in mg for all four drugs.",
    sep = " "
  )
  reference <- paste(
    "Volkova A, Sokolov V, Tettamanti F, Verma M, Ugolkov Y, Peskov K,",
    "Tang W, Kimko H (2025). An Integrative Mechanistic Model of Type 1",
    "IFN-Mediated Inflammation in Systemic Lupus Erythematosus. CPT",
    "Pharmacometrics Syst Pharmacol 14(7):1225-1236.",
    "doi:10.1002/psp4.70043. PMID: 40364448; PMCID: PMC12256560.",
    "Equations from Data S3 (Supplementary Files 1-2, Equations 1-25);",
    "parameter values from Table S2, the Monolix model code in Data S2 and",
    "the typical baseline IFN-alpha from Table S3.",
    sep = " "
  )
  vignette <- "Volkova_2025_SLE_ifn_qsp"

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Every state is either a paper-mechanistic biological species or a PK
  # compartment of one of four co-modelled drugs, none of which has a
  # registered drug suffix. Names follow the canonical role names with the
  # spelled-out INN appended.
  paper_specific_compartments <- c(
    "depot_anifrolumab",
    "central_anifrolumab",
    "peripheral1_anifrolumab",
    "depot_sifalimumab",
    "central_sifalimumab",
    "peripheral1_sifalimumab",
    "depot_daxdilimab",
    "central_daxdilimab",
    "depot_litifilimab",
    "central_litifilimab",
    "peripheral1_litifilimab",
    "ifna",
    "ifnb",
    "ifnar1",
    "ifna_ifnar1",
    "ifnb_ifnar1",
    "anifrolumab_ifnar1",
    "sifalimumab_ifna",
    "ifngs",
    "ic"
  )

  covariateData <- list(
    BGENE21_HIGH = list(
      description = "Baseline 21-gene IFN gene signature status: 1 = high (baseline IFNGS >= 2.75), 0 = low (< 2.75)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low)",
      notes = paste(
        "Volkova 2025 Data S3 Supplementary File 3: cut-off 2.75 chosen from",
        "the bimodal distribution of baseline IFNGS in the MUSE study",
        "(Figure S9). The source codes the opposite indicator (GRNUM = 1 for",
        "IFNGS-low; Table S2 'IFNGS-low status effect on kin_IFNGS' = -1.57),",
        "so the model uses (1 - BGENE21_HIGH) and the coefficient keeps the",
        "published sign. The paper's virtual populations assume every patient",
        "is IFNGS-high (Data S3 Supplementary File 4).",
        sep = " "
      ),
      source_name = "GRNUM (1 = IFNGS-low)"
    )
  )

  compartmentData <- list(
    depot_anifrolumab = list(
      analyte = "anifrolumab",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_anifrolumab = list(analyte = "anifrolumab", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_anifrolumab = list(analyte = "anifrolumab", units = "nmol", specimen = "tissue", verified = TRUE),
    depot_sifalimumab = list(
      analyte = "sifalimumab",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_sifalimumab = list(analyte = "sifalimumab", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_sifalimumab = list(analyte = "sifalimumab", units = "nmol", specimen = "tissue", verified = TRUE),
    depot_daxdilimab = list(analyte = "daxdilimab", units = "mg", specimen = "administration site", verified = TRUE),
    central_daxdilimab = list(analyte = "daxdilimab", units = "mg", specimen = "plasma", verified = TRUE),
    depot_litifilimab = list(analyte = "litifilimab", units = "ug", specimen = "administration site", verified = TRUE),
    central_litifilimab = list(analyte = "litifilimab", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_litifilimab = list(analyte = "litifilimab", units = "ug", specimen = "tissue", verified = TRUE),
    ifna = list(analyte = "free interferon alpha", units = "nmol", specimen = "plasma", verified = TRUE),
    ifnb = list(analyte = "free interferon beta", units = "nmol", specimen = "plasma", verified = TRUE),
    ifnar1 = list(analyte = "free IFNAR1 receptor subunit", units = "nmol", specimen = "plasma", verified = TRUE),
    ifna_ifnar1 = list(analyte = "IFNAR1-IFN-alpha complex", units = "nmol", specimen = "plasma", verified = TRUE),
    ifnb_ifnar1 = list(analyte = "IFNAR1-IFN-beta complex", units = "nmol", specimen = "plasma", verified = TRUE),
    anifrolumab_ifnar1 = list(
      analyte = "anifrolumab-IFNAR1 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    sifalimumab_ifna = list(
      analyte = "sifalimumab-IFN-alpha complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    ifngs = list(
      analyte = "21-gene IFN gene signature score",
      units = "unitless fold change",
      specimen = "whole blood",
      verified = TRUE
    ),
    ic = list(
      analyte = "IFN1-producing immune cells (pDC surrogate)",
      units = "cell",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 254L,
    n_studies = 11L,
    disease_state = paste(
      "Adults with moderate-to-severe systemic lupus erythematosus (the",
      "litifilimab and daxdilimab validation trials also enrolled cutaneous",
      "lupus / autoimmune-disease patients).",
      sep = " "
    ),
    dose_range = paste(
      "Anifrolumab 100-1000 mg IV Q4W and 150-300 mg SC Q2W; sifalimumab up",
      "to 30 mg/kg IV and SC regimens; daxdilimab 5-150 mg SC; litifilimab",
      "450 mg SC Q4W with an extra dose at week 2 (Table S1, Figure S1,",
      "Data S3 Supplementary File 4).",
      sep = " "
    ),
    notes = paste(
      "Patient-level calibration data: 254 patients of the Phase IIb MUSE",
      "anifrolumab trial (NCT01438489; anifrolumab PK, IFN-alpha and IFNGS).",
      "Daxdilimab PK/PD parameters were fitted to aggregated Phase I data",
      "(NCT03817424). Study-level validation data from five anifrolumab, four",
      "sifalimumab, one daxdilimab and one litifilimab trials (Table S1).",
      "Demographics are not reported.",
      sep = " "
    )
  )

  ini({
    # ---- Systemic IFN1 / IFNGS parameters (Table S2, first block) ----------
    bl_ifna <- fixed(6.35); label("Baseline plasma IFN-alpha (pg/mL)") # Table S3 typical IFNa_BL = 6.35 pg/mL; Table S2 lists it as an individual per-subject covariate (code: log-normal individual parameter, variance not reported)
    ifna_fraction <- fixed(0.99); label("Fraction of baseline IFN1 that is IFN-alpha (unitless)") # Table S2 IFNa_fraction = 0.99 (assumed); code A_to_B = 0.99
    mw_ifna <- fixed(18000); label("IFN-alpha molecular weight (g/mol)") # Table S2 MW_IFNa = 18 kDa; code MW_IFNa = 18000
    mw_ifnb <- fixed(22000); label("IFN-beta molecular weight (g/mol)") # Table S2 MW_IFNb = 22 kDa; code MW_IFNb = 22000
    bl_ifnar1 <- fixed(0.26); label("Baseline free IFNAR1 concentration (nmol/L)") # Table S2 IFNAR1_BL = 0.26 nmol/L (ref 1, Wang 2013)
    kout_ifnar1 <- fixed(24.2); label("IFNAR1 elimination rate constant (1/day)") # Table S2 kout_IFNAR1 = 24.2; code kdeg_IFNAR1 = 24.2
    kout_ifna_ifnar1 <- fixed(16.6); label("IFNAR1-IFN-alpha complex elimination rate constant (1/day)") # Table S2 = 16.6; code kdeg_IFNa_IFNAR1 = 16.6
    kout_ifnb_ifnar1 <- fixed(16.6); label("IFNAR1-IFN-beta complex elimination rate constant (1/day)") # Table S2 = 16.6; code kdeg_IFNb_IFNAR1 = 16.6
    kon_ifna <- fixed(43.2); label("IFN-alpha/IFNAR1 association rate constant (1/(nmol/L*day))") # Table S2 kon_IFNa = 43.2; code kon_IFNAR1 = 43.2
    kon_ifnb <- fixed(43.2); label("IFN-beta/IFNAR1 association rate constant (1/(nmol/L*day))") # Table S2 kon_IFNb = 43.2; code uses the same kon_IFNAR1
    kd_ifna <- fixed(500); label("IFN-alpha/IFNAR1 equilibrium dissociation constant (nmol/L)") # Table S2 KD_IFNa = 500; code Kd_IFNa_IFNAR1 = 500
    kd_ifnb <- fixed(100); label("IFN-beta/IFNAR1 equilibrium dissociation constant (nmol/L)") # Table S2 KD_IFNb = 100; code Kd_IFNb_IFNAR1 = 100
    kout_ifna <- fixed(8.3); label("Free IFN-alpha elimination rate constant (1/day)") # Table S2 kout_IFNa = 8.3; code kdeg_IFNa = 8.3
    kout_ifnb <- fixed(4.16); label("Free IFN-beta elimination rate constant (1/day)") # Table S2 kout_IFNb = 4.16; code kdeg_IFNb = 4.16
    lkin_ifngs <- log(2.7); label("IFNGS synthesis rate constant, IFNGS-high patients (1/day)") # Table S2 kin_IFNGS = 2.7 (CV 20.7), RSE 5.65, estimated
    e_bgene21_low_kin_ifngs <- -1.57; label("IFNGS-low status effect on log kin_IFNGS (unitless)") # Table S2 beta_kin_IFNGS = -1.57, RSE 8.23, estimated; Results: exp(1.57) = 4.8-fold lower response in IFNGS-low patients
    kout_ifngs <- fixed(1.92); label("IFNGS elimination rate constant (1/day)") # Table S2 kout_IFNGS = 1.92 (ref 1)
    kn <- fixed(1); label("Slope of IFNGS stimulation by IFN1-IFNAR1 complexes (1/(nmol*1e8))") # Table S2 k_n = 1 (assumed); Equation 7 and code multiply the complex amount by 10^8
    ifngs_ll <- fixed(0.764); label("Residual IFNGS in the absence of IFN1 signalling (unitless)") # Table S2 IFNGS_LL = min(IFNGS, 0.764) per subject; code lowIFNGS = 0.764
    lvc <- fixed(log(2.75)); label("Plasma volume shared by all four drugs and the IFN1 species (L)") # Table S2 Vc = 2.75 L (ref 10); code Vc = 2.75
    bl_ic <- fixed(1e7); label("Baseline IFN1-producing immune cells (cell)") # Table S2 IC_BL = 10^7 cells (ref 11); code pDC_BL = 10000000
    kout_ic <- fixed(0.1); label("Immune cell elimination rate constant (1/day)") # Table S2 kout_IC = 0.1 (ref 12); code kout_pDC = 0.1

    # ---- Anifrolumab PK and IFNAR1 binding (Table S2, second block) -------
    lka_anifrolumab <- fixed(log(0.2)); label("Anifrolumab SC absorption rate constant (1/day)") # Table S2 ka_ANF = 0.2 (assumed); Results: fixed from internal data; code ka_ani = 0.2
    lk12_anifrolumab <- log(4.69); label("Anifrolumab central-to-peripheral rate constant (1/day)") # Table S2 k12 = 4.69, RSE 2.35 (the table's description column swaps the k12/k21 directions; Equation 9 and the code use k12 for central -> peripheral)
    lk21_anifrolumab <- log(3.98); label("Anifrolumab peripheral-to-central rate constant (1/day)") # Table S2 k21 = 3.98, RSE 7.05
    logitfbio_anifrolumab <- log(0.79 / (1 - 0.79)); label("Logit of anifrolumab IV bioavailability (unitless)") # Table S2 Fbio_ANF_IV = 0.79 (omega 1.05), RSE 4.59; code: logitNormal distribution
    fsc_anifrolumab <- fixed(0.87); label("Anifrolumab SC bioavailability relative to IV (unitless)") # Table S2 Fbio_ANF_SC = 0.87 * Fbio_ANF_IV (internal data); code Fbio_ani_sc = 0.87
    lkel_anifrolumab <- log(0.0486); label("Anifrolumab linear elimination rate constant (1/day)") # Table S2 kel_ANF = 0.0486 (CV 28), RSE 5.81
    lvm_anifrolumab <- log(21.4); label("Anifrolumab maximum saturable elimination rate (nmol/day)") # Table S2 Vm = 21.4 (CV 32), RSE 17.4
    lkm_anifrolumab <- log(38.9); label("Anifrolumab amount at half-maximal saturable elimination (nmol)") # Table S2 Km = 38.9, RSE 18.3
    mw_anifrolumab <- fixed(145119.2); label("Anifrolumab molecular weight (g/mol)") # Data S2 code MW = 145119.20; Table S2 and Data S3 print 148 kDa. Code value reproduces the published anifrolumab IC50 of 132 ng/mL (see vignette)
    kd_anifrolumab <- fixed(0.067); label("Anifrolumab/IFNAR1 equilibrium dissociation constant (nmol/L)") # Table S2 KD_ANF = 0.067 (ref 14); code Kd = 0.067
    kon_anifrolumab <- fixed(233); label("Anifrolumab/IFNAR1 association rate constant (1/(nmol/L*day))") # Table S2 kon_ANF = 233 (ref 14); code kon = 233
    kint_anifrolumab <- fixed(77.1); label("Anifrolumab-IFNAR1 complex internalisation rate constant (1/day)") # Data S2 code kint = 77.1; Table S2 prints 77.4 (ref 1); see vignette

    # ---- Sifalimumab PK and IFN-alpha binding (Table S2, third block) -----
    lka_sifalimumab <- fixed(log(0.2)); label("Sifalimumab SC absorption rate constant (1/day)") # Table S2 ka_SIF = 0.2 (assumed); code ka_Sif = 0.2
    lfdepot_sifalimumab <- fixed(log(0.7)); label("Sifalimumab SC bioavailability (unitless)") # Table S2 Fbio_SIF = 0.7 (internal data); code Fbio_Sif = 0.7
    lcl_sifalimumab <- fixed(log(0.19)); label("Sifalimumab linear clearance (L/day)") # Table S2 CL_SIF = 0.19 (ref 15, Zheng 2016); code CL_Sif = 0.19
    lvp_sifalimumab <- fixed(log(1.89)); label("Sifalimumab peripheral volume (L)") # Table S2 Vp_SIF = 1.89 (ref 15); code Vp_Sif = 1.89
    lq_sifalimumab <- fixed(log(0.32)); label("Sifalimumab intercompartmental clearance (L/day)") # Table S2 Q_SIF = 0.32 (ref 15); code Q_Sif = 0.32
    kd_sifalimumab <- fixed(0.18); label("Sifalimumab/IFN-alpha equilibrium dissociation constant (nmol/L)") # Table S2 KD_SIF = 0.18 (internal in vitro data; the table's description cell is mislabelled); code Kd_Sif = 0.18
    koff_sifalimumab <- fixed(11.232); label("Sifalimumab/IFN-alpha dissociation rate constant (1/day)") # Table S2 koff_SIF = 11.232 (internal in vitro data; description cell mislabelled); code koff_Sif = 11.232
    mw_sifalimumab <- fixed(147000); label("Sifalimumab molecular weight (g/mol)") # Table S2 MW_SIF = 147 kDa; code MW_Sif = 147000
    kout_sifalimumab_ifna <- fixed(0.033); label("Sifalimumab-IFN-alpha complex elimination rate constant (1/day)") # Table S2 kout_SIF-IFNa = 0.033 (assumed); code kdeg_Sif = 0.033

    # ---- Daxdilimab PK and effect on IC turnover (Table S2, fourth block) -
    lka_daxdilimab <- log(0.573); label("Daxdilimab SC absorption rate constant (1/day)") # Table S2 ka_DAX = 0.573, RSE 22.2
    lfdepot_daxdilimab <- log(0.147); label("Daxdilimab SC bioavailability (unitless)") # Table S2 Fbio_DAX = 0.147, RSE 10.9
    lkel_daxdilimab <- log(0.0485); label("Daxdilimab elimination rate constant (1/day)") # Table S2 kel_DAX = 0.0485, RSE 8.52
    mw_daxdilimab <- fixed(145000); label("Daxdilimab molecular weight (g/mol)") # Table S2 MW_DAX = 145 kDa; code MW_VIB7734 = 145000
    lec50_daxdilimab <- log(0.52); label("Daxdilimab concentration at half-maximal stimulation of IC loss (nmol/L)") # Table S2 EC50_DAX = 0.52, RSE 38.2
    emax_daxdilimab <- fixed(4); label("Daxdilimab maximal stimulation of IC loss (unitless)") # Table S2 Emax_DAX = 4 (calculated from observed data, ref 17); code Emax_VIB = 4

    # ---- Litifilimab PK and effect on IFN1 production (Table S2, fifth) ---
    lka_litifilimab <- fixed(log(0.234)); label("Litifilimab SC absorption rate constant (1/day)") # Table S2 ka_LIT = 0.234 (ref 18, Hartmann 2020); code ka_Lit = 0.234
    lfdepot_litifilimab <- fixed(log(0.713)); label("Litifilimab SC bioavailability (unitless)") # Table S2 Fbio_LIT = 0.713 (ref 18); code Fbio_Lit = 0.713
    tlag_litifilimab <- fixed(0.09); label("Litifilimab absorption lag time (day)") # Table S2 Tlag_LIT = 0.09 (ref 18); code Lag_Lit = 0.09
    lkel_litifilimab <- fixed(log(0.071)); label("Litifilimab linear elimination rate constant (1/day)") # Table S2 kel_LIT = 0.071 1/d (ref 18); code kel_Lit = 0.071 with CL_Lit = kel_Lit * Vc
    lvp_litifilimab <- fixed(log(2.31)); label("Litifilimab peripheral volume (L)") # Table S2 Vp_LIT = 2.31 (ref 18); code Vp_Lit = 2.31
    lq_litifilimab <- fixed(log(0.4896)); label("Litifilimab intercompartmental clearance (L/day)") # Table S2 Q_LIT = 0.4896 (ref 18); code Q_Lit = 0.4896
    lkm_litifilimab <- fixed(log(547)); label("Litifilimab concentration at half-maximal saturable elimination (ug/L)") # Table S2 Km_LIT = 547 (ref 18); code Km_Lit = 547
    lvm_litifilimab <- fixed(log(256.8)); label("Litifilimab maximal saturable elimination rate (ug/day)") # Table S2 Vm_LIT = 256.8 (ref 18); code Vm_Lit = 256.8
    mw_litifilimab <- fixed(146000); label("Litifilimab molecular weight (g/mol)") # Table S2 MW_LIT = 146 kDa; code MW_Lit = 146000
    ic50_litifilimab <- fixed(350); label("Litifilimab concentration at half-maximal inhibition of IFN1 production (ug/L)") # Table S2 IC50_LIT = 350 (ref 18); code IC50_Lit = 350
    imax_litifilimab <- fixed(0.875); label("Litifilimab maximal inhibition of IFN1 production (fraction)") # Table S2 Imax_LIT = 0.875 (ref 18); code Imax_Lit = 0.875
    hill_litifilimab <- fixed(3.63); label("Litifilimab Hill coefficient (unitless)") # Table S2 Hill_LIT = 3.63 (ref 18); code Hill_Lit = 3.63

    # ---- IIV (Table S2 CV; omega^2 = log(1 + CV^2) for log-normal) --------
    etalkin_ifngs ~ log(1 + 0.207^2) # Table S2 kin_IFNGS CV 20.7
    etalogitfbio_anifrolumab ~ 1.05^2 # Table S2 Fbio_ANF_IV omega = 1.05 (logit scale SD)
    etalkel_anifrolumab ~ log(1 + 0.28^2) # Table S2 kel_ANF CV 28
    etalvm_anifrolumab ~ log(1 + 0.32^2) # Table S2 Vm CV 32
  })
  model({
    # ---- Individual parameters --------------------------------------------
    kin_ifngs <- exp(lkin_ifngs + etalkin_ifngs + e_bgene21_low_kin_ifngs * (1 - BGENE21_HIGH)) # code kin_IFNGS_t = kin_IFNGS * exp(beta_gr * GRNUM), GRNUM = 1 for IFNGS-low
    vc <- exp(lvc)

    ka_anifrolumab <- exp(lka_anifrolumab)
    k12_anifrolumab <- exp(lk12_anifrolumab)
    k21_anifrolumab <- exp(lk21_anifrolumab)
    fbio_anifrolumab <- expit(logitfbio_anifrolumab + etalogitfbio_anifrolumab)
    kel_anifrolumab <- exp(lkel_anifrolumab + etalkel_anifrolumab)
    vm_anifrolumab <- exp(lvm_anifrolumab + etalvm_anifrolumab)
    km_anifrolumab <- exp(lkm_anifrolumab)

    ka_sifalimumab <- exp(lka_sifalimumab)
    fdepot_sifalimumab <- exp(lfdepot_sifalimumab)
    cl_sifalimumab <- exp(lcl_sifalimumab)
    vp_sifalimumab <- exp(lvp_sifalimumab)
    q_sifalimumab <- exp(lq_sifalimumab)

    ka_daxdilimab <- exp(lka_daxdilimab)
    fdepot_daxdilimab <- exp(lfdepot_daxdilimab)
    kel_daxdilimab <- exp(lkel_daxdilimab)
    ec50_daxdilimab <- exp(lec50_daxdilimab)

    ka_litifilimab <- exp(lka_litifilimab)
    fdepot_litifilimab <- exp(lfdepot_litifilimab)
    kel_litifilimab <- exp(lkel_litifilimab)
    vp_litifilimab <- exp(lvp_litifilimab)
    q_litifilimab <- exp(lq_litifilimab)
    km_litifilimab <- exp(lkm_litifilimab)
    vm_litifilimab <- exp(lvm_litifilimab)

    # ---- Derived binding constants (Data S3 Equations 11 and 15; code) ----
    koff_ifna <- kd_ifna * kon_ifna
    koff_ifnb <- kd_ifnb * kon_ifnb
    koff_anifrolumab <- kd_anifrolumab * kon_anifrolumab
    kon_sifalimumab <- koff_sifalimumab / kd_sifalimumab

    # ---- Steady-state baseline (Data S3 Supplementary File 2) -------------
    # IFN amounts in nmol: pg/mL * L / (g/mol) = nmol.
    bl_ifnb <- bl_ifna / ifna_fraction * (1 - ifna_fraction) # pg/mL; code IFNb_BL_pgmL
    ifna_bl_amt <- bl_ifna * vc / mw_ifna
    ifnb_bl_amt <- bl_ifnb * vc / mw_ifnb
    ifnar1_bl_amt <- bl_ifnar1 * vc

    ksyn_ifna <- ifna_bl_amt *
      (kout_ifna + kon_ifna * ifnar1_bl_amt / vc -
        (koff_ifna * kon_ifna * ifnar1_bl_amt / vc) / (koff_ifna + kout_ifna_ifnar1)) /
      bl_ic # Equation 21
    ksyn_ifnb <- ifnb_bl_amt *
      (kout_ifnb + kon_ifnb * ifnar1_bl_amt / vc -
        (koff_ifnb * kon_ifnb * ifnar1_bl_amt / vc) / (koff_ifnb + kout_ifnb_ifnar1)) /
      bl_ic # Equation 21 mirrored for IFN-beta
    ifna_ifnar1_bl <- kon_ifna * ifna_bl_amt * ifnar1_bl_amt / vc / (koff_ifna + kout_ifna_ifnar1) # Equation 24
    ifnb_ifnar1_bl <- kon_ifnb * ifnb_bl_amt * ifnar1_bl_amt / vc / (koff_ifnb + kout_ifnb_ifnar1) # Equation 24 mirrored
    ksyn_ifnar1 <- kout_ifnar1 * ifnar1_bl_amt +
      kon_ifna * ifna_bl_amt * ifnar1_bl_amt / vc - koff_ifna * ifna_ifnar1_bl +
      kon_ifnb * ifnb_bl_amt * ifnar1_bl_amt / vc - koff_ifnb * ifnb_ifnar1_bl # Equation 22
    ifngs_bl <- (kin_ifngs * log(1 + kn * (ifna_ifnar1_bl + ifnb_ifnar1_bl) * 1e8) + kout_ifngs * ifngs_ll) /
      kout_ifngs # Equation 25
    kin_ic <- kout_ic * bl_ic # Equation 23

    # ---- Drug concentrations driving the system ---------------------------
    cc_anifrolumab_nm <- central_anifrolumab / vc
    cc_sifalimumab_nm <- central_sifalimumab / vc
    cp_sifalimumab_nm <- peripheral1_sifalimumab / vp_sifalimumab
    cc_daxdilimab_nm <- central_daxdilimab / vc * 1e6 / mw_daxdilimab # mg/L -> nmol/L
    cc_litifilimab_ugl <- central_litifilimab / vc
    cp_litifilimab_ugl <- peripheral1_litifilimab / vp_litifilimab

    # Feedback from IFNGS on IFN1 production and IC loss; drug effects
    fb_ifngs <- ifngs / ifngs_bl
    inh_litifilimab <- imax_litifilimab * cc_litifilimab_ugl^hill_litifilimab /
      (ic50_litifilimab^hill_litifilimab + cc_litifilimab_ugl^hill_litifilimab)
    stim_daxdilimab <- emax_daxdilimab * cc_daxdilimab_nm / (ec50_daxdilimab + cc_daxdilimab_nm)

    # ---- Drug PK (Equations 8-10, 12-14, 16-17, 18-20) --------------------
    d/dt(depot_anifrolumab) <- -ka_anifrolumab * depot_anifrolumab
    d/dt(central_anifrolumab) <- ka_anifrolumab * depot_anifrolumab -
      k12_anifrolumab * central_anifrolumab + k21_anifrolumab * peripheral1_anifrolumab -
      kel_anifrolumab * central_anifrolumab -
      vm_anifrolumab * central_anifrolumab / (km_anifrolumab + central_anifrolumab)
    d/dt(peripheral1_anifrolumab) <- k12_anifrolumab * central_anifrolumab - k21_anifrolumab * peripheral1_anifrolumab

    d/dt(depot_sifalimumab) <- -ka_sifalimumab * depot_sifalimumab
    d/dt(central_sifalimumab) <- ka_sifalimumab * depot_sifalimumab -
      cl_sifalimumab * cc_sifalimumab_nm -
      q_sifalimumab * (cc_sifalimumab_nm - cp_sifalimumab_nm) -
      kon_sifalimumab * cc_sifalimumab_nm * ifna + koff_sifalimumab * sifalimumab_ifna
    d/dt(peripheral1_sifalimumab) <- q_sifalimumab * (cc_sifalimumab_nm - cp_sifalimumab_nm)

    d/dt(depot_daxdilimab) <- -ka_daxdilimab * depot_daxdilimab
    d/dt(central_daxdilimab) <- ka_daxdilimab * depot_daxdilimab - kel_daxdilimab * central_daxdilimab

    d/dt(depot_litifilimab) <- -ka_litifilimab * depot_litifilimab
    d/dt(central_litifilimab) <- ka_litifilimab * depot_litifilimab -
      q_litifilimab * (cc_litifilimab_ugl - cp_litifilimab_ugl) -
      kel_litifilimab * central_litifilimab -
      vm_litifilimab * cc_litifilimab_ugl / (km_litifilimab + cc_litifilimab_ugl)
    d/dt(peripheral1_litifilimab) <- q_litifilimab * (cc_litifilimab_ugl - cp_litifilimab_ugl)

    # ---- IFN1 signalling core (Equations 1-7, 11, 15) ---------------------
    d/dt(ifna) <- ksyn_ifna * ic * fb_ifngs * (1 - inh_litifilimab) -
      kout_ifna * ifna -
      kon_ifna * ifna * ifnar1 / vc + koff_ifna * ifna_ifnar1 -
      kon_sifalimumab * cc_sifalimumab_nm * ifna + koff_sifalimumab * sifalimumab_ifna # Equation 2
    d/dt(ifnb) <- ksyn_ifnb * ic * fb_ifngs * (1 - inh_litifilimab) -
      kout_ifnb * ifnb -
      kon_ifnb * ifnb * ifnar1 / vc + koff_ifnb * ifnb_ifnar1 # Equation 3
    d/dt(ifnar1) <- ksyn_ifnar1 - kout_ifnar1 * ifnar1 -
      kon_ifna * ifna * ifnar1 / vc + koff_ifna * ifna_ifnar1 -
      kon_ifnb * ifnb * ifnar1 / vc + koff_ifnb * ifnb_ifnar1 -
      kon_anifrolumab * cc_anifrolumab_nm * ifnar1 + koff_anifrolumab * anifrolumab_ifnar1 # Equation 4
    d/dt(ifna_ifnar1) <- kon_ifna * ifna * ifnar1 / vc - koff_ifna * ifna_ifnar1 - kout_ifna_ifnar1 * ifna_ifnar1 # Equation 5
    d/dt(ifnb_ifnar1) <- kon_ifnb * ifnb * ifnar1 / vc - koff_ifnb * ifnb_ifnar1 - kout_ifnb_ifnar1 * ifnb_ifnar1 # Equation 6
    d/dt(anifrolumab_ifnar1) <- kon_anifrolumab * cc_anifrolumab_nm * ifnar1 -
      koff_anifrolumab * anifrolumab_ifnar1 - kint_anifrolumab * anifrolumab_ifnar1 # Equation 11
    d/dt(sifalimumab_ifna) <- kon_sifalimumab * cc_sifalimumab_nm * ifna -
      koff_sifalimumab * sifalimumab_ifna - kout_sifalimumab_ifna * sifalimumab_ifna # Equation 15
    d/dt(ifngs) <- kin_ifngs * log(1 + kn * (ifna_ifnar1 + ifnb_ifnar1) * 1e8) -
      kout_ifngs * (ifngs - ifngs_ll) # Equation 7
    d/dt(ic) <- kin_ic - kout_ic * fb_ifngs * (1 + stim_daxdilimab) * ic # Equation 1

    ifna(0) <- ifna_bl_amt
    ifnb(0) <- ifnb_bl_amt
    ifnar1(0) <- ifnar1_bl_amt
    ifna_ifnar1(0) <- ifna_ifnar1_bl
    ifnb_ifnar1(0) <- ifnb_ifnar1_bl
    ifngs(0) <- ifngs_bl
    ic(0) <- bl_ic

    # ---- Dosing: mg doses converted to the state amount units (code iv()) -
    f(central_anifrolumab) <- fbio_anifrolumab * 1e6 / mw_anifrolumab # IV, mg -> nmol
    f(depot_anifrolumab) <- fsc_anifrolumab * fbio_anifrolumab * 1e6 / mw_anifrolumab # SC, mg -> nmol
    f(central_sifalimumab) <- 1e6 / mw_sifalimumab # IV, mg -> nmol
    f(depot_sifalimumab) <- fdepot_sifalimumab * 1e6 / mw_sifalimumab # SC, mg -> nmol
    f(depot_daxdilimab) <- fdepot_daxdilimab # SC, mg
    f(depot_litifilimab) <- fdepot_litifilimab * 1000 # SC, mg -> ug
    alag(depot_litifilimab) <- tlag_litifilimab

    # ---- Outputs ----------------------------------------------------------
    Cc_anifrolumab <- cc_anifrolumab_nm * mw_anifrolumab / 1000 # ng/mL
    Cc_sifalimumab <- cc_sifalimumab_nm * mw_sifalimumab / 1000 # ng/mL
    Cc_daxdilimab <- central_daxdilimab / vc * 1000 # ng/mL
    Cc_litifilimab <- cc_litifilimab_ugl # ug/L = ng/mL
    cc_litifilimab_nm <- cc_litifilimab_ugl / mw_litifilimab * 1000 # nmol/L; code Cc_Lit_nM
    ifna_pgml <- ifna / vc * mw_ifna # free IFN-alpha, pg/mL
    ifngs_cfb <- 100 * (ifngs - ifngs_bl) / ifngs_bl # percent change from baseline
    ifna_cfb <- 100 * (ifna - ifna_bl_amt) / ifna_bl_amt # percent change from baseline
    ic_cfb <- 100 * (ic - bl_ic) / bl_ic # percent change from baseline (pDC surrogate)
    ifnar1_cfb <- 100 * (ifnar1 - ifnar1_bl_amt) / ifnar1_bl_amt
  })
}
