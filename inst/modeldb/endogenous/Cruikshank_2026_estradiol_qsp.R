Cruikshank_2026_estradiol_qsp <- function() {
  description <- paste(
    "QSP. Mechanistic model of hepatocyte hydrogen peroxide (H2O2) and",
    "glutathione (GSH) metabolism with estradiol and progesterone regulation",
    "(Cruikshank 2026): 17 mass-balance ODEs for the transsulfuration",
    "pathway (homocysteine, cystathionine, cysteine, glutamyl-cysteine), GSH",
    "and GSSG in cytosol, mitochondria and blood, cytosolic and mitochondrial",
    "H2O2, and reduced/oxidised thioredoxin in both compartments. H2O2 is",
    "scavenged by glutathione peroxidase (GPX), thioredoxin peroxidase (TPX)",
    "and catalase (CAT). Estradiol stimulates CBS, GCL and GPX (saturable",
    "Michaelis-Menten activation relative to 0.09 nM); progesterone inhibits",
    "GPX (relative to 0.318 nM); homocysteine inhibits cytosolic GPX and",
    "CAT. Serum estradiol and progesterone (nM) are supplied as",
    "(time-varying) covariates. Cystathionine beta-synthase deficiency",
    "(CBSD) is simulated with f_cbs = 0.025. Deterministic: the authors'",
    "virtual population (20,000 plausible patients varying 20 parameters",
    "within 50-150% of nominal) is not a parametric distribution, so no IIV",
    "or residual error is included. Time in hours, concentrations in uM.",
    sep = " "
  )
  reference <- paste(
    "Cruikshank A, Nijhout HF, Reed MC (2026).",
    "The Role of Estradiol in Oxidative Stress in Health and Cystathionine",
    "beta-Synthase Deficiency: A Mechanistic Model.",
    "CPT Pharmacometrics Syst Pharmacol 15(4):e70243. doi:10.1002/psp4.70243.",
    "Parameter values and the 38 explicit rate laws are from the authors'",
    "SimBiology SBML export (Data S2, model 'H2O2Dynamics'); the five rate",
    "laws that the export leaves empty (GCL, cytosolic and mitochondrial",
    "GPX and TPX) are from the Data S1 Supplementary Methods equations.",
    sep = " "
  )
  vignette <- "Cruikshank_2026_estradiol_qsp"

  units <- list(time = "h", dosing = "umol", concentration = "uM")

  paper_specific_compartments <- c(
    "hcy",
    "cyst",
    "ccys",
    "glutcys",
    "cgsh",
    "cgssg",
    "bgsh",
    "bgssg",
    "bcys",
    "ch2o2",
    "mh2o2",
    "mgsh",
    "mgssg",
    "ctrxred",
    "ctrxox",
    "mtrxred",
    "mtrxox"
  )

  compartmentData <- list(
    hcy = list(analyte = "homocysteine (hepatocyte cytosol)", units = "uM", specimen = "tissue", verified = TRUE),
    cyst = list(analyte = "cystathionine (hepatocyte cytosol)", units = "uM", specimen = "tissue", verified = TRUE),
    ccys = list(analyte = "cysteine (hepatocyte cytosol)", units = "uM", specimen = "tissue", verified = TRUE),
    glutcys = list(
      analyte = "gamma-glutamylcysteine (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    cgsh = list(
      analyte = "reduced glutathione (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    cgssg = list(
      analyte = "glutathione disulfide (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    bgsh = list(analyte = "reduced glutathione (blood)", units = "uM", specimen = "whole blood", verified = TRUE),
    bgssg = list(analyte = "glutathione disulfide (blood)", units = "uM", specimen = "whole blood", verified = TRUE),
    bcys = list(analyte = "cysteine (blood)", units = "uM", specimen = "whole blood", verified = TRUE),
    ch2o2 = list(
      analyte = "hydrogen peroxide (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    mh2o2 = list(
      analyte = "hydrogen peroxide (hepatocyte mitochondria)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    mgsh = list(
      analyte = "reduced glutathione (hepatocyte mitochondria)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    mgssg = list(
      analyte = "glutathione disulfide (hepatocyte mitochondria)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    ctrxred = list(
      analyte = "reduced thioredoxin (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    ctrxox = list(
      analyte = "oxidised thioredoxin (hepatocyte cytosol)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    mtrxred = list(
      analyte = "reduced thioredoxin (hepatocyte mitochondria)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    ),
    mtrxox = list(
      analyte = "oxidised thioredoxin (hepatocyte mitochondria)",
      units = "uM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    ESTRADIOL = list(
      description = "Serum estradiol (17-beta-estradiol) concentration driving the estradiol activation of CBS, GCL and GPX.",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Constant or time-varying. Cruikshank 2026 Methods (Hormone",
        "Concentrations): 0.09 nM for males and postmenopausal females;",
        "0.4368 nM as the 28-day average for premenopausal females; a",
        "spline-interpolated 28-day menstrual-cycle profile (Figure 5A) for",
        "the cycle simulations. Every activation term is zero at 0.09 nM.",
        sep = " "
      ),
      source_name = "E2"
    ),
    PROGESTERONE = list(
      description = "Serum progesterone concentration driving the progesterone inhibition of cytosolic and mitochondrial GPX.",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Constant or time-varying. Cruikshank 2026 Methods (Hormone",
        "Concentrations): 0.318 nM for males and postmenopausal females;",
        "12.4859 nM as the 28-day average for premenopausal females; the",
        "menstrual-cycle profile in Figure 5B. The inhibition term is zero",
        "at 0.318 nM.",
        sep = " "
      ),
      source_name = "P4"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20000,
    disease_state = paste(
      "Healthy adult liver (males, postmenopausal and premenopausal",
      "females) and classical homocystinuria (cystathionine beta-synthase",
      "deficiency, CBSD).",
      sep = " "
    ),
    dose_range = "No drug dosing; serum estradiol 0.09-2.5 nM at steady state and the 28-day menstrual-cycle estradiol and progesterone profiles.",
    notes = paste(
      "Not fitted to individual data. Km values come from the enzyme",
      "literature; Vmax values and the mitochondrial H2O2 transport rate",
      "were calibrated so that substrates stay in physiological ranges and",
      "to reproduce rat liver mitochondria H2O2-removal experiments (Data",
      "S1 section 2). The 20,000 virtual individuals (Table S3) were",
      "generated with a plausible-patient algorithm after a Sobol",
      "sensitivity analysis on cytosolic H2O2; the 15 most influential",
      "kinetic parameters plus the five estradiol-activation parameters",
      "were varied between 50% and 150% of nominal (Table S2). The",
      "progesterone parameters were held at nominal.",
      sep = " "
    )
  )

  ini({
    # All values are the nominal parameters of the authors' SBML export
    # (Data S2, model H2O2Dynamics); the SBML parameter name is quoted after
    # each line and the matching Data S1 table where one prints the value.
    # Nothing was estimated, so every value is fixed. Concentrations uM,
    # velocities uM/h, first-order rates 1/h.

    # --- Compartment volumes (relative) ---------------------------------
    vol_cytosol <- fixed(0.85); label("Relative cytosol volume (fraction of hepatocyte)") # SBML compartment Cytosol; Data S1: mitochondria are 15% of the cell
    vol_mitochondria <- fixed(0.15); label("Relative mitochondrial volume (fraction of hepatocyte)") # SBML compartment Mitochondria; Data S1 section 1
    vol_blood <- fixed(0.85); label("Relative blood compartment volume") # SBML compartment Blood

    # --- Constant inputs and co-substrates (Data S1 Table B) ------------
    metin <- fixed(100); label("Methionine-cycle input to homocysteine (uM/h)") # SBML METIN; Table B METIN 100
    cysin <- fixed(105); label("Cysteine input into blood (uM/h)") # SBML K_CysIntoBlood; Data S1 text CYSIN = 105
    cglu <- fixed(3219); label("Cytosolic glutamate, held constant (uM)") # SBML cglu; Table B cGlu 3219
    cgly <- fixed(1850); label("Cytosolic glycine, held constant (uM)") # SBML cgly; Table B cGly 1850
    cnadph <- fixed(50); label("Cytosolic NADPH, held constant (uM)") # SBML cNADPH; Table B c(m)NADPH 50
    mnadph <- fixed(50); label("Mitochondrial NADPH, held constant (uM)") # SBML mNADPH; Table B c(m)NADPH 50

    # --- Hormone reference levels (males / postmenopausal females) -------
    e2_ref <- fixed(0.09); label("Reference serum estradiol at which estradiol activation is zero (nM)") # Data S1 CBS, GCL, GPX equations (E2 - 0.09); Methods 0.09 nM
    p4_ref <- fixed(0.318); label("Reference serum progesterone at which progesterone inhibition is zero (nM)") # Data S1 GPX equations (P4 - 0.318); Methods 0.318 nM

    # --- Steady-state reference concentrations used in activation terms -
    ss_cgsh <- fixed(5000); label("Reference cytosolic GSH in the CBS GSH-activation term (uM)") # SBML sscGSH; Data S1 CBS table sscGSH 5000
    ss_ch2o2 <- fixed(0.01); label("Reference cytosolic H2O2 in the CBS, CTGL and GCL H2O2 terms (uM)") # SBML sscH2O2; Data S1 sscH2O2 0.01
    ss_hcy <- fixed(1.8); label("Reference homocysteine in the cytosolic GPX inhibition term (uM)") # SBML ssHcy; Data S1 GPXcyto equation (Ki + 1.8)

    # --- CBS (homocysteine -> cystathionine) ----------------------------
    f_cbs <- fixed(1); label("Fraction of healthy CBS activity (1 = healthy; 0.025 = CBS deficiency)") # Results 'CBS Deficiency': CBSD = 2.5% of normal CBS activity
    vmax_cbs <- fixed(7042.85714285714); label("CBS Vmax (uM/h)") # SBML Vmax_CBS; Data S1 CBS table Vmax 7043
    km_cbs <- fixed(170); label("CBS Km for homocysteine (uM)") # SBML Km_CBS; Data S1 CBS table Km 170
    ka_cbs_gsh <- fixed(6000); label("CBS activation constant for cytosolic GSH (uM)") # SBML Km_CBS_GSH; Data S1 CBS table Ka,GSH 6000
    ka_cbs_h2o2 <- fixed(0.01); label("CBS activation constant for cytosolic H2O2 (uM)") # SBML Km_CBS_h2o2; Data S1 CBS table Ka,H2O2 0.01
    ve2_cbs <- fixed(1); label("CBS maximal fold-activation by estradiol (unitless)") # Data S1 CBS table V_E2 1 (implicit 1 in SBML)
    ke2_cbs <- fixed(0.67); label("Estradiol excess at half-maximal CBS activation (nM)") # SBML K_E2_CBS; Data S1 CBS table K_E2 0.67

    # --- CTGL (cystathionine -> cysteine) -------------------------------
    vmax_ctgl <- fixed(1200); label("CTGL Vmax (uM/h)") # SBML Vmax_CTGL; Data S1 CTGL table Vmax 1200
    km_ctgl <- fixed(500); label("CTGL Km for cystathionine (uM)") # SBML Km_CTGL; Data S1 CTGL table Km 500
    ka_ctgl_h2o2 <- fixed(0.01); label("CTGL activation constant for cytosolic H2O2 (uM)") # SBML K_CTGL_H2O2; Data S1 CTGL table Ka 0.01

    # --- GCL (cysteine + glutamate -> glutamyl-cysteine) ----------------
    vmax_gcl <- fixed(3600); label("GCL Vmax (uM/h)") # SBML Vmax_GCL; Data S1 GCL table Vmax 3600
    km_gcl_ccys <- fixed(100); label("GCL Km for cytosolic cysteine (uM)") # SBML Km_GCL_cys; Data S1 GCL table 100
    km_gcl_cglu <- fixed(1900); label("GCL Km for glutamate (uM)") # SBML Km_GCL_cglu; Data S1 GCL table 1900
    ki_gcl_cgsh <- fixed(8200); label("GCL inhibition constant for cytosolic GSH (uM)") # SBML Ki_GCL_GSH; Data S1 GCL table Ki 8200
    ke_gcl <- fixed(5597); label("GCL equilibrium constant (1/uM)") # SBML Ke_GCL_glutcys; Data S1 GCL table Ke 5597
    kp_gcl <- fixed(300); label("GCL product constant for glutamyl-cysteine (1/uM)") # SBML Kp_GCL_glutcys; Data S1 GCL table Kp 300
    ka_gcl_h2o2 <- fixed(0.01); label("GCL activation constant for cytosolic H2O2 (uM)") # SBML Ka_GCL_H2O2; Data S1 GCL table Ka 0.01
    ve2_gcl <- fixed(25.7716); label("GCL maximal fold-activation by estradiol (unitless)") # SBML V_E2_GCL; Data S1 GCL table V_E2 25.77
    ke2_gcl <- fixed(10.8); label("Estradiol excess at half-maximal GCL activation (nM)") # SBML K_E2_GCL; Data S1 GCL table K_E2 10.8

    # --- GS (glutamyl-cysteine + glycine -> GSH) ------------------------
    vmax_gs <- fixed(5400); label("GS Vmax (uM/h)") # SBML Vmax_GS; Data S1 GS table Vmax 5400
    km_gs_glutcys <- fixed(22); label("GS Km for glutamyl-cysteine (uM)") # SBML Km_GS_glutcys; Data S1 GS table 22
    km_gs_cgly <- fixed(300); label("GS Km for glycine (uM)") # SBML Km_GS_cgly; Data S1 GS table 300
    ke_gs <- fixed(5600); label("GS equilibrium constant (1/uM)") # SBML Ke_GS_GSH; Data S1 GS table Ke 5600
    kp_gs <- fixed(30); label("GS product constant for GSH (1/uM)") # SBML Kp_GS_GSH; Data S1 GS table Kp 30

    # --- GPX (2 GSH + H2O2 -> GSSG), shared hormone regulation ----------
    ve2_gpx <- fixed(25.7716); label("GPX maximal fold-activation by estradiol (unitless)") # SBML V_E2_GPX; Table S2 nominal 25.7716 (Data S1 GPX tables print 27)
    ke2_gpx <- fixed(10.8); label("Estradiol excess at half-maximal GPX activation (nM)") # SBML K_E2_GPX; Table S2 nominal 10.8 (Data S1 GPX tables print 11)
    vp4_gpx <- fixed(0.3); label("GPX maximal fractional inhibition by progesterone (unitless)") # SBML V_P4_GPX; Data S1 GPX tables V_P4 0.3
    kp4_gpx <- fixed(70); label("Progesterone excess at half-maximal GPX inhibition (nM)") # SBML K_P4_GPX; Data S1 GPX tables K_P4 70
    vmax_gpx_cyto <- fixed(1.8e6); label("Cytosolic GPX Vmax (uM/h)") # SBML Vmax_GPXcyto; Data S1 GPXcyto table 1.8e6
    km_gpx_cyto_gsh <- fixed(1330); label("Cytosolic GPX Km for GSH (uM)") # SBML Km_GPXcyto_GSH; Data S1 GPXcyto table 1330
    km_gpx_cyto_h2o2 <- fixed(6.8); label("Cytosolic GPX Km for H2O2 (uM)") # SBML Km_GPXcyto_H2O2; Data S1 GPXcyto table 6.8
    ki_gpx_cyto_hcy <- fixed(27); label("Cytosolic GPX inhibition constant for homocysteine (uM)") # SBML Ki_GPXcyto_Hcy; Data S1 GPXcyto table Ki 27
    vmax_gpx_mito <- fixed(2.5e6); label("Mitochondrial GPX Vmax (uM/h)") # SBML Vmax_GPXmito; Data S1 GPXmito table 2.5e6
    km_gpx_mito_gsh <- fixed(1330); label("Mitochondrial GPX Km for GSH (uM)") # SBML Km_GPXmito_GSH; Data S1 GPXmito table 1330
    km_gpx_mito_h2o2 <- fixed(6.8); label("Mitochondrial GPX Km for H2O2 (uM)") # SBML Km_GPXmito_H2O2; Data S1 GPXmito table 6.8

    # --- GR (GSSG -> 2 GSH) ----------------------------------------------
    vmax_gr_cyto <- fixed(7000); label("Cytosolic GR Vmax (uM/h)") # SBML Vmax_GRcyto; Table C GRcyto Vmax 7000
    km_gr_cyto_gssg <- fixed(107); label("Cytosolic GR Km for GSSG (uM)") # SBML Km_GRcyto_GSSG; Table C 107
    km_gr_cyto_nadph <- fixed(10.4); label("Cytosolic GR Km for NADPH (uM)") # SBML Km_GRcyto_cNADPH; Table C 10.4
    vmax_gr_mito <- fixed(5300); label("Mitochondrial GR Vmax (uM/h)") # SBML Vmax_GRmito; Table C GRmito Vmax 5300
    km_gr_mito_gssg <- fixed(107); label("Mitochondrial GR Km for GSSG (uM)") # SBML Km_GRmito_GSSG; Table C 107
    km_gr_mito_nadph <- fixed(10.4); label("Mitochondrial GR Km for NADPH (uM)") # SBML Km_GRmito_mNADPH; Table C 10.4

    # --- TPX (H2O2 + reduced Trx -> oxidised Trx) ------------------------
    vmax_tpx_cyto <- fixed(3e5); label("Cytosolic TPX Vmax (uM/h)") # SBML Vmax_TPXcyto; Data S1 TPXcyto table 3e5
    km_tpx_cyto_trxred <- fixed(2.4); label("Cytosolic TPX Km for reduced thioredoxin (uM)") # SBML Km_TPXcyto_Trxr; Data S1 TPXcyto table 2.4
    km_tpx_cyto_h2o2 <- fixed(0.7); label("Cytosolic TPX Km for H2O2 (uM)") # SBML Km_TPXcyto_H2O2; Data S1 TPXcyto table 0.7
    ic50_tpx_cyto_h2o2 <- fixed(0.17); label("H2O2 at half-maximal cytosolic TPX overoxidation (uM)") # SBML IC50_TPXcyto_H2O2; Data S1 TPXcyto table IC50 0.17
    vmax_tpx_mito <- fixed(8e5); label("Mitochondrial TPX Vmax (uM/h)") # SBML Vmax_TPXmito; Data S1 TPXmito table 8e5
    km_tpx_mito_trxred <- fixed(4.3); label("Mitochondrial TPX Km for reduced thioredoxin (uM)") # SBML Km_TPXmito_Trxr; Data S1 TPXmito table 4.3
    km_tpx_mito_h2o2 <- fixed(0.7); label("Mitochondrial TPX Km for H2O2 (uM)") # SBML Km_TPXmito_H2O2; Data S1 TPXmito table 0.7
    ic50_tpx_mito_h2o2 <- fixed(1.8); label("H2O2 at half-maximal mitochondrial TPX overoxidation (uM)") # SBML IC50_TPXmito_H2O2; Data S1 TPXmito table IC50 1.8

    # --- TR (oxidised Trx -> reduced Trx) --------------------------------
    vmax_tr_cyto <- fixed(40000); label("Cytosolic TR Vmax (uM/h)") # SBML Vmax_TRcyto; Table C TRcyto Vmax 40000
    km_tr_cyto_trxox <- fixed(2.5); label("Cytosolic TR Km for oxidised thioredoxin (uM)") # SBML Km_TRcyto_Trxo; Table C 2.5
    km_tr_cyto_nadph <- fixed(6); label("Cytosolic TR Km for NADPH (uM)") # SBML Km_TRcyto_cNADPH; Table C 6
    vmax_tr_mito <- fixed(50000); label("Mitochondrial TR Vmax (uM/h)") # SBML Vmax_TRmito; Table C TRmito Vmax 50000
    km_tr_mito_trxox <- fixed(2.5); label("Mitochondrial TR Km for oxidised thioredoxin (uM)") # SBML Km_TRmito_Trxo; Table C 2.5
    km_tr_mito_nadph <- fixed(6); label("Mitochondrial TR Km for NADPH (uM)") # SBML Km_TRmito_mNADPH; Table C 6

    # --- CAT (H2O2 removal) ----------------------------------------------
    vmax_cat_cyto <- fixed(6.5e8); label("Cytosolic catalase Vmax (uM/h)") # SBML Vmax_CATcyto; Data S1 CATcyto table 6.5e8
    km_cat_cyto_h2o2 <- fixed(32300); label("Cytosolic catalase Km for H2O2 (uM)") # SBML Km_CATcyto_H2O2; Data S1 CATcyto table 32300
    ic50_cat_cyto_hcy <- fixed(25); label("Homocysteine at half-maximal cytosolic catalase inhibition (uM)") # SBML IC50_CATcyto_Hcy; Data S1 CATcyto table IC50 25
    vmax_cat_mito <- fixed(5.5e10); label("Mitochondrial catalase Vmax (uM/h)") # SBML Vmax_CATmito; Table C CATmito Vmax 5.5e10
    km_cat_mito_h2o2 <- fixed(32300); label("Mitochondrial catalase Km for H2O2 (uM)") # SBML Km_CATmito_H2O2; Table C 32300

    # --- H2O2 production and mitochondrial exchange ----------------------
    h2o2prod_cyto <- fixed(3026); label("Cytosolic H2O2 production (uM/h)") # SBML H2O2prod_cyto; Data S1 Hydrogen Peroxide Production 3026
    h2o2prod_mito <- fixed(12000); label("Mitochondrial H2O2 production (uM/h)") # SBML H2O2prod_mito; Data S1 Hydrogen Peroxide Production 12000
    k_h2o2_transport <- fixed(120); label("H2O2 diffusion rate across the mitochondrial membrane (1/h)") # SBML H2O2Transport; Data S1 H2O2 Transport table 120

    # --- GSH exchange between cytosol and mitochondria --------------------
    vmax_gsh_tomito_h <- fixed(388); label("High-affinity cytosol-to-mitochondria GSH transporter Vmax (uM/h)") # SBML Vmax_GSHtoMitoH; Data S1 GSH Transport table 388
    km_gsh_tomito_h <- fixed(36.4286); label("High-affinity cytosol-to-mitochondria GSH transporter Km (uM)") # SBML Km_GSHtoMitoH; Data S1 GSH Transport table 36.4
    vmax_gsh_tomito_l <- fixed(2711.42857142857); label("Low-affinity cytosol-to-mitochondria GSH transporter Vmax (uM/h)") # SBML Vmax_GSHtoMitoL; Data S1 GSH Transport table 2711
    km_gsh_tomito_l <- fixed(9000); label("Low-affinity cytosol-to-mitochondria GSH transporter Km (uM)") # SBML Km_GSHtoMitoL; Data S1 GSH Transport table 9000
    k_gsh_tocyto <- fixed(1.03); label("Mitochondria-to-cytosol GSH transport rate (1/h)") # SBML K_GSHtoCyto; Data S1 GSH Transport table K_GSH 1.03

    # --- Cytosol <-> blood transport -------------------------------------
    vmax_gsh_toblood_h <- fixed(330); label("High-affinity GSH export Vmax (uM/h)") # SBML Vmax_GSHtoBloodH; Data S1 GSH/GSSG Transport table 330
    km_gsh_toblood_h <- fixed(150); label("High-affinity GSH export Km (uM)") # SBML Km_GSHtoBloodH; Data S1 GSH/GSSG Transport table 150
    vmax_gsh_toblood_l <- fixed(1100); label("Low-affinity GSH export Vmax (uM/h)") # SBML Vmax_GSHtoBloodL; Data S1 GSH/GSSG Transport table 1100
    km_gsh_toblood_l <- fixed(3000); label("Low-affinity GSH export half-saturation (uM)") # SBML Km_GSHtoBloodL; Data S1 GSH/GSSG Transport table 3000
    hill_gsh_toblood_l <- fixed(3); label("Low-affinity GSH export Hill coefficient (unitless)") # SBML Hill_GSHtoBloodL; Data S1 GSH/GSSG Transport table Hill 3
    vmax_gssg_toblood_h <- fixed(40); label("High-affinity GSSG export Vmax (uM/h)") # SBML Vmax_GSSGtoBloodH; Data S1 GSH/GSSG Transport table 40
    km_gssg_toblood_h <- fixed(1250); label("High-affinity GSSG export Km (uM)") # SBML Km_GSSGtoBloodH; Data S1 GSH/GSSG Transport table 1250
    vmax_gssg_toblood_l <- fixed(4025); label("Low-affinity GSSG export Vmax (uM/h)") # SBML Vmax_GSSGtoBloodL; Data S1 GSH/GSSG Transport table 4025
    km_gssg_toblood_l <- fixed(7100); label("Low-affinity GSSG export Km (uM)") # SBML Km_GSSGtoBloodL; Data S1 GSH/GSSG Transport table 7100
    vmax_cys_fromblood <- fixed(14950); label("Blood-to-cytosol cysteine uptake Vmax (uM/h)") # SBML Vmax_CysfromBlood; Data S1 Cysteine Transport table 14950
    km_cys_fromblood <- fixed(2100); label("Blood-to-cytosol cysteine uptake Km (uM)") # SBML Km_CysfromBlood; Data S1 Cysteine Transport table 2100
    k_cys_toblood <- fixed(1); label("Cytosol-to-blood cysteine efflux rate (1/h)") # SBML K_CystoBlood; Data S1 Cysteine Transport table K_out 1

    # --- First-order removal and breakdown --------------------------------
    k_hcy_removal <- fixed(3); label("Homocysteine removal (remethylation and export) rate (1/h)") # SBML K_HcyRemoval; Data S1 Hcy equation (3)(Hcy)
    k_ccys_removal <- fixed(0.0017); label("Second-order cytosolic cysteine loss to other pathways (1/(uM*h))") # SBML K_cCysRemoval; Data S1 cCys equation prints (0.35) cCys^2 / 200
    k_cgsh_removal <- fixed(0.002); label("Cytosolic GSH removal rate (1/h)") # SBML K_cGSHRemoval; Data S1 cGSH equation (0.002)
    k_cgssg_removal <- fixed(0.1); label("Cytosolic GSSG removal rate (1/h)") # SBML K_cGSSGRemoval; Data S1 cGSSG equation (.1)
    k_mgssg_removal <- fixed(1); label("Mitochondrial GSSG removal rate (1/h)") # SBML K_mGSSGRemoval; Data S1 mGSSG equation (1)
    k_bgsh_breakdown <- fixed(90); label("Blood GSH breakdown to cysteine rate (1/h)") # SBML K_GSHBreakdownBlood; Data S1 bGSH equation (90)
    k_bgssg_breakdown <- fixed(67.5); label("Blood GSSG breakdown to cysteine rate (1/h)") # SBML K_GSSGBreakdownBlood; Data S1 bGSSG equation (67.5)
    k_bgsh_removal <- fixed(0.7); label("Blood GSH removal rate (1/h)") # SBML K_bGSHRemoval; Data S1 bGSH equation (0.7)
    k_bgssg_removal <- fixed(7.5); label("Blood GSSG removal rate (1/h)") # SBML K_bGSSGRemoval; Data S1 bGSSG equation (7.5)
    k_bcys_removal <- fixed(0.35); label("Blood cysteine removal rate (1/h)") # SBML K_bCysRemoval; Data S1 bCys equation (0.35)

    # --- Thioredoxin synthesis and removal ----------------------------------
    trxred_prod_cyto <- fixed(0.1); label("Cytosolic reduced thioredoxin synthesis (uM/h)") # SBML Trxr_cytoProd; Data S1 Thioredoxin Synthesis 0.1
    trxred_prod_mito <- fixed(0.1); label("Mitochondrial reduced thioredoxin synthesis (uM/h)") # SBML Trxr_mitoProd; Data S1 Thioredoxin Synthesis 0.1
    k_trxred_removal_cyto <- fixed(0.08); label("Cytosolic reduced thioredoxin removal rate (1/h)") # SBML K_Trxr_cytoRemoval; Data S1 Thioredoxin Synthesis 0.08
    k_trxox_removal_cyto <- fixed(0.05); label("Cytosolic oxidised thioredoxin removal rate (1/h)") # SBML K_Trxo_cytoRemoval; Data S1 Thioredoxin Synthesis 0.05
    k_trxred_removal_mito <- fixed(0.08); label("Mitochondrial reduced thioredoxin removal rate (1/h)") # SBML K_Trxr_mitoRemoval; Data S1 Thioredoxin Synthesis 0.08
    k_trxox_removal_mito <- fixed(0.05); label("Mitochondrial oxidised thioredoxin removal rate (1/h)") # SBML K_Trxo_mitoRemoval; Data S1 Thioredoxin Synthesis 0.05

    # --- Initial conditions (SBML initial concentrations) -----------------
    bl_hcy <- fixed(2.33); label("Initial cytosolic homocysteine (uM)") # SBML species Hcy
    bl_cyst <- fixed(42.67); label("Initial cystathionine (uM)") # SBML species Cyst
    bl_ccys <- fixed(198.9); label("Initial cytosolic cysteine (uM)") # SBML species Cys_cytosol
    bl_glutcys <- fixed(8.52); label("Initial glutamyl-cysteine (uM)") # SBML species Glutcys
    bl_cgsh <- fixed(4889); label("Initial cytosolic GSH (uM)") # SBML species GSH_cytosol
    bl_cgssg <- fixed(38.42); label("Initial cytosolic GSSG (uM)") # SBML species GSSG_cytosol
    bl_bgsh <- fixed(13.38); label("Initial blood GSH (uM)") # SBML species GSH_blood
    bl_bgssg <- fixed(0.30478); label("Initial blood GSSG (uM)") # SBML species GSSG_blood
    bl_bcys <- fixed(228.89); label("Initial blood cysteine (uM)") # SBML species Cys_blood
    bl_ch2o2 <- fixed(0.00971); label("Initial cytosolic H2O2 (uM)") # SBML species H2O2_cytosol
    bl_mh2o2 <- fixed(0.0054); label("Initial mitochondrial H2O2 (uM)") # SBML species H2O2_mitochondria
    bl_mgsh <- fixed(7274); label("Initial mitochondrial GSH (uM)") # SBML species GSH_mitochondria
    bl_mgssg <- fixed(49); label("Initial mitochondrial GSSG (uM)") # SBML species GSSG_mitochondria
    bl_ctrxred <- fixed(1.19); label("Initial cytosolic reduced thioredoxin (uM)") # SBML species Trxr_cytosol
    bl_ctrxox <- fixed(0.09); label("Initial cytosolic oxidised thioredoxin (uM)") # SBML species Trxo_cytosol
    bl_mtrxred <- fixed(1.2); label("Initial mitochondrial reduced thioredoxin (uM)") # SBML species Trxr_mitochondria
    bl_mtrxox <- fixed(0.07); label("Initial mitochondrial oxidised thioredoxin (uM)") # SBML species Trxo_mitochondria
  })

  model({
    # Hormone excess over the male / postmenopausal reference levels (nM).
    de2 <- ESTRADIOL - e2_ref
    dp4 <- PROGESTERONE - p4_ref

    # Velocities are in uM/h of the compartment where the reaction runs
    # (SBML kinetic law divided by that compartment's volume).

    # CBS: Data S1 CBS equation; SBML Reaction_2.
    v_cbs <- f_cbs * vmax_cbs * hcy / (km_cbs + hcy) *
      (1 + ve2_cbs * de2 / (ke2_cbs + de2)) *
      (ka_cbs_gsh + cgsh) / (ka_cbs_gsh + ss_cgsh) *
      (ka_cbs_h2o2 + ch2o2) / (ka_cbs_h2o2 + ss_ch2o2)

    # CTGL: Data S1 CTGL equation; SBML Reaction_3.
    v_ctgl <- vmax_ctgl * cyst / (km_ctgl + cyst) *
      (ka_ctgl_h2o2 + ch2o2) / (ka_ctgl_h2o2 + ss_ch2o2)

    # GCL: Data S1 GCL equation (empty kinetic law in the SBML, Reaction_4).
    # The H2O2 factor is the activation form (Ka + cH2O2)/(Ka + sscH2O2),
    # as in CBS and CTGL; Data S1 prints the reciprocal, which does not
    # reproduce the Figure 2 steady states (see the vignette).
    v_gcl <- vmax_gcl * (ccys * cglu - glutcys / ke_gcl) /
      (km_gcl_ccys * km_gcl_cglu + km_gcl_ccys * cglu +
        km_gcl_cglu * ccys * (1 + cgsh / ki_gcl_cgsh + cglu / km_gcl_cglu) +
        glutcys / kp_gcl + cgsh / ki_gcl_cgsh) *
      (ka_gcl_h2o2 + ch2o2) / (ka_gcl_h2o2 + ss_ch2o2) *
      (1 + ve2_gcl * de2 / (ke2_gcl + de2))

    # GS: SBML Reaction_5 (Data S1 prints a different denominator).
    v_gs <- vmax_gs * (cgly * glutcys - cgsh / ke_gs) /
      (km_gs_glutcys * km_gs_cgly + km_gs_cgly * glutcys +
        km_gs_glutcys * cgly * (1 + glutcys / km_gs_glutcys) + cgsh / kp_gs)

    # GPX: Data S1 GPXcyto and GPXmito equations (empty in the SBML,
    # Reactions 6 and 17). Velocity in uM H2O2/h; 2 GSH consumed per H2O2.
    f_p4_gpx <- 1 - vp4_gpx * dp4 / (kp4_gpx + dp4)
    f_e2_gpx <- 1 + ve2_gpx * de2 / (ke2_gpx + de2)
    v_gpx_cyto <- vmax_gpx_cyto * cgsh^2 * ch2o2 /
      ((km_gpx_cyto_gsh + cgsh)^2 * (km_gpx_cyto_h2o2 + ch2o2)) *
      f_p4_gpx * f_e2_gpx * (ki_gpx_cyto_hcy + ss_hcy) / (ki_gpx_cyto_hcy + hcy)
    v_gpx_mito <- vmax_gpx_mito * mgsh^2 * mh2o2 /
      ((km_gpx_mito_gsh + mgsh)^2 * (km_gpx_mito_h2o2 + mh2o2)) *
      f_p4_gpx * f_e2_gpx

    # GR: Table C bidirectional-form constants; SBML Reactions 7 and 18.
    gr_cyto_s1 <- cgssg / km_gr_cyto_gssg
    gr_cyto_s2 <- cnadph / km_gr_cyto_nadph
    v_gr_cyto <- vmax_gr_cyto * gr_cyto_s1 * gr_cyto_s2 / (1 + gr_cyto_s1 + gr_cyto_s2 + gr_cyto_s1 * gr_cyto_s2)
    gr_mito_s1 <- mgssg / km_gr_mito_gssg
    gr_mito_s2 <- mnadph / km_gr_mito_nadph
    v_gr_mito <- vmax_gr_mito * gr_mito_s1 * gr_mito_s2 / (1 + gr_mito_s1 + gr_mito_s2 + gr_mito_s1 * gr_mito_s2)

    # TPX with H2O2 overoxidation: Data S1 TPXcyto / TPXmito equations
    # (empty in the SBML, Reactions 8 and 19). The mitochondrial
    # overoxidation term uses mitochondrial H2O2 (Data S1 writes cH2O2 there
    # but lists V_TPXmito(mtrx_red, mH2O2) in the ODEs).
    v_tpx_cyto <- vmax_tpx_cyto * ctrxred * ch2o2 /
      ((km_tpx_cyto_trxred + ctrxred) * (km_tpx_cyto_h2o2 + ch2o2)) *
      (1 - ch2o2 / (ic50_tpx_cyto_h2o2 + ch2o2))
    v_tpx_mito <- vmax_tpx_mito * mtrxred * mh2o2 /
      ((km_tpx_mito_trxred + mtrxred) * (km_tpx_mito_h2o2 + mh2o2)) *
      (1 - mh2o2 / (ic50_tpx_mito_h2o2 + mh2o2))

    # TR: SBML Reactions 9 and 20.
    v_tr_cyto <- vmax_tr_cyto * ctrxox / (km_tr_cyto_trxox + ctrxox) * cnadph / (km_tr_cyto_nadph + cnadph)
    v_tr_mito <- vmax_tr_mito * mtrxox / (km_tr_mito_trxox + mtrxox) * mnadph / (km_tr_mito_nadph + mnadph)

    # CAT: SBML Reactions 10 and 21; cytosolic CAT is inhibited by Hcy.
    v_cat_cyto <- vmax_cat_cyto * ch2o2 / (km_cat_cyto_h2o2 + ch2o2) * (1 - hcy / (ic50_cat_cyto_hcy + hcy))
    v_cat_mito <- vmax_cat_mito * mh2o2 / (km_cat_mito_h2o2 + mh2o2)

    # Transport (SBML Reactions 13-16, 33-35, 42-43).
    v_h2o2_tomito <- k_h2o2_transport * ch2o2 # uM/h of cytosol
    v_h2o2_tocyto <- k_h2o2_transport * mh2o2 # uM/h of mitochondria
    v_gsh_tomito <- vmax_gsh_tomito_h * cgsh / (cgsh + km_gsh_tomito_h) +
      vmax_gsh_tomito_l * cgsh / (cgsh + km_gsh_tomito_l) # uM/h of cytosol
    v_gsh_tocyto <- k_gsh_tocyto * mgsh # uM/h of mitochondria
    v_gsh_toblood <- vmax_gsh_toblood_h * cgsh / (km_gsh_toblood_h + cgsh) +
      vmax_gsh_toblood_l * cgsh^hill_gsh_toblood_l / (km_gsh_toblood_l^hill_gsh_toblood_l + cgsh^hill_gsh_toblood_l) # uM/h of cytosol
    v_gssg_toblood <- vmax_gssg_toblood_h * cgssg / (km_gssg_toblood_h + cgssg) +
      vmax_gssg_toblood_l * cgssg / (km_gssg_toblood_l + cgssg) # uM/h of cytosol
    v_cys_fromblood <- vmax_cys_fromblood * bcys / (km_cys_fromblood + bcys) - k_cys_toblood * ccys # uM/h of blood

    r_cm <- vol_cytosol / vol_mitochondria
    r_mc <- vol_mitochondria / vol_cytosol
    r_cb <- vol_cytosol / vol_blood
    r_bc <- vol_blood / vol_cytosol

    # Mass balances (Data S1 section 1 differential equations).
    d/dt(hcy) <- metin - v_cbs - k_hcy_removal * hcy
    d/dt(cyst) <- v_cbs - v_ctgl
    d/dt(ccys) <- v_ctgl - v_gcl + r_bc * v_cys_fromblood - k_ccys_removal * ccys^2
    d/dt(glutcys) <- v_gcl - v_gs
    d/dt(cgsh) <- v_gs - v_gsh_toblood - 2 * v_gpx_cyto + 2 * v_gr_cyto -
      k_cgsh_removal * cgsh - v_gsh_tomito + r_mc * v_gsh_tocyto
    d/dt(cgssg) <- v_gpx_cyto - v_gr_cyto - v_gssg_toblood - k_cgssg_removal * cgssg
    d/dt(bgsh) <- r_cb * v_gsh_toblood - k_bgsh_removal * bgsh - k_bgsh_breakdown * bgsh
    d/dt(bgssg) <- r_cb * v_gssg_toblood - k_bgssg_removal * bgssg - k_bgssg_breakdown * bgssg
    d/dt(bcys) <- k_bgsh_breakdown * bgsh + 2 * k_bgssg_breakdown * bgssg + cysin -
      v_cys_fromblood - k_bcys_removal * bcys
    d/dt(ch2o2) <- h2o2prod_cyto + r_mc * v_h2o2_tocyto - v_h2o2_tomito -
      v_gpx_cyto - v_cat_cyto - v_tpx_cyto
    d/dt(mh2o2) <- h2o2prod_mito - v_h2o2_tocyto + r_cm * v_h2o2_tomito -
      v_gpx_mito - v_cat_mito - v_tpx_mito
    d/dt(mgsh) <- -2 * v_gpx_mito + 2 * v_gr_mito - v_gsh_tocyto + r_cm * v_gsh_tomito
    d/dt(mgssg) <- v_gpx_mito - v_gr_mito - k_mgssg_removal * mgssg
    d/dt(ctrxred) <- trxred_prod_cyto - v_tpx_cyto + v_tr_cyto - k_trxred_removal_cyto * ctrxred
    d/dt(ctrxox) <- v_tpx_cyto - v_tr_cyto - k_trxox_removal_cyto * ctrxox
    d/dt(mtrxred) <- trxred_prod_mito - v_tpx_mito + v_tr_mito - k_trxred_removal_mito * mtrxred
    d/dt(mtrxox) <- v_tpx_mito - v_tr_mito - k_trxox_removal_mito * mtrxox

    hcy(0) <- bl_hcy
    cyst(0) <- bl_cyst
    ccys(0) <- bl_ccys
    glutcys(0) <- bl_glutcys
    cgsh(0) <- bl_cgsh
    cgssg(0) <- bl_cgssg
    bgsh(0) <- bl_bgsh
    bgssg(0) <- bl_bgssg
    bcys(0) <- bl_bcys
    ch2o2(0) <- bl_ch2o2
    mh2o2(0) <- bl_mh2o2
    mgsh(0) <- bl_mgsh
    mgssg(0) <- bl_mgssg
    ctrxred(0) <- bl_ctrxred
    ctrxox(0) <- bl_ctrxox
    mtrxred(0) <- bl_mtrxred
    mtrxox(0) <- bl_mtrxox

    # Cytosolic GSH/GSSG ratio, the redox index plotted in Figures 3D and 5F.
    cgsh_cgssg_ratio <- cgsh / cgssg
  })
}
