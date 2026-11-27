Sluka_2016_acetaminophen_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, perfusion-limited, SBML). The standalone whole-body",
    "acetaminophen (APAP) module of the Sluka 2016 liver-centric multiscale",
    "framework, run on its own for a 1.4 g oral dose in a 70 kg man (paper",
    "Fig 4). Gut lumen, gut, liver, kidney, lung, arterial blood, venous",
    "blood and a lumped rest-of-body compartment, translated by the authors",
    "from the Wambaugh and Shah 2010 PBPK model into SBML. First-order",
    "absorption from the gut lumen into the gut tissue, which drains with",
    "the hepatic artery into the liver; hepatic metabolism is a single",
    "first-order removal of liver drug and renal excretion is filtration",
    "into the kidney tubules at a reduced glomerular filtration rate that",
    "stands in for APAP tubular reabsorption. Organ volumes and blood flows",
    "scale allometrically from body weight. Deterministic typical-value",
    "model: the paper reports no between-subject variability and no",
    "residual-error model. The coupled multiscale model (agent-based",
    "sinusoid in CompuCell3D linking this module to per-hepatocyte",
    "metabolism) is not an ODE system and is not included; the",
    "per-hepatocyte reaction network is Sluka_2016_acetaminophen_hepatocyte_qsp.",
    sep = " "
  )
  reference <- paste(
    "Sluka JP, Fu X, Swat M, Belmonte JM, Cosmanescu A, Clendenon SG,",
    "Wambaugh JF, Glazier JA. A Liver-Centric Multiscale Modeling Framework",
    "for Xenobiotics. PLoS One. 2016;11(9):e0162428.",
    "doi:10.1371/journal.pone.0162428. PMID 27636091. PMCID PMC5026379.",
    "Structure from Table 1 and the deposited SBML file",
    "Wambaugh_vLiver_PBPK_based_v4.xml (S2 File, BioModels MODEL1509230010);",
    "parameter values from Table 2 (parameter set REFSIM), S1 File Table 2",
    "(allometric scaling factors) and the deposited SBML.",
    sep = " "
  )
  vignette <- "Sluka_2016_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The oral dose is the initial amount of the gut lumen (SBML species
  # AGutlumen); the lumen is modelled as the depot.
  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "acetaminophen", units = "mg", specimen = "administration site", verified = TRUE),
    gut = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "acetaminophen", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "acetaminophen", units = "mg", specimen = "whole blood", verified = TRUE),
    a_urine = list(analyte = "acetaminophen", units = "mg", specimen = "urine", verified = TRUE),
    a_metabolized = list(
      analyte = "acetaminophen metabolites (cumulative amount of parent metabolized)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # No covariates: one 70 kg reference man. Body weight is a fixed system
  # parameter, as in the deposited SBML, rather than a data column.
  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 9L,
    n_studies = 1L,
    age_range = "23-44 years (calibration data)",
    weight_range = "55-97 kg (calibration data); simulated reference 70 kg",
    sex_female_pct = 33,
    race_ethnicity = c(White = 100),
    disease_state = "Healthy adult volunteers (calibration data); the model simulates a single 70 kg male.",
    dose_range = "1.4 g (20 mg/kg) single oral dose of APAP syrup.",
    notes = paste(
      "No population model was fitted. The parameter set REFSIM was hand-tuned",
      "(Model Parameterization section; S1 File section 6.1) so that the complete",
      "multiscale model reproduces mean plasma APAP, APAP-glucuronide and",
      "APAP-sulfate concentrations from the nine Caucasian subjects of Critchley",
      "et al. (J Clin Pharm Ther 2005;30(2):179-84, doi:10.1111/j.1365-2710.2004.00626.x;",
      "paper reference 58). Ages, weights and the 6:3 male:female ratio are from",
      "the Results 'Population Variability' section; the digitised data are S1",
      "File Table 15. The standalone PBPK module packaged here is the one shown",
      "against the APAP data in Fig 4.",
      sep = " "
    )
  )

  ini({
    # ============ Compound-independent (system) parameters ============
    bw <- fixed(70); label("Body weight (kg)") # Table 2 'pbpk_bw' 70 Kg; SBML parameter BW = 70

    # Cardiac output scales as bw^0.75 (S1 File section 3.1.1).
    qcc <- fixed(15); label("Cardiac output per kg^0.75 body weight (L/h/kg^0.75)") # S1 File Table 2 'QArt' factor 15.0 per Kg 3/4; SBML Cardiac_flow_ref = 15

    fq_gut <- fixed(0.205); label("Gut blood flow as fraction of cardiac output") # S1 File Table 2 'QGut' 0.2050 F of QArt; SBML QGut_fraction_QCardiac
    fq_liver <- fixed(0.0535); label("Hepatic arterial blood flow as fraction of cardiac output") # S1 File Table 2 'QLiver' 0.0535 F of QArt; SBML QLiver_fraction_QCardiac
    fq_kidney <- fixed(0.2214); label("Kidney blood flow as fraction of cardiac output") # S1 File Table 2 'QKidney' 0.2214 F of QArt; SBML QKidney_fraction_QCardiac

    # Total perfusable volume scales linearly with body weight in the
    # deposited SBML (VTotal = VTotal_ref * BW); see the vignette for the
    # bw^0.75 form used by the multiscale driver script.
    fvol_total <- fixed(0.6); label("Total perfusable volume per kg body weight (L/kg)") # SBML VTotal_ref = 0.6 with VTotal = VTotal_ref * BW; 0.6 * 70 = 42.0 L = S1 File Table 2 'VTotal' 42.0
    fvol_arterial <- fixed(0.0357); label("Arterial blood volume as fraction of total perfusable volume") # S1 File Table 2 'VArt' 0.0357 F of VTotal
    fvol_gut <- fixed(0.0263); label("Gut volume as fraction of total perfusable volume") # S1 File Table 2 'VGut' 0.0263 F of VTotal
    fvol_kidney <- fixed(0.007); label("Kidney volume as fraction of total perfusable volume") # S1 File Table 2 'VKidney' 0.0070 F of VTotal
    fvol_liver <- fixed(0.0408); label("Liver volume as fraction of total perfusable volume") # S1 File Table 2 'VLiver' 0.0408 F of VTotal
    fvol_lung <- fixed(0.0121); label("Lung volume as fraction of total perfusable volume") # S1 File Table 2 'VLung' 0.0121 F of VTotal
    fvol_venous <- fixed(0.0812); label("Venous blood volume as fraction of total perfusable volume") # S1 File Table 2 'VVen' 0.0812 F of VTotal

    # ============ Compound-dependent (APAP) parameters, set REFSIM ============
    lka <- fixed(log(1.5)); label("First-order absorption rate constant from the gut lumen (1/h)") # Table 2 'pbpk_kGutabs' 1.5 h-1; S1 File Table 10 REFSIM
    fup <- fixed(0.8); label("Fraction of APAP unbound in plasma") # Table 2 'pbpk_Fup' 0.8; S1 File Table 10 REFSIM
    bp <- fixed(1.09); label("Blood-to-plasma concentration ratio of APAP") # Table 2 'pbpk_Rb2p' 1.09; S1 File Table 10 REFSIM
    lkp_liver <- fixed(log(1)); label("Liver-to-plasma partition coefficient of APAP (unitless)") # Table 2 'pbpk_Kl2p' 1 (standalone value; replaced by the sinusoid model in the multiscale model)
    lkp_kidney <- fixed(log(1)); label("Kidney-to-plasma partition coefficient of APAP (unitless)") # Table 2 'pbpk_Kk2p' 1; S1 File Table 10 REFSIM
    lkp_other <- fixed(log(1.6)); label("Rest-of-body-to-plasma partition coefficient of APAP (unitless)") # Table 2 'pbpk_Kr2p' 1.6; S1 File Table 10 REFSIM
    qgfr <- fixed(0.714); label("Effective glomerular filtration rate of APAP, net of tubular reabsorption (L/h)") # Table 2 'pbpk_Qgfr' 0.714 L/h; S1 File Table 10 REFSIM 7.14E-1
    kmet <- fixed(9.5); label("First-order hepatic metabolism rate constant applied to liver amount / (kp_liver * fup) (1/h)") # deposited SBML parameter CLmetabolism = 9.5 (units first_order_rate_constant, per hour); Table 1 'CLiver -> CMetabolized' rate law
  })

  model({
    # ================= Allometric physiology =================
    # S1 File section 3.1.1 / SBML assignment rules.
    qc <- qcc * bw^0.75
    q_gut <- qc * fq_gut
    q_liver <- qc * fq_liver
    q_kidney <- qc * fq_kidney
    q_other <- qc - q_gut - q_kidney - q_liver # SBML QRest = QCardiac - QGut - QKidney - QLiver

    v_total <- fvol_total * bw
    v_arterial <- v_total * fvol_arterial
    v_gut <- v_total * fvol_gut
    v_kidney <- v_total * fvol_kidney
    v_liver <- v_total * fvol_liver
    v_lung <- v_total * fvol_lung
    v_venous <- v_total * fvol_venous
    v_other <- v_total - v_arterial - v_gut - v_kidney - v_liver - v_lung - v_venous # SBML VRest initial assignment

    ka <- exp(lka)
    kp_liver <- exp(lkp_liver)
    kp_kidney <- exp(lkp_kidney)
    kp_other <- exp(lkp_other)

    # ================= Transfer rates (amount/h), Table 1 =================
    # Every state holds an amount; Table 1 writes each rate as flow * amount
    # / volume, so amount / volume is the compartment concentration.
    r_absorb <- ka * depot # AGutlumen -> CGut
    r_art_gut <- q_gut * arterial / v_arterial # CArt -> CGut
    r_art_liver <- q_liver * arterial / v_arterial # CArt -> CLiver
    r_art_kidney <- q_kidney * arterial / v_arterial # CArt -> CKidney
    r_art_other <- q_other * arterial / v_arterial # CArt -> CRest
    r_gut_liver <- q_gut * gut / v_gut # CGut -> CLiver
    r_liver_ven <- (q_liver + q_gut) * liver * bp / (kp_liver * fup * v_liver) # CLiver -> CVen
    r_metab <- kmet * liver / (kp_liver * fup) # CLiver -> CMetabolized
    r_kidney_ven <- q_kidney * kidney * bp / (kp_kidney * fup * v_kidney) # CKidney -> CVen
    r_filtration <- qgfr * kidney / (kp_kidney * v_kidney) # CKidney -> CTubules
    r_other_ven <- q_other * other * bp / (kp_other * fup * v_other) # CRest -> CVen
    r_ven_lung <- qc * venous / v_venous # CVen -> CLung
    r_lung_art <- qc * lung / v_lung # CLung -> CArt

    d/dt(depot) <- -r_absorb
    d/dt(gut) <- r_absorb + r_art_gut - r_gut_liver
    d/dt(liver) <- r_gut_liver + r_art_liver - r_liver_ven - r_metab
    d/dt(kidney) <- r_art_kidney - r_kidney_ven - r_filtration
    d/dt(other) <- r_art_other - r_other_ven
    d/dt(venous) <- r_liver_ven + r_kidney_ven + r_other_ven - r_ven_lung
    d/dt(lung) <- r_ven_lung - r_lung_art
    d/dt(arterial) <- r_lung_art - r_art_gut - r_art_liver - r_art_kidney - r_art_other
    d/dt(a_urine) <- r_filtration
    d/dt(a_metabolized) <- r_metab

    # ================= Observations =================
    # Venous whole-blood concentration, the quantity plotted in Fig 4.
    Cvenous <- venous / v_venous
    # Venous plasma concentration, computed as in the authors' multiscale
    # driver script (CVen / VVen / Rb2p) and compared there with the plasma
    # data of Critchley et al.
    Cc <- Cvenous / bp
  })
}
