Yang_2026_deha_rat_pbpk <- function() {
  description <- paste0(
    "Preclinical (rat). PBPK (whole-body, 11 tissue compartments plus ",
    "arterial and venous blood pools and a lymph compartment, parent plus ",
    "metabolite, coded in Berkeley Madonna). Di(2-ethylhexyl) adipate (DEHA), ",
    "a lipophilic ester plasticizer and candidate pharmaceutical excipient, ",
    "and its primary hydrolysis metabolite mono(2-ethylhexyl) adipate (MEHA) ",
    "in Sprague-Dawley rats after intravenous and oral administration (Yang ",
    "et al. 2026, Pharmaceutics). Distribution is capillary-permeability ",
    "limited (fractional-distribution parameter fd per tissue from a PAMPA ",
    "permeability-surface product). Both adipates undergo rapid multi-tissue ",
    "hydrolysis: carboxylesterase-mediated parenchymal clearance (the ",
    "BNPP-inhibited fraction) plus, for DEHA only, a BNPP-insensitive ",
    "vascular (lipase-mediated) clearance acting in arterial and venous ",
    "plasma; hepatic clearance is scaled from liver-microsome stability. DEHA ",
    "metabolism in every compartment forms MEHA (complete conversion, fm=1). ",
    "Oral DEHA absorption uses a Weibull-dissolution GI-transit model with ",
    "parallel portal-hepatic and intestinal-lymphatic routes; lymphatic DEHA ",
    "returns to venous blood after a thoracic-duct lag and is hydrolysed to ",
    "MEHA by saturable (Michaelis-Menten) kinetics within the lymph. The ODE ",
    "system is transcribed from the authors' deposited Berkeley Madonna code; ",
    "where the printed tables disagree with the code, the code governs (see ",
    "the vignette Errata). Deterministic: no inter-individual variability or ",
    "residual-error magnitude is reported, so the model is for typical-value ",
    "simulation. The lymphatic return uses delay(), so every rxSolve needs ",
    "method = 'dop853'."
  )
  reference <- paste0(
    "Yang E, Jeong YS, Kim M, Kim S, Chung SJ. Physiologically Based ",
    "Pharmacokinetic Modeling of Di(2-ethylhexyl) Adipate and Its Primary ",
    "Metabolite, Mono(2-ethylhexyl) Adipate, in Rats, Incorporating ",
    "Circulatory Topology-Based Multi-Tissue Metabolism and Lymphatic ",
    "Absorption. Pharmaceutics. 2026;18(8):1001. ",
    "doi:10.3390/pharmaceutics18081001. ",
    "Structural parameters and ODEs are transcribed from the deposited ",
    "Berkeley Madonna model code (Supplementary Material File S1); measured ",
    "inputs are Tables 1, 2 and 4 and Supplementary Tables S1-S2."
  )
  vignette <- "Yang_2026_deha"
  # Intravenous DEHA and MEHA are infused into venous blood; oral DEHA is
  # dosed into the stomach of the GI-transit chain.
  dosing <- c("ven_deha", "ven_meha", "gi_stomach")
  units <- list(
    time = "min",
    dosing = "ng",
    concentration = "ng/mL",
    weight = "kg"
  )

  covariateData <- list()

  # Whole-body PBPK amount states (ng) for DEHA, its gastrointestinal transit
  # chain, and MEHA. Named after the tissue + adipate so the ODE block can be
  # audited line-by-line against the deposited Madonna code. Declared here so
  # checkModelConventions() does not flag them as non-canonical compartments.
  paper_specific_compartments <- c(
    "art_deha",
    "ven_deha",
    "adipose_deha",
    "bone_deha",
    "brain_deha",
    "gut_deha",
    "heart_deha",
    "kidney_deha",
    "liver_deha",
    "lung_deha",
    "muscle_deha",
    "skin_deha",
    "spleen_deha",
    "lymph_deha",
    "gi_stomach",
    "gi_duodenum",
    "gi_jejunum",
    "gi_ileum",
    "gi_colon",
    "gi_feces",
    "art_meha",
    "ven_meha",
    "adipose_meha",
    "bone_meha",
    "brain_meha",
    "gut_meha",
    "heart_meha",
    "kidney_meha",
    "liver_meha",
    "lung_meha",
    "muscle_meha",
    "skin_meha",
    "spleen_meha",
    "lymph_meha"
  )

  compartmentData <- list(
    art_deha = list(analyte = "DEHA", units = "ng", specimen = "whole blood", verified = TRUE),
    ven_deha = list(analyte = "DEHA", units = "ng", specimen = "whole blood", verified = TRUE),
    adipose_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    bone_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    brain_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    gut_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    heart_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    kidney_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    liver_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    lung_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    muscle_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    skin_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    spleen_deha = list(analyte = "DEHA", units = "ng", specimen = "tissue", verified = TRUE),
    lymph_deha = list(analyte = "DEHA", units = "ng", specimen = "lymph", verified = TRUE),
    gi_stomach = list(analyte = "DEHA", units = "ng", specimen = "administration site", verified = TRUE),
    gi_duodenum = list(analyte = "DEHA", units = "ng", specimen = "administration site", verified = TRUE),
    gi_jejunum = list(analyte = "DEHA", units = "ng", specimen = "administration site", verified = TRUE),
    gi_ileum = list(analyte = "DEHA", units = "ng", specimen = "administration site", verified = TRUE),
    gi_colon = list(analyte = "DEHA", units = "ng", specimen = "administration site", verified = TRUE),
    gi_feces = list(analyte = "DEHA", units = "ng", specimen = "faeces", verified = TRUE),
    art_meha = list(analyte = "MEHA", units = "ng", specimen = "whole blood", verified = TRUE),
    ven_meha = list(analyte = "MEHA", units = "ng", specimen = "whole blood", verified = TRUE),
    adipose_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    bone_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    brain_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    gut_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    heart_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    kidney_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    liver_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    lung_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    muscle_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    skin_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    spleen_meha = list(analyte = "MEHA", units = "ng", specimen = "tissue", verified = TRUE),
    lymph_meha = list(analyte = "MEHA", units = "ng", specimen = "lymph", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley, male)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = NULL,
    weight_range = "240-260 g (PK studies); 300-330 g (lymph cannulation); the physiology is that of a 250 g rat",
    sex_female_pct = 0,
    disease_state = "Healthy rats.",
    dose_range = "Intravenous bolus DEHA 5.94 mg/kg and MEHA 5 mg/kg; oral DEHA 99.8 mg/kg (in soybean oil).",
    regions = "Seoul National University, Seoul, Republic of Korea",
    notes = "Per-cohort n = 3 (DEHA IV, MEHA IV, DEHA oral); MEHA after DEHA IV n = 6. DEHA was dosed as a 1:1 mix of unlabelled and deuterated (DEHA-d8) material and reported as DEHA equivalents. Physiological constants are a 250 g rat from the Simcyp V23 rat module and Supplementary Tables S1-S2."
  )

  ini({
    # ----------------------------------------------------------------
    # Fitted parameters (ADAPT 5 maximum likelihood; Table 4, CV% in text).
    # Disposition (PPAMPA,DEHA and Kp,scalar) fitted to IV DEHA plasma DEHA +
    # MEHA; oral absorption (CLLY, Tlag,LY, Km,LY, Fabs) fitted to oral DEHA
    # plasma + mesenteric lymph. Values are the deposited-code point estimates.
    # ----------------------------------------------------------------
    ppampa_deha <- fixed(8.050915946e-07)
    label("DEHA PAMPA permeability (cm/s), fitted") # Table 4 PPAMPA,DEHA = 0.805e-6 (CV 36.1%)
    kp_scaler_deha <- fixed(0.5419000121)
    label("DEHA tissue-partition scalar (unitless), fitted") # Table 4 footnote Kp,scalar = 0.542 (CV 60.0%)
    cl_lymph <- fixed(33.79193675)
    label("Gut-to-lymph uptake clearance CLLY (mL/min), fitted") # Table 4 CLLY = 33.8 (CV 21.4%)
    tlag_lymph <- fixed(66.48311772)
    label("Lymphatic (thoracic-duct) lag time Tlag,LY (min), fitted") # Table 4 Tlag,LY = 66.5 (CV 14.2%)
    km_lymph_deha <- fixed(4643.08215)
    label("DEHA lymphatic-hydrolysis Michaelis constant Km,LY (ng/mL), fitted") # Table 4 Km,LY = 4640 (CV 19.3%)
    fabs_deha <- fixed(0.1391578457)
    label("Absorbed fraction Fabs scaling GI uptake (unitless), fitted") # Table 4 Fabs = 0.139 (CV 10.2%)

    # ----------------------------------------------------------------
    # DEHA physicochemical and binding constants.
    # ----------------------------------------------------------------
    fup_deha <- fixed(0.0208)
    label("DEHA fraction unbound in plasma (unitless)") # Table 4 fup (SPME)
    rbp_deha <- fixed(0.581)
    label("DEHA blood-to-plasma ratio R (unitless)") # Table 4 R; fixed to 1 - hematocrit (theoretical minimum)
    mw_deha <- fixed(370.6)
    label("DEHA molecular weight (g/mol)") # Table 4
    kplasma_deha <- fixed(0.1217)
    label("DEHA plasma disappearance rate constant (1/min)") # Table 1 plasma control k = 0.122
    kmic_deha <- fixed(0.08915)
    label("DEHA liver-microsome disappearance rate constant (1/min)") # Table 1 rat liver microsome control k = 0.0891
    klymph_deha <- fixed(0.0122)
    label("DEHA mesenteric-lymph disappearance rate constant (1/min)") # Table 1 mesenteric lymph k = 0.0122

    # DEHA tissue-homogenate disappearance rate constants (1/min); deposited
    # code values. All match Table 1 control k except adipose (code 0.02436 vs
    # Table 1 0.0312); the code value governs (see Errata). Adipose/brain have
    # no parenchymal metabolism (fBNPP = 0), so their k feeds vascular clearance.
    ka_deha <- fixed(0.02436)
    label("DEHA adipose-homogenate k (1/min)") # deposited code (Table 1 adipose control 0.0312; code value used)
    kb_deha <- fixed(0.07674)
    label("DEHA brain-homogenate k (1/min)") # Table 1 brain control k = 0.0767
    kg_deha <- fixed(0.1392)
    label("DEHA gut-homogenate k (1/min)") # Table 1 gut control k = 0.139
    kh_deha <- fixed(0.03162)
    label("DEHA heart-homogenate k (1/min)") # Table 1 heart control k = 0.0316
    kk_deha <- fixed(0.09506)
    label("DEHA kidney-homogenate k (1/min)") # Table 1 kidney control k = 0.0951
    klg_deha <- fixed(0.03799)
    label("DEHA lung-homogenate k (1/min)") # Table 1 lung control k = 0.0380
    km_deha <- fixed(0.03216)
    label("DEHA muscle-homogenate k (1/min)") # Table 1 muscle control k = 0.0322
    ksp_deha <- fixed(0.06294)
    label("DEHA spleen-homogenate k (1/min)") # Table 1 spleen control k = 0.0629

    # DEHA homogenate unbound fractions fu,T,homogenate (unitless); Table 2.
    futa_deha <- fixed(0.02494)
    label("DEHA adipose fu,homogenate (unitless)") # Table 2 adipose DEHA = 0.0249
    futb_deha <- fixed(0.01704)
    label("DEHA brain fu,homogenate (unitless)") # Table 2 brain DEHA = 0.0170
    futg_deha <- fixed(0.1644)
    label("DEHA gut fu,homogenate (unitless)") # Table 2 gut DEHA = 0.164
    futh_deha <- fixed(0.01668)
    label("DEHA heart fu,homogenate (unitless)") # Table 2 heart DEHA = 0.0167
    futk_deha <- fixed(0.03199)
    label("DEHA kidney fu,homogenate (unitless)") # Table 2 kidneys DEHA = 0.0320
    futli_deha <- fixed(0.02036)
    label("DEHA liver fu,homogenate (unitless)") # Table 2 liver DEHA = 0.0204
    futlg_deha <- fixed(0.01063)
    label("DEHA lung fu,homogenate (unitless)") # Table 2 lungs DEHA = 0.0106
    futm_deha <- fixed(0.01601)
    label("DEHA muscle fu,homogenate (unitless)") # Table 2 muscle DEHA = 0.0160
    futsp_deha <- fixed(0.01172)
    label("DEHA spleen fu,homogenate (unitless)") # Table 2 spleen DEHA = 0.0117

    # DEHA BNPP-inhibited (carboxylesterase) fraction fBNPP = fraction of the
    # tissue disappearance that is parenchymal metabolism; Table 1 % reduction.
    # The complement (1 - fBNPP) is the BNPP-insensitive vascular (lipase) route.
    fbnpp_a_deha <- fixed(0)
    label("DEHA adipose BNPP-inhibited fraction (unitless)") # Table 1 adipose reduction -26.7% -> 0
    fbnpp_b_deha <- fixed(0)
    label("DEHA brain BNPP-inhibited fraction (unitless)") # Table 1 brain reduction -0.995% -> 0
    fbnpp_g_deha <- fixed(0.312)
    label("DEHA gut BNPP-inhibited fraction (unitless)") # Table 1 gut reduction 31.2%
    fbnpp_h_deha <- fixed(0.198)
    label("DEHA heart BNPP-inhibited fraction (unitless)") # Table 1 heart reduction 19.8%
    fbnpp_k_deha <- fixed(0.702)
    label("DEHA kidney BNPP-inhibited fraction (unitless)") # Table 1 kidney reduction 70.2%
    fbnpp_lg_deha <- fixed(0.769)
    label("DEHA lung BNPP-inhibited fraction (unitless)") # Table 1 lung reduction 76.9%
    fbnpp_m_deha <- fixed(0.390)
    label("DEHA muscle BNPP-inhibited fraction (unitless)") # Table 1 muscle reduction 39.0%
    fbnpp_sp_deha <- fixed(0.272)
    label("DEHA spleen BNPP-inhibited fraction (unitless)") # Table 1 spleen reduction 27.2%

    # DEHA in-silico tissue partition coefficients for bone and skin (no
    # homogenate measurement); NOT scaled by kp_scaler (deposited code).
    kp_bone_deha <- fixed(6.508)
    label("DEHA bone-to-plasma partition coefficient (unitless)") # Table 4 bone 6.51 (Berezhkovskiy)
    kp_skin_deha <- fixed(22.79)
    label("DEHA skin-to-plasma partition coefficient (unitless)") # Table 4 skin 22.8 (Schmitt)

    # ----------------------------------------------------------------
    # MEHA physicochemical and binding constants. MEHA distribution is
    # perfusion-limited (measured PAMPA permeability) and its metabolism is
    # essentially fully carboxylesterase-mediated (fBNPP ~ 1), so every tissue
    # disappearance is parenchymal (no MEHA vascular clearance term).
    # ----------------------------------------------------------------
    fup_meha <- fixed(0.0278)
    label("MEHA fraction unbound in plasma (unitless)") # Table 4 fup (RED)
    rbp_meha <- fixed(0.592)
    label("MEHA blood-to-plasma ratio R (unitless)") # Table 4 R
    mw_meha <- fixed(258.35)
    label("MEHA molecular weight (g/mol)") # deposited code (Table 4 reports 258.4)
    ppampa_meha <- fixed(4.88e-5)
    label("MEHA PAMPA permeability (cm/s)") # Table 4 PPAMPA,MEHA = 48.8e-6
    kplasma_meha <- fixed(0.03138)
    label("MEHA plasma disappearance rate constant (1/min)") # Table 1 plasma control k = 0.0316
    kmic_meha <- fixed(0.8113)
    label("MEHA liver-microsome disappearance rate constant (1/min)") # Table 1 rat liver microsome control k = 0.811
    klymph_meha <- fixed(0.146)
    label("MEHA mesenteric-lymph disappearance rate constant (1/min)") # Table 1 mesenteric lymph k = 0.146

    # MEHA tissue-homogenate disappearance rate constants (1/min); Table 1.
    ka_meha <- fixed(2.009)
    label("MEHA adipose-homogenate k (1/min)") # Table 1 adipose control k = 2.01
    kb_meha <- fixed(0.05720)
    label("MEHA brain-homogenate k (1/min)") # Table 1 brain control k = 0.0572
    kg_meha <- fixed(8.066)
    label("MEHA gut-homogenate k (1/min)") # Table 1 gut control k = 8.07
    kh_meha <- fixed(0.5530)
    label("MEHA heart-homogenate k (1/min)") # Table 1 heart control k = 0.553
    kk_meha <- fixed(4.580)
    label("MEHA kidney-homogenate k (1/min)") # Table 1 kidney control k = 4.58
    klg_meha <- fixed(2.011)
    label("MEHA lung-homogenate k (1/min)") # Table 1 lung control k = 2.01
    km_meha <- fixed(0.2287)
    label("MEHA muscle-homogenate k (1/min)") # Table 1 muscle control k = 0.229
    ksp_meha <- fixed(0.2748)
    label("MEHA spleen-homogenate k (1/min)") # Table 1 spleen control k = 0.275

    # MEHA homogenate unbound fractions fu,T,homogenate (unitless). Deposited
    # code values; all match Table 2 except brain (code 0.04432 vs Table 2
    # 0.443). The code value governs (see Errata); brain contributes
    # negligibly to plasma.
    futa_meha <- fixed(0.4306)
    label("MEHA adipose fu,homogenate (unitless)") # Table 2 adipose MEHA = 0.431
    futb_meha <- fixed(0.04432)
    label("MEHA brain fu,homogenate (unitless)") # deposited code (Table 2 brain MEHA 0.443; code value used)
    futg_meha <- fixed(0.1348)
    label("MEHA gut fu,homogenate (unitless)") # Table 2 gut MEHA = 0.135
    futh_meha <- fixed(0.3546)
    label("MEHA heart fu,homogenate (unitless)") # Table 2 heart MEHA = 0.355
    futk_meha <- fixed(0.6838)
    label("MEHA kidney fu,homogenate (unitless)") # Table 2 kidneys MEHA = 0.684
    futli_meha <- fixed(0.5212)
    label("MEHA liver fu,homogenate (unitless)") # Table 2 liver MEHA = 0.521
    futlg_meha <- fixed(0.1305)
    label("MEHA lung fu,homogenate (unitless)") # Table 2 lungs MEHA = 0.131
    futm_meha <- fixed(0.4458)
    label("MEHA muscle fu,homogenate (unitless)") # Table 2 muscle MEHA = 0.446
    futsp_meha <- fixed(0.2684)
    label("MEHA spleen fu,homogenate (unitless)") # Table 2 spleen MEHA = 0.268

    # MEHA in-silico partition coefficients for bone and skin (Kp_scaler = 1).
    kp_bone_meha <- fixed(0.6388)
    label("MEHA bone-to-plasma partition coefficient (unitless)") # Table 4 bone 0.639 (Schmitt)
    kp_skin_meha <- fixed(3.276)
    label("MEHA skin-to-plasma partition coefficient (unitless)") # Table 4 skin 3.28 (Schmitt)

    # ----------------------------------------------------------------
    # System physiology of a 250 g rat (Simcyp V23 rat module; Supplementary
    # Tables S1-S2 and deposited code). Blood flows (mL/min).
    # ----------------------------------------------------------------
    qc <- fixed(80)
    label("Cardiac output / lung blood flow (mL/min)") # code Q = 80; Table S2 lung 80
    q_adipose <- fixed(4.72)
    label("Adipose blood flow (mL/min)") # Table S2 adipose
    q_bone <- fixed(8.08)
    label("Bone blood flow (mL/min)") # Table S2 bone
    q_brain <- fixed(1.12)
    label("Brain blood flow (mL/min)") # Table S2 brain
    q_gut <- fixed(8.08)
    label("Gut blood flow (mL/min)") # Table S2 gut
    q_heart <- fixed(3.2)
    label("Heart blood flow (mL/min)") # Table S2 heart
    q_kidney <- fixed(11.6)
    label("Kidney blood flow (mL/min)") # Table S2 kidneys
    q_liver <- fixed(19.36)
    label("Total hepatic blood flow (mL/min)") # Table S2 liver (hepatic artery + gut + spleen)
    q_muscle <- fixed(18.96)
    label("Muscle blood flow (mL/min)") # Table S2 muscle
    q_skin <- fixed(4.08)
    label("Skin blood flow (mL/min)") # Table S2 skin
    q_spleen <- fixed(0.88)
    label("Spleen blood flow (mL/min)") # Table S2 spleen
    q_lymph <- fixed(0.0416667)
    label("Lymph flow (mL/min)") # Table S1 lymph flow 2.5 mL/h

    # Tissue volumes (mL).
    v_venous <- fixed(11.49333)
    label("Venous blood volume (mL)") # code (10.04+7.2)*2/3
    v_arterial <- fixed(5.74667)
    label("Arterial blood volume (mL)") # code (10.04+7.2)*1/3
    v_adipose <- fixed(16.66)
    label("Adipose volume (mL)") # Table S2 adipose
    v_bone <- fixed(15.68)
    label("Bone volume (mL)") # Table S2 bone
    v_brain <- fixed(1.24)
    label("Brain volume (mL)") # Table S2 brain
    v_gut <- fixed(6.19)
    label("Gut volume (mL)") # Table S2 gut
    v_heart <- fixed(1.05)
    label("Heart volume (mL)") # Table S2 heart
    v_kidney <- fixed(2.19)
    label("Kidney volume (mL)") # Table S2 kidneys
    v_liver <- fixed(8.57)
    label("Liver volume (mL)") # Table S2 liver
    v_lung <- fixed(1.24)
    label("Lung volume (mL)") # Table S2 lungs
    v_muscle <- fixed(116.13)
    label("Muscle volume (mL)") # Table S2 muscle
    v_skin <- fixed(39.41)
    label("Skin volume (mL)") # Table S2 skin
    v_spleen <- fixed(0.57)
    label("Spleen volume (mL)") # Table S2 spleen
    v_lymph <- fixed(17.24)
    label("Lymph volume (mL)") # code Vlymph = blood volume

    hct <- fixed(0.417633)
    label("Hematocrit (unitless)") # code 7.2/(10.04+7.2)
    mppgl <- fixed(44.5)
    label("Microsomal protein per gram liver (mg/g)") # code / doi:10.1002/bdd.1959
    vmic <- fixed(2)
    label("Microsomal incubation volume per protein (mL/mg)") # code 1/0.5

    # Tissue effective surface-area coefficients (m^2/g tissue); Table S2.
    # Seff = coefficient * tissue weight; used with PAMPA permeability for fd.
    seffc_adipose <- fixed(3.02)
    label("Adipose effective-surface coefficient (m2/g)") # Table S2 adipose
    seffc_bone <- fixed(2.34)
    label("Bone effective-surface coefficient (m2/g)") # Table S2 bone
    seffc_brain <- fixed(1.67)
    label("Brain effective-surface coefficient (m2/g)") # Table S2 brain
    seffc_gut <- fixed(16.1)
    label("Gut effective-surface coefficient (m2/g)") # Table S2 gut
    seffc_heart <- fixed(18.1)
    label("Heart effective-surface coefficient (m2/g)") # Table S2 heart
    seffc_kidney <- fixed(37.4)
    label("Kidney effective-surface coefficient (m2/g)") # Table S2 kidneys
    seffc_liver <- fixed(14.8)
    label("Liver effective-surface coefficient (m2/g)") # Table S2 liver
    seffc_lung <- fixed(14.4)
    label("Lung effective-surface coefficient (m2/g)") # Table S2 lungs
    seffc_muscle <- fixed(5.71)
    label("Muscle effective-surface coefficient (m2/g)") # Table S2 muscle
    seffc_skin <- fixed(0.843)
    label("Skin effective-surface coefficient (m2/g)") # Table S2 skin
    seffc_spleen <- fixed(4.56)
    label("Spleen effective-surface coefficient (m2/g)") # Table S2 spleen

    # Oral DEHA absorption (deposited code; Supplementary Eqs S15-S27).
    peff_human_deha <- fixed(4.885)
    label("DEHA human effective permeability (1e-4 cm/s)") # code Peff_human_DEHA
    radius_intestine <- fixed(0.18)
    label("Intestinal radius (cm)") # code / Issa et al. 2003
    k_st <- fixed(0.024)
    label("Stomach transit rate constant (1/min)") # Table S1 / code 1.44/60
    k_du <- fixed(0.041)
    label("Duodenum transit rate constant (1/min)") # Table S1 / code 2.46/60
    k_je <- fixed(0.0136667)
    label("Jejunum transit rate constant (1/min)") # Table S1 / code 0.82/60
    k_il <- fixed(0.0136667)
    label("Ileum transit rate constant (1/min)") # Table S1 / code 0.82/60
    k_co <- fixed(0.0026667)
    label("Colon transit rate constant (1/min)") # Table S1 / code 0.16/60

    # ----------------------------------------------------------------
    # Residual error. The model is a deterministic PBPK fit (ADAPT 5 maximum
    # likelihood on pooled data); no between-subject variance or residual-error
    # magnitude is reported, so propSd is held at zero rather than invented.
    # ----------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # ================================================================
    # Derived physiology (deposited code).
    # ================================================================
    # Hepatic-artery flow = total hepatic flow minus splanchnic inflows; the
    # rest-of-body flow carries the residual cardiac output with no tissue.
    q_hepatic_artery <- q_liver - q_gut - q_spleen
    q_rest <- qc - (q_adipose + q_bone + q_brain + q_heart + q_kidney + q_liver + q_muscle + q_skin)

    # Tissue weights (g): adipose 0.9 g/mL, bone 1.85 g/mL, else 1 g/mL.
    w_adipose <- v_adipose * 0.9
    w_bone <- v_bone * 1.85
    w_brain <- v_brain
    w_gut <- v_gut
    w_heart <- v_heart
    w_kidney <- v_kidney
    w_liver <- v_liver
    w_lung <- v_lung
    w_muscle <- v_muscle
    w_skin <- v_skin
    w_spleen <- v_spleen

    # Permeability-surface products (mL/min) = PAMPA perm * effective area, unit
    # bridge cm/s * m^2 -> mL/min via * 60 * 10000 (deposited code Eq S1).
    ps_adipose_de <- ppampa_deha * seffc_adipose * w_adipose * 60 * 10000
    ps_bone_de <- ppampa_deha * seffc_bone * w_bone * 60 * 10000
    ps_brain_de <- ppampa_deha * seffc_brain * w_brain * 60 * 10000
    ps_gut_de <- ppampa_deha * seffc_gut * w_gut * 60 * 10000
    ps_heart_de <- ppampa_deha * seffc_heart * w_heart * 60 * 10000
    ps_kidney_de <- ppampa_deha * seffc_kidney * w_kidney * 60 * 10000
    ps_liver_de <- ppampa_deha * seffc_liver * w_liver * 60 * 10000
    ps_lung_de <- ppampa_deha * seffc_lung * w_lung * 60 * 10000
    ps_muscle_de <- ppampa_deha * seffc_muscle * w_muscle * 60 * 10000
    ps_skin_de <- ppampa_deha * seffc_skin * w_skin * 60 * 10000
    ps_spleen_de <- ppampa_deha * seffc_spleen * w_spleen * 60 * 10000

    # Fractional distribution fd (Eq S2): fraction permeating per capillary pass.
    fd_adipose_de <- 1 - exp(-fup_deha * ps_adipose_de / rbp_deha / q_adipose)
    fd_bone_de <- 1 - exp(-fup_deha * ps_bone_de / rbp_deha / q_bone)
    fd_brain_de <- 1 - exp(-fup_deha * ps_brain_de / rbp_deha / q_brain)
    fd_gut_de <- 1 - exp(-fup_deha * ps_gut_de / rbp_deha / q_gut)
    fd_heart_de <- 1 - exp(-fup_deha * ps_heart_de / rbp_deha / q_heart)
    fd_kidney_de <- 1 - exp(-fup_deha * ps_kidney_de / rbp_deha / q_kidney)
    fd_liver_de <- 1 - exp(-fup_deha * ps_liver_de / rbp_deha / q_liver)
    fd_lung_de <- 1 - exp(-fup_deha * ps_lung_de / rbp_deha / qc)
    fd_muscle_de <- 1 - exp(-fup_deha * ps_muscle_de / rbp_deha / q_muscle)
    fd_skin_de <- 1 - exp(-fup_deha * ps_skin_de / rbp_deha / q_skin)
    fd_spleen_de <- 1 - exp(-fup_deha * ps_spleen_de / rbp_deha / q_spleen)

    ps_adipose_me <- ppampa_meha * seffc_adipose * w_adipose * 60 * 10000
    ps_bone_me <- ppampa_meha * seffc_bone * w_bone * 60 * 10000
    ps_brain_me <- ppampa_meha * seffc_brain * w_brain * 60 * 10000
    ps_gut_me <- ppampa_meha * seffc_gut * w_gut * 60 * 10000
    ps_heart_me <- ppampa_meha * seffc_heart * w_heart * 60 * 10000
    ps_kidney_me <- ppampa_meha * seffc_kidney * w_kidney * 60 * 10000
    ps_liver_me <- ppampa_meha * seffc_liver * w_liver * 60 * 10000
    ps_lung_me <- ppampa_meha * seffc_lung * w_lung * 60 * 10000
    ps_muscle_me <- ppampa_meha * seffc_muscle * w_muscle * 60 * 10000
    ps_skin_me <- ppampa_meha * seffc_skin * w_skin * 60 * 10000
    ps_spleen_me <- ppampa_meha * seffc_spleen * w_spleen * 60 * 10000

    fd_adipose_me <- 1 - exp(-fup_meha * ps_adipose_me / rbp_meha / q_adipose)
    fd_bone_me <- 1 - exp(-fup_meha * ps_bone_me / rbp_meha / q_bone)
    fd_brain_me <- 1 - exp(-fup_meha * ps_brain_me / rbp_meha / q_brain)
    fd_gut_me <- 1 - exp(-fup_meha * ps_gut_me / rbp_meha / q_gut)
    fd_heart_me <- 1 - exp(-fup_meha * ps_heart_me / rbp_meha / q_heart)
    fd_kidney_me <- 1 - exp(-fup_meha * ps_kidney_me / rbp_meha / q_kidney)
    fd_liver_me <- 1 - exp(-fup_meha * ps_liver_me / rbp_meha / q_liver)
    fd_lung_me <- 1 - exp(-fup_meha * ps_lung_me / rbp_meha / qc)
    fd_muscle_me <- 1 - exp(-fup_meha * ps_muscle_me / rbp_meha / q_muscle)
    fd_skin_me <- 1 - exp(-fup_meha * ps_skin_me / rbp_meha / q_skin)
    fd_spleen_me <- 1 - exp(-fup_meha * ps_spleen_me / rbp_meha / q_spleen)

    # Dilution-corrected tissue unbound fractions (Eq 1, df = 4) and
    # tissue-to-plasma partition coefficients Kp = fup / fu,T * scalar (code).
    futa_de_c <- 0.25 / (1 / futa_deha - 1 + 0.25)
    futb_de_c <- 0.25 / (1 / futb_deha - 1 + 0.25)
    futg_de_c <- 0.25 / (1 / futg_deha - 1 + 0.25)
    futh_de_c <- 0.25 / (1 / futh_deha - 1 + 0.25)
    futk_de_c <- 0.25 / (1 / futk_deha - 1 + 0.25)
    futli_de_c <- 0.25 / (1 / futli_deha - 1 + 0.25)
    futlg_de_c <- 0.25 / (1 / futlg_deha - 1 + 0.25)
    futm_de_c <- 0.25 / (1 / futm_deha - 1 + 0.25)
    futsp_de_c <- 0.25 / (1 / futsp_deha - 1 + 0.25)

    kp_adipose_de <- fup_deha / futa_de_c * kp_scaler_deha
    kp_brain_de <- fup_deha / futb_de_c * kp_scaler_deha
    kp_gut_de <- fup_deha / futg_de_c * kp_scaler_deha
    kp_heart_de <- fup_deha / futh_de_c * kp_scaler_deha
    kp_kidney_de <- fup_deha / futk_de_c * kp_scaler_deha
    kp_liver_de <- fup_deha / futli_de_c * kp_scaler_deha
    kp_lung_de <- fup_deha / futlg_de_c * kp_scaler_deha
    kp_muscle_de <- fup_deha / futm_de_c * kp_scaler_deha
    kp_spleen_de <- fup_deha / futsp_de_c * kp_scaler_deha

    futa_me_c <- 0.25 / (1 / futa_meha - 1 + 0.25)
    futb_me_c <- 0.25 / (1 / futb_meha - 1 + 0.25)
    futg_me_c <- 0.25 / (1 / futg_meha - 1 + 0.25)
    futh_me_c <- 0.25 / (1 / futh_meha - 1 + 0.25)
    futk_me_c <- 0.25 / (1 / futk_meha - 1 + 0.25)
    futli_me_c <- 0.25 / (1 / futli_meha - 1 + 0.25)
    futlg_me_c <- 0.25 / (1 / futlg_meha - 1 + 0.25)
    futm_me_c <- 0.25 / (1 / futm_meha - 1 + 0.25)
    futsp_me_c <- 0.25 / (1 / futsp_meha - 1 + 0.25)

    kp_adipose_me <- fup_meha / futa_me_c
    kp_brain_me <- fup_meha / futb_me_c
    kp_gut_me <- fup_meha / futg_me_c
    kp_heart_me <- fup_meha / futh_me_c
    kp_kidney_me <- fup_meha / futk_me_c
    kp_liver_me <- fup_meha / futli_me_c
    kp_lung_me <- fup_meha / futlg_me_c
    kp_muscle_me <- fup_meha / futm_me_c
    kp_spleen_me <- fup_meha / futsp_me_c

    # Unbound tissue disappearance rate constants ku,T = k,homogenate/fu * df.
    kua_de <- ka_deha / futa_deha * 4
    kub_de <- kb_deha / futb_deha * 4
    kug_de <- kg_deha / futg_deha * 4
    kuh_de <- kh_deha / futh_deha * 4
    kuk_de <- kk_deha / futk_deha * 4
    kulg_de <- klg_deha / futlg_deha * 4
    kum_de <- km_deha / futm_deha * 4
    kusp_de <- ksp_deha / futsp_deha * 4

    kua_me <- ka_meha / futa_meha * 4
    kub_me <- kb_meha / futb_meha * 4
    kug_me <- kg_meha / futg_meha * 4
    kuh_me <- kh_meha / futh_meha * 4
    kuk_me <- kk_meha / futk_meha * 4
    kulg_me <- klg_meha / futlg_meha * 4
    kum_me <- km_meha / futm_meha * 4
    kusp_me <- ksp_meha / futsp_meha * 4

    # MEHA parenchymal rate constants acting on tissue concentration
    # (k,MEHA = ku * fu,T,corrected); full carboxylesterase metabolism.
    ka_me_t <- kua_me * futa_me_c
    kb_me_t <- kub_me * futb_me_c
    kg_me_t <- kug_me * futg_me_c
    kh_me_t <- kuh_me * futh_me_c
    kk_me_t <- kuk_me * futk_me_c
    klg_me_t <- kulg_me * futlg_me_c
    km_me_t <- kum_me * futm_me_c
    ksp_me_t <- kusp_me * futsp_me_c

    # Hepatic intrinsic clearance (mL/min) scaled from microsomal stability.
    clint_mic_de <- kmic_deha * vmic / fup_deha * mppgl * w_liver
    clint_mic_me <- kmic_meha * vmic / fup_meha * mppgl * w_liver

    # Arterial (sum over arterial-side tissues incl. gut/spleen) and venous
    # (lung) vascular lipase clearances for DEHA (Eqs 6, 8).
    cl_add_ar_de <- (kua_de * (1 - fbnpp_a_deha) * v_adipose +
      kub_de * (1 - fbnpp_b_deha) * v_brain +
      kug_de * (1 - fbnpp_g_deha) * v_gut +
      kuh_de * (1 - fbnpp_h_deha) * v_heart +
      kuk_de * (1 - fbnpp_k_deha) * v_kidney +
      kum_de * (1 - fbnpp_m_deha) * v_muscle +
      kusp_de * (1 - fbnpp_sp_deha) * v_spleen) * fup_deha
    cl_add_ve_de <- kulg_de * (1 - fbnpp_lg_deha) * v_lung * fup_deha

    # First-order absorption rate constants per intestinal segment (Eqs S17-S21).
    peff_rat_de <- (peff_human_deha - 0.03) / 3.6 * 60 / 10000
    ka_du <- 2 * 1.7 * peff_rat_de * 0.49 / radius_intestine
    ka_je <- 2 * 1.7 * peff_rat_de * 1 / radius_intestine
    ka_il <- 2 * 1.7 * peff_rat_de * 0.58 / radius_intestine
    ka_co <- 2 * 1.7 * peff_rat_de * 0.03 / radius_intestine
    # Weibull dissolution (alpha = beta = 1 min, lag 0.1 min): immediate release.
    frac_dis <- (time > 0.1) * (1 - exp(-(time - 0.1)))

    mwr <- mw_meha / mw_deha

    # ================================================================
    # Concentrations (ng/mL). Blood pools carry drug as blood concentration
    # (amount / volume); plasma concentration = that / R.
    # ================================================================
    cb_art_de <- art_deha / v_arterial
    cb_ven_de <- ven_deha / v_venous
    cp_art_de <- art_deha / (rbp_deha * v_arterial)
    cp_ven_de <- ven_deha / (rbp_deha * v_venous)
    c_adipose_de <- adipose_deha / v_adipose
    c_bone_de <- bone_deha / v_bone
    c_brain_de <- brain_deha / v_brain
    c_gut_de <- gut_deha / v_gut
    c_heart_de <- heart_deha / v_heart
    c_kidney_de <- kidney_deha / v_kidney
    c_liver_de <- liver_deha / v_liver
    c_lung_de <- lung_deha / v_lung
    c_muscle_de <- muscle_deha / v_muscle
    c_skin_de <- skin_deha / v_skin
    c_spleen_de <- spleen_deha / v_spleen
    c_lymph_de <- lymph_deha / v_lymph

    cb_art_me <- art_meha / v_arterial
    cb_ven_me <- ven_meha / v_venous
    cp_art_me <- art_meha / (rbp_meha * v_arterial)
    cp_ven_me <- ven_meha / (rbp_meha * v_venous)
    c_adipose_me <- adipose_meha / v_adipose
    c_bone_me <- bone_meha / v_bone
    c_brain_me <- brain_meha / v_brain
    c_gut_me <- gut_meha / v_gut
    c_heart_me <- heart_meha / v_heart
    c_kidney_me <- kidney_meha / v_kidney
    c_liver_me <- liver_meha / v_liver
    c_lung_me <- lung_meha / v_lung
    c_muscle_me <- muscle_meha / v_muscle
    c_skin_me <- skin_meha / v_skin
    c_spleen_me <- spleen_meha / v_spleen
    c_lymph_me <- lymph_meha / v_lymph

    # ================================================================
    # DEHA elimination rates per compartment (ng/min). Parenchymal metabolism
    # acts on unbound tissue concentration (ku * fBNPP * fu,T * C * V); the
    # complement forms the vascular route accounted for above.
    # ================================================================
    xmet_a_de <- kua_de * fbnpp_a_deha * c_adipose_de * futa_de_c * v_adipose
    xmet_b_de <- kub_de * fbnpp_b_deha * c_brain_de * futb_de_c * v_brain
    xmet_g_de <- kug_de * fbnpp_g_deha * c_gut_de * futg_de_c * v_gut
    xmet_h_de <- kuh_de * fbnpp_h_deha * c_heart_de * futh_de_c * v_heart
    xmet_k_de <- kuk_de * fbnpp_k_deha * c_kidney_de * futk_de_c * v_kidney
    xmet_lg_de <- kulg_de * fbnpp_lg_deha * c_lung_de * futlg_de_c * v_lung
    xmet_m_de <- kum_de * fbnpp_m_deha * c_muscle_de * futm_de_c * v_muscle
    xmet_sp_de <- kusp_de * fbnpp_sp_deha * c_spleen_de * futsp_de_c * v_spleen
    xmet_li_de <- clint_mic_de * c_liver_de * fup_deha / kp_liver_de
    xmet_ar_de <- kplasma_deha * cp_art_de * v_arterial * (1 - hct) + cl_add_ar_de * cp_art_de
    xmet_ve_de <- kplasma_deha * cp_ven_de * v_venous * (1 - hct) + cl_add_ve_de * cp_ven_de
    xmet_lymph_de <- km_lymph_deha * klymph_deha / (km_lymph_deha + c_lymph_de) * c_lymph_de * v_lymph

    # Tissue in/out mass-transfer rates (ng/min). Inflow = fd * Q * C_in,blood;
    # the (1 - fd) bypass fraction returns unchanged to the downstream pool.
    xin_a_de <- fd_adipose_de * q_adipose * cb_art_de
    xin_bo_de <- fd_bone_de * q_bone * cb_art_de
    xin_b_de <- fd_brain_de * q_brain * cb_art_de
    xin_g_de <- fd_gut_de * q_gut * cb_art_de
    xin_h_de <- fd_heart_de * q_heart * cb_art_de
    xin_k_de <- fd_kidney_de * q_kidney * cb_art_de
    xin_m_de <- fd_muscle_de * q_muscle * cb_art_de
    xin_sk_de <- fd_skin_de * q_skin * cb_art_de
    xin_sp_de <- fd_spleen_de * q_spleen * cb_art_de
    xin_lg_de <- fd_lung_de * qc * cb_ven_de

    xout_a_de <- fd_adipose_de * q_adipose * c_adipose_de / kp_adipose_de * rbp_deha
    xout_bo_de <- fd_bone_de * q_bone * c_bone_de / kp_bone_deha * rbp_deha
    xout_b_de <- fd_brain_de * q_brain * c_brain_de / kp_brain_de * rbp_deha
    xout_g_de <- fd_gut_de * q_gut * c_gut_de / kp_gut_de * rbp_deha + cl_lymph * c_gut_de
    xout_h_de <- fd_heart_de * q_heart * c_heart_de / kp_heart_de * rbp_deha
    xout_k_de <- fd_kidney_de * q_kidney * c_kidney_de / kp_kidney_de * rbp_deha
    xout_li_de <- fd_liver_de * q_liver * c_liver_de / kp_liver_de * rbp_deha
    xout_m_de <- fd_muscle_de * q_muscle * c_muscle_de / kp_muscle_de * rbp_deha
    xout_sk_de <- fd_skin_de * q_skin * c_skin_de / kp_skin_deha * rbp_deha
    xout_sp_de <- fd_spleen_de * q_spleen * c_spleen_de / kp_spleen_de * rbp_deha
    xout_lg_de <- fd_lung_de * qc * c_lung_de / kp_lung_de * rbp_deha

    # Liver inflow: hepatic artery + gut outflow (equilibrated + bypass) +
    # spleen outflow (equilibrated + bypass), all times fd_liver (Eq S7-S8).
    xin_li_de <- fd_liver_de * (q_hepatic_artery * cb_art_de +
      fd_gut_de * q_gut * c_gut_de / kp_gut_de * rbp_deha +
      xin_g_de * (1 - fd_gut_de) / fd_gut_de +
      xout_sp_de + xin_sp_de * (1 - fd_spleen_de) / fd_spleen_de)

    # Oral absorption input to the gut tissue (ng/min).
    xin_oral_g_de <- fabs_deha * frac_dis *
      (ka_du * gi_duodenum + ka_je * gi_jejunum + ka_il * gi_ileum + ka_co * gi_colon)

    # Lymph: gut uptake in, flow + saturable hydrolysis out; return to venous
    # blood after the thoracic-duct lag via delay() of the lymph amount.
    xin_lymph_de <- cl_lymph * c_gut_de
    xout_lymph_de <- q_lymph * c_lymph_de
    xreturn_lymph_de <- q_lymph * delay(lymph_deha, tlag_lymph) / v_lymph

    # Arterial / venous blood inflow (Eqs S10-S11). Venous collects the
    # non-splanchnic tissues and the liver (each equilibrated + bypass), the
    # residual cardiac output, and the delayed lymphatic return; arterial
    # collects the lung outflow.
    xin_ar_de <- xout_lg_de + xin_lg_de * (1 - fd_lung_de) / fd_lung_de
    xin_ve_de <- xout_a_de + xin_a_de * (1 - fd_adipose_de) / fd_adipose_de +
      xout_bo_de + xin_bo_de * (1 - fd_bone_de) / fd_bone_de +
      xout_b_de + xin_b_de * (1 - fd_brain_de) / fd_brain_de +
      xout_h_de + xin_h_de * (1 - fd_heart_de) / fd_heart_de +
      xout_k_de + xin_k_de * (1 - fd_kidney_de) / fd_kidney_de +
      xout_li_de + xin_li_de * (1 - fd_liver_de) / fd_liver_de +
      xout_m_de + xin_m_de * (1 - fd_muscle_de) / fd_muscle_de +
      xout_sk_de + xin_sk_de * (1 - fd_skin_de) / fd_skin_de +
      q_rest * cb_art_de + xreturn_lymph_de
    xout_ar_de <- qc * cb_art_de
    xout_ve_de <- qc * cb_ven_de

    # ================================================================
    # MEHA source terms = DEHA elimination in each compartment * MW ratio
    # (complete conversion, fm = 1; Eqs S12-S14).
    # ================================================================
    rmeha_a <- xmet_a_de * mwr
    rmeha_b <- xmet_b_de * mwr
    rmeha_g <- xmet_g_de * mwr
    rmeha_h <- xmet_h_de * mwr
    rmeha_k <- xmet_k_de * mwr
    rmeha_lg <- xmet_lg_de * mwr
    rmeha_m <- xmet_m_de * mwr
    rmeha_sp <- xmet_sp_de * mwr
    rmeha_li <- xmet_li_de * mwr
    rmeha_ar <- xmet_ar_de * mwr
    rmeha_ve <- xmet_ve_de * mwr
    rmeha_lymph <- xmet_lymph_de * mwr

    # MEHA elimination (parenchymal + plasma only; no vascular term).
    xmet_a_me <- ka_me_t * c_adipose_me * v_adipose
    xmet_b_me <- kb_me_t * c_brain_me * v_brain
    xmet_g_me <- kg_me_t * c_gut_me * v_gut
    xmet_h_me <- kh_me_t * c_heart_me * v_heart
    xmet_k_me <- kk_me_t * c_kidney_me * v_kidney
    xmet_lg_me <- klg_me_t * c_lung_me * v_lung
    xmet_m_me <- km_me_t * c_muscle_me * v_muscle
    xmet_sp_me <- ksp_me_t * c_spleen_me * v_spleen
    xmet_li_me <- clint_mic_me * c_liver_me * fup_meha / kp_liver_me
    xmet_ar_me <- kplasma_meha * cp_art_me * v_arterial * (1 - hct)
    xmet_ve_me <- kplasma_meha * cp_ven_me * v_venous * (1 - hct)
    xmet_lymph_me <- klymph_meha * c_lymph_me * v_lymph

    xin_a_me <- fd_adipose_me * q_adipose * cb_art_me + rmeha_a
    xin_bo_me <- fd_bone_me * q_bone * cb_art_me
    xin_b_me <- fd_brain_me * q_brain * cb_art_me + rmeha_b
    xin_g_me <- fd_gut_me * q_gut * cb_art_me + rmeha_g
    xin_h_me <- fd_heart_me * q_heart * cb_art_me + rmeha_h
    xin_k_me <- fd_kidney_me * q_kidney * cb_art_me + rmeha_k
    xin_m_me <- fd_muscle_me * q_muscle * cb_art_me + rmeha_m
    xin_sk_me <- fd_skin_me * q_skin * cb_art_me
    xin_sp_me <- fd_spleen_me * q_spleen * cb_art_me + rmeha_sp
    xin_lg_me <- fd_lung_me * qc * cb_ven_me + rmeha_lg

    xout_a_me <- fd_adipose_me * q_adipose * c_adipose_me / kp_adipose_me * rbp_meha
    xout_bo_me <- fd_bone_me * q_bone * c_bone_me / kp_bone_meha * rbp_meha
    xout_b_me <- fd_brain_me * q_brain * c_brain_me / kp_brain_me * rbp_meha
    xout_g_me <- fd_gut_me * q_gut * c_gut_me / kp_gut_me * rbp_meha
    xout_h_me <- fd_heart_me * q_heart * c_heart_me / kp_heart_me * rbp_meha
    xout_k_me <- fd_kidney_me * q_kidney * c_kidney_me / kp_kidney_me * rbp_meha
    xout_li_me <- fd_liver_me * q_liver * c_liver_me / kp_liver_me * rbp_meha
    xout_m_me <- fd_muscle_me * q_muscle * c_muscle_me / kp_muscle_me * rbp_meha
    xout_sk_me <- fd_skin_me * q_skin * c_skin_me / kp_skin_meha * rbp_meha
    xout_sp_me <- fd_spleen_me * q_spleen * c_spleen_me / kp_spleen_me * rbp_meha
    xout_lg_me <- fd_lung_me * qc * c_lung_me / kp_lung_me * rbp_meha

    xin_li_me <- fd_liver_me * (q_hepatic_artery * cb_art_me +
      fd_gut_me * q_gut * c_gut_me / kp_gut_me * rbp_meha +
      xin_g_me * (1 - fd_gut_me) / fd_gut_me +
      xout_sp_me + xin_sp_me * (1 - fd_spleen_me) / fd_spleen_me) + rmeha_li

    xout_lymph_me <- q_lymph * c_lymph_me
    xin_ar_me <- xout_lg_me + xin_lg_me * (1 - fd_lung_me) / fd_lung_me + rmeha_ar
    xin_ve_me <- xout_a_me + xin_a_me * (1 - fd_adipose_me) / fd_adipose_me +
      xout_bo_me + xin_bo_me * (1 - fd_bone_me) / fd_bone_me +
      xout_b_me + xin_b_me * (1 - fd_brain_me) / fd_brain_me +
      xout_h_me + xin_h_me * (1 - fd_heart_me) / fd_heart_me +
      xout_k_me + xin_k_me * (1 - fd_kidney_me) / fd_kidney_me +
      xout_li_me + xin_li_me * (1 - fd_liver_me) / fd_liver_me +
      xout_m_me + xin_m_me * (1 - fd_muscle_me) / fd_muscle_me +
      xout_sk_me + xin_sk_me * (1 - fd_skin_me) / fd_skin_me +
      q_rest * cb_art_me + rmeha_ve + xout_lymph_me
    xout_ar_me <- qc * cb_art_me
    xout_ve_me <- qc * cb_ven_me

    # ================================================================
    # Differential equations (ng/min). Blood pools dose by compartment:
    # IV DEHA -> ven_deha, IV MEHA -> ven_meha, oral DEHA -> gi_stomach.
    # ================================================================
    d/dt(art_deha) <- xin_ar_de - xout_ar_de - xmet_ar_de
    d/dt(ven_deha) <- xin_ve_de - xout_ve_de - xmet_ve_de
    d/dt(adipose_deha) <- xin_a_de - xout_a_de - xmet_a_de
    d/dt(bone_deha) <- xin_bo_de - xout_bo_de
    d/dt(brain_deha) <- xin_b_de - xout_b_de - xmet_b_de
    d/dt(gut_deha) <- xin_g_de + xin_oral_g_de - xout_g_de - xmet_g_de
    d/dt(heart_deha) <- xin_h_de - xout_h_de - xmet_h_de
    d/dt(kidney_deha) <- xin_k_de - xout_k_de - xmet_k_de
    d/dt(liver_deha) <- xin_li_de - xout_li_de - xmet_li_de
    d/dt(lung_deha) <- xin_lg_de - xout_lg_de - xmet_lg_de
    d/dt(muscle_deha) <- xin_m_de - xout_m_de - xmet_m_de
    d/dt(skin_deha) <- xin_sk_de - xout_sk_de
    d/dt(spleen_deha) <- xin_sp_de - xout_sp_de - xmet_sp_de
    d/dt(lymph_deha) <- xin_lymph_de - xout_lymph_de - xmet_lymph_de

    d/dt(gi_stomach) <- -k_st * gi_stomach
    d/dt(gi_duodenum) <- k_st * gi_stomach - frac_dis * ka_du * gi_duodenum - k_du * gi_duodenum
    d/dt(gi_jejunum) <- k_du * gi_duodenum - frac_dis * ka_je * gi_jejunum - k_je * gi_jejunum
    d/dt(gi_ileum) <- k_je * gi_jejunum - frac_dis * ka_il * gi_ileum - k_il * gi_ileum
    d/dt(gi_colon) <- k_il * gi_ileum - frac_dis * ka_co * gi_colon - k_co * gi_colon
    d/dt(gi_feces) <- k_co * gi_colon

    d/dt(art_meha) <- xin_ar_me - xout_ar_me - xmet_ar_me
    d/dt(ven_meha) <- xin_ve_me - xout_ve_me - xmet_ve_me
    d/dt(adipose_meha) <- xin_a_me - xout_a_me - xmet_a_me
    d/dt(bone_meha) <- xin_bo_me - xout_bo_me
    d/dt(brain_meha) <- xin_b_me - xout_b_me - xmet_b_me
    d/dt(gut_meha) <- xin_g_me - xout_g_me - xmet_g_me
    d/dt(heart_meha) <- xin_h_me - xout_h_me - xmet_h_me
    d/dt(kidney_meha) <- xin_k_me - xout_k_me - xmet_k_me
    d/dt(liver_meha) <- xin_li_me - xout_li_me - xmet_li_me
    d/dt(lung_meha) <- xin_lg_me - xout_lg_me - xmet_lg_me
    d/dt(muscle_meha) <- xin_m_me - xout_m_me - xmet_m_me
    d/dt(skin_meha) <- xin_sk_me - xout_sk_me
    d/dt(spleen_meha) <- xin_sp_me - xout_sp_me - xmet_sp_me
    d/dt(lymph_meha) <- rmeha_lymph - xout_lymph_me - xmet_lymph_me

    # ================================================================
    # Observations. The paper reports arterial plasma concentrations (IV NCA)
    # and, for oral dosing, plasma and mesenteric-lymph concentrations. The
    # metabolite and lymph concentrations are returned as output columns; the
    # primary endpoint is arterial plasma DEHA.
    # ================================================================
    Cc <- cp_art_de
    Cc ~ prop(propSd)
  })
}
