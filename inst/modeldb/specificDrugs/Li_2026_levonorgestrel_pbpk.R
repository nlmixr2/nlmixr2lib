Li_2026_levonorgestrel_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body with a female-reproductive-tract sub-model; original",
    "implementation in mrgsolve). Levonorgestrel (LNG) given orally, vaginally",
    "or by intrauterine system (IUS), from Li et al. 2026. Twenty-two",
    "mass-balance states: an oral gut depot and a vaginal-tablet depot;",
    "venous and arterial blood, lungs, liver, kidneys, adipose, GI tract,",
    "spleen and a lumped rest-of-body, all perfusion-limited; and a pelvic",
    "sub-model in which the vagina, cervix and uterus each have a luminal",
    "fluid layer, an outer tissue layer (vaginal epithelium, cervical",
    "epithelium, uterine endometrium) and an inner tissue layer (vaginal",
    "stroma, cervical stroma, uterine myometrium). The inner layers exchange",
    "with a shared pelvic blood supply fed from arterial and drained to",
    "venous blood; neighbouring layers exchange by first-order passive",
    "permeation. The vaginal lumen carries a dissolved and a suspended",
    "(precipitated) pool split by a Hill-type solubility switch, and is",
    "cleared by fluid turnover. Elimination is hepatic (plus the vaginal",
    "luminal fluid turnover). Tissue volumes scale with body weight (WT).",
    "Oral and vaginal doses go to their depots and are released first-order",
    "with a common rate constant; IUS release is a zero-order infusion into",
    "the uterine luminal fluid (rate in the dose record). States are AMOUNTS",
    "(mg); outputs are concentrations (mg/L). CAUTION, vaginal route: the",
    "printed equations have no redissolution of suspended drug, so a vaginal",
    "dose above the luminal solubility is almost entirely lost to fluid",
    "turnover and the model does not reproduce the paper's own vaginal",
    "simulation (Figure S4); see the vignette. Deterministic typical-value",
    "model (no IIV, no residual error): the paper's population predictions",
    "come from Monte Carlo sampling of the truncated-normal distributions in",
    "Supplemental Table 7, which the vignette reproduces."
  )
  reference <- paste(
    "Li X, Straubinger T, Rohan LC, Achilles SL, Chen BA, Valicherla G,",
    "Zhang Z, Donnelly M, Tsakalozou E, Zhao L, Bies R. (2026).",
    "Mechanism-Based Predictions of Local Tissue and Systemic Exposure for",
    "Drug Products Delivered Through the Female Reproductive Tract.",
    "J Clin Pharmacol 66(2):e70159. doi:10.1002/jcph.70159. PMCID",
    "PMC12910638. Structure from Supplemental Information S1 (ODEs) and S2",
    "(mrgsolve code); physiological parameters from Supplemental Table S1;",
    "drug-specific parameters from Supplemental Table S3; population",
    "variability from Supplemental Table S7."
  )
  vignette <- "Li_2026_levonorgestrel"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    amount = "mg"
  )

  # The female-reproductive-tract states have no canonical compartment names.
  # Kay_2018_dapivirine_pbpk.R, from the same group and the code base this
  # model's mrgsolve file extends, declares its vaginal states the same way.
  paper_specific_compartments <- c(
    "depot_vaginal",
    "vaginal_lumen",
    "vaginal_lumen_precip",
    "vaginal_epithelium",
    "vaginal_stroma",
    "pelvic_blood",
    "cervical_fluid",
    "cervical_epithelium",
    "cervical_stroma",
    "uterine_fluid",
    "uterine_endometrium",
    "uterine_myometrium"
  )

  covariateData <- list(
    WT = list(
      description = paste(
        "Body weight. Every tissue volume is a fixed fraction of body weight",
        "(Supplemental Table S1); blood flows do not scale with weight."
      ),
      units = "kg",
      type = "continuous",
      reference_category = "58 kg (reference woman, Supplemental Table S1)",
      notes = paste(
        "Linear (not allometric) scaling of volumes only. Supplemental Table",
        "S7 samples WT from a normal distribution with mean 58 kg and SD",
        "8.6 kg truncated to 35-100 kg for the population predictions."
      ),
      source_name = "weight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    depot_vaginal = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    # Vaginal, cervical and uterine luminal fluids: the specimen vocabulary has
    # no genital-tract fluid entry; "epithelial lining fluid" is literally the
    # matrix (as in Kay_2018_dapivirine_pbpk.R).
    vaginal_lumen = list(
      analyte = "levonorgestrel (dissolved)",
      units = "mg",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    vaginal_lumen_precip = list(
      analyte = "levonorgestrel (suspended)",
      units = "mg",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    vaginal_epithelium = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    vaginal_stroma = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    pelvic_blood = list(analyte = "levonorgestrel", units = "mg", specimen = "plasma", verified = TRUE),
    cervical_fluid = list(
      analyte = "levonorgestrel",
      units = "mg",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    cervical_epithelium = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    cervical_stroma = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    uterine_fluid = list(
      analyte = "levonorgestrel",
      units = "mg",
      specimen = "epithelial lining fluid",
      verified = TRUE
    ),
    uterine_endometrium = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    uterine_myometrium = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    # The venous pool is what the paper plots as "plasma concentration".
    venous = list(analyte = "levonorgestrel", units = "mg", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "levonorgestrel", units = "mg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "levonorgestrel", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 11,
    sex_female_pct = 100,
    weight_range = "58 kg reference woman; 35-100 kg in the population predictions",
    disease_state = "healthy women using levonorgestrel contraception or emergency contraception",
    dose_range = paste(
      "Oral 0.15-1.5 mg tablets (single, twice 12 or 24 h apart, or once",
      "daily); vaginal 1.5 mg (2 x 0.75 mg tablets); intrauterine systems",
      "holding 13.5 mg (Skyla), 19.5 mg (Kyleena), 52 mg (Mirena, CHIC",
      "LNG-IUS) and 63 mg (Nilsson)"
    ),
    notes = paste(
      "Not fitted to individual data: a literature-parameterised PBPK model",
      "calibrated against published mean profiles (Supplemental Table S6).",
      "Model development: Nilsson 1982 (oral and IUS; plasma, endometrium,",
      "myometrium, adipose), Devoto 2005 (oral and vaginal; plasma,",
      "endometrium), Landgren 1986, He 1990 (two products), Kook 2002,",
      "Tremblay 2001 group B, and the Zim CHIC cohort (oral COC and LNG-IUS;",
      "plasma and vaginal fluid). External validation: Tremblay 2001 and",
      "Johansson 2002 (two 0.75 mg doses 12 or 24 h apart) and the FDA",
      "clinical-pharmacology reviews of Mirena, Skyla and Kyleena. Uterine",
      "volume and uterine artery flow were measured by ultrasound in 12 women",
      "(Supplemental Table S2); uterine and cervical permeation rate",
      "constants come from ex vivo permeation of excised human tissue",
      "(Supplemental Table S4)."
    )
  )

  ini({
    # ---- Blood flows (Supplemental Table S1; mrgsolve code S2 $MAIN) --------
    q_cardiac <- fixed(294); label("Cardiac output (L/h)") # Table S1, QCO = 294 (70 bpm x 70 mL/beat); S2 Q_CO = 4.9 * 60
    qfrac_liver <- fixed(0.1); label("Hepatic arterial blood flow as a fraction of cardiac output (fraction)") # Table S1, QH = 0.1 * QCO
    qfrac_kidney <- fixed(0.15); label("Renal blood flow as a fraction of cardiac output (fraction)") # Table S1, QR = 0.15 * QCO
    qfrac_adipose <- fixed(0.11); label("Adipose blood flow as a fraction of cardiac output (fraction)") # Table S1, QAD = 0.11 * QCO
    qfrac_gut <- fixed(0.14); label("GI-tract blood flow as a fraction of cardiac output (fraction)") # Table S1, QGI = 0.14 * QCO
    qfrac_spleen <- fixed(0.03); label("Splenic blood flow as a fraction of cardiac output (fraction)") # Table S1, QS = 0.03 * QCO
    q_pelvic <- fixed(2.44); label("Pelvic blood supply flow (L/h)") # Table S1, QVS = 2.44, from measured uterine artery flow (Table S2)

    # ---- Tissue mass fractions of body weight (Supplemental Table S1; S2) ---
    bwfrac_lung <- fixed(0.008); label("Lung mass fraction of body weight (fraction)") # Table S1, VLungs = 0.008 * BM
    bwfrac_liver <- fixed(0.026); label("Liver mass fraction of body weight (fraction)") # Table S1, VLiver = 0.026 * BM
    bwfrac_kidney <- fixed(0.004); label("Kidney mass fraction of body weight (fraction)") # Table S1, VKidneys = 0.004 * BM
    bwfrac_gut <- fixed(0.0165); label("GI-tract mass fraction of body weight (fraction)") # S2 F_Gut = 0.0165; Table S1 prints the rounded 0.017 * BM
    bwfrac_spleen <- fixed(0.0026); label("Spleen mass fraction of body weight (fraction)") # Table S1, VSpleen = 2.6e-3 * BM
    bwfrac_adipose <- fixed(0.327); label("Adipose mass fraction of body weight (fraction)") # Table S1, VAdipose = 0.327 * BM
    bwfrac_arterial <- fixed(0.019); label("Arterial blood mass fraction of body weight (fraction)") # Table S1, VArteries = 0.019 * BM
    bwfrac_venous <- fixed(0.056); label("Venous blood mass fraction of body weight (fraction)") # Table S1, VVeins = 0.056 * BM
    bwfrac_uterus <- fixed(0.00138); label("Uterus mass fraction of body weight (fraction)") # Table S1, VUterus = 1.38e-3 * BM
    bwfrac_uterine_fluid <- fixed(1.72e-6); label("Uterine luminal fluid mass fraction of body weight (fraction)") # Table S1, VUfluid = 1.72e-6 * BM
    bwfrac_cervix <- fixed(0.000431); label("Cervix mass fraction of body weight (fraction)") # Table S1, VCervix = 4.31e-4 * BM
    bwfrac_cervical_fluid <- fixed(8.62e-6); label("Cervical luminal fluid mass fraction of body weight (fraction)") # Table S1, VCvfluid = 8.62e-6 * BM
    bwfrac_vaginal_lumen <- fixed(8.62e-6); label("Vaginal luminal fluid mass fraction of body weight (fraction)") # Table S1, VLumen = 8.62e-6 * BM
    bwfrac_vaginal_epithelium <- fixed(0.000151); label("Vaginal epithelium mass fraction of body weight (fraction)") # Table S1, VEpithelium = 1.51e-4 * BM
    bwfrac_vaginal_stroma <- fixed(0.000181); label("Vaginal stroma mass fraction of body weight (fraction)") # Table S1, VStroma_Tissue = 1.81e-4 * BM
    bwfrac_pelvic_blood <- fixed(0.00259); label("Pelvic blood supply mass fraction of body weight (fraction)") # Table S1, VStroma_Blood = 2.59e-3 * BM

    # ---- Drug-specific physicochemistry (Supplemental Table S3) -------------
    fu <- fixed(0.015); label("Fraction unbound in the genital-tract luminal fluids (fraction)") # Table S3, funbound = 0.015 (S2 code carries 0.02)
    bpr <- fixed(1.07); label("Blood-to-plasma concentration ratio (ratio)") # Table S3, BP = 1.07, predicted
    funi_lumen <- fixed(1); label("Un-ionised fraction in the vaginal luminal fluid (fraction)") # Table S3, funi = 1 (LNG neutral)
    funi_tissue <- fixed(1); label("Un-ionised fraction in the vaginal epithelium and stroma (fraction)") # Table S3, funibo = 1 (LNG neutral)
    sol_lumen <- fixed(2.05); label("Aqueous solubility of LNG in the vaginal luminal fluid (mg/L)") # Table S3, solubility = 2.05 mg/L
    hill_sol <- fixed(5); label("Hill exponent of the vaginal-lumen solubility switch (unitless)") # Main text Eq. 9 text, gamma = 5; S2 gamma = 5.0

    # ---- Clearance (Supplemental Table S3) -----------------------------------
    lcl <- fixed(log(7.06)); label("Systemic (hepatic) clearance (L/h)") # Table S3, CLH = 7.06 L/h
    lcl_lumen <- fixed(log(2.5e-5)); label("Vaginal luminal fluid-turnover clearance (L/h)") # Table S3, CLL = 2.5e-5 L/h

    # ---- Tissue-to-plasma partition coefficients (Supplemental Table S3) ----
    lkp_liver <- fixed(log(2.1)); label("Liver-to-plasma partition coefficient (unitless)") # Table S3, PH = 2.1
    lkp_kidney <- fixed(log(1.9)); label("Kidney-to-plasma partition coefficient (unitless)") # Table S3, PR = 1.9
    lkp_lung <- fixed(log(2.5)); label("Lung-to-plasma partition coefficient (unitless)") # Table S3, PLungs = 2.5
    lkp_adipose <- fixed(log(5)); label("Adipose-to-plasma partition coefficient (unitless)") # Table S3, PAD = 5
    lkp_gut <- fixed(log(4)); label("GI-tract-to-plasma partition coefficient (unitless)") # Table S3, PGI = 4.0
    lkp_spleen <- fixed(log(1.1)); label("Spleen-to-plasma partition coefficient (unitless)") # S2 code P_S = 1.1 (not in Table S3)
    lkp_other <- fixed(log(0.6)); label("Rest-of-body-to-plasma partition coefficient (unitless)") # Table S3, PB = 0.6
    lkp_vaginal_stroma <- fixed(log(0.01)); label("Vaginal-stroma-to-plasma partition coefficient (unitless)") # Table S3, PVS = 0.01
    lkp_cervical_stroma <- fixed(log(0.01)); label("Cervical-stroma-to-plasma partition coefficient (unitless)") # Table S3, PCM = 0.01
    lkp_uterine_myometrium <- fixed(log(0.01)); label("Uterine-myometrium-to-plasma partition coefficient (unitless)") # Table S3, PUM = 0.01

    # ---- First-order permeation rate constants (1/h; Supplemental Table S3) -
    # Naming: k_<from>_<to>; vl/ve/vs = vaginal lumen/epithelium/stroma,
    # cf/ce/cm = cervical fluid/epithelium/stroma, uf/ue/um = uterine
    # fluid/endometrium/myometrium, pb = pelvic blood supply.
    k_vl_ve <- fixed(5); label("Vaginal luminal fluid to epithelium (1/h)") # Table S3, kLE = 5
    k_ve_vl <- fixed(0.01); label("Vaginal epithelium to luminal fluid (1/h)") # Table S3, kEL = 0.01
    k_ve_vs <- fixed(0.43); label("Vaginal epithelium to stroma (1/h)") # Table S3, kEVS = 0.43
    k_ve_vs_active <- fixed(0); label("Active transport of ionised drug, vaginal epithelium to stroma (1/h)") # Table S3, kactiveEtoS = 0
    k_vs_ve <- fixed(0.01); label("Vaginal stroma to epithelium (1/h)") # Table S3, kVSE = 0.01
    k_vs_pb_active <- fixed(0); label("Active transport of ionised drug, vaginal stroma to pelvic blood (1/h)") # Table S3, kactiveStrTtoB = 0
    k_vs_pb <- fixed(0.15); label("Vaginal stroma to pelvic blood (1/h)") # Table S3, kVS_TtoB = 0.15
    k_pb_vs <- fixed(0.2); label("Pelvic blood to vaginal stroma (1/h)") # Table S3, kVS_BtoT = 0.2
    k_vl_cf <- fixed(1e-4); label("Vaginal to cervical luminal fluid (1/h)") # Table S3, kL_CF = 1e-4
    k_cf_vl <- fixed(1e-4); label("Cervical to vaginal luminal fluid (1/h)") # Table S3, kCF_L = 1e-4
    k_ve_ce <- fixed(1e-3); label("Vaginal to cervical epithelium (1/h)") # Table S3, kVE_CE = 1e-3
    k_ce_ve <- fixed(1e-4); label("Cervical to vaginal epithelium (1/h)") # Table S3, kCE_VE = 1e-4
    k_vs_cm <- fixed(1e-4); label("Vaginal stroma to cervical stroma (1/h)") # Table S3, kVS_CM = 1e-4
    k_cm_vs <- fixed(1e-4); label("Cervical stroma to vaginal stroma (1/h)") # Table S3, kCM_VS = 1e-4
    k_cf_ce <- fixed(2); label("Cervical luminal fluid to epithelium (1/h)") # Table S3, kCF_CE = 2
    k_ce_cf <- fixed(0.212); label("Cervical epithelium to luminal fluid (1/h)") # Table S3, kCE_CF = 0.212
    k_ce_cm <- fixed(0.043); label("Cervical epithelium to stroma (1/h)") # Table S3, kCE_CM = 0.043
    k_cm_ce <- fixed(0.212); label("Cervical stroma to epithelium (1/h)") # Table S3, kCM_CE = 0.212
    k_cm_pb <- fixed(0.043); label("Cervical stroma to pelvic blood (1/h)") # Table S3, kCM_TtoB = 0.043
    k_pb_cm <- fixed(1); label("Pelvic blood to cervical stroma (1/h)") # Table S3, kCM_BtoT = 1
    k_cf_uf <- fixed(1e-4); label("Cervical to uterine luminal fluid (1/h)") # Table S3, kCF_UF = 1e-4
    k_uf_cf <- fixed(1e-4); label("Uterine to cervical luminal fluid (1/h)") # Table S3, kUF_CF = 1e-4
    k_ce_ue <- fixed(1e-2); label("Cervical epithelium to uterine endometrium (1/h)") # Table S3, kCE_UE = 1e-2
    k_ue_ce <- fixed(1e-4); label("Uterine endometrium to cervical epithelium (1/h)") # Table S3, kUE_CE = 1e-4
    k_cm_um <- fixed(1e-4); label("Cervical stroma to uterine myometrium (1/h)") # Table S3, kCM_UM = 1e-4
    k_um_cm <- fixed(1e-4); label("Uterine myometrium to cervical stroma (1/h)") # Table S3, kUM_CM = 1e-4
    k_uf_ue <- fixed(2); label("Uterine luminal fluid to endometrium (1/h)") # Table S3, kUF_UE = 2
    k_ue_uf <- fixed(0.148); label("Uterine endometrium to luminal fluid (1/h)") # Table S3, kUE_UF = 0.148
    k_ue_um <- fixed(0.076); label("Uterine endometrium to myometrium (1/h)") # Table S3, kUE_UM = 0.076
    k_um_ue <- fixed(0.097); label("Uterine myometrium to endometrium (1/h)") # Table S3, kUM_UE = 0.097
    k_um_pb <- fixed(0.076); label("Uterine myometrium to pelvic blood (1/h)") # Table S3, kUM_TtoB = 0.076
    k_pb_um <- fixed(3); label("Pelvic blood to uterine myometrium (1/h)") # Table S3, kUM_BtoT = 3

    # ---- Dosing function (mrgsolve code S2 $PARAM) --------------------------
    lka <- fixed(log(1.5)); label("First-order release rate constant from the oral and vaginal depots (1/h)") # S2 TVka = 1.5
    lfdepot <- fixed(log(1)); label("Bioavailability of the oral and vaginal doses (fraction)") # S2 F = 1.0
  })

  model({
    # ---- 1. Derived terms --------------------------------------------------
    cl <- exp(lcl)
    cl_lumen <- exp(lcl_lumen)
    ka <- exp(lka)
    kp_liver <- exp(lkp_liver)
    kp_kidney <- exp(lkp_kidney)
    kp_lung <- exp(lkp_lung)
    kp_adipose <- exp(lkp_adipose)
    kp_gut <- exp(lkp_gut)
    kp_spleen <- exp(lkp_spleen)
    kp_other <- exp(lkp_other)
    kp_vaginal_stroma <- exp(lkp_vaginal_stroma)
    kp_cervical_stroma <- exp(lkp_cervical_stroma)
    kp_uterine_myometrium <- exp(lkp_uterine_myometrium)

    # Blood flows (L/h). Rest-of-body flow closes the cardiac-output balance,
    # the formula Table S1 prints beside its value (see the vignette Errata).
    q_liver <- qfrac_liver * q_cardiac
    q_kidney <- qfrac_kidney * q_cardiac
    q_adipose <- qfrac_adipose * q_cardiac
    q_gut <- qfrac_gut * q_cardiac
    q_spleen <- qfrac_spleen * q_cardiac
    q_other <- q_cardiac - q_pelvic - q_liver - q_kidney - q_adipose -
      q_gut - q_spleen

    # Tissue volumes (L), density 1 kg/L except adipose at 0.916 kg/L. The
    # endometrium:myometrium split 0.6:1.0 follows their relative thickness,
    # and the cervical epithelium is 1/2.2 of the cervix (Table S1).
    v_lung <- bwfrac_lung * WT
    v_liver <- bwfrac_liver * WT
    v_kidney <- bwfrac_kidney * WT
    v_gut <- bwfrac_gut * WT
    v_spleen <- bwfrac_spleen * WT
    v_adipose <- bwfrac_adipose * WT / 0.916
    v_arterial <- bwfrac_arterial * WT
    v_venous <- bwfrac_venous * WT
    v_uterus <- bwfrac_uterus * WT
    v_uterine_fluid <- bwfrac_uterine_fluid * WT
    v_uterine_endometrium <- 0.6 / 1.6 * v_uterus
    v_uterine_myometrium <- 1.0 / 1.6 * v_uterus
    v_cervix <- bwfrac_cervix * WT
    v_cervical_fluid <- bwfrac_cervical_fluid * WT
    v_cervical_epithelium <- v_cervix / 2.2
    v_cervical_stroma <- v_cervix - v_cervical_epithelium
    v_vaginal_lumen <- bwfrac_vaginal_lumen * WT
    v_vaginal_epithelium <- bwfrac_vaginal_epithelium * WT
    v_vaginal_stroma <- bwfrac_vaginal_stroma * WT
    v_pelvic_blood <- bwfrac_pelvic_blood * WT
    # S2 F_Body: the mass-fraction balance uses the adipose MASS fraction.
    bwfrac_other <- 1 - bwfrac_lung - bwfrac_liver - bwfrac_kidney -
      bwfrac_adipose - bwfrac_spleen - bwfrac_arterial - bwfrac_venous -
      bwfrac_uterus - bwfrac_cervix - bwfrac_vaginal_lumen -
      bwfrac_vaginal_epithelium - bwfrac_vaginal_stroma - bwfrac_pelvic_blood -
      bwfrac_uterine_fluid - bwfrac_cervical_fluid - bwfrac_gut
    v_other <- bwfrac_other * WT

    # Concentrations (mg/L) from the amount states.
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_lung <- lung / v_lung
    c_liver <- liver / v_liver
    c_kidney <- kidney / v_kidney
    c_adipose <- adipose / v_adipose
    c_gut <- gut / v_gut
    c_spleen <- spleen / v_spleen
    c_other <- other / v_other
    c_pelvic_blood <- pelvic_blood / v_pelvic_blood
    c_vaginal_total <- (vaginal_lumen + vaginal_lumen_precip) / v_vaginal_lumen

    # Main-text Eq. 9: Hill-type switch sending released vaginal drug to the
    # suspended pool once the total luminal concentration passes the
    # solubility (half-maximal at solubility / 2, as in S2).
    sol_switch <- c_vaginal_total^hill_sol /
      ((sol_lumen / 2)^hill_sol + c_vaginal_total^hill_sol)

    # ---- 2. Mass fluxes (mg/h), S2 $ODE Flux1-Flux43 ----------------------
    # Each S2 flux is k * C * V; with amount states the volume cancels.
    rel_oral <- ka * depot
    rel_vaginal <- ka * depot_vaginal

    # Vaginal tissues (S1 dC_VL, dC_VE, dC_VS; S2 Flux1-3, 16, 17, 37, 38)
    j_vl_ve <- k_vl_ve * vaginal_lumen * fu * funi_lumen -
      k_ve_vl * vaginal_epithelium * funi_tissue
    j_ve_vs <- k_ve_vs * vaginal_epithelium * funi_tissue -
      k_vs_ve * vaginal_stroma * funi_tissue
    j_vs_pb <- k_vs_pb * bpr * vaginal_stroma * funi_tissue / kp_vaginal_stroma -
      k_pb_vs * pelvic_blood
    j_ve_vs_active <- k_ve_vs_active * vaginal_epithelium * (1 - funi_tissue)
    j_vs_pb_active <- k_vs_pb_active * vaginal_stroma * (1 - funi_tissue)
    j_vl_out <- cl_lumen * vaginal_lumen / v_vaginal_lumen
    j_vlp_out <- cl_lumen * vaginal_lumen_precip / v_vaginal_lumen

    # Cervical exchanges (S1 dC_CF, dC_CE, dC_CM; S2 Flux19-24). S1 prints
    # V_VE where the cervical-epithelium volume belongs in the CE back-flux
    # terms; S2 and mass balance use the cervical-epithelium volume.
    j_vl_cf <- k_vl_cf * vaginal_lumen * fu - k_cf_vl * cervical_fluid
    j_ve_ce <- k_ve_ce * vaginal_epithelium - k_ce_ve * cervical_epithelium
    j_vs_cm <- k_vs_cm * vaginal_stroma - k_cm_vs * cervical_stroma
    j_cf_ce <- k_cf_ce * cervical_fluid - k_ce_cf * cervical_epithelium
    j_ce_cm <- k_ce_cm * cervical_epithelium - k_cm_ce * cervical_stroma
    j_cm_pb <- k_cm_pb * bpr * cervical_stroma / kp_cervical_stroma -
      k_pb_cm * pelvic_blood

    # Uterine exchanges (S1 dC_UF, dC_UE, dC_UM; S2 Flux25-30)
    j_cf_uf <- k_cf_uf * cervical_fluid - k_uf_cf * uterine_fluid * fu
    j_ce_ue <- k_ce_ue * cervical_epithelium - k_ue_ce * uterine_endometrium
    j_cm_um <- k_cm_um * cervical_stroma - k_um_cm * uterine_myometrium
    j_uf_ue <- k_uf_ue * uterine_fluid * fu - k_ue_uf * uterine_endometrium
    j_ue_um <- k_ue_um * uterine_endometrium - k_um_ue * uterine_myometrium
    j_um_pb <- k_um_pb * bpr * uterine_myometrium / kp_uterine_myometrium -
      k_pb_um * pelvic_blood

    # Perfusion-limited outflows: Q * BP * C_tissue / Kp (main-text Eq. 1).
    out_lung <- q_cardiac * bpr * c_lung / kp_lung
    out_liver <- (q_liver + q_gut + q_spleen) * bpr * c_liver / kp_liver
    out_kidney <- q_kidney * bpr * c_kidney / kp_kidney
    out_adipose <- q_adipose * bpr * c_adipose / kp_adipose
    out_gut <- q_gut * bpr * c_gut / kp_gut
    out_spleen <- q_spleen * bpr * c_spleen / kp_spleen
    out_other <- q_other * bpr * c_other / kp_other
    elim_liver <- cl * bpr * c_liver / kp_liver

    # ---- 3. ODE system (amounts, mg) ---------------------------------------
    d/dt(depot) <- -rel_oral
    d/dt(depot_vaginal) <- -rel_vaginal

    # Main-text Eqs. 7-8: released vaginal drug splits between solution and
    # suspension; only dissolved, un-ionised, unbound drug permeates.
    d/dt(vaginal_lumen) <- rel_vaginal * (1 - sol_switch) - j_vl_ve -
      j_vl_cf - j_vl_out
    d/dt(vaginal_lumen_precip) <- rel_vaginal * sol_switch - j_vlp_out
    d/dt(vaginal_epithelium) <- j_vl_ve - j_ve_vs - j_ve_vs_active - j_ve_ce
    d/dt(vaginal_stroma) <- j_ve_vs + j_ve_vs_active - j_vs_pb -
      j_vs_pb_active - j_vs_cm
    d/dt(pelvic_blood) <- j_vs_pb - q_pelvic * c_pelvic_blood +
      q_pelvic * c_arterial + j_vs_pb_active + j_cm_pb + j_um_pb

    d/dt(cervical_fluid) <- j_vl_cf - j_cf_ce - j_cf_uf
    d/dt(cervical_epithelium) <- j_ve_ce + j_cf_ce - j_ce_cm - j_ce_ue
    d/dt(cervical_stroma) <- j_vs_cm + j_ce_cm - j_cm_pb - j_cm_um

    # IUS release enters here as a zero-order infusion (D_formulation in S1).
    d/dt(uterine_fluid) <- j_cf_uf - j_uf_ue
    d/dt(uterine_endometrium) <- j_ce_ue + j_uf_ue - j_ue_um
    d/dt(uterine_myometrium) <- j_cm_um + j_ue_um - j_um_pb

    # Main-text Eqs. 1-5 and S1: systemic circulation. The GI tract and
    # spleen drain into the liver; the lungs see the whole cardiac output.
    d/dt(venous) <- q_pelvic * c_pelvic_blood + out_other + out_liver +
      out_kidney + out_adipose - q_cardiac * c_venous
    d/dt(lung) <- q_cardiac * c_venous - out_lung
    d/dt(arterial) <- out_lung - (q_pelvic + q_other + q_liver + q_kidney +
      q_adipose + q_gut + q_spleen) * c_arterial
    d/dt(liver) <- q_liver * c_arterial + out_gut + out_spleen - out_liver -
      elim_liver
    d/dt(kidney) <- q_kidney * c_arterial - out_kidney
    d/dt(adipose) <- q_adipose * c_arterial - out_adipose
    d/dt(gut) <- rel_oral + q_gut * c_arterial - out_gut
    d/dt(spleen) <- q_spleen * c_arterial - out_spleen
    d/dt(other) <- q_other * c_arterial - out_other

    # ---- 4. Bioavailability (S2 F; dose amount x F enters the depot) -------
    f(depot) <- exp(lfdepot)
    f(depot_vaginal) <- exp(lfdepot)

    # ---- 5. Observations (mg/L; tissue mg/kg at density 1) -----------------
    Cc <- c_venous
    Cvaginal_fluid <- c_vaginal_total
    Cvaginal_epithelium <- vaginal_epithelium / v_vaginal_epithelium
    Cvaginal_stroma <- vaginal_stroma / v_vaginal_stroma
    Ccervical_fluid <- cervical_fluid / v_cervical_fluid
    Ccervical_epithelium <- cervical_epithelium / v_cervical_epithelium
    Ccervical_stroma <- cervical_stroma / v_cervical_stroma
    Cuterine_fluid <- uterine_fluid / v_uterine_fluid
    Cendometrium <- uterine_endometrium / v_uterine_endometrium
    Cmyometrium <- uterine_myometrium / v_uterine_myometrium
    Cadipose <- c_adipose
  })
}
