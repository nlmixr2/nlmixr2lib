Berton_2024_pregnancy_physiology <- function() {
  description <- paste(
    "PBPK system-parameter (physiology) model of healthy pregnancy and",
    "the postpartum period (Berton 2024 Open Forum Infect Dis), the",
    "pregnant / postpartum virtual population of the Basel in-house",
    "Matlab whole-body PBPK framework. This is the SYSTEM layer only: it",
    "has no drug, no dosing, no compartments and no ODEs, and every",
    "output is an algebraic function of the rxode2 time variable, which",
    "the model interprets as gestational age in WEEKS, with delivery at",
    "40 weeks and the postpartum period beyond it (valid 0 to 66 weeks).",
    "Each woman's non-pregnant anatomy and physiology come from the",
    "framework's age- and height-dependent regressions with their",
    "population CVs (Supplementary Table 2, citing Stader 2019); the",
    "paper's own gestational-age fold-change regressions then scale body",
    "weight, adipose / brain / heart / kidney / uterus / mammary /",
    "plasma / red-cell weights, cardiac output, regional blood flows,",
    "GFR, haematocrit, albumin, alpha-1-acid glycoprotein and the",
    "hepatic and intestinal abundance of CYP1A2, CYP2B6, CYP2C9,",
    "CYP2C19, CYP2D6, CYP3A4, CYP3A5, UGT1A1 and UGT1A4. Fetal,",
    "placental and amniotic-fluid weights and the fetoplacental blood",
    "flow are added for the pregnancy. The seven antiretroviral drug",
    "PBPK models the paper verifies against this population (ritonavir,",
    "efavirenz, rilpivirine, doravirine, dolutegravir, bictegravir,",
    "raltegravir) are NOT included: their drug parameters are inherited",
    "from earlier publications and not printed here, and the way the",
    "fetoplacental unit is wired into the PBPK ODEs is not described.",
    sep = " "
  )
  reference <- paste(
    "Berton M, Stader F, Bettonte S, Battegay M, Marzolini C.",
    "Ideal Time to Conduct a Pharmacokinetic Investigation After Delivery",
    "to Fully Capture the Effect of Pregnancy on Drug Exposure.",
    "Open Forum Infect Dis. 2024;11(10):ofae585.",
    "doi:10.1093/ofid/ofae585.",
    "Equations: Supplementary Table 2 (supplementary data file 1).",
    sep = " "
  )
  vignette <- "Berton_2024_pregnancy_physiology"
  units <- list(
    time = "week (gestational age; delivery at 40; valid 0 to 66)",
    dosing = "n/a (no exogenous dosing; PBPK system-physiology model)",
    concentration = paste(
      "n/a (no drug concentration). Outputs carry their own units:",
      "height cm; body and organ weights kg; BMI kg/m^2; BSA m^2;",
      "cardiac output and absolute blood flows L/h; regional flows",
      "fq_* percent of cardiac output; GFR mL/min; haematocrit as a",
      "fraction; albumin (hsa) and alpha-1-acid glycoprotein (agp) g/L;",
      "enzyme-abundance fold changes fc_* unitless relative to the",
      "non-pregnant woman."
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age of the woman",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the non-pregnant regressions of Supplementary Table 2:",
        "height, non-pregnant body weight, lung, brain, gonad, heart,",
        "kidney, muscle, skin and liver weights, blood weight, cardiac",
        "output, adipose / liver / kidney blood-flow fractions, GFR and",
        "albumin. The generator draws age from a uniform distribution",
        "(Supplementary Table 2, row 'Age'); the range is set per",
        "simulated trial (Supplementary Table 3 lists 18 to 44 years for",
        "the verification trials). Age is held constant over the",
        "simulated gestational / postpartum window."
      ),
      source_name = "Age"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1000L,
    n_studies = 224L,
    age_range = "Uniform draw per simulated trial; verification trials 18 to 44 years (Supplementary Table 3)",
    weight_range = "Generated from the non-pregnant regression (CV 15.2%) times the gestational-age fold change; Figure 1A spans about 40 to 120 kg",
    sex_female_pct = 100,
    race_ethnicity = "Not reported; the non-pregnant regressions come from the Stader 2019 ageing repository (White adults)",
    disease_state = "Healthy pregnant and postpartum women (physiology development data); the verification trials enrolled pregnant women living with HIV",
    dose_range = "n/a (no drug administered)",
    gestational_age_range = "0 to 40 weeks pregnancy and 40 to 66 weeks gestational age postpartum (0 to 26 weeks after delivery)",
    regions = "Not restricted (PubMed and Google Scholar literature search)",
    notes = paste(
      "224 articles (Supplementary Table 1) were split into development",
      "and verification data sets; parameter values were normalised to",
      "non-pregnant values and regressed on gestational age (linear,",
      "polynomial or exponential, selected by R^2 and AIC) in Matlab",
      "2020a. The fitted equations were implemented in the in-house PBPK",
      "framework and 1000 virtual women were generated across 0 to 66",
      "weeks (Methods; Figure 1; Supplementary Figures 1 and 2). All",
      "women are assumed to deliver at exactly 40 weeks."
    )
  )

  # The etas are the per-subject random draws of the virtual-population
  # generator (the CV column of Supplementary Table 2), not estimated
  # between-subject variances.
  paper_specific_etas <- c(
    "etaht",
    "etawt",
    "etaw_adipose",
    "etaw_bone",
    "etaw_brain",
    "etaw_gonads",
    "etaw_heart",
    "etaw_kidney",
    "etaw_uterus",
    "etaw_mammary",
    "etaw_muscle",
    "etaw_skin",
    "etaw_thymus",
    "etaw_gut",
    "etaw_spleen",
    "etaw_pancreas",
    "etaw_liver",
    "etaw_blood",
    "etaw_fetus",
    "etaw_placenta",
    "etaw_amniotic",
    "etaco",
    "etagfr",
    "etahct",
    "etahsa",
    "etaagp"
  )

  ini({
    # ---------------------------------------------------------------
    # Gestational-age (GA, weeks) fold-change (FC) regressions fitted by
    # the authors (Supplementary Table 2 rows without a literature
    # citation). Each FC multiplies the non-pregnant value. Pregnancy
    # rows (GA <= 40 unless stated) have an intercept of exactly 1, so
    # only the GA coefficients are parameters; postpartum rows (GA > 40)
    # carry their own intercept. No uncertainty is reported for any
    # coefficient, so all are fixed().
    # ---------------------------------------------------------------

    # ---- Body weight FC ----
    bw_preg_c2 <- fixed(4.85e-05)
    label("Body weight FC, pregnancy, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Body weight FC, Pregnant, <= 40': 4.85E-05*GesAge^2 + 4.05E-03*GesAge + 1
    bw_preg_c1 <- fixed(4.05e-03)
    label("Body weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Body weight FC, Pregnant, <= 40'
    bw_pp_c2 <- fixed(3.38e-06)
    label("Body weight FC, postpartum, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Body weight FC, Postpartum, > 40': 3.38E-06*GesAge^2 - 1.28E-03*GesAge + 1.11
    bw_pp_c1 <- fixed(-1.28e-03)
    label("Body weight FC, postpartum, linear GA coefficient (1/week)") # Suppl. Table 2 'Body weight FC, Postpartum, > 40'
    bw_pp_c0 <- fixed(1.11)
    label("Body weight FC, postpartum, intercept (unitless)") # Suppl. Table 2 'Body weight FC, Postpartum, > 40'

    # ---- Organ weight FCs ----
    adipose_preg_c1 <- fixed(6.75e-03)
    label("Adipose weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Adipose weight FC, Pregnant, <= 40': 6.75E-03*GesAge + 1
    adipose_pp_c1 <- fixed(-1.88e-03)
    label("Adipose weight FC, postpartum, linear GA coefficient (1/week)") # Suppl. Table 2 'Adipose weight FC, Postpartum, > 40': -1.88E-03*GesAge + 1.24
    adipose_pp_c0 <- fixed(1.24)
    label("Adipose weight FC, postpartum, intercept (unitless)") # Suppl. Table 2 'Adipose weight FC, Postpartum, > 40'
    brain_preg_c1 <- fixed(-1.49e-03)
    label("Brain weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Brain weight FC, Pregnant, <= 40': -1.49E-03*GesAge + 1
    brain_pp_c1 <- fixed(1.22e-03)
    label("Brain weight FC, postpartum (40 to 42 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Brain weight FC, Postpartum, > 40': 1.22E-03*GesAge + 0.95
    brain_pp_c0 <- fixed(0.95)
    label("Brain weight FC, postpartum (40 to 42 weeks), intercept (unitless)") # Suppl. Table 2 'Brain weight FC, Postpartum, > 40'
    heart_preg_c1 <- fixed(6.74e-03)
    label("Heart weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Heart weight FC, Pregnant, <= 40': 6.74E-03*GesAge + 1
    heart_pp_c1 <- fixed(-1.35e-03)
    label("Heart weight FC, postpartum (40 to 86 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Heart weight FC, Postpartum, > 40': -1.35E-03*GesAge + 1.12
    heart_pp_c0 <- fixed(1.12)
    label("Heart weight FC, postpartum (40 to 86 weeks), intercept (unitless)") # Suppl. Table 2 'Heart weight FC, Postpartum, > 40'
    kidney_preg_c1 <- fixed(7.32e-03)
    label("Kidney weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Kidney weight FC, Pregnant, <= 40': 7.32E-03*GesAge + 1
    kidney_pp_c1 <- fixed(-2.93e-02)
    label("Kidney weight FC, postpartum (40 to 50 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Kidney weight FC, Postpartum, > 40': -2.93E-02*GesAge + 2.46
    kidney_pp_c0 <- fixed(2.46)
    label("Kidney weight FC, postpartum (40 to 50 weeks), intercept (unitless)") # Suppl. Table 2 'Kidney weight FC, Postpartum, > 40'

    # ---- Plasma and red-blood-cell weight FCs ----
    plasma_preg_c1 <- fixed(1.43e-02)
    label("Plasma weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Plasma weight FC, Pregnant, <= 40': 1.43E-02*GesAge + 1
    plasma_pp_c1 <- fixed(-4.16e-03)
    label("Plasma weight FC, postpartum (40 to 57.7 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Plasma weight FC, Postpartum, > 40': -4.16E-03*GesAge + 1.24
    plasma_pp_c0 <- fixed(1.24)
    label("Plasma weight FC, postpartum (40 to 57.7 weeks), intercept (unitless)") # Suppl. Table 2 'Plasma weight FC, Postpartum, > 40'
    rbc_preg_c2 <- fixed(2.47e-04)
    label("Red-blood-cell weight FC, pregnancy, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Red blood cells weight FC, Pregnant, <= 40': 2.47E-04*GesAge^2 + 2.49E-04*GesAge + 1
    rbc_preg_c1 <- fixed(2.49e-04)
    label("Red-blood-cell weight FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Red blood cells weight FC, Pregnant, <= 40'
    rbc_pp_c2 <- fixed(2.07e-03)
    label("Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Red blood cells weight FC, Postpartum, > 40': 2.07E-03*GesAge^2 - 1.95E-01*GesAge + 5.64
    rbc_pp_c1 <- fixed(-1.95e-01)
    label("Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Red blood cells weight FC, Postpartum, > 40'
    rbc_pp_c0 <- fixed(5.64)
    label("Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), intercept (unitless)") # Suppl. Table 2 'Red blood cells weight FC, Postpartum, > 40'

    # ---- Cardiac output and blood-flow FCs ----
    co_c2 <- fixed(-4.77e-04)
    label("Cardiac output FC (0 to 50.9 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Cardiac output FC, Pregnant and postpartum, <= 50.9': -4.77E-04*GesAge^2 + 2.43E-02*GesAge + 1
    co_c1 <- fixed(2.43e-02)
    label("Cardiac output FC (0 to 50.9 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Cardiac output FC, Pregnant and postpartum, <= 50.9'
    qfrac_preg_c1 <- fixed(-0.0099)
    label("FC on the percent-of-cardiac-output flows of organs with unchanged absolute flow, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Blood flows FC, Pregnant, <= 40': -0.0099*GesAge + 1
    qfrac_pp_c1 <- fixed(0.0073)
    label("FC on the percent-of-cardiac-output flows, postpartum (40 to 55 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Blood flows FC, Postpartum, > 40': 0.0073*GesAge + 0.595
    qfrac_pp_c0 <- fixed(0.595)
    label("FC on the percent-of-cardiac-output flows, postpartum (40 to 55 weeks), intercept (unitless)") # Suppl. Table 2 'Blood flows FC, Postpartum, > 40'
    qbrain_c2 <- fixed(-1.24e-04)
    label("Absolute brain blood flow FC (0 to 82.4 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Brain blood flow FC, <= 82.4': (-1.24E-04*GesAge^2 + 1.13E-02*GesAge + 1) / Cardiac Output FC
    qbrain_c1 <- fixed(1.13e-02)
    label("Absolute brain blood flow FC (0 to 82.4 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Brain blood flow FC, <= 82.4'
    qkidney_c2 <- fixed(-1.03e-03)
    label("Absolute kidney blood flow FC (0 to 40.7 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Kidney blood flow FC, Pregnant, <= 40.7': (-1.03E-03*GesAge^2 + 4.19E-02*GesAge + 1) / Cardiac Output FC
    qkidney_c1 <- fixed(4.19e-02)
    label("Absolute kidney blood flow FC (0 to 40.7 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Kidney blood flow FC, Pregnant, <= 40.7'
    quterus_preg_c2 <- fixed(2.62e-03)
    label("Uterus blood flow FC, pregnancy, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Uterus blood flow FC, Pregnant, <= 40': 2.62E-03*GesAge^2 + 1.83E-01*GesAge + 1
    quterus_preg_c1 <- fixed(1.83e-01)
    label("Uterus blood flow FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Uterus blood flow FC, Pregnant, <= 40'
    quterus_pp_c1 <- fixed(-0.66)
    label("Uterus blood flow FC, postpartum (40 to 53 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Uterus blood flow FC, Postpartum, > 40': -0.66*GesAge + 36
    quterus_pp_c0 <- fixed(36)
    label("Uterus blood flow FC, postpartum (40 to 53 weeks), intercept (unitless)") # Suppl. Table 2 'Uterus blood flow FC, Postpartum, > 40'

    # ---- GFR, haematocrit and plasma-protein FCs ----
    gfr_c2 <- fixed(-2.33e-04)
    label("GFR FC (0 to 67.5 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'GFR FC, <= 67.5': -2.33E-04*GesAge^2 + 1.57E-02*GesAge + 1
    gfr_c1 <- fixed(1.57e-02)
    label("GFR FC (0 to 67.5 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'GFR FC, <= 67.5'
    hct_c2 <- fixed(8.31e-05)
    label("Haematocrit FC (0 to 65 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'HCT FC, <= 65': 8.31E-05*GesAge^2 - 5.37E-03*GesAge + 1
    hct_c1 <- fixed(-5.37e-03)
    label("Haematocrit FC (0 to 65 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'HCT FC, <= 65'
    hsa_preg_c2 <- fixed(1.75e-04)
    label("Albumin FC, pregnancy, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'Albumin FC, Pregnant, <= 40': 1.75E-04*GesAge^2 - 1.32E-02*GesAge + 1
    hsa_preg_c1 <- fixed(-1.32e-02)
    label("Albumin FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'Albumin FC, Pregnant, <= 40'
    hsa_pp_c1 <- fixed(0.05)
    label("Albumin FC, postpartum (40 to 45 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'Albumin FC, Postpartum, > 40': 0.05*GesAge - 1.23
    hsa_pp_c0 <- fixed(-1.23)
    label("Albumin FC, postpartum (40 to 45 weeks), intercept (unitless)") # Suppl. Table 2 'Albumin FC, Postpartum, > 40'
    agp_c2 <- fixed(2.07e-04)
    label("Alpha-1-acid glycoprotein FC (0 to 63.7 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'AAG FC, <= 63.7': 2.07E-04*GesAge^2 - 1.32E-02*GesAge + 1
    agp_c1 <- fixed(-1.32e-02)
    label("Alpha-1-acid glycoprotein FC (0 to 63.7 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'AAG FC, <= 63.7'

    # ---- Hepatic (and intestinal) enzyme-abundance FCs ----
    cyp1a2_c2 <- fixed(9.54e-04)
    label("CYP1A2 abundance FC (0 to 51.8 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP1A2 FC, <= 51.8': 9.54E-04*GesAge^2 - 4.94E-02*GesAge + 1
    cyp1a2_c1 <- fixed(-4.94e-02)
    label("CYP1A2 abundance FC (0 to 51.8 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP1A2 FC, <= 51.8'
    cyp2b6_c2 <- fixed(-7.36e-04)
    label("CYP2B6 abundance FC (0 to 50.8 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP2B6 FC, <= 50.8': -7.36E-04*GesAge^2 + 3.73E-02*GesAge + 1
    cyp2b6_c1 <- fixed(3.73e-02)
    label("CYP2B6 abundance FC (0 to 50.8 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP2B6 FC, <= 50.8'
    cyp2c9_c2 <- fixed(-5.26e-03)
    label("CYP2C9 abundance FC (0 to 42 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP2C9 FC, <= 42': -5.26E-03*GesAge^2 + 2.18E-01*GesAge + 1
    cyp2c9_c1 <- fixed(2.18e-01)
    label("CYP2C9 abundance FC (0 to 42 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP2C9 FC, <= 42'
    cyp2c19_c2 <- fixed(-1.75e-03)
    label("CYP2C19 abundance FC (0 to 43.5 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP2C19 FC, <= 43.5': -1.75E-03*GesAge^2 + 7.61E-02*GesAge + 1
    cyp2c19_c1 <- fixed(7.61e-02)
    label("CYP2C19 abundance FC (0 to 43.5 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP2C19 FC, <= 43.5'
    cyp2d6_c2 <- fixed(-1.17e-03)
    label("CYP2D6 abundance FC (0 to 51.4 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP2D6 FC, <= 51.4': -1.17E-03*GesAge^2 + 5.99E-02*GesAge + 1
    cyp2d6_c1 <- fixed(5.99e-02)
    label("CYP2D6 abundance FC (0 to 51.4 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP2D6 FC, <= 51.4'
    cyp3a4_preg_c2 <- fixed(-1.11e-03)
    label("CYP3A4 abundance FC, pregnancy, quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'CYP3A4 FC, Pregnant, <= 40': -1.11E-03*GesAge^2 + 6.23E-02*GesAge + 1
    cyp3a4_preg_c1 <- fixed(6.23e-02)
    label("CYP3A4 abundance FC, pregnancy, linear GA coefficient (1/week)") # Suppl. Table 2 'CYP3A4 FC, Pregnant, <= 40'
    cyp3a4_pp_c1 <- fixed(-0.143)
    label("CYP3A4 abundance FC, postpartum (40 to 45 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'CYP3A4 FC, Postpartum, > 40': -0.143*GesAge + 7.44
    cyp3a4_pp_c0 <- fixed(7.44)
    label("CYP3A4 abundance FC, postpartum (40 to 45 weeks), intercept (unitless)") # Suppl. Table 2 'CYP3A4 FC, Postpartum, > 40'
    ugt1a1_c2 <- fixed(-4.14e-04)
    label("UGT1A1 abundance FC (0 to 46.5 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'UGT1A1 FC, <= 46.5': -4.14E-04*GesAge^2 + 1.93E-02*GesAge + 1
    ugt1a1_c1 <- fixed(1.93e-02)
    label("UGT1A1 abundance FC (0 to 46.5 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'UGT1A1 FC, <= 46.5'
    ugt1a4_c2 <- fixed(-2.10e-03)
    label("UGT1A4 abundance FC (0 to 52.1 weeks), quadratic GA coefficient (1/week^2)") # Suppl. Table 2 'UGT1A4 FC, <= 52.1': -2.10E-03*GesAge^2 + 1.09E-01*GesAge + 1
    ugt1a4_c1 <- fixed(1.09e-01)
    label("UGT1A4 abundance FC (0 to 52.1 weeks), linear GA coefficient (1/week)") # Suppl. Table 2 'UGT1A4 FC, <= 52.1'

    # ---------------------------------------------------------------
    # Virtual-population draws: the CV column of Supplementary Table 2,
    # attached to the non-pregnant regressions and to the fetal,
    # placental and amniotic-fluid weights. Each quantity is drawn as
    # Mean * (1 + eta) with eta ~ N(0, (CV/100)^2), the normrnd(Mean,
    # (CV/100)*Mean) form of the framework's population generator (Stader
    # 2019), so each variance is (CV/100)^2. Lung weight has CV 0.
    # ---------------------------------------------------------------
    etaht ~ fixed(0.001444) # Body height CV 3.8 (Suppl. Table 2)
    etawt ~ fixed(0.023104) # Body weight (non-pregnant) CV 15.2 (Suppl. Table 2)
    etaw_adipose ~ fixed(0.087616) # Adipose weight CV 29.6 (Suppl. Table 2)
    etaw_bone ~ fixed(0.017424) # Bone weight CV 13.2 (Suppl. Table 2)
    etaw_brain ~ fixed(0.0081) # Brain weight CV 9.0 (Suppl. Table 2)
    etaw_gonads ~ fixed(0.121104) # Gonad weight CV 34.8 (Suppl. Table 2)
    etaw_heart ~ fixed(0.041209) # Heart weight CV 20.3 (Suppl. Table 2)
    etaw_kidney ~ fixed(0.053824) # Kidney weight CV 23.2 (Suppl. Table 2)
    etaw_uterus ~ fixed(0.021904) # Uterus weight CV 14.8 (Suppl. Table 2)
    etaw_mammary ~ fixed(0.26677225) # Mammary gland CV 51.65 (Suppl. Table 2)
    etaw_muscle ~ fixed(0.013924) # Muscle weight CV 11.8 (Suppl. Table 2)
    etaw_skin ~ fixed(0.006889) # Skin weight CV 8.3 (Suppl. Table 2)
    etaw_thymus ~ fixed(0.200704) # Thymus weight CV 44.8 (Suppl. Table 2)
    etaw_gut ~ fixed(0.005329) # Gut weight CV 7.3 (Suppl. Table 2)
    etaw_spleen ~ fixed(0.267289) # Spleen weight CV 51.7 (Suppl. Table 2)
    etaw_pancreas ~ fixed(0.077284) # Pancreas weight CV 27.8 (Suppl. Table 2)
    etaw_liver ~ fixed(0.056169) # Liver weight CV 23.7 (Suppl. Table 2)
    etaw_blood ~ fixed(0.010816) # Blood weight CV 10.4 (Suppl. Table 2)
    etaw_fetus ~ fixed(0.07568001) # Foetal weight CV 27.51 (Suppl. Table 2)
    etaw_placenta ~ fixed(0.032041) # Placental weight CV 17.9 (Suppl. Table 2)
    etaw_amniotic ~ fixed(0.08952064) # Amniotic fluid weight CV 29.92 (Suppl. Table 2)
    etaco ~ fixed(0.044521) # Cardiac output CV 21.1 (Suppl. Table 2)
    etagfr ~ fixed(0.021609) # GFR CV 14.7 (Suppl. Table 2)
    etahct ~ fixed(0.020736) # HCT CV 14.4 (Suppl. Table 2)
    etahsa ~ fixed(0.006241) # Albumin CV 7.9 (Suppl. Table 2)
    etaagp ~ fixed(0.059049) # AAG CV 24.3 (Suppl. Table 2)
  })

  model({
    # The rxode2 time variable IS gestational age in weeks; every woman
    # delivers at exactly 40 weeks (Discussion), so GA > 40 is postpartum
    # and GA - 40 is the week after delivery. The equations were fitted
    # and the population generated over 0 to 66 weeks (Methods).
    ga <- time
    preg <- 1 * (ga <= 40)
    pp <- 1 - preg

    # All women: the 'Sex' term of the framework regressions is 1 for
    # female (Stader 2019 generator convention, 1 = female, 0 = male).
    sexf <- 1

    # ==================================================================
    # 1. Anthropometry (Suppl. Table 2, 'Anthropometric parameters')
    # ==================================================================
    ht <- (-0.0039 * AGE^2 + 0.238 * AGE - 12.5 * sexf + 176) * (1 + etaht) # Body height (cm)
    wt_np <- (-0.0039 * AGE^2 + 1.12 * ht + 0.611 * AGE - 0.424 * sexf - 137) * (1 + etawt) # Body weight, non-pregnant (kg)
    fc_bw <- preg * (bw_preg_c2 * ga^2 + bw_preg_c1 * ga + 1) +
      pp * (bw_pp_c2 * ga^2 + bw_pp_c1 * ga + bw_pp_c0)
    wt <- wt_np * fc_bw # Body weight at GA (kg)
    bmi <- wt / (ht / 100)^2 # BMI (kg/m^2)
    bsa <- 0.007184 * ht^0.725 * wt^0.425 # BSA at GA (m^2), DuBois & DuBois
    # The non-pregnant regressions below are 'Non pregnant' rows, so they
    # take the non-pregnant body surface area.
    bsa_np <- 0.007184 * ht^0.725 * wt_np^0.425

    # ==================================================================
    # 2. Organ weights, kg (Suppl. Table 2, 'Organ weights')
    # ==================================================================
    w_lung <- exp(0.028 * ht + 0.0077 * AGE - 5.6) # Lung weight; CV 0; FC 1

    fc_adipose <- preg * (adipose_preg_c1 * ga + 1) +
      pp * (adipose_pp_c1 * ga + adipose_pp_c0)
    w_adipose <- (0.68 * wt_np - 0.56 * ht + 6.1 * sexf + 65) * (1 + etaw_adipose) * fc_adipose

    w_bone <- exp(0.024 * ht - 1.9) * (1 + etaw_bone) # FC 1

    fc_brain <- preg * (brain_preg_c1 * ga + 1) +
      pp * ((ga <= 42) * (brain_pp_c1 * ga + brain_pp_c0) + (ga > 42))
    w_brain <- exp(-0.0075 * AGE + 0.0078 * ht - 0.97) * (1 + etaw_brain) * fc_brain

    w_gonads <- (-0.00034 * wt_np - 0.00022 * AGE - 0.03 * sexf + 0.072) * (1 + etaw_gonads) # FC 1

    fc_heart <- preg * (heart_preg_c1 * ga + 1) +
      pp * ((ga <= 86) * (heart_pp_c1 * ga + heart_pp_c0) + (ga > 86))
    w_heart <- (0.34 * bsa_np + 0.0018 * AGE - 0.36) * (1 + etaw_heart) * fc_heart

    fc_kidney <- preg * (kidney_preg_c1 * ga + 1) +
      pp * ((ga <= 50) * (kidney_pp_c1 * ga + kidney_pp_c0) + (ga > 50))
    w_kidney <- (-0.00038 * AGE - 0.056 * sexf + 0.33) * (1 + etaw_kidney) * fc_kidney

    # Uterus and mammary gland (Abduljalil 2012 equations) sit inside the
    # muscle compartment (Results, after Gaohua 2012).
    w_uterus_np <- 0.079 * (1 + etaw_uterus)
    fc_uterus <- preg * (2.74e-01 * ga + 1) +
      pp * ((ga <= 49.7) * (71822 * exp(-0.225 * ga)) + (ga > 49.7))
    w_uterus <- w_uterus_np * fc_uterus
    w_mammary_np <- 0.720 * (1 + etaw_mammary)
    fc_mammary <- preg * (3.67e-02 * ga + 1) +
      pp * ((ga <= 41.2) * (-1.29 * ga + 54) + (ga > 41.2))
    w_mammary <- w_mammary_np * fc_mammary
    w_muscle_np <- (17.9 * bsa_np - 0.0667 * AGE - 5.68 * sexf - 1.22) * (1 + etaw_muscle)
    # 'Total muscle weight' row; the muscle-only FC is 1.
    w_muscle <- (w_muscle_np - w_uterus_np - w_mammary_np) + w_uterus + w_mammary

    w_skin <- exp(-0.0058 * AGE - 0.37 * sexf + 1.13) * (1 + etaw_skin) # FC 1
    w_thymus <- 0.0221 * (1 + etaw_thymus) # FC 1
    w_gut <- 3e-06 * ht^2.49 * (1 + etaw_gut) # FC 1
    w_spleen <- exp(1.13 * bsa_np - 3.93) * (1 + etaw_spleen) # FC 1
    w_pancreas <- 0.103 * (1 + etaw_pancreas) # FC 1
    w_liver <- exp(0.87 * bsa_np - 0.0014 * AGE - 1.0) * (1 + etaw_liver) # FC 1

    # Blood: non-pregnant blood weight split into plasma and red cells by
    # the non-pregnant haematocrit, each scaled by its own FC. The printed
    # 'Blood weight, Pregnant and postpartum' row also divides the plasma
    # term by (1 - 0.91*HCT), which counts the red cells twice (+105% at
    # 38 weeks); the sum of the two components gives the +47% the Results
    # report for blood at the third trimester and is used here.
    hct_np <- (0.443 - 0.033 * sexf) * (1 + etahct)
    w_blood_np <- exp(0.067 * bsa_np - 0.0025 * AGE - 0.38 * sexf + 1.7) * (1 + etaw_blood)
    w_plasma_np <- w_blood_np * (1 - 0.91 * hct_np)
    w_rbc_np <- w_blood_np - w_plasma_np
    fc_plasma <- preg * (plasma_preg_c1 * ga + 1) +
      pp * ((ga <= 57.7) * (plasma_pp_c1 * ga + plasma_pp_c0) + (ga > 57.7))
    fc_rbc <- preg * (rbc_preg_c2 * ga^2 + rbc_preg_c1 * ga + 1) +
      pp * ((ga <= 46.5) * (rbc_pp_c2 * ga^2 + rbc_pp_c1 * ga + rbc_pp_c0) + (ga > 46.5))
    w_plasma <- w_plasma_np * fc_plasma
    w_rbc <- w_rbc_np * fc_rbc
    w_blood <- w_plasma + w_rbc

    # Fetoplacental unit (Abduljalil 2012 equations as printed in Suppl.
    # Table 2): volumes in mL / 1000 times tissue density; zero after
    # delivery.
    w_fetus <- preg * ((0.01 * exp((0.945 / 0.0702) * (1 - exp(-0.0702 * ga)))) / 1000) * 1.03 * (1 + etaw_fetus)
    w_placenta <- preg * ((-0.0122 * ga^3 + 0.9149 * ga^2 - 0.716 * ga) / 1000) * 1.048 * (1 + etaw_placenta)
    w_amniotic <- preg * ((0.00005 * ga^5 - 0.0061 * ga^4 + 0.2064 * ga^3 - 1.2056 * ga^2 + 1.9648 * ga) / 1000) * (1 + etaw_amniotic)
    w_fetoplacental <- w_fetus + w_placenta + w_amniotic

    # ==================================================================
    # 3. Cardiac output (L/h) and regional blood flows (Suppl. Table 2,
    # 'Blood flow parameters'). Regional flows are percent of cardiac
    # output; absolute flows q_* are fq_* / 100 * co.
    # ==================================================================
    fc_co <- (ga <= 50.9) * (co_c2 * ga^2 + co_c1 * ga + 1) + (ga > 50.9)
    co_np <- (159 * bsa_np - 1.56 * AGE + 114) * (1 + etaco)
    co <- co_np * fc_co

    # Organs whose absolute flow the authors held unchanged share one FC
    # on their percent of cardiac output. As printed it jumps at delivery
    # (0.604 to 0.887 at 40 weeks).
    fc_qfrac <- preg * (qfrac_preg_c1 * ga + 1) +
      pp * ((ga <= 55) * (qfrac_pp_c1 * ga + qfrac_pp_c0) + (ga > 55))
    fq_adipose <- ((0.044 + 0.027 * sexf) * AGE + 2.4 * sexf + 3.9) * fc_qfrac
    fq_bone <- 5 * fc_qfrac
    fq_gonads <- (-0.03 * sexf + 0.05) * fc_qfrac
    fq_heart <- (-0.72 * ht - 10 * sexf + 134) * fc_qfrac
    fq_skin <- 5 * fc_qfrac
    fq_thymus <- 1.5 * fc_qfrac
    fq_liver <- (-0.108 * AGE + 1.04 * sexf + 27.9) * fc_qfrac
    fq_gut <- (2 * sexf + 14) * fc_qfrac
    fq_spleen <- 3 * fc_qfrac
    fq_pancreas <- 1 * fc_qfrac

    # Brain and kidney: the regressions give the change in ABSOLUTE flow,
    # so the percent is divided by the cardiac-output FC.
    fc_qbrain <- (ga <= 82.4) * (qbrain_c2 * ga^2 + qbrain_c1 * ga + 1) / fc_co + (ga > 82.4)
    fq_brain <- exp(-0.48 * bsa_np + 0.04 * sexf + 3.5) * fc_qbrain
    fc_qkidney <- (ga <= 40.7) * (qkidney_c2 * ga^2 + qkidney_c1 * ga + 1) / fc_co + (ga > 40.7)
    fq_kidney <- (-8.7 * bsa_np + 0.29 * ht - 0.081 * AGE - 13) * fc_qkidney

    # Uterus flow is part of the muscle flow.
    fq_uterus_np <- 0.947
    fc_quterus <- preg * (quterus_preg_c2 * ga^2 + quterus_preg_c1 * ga + 1) +
      pp * ((ga <= 53) * (quterus_pp_c1 * ga + quterus_pp_c0) + (ga > 53))
    fq_uterus <- fq_uterus_np * fc_quterus
    fq_muscle <- ((-6.4 * sexf + 17.5) - fq_uterus_np) * fc_qfrac + fq_uterus

    # Fetoplacental flow (Gaohua 2012 form, refitted coefficients as
    # printed), listed under the percent-of-cardiac-output flows; zero
    # before 2.23 weeks, where the cubic is negative, and after delivery.
    # The row's range is '> 2.23 and < 40'; it is closed at 40 like every
    # other pregnancy row so that the flow and the fetoplacental weight
    # end together at delivery.
    fq_fetoplacental <- (ga > 2.23) * preg * (-4.00e-04 * ga^3 + 2.76e-02 * ga^2 - 5.95e-02 * ga)

    q_adipose <- fq_adipose / 100 * co
    q_bone <- fq_bone / 100 * co
    q_gonads <- fq_gonads / 100 * co
    q_heart <- fq_heart / 100 * co
    q_skin <- fq_skin / 100 * co
    q_thymus <- fq_thymus / 100 * co
    q_liver <- fq_liver / 100 * co
    q_gut <- fq_gut / 100 * co
    q_spleen <- fq_spleen / 100 * co
    q_pancreas <- fq_pancreas / 100 * co
    q_brain <- fq_brain / 100 * co
    q_kidney <- fq_kidney / 100 * co
    q_uterus <- fq_uterus / 100 * co
    q_muscle <- fq_muscle / 100 * co
    q_fetoplacental <- fq_fetoplacental / 100 * co

    # ==================================================================
    # 4. GFR (mL/min), haematocrit (fraction), albumin and AAG (g/L)
    # ==================================================================
    fc_gfr <- (ga <= 67.5) * (gfr_c2 * ga^2 + gfr_c1 * ga + 1) + (ga > 67.5)
    gfr <- exp(-0.0079 * AGE + 0.5 * bsa_np + 4.2) * (1 + etagfr) * fc_gfr

    fc_hct <- (ga <= 65) * (hct_c2 * ga^2 + hct_c1 * ga + 1) + (ga > 65)
    hct <- hct_np * fc_hct

    fc_hsa <- preg * (hsa_preg_c2 * ga^2 + hsa_preg_c1 * ga + 1) +
      pp * ((ga <= 45) * (hsa_pp_c1 * ga + hsa_pp_c0) + (ga > 45))
    hsa <- (-0.0709 * AGE + 47.7) * (1 + etahsa) * fc_hsa

    fc_agp <- (ga <= 63.7) * (agp_c2 * ga^2 + agp_c1 * ga + 1) + (ga > 63.7)
    agp <- 0.798 * (1 + etaagp) * fc_agp

    # ==================================================================
    # 5. Liver enzyme-abundance FCs, applied unchanged to the intestinal
    # enzymes (last row of Suppl. Table 2). Absolute abundances are the
    # framework's and are not printed in the paper.
    # ==================================================================
    fc_cyp1a2 <- (ga <= 51.8) * (cyp1a2_c2 * ga^2 + cyp1a2_c1 * ga + 1) + (ga > 51.8)
    fc_cyp2b6 <- (ga <= 50.8) * (cyp2b6_c2 * ga^2 + cyp2b6_c1 * ga + 1) + (ga > 50.8)
    fc_cyp2c9 <- (ga <= 42) * (cyp2c9_c2 * ga^2 + cyp2c9_c1 * ga + 1) + (ga > 42)
    fc_cyp2c19 <- (ga <= 43.5) * (cyp2c19_c2 * ga^2 + cyp2c19_c1 * ga + 1) + (ga > 43.5)
    fc_cyp2d6 <- (ga <= 51.4) * (cyp2d6_c2 * ga^2 + cyp2d6_c1 * ga + 1) + (ga > 51.4)
    fc_cyp3a4 <- preg * (cyp3a4_preg_c2 * ga^2 + cyp3a4_preg_c1 * ga + 1) +
      pp * ((ga <= 45) * (cyp3a4_pp_c1 * ga + cyp3a4_pp_c0) + (ga > 45))
    # 'Assumed similar fold change as for the CYP3A4' up to 53.4 weeks,
    # 1 after; CYP3A4 is already back to 1 from 45 weeks.
    fc_cyp3a5 <- (ga <= 53.4) * fc_cyp3a4 + (ga > 53.4)
    fc_ugt1a1 <- (ga <= 46.5) * (ugt1a1_c2 * ga^2 + ugt1a1_c1 * ga + 1) + (ga > 46.5)
    fc_ugt1a4 <- (ga <= 52.1) * (ugt1a4_c2 * ga^2 + ugt1a4_c1 * ga + 1) + (ga > 52.1)
  })
}
