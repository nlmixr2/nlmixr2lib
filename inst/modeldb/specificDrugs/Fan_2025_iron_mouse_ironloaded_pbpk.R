Fan_2025_iron_mouse_ironloaded_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 14-compartment, perfusion-limited; NONMEM ADVAN15).",
    "Preclinical (mouse, iron-loaded diet).",
    "Mechanistic whole-body PBPK model of iron disposition after a single",
    "intravenous bolus of 0.2 umol/kg 59Fe. Thirteen measured compartments",
    "(plasma, bone, gut, heart, kidney, liver, muscle, skin, lung, spleen,",
    "red-cell haemoglobin iron, brain, fat) plus a lumped remainder were",
    "fitted simultaneously by naive pooled analysis of digitised mean data.",
    "Every organ is perfusion-limited with a tissue-to-plasma partition",
    "coefficient kp; the mechanistic limbs are erythropoietic iron",
    "utilisation from bone into the red-cell pool (q_bone_rbc), release of",
    "that iron back into the spleen by macrophage-mediated",
    "erythrophagocytosis over the red-cell lifespan (mtt_rbc), portal",
    "drainage of gut and spleen into liver, and an unavoidable physiologic",
    "iron loss clearance from plasma (cl). Parameters are the",
    "iron-loaded diet column of Table 2; organ volumes and blood flows are",
    "fixed to mouse physiology. No between-subject variability was",
    "estimated (all seven $OMEGA elements are 0 FIX), so the model carries",
    "no etas. States hold AMOUNTS (nmol), matching the published $ERROR",
    "block, and the red-cell pool has no volume."
  )
  reference <- paste(
    "Fan X, Cao K, Wong RSM, Yan X. A whole-body mechanistic",
    "physiologically-based pharmacokinetic modeling of intravenous iron.",
    "Drug Deliv Transl Res. 2025;15(3):1109-1120.",
    "doi:10.1007/s13346-024-01675-x (PMCID: PMC11870943).",
    "The complete NONMEM control stream ($MODEL / $PK / $DES / $ERROR /",
    "$THETA / $OMEGA / $SIGMA), the rat and human physiology tables and",
    "every digitised observed dataset are in the Electronic Supplementary",
    "Material (MOESM1)."
  )
  vignette <- "Fan_2025_iron_wholebody_pbpk"

  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  dosing <- "plasma"

  compartmentData <- list(
    plasma   = list(analyte = "iron-59", units = "nmol", specimen = "plasma", verified = TRUE),
    bone     = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    gut      = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    heart    = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    skin     = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    lung     = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    spleen   = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    rbc_iron = list(analyte = "iron-59", units = "nmol", specimen = "blood cell", verified = TRUE),
    brain    = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    adipose  = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = TRUE),
    other    = list(analyte = "iron-59", units = "nmol", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species       = "mouse (iron-loaded diet)",
    n_subjects    = "not reported; digitised group mean data only",
    n_studies     = 1,
    weight_range  = "25 g (single body weight used for all physiology)",
    disease_state = "Iron-loaded dietary iron status induced by diet",
    dose_range    = "single 0.2 umol/kg 59Fe intravenous bolus (5 nmol in a 25 g mouse)",
    notes         = paste(
      "Source PK / biodistribution study is reference [23] of the paper;",
      "mice were killed and dissected at 15 min, 12 h, 24 h and days 4, 7,",
      "14 and 28. Only group means were available, so a naive pooled data",
      "method was used, treating all observations as one individual.",
      "Fig. 3's caption states 0.02 umol/kg but Methods states 0.2 umol/kg;",
      "the supplement's data tables start plasma at 5 nmol, which confirms",
      "0.2 umol/kg x 0.025 kg = 5 nmol. See the vignette Errata."
    )
  )

  ini({
    # ---- Estimated iron-disposition parameters: Table 2, iron-loaded diet column ----
    # RSE% from Table 2 is given after each value.
    lcl <- log(1.184e-4); label("Unavoidable physiologic iron loss clearance from plasma (L/h)")  # Table 2 CLloss = 1.184 x 10^-4 L/h (RSE 0.035%)
    lkp_adipose <- log(0.328); label("Fat (adipose)-to-plasma partition coefficient (unitless)")  # Table 2 KPfat (RSE 0.025%)
    lkp_bone    <- log(7.753); label("Bone-to-plasma partition coefficient (unitless)")  # Table 2 KPbon (RSE 0.038%)
    lkp_brain   <- log(0.9304); label("Brain-to-plasma partition coefficient (unitless)")  # Table 2 KPbra (RSE 0.048%)
    lkp_gut     <- log(6.979); label("Gut-to-plasma partition coefficient (unitless)")  # Table 2 KPgut (RSE 0.033%)
    lkp_heart   <- log(25.88); label("Heart-to-plasma partition coefficient (unitless)")  # Table 2 KPhea (RSE 0.043%)
    lkp_kidney  <- log(29.06); label("Kidney-to-plasma partition coefficient (unitless)")  # Table 2 KPkid (RSE 0.042%)
    lkp_liver   <- log(38.53); label("Liver-to-plasma partition coefficient (unitless)")  # Table 2 KPliv (RSE 0.036%)
    lkp_lung    <- log(45.24); label("Lung-to-plasma partition coefficient (unitless)")  # Table 2 KPlun (RSE 0.027%)
    lkp_muscle  <- log(2.394); label("Muscle-to-plasma partition coefficient (unitless)")  # Table 2 KPmus (RSE 0.034%)
    lkp_other   <- log(3.06e-8); label("Remainder (rest of body)-to-plasma partition coefficient (unitless)")  # Table 2 KPrem (RSE 8.16%)
    lkp_skin    <- log(4.252); label("Skin-to-plasma partition coefficient (unitless)")  # Table 2 KPski (RSE 0.027%)
    lkp_spleen  <- log(23.93); label("Spleen-to-plasma partition coefficient (unitless)")  # Table 2 KPspl (RSE 0.046%)
    lq_bone_rbc <- log(0.0388e-3); label("Erythropoietic iron utilisation flow, bone to red cells (L/h)")  # Table 2 QE (RSE 0.089%)
    lmtt_rbc <- log(197.6); label("Mean red-blood-cell lifespan (h)")  # Table 2 TRBC (RSE 0.016%)

    # ---- Mouse physiology, fixed. Taken from the supplement's NONMEM $PK
    # block, which reproduces every entry of Table 1 to 3 significant
    # figures. Volumes keep the source arithmetic visible: organ weights are
    # reported for a 30 g mouse and rescaled to the 25 g animal studied.
    bw       <- fixed(0.025); label("Body weight (kg)")  # $PK: BW = 0.025
    co_allom <- fixed(0.275); label("Cardiac output coefficient (L/min per kg^0.75)")  # $PK: CO = 0.275*BW**0.75, Brown 1988
    v_body   <- fixed(0.0326); label("Total body volume (L)")  # $PK: VBODY = 0.0326

    # Organ blood flows as fractions of cardiac output ($PK block)
    fq_adipose <- fixed(0.07);   label("Fat blood flow (fraction of cardiac output)")  # $PK QADI
    fq_bone    <- fixed(0.0407); label("Bone blood flow (fraction of cardiac output)")  # $PK QBON
    fq_brain   <- fixed(0.033);  label("Brain blood flow (fraction of cardiac output)")  # $PK QBRA
    fq_gut     <- fixed(0.1408); label("Gut blood flow (fraction of cardiac output)")  # $PK QGUT
    fq_heart   <- fixed(0.066);  label("Heart blood flow (fraction of cardiac output)")  # $PK QHEA
    fq_kidney  <- fixed(0.091);  label("Kidney blood flow (fraction of cardiac output)")  # $PK QKID
    fq_liver   <- fixed(0.161);  label("Total liver blood flow (fraction of cardiac output)")  # $PK QLIV
    fq_muscle  <- fixed(0.159);  label("Muscle blood flow (fraction of cardiac output)")  # $PK QMUS
    fq_skin    <- fixed(0.058);  label("Skin blood flow (fraction of cardiac output)")  # $PK QSKI
    fq_spleen  <- fixed(0.0112); label("Spleen blood flow (fraction of cardiac output)")  # $PK QSPL

    # Absolute organ volumes (L) ($PK block; reproduce Table 1 column 2 in mL)
    v_adipose <- fixed(2.59 * 25 / (30 * 1000));  label("Fat volume (L)")  # $PK VADI; Table 1: 2.16 mL
    v_bone    <- fixed(0.1073 * 0.0326);          label("Bone volume (L)")  # $PK VBON; Table 1: 3.50 mL
    v_brain   <- fixed(0.50 * 25 / (30 * 1000));  label("Brain volume (L)")  # $PK VBRA; Table 1: 0.42 mL
    v_gut     <- fixed(1.27 * 25 / (30 * 1000));  label("Gut volume (L)")  # $PK VGUT; Table 1: 1.06 mL
    v_heart   <- fixed(0.15 * 25 / (30 * 1000));  label("Heart volume (L)")  # $PK VHEA; Table 1: 0.125 mL
    v_kidney  <- fixed(0.5 * 25 / (30 * 1000));   label("Kidney volume (L)")  # $PK VKID; Table 1: 0.42 mL
    v_liver   <- fixed(1.65 * 25 / (30 * 1000));  label("Liver volume (L)")  # $PK VLIV; Table 1: 1.375 mL
    v_lung    <- fixed(0.22 * 25 / (30 * 1000));  label("Lung volume (L)")  # $PK VLUN; Table 1: 0.18 mL
    v_muscle  <- fixed(11.5 * 25 / (30 * 1000));  label("Muscle volume (L)")  # $PK VMUS; Table 1: 9.58 mL
    v_skin    <- fixed(4.95 * 25 / (30 * 1000));  label("Skin volume (L)")  # $PK VSKI; Table 1: 4.125 mL
    v_spleen  <- fixed(0.11 * 25 / (30 * 1000));  label("Spleen volume (L)")  # $PK VSPL; Table 1: 0.092 mL
    v_plasma  <- fixed(1.2 * 25 / (27.4 * 1000)); label("Plasma volume (L)")  # $PK VBLD; Table 1: 1.09 mL

    # ---- Residual error. Structure is fully determined by the $ERROR block
    # (per-compartment proportional or additive); the VALUES below are the
    # $SIGMA entries of the control stream, which are INITIAL estimates --
    # the paper reports no final residual estimates anywhere. Encoded as SDs
    # (nlmixr2) from the NONMEM variances. See the vignette Errata.
    propSd_plasma   <- sqrt(0.1); label("Proportional residual error, plasma (fraction)")  # $SIGMA(1) initial = 0.1 (variance)
    propSd_bone     <- sqrt(0.5); label("Proportional residual error, bone (fraction)")  # $SIGMA(2) initial = 0.5 (variance)
    propSd_gut      <- sqrt(0.1); label("Proportional residual error, gut (fraction)")  # $SIGMA(3) initial = 0.1
    propSd_skin     <- sqrt(0.1); label("Proportional residual error, skin (fraction)")  # $SIGMA(8) initial = 0.1
    propSd_lung     <- sqrt(0.1); label("Proportional residual error, lung (fraction)")  # $SIGMA(9) initial = 0.1
    propSd_spleen   <- sqrt(0.1); label("Proportional residual error, spleen (fraction)")  # $SIGMA(10) initial = 0.1
    propSd_brain    <- sqrt(0.1); label("Proportional residual error, brain (fraction)")  # $SIGMA(12) initial = 0.1
    addSd_heart     <- sqrt(0.1); label("Additive residual error, heart (nmol)")  # $SIGMA(17) initial = 0.1; $SIGMA(4) prop arm is 0 FIX
    addSd_kidney    <- sqrt(0.1); label("Additive residual error, kidney (nmol)")  # $SIGMA(18) initial = 0.1; $SIGMA(5) prop arm is 0 FIX
    addSd_liver     <- sqrt(0.1); label("Additive residual error, liver (nmol)")  # $SIGMA(19) initial = 0.1; $SIGMA(6) prop arm is 0 FIX
    addSd_muscle    <- sqrt(0.1); label("Additive residual error, muscle (nmol)")  # $SIGMA(20) initial = 0.1; $SIGMA(7) prop arm is 0 FIX
    addSd_rbc_iron  <- sqrt(0.1); label("Additive residual error, red-cell iron (nmol)")  # $SIGMA(24) initial = 0.1; $SIGMA(11) prop arm is 0 FIX
    addSd_adipose   <- sqrt(0.1); label("Additive residual error, fat (nmol)")  # $SIGMA(26) initial = 0.1; $SIGMA(13) prop arm is 0 FIX
  })

  model({
    # ---- 1. Estimated iron-disposition parameters back-transformed ----
    cl <- exp(lcl)
    kp_adipose <- exp(lkp_adipose)
    kp_bone <- exp(lkp_bone)
    kp_brain <- exp(lkp_brain)
    kp_gut <- exp(lkp_gut)
    kp_heart <- exp(lkp_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_lung <- exp(lkp_lung)
    kp_muscle <- exp(lkp_muscle)
    kp_other <- exp(lkp_other)
    kp_skin <- exp(lkp_skin)
    kp_spleen <- exp(lkp_spleen)
    q_bone_rbc <- exp(lq_bone_rbc)
    mtt_rbc <- exp(lmtt_rbc)

    # ---- 2. System physiology: cardiac output, organ blood flows, volumes ----
    qc <- 60 * co_allom * bw^0.75
    q_adipose <- fq_adipose * qc
    q_bone <- fq_bone * qc
    q_brain <- fq_brain * qc
    q_gut <- fq_gut * qc
    q_heart <- fq_heart * qc
    q_kidney <- fq_kidney * qc
    q_liver <- fq_liver * qc
    q_muscle <- fq_muscle * qc
    q_skin <- fq_skin * qc
    q_spleen <- fq_spleen * qc
    q_lung <- qc
    # Rest-of-body flow, computed as in the fitted $PK block. NOTE: gut and
    # spleen are NOT subtracted again here (they are inside q_liver), which
    # is what makes the plasma efflux sum to qc exactly. Table 1's printed
    # Remainder flow of 0.176 L/h double-subtracts them; see vignette Errata.
    q_other <- qc - q_bone - q_heart - q_kidney - q_liver - q_muscle -
      q_skin - q_adipose - q_brain
    v_other <- v_body - v_bone - v_gut - v_heart - v_kidney - v_liver -
      v_muscle - v_skin - v_lung - v_spleen - v_plasma - v_adipose - v_brain

    # ---- 3. Compartment concentrations (every state holds an AMOUNT) ----
    Cp <- plasma / v_plasma
    Cadipose <- adipose / v_adipose
    Cbone <- bone / v_bone
    Cbrain <- brain / v_brain
    Cgut <- gut / v_gut
    Cheart <- heart / v_heart
    Ckidney <- kidney / v_kidney
    Cliver <- liver / v_liver
    Clung <- lung / v_lung
    Cmuscle <- muscle / v_muscle
    Cother <- other / v_other
    Cskin <- skin / v_skin
    Cspleen <- spleen / v_spleen

    # ---- 4. Whole-body mass balance ($DES) ----
    # Arterial / plasma compartment (Eq. 6). The efflux coefficient below sums
    # EXACTLY to the cardiac output, because gut and spleen arterial inflow are
    # already inside q_liver (the liver receives only the hepatic-artery share).
    d/dt(plasma) <- -(q_other + q_bone + q_heart + q_kidney + q_liver + q_muscle +
                      q_skin + cl + q_brain + q_adipose) * Cp +
      (q_bone / kp_bone) * Cbone + (q_heart / kp_heart) * Cheart +
      (q_kidney / kp_kidney) * Ckidney + (q_liver / kp_liver) * Cliver +
      (q_muscle / kp_muscle) * Cmuscle + (q_skin / kp_skin) * Cskin +
      (q_other / kp_other) * Cother + (q_brain / kp_brain) * Cbrain +
      (q_adipose / kp_adipose) * Cadipose
    # Bone: perfusion-limited uptake minus erythropoietic iron utilisation (Eq. 2)
    d/dt(bone) <- q_bone * Cp - (q_bone / kp_bone + q_bone_rbc) * Cbone
    # Perfusion-limited organs (Eq. 1)
    d/dt(gut)     <- q_gut * Cp - (q_gut / kp_gut) * Cgut
    d/dt(heart)   <- q_heart * Cp - (q_heart / kp_heart) * Cheart
    d/dt(kidney)  <- q_kidney * Cp - (q_kidney / kp_kidney) * Ckidney
    d/dt(muscle)  <- q_muscle * Cp - (q_muscle / kp_muscle) * Cmuscle
    d/dt(skin)    <- q_skin * Cp - (q_skin / kp_skin) * Cskin
    d/dt(brain)   <- q_brain * Cp - (q_brain / kp_brain) * Cbrain
    d/dt(adipose) <- q_adipose * Cp - (q_adipose / kp_adipose) * Cadipose
    d/dt(other)   <- q_other * Cp - (q_other / kp_other) * Cother
    # Lung. Reproduced verbatim from the published $DES: the lung draws the full
    # cardiac output out of plasma and returns nothing, yet q_lung appears in
    # NEITHER the plasma efflux sum NOR the plasma influx sum. The lung is
    # therefore a non-mass-balanced observer compartment in the published model.
    # See the vignette Errata; kept as fitted rather than silently corrected.
    d/dt(lung) <- q_lung * Cp - (q_lung / kp_lung) * Clung
    # Liver: hepatic-artery inflow plus portal return from gut and spleen (Eq. 5)
    d/dt(liver) <- (q_liver - q_gut - q_spleen) * Cp + (q_gut / kp_gut) * Cgut +
      (q_spleen / kp_spleen) * Cspleen - (q_liver / kp_liver) * Cliver
    # Spleen: erythrophagocytic recycling of senescent RBC iron (Eq. 4)
    d/dt(spleen) <- q_spleen * Cp + rbc_iron / mtt_rbc -
      (q_spleen / kp_spleen) * Cspleen
    # Red-cell haemoglobin iron pool; amount-only state, no volume (Eq. 3)
    d/dt(rbc_iron) <- q_bone_rbc * Cbone - rbc_iron / mtt_rbc

    # ---- 5. Observations ----
    # The published $ERROR block observes AMOUNTS: IPRED = A(n) for every
    # compartment, in nmol. Concentrations are exposed alongside for
    # convenience; Cc is the plasma concentration.
    Cc <- Cp
    plasma   ~ prop(propSd_plasma)
    bone     ~ prop(propSd_bone)
    gut      ~ prop(propSd_gut)
    skin     ~ prop(propSd_skin)
    lung     ~ prop(propSd_lung)
    spleen   ~ prop(propSd_spleen)
    brain    ~ prop(propSd_brain)
    heart    ~ add(addSd_heart)
    kidney   ~ add(addSd_kidney)
    liver    ~ add(addSd_liver)
    muscle   ~ add(addSd_muscle)
    rbc_iron ~ add(addSd_rbc_iron)
    adipose  ~ add(addSd_adipose)
  })
}
