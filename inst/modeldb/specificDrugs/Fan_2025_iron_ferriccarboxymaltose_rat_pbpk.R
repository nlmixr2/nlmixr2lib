Fan_2025_iron_ferriccarboxymaltose_rat_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 15-compartment, perfusion-limited; NONMEM ADVAN15).",
    "Rat (Sprague-Dawley, iron-deficiency anaemia).",
    "Interspecies extrapolation of the Fan 2025 mouse whole-body iron PBPK",
    "model to rat, for the iron-carbohydrate complex preparation (IVIP)",
    "ferric carboxymaltose (FCM). The 14-compartment mouse structure is",
    "carried over unchanged and rat physiology is fixed to the literature",
    "(Supplementary Table 1); one extra circulating `ivip` compartment receives the",
    "whole intravenous dose and is taken up by splenic macrophages at a rate",
    "set equal to the spleen blood flow (KA := q_spleen), the authors having",
    "removed the direct IVIP-to-plasma release limb as ~0.1% of dose. The",
    "mechanistic limbs are unchanged: erythropoietic iron utilisation from",
    "bone into the red-cell pool, erythrophagocytic return of that iron to",
    "the spleen over the red-cell lifespan, portal drainage of gut and",
    "spleen into liver, and a physiologic iron loss clearance from plasma.",
    "Only bone, liver, spleen, heart, muscle and kidney kp values were estimated from the rat data; every other kp, and the loss clearance, were inherited from the mouse fit. The paper does not say WHICH mouse iron-status column it inherited from; the iron-adequate column is used here per the operator ruling (sidecar q1 = B). See the vignette Errata.",
    "The routing of the released iron into the SPLEEN (rather than directly",
    "to plasma) is a structural assumption -- the paper describes the limb",
    "only in prose and the supplied control stream is the mouse, pure-iron",
    "run. No between-subject variability is estimated, so the model carries",
    "no etas."
  )
  reference <- paste(
    "Fan X, Cao K, Wong RSM, Yan X. A whole-body mechanistic",
    "physiologically-based pharmacokinetic modeling of intravenous iron.",
    "Drug Deliv Transl Res. 2025;15(3):1109-1120.",
    "doi:10.1007/s13346-024-01675-x (PMCID: PMC11870943).",
    "The complete NONMEM control stream for the MOUSE run, the rat and human",
    "physiology tables and every digitised observed dataset are in the",
    "Electronic Supplementary Material (MOESM1)."
  )
  vignette <- "Fan_2025_iron_wholebody_pbpk"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- "ivip"

  # The circulating iron-carbohydrate complex preparation. Paper-mechanistic:
  # it is the undissociated FCM colloid in blood, not a canonical PK state.
  paper_specific_compartments <- "ivip"

  compartmentData <- list(
    ivip     = list(analyte = "ferric carboxymaltose", units = "mg", specimen = "plasma", verified = TRUE),
    plasma   = list(analyte = "iron", units = "mg", specimen = "serum", verified = TRUE),
    bone     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    gut      = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    heart    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    skin     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    lung     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    spleen   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = TRUE),
    rbc_iron = list(analyte = "iron", units = "mg", specimen = "blood cell", verified = FALSE),
    brain    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    adipose  = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    other    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species       = "rat (Sprague-Dawley, iron-deficiency anaemia)",
    n_subjects    = "not reported; digitised group mean data only",
    n_studies     = "1 (reference [24] of the paper)",
    weight_range  = "0.345 kg (single body weight used for all physiology)",
    disease_state = "iron-deficiency anaemia",
    dose_range    = "single 30 mg Fe/kg intravenous ferric carboxymaltose (10.35 mg in a 345 g rat)",
    notes         = paste(
      "Only bone, liver, spleen, heart, muscle and kidney kp values were",
      "estimated from the rat data; every other kp, and the loss clearance,",
      "were inherited from the mouse fit. The paper does not say WHICH mouse",
      "iron-status column it inherited from; the iron-adequate column is used",
      "here per the operator ruling (sidecar q1 = B). See the vignette Errata."
    )
  )

  ini({
    # ---- Iron-disposition parameters ----
    # Loss clearance. NOT estimated on the rat data; inherited from the mouse
    # iron-adequate column (Table 2) per the operator ruling (sidecar q1 = B).
    lcl <- log(1.647e-4); label("Unavoidable physiologic iron loss clearance from plasma (L/h)")  # Table 2 CLloss, iron-adequate mouse column

    # Partition coefficients ESTIMATED on the rat tissue data (Results text:
    # "bone, liver, spleen, heart, muscle and kidney ... were 8.64, 59.4,
    # 147.9, 17.64, 2.28 and 6.41, respectively"). The human model re-uses
    # these unchanged ("assumed to be identical among tissues in rats").
    lkp_bone    <- log(8.64); label("Bone-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate
    lkp_heart   <- log(17.64); label("Heart-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate
    lkp_kidney  <- log(6.41); label("Kidney-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate
    lkp_liver   <- log(59.4); label("Liver-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate
    lkp_muscle  <- log(2.28); label("Muscle-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate
    lkp_spleen  <- log(147.9); label("Spleen-to-plasma partition coefficient (unitless)")  # Results text, IDA rat estimate

    # Partition coefficients for the tissues that were NOT measured. Inherited
    # from the mouse iron-adequate column (sidecar q1 = B); the paper says only
    # that they "were assumed to be identical ... in mice".
    lkp_adipose <- log(0.258); label("Fat (adipose)-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column
    lkp_brain   <- log(1.087); label("Brain-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column
    lkp_gut     <- log(5.349); label("Gut-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column
    lkp_lung    <- log(13.7); label("Lung-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column
    lkp_other   <- log(2.30e-9); label("Remainder (rest of body)-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column
    lkp_skin    <- log(4.651); label("Skin-to-plasma partition coefficient (unitless)")  # Table 2, iron-adequate mouse column

    # NOT REPORTED for the rat anywhere on disk. Derived by the operator ruling
    # (sidecar q2 = B) from the paper's OWN two lifespan anchors: mouse
    # iron-adequate TRBC = 34.44 h at 0.025 kg and human TRBC = 120 d = 2880 h
    # at 73 kg fix an allometric exponent b_T = log(2880/34.44)/log(73/0.025) =
    # 0.5547, evaluated at the paper's rat body weight of 0.345 kg -> 147.7 h.
    lmtt_rbc <- log(34.44 * (0.345 / 0.025)^(log(2880 / 34.44) / log(73 / 0.025))); label("Mean red-blood-cell lifespan (h)")
    # NOT REPORTED for the rat. Scaled from the mouse iron-adequate QE by the
    # paper's own RBC-lifespan allometry, Eq. (9) with b = 0.75, using the
    # derived rat TRBC above -> 6.47e-4 L/h. Operator sidecar q2 = B.
    lq_bone_rbc <- log(0.217e-3 * (0.345 / 0.025)^(0.75 * log(2880 / 34.44) / log(73 / 0.025))); label("Erythropoietic iron utilisation flow, bone to red cells (L/h)")

    # ---- Rat physiology, fixed to Supplementary Table 1 ----
    qc <- fixed(6.624); label("Cardiac output (L/h)")  # Methods; equals the tabulated lung blood flow
    q_adipose <- fixed(0.0265); label("Fat (adipose) blood flow (L/h)")  # Supplementary Table 1
    q_bone    <- fixed(0.168); label("Bone blood flow (L/h)")  # Supplementary Table 1
    q_brain   <- fixed(0.132); label("Brain blood flow (L/h)")  # Supplementary Table 1
    q_gut     <- fixed(0.498); label("Gut blood flow (L/h)")  # Supplementary Table 1
    q_heart   <- fixed(0.260); label("Heart blood flow (L/h)")  # Supplementary Table 1
    q_kidney  <- fixed(0.611); label("Kidney blood flow (L/h)")  # Supplementary Table 1
    q_liver   <- fixed(0.782); label("Liver blood flow (L/h)")  # Supplementary Table 1
    q_muscle  <- fixed(0.134); label("Muscle blood flow (L/h)")  # Supplementary Table 1
    q_skin    <- fixed(0.386); label("Skin blood flow (L/h)")  # Supplementary Table 1
    q_spleen  <- fixed(0.0417); label("Spleen blood flow (L/h)")  # Supplementary Table 1

    v_adipose <- fixed(10.84 / 1000); label("Fat (adipose) volume (L)")  # Supplementary Table 1
    v_bone    <- fixed(17.13 / 1000); label("Bone volume (L)")  # Supplementary Table 1
    v_brain   <- fixed(11.17 / 1000); label("Brain volume (L)")  # Supplementary Table 1
    v_gut     <- fixed(10.84 / 1000); label("Gut volume (L)")  # Supplementary Table 1
    v_heart   <- fixed(0.87 / 1000); label("Heart volume (L)")  # Supplementary Table 1
    v_kidney  <- fixed(2.49 / 1000); label("Kidney volume (L)")  # Supplementary Table 1
    v_liver   <- fixed(11.17 / 1000); label("Liver volume (L)")  # Supplementary Table 1
    v_lung    <- fixed(1.08 / 1000); label("Lung volume (L)")  # Supplementary Table 1
    v_muscle  <- fixed(132.28 / 1000); label("Muscle volume (L)")  # Supplementary Table 1
    v_other   <- fixed(27.27 / 1000); label("Remainder (rest of body) volume (L)")  # Supplementary Table 1
    v_skin    <- fixed(43.37 / 1000); label("Skin volume (L)")  # Supplementary Table 1
    v_spleen  <- fixed(0.65 / 1000); label("Spleen volume (L)")  # Supplementary Table 1
    v_plasma  <- fixed(7.8 / 1000); label("Plasma volume (L)")  # Supplementary Table 1

    # ---- Residual error. The paper reports NO residual estimate for the rat
    # model; the control stream in the supplement is the mouse run. The mouse
    # $SIGMA(1) proportional initial (variance 0.1) is carried forward. The
    # mouse additive arms are on a nmol AMOUNT scale and are not transferable
    # to this model's mg scale, so proportional error is used throughout.
    propSd <- sqrt(0.1); label("Proportional residual error, serum iron (fraction)")  # carried from mouse $SIGMA(1) initial
    propSd_Cbone    <- sqrt(0.1); label("Proportional residual error, bone (fraction)")  # carried from mouse $SIGMA initial
    propSd_Cheart   <- sqrt(0.1); label("Proportional residual error, heart (fraction)")  # carried from mouse $SIGMA initial
    propSd_Ckidney  <- sqrt(0.1); label("Proportional residual error, kidney (fraction)")  # carried from mouse $SIGMA initial
    propSd_Cliver   <- sqrt(0.1); label("Proportional residual error, liver (fraction)")  # carried from mouse $SIGMA initial
    propSd_Cmuscle  <- sqrt(0.1); label("Proportional residual error, muscle (fraction)")  # carried from mouse $SIGMA initial
    propSd_Cspleen  <- sqrt(0.1); label("Proportional residual error, spleen (fraction)")  # carried from mouse $SIGMA initial
  })

  model({
    # ---- 1. Iron-disposition parameters back-transformed ----
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

    # ---- 2. Rest-of-body blood flow ----
    # Computed as in the fitted mouse $PK block, so that the plasma efflux sum
    # below equals the cardiac output exactly. Gut and spleen are NOT subtracted
    # again (they are inside q_liver). Supplementary Table 1's printed Remainder flow of
    # 3.584 L/h double-subtracts them and breaks the plasma mass balance;
    # see the vignette Errata.
    q_lung  <- qc
    q_other <- qc - q_bone - q_heart - q_kidney - q_liver - q_muscle -
      q_skin - q_adipose - q_brain

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
    # Circulating iron-carbohydrate complex: the whole IV dose lands here and is
    # taken up by splenic macrophages at KA := q_spleen (Methods; the direct
    # IVIP-to-plasma release limb P1 was removed by the authors as ~0.1% of
    # dose). Its concentration is referenced to the plasma volume.
    d/dt(ivip) <- -q_spleen * (ivip / v_plasma)
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
    # Spleen: macrophage uptake of the FCM colloid (the IVIP release limb) plus
    # erythrophagocytic recycling of senescent RBC iron (Eq. 4)
    d/dt(spleen) <- q_spleen * Cp + q_spleen * (ivip / v_plasma) +
      rbc_iron / mtt_rbc - (q_spleen / kp_spleen) * Cspleen
    # Red-cell haemoglobin iron pool; amount-only state, no volume (Eq. 3)
    d/dt(rbc_iron) <- q_bone_rbc * Cbone - rbc_iron / mtt_rbc

    # ---- 5. Observations ----
    # Total serum iron, the quantity a clinical iron assay reports: released
    # iron PLUS iron still bound in the circulating FCM colloid. Cc is the
    # primary observable; Cserum_released exposes the released fraction alone.
    Cc <- (plasma + ivip) / v_plasma
    Cserum_released <- Cp
    Cc ~ prop(propSd)
    Cbone    ~ prop(propSd_Cbone)
    Cheart   ~ prop(propSd_Cheart)
    Ckidney  ~ prop(propSd_Ckidney)
    Cliver   ~ prop(propSd_Cliver)
    Cmuscle  ~ prop(propSd_Cmuscle)
    Cspleen  ~ prop(propSd_Cspleen)
  })
}
