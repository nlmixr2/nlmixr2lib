Fan_2025_iron_ferriccarboxymaltose_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 15-compartment, perfusion-limited; NONMEM ADVAN15).",
    "Human (adults with iron-deficiency anaemia).",
    "Interspecies extrapolation of the Fan 2025 mouse whole-body iron PBPK",
    "model to human, for the iron-carbohydrate complex preparation (IVIP)",
    "ferric carboxymaltose (FCM). The 14-compartment mouse structure is",
    "carried over unchanged and human physiology is fixed to the literature",
    "(Supplementary Table 2); one extra circulating `ivip` compartment receives the",
    "whole intravenous dose and is taken up by splenic macrophages at a rate",
    "set equal to the spleen blood flow (KA := q_spleen), the authors having",
    "removed the direct IVIP-to-plasma release limb as ~0.1% of dose. The",
    "mechanistic limbs are unchanged: erythropoietic iron utilisation from",
    "bone into the red-cell pool, erythrophagocytic return of that iron to",
    "the spleen over the red-cell lifespan, portal drainage of gut and",
    "spleen into liver, and a physiologic iron loss clearance from plasma.",
    "No parameter was estimated on the human data: the six kp values estimated in the rat are re-used unchanged, the remaining kp values and the loss clearance are the inherited mouse iron-adequate values, and QE follows Eq. (9). The human step is a pure forward prediction.",
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
    bone     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    gut      = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    heart    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    kidney   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    liver    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    muscle   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    skin     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    lung     = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    spleen   = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    rbc_iron = list(analyte = "iron", units = "mg", specimen = "blood cell", verified = FALSE),
    brain    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    adipose  = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE),
    other    = list(analyte = "iron", units = "mg", specimen = "tissue", verified = FALSE)
  )

  population <- list(
    species       = "human (adults with iron-deficiency anaemia)",
    n_subjects    = "not reported; digitised group mean data only",
    n_studies     = "1 (reference [25] of the paper)",
    weight_range  = "73 kg (single body weight used for all physiology)",
    disease_state = "iron-deficiency anaemia",
    dose_range    = "single ascending intravenous ferric carboxymaltose doses of 100, 500, 800 and 1000 mg Fe",
    notes         = paste(
      "No parameter was estimated on the human data: the six kp values estimated",
      "in the rat are re-used unchanged, the remaining kp values and the loss",
      "clearance are the inherited mouse iron-adequate values, and QE follows",
      "Eq. (9). The human step is a pure forward prediction."
    )
  )

  ini({
    # ---- Iron-disposition parameters ----
    # Loss clearance. NOT estimated on the human data; inherited from the mouse
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

    # Methods: "TRBC used in the human simulation was 120 days" -> 2880 h.
    lmtt_rbc <- log(120 * 24); label("Mean red-blood-cell lifespan (h)")
    # Eq. (9): QE(human) = QE(mouse) * (TRBC_human / TRBC_mouse)^0.75 with b = 0.75
    # fixed. QE(mouse) and TRBC(mouse) are the iron-adequate column per the
    # operator ruling (sidecar q1 = B) -> 6.00e-3 L/h.
    lq_bone_rbc <- log(0.217e-3 * (2880 / 34.44)^0.75); label("Erythropoietic iron utilisation flow, bone to red cells (L/h)")

    # ---- Human physiology, fixed to Supplementary Table 2 ----
    qc <- fixed(336); label("Cardiac output (L/h)")  # Methods; equals the tabulated lung blood flow
    q_adipose <- fixed(16.8); label("Fat (adipose) blood flow (L/h)")  # Supplementary Table 2
    q_bone    <- fixed(16.8); label("Bone blood flow (L/h)")  # Supplementary Table 2
    q_brain   <- fixed(40.32); label("Brain blood flow (L/h)")  # Supplementary Table 2
    q_gut     <- fixed(50.4); label("Gut blood flow (L/h)")  # Supplementary Table 2
    q_heart   <- fixed(13.44); label("Heart blood flow (L/h)")  # Supplementary Table 2
    q_kidney  <- fixed(63.84); label("Kidney blood flow (L/h)")  # Supplementary Table 2
    q_liver   <- fixed(21.84); label("Liver blood flow (L/h)")  # Supplementary Table 2
    q_muscle  <- fixed(56.99); label("Muscle blood flow (L/h)")  # Supplementary Table 2
    q_skin    <- fixed(16.8); label("Skin blood flow (L/h)")  # Supplementary Table 2
    q_spleen  <- fixed(10.08); label("Spleen blood flow (L/h)")  # Supplementary Table 2

    v_adipose <- fixed(12.5); label("Fat (adipose) volume (L)")  # Supplementary Table 2
    v_bone    <- fixed(10.5); label("Bone volume (L)")  # Supplementary Table 2
    v_brain   <- fixed(1.4); label("Brain volume (L)")  # Supplementary Table 2
    v_gut     <- fixed(1.2); label("Gut volume (L)")  # Supplementary Table 2
    v_heart   <- fixed(0.33); label("Heart volume (L)")  # Supplementary Table 2
    v_kidney  <- fixed(0.31); label("Kidney volume (L)")  # Supplementary Table 2
    v_liver   <- fixed(1.8); label("Liver volume (L)")  # Supplementary Table 2
    v_lung    <- fixed(0.47); label("Lung volume (L)")  # Supplementary Table 2
    v_muscle  <- fixed(30); label("Muscle volume (L)")  # Supplementary Table 2
    v_other   <- fixed(8.01); label("Remainder (rest of body) volume (L)")  # Supplementary Table 2
    v_skin    <- fixed(3.3); label("Skin volume (L)")  # Supplementary Table 2
    v_spleen  <- fixed(0.18); label("Spleen volume (L)")  # Supplementary Table 2
    v_plasma  <- fixed(3); label("Plasma volume (L)")  # Supplementary Table 2

    # ---- Residual error. The paper reports NO residual estimate for the human
    # model; the control stream in the supplement is the mouse run. The mouse
    # $SIGMA(1) proportional initial (variance 0.1) is carried forward. The
    # mouse additive arms are on a nmol AMOUNT scale and are not transferable
    # to this model's mg scale, so proportional error is used throughout.
    propSd <- sqrt(0.1); label("Proportional residual error, serum iron (fraction)")  # carried from mouse $SIGMA(1) initial
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
    # again (they are inside q_liver). Supplementary Table 2's printed Remainder flow of
    # 28.69 L/h double-subtracts them and breaks the plasma mass balance;
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
  })
}
