Chen_2016_deoxypodophyllotoxin_dog_pbpk <- function() {
  description <- paste0(
    "Preclinical (beagle dog). PBPK (whole-body, 13 tissue compartments, ",
    "coded in Phoenix WinNonlin). Deoxypodophyllotoxin (DPT) disposition in ",
    "the beagle dog after intravenous bolus dosing (Chen et al. 2016, Front ",
    "Pharmacol). Adipose, liver, muscle, lung, kidney, brain, heart, spleen, ",
    "skin, gastrointestinal tract and rest of body, plus arterial and venous ",
    "blood pools. All tissues are perfusion-rate limited except the brain, ",
    "which is permeability-limited and split into vascular (3% of brain ",
    "volume) and extravascular (97%) sub-compartments coupled by a ",
    "permeability-surface product PS. The gastrointestinal tract and spleen ",
    "drain into the liver rather than into venous blood. Elimination is ",
    "hepatic only, as the sum of two in-vitro metabolic pathways: ",
    "Michaelis-Menten formation of M2 and auto-activating Hill-type formation ",
    "of M7, both driven by the unbound hepatic outflow concentration and ",
    "scaled to the whole liver by the microsomal protein amount PBSF. The ",
    "blood-to-plasma ratio is assumed to be 1. The dog is the outlier ",
    "species: its very low unbound fraction (0.16%) and low hepatic metabolic ",
    "velocity give a clearance only about 3% of hepatic blood flow, and the ",
    "dog was excluded from the paper's interspecies allometric scaling. Dog ",
    "M2 formation needed two Michaelis-Menten sites. Deterministic: the ",
    "publication reports no inter-individual variance estimates and no ",
    "residual-error magnitude, so the model is intended for typical-value ",
    "simulation."
  )
  reference <- paste0(
    "Chen Y, Zhao K, Liu F, Xie Q, Zhong Z, Miao M, Liu X, Liu L. ",
    "Prediction of Deoxypodophyllotoxin Disposition in Mouse, Rat, ",
    "Monkey, and Dog by Physiologically Based Pharmacokinetic Model ",
    "and the Extrapolation to Human. Front Pharmacol. 2016;7:488. ",
    "doi:10.3389/fphar.2016.00488"
  )
  vignette <- "Chen_2016_deoxypodophyllotoxin"
  units <- list(
    time = "min",
    dosing = "ug",
    concentration = "ug/mL",
    weight = "kg"
  )

  # Issue #482: what each ODE state holds. Every state is an AMOUNT of DPT in
  # ug; the model() block divides by the Table 1 volume to obtain the tissue
  # concentration.
  compartmentData <- list(
    venous = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "whole blood", verified = TRUE),
    adipose = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "whole blood", verified = TRUE),
    brain_extravascular = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE),
    other = list(analyte = "deoxypodophyllotoxin", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "dog (beagle, male and female)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "about 10 months",
    weight_range = "8.65 +/- 0.64 kg",
    sex_female_pct = 50,
    disease_state = "Healthy beagle dogs, fasted 12 h before dosing.",
    dose_range = "Intravenous bolus 0.3 mg/kg via the cephalic vein.",
    regions = "China Pharmaceutical University, Nanjing, China",
    notes = "Six beagles (three male, three female). Blood sampled pre-dose and at 5, 15, 30, 60, 120, 180, 240, 300, 360, 420, 480, 600, 720 and 840 min post-dose -- the longest sampling window of the four animal species, because dog clearance is roughly 30-fold lower than the other species. The dog was excluded from the paper's interspecies allometric scaling for this reason."
  )

  ini({
    # ----------------------------------------------------------------
    # Organ volumes (mL) -- Table 1, "Dog (8.5 kg)" column, "Volume" sub-column.
    # Every physiological constant is FIXED: the paper takes them from the
    # cited physiology literature and estimates none of them.
    # ----------------------------------------------------------------
    v_adipose <- fixed(1500)
    label("Volume of adipose (mL)") # Table 1 row 'Adipose'
    v_liver <- fixed(213)
    label("Volume of liver (mL)") # Table 1 row 'Liver'
    v_muscle <- fixed(4250)
    label("Volume of muscle (mL)") # Table 1 row 'Muscle'
    v_lung <- fixed(85)
    label("Volume of lungs (mL)") # Table 1 row 'Lungs'
    v_kidney <- fixed(97)
    label("Volume of kidneys (mL)") # Table 1 row 'Kidneys'
    v_brain <- fixed(50)
    label("Volume of brain (mL)") # Table 1 row 'Brain'
    v_heart <- fixed(43)
    label("Volume of heart (mL)") # Table 1 row 'Heart'
    v_spleen <- fixed(22)
    label("Volume of spleen (mL)") # Table 1 row 'Spleen'
    v_skin <- fixed(364)
    label("Volume of skin (mL)") # Table 1 row 'Skin'
    v_gut <- fixed(228)
    label("Volume of gastrointestinal tract (mL)") # Table 1 row 'Gastrointestinal tract'
    v_other <- fixed(1223)
    label("Volume of rest of body (mL)") # Table 1 row 'Rest of body'
    v_venous <- fixed(284)
    label("Volume of venous blood (mL)") # Table 1 row 'Vein'
    v_arterial <- fixed(141)
    label("Volume of arterial blood (mL)") # Table 1 row 'Artery'

    # ----------------------------------------------------------------
    # Organ blood flows (mL/min) -- Table 1, "Dog (8.5 kg)" column, "Blood flow"
    # sub-column. q_lung is the cardiac output (see qc below); q_liver is
    # the TOTAL hepatic flow, which Table 1 footnote h defines as the sum
    # of hepatic-artery, gastrointestinal and splenic flow -- so the
    # hepatic-artery flow is derived in model() by subtraction.
    # ----------------------------------------------------------------
    q_adipose <- fixed(50)
    label("Blood flow to adipose (mL/min)") # Table 1 row 'Adipose'
    q_liver <- fixed(323.33)
    label("Blood flow to liver (mL/min)") # Table 1 row 'Liver'
    q_muscle <- fixed(170)
    label("Blood flow to muscle (mL/min)") # Table 1 row 'Muscle'
    q_lung <- fixed(968.33)
    label("Blood flow to lungs (mL/min)") # Table 1 row 'Lungs'
    q_kidney <- fixed(170)
    label("Blood flow to kidneys (mL/min)") # Table 1 row 'Kidneys'
    q_brain <- fixed(145)
    label("Blood flow to brain (mL/min)") # Table 1 row 'Brain'
    q_heart <- fixed(43.33)
    label("Blood flow to heart (mL/min)") # Table 1 row 'Heart'
    q_spleen <- fixed(13.33)
    label("Blood flow to spleen (mL/min)") # Table 1 row 'Spleen'
    q_skin <- fixed(18.33)
    label("Blood flow to skin (mL/min)") # Table 1 row 'Skin'
    q_gut <- fixed(265)
    label("Blood flow to gastrointestinal tract (mL/min)") # Table 1 row 'Gastrointestinal tract'
    q_other <- fixed(48.33)
    label("Blood flow to rest of body (mL/min)") # Table 1 row 'Rest of body'

    # ----------------------------------------------------------------
    # Tissue-to-plasma concentration ratios Kt:pl -- Table 1, "Dog (8.5 kg)"
    # column, "Kt:pl" sub-column. The rat values were computed from rat
    # tissue composition by the Ruark 2014 method and the beagle dog values
    # follow from the paper's identical-Kt:pl,u assumption,
    # Kt:pl = Kt:pl,rat * fu / fu_rat (Results, "PBPK Model Development
    # and Validation").
    # ----------------------------------------------------------------
    lkp_adipose <- fixed(log(0.64))
    label("Log adipose-to-plasma concentration ratio (unitless)") # Table 1 row 'Adipose' Kt:pl = 0.64
    lkp_liver <- fixed(log(0.05))
    label("Log liver-to-plasma concentration ratio (unitless)") # Table 1 row 'Liver' Kt:pl = 0.05
    lkp_muscle <- fixed(log(0.02))
    label("Log muscle-to-plasma concentration ratio (unitless)") # Table 1 row 'Muscle' Kt:pl = 0.02
    lkp_lung <- fixed(log(0.05))
    label("Log lungs-to-plasma concentration ratio (unitless)") # Table 1 row 'Lungs' Kt:pl = 0.05
    lkp_kidney <- fixed(log(0.04))
    label("Log kidneys-to-plasma concentration ratio (unitless)") # Table 1 row 'Kidneys' Kt:pl = 0.04
    lkp_brain <- fixed(log(0.08))
    label("Log brain-to-plasma concentration ratio (unitless)") # Table 1 row 'Brain' Kt:pl = 0.08
    lkp_heart <- fixed(log(0.03))
    label("Log heart-to-plasma concentration ratio (unitless)") # Table 1 row 'Heart' Kt:pl = 0.03
    lkp_spleen <- fixed(log(0.03))
    label("Log spleen-to-plasma concentration ratio (unitless)") # Table 1 row 'Spleen' Kt:pl = 0.03
    lkp_skin <- fixed(log(0.04))
    label("Log skin-to-plasma concentration ratio (unitless)") # Table 1 row 'Skin' Kt:pl = 0.04
    lkp_gut <- fixed(log(0.03))
    label("Log gastrointestinal tract-to-plasma concentration ratio (unitless)") # Table 1 row 'Gastrointestinal tract' Kt:pl = 0.03
    lkp_other <- fixed(log(0.0003))
    label("Log rest of body-to-plasma concentration ratio (unitless)") # Table 1 row 'Rest of body' Kt:pl = 0.0003

    # ----------------------------------------------------------------
    # Brain permeability-surface area product.
    # ----------------------------------------------------------------
    ps_brain <- fixed(0.1407)
    label("Brain permeability-surface area product PS (mL/min)")
    # Table 1 footnote row 'PS (mL/min)', Dog column. Scaled from the fitted mouse value as PS_i = PS_mouse * (W_i / W_mouse)^0.67.

    # ----------------------------------------------------------------
    # Drug-specific constants.
    # ----------------------------------------------------------------
    fu <- fixed(0.0016)
    label("Fraction of DPT unbound in beagle dog plasma (unitless)")
    # Results 'Plasma Protein Binding of DPT in Five Species': dog plasma protein binding 99.84 +/- 0.03%, so fu = 1 - 0.9984; the Abstract and Results quote the same value as '0.16% in dog'
    rbp <- fixed(1)
    label("Blood-to-plasma concentration ratio Rbp (unitless)")
    # Methods 'PBPK Model Development': Rbp "was assumed to be unity in the study, according to the reasons provided in Poulin and Theil (2002)"
    mw <- fixed(398.4)
    label("Deoxypodophyllotoxin molecular weight (g/mol)")
    # NOT printed in Chen 2016. Chemical constant for the DPT molecular formula C22H22O7 (PubChem CID 73435, 398.4 g/mol), needed only as the unit bridge between the ug/mL concentration scale of the ODEs and the uM scale of the Table 2 Km values. Not a fitted model parameter.

    # ----------------------------------------------------------------
    # Hepatic metabolism -- Table 2. Elimination is assumed
    # to occur only in liver, only via formation of metabolites M2 and M7
    # (Methods "PBPK Model Development"). The M2 kinetic parameters are
    # cited from the authors' companion in vitro paper (Xie et al. 2016);
    # the M7 kinetics were measured in the present study and are
    # auto-activating, hence the Hill form.
    # ----------------------------------------------------------------
    pbsf <- fixed(16592.7)
    label("Total hepatic microsomal protein per body (mg protein)")
    # Table 2 row 'PBSF'; Table 2 footnote: microsomal protein yield x liver weight
    km1_m2 <- fixed(0.09)
    label("Michaelis constant of the first M2-formation site (uM)") # Table 2 row 'Km1,M2'
    vmax1_m2 <- fixed(0.03)
    label("Vmax of the first M2-formation site (nmol/min/mg protein)") # Table 2 row 'Vmax1,M2'
    km2_m2 <- fixed(9.47)
    label("Michaelis constant of the second M2-formation site (uM)") # Table 2 row 'Km2,M2'
    vmax2_m2 <- fixed(0.1)
    label("Vmax of the second M2-formation site (nmol/min/mg protein)") # Table 2 row 'Vmax2,M2'
    km_m7 <- fixed(0.52)
    label("Michaelis constant of M7 formation (uM)") # Table 2 row 'Km,M7'
    vmax_m7 <- fixed(2.64)
    label("Vmax of M7 formation (pmol/min/mg protein)")
    # Table 2 row 'Vmax,M7'. NOTE the unit: Vmax,M7 is printed in pmol/min/mg protein while Vmax,M2 is printed in nmol/min/mg protein, so model() divides this by 1000 before adding the two pathways. Keeping the printed pmol value here preserves the source trace; the Discussion sanity-checks the conversion by stating CL_M2 is "at least fifty times larger than CL_M7".
    hill_m7 <- fixed(1.45)
    label("Hill coefficient of the auto-activating M7 formation (unitless)") # Table 2 row 'gamma'

    # ----------------------------------------------------------------
    # Residual error. The "Visual Predictive Checks of the Model" section
    # states that a multiplicative (proportional) residual-error model was
    # fitted in Phoenix NLME, but neither its magnitude nor the estimated
    # inter-individual variances of hepatic blood flow and metabolic
    # velocity are reported anywhere in the paper or its figures. They are
    # therefore encoded as zero rather than invented; see the vignette
    # Errata. The PBPK predictions the paper validates against Table 3 are
    # deterministic typical-value simulations in any case.
    # ----------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - magnitude not reported in source)")
  })

  model({
    # ----------------------------------------------------------------
    # Cardiac output. Table 1 gives the lung blood flow, which the
    # Methods equation for the lung compartment calls Qtotal.
    # ----------------------------------------------------------------
    qc <- q_lung

    # Table 1 footnote h: "The blood flow rate of liver was assumed to be
    # the sum of blood flow rates of hepatic artery, gastrointestinal
    # tract, and spleen." q_liver is that sum, so the hepatic-artery
    # share is the remainder.
    q_hepatic_artery <- q_liver - q_gut - q_spleen

    # Table 1 footnote i: "The vascular and extravascular volumes of brain
    # were assumed to account for 3 and 97% of the total brain volume".
    v_brain_vascular <- 0.03 * v_brain
    v_brain_extravascular <- 0.97 * v_brain
    ps <- ps_brain

    # ----------------------------------------------------------------
    # Tissue-to-BLOOD ratios. The Methods ODEs are written with the
    # grouping Kt:pl / Rbp throughout, so the venous concentration
    # leaving a tissue is Ct / (Kt:pl / Rbp). Rbp is 1 in this paper, but
    # the division is kept explicit so the encoded equations match the
    # printed ones term for term.
    # ----------------------------------------------------------------
    kp_adipose <- exp(lkp_adipose)
    kp_liver <- exp(lkp_liver)
    kp_muscle <- exp(lkp_muscle)
    kp_lung <- exp(lkp_lung)
    kp_kidney <- exp(lkp_kidney)
    kp_brain <- exp(lkp_brain)
    kp_heart <- exp(lkp_heart)
    kp_spleen <- exp(lkp_spleen)
    kp_skin <- exp(lkp_skin)
    kp_gut <- exp(lkp_gut)
    kp_other <- exp(lkp_other)

    # Tissue-to-BLOOD ratios, Kt:pl / Rbp.
    kb_adipose <- kp_adipose / rbp
    kb_liver <- kp_liver / rbp
    kb_muscle <- kp_muscle / rbp
    kb_lung <- kp_lung / rbp
    kb_kidney <- kp_kidney / rbp
    kb_brain <- kp_brain / rbp
    kb_heart <- kp_heart / rbp
    kb_spleen <- kp_spleen / rbp
    kb_skin <- kp_skin / rbp
    kb_gut <- kp_gut / rbp
    kb_other <- kp_other / rbp
    kb_brain_ev <- kp_brain # the brain ODEs are printed with a bare Kbra:pl, not Kbra:pl / Rbp

    # ----------------------------------------------------------------
    # Concentrations (ug/mL) = amount (ug) / volume (mL).
    # ----------------------------------------------------------------
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_lung <- lung / v_lung
    c_adipose <- adipose / v_adipose
    c_liver <- liver / v_liver
    c_muscle <- muscle / v_muscle
    c_kidney <- kidney / v_kidney
    c_heart <- heart / v_heart
    c_spleen <- spleen / v_spleen
    c_skin <- skin / v_skin
    c_gut <- gut / v_gut
    c_other <- other / v_other
    c_brain_vascular <- brain_vascular / v_brain_vascular
    c_brain_extravascular <- brain_extravascular / v_brain_extravascular

    # Total brain concentration, the quantity a brain homogenate assay
    # measures: vascular plus extravascular amount over whole brain volume.
    c_brain <- (brain_vascular + brain_extravascular) / v_brain

    # ----------------------------------------------------------------
    # Hepatic elimination. Both clearance terms are driven by the unbound
    # plasma-equivalent concentration at the liver, fu * Cliv / Kliv:pl,
    # expressed in uM to match the Table 2 Km values:
    #   Cu [uM] = Cu [ug/mL] * 1000 / MW [g/mol]
    #
    # CL_M2 = sum_i Vmax,i,M2 * Cu / (Km,i,M2 + Cu)   (Results, M2 equation)
    # CL_M7 = Vmax,M7 * Cu^gamma / (Km,M7^gamma + Cu^gamma)  (Results, M7
    #   equation; a Hill form because "the kinetics exhibited
    #   auto-activation features")
    #
    # Vmax,M7 is printed in pmol/min/mg protein and Vmax,M2 in
    # nmol/min/mg protein, hence the 1000-fold conversion before summing.
    # The summed velocity is per mg microsomal protein, so it is scaled to
    # the whole liver by PBSF, giving nmol/min, then converted to the
    # amount unit of the ODEs: nmol/min -> ug/min (1 nmol = mw * 1e-3 ug).
    # ----------------------------------------------------------------
    cu_liver <- fu * (c_liver / kb_liver) * 1000 / mw
    cl_m2 <- vmax1_m2 * cu_liver / (km1_m2 + cu_liver) +
      vmax2_m2 * cu_liver / (km2_m2 + cu_liver)
    cl_m7 <- (vmax_m7 / 1000) * cu_liver^hill_m7 /
      (km_m7^hill_m7 + cu_liver^hill_m7)
    rate_metabolism <- pbsf * (cl_m2 + cl_m7) * mw / 1000

    # ----------------------------------------------------------------
    # Venous blood. Methods: Vven dCven/dt = sum_t (Qt * Ct / (Kt:pl/Rbp))
    # - Qtotal * Cven. The gastrointestinal tract and spleen are NOT in
    # this sum: their outflow feeds the liver (Figure 2 and the liver
    # equation), and only the liver outflow reaches venous blood. The
    # brain contributes Qbra * C1, its VASCULAR concentration, because the
    # permeability-limited brain equation writes its perfusion term as
    # Qbra * (Cart - C1); that is the mass balance the refined brain model
    # forces, replacing the generic Ct / (Kt:pl/Rbp) term.
    # Dose the IV bolus into this compartment.
    # ----------------------------------------------------------------
    d/dt(venous) <- q_adipose * (c_adipose / kb_adipose) +
      q_muscle * (c_muscle / kb_muscle) +
      q_kidney * (c_kidney / kb_kidney) +
      q_heart * (c_heart / kb_heart) +
      q_skin * (c_skin / kb_skin) +
      q_other * (c_other / kb_other) +
      q_liver * (c_liver / kb_liver) +
      q_brain * c_brain_vascular -
      qc * c_venous

    # Lung. Methods: Vlun dClun/dt = Qtotal * (Cven - Clun / (Klun:pl/Rbp)).
    d/dt(lung) <- qc * (c_venous - c_lung / kb_lung)

    # Arterial blood. Methods: Vart dCart/dt = Qtotal * (Clun /
    # (Klun:pl/Rbp) - Cart).
    d/dt(arterial) <- qc * (c_lung / kb_lung - c_arterial)

    # Perfusion-rate limited tissues. Methods: Vt dCt/dt = Qt * (Cart -
    # Ct / (Kt:pl/Rbp)).
    d/dt(adipose) <- q_adipose * (c_arterial - c_adipose / kb_adipose)
    d/dt(muscle) <- q_muscle * (c_arterial - c_muscle / kb_muscle)
    d/dt(kidney) <- q_kidney * (c_arterial - c_kidney / kb_kidney)
    d/dt(heart) <- q_heart * (c_arterial - c_heart / kb_heart)
    d/dt(spleen) <- q_spleen * (c_arterial - c_spleen / kb_spleen)
    d/dt(skin) <- q_skin * (c_arterial - c_skin / kb_skin)
    d/dt(gut) <- q_gut * (c_arterial - c_gut / kb_gut)
    d/dt(other) <- q_other * (c_arterial - c_other / kb_other)

    # Brain, permeability-limited. Results: "Preliminary prediction of DPT
    # concentrations in brain of mice demonstrated a significant
    # overestimation using the perfusion-rate limited model, therefore
    # brain model was refined to be a permeability-limited model":
    #   V1 dC1/dt = Qbra (Cart - C1) - PS (C1 - C2 / Kbra:pl)
    #   V2 dC2/dt = PS (C1 - C2 / Kbra:pl)
    d/dt(brain_vascular) <- q_brain * (c_arterial - c_brain_vascular) -
      ps * (c_brain_vascular - c_brain_extravascular / kb_brain_ev)
    d/dt(brain_extravascular) <-
      ps * (c_brain_vascular - c_brain_extravascular / kb_brain_ev)

    # Liver, the only eliminating organ. Methods:
    #   Vliv dCliv/dt = Qhep * Cart + Qgut * Cgut / (Kgut:pl/Rbp)
    #     + Qspl * Cspl / (Kspl:pl/Rbp)
    #     - (Qhep + Qgut + Qspl) * Cliv / (Kliv:pl/Rbp)
    #     - PBSF * (CL_M2 + CL_M7)
    # with (Qhep + Qgut + Qspl) = q_liver by Table 1 footnote h.
    d/dt(liver) <- q_hepatic_artery * c_arterial +
      q_gut * (c_gut / kb_gut) +
      q_spleen * (c_spleen / kb_spleen) -
      q_liver * (c_liver / kb_liver) -
      rate_metabolism

    # ----------------------------------------------------------------
    # Observation. The paper measures and reports PLASMA concentrations;
    # c_venous is a venous BLOOD concentration, so plasma = blood / Rbp.
    # With the paper's Rbp = 1 the two coincide.
    # ----------------------------------------------------------------
    Cc <- c_venous / rbp
    Cc ~ prop(propSd)
  })
}
