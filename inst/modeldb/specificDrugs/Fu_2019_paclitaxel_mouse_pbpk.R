Fu_2019_paclitaxel_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, Phoenix WinNonlin 8.3 PML).",
    "Paclitaxel (Taxol, Cremophor EL formulation) disposition in female FVB",
    "mice after a single 20 mg/kg IV bolus (Fu 2019). Separate venous and",
    "arterial blood pools with the lung in series between them, eight",
    "perfusion-limited tissues (spleen, liver, kidney, heart, gut, muscle,",
    "adipose, brain) and a permeability-limited remainder split into a",
    "vascular/interstitial (is_remainder, 33.3%) and an intracellular",
    "(int_remainder, 66.7%) subspace. Spleen and gut drain into the liver;",
    "hepatic metabolism (fu * CLint acting on the liver venous concentration)",
    "is the only elimination route. Plasma unbound fraction fixed to 0.05.",
    "Partition coefficients, CLint and PS are the Table 1 final mouse",
    "estimates; ODE structure and 30-g physiology from Supplementary Code S1.",
    "Typical-value model; no IIV was reported."
  )
  reference <- paste(
    "Fu Q, Sun X, Lustberg MB, Sparreboom A, Hu S. Predicting Paclitaxel",
    "Disposition in Humans With Whole-Body Physiologically-Based",
    "Pharmacokinetic Modeling. CPT Pharmacometrics Syst Pharmacol.",
    "2019;8(12):931-939. doi:10.1002/psp4.12472.",
    "Supplementary Table S1 (physiology) and Code S1 (Phoenix PML code for",
    "the final PBPK model)."
  )
  vignette <- "Fu_2019_paclitaxel_pbpk"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    venous = list(analyte = "paclitaxel", units = "ug", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "paclitaxel", units = "ug", specimen = "plasma", verified = TRUE),
    spleen = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    is_remainder = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    int_remainder = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "paclitaxel", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (FVB, female)",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "10-14 weeks",
    weight_range = "23-29 g (model physiology fixed at a 30-g reference mouse)",
    sex_female_pct = 100,
    disease_state = "Healthy (non-tumor-bearing) mice",
    dose_range = "20 mg/kg paclitaxel (Taxol) single IV bolus via the tail vein",
    regions = "USA",
    notes = paste(
      "Four animals per sacrifice time (0.5, 1, 4, 8 and 24 h); plasma plus",
      "brain, fat, colon, cecum, small intestine, stomach, liver, kidneys,",
      "lungs, spleen and heart were assayed by HPLC-UV (Fu 2019 Methods).",
      "n_subjects = 4 animals x 5 time points. Physiology (organ volume",
      "fractions of body weight and blood-flow fractions of cardiac output)",
      "is taken from Code S1 for a 30-g mouse (Table S1)."
    )
  )

  ini({
    # Perfusion-limited tissue:plasma partition coefficients -- Fu 2019 Table 1,
    # column 'PBPK estimated' (final mouse estimates, CV% 13-14%).
    lkp_spleen <- log(0.83); label("Spleen tissue:plasma partition coefficient Ks (unitless)") # Table 1 Spleen Ks = 0.83
    lkp_liver <- log(2.8); label("Liver tissue:plasma partition coefficient Kl (unitless)") # Table 1 Liver Kl = 2.8
    lkp_kidney <- log(1); label("Kidney tissue:plasma partition coefficient Kk (unitless)") # Table 1 Kidney Kk = 1
    lkp_heart <- log(0.54); label("Heart tissue:plasma partition coefficient Kh (unitless)") # Table 1 Heart Kh = 0.54
    lkp_lung <- log(0.77); label("Lung tissue:plasma partition coefficient Klu (unitless)") # Table 1 Lung Klu = 0.77
    lkp_gut <- log(1.61); label("Gut tissue:plasma partition coefficient Kg (unitless)") # Table 1 Gut Kg = 1.61
    lkp_muscle <- log(0.37); label("Muscle tissue:plasma partition coefficient Km (unitless)") # Table 1 Muscle Km = 0.37
    lkp_adipose <- log(0.58); label("Fat tissue:plasma partition coefficient Kf (unitless)") # Table 1 Fat Kf = 0.58
    lkp_brain <- log(0.03); label("Brain tissue:plasma partition coefficient Kbr (unitless)") # Table 1 Brain Kbr = 0.03

    # Permeability-limited remainder -- Table 1 'Remainder' rows. Kr is the
    # vascular/interstitial-side coefficient (Code S1 'Kr'); K_ISF is the
    # intracellular-side coefficient (Code S1 'Krt'). Mapping by value match
    # with the Code S1 human block (Kr 0.528, Krt 2.606); see vignette.
    lkp_remainder <- log(0.52); label("Remainder vascular/interstitial partition coefficient Kr (unitless)") # Table 1 Remainder Kr = 0.52
    lkp_int_remainder <- log(2.63); label("Remainder intracellular partition coefficient K_ISF / Krt (unitless)") # Table 1 Remainder K_ISF = 2.63
    lps_remainder <- log(43); label("Remainder permeability-surface area product PSr (mL/h)") # Table 1 PS r = 43 mL/hour

    # Hepatic intrinsic clearance
    lclint <- log(79.72); label("Hepatic intrinsic clearance CLint,H (mL/h)") # Table 1 Cl int,H = 79.72 mL/hour

    # Plasma unbound fraction, fixed (Methods; Brouwer 2000 ref 18)
    fu <- fixed(0.05); label("Plasma unbound fraction fu (unitless)") # Methods 'the fu of paclitaxel was fixed to a value of 0.05'; Code S1 fu = 0.05

    # Residual error: Code S1 mouse block, error(CEps2 = 0.218255) with
    # observe(C*Obs = C * (1 + CEps2)) -- proportional, shared across plasma
    # and all tissues.
    propSd <- 0.218255; label("Proportional residual error (fraction)") # Code S1 mouse block CEps2 = 0.218255
  })

  model({
    # ===== Mouse physiology (Code S1 mouse block; Table S1, 30-g mouse) =====
    bw <- 0.03 # kg, Code S1 BW = 0.03

    # Organ volume fractions of body weight (Code S1 FV_*)
    fv_spleen <- 0.00275
    fv_liver <- 0.0413
    fv_kidney <- 0.0125
    fv_gut <- 0.0318
    fv_heart <- 0.00375
    fv_lung <- 0.0055
    fv_blood <- 0.037
    fv_muscle <- 0.2875
    fv_adipose <- 0.06475
    fv_brain <- 0.0125
    fv_remainder <- 1 - fv_spleen - fv_liver - fv_kidney - fv_gut - fv_heart -
      fv_lung - fv_blood - fv_muscle - fv_adipose - fv_brain

    # Organ volumes (mL) = BW (kg) * FV * 1000
    v_spleen <- bw * fv_spleen * 1000
    v_liver <- bw * fv_liver * 1000
    v_kidney <- bw * fv_kidney * 1000
    v_gut <- bw * fv_gut * 1000
    v_heart <- bw * fv_heart * 1000
    v_lung <- bw * fv_lung * 1000
    v_blood <- bw * fv_blood * 1000
    v_muscle <- bw * fv_muscle * 1000
    v_adipose <- bw * fv_adipose * 1000
    v_brain <- bw * fv_brain * 1000
    v_remainder <- bw * fv_remainder * 1000

    # Blood split 75% venous / 25% arterial; remainder split 33.3% vascular +
    # interstitial / 66.7% intracellular (Code S1 Vbla, Vblb, Vrv, Vre)
    v_venous <- 0.75 * v_blood
    v_arterial <- 0.25 * v_blood
    v_is_remainder <- 0.333 * v_remainder
    v_int_remainder <- 0.667 * v_remainder

    # Cardiac output (mL/h) and blood-flow fractions (Code S1 CO, FQ_*)
    co <- 885
    fq_spleen <- 0.0048
    fq_liver <- 0.131
    fq_gut <- 0.114
    fq_kidney <- 0.0739
    fq_heart <- 0.0536
    fq_muscle <- 0.130
    fq_lung <- 0.814
    fq_adipose <- 0.057
    fq_brain <- 0.027
    fq_remainder <- 1 - fq_liver - fq_spleen - fq_gut - fq_kidney - fq_heart -
      fq_muscle - fq_adipose - fq_brain

    q_spleen <- co * fq_spleen
    q_liver <- co * fq_liver
    q_gut <- co * fq_gut
    q_kidney <- co * fq_kidney
    q_heart <- co * fq_heart
    q_lung <- co * fq_lung
    q_muscle <- co * fq_muscle
    q_adipose <- co * fq_adipose
    q_remainder <- co * fq_remainder
    q_brain <- co * fq_brain

    # ===== Individual parameters =====
    kp_spleen <- exp(lkp_spleen)
    kp_liver <- exp(lkp_liver)
    kp_kidney <- exp(lkp_kidney)
    kp_heart <- exp(lkp_heart)
    kp_lung <- exp(lkp_lung)
    kp_gut <- exp(lkp_gut)
    kp_muscle <- exp(lkp_muscle)
    kp_adipose <- exp(lkp_adipose)
    kp_brain <- exp(lkp_brain)
    kp_remainder <- exp(lkp_remainder)
    kp_int_remainder <- exp(lkp_int_remainder)
    ps_remainder <- exp(lps_remainder)
    clint <- exp(lclint)

    # ===== Concentrations (ug/mL) =====
    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_lung <- lung / v_lung
    c_spleen <- spleen / v_spleen
    c_liver <- liver / v_liver
    c_kidney <- kidney / v_kidney
    c_heart <- heart / v_heart
    c_gut <- gut / v_gut
    c_muscle <- muscle / v_muscle
    c_adipose <- adipose / v_adipose
    c_brain <- brain / v_brain
    c_is_remainder <- is_remainder / v_is_remainder
    c_int_remainder <- int_remainder / v_int_remainder

    # ===== ODEs (Code S1; paper Eqs. 1-7) =====
    # Venous blood (Code S1 Ab; dose point)
    d/dt(venous) <- -q_lung * c_venous + q_liver * c_liver / kp_liver +
      q_remainder * c_is_remainder / kp_remainder +
      q_kidney * c_kidney / kp_kidney + q_heart * c_heart / kp_heart +
      q_muscle * c_muscle / kp_muscle + q_adipose * c_adipose / kp_adipose +
      q_brain * c_brain / kp_brain
    # Lung, in series between venous and arterial blood (Eq. 5)
    d/dt(lung) <- q_lung * c_venous - q_lung * c_lung / kp_lung
    # Arterial blood (Eq. 6). The flows written out of the arterial pool are
    # reproduced term by term from Code S1.
    d/dt(arterial) <- q_lung * c_lung / kp_lung -
      (q_liver - q_spleen - q_gut) * c_arterial - q_remainder * c_arterial -
      q_kidney * c_arterial - q_heart * c_arterial - q_muscle * c_arterial -
      q_adipose * c_arterial - q_gut * c_arterial - q_spleen * c_arterial -
      q_brain * c_arterial
    # Spleen and gut drain into the liver
    d/dt(spleen) <- q_spleen * c_arterial - q_spleen * c_spleen / kp_spleen
    d/dt(gut) <- q_gut * c_arterial - q_gut * c_gut / kp_gut
    # Liver, the eliminating organ (Eq. 2 with the Code S1 signs)
    d/dt(liver) <- (q_liver - q_spleen - q_gut) * c_arterial -
      q_liver * c_liver / kp_liver - fu * clint * c_liver / kp_liver +
      q_spleen * c_spleen / kp_spleen + q_gut * c_gut / kp_gut
    # Permeability-limited remainder (Eqs. 3-4)
    d/dt(is_remainder) <- q_remainder *
      (c_arterial - c_is_remainder / kp_remainder) -
      ps_remainder * fu *
        (c_is_remainder / kp_remainder - c_int_remainder / kp_int_remainder)
    d/dt(int_remainder) <- ps_remainder * fu *
      (c_is_remainder / kp_remainder - c_int_remainder / kp_int_remainder)
    # Perfusion-limited tissues (Eq. 1)
    d/dt(kidney) <- q_kidney * c_arterial - q_kidney * c_kidney / kp_kidney
    d/dt(heart) <- q_heart * c_arterial - q_heart * c_heart / kp_heart
    d/dt(muscle) <- q_muscle * c_arterial - q_muscle * c_muscle / kp_muscle
    d/dt(adipose) <- q_adipose * c_arterial - q_adipose * c_adipose / kp_adipose
    d/dt(brain) <- q_brain * (c_arterial - c_brain / kp_brain)

    # Observed plasma concentration is the arterial pool (Code S1 Cp = Ca)
    Cc <- c_arterial
    Cc ~ prop(propSd)
  })
}
