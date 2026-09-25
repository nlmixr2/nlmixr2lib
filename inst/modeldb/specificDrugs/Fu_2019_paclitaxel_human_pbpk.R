Fu_2019_paclitaxel_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, Phoenix WinNonlin 8.3 PML). Paclitaxel (Taxol,",
    "Cremophor EL formulation) disposition in adult cancer patients after",
    "175 mg/m^2 as a 3-h IV infusion, scaled up from the mouse model in",
    "Fu_2019_paclitaxel_mouse_pbpk (Fu 2019). Same topology: separate venous",
    "and arterial blood pools with the lung in series between them, eight",
    "perfusion-limited tissues (spleen, liver, kidney, heart, gut, muscle,",
    "adipose, brain) and a permeability-limited remainder split into a",
    "vascular/interstitial (is_remainder, 33.3%) and an intracellular",
    "(int_remainder, 66.7%) subspace. Spleen and gut drain into the liver;",
    "hepatic metabolism (fu * CLint acting on the liver venous concentration)",
    "is the only elimination route; fu fixed to 0.05. Physiology is a",
    "70-kg reference adult (Table S1). All parameter values are from the",
    "human block of Supplementary Code S1 ('Code for the final PBPK model').",
    "The 1,410 L/h CLint quoted in the Methods text is not used: it does not",
    "reproduce the observed plasma profile (Figure 3; see vignette).",
    "Typical-value model; no IIV was reported."
  )
  reference <- paste(
    "Fu Q, Sun X, Lustberg MB, Sparreboom A, Hu S. Predicting Paclitaxel",
    "Disposition in Humans With Whole-Body Physiologically-Based",
    "Pharmacokinetic Modeling. CPT Pharmacometrics Syst Pharmacol.",
    "2019;8(12):931-939. doi:10.1002/psp4.12472.",
    "Supplementary Table S1 (physiology) and Code S1 (Phoenix PML code for",
    "the final PBPK model, human block)."
  )
  vignette <- "Fu_2019_paclitaxel_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    venous = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    spleen = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    is_remainder = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    int_remainder = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "paclitaxel", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  # BSA only converts the 175 mg/m^2 dose to an amount; the physiology is
  # fixed at the 70-kg reference adult and BSA does not enter model().
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "External dose-amount multiplier only (dose mg = 175 mg/m^2 * BSA).",
        "The Fu 2019 human PBPK model fixes organ volumes and blood flows at",
        "a 70-kg reference adult (Table S1, Code S1); BSA does not enter the",
        "ODE system."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    age_range = "> 18 years",
    weight_range = "not reported (model physiology fixed at a 70-kg reference adult)",
    sex_female_pct = NA_real_,
    disease_state = "Adults with confirmed solid tumors eligible for paclitaxel chemotherapy",
    dose_range = "175 mg/m^2 single-agent paclitaxel (Taxol) as a 3-h IV infusion",
    regions = "Netherlands (Rotterdam Cancer Institute)",
    notes = paste(
      "Plasma sampled pre-dose, at 1 and 2 h after the start of infusion,",
      "at the end of infusion and at 5, 15, 30, 45 min and 1, 2, 4, 6, 8,",
      "12 and 21 h after the end of infusion; patients previously reported",
      "by van Zuylen 2001 (Fu 2019 ref 13). Only plasma was observed in",
      "humans; tissue concentrations are model predictions (Figure 5a)."
    )
  )

  ini({
    # Code S1 human block fixef() values (the block below the dashed line).
    lkp_spleen <- log(0.833372); label("Spleen tissue:plasma partition coefficient Ks (unitless)") # Code S1 human tvKs = 0.833372
    lkp_liver <- log(2.82516); label("Liver tissue:plasma partition coefficient Kl (unitless)") # Code S1 human tvKl = 2.82516
    lkp_kidney <- log(1.01274); label("Kidney tissue:plasma partition coefficient Kk (unitless)") # Code S1 human tvKk = 1.01274
    lkp_heart <- log(0.541566); label("Heart tissue:plasma partition coefficient Kh (unitless)") # Code S1 human tvKh = 0.541566
    lkp_lung <- log(0.778136); label("Lung tissue:plasma partition coefficient Klu (unitless)") # Code S1 human tvKlu = 0.778136
    lkp_gut <- log(1.62708); label("Gut tissue:plasma partition coefficient Kg (unitless)") # Code S1 human tvKg = 1.62708
    lkp_muscle <- log(0.37604); label("Muscle tissue:plasma partition coefficient Km (unitless)") # Code S1 human tvKm = 0.37604
    lkp_adipose <- log(0.586303); label("Fat tissue:plasma partition coefficient Kf (unitless)") # Code S1 human tvKf = 0.586303
    lkp_brain <- log(0.035072); label("Brain tissue:plasma partition coefficient Kbr (unitless)") # Code S1 human tvKbr = 0.035072

    lkp_remainder <- log(0.527921); label("Remainder vascular/interstitial partition coefficient Kr (unitless)") # Code S1 human tvKr = 0.527921
    lkp_int_remainder <- log(2.60612); label("Remainder intracellular partition coefficient Krt (unitless)") # Code S1 human tvKrt = 2.60612
    lps_remainder <- log(61.7295); label("Remainder permeability-surface area product PSr (L/h)") # Code S1 human tvPSr = 61729.5 mL/h

    lclint <- log(449.886); label("Hepatic intrinsic clearance CLint,H (L/h)") # Code S1 human tvCll = 449886 mL/h

    fu <- fixed(0.05); label("Plasma unbound fraction fu (unitless)") # Methods 'the fu of paclitaxel was fixed to a value of 0.05'; Code S1 fu = 0.05

    propSd <- 0.252886; label("Proportional residual error (fraction)") # Code S1 human block CEps2 = 0.252886
  })

  model({
    # ===== Human physiology (Code S1 human block; Table S1, 70-kg adult) =====
    bw <- 70 # kg, Code S1 BW = 70

    # Organ volume fractions of body weight (Code S1 FV_*)
    fv_spleen <- 0.0027
    fv_liver <- 0.026
    fv_kidney <- 0.004
    fv_gut <- 0.017
    fv_heart <- 0.005
    fv_lung <- 0.008
    fv_blood <- 0.079
    fv_muscle <- 0.4
    fv_adipose <- 0.21
    fv_brain <- 0.02
    fv_remainder <- 1 - fv_spleen - fv_liver - fv_kidney - fv_gut - fv_heart -
      fv_lung - fv_blood - fv_muscle - fv_adipose - fv_brain

    # Organ volumes (L) = BW (kg) * FV (Code S1 multiplies by 1000 to mL)
    v_spleen <- bw * fv_spleen
    v_liver <- bw * fv_liver
    v_kidney <- bw * fv_kidney
    v_gut <- bw * fv_gut
    v_heart <- bw * fv_heart
    v_lung <- bw * fv_lung
    v_blood <- bw * fv_blood
    v_muscle <- bw * fv_muscle
    v_adipose <- bw * fv_adipose
    v_brain <- bw * fv_brain
    v_remainder <- bw * fv_remainder

    # Blood split 75% venous / 25% arterial; remainder split 33.3% vascular +
    # interstitial / 66.7% intracellular (Code S1 Vbla, Vblb, Vrv, Vre)
    v_venous <- 0.75 * v_blood
    v_arterial <- 0.25 * v_blood
    v_is_remainder <- 0.333 * v_remainder
    v_int_remainder <- 0.667 * v_remainder

    # Cardiac output (L/h; Code S1 CO = 394600 mL/h) and flow fractions
    co <- 394.6
    fq_spleen <- 0.0117
    fq_liver <- 0.221
    fq_gut <- 0.167
    fq_kidney <- 0.189
    fq_heart <- 0.0365
    fq_muscle <- 0.114
    fq_lung <- 0.846
    fq_adipose <- 0.0395
    fq_brain <- 0.106
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

    # ===== Concentrations (mg/L) =====
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
