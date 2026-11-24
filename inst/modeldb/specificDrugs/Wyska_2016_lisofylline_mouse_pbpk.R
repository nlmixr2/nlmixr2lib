Wyska_2016_lisofylline_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, flow-limited, ADAPT II) model of",
    "lisofylline in male CD-1 mice after a single 50 mg/kg intravenous bolus",
    "(Wyska 2016). Brain, muscle, kidneys, liver, gut and lungs plus a lumped",
    "remainder, arterial and venous blood compartments; the gut drains into",
    "the liver. Elimination is Michaelis-Menten hepatic metabolism of the",
    "unbound liver outflow concentration (Vmax and KM estimated); renal",
    "clearance was fixed to zero. Organ volumes, blood flows and the",
    "AUC-ratio tissue-to-serum partition coefficients are fixed (Table 1).",
    "The published remainder blood flow subtracts the gut flow twice, so",
    "Q_gut x C_arterial leaves the arterial blood without returning; the",
    "model keeps the equations and Table 1 exactly as printed. Fitted to",
    "naive-pooled mean profiles, so there is no between-animal variability;",
    "the residual-error model is the ADAPT variance model",
    "SD = sigma1 + sigma2 x Y on the serum concentration."
  )
  reference <- paste(
    "Wyska E, Swierczek A, Pociecha K, Przejczowska-Pomierny K. Physiologically",
    "based modeling of lisofylline pharmacokinetics following intravenous",
    "administration in mice. Eur J Drug Metab Pharmacokinet.",
    "2016;41(4):403-412. doi:10.1007/s13318-015-0260-y"
  )
  vignette <- "Wyska_2016_lisofylline_pbpk"
  dosing <- "venous"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    venous = list(analyte = "lisofylline", units = "ug", specimen = "serum", verified = TRUE),
    lung = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "lisofylline", units = "ug", specimen = "serum", verified = TRUE),
    brain = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE),
    other = list(analyte = "lisofylline", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (male Crl:CD-1)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "8-10 weeks",
    weight_range = "30.25 +/- 2.11 g (mean +/- SD); model physiology for a 30 g mouse",
    sex_female_pct = 0,
    disease_state = "Healthy male CD-1 mice, fasted overnight before dosing.",
    dose_range = "50 mg/kg single intravenous bolus into the tail vein.",
    regions = "Jagiellonian University Medical College, Krakow, Poland",
    notes = paste(
      "Destructive sampling with 3-4 mice exsanguinated per time point up",
      "to 60 min; serum, brain, muscle, kidneys, liver, gut and lungs",
      "assayed by chiral HPLC-UV (LLOQ 0.05 ug/mL). The total number of",
      "animals is not reported. The 7 mean concentration-time profiles were",
      "fitted simultaneously by maximum likelihood in ADAPT II. A separate",
      "metabolic-cage study (3 mice, 50 mg/kg IV) recovered 0.71 +/- 0.13%",
      "of the dose unchanged in urine over 4 h, supporting CLR = 0."
    )
  )

  ini({
    # Organ weights or volumes (g or mL) for a 30 g mouse, Table 1
    v_brain <- fixed(0.50)
    label("Brain volume (mL)") # Table 1 row 'Brain' weight 0.50
    v_muscle <- fixed(11.40)
    label("Muscle volume (mL)") # Table 1 row 'Muscle' weight 11.40
    v_kidney <- fixed(0.50)
    label("Kidneys volume (mL)") # Table 1 row 'Kidneys' weight 0.50
    v_liver <- fixed(1.65)
    label("Liver volume (mL)") # Table 1 row 'Liver' weight 1.65
    v_gut <- fixed(1.50)
    label("Gut volume (mL)") # Table 1 row 'Gut' weight 1.50
    v_lung <- fixed(0.22)
    label("Lungs volume (mL)") # Table 1 row 'Lungs' weight 0.22
    v_other <- fixed(13.05)
    label("Remainder volume (mL)") # Table 1 row 'Remainder' weight 13.05 (30 g minus the other tissues, footnote b)
    v_arterial <- fixed(0.24)
    label("Arterial blood volume (mL)") # Table 1 row 'Arterial blood' 0.24
    v_venous <- fixed(0.94)
    label("Venous blood volume (mL)") # Table 1 row 'Venous blood' 0.94

    # Blood flows (mL/min), Table 1
    q_brain <- fixed(0.47)
    label("Brain blood flow (mL/min)") # Table 1 row 'Brain' blood flow 0.47
    q_muscle <- fixed(0.91)
    label("Muscle blood flow (mL/min)") # Table 1 row 'Muscle' blood flow 0.91
    q_kidney <- fixed(1.27)
    label("Kidneys blood flow (mL/min)") # Table 1 row 'Kidneys' blood flow 1.27
    q_liver <- fixed(2.25)
    label("Total liver blood flow, hepatic artery plus portal (mL/min)") # Table 1 row 'Liver' blood flow 2.25; Eq. 4 feeds the liver with C_A x (Q_Li - Q_Gu) plus the gut outflow
    q_gut <- fixed(1.50)
    label("Gut blood flow (mL/min)") # Table 1 row 'Gut' blood flow 1.50
    q_lung <- fixed(13.98)
    label("Lung blood flow, cardiac output (mL/min)") # Table 1 row 'Lungs' blood flow 13.98
    q_other <- fixed(7.58)
    label("Remainder blood flow (mL/min)") # Table 1 row 'Remainder' blood flow 7.58 (footnote a: Q_Lu minus the sum of all other flows, gut included)

    # Tissue-to-serum partition coefficients (AUC-ratio method), Table 1
    lkp_brain <- fixed(log(0.345))
    label("Log brain-to-serum partition coefficient (unitless)") # Table 1 row 'Brain' P 0.345
    lkp_muscle <- fixed(log(0.640))
    label("Log muscle-to-serum partition coefficient (unitless)") # Table 1 row 'Muscle' P 0.640
    lkp_kidney <- fixed(log(0.777))
    label("Log kidneys-to-serum partition coefficient (unitless)") # Table 1 row 'Kidneys' P 0.777
    lkp_liver <- fixed(log(0.859))
    label("Log liver-to-serum partition coefficient (unitless)") # Table 1 row 'Liver' P 0.859
    lkp_gut <- fixed(log(0.892))
    label("Log gut-to-serum partition coefficient (unitless)") # Table 1 row 'Gut' P 0.892
    lkp_lung <- fixed(log(0.758))
    label("Log lungs-to-serum partition coefficient (unitless)") # Table 1 row 'Lungs' P 0.758

    # Elimination, Table 3
    lvmax <- log(9.644)
    label("Maximum hepatic metabolic rate Vmax (ug/min)") # Table 3 'V max (ug/min)' 9.644, CV 30.34%
    lkm <- log(1.777)
    label("Michaelis-Menten constant KM (ug/mL)") # Table 3 'K M (ug/ml)' 1.777, CV 89.52%
    cl_renal <- fixed(0)
    label("Renal clearance CLR (mL/min)") # Table 3 'CLR (ml/min)' 0 (fixed)
    fu <- fixed(0.76)
    label("Unbound fraction applied to the liver outflow concentration (fuT, unitless)") # Methods text: fraction unbound in murine serum 0.76; unbound liver-water and serum concentrations assumed equal

    # Residual error: var(C) = (sigma1 + sigma2 x Y)^2, Eq. 10
    addSd <- 0.311
    label("Additive residual error sigma1 (ug/mL)") # Table 3 'sigma 1' 0.311, CV 15.50%
    propSd <- 0.250
    label("Proportional residual error sigma2 (fraction)") # Table 3 'sigma 2' 0.250, CV 10.86%
  })

  model({
    vmax <- exp(lvmax)
    km <- exp(lkm)
    kp_brain <- exp(lkp_brain)
    kp_muscle <- exp(lkp_muscle)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_gut <- exp(lkp_gut)
    kp_lung <- exp(lkp_lung)

    c_venous <- venous / v_venous
    c_arterial <- arterial / v_arterial
    c_brain <- brain / v_brain
    c_muscle <- muscle / v_muscle
    c_kidney <- kidney / v_kidney
    c_liver <- liver / v_liver
    c_gut <- gut / v_gut
    c_lung <- lung / v_lung
    c_other <- other / v_other

    # Eq. 4 metabolic term, Vmax x fuT x C_Li / (fuT x C_Li + KM x P_Li)
    r_metab <- vmax * fu * c_liver / (fu * c_liver + km * kp_liver)

    # Eqs. 1-9 in amounts (V x dC/dt = dA/dt); the remainder (Eq. 7) carries no
    # partition coefficient. Drug input In(t) enters the venous blood (Eq. 9).
    d/dt(venous) <- q_brain * c_brain / kp_brain +
      q_muscle * c_muscle / kp_muscle +
      q_kidney * c_kidney / kp_kidney +
      q_liver * c_liver / kp_liver +
      q_other * c_other -
      q_lung * c_venous
    d/dt(lung) <- q_lung * (c_venous - c_lung / kp_lung)
    d/dt(arterial) <- q_lung * (c_lung / kp_lung - c_arterial)
    d/dt(brain) <- q_brain * (c_arterial - c_brain / kp_brain)
    d/dt(muscle) <- q_muscle * (c_arterial - c_muscle / kp_muscle)
    d/dt(kidney) <- q_kidney * (c_arterial - c_kidney / kp_kidney) -
      cl_renal * c_arterial
    d/dt(liver) <- c_arterial * (q_liver - q_gut) +
      q_gut * c_gut / kp_gut -
      q_liver * c_liver / kp_liver -
      r_metab
    d/dt(gut) <- q_gut * (c_arterial - c_gut / kp_gut)
    d/dt(other) <- q_other * (c_arterial - c_other)

    # Tissue concentrations in ug/g (1 g = 1 mL); serum is the venous plasma.
    Cbrain <- c_brain
    Cmuscle <- c_muscle
    Ckidney <- c_kidney
    Cliver <- c_liver
    Cgut <- c_gut
    Clung <- c_lung
    Carterial <- c_arterial
    Cother <- c_other

    Cc <- c_venous
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
