Wyska_2016_lisofylline_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited). Human projection of the Wyska 2016",
    "mouse lisofylline PBPK model to a 70 kg adult: the mouse structure and",
    "tissue-to-serum partition coefficients with human organ volumes and",
    "blood flows (Table 2), and linear hepatic metabolism driven by an in",
    "vivo intrinsic clearance of 1545.24 mL/min scaled from human",
    "liver-microsome Vmax/KM (Eq. 11). Intravenous doses enter the venous",
    "blood; oral doses are absorbed first-order from a depot into the venous",
    "blood using the mouse-derived ka 0.024 1/min and F 0.16. As in the",
    "mouse model, the remainder blood flow subtracts the gut flow twice, so",
    "Q_gut x C_arterial leaves the arterial blood without returning; the",
    "equations and Table 2 are kept as printed. A simulation-only",
    "projection: no variability and no residual error."
  )
  reference <- paste(
    "Wyska E, Swierczek A, Pociecha K, Przejczowska-Pomierny K. Physiologically",
    "based modeling of lisofylline pharmacokinetics following intravenous",
    "administration in mice. Eur J Drug Metab Pharmacokinet.",
    "2016;41(4):403-412. doi:10.1007/s13318-015-0260-y"
  )
  vignette <- "Wyska_2016_lisofylline_pbpk"
  dosing <- c("depot", "venous")
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "lisofylline", units = "ug", specimen = "administration site", verified = TRUE),
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
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    weight_range = "70 kg reference adult (Table 2)",
    disease_state = paste(
      "Simulated reference adult; no human data were fitted. The",
      "intravenous projection was compared with digitised serum",
      "concentrations after a 3 mg/kg 10-min infusion from Bursten et al.",
      "1998."
    ),
    dose_range = "3 mg/kg 10-min intravenous infusion; 6 mg/kg oral (simulated).",
    notes = paste(
      "Projection from the mouse PBPK model of the same paper. Partition",
      "coefficients were assumed not to differ between species. The",
      "intrinsic clearance sums the CYP3A4, CYP2A6 (lisofylline diol) and",
      "CYP1A2 plus non-P450 (pentoxifylline) microsomal pathways scaled with",
      "45 mg microsomal protein per g liver and a 1800 g liver."
    )
  )

  ini({
    # Organ weights or volumes (g or mL) for a 70 kg adult, Table 2
    v_brain <- fixed(1400)
    label("Brain volume (mL)") # Table 2 row 'Brain' weight 1,400
    v_muscle <- fixed(35000)
    label("Muscle volume (mL)") # Table 2 row 'Muscle' weight 35,000
    v_kidney <- fixed(310)
    label("Kidneys volume (mL)") # Table 2 row 'Kidneys' weight 310
    v_liver <- fixed(1800)
    label("Liver volume (mL)") # Table 2 row 'Liver' weight 1,800
    v_gut <- fixed(1650)
    label("Gut volume (mL)") # Table 2 row 'Gut' weight 1,650
    v_lung <- fixed(1000)
    label("Lungs volume (mL)") # Table 2 row 'Lungs' weight 1,000
    v_other <- fixed(24290)
    label("Remainder volume (mL)") # Table 2 row 'Remainder' weight 24,290
    v_arterial <- fixed(1000)
    label("Arterial blood volume (mL)") # Table 2 row 'Arterial blood' 1,000
    v_venous <- fixed(3550)
    label("Venous blood volume (mL)") # Table 2 row 'Venous blood' 3,550

    # Blood flows (mL/min), Table 2
    q_brain <- fixed(700)
    label("Brain blood flow (mL/min)") # Table 2 row 'Brain' blood flow 700
    q_muscle <- fixed(750)
    label("Muscle blood flow (mL/min)") # Table 2 row 'Muscle' blood flow 750
    q_kidney <- fixed(1240)
    label("Kidneys blood flow (mL/min)") # Table 2 row 'Kidneys' blood flow 1,240
    q_liver <- fixed(1450)
    label("Total liver blood flow, hepatic artery plus portal (mL/min)") # Table 2 row 'Liver' blood flow 1,450; Eq. 4 feeds the liver with C_A x (Q_Li - Q_Gu) plus the gut outflow
    q_gut <- fixed(1100)
    label("Gut blood flow (mL/min)") # Table 2 row 'Gut' blood flow 1,100
    q_lung <- fixed(5600)
    label("Lung blood flow, cardiac output (mL/min)") # Table 2 row 'Lungs' blood flow 5,600
    q_other <- fixed(360)
    label("Remainder blood flow (mL/min)") # Table 2 row 'Remainder' blood flow 360 (= 5600 minus all other flows, gut included)

    # Tissue-to-serum partition coefficients, taken unchanged from the mouse
    # (Methods: 'It was assumed that P values did not differ between species')
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

    # Elimination
    lclint <- fixed(log(1545.24))
    label("In vivo whole-liver unbound intrinsic clearance CLint (mL/min)") # Results text: 'calculated based on these data and Eq. 11 was 1,545.24 ml/min'
    cl_renal <- fixed(0)
    label("Renal clearance CLR (mL/min)") # Methods: only hepatic CLint used for the human simulations; CLR fixed to 0 as in the mouse (Table 3)
    fu <- fixed(0.54)
    label("Unbound fraction in human plasma applied to the liver outflow concentration (unitless)") # Methods text after Eq. 12: 'f u in human plasma is 0.54'

    # Oral absorption (Eq. 13), mouse values used for the dashed lines of Fig. 5
    lka <- fixed(log(0.024))
    label("First-order absorption rate constant ka (1/min)") # Methods after Eq. 13: ka 0.024 min-1 (Wyska et al. 2007, mice)
    lfdepot <- fixed(log(0.16))
    label("Oral bioavailability F (fraction)") # Methods after Eq. 13: F 0.16 (Wyska et al. 2007, mice)

    propSd <- fixed(0)
    label("Proportional residual error (fraction); none in the simulation-only projection")
  })

  model({
    clint <- exp(lclint)
    ka <- exp(lka)
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

    # Eq. 4 with the saturable term replaced by the single linear CLint
    # (Methods: 'a single parameter, CLint in vivo was used in Eq. 4')
    r_metab <- clint * fu * c_liver / kp_liver

    # Eq. 13: In(t) = F x Dose x ka x exp(-ka t), entering the venous blood
    d/dt(depot) <- -ka * depot
    d/dt(venous) <- ka * depot +
      q_brain * c_brain / kp_brain +
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

    f(depot) <- exp(lfdepot)

    Cbrain <- c_brain
    Cmuscle <- c_muscle
    Ckidney <- c_kidney
    Cliver <- c_liver
    Cgut <- c_gut
    Clung <- c_lung
    Carterial <- c_arterial
    Cother <- c_other

    Cc <- c_venous
    Cc ~ prop(propSd)
  })
}
