Tan_2021_cefadroxil_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, perfusion-limited, GastroPlus 9.7). Cefadroxil",
    "disposition in the healthy adult human (Tan 2021), encoded as the",
    "disposition layer only. Thirteen well-stirred, perfusion-rate-limited",
    "tissues (lung, adipose, muscle, liver, spleen, heart, brain, kidney, skin,",
    "reproductive organ, red marrow, yellow marrow and rest of body) are strung",
    "between arterial and venous blood pools; the spleen and the gut drain into",
    "the liver through the portal vein, and the liver additionally receives the",
    "hepatic artery. Every tissue volume and blood flow is the paper's Table 1",
    "human column (the GastroPlus default fasted human), the tissue-to-plasma",
    "partition coefficients are the species-invariant Table 3 values, the",
    "blood-to-plasma ratio is 1 and the only elimination is a renal clearance",
    "of 8.50 L/h acting on the arterial plasma that perfuses the kidney. The",
    "paper's oral absorption layer (the GastroPlus ACAT intestinal model with",
    "hPEPT1-mediated transport) is NOT reproduced: its segment volumes, transit",
    "times and absorption scale factors are platform-internal and unprinted.",
    "Doses enter the venous blood directly, which corresponds to the paper's",
    "simulated 99.9% fraction absorbed but not to its absorption time course.",
    "The ACAT gut tissue has a printed blood flow but no printed volume or",
    "partition coefficient, so it is encoded as a zero-volume pass-through.",
    "Deterministic typical-value model: the paper reports no between-subject",
    "variability and no residual-error model. Companion models:",
    "Tan_2021_cefadroxil_mouse_pbpk, Tan_2021_cefadroxil_rat_pbpk.",
    sep = " "
  )
  reference <- paste(
    "Tan Z, Zhang Y, Wang C, Sun L. Physiologically Based Pharmacokinetic",
    "Modeling of Cefadroxil in Mouse, Rat, and Human to Predict",
    "Concentration-Time Profile at Infected Tissue. Front Pharmacol.",
    "2021;12:692741. doi:10.3389/fphar.2021.692741. PMCID PMC8733657.",
    "Tissue volumes and blood flows are Table 1 (human columns); partition",
    "coefficients are Table 3 ('The K p values used in PBPK models' column);",
    "the blood-to-plasma ratio, doses, renal clearance and Vss are Table 2;",
    "predicted Cmax and AUC0-inf are Table 4; the predicted tissue Cmax values",
    "after 15 mg/kg are Results ('The Physiologically Based Pharmacokinetic",
    "Models in Mouse, Rat, and Human'); the fraction absorbed is Results",
    "('Effect of hPEPT1 on the Oral Absorption of Cefadroxil').",
    sep = " "
  )
  vignette <- "Tan_2021_cefadroxil_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Intravenous bolus into venous blood, the GastroPlus IV-bolus entry point.
  dosing <- "venous"

  # Every state is an AMOUNT of cefadroxil (mg); model() divides by the
  # Table 1 volume to obtain each concentration.
  compartmentData <- list(
    venous = list(analyte = "cefadroxil", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "cefadroxil", units = "mg", specimen = "whole blood", verified = TRUE),
    adipose = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    reproductive_organ = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    red_marrow = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    yellow_marrow = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "cefadroxil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    weight_range = "about 75 kg (implied: 375.5 mg is 5 mg/kg)",
    disease_state = "Healthy volunteers; the physiology is the GastroPlus default fasted human.",
    dose_range = "Oral 5, 15 and 30 mg/kg (375.5, 1126.5 and 2253 mg; Table 2) as an immediate-release suspension in the paper; encoded here as the absorbed dose entering venous blood.",
    regions = "Not stated; observed data from Garrigues et al. 1991.",
    notes = paste(
      "The model is not a population fit: it is a GastroPlus default-physiology",
      "forward simulation whose renal clearance was optimised against one dose",
      "level and validated on the other two (Discussion, Model Calibration).",
      "Six volunteers received 5 and 15 mg/kg orally and three also received 30",
      "mg/kg (Methods, Model Validation). Partition coefficients were measured",
      "in mice and rats and assumed identical in humans. The plasma unbound",
      "fraction (71.9%, Table 2) does not enter the perfusion-limited equations",
      "because both the partition coefficients and the clearance are referenced",
      "to total plasma. The hPEPT1 kinetic constants (Km 860.31 mg/L, Vmax",
      "0.0025 ug/s/mg hPEPT1) and the per-segment hPEPT1 contents belong to the",
      "unreproduced absorption layer.",
      sep = " "
    )
  )

  ini({
    # Tissue volumes (mL) -- Table 1, 'Human' 'Volume (ml)' column.
    v_lung <- fixed(1125.5); label("Volume of lung (mL)") # Table 1 row 'Lung'
    v_arterial <- fixed(2148.17); label("Volume of arterial blood (mL)") # Table 1 row 'Arterial supply'
    v_venous <- fixed(4296.34); label("Volume of venous blood (mL)") # Table 1 row 'Venous return'
    v_adipose <- fixed(24307.3); label("Volume of adipose (mL)") # Table 1 row 'Adipose'
    v_muscle <- fixed(25174.7); label("Volume of muscle (mL)") # Table 1 row 'Muscle'
    v_liver <- fixed(1590.1); label("Volume of liver (mL)") # Table 1 row 'Liver'
    v_spleen <- fixed(169.973); label("Volume of spleen (mL)") # Table 1 row 'Spleen'
    v_heart <- fixed(337.26); label("Volume of heart (mL)") # Table 1 row 'Heart'
    v_brain <- fixed(1493.16); label("Volume of brain (mL)") # Table 1 row 'Brain'
    v_kidney <- fixed(360.411); label("Volume of kidney (mL)") # Table 1 row 'Kidney'
    v_skin <- fixed(2831.73); label("Volume of skin (mL)") # Table 1 row 'Skin'
    v_reproductive_organ <- fixed(50.6174); label("Volume of reproductive organ (mL)") # Table 1 row 'ReproOrg'
    v_red_marrow <- fixed(1106.48); label("Volume of red marrow (mL)") # Table 1 row 'Red marrow'
    v_yellow_marrow <- fixed(3075.64); label("Volume of yellow marrow (mL)") # Table 1 row 'Yellow marrow'
    v_other <- fixed(2681.07); label("Volume of rest of body (mL)") # Table 1 row 'RestOfBody'

    # Tissue blood flows (mL/s) -- Table 1, 'Human' 'The blood flow (ml/s)'
    # column. q_lung is the cardiac output; q_liver is the TOTAL hepatic
    # outflow (hepatic artery + gut + spleen), which is what makes the flows
    # balance exactly against cardiac output, so the hepatic artery is derived
    # in model() by subtraction. q_gut is the ACAT gut row.
    q_lung <- fixed(98.1779); label("Blood flow to lung, cardiac output (mL/s)") # Table 1 row 'Lung'
    q_adipose <- fixed(8.09433); label("Blood flow to adipose (mL/s)") # Table 1 row 'Adipose'
    q_muscle <- fixed(12.5874); label("Blood flow to muscle (mL/s)") # Table 1 row 'Muscle'
    q_liver <- fixed(24.3444); label("Total hepatic blood outflow (mL/s)") # Table 1 row 'Liver'
    q_gut <- fixed(12.263); label("Blood flow to the gut (portal tributary) (mL/s)") # Table 1 row 'ACAT gut'
    q_spleen <- fixed(2.83294); label("Blood flow to spleen (mL/s)") # Table 1 row 'Spleen'
    q_heart <- fixed(4.10344); label("Blood flow to heart (mL/s)") # Table 1 row 'Heart'
    q_brain <- fixed(12.6919); label("Blood flow to brain (mL/s)") # Table 1 row 'Brain'
    q_kidney <- fixed(22.1051); label("Blood flow to kidney (mL/s)") # Table 1 row 'Kidney'
    q_skin <- fixed(5.66346); label("Blood flow to skin (mL/s)") # Table 1 row 'Skin'
    q_reproductive_organ <- fixed(0.17716); label("Blood flow to reproductive organ (mL/s)") # Table 1 row 'ReproOrg'
    q_red_marrow <- fixed(5.53242); label("Blood flow to red marrow (mL/s)") # Table 1 row 'Red marrow'
    q_yellow_marrow <- fixed(1.53782); label("Blood flow to yellow marrow (mL/s)") # Table 1 row 'Yellow marrow'
    q_other <- fixed(1.34054); label("Blood flow to rest of body (mL/s)") # Table 1 row 'RestOfBody'

    # Tissue-to-plasma partition coefficients -- Table 3, column 'The K p
    # values used in PBPK models'; the same values are used in all three
    # species (Methods).
    lkp_lung <- fixed(log(0.52)); label("Log lung-to-plasma partition coefficient (unitless)") # Table 3 row 'Lung' = 0.52
    lkp_adipose <- fixed(log(0.1)); label("Log adipose-to-plasma partition coefficient (unitless)") # Table 3 row 'Adipose' = 0.1
    lkp_muscle <- fixed(log(0.25)); label("Log muscle-to-plasma partition coefficient (unitless)") # Table 3 row 'Muscle' = 0.25
    lkp_liver <- fixed(log(0.81)); label("Log liver-to-plasma partition coefficient (unitless)") # Table 3 row 'Liver' = 0.81
    lkp_spleen <- fixed(log(0.41)); label("Log spleen-to-plasma partition coefficient (unitless)") # Table 3 row 'Spleen' = 0.41
    lkp_heart <- fixed(log(0.23)); label("Log heart-to-plasma partition coefficient (unitless)") # Table 3 row 'Heart' = 0.23
    lkp_brain <- fixed(log(0.11)); label("Log brain-to-plasma partition coefficient (unitless)") # Table 3 row 'Brain' = 0.11
    lkp_kidney <- fixed(log(6.87)); label("Log kidney-to-plasma partition coefficient (unitless)") # Table 3 row 'Kidney' = 6.87
    lkp_skin <- fixed(log(0.89)); label("Log skin-to-plasma partition coefficient (unitless)") # Table 3 row 'Skin' = 0.89
    lkp_reproductive_organ <- fixed(log(0.5)); label("Log reproductive-organ-to-plasma partition coefficient (unitless)") # Table 3 row 'Reproductive organ' = 0.5
    lkp_red_marrow <- fixed(log(0.47)); label("Log red-marrow-to-plasma partition coefficient (unitless)") # Table 3 row 'Red marrow' = 0.47
    lkp_yellow_marrow <- fixed(log(0.25)); label("Log yellow-marrow-to-plasma partition coefficient (unitless)") # Table 3 row 'Yellow marrow' = 0.25
    lkp_other <- fixed(log(0.5)); label("Log rest-of-body-to-plasma partition coefficient (unitless)") # Table 3 row 'Rest of body' = 0.5

    bp <- fixed(1); label("Blood-to-plasma concentration ratio (unitless)") # Table 2 row 'Blood/plasma conc. ratio' = 1

    # Renal clearance optimised against the observed profiles (Discussion,
    # Model Calibration); it acts on the arterial plasma perfusing the kidney.
    lcl_renal <- fixed(log(8.50)); label("Renal clearance (L/h)") # Table 2 row 'Clearance (L/h)', human 'Kidney: 8.50'
  })

  model({
    # Unit conversion: Table 1 volumes mL -> L, flows mL/s -> L/h.
    vl_lung <- v_lung / 1000
    vl_arterial <- v_arterial / 1000
    vl_venous <- v_venous / 1000
    vl_adipose <- v_adipose / 1000
    vl_muscle <- v_muscle / 1000
    vl_liver <- v_liver / 1000
    vl_spleen <- v_spleen / 1000
    vl_heart <- v_heart / 1000
    vl_brain <- v_brain / 1000
    vl_kidney <- v_kidney / 1000
    vl_skin <- v_skin / 1000
    vl_reproductive_organ <- v_reproductive_organ / 1000
    vl_red_marrow <- v_red_marrow / 1000
    vl_yellow_marrow <- v_yellow_marrow / 1000
    vl_other <- v_other / 1000

    ql_lung <- q_lung * 3.6
    ql_adipose <- q_adipose * 3.6
    ql_muscle <- q_muscle * 3.6
    ql_liver <- q_liver * 3.6
    ql_gut <- q_gut * 3.6
    ql_spleen <- q_spleen * 3.6
    ql_heart <- q_heart * 3.6
    ql_brain <- q_brain * 3.6
    ql_kidney <- q_kidney * 3.6
    ql_skin <- q_skin * 3.6
    ql_reproductive_organ <- q_reproductive_organ * 3.6
    ql_red_marrow <- q_red_marrow * 3.6
    ql_yellow_marrow <- q_yellow_marrow * 3.6
    ql_other <- q_other * 3.6

    # Hepatic artery = total hepatic outflow minus the two portal tributaries.
    ql_ha <- ql_liver - ql_gut - ql_spleen

    kp_lung <- exp(lkp_lung)
    kp_adipose <- exp(lkp_adipose)
    kp_muscle <- exp(lkp_muscle)
    kp_liver <- exp(lkp_liver)
    kp_spleen <- exp(lkp_spleen)
    kp_heart <- exp(lkp_heart)
    kp_brain <- exp(lkp_brain)
    kp_kidney <- exp(lkp_kidney)
    kp_skin <- exp(lkp_skin)
    kp_reproductive_organ <- exp(lkp_reproductive_organ)
    kp_red_marrow <- exp(lkp_red_marrow)
    kp_yellow_marrow <- exp(lkp_yellow_marrow)
    kp_other <- exp(lkp_other)
    cl_renal <- exp(lcl_renal)

    # Concentrations (mg/L). Blood pools hold whole blood; with bp = 1 whole
    # blood and plasma concentrations coincide.
    c_venous <- venous / vl_venous
    c_arterial <- arterial / vl_arterial
    c_lung <- lung / vl_lung
    c_adipose <- adipose / vl_adipose
    c_muscle <- muscle / vl_muscle
    c_liver <- liver / vl_liver
    c_spleen <- spleen / vl_spleen
    c_heart <- heart / vl_heart
    c_brain <- brain / vl_brain
    c_kidney <- kidney / vl_kidney
    c_skin <- skin / vl_skin
    c_reproductive_organ <- reproductive_organ / vl_reproductive_organ
    c_red_marrow <- red_marrow / vl_red_marrow
    c_yellow_marrow <- yellow_marrow / vl_yellow_marrow
    c_other <- other / vl_other

    # Plasma concentration of the blood leaving the kidney, the renal clearance
    # driver: arterial plasma.
    cp_arterial <- c_arterial / bp

    # Venous-blood concentration leaving each perfusion-limited tissue.
    cout_lung <- c_lung / kp_lung * bp
    cout_adipose <- c_adipose / kp_adipose * bp
    cout_muscle <- c_muscle / kp_muscle * bp
    cout_liver <- c_liver / kp_liver * bp
    cout_spleen <- c_spleen / kp_spleen * bp
    cout_heart <- c_heart / kp_heart * bp
    cout_brain <- c_brain / kp_brain * bp
    cout_kidney <- c_kidney / kp_kidney * bp
    cout_skin <- c_skin / kp_skin * bp
    cout_reproductive_organ <- c_reproductive_organ / kp_reproductive_organ * bp
    cout_red_marrow <- c_red_marrow / kp_red_marrow * bp
    cout_yellow_marrow <- c_yellow_marrow / kp_yellow_marrow * bp
    cout_other <- c_other / kp_other * bp

    d/dt(venous) <- ql_adipose * cout_adipose + ql_muscle * cout_muscle +
      ql_liver * cout_liver + ql_heart * cout_heart + ql_brain * cout_brain +
      ql_kidney * cout_kidney + ql_skin * cout_skin +
      ql_reproductive_organ * cout_reproductive_organ +
      ql_red_marrow * cout_red_marrow + ql_yellow_marrow * cout_yellow_marrow +
      ql_other * cout_other - ql_lung * c_venous
    d/dt(lung) <- ql_lung * (c_venous - cout_lung)
    d/dt(arterial) <- ql_lung * (cout_lung - c_arterial)
    d/dt(adipose) <- ql_adipose * (c_arterial - cout_adipose)
    d/dt(muscle) <- ql_muscle * (c_arterial - cout_muscle)
    # The ACAT gut is a zero-volume pass-through (Table 1 prints no volume and
    # Table 3 no partition coefficient), so its portal flow enters the liver at
    # arterial concentration.
    d/dt(liver) <- ql_ha * c_arterial + ql_gut * c_arterial +
      ql_spleen * cout_spleen - ql_liver * cout_liver
    d/dt(spleen) <- ql_spleen * (c_arterial - cout_spleen)
    d/dt(heart) <- ql_heart * (c_arterial - cout_heart)
    d/dt(brain) <- ql_brain * (c_arterial - cout_brain)
    d/dt(kidney) <- ql_kidney * (c_arterial - cout_kidney) - cl_renal * cp_arterial
    d/dt(skin) <- ql_skin * (c_arterial - cout_skin)
    d/dt(reproductive_organ) <- ql_reproductive_organ * (c_arterial - cout_reproductive_organ)
    d/dt(red_marrow) <- ql_red_marrow * (c_arterial - cout_red_marrow)
    d/dt(yellow_marrow) <- ql_yellow_marrow * (c_arterial - cout_yellow_marrow)
    d/dt(other) <- ql_other * (c_arterial - cout_other)

    # Venous plasma concentration (GastroPlus 'plasma'). Deterministic: the
    # paper reports no residual-error model.
    Cc <- c_venous / bp
  })
}
