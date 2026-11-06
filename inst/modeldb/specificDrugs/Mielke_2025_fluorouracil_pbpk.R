Mielke_2025_fluorouracil_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited, eight tissue compartments; rxode2 script deposited",
    "as Supplementary Material S1). Physiologically based pharmacokinetic model of",
    "intravenous 5-fluorouracil (5-FU) for a 73 kg ICRP reference adult (Mielke 2025),",
    "used to compare Cmax, the infusion steady-state concentration and AUC across",
    "infusion rates. Tissues: fat, brain, kidney, lumped muscle/skin, vessel-rich tissue",
    "(excluding red marrow), red marrow, skeleton (excluding marrow) and liver, each",
    "perfusion-limited with a Schmitt-algorithm tissue:blood partition coefficient.",
    "Blood is not a compartment: mixed venous and arterial blood concentrations are",
    "algebraic flow-weighted averages, and the infused drug enters arterial blood in a",
    "hypothetical 1 L/h infusate stream. Elimination is Michaelis-Menten metabolism of",
    "the liver tissue concentration only (Vmax and Km are the means of three published",
    "popPK estimates). The infusion is driven, as in the authors' script, by the",
    "`infusion_rate` state, which holds the current infusion rate in mg/h: dose",
    "amt = rate (mg/h) into `infusion_rate` at the start of an infusion and",
    "amt = -rate at its end. Deterministic: no IIV, no residual error. Cc is the",
    "arterial blood concentration (mg/L)."
  )
  reference <- paste(
    "Mielke H, Algharably EAE, Gundert-Remy U, Enioutina EY (2025).",
    "5-Fluorouracil Toxicity: Revisiting the Relevance of Pharmacokinetic Parameters.",
    "Pharmaceuticals (Basel) 18(5):653. doi:10.3390/ph18050653. PMCID PMC12114802.",
    "Model code: Supplementary Material S1, 'R script for 5-FU PBPK model code'.",
    "Vmax and Km are the means of Terret et al. (2000) Clin Pharmacol Ther 68:270-279,",
    "van Kuilenburg et al. (2012) Clin Pharmacokinet 51:163-174 and Woloch et al.",
    "(2012) Curr Top Med Chem 12:1713-1719.",
    sep = " "
  )
  vignette <- "Mielke_2025_fluorouracil_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Lumped muscle/skin, red marrow and skeleton have no registered amount-state
  # name; infusion_rate is the authors' rate-holder state; the two auc_* states
  # integrate the arterial and red-marrow concentrations (script AUC, AUCRM).
  paper_specific_compartments <- c(
    "infusion_rate",
    "a_muscle_skin",
    "a_red_marrow",
    "a_bone",
    "auc_arterial",
    "auc_red_marrow"
  )

  compartmentData <- list(
    infusion_rate = list(analyte = "fluorouracil", units = "mg/h", specimen = "administration site", verified = TRUE),
    a_fat = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle_skin = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_red_marrow = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_bone = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "fluorouracil", units = "mg", specimen = "tissue", verified = TRUE),
    a_metabolized = list(analyte = "fluorouracil", units = "mg", specimen = "not applicable", verified = TRUE),
    auc_arterial = list(analyte = "fluorouracil", units = "mg*h/L", specimen = "not applicable", verified = TRUE),
    auc_red_marrow = list(analyte = "fluorouracil", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "adult (single ICRP reference adult; no age term in the model)",
    weight_range = "73 kg (Table 1 'Body weight (bw) (kg) 73'); the script carries absolute volumes and flows, so body weight does not enter the model",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "adults receiving intravenous 5-FU for solid tumours (simulation scenarios); parameters are literature constants, not a fit",
    dose_range = paste(
      "Simulated: 3052 mg over 2 h (the patient's erroneous infusion) and 4612 mg over 24 h",
      "(the planned infusion); 650 and 1000 mg/m^2/day as 24 h infusions (evaluation against",
      "Arshad 2020); 1.0 g IV bolus (evaluation against Finch 1979); and the 138 individual",
      "dosing records of Ma et al. 2017 (Cancer).",
      sep = " "
    ),
    regions = "Germany (Federal Institute for Risk Assessment) and United States (University of Utah)",
    notes = paste(
      "A forward-simulation PBPK model, not fitted to data, so n_subjects = 0. Physiology",
      "from ICRP Publication 89 (Table 1); partition coefficients from the Schmitt",
      "algorithm; hepatic Michaelis-Menten metabolism only (the paper attributes 85% of",
      "5-FU elimination to hepatic DPD); DPD-deficient patients are not represented. The",
      "model was evaluated against digitised concentration-time data from Arshad 2020",
      "(continuous infusion) and Finch 1979 (1.0 g bolus), Figure 1.",
      sep = " "
    )
  )

  ini({
    # Every value is a literature or physiological constant from the deposited
    # script (Supplementary Material S1); none was estimated from data.

    # Tissue volumes (L), script 'Tissue volumes'
    v_fat <- fixed(18.2); label("Fat tissue volume (L)") # S1 VF = 18.200; Table 1 Fat 18.2
    v_liver <- fixed(1.8); label("Liver volume (L)") # S1 VL = 1.800; Table 1 Liver 1.8
    v_brain <- fixed(1.45); label("Brain volume (L)") # S1 VBr = 1.450; Table 1 Brain 1.45
    v_kidney <- fixed(0.31); label("Kidney volume (L)") # S1 VK = 0.310; Table 1 Kidney 0.31
    v_muscle_skin <- fixed(32.3); label("Lumped muscle and skin volume (L)") # S1 VMS = 32.300 (Table 1 lists muscle 0.40 x bw and skin 0.0371 x bw separately)
    v_rapidly_perfused <- fixed(2.598); label("Vessel-rich tissue volume excluding red marrow (L)") # S1 VVR = 3.768 - 1.170; Table 1 2.6
    v_red_marrow <- fixed(1.17); label("Red marrow volume (L)") # S1 VRM = 1.170; Table 1 1.2
    v_bone <- fixed(9.33); label("Skeleton volume excluding marrow (L)") # S1 VSkel = 9.330 (Table 1 prints 6.9)

    # Blood flows (L/h), script 'Blood flows'
    q_fat <- fixed(19.5); label("Fat blood flow (L/h)") # S1 QF = 19.5; Table 1 Fat 19.5
    q_liver <- fixed(99.5); label("Liver blood flow (L/h)") # S1 QL = 99.5; Table 1 Liver 99.5
    q_brain <- fixed(46.8); label("Brain blood flow (L/h)") # S1 QBr = 46.8; Table 1 Brain 46.8
    q_kidney <- fixed(74.1); label("Kidney blood flow (L/h)") # S1 QK = 74.1; Table 1 Kidney 74.1
    q_muscle_skin <- fixed(85.8); label("Lumped muscle and skin blood flow (L/h)") # S1 QMS = 85.8; Table 1 muscle 65.8 + skin 20
    q_rapidly_perfused <- fixed(44.8); label("Vessel-rich tissue blood flow excluding red marrow (L/h)") # S1 QVR = 56.5 - 11.7; Table 1 44.8
    q_red_marrow <- fixed(11.7); label("Red marrow blood flow (L/h)") # S1 QRM = 11.7; Table 1 11.7
    q_bone <- fixed(7.8); label("Skeleton blood flow (L/h)") # S1 QSkel = 7.8; Table 1 7.8
    q_infusion <- fixed(1); label("Hypothetical infusate flow into arterial blood (L/h)") # S1 Qinf = 1, 'hypothetic infusion 1 l/h'

    # Tissue:blood partition coefficients, script 'according to Schmitt'
    lkp_fat <- fixed(log(0.2)); label("Fat:blood partition coefficient (unitless)") # S1 PF = 0.2; Table 1 0.20
    lkp_liver <- fixed(log(0.92)); label("Liver:blood partition coefficient (unitless)") # S1 PL = 0.92; Table 1 0.92
    lkp_brain <- fixed(log(1.35)); label("Brain:blood partition coefficient (unitless)") # S1 PBr = 1.35; Table 1 1.35
    lkp_kidney <- fixed(log(0.98)); label("Kidney:blood partition coefficient (unitless)") # S1 PK = 0.98; Table 1 0.98
    lkp_muscle_skin <- fixed(log(2)); label("Muscle/skin:blood partition coefficient (unitless)") # S1 PMS = 2 (Table 1 muscle 2.19, skin 10.04; S1 comment 'Schmitt: muscle 2.19, skin 1.04')
    lkp_rapidly_perfused <- fixed(log(0.92)); label("Vessel-rich tissue:blood partition coefficient (unitless)") # S1 PVR = 0.92, set equal to liver (Table 1 prints 0.98)
    lkp_red_marrow <- fixed(log(0.9)); label("Red marrow:blood partition coefficient (unitless)") # S1 PRM = 0.9, from the Schmitt red-blood-cell value (Table 1 prints 0.92)
    lkp_bone <- fixed(log(0.63)); label("Skeleton:blood partition coefficient (unitless)") # S1 PSkel = 0.63; Table 1 0.63

    # Hepatic Michaelis-Menten metabolism: means of three popPK estimates
    # (Terret 2000: 1390 mg/h, 5.57 mg/L; van Kuilenburg 2012: 1512, 4.55;
    # Woloch 2012: 763, 25).
    lvmax <- fixed(log(1221.667)); label("Maximum rate of hepatic metabolism (mg/h)") # S1 VMAX = mean(c(1390, 1512, 763)); Table 1 1221.7
    lkm <- fixed(log(11.70667)); label("Michaelis-Menten constant on the liver tissue concentration (mg/L)") # S1 KM = mean(c(5.57, 4.55, 25)); Table 1 11.7 (Section 4.3.1 text prints 11.1)
  })

  model({
    kp_fat <- exp(lkp_fat)
    kp_liver <- exp(lkp_liver)
    kp_brain <- exp(lkp_brain)
    kp_kidney <- exp(lkp_kidney)
    kp_muscle_skin <- exp(lkp_muscle_skin)
    kp_rapidly_perfused <- exp(lkp_rapidly_perfused)
    kp_red_marrow <- exp(lkp_red_marrow)
    kp_bone <- exp(lkp_bone)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # Cardiac output is the sum of the tissue flows (S1 QC = 390 L/h; Table 1
    # prints 401.7).
    q_cardiac <- q_fat + q_brain + q_kidney + q_muscle_skin + q_rapidly_perfused +
      q_red_marrow + q_bone + q_liver

    # Tissue concentrations (mg/L) and the venous concentrations leaving each
    # tissue, C_VT = C_T / P_T (paper Eq. 4).
    c_fat <- a_fat / v_fat
    c_brain <- a_brain / v_brain
    c_kidney <- a_kidney / v_kidney
    c_muscle_skin <- a_muscle_skin / v_muscle_skin
    c_rapidly_perfused <- a_rapidly_perfused / v_rapidly_perfused
    c_red_marrow <- a_red_marrow / v_red_marrow
    c_bone <- a_bone / v_bone
    c_liver <- a_liver / v_liver

    cv_fat <- c_fat / kp_fat
    cv_brain <- c_brain / kp_brain
    cv_kidney <- c_kidney / kp_kidney
    cv_muscle_skin <- c_muscle_skin / kp_muscle_skin
    cv_rapidly_perfused <- c_rapidly_perfused / kp_rapidly_perfused
    cv_red_marrow <- c_red_marrow / kp_red_marrow
    cv_bone <- c_bone / kp_bone
    cv_liver <- c_liver / kp_liver

    # Mixed venous blood is the flow-weighted mean of the tissue effluents (S1 CV).
    c_venous <- (q_fat * cv_fat + q_brain * cv_brain + q_kidney * cv_kidney +
      q_muscle_skin * cv_muscle_skin + q_rapidly_perfused * cv_rapidly_perfused +
      q_red_marrow * cv_red_marrow + q_bone * cv_bone + q_liver * cv_liver) / q_cardiac

    # Arterial blood mixes the venous return with the infusate stream (S1 CA).
    # The script writes the infusate term as Qinf * iv.dose with iv.dose holding
    # the rate in mg/h and Qinf = 1 L/h; it is written here as the mass input
    # rate (mg/h), which is identical at the published q_infusion = 1 L/h.
    c_arterial <- (q_cardiac * c_venous + infusion_rate) / (q_cardiac + q_infusion)

    # Michaelis-Menten metabolism of the liver tissue concentration (S1 RAM).
    r_metab <- vmax * c_liver / (km + c_liver)

    # Holds the running infusion rate (mg/h); changed only by dose records
    # (S1 d/dt(iv.dose) = 0).
    d/dt(infusion_rate) <- 0

    # Perfusion-limited tissues (paper Eqs. 2-3; S1 ODEs).
    d/dt(a_fat) <- q_fat * (c_arterial - cv_fat)
    d/dt(a_brain) <- q_brain * (c_arterial - cv_brain)
    d/dt(a_kidney) <- q_kidney * (c_arterial - cv_kidney)
    d/dt(a_muscle_skin) <- q_muscle_skin * (c_arterial - cv_muscle_skin)
    d/dt(a_rapidly_perfused) <- q_rapidly_perfused * (c_arterial - cv_rapidly_perfused)
    d/dt(a_red_marrow) <- q_red_marrow * (c_arterial - cv_red_marrow)
    d/dt(a_bone) <- q_bone * (c_arterial - cv_bone)
    d/dt(a_liver) <- q_liver * (c_arterial - cv_liver) - r_metab
    d/dt(a_metabolized) <- r_metab

    # Cumulative AUCs of the arterial and red-marrow concentrations (S1 AUC, AUCRM).
    d/dt(auc_arterial) <- c_arterial
    d/dt(auc_red_marrow) <- c_red_marrow

    Cc <- c_arterial
    Cvenous <- c_venous
    Cred_marrow <- c_red_marrow
  })
}
