Amice_2021_nicotine_pbpk <- function() {
  description <- "PBPK (whole-body, in-house MATLAB). Nicotine and cotinine disposition in a non-pregnant adult woman after intravenous dosing (Amice 2021, 'Baseline PBPK Model'). Nine flow-limited nicotine compartments (arterial and venous blood, lung, brain, liver, kidney, rapidly perfused, slowly perfused and fat) and nine matching cotinine compartments, coupled in the liver where 80 percent of hepatic nicotine clearance forms cotinine; renal and hepatic clearances act on the unbound tissue-outflow concentration, with urine and metabolised-mass accumulators for each analyte. Deterministic typical-value model: the source reports no between-subject variability and no residual-error model. The paper's pregnancy model (mammary, uterus, placenta and a four-compartment fetus) is NOT implemented, because its placental transfer equations contradict the supplement tables and several of its parameters are unpublished; see the vignette Errata."
  reference <- paste(
    "Amice B, Ho H, Zhang E, Bullen C. (2021).",
    "Physiologically Based Pharmacokinetic Modelling for Nicotine and Cotinine Clearance in Pregnant Women.",
    "Front Pharmacol 12:688597. doi:10.3389/fphar.2021.688597. PMCID PMC8329445.",
    "Equations from Supplementary DataSheet 1 (Eqs 1-22); physiological and drug parameters from",
    "Supplementary DataSheet 2 (Tables S.1 and S.3) and main-text Table 1.",
    "Framework inherited from Robinson DE, Balter NJ, Schwartz SL (1992) J Pharmacokinet Biopharm 20:591-609,",
    "doi:10.1007/BF01064421.",
    sep = " "
  )
  vignette <- "Amice_2021_nicotine_pbpk"

  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL",
    amount = "ug",
    weight = "kg"
  )

  # Nicotine is infused into venous blood (Eq 11 'DoseI.V'); the cotinine
  # infusion validation regimen (Figure 2B) doses cotinine venous blood.
  dosing <- c("a_venous", "a_venous_cot")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every organ volume (fraction of body weight), the cardiac output",
        "(15.126 x WT^0.75 L/h) and the hepatic and renal clearances (mL/min/kg x WT).",
        "Amice 2021 Results 'Baseline PBPK Model' sets the adult woman at 70 kg."
      ),
      source_name = "BW"
    )
  )

  # Every ODE state. Analyte / specimen checked against Amice 2021
  # Supplementary DataSheet 1 Eqs 1-22 and its Table 1 compartment designations.
  compartmentData <- list(
    a_arterial = list(analyte = "nicotine", units = "ug", specimen = "whole blood", verified = TRUE),
    a_fat = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_brain = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_lung = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    a_venous = list(analyte = "nicotine", units = "ug", specimen = "whole blood", verified = TRUE),
    a_urine = list(analyte = "nicotine", units = "ug", specimen = "urine", verified = TRUE),
    a_hepatic = list(analyte = "nicotine", units = "ug", specimen = "not applicable", verified = TRUE),
    a_arterial_cot = list(analyte = "cotinine", units = "ug", specimen = "whole blood", verified = TRUE),
    a_fat_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_brain_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_kidney_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_liver_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_lung_cot = list(analyte = "cotinine", units = "ug", specimen = "tissue", verified = TRUE),
    a_venous_cot = list(analyte = "cotinine", units = "ug", specimen = "whole blood", verified = TRUE),
    a_urine_cot = list(analyte = "cotinine", units = "ug", specimen = "urine", verified = TRUE),
    a_hepatic_cot = list(analyte = "cotinine", units = "ug", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "adult",
    weight_median = "70 kg reference adult woman (Amice 2021 Results 'Baseline PBPK Model', after Robinson 1992)",
    sex_female_pct = 100,
    disease_state = "healthy non-pregnant adult women (non-smokers in the validation datasets)",
    dose_range = paste(
      "intravenous nicotine 10 ug/min for 60 min (Curvall 1990 regimen);",
      "intravenous cotinine 0.67 mg/min for 30 min (De Schepper 1987 regimen);",
      "single intravenous nicotine injection 1 ug/kg (Figure 3A)"
    ),
    notes = paste(
      "A PBPK model, not fitted to individual-level data. The hepatic and renal",
      "clearances for women are the Dempsey 2002 values (Table 1, column 2).",
      "Model outputs were compared with the mean plasma nicotine / cotinine data of",
      "six non-smokers (Curvall 1990) and the plasma cotinine data of four",
      "non-smokers (De Schepper 1987), as in Robinson 1992."
    )
  )

  ini({
    # =====================================================================
    # PHYSIOLOGY -- Supplementary Table S.1 ('Compartment specific
    # physiological parameters for adults'). Volumes are fractions of body
    # weight (L/kg, unit tissue density); flows are fractions of cardiac
    # output. The arterial, venous and lung blood flows are the whole
    # cardiac output (FQ = 1 in Table S.1).
    # =====================================================================
    qcc <- fixed(15.126); label("Cardiac output constant (L/h/kg^0.75)") # Table S.1 'Cardiac output fraction to body weight' FCO = 15.126 L/h/kg0.75

    fvol_arterial <- fixed(0.02); label("Arterial blood volume (fraction of body weight)") # Table S.1 Arterial Ftissue = 0.02
    fvol_fat <- fixed(0.187); label("Fat volume (fraction of body weight)") # Table S.1 Fat Ftissue = 0.187
    fvol_slowly_perfused <- fixed(0.536); label("Slowly perfused tissue volume (fraction of body weight)") # Table S.1 Slow Ftissue = 0.536
    fvol_rapidly_perfused <- fixed(0.032); label("Rapidly perfused tissue volume (fraction of body weight)") # Table S.1 Rapid Ftissue = 0.032
    fvol_brain <- fixed(0.021); label("Brain volume (fraction of body weight)") # Table S.1 Brain Ftissue = 0.021
    fvol_kidney <- fixed(0.004); label("Kidney volume (fraction of body weight)") # Table S.1 Kidney Ftissue = 0.004
    fvol_liver <- fixed(0.026); label("Liver volume (fraction of body weight)") # Table S.1 Liver Ftissue = 0.026
    fvol_lung <- fixed(0.008); label("Lung volume (fraction of body weight)") # Table S.1 Lung Ftissue = 0.008
    fvol_venous <- fixed(0.0571); label("Venous blood volume (fraction of body weight)") # Table S.1 Venous Ftissue = 0.0571

    fq_fat <- fixed(0.052); label("Fat blood flow (fraction of cardiac output)") # Table S.1 Fat FQtissue = 0.052
    fq_slowly_perfused <- fixed(0.188); label("Slowly perfused tissue blood flow (fraction of cardiac output)") # Table S.1 Slow FQtissue = 0.188
    fq_rapidly_perfused <- fixed(0.216); label("Rapidly perfused tissue blood flow (fraction of cardiac output)") # Table S.1 Rapid FQtissue = 0.216
    fq_brain <- fixed(0.117); label("Brain blood flow (fraction of cardiac output)") # Table S.1 Brain FQtissue = 0.117
    fq_kidney <- fixed(0.177); label("Kidney blood flow (fraction of cardiac output)") # Table S.1 Kidney FQtissue = 0.177
    fq_liver <- fixed(0.25); label("Liver blood flow (fraction of cardiac output)") # Table S.1 Liver FQtissue = 0.25

    # =====================================================================
    # NICOTINE -- Table 1 column 'Parameters for women' (Dempsey 2002,
    # mL/min/kg) and Supplementary Table S.3
    # =====================================================================
    lclint <- fixed(log(16.2)); label("Nicotine hepatic clearance, applied to the unbound liver-outflow concentration (log mL/min/kg)") # Table 1 'Nicotine hepatic clearance', women = 16.2 ml/min/kg; Table S.3 Hepatic 16.2
    lcl_renal <- fixed(log(0.7)); label("Nicotine renal clearance, applied to the unbound kidney-outflow concentration (log mL/min/kg)") # Table 1 'Nicotine renal clearance', women = 0.7 ml/min/kg; Table S.3 Renal 0.7
    fu <- fixed(0.95); label("Nicotine unbound fraction (unitless)") # Table S.3 'Unbound fraction' nicotine = 0.95
    bpr <- fixed(0.8); label("Nicotine blood:plasma concentration ratio (unitless)") # Table S.3 'Blood/plasma ratio' nicotine = 0.8
    fm <- fixed(0.8); label("Fraction of hepatic nicotine clearance forming cotinine (unitless)") # Supplementary DataSheet 1 Eq 20 term '0.8 fuN CLN Li CN Li/KN Li'; main text Methods 'approximately 80% of nicotine is metabolized into cotinine'
    lkp_fat <- fixed(log(1)); label("Nicotine fat:blood partition coefficient (log unitless)") # Table S.3 Fat/Blood nicotine = 1
    lkp_slowly_perfused <- fixed(log(2.5)); label("Nicotine slowly perfused:blood partition coefficient (log unitless)") # Table S.3 Slow/Blood nicotine = 2.5
    lkp_rapidly_perfused <- fixed(log(3)); label("Nicotine rapidly perfused:blood partition coefficient (log unitless)") # Table S.3 Rapid/Blood nicotine = 3
    lkp_brain <- fixed(log(3)); label("Nicotine brain:blood partition coefficient (log unitless)") # Table S.3 Brain/Blood nicotine = 3
    lkp_kidney <- fixed(log(15)); label("Nicotine kidney:blood partition coefficient (log unitless)") # Table S.3 Kidney/Blood nicotine = 15
    lkp_liver <- fixed(log(9)); label("Nicotine liver:blood partition coefficient (log unitless)") # Table S.3 Liver/Blood nicotine = 9
    lkp_lung <- fixed(log(2)); label("Nicotine lung:blood partition coefficient (log unitless)") # Table S.3 Lung/Blood nicotine = 2

    # =====================================================================
    # COTININE -- Table 1 column 'Parameters for women' and Table S.3
    # =====================================================================
    lclint_cot <- fixed(log(0.5)); label("Cotinine hepatic clearance, applied to the unbound liver-outflow concentration (log mL/min/kg)") # Table 1 'Cotinine hepatic clearance', women = 0.5 ml/min/kg; Table S.3 Hepatic 0.5
    lcl_renal_cot <- fixed(log(0.2)); label("Cotinine renal clearance, applied to the unbound kidney-outflow concentration (log mL/min/kg)") # Table 1 'Cotinine renal clearance', women = 0.2 ml/min/kg; Table S.3 Renal 0.2
    fu_cot <- fixed(0.974); label("Cotinine unbound fraction (unitless)") # Table S.3 'Unbound fraction' cotinine = 0.974
    bpr_cot <- fixed(0.88); label("Cotinine blood:plasma concentration ratio (unitless)") # Table S.3 'Blood/plasma ratio' cotinine = 0.88
    lkp_fat_cot <- fixed(log(0.5)); label("Cotinine fat:blood partition coefficient (log unitless)") # Table S.3 Fat/Blood cotinine = 0.5
    lkp_slowly_perfused_cot <- fixed(log(0.85)); label("Cotinine slowly perfused:blood partition coefficient (log unitless)") # Table S.3 Slow/Blood cotinine = 0.85
    lkp_rapidly_perfused_cot <- fixed(log(1.5)); label("Cotinine rapidly perfused:blood partition coefficient (log unitless)") # Table S.3 Rapid/Blood cotinine = 1.5
    lkp_brain_cot <- fixed(log(1.5)); label("Cotinine brain:blood partition coefficient (log unitless)") # Table S.3 Brain/Blood cotinine = 1.5
    lkp_kidney_cot <- fixed(log(2)); label("Cotinine kidney:blood partition coefficient (log unitless)") # Table S.3 Kidney/Blood cotinine = 2
    lkp_liver_cot <- fixed(log(2)); label("Cotinine liver:blood partition coefficient (log unitless)") # Table S.3 Liver/Blood cotinine = 2
    lkp_lung_cot <- fixed(log(1)); label("Cotinine lung:blood partition coefficient (log unitless)") # Table S.3 Lung/Blood cotinine = 1
  })

  model({
    # =====================================================================
    # 1. VOLUMES (L), FLOWS (L/h) AND CLEARANCES (L/h)
    #    Table 1 clearances are mL/min/kg: x WT x 60 / 1000 gives L/h.
    # =====================================================================
    qc <- qcc * WT^0.75
    v_arterial <- fvol_arterial * WT
    v_fat <- fvol_fat * WT
    v_slowly_perfused <- fvol_slowly_perfused * WT
    v_rapidly_perfused <- fvol_rapidly_perfused * WT
    v_brain <- fvol_brain * WT
    v_kidney <- fvol_kidney * WT
    v_liver <- fvol_liver * WT
    v_lung <- fvol_lung * WT
    v_venous <- fvol_venous * WT

    q_fat <- fq_fat * qc
    q_slowly_perfused <- fq_slowly_perfused * qc
    q_rapidly_perfused <- fq_rapidly_perfused * qc
    q_brain <- fq_brain * qc
    q_kidney <- fq_kidney * qc
    q_liver <- fq_liver * qc

    clint <- exp(lclint) * WT * 0.06
    cl_renal <- exp(lcl_renal) * WT * 0.06
    clint_cot <- exp(lclint_cot) * WT * 0.06
    cl_renal_cot <- exp(lcl_renal_cot) * WT * 0.06

    kp_fat <- exp(lkp_fat)
    kp_slowly_perfused <- exp(lkp_slowly_perfused)
    kp_rapidly_perfused <- exp(lkp_rapidly_perfused)
    kp_brain <- exp(lkp_brain)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_lung <- exp(lkp_lung)
    kp_fat_cot <- exp(lkp_fat_cot)
    kp_slowly_perfused_cot <- exp(lkp_slowly_perfused_cot)
    kp_rapidly_perfused_cot <- exp(lkp_rapidly_perfused_cot)
    kp_brain_cot <- exp(lkp_brain_cot)
    kp_kidney_cot <- exp(lkp_kidney_cot)
    kp_liver_cot <- exp(lkp_liver_cot)
    kp_lung_cot <- exp(lkp_lung_cot)

    # =====================================================================
    # 2. NICOTINE (Supplementary DataSheet 1 Eqs 1-5 and 8-11, written as
    #    amounts; the mammary and uterus terms of Eqs 6, 7 and 11 belong to
    #    the pregnancy model and are absent in the non-pregnant woman).
    #    Each tissue's venous outflow is at concentration C_tissue / Kp.
    # =====================================================================
    c_arterial <- a_arterial / v_arterial
    c_venous <- a_venous / v_venous
    cv_fat <- a_fat / v_fat / kp_fat
    cv_slowly_perfused <- a_slowly_perfused / v_slowly_perfused / kp_slowly_perfused
    cv_rapidly_perfused <- a_rapidly_perfused / v_rapidly_perfused / kp_rapidly_perfused
    cv_brain <- a_brain / v_brain / kp_brain
    cv_kidney <- a_kidney / v_kidney / kp_kidney
    cv_liver <- a_liver / v_liver / kp_liver
    cv_lung <- a_lung / v_lung / kp_lung

    r_renal <- fu * cl_renal * cv_kidney
    r_hepatic <- fu * clint * cv_liver

    d/dt(a_arterial) <- qc * (cv_lung - c_arterial)
    d/dt(a_fat) <- q_fat * (c_arterial - cv_fat)
    d/dt(a_slowly_perfused) <- q_slowly_perfused * (c_arterial - cv_slowly_perfused)
    d/dt(a_rapidly_perfused) <- q_rapidly_perfused * (c_arterial - cv_rapidly_perfused)
    d/dt(a_brain) <- q_brain * (c_arterial - cv_brain)
    d/dt(a_kidney) <- q_kidney * (c_arterial - cv_kidney) - r_renal
    d/dt(a_liver) <- q_liver * (c_arterial - cv_liver) - r_hepatic
    d/dt(a_lung) <- qc * (c_venous - cv_lung)
    d/dt(a_venous) <- q_fat * cv_fat + q_slowly_perfused * cv_slowly_perfused +
      q_rapidly_perfused * cv_rapidly_perfused + q_brain * cv_brain +
      q_kidney * cv_kidney + q_liver * cv_liver - qc * c_venous
    d/dt(a_urine) <- r_renal
    d/dt(a_hepatic) <- r_hepatic

    # =====================================================================
    # 3. COTININE (Eqs 12-16 and 19-22). Cotinine is formed in the liver at
    #    fm x the hepatic nicotine elimination rate, on a mass basis (Eq 20
    #    carries no nicotine / cotinine molecular-weight ratio).
    # =====================================================================
    c_arterial_cot <- a_arterial_cot / v_arterial
    c_venous_cot <- a_venous_cot / v_venous
    cv_fat_cot <- a_fat_cot / v_fat / kp_fat_cot
    cv_slowly_perfused_cot <- a_slowly_perfused_cot / v_slowly_perfused / kp_slowly_perfused_cot
    cv_rapidly_perfused_cot <- a_rapidly_perfused_cot / v_rapidly_perfused / kp_rapidly_perfused_cot
    cv_brain_cot <- a_brain_cot / v_brain / kp_brain_cot
    cv_kidney_cot <- a_kidney_cot / v_kidney / kp_kidney_cot
    cv_liver_cot <- a_liver_cot / v_liver / kp_liver_cot
    cv_lung_cot <- a_lung_cot / v_lung / kp_lung_cot

    r_renal_cot <- fu_cot * cl_renal_cot * cv_kidney_cot
    r_hepatic_cot <- fu_cot * clint_cot * cv_liver_cot

    d/dt(a_arterial_cot) <- qc * (cv_lung_cot - c_arterial_cot)
    d/dt(a_fat_cot) <- q_fat * (c_arterial_cot - cv_fat_cot)
    d/dt(a_slowly_perfused_cot) <- q_slowly_perfused * (c_arterial_cot - cv_slowly_perfused_cot)
    d/dt(a_rapidly_perfused_cot) <- q_rapidly_perfused * (c_arterial_cot - cv_rapidly_perfused_cot)
    d/dt(a_brain_cot) <- q_brain * (c_arterial_cot - cv_brain_cot)
    d/dt(a_kidney_cot) <- q_kidney * (c_arterial_cot - cv_kidney_cot) - r_renal_cot
    d/dt(a_liver_cot) <- q_liver * (c_arterial_cot - cv_liver_cot) + fm * r_hepatic - r_hepatic_cot
    d/dt(a_lung_cot) <- qc * (c_venous_cot - cv_lung_cot)
    d/dt(a_venous_cot) <- q_fat * cv_fat_cot + q_slowly_perfused * cv_slowly_perfused_cot +
      q_rapidly_perfused * cv_rapidly_perfused_cot + q_brain * cv_brain_cot +
      q_kidney * cv_kidney_cot + q_liver * cv_liver_cot - qc * c_venous_cot
    d/dt(a_urine_cot) <- r_renal_cot
    d/dt(a_hepatic_cot) <- r_hepatic_cot

    # =====================================================================
    # 4. OUTPUTS (ug/L = ng/mL). The equations track blood; plasma is blood
    #    divided by the Table S.3 blood:plasma ratio, which is the reading
    #    that matches the post-infusion plasma curves of Figure 2A. Brain
    #    outputs are tissue concentrations (Figure 3A). No residual-error
    #    model is reported, so none is declared.
    # =====================================================================
    Cc <- c_venous / bpr
    Cc_cot <- c_venous_cot / bpr_cot
    Cblood <- c_venous
    Cblood_cot <- c_venous_cot
    Cbrain <- a_brain / v_brain
    Cbrain_cot <- a_brain_cot / v_brain
  })
}
