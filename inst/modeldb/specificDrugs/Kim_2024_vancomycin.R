Kim_2024_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in healthy Korean adult volunteers (Kim 2024). Clearance scales as power functions of body weight (exponent 1.28, reference 69.45 kg) and age (exponent -0.377, reference 38.5 years); intercompartmental clearance scales as a power function of body weight (exponent 2.38) and peripheral volume as a power function of body surface area (exponent 3.09, reference 1.80 m^2). IIV on CL and on the central volume (the latter fixed); proportional residual error."
  reference <- "Kim YK, Kim D, Kang G, Zang DY, Lee DH.Pharmacokinetics of Vancomycin in Healthy Korean Volunteers and Monte Carlo Simulations to Explore Optimal Dosage Regimens in Patients with Normal Renal Function. Antibiotics (Basel). 2024;13(10):993. doi:10.3390/antibiotics13100993"
  vignette <- "Kim_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin is given as a 2-h IV infusion, so the dose enters `central`
  # directly and there is no depot state. Kim 2024 Section 4.2-4.3 states that
  # venous blood was drawn into heparinised tubes and vancomycin was assayed in
  # PLASMA by LC-MS/MS. The paper never states what matrix the peripheral
  # distribution compartment represents; "plasma" follows the repository
  # default for a mathematical peripheral compartment.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (single-dose study). Reference weight 69.45 kg is the normalising constant printed in the Kim 2024 Table 2 CL and Q equations (Table 1 median 69.5 kg, IQR 60.7-78.2; range 53.2-89.4 kg). Enters CL with exponent 1.28 and Q with exponent 2.38.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Reference age 38.5 years is the normalising constant printed in the Kim 2024 Table 2 CL equation and equals the Table 1 median (IQR 31.5-43.3; range 19-49 years). Enters CL with exponent -0.377.",
      source_name = "AGE"
    ),
    BSA = list(
      description = "Body surface area by the DuBois formula, BSA (m^2) = 0.007184 * height(cm)^0.725 * weight(kg)^0.425",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Kim 2024 Section 2.1 states the DuBois formula (printed coefficient 0.007185). Reference 1.80 m^2 is the normalising constant printed in the Kim 2024 Table 2 V2 equation and equals the Table 1 median (IQR 1.70-1.92). Enters the peripheral volume with exponent 3.09.",
      source_name = "BSA"
    )
  )

  # Screened in the Kim 2024 stepwise covariate search (Section 4.4; Table S2)
  # but not retained in the final model. Kept here for provenance only; none is
  # referenced in model().
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL, Q and V2 in forward steps 1-3 of Table S2 but displaced by weight (CL, Q) and BSA (V2). Not retained.",
      source_name = "BMI"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on CL and V2 in forward step 1 of Table S2 but displaced by weight and BSA. Not retained.",
      source_name = "Height"
    ),
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault and eGFR by MDRD and CKD-EPI (creatinine, cystatin C, and combined) equations, each de-normalised to mL/min by BSA",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Each renal-function estimate was significant on CL in forward steps 1-3 of Table S2 (dOFV 6.9-9.97) but none beat weight and age; the paper attributes this to the narrow range of renal function in healthy volunteers (Discussion). Not retained.",
      source_name = "CLCR_CG, eGFRMDRD2, eGFRCE_CR2, eGFRCE_CYS2, eGFRCE_CRCYS2"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Gender was among the demographic covariates considered (Section 4.4); no significant effect is listed in Table S2. Not retained.",
      source_name = "gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    n_observations = 60L,
    age_range = "19-49 years (Table 1 median 38.5, IQR 31.5-43.3)",
    weight_range = "53.2-89.4 kg (Table 1 median 69.5, IQR 60.7-78.2)",
    height_range = "159-180 cm",
    bsa_median = "1.80 m^2 (IQR 1.70-1.92)",
    sex_female_pct = 25,
    race_ethnicity = c(Korean = 100),
    disease_state = "Healthy adult volunteers with no congenital or chronic disease and normal renal function (Cockcroft-Gault CLCR median 113 mL/min, IQR 102-124; Table 1).",
    dose_range = "Single 1 g IV dose infused over 2 h in 200 mL normal saline",
    regions = "Republic of Korea (Hallym University Sacred Heart Hospital, Anyang)",
    renal_function = "CLCR (Cockcroft-Gault) mean 117 mL/min (CV 29.1%); eGFR CKD-EPI creatinine median 110 mL/min/1.73 m^2",
    notes = "Prospective single-centre study conducted August 2021 to July 2022 (Kim 2024 Section 4.1). Three females and nine males. Plasma was sampled at 2 (end of infusion), 3, 4, 6 and 8 h after the start of infusion, giving 60 concentrations. Model fit in NONMEM 7.5 with FOCE-I; covariates selected by stepwise forward inclusion (p < 0.01) and backward elimination (p < 0.001)."
  )

  ini({
    # Structural parameters (Kim 2024 Table 2). The reference subject weighs
    # 69.45 kg, is 38.5 years old and has a BSA of 1.80 m^2.
    lcl <- log(8.83); label("Clearance at WT=69.45 kg and AGE=38.5 years (L/h)") # Kim 2024 Table 2: theta1 = 8.83 L/h (RSE 3.16%)
    lvc <- log(12.9); label("Central volume of distribution V1 (L)") # Kim 2024 Table 2: theta4 = 12.9 L (RSE 26.4%)
    lq <- log(20.2); label("Intercompartmental clearance at WT=69.45 kg (L/h)") # Kim 2024 Table 2: theta5 = 20.2 L/h (RSE 10.0%)
    lvp <- log(27.4); label("Peripheral volume V2 at BSA=1.80 m^2 (L)") # Kim 2024 Table 2: theta7 = 27.4 L (RSE 4.17%)

    # Covariate exponents (Kim 2024 Table 2), all estimated.
    e_wt_cl <- 1.28; label("Power exponent on (WT/69.45) for CL (unitless)") # Kim 2024 Table 2: theta2 = 1.28 (RSE 8.13%)
    e_age_cl <- -0.377; label("Power exponent on (AGE/38.5) for CL (unitless)") # Kim 2024 Table 2: theta3 = -0.377 (RSE 19.2%)
    e_wt_q <- 2.38; label("Power exponent on (WT/69.45) for Q (unitless)") # Kim 2024 Table 2: theta6 = 2.38 (RSE 17.6%)
    e_bsa_vp <- 3.09; label("Power exponent on (BSA/1.80) for V2 (unitless)") # Kim 2024 Table 2: theta8 = 3.09 (RSE 11.7%)

    # Interindividual variability. Section 4.4 specifies the exponential model
    # theta_i = theta * exp(eta_i); Table 2 reports the IIV as %CV, converted
    # here with omega^2 = log(CV^2 + 1). IIV on Q and V2 shrank to nearly 0
    # after the covariates entered (Section 2.2) and is absent from Table 2.
    etalcl ~ 0.006348 # log(0.0798^2 + 1); Kim 2024 Table 2: CL 7.98% (RSE 18.3%, shrinkage 9.92%)
    etalvc ~ fixed(0.202847) # log(0.474^2 + 1); Kim 2024 Table 2: V1 47.4%, footnote f (shrinkage 16.6%)

    # Residual variability (Kim 2024 Table 2).
    propSd <- 0.0694; label("Proportional residual error (fraction)") # Kim 2024 Table 2: 6.94% (RSE 9.21%, shrinkage 13.3%)
  })
  model({
    # Individual PK parameters (Kim 2024 Table 2):
    #   CL = theta1 * (WT/69.45)^theta2 * (AGE/38.5)^theta3
    #   V1 = theta4
    #   Q  = theta5 * (WT/69.45)^theta6
    #   V2 = theta7 * (BSA/1.80)^theta8
    cl <- exp(lcl + etalcl) * (WT / 69.45)^e_wt_cl * (AGE / 38.5)^e_age_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq) * (WT / 69.45)^e_wt_q
    vp <- exp(lvp) * (BSA / 1.80)^e_bsa_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L, so central/vc is in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
