Fu_2022_nifedipine <- function() {
  description <- paste(
    "One-compartment population PK model for oral controlled-release",
    "nifedipine in Chinese adults with chronic kidney disease (stages 2-5, not",
    "on dialysis) and hypertension (Fu 2022). Zero-order input into the",
    "central compartment at an estimated constant rate K0 and first-order",
    "elimination; apparent clearance falls with age, hemoglobin and serum",
    "parathyroid hormone (power effects on the cohort medians), apparent volume",
    "rises with ALT; additive residual error."
  )
  reference <- paste(
    "Fu C, Pei Q, Liang W, Yang B, Li W, Liu J, Tan H, Guo C, Zhang H, Yang G.",
    "Population Pharmacokinetic Modelling for Nifedipine to Evaluate the Effect",
    "of Parathyroid Hormone on CYP3A in Patients with Chronic Kidney Disease.",
    "Drug Des Devel Ther. 2022;16:2261-2274. doi:10.2147/DDDT.S362607.",
    sep = " "
  )
  vignette <- "Fu_2022_nifedipine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "nifedipine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE/54)^-0.22 on CL/F; 54 years is the cohort median (Table 1, range 24-84).",
      source_name = "Age"
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (HGB/94)^-0.35 on CL/F; 94 g/L is the cohort median (Table 1, range 40-177 g/L). Values in g/dL must be multiplied by 10.",
      source_name = "Hb"
    ),
    PARATHYROID_HORMONE = list(
      description = "Serum parathyroid hormone concentration",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (PARATHYROID_HORMONE/156.84)^-0.22 on CL/F; 156.84 pg/mL is the cohort median (Table 1, range 9.34-687.27). The paper does not state whether the assay measured intact PTH.",
      source_name = "PTH"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (ALT/14)^0.48 on V/F; 14 U/L is the cohort median (Table 1, range 1-81).",
      source_name = "ALT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 157,
    n_studies = 1,
    n_observations = 612,
    age_range = "24-84 years",
    age_median = "54 years",
    bmi_range = "15.81-36.70 kg/m^2",
    bmi_median = "23.4 kg/m^2",
    sex_female_pct = 59.9,
    race_ethnicity = "Chinese (100%)",
    disease_state = paste(
      "Chronic kidney disease with hypertension, not on dialysis; KDIGO stage",
      "2/3/4/5 in 11/20/62/64 patients; median GFR 16.85 mL/min/1.73 m^2",
      "(4.84-73.93)."
    ),
    dose_range = paste(
      "Nifedipine controlled-release 30 mg tablets (Bayer) for 7 days; the",
      "daily dose per patient is not reported (the paper's simulations use",
      "30 mg once daily)."
    ),
    regions = "China (Third Xiangya Hospital, Changsha)",
    renal_function = "Median serum creatinine 478 umol/L (35-2017); median cystatin C 3.12 mg/L (0.66-8.22).",
    parathyroid_hormone = "Median 156.84 pg/mL (9.34-687.27).",
    notes = paste(
      "Prospective single-centre study (NCT03695107), August 2017 to April 2019.",
      "Sparse sampling on day 7 at predose, 3-4 h, 6-8 h and 24 h post-dose",
      "(2-4 samples per patient). Demographics from Table 1 (male:female",
      "63:94); body weight was not tabulated."
    )
  )

  ini({
    lcl <- log(49.61)
    label("Apparent clearance CL/F at AGE 54 y, HGB 94 g/L, PTH 156.84 pg/mL (L/h)") # Table 3 tvCL = 49.61 L/h
    lvc <- log(2300.26)
    label("Apparent volume of distribution V/F at ALT 14 U/L (L)") # Table 3 tvV = 2300.26 L
    lr1 <- log(4.21)
    label("Zero-order input rate K0 of the oral dose into central (mg/h)") # Table 3 tvK0 = 4.21 mg/h

    e_age_cl <- -0.22
    label("Power exponent of age on CL/F (unitless)") # Table 3 dCLdAge = -0.22
    e_hgb_cl <- -0.35
    label("Power exponent of hemoglobin on CL/F (unitless)") # Table 3 dCLdHb = -0.35
    e_parathyroid_hormone_cl <- -0.22
    label("Power exponent of serum parathyroid hormone on CL/F (unitless)") # Table 3 dCLdPTH = -0.22
    e_alt_vc <- 0.48
    label("Power exponent of ALT on V/F (unitless)") # Table 3 dVdALT = 0.48

    # IIV (%) of Table 3 read as a log-normal CV: omega^2 = log(CV^2 + 1)
    etalcl ~ 0.29276 # Table 3 IIV CL 58.32%
    etalvc ~ 0.18925 # Table 3 IIV V 45.644%
    etalr1 ~ 0.27107 # Table 3 IIV K0 55.80%

    addSd <- 8.12
    label("Additive residual error (ng/mL)") # Table 3 stdev0 = 8.12
  })
  model({
    # Final-model equations, p.4: power effects on the Table 1 medians
    cl <- exp(lcl + etalcl) *
      (AGE / 54)^e_age_cl *
      (HGB / 94)^e_hgb_cl *
      (PARATHYROID_HORMONE / 156.84)^e_parathyroid_hormone_cl
    vc <- exp(lvc + etalvc) * (ALT / 14)^e_alt_vc
    r1 <- exp(lr1 + etalr1)

    kel <- cl / vc

    d / dt(central) <- -kel * central
    # Zero-order input at rate K0 (duration = dose / K0); dose records
    # need rate = -1 so rxode2 reads the modelled rate.
    rate(central) <- r1

    # Dose in mg and volume in L give mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
