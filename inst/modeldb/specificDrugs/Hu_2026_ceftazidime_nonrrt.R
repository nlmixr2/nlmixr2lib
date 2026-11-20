Hu_2026_ceftazidime_nonrrt <- function() {
  description <- "One-compartment IV population PK model for the ceftazidime component of ceftazidime-avibactam in critically ill Chinese adults with carbapenem-resistant Gram-negative infection who were NOT receiving renal replacement therapy (Hu 2026), with a power-form Cockcroft-Gault creatinine-clearance effect on clearance."
  reference <- paste(
    "Hu J, Wang CH, Wang XJ, Yu YW, Ji JR, Zhu W, Shen P, Chen YB, Zheng X,",
    "Xiao TT, Cai HL, Dong YL, Xiao YH. Model-Informed Ceftazidime-Avibactam",
    "Dosing Stratified by Renal Function and Renal Replacement Therapy Status",
    "in Critically Ill Patients: A Prospective Population Pharmacokinetic and",
    "Monte Carlo Simulation Study. Drug Des Devel Ther. 2026;20:632562.",
    "doi:10.2147/DDDT.S632562.",
    sep = " "
  )
  vignette <- "Hu_2026_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance calculated with the Cockcroft-Gault equation (Hu 2026 Table 1 note); the body-weight convention inside the equation is not stated",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw Cockcroft-Gault mL/min, NOT BSA-normalized; stored under the canonical CRCL column, which accepts raw mL/min when the source applies no BSA normalization (same convention as Chen_2025_ceftazidime.R and Li_2019_ceftazidime.R). Reference value 100.45 mL/min as printed in the covariate equation (Hu 2026 Results); the Table 1 non-RRT median is 100.48 mL/min (IQR 43.89-122.32). Cohort composition: 2 subjects with augmented renal clearance (> 130 mL/min), 9 with 81-130 mL/min and 4 with 31-50 mL/min. The authors' Monte Carlo simulations extrapolated the model to 6-160 mL/min.",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 59 kg (IQR 55.0-63.0), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 67 years (IQR 62.0-76.0), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "2 of 15 non-RRT subjects female, Hu 2026 Table 1 (13 male). Screened; not retained.",
      source_name = "Sex"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 20.5 (IQR 18.40-26.0), Hu 2026 Table 1. Screened; not retained.",
      source_name = "BMI"
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 4 (IQR 2.0-5.0), Hu 2026 Table 1. Screened; not retained.",
      source_name = "SOFA"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 15 (IQR 10.0-19.0), Hu 2026 Table 1. Screened; not retained.",
      source_name = "APACHE II"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 33.70 g/L (IQR 28.70-38.40), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Albumin"
    ),
    LACT = list(
      description = "Serum lactate",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Non-RRT median 1.10 mmol/L (IQR 0.65-1.70), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Lactate"
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not ventilated",
      notes = "10 of 15 non-RRT subjects (66.67%) ventilated, Hu 2026 Table 1. Screened; not retained.",
      source_name = "Mechanical ventilation"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    n_observations = 91L,
    age_median = "67 years (IQR 62.0-76.0)",
    weight_median = "59 kg (IQR 55.0-63.0)",
    sex_female_pct = 13.3,
    race_ethnicity = "Not reported by category; single-centre Chinese cohort (First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    disease_state = "Critically ill ICU adults with culture-confirmed carbapenem-resistant organism infection susceptible to ceftazidime-avibactam (carbapenem-resistant Klebsiella pneumoniae 12/15, carbapenem-resistant Pseudomonas aeruginosa 3/15), NOT receiving renal replacement therapy. Infection sites: lungs 11, bloodstream 1, lungs + bloodstream 2, lungs + urinary tract 1.",
    dose_range = "Ceftazidime-avibactam 2.5 g (2000/500 mg) q8h in 13 subjects and 1.25 g (1000/250 mg) q8h in 2 subjects, each as a 2-hour IV infusion.",
    regions = "China (single centre, Hangzhou)",
    renal_function = "Cockcroft-Gault creatinine clearance median 100.48 mL/min (IQR 43.89-122.32): augmented renal clearance (> 130 mL/min) 2, 81-130 mL/min 9, 31-50 mL/min 4 (Hu 2026 Table 1).",
    notes = "Prospective single-centre cohort, June 2021 to January 2024 (NCT05413343). Arterial plasma sampled at 0, 1, 2, 4, 6 and 8 h after the first and a steady-state infusion; first-dose and steady-state data fitted jointly. Assay UHPLC-Q-TOF, calibration 10-640 mg/L for ceftazidime; no BLQ samples. NONMEM 7.2, FOCE-I. Separate models were fitted for each analyte and for the RRT and non-RRT strata (sibling files Hu_2026_ceftazidime_rrt, Hu_2026_avibactam_nonrrt, Hu_2026_avibactam_rrt). Also screened and not retained: cumulative fluid balance, vasopressor use and infection site."
  )

  ini({
    lcl <- log(3.36); label("Clearance at CRCL 100.45 mL/min (L/h)") # Hu 2026 Table 2, non-RRT CAZ CLPOP = 3.36 L/h (RSE 7%; bootstrap median 3.33, 95% CI 2.82-3.99)
    lvc <- log(20.5); label("Volume of distribution (L)")             # Hu 2026 Table 2, non-RRT CAZ VPOP = 20.5 L (RSE 20%; bootstrap median 20.71, 95% CI 13.81-30.26)

    e_crcl_cl <- 0.792; label("Exponent of CRCL/100.45 on CL (unitless)") # Hu 2026 Table 2, non-RRT CAZ alpha = 0.792 (RSE 17%); also Results CL = CLPOP x (CrCL/100.45)^0.792

    # Table 2 IIV entries are omega^2 VARIANCES (Methods Eq. 1 defines eta with
    # variance omega^2); sqrt(0.368) = 61% matches the 62% CV of the post-hoc
    # Vd in Table 3 (23.77 +/- 14.70 L).
    etalcl ~ 0.0415 # Hu 2026 Table 2, non-RRT CAZ IIV CL = 0.0415 (RSE 40%)
    etalvc ~ 0.368  # Hu 2026 Table 2, non-RRT CAZ IIV V = 0.368 (RSE 35%)

    # Table 2 RV is the proportional-error VARIANCE sigma^2 (Methods Eq. 3 defines
    # epsilon with variance sigma^2); propSd is its square root.
    propSd <- sqrt(0.0559); label("Proportional residual error (fraction)") # Hu 2026 Table 2, non-RRT CAZ RV (proportional) = 0.0559 (RSE 21%)
  })
  model({
    cl <- exp(lcl + etalcl) * (CRCL / 100.45)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
