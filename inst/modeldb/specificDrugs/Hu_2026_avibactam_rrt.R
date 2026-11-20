Hu_2026_avibactam_rrt <- function() {
  description <- "One-compartment IV population PK model for the avibactam component of ceftazidime-avibactam in critically ill Chinese adults with carbapenem-resistant Gram-negative infection receiving continuous renal replacement therapy (CVVH or CVVHDF) (Hu 2026); no covariates retained."
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
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 68.62 mL/min (IQR 48.46-103.52), Hu 2026 Table 1. Screened; not retained in the RRT model (it is retained in the sibling non-RRT model).",
      source_name = "CrCL"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 64 kg (IQR 57.25-77.75), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 72 years (IQR 59.50-75.25), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "5 of 18 RRT subjects female, Hu 2026 Table 1 (13 male). Screened; not retained.",
      source_name = "Sex"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 22.6 (IQR 20.05-28.60), Hu 2026 Table 1. Screened; not retained.",
      source_name = "BMI"
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 8.5 (IQR 6.75-12.50), Hu 2026 Table 1. Screened; not retained.",
      source_name = "SOFA"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 17.5 (IQR 13.50-22.50), Hu 2026 Table 1. Screened; not retained.",
      source_name = "APACHE II"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 34.10 g/L (IQR 30.7-37.5), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Albumin"
    ),
    LACT = list(
      description = "Serum lactate",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "RRT median 1.00 mmol/L (IQR 0.83-1.77), Hu 2026 Table 1. Screened; not retained.",
      source_name = "Lactate"
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not ventilated",
      notes = "16 of 18 RRT subjects (88.89%) ventilated, Hu 2026 Table 1. Screened; not retained.",
      source_name = "Mechanical ventilation"
    ),
    RRT_CVVHDF_STATUS = list(
      description = "Continuous venovenous hemodiafiltration (vs hemofiltration) modality indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = CVVH",
      notes = "CVVHDF in 14 of 18 RRT subjects (Hu 2026 Discussion); per-patient modes in Supplementary Tables 1-2. Screened as 'RRT mode'; not retained.",
      source_name = "RRT mode"
    ),
    URINE_VOL_24H = list(
      description = "Daily residual urine output",
      units = "mL/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Recorded per patient in Supplementary Tables 1-2. Screened; not retained. Residual renal clearance was not quantified and urine drug concentrations were not measured.",
      source_name = "Daily residual urine volume"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    n_observations = 120L,
    age_median = "72 years (IQR 59.50-75.25)",
    weight_median = "64 kg (IQR 57.25-77.75)",
    sex_female_pct = 27.8,
    race_ethnicity = "Not reported by category; single-centre Chinese cohort (First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    disease_state = "Critically ill ICU adults with culture-confirmed carbapenem-resistant organism infection susceptible to ceftazidime-avibactam (carbapenem-resistant Klebsiella pneumoniae 9/18, carbapenem-resistant Pseudomonas aeruginosa 9/18), receiving continuous renal replacement therapy. Infection sites: lungs 12, bloodstream 2, urinary tract 2, lungs + bloodstream 2.",
    dose_range = "Ceftazidime-avibactam 2.5 g (2000/500 mg) q8h in 17 subjects and 2.5 g q12h in 1 subject, each as a 2-hour IV infusion.",
    regions = "China (single centre, Hangzhou)",
    renal_function = "Continuous renal replacement therapy: CVVHDF in 14 and CVVH in 4 subjects, Prismaflex machine with the Prismaflex M150 set, blood flow 160-220 mL/min; actual operating duration on the sampling day 12-24 h (Hu 2026 Supplementary Tables 1-2). Cockcroft-Gault creatinine clearance median 68.62 mL/min (IQR 48.46-103.52), Table 1.",
    notes = "Prospective single-centre cohort, June 2021 to January 2024 (NCT05413343). Arterial plasma sampled at 0, 1, 2, 4, 6 and 8 h after the first and a steady-state infusion; first-dose and steady-state data fitted jointly. Assay UHPLC-Q-TOF, calibration 2-160 mg/L for avibactam; no BLQ samples. NONMEM 7.2, FOCE-I. The model describes TOTAL clearance (residual renal plus extracorporeal) with no time-varying on/off RRT term. Separate models were fitted for each analyte and for the RRT and non-RRT strata (sibling files Hu_2026_ceftazidime_nonrrt, Hu_2026_ceftazidime_rrt, Hu_2026_avibactam_nonrrt). Also screened and not retained: cumulative fluid balance, vasopressor use, infection site, actual RRT operating duration, weight-normalized total effluent rate, filtration fraction and mean net ultrafiltration rate."
  )

  ini({
    lcl <- log(1.77); label("Clearance (L/h)")            # Hu 2026 Table 2, RRT AVI CLPOP = 1.77 L/h (RSE 11%; bootstrap median 1.75, 95% CI 1.39-2.16)
    lvc <- log(14); label("Volume of distribution (L)")   # Hu 2026 Table 2, RRT AVI VPOP = 14 L (RSE 10%; bootstrap median 13.97, 95% CI 11.29-16.67)

    # Table 2 IIV entries are omega^2 VARIANCES (Methods Eq. 1); sqrt(0.196) =
    # 44% matches the 41% CV of the post-hoc CL in Table 3 (1.95 +/- 0.79 L/h).
    etalcl ~ 0.196  # Hu 2026 Table 2, RRT AVI IIV CL = 0.196 (RSE 30%)
    etalvc ~ 0.0814 # Hu 2026 Table 2, RRT AVI IIV V = 0.0814 (RSE 60%)

    # Table 2 RV is the proportional-error VARIANCE sigma^2 (Methods Eq. 3).
    propSd <- sqrt(0.0857); label("Proportional residual error (fraction)") # Hu 2026 Table 2, RRT AVI RV (proportional) = 0.0857 (RSE 22%)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
