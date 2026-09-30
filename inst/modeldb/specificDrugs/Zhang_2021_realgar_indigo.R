Zhang_2021_realgar_indigo <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for total plasma arsenic after oral tetra-arsenic tetra-sulfide (As4S4) formula (Realgar-Indigo Naturalis Formula, RIF) given three times daily to Chinese children with acute promyelocytic leukemia aged 4-14 years (Zhang 2021). Dose is the mass of RIF formula (mg) and the observation is total arsenic concentration (ug/L), so CL/F and V/F are apparent values relative to the formula dose. Apparent oral clearance CL/F scales with body weight as a power function with an estimated exponent of 0.629 referenced to the 27 kg cohort median; V/F and ka carry no covariates. Very slow absorption (ka 0.013 1/h, flip-flop kinetics). Diagonal exponential inter-individual variability on ka, CL/F and V/F; additive residual error."
  reference <- "Zhang L, Yang XM, Chen J, Hu L, Yang F, Zhou Y, Zhao BB, Zhao W, Zhu XF. Population Pharmacokinetics and Safety of Oral Tetra-Arsenic Tetra-Sulfide Formula in Pediatric Acute Promyelocytic Leukemia. Drug Des Devel Ther. 2021;15:1633-1640. doi:10.2147/DDDT.S305244"
  vignette <- "Zhang_2021_realgar_indigo"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(
      analyte = "Realgar-Indigo Naturalis Formula (As4S4 formula)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "Realgar-Indigo Naturalis Formula (As4S4 formula) dose-equivalent, observed as total arsenic",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F, (WT/27)^0.629, referenced to the cohort median weight of 27 kg (Table 2 'F WT-CL = (CW/27)^theta4'). Cohort mean (SD) 30.03 (14.06) kg, median 27 kg, range 16.0-63.0 kg (Table 1). Treated as time-fixed.",
      source_name = "CW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the PK parameters but not retained in the final model (Results, 'Population Pharmacokinetic Analysis'). Cohort median 8 years, range 4-14 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 32.6 umol/L, range 18.4-58.4 (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 43.8 g/L, range 40.4-49.6 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 26.2 U/L, range 17.7-38.0 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 15.2 U/L (Table 1 prints the range as '9.0-4.4', an evident typo)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    n_observations = 107L,
    age_range = "4-14 years",
    age_median = "8 years",
    age_mean_sd = "7.73 (3.169) years",
    weight_range = "16.0-63.0 kg",
    weight_median = "27 kg",
    weight_mean_sd = "30.03 (14.06) kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Chinese (12 of 12; Table 1).",
    disease_state = "Pediatric acute promyelocytic leukemia (PML/RARa-positive) in complete remission: 11 newly diagnosed patients and 1 relapsed patient. Seven received RIF in maintenance therapy (regimen A) and five in consolidation therapy (regimen B), together with ATRA; patients with abnormal renal or liver function were excluded.",
    dose_range = "RIF 60 mg/kg/day divided three times daily (20 mg/kg per administration); actual per-administration dose median 540 mg, range 320-1350 mg (Table 1).",
    regions = "China (Institute of Hematology and Blood Diseases Hospital, Tianjin).",
    sampling = "Pre-dose and 1, 2, 4, 6, 7, 8, 9, 10, 16, 20, 24 and 28 h after administration on day 1; weekly pre-dose samples from day 8 to day 28; four samples over the 2 weeks after stopping RIF.",
    assay = "Total plasma arsenic by ICP-MS (Agilent 7700x); calibration range 0.015-50 ug/L; lower limit of detection 0.015 ug/L. Observed concentrations 0.1-75.0 ug/L.",
    notes = "Single-centre prospective open-label study, July 2016 to July 2019 (ChiCTR-OIC-16010014). NONMEM 7.2, FOCE with interaction. Sex distribution not reported."
  )

  ini({
    lka <- log(0.013); label("Absorption rate constant ka (1/h)") # Table 2 theta1 = 0.013 h-1 (RSE 17.3%; bootstrap median 0.013, 5th-95th 0.0098-0.0180)
    lcl <- log(1380); label("Apparent oral clearance CL/F at WT = 27 kg (L/h)") # Table 2 theta2 = 1380 L/h (RSE 7.0%; bootstrap median 1350, 5th-95th 1220-1530)
    lvc <- log(7080); label("Apparent volume of distribution V/F (L)") # Table 2 theta3 = 7080 L (RSE 44.6%; bootstrap median 6450, 5th-95th 3370-13600)
    e_wt_cl <- 0.629; label("Power exponent of (WT/27) on CL/F (unitless)") # Table 2 theta4 = 0.629 (RSE 23.2%; bootstrap median 0.622, 5th-95th 0.361-0.887)

    # Table 2 prints the IIV rows under 'Inter-individual variability (%)' as
    # 0.353 (Ka), 0.167 (CL) and 0.787 (V). They are read as omega, the SD of
    # eta on the log scale, and squared here. The variance reading is ruled
    # out by the row RSEs: a variance estimated from N = 12 subjects cannot
    # have a relative SE below sqrt(2/12) = 40.8%, yet the printed CL RSE is
    # 39.3% and the bootstrap 5th-95th intervals imply RSEs of about 32% (Ka)
    # and 22% (CL). The SD and lognormal-CV readings cannot be separated by
    # the table; the CV reading would give omega^2 = log(1 + CV^2) = 0.117,
    # 0.0275 and 0.486.
    etalka ~ 0.124609 # Table 2 IIV Ka 0.353 (RSE 48.2%; bootstrap median 0.339, 5th-95th 0.075-0.445) -> 0.353^2
    etalcl ~ 0.027889 # Table 2 IIV CL 0.167 (RSE 39.3%; bootstrap median 0.139, 5th-95th 0.078-0.196) -> 0.167^2
    etalvc ~ 0.619369 # Table 2 IIV V 0.787 (RSE 52.6%; bootstrap median 0.781, 5th-95th 0.023-1.148) -> 0.787^2

    # Table 2 prints 'Residual variability (%)' = 3.619 and the Results text
    # says a proportional model was used. A 3.6% proportional error is not
    # compatible with the paper's own Figure 1B (DV vs IPRED): digitised, the
    # individual-fit residuals are about 10% at 40-70 ug/L and 50-250% at
    # 1-15 ug/L, i.e. roughly constant in absolute size (about 3-4 ug/L), which
    # is the signature of an additive error. The value is therefore encoded as
    # an additive SD in ug/L. See the vignette for the digitisation.
    addSd <- 3.619; label("Additive residual error (ug/L)") # Table 2 'Residual variability' 3.619 (RSE 31.9%; bootstrap median 3.688, 5th-95th 2.747-4.609)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 27)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg of RIF formula, vc in L -> mg/L; x1000 for ug/L arsenic.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
