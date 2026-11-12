Tan_2025_pazopanib <- function() {
  description <- paste(
    "One-compartment population PK model for oral pazopanib in real-world adults with metastatic",
    "renal cell carcinoma or soft tissue sarcoma under routine therapeutic drug monitoring (Tan 2025).",
    "First-order absorption with ka fixed to the Ozbey 2021 value, linear elimination, and a",
    "dose-nonlinear relative bioavailability F1 = (200 / DOSE)^0.42 with the 200 mg dose as the",
    "reference. Exponential inter-individual variability on V/F and F1 (none on CL/F) and a",
    "combined proportional-plus-additive residual error. The companion exposure-liver toxicity",
    "(Tan_2025_pazopanib_livertox) and tumor size models (Tan_2025_pazopanib_tumor_mrcc,",
    "Tan_2025_pazopanib_tumor_sts) are packaged separately.",
    sep = " "
  )
  reference <- paste(
    "Tan Z, Voller S, Yin A, Rieborn A, Gelderblom H, van der Hulle T, Knibbe CAJ, Moes DJAR.",
    "Model-Informed Dose Optimization of Pazopanib in Real-World Patients with Cancer.",
    "Clin Pharmacokinet. 2025;64(5):715-728. doi:10.1007/s40262-025-01504-5.",
    sep = " "
  )
  vignette <- "Tan_2025_pazopanib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "pazopanib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pazopanib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description = "Pazopanib dose amount of the current dosing record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying, record-level (register use case b). Drives the dose-nonlinear bioavailability",
        "F1 = TVF1 * (200 / DOSE)^EXP (Table 2 equation row). In the supplement control stream DOS is",
        "set to AMT on each dose record (`IF (AMT.GT.0) DOS = AMT`) and F1 uses it, so every dose",
        "record must carry its own amount in DOSE; observation records may carry the last dose. Must",
        "be > 0 on dose records. Doses in the cohort ranged 200-1000 mg once daily (Table 1), with",
        "some alternating-day regimens (Supplementary Fig. S2).",
        sep = " "
      ),
      source_name = "DOS (= AMT of the dose record)"
    )
  )

  # Covariates screened by the stepwise covariate search on CL/F and V/F and not retained
  # (Results, 'POPPK Model of Pazopanib'). CYP3A4 inhibitor / inducer use did not occur.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table 1: median 67 (range 21-89) years. Screened on CL/F and V/F, not retained."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Table 1: 85 male / 50 female. Screened on CL/F and V/F, not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Table 1: median 79 (range 45-130) kg. Screened on CL/F and V/F, not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Table 1: median 25.4 (range 17.6-42.5) kg/m^2. Screened on CL/F and V/F, not retained."
    ),
    TUMTP_SARC = list(
      description = "Sarcoma tumor-type indicator (1 = soft tissue sarcoma, 0 = metastatic renal cell carcinoma)",
      units = "(binary)",
      type = "binary",
      notes = "Table 1: 96 mRCC / 39 STS. Screened on CL/F (Methods), not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant gastric acid-suppressive agent",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Table 1: 25 / 135 (18%) used gastric acid-suppressive agents. Screened, not retained. The paper",
        "pools proton-pump inhibitors and other acid-suppressive agents in one group.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 135L,
    n_studies = 1L,
    n_concentrations = 460L,
    age_range = "21-89 years (median 67)",
    weight_range = "45-130 kg (median 79)",
    sex_female_pct = 37.0,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Adults with metastatic renal cell carcinoma (n = 96) or soft tissue sarcoma (n = 39) treated",
      "with pazopanib at Leiden University Medical Center, 2014-2023. 132 unique patients; patients",
      "who interrupted pazopanib for more than 6 months were entered as two individuals (135).",
      sep = " "
    ),
    dose_range = paste(
      "Oral pazopanib once daily, fasted; initial dose median 800 mg (range 200-800), median dose",
      "during treatment 600 mg (range 200-1000). Doses adjusted by therapeutic drug monitoring.",
      sep = " "
    ),
    regions = "The Netherlands (Leiden University Medical Center)",
    notes = paste(
      "460 plasma concentrations (326 troughs, 20-28 h post-dose), median 3 per patient; time after",
      "dose 0.5-166 h. LC-MS/MS assay; no sample below the LLOQ. NONMEM 7.4.4, FOCE-I (Table 1,",
      "Methods).",
      sep = " "
    )
  )

  ini({
    lka <- fixed(log(0.976)); label("Absorption rate constant (1/h)") # Table 2 row 'Ka (h)' 0.976 FIX, taken from Ozbey 2021; control stream THETA(3) (0.976) FIX
    lcl <- log(0.497); label("Apparent clearance CL/F (L/h)") # Table 2 row 'CL/F (L/h)' 0.497 (RSE 7%); control stream THETA(1)
    lvc <- log(46.1); label("Apparent volume of distribution V/F (L)") # Table 2 row 'V/F (L)' 46.1 (RSE 19%); control stream THETA(2)
    lfdepot <- fixed(log(1)); label("Relative bioavailability at the 200 mg reference dose (unitless)") # Table 2 row 'TVF1' 1 FIX; control stream THETA(6) (1) FIX
    e_dose_fdepot <- 0.42; label("Power exponent on (200 / DOSE) for relative bioavailability (unitless)") # Table 2 row 'EXP' 0.42 (RSE 12%); control stream THETA(7)

    # Table 2 reports IIV as CV% = sqrt(omega): the control stream $OMEGA
    # finals 0.58 (V) and 0.12 (F1) give sqrt = 76.2% and 34.6%, matching the
    # printed 76.5% and 34.2%. Variances below are the squared Table 2 values.
    etalvc ~ 0.585225 # Table 2 row 'IIV_V/F (CV%)' 76.5% -> 0.765^2; control stream $OMEGA 0.58
    etalfdepot ~ 0.116964 # Table 2 row 'IIV_F1 (CV%)' 34.2% -> 0.342^2; control stream $OMEGA 0.12

    propSd <- 0.24; label("Proportional residual error (fraction)") # Table 2 row 'Proportional error (CV%)' 24%; control stream THETA(4) (0.24)
    addSd <- 4.71; label("Additive residual error (mg/L)") # Table 2 row 'Additive error (mg/L)' 4.71 (RSE 20%); control stream THETA(5)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose-nonlinear relative bioavailability (Table 2 equation row;
    # control stream F1 = THETA(6) * EXP(ETA(2)) * (200 / DOS)**THETA(7)).
    f(depot) <- exp(lfdepot + etalfdepot) * (200 / DOSE)^e_dose_fdepot

    # Dose in mg over V/F in L gives mg/L.
    Cc <- central / vc
    # Control stream W = SQRT(THETA(4)**2 * IPRED**2 + THETA(5)**2): combined2.
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
