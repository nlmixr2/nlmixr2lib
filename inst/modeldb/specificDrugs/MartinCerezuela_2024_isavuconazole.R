MartinCerezuela_2024_isavuconazole <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous isavuconazole (dosed as the",
    "prodrug isavuconazonium sulfate, doses expressed as isavuconazole) in critically",
    "ill adults in a medical intensive care unit, a third of them on continuous renal",
    "replacement therapy and a third on extracorporeal membrane oxygenation, with a 1-h",
    "infusion into the central compartment and first-order elimination",
    "(Martin-Cerezuela 2024). No covariates were retained. The residual error model is",
    "not reported and is carried as zero.",
    sep = " "
  )
  reference <- paste(
    "Martin-Cerezuela M, Maya Gallegos C, Marques-Minana MR, Broch Porcar MJ,",
    "Cruz-Sanchez A, Mateo-Pardo JC, Peris Ribera JE, Gimeno R, Castellanos-Ortega A,",
    "Poveda Andres JL, Ramirez Galleymore P, Lipman J, Di Paolo A, Breilh D.",
    "Isavuconazole Pharmacokinetics in Critically Ill Patients: Relationship with",
    "Clinical Effectiveness and Patient Safety. Antibiotics (Basel). 2024;13(8):706.",
    "doi:10.3390/antibiotics13080706. PMCID PMC11350865.",
    sep = " "
  )
  vignette <- "MartinCerezuela_2024_isavuconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Checked against Martin-Cerezuela 2024 Methods 2.1-2.2.
  compartmentData <- list(
    central = list(
      analyte = "isavuconazole",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Methods 2.2: blood collected in EDTA tubes and centrifuged to plasma; total plasma isavuconazole measured by HPLC with fluorescence detection. Doses are 1-h intravenous infusions into central (Methods 2.1)."
    )
  )

  covariateData <- list()

  # Methods 2.2 lists the covariates screened by stepwise forward selection;
  # Results 3.2: 'None of the studied covariates improved the predictive
  # ability of the model, so no covariates were introduced in the final model.'
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort mean 53.2 years (SD 16.4), Table 1."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). 17 of 24 patients male (70.8 percent), Table 1."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort mean 75.1 kg (SD 12.6), Table 1."
    ),
    SAPS3 = list(
      description = "Simplified Acute Physiology Score 3 at ICU admission",
      units = "points",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort median 63 (IQR 57-71.5), Table 1. SAPS 3 is a different score from SAPS II (canonical SAPS_II), so it is not mapped onto that column; same non-canonical name as Swartling_2022_cefotaxime.R."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort median 0.61 (IQR 0.47-0.79), Table 1; the unit is not printed and mg/dL is inferred from the magnitude."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort mean 1.1 (SD 1.2), Table 1. Table 1 prints the unit as umol/L, but values of this size are mg/dL (1.1 umol/L would be far below the normal range), so mg/dL is recorded here."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort median 60 U/L (IQR 28-70), Table 1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort median 51.5 U/L, Table 1 (IQR printed as '22.102')."
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort median 270 U/L (IQR 211-309), Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2). Cohort mean 2.4 (SD 0.4), Table 1. Table 1 prints the unit as g/L, but 2.4 is a g/dL value (24 g/L)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2; Figure 2 shows the individual CL and V by CRRT status). 8 of 24 patients (33.3 percent), Table 1."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V (Methods 2.2), not retained (Results 3.2; Figure 2 shows the individual CL and V by ECMO status). 7 of 24 patients (29.1 percent), Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 71L,
    age_mean = "53.2 years (SD 16.4)",
    weight_mean = "75.1 kg (SD 12.6)",
    sex_female_pct = 29.2,
    race_ethnicity = "not reported in the source paper",
    disease_state = paste(
      "Critically ill adults in a medical ICU treated with intravenous isavuconazole for",
      "suspected or confirmed invasive fungal infection. Admission causes: acute",
      "respiratory failure 62.5 percent, cardiogenic shock 16.7 percent, septic shock",
      "12.5 percent, acute stroke 8.3 percent. 23 of 24 mechanically ventilated; CRRT",
      "33.3 percent; ECMO 29.1 percent. Fungal infection confirmed in 11 (candidemia 7,",
      "mostly Candida auris; Candida parapsilosis ventriculitis 1; invasive pulmonary",
      "aspergillosis 2; mucormycosis 1). SAPS3 median 63. ICU mortality 58.3 percent.",
      sep = " "
    ),
    dose_range = "200 mg every 8 h for six doses (loading), then 200 mg every 24 h, each as a 60-min intravenous infusion",
    regions = "Spain (24-bed medical ICU, tertiary university hospital, Valencia)",
    notes = paste(
      "Prospective observational study, January to December 2023 (Methods 2.1).",
      "Samples drawn after 48-96 h of treatment at predose (Cmin), 1 h (Cmax) and 12 h",
      "(C50) after the infusion; 71 plasma concentrations from 24 patients. Observed",
      "means: Cmin 1.76 (SD 1.02), Cmax 3.99 (SD 1.96), C50 2.26 (SD 1.20) mg/L",
      "(Results 3.2). NONMEM ADVAN1 TRANS2, FOCE with interaction. Baseline",
      "characteristics from Table 1.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Martin-Cerezuela 2024 Table 2 ('Pharmacokinetic
    # parameters estimated from the isavuconazole model'), column 'Mean (SE)'.
    lcl <- log(5.36); label("Clearance (L/h)")
    # Table 2, 'Clearance (L/h)' = 5.36 (SE 0.576; CI95% 4.23-6.49)
    lvc <- log(122); label("Volume of distribution (L)")
    # Table 2, 'Volume of distribution (L)' = 122 (SE 20.7; CI95% 81.43-162.57)

    # Inter-individual variability: Table 2 column '%CV', which the preprint
    # footnote defines as the 'percent covariance of intersubject-variability'.
    # Exponential etas (ADVAN1 TRANS2); converted with omega^2 = log(CV^2 + 1).
    # The paper does not say whether its CV is sqrt(omega^2) instead; see the
    # vignette Assumptions.
    etalcl ~ 0.186653 # Table 2 Clearance %CV 45.3 -> log(0.453^2 + 1)
    etalvc ~ 0.370805 # Table 2 Volume of distribution %CV 67 -> log(0.67^2 + 1)

    # Residual unexplained variability is not reported anywhere in the paper
    # (main text or the Research Square preprint): no error model and no
    # magnitude. Carried as fixed(0) rather than invented; see vignette Errata.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd <- fixed(0); label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
