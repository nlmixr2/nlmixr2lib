Ma_2025_lenalidomide <- function() {
  description <- "One-compartment oral population PK model with first-order absorption (Ka fixed) and first-order elimination for lenalidomide in Chinese adults with multiple myeloma or B-cell non-Hodgkin lymphoma, with Cockcroft-Gault creatinine clearance (power, centred at 94.2 mL/min) and fixed-exponent allometric body weight (centred at 63 kg) on CL/F and allometric body weight on V/F (Ma 2025)"
  reference <- "Ma Y, Song Z, Bing H, He H, Zhao L, Zhao R. Push Forward Clinical Management of Hematological Toxicity due to Lenalidomide Overexposure: Model-Informed Precision Dosing for Chinese Population With Renal Insufficiency. CPT Pharmacometrics Syst Pharmacol. 2025;14(7):1201-1210. doi:10.1002/psp4.70040"
  vignette <- "Ma_2025_lenalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "lenalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lenalidomide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed-exponent allometric scaling (0.75 on CL/F, 1.0 on V/F) centred at",
        "63 kg, the cohort median (Ma 2025 Table 1: mean 64.5, SD 11.7, median 63.0,",
        "range 37-101 kg). Included in the base model rather than selected by the",
        "covariate search (Ma 2025 Results 'PopPK Modeling'; supplement Material IV",
        "'stparm(V = tvV * (wt/63)^dVdwt)' and 'stparm(Cl = tvCl * (wt/63)^dCldwt ...)')."
      ),
      source_name = "wt / WT"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault formula",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalised",
        "(Ma 2025 Methods 'Covariate Model Development': 'CCr (Cockcroft-Gault",
        "formula)'); the body weight used inside the Cockcroft-Gault equation is not",
        "stated. Power effect on CL/F centred at 94.2 mL/min, the cohort median",
        "(Ma 2025 Table 1: mean 93.0, SD 31.6, median 94.2, range 34.8-227.7 mL/min).",
        "The fitted cohort had no subject below 34.8 mL/min, so the paper's dosing",
        "simulations at 10-30 mL/min extrapolate the power relationship."
      ),
      source_name = "CCr"
    )
  )

  covariatesDataExcluded <- list(
    DIS_MM = list(
      description = "Multiple myeloma indicator (1 = MM, 0 = B-cell NHL)",
      units = "(binary)",
      type = "binary",
      notes = "Disease type on CL/F passed forward inclusion (supplement Table S1 model 201, dOFV -9.593) but was removed at backward elimination because removing it raised the OFV by less than 10.83 (model 401). Cohort 55 MM (56.7%), 42 NHL (43.3%) (Ma 2025 Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F; not retained (supplement Table S1 models 102, 202, 301). Cohort 53 female (54.6%), 44 male (45.4%) (Ma 2025 Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 61.7, SD 10.9, median 62.0, range 33-87 years (Ma 2025 Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "(not reported)",
      type = "continuous",
      notes = "Screened on CL/F; significant alone (supplement Table S1 model 104, dOFV -8.696) but not retained once CCr was in the model (model 204, dOFV -0.635). Units and cohort summary are not reported."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 66.0, SD 10.9, median 65.8, range 40.8-116.6 g/L (Ma 2025 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 38.5, SD 4.7, median 38.0, range 25.3-49.3 g/L (Ma 2025 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 24.8, SD 21.2, median 20.0, range 2-148 U/L (Ma 2025 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 22.2, SD 12.7, median 20.2, range 1.1-89 U/L (Ma 2025 Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F; not retained (supplement Table S1). Cohort mean 82.1, SD 28.6, median 77.0, range 32-204 U/L (Ma 2025 Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 97,
    n_observations = 164,
    n_studies = 1,
    age_range = "33-87 years",
    age_median = "62 years",
    weight_range = "37-101 kg",
    weight_median = "63 kg",
    sex_female_pct = 54.6,
    race_ethnicity = c(Asian = 100),
    disease_state = "B-cell malignancy: multiple myeloma 55 (56.7%), B-cell non-Hodgkin lymphoma 42 (43.3%; DLBCL, FL, MCL, MZL, HGBL)",
    dose_range = "10-25 mg orally once daily for 21 days of a 28-day cycle (10 mg 23.7%, 12.5 mg 6.2%, 20 mg 1.0%, 25 mg 69.1%)",
    renal_function = "Cockcroft-Gault CCr 34.8-227.7 mL/min (median 94.2)",
    regions = "China (four hospitals in Beijing)",
    notes = paste(
      "Multicentre prospective real-world study, October 2021 to February 2022, at Peking",
      "University Third Hospital, Beijing Jishuitan Hospital, Beijing Chaoyang Hospital and",
      "Peking Union Medical College Hospital. Patients had taken lenalidomide for at least",
      "3 days (steady state); one trough sample 0.5 h before the dose and one sample 1 h",
      "after it were drawn on the same day, giving 164 concentrations in 97 patients.",
      "Lenalidomide given as monotherapy or in combination. Baseline demographics:",
      "Ma 2025 Table 1; by disease: Table 2."
    )
  )

  ini({
    # Structural parameters -- Ma 2025 Table 3 (final model).
    lka <- fixed(log(8.34)); label("Absorption rate constant (1/h)") # Ma 2025 Table 3 'tvKa (h-1) [17] = 8.34 (fixed)'; taken from Connarn 2018 (Ma 2025 reference 17) because the sparse absorption-phase data could not support estimation (Results 'PopPK Modeling')
    lcl <- log(7.4); label("Apparent clearance CL/F for a 63 kg subject with CCr 94.2 mL/min (L/h)") # Ma 2025 Table 3 'tvCL/F = 7.4 (5.66)'; bootstrap 7.4 (6.74-8.24)
    lvc <- log(37.4); label("Apparent volume of distribution V/F for a 63 kg subject (L)") # Ma 2025 Table 3 'tvV/F = 37.4 (6.37)'; bootstrap 37.4 (33.40-42.13)

    # Covariate effects -- Ma 2025 Table 3 equations
    #   CL/F (L/h) = tvCL/F x (CCr/94.2)^thetaCCr x (WT/63)^0.75
    #   V/F (L)    = tvV/F x (WT/63)^1.0
    e_crcl_cl <- 0.25; label("Power exponent of CCr/94.2 on CL/F (unitless)") # Ma 2025 Table 3 'thetaCCr = 0.25 (30.88)'; bootstrap 0.25 (0.14-0.34)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/63 on CL/F (unitless)") # Ma 2025 Table 3 CL/F equation '(WT/63)^0.75'; supplement Material IV 'fixef(dCldwt(freeze) = c(, 0.75, ))'
    e_wt_vc <- fixed(1); label("Allometric exponent of WT/63 on V/F (unitless)") # Ma 2025 Table 3 V/F equation '(WT/63)^1.0'; supplement Material IV 'fixef(dVdwt(freeze) = c(, 1, ))'

    # Between-subject variability on CL/F only (no BSV on Ka or V/F; Results
    # 'PopPK Modeling'). The text gives the eta variance as 0.025, which equals
    # the Table 3 15.8% read as omega (0.158^2 = 0.02496); the supplement
    # control stream likewise pairs the base-model 17.9% with omega^2 = 0.032.
    etalcl ~ 0.025 # Ma 2025 Results 'etaCL/F ... a variance of 0.025'; Table 3 'etaCL/F, %CV = 15.8 (22.2)'

    # Residual error -- proportional (Results 'PopPK Modeling'; supplement
    # Material IV 'observe(CObs = C * (1 + CEps))').
    propSd <- 0.495; label("Proportional residual error (fraction)") # Ma 2025 Table 3 'RV, %CV = 49.5 (13.0)'; bootstrap 48.9 (41.7-56.1)
  })

  model({
    # Individual parameters. Reference subject: WT = 63 kg, CRCL = 94.2 mL/min.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 94.2)^e_crcl_cl * (WT / 63)^e_wt_cl
    vc <- exp(lvc) * (WT / 63)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x 1000 converts to the ng/mL used
    # throughout Ma 2025.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
