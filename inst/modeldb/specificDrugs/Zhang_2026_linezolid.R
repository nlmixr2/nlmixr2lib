Zhang_2026_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination for",
    "intravenous linezolid in critically ill Chinese children treated in a",
    "paediatric intensive care unit (Zhang 2026).",
    "Clearance carries two power covariates referenced to the cohort",
    "medians, CL = 2.80 * (WT/20.00)^0.69 * (eGFR/126.39)^0.34 (Table 2,",
    "Equation 6), where eGFR is the Schwartz-formula estimated glomerular",
    "filtration rate; central volume carries a single weight power term,",
    "V = 14.48 * (WT/20.00)^0.89 (Table 2, Equation 7). Weight and eGFR",
    "were the only covariates retained by stepwise forward inclusion and",
    "backward elimination (see covariatesDataExcluded). Inter-individual",
    "variability is exponential and was retained on CL only; the sparse,",
    "trough-dominant therapeutic-drug-monitoring design did not support an",
    "IIV term on V. Residual variability is proportional.",
    "The paper's second half trains a LightGBM machine-learning model on",
    "empirical-Bayes CL and V plus clinical features; that arm is a",
    "gradient-boosted tree ensemble rather than a structural model, and",
    "its fitted trees are not published, so only the population PK model",
    "is encoded here.",
    sep = " "
  )
  reference <- paste(
    "Zhang Y, Zhu L, Shen L, Wang G, Zhang J, Chen Y, Wang Y, Li Z.",
    "Machine learning combined with population pharmacokinetics: a hybrid",
    "model for predicting the plasma concentration of linezolid in",
    "critically ill pediatric patients.",
    "Front Pharmacol. 2026;17:1817282.",
    "doi:10.3389/fphar.2026.1817282.",
    "PMCID PMC13269215.",
    sep = " "
  )
  vignette <- "Zhang_2026_linezolid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Zhang 2026 Table 2 centres both",
        "the CL and the V power terms on 20.00 kg. Table 1 reports a",
        "per-sample median weight of 21.00 kg (IQR 9.10-33.50) in the",
        "training set and 20.00 kg (IQR 12.03-30.50) in the testing set, so",
        "20.00 kg is the cohort median rather than a rounded standard",
        "reference; the paper does not state the centring value in prose.",
        "Weight is the dosing descriptor as well as a covariate: children",
        "under 12 years were dosed 10 mg/kg q8h or q12h.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (Schwartz formula)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline). Computed by the authors with the",
        "Schwartz formula (Methods 2.1, Equation 1) using k = 0.45 for",
        "children under 1 year and k = 0.413 for children aged 1 year or",
        "older, so this is a creatinine-based, body-surface-area-normalised",
        "estimate rather than a measured clearance. Zhang 2026 Table 2",
        "centres the CL power term on 126.39 mL/min/1.73 m^2; Table 1",
        "reports a per-sample median of 124.98 (IQR 90.43-171.96) in the",
        "training set and 116.06 (IQR 84.02-153.46) in the testing set, so",
        "126.39 is the cohort median. That reference sits at the upper end",
        "of normal paediatric renal function, reflecting the mix of acute",
        "kidney injury and augmented renal clearance seen in the PICU",
        "(Discussion); the Monte Carlo dosing simulations span eGFR < 30 to",
        "200-400 mL/min/1.73 m^2.",
        sep = " "
      ),
      source_name        = "eGFR"
    )
  )

  # Variables collected for every patient (Zhang 2026 Methods 2.1) and carried
  # into the 32-feature machine-learning dataset, but NOT retained as
  # covariates in the final population PK model: "The covariate analysis
  # identified that weight and eGFR were significant covariates for clearance
  # (CL). Weight was the significant covariate for volume of distribution (V)"
  # (Results 3.2). The paper does not enumerate which of the collected
  # variables entered the stepwise popPK screen, nor report a point estimate
  # for any rejected covariate, so no effect size is recoverable for these.
  # Documented here to preserve the provenance of the paper's covariate set
  # without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Collected, not retained. Supplementary Table S3 reports a cohort median of 6.3 years (range 0.2-15.2); Table 1 per-sample medians 6.60 (training) and 6.20 (testing) years. Age also sets the dosing rule (10 mg/kg under 12 years, 600 mg q12h at 12 years and above)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected, not retained. Table 1 reports counts per sample rather than per patient: 63/149 female in the training set (42.28%) and 32/64 in the testing set (50.00%)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 15.72 (training) and 15.52 (testing) kg/m^2. Body size entered the model through WT rather than BMI."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 28.90 (training) and 32.80 (testing) umol/L. Serum creatinine is the input to the Schwartz eGFR that WAS retained, so the two are strongly correlated and only the derived eGFR survives the covariate screen. SCR is nonetheless the fifth-ranked feature of the LightGBM model by mean absolute SHAP value (Results 3.4)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 36.16 (training) and 35.93 (testing) g/L."
    ),
    TP = list(
      description = "Total serum protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 59.80 (training) and 59.47 (testing) g/L."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 21.80 (training) and 22.65 (testing) U/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 39.76 (training) and 38.92 (testing) U/L."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 8.40 (training) and 11.40 (testing) umol/L. Abbreviated TBIL in Table 1."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 3.80 (training) and 4.60 (testing) umol/L."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 212.00 (training) and 198.50 (testing) x10^9/L. Clinically relevant because linezolid overexposure is associated with thrombocytopenia (Introduction), which is why the paper carries a Cmin = 7 mg/L safety threshold."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 8.26 (training) and 7.78 (testing) x10^9/L."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^12 cells/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 3.05 (training) and 2.94 (testing) x10^12/L."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 85.00 (training) and 84.00 (testing) g/L."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Collected, not retained. Table 1 per-sample medians 14.00 (training) and 15.32 (testing) mg/L. The Discussion argues that systemic inflammation is a 'hidden covariate' acting through augmented renal clearance."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected as a yes/no flag, not retained. Table 1: 17/149 training samples (11.41%) and 3/64 testing samples (4.69%). The Discussion states explicitly that ECMO and CRRT entered the analysis 'only as binary covariates (yes/no) in the model, without integrating specific treatment parameters (e.g. ... replacement fluid rates, and operational modes for CRRT)', and names this simplification as a limitation that may underestimate PK variability."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected as a yes/no flag, not retained. Table 1: 9/149 training samples (6.04%) and 3/64 testing samples (4.69%). Same Discussion limitation as RRT_CRRT_STATUS: flow rates and oxygenator types were not modelled."
    ),
    CONMED_MEROPENEM = list(
      description = "Concomitant meropenem",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected, not retained. Table 1: 77/149 training samples (51.68%) and 37/64 testing samples (57.81%)."
    ),
    CONMED_OMEPRAZOLE = list(
      description = "Concomitant omeprazole",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected, not retained. Table 1: 60/149 training samples (40.27%) and 28/64 testing samples (43.75%)."
    ),
    CONMED_VORICONAZOLE = list(
      description = "Concomitant voriconazole",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected, not retained. Table 1: 21/149 training samples (14.09%) and 8/64 testing samples (12.50%)."
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole",
      units       = "(binary)",
      type        = "binary",
      notes       = "Collected, not retained. Table 1: 29/149 training samples (19.46%) and 12/64 testing samples (18.75%)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 145L,
    n_studies      = 1L,
    age_range      = "under 18 years by inclusion criterion; observed median 6.3 years (range 0.2-15.2) per Supplementary Table S3",
    weight_range   = "per-sample median 21.00 kg (IQR 9.10-33.50) in the training set; covariate model centred on 20.00 kg",
    sex_female_pct = 44.6,
    race_ethnicity = c(Asian = 100),
    renal_function = paste(
      "Schwartz eGFR per-sample median 124.98 mL/min/1.73 m^2 (IQR",
      "90.43-171.96) in the training set; covariate model centred on",
      "126.39. The PICU cohort spans acute kidney injury through augmented",
      "renal clearance, and the paper's Monte Carlo dosing simulations are",
      "stratified into eGFR bands of < 30, 30-59, 60-89, 90-129, 130-199",
      "and 200-400 mL/min/1.73 m^2.",
      sep = " "
    ),
    disease_state  = paste(
      "Critically ill children admitted to the paediatric intensive care",
      "unit and treated with intravenous linezolid for more than 3 days.",
      "Of the 145 patients, 69 (47.6%) had severe pneumonia, 43 (29.7%)",
      "sepsis, 19 (13.1%) central nervous system infection and 14 (9.6%)",
      "other conditions. Linezolid is used here against multidrug-resistant",
      "Gram-positive organisms including MRSA and vancomycin-resistant",
      "enterococci.",
      sep = " "
    ),
    dose_range     = paste(
      "Intravenous linezolid 10 mg/kg every 8 or 12 h for children under",
      "12 years, and 600 mg every 12 h for children aged 12 years and",
      "above, each given as an infusion over 1-2 h. Table 1 reports a",
      "per-sample median daily dose of 600.00 mg (IQR 240.00-930.00) in",
      "the training set.",
      sep = " "
    ),
    regions        = "China (single centre: Children's Hospital of Fudan University, National Children's Medical Center, Shanghai)",
    notes          = paste(
      "Retrospective single-centre study, January 2022 to October 2025,",
      "approved by the Ethics Committee of the Children's Hospital of Fudan",
      "University (No. 2025-506). 145 children contributed 213 steady-state",
      "linezolid concentrations (39 peak, 174 trough) after 5 samples below",
      "the limit of quantification were discarded by the M1 method. Routine",
      "therapeutic drug monitoring sampled half an hour after the end of",
      "infusion (peak) and half an hour before the next dose (trough), so",
      "the design is sparse and trough-dominant: time after dose is 7.5 h",
      "for roughly 68% of the training set. Concentrations measured by",
      "HPLC-UV at 253 nm, calibration range 0.25-50.00 mg/L, lower limit of",
      "quantification 0.25 mg/L. Model fitted in NONMEM VII with FOCE-I and",
      "evaluated by goodness-of-fit plots, a 1000-replicate nonparametric",
      "bootstrap (99.9% minimisation success), pvcVPC and NPDE (t-test",
      "P = 0.285, Shapiro-Wilk P = 0.408, Fisher variance test P = 0.208,",
      "global normality P = 0.623). A separate 10-patient, 13-sample",
      "external validation set (November-December 2025) was used only for",
      "the machine-learning arm. Supplementary Table S3 reports the",
      "predictive performance of this model on the study data as",
      "MDPE -14.41%, MAPE 37.48%, F20 33.62% and F30 49.10%, against",
      "18.17% / 56.91% / 15.05% / 22.58% for the two-compartment Yang 2021",
      "linezolid model evaluated on the same sparse dataset.",
      sep = " "
    )
  )

  ini({
    # Final-model fixed-effect estimates from Zhang 2026 Table 2 ("Final
    # model", Estimation column), which prints the covariate equations in its
    # row headers; the same two equations are repeated as Equations 6 and 7 in
    # Results 3.2. Reference subject: WT = 20.00 kg and eGFR = 126.39
    # mL/min/1.73 m^2, the centring constants named in those equations.
    lcl <- log(2.80); label("Clearance at WT = 20.00 kg and eGFR = 126.39 mL/min/1.73 m^2 (L/h)")  # Zhang 2026 Table 2 theta1 = 2.80 L/h (RSE 8.40%; bootstrap median 2.74, 95% CI 2.32-3.27; bias -2.14%)
    lvc <- log(14.48); label("Central volume of distribution at WT = 20.00 kg (L)")                 # Zhang 2026 Table 2 theta2 = 14.48 L (RSE 17.60%; bootstrap median 14.08, 95% CI 9.73-20.30; bias -2.76%)

    # Power exponents of the retained covariates.
    # Table 2 / Eq. 6: CL = theta1 * (WT/20.00)^theta3 * (eGFR/126.39)^theta4.
    # Table 2 / Eq. 7: V  = theta2 * (WT/20.00)^theta5.
    e_wt_cl   <- 0.69; label("Power exponent of (WT/20.00) on CL (unitless)")                # Zhang 2026 Table 2 theta3 = 0.69 (RSE 12.60%; bootstrap median 0.69, 95% CI 0.48-0.90). Estimated, not fixed at an allometric 0.75.
    e_crcl_cl <- 0.34; label("Power exponent of (eGFR/126.39) on CL (unitless)")             # Zhang 2026 Table 2 theta4 = 0.34 (RSE 23.10%; bootstrap median 0.34, 95% CI 0.20-0.54)
    e_wt_vc   <- 0.89; label("Power exponent of (WT/20.00) on central volume (unitless)")    # Zhang 2026 Table 2 theta5 = 0.89 (RSE 17.50%; bootstrap median 0.89, 95% CI 0.43-1.30). Estimated, not fixed at an allometric 1.

    # Inter-individual variability. Zhang 2026 Methods 2.4: "Inter-individual
    # variability was described using an exponential error model", i.e.
    # theta_i = theta_typical * exp(eta_i). Table 2 reports a single IIV term,
    # on CL, in a block headed 'Inter-individual variability (%)' as 48.79.
    # That percentage is read here as omega on the log scale (100 *
    # sqrt(omega^2)), not as an exact log-normal CV, because the SAME table
    # prints the proportional residual term in the same percent style, and
    # NONMEM / PsN report that quantity as 100 * sqrt(sigma^2). Reading one row
    # of the table on the SD scale and its neighbour on a CV-transformed scale
    # is not a convention that toolchain uses. The competing reading
    # (omega^2 = log(1 + CV^2)) would give 0.21357 instead of 0.23805, a 5%
    # difference in the eta SD; see the vignette Assumptions and deviations
    # section. No IIV on V was retained -- the trough-dominant sampling
    # carries little distribution-phase information (Discussion, limitation 2).
    etalcl ~ 0.23805  # Zhang 2026 Table 2 omega CL = 48.79% (RSE 11.20%; eta shrinkage 6%; bootstrap median 47.43, 95% CI 37.28-57.96; bias -2.79%) -> omega = 0.4879, omega^2 = 0.4879^2 = 0.23805; lognormal CV sqrt(exp(0.23805) - 1) = 51.8%

    # Residual error. Zhang 2026 Methods 2.4: "residual variability was
    # modeled using a proportional error structure"; Results 3.2 repeats "A
    # proportional model was used to explain the residual variability".
    propSd <- 0.4012; label("Proportional residual SD (unitless fraction of the prediction)")  # Zhang 2026 Table 2 sigma PROP = 40.12% (RSE 9.40%; epsilon shrinkage 24%; bootstrap median 39.37, 95% CI 31.46-47.75; bias -1.87%) -> 0.4012
  })

  model({
    # Covariate centring constants, Zhang 2026 Table 2 row headers and
    # Equations 6-7: CL and V are both centred on WT = 20.00 kg, and CL is
    # additionally centred on eGFR = 126.39 mL/min/1.73 m^2.
    wt_norm   <- WT / 20.00
    crcl_norm <- CRCL / 126.39

    # Individual PK parameters.
    # CL = 2.80 * (WT/20.00)^0.69 * (eGFR/126.39)^0.34   (Eq. 6)
    # V  = 14.48 * (WT/20.00)^0.89                        (Eq. 7)
    cl <- exp(lcl + etalcl) * wt_norm^e_wt_cl * crcl_norm^e_crcl_cl
    vc <- exp(lvc) * wt_norm^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition with first-order elimination. Linezolid was
    # given as an intravenous infusion over 1-2 h, so the dose enters `central`
    # directly; there is no absorption compartment and no bioavailability term.
    d/dt(central) <- -kel * central

    # Doses in mg with vc in L give central/vc in mg/L, the unit used
    # throughout Zhang 2026 (assay calibration range 0.25-50.00 mg/L, safety
    # threshold Cmin = 7 mg/L, MIC values 0.5-4 mg/L). mg/L is numerically
    # identical to ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
