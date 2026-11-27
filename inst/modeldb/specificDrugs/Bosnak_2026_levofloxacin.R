Bosnak_2026_levofloxacin <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model with first-order",
    "absorption and first-order elimination for oral levofloxacin",
    "prophylaxis (750 mg once daily) in adults with acute leukaemia or",
    "undergoing allogeneic haematopoietic stem cell transplantation during",
    "neutropenia, from Bosnak 2026. Cockcroft-Gault creatinine clearance",
    "enters apparent clearance as a power function centred on 94 mL/min,",
    "and body weight scales the apparent volume linearly (exponent fixed at",
    "1) centred on 66.6 kg. The absorption rate constant is fixed. Log-normal",
    "between-subject variability on CL/F and V/F; proportional residual",
    "error. The unbound concentration Cu = 0.7 * Cc carries the 30% plasma",
    "protein binding the authors assumed for the fAUC0-24/MIC target.",
    sep = " "
  )
  reference <- paste(
    "Bosnak C, Metan G, Zarakolu P, Aladag Karakulak E, Pinar A, Goker H,",
    "Roberts JA, Hernandez-Mitre MP, Akova M.",
    "Determinants of levofloxacin prophylaxis failure in high-risk",
    "haematology patients: a prospective study integrating resistant",
    "colonization and pharmacokinetic/pharmacodynamic target attainment.",
    "J Antimicrob Chemother. 2026;81(9):dkag284. doi:10.1093/jac/dkag284.",
    "PMC13464500.",
    "Final estimates are Table 3; the Monolix model code is Supplementary",
    "Methods Section S2; the PK-analysis demographics are Supplementary",
    "Table S2.",
    sep = " "
  )
  vignette <- "Bosnak_2026_levofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear (power exponent fixed at 1) effect on V/F centred on 66.6 kg,",
        "the PK-analysis median (Table 3 footnote 'Continuous covariates were",
        "normalized to their median values (WT = 66.6 kg ...)'; Supplementary",
        "Table S2 median 66.6 kg, IQR 59.5-72.0). Supplementary S2 code:",
        "logtWT = log(WT/66.6) with coefficient beta_V_logtWT.",
        sep = " "
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated by the Cockcroft-Gault equation, in",
        "raw mL/min (NOT BSA-normalized).",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent 0.81) centred on 94 mL/min, the",
        "PK-analysis median (Table 3 footnote; Supplementary Table S2 median",
        "94 mL/min, IQR 75-143). Supplementary S2 code: logtCRCL =",
        "log(CRCL/94). The body-weight convention inside the Cockcroft-Gault",
        "equation is not stated by the source. Measured once per neutropenic",
        "episode, at the PK sampling day.",
        sep = " "
      ),
      source_name = "CRCL"
    )
  )

  # Covariates the source screened and did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search (Methods, 'Population pharmacokinetic modeling') and not retained. No coefficient is reported."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained (Cockcroft-Gault CRCL was retained instead). No coefficient is reported."
    ),
    TPRO = list(
      description = "Serum total protein.",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score.",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained. No coefficient is reported."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score.",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search and not retained (Supplementary Table S2 median 14, IQR 12-16). No coefficient is reported. The MDRD and CKD-EPI eGFR variants were also screened against Cockcroft-Gault CRCL and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 103,
    n_studies = 1,
    n_patients = 108,
    n_observations = "plasma samples at 2, 4, 8 and 24 h after one 750 mg dose in each of 103 neutropenic episodes",
    age_median = "58 years (IQR 34-68)",
    weight_median = "66.6 kg (IQR 59.5-72.0)",
    sex_female_pct = 42.7,
    disease_state = paste(
      "adults with acute leukaemia (88% acute myeloid leukaemia) or",
      "undergoing allogeneic haematopoietic stem cell transplantation, with",
      "an absolute neutrophil count below 500 cells/mm^3 for at least one",
      "week (median 70 cells/mm^3, IQR 10-290)"
    ),
    renal_function = "Cockcroft-Gault creatinine clearance median 94 mL/min (IQR 75-143); serum creatinine median 0.77 mg/dL (IQR 0.58-0.96)",
    dose_range = "oral levofloxacin 750 mg once daily as febrile-neutropenia prophylaxis for at least 7 days",
    regions = "Turkiye (Hacettepe University Hospitals, Ankara; single centre)",
    notes = paste(
      "Prospective case-control study of 120 neutropenic episodes in 108",
      "patients (60 prophylaxis successes, 60 failures); PK sampling was",
      "complete in 103 episodes (60 successes, 43 failures), which form the",
      "PK analysis set (Results; Supplementary Table S2). A patient could",
      "contribute more than one episode if they were at least a week apart,",
      "so n_subjects counts episodes. Failures were sampled on the day after",
      "febrile-neutropenia onset, successes after a scheduled prophylactic",
      "dose. Patients with fluoroquinolone-resistant colonization or",
      "infection in the preceding three months were excluded. Total plasma",
      "levofloxacin by LC-MS/MS (LLOQ 0.13 mg/L, Supplementary Table S1).",
      "Fitted in Monolix 2024R1 (SAEM)."
    )
  )

  ini({
    # Structural parameters. Table 3 'Final estimate' column; parameters are
    # APPARENT (CL/F, V/F) because only oral data were collected (Discussion,
    # limitations). Supplementary S2 code: Cc = pkmodel(ka, V, Cl), F = 1.
    lka <- fixed(log(1.12)); label("First-order absorption rate constant, ka (1/h)") # Table 3 ka = 1.12 (fixed); Results: sampling began after the expected Tmax, so ka was fixed to the preliminary estimate from model development
    lcl <- log(4.77); label("Apparent clearance at CRCL 94 mL/min, CL/F (L/h)") # Table 3 CL/F = 4.77 L/h (RSE 4.06%; bootstrap median 4.73, 95% CI 4.34-5.11)
    lvc <- log(97.57); label("Apparent volume of distribution at WT 66.6 kg, V/F (L)") # Table 3 V/F = 97.57 L (RSE 4.78%; bootstrap median 97.25, 95% CI 88.76-106.95)

    # Covariate effects. Supplementary S2 code: Monolix 'covariate=logtCRCL'
    # and 'covariate=logtWT' on log-normal parameters, i.e. power functions of
    # the median-normalised covariate; Table 3 'Structural model equations'.
    e_crcl_cl <- 0.81; label("Power exponent of CRCL/94 on CL/F (unitless)") # Table 3 'Effect of CRCL on CL/F' = 0.81 (RSE 9.46%; bootstrap median 0.84, 95% CI 0.69-1.01)
    e_wt_vc <- fixed(1); label("Power exponent of WT/66.6 on V/F (unitless)") # Table 3 'Effect of WT on V/F' = 1 (fixed)

    # Between-subject variability. Table 3 prints Monolix C.V.(%) values; for
    # log-normal parameters Monolix computes CV = sqrt(exp(omega^2) - 1), so
    # omega^2 = log(CV^2 + 1). Diagonal: the S2 code declares no correlation.
    # ka carries no variability ('no-variability' in the S2 code).
    etalcl ~ 0.150832 # Table 3 IIV_CL/F = 40.35 CV% -> log(0.4035^2 + 1) (RSE 7.83%, shrinkage 5.61%)
    etalvc ~ 0.117875 # Table 3 IIV_V/F = 35.37 CV% -> log(0.3537^2 + 1) (RSE 11.7%, shrinkage 24.7%)

    # Plasma protein binding: a literature value the authors applied to
    # compute the fAUC0-24/MIC target, not an estimate. Supplementary S2 code:
    # 'Cfree = Cc * 0.7'.
    fu <- fixed(0.7); label("Fraction of levofloxacin unbound in plasma (unitless)") # Methods 'PK/PD target attainment': protein binding assumed 30%; Supplementary S2 Cfree = Cc * 0.7

    # Residual error. Supplementary S2 code: errorModel = proportional(b).
    propSd <- 0.25; label("Proportional residual SD (fraction)") # Table 3 proportional error = 25.0 CV% (RSE 4.72%; bootstrap median 26.0, 95% CI 22.0-29.0)
  })

  model({
    # 1. Individual parameters (Table 3 'Structural model equations':
    #    V/F = 97.57 x (WT/66.6); CL/F = 4.77 x (CRCL/94)^0.81).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 94)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 66.6)^e_wt_vc

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment model with first-order absorption (S2 'pkmodel(ka,
    #    V, Cl)'); oral doses go into depot.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Observation. The assay measured TOTAL plasma levofloxacin (Cc). Cu
    #    is the unbound concentration whose AUC0-24 over the MIC defines the
    #    fAUC0-24/MIC target; the S2 code integrates it in its own fAUCtau
    #    table output.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
