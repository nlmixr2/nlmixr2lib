Wu_2026_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model with zero-order (intravenous",
    "infusion) input and first-order elimination for linezolid in",
    "critically ill adults in an emergency intensive care unit (Wu 2026).",
    "Central volume carries a body-weight power term, V = 30.7 *",
    "(WT/61.5)^2.16, and clearance a Cockcroft-Gault creatinine-clearance",
    "power term, CL = 7.37 * (CLcr/106.5)^1.35 (Results 3.2.2, Table 2).",
    "Inter-individual variability is exponential on V and CL; residual",
    "variability is exponential (log-normal). The paper uses the model in",
    "Monte Carlo simulations of infusion duration (0.5-12 h) against a",
    "composite efficacy (AUC24,ss/MIC >= 100 and %T>MIC >= 85%) and",
    "safety (Cmin,ss < 8 mg/L) target.",
    sep = " "
  )
  reference <- paste(
    "Wu C, Wang J, Sun W, Zhao W, Zhao W, Yang Q, Xue L, Kong L.",
    "Optimization of linezolid infusion duration in critically ill",
    "patients: a population pharmacokinetic analysis of the trade-off",
    "between efficacy and toxicity.",
    "Front Pharmacol. 2026;17:1766346.",
    "doi:10.3389/fphar.2026.1766346.",
    "PMCID PMC13364645.",
    sep = " "
  )
  vignette <- "Wu_2026_linezolid"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject (baseline). Wu 2026 Results 3.2.2 prints the",
        "final volume equation as V(L) = 30.7 x (WT/61.5)^2.16. Methods",
        "Equation 3 describes the centring value as the covariate median,",
        "while the Monte Carlo section calls 61.5 kg the study-population",
        "mean; Table 1 reports weight only per infusion group (1-h group",
        "65.2 +/- 11.7 kg, median 62, range 55-96; 3-h group 58.9 +/- 4.7 kg,",
        "median 59, range 52-65). The printed 61.5 kg is used as written.",
        sep = " "
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject (baseline). Wu 2026 Methods 2.1: 'CLcr was",
        "estimated using the Cockcroft-Gault equation'; no body-surface-area",
        "normalisation is mentioned and Table 1 reports the values in mL/min,",
        "so this is raw mL/min (the CRCL register entry accepts raw",
        "Cockcroft-Gault when the source does not BSA-normalise). Results",
        "3.2.2 prints CL (L/h) = 7.37 x (CLcr/106.5)^1.35; the Monte Carlo",
        "section calls 106.5 mL/min the study-population mean. Table 1: 1-h",
        "group 107.2 +/- 14 mL/min (median 110, range 71-124); 3-h group",
        "99.1 +/- 17.6 mL/min (median 97.5, range 72-128). Patients on renal",
        "replacement therapy were excluded.",
        sep = " "
      ),
      source_name = "CLcr"
    )
  )

  # Candidate covariates screened by stepwise covariate modelling (Methods
  # 2.2.2; Supplementary Table S1) but not retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained. Supplementary Table S1: age on V dOFV -1.345, age on CL dOFV -0.041 (both P > 0.05). Table 1: 60.8 +/- 6.1 years (1-h group), 55.7 +/- 8.2 years (3-h group)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained. Supplementary Table S1: sex on V passed the forward step alone (dOFV -4.605) but not once WT was on V (model 12, dOFV -2.404, P > 0.05); sex on CL dOFV -0.177. 9 of 20 patients female (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. Supplementary Table S1: albumin on V dOFV -0.026, on CL dOFV -0.573. Table 1: 44.1 +/- 6.0 g/L (1-h group), 43.0 +/- 6.9 g/L (3-h group)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "40-68 years (Table 1 ranges: 51-68 in the 1-h group, 40-67 in the 3-h group); adults 18 years or older by inclusion",
    weight_range = "52-96 kg (Table 1 ranges: 55-96 in the 1-h group, 52-65 in the 3-h group)",
    sex_female_pct = 45,
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 71-128 mL/min (Table 1),",
      "described as spanning mild impairment to the threshold of augmented",
      "renal clearance. Patients on renal replacement therapy were excluded.",
      sep = " "
    ),
    disease_state = paste(
      "Critically ill adults in an emergency intensive care unit with",
      "confirmed or suspected Gram-positive infection treated with",
      "intravenous linezolid.",
      sep = " "
    ),
    dose_range = paste(
      "Intravenous linezolid 600 mg every 12 h, infused over 1 h (10",
      "patients, standard infusion) or 3 h (10 patients, prolonged",
      "intermittent infusion).",
      sep = " "
    ),
    regions = "China (single centre: Emergency ICU, First Affiliated Hospital of Bengbu Medical University)",
    notes = paste(
      "Retrospective analysis of prospectively collected therapeutic drug",
      "monitoring data, January 2017 to January 2021. Of 42 ICU patients",
      "receiving linezolid, 20 were included (Supplementary Figure S1).",
      "Rich sampling after the first dose: pre-dose and 0.5, 1, 2, 3, 4, 6,",
      "8 and 12 h after the start of infusion, 9 samples per patient, 180",
      "in total (Table 1). Plasma linezolid measured by validated HPLC.",
      "Model fitted in NONMEM 7.5 with FOCE-I and PsN; evaluated by",
      "goodness-of-fit plots, a 1000-replicate bootstrap and a VPC. Race and",
      "ethnicity are not reported.",
      sep = " "
    )
  )

  ini({
    # Final-model estimates, Wu 2026 Table 2 ('Final model' column, values
    # with %RSE in parentheses) and the covariate equations printed in
    # Results 3.2.2. Reference subject: WT = 61.5 kg, CLcr = 106.5 mL/min.
    lcl <- log(7.37); label("Clearance at CLcr = 106.5 mL/min (L/h)")    # Wu 2026 Table 2 CL = 7.37 L/h (RSE 5.2%; bootstrap median 7.36, 95% CI 6.53-8.25)
    lvc <- log(30.7); label("Central volume of distribution at WT = 61.5 kg (L)") # Wu 2026 Table 2 V = 30.7 L (RSE 6%; bootstrap median 30.67, 95% CI 26.96-34.98)

    # Power exponents, Results 3.2.2: V = 30.7 x (WT/61.5)^2.16 and
    # CL = 7.37 x (CLcr/106.5)^1.35.
    e_wt_vc <- 2.16; label("Power exponent of (WT/61.5) on central volume (unitless)") # Wu 2026 Table 2 WT_V = 2.16 (RSE 21.5%; bootstrap median 2.19, 95% CI 0.51-2.99)
    e_crcl_cl <- 1.35; label("Power exponent of (CLcr/106.5) on CL (unitless)")       # Wu 2026 Table 2 CLcr_CL = 1.35 (RSE 22.9%; bootstrap median 1.35, 95% CI 0.72-1.93)

    # Inter-individual variability, exponential (Methods Equation 1). Methods
    # 2.2.1 defines the reported percentage as CV% = sqrt(omega^2) x 100, so
    # omega = CV%/100 and omega^2 = (CV%/100)^2. The two etas are reported
    # without a covariance, so they are independent.
    etalvc ~ 0.058081 # Wu 2026 Table 2 omega V = 24.1% (RSE 19.7%; shrinkage 6.6%) -> omega^2 = 0.241^2 = 0.058081
    etalcl ~ 0.043264 # Wu 2026 Table 2 omega CL = 20.8% (RSE 18.2%; shrinkage 4.8%) -> omega^2 = 0.208^2 = 0.043264

    # Residual error, Methods Equation 2: Cobs = Cpred x exp(eps), the
    # exponential model chosen over additive and combined forms. That is
    # nlmixr2's lnorm() residual, so 21.7% is the SD on the log-concentration
    # scale, read on the same sqrt(variance) x 100 scale as the IIV rows.
    expSd <- 0.217; label("Residual SD on the log-concentration scale (log-normal)") # Wu 2026 Table 2 epsilon = 21.7% (RSE 13.4%; shrinkage 10.1%; bootstrap median 21.79, 95% CI 19.71-23.77)
  })

  model({
    # Covariate centring values printed in the Results 3.2.2 equations.
    wt_norm <- WT / 61.5
    crcl_norm <- CRCL / 106.5

    # Individual PK parameters, Results 3.2.2 with exponential IIV
    # (Methods Equation 1).
    vc <- exp(lvc + etalvc) * wt_norm^e_wt_vc
    cl <- exp(lcl + etalcl) * crcl_norm^e_crcl_cl

    kel <- cl / vc

    # One-compartment model with zero-order input (Results 3.2.1): linezolid
    # is given as an intravenous infusion into `central`, so supply the
    # infusion duration (or rate) on the dosing record.
    d/dt(central) <- -kel * central

    # mg / L = mg/L, numerically identical to the ug/mL used in the paper's
    # simulation figures.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
