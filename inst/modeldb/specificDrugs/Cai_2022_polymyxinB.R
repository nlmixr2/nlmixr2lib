Cai_2022_polymyxinB <- function() {
  description <- "One-compartment intravenous population PK model for polymyxin B in Chinese adult lung transplant recipients with carbapenem-resistant Gram-negative pneumonia (Cai 2022). Cockcroft-Gault creatinine clearance is the sole retained covariate, entering clearance as a power term normalized to the cohort median of 78.49 mL/min with exponent 0.681. Independent exponential inter-individual variability on CL and V; proportional residual error."
  reference <- paste(
    "Cai X-J, Chen Y, Zhang X-S, Wang Y-Z, Zhou W-B, Zhang C-H, Wu B, Song H-Z, Yang H, Yu X-B.",
    "Population pharmacokinetic analysis, renal safety, and dosing optimization",
    "of polymyxin B in lung transplant recipients with pneumonia: A prospective study.",
    "Front Pharmacol. 2022;13:1019411.",
    "doi:10.3389/fphar.2022.1019411. PMCID PMC9608142.",
    sep = " "
  )
  vignette <- "Cai_2022_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Cai 2022: polymyxin B in plasma by LC-MS/MS
  # (Methods 2.3), described by a one-compartment model with first-order
  # elimination (Results 3.3).
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance calculated with the Cockcroft-Gault equation. NOT body-surface-area normalised: the source reports raw Cockcroft-Gault mL/min.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Cai 2022 Table 4 footnote: 'CrCL, creatinine clearance calculated using the Cockcroft-Gault equation'. The body-weight convention inside the Cockcroft-Gault equation is not stated. Applied as a power covariate on CL, CL = 1.72 * (CrCL / 78.49)^0.681 (Results Eq. 5, read as multiplicative; see the vignette), where 78.49 mL/min is 'the median value of CrCL for the included patients'. Table 1 cohort CrCL 80.81 +/- 29.97 mL/min. Stored under the canonical CRCL column with raw mL/min units, following the raw Cockcroft-Gault precedents in the CRCL register entry (e.g. Wang_2020_polymyxinB.R).",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    SCREENED_NOT_RETAINED = list(
      description = "Candidate covariates screened by stepwise forward selection / backward deletion on CL and V and not retained",
      units = "(various)",
      type = "continuous",
      notes = "Cai 2022 Methods 2.4.2: age, gender, body weight, height, haemoglobin, white blood cell count, percentage of neutrophils, platelet count, ALT, AST, total bilirubin, serum albumin, serum total protein, CRP, procalcitonin, blood urea nitrogen, serum creatinine, CrCL, and concomitant furosemide or albumin were screened. Only CrCL on CL was retained (Results 3.3). The Discussion states that body weight (including an allometric form) did not improve the fit and that concomitant furosemide or albumin had no significant effect. Co-administered inhaled polymyxin B (25 mg q12h, given to every patient) could not be assessed as a covariate."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "adults >= 18 years; mean 56 +/- 12.76 years (Table 1)",
    weight_range = "mean 52.15 +/- 10.00 kg (Table 1)",
    height_range = "mean 166.93 +/- 7.20 cm (Table 1)",
    sex_female_pct = 26.47,
    race_ethnicity = "Chinese (single-centre cohort in Wuxi, China)",
    disease_state = "Adult lung transplant recipients with pneumonia caused by carbapenem-resistant Gram-negative organisms (A. baumannii 47.5%, P. aeruginosa 30%, K. pneumoniae 17.5%, E. cloacae 5% of isolates), treated with polymyxin B for >= 3 days. No concentrations were obtained during extracorporeal membrane oxygenation or renal replacement therapy.",
    dose_range = "Intravenous polymyxin B sulfate q12h infused over 1 h; guideline loading dose 2.0-2.5 mg/kg and maintenance 1.25-1.5 mg/kg q12h; daily IV dose median 100 (range 100-150) per Table 1. Every patient also received inhaled polymyxin B sulfate 25 mg q12h.",
    regions = "China (Affiliated Wuxi People's Hospital of Nanjing Medical University), January 2020 to December 2021",
    renal_function = "Cockcroft-Gault CrCL 80.81 +/- 29.97 mL/min (Table 1), median 78.49 mL/min (Results 3.3); serum creatinine 72.40 +/- 33.45 umol/L.",
    n_observations = "164 plasma polymyxin B concentrations (0.56-11.66 mg/L); four samples per patient at least 48 h after starting therapy: 0.5 h before an infusion and 1, 2 and 6 h after the end of the infusion.",
    notes = "Single-centre prospective study. NONMEM 7.4 with FOCE-I. Evaluation by goodness-of-fit plots, 1000-sample nonparametric bootstrap (Table 3), prediction- and variability-corrected VPC and NPDE. Monte Carlo PTA (fAUC/MIC >= 20, fu = 0.42) and AUCss,24h safety (< 100 mg*h/L) simulations by CrCL (Table 4, Figures 4-5)."
  )

  ini({
    # Structural parameters -- Cai 2022 Table 3 (final model) and Results
    # Eq. 5-6. The typical clearance is the value at CrCL = 78.49 mL/min.
    lcl <- log(1.72) ; label("Clearance CL (L/h) at CRCL = 78.49 mL/min") # Table 3: theta CL = 1.72 L/h (RSE 8%; bootstrap median 1.72, 95% CI 1.44-1.99); Eq. 5
    lvc <- log(14.4) ; label("Volume of distribution V (L)")              # Table 3: theta V = 14.4 L (RSE 11%; bootstrap median 14.3, 95% CI 11.4-17.4); Eq. 6

    # Covariate effect on CL -- Cai 2022 Eq. 5, CL = 1.72 * (CrCL / 78.49)^0.681
    # (the typeset equation prints '+' between the two factors; see the
    # vignette for why the multiplicative reading is used).
    e_crcl_cl <- 0.681 ; label("Power exponent on (CRCL / 78.49 mL/min) for CL (unitless)") # Table 3: CrCL on CL theta1 = 0.681 (RSE 20%; bootstrap median 0.686, 95% CI 0.360-1.002)

    # Inter-individual variability -- Cai 2022 Methods Eq. 1,
    # P_i = P * exp(eta_i). Table 3 footnote defines omega as the 'square root
    # of between-subject variability', so the reported 32.6% and 40.6% are
    # SD(eta) x 100; variances are 0.326^2 and 0.406^2. No correlation reported.
    etalcl ~ 0.106276 # Table 3: omega CL = 32.6% (RSE 12%), variance = 0.326^2
    etalvc ~ 0.164836 # Table 3: omega V = 40.6% (RSE 17%), variance = 0.406^2

    # Residual variability -- Cai 2022 Results 3.3: 'Proportional error model
    # was selected'; Methods Eq. 3, Y = F + F x EPS(1).
    propSd <- 0.379 ; label("Proportional residual error (fraction)") # Table 3: sigma pro = 37.9% (RSE 15%; bootstrap median 37.8%)
  })

  model({
    # Individual parameters -- Cai 2022 Eq. 1, 5 and 6.
    cl <- exp(lcl + etalcl) * (CRCL / 78.49)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination. Polymyxin B
    # is given as a 1-hour intravenous infusion into the central compartment.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
