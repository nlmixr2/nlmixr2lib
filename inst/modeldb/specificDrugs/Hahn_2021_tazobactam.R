Hahn_2021_tazobactam <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam in 26 critically ill",
    "Korean adults on venoarterial extracorporeal membrane oxygenation (VA-ECMO),",
    "13 of whom also received continuous venovenous hemodiafiltration (CVVHDF)",
    "(Hahn 2021). Zero-order IV input into the central compartment and",
    "first-order elimination. Clearance is an exponential ECMO reduction of the",
    "typical value plus an additive term linear in Cockcroft-Gault creatinine",
    "clearance centred on 54.7 mL/min. Log-normal IIV on CL and V1 and a",
    "combined proportional + additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Hahn J, Min KL, Kang S, Yang S, Park MS, Wi J, Chang MJ.",
    "Population Pharmacokinetics and Dosing Optimization of",
    "Piperacillin-Tazobactam in Critically Ill Patients on Extracorporeal",
    "Membrane Oxygenation and the Influence of Concomitant Renal Replacement",
    "Therapy.",
    "Microbiol Spectr. 2021;9(3):e00633-21.",
    "doi:10.1128/spectrum.00633-21.",
    "Structural and covariate equations: Results, 'Population PK analysis'",
    "(final tazobactam model equations for CL, V1, V2 and Q).",
    "All parameter estimates: Table 3, 'Final model' column.",
    "The piperacillin counterpart fitted in the same paper is",
    "modellib('Hahn_2021_piperacillin').",
    sep = " "
  )
  vignette <- "Hahn_2021_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Hahn 2021 puts one exponential eta on total CL (Table 3 'omega CL^2'), but
  # the typical CL is a sum of an ECMO-scaled intercept and an additive CrCL
  # term, so the eta multiplies the whole sum rather than pairing with lcl
  # alone. Declared here so checkModelConventions() recognises the pairing.
  paper_specific_etas <- "etalcl"

  # What each ODE state holds, in what amount units, in what biological matrix.
  # Verified against Hahn 2021 Methods, 'Plasma concentration assay' (total
  # tazobactam in plasma by LC-MS/MS) and Results, 'Population PK analysis'
  # (two-compartment model with first-order linear elimination).
  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation",
        "(Hahn 2021 Table 1 footnote b and Results). Absolute clearance in",
        "mL/min, NOT normalised to 1.73 m^2 body surface area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL additively and linearly: CL = 7.93 * exp(-0.0723 * ECMO) +",
        "0.104 * (CrCL - 54.7) (Hahn 2021 Results equation; Table 3",
        "Theta_CrCL = 0.104 L/h per mL/min). 54.7 mL/min is the cohort median",
        "(Table 1; Methods: 'Continuous covariates were centered on the median",
        "population value'). Cohort median 54.7 mL/min, range 16.2-157 (Table",
        "1). Stored under the canonical CRCL column following the raw",
        "Cockcroft-Gault precedent of Kim_2016_tazobactam.R and",
        "Conil_2010_tobramycin.R. The typical CL stays positive for CrCL >= 0",
        "(1.7 L/h on ECMO at CrCL 0)."
      ),
      source_name = "CrCL"
    ),
    ECMO_STATUS = list(
      description = "Venoarterial ECMO support indicator (1 = on ECMO, 0 = weaned from ECMO)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (weaned from ECMO)",
      notes = paste(
        "Exponential multiplier on the CL intercept: 7.93 * exp(-0.0723 * ECMO)",
        "(Table 3 Theta_ECMO = -0.0723), i.e. a 7.0% lower intercept on ECMO.",
        "The printed piperacillin equation uses the linear form",
        "(1 + theta * ECMO) while the tazobactam equation uses the exponential",
        "form; each is encoded as printed. Time-varying within subject: all 26",
        "patients were sampled on ECMO (days 2-4) and the 14 who were weaned",
        "were sampled again on day 2 after ECMO discontinuation (Methods,",
        "'Dosing and sampling procedure')."
      ),
      source_name = "ECMO"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Cohort median 57 years, range 20-89 (Table 1).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Methods) but not retained. 19 of 26 (73.1%) male (Table 1).",
      source_name = "sex"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Cohort median 70 kg, range 40.8-92.5 (Table 1).",
      source_name = "weight"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous venovenous hemodiafiltration indicator (1 = on CVVHDF, 0 = not on CVVHDF)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHDF)",
      notes = paste(
        "Screened (Methods, 'use of CRRT') but not retained for tazobactam,",
        "unlike piperacillin where it enlarges V1 (Hahn_2021_piperacillin)."
      ),
      source_name = "CVVHDF"
    ),
    TPRO = list(
      description = "Total plasma protein",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Cohort median 4.9 g/dL, range 2.7-6.8 (Table 1).",
      source_name = "total plasma protein"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Table S2, univariate dOFV -5.3) but not retained. Cohort median 23.1 mg/dL, range 7-64.2 (Table 1).",
      source_name = "BUN"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Table S2, univariate dOFV -9) but not retained; CrCL was preferred. Cohort median 1.4 mg/dL, range 0.37-5.22 (Table 1).",
      source_name = "sCr"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Cohort median 2.1 mg/dL, range 0.5-8.1 (Table 1).",
      source_name = "T.bili"
    ),
    BFR = list(
      description = "CVVHDF blood flow rate",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods; Discussion) but not retained. Median 150 mL/min, range 100-160 (Table 1).",
      source_name = "blood flow rate"
    ),
    DFR = list(
      description = "CVVHDF dialysate flow rate",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods; Discussion) but not retained. Median 1200 mL/h, range 800-1400 (Table 1).",
      source_name = "dialysate flow rate"
    ),
    ECMO_PUMP_SPEED = list(
      description = "ECMO centrifugal pump speed",
      units = "RPM",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Not summarised in Table 1.",
      source_name = "ECMO pump speed"
    ),
    Q_ECMO = list(
      description = "ECMO circuit blood flow rate",
      units = "L/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Median 3.0 L/min, range 0.3-4.1 (Table 1).",
      source_name = "ECMO flow rate"
    ),
    T_ECMO = list(
      description = "Duration of ECMO support",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods) but not retained. Median 135 h, range 52.1-264 (Table 1).",
      source_name = "ECMO duration"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26,
    n_studies = 1,
    n_observations = 244,
    age_range = "20-89 years",
    age_median = "57 years",
    weight_range = "40.8-92.5 kg",
    weight_median = "70 kg",
    sex_female_pct = 26.9,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "critically ill adults on venoarterial ECMO in a cardiac intensive care",
      "unit (indications: ST-elevation MI 16, valvular heart disease 4,",
      "cardiomyopathy 3, non-ST-elevation MI 3) receiving",
      "piperacillin-tazobactam; 13 of 26 also on CVVHDF; APACHE II median 32",
      "(range 6-46)"
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance median 54.7 mL/min (range",
      "16.2-157); serum creatinine median 1.4 mg/dL (range 0.37-5.22)"
    ),
    dose_range = paste(
      "piperacillin-tazobactam 2/0.25 g q6h (9 patients), 3/0.375 g q6h (4),",
      "4/0.5 g q6h (9) or 4/0.5 g q8h (4), each infused intravenously over",
      "40 min"
    ),
    regions = "South Korea (Severance Cardiovascular Hospital, Seoul; single centre)",
    notes = paste(
      "Baseline demographics: Table 1. Prospective observational study,",
      "November 2015 to January 2019 (NCT02581280). Samples pre-dose and at",
      "0-0.5, 0.5-1, 1-2, 2-4, 4-6 and 6-8 h after a dose on ECMO days 2-4,",
      "repeated on day 2 after ECMO weaning in 14 patients; 244 samples in",
      "all (67 on ECMO + CVVHDF, 96 on ECMO only, 27 off ECMO on CVVHDF, 54",
      "off both). Total tazobactam by LC-MS/MS, LLOQ 0.5 mg/L.",
      "Glomerular filtration rate and CVVHDF duration were also screened as",
      "covariates and not retained. Table 3 lists no IIV on V2 or Q."
    )
  )

  ini({
    # ===== Structural parameters (Hahn 2021 Table 3, 'Final model') =====
    lcl <- log(7.93)
    label("CL intercept off ECMO at CrCL 54.7 mL/min (L/h)")
    # Table 3 Theta_CL = 7.93 L/h (RSE 6%; bootstrap median 7.71, 95% CI 6.50-8.91)
    lvc <- log(8.58)
    label("Central volume V1 (L)")
    # Table 3 Theta_V1 = 8.58 L (RSE 12%; bootstrap median 8.45, 95% CI 5.26-13.92)
    lvp <- log(10.5)
    label("Peripheral volume V2 (L)")
    # Table 3 Theta_V2 = 10.5 L (RSE 7%; bootstrap median 10.3, 95% CI 8.75-12.40)
    lq <- log(17.1)
    label("Intercompartmental clearance Q (L/h)")
    # Table 3 Theta_Q = 17.1 L/h (RSE 14%; bootstrap median 15.9, 95% CI 10.42-28.04)

    # ===== Covariate effects (Hahn 2021 Table 3 and Results equations) =====
    e_ecmo_cl <- -0.0723
    label("Exponential coefficient of ECMO on the CL intercept (unitless)")
    # Table 3 Theta_ECMO = -0.0723 (RSE 43%; bootstrap median -0.0722, 95% CI -0.175 to -0.010)
    e_crcl_cl <- 0.104
    label("Additive CL slope per mL/min of CrCL above 54.7 mL/min (L/h per mL/min)")
    # Table 3 Theta_CrCL = 0.104 (RSE 13%; bootstrap median 0.0979, 95% CI 0.018-0.121)

    # ===== IIV (Hahn 2021 Table 3, printed as omega^2 variances) =====
    # Results: 'IIV was included for CL and V1'; Methods: modeled exponentially.
    # No covariance is reported, so the etas are independent.
    etalcl ~ 0.0724 # Table 3 'omega_CL^2' = 0.0724 (RSE 40%, shrinkage 5%)
    etalvc ~ 0.705 # Table 3 'omega_V1^2' = 0.705 (RSE 61%, shrinkage 18%)

    # ===== Residual error (Hahn 2021 Table 3: combined) =====
    # Results: 'combined residual variability'. Table 3 prints separate
    # proportional and additive sigma^2 variances; with no control stream the
    # two are encoded as independent components (nlmixr2 combined2, the
    # NONMEM Y = F*(1 + EPS1) + EPS2 form with a diagonal $SIGMA).
    propSd <- 0.2598
    label("Proportional residual error (fraction)")
    # Table 3 'sigma^2 proportional' = 0.0675 (RSE 35%); SD = sqrt(0.0675) = 0.2598
    addSd <- 0.7190
    label("Additive residual error (mg/L)")
    # Table 3 'sigma^2 additional' = 0.517 (RSE 56%); SD = sqrt(0.517) = 0.7190 mg/L
  })

  model({
    # ----- Individual PK parameters (Hahn 2021 Results equations) -----
    # CL (L/h) = 7.93 * e^(-0.0723 * ECMO) + 0.104 * (CrCL - 54.7), with the
    # exponential eta applied to the whole typical value.
    cl <- (exp(lcl) * exp(e_ecmo_cl * ECMO_STATUS) + e_crcl_cl * (CRCL - 54.7)) * exp(etalcl)
    # V1 (L) = 8.58; V2 (L) = 10.5; Q (L/h) = 17.1
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Tazobactam is given as a zero-order IV infusion into the central
    # compartment (rate or duration from the data); there is no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma tazobactam: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
