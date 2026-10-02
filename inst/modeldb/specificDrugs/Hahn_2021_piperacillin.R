Hahn_2021_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for piperacillin in 26 critically ill",
    "Korean adults on venoarterial extracorporeal membrane oxygenation (VA-ECMO),",
    "13 of whom also received continuous venovenous hemodiafiltration (CVVHDF)",
    "(Hahn 2021). Zero-order IV input into the central compartment and",
    "first-order elimination. Clearance is a fractional ECMO reduction of the",
    "typical value plus an additive term linear in Cockcroft-Gault creatinine",
    "clearance centred on 54.7 mL/min; the central volume is increased",
    "2.46-fold during CVVHDF. Log-normal IIV on CL, V1 and V2 and a",
    "proportional residual error.",
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
    "(final piperacillin model equations for CL, V1, V2 and Q).",
    "All parameter estimates: Table 2, 'Final model' column.",
    "The tazobactam counterpart fitted in the same paper is",
    "modellib('Hahn_2021_tazobactam').",
    sep = " "
  )
  vignette <- "Hahn_2021_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Hahn 2021 puts one exponential eta on total CL (Table 2 'omega CL^2'), but
  # the typical CL is a sum of an ECMO-scaled intercept and an additive CrCL
  # term, so the eta multiplies the whole sum rather than pairing with lcl
  # alone. Declared here so checkModelConventions() recognises the pairing.
  paper_specific_etas <- "etalcl"

  # What each ODE state holds, in what amount units, in what biological matrix.
  # Verified against Hahn 2021 Methods, 'Plasma concentration assay' (total
  # piperacillin in plasma by LC-MS/MS) and Results, 'Population PK analysis'
  # (two-compartment model with first-order linear elimination).
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
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
        "Enters CL additively and linearly: CL = 9.4 * (1 - 0.092 * ECMO) +",
        "0.115 * (CrCL - 54.7) (Hahn 2021 Results equation; Table 2",
        "Theta_CrCL = 0.115 L/h per mL/min). 54.7 mL/min is the cohort median",
        "(Table 1; Methods: 'Continuous covariates were centered on the median",
        "population value'). Cohort median 54.7 mL/min, range 16.2-157; 40.5",
        "(18.0-111) while on CVVHDF and 63.2 (16.2-157) while off (Table 1).",
        "Stored under the canonical CRCL column following the raw",
        "Cockcroft-Gault precedent of Kim_2016_piperacillin.R and",
        "Conil_2010_tobramycin.R. Values below about 20 mL/min make the",
        "additive term large relative to the intercept (typical CL 4.1 L/h at",
        "16.2 mL/min on ECMO); the typical CL stays positive for CrCL >= 0."
      ),
      source_name = "CrCL"
    ),
    ECMO_STATUS = list(
      description = "Venoarterial ECMO support indicator (1 = on ECMO, 0 = weaned from ECMO)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (weaned from ECMO)",
      notes = paste(
        "Fractional reduction of the CL intercept: 9.4 * (1 - 0.092 * ECMO)",
        "(Table 2 Theta_ECMO = -0.092), i.e. typical CL 8.54 L/h on ECMO vs",
        "9.4 L/h after weaning at CrCL 54.7 mL/min (Discussion). Time-varying",
        "within subject: all 26 patients were sampled on ECMO (days 2-4) and",
        "the 14 who were weaned were sampled again on day 2 after ECMO",
        "discontinuation (Methods, 'Dosing and sampling procedure')."
      ),
      source_name = "ECMO"
    ),
    RRT_CRRT_STATUS = list(
      description = paste(
        "Continuous venovenous hemodiafiltration indicator (1 = on CVVHDF,",
        "0 = not on CVVHDF)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHDF)",
      notes = paste(
        "Linear multiplier on V1: V1 = 6.56 * (1 + 1.46 * CVVHDF) (Table 2",
        "Theta_CVVHDF = 1.46), i.e. 16.14 L on CVVHDF vs 6.56 L off",
        "(Discussion). CVVHDF was the only continuous RRT modality used",
        "(Prismaflex, AN69 membrane; Methods), so the modality-agnostic",
        "RRT_CRRT_STATUS column is used, following Wi_2017_teicoplanin.R from",
        "the same centre. Time-varying within subject: 67 samples were drawn",
        "on ECMO + CVVHDF and 27 after ECMO weaning while still on CVVHDF",
        "(Results, 'Study population')."
      ),
      source_name = "CVVHDF"
    )
  )

  # Screened in the covariate search (Hahn 2021 Methods, 'Population PK
  # analyses'; Table S1) but not retained in the final piperacillin model.
  # Documented for provenance; none is referenced in model(). Glomerular
  # filtration rate and CVVHDF duration were also screened; they have no
  # distinct canonical column here and are recorded in population$notes.
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
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 (Table S1, univariate dOFV -4.6) but not retained. Cohort median 26, range 18.1-31.3 (Table 1).",
      source_name = "BMI"
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
      notes = "Screened on CL (Table S1, univariate dOFV -8.3) but not retained. Cohort median 23.1 mg/dL, range 7-64.2 (Table 1).",
      source_name = "BUN"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Table S1, univariate dOFV -9.1) but not retained; CrCL was preferred. Cohort median 1.4 mg/dL, range 0.37-5.22 (Table 1).",
      source_name = "sCr"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Entered the full model on V1 but was removed in backward elimination",
        "(OFV increase 6.2 < 6.64; Table S1). Cohort median 2.1 mg/dL, range",
        "0.5-8.1 (Table 1)."
      ),
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
      "off both). Total piperacillin by LC-MS/MS, LLOQ 0.5 mg/L.",
      "Glomerular filtration rate and CVVHDF duration were also screened as",
      "covariates and not retained."
    )
  )

  ini({
    # ===== Structural parameters (Hahn 2021 Table 2, 'Final model') =====
    lcl <- log(9.4)
    label("CL intercept off ECMO at CrCL 54.7 mL/min (L/h)")
    # Table 2 Theta_CL = 9.4 L/h (RSE 7%; bootstrap median 9.601, 95% CI 8.307-11.43)
    lvc <- log(6.56)
    label("Central volume V1 off CVVHDF (L)")
    # Table 2 Theta_V1 = 6.56 L (RSE 18%; bootstrap median 7.824, 95% CI 3.038-14.52)
    lvp <- log(14.2)
    label("Peripheral volume V2 (L)")
    # Table 2 Theta_V2 = 14.2 L (RSE 14%; bootstrap median 12.70, 95% CI 9.095-17.04)
    lq <- log(17.2)
    label("Intercompartmental clearance Q (L/h)")
    # Table 2 Theta_Q = 17.2 L/h (RSE 22%; bootstrap median 14.53, 95% CI 7.318-23.92)

    # ===== Covariate effects (Hahn 2021 Table 2 and Results equations) =====
    e_ecmo_cl <- -0.092
    label("Fractional change in the CL intercept on ECMO (unitless)")
    # Table 2 Theta_ECMO = -0.092 (RSE 33%; bootstrap median -0.099, 95% CI -0.170 to -0.020)
    e_crcl_cl <- 0.115
    label("Additive CL slope per mL/min of CrCL above 54.7 mL/min (L/h per mL/min)")
    # Table 2 Theta_CrCL = 0.115 (RSE 17%; bootstrap median 0.120, 95% CI 0.051-0.166)
    e_crrt_vc <- 1.46
    label("Fractional increase in V1 on CVVHDF (unitless)")
    # Table 2 Theta_CVVHDF = 1.46 (RSE 40%; bootstrap median 1.837, 95% CI 0.076-5.257)

    # ===== IIV (Hahn 2021 Table 2, printed as omega^2 variances) =====
    # Methods: 'Interindividual variability (IIV) was modeled exponentially'.
    # No covariance is reported, so the etas are independent.
    etalcl ~ 0.0523 # Table 2 'omega_CL^2' = 0.0523 (RSE 33%, shrinkage 8%)
    etalvc ~ 0.291 # Table 2 'omega_V1^2' = 0.291 (RSE 65%, shrinkage 31%)
    etalvp ~ 0.138 # Table 2 'omega_V2^2' = 0.138 (RSE 34%, shrinkage 28%)

    # ===== Residual error (Hahn 2021 Table 2: proportional only) =====
    propSd <- 0.3129
    label("Proportional residual error (fraction)")
    # Table 2 'sigma^2 proportional' = 0.0979 (RSE 20%); SD = sqrt(0.0979) = 0.3129
  })

  model({
    # ----- Individual PK parameters (Hahn 2021 Results equations) -----
    # CL = 9.4 * (1 - 0.092 * ECMO) + [0.115 * (CrCL - 54.7)], with the
    # exponential eta applied to the whole typical value.
    cl <- (exp(lcl) * (1 + e_ecmo_cl * ECMO_STATUS) + e_crcl_cl * (CRCL - 54.7)) * exp(etalcl)
    # V1 = 6.56 * (1 + 1.46 * CVVHDF)
    vc <- exp(lvc + etalvc) * (1 + e_crrt_vc * RRT_CRRT_STATUS)
    # V2 = 14.2; Q = 17.2
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Piperacillin is given as a zero-order IV infusion into the central
    # compartment (rate or duration from the data); there is no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma piperacillin: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
