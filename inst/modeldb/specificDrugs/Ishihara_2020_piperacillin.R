Ishihara_2020_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for piperacillin in 18 Japanese late",
    "elderly (over 75 years) inpatients with pneumonia (Ishihara 2020), given",
    "piperacillin/tazobactam 4/0.5 g or 2/0.25 g as a 1-h IV infusion three",
    "times daily. Zero-order IV input into the central compartment,",
    "first-order elimination, a total clearance linear in Cockcroft-Gault",
    "creatinine clearance centred at the cohort median of 37.4 mL/min,",
    "exponential IIV on CL, Vc and Q (none on Vp) and a combined",
    "proportional-plus-additive residual error.",
    sep = " "
  )
  reference <- paste(
    "Ishihara N, Nishimura N, Ikawa K, Karino F, Miura K, Tamaki H, Yano T,",
    "Isobe T, Morikawa N, Naora K.",
    "Population pharmacokinetic modeling and pharmacodynamic target attainment",
    "simulation of piperacillin/tazobactam for dosing optimization in late",
    "elderly patients with pneumonia.",
    "Antibiotics (Basel). 2020;9(3):113.",
    "doi:10.3390/antibiotics9030113.",
    "All parameter estimates and the clearance equation: Table 2.",
    "The tazobactam model fitted in the same paper is",
    "modellib('Ishihara_2020_tazobactam').",
    sep = " "
  )
  vignette <- "Ishihara_2020_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Verified against Ishihara 2020 section 4.2 (total plasma
  # concentrations by HPLC-UV) and section 2.2.1 (two-compartment model).
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault formula",
        "(Ishihara 2020 Table 1 footnote a). Absolute clearance in mL/min, NOT",
        "normalised to 1.73 m^2 body surface area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters total clearance additively and linearly:",
        "CL (L/h) = 4.58 + 0.061 * (CLcr - 37.4) (Table 2), where 37.4 mL/min",
        "is the cohort median (Table 2 footnote). Cohort CLcr",
        "38.0 +/- 11.1 mL/min, range 21.5-59.1 (Table 1). The typical CL stays",
        "positive for any CLcr above -37.7 mL/min, so no guard is needed.",
        "Follows the additive centred-linear precedent of",
        "Conil_2010_tobramycin.R."
      ),
      source_name = "CLcr"
    )
  )

  # Screened in the forward-inclusion covariate search (section 2.2.1) but
  # not retained. Documented for provenance; none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at enrollment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Significantly affected CL in univariate screening but was not",
        "incorporated because of its high correlation with CLcr (section",
        "2.2.1). Cohort 86.5 +/- 6.0 years, range 75-101 (Table 1)."
      ),
      source_name = "age"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Significantly affected CL in univariate screening but was not",
        "incorporated because of its high correlation with CLcr; not",
        "significant on Vc, Q or Vp (section 2.2.1). Cohort 45.5 +/- 10.0 kg,",
        "range 32.0-68.7 (Table 1)."
      ),
      source_name = "body weight"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on Vc, Q and Vp (section 2.2.1), not retained. Cohort",
        "19.1 +/- 3.5, range 13.9-27.3 (Table 1)."
      ),
      source_name = "body mass index"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on Vc, Q and Vp (section 2.2.1), not retained. Cohort",
        "2.9 +/- 0.6 g/dL, range 2.1-3.7 (Table 1)."
      ),
      source_name = "serum albumin"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on Vc, Q and Vp (section 2.2.1), not retained. Cohort",
        "0.91 +/- 0.31 mg/dL, range 0.60-1.55 (Table 1)."
      ),
      source_name = "serum creatinine"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on Vc, Q and Vp (section 2.2.1), not retained. No summary",
        "statistic is reported in Table 1."
      ),
      source_name = "total protein"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18,
    n_studies = 1,
    age_range = "75-101 years",
    age_mean = "86.5 years (SD 6.0)",
    weight_range = "32.0-68.7 kg",
    weight_mean = "45.5 kg (SD 10.0)",
    sex_female_pct = 22.2,
    race_ethnicity = c(Asian = 100),
    disease_state = "pneumonia in late elderly (over 75 years) inpatients",
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 38.0 +/- 11.1 mL/min, range",
      "21.5-59.1 (median 37.4); serum creatinine 0.91 +/- 0.31 mg/dL, range",
      "0.60-1.55"
    ),
    dose_range = paste(
      "piperacillin/tazobactam 4/0.5 g, or 2/0.25 g when eGFR < 50 mL/min,",
      "each infused intravenously over 1 h three times daily"
    ),
    regions = "Japan (Shimane University Hospital)",
    notes = paste(
      "Baseline demographics: Table 1 (14 men, 4 women; height",
      "154.1 +/- 7.8 cm; body mass index 19.1 +/- 3.5; serum albumin",
      "2.9 +/- 0.6 g/dL). Patients with serious heart, liver or renal failure,",
      "suspected atypical pneumonia or beta-lactam allergy were excluded",
      "(section 4.1). 100 piperacillin concentrations were sampled at 0, 1,",
      "1.5, 2, 3 and 5 h after the start of an infusion (section 4.2)."
    )
  )

  ini({
    # ===== Structural parameters (Ishihara 2020 Table 2) =====
    # CL (L/h) = theta1 + theta2 * (CLcr - 37.4), additive linear centred at
    # the cohort median CLcr (Table 2 footnote: 'The median value of creatinine
    # clearance was 37.4').
    lcl <- log(4.58)
    label("Clearance at CLcr = 37.4 mL/min (L/h)")
    # Table 2 theta1 = 4.58 (SE 0.289; bootstrap 95% CI 4.04-5.22)
    e_crcl_cl <- 0.061
    label("Slope of CL per (CRCL - 37.4) (L/h per mL/min)")
    # Table 2 theta2 = 0.061 (SE 0.0293; bootstrap 95% CI 0.0107-0.114)
    lvc <- log(5.39)
    label("Central volume of distribution (L)")
    # Table 2 theta3 = 5.39 (SE 0.969; bootstrap 95% CI 3.81-8.33)
    lq <- log(20.7)
    label("Intercompartmental clearance (L/h)")
    # Table 2 theta4 = 20.7 (SE 6.11; bootstrap 95% CI 12.0-37.1)
    lvp <- log(6.96)
    label("Peripheral volume of distribution (L)")
    # Table 2 theta5 = 6.96 (SE 0.314; bootstrap 95% CI 5.09-7.84)

    # ===== IIV (Ishihara 2020 Table 2, variances) =====
    # Section 4.3: exponential IIV, theta_i = theta * exp(eta_i). Table 2 prints
    # omega^2 with the matching CV = sqrt(exp(omega^2) - 1) (0.0705 -> 27.0%),
    # so the values are variances. No covariance terms are reported. The IIV on
    # Vp was fixed to zero (section 2.2.1) and is therefore omitted.
    etalcl ~ 0.0705 # Table 2 'omega2 CL' = 0.0705 (CV 27.0%; SE 0.0177)
    etalvc ~ 0.389 # Table 2 'omega2 Vc' = 0.389 (CV 69.0%; SE 0.153)
    etalq ~ 0.311 # Table 2 'omega2 Q' = 0.311 (CV 60.4%; SE 0.327)

    # ===== Residual error (Ishihara 2020 Table 2; section 4.3) =====
    # Cobs = Cpred * (1 + eps_prop) + eps_add; Table 2 prints the sigma^2
    # values, so the SDs are their square roots.
    propSd <- 0.030447
    label("Proportional residual error (fraction)")
    # Table 2 'sigma2 proportional' = 0.000927 (SE 0.000798); sqrt = 0.030447
    addSd <- 5.00999
    label("Additive residual error (mg/L)")
    # Table 2 'sigma2 additive' = 25.1 (SE 12.8); sqrt = 5.00999
  })

  model({
    # ----- Individual PK parameters (Ishihara 2020 Table 2) -----
    # CL (L/h) = 4.58 + 0.061 * (CLcr - 37.4), with the exponential eta of
    # section 4.3 multiplying the covariate-adjusted typical value.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 37.4)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Piperacillin is given as a 1-h zero-order IV infusion into the central
    # compartment (rate or duration supplied through the data); no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma piperacillin: dose in mg, vc in L -> mg/L (= ug/mL, the
    # HPLC assay scale of section 4.2).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
