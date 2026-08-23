Lee_2025_tazobactam <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam in 12 healthy Korean adults",
    "given a single 4 g / 0.5 g piperacillin-tazobactam IV infusion over 30 min",
    "(Lee 2025); zero-order IV input into the central compartment, first-order",
    "elimination, a power effect of BSA-adjusted CKD-EPI creatinine eGFR on CL, a fixed",
    "typical intercompartmental clearance Q, and an exponential effect of body weight on",
    "the peripheral volume V2.",
    sep = " "
  )
  reference <- paste(
    "Lee YJ, Kang G, Zang DY, Lee DH.",
    "Population pharmacokinetic modeling of piperacillin/tazobactam in healthy adults",
    "and exploration of optimal dosing strategies.",
    "Pharmaceuticals (Basel). 2025;18(8):1124.",
    "doi:10.3390/ph18081124.",
    "Structural model, covariate equations and all parameter estimates: Table 3.",
    "Covariate-selection path: Table A2.",
    sep = " "
  )
  vignette <- "Lee_2025_piperacillin_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Lee 2025 Section 4.2 (plasma samples
  # assayed by LC-MS/MS) and Section 2.2 (two-compartment structural model).
  compartmentData <- list(
    central     = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the 2021 CKD-EPI creatinine equation, adjusted to the individual's body surface area (the source paper's 'CE' / 'BSA adjusted eGFR CKD-EPI_CR')",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL: (CRCL / 108.25)^theta2 with theta2 = 0.857",
        "(Lee 2025 Table 3). The reference 108.25 mL/min is the cohort median",
        "(Table 1 reports 108 mL/min, range 86.2-136). NOT normalised to",
        "1.73 m^2: Lee 2025 Table 1 footnote c computes the standard CKD-EPI",
        "value and then multiplies by BSA / 1.73, so the covariate is an",
        "absolute individual-BSA eGFR in mL/min. Section 2.2 describes this",
        "covariate as 'CrCl' in the tazobactam sentence, but Table 3 and",
        "Table A2 (backward step 2, delta OFV = 13.318) both name the CKD-EPI",
        "creatinine eGFR; the printed equation is taken as authoritative.",
        "Removing it inflated IIV on CL from 17.6% to 6.95%."
      ),
      source_name        = "CE (CKD-EPI_CR eGFR, BSA adjusted)"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential effect on the peripheral volume V2:",
        "V2 = theta5 * exp(theta6 * (WT - 61.7)) with theta6 = 0.0145 per kg",
        "(Lee 2025 Table 3). The centering constant 61.7 kg is the cohort median",
        "weight (Table 1). Lee 2025 tested WT on V1 as well and found no",
        "significant improvement (Section 2.2). Note the exponential (not power)",
        "form: it beat the power form by delta OFV = -0.539 at forward step 3",
        "(Table A2)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "26-50 years (inclusion criteria 19-55 years)",
    age_median     = "36.0 years",
    weight_range   = "45.8-88.5 kg",
    weight_median  = "61.7 kg",
    sex_female_pct = 33.3,
    race_ethnicity = "Korean (all participants)",
    disease_state  = "Healthy adults with no congenital or chronic health conditions; all baseline laboratory values within normal clinical ranges",
    dose_range     = "Single 4 g piperacillin / 0.5 g tazobactam intravenous dose in 100 mL saline, infused over 30 min",
    regions        = "Republic of Korea (Clinical Trial Center, Hallym University Sacred Heart Hospital, Anyang)",
    renal_function = "Normal; CrCl (Cockcroft-Gault) 105 mL/min (76.2-146), BSA-adjusted CKD-EPI creatinine eGFR 108 mL/min (86.2-136)",
    height_range   = "158-182 cm (median 168)",
    lbm_range      = "36.6-65.9 kg (median 50.1)",
    bsa_range      = "1.44-2.07 m^2 (median 1.71)",
    notes          = paste(
      "12 healthy Korean adults (8 male, 4 female) studied in January 2023;",
      "IRB 2022-08-006, trial registration KCT0009855. Rich sampling: pre-dose",
      "and 0.5, 0.75, 1, 2, 3 and 6 h after the start of the infusion, giving 84",
      "plasma samples across both analytes' models. Demographics in Lee 2025",
      "Table 1; study design in Section 4.2."
    )
  )

  ini({
    # ===== Structural PK (Lee 2025 Table 3, final tazobactam model) =====
    # Reference subject: CE = 108.25 mL/min, WT = 61.7 kg (cohort medians used
    # as centering constants in Table 3).
    lcl <- log(12.4); label("Typical CL at CRCL = 108.25 mL/min (L/h)")  # Table 3: theta1 = 12.4 L/h (RSE 3.26%; bootstrap 12.3, 95% CI 11.6-13.3)
    lvc <- log(9.03); label("Typical central volume V1 (L)")             # Table 3: theta3 = 9.03 L (RSE 6.40%; bootstrap 9.02, 95% CI 8.05-10.4)
    # theta4 carries a footnote "a, fixed" in Table 3 (no RSE, no bootstrap CI).
    # Unlike the piperacillin model, no covariate was retained on Q.
    lq  <- fixed(log(4.39)); label("Typical intercompartmental clearance Q (L/h)")  # Table 3: theta4 = 4.39 L/h, FIXED
    lvp <- log(3.21); label("Typical peripheral volume V2 at WT = 61.7 kg (L)")     # Table 3: theta5 = 3.21 L (RSE 5.48%; bootstrap 3.23, 95% CI 2.68-3.48)

    # ===== Covariate effects (Lee 2025 Table 3 structural-model equations) =====
    # CL = theta1 * (CE / 108.25)^theta2
    # V2 = theta5 * exp(theta6 * (WT - 61.7))
    e_crcl_cl <- 0.857;  label("Power exponent on (CRCL / 108.25) for CL (unitless)")     # Table 3: theta2 = 0.857 (RSE 13.1%; bootstrap 0.858, 95% CI 0.602-1.21)
    e_wt_vp   <- 0.0145; label("Exponential coefficient on (WT - 61.7) for V2 (per kg)")  # Table 3: theta6 = 0.0145 (RSE 16.9%; bootstrap 0.0142, 95% CI 0.0106-0.0238)

    # Reference covariate values (Lee 2025 Table 3 equations; cohort medians in Table 1)
    crcl_ref <- 108.25; label("Reference BSA-adjusted CKD-EPI creatinine eGFR (mL/min, cohort median)")
    wt_ref   <- 61.7;   label("Reference body weight (kg, cohort median)")

    # ===== IIV (Lee 2025 Table 3) =====
    # Lee 2025 Section 4.3: theta_i = theta * exp(eta_i), eta ~ N(0, omega^2).
    # Table 3 reports IIV as a percentage on the SD scale (the 29.0% RSE is
    # below the sqrt(2/n) = 40.8% floor a variance-scale column could not beat
    # for n = 12); converted here with omega^2 = log(1 + CV^2):
    #   6.95% -> log(1 + 0.0695^2) = 0.00481862
    # IIV on V1 was tested but dropped from the final tazobactam model for poor
    # precision (RSE > 70%) and negligible fit improvement (Section 2.2), so no
    # etalvc term is present here. Covariance (block) terms were tested and did
    # not improve the fit (Section 4.3).
    etalcl ~ 0.00481862  # Table 3: IIV CL 6.95% (RSE 29.0%, shrinkage 7.37%; bootstrap 6.17, 95% CI 0.403-9.94)

    # ===== Residual error (Lee 2025 Table 3: proportional only) =====
    propSd <- 0.135; label("Proportional residual error (fraction)")  # Table 3: proportional error 13.5% (RSE 9.57%, shrinkage 6.06%; bootstrap 13.2, 95% CI 10.4-15.6)
  })

  model({
    # ----- Individual PK parameters (Lee 2025 Table 3 structural equations) -----
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp + e_wt_vp * (WT - wt_ref))

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Tazobactam given as a zero-order IV infusion into the central compartment
    # (rate supplied through the data-level RATE / DUR column); no depot
    # compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma tazobactam concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
