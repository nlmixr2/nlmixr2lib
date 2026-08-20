Lee_2025_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for piperacillin in 12 healthy Korean adults",
    "given a single 4 g / 0.5 g piperacillin-tazobactam IV infusion over 30 min",
    "(Lee 2025); zero-order IV input into the central compartment, first-order",
    "elimination, a power effect of BSA-adjusted CKD-EPI creatinine eGFR on CL, a power",
    "effect of lean body mass on Q (with the typical Q fixed), and an exponential effect",
    "of body weight on the peripheral volume V2.",
    sep = " "
  )
  reference <- paste(
    "Lee YJ, Kang G, Zang DY, Lee DH.",
    "Population pharmacokinetic modeling of piperacillin/tazobactam in healthy adults",
    "and exploration of optimal dosing strategies.",
    "Pharmaceuticals (Basel). 2025;18(8):1124.",
    "doi:10.3390/ph18081124.",
    "Structural model, covariate equations and all parameter estimates: Table 2.",
    "Covariate-selection path: Table A1.",
    sep = " "
  )
  vignette <- "Lee_2025_piperacillin_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Lee 2025 Section 4.2 (plasma samples
  # assayed by LC-MS/MS) and Section 2.2 (two-compartment structural model).
  compartmentData <- list(
    central     = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the 2021 CKD-EPI creatinine equation, adjusted to the individual's body surface area (the source paper's 'CE' / 'BSA adjusted eGFR CKD-EPI_CR')",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL: (CRCL / 108.25)^theta2 with theta2 = 1.16",
        "(Lee 2025 Table 2). The reference 108.25 mL/min is the cohort median",
        "(Table 1 reports 108 mL/min, range 86.2-136). NOT normalised to",
        "1.73 m^2: Lee 2025 Table 1 footnote c computes the standard CKD-EPI",
        "value and then multiplies by BSA / 1.73, so the covariate is an",
        "absolute individual-BSA eGFR in mL/min. Removing this covariate raised",
        "the OFV by 16.414 and inflated IIV on CL from 7.17% to 13.2%",
        "(Section 2.2, Table A1 backward step 1)."
      ),
      source_name        = "CE (CKD-EPI_CR eGFR, BSA adjusted)"
    ),
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on Q: (LBM / 50.08)^theta5 with theta5 = 2.50",
        "(Lee 2025 Table 2). The reference 50.08 kg is the cohort median",
        "(Table 1 reports 50.1 kg, range 36.6-65.9). The typical Q (theta4 =",
        "4.32 L/h) was FIXED, so only the exponent was estimated. The paper",
        "flags this LBM-on-Q association as atypical but statistically robust",
        "(delta OFV = 14.798 on backward elimination; Table A1) and explicitly",
        "cautions that it needs external validation (Discussion). Lee 2025 does",
        "not report the formula used to compute LBM."
      ),
      source_name        = "LBM"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential effect on the peripheral volume V2:",
        "V2 = theta6 * exp(theta7 * (WT - 61.7)) with theta7 = 0.0288 per kg",
        "(Lee 2025 Table 2). The centering constant 61.7 kg is the cohort median",
        "weight (Table 1). Lee 2025 tested WT on V1 as well and found no",
        "significant improvement (Section 2.2), so no weight effect is applied",
        "to the central volume. Note the exponential (not power) form: the",
        "exponential form beat the power form by delta OFV = -0.619 at forward",
        "step 3 (Table A1)."
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
    # ===== Structural PK (Lee 2025 Table 2, final piperacillin model) =====
    # Reference subject: CE = 108.25 mL/min, LBM = 50.08 kg, WT = 61.7 kg
    # (the cohort medians used as centering constants in Table 2).
    lcl <- log(11.2); label("Typical CL at CRCL = 108.25 mL/min (L/h)")  # Table 2: theta1 = 11.2 L/h (RSE 3.40%; bootstrap 11.2, 95% CI 10.5-12.1)
    lvc <- log(6.24); label("Typical central volume V1 (L)")             # Table 2: theta3 = 6.24 L (RSE 8.99%; bootstrap 6.19, 95% CI 5.27-7.57)
    # theta4 carries a footnote "a, fixed" in Table 2 (no RSE, no bootstrap CI).
    lq  <- fixed(log(4.32)); label("Typical intercompartmental clearance Q at LBM = 50.08 kg (L/h)")  # Table 2: theta4 = 4.32 L/h, FIXED
    lvp <- log(2.59); label("Typical peripheral volume V2 at WT = 61.7 kg (L)")  # Table 2: theta6 = 2.59 L (RSE 3.11%; bootstrap 2.59, 95% CI 2.28-2.74)

    # ===== Covariate effects (Lee 2025 Table 2 structural-model equations) =====
    # CL  = theta1 * (CE / 108.25)^theta2
    # Q   = theta4 * (LBM / 50.08)^theta5
    # V2  = theta6 * exp(theta7 * (WT - 61.7))
    e_crcl_cl <- 1.16;   label("Power exponent on (CRCL / 108.25) for CL (unitless)")           # Table 2: theta2 = 1.16 (RSE 13.1%; bootstrap 1.15, 95% CI 0.811-1.59)
    e_lbm_q   <- 2.50;   label("Power exponent on (LBM / 50.08) for Q (unitless)")              # Table 2: theta5 = 2.50 (RSE 13.9%; bootstrap 2.45, 95% CI 1.39-3.56)
    e_wt_vp   <- 0.0288; label("Exponential coefficient on (WT - 61.7) for V2 (per kg)")        # Table 2: theta7 = 0.0288 (RSE 8.38%; bootstrap 0.0284, 95% CI 0.0208-0.0371)

    # Reference covariate values (Lee 2025 Table 2 equations; cohort medians in Table 1)
    crcl_ref <- 108.25; label("Reference BSA-adjusted CKD-EPI creatinine eGFR (mL/min, cohort median)")
    lbm_ref  <- 50.08;  label("Reference lean body mass (kg, cohort median)")
    wt_ref   <- 61.7;   label("Reference body weight (kg, cohort median)")

    # ===== IIV (Lee 2025 Table 2) =====
    # Lee 2025 Section 4.3: theta_i = theta * exp(eta_i), eta ~ N(0, omega^2).
    # Table 2 reports IIV as a percentage on the SD scale (each %RSE is below
    # the sqrt(2/n) = 40.8% floor that a variance-scale column could not beat
    # for n = 12); converted here with omega^2 = log(1 + CV^2):
    #   7.17% -> log(1 + 0.0717^2) = 0.00512772
    #  18.4%  -> log(1 + 0.184^2)  = 0.03329550
    # Lee 2025 Section 4.3 states that covariance (block) terms between CL, V1,
    # V2 and Q were tested and did not improve the fit, so the etas are
    # independent.
    etalcl ~ 0.00512772  # Table 2: IIV CL 7.17% (RSE 30.3%, shrinkage 18.7%; bootstrap 6.08, 95% CI 0-10.5)
    etalvc ~ 0.03329550  # Table 2: IIV V1 18.4% (RSE 28.7%, shrinkage 19.1%; bootstrap 17.5, 95% CI 0-29.9)

    # ===== Residual error (Lee 2025 Table 2: proportional only) =====
    propSd <- 0.134; label("Proportional residual error (fraction)")  # Table 2: proportional error 13.4% (RSE 12.2%, shrinkage 9.48%; bootstrap 13.1, 95% CI 9.39-16.0)
  })

  model({
    # ----- Individual PK parameters (Lee 2025 Table 2 structural equations) -----
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq) * (LBM / lbm_ref)^e_lbm_q
    vp <- exp(lvp + e_wt_vp * (WT - wt_ref))

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Piperacillin given as a zero-order IV infusion into the central
    # compartment (rate supplied through the data-level RATE / DUR column);
    # no depot compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total (protein-bound plus unbound) plasma piperacillin concentration:
    # dose in mg, vc in L -> mg/L. Lee 2025 assumed a fixed unbound fraction of
    # 0.7 when converting to free concentrations for the fT>MIC target-attainment
    # simulations (Discussion, limitation 3); that factor is a property of the
    # PD analysis, not of this PK model, so it is applied in the vignette rather
    # than here.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
