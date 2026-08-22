Lee_2025_levofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous levofloxacin in",
    "healthy Korean adults (Lee 2025; n = 12; 84 plasma samples; single",
    "500 mg dose given as a 1-hour infusion). First-order disposition with",
    "a zero-order infusion input. Covariate effects: Cockcroft-Gault",
    "creatinine clearance on CL as a power function centred on the cohort",
    "median 105.71 mL/min (exponent 0.901), and James-formula lean body",
    "mass on the peripheral volume as a power function centred on the",
    "cohort median 47.91 kg (exponent 1.75). Inter-individual variability",
    "was estimated on CL and Q only; residual error is proportional."
  )
  reference <- paste(
    "Lee YJ, Kang G, Zang DY, Lee DH (2025).",
    "Development of a population pharmacokinetic model of levofloxacin in",
    "healthy adults and identification of optimal dosing regimens.",
    "Pharmaceuticals 18(5):621.",
    "doi:10.3390/ph18050621.",
    sep = " "
  )
  vignette <- "Lee_2025_levofloxacin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Plasma levofloxacin was the measured analyte (Methods 4.3 "Drug Assay":
    # LC-MS/MS on plasma, calibration range 0.05-50 mg/L).
    central     = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    # V2 is described only as "volume of distribution for the first peripheral
    # compartment" (Table 2 footnote); the source names no biological matrix
    # for it, so the specimen assignment is conventional, not source-verified.
    peripheral1 = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation,",
        "reported as raw mL/min and NOT normalised to 1.73 m^2 body surface",
        "area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on CL:",
        "CL = theta1 * (CRCL / 105.71)^theta2 with theta2 = 0.901 per",
        "Lee 2025 Table 2. The reference 105.71 mL/min is the cohort",
        "median stated in the Table 2 footnote ('105.71 mL/min used as",
        "the median reference value'); Table 1 reports the same median",
        "rounded to 106 mL/min (range 74.8-113). The Cockcroft-Gault",
        "form used is given in Table 1 footnote c:",
        "CrCl = (140 - Age) * TBW / (CR * 72), multiplied by 0.85 for",
        "female subjects. Note that the model was fit only over the",
        "narrow healthy range 74.8-113 mL/min; the paper's own Monte",
        "Carlo simulations (Methods 4.5) extrapolate it to 10-170",
        "mL/min."
      ),
      source_name        = "CrCl"
    ),
    LBM = list(
      description        = "Lean body mass derived from sex, total body weight, and height by the James equation.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on the peripheral volume:",
        "V2 = theta3 * (LBM / 47.91)^theta4 with theta4 = 1.75 per",
        "Lee 2025 Table 2. The reference 47.91 kg is the cohort median",
        "stated in the Table 2 footnote; Table 1 reports the same median",
        "rounded to 47.9 kg (range 37.7-60.3). The James formulae are",
        "given in Table 1 footnote a:",
        "LBM (female) = 1.07 * TBW - 148 * (TBW / height)^2 and",
        "LBM (male) = 1.1 * TBW - 128 * (TBW / height)^2, with TBW in kg",
        "and height in cm. Fat-free mass (Janmahasatian) and normal fat",
        "mass (Anderson) were also screened but not retained",
        "(Methods 4.4)."
      ),
      source_name        = "LBM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "29.0 to 44.0 years",
    age_median     = "35.5 years",
    weight_range   = "47.3 to 77.2 kg (total body weight)",
    weight_median  = "68.0 kg (total body weight)",
    sex_female_pct = 100 * 8 / 12,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Healthy adult volunteers aged 19-55 years with no clinically",
      "significant medical history, screened to exclude hepatic or renal",
      "disease. Serum protein, albumin, creatinine, and cystatin C were",
      "all within the normal range (Results 2.1). No adverse events",
      "occurred."
    ),
    dose_range     = paste(
      "Single 500 mg dose of levofloxacin dissolved in 100 mL normal",
      "saline and given as a 1-hour intravenous infusion by infusion pump",
      "(Methods 4.2). No multiple-dose data were collected."
    ),
    sampling       = paste(
      "Immediately before the start of the infusion (0 min) and at 61 min,",
      "75 min, 90 min, 4 h, 8 h, and 24 h after the start of the infusion.",
      "84 plasma samples were analysed (Results 2.2)."
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCl median 106 mL/min (range 74.8-113); MDRD eGFR",
      "median 80.1 mL/min/1.73 m^2 (69.8-121); CKD-EPI creatinine eGFR",
      "median 95.8 (83.1-120) (Table 1). No subject had impaired or",
      "augmented renal clearance."
    ),
    regions        = "Republic of Korea (Anyang; Hallym University Sacred Heart Hospital Clinical Trial Center).",
    notes          = paste(
      "Baseline demographics from Lee 2025 Table 1 and Results 2.1.",
      "The study ran August-September 2024 (IRB 2024-05-015, protocol",
      "LVF-001). The final model was fit in NONMEM 7.5 with FOCE-I",
      "(ADVAN3 TRANS4) and evaluated with PsN 5.5.0 bootstrap (2000",
      "samples) and VPC (1000 simulations). The unbound fraction used for",
      "the paper's fAUC/MIC target-attainment simulations was fixed at 0.7",
      "(Methods 4.5); protein binding is not part of the PK model itself",
      "and is therefore not encoded here."
    )
  )

  ini({
    # Structural parameters. Typical values refer to the covariate reference
    # subject: CrCl 105.71 mL/min and LBM 47.91 kg (Lee 2025 Table 2 footnote).

    lcl <- log(13.4); label("Clearance (CL, L/h) at CrCl 105.71 mL/min")               # Table 2, theta1 (RSE 3.36%)
    lvc <- log(34.3); label("Central volume of distribution (V1, L)")                  # Table 2, V1 (RSE 8.93%)
    lq  <- log(72.8); label("Intercompartmental clearance (Q, L/h)")                   # Table 2, Q (RSE 10.9%)
    lvp <- log(67.7); label("Peripheral volume of distribution (V2, L) at LBM 47.91 kg") # Table 2, theta3 (RSE 3.42%)

    # Covariate effects, both entered as power functions of the covariate
    # divided by its cohort median (Methods 4.4: "(individual value/median
    # value)^theta ... The median value of each covariate in the study
    # population was used as the reference value for scaling").
    e_crcl_cl <- 0.901; label("Power exponent of CrCl on CL (unitless)")    # Table 2, theta2 (RSE 16.8%)
    e_lbm_vp  <- 1.75;  label("Power exponent of LBM on V2 (unitless)")     # Table 2, theta4 (RSE 12.5%)

    # Inter-individual variability. Lee 2025 Table 2 reports IIV on CL and Q
    # only (no IIV on V1 or V2). The table gives bare percentages under
    # "Interindividual variability ... (%)" with no stated convention and the
    # paper prints no exponential-IIV equation and no OMEGA covariance block,
    # so the two candidate readings cannot be arbitrated from the source.
    # They are numerically almost identical here, so the NONMEM
    # approximate-CV convention omega = pct / 100 is used. Under the
    # alternative log-normal-CV reading omega = sqrt(log(1 + CV^2)):
    #   CL: 0.0899 here vs 0.0897 (0.2% lower)
    #   Q : 0.360  here vs 0.349  (3.1% lower)
    # The third possible reading -- that the percentage is a variance, giving
    # omega = sqrt(0.0899) = 0.300 on CL -- IS material, and it is refuted by
    # the paper's own Table 3: at that spread the published target-attainment
    # values are unreachable for any within-stratum CrCl distribution. The
    # vignette runs that check ("iiv-scale-check").
    etalcl ~ 0.008082  # IIV CL = 8.99% -> 0.0899^2 (Table 2; RSE 15.3%, shrinkage 3.58%)
    etalq  ~ 0.129600  # IIV Q  = 36.0% -> 0.360^2  (Table 2; RSE 30.6%, shrinkage 10.2%)

    # Residual error: proportional only (Table 2 "Residual variability").
    propSd <- 0.0699; label("Proportional residual error (fraction)")  # Table 2, proportional error 6.99% (RSE 13.8%, shrinkage 7.09%)
  })

  model({
    # Individual PK parameters. Covariates enter as power functions centred on
    # the cohort medians (Lee 2025 Table 2 structural-model rows:
    # CL = theta1 * (CrCl/105.71)^theta2 and V2 = theta3 * (LBM/47.91)^theta4).
    cl <- exp(lcl + etalcl) * (CRCL / 105.71)^e_crcl_cl
    vc <- exp(lvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp) * (LBM / 47.91)^e_lbm_vp

    # Levofloxacin is given intravenously; dose records target `central`
    # directly (a 1-hour infusion in the source study), so there is no depot.
    #
    # The transfer terms are written in MACRO form (q / vc, q / vp) rather than
    # as k12 / k21 micro-constants. rxSolve()'s default useLinCmt = TRUE
    # detector matches on canonical macro names; a two-compartment body written
    # with k12 / k21 is silently collapsed to a one-compartment closed form
    # (peripheral1 dropped, terminal half-life wrong) while total AUC still
    # equals Dose / CL, so the usual sanity check does not catch it. The macro
    # form keeps the default solve correct for every downstream user.
    d/dt(central)     <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # Concentration in mg/L (= ug/mL). Proportional residual error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
