Kim_2022_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for piperacillin in 38 critically ill",
    "Korean adults (19 on extracorporeal membrane oxygenation, ECMO) given",
    "piperacillin/tazobactam as 30-min IV infusions every 6 or 8 h (Kim 2022);",
    "zero-order IV input into the central compartment, first-order elimination,",
    "fixed allometric weight scaling (exponent 0.75 on CL and Q, 1 on VC and VP,",
    "reference 70 kg), an exponential effect of cystatin-C CKD-EPI eGFR on CL,",
    "and separate typical central volumes (each with its own IIV) for patients",
    "on and off ECMO.",
    sep = " "
  )
  reference <- paste(
    "Kim YK, Kim HS, Park S, Kim HI, Lee SH, Lee DH.",
    "Population pharmacokinetics of piperacillin/tazobactam in critically ill",
    "Korean patients and the effects of extracorporeal membrane oxygenation.",
    "J Antimicrob Chemother. 2022;77(5):1353-1364.",
    "doi:10.1093/jac/dkac059.",
    "Parameter estimates: Table 2. Structural equations: NONMEM control stream in",
    "the Supplementary data (ADVAN3 TRANS4).",
    sep = " "
  )
  vignette <- "Kim_2022_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Kim 2022 Methods: total plasma piperacillin assayed by
  # LC-MS/MS; two-compartment structural model (Results, Table 2).
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling (Kim 2022 Results; Supplementary NONMEM code):",
        "(WT / 70)^0.75 on CL and Q, (WT / 70)^1 on VC (both ECMO strata) and VP.",
        "Cohort median 60 kg (IQR 50-70); ECMO 70 kg, non-ECMO 54 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate from the CKD-EPI cystatin C",
        "equation, BSA-normalised (mL/min/1.73 m^2); the source paper's 'CECYS'"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Exponential effect on CL: exp(theta2 * (CECYS - 52.77)) with theta2 =",
        "0.00932 per mL/min/1.73 m^2 (Kim 2022 Table 2). 52.77 is the centring",
        "value printed in the Table 2 equation and the control stream (the",
        "Table 1 cohort median CKD-EPI CYS is 52.8, IQR 38.7-81.0). Selected over",
        "Cockcroft-Gault, MDRD, modified MDRD, creatinine CKD-EPI and the",
        "BSA-de-normalised CKD-EPI variants (Methods). The PTA simulations sampled",
        "it uniformly over 0-170 mL/min/1.73 m^2."
      ),
      source_name = "CECYS"
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation support (1 = on ECMO, 0 = not on ECMO)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the central volume: V1 = ECMO * VE1 + (1 - ECMO) * VE0",
        "(Supplementary NONMEM code), with a separate typical value and a",
        "separate IIV for each group (Table 2). 18 of the 19 ECMO patients were on",
        "veno-arterial ECMO (Table 1); ECMO type, pump speed and flow rate were",
        "tested but not retained (Methods)."
      ),
      source_name = "ECMO"
    )
  )

  # Screened in the stepwise covariate search (Kim 2022 Methods) but not
  # retained in the final model.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Tested on CL and the other PK parameters; not retained.",
      source_name = "sex"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and the other PK parameters; not retained.",
      source_name = "age"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and the other PK parameters; not retained.",
      source_name = "height"
    ),
    BSA = list(
      description = "Body surface area (Du Bois)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and the other PK parameters; not retained.",
      source_name = "BSA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported in g/dL in Kim 2022 Table 1; tested on CL and the other PK parameters; not retained.",
      source_name = "serum albumin level"
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported in g/dL in Kim 2022 Table 1; tested on CL and the other PK parameters; not retained.",
      source_name = "serum protein level"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL; not retained (the cystatin-C eGFR was).",
      source_name = "serum creatinine level"
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL as a raw level; not retained (it enters only through the CKD-EPI cystatin C eGFR).",
      source_name = "serum cystatin C level"
    ),
    ECMO_PUMP_SPEED = list(
      description = "ECMO centrifugal-pump speed",
      units = "RPM (revolutions per minute)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on the ECMO-affected parameter (VC); not retained.",
      source_name = "rpm"
    ),
    Q_ECMO = list(
      description = "ECMO circuit blood flow rate",
      units = "L/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on the ECMO-affected parameter (VC); not retained.",
      source_name = "flow rate"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38,
    n_studies = 1,
    n_observations = 226,
    age_range = "Median 66.5 years (IQR 53.3-78.8); ECMO 58 (46.5-64.5), non-ECMO 79 (67.5-83)",
    weight_range = "Median 60 kg (IQR 50-70); ECMO 70 (55-72.4), non-ECMO 54 (46-61)",
    sex_female_pct = 34.2,
    race_ethnicity = "Korean (all participants)",
    disease_state = paste(
      "Critically ill adults in the ICU receiving piperacillin/tazobactam for",
      "nosocomial infection, empirical management of septic shock, or",
      "prophylaxis during ECMO; 19 on ECMO (18 veno-arterial, 1 veno-venous)",
      "and 19 not on ECMO"
    ),
    dose_range = paste(
      "Piperacillin/tazobactam 2000/250, 3000/375 or 4000/500 mg as 30-min IV",
      "infusions every 6 or 8 h"
    ),
    regions = "Republic of Korea (Hallym University Sacred Heart Hospital, Anyang)",
    renal_function = paste(
      "CKD-EPI cystatin C eGFR median 52.8 mL/min/1.73 m^2 (IQR 38.7-81.0);",
      "8 patients on CRRT"
    ),
    notes = paste(
      "Prospective study, September 2020 to April 2021. Six samples per patient",
      "over the first dosing interval after enrolment (pre-dose and 0.5, 1, 2,",
      "3 [q6h] or 4 [q8h], and 6 or 8 h). Demographics: Kim 2022 Table 1."
    )
  )

  ini({
    # ===== Structural PK (Kim 2022 Table 2, final piperacillin model) =====
    # Reference subject: WT = 70 kg, CECYS = 52.77 mL/min/1.73 m^2.
    lcl <- log(5.05); label("Typical CL at WT = 70 kg and CRCL = 52.77 mL/min/1.73 m^2 (L/h)") # Table 2: theta1 = 5.05 L/h (RSE 6.67%)
    lvc_ecmo <- log(7.38); label("Typical central volume, patients on ECMO, at WT = 70 kg (L)") # Table 2: VC_ECMO = 7.38 L (RSE 15.0%)
    lvc_nonecmo <- log(16.5); label("Typical central volume, patients not on ECMO, at WT = 70 kg (L)") # Table 2: VC_nonECMO = 16.5 L (RSE 18.0%)
    lq <- log(6.28); label("Typical intercompartmental clearance at WT = 70 kg (L/h)") # Table 2: Q = 6.28 L/h (RSE 23.9%)
    lvp <- log(6.27); label("Typical peripheral volume at WT = 70 kg (L)") # Table 2: VP = 6.27 L (RSE 16.2%)

    # ===== Covariate effects =====
    e_crcl_cl <- 0.00932; label("Exponential coefficient on (CRCL - 52.77) for CL (per mL/min/1.73 m^2)") # Table 2: theta2 = 0.00932 (RSE 13.5%)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless)") # Results: k = 0.75 for clearance terms; control stream (WT/70)**0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/70) for VC and VP (unitless)") # Results: k = 1 for volume terms; control stream (WT/70)

    # ===== IIV (Kim 2022 Table 2) =====
    # theta_i = theta * exp(eta_i) (Methods). Table 2 prints IIV as %CV only;
    # converted with omega^2 = log(1 + CV^2). The control stream gives each
    # ECMO stratum its own eta on VC (ETA(2) for ECMO, ETA(3) for non-ECMO); no
    # IIV on Q or VP. No covariance terms are reported.
    etalcl ~ 0.107570 # Table 2: IIV CL 33.7% (RSE 11.8%)
    etalvc_ecmo ~ 0.191183 # Table 2: IIV VC_ECMO 45.9% (RSE 17.4%)
    etalvc_nonecmo ~ 0.355160 # Table 2: IIV VC_nonECMO 65.3% (RSE 31.5%)

    # ===== Residual error (Kim 2022 Table 2: proportional only) =====
    # Control stream: W = SQRT(THETA(5)**2 * IPRED**2) with $SIGMA 1 FIX, so
    # THETA(5) is the proportional SD.
    propSd <- 0.269; label("Proportional residual error (fraction)") # Table 2: proportional error 26.9% (RSE 11.3%)
  })

  model({
    # ----- Individual PK parameters (Supplementary NONMEM $PK) -----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * exp(e_crcl_cl * (CRCL - 52.77))
    vc_ecmo <- exp(lvc_ecmo + etalvc_ecmo) * (WT / 70)^e_wt_vc_vp
    vc_nonecmo <- exp(lvc_nonecmo + etalvc_nonecmo) * (WT / 70)^e_wt_vc_vp
    vc <- ECMO_STATUS * vc_ecmo + (1 - ECMO_STATUS) * vc_nonecmo
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Zero-order IV infusion into the central compartment (rate supplied
    # through the data-level RATE / DUR column).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma piperacillin: dose in mg, vc in L -> mg/L. The free fraction
    # of 0.91 used in the fT>MIC simulations is applied in the vignette.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
