Tsai_2024_vancomycin <- function() {
  description <- "Two-compartment IV-infusion population PK model for vancomycin in critically ill adults in a cardiac ICU, with and without temporary mechanical circulatory support (ECMO or extracorporeal VAD) or continuous renal replacement therapy (Tsai 2024). Separate typical clearances off and on CRRT sharing one eta, with serum creatinine on the off-CRRT arm and body weight on both; CRRT expands the central volume and temporary mechanical circulatory support lowers intercompartmental clearance. All three status covariates are time-varying."
  reference <- paste(
    "Tsai M-T, Wang W-C, Roan J-N, Luo C-Y, Chou C-H.",
    "Population Pharmacokinetics of Vancomycin in Intensive Care Patients with",
    "the Time-Varying Status of Temporary Mechanical Circulatory Support or",
    "Continuous Renal Replacement Therapy.",
    "Infect Dis Ther. 2024;13(12):2617-2635.",
    "doi:10.1007/s40121-024-01071-5.",
    sep = " "
  )
  vignette <- "Tsai_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column BW. Power form (WT / 70)^0.73 on both clearance arms (Table 2). Recorded time-varying from the first dose to the index dose, but essentially constant in this cohort (median change 0.0% [0.0, 0.0], Table 1).",
      source_name = "BW"
    ),
    CREAT = list(
      description = "Serum creatinine (time-varying)",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column Cr. Enters UNCENTRED as CREAT^-0.69 on the off-CRRT clearance arm only (Table 2: CL = theta3 x (Cr)^theta5 x (BW/70)^theta6; Supplementary Table S1 Model 3, 'Cr as a covariate only on CL for non-dialysis'), so theta3 = 2.85 L/h is the clearance at Cr = 1 mg/dL and 70 kg. The Results prose ('an independent CL for CRRT with serum Cr and BW as covariates') swaps the two arms; the printed equation, Supplementary Table S1 and the Discussion ('Cr level on decreasing the CL for non-CRRT') agree and are followed. Must be a positive number on every record, including CRRT records where it does not affect the solution. Cohort median 1.1 mg/dL (IQR 0.8-1.7) at the index dose (Table 1).",
      source_name = "Cr"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Continuous renal replacement therapy running (1) or not (0), time-varying",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source label CRRT (CVVH or CVVHDF; median effluent flow 29.7 mL/kg/h, Table 1). Selects the clearance arm (theta3 x Cr^theta5 off CRRT, theta4 on CRRT, both x (BW/70)^theta6) and multiplies Vc by exp(0.54). Time-varying: one patient crossed over from CRRT to non-CRRT within the vancomycin course (Results, Patient Characteristics). No patient received intermittent hemodialysis. CRRT was connected to the post-oxygenator ECMO circuit when both were running.",
      source_name = "CRRT"
    ),
    MCS_TEMP_STATUS = list(
      description = "Temporary mechanical circulatory support (VA-ECMO, VV-ECMO or extracorporeal VAD) running (1) or not (0), time-varying",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source label tMCS, defined as 'temporary mechanical circulatory support, including ECMO or extracorporeal VAD' (Supplementary Table S2 footnote). Index-dose cohort: 11 VA-ECMO, 1 VV-ECMO, 1 BiVAD and 1 LVAD (Levitronix CentriMag) of 25 (Table 1); the modes were pooled because separate ECMO / VAD covariates did not improve the model (Discussion, limitations). Multiplies Q by exp(-0.66) = 0.52. Time-varying: five patients started vancomycin on VA-ECMO that was explanted before PK sampling.",
      source_name = "tMCS"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    n_observations = 184L,
    age_range = "median 61.0 years (IQR 49.0-66.0); adults >= 18 years",
    weight_range = "median 68.8 kg (IQR 61.0-82.1)",
    sex_female_pct = 4,
    race_ethnicity = "Not reported (single-centre Taiwanese cohort)",
    disease_state = "Critically ill adults in a cardiac ICU receiving vancomycin (prophylactic, empirical or therapeutic), with or without temporary mechanical circulatory support (VA-ECMO, VV-ECMO or extracorporeal VAD) and with or without CRRT. Indications for support: myocardial infarction, cardiogenic shock, post-cardiotomy, ARDS. Median APACHE II 24.0.",
    dose_range = "Index dose median 900 mg (IQR 750-1000) infused over median 1.8 h (IQR 1.2-2.2); any dose and interval per clinical requirement",
    regions = "Taiwan (National Cheng Kung University Hospital, Tainan; January 2018 - March 2022)",
    renal_function = "Serum creatinine median 1.1 mg/dL (IQR 0.8-1.7); CrCL (Cockcroft-Gault) median 67 mL/min; 6/25 on CRRT at the index dose (7/25 at the first dose); no hemodialysis",
    mcs_status = "At the index dose: 11 none, 11 VA-ECMO, 1 VV-ECMO, 1 BiVAD, 1 LVAD (Table 1)",
    notes = "Prospective single-centre PK study; rich sampling (pre-dose and 0.5, 1, 2, 3, 6, 9 and 12 h after the end of infusion) after the index dose; HPLC-UV assay validated 0.5-80 ug/mL. Phoenix NLME 8.2; bootstrap (10,000 replicates) for confidence intervals. External validation in 130 retrospective patients (200 sparse levels). Baseline characteristics in Table 1."
  )

  ini({
    # Table 2 (final model). Vc = theta1, or theta1 x exp(theta2) if CRRT.
    lvc <- log(13.49); label("Central volume off CRRT (L)") # Table 2: theta1 = 13.49 L (RSE 6.70%)
    e_crrt_vc <- 0.54; label("Log-scale CRRT effect on Vc (unitless)") # Table 2: theta2 = 0.54 (RSE 15.12%)

    # CL = theta3 x Cr^theta5 x (BW/70)^theta6, or theta4 x (BW/70)^theta6 if CRRT.
    lcl <- log(2.85); label("Clearance off CRRT at CREAT = 1 mg/dL and 70 kg (L/h)") # Table 2: theta3 = 2.85 L/h (RSE 7.24%)
    lcl_crrt <- log(2.30); label("Clearance on CRRT at 70 kg (L/h)") # Table 2: theta4 = 2.30 L/h (RSE 15.43%)
    e_creat_cl <- -0.69; label("Power exponent of serum creatinine on off-CRRT CL (unitless)") # Table 2: theta5 = -0.69 (RSE 13.48%)
    e_wt_cl <- 0.73; label("Power exponent of body weight on CL (unitless)") # Table 2: theta6 = 0.73 (RSE 32.76%)

    lvp <- log(32.68); label("Peripheral volume (L)") # Table 2: theta7 = 32.68 L (RSE 9.42%)

    # Q = theta8, or theta8 x exp(theta9) if tMCS.
    lq <- log(18.81); label("Intercompartmental clearance without tMCS (L/h)") # Table 2: theta8 = 18.81 L/h (RSE 11.58%)
    e_mcs_temp_q <- -0.66; label("Log-scale tMCS effect on Q (unitless)") # Table 2: theta9 = -0.66 (RSE 22.29%)

    # IIV: exponential; Table 2 prints omega as the VARIANCE (footnote IIV% =
    # sqrt(exp(omega) - 1), reproduced to within the rounding of the printed
    # omegas). No IIV on Vc (removed for high shrinkage). One eta
    # is shared by both clearance arms (Supplementary Table S1, Model 3).
    etalcl ~ 0.100 # Table 2: omega CL = 0.100 (32.4%)
    etalvp ~ 0.095 # Table 2: omega Vp = 0.095 (31.7%)
    etalq ~ 0.081 # Table 2: omega Q = 0.081 (29.1%)

    addSd <- 1.39; label("Additive residual error (mg/L)") # Table 2: residual variability (additive) = 1.39
  })
  model({
    cl_noncrrt <- exp(lcl + etalcl) * CREAT^e_creat_cl * (WT / 70)^e_wt_cl
    cl_crrt <- exp(lcl_crrt + etalcl) * (WT / 70)^e_wt_cl
    cl <- (1 - RRT_CRRT_ACTIVE) * cl_noncrrt + RRT_CRRT_ACTIVE * cl_crrt
    vc <- exp(lvc + e_crrt_vc * RRT_CRRT_ACTIVE)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + e_mcs_temp_q * MCS_TEMP_STATUS + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
