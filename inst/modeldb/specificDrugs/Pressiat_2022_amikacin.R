Pressiat_2022_amikacin <- function() {
  description <- "Two-compartment population PK model for a single 30-minute intravenous infusion of amikacin in critically ill adults with nosocomial sepsis, with or without veno-arterial extracorporeal membrane oxygenation (V-A ECMO) (Pressiat 2022). KDIGO acute-kidney-injury stage (1, 2, 3 vs 0) lowers clearance, total body weight raises the central volume, and V-A ECMO support enlarges the peripheral volume; all three enter as Monolix exponential covariate effects. Proportional residual error."
  reference <- "Pressiat C, Kudela A, De Roux Q, Khoudour N, Alessandri C, Haouache H, Vodovar D, Woerther PL, Hutin A, Ghaleh B, Hulin A, Mongardon N. Population Pharmacokinetics of Amikacin in Patients on Veno-Arterial Extracorporeal Membrane Oxygenation. Pharmaceutics. 2022;14(2):289. doi:10.3390/pharmaceutics14020289. PMC8879580."
  vignette <- "Pressiat_2022_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Amikacin is infused into plasma and total plasma amikacin was assayed by
  # enzyme-linked immunoturbidimetry (Pressiat 2022 Section 2.3).
  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight (weight of the day)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Pressiat 2022 Section 2.2 defines TBW as 'total body weight: weight of the day'; the 30 mg/kg dose was computed on it. Enters the central volume as the UNCENTRED Monolix exponential form V1 = 8.77 * exp(0.015 * TBW), so 8.77 L is the intercept at TBW = 0 and V1 = 25.8 L at 72 kg. The printed centred power form V1 = 8.77 * (TBW/72)^0.014 is contradicted by the paper's own PC-VPC (Figure 2): see the lvc comment in ini() and the vignette. Cohort TBW median (IQR) 70 (65-84) kg control and 75 (60-87) kg V-A ECMO (Table 1).",
      source_name = "TBW"
    ),
    KDIGO_AKI_1 = list(
      description = "KDIGO acute-kidney-injury stage 1 indicator (1 = stage 1, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (KDIGO stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Level indicator for the paper's KDIGO covariate, stage 1 (Table 3 row 'KDIGO stage 1/CL'). KDIGO stages AKI from serum creatinine and urine output (Section 2.1). Mutually exclusive with KDIGO_AKI_2 and KDIGO_AKI_3. 5 of 39 patients (Table 1).",
      source_name = "KDIGO"
    ),
    KDIGO_AKI_2 = list(
      description = "KDIGO acute-kidney-injury stage 2 indicator (1 = stage 2, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (KDIGO stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Level indicator for the paper's KDIGO covariate, stage 2 (Table 3 row 'KDIGO stage 2/CL'). Mutually exclusive with KDIGO_AKI_1 and KDIGO_AKI_3. 8 of 39 patients (Table 1).",
      source_name = "KDIGO"
    ),
    KDIGO_AKI_3 = list(
      description = "KDIGO acute-kidney-injury stage 3 indicator (1 = stage 3, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (KDIGO stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Level indicator for the paper's KDIGO covariate, stage 3 (Table 3 row 'KDIGO stage 3/CL'). Mutually exclusive with KDIGO_AKI_1 and KDIGO_AKI_2. Patients on renal replacement therapy were excluded from the study (Section 2.1). 6 of 39 patients (Table 1).",
      source_name = "KDIGO"
    ),
    ECMO_STATUS = list(
      description = "Veno-arterial extracorporeal membrane oxygenation support indicator (1 = V-A ECMO, 0 = control)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control group, no ECMO)",
      notes = "The paper's 'ECMO group' covariate (Table 3 row 'ECMO group/V2'). Group-level, time-fixed: all 24 ECMO-group patients were on V-A ECMO during the single studied amikacin dose (flow rate 4.2 (4-4.8) L/min, Table 1). Enters the peripheral volume as exp(0.76 * ECMO_STATUS).",
      source_name = "ECMO"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 39L,
    n_studies = 1L,
    n_observations = 215L,
    age_median = "62 (52-75) years control; 60 (51-64) years V-A ECMO (median (IQR), Table 1)",
    weight_median = "70 (65-84) kg control; 75 (60-87) kg V-A ECMO (median (IQR), Table 1)",
    sex_female_pct = 33,
    race_ethnicity = "Not reported (single-centre French cohort)",
    disease_state = "Critically ill adults in a surgical ICU requiring empirical antimicrobial therapy including amikacin for nosocomial sepsis: 15 controls (87% post-cardiac-surgery) and 24 on V-A ECMO (cardiogenic shock 54%, right heart failure 25%, refractory cardiac arrest 13%, graft dysfunction 8%). SAPS II 46 (39-58) vs 60 (52-69); SOFA 11 (9-12) vs 13 (11-16).",
    renal_function = "Measured (UV/P) creatinine clearance 120 (86-191) mL/min control and 18 (7-54) mL/min V-A ECMO. KDIGO stage 0/1/2/3: 10/2/2/1 control and 11/3/6/5 V-A ECMO as printed in Table 1 (the ECMO counts sum to 25 for 24 patients). Chronic dialysis, renal replacement therapy during the study and cirrhosis were exclusion criteria.",
    dose_range = "Single dose, 30 mg/kg total body weight in 50 mL 5% glucose over 30 minutes, rounded up to a multiple of 125 mg; received 29 (24-33) mg/kg control and 32 (30-35) mg/kg V-A ECMO. Only the first administration was studied.",
    regions = "France (Henri Mondor Hospital, Creteil)",
    notes = "Prospective study July 2013 - September 2015 (Section 2.1). Samples at end of infusion, 3-6 h, 7-9 h, 9-12 h and 24 h, then every 12 h while the concentration exceeded 2.5 mg/L; 5.5 samples per patient on average. LLOQ 0.8 mg/L, ULOQ 40 mg/L (Section 2.3). Fitted with Monolix 2020R1 (SAEM)."
  )

  ini({
    lcl <- log(4.45)
    label("Clearance at KDIGO stage 0 (L/h)") # Pressiat 2022 Table 3: CL = 4.45 L/h (RSE 6.33%)

    # Pressiat 2022 V1 = 8.77 L is the Monolix intercept of an UNCENTRED
    # exponential weight effect, log(V1) = log(8.77) + 0.015 * TBW, so the
    # typical V1 is 25.8 L at 72 kg. The printed equation
    # 'V1 = 8.77 x (TBW/72)^0.014' would make 8.77 L the typical V1 at
    # 72 kg; that reading predicts a median end-of-infusion concentration
    # of about 190 mg/L after 30 mg/kg, against about 85 mg/L for the
    # predicted median of the paper's PC-VPC (Figure 2) and a highest
    # observed concentration of about 122 mg/L. The uncentred reading
    # predicts about 80 mg/L. See the vignette.
    lvc <- log(8.77)
    label("Central volume intercept at TBW = 0 of the exponential weight effect (L)") # Pressiat 2022 Table 3: V1 = 8.77 L (RSE 40.2%)
    lvp <- log(15.90)
    label("Peripheral volume in the control group (L)") # Pressiat 2022 Table 3: V2 = 15.90 L (RSE 20.0%)
    lq <- log(6.96)
    label("Inter-compartmental clearance (L/h)") # Pressiat 2022 Table 3: Q = 6.96 L/h (RSE 17.5%)

    # Monolix categorical-covariate coefficients: CL = 4.45 * exp(beta) for
    # each KDIGO stage versus stage 0. The printed equation writes them as
    # bases, 'Cl = 4.45 x (-0.41)^(KDIGO=1) x ...', which cannot hold for
    # negative numbers; exp(beta) gives 0.66, 0.55 and 0.39.
    e_kdigo_aki_1_cl <- -0.41
    label("Log-scale effect of KDIGO AKI stage 1 on CL (unitless)") # Pressiat 2022 Table 3: 'KDIGO stage 1/CL' = -0.41 (RSE 35.4%)
    e_kdigo_aki_2_cl <- -0.59
    label("Log-scale effect of KDIGO AKI stage 2 on CL (unitless)") # Pressiat 2022 Table 3: 'KDIGO stage 2/CL' = -0.59 (RSE 18.8%)
    e_kdigo_aki_3_cl <- -0.93
    label("Log-scale effect of KDIGO AKI stage 3 on CL (unitless)") # Pressiat 2022 Table 3: 'KDIGO stage 3/CL' = -0.93 (RSE 14.9%)

    e_wt_vc <- 0.015
    label("Exponential coefficient of TBW on V1 (1/kg)") # Pressiat 2022 Table 3: 'TBW/V1' = 0.015 (RSE 34.0%); the printed equation shows 0.014

    e_ecmo_status_vp <- 0.76
    label("Log-scale effect of V-A ECMO support on V2 (unitless)") # Pressiat 2022 Table 3: 'ECMO group/V2' = 0.76 (RSE 35.8%); the printed equation shows 0.79 as a base

    # Monolix omega = SD of the log-normal random effect; variance = omega^2.
    etalcl ~ 0.0625 # Pressiat 2022 Table 3: omega CL = 0.25 (RSE 13.0%); 0.25^2
    etalvc ~ 0.1296 # Pressiat 2022 Table 3: omega V1 = 0.36 (RSE 19.5%); 0.36^2
    etalvp ~ 0.3969 # Pressiat 2022 Table 3: omega V2 = 0.63 (RSE 19.1%); 0.63^2

    propSd <- 0.16
    label("Proportional residual error (fraction)") # Pressiat 2022 Table 3: Sigma prop = 0.16 (RSE 7.1%)
  })

  model({
    # Monolix categorical effects on CL, stage 0 as reference.
    cl <- exp(lcl + e_kdigo_aki_1_cl * KDIGO_AKI_1 + e_kdigo_aki_2_cl * KDIGO_AKI_2 +
      e_kdigo_aki_3_cl * KDIGO_AKI_3 + etalcl)
    # Uncentred exponential weight effect on V1 (see the ini() note).
    vc <- exp(lvc + e_wt_vc * WT + etalvc)
    # V-A ECMO enlarges V2 by exp(0.76) = 2.14-fold.
    vp <- exp(lvp + e_ecmo_status_vp * ECMO_STATUS + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
