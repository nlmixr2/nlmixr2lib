Amajor_2026_cefepime_plasma <- function() {
  description <- "Two-compartment IV population PK model for cefepime in PLASMA-EQUIVALENT concentrations in critically ill children with multiple organ dysfunction syndrome (MODS), 1 month to 17 years (Amajor 2026, plasma re-fit). The authors re-ran their final whole-blood volumetric absorptive microsampling (VAMS) model after dividing every VAMS concentration by the in vitro whole-blood:plasma ratio of 0.58, giving estimated plasma PK. Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 36 kg) and volumes linearly; clearance scales as a power of CKiD Under-25 eGFR (reference 100 mL/min/1.73 m^2) and central volume as a power of age (reference 10 years). Exponential between-subject variability on CL and V1 only, proportional residual error (magnitude not reported for this fit; carried from the whole-blood model). Companion model Amajor_2026_cefepime_wholeblood is the primary VAMS fit."
  reference <- paste(
    "Amajor V, Bwint A, Shein SL, Hume J, Rowan CM, Flygt LP, Malone MP,",
    "Stenson EK, Weiss SL, Zinter MS, Kannan S, Dahmer MK, Nofziger RA, Pai MP,",
    "Zuppa AF, Rhodes NJ, Hall MW, Scheetz MH, Downes KJ; Pediatric Acute Lung",
    "Injury and Sepsis Investigators (PALISI) Network. Cefepime pharmacokinetics",
    "in critically ill children with multiple organ dysfunction syndrome using",
    "volumetric absorptive microsampling. Antimicrob Agents Chemother.",
    "2026;70(7):e01736-25. doi:10.1128/aac.01736-25. PMCID PMC13321836.",
    sep = " "
  )
  vignette <- "Amajor_2026_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefepime is given intravenously, so the dose enters
  # `central` directly. The fitted concentrations are PLASMA-EQUIVALENT:
  # measured whole-blood VAMS concentrations divided by the in vitro
  # whole-blood:plasma ratio 0.58 (Methods, 'Population PK model
  # development'; Table S2). No plasma sample was assayed in the patients.
  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 36 kg is the cohort median (Amajor 2026 Table 1, 'Weight (kg)' 36 (4.3-213.6)) and the normalizing constant printed in all four Table 2 individual-parameter equations ('allometric scaling (centered on a weight of 36 kg)', Results 'PK model').",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate from the creatinine-only CKiD Under-25 (U25) equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Amajor 2026 Methods ('Data collection'): eGFR was calculated for each study day with the Chronic Kidney Disease in Children Under 25 (U25) formula from the HIGHEST daily serum creatinine, using serum creatinine only (cystatin C was not consistently available). Time-varying: handled as a Monolix regressor with interpolation between measured values (Methods, 'Population PK model development'). Reference 100 mL/min/1.73 m^2 is the normalizing constant of the Table 2 clearance equation, (GFR/100)^theta_GFR. Cohort median 87.37 (range 19.35-278.29) at the start of sampling (Table 1). The U25 equation is applied when preparing this column, not inside the model; an eGFR from a different estimator (bedside Schwartz was also screened, Table S4) rescales the renal term. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates.",
      source_name = "GFR (U25)"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 10 years is the cohort median (Table 1, 'Age' 10y (1m-17y)) and the normalizing constant of the Table 2 central-volume equation, (Age/10)^theta_Age. Enter infant ages as fractional years (1 month = 1/12); the power form requires AGE > 0. The three infants in the cohort were 1, 5 and 5 months old, all born at 37+ weeks gestation and > 40 weeks postmenstrual age (Table 1 footnote).",
      source_name = "Age"
    )
  )

  covariatesDataExcluded <- list(
    PELOD2 = list(
      description = "Daily Pediatric Logistic Organ Dysfunction-2 (PELOD-2) score",
      units = "(score)",
      type = "continuous",
      notes = "Screened on CL and V1 (Table S4). PELOD-2 on CL and on V1 met the OFV/BICc criteria in the third forward step but were not retained: PELOD-2 on CL did not reduce BSV, and PELOD-2 on V1 had a small effect and worsened goodness-of-fit plots (Results, 'PK model')."
    ),
    PROULX = list(
      description = "Daily Proulx MODS score",
      units = "(score)",
      type = "continuous",
      notes = "Screened on CL and V1 (Table S4). Proulx on V1 met the statistical criteria but did not reduce BSV and was not retained (Results, 'PK model')."
    ),
    FLUID_STATUS = list(
      description = "Daily fluid status (neutral, positive or negative volume)",
      units = "(category)",
      type = "categorical",
      notes = "Screened on CL and V1 (Table S4). Reduced BSV on both but the model became unstable (Monolix returned NaN eigenvalues) and other precisions worsened, so not retained (Results, 'PK model'; Table S4 phase 3)."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation in use",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V1 (Table S4), not significant. The empirical Bayes CL and V1 of the 4 ECMO patients (5.42 L/h and 21.2 L) were similar to the population (Results, 'PK model')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    n_observations = 377L,
    age_range = "1 month to 17 years",
    age_median = "10 years",
    weight_range = "4.3 - 213.6 kg",
    weight_median = "36 kg",
    height_median = "135 cm (range 56 - 198.1)",
    sex_female_pct = 54,
    disease_state = "Critically ill children with multiple organ dysfunction syndrome (MODS; 35 of 37 had MODS when sampled, 2 had single-organ dysfunction) receiving intravenous cefepime as standard of care. Median PELOD-2 6 (0-21), PRISM III 11 (0-30), Proulx 2 (0-5); 4 (11%) on ECMO; 22 (59%) immunoparalysed (TNF-alpha response < 200 pg/mL) at the start of dosing. Patients on continuous renal replacement therapy or plasmapheresis were excluded.",
    renal_function = "CKiD U25 creatinine-only eGFR median 87.37 mL/min/1.73 m^2 (range 19.35 - 278.29) at the start of sampling; serum creatinine median 0.48 mg/dL (0.15 - 3.08). Albumin median 3 g/dL (1.6 - 3.9). Daily values are in Table S3.",
    dose_range = "Median 47.6 mg/kg per dose (range 9.4 - 54.3) intravenously; dosing interval and infusion duration were clinician-chosen and are not reported.",
    regions = "Multicentre: 14 PARADIGM study sites of the Pediatric Acute Lung Injury and Sepsis Investigators (PALISI) network; countries not stated in the article.",
    sampling = "Up to 15 whole-blood VAMS samples (2 x 20 uL) per child over 3 consecutive calendar days, median 11 per child (range 2 - 15); 426 collected, 49 excluded (32 drawn during infusion, 13 collection errors, 4 implausible); LC-MS/MS, LLOQ 0.1 mg/L, none of the 377 analysed samples below LLOQ. Median time on drug before sampling 24 h (1 - 199).",
    notes = "AMPLE ancillary study of the prospective observational PARADIGM study. Fitted in Monolix 2024R1. One- and two-compartment models tested; the two-compartment model with allometric scaling and proportional error was retained. BSV on Q and V2 was poorly estimated and fixed to zero. Covariates: stepwise forward addition (dOFV >= 3.84) and backward elimination (dOFV >= 6.635), plus BICc, precision and GOF. Also screened and not retained (Table S4): bedside-Schwartz eGFR on CL (U25 was preferred), days since MODS onset, PRISM III, sex, mechanical ventilation, vasopressor use, immunoparalysis, TNF-alpha response and age on CL; PELOD-2, Proulx, fluid status and ECMO are listed in covariatesDataExcluded. This plasma model is a re-fit of the final VAMS model to VAMS concentrations divided by 0.58, the mean whole-blood:plasma ratio of 30 paired in vitro spiked samples (5 QC levels 0.2-150 mg/L, run-level ratios 0.52-0.62; Table S2); the paper calls it 'estimated plasma PK', and the ratio has not been validated clinically (Discussion). See Amajor_2026_cefepime_wholeblood for the primary whole-blood fit. The paper assumed 20% protein binding to compute free concentrations for fT>MIC."
  )

  ini({
    # Structural parameters, Amajor 2026 Table 3 'Estimated plasma PK
    # parameters using adjusted concentrations based on whole blood
    # VAMS:plasma ratio'. Reference subject: WT 36 kg, U25 eGFR
    # 100 mL/min/1.73 m^2, age 10 years.
    lcl <- log(2.97); label("Clearance at WT=36 kg, CRCL=100 mL/min/1.73 m^2 (L/h)") # Table 3 'CL (L/h)': 2.97 (SE 0.19, RSE 6.48%, 95% CI 2.62-3.37)
    lvc <- log(11.51); label("Central volume at WT=36 kg, AGE=10 years (L)") # Table 3 'V1 (L)': 11.51 (SE 0.93, RSE 8.12%, 95% CI 9.83-13.49)
    lq <- log(0.46); label("Intercompartmental clearance at WT=36 kg (L/h)") # Table 3 'Q (L/h)': 0.46 (SE 0.27, RSE 57.8%, 95% CI 0.18-1.20)
    lvp <- log(3.21); label("Peripheral volume at WT=36 kg (L)") # Table 3 'V2 (L)': 3.21 (SE 0.89, RSE 27.9%, 95% CI 1.91-5.38)

    # Allometric exponents are printed inside the Table 3 equations with no SE
    # row, as in the whole-blood model -- fixed, not estimated.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/36) for CL and Q (unitless)") # Table 3 equations: (WT/36)^0.75 on CL and Q
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/36) for V1 and V2 (unitless)") # Table 3 equations: (WT/36) on V1 and V2

    # Covariate exponents (estimated; each has its own SE row).
    e_crcl_cl <- 0.49; label("Power exponent on (CRCL/100) for CL (unitless)") # Table 3 'theta_GFR': 0.49 (SE 0.13, RSE 26.0%, 95% CI 0.24-0.74)
    e_age_vc <- -0.20; label("Power exponent on (AGE/10) for V1 (unitless)") # Table 3 'theta_Age': -0.20 (SE 0.06, RSE 27.7%, 95% CI -0.32 to -0.09)

    # Between-subject variability. Table 3 footnote: 'BSV reported as SD and
    # CV %' -- the first number is the log-scale SD omega, squared here to the
    # variance nlmixr2 expects. No BSV on Q or V2, as in the whole-blood model.
    etalcl ~ 0.1369 # 0.37^2; Table 3 'BSV CL (CV%)' = 0.37 (38.7%)
    etalvc ~ 0.0841 # 0.29^2; Table 3 'BSV V1 (CV%)' = 0.29 (29.5%)

    # Proportional residual error. Table 3 does not print an RV row (its
    # abbreviation footnote defines 'RV' but the table has no such row). The
    # value is carried from the whole-blood fit, Table 2 'RV' = 0.25: dividing
    # every observation by the same constant 0.58 leaves a proportional
    # residual unchanged, so the re-fit's residual magnitude is the same
    # quantity. See the vignette Assumptions section.
    propSd <- 0.25; label("Proportional residual error (fraction)") # carried from Table 2 'RV' 0.25 (whole-blood fit); not printed in Table 3
  })
  model({
    # Individual parameters, Amajor 2026 Table 3 equations (same form as
    # Table 2):
    #   CLi = CL * (WT/36)^0.75 * (GFR/100)^theta_GFR * exp(eta_CL)
    #   V1i = V1 * (WT/36)      * (Age/10)^theta_Age  * exp(eta_V1)
    #   Qi  = Q  * (WT/36)^0.75
    #   V2i = V2 * (WT/36)
    cl <- exp(lcl + etalcl) * (WT / 36)^e_wt_cl_q * (CRCL / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 36)^e_wt_vc_vp * (AGE / 10)^e_age_vc
    q <- exp(lq) * (WT / 36)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 36)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L: central / vc is plasma-equivalent cefepime in
    # mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
