Amajor_2026_cefepime_wholeblood <- function() {
  description <- "Two-compartment IV population PK model for cefepime in WHOLE BLOOD collected by volumetric absorptive microsampling (VAMS) in critically ill children with multiple organ dysfunction syndrome (MODS), 1 month to 17 years (Amajor 2026, final VAMS model). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 36 kg) and central and peripheral volumes linearly (exponent 1); clearance additionally scales as a power function of CKiD Under-25 (U25) creatinine-based eGFR (reference 100 mL/min/1.73 m^2) and central volume as a power function of age (reference 10 years). Exponential between-subject variability on CL and V1 only, proportional residual error. Companion model Amajor_2026_cefepime_plasma is the authors' re-fit to plasma-equivalent concentrations."
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
  # `central` directly. The concentrations this model was fitted to are
  # total WHOLE-BLOOD cefepime from VAMS tips (Methods, 'Cefepime sampling':
  # 'total, whole blood cefepime concentrations were measured'), not plasma.
  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "whole blood", verified = TRUE)
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
    notes = "AMPLE ancillary study of the prospective observational PARADIGM study. Fitted in Monolix 2024R1. One- and two-compartment models tested; the two-compartment model with allometric scaling and proportional error was retained. BSV on Q and V2 was poorly estimated and fixed to zero. Covariates: stepwise forward addition (dOFV >= 3.84) and backward elimination (dOFV >= 6.635), plus BICc, precision and GOF. Also screened and not retained (Table S4): bedside-Schwartz eGFR on CL (U25 was preferred), days since MODS onset, PRISM III, sex, mechanical ventilation, vasopressor use, immunoparalysis, TNF-alpha response and age on CL; PELOD-2, Proulx, fluid status and ECMO are listed in covariatesDataExcluded. Whole-blood VAMS concentrations are about 0.58 times plasma (in vitro, Table S2), so CL and V1 here are about 1.7-fold the plasma-referenced values; see Amajor_2026_cefepime_plasma. The paper assumed 20% protein binding to compute free concentrations for fT>MIC."
  )

  ini({
    # Structural parameters, Amajor 2026 Table 2 'Final VAMS population
    # pharmacokinetic parameter estimates'. Reference subject: WT 36 kg,
    # U25 eGFR 100 mL/min/1.73 m^2, age 10 years.
    lcl <- log(5.14); label("Clearance at WT=36 kg, CRCL=100 mL/min/1.73 m^2 (L/h)") # Table 2 'CL (L/h)': 5.14 (SE 0.33, RSE 6.55%, 95% CI 4.53-5.84)
    lvc <- log(20.53); label("Central volume at WT=36 kg, AGE=10 years (L)") # Table 2 'V1 (L)': 20.53 (SE 1.38, RSE 6.72%, 95% CI 18.01-23.41)
    lq <- log(0.64); label("Intercompartmental clearance at WT=36 kg (L/h)") # Table 2 'Q (L/h)': 0.64 (SE 0.31, RSE 48.4%, 95% CI 0.28-1.47)
    lvp <- log(4.76); label("Peripheral volume at WT=36 kg (L)") # Table 2 'V2 (L)': 4.76 (SE 1.5, RSE 31.5%, 95% CI 2.67-8.48)

    # Allometric exponents are printed inside the Table 2 equations with no SE
    # row, and Methods states weight was tested 'with powers of 0.75 and 1 for
    # clearance and volume terms, respectively' -- fixed, not estimated.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/36) for CL and Q (unitless)") # Table 2 equations: (WT/36)^0.75 on CL and Q
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/36) for V1 and V2 (unitless)") # Table 2 equations: (WT/36) on V1 and V2

    # Covariate exponents (estimated; each has its own SE row).
    e_crcl_cl <- 0.47; label("Power exponent on (CRCL/100) for CL (unitless)") # Table 2 'theta_GFR': 0.47 (SE 0.043, RSE 9.23%, 95% CI 0.38-0.55)
    e_age_vc <- -0.20; label("Power exponent on (AGE/10) for V1 (unitless)") # Table 2 'theta_Age': -0.20 (SE 0.058, RSE 28.5%, 95% CI -0.32 to -0.09)

    # Between-subject variability. Table 2 footnote: 'BSV reported as SD and
    # CV %' -- the first number is the log-scale SD omega, squared here to the
    # variance nlmixr2 expects. No BSV on Q or V2 (fixed to zero, Results).
    etalcl ~ 0.1444 # 0.38^2; Table 2 'BSV CL (CV%)' = 0.38 (39.0%)
    etalvc ~ 0.0784 # 0.28^2; Table 2 'BSV V1 (CV%)' = 0.28 (29.4%)

    # Proportional residual error (Results: 'a proportional error structure
    # best described the data').
    propSd <- 0.25; label("Proportional residual error (fraction)") # Table 2 'RV': 0.25 (SE 0.011, RSE 4.21%, 95% CI 0.23-0.27)
  })
  model({
    # Individual parameters, Amajor 2026 Table 2 equations:
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

    # Dose in mg, volumes in L: central / vc is whole-blood cefepime in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
