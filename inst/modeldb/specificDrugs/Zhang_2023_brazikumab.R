Zhang_2023_brazikumab <- function() {
  description <- "Two-compartment population PK model with first-order lagged SC absorption and linear clearance for brazikumab (anti-IL-23 p19 mAb) in healthy adults and adults with Crohn's disease (Zhang 2023)"
  reference <- "Zhang N, Chan ML, Li J, Brohawn PZ, Sun B, Vainshtein I, Roskos LK, Faggioni R, Savic RM. Combining pharmacometric models with predictive and prognostic biomarkers for precision therapy in Crohn's disease: A case study of brazikumab. CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1945-1959. doi:10.1002/psp4.13044"
  vignette <- "Zhang_2023_brazikumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "brazikumab", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "brazikumab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "brazikumab", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL centred at 39 g/L (the pooled-cohort median):",
        "(ALB / 39)^-1.32 per Table 2 footnote a. Baseline value (time-fixed).",
        "UNITS CAVEAT: Zhang 2023 Table 1 labels this row 'BALB, mg/dL' with median 39.0 and",
        "range 32.0-46.0 (treatment arm) / 33.0-50.0 (placebo arm). Those magnitudes are",
        "unambiguously g/L (equivalently 3.2-5.0 g/dL) -- normal human serum albumin is",
        "approximately 35-50 g/L, and 39 mg/dL would be a lethal hypoalbuminaemia. The",
        "supplement analysis dataset column BALBU carries the same 32-50 range. The published",
        "'mg/dL' label is therefore a publication typo; the canonical g/L is used here and the",
        "reference value 39 is unchanged, so the published exponent applies as printed.",
        sep = " "
      ),
      source_name        = "BALBU"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with Crohn's disease)",
      notes              = paste(
        "Linear fractional effect on CL: (1 - 0.362 * DIS_HEALTHY), i.e. healthy participants",
        "have 36.2% lower CL than patients with CD at the same albumin (Table 2 footnote a:",
        "'for health status, patients with CD = 0, and healthy subjects = 1'). Source NONMEM",
        "coding in PSP-2023-0072-T-s02.mod is the GRP column with GRP = 1 (most common,",
        "reference) for CD patients and GRP = 2 for healthy participants; verified against the",
        "supplement analysis dataset, where GRP = 2 selects exactly the 30 phase Ib healthy",
        "participants. Zhang 2023 attributes the higher patient CL to faecal protein loss",
        "through the inflamed gut wall ('leaky gut').",
        sep = " "
      ),
      source_name        = "GRP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the reference category of the published Vc typical value)",
      notes              = paste(
        "Linear fractional effect on Vc applied to the MALE category:",
        "(1 + 0.214 * (1 - SEXF)), giving Vc = 3.27 L in females and 3.97 L in males",
        "(Zhang 2023 Results and Table 2 footnote b). Source column GNDR is coded male = 1,",
        "female = 0, so the canonical SEXF is the inversion SEXF = 1 - GNDR and the published",
        "coefficient 0.214 is carried verbatim on (1 - SEXF); this preserves the paper's",
        "female-reference typical value rather than re-deriving a male-reference one.",
        "Verified against the supplement analysis dataset: GNDR = 0 for 74 of the 119 phase IIa",
        "subjects, matching the paper's '74 women and 45 men'.",
        sep = " "
      ),
      source_name        = "GNDR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate model (SCM) on CL as linear and hockey-stick relationships; not retained in the final model. The s02.mod CLWGT1 / CLWGT2 thetas are commented out and fixed to 0."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened in the SCM on CL and on Vc (hockey-stick, centred at 23.44 kg/m^2); not retained in the final model. The s02.mod CLBMI1 / V2BMI1 / V2BMI2 thetas are commented out and fixed to 0."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the SCM (Methods: 'Covariates, including subject population, sex, race, age, body weight, baseline albumin, and body mass index, were screened'); not retained in the final model."
    ),
    RACE = list(
      description = "Race",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened in the SCM; not retained in the final model. The phase IIa cohort was 93% White (Table 1), so the non-White strata were too small to support an effect."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 153L,
    n_studies      = 2L,
    age_range      = "18-61 years (phase IIa; not reported for phase Ib)",
    age_median     = "35 years (phase IIa pooled arms)",
    weight_range   = "44.0-158.8 kg (phase IIa; not reported for phase Ib)",
    weight_median  = "66.9 kg (phase IIa pooled arms)",
    sex_female_pct = 51.0,
    race_ethnicity = c(White = 93.3, Black = 5.0, Other = 1.7),
    disease_state  = "Pooled: 30 healthy adults (phase Ib) plus 123 adults with mild-to-severe or moderate-to-severe active Crohn's disease who had failed or were intolerant to anti-TNF-alpha therapy (4 in phase Ib, 119 in phase IIa)",
    dose_range     = "Phase Ib: 70, 210, 420, or 700 mg IV, or 210 mg SC, Q4W (3 doses). Phase IIa double-blind induction: 700 mg IV over at least 60 min on day 1 and day 29; open-label period: 210 mg SC Q4W from week 12 through week 112.",
    regions        = "Multinational (NCT01258205 phase Ib; NCT01714726 phase IIa)",
    trials         = c("NCT01258205 (phase Ib)", "NCT01714726 (phase IIa)"),
    notes          = paste(
      "Baseline demographics from Zhang 2023 Table 1. The population PK analysis pooled all PK",
      "data from both studies including the phase IIa open-label period: 34 phase Ib subjects",
      "(30 healthy participants, 4 with CD; 4 women / 30 men; median 21 PK samples per subject,",
      "range 14-21) and 119 phase IIa patients with CD (74 women / 45 men; median 9 PK samples",
      "per subject, range 1-10). Race percentages above are the phase IIa cohort (Table 1;",
      "race was not reported for phase Ib) and sex_female_pct is the pooled 78/153.",
      "Reference covariate values for typical-subject predictions: ALB = 39 g/L, patient with CD",
      "(DIS_HEALTHY = 0), female (SEXF = 1). Serum brazikumab was measured by a validated ECL",
      "immunoassay with LLOQ 12.5 ng/mL and quantitative range 12.5-12,800 ng/mL",
      "(Appendix S1, 'PK Assay for Brazikumab').",
      sep = " "
    )
  )

  ini({
    # Structural parameters (Zhang 2023 Table 2, final population PK model).
    # Reference subject: patient with Crohn's disease (DIS_HEALTHY = 0), female
    # (SEXF = 1), baseline albumin 39 g/L.
    lcl     <- log(0.26);    label("Clearance (CL, L/day)")                                # Table 2: CL in female patients with CD = 0.26 L/d (RSE 5%)
    lvc     <- log(3.27);    label("Central volume of distribution (Vc, L)")               # Table 2: Vc in female subjects = 3.27 L (RSE 5%)
    lvp     <- log(2.64);    label("Peripheral volume of distribution (Vp, L)")            # Table 2: Vp = 2.64 L (RSE 8%)
    lq      <- log(0.412);   label("Intercompartmental clearance (Q, L/day)")              # Table 2: Q = 0.412 L/d (RSE 19%)
    lka     <- log(0.286);   label("First-order SC absorption rate constant (ka, 1/day)")  # Table 2: ka = 0.286 1/d (RSE 11%)
    ltlag   <- log(0.0296);  label("SC absorption lag time (Tlag, day)")                   # Table 2: Tlag = 0.0296 d (RSE 14%)
    lfdepot <- log(0.88);    label("SC bioavailability (F, fraction)")                     # Table 2: F = 0.88 (RSE 5%)

    # Covariate effects (Zhang 2023 Table 2 with the parameterisation given in
    # footnotes a and b).
    e_alb_cl         <- -1.32;  label("Power exponent of baseline albumin on CL (centred at 39 g/L)")  # Table 2: effect of baseline albumin on CL = -1.32 (RSE 40%)
    e_dis_healthy_cl <- -0.362; label("Healthy-participant fractional change in CL")                   # Table 2: effect of health status on CL = -0.362 (RSE 10%)
    e_male_vc        <-  0.214; label("Male fractional change in Vc vs. the female reference")         # Table 2: effect of male gender on Vc = 0.214 (RSE 37%)

    # Inter-individual variability. Zhang 2023 Methods state the PK IIVs are
    # log-normal (theta_ij = theta_TVj * exp(eta_ij)), so the Table 2 IIV column
    # is read as the log-scale variance omega^2 directly. Vp, Q, ka, and F carry
    # no IIV: the corresponding $OMEGA elements are 0 FIX in
    # PSP-2023-0072-T-s02.mod and Table 2 leaves their IIV cells blank.
    etalcl ~ 0.131   # Table 2: IIV on CL = 0.131 (RSE 33%)
    etalvc ~ 0.0502  # Table 2: IIV on Vc = 0.0502 (RSE 18%)

    # Residual error. Zhang 2023 Methods: "A proportional residual error model was
    # found to best fit the PK data ... Cobs(t)i = Cpred(t)i * (1 + eps_i)". The
    # additive term of the s02.mod combined-error expression
    # W = sqrt(PROP^2 * IPRED^2 + ADD^2) is 0 FIX, leaving pure proportional error.
    propSd <- 0.249; label("Proportional residual error (fraction)")  # Table 2: proportional error = 24.9% (RSE 9%)
  })

  model({
    # 1. Covariate-effect terms.
    # Power-form albumin effect on CL, centred at the 39 g/L cohort median.
    alb_cl <- (ALB / 39)^e_alb_cl
    # Linear fractional health-status effect on CL (reference: patient with CD).
    dis_cl <- 1 + e_dis_healthy_cl * DIS_HEALTHY
    # Linear fractional sex effect on Vc. The published coefficient is the MALE
    # effect and the canonical column is SEXF, so it is applied to (1 - SEXF).
    sex_vc <- 1 + e_male_vc * (1 - SEXF)

    # 2. Individual PK parameters.
    cl <- exp(lcl + etalcl) * alb_cl * dis_cl
    vc <- exp(lvc + etalvc) * sex_vc
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (Appendix S1 Equations 1-6). SC doses enter depot; IV doses
    # enter central directly.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability and lag time apply to the SC depot only.
    f(depot)    <- exp(lfdepot)
    alag(depot) <- exp(ltlag)

    # 6. Observation. Dose in mg / volume in L gives mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
