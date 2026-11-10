Hui_2024_linezolid_continuous <- function() {
  description <- paste(
    "Two-compartment pharmacokinetic model with first-order elimination for",
    "intravenous linezolid given as a continuous infusion (600 mg loading",
    "infusion over 1 h, then 1200 mg/day at 50 mg/h) to critically ill adults",
    "in an intensive care unit (Hui 2024, continuous-infusion arm, n = 14).",
    "Hui 2024 fitted each patient individually in Phoenix WinNonlin, selected",
    "the two-compartment model over the one-compartment model by AIC, and",
    "reported only summary statistics (mean, SD, CV% and median) of the",
    "individual estimates per infusion arm (Table 2). This is therefore a",
    "typical-value model, not a population model: every parameter is fixed at",
    "the arm's median individual estimate of CL, V1, V2 and the distribution",
    "clearance CLD2, and there is no between-subject variability or",
    "residual-error model (residual SD fixed at zero). Covariates were only",
    "screened by ANOVA on the individual estimates (Table 4) and none enters",
    "the model. The intermittent-infusion arm of the same study is",
    "modellib('Hui_2024_linezolid_intermittent')."
  )
  reference <- paste(
    "Hui LA, Bodolea C, Popa A, Vlase AM, Hiriscau EI, Vlase L (2024).",
    "Linezolid Pharmacokinetics in Critically Ill Patients: Continuous Versus",
    "Intermittent Infusion. Antibiotics (Basel) 13(10):961.",
    "doi:10.3390/antibiotics13100961"
  )
  vignette <- "Hui_2024_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "linezolid", units = "mg", specimen = "tissue", verified = TRUE)
  )

  # Hui 2024 tested each of these against every individual PK parameter by
  # ANOVA in the Phoenix linear mixed-effects module (Table 4), pooling both
  # infusion arms, and reports p-values only -- no effect size, slope or
  # ratio -- so nothing is recoverable as a covariate effect. Documented here
  # to preserve the screen.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not modelled. Table 4 p-values were significant for CLD2 (0.018), K12 (0.025), V2 (0.023) and Cmax (0.024) but not for CL (0.419), V1 (0.604) or K10 (0.618); no effect size is reported. 2 female / 12 male in this arm (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not modelled. Table 4 p = 0.084 for CLD2 and 0.097 for K12, not significant elsewhere. Median 72.0 years (IQR 61.0-80.0) in this arm, younger than the intermittent-infusion arm (p = 0.043; Table 1)."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh class B)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not modelled. Hui 2024 grouped patients as Child-Pugh A (none or mild impairment) versus Child-Pugh B; Child-Pugh C was an exclusion criterion. Table 4 p = 0.054 for K10, not significant elsewhere. 8 of the 14 patients in this arm were Child-Pugh B (Results 2.1)."
    ),
    CRCL = list(
      description = "Creatinine clearance, raw mL/min (estimating equation not stated)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened, not modelled, as a five-band category (> 130, 61-130, 31-60, 15-30 and 0-14 mL/min; Methods 4.2). Table 4 p = 0.095 for K10, not significant elsewhere. This arm had 3 / 2 / 4 / 3 / 2 patients in those bands (Results 2.1). The paper does not say whether CrCl was measured or estimated, nor whether it was BSA-normalised; it is recorded here as raw mL/min, the unit printed with the bands."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not modelled. Table 4 p >= 0.241 for every parameter. Albumin values are not tabulated."
    ),
    RRT_CRRT_STATUS = list(
      description = "Extracorporeal renal replacement therapy (continuous venovenous haemodiafiltration or intermittent haemodialysis 3 times per week)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not modelled, as a single 'HD' column that pools CVVHDF and intermittent haemodialysis (Results 2.2). Table 4 p >= 0.120 for every parameter. The number of patients on each modality is not reported."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not modelled. Table 4 p >= 0.098 for every parameter. Mean 25.6 (SD 6.4) kg/m^2 in this arm (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14,
    n_studies = 1,
    age_median = "72.0 years (IQR 61.0-80.0)",
    weight_mean = "77.7 kg (SD 22.7)",
    sex_female_pct = 14.3,
    race_ethnicity = "Caucasian (100%)",
    disease_state = paste(
      "Critically ill adults in an intensive care unit (SOFA >= 2 with multiple",
      "comorbidities; mean SOFA 10.7, mean SAPS II 50.6) with a suspected or",
      "documented Gram-positive infection (urinary 8/14, septic shock 7/14,",
      "cutaneous 6/14, pulmonary 3/14, abdominal 1/14, bacteraemia 1/14).",
      "Child-Pugh C excluded."
    ),
    dose_range = "Linezolid 600 mg IV loading infusion over 1 h (600 mg/h), then 1200 mg/day by continuous IV infusion (50 mg/h) for 2-13 days (mean 7.2 days).",
    regions = "Romania (single centre, Cluj-Napoca)",
    notes = paste(
      "Hui 2024 Table 1, Results 2.1 and Methods 4.2-4.5. Prospective,",
      "open-label, two-arm study (NCT05801484), July 2022 - March 2023, with",
      "convenience allocation: the first 14 eligible patients received",
      "continuous infusion (the ward's standard practice) and the next 15",
      "intermittent infusion. Serum was sampled before the loading infusion and",
      "at 1, 2, 4, 8, 12, 16, 24, 36 and 48 h from its start, then once daily",
      "to the end of treatment (396 samples across both arms); linezolid was",
      "assayed by LC-MS (LOQ 0.256 ug/mL). The companion intermittent-infusion",
      "arm is Hui_2024_linezolid_intermittent."
    )
  )

  ini({
    # Hui 2024 Table 2, continuous-infusion (CI) group, MEDIAN of the 14
    # individual two-compartment WinNonlin estimates. The median is used rather
    # than the arithmetic mean because the individual estimates are strongly
    # right-skewed (CV 46-114% for these four parameters; the mean CLD2 of
    # 15.43 L/h is 1.5-fold the median) and a median is invariant under the
    # monotone transforms (log, 1/x) that relate a parameter to its half-life
    # or rate constant. The mean set is CL 6.95 L/h, V1 36.99 L, V2 26.17 L,
    # CLD2 15.43 L/h. The four medians give alpha = 0.77 /h (alpha half-life
    # 0.90 h) and beta = 0.086 /h (beta half-life 8.1 h) against the tabulated
    # medians 0.56 /h, 1.47 h, 0.11 /h and 6.20 h; the continuous profile
    # carries little information on the distribution phase, which is where
    # the per-parameter medians disagree most.
    lcl <- fixed(log(5.66)); label("Clearance CL (L/h)") # Table 2, CI group, CL median
    lvc <- fixed(log(35.91)); label("Central volume of distribution V1 (L)") # Table 2, CI group, V1 median
    lvp <- fixed(log(23.75)); label("Peripheral volume of distribution V2 (L)") # Table 2, CI group, V2 median
    lq <- fixed(log(10.01)); label("Intercompartmental distribution clearance CLD2 (L/h)") # Table 2, CI group, CLD2 median

    addSd <- fixed(0); label("Additive residual SD on linezolid serum concentration, not reported (mg/L)") # individual WinNonlin fits; no residual-error model reported
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
