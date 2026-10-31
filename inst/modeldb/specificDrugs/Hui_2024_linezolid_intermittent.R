Hui_2024_linezolid_intermittent <- function() {
  description <- paste(
    "Two-compartment pharmacokinetic model with first-order elimination for",
    "intravenous linezolid given as an intermittent infusion (600 mg over 1-2 h",
    "every 12 h) to critically ill adults in an intensive care unit (Hui 2024,",
    "intermittent-infusion arm, n = 15). Hui 2024 fitted each patient",
    "individually in Phoenix WinNonlin, selected the two-compartment model over",
    "the one-compartment model by AIC, and reported only summary statistics",
    "(mean, SD, CV% and median) of the individual estimates per infusion arm",
    "(Table 2). This is therefore a typical-value model, not a population",
    "model: every parameter is fixed at the arm's median individual estimate",
    "of CL, V1, V2 and the distribution clearance CLD2, and there is no",
    "between-subject variability or residual-error model (residual SD fixed",
    "at zero). Covariates were only screened by ANOVA on the individual",
    "estimates (Table 4) and none enters the model. The continuous-infusion",
    "arm of the same study is modellib('Hui_2024_linezolid_continuous')."
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
  # ANOVA in the Phoenix linear mixed-effects module (Table 4) and reports
  # p-values only -- no effect size, slope or ratio -- so nothing is
  # recoverable as a covariate effect. Documented here to preserve the screen.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not modelled. Table 4 p-values were significant for CLD2 (0.018), K12 (0.025), V2 (0.023) and Cmax (0.024) but not for CL (0.419), V1 (0.604) or K10 (0.618); no effect size is reported. 5 female / 10 male in this arm (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not modelled. Table 4 p = 0.084 for CLD2 and 0.097 for K12, not significant elsewhere. Median 75.0 years (IQR 75.0-81.0) in this arm, older than the continuous-infusion arm (p = 0.043; Table 1)."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh class B)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not modelled. Hui 2024 grouped patients as Child-Pugh A (none or mild impairment) versus Child-Pugh B; Child-Pugh C was an exclusion criterion. Table 4 p = 0.054 for K10, not significant elsewhere. 7 of the 15 patients in this arm were Child-Pugh B (Results 2.1)."
    ),
    CRCL = list(
      description = "Creatinine clearance, raw mL/min (estimating equation not stated)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened, not modelled, as a five-band category (> 130, 61-130, 31-60, 15-30 and 0-14 mL/min; Methods 4.2). Table 4 p = 0.095 for K10, not significant elsewhere. This arm had 0 / 6 / 6 / 2 / 1 patients in those bands (Results 2.1). The paper does not say whether CrCl was measured or estimated, nor whether it was BSA-normalised; it is recorded here as raw mL/min, the unit printed with the bands."
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
      notes = "Screened, not modelled. Table 4 p >= 0.098 for every parameter. Mean 29.5 (SD 5.2) kg/m^2 in this arm (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15,
    n_studies = 1,
    age_median = "75.0 years (IQR 75.0-81.0)",
    weight_mean = "85.8 kg (SD 14.5)",
    sex_female_pct = 33.3,
    race_ethnicity = "Caucasian (100%)",
    disease_state = paste(
      "Critically ill adults in an intensive care unit (SOFA >= 2 with multiple",
      "comorbidities; mean SOFA 8.7, mean SAPS II 46.1) with a suspected or",
      "documented Gram-positive infection (septic shock 11/15, pulmonary 6/15,",
      "urinary 4/15, cutaneous 4/15, abdominal 2/15). Child-Pugh C excluded."
    ),
    dose_range = "Linezolid 600 mg IV infused over 1 or 2 h (300-600 mg/h) every 12 h for 2-10 days (mean 5.8 days).",
    regions = "Romania (single centre, Cluj-Napoca)",
    notes = paste(
      "Hui 2024 Table 1, Results 2.1 and Methods 4.2-4.5. Prospective,",
      "open-label, two-arm study (NCT05801484), July 2022 - March 2023, with",
      "convenience allocation: the first 14 eligible patients received",
      "continuous infusion and the next 15 intermittent infusion. Serum was",
      "sampled before the first infusion and at 1, 2, 4, 8, 12, 16, 24, 36 and",
      "48 h from its start, then once daily to the end of treatment (396",
      "samples across both arms); linezolid was assayed by LC-MS (LOQ 0.256",
      "ug/mL). The companion continuous-infusion arm is",
      "Hui_2024_linezolid_continuous."
    )
  )

  ini({
    # Hui 2024 Table 2, intermittent-infusion (II) group, MEDIAN of the 15
    # individual two-compartment WinNonlin estimates. The median is used rather
    # than the arithmetic mean because the individual estimates are strongly
    # right-skewed (CV 53-118% for these four parameters; the mean CL of 7.88
    # L/h is 1.5-fold the median) and a median is invariant under the monotone
    # transforms (log, 1/x) that relate a parameter to its half-life or rate
    # constant. The mean set is CL 7.88 L/h, V1 47.46 L, V2 26.50 L, CLD2 29.47
    # L/h. The four medians are mutually consistent: they give alpha = 1.25 /h
    # (alpha half-life 0.55 h) and beta = 0.079 /h (beta half-life 8.8 h)
    # against the tabulated medians 1.21 /h, 0.57 h, 0.07 /h and 9.25 h.
    lcl <- fixed(log(5.36)); label("Clearance CL (L/h)") # Table 2, II group, CL median
    lvc <- fixed(log(40.99)); label("Central volume of distribution V1 (L)") # Table 2, II group, V1 median
    lvp <- fixed(log(24.01)); label("Peripheral volume of distribution V2 (L)") # Table 2, II group, V2 median
    lq <- fixed(log(18.16)); label("Intercompartmental distribution clearance CLD2 (L/h)") # Table 2, II group, CLD2 median

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
