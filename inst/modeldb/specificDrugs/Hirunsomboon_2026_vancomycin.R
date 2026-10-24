Hirunsomboon_2026_vancomycin <- function() {
  description <- "One-compartment population PK model for intravenous vancomycin in critically ill adults with acute kidney injury on continuous renal replacement therapy (CVVHDF), with an uncentered log-linear age effect and an adsorptive-membrane (oXiris) effect on clearance (Hirunsomboon 2026)."
  reference <- "Hirunsomboon C, Pathumarak A, Nosoongnoen W, Tharanon V, Sukkha S. Population Pharmacokinetics of Vancomycin in Critically Ill Patients Undergoing Continuous Renal Replacement Therapy: Influence of Adsorptive Membranes. Clin Transl Sci. 2026;19:e70511. doi:10.1111/cts.70511"
  vignette <- "Hirunsomboon_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin is given by intravenous infusion straight into the sampled
  # serum compartment (Methods 2.2: concentrations from routine TDM by CMIA).
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Collected on the day of the vancomycin level measurement (Methods 2.2).",
        "Enters log-linearly and UNCENTERED on clearance, per the printed final CL equation",
        "'CL (L/h) = 2.82 x e^(-0.007 x Age (years) + 0.053 x MM)' (Results 3.2), so",
        "lcl = log(2.82) is the intercept at AGE = 0 and not the clearance of any study subject.",
        "At the cohort mean age of 61.5 years the typical CL is 2.82 * exp(-0.007 * 61.5) = 1.84 L/h",
        "without an adsorptive membrane. The uncentered reading is confirmed by the paper's own",
        "Monte Carlo PTA tables (Tables S3-S4); see the vignette. Cohort 61.51 +/- 16.01 years (Table 1).",
        sep = " "
      ),
      source_name = "Age"
    ),
    FILT_ADSORPTIVE = list(
      description = "Adsorptive CRRT hemofilter membrane (oXiris, AN69ST-derived) in use",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = oXiris adsorptive membrane, 0 = non-adsorptive membrane such as AN69 (Methods 2.2;",
        "Results 3.2 'use = 1, not use = 0'). Subject-level in the source (26 of 55 patients, 47.27%,",
        "Table 1). Exponential effect on total CL: exp(0.053) = 1.054, the 5.4% increase quoted in",
        "the Discussion. All patients were on CRRT (CVVHDF, Table S1), so the covariate is defined",
        "only while the circuit runs; the model has no separate off-CRRT clearance.",
        sep = " "
      ),
      source_name = "MM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and V (Table S2, CL models 2 and V model 4) but not retained. Cohort 60.49 +/- 17.33 kg (Table 1). Correlated with CRRT effluent flow (r = 0.5, Figure S3)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and V (Table S2) but not retained. Cohort 25.82 +/- 4.80 g/L (Table 1)."
    ),
    URINE_VOL_24H = list(
      description = "Residual urine volume",
      units = "mL/24h",
      type = "continuous",
      notes = "Screened on CL and V (Table S2) but not retained. Cohort 156.82 +/- 152.03 mL/day (Table 1)."
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "CRRT effluent flow rate",
      units = "mL/h",
      type = "continuous",
      notes = "Listed among the evaluated covariates (Methods 2.3) but absent from the Table S2 candidate models, having been collinear with body weight (r = 0.5, Figure S3). Cohort 1.99 +/- 0.59 L/h, i.e. 33.20 +/- 5.25 mL/kg/h (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    n_observations = 126L,
    age_range = "Adults >= 18 years; range not reported",
    age_median = "61.51 years (mean, SD 16.01)",
    weight_range = "Not reported",
    weight_median = "60.49 kg (mean, SD 17.33)",
    sex_female_pct = 45.5,
    race_ethnicity = "Not reported; single-centre Thai cohort",
    disease_state = "Critically ill ICU adults with acute kidney injury treated with CRRT in CVVHDF mode, receiving intravenous vancomycin. Mean SOFA 15.63 +/- 2.16; creatinine clearance 27.35 +/- 20.16 mL/min; serum albumin 25.82 +/- 4.80 g/L. End-stage renal disease (dialysis or kidney transplant history) excluded. 26 of 55 (47.27%) treated with the adsorptive oXiris membrane.",
    dose_range = "All patients received a loading dose (mean 1361.81 +/- 377.36 mg) followed by a maintenance dose of 1121.82 +/- 331.49 mg/day, adjusted by pharmacist-led TDM",
    regions = "Thailand (Ramathibodi Hospital, Bangkok)",
    renal_function = "AKI on CVVHDF; effluent flow 1.99 +/- 0.59 L/h (33.20 +/- 5.25 mL/kg/h); residual urine 156.82 +/- 152.03 mL/day",
    notes = "Retrospective single-centre study of ICU admissions 2014-2024 (Methods 2.1-2.2; Figure S1). Sparse routine-TDM sampling; vancomycin assayed by CMIA (Abbott ARCHITECT i), LLOQ 0.24 ug/mL. Fitted in MonolixSuite 2024R1 with SAEM. Baseline characteristics per Table 1."
  )

  ini({
    # Structural parameters -- Hirunsomboon 2026 Table 3 and the Results 3.2
    # final CL equation, CL = 2.82 * exp(-0.007 * Age + 0.053 * MM). Age is
    # uncentered, so 2.82 L/h is the intercept at AGE = 0.
    lcl <- log(2.82); label("Clearance intercept at AGE = 0 with a non-adsorptive membrane (L/h)") # Table 3: CL pop = 2.82 L/h (RSE 19.0%)
    lvc <- log(94.97); label("Volume of distribution (L)") # Table 3: V pop = 94.97 L (RSE 10.5%)

    e_age_cl <- -0.007; label("Log-linear effect of age on CL (1/year; CL ~ exp(e_age_cl * AGE))") # Table 3: theta CL,age = -0.007 (RSE 44.0%)
    e_filt_adsorptive_cl <- 0.053; label("Log-linear effect of the adsorptive oXiris membrane on CL (unitless)") # Table 3: theta CL,MM=1 = 0.053 (RSE 49.5%)

    # IIV. Table 3 prints 'IIV V (%) 21.63' and 'IIV CL (%) 20.01' beside
    # bootstrap medians of 0.19 and 0.17 on the raw scale, so the printed
    # values are the Monolix omega (SD of eta) x 100; variance = omega^2. No
    # CL-V correlation is reported.
    etalcl ~ 0.04004 # Table 3: IIV CL = 20.01 (omega 0.2001, RSE 22.4%); 0.2001^2
    etalvc ~ 0.04679 # Table 3: IIV V = 21.63 (omega 0.2163, RSE 46.6%); 0.2163^2

    # Residual error: proportional (Results 3.2), Monolix y = f * (1 + b * eps).
    propSd <- 0.26; label("Proportional residual error (fraction)") # Table 3: sigma = 0.26 (RSE 8.60%)
  })
  model({
    # Individual parameters -- Results 3.2 final CL equation; no covariates on V.
    cl <- exp(lcl + etalcl) * exp(e_age_cl * AGE + e_filt_adsorptive_cl * FILT_ADSORPTIVE)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One compartment, first-order elimination (Results 3.2), IV infusion.
    d/dt(central) <- -kel * central

    # Dose mg, V L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
