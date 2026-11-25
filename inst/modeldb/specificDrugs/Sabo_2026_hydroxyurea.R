Sabo_2026_hydroxyurea <- function() {
  description <- paste(
    "One-compartment population PK model for once-daily oral hydroxyurea",
    "in patients with sickle cell disease (Sabo 2026, OPTIMDREP trial",
    "NCT06464458; 22 patients, 20 children aged 2.5-16.4 years and 2 adult",
    "women). First-order absorption whose rate constant ka depends on the",
    "observed absorption phenotype (rapid absorber, NCA Tmax < 1 h: ka 9.93",
    "1/h; slow absorber: ka 1.36 1/h) and linear elimination. Apparent",
    "clearance and volume scale with body weight normalised to the cohort",
    "median of 34.9 kg (estimated exponents 0.64 and 0.80). Log-normal IIV",
    "on ka, CL/F and V/F; combined additive and proportional residual error."
  )
  reference <- paste(
    "Sabo AN, Nazon C, Paillard C, Kemmel V (2026). The Hydroxyurea",
    "Absorption Phenotype: A Key PK/PD Determinant in Sickle Cell Disease",
    "Treatment. Pharmaceutics 18(6):654.",
    "doi:10.3390/pharmaceutics18060654.",
    sep = " "
  )
  vignette <- "Sabo_2026_hydroxyurea"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "hydroxyurea", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "hydroxyurea", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model on CL/F (exponent 0.64) and V/F (exponent 0.80) with",
        "body weight normalised to 34.9 kg, the cohort median (Sabo 2026",
        "Eq. 3-4 and Table 4 'Present model' row: '34.9 kg (median BW,",
        "used for normalization)'). Monolix encodes it as a coefficient on",
        "log(BW/34.9) ('log-transformed median-normalized BW'), which is",
        "the same power model. Cohort mean 35.8 +/- 16.7 kg at inclusion",
        "(Table 1); paediatric range 11.4-60.0 kg (Results 3.1)."
      ),
      source_name = "BW"
    ),
    ABS_RAPID = list(
      description = "Observed absorption phenotype: 1 = rapid absorber, 0 = slow absorber",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (slow absorber, NCA Tmax >= 1 h)",
      notes = paste(
        "Sabo 2026 covariate KINETIC ('rapid = 1' or 'slow = 0', Results",
        "3.2.2). Assigned per pharmacokinetic profile from the observed",
        "non-compartmental Tmax: Tmax < 1 h (observed at 10 or 20 min) =",
        "rapid, Tmax >= 1 h = slow (Methods 2.4.1). 16 of 44 profiles",
        "(36.4%) were rapid (Results 3.2.1); two patients rapid at the",
        "first visit were slow at the fifth, so the phenotype can change",
        "between occasions. Enters ka as exp(1.99 * ABS_RAPID) (Eq. 2). For",
        "population simulation draw ABS_RAPID ~ Bernoulli(0.364)."
      ),
      source_name = "KINETIC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested as log(median-normalised age) on V/F (dOFV 0, Supplementary",
        "Table S3) and as a paediatric (0) / adult (1) category on ka, V/F",
        "and CL/F; not retained (Results 3.2.2, Supplementary Table S3)."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Tested on V/F (dOFV -1.34) and CL/F (dOFV -1.18); not retained",
        "(Supplementary Table S3). Cohort 14 male / 8 female (Table 1)."
      ),
      source_name = "sex"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on CL/F (dOFV -6.49) and V/F (dOFV -5.53); not retained",
        "against the backward-deletion threshold of 7.88 (Results 3.2.2,",
        "Supplementary Table S3). Cohort mean 34.2 +/- 7.7 umol/L (Table 1)."
      ),
      source_name = "creatinine"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    age_range = "2.5-35.7 years",
    age_median = "11.2 years (paediatric subgroup 9.7 years)",
    weight_range = "11.4-60.0 kg in the 20 children (adult weights not reported)",
    weight_median = "34.9 kg (31.0 kg in the 20 children)",
    sex_female_pct = 100 * 8 / 22,
    disease_state = paste(
      "Sickle cell disease requiring intensification of hydroxyurea",
      "therapy or with poorly controlled hydroxyurea treatment (failure to",
      "reach the maximum tolerated dose); normal renal function (cystatin",
      "C-based eGFR 126.7 +/- 14.2 mL/min/1.73 m2)."
    ),
    dose_range = paste(
      "Once-daily oral hydroxyurea, 20.3 +/- 2.4 mg/kg/day at inclusion;",
      "17 of 22 patients already at steady state, 5 initiating treatment."
    ),
    regions = "France (Strasbourg)",
    notes = paste(
      "OPTIMDREP randomised trial (NCT06464458). 391 hydroxyurea plasma",
      "concentrations: 44 full profiles (pre-dose, 10 and 20 min, 1, 2, 4",
      "and 6 h post-dose) at visits V0 and V5 plus 88 two-hour",
      "concentrations at V1-V4, assayed by GC-MS (LLOQ 0.79 mg/L).",
      "Demographics from Sabo 2026 Table 1 and Results 3.1."
    )
  )

  ini({
    # Structural parameters for a slow absorber (ABS_RAPID = 0) weighing
    # 34.9 kg (Sabo 2026 Table 2 'Model Estimation' column; Eq. 2-4).
    lka <- log(1.36); label("Absorption rate constant for a slow absorber (1/h)")            # Table 2 'ka (h-1)' = 1.36 (RSE 14.80%)
    lcl <- log(9.06); label("Apparent clearance CL/F at 34.9 kg (L/h)")                     # Table 2 'Cl/F (L/h)' = 9.06 (RSE 4.71%)
    lvc <- log(25.23); label("Apparent volume of distribution V/F at 34.9 kg (L)")          # Table 2 'Vd/F (L)' = 25.23 (RSE 7.69%)

    # Covariate effects (Table 2 'Covariates' rows; Eq. 2-4).
    e_abs_rapid_ka <- 1.99; label("Log-scale effect of the rapid absorption phenotype on ka (unitless)") # Table 2 'beta_ka_KINETIC' = 1.99 (RSE 21.20%); ka 9.93 vs 1.36 1/h (Results 3.2.2)
    e_wt_cl <- 0.64; label("Power exponent of body weight on CL/F (unitless)")              # Table 2 'beta_Cl_NBW' = 0.64 (RSE 13.90%)
    e_wt_vc <- 0.80; label("Power exponent of body weight on V/F (unitless)")               # Table 2 'beta_Vd_NBW' = 0.80 (RSE 14.80%)

    # IIV: Table 2 'IIV-... (%)' and Supplementary Table S3 'IIV final model
    # (CV%)'; omega^2 = log(1 + CV^2) for a log-normal parameter.
    etalka ~ log(1 + 0.6407^2) # Table 2 'IIV-ka (%)' = 64.07 CV%
    etalcl ~ log(1 + 0.3299^2) # Table 2 'IIV-Cl/F (%)' = 32.99 CV%
    etalvc ~ log(1 + 0.3473^2) # Table 2 'IIV-Vd/F (%)' = 34.73 CV%

    # Residual error: combined additive + proportional (Results 3.2.2).
    addSd <- 1.60; label("Additive residual error (mg/L)")                                  # Table 2 residual error 'a' = 1.60 mg/L (RSE 8.88%)
    propSd <- 0.35; label("Proportional residual error (fraction)")                         # Table 2 residual error 'b' = 0.35 (RSE 6.33%)
  })

  model({
    # Individual parameters (Sabo 2026 Eq. 2-4).
    ka <- exp(lka + e_abs_rapid_ka * ABS_RAPID + etalka)
    cl <- exp(lcl + etalcl) * (WT / 34.9)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 34.9)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give Cc in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
