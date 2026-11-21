GandaraMireles_2026_sertraline <- function() {
  description <- paste(
    "Two-compartment population PK model for oral sertraline in Mexican pediatric",
    "patients (mean age 12.2 years) with anxiety disorders (Gandara Mireles 2026), fitted in",
    "Monolix (SAEM) to 218 sparse steady-state concentrations from 85 patients.",
    "First-order absorption with a lag time, first-order elimination, log-normal IIV on",
    "all six parameters and a combined additive + proportional residual error.",
    "Body weight and BMI act on CL/F and V1/F as power functions centred on the paper's",
    "typical patient (47.8 kg, 22.0 kg/m^2). The paper prints neither the covariate",
    "equation nor its coefficients; the four exponents were back-solved by the",
    "maintainers from the Supplementary Figure S1 forest plot (see the vignette), and",
    "the BMI exponent on V1/F is the least certain of the four.",
    sep = " "
  )
  reference <- paste(
    "Gandara Mireles JA, Cesar Grijalva Avila J, Villanueva Fierro I, Mora Munoz C,",
    "Sanchez Campos P, Salas Name S, Lazcano Franco MA, Arceo Duran A, Rosales Silva KP,",
    "Ontiveros Sanchez H, Nieto Pescador MG, Romero Soto O, Loera Castaneda V,",
    "Patron Romero L, Almanza Reyes H. Exposure-response relationship of sertraline in",
    "pediatric patients with anxiety disorders: a population pharmacokinetic analysis.",
    "Front Pharmacol. 2026;17:1862275. doi:10.3389/fphar.2026.1862275.",
    "PMID 42422083. PMCID PMC13342236.",
    sep = " "
  )
  vignette <- "GandaraMireles_2026_sertraline"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F and V1/F, centred on the 47.8 kg typical patient of",
        "Supplementary Figure S1. Exponents back-solved from that forest plot by the",
        "maintainers (no coefficient is printed in the paper). Cohort mean (SD) 47.76",
        "(15.02) kg without AEs and 51.87 (13.22) kg with AEs (Table 1).",
        sep = " "
      ),
      source_name = "Body weight (Table 1; Supplementary Figure S1)"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F and V1/F, centred on the 22.0 kg/m^2 typical patient of",
        "Supplementary Figure S1. Exponents back-solved from that forest plot by the",
        "maintainers. Cohort mean (SD) 21.2 (4.5) kg/m^2 without AEs and 23.9 (4.2) with",
        "AEs (Table 1).",
        sep = " "
      ),
      source_name = "BMI (Table 1; Supplementary Figure S1)"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods 2.7 states that 'PK parameters were scaled according to body surface",
        "area using linear normalization, with a reference value of 1.24 m2', but BSA is",
        "not mentioned again: the Results, Table 2 and Supplementary Figure S1 describe",
        "a weight + BMI covariate model whose typical patient (47.8 kg, BMI 22.0) has a",
        "Mosteller BSA near 1.41 m^2 at the 1.50 m mean height, not 1.24 m^2. Not",
        "encoded; see the vignette's Assumptions and deviations.",
        sep = " "
      ),
      source_name = "Body surface area (Methods 2.7)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened in the stepwise covariate search (Methods 2.8) and not retained.",
        "Cohort mean (SD) 12.22 (4.01) years without AEs and 12.08 (3.99) with AEs.",
        sep = " "
      ),
      source_name = "Age (Table 1)"
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened (Methods 2.8) and not retained. 48 male / 37 female (Table 1).",
        sep = " "
      ),
      source_name = "Sex (Table 1)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sertraline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sertraline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sertraline", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 85L,
    n_studies = 1L,
    n_concentrations = 218L,
    age_mean = "12.2 years (SD about 4.0; Table 1 by AE group: 12.22 +/- 4.01 and 12.08 +/- 3.99)",
    weight_mean = "about 49.0 kg (Table 1 by AE group: 47.76 +/- 15.02 and 51.87 +/- 13.22 kg)",
    height_mean = "1.50 m (Table 1: 1.50 +/- 0.10 and 1.49 +/- 0.10 m)",
    bmi_mean = "about 22.0 kg/m^2 (Table 1 by AE group: 21.2 +/- 4.5 and 23.9 +/- 4.2)",
    sex_female_pct = 43.5,
    race_ethnicity = "Mexican (two centres in Durango, Mexico); not further stratified",
    disease_state = "Pediatric anxiety disorders (DSM-5) on sertraline monotherapy",
    dose_range = "Started at 25 mg/day and titrated; mean 35 +/- 15 mg/day once daily",
    regions = "Mexico (Mental Health Hospital and AMCII Specialty Hospital, Durango)",
    sampling = paste(
      "Sparse steady-state sampling, 1-3 samples per patient at predose and 1-25 h",
      "post-dose after at least about 5 days of treatment; HPLC with fluorescence",
      "detection, LLOQ 2 ng/mL, values below LLOQ handled as censored.",
      sep = " "
    ),
    notes = paste(
      "Baseline demographics from Table 1 (59 patients without and 26 with adverse",
      "events). Model fitted in Monolix with SAEM; evaluated by goodness-of-fit, VPC",
      "and a 1000-replicate bootstrap. The paper's exposure-response analyses (t-tests",
      "and logistic regressions on post hoc Cavg, Cmax and CL/F) are descriptive and",
      "carry no PD model.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Table 2, 'Final model' column. They are the values for
    # the typical patient of Supplementary Figure S1 (47.8 kg, BMI 22.0 kg/m^2).
    lka <- log(0.7)
    label("Absorption rate constant ka (1/h)") # Table 2 'Ka (h-1)' final model 0.7 (RSE 18.0%)
    ltlag <- log(0.30)
    label("Absorption lag time (h)") # Table 2 'Tlag (h)' final model 0.30 (RSE 22.0%)
    lcl <- log(59.7)
    label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' final model 59.7 (RSE 9.0%)
    lvc <- log(1200)
    label("Apparent central volume V1/F (L)") # Table 2 'V1/F (L)' final model 1200 (RSE 12.0%)
    lq <- log(161)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h)' final model 161 (RSE 15.0%)
    lvp <- log(1350)
    label("Apparent peripheral volume V2/F (L)") # Table 2 'V2/F (L)' final model 1350 (RSE 14.0%)

    # Covariate exponents. The paper prints no equation and no coefficient for the
    # retained weight and BMI effects; these are back-solved from the Supplementary
    # Figure S1 forest plot. Each effect there is a pair of percentage changes
    # relative to the typical patient, at covariate values the figure does not
    # state. A power model (COV / ref)^e evaluated at ref * (1 - x) and ref * (1 + x)
    # solves both points of an effect exactly for one (x, e) pair: both weight
    # effects give the same x = 0.3079 (47.8 +/- 14.72 kg, about one pooled SD of
    # Table 1 weight), which supports both the power form and the symmetric-endpoint
    # reading. BMI on CL/F gives x = 0.1804 (22.0 +/- 3.97 kg/m^2). The vignette
    # reproduces the arithmetic.
    e_wt_cl <- 0.635
    label("Power exponent of (WT/47.8) on CL/F (unitless)") # Back-solved from Supplementary Figure S1 'Weight -> CL/F' -20.85% / +18.60% (digitised)
    e_bmi_cl <- 0.303
    label("Power exponent of (BMI/22) on CL/F (unitless)") # Back-solved from Supplementary Figure S1 'BMI -> CL/F' -5.85% / +5.15% (digitised)
    e_wt_vc <- 0.847
    label("Power exponent of (WT/47.8) on V1/F (unitless)") # Back-solved from Supplementary Figure S1 'Weight -> V1/F' -26.78% / +25.53% (digitised)
    # BMI on V1/F (-6.79% / +8.81%) has a larger upward than downward change, which
    # no power or exponential model reproduces at symmetric endpoints. The exponent
    # is the log-scale least-squares fit of both points at the BMI endpoints found
    # for CL/F (22.0 +/- 3.97); the single-point values are 0.353 (low) and 0.509
    # (high). The least certain of the four exponents.
    e_bmi_vc <- 0.417
    label("Power exponent of (BMI/22) on V1/F (unitless)") # Least-squares back-solve from Supplementary Figure S1 'BMI -> V1/F' -6.79% / +8.81% (digitised)

    # IIV: Table 2 'Interindividual variability (CV, %)' column, log-normal (Methods
    # 2.7), converted with omega^2 = log(1 + CV^2). Diagonal: the paper reports the
    # ETA_CL-ETA_V1 r = 0.73 only as a post hoc eta-correlation diagnostic.
    etalka ~ log(1 + 0.25^2) # Table 2 Ka IIV 25 CV%
    etaltlag ~ log(1 + 0.30^2) # Table 2 Tlag IIV 30 CV%
    etalcl ~ log(1 + 0.28^2) # Table 2 CL/F IIV 28 CV%
    etalvc ~ log(1 + 0.35^2) # Table 2 V1/F IIV 35 CV%
    etalq ~ log(1 + 0.40^2) # Table 2 Q/F IIV 40 CV%
    etalvp ~ log(1 + 0.38^2) # Table 2 V2/F IIV 38 CV%

    # Residual error: Table 2 'Model error parameters' a and B, the Monolix
    # combined error with additive a (ng/mL) and proportional b.
    addSd <- 0.882
    label("Additive residual error (ng/mL)") # Table 2 'Model error parameters' a = 0.882
    propSd <- 0.723
    label("Proportional residual error (fraction)") # Table 2 'Model error parameters' B = 0.723
  })

  model({
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl <- exp(lcl + etalcl) * (WT / 47.8)^e_wt_cl * (BMI / 22)^e_bmi_cl
    vc <- exp(lvc + etalvc) * (WT / 47.8)^e_wt_vc * (BMI / 22)^e_bmi_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
