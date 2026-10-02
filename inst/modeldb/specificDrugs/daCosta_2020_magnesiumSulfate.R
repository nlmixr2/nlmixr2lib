daCosta_2020_magnesiumSulfate <- function() {
  description <- "One-compartment population PK model of intravenous magnesium sulfate (MgSO4-7H2O) in pregnant women with preeclampsia, with an exponential serum-creatinine effect on clearance and an exponential body-weight effect on volume; no endogenous magnesium baseline (da Costa 2020)."
  reference <- "da Costa TX, Azeredo FJ, Ururahy MAG, da Silva Filho MA, Martins RR, Oliveira AG. Population Pharmacokinetics of Magnesium Sulfate in Preeclampsia and Associated Factors. Drugs R D 2020;20:257-266. doi:10.1007/s40268-020-00315-2"
  vignette <- "daCosta_2020_magnesiumSulfate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/dL")

  compartmentData <- list(
    central = list(analyte = "magnesium", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Maternal body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on V: V = V_pop * exp(beta * (WT - 79.2)). da Costa 2020 prints the Results equation uncentered ('V = exp(2.5878 + 0.0752 x weight (kg) + 0.404)'), which gives V of about 5000 L at the cohort-mean weight; the 79.2 kg centring is the cohort mean from Table 1 and is supported by the Figure 3 visual predictive check (see the vignette's Assumptions and deviations).",
      source_name = "weight"
    ),
    CREAT = list(
      description = "Serum creatinine at baseline",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on CL: CL = CL_pop * exp(beta * (CREAT - 0.7)). Printed uncentered in the Results equation; centred here on the Table 1 cohort mean (0.7 mg/dL) consistently with the weight effect. Higher creatinine lowers CL.",
      source_name = "creatinine"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL (p = 0.33) and V (p = 0.35) in da Costa 2020 Table 4; not retained."
    ),
    TPROT = list(
      description = "Serum total protein",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL (p = 0.27) and V (p = 0.22) in da Costa 2020 Table 4; not retained."
    ),
    AGE = list(
      description = "Maternal age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL (p = 0.071) and V (p = 0.13) in da Costa 2020 Table 4; not retained."
    ),
    EGA = list(
      description = "Gestational age",
      units = "weeks",
      type = "continuous",
      notes = "Listed among the covariates investigated (Methods 2.5); no estimate reported; not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed among the covariates investigated (Methods 2.5); no estimate reported; not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 109L,
    n_studies = 1L,
    n_observations = 347L,
    age_range = "mean 25.8 +/- 7.4 years",
    weight_range = "mean 79.2 +/- 14.7 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported (Brazil)",
    disease_state = "Preeclampsia, third trimester (gestational age 35.2 +/- 4.2 weeks); maternal ICU. Baseline serum creatinine 0.7 +/- 0.3 mg/dL; baseline serum magnesium 1.9 +/- 0.6 mg/dL.",
    dose_range = "Zuspan regimen: 4 g MgSO4-7H2O IV over 30 min, then 1 g/h continuous IV infusion.",
    regions = "Brazil (Maternity School Januario Cicco, Natal)",
    notes = "Prospective observational cohort, June 2016 - February 2018. Serum magnesium (colorimetric, Mann-Yoe) sampled before the loading dose and 2, 6, 12 and 18 h after; 3-4 samples per woman. Estimation in Monolix 2018R1 (SAEM); 1000-replicate bootstrap with Rsmlx. UNITS: dose is mg of elemental Mg (multiply g MgSO4-7H2O by 24.305/246.47 = 0.0986; 1 g/h = 98.6 mg Mg/h) and Cc is mg/dL of elemental Mg, the axis unit of the paper's VPC (Figure 3). NO BASELINE: Methods state 'No endogenous magnesium baseline adjustment was made to the model', and the population predictions in Figure 1 are zero at the pre-dose samples, so Cc is the drug-derived magnesium only even though the fitted observations were total serum magnesium."
  )

  ini({
    # Structural parameters: da Costa 2020 Table 3 ('Estimated magnesium sulfate
    # population parameters'). The Results equation prints the same values on
    # the log scale: exp(0.3221) = 1.38 L/h and exp(2.5878) = 13.3 L.
    lcl <- log(1.38); label("Clearance for the reference subject (L/h)") # Table 3, CL_pop = 1.38
    lvc <- log(13.3); label("Volume of distribution for the reference subject (L)") # Table 3, V_pop = 13.3

    # Covariate coefficients (Monolix beta, log-linear): Table 3.
    e_creat_cl <- -0.0814; label("Exponential coefficient of serum creatinine on CL (per mg/dL)") # Table 3, Beta_CL_creatinine_mg_dL = -0.0814
    e_wt_vc <- 0.0752; label("Exponential coefficient of body weight on V (per kg)") # Table 3, Beta_V_Weith_kg = +0.0752

    # Random effects: Monolix reports omega as the SD of the log-normal eta, so
    # the variance is omega^2. Table 3 'Omega_V' = 0.404 and 'Omega_CL' = 0.015.
    etalcl ~ 0.000225 # Table 3, Omega_CL = 0.015 (SD); variance 0.015^2
    etalvc ~ 0.163216 # Table 3, Omega_V = 0.404 (SD); variance 0.404^2

    # Residual error: Monolix combined model with constant a and proportional b
    # (Table 3). Monolix's combined1 form, sd = a + b * f.
    addSd <- 1.76; label("Additive residual error (mg/dL)") # Table 3, a = 1.76
    propSd <- 0.000552; label("Proportional residual error (fraction)") # Table 3, b = 0.000552
  })
  model({
    # Individual parameters (Results equations, centred on the Table 1 cohort
    # means; see covariateData notes).
    cl <- exp(lcl + etalcl + e_creat_cl * (CREAT - 0.7))
    vc <- exp(lvc + etalvc + e_wt_vc * (WT - 79.2))
    kel <- cl / vc

    # One-compartment model with first-order elimination (Methods 2.4); IV
    # loading and maintenance infusions go to central.
    d/dt(central) <- -kel * central

    # mg/L -> mg/dL of elemental magnesium; drug-derived only (no baseline).
    Cc <- central / vc / 10

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
