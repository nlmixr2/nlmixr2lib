Cojutti_2020_darunavir_hiv <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for cobicistat-boosted",
    "oral darunavir (darunavir/cobicistat/emtricitabine/tenofovir alafenamide 800/150/200/10 mg",
    "once daily) in adults with chronic HIV infection at steady state. No covariate was retained.",
    "Fitted as the comparator for the SARS-CoV-2 model Cojutti_2020_darunavir_covid19.",
    sep = " "
  )
  reference <- paste(
    "Cojutti PG, Londero A, Della Siega P, Givone F, Fabris M, Biasizzo J, Tascini C, Pea F.",
    "Comparative Population Pharmacokinetics of Darunavir in SARS-CoV-2 Patients vs. HIV",
    "Patients: The Role of Interleukin-6. Clin Pharmacokinet. 2020;59(10):1251-1260.",
    "doi:10.1007/s40262-020-00933-8. PMCID: PMC7453069.",
    sep = " "
  )
  vignette <- "Cojutti_2020_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariatesDataExcluded <- list(
    IL6 = list(
      description = "Serum interleukin-6 concentration",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested in the stepwise covariate search but not retained in the HIV model (Results",
        "3.1: 'no covariate was found to be significantly associated with any of the darunavir",
        "pharmacokinetic parameters in the HIV final pharmacokinetic model'). HIV cohort median",
        "2.0 pg/mL (IQR 2.0-2.75; Table 1)."
      ),
      source_name = "IL-6"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested but not retained in the HIV model (Results 3.1). HIV cohort median 1.86 m^2",
        "(IQR 1.75-2.02; Table 1)."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25,
    n_studies = 1,
    n_observations = 50,
    age_median = "47 years (IQR 40-51)",
    weight_median = "75.0 kg (IQR 66.0-84.0)",
    bsa_median = "1.86 m^2 (IQR 1.75-2.02)",
    sex_female_pct = 28,
    race_ethnicity = "Not reported.",
    disease_state = "Chronic HIV infection on stable antiretroviral therapy; serum IL-6 median 2.0 pg/mL (IQR 2.0-2.75).",
    dose_range = "Darunavir/cobicistat/emtricitabine/tenofovir alafenamide 800/150/200/10 mg once daily (all 25).",
    regions = "Italy (Santa Maria della Misericordia University Hospital, Udine)",
    notes = paste(
      "Retrospective TDM comparator group sampled over the same period as the SARS-CoV-2",
      "cohort (15 March - 15 May 2020). Two plasma samples per patient (trough and 2 h",
      "post-dose). Because treatment was chronic, six once-daily doses were added virtually",
      "before the TDM samples to represent steady state (Methods 2.2). Fitted in Monolix",
      "2019R1 (SAEM)."
    )
  )

  ini({
    # Table 2, 'HIV patients' 'Value (RSE%)' column (final model = base model;
    # no covariate retained).
    lka <- log(0.58)
    label("Absorption rate constant ka (1/h)") # Table 2, ka = 0.58 1/h (RSE 36.6%)
    lcl <- log(10.3)
    label("Apparent clearance CL/F (L/h)") # Table 2, CL/F = 10.3 L/h (RSE 10.4%)
    lvc <- log(96.9)
    label("Apparent volume of distribution Vd/F (L)") # Table 2, Vd = 96.9 L (RSE 11.5%)

    # Monolix omegas (SD of the log-scale random effect); variance = omega^2.
    etalka ~ 1.2996 # 1.14^2; Table 2, omega ka = 1.14 (RSE 24.6%)
    etalcl ~ 0.1936 # 0.44^2; Table 2, omega CL = 0.44 (RSE 19.5%)
    etalvc ~ 0.0169 # 0.13^2; Table 2, omega Vd = 0.13 (RSE 61.5%)

    propSd <- 0.155
    label("Proportional residual error (fraction)") # Table 2, b (proportional) = 0.155 (RSE 64.2%)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
