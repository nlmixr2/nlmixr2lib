Cojutti_2020_darunavir_covid19 <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for boosted oral darunavir",
    "(darunavir/cobicistat 800/150 mg or darunavir/ritonavir 800/100 mg once daily) in",
    "hospitalised adults with SARS-CoV-2 disease (COVID-19). Apparent clearance CL/F falls with",
    "serum interleukin-6 by a power function (exponent -0.23) and the apparent volume rises with",
    "body surface area by a power function (exponent 1.44). Both covariates are centred on the",
    "cohort medians (IL-6 31 pg/mL, BSA 1.86 m^2); the centring values are not printed in the",
    "paper and were adopted from Table 1 after the paper's own Figure 4 simulations confirmed the",
    "IL-6 centring. The companion HIV-patient model is Cojutti_2020_darunavir_hiv.",
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

  covariateData <- list(
    IL6 = list(
      description = "Serum interleukin-6 concentration",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured on the day of the therapeutic-drug-monitoring assessment (Methods 2.1), i.e.",
        "one value per patient that the source treated as time-fixed over the sampled interval.",
        "Power effect on CL/F, cl = exp(lcl) * (IL6 / 31)^e_il6_cl. The centring value is NOT",
        "printed in the paper; 31 pg/mL is the SARS-CoV-2 cohort median (Table 1, IQR",
        "10-114.75). The paper's Figure 4 steady-state simulations (Results 3.2: median Cmin",
        "1460 / 7020 / 14140 ng/mL at IL-6 1 / 100 / 1000 pg/mL) are reproduced within a few",
        "percent at the upper two levels only when IL-6 is centred near 31 pg/mL; an uncentred",
        "form (reference 1 pg/mL) over-predicts every trough by 2.5-fold or more. See the",
        "vignette. Cohort range not reported (IQR 10-114.75 pg/mL; Figure 4 simulates up to",
        "1000 pg/mL)."
      ),
      source_name = "IL-6"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the apparent volume, vc = exp(lvc) * (BSA / 1.86)^e_bsa_vc. The",
        "centring value is NOT printed in the paper; 1.86 m^2 is the SARS-CoV-2 cohort median",
        "(Table 1, IQR 1.77-1.96), which is also the HIV cohort median. The BSA formula",
        "(DuBois / Mosteller / other) is unspecified in the source."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    n_observations = 60,
    age_median = "63 years (IQR 55-70.5)",
    weight_median = "75.0 kg (IQR 69.25-81.50)",
    bsa_median = "1.86 m^2 (IQR 1.77-1.96)",
    sex_female_pct = 40,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Hospitalised SARS-CoV-2 disease (COVID-19); Siddiqi-Mehra stage I n = 3, IIa n = 8,",
      "IIb n = 10, III n = 9. Serum IL-6 median 31.0 pg/mL (IQR 10-114.75)."
    ),
    dose_range = paste(
      "Darunavir 800 mg once daily boosted with cobicistat 150 mg (22 of 30) or ritonavir",
      "100 mg (8 of 30); TDM at least 48 h after treatment start (median 3 days)."
    ),
    regions = "Italy (Santa Maria della Misericordia University Hospital, Udine)",
    notes = paste(
      "Retrospective single-centre TDM study, 15 March - 15 May 2020. Each patient",
      "contributed two plasma samples, a trough just before a daily dose and a sample 2 h",
      "post-dose (Methods 2.1); LC-MS/MS assay with LLOQ 6 ng/mL. Demographics from Table 1.",
      "Fitted in Monolix 2019R1 (SAEM) separately from the HIV comparator cohort."
    )
  )

  ini({
    # Table 2, 'SARS-CoV-2 patients' 'Value (RSE%)' column (final model).
    # Monolix: all individual parameters log-normally distributed (Methods 2.2).
    lka <- log(0.74)
    label("Absorption rate constant ka (1/h)") # Table 2, ka = 0.74 1/h (RSE 30.6%)
    lcl <- log(4.10)
    label("Apparent clearance CL/F at IL-6 = 31 pg/mL (L/h)") # Table 2, CL/F = 4.10 L/h (RSE 10.1%)
    lvc <- log(88.41)
    label("Apparent volume of distribution Vd/F at BSA = 1.86 m^2 (L)") # Table 2, Vd = 88.41 L (RSE 7.8%)

    e_il6_cl <- -0.23
    label("Power exponent of IL-6 on CL/F (unitless)") # Table 2, beta IL6-CL/F = -0.23 (RSE 24.9%)
    e_bsa_vc <- 1.44
    label("Power exponent of BSA on Vd (unitless)") # Table 2, beta BSA-Vd = 1.44 (RSE 35.9%)

    # Table 2 'Between-subject variability' rows are Monolix omegas, i.e. SDs
    # of the normal random effect on the log scale (the '(%)' in the row label
    # cannot be a percent: 0.53% CV would be negligible). Variance = omega^2.
    etalka ~ 0.6724 # 0.82^2; Table 2, omega ka = 0.82 (RSE 28.3%)
    etalcl ~ 0.2809 # 0.53^2; Table 2, omega CL = 0.53 (RSE 15.2%)
    etalvc ~ 0.0225 # 0.15^2; Table 2, omega Vd = 0.15 (RSE 41.0%)

    propSd <- 0.09
    label("Proportional residual error (fraction)") # Table 2, b (proportional) = 0.09 (RSE 22.3%)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (IL6 / 31)^e_il6_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.86)^e_bsa_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
