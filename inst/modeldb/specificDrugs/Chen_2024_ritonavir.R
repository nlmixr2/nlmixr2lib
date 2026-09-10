Chen_2024_ritonavir <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral ritonavir 100 mg twice daily, given as the pharmacokinetic booster of the nirmatrelvir/ritonavir fixed combination, in critically ill adults with COVID-19; no covariate effects were retained, and the correlated interindividual variability on CL/F and V/F reproduces the ritonavir exposure distribution that drives the companion nirmatrelvir model (Chen 2024)."
  reference <- "Chen N, Yu X, Li L, Yang P, Dong R, Huang Y, Ling X, Shentu Q, Yu W, Jiang S. Target Attainment and Population Pharmacokinetics of Nirmatrelvir/Ritonavir in Critically Ill Adult Patients. Infect Drug Resist. 2024;17:4055-4065. doi:10.2147/IDR.S471918"
  vignette <- "Chen_2024_nirmatrelvir_ritonavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariate effects were retained in the final ritonavir model: "During the
  # covariate model building step, no significant covariate effects were
  # identified." (Chen 2024 Results, Population PK Model Development). The
  # demographics and laboratory values that were screened are documented in
  # covariatesDataExcluded below.
  covariateData <- list()

  # Screened during covariate model building but not retained in the final
  # ritonavir model; the paper reports no point estimate for any of them.
  # Source: Chen 2024 Methods (Data Collection) for the screened list and
  # Results (Population PK Model Development) for the negative outcome.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model (no significant covariate effects identified)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model. Cohort was 65% male (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the CKD-EPI equation",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model. Retained on CL/F in the companion nirmatrelvir model (see Chen_2024_nirmatrelvir)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    TPRO = list(
      description = "Total protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation (APACHE) II score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model. Cohort median 14 (IQR 11-21) per Table 1."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained in the final ritonavir model. Cohort median 6 (IQR 3-9) per Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    n_observations = 89L,
    age_median     = "69 years (IQR 63-77)",
    weight_median  = "60.6 kg (SD 12.8)",
    sex_female_pct = 35,
    disease_state  = "Critically ill adults (>=18 years) admitted to the intensive care unit with RT-PCR-confirmed SARS-CoV-2 infection and compatible pulmonary CT findings. Patients receiving continuous renal replacement therapy or extracorporeal membrane oxygenation were excluded. Median APACHE II score 14 (IQR 11-21); median SOFA score 6 (IQR 3-9).",
    renal_function = "Creatinine clearance (CKD-EPI) 78.5 +/- 34.4 mL/min; serum creatinine 100.2 +/- 79.9 umol/L (Table 1).",
    dose_range     = "Oral nirmatrelvir/ritonavir 300 mg/100 mg twice daily in 27 of 31 patients (87%) and 150 mg/100 mg twice daily in 4 of 31 (13%); the ritonavir dose was 100 mg twice daily in every patient. 29 of 31 (94%) received the tablets by nasal feeding and 2 (6%) orally, all for more than 5 days.",
    co_medication  = "Concomitant CYP3A4/5 weak inhibitor in 7 of 31 (22.6%) and weak inducer in 10 of 31 (32.3%) per Table 1; named co-medications in the Discussion include voriconazole, dexamethasone, methylprednisolone and omeprazole.",
    regions        = "China (The First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    notes          = "Prospective observational study, January-June 2023. Two to three serial plasma samples per patient collected after the second dose; 89 plasma samples from 31 patients. Ritonavir quantified by validated LC-MS/MS over 2.03-961.37 ng/mL. Measured ritonavir concentrations ranged 21.71-4966.72 ng/mL. Model fit in Phoenix NLME 8.1 with FOCE; final OFV 174.5. Baseline demographics in Table 1; final parameter estimates and 1000-sample bootstrap in Table 2."
  )

  ini({
    # Structural parameters, Chen 2024 Table 2 (ritonavir), "Final Model /
    # Estimate, Mean" column. The "CV, %" column of Table 2 is the relative
    # standard error of the estimate, not interindividual variability: the
    # Abstract reports the companion nirmatrelvir estimates as "Mean (SD)
    # 0.42 (0.10) h-1, 36.5 (8.5) L, 3.6 (0.26) L/h", which are exactly the
    # Table 3 estimates multiplied by the Table 3 "CV, %" values.
    lka <- log(0.64); label("Apparent first-order absorption rate constant (1/h)")  # Table 2: tvKa = 0.64 1/h (RSE 37.4%; bootstrap median 1.10, 95% CI 0.13-9.99)
    lvc <- log(84.9); label("Apparent central volume of distribution (L)")         # Table 2: tvV = 84.9 L (RSE 29.2%; bootstrap median 82.4, 95% CI 10.1-158.4)
    lcl <- log(10.3); label("Apparent oral clearance (L/h)")                      # Table 2: tvCL = 10.3 L/h (RSE 14.7%; bootstrap median 10.3, 95% CI 7.5-13.7)

    # Interindividual variability, exponential model on CL/F and V/F, with the
    # CL-V correlation carried in the off-diagonal element of the
    # variance-covariance matrix (Chen 2024 Results: "a correlation between CL
    # and V was observed and incorporated into the off-diagonal elements of the
    # variance-covariance matrix, which resulted in a significantly improved OFV
    # of 174.5"). Table 2 reports variances directly (omega^2), so no CV%
    # back-transformation is needed. The off-diagonal covariance is
    # CorrV-CL * sqrt(omega2CL * omega2V) = 0.73 * sqrt(0.57 * 1.27) = 0.621.
    # No interindividual variability was reported for ka.
    etalcl + etalvc ~ c(0.57,
                        0.621, 1.27)  # Table 2: omega^2 CL = 0.57, omega^2 V = 1.27, CorrV-CL = 0.73 (bootstrap 0.53, 1.64, 0.71)

    # Residual variability. The paper selected "a one-compartment model with a
    # log-additive error option" (Results, Population PK Model Development); in
    # Phoenix NLME the log-additive error model is C = Cpred * exp(eps), which
    # is nlmixr2's lnorm() residual, so the reported stdev is the additive SD on
    # the log-transformed concentration scale.
    expSd <- 0.34; label("Additive residual SD on the log-transformed concentration scale (log-normal)")  # Table 2: Stdev0 = 0.34 (RSE 12.7%; bootstrap median 0.34, 95% CI 0.24-0.44)
  })

  model({
    # Individual parameters (apparent, oral administration; F is not identifiable
    # separately and is subsumed into CL/F and V/F)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # One-compartment model with first-order absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    Cc ~ lnorm(expSd)
  })
}
