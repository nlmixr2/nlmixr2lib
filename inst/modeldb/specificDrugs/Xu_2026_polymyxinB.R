Xu_2026_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B in critically ill Chinese adults in the intensive care unit, fitted to routine therapeutic drug monitoring samples (Xu 2026). No covariate was retained from the 22 screened (demographics, haematology, liver and renal markers, creatinine clearance, CRRT). Log-normal inter-individual variability on CL and V1 (none on V2 and Q); combined additive and proportional residual error."
  reference <- paste(
    "Xu Q, Li X, Huang S, Ding Q, Li Y, Yang N, Deng C, Tang B, Yang G, Pei Q.",
    "Exposure prediction and dose optimization of polymyxin B based on",
    "bayesian and machine learning.",
    "Front Pharmacol. 2026;17:1859709.",
    "doi:10.3389/fphar.2026.1859709. PMCID PMC13349818.",
    sep = " "
  )
  vignette <- "Xu_2026_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Plasma polymyxin B1 and B2 were quantified by UPLC-MS/MS
  # (Methods, 'Measurement of PMB concentrations'); the model describes a
  # single polymyxin B (PMB) concentration in plasma. The final model is
  # two-compartment with central volume V1 and peripheral volume V2
  # (Results, 'Development and evaluation of the initial PopPK model'; Table 3).
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SCREENED_NOT_RETAINED = list(
      description = "Candidate covariates screened by two-stage stepwise covariate modelling and not retained",
      units = "(various)",
      type = "continuous",
      notes = "Xu 2026 Methods 'PopPK modeling': gender, age, height, body weight, white blood cells, hemoglobin, platelets, neutrophils, C-reactive protein, procalcitonin, alanine-transaminase, aspartate-transaminase, total bilirubin, direct bilirubin, albumin, globulin, total bile acids, blood urea nitrogen, serum creatinine, urea, creatinine clearance and receipt of continuous renal replacement therapy were tested (forward p < 0.05, backward p < 0.01). Results: 'No covariates were included in the final model in the two-stage stepwise covariate modeling procedure.'"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = "92 plasma polymyxin B concentrations from routine therapeutic drug monitoring",
    age_range = "27-90 years",
    age_median = "53.5 years",
    weight_range = "40-120 kg",
    weight_median = "64 kg",
    height_range = "155-175 cm (median 169)",
    sex_female_pct = 15,
    race_ethnicity = "Chinese (single-centre cohort in Changsha, China)",
    disease_state = "Critically ill adults in the intensive care unit receiving intravenous polymyxin B with at least one therapeutic drug monitoring sample; 7 of 20 (35%) received continuous renal replacement therapy.",
    dose_range = "Intravenous polymyxin B 50-200 mg twice daily (median 100 mg), infusion rate 33.33-200 mg/h (median 100 mg/h); dose chosen by the physician.",
    regions = "China (Third Xiangya Hospital of Central South University), October 2019 to December 2020",
    renal_function = "Creatinine clearance median 47.38 (range 20.07-422.44) mL/min; serum creatinine 2.06 (0.17-4.65) mg/dL; blood urea nitrogen 17.6 (5.9-53.9) (Table 2).",
    notes = "Retrospective single-centre study. NONMEM 7.5, FOCE with interaction. Evaluation by NPDE, goodness-of-fit plots, a 1,000-sample VPC and a 1,000-sample nonparametric bootstrap (Table 3). Assay: UPLC-MS/MS, polymyxin B1 0.05-10 ug/mL and B2 0.011-1.097 ug/mL. Demographics from Table 2. The paper then used this model to simulate dense steady-state profiles (100 mg twice daily at 100 mg/h) for comparing Bayesian and XGBoost AUC0-12h predictions."
  )

  ini({
    # Structural parameters -- Xu 2026 Table 3 (final model). No covariate
    # was retained, so these are the population values for every patient.
    lcl <- log(3.05) ; label("Clearance CL (L/h)")                         # Table 3: CL = 3.05 L/h (RSE 11%; bootstrap median 3.10)
    lvc <- log(40.3) ; label("Central volume of distribution V1 (L)")      # Table 3: V1 = 40.3 L (RSE 11%; bootstrap median 41.93)
    lvp <- log(7.52) ; label("Peripheral volume of distribution V2 (L)")   # Table 3: V2 = 7.52 L (RSE 9%; bootstrap median 6.91)
    lq  <- log(29.7) ; label("Intercompartmental clearance Q (L/h)")       # Table 3: Q = 29.7 L/h (RSE 14%; bootstrap median 27.41)

    # Inter-individual variability -- Table 3 rows 'IIV on CL (%)' = 29.1 and
    # 'IIV on V1 (%)' = 33.3. The paper prints omega^2 x 100 in this column:
    # the same layout in Supplementary Table 2 ('IIV on CL (%)' 16.1, 'IIV on
    # V1 (%)' 31.9) is restated in the Methods (Scenario 1) as 'between-subject
    # variances (omega2) for CL and V1 ... 0.161 and 0.319'. The vignette
    # checks this reading against the simulated AUC0-12h range of Supplementary
    # Table 4. Results: IIV on CL and V1 only, no correlation reported.
    etalcl ~ 0.291 # Table 3: IIV on CL = 29.1, read as omega^2 = 0.291 (RSE 35%)
    etalvc ~ 0.333 # Table 3: IIV on V1 = 33.3, read as omega^2 = 0.333 (RSE 44%)

    # Residual variability -- Table 3, mixed additive-proportional error
    # (Results). The Scenario 1 residual model of the Methods,
    # sigma = sqrt(0.01^2 + (0.01 * Cpred)^2) for Supplementary Table 2
    # values 1% and 0.01 mg/L, shows both rows are standard deviations.
    propSd <- 0.192 ; label("Proportional residual error (fraction)")   # Table 3: proportional error = 19.2% (RSE 46%)
    addSd  <- 0.06  ; label("Additive residual error (mg/L)")           # Table 3: additive error = 0.06 mg/L (RSE 7%)
  })

  model({
    # Individual parameters -- Xu 2026 Table 3
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment. Polymyxin B is given as an intravenous infusion
    # (median 100 mg/h in the cohort; 100 mg/h in the paper's simulations),
    # so drug enters the central compartment directly.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
