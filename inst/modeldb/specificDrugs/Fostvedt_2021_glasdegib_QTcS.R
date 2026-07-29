Fostvedt_2021_glasdegib_QTcS <- function() {
  description <- paste(
    "Population PD model for glasdegib concentration-driven prolongation",
    "of the QT interval corrected for heart rate using a population-",
    "specific Fridericia-style correction factor (QTcS; beta estimated at",
    "0.312 in this cohort versus the fixed beta = 1/3 of QTcF) in 70 adult",
    "patients with advanced cancer pooled from two phase 1 dose-",
    "escalation trials (B1371001 in hematologic malignancies; B1371002",
    "in solid tumors). The exposure-response form is the same linear",
    "mixed-effects model as the companion QTcF extraction:",
    "QTcS = theta1 + theta2 * (CP_GLASDEGIB_NGML / 1000) + eta1 + W * eps,",
    "with additive random intercept (eta1 ~ N(0, omega2_1)) and a",
    "'thetarized' additive residual error (W is a fitted scalar;",
    "eps ~ N(0, 1)). Covariate analysis (age, sex, study) retained no",
    "covariates; a random effect on the slope was considered and dropped",
    "because shrinkage exceeded 20%. PD-only model: plasma glasdegib",
    "concentration is supplied as a time-varying covariate",
    "CP_GLASDEGIB_NGML (ng/mL). The slope is reported in the source",
    "publication on the microgram-per-mL scale (4.31 msec per microgram",
    "per mL); this file keeps that scaling and applies the unit",
    "conversion `CP_GLASDEGIB_NGML / 1000` inside model(). The source",
    "publication does not fit a population PK model; users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "concentration trajectory (no glasdegib popPK model exists in the",
    "nlmixr2lib registry).",
    sep = " "
  )

  reference <- paste(
    "Fostvedt LK, Shaik N, Martinelli G, Wagner AJ, Ruiz-Garcia A.",
    "Exposure-response modeling of the effect of glasdegib on cardiac",
    "repolarization in patients with cancer.",
    "Expert Review of Clinical Pharmacology 2021;14(7):927-935.",
    "doi:10.1080/17512433.2021.1925538.",
    sep = " "
  )

  vignette <- "Fostvedt_2021_glasdegib"

  # Paper-specific etas: as in the companion QTcF file, the additive
  # linear-scale random intercept (omega^2_1 = 291.6 msec^2) is translated
  # to a log-normal-equivalent omega^2_log = log(1 + 291.6 / 422.0^2) =
  # 0.0016379 (CV ~4.05% on the baseline). Approximation accuracy at
  # this CV is within ~0.8%. The eta name `etale0` pairs canonically
  # with the `le0` parameter (no paper_specific_etas declaration needed).

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external glasdegib plasma-concentration covariate)",
    concentration = "(observation QTcS is the QT interval corrected for heart rate using a population-specific correction factor beta estimated at 0.312 in this cohort, msec; driving covariate CP_GLASDEGIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_GLASDEGIB_NGML = list(
      description        = "Instantaneous glasdegib plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-QTcS expression QTcS = e0 + slope * (CP_GLASDEGIB_NGML / 1000).",
        "In Fostvedt 2021 this was the observed glasdegib plasma concentration around the time of the ECG collection (paired within +/-15 minutes of the PK sample; Methods 2.2 and 2.3). The source NM-TRAN column name in the model is CONC.",
        "The slope is reported on the microgram-per-mL scale (Table 3 'theta_2; msec/1000 ng/mL'); the canonical covariate is in ng/mL so the `/ 1000` rescaling lives inside model().",
        "Reference values observed: geometric-mean steady-state Cmax was 1137 ng/mL at 100 mg QD (therapeutic dose) and 2445 ng/mL at 200 mg QD (Fostvedt 2021 Methods 2.3 and Table 4).",
        "Set to 0 outside the drug-exposure window (the concentration-slope term then collapses to 0)."
      ),
      source_name        = "CONC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by forward-selection (alpha = 0.05) / backward-elimination (alpha = 0.001) stepwise covariate procedure on the QTcF and QTcS model parameters and not retained (Fostvedt 2021 Methods 2.3 and Results 3.3)."
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'sex') and not retained (Fostvedt 2021 Results 3.3). Of 70 patients in the pooled analysis, 28 (40%) were female (Fostvedt 2021 Table 1)."
    ),
    STUDY = list(
      description = "Source phase 1 study (hematologic malignancies B1371001 vs solid tumors B1371002)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'study (hematologic versus solid tumors)') and not retained (Fostvedt 2021 Results 3.3 and Discussion: 'no statistically significant differences were found based on the pre-specified backwards elimination significant criteria of 0.001'). Electrolyte imbalances and CYP3A4 inhibitor / inducer comedications were prohibited by protocol and therefore could not be tested as covariates (Fostvedt 2021 Methods 2.3)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 70L,
    n_studies        = 2L,
    n_observations   = 747L,
    age_range        = "25-89 years (Study B1371001 25-89, median 69; Study B1371002 27-76, median 61; pooled median 67)",
    weight_range     = NA_character_,
    sex_female_pct   = 40,
    race_ethnicity   = c(White = 80, Asian = 6, Black = 6, Other = 9),
    disease_state    = "Patients with advanced cancer pooled from two phase 1 dose-escalation trials: hematologic malignancies (Study B1371001, NCT00953758, n = 47) and solid tumor malignancies (Study B1371002, NCT01286467, n = 23). Eligibility excluded patients with screening QTcF > 470 msec.",
    dose_range       = paste(
      "Glasdegib QD oral monotherapy across a wide dose range (Fostvedt 2021 Section 2.1).",
      "Study B1371001: 5, 10, 20, 40, 80, 120, 180, 270, 400, and 600 mg QD; single lead-in dose for single-dose PK characterisation, continuous once-daily dosing in 28-day cycles thereafter.",
      "Study B1371002: 80, 160, 320, and 640 mg QD. Maximum dose tested (640 mg QD) was more than 6-fold above the selected clinical dose of 100 mg QD."
    ),
    regions          = NA_character_,
    notes            = paste(
      "Concentration-QTc analysis data set: 747 PK-ECG pairs collected within +/-15 minutes of paired PK draws, from 70 patients (Fostvedt 2021 Methods 2.2 and Table 1).",
      "Both studies were open-label, multicentre, phase 1 first-in-patient (B1371001) or solid-tumor phase 1 (B1371002) dose-escalation designs.",
      "QTcS is QT corrected for heart rate via the standard QTc = QT / (RR / 1000)^beta formula with beta estimated from the cohort by linear-mixed-effects regression of ln(QT) on ln(RR). Fostvedt 2021 Results 3.2 reports the QTcS correction factor beta = 0.312, which is very close to the Fridericia fixed beta = 1/3 = 0.333 (Fostvedt 2021 Results 3.2). Baseline QTcS mean (SD) = 419.1 (18.3) msec; on-treatment QTcS mean (SD) = 428.9 (24.4) msec (Fostvedt 2021 Table 2).",
      "Both QTcF (the companion model file) and QTcS were judged to adequately remove the correlation between heart rate and QTc (Fostvedt 2021 Results 3.2), and the authors developed parallel final E-R models for both endpoints (Fostvedt 2021 Results 3.3).",
      "A subsequent ICH E14-compliant thorough QT (TQT) study confirmed the conclusions: estimated slope 0.004 msec/(ng/mL) (TQT) vs. 0.00431 msec/(ng/mL) (this E-R analysis) for QTcS; estimated mean prolongation at 2445 ng/mL was 12.23 msec (95% CI 9.85-14.26) in this E-R analysis vs. the TQT 90% CI 10.14-13.64 msec (Fostvedt 2021 Discussion)."
    )
  )

  ini({
    # ==================================================================
    # Linear concentration-QTcS model (Fostvedt 2021 Section 2.3 and
    # Table 3 QTcS row). Same equation form as the companion QTcF model:
    #   QTcS_ij = theta1 + theta2 * CONC_ij + eta1_i + W * eps_ij
    # ==================================================================

    le0 <- log(422.0)
    label("Baseline QTcS intercept theta1 (msec)")
    # Fostvedt 2021 Table 3 QTcS 'theta_1, msec' = 422.0 (SE 2.2; RSE 0.5%).

    lslope <- log(4.31)
    label("Concentration-QTcS slope theta2 (msec per microgram per mL)")
    # Fostvedt 2021 Table 3 QTcS 'theta_2; msec/1000 ng/mL' = 4.31
    # (SE 0.37; RSE 8.6%). Equivalent to 0.00431 msec per ng/mL
    # (Fostvedt 2021 Discussion). Paper-specific log-transformed linear
    # concentration-effect slope; precedent in the registry includes
    # Liesenfeld_2006_dabigatran_aPTT.R and
    # deVriesSchultink_2018_anthracycline_troponinT.R.

    # ==================================================================
    # Random effect on the intercept (Fostvedt 2021 Table 3 QTcS
    # 'Variance of the random effect on the intercept, omega^2_1' =
    # 291.6 msec^2). No random effect on the slope (>20% shrinkage,
    # Section 2.3).
    #
    # Translation to log-normal-equivalent omega^2_log:
    #   omega^2_log = log(1 + 291.6 / 422.0^2)
    #               = log(1.001638) ~= 0.0016379
    # CV at the linear-scale typical value: sqrt(291.6) / 422.0 ~= 4.05%.
    # The vignette Errata documents the conversion.
    # ==================================================================
    etale0 ~ 0.0016379
    # Fostvedt 2021 Table 3 QTcS omega^2_1 = 291.6 (SE 51.1; RSE 17.5%;
    # shrinkage 4.4%); converted to log-normal equivalent for canonical
    # log-transformed intercept encoding.

    # ==================================================================
    # Residual error (Fostvedt 2021 Section 2.3 / Table 3 QTcS 'Residual
    # standard deviation, W' = 14.4 msec; additive on the QTcS scale).
    # ==================================================================
    addSd <- 14.4
    label("Additive residual error standard deviation on QTcS (msec)")
    # Fostvedt 2021 Table 3 QTcS W = 14.4 (SE 1.1; RSE 7.8%; shrinkage 3.9%).
  })

  model({
    # ==================================================================
    # Linear concentration-QTcS model. CP_GLASDEGIB_NGML is supplied as
    # a time-varying covariate (ng/mL); the slope is on the microgram-
    # per-mL scale (per Fostvedt 2021 Table 3 'theta_2; msec/1000 ng/mL'),
    # so the in-model `/ 1000` performs the ng/mL -> microgram/mL
    # conversion.
    # ==================================================================
    e0    <- exp(le0 + etale0)
    slope <- exp(lslope)

    QTcS <- e0 + slope * (CP_GLASDEGIB_NGML / 1000)
    QTcS ~ add(addSd)
  })
}
