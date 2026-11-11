Tachibana_2025_valemetostat_QTcP <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-QTc model for the",
    "EZH1/EZH2 dual inhibitor valemetostat in 100 adults with",
    "relapsed/refractory non-Hodgkin lymphoma (phase 1 study DS3201-A-J101)",
    "or acute myeloid / acute lymphoblastic leukemia (phase 1 study",
    "DS3201-A-U102) receiving 100-700 mg once daily. The endpoint is the",
    "change from baseline in the population-specific corrected QT interval",
    "(DeltaQTcP, ms; QTcP = QT / RR^0.443): DeltaQTcP = (e0 + etae0) +",
    "(slope + etaslope) * CP_VALEMETOSTAT_NGML + e_qtc_bl_e0 * (QTC_BL -",
    "429) + eps, with intercept -0.00395 ms, slope 0.000938 ms/(ng/mL)",
    "(95% CI excludes 0) and baseline-QTcP coefficient -0.287 ms/ms.",
    "Additive, correlated random effects on the intercept (SD 5.35 ms) and",
    "slope (SD 0.00192 ms/(ng/mL), correlation 0.451). Additive residual",
    "SD 9.94 ms in J101, multiplied by 1.20 in U102 (STUDY_U102). No",
    "placebo arm and no nominal-time term, so the typical prediction at",
    "the cohort-mean baseline is e0 + slope * C (1.71 ms at the 1830 ng/mL",
    "steady-state Cmax of 200 mg once daily). The 429 ms centering value is",
    "derived (the paper does not print the mean baseline QTcP). PD-only",
    "model: the valemetostat plasma concentration is supplied as the",
    "time-varying covariate CP_VALEMETOSTAT_NGML (ng/mL). The companion",
    "Fridericia model is Tachibana_2025_valemetostat_QTcF.R.",
    sep = " "
  )

  reference <- paste(
    "Tachibana M, Poland B, Kakurai Y, Chen Y, Li C, Lau Y.",
    "Concentration-QTc Analysis of Valemetostat in Patients With",
    "Hematologic Malignancies. Clin Transl Sci. 2025;18(11):e70391.",
    "doi:10.1111/cts.70391.",
    sep = " "
  )

  vignette <- "Tachibana_2025_valemetostat_QTc"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external valemetostat plasma-concentration covariate)",
    concentration = "(observation QTcP is the CHANGE FROM BASELINE in the population-specific corrected QT interval, DeltaQTcP, in ms; driving covariate CP_VALEMETOSTAT_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_VALEMETOSTAT_NGML = list(
      description = "Instantaneous total valemetostat (parent) plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the linear term (slope + etaslope) * CP_VALEMETOSTAT_NGML; Tachibana 2025 Table 2 reports the slope in ms per ng/mL, so no in-model unit rescaling is required.",
        "In the source analysis this was the observed concentration of a PK sample drawn 0-15 min after the last ECG of a triplicate; concentrations below the lower limit of quantification (1 ng/mL in J101, 20 ng/mL in U102) and missing pre-dose concentrations on C1D1 were set to 0 (Methods 'Concentration-QTc Model Data Collection and Assembly').",
        "Observed range 0-8310 ng/mL, overall mean 1021 ng/mL (Table S1). Steady-state geometric mean Cmax by dose: 1460-3690 ng/mL in J101 (150-300 mg) and 330-4850 ng/mL in U102 (100-700 mg); 1830 ng/mL at the approved 200 mg once-daily dose (Table 3)."
      ),
      source_name = "Conc_i,t (observed valemetostat concentration)"
    ),
    QTC_BL = list(
      description = "Patient's pre-dose baseline population-specific corrected QT interval (QTcP_i,0 = QT / RR^0.443), the triplicate-averaged ECG taken pre-dose on Cycle 1 Day 1. Time-fixed per patient. Enters the intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 429); set QTC_BL = 429 ms for the typical patient.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Population-specific correction: the exponent 0.443 minimised the absolute Pearson correlation between baseline QTc and RR in the pooled baseline ECGs (Methods 'Data Analysis'; Results 'QT Correction and Exploratory Analysis'). Supply QTcP computed with the same exponent.",
        "Centering reference: the model centers on the mean baseline QTcP over all patients, which the paper does not print. 429 ms is derived from the Table 1 means: QTcP = QTcF * RR^(1/3 - 0.443), so 417 ms * (60 / 77.2 bpm)^(-0.1097) = 428.7 ms. Using the mean heart rate in place of the per-patient RR leaves an error of the order of 1 ms in the centering constant, which shifts a prediction by about 0.3 ms (1 * 0.287) and does not affect the typical prediction at QTC_BL = 429.",
        "Eligibility excluded QTcF > 450 ms at screening in both studies."
      ),
      source_name = "QTcP_i,0 (baseline QTcP)"
    ),
    STUDY_U102 = list(
      description = "Study indicator: 1 = record from DS3201-A-U102 (relapsed/refractory AML or ALL), 0 = record from DS3201-A-J101 (relapsed/refractory B- or T-cell NHL).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (DS3201-A-J101)",
      notes = paste(
        "Selects the residual SD only: addSd (9.94 ms) in J101 and addSd * e_study_u102_addSd (9.94 * 1.20 = 11.93 ms) in U102 (Table 2 'SD multiplier for study U102').",
        "'Different residual standard deviations (SDs) for the two studies were found to significantly improve the fit'; study population (NHL versus leukemia) was not retained on any fixed effect (Results 'Modeling Results').",
        "Does not change the typical prediction."
      ),
      source_name = "study (J101 / U102)"
    )
  )

  covariatesDataExcluded <- list(
    TUMTP_LEUK = list(
      description = "Study population: AML/ALL (U102) versus NHL (J101).",
      units = "(binary)",
      type = "binary",
      notes = "Tested individually in the stepwise forward-selection / backward-elimination search (p < 0.01) and not retained on any fixed effect (Results 'Modeling Results'); no estimate is reported. The study does enter the residual SD, as STUDY_U102."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained; no estimate reported."
    ),
    SEXF = list(
      description = "Sex (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained; no estimate reported."
    ),
    RACE_ASIAN = list(
      description = "Race (Asian, White, Black, other were recorded).",
      units = "(binary)",
      type = "binary",
      notes = "Race was screened and not retained; no estimate reported."
    ),
    POT = list(
      description = "Baseline serum potassium.",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained; no estimate reported."
    ),
    CALCIUM = list(
      description = "Baseline serum calcium.",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained; no estimate reported. Documentation only; no canonical calcium column is registered."
    ),
    MAGNESIUM = list(
      description = "Baseline serum magnesium.",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained; no estimate reported. Documentation only; no canonical magnesium column is registered."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 2L,
    n_observations = 769L,
    age_range = "mean 63.8 years (SD 13.4); J101 66.1 (10.2), U102 57.0 (18.2)",
    weight_range = "not reported",
    sex_female_pct = 50,
    race_ethnicity = c(Asian = 49, White = 41, Black = 8, Other = 2),
    disease_state = paste(
      "Relapsed/refractory B- or T-cell non-Hodgkin lymphoma (J101, n = 72)",
      "or relapsed/refractory acute myeloid / acute lymphoblastic leukemia",
      "(U102, n = 28)."
    ),
    dose_range = paste(
      "Oral valemetostat once daily in continuous 28-day cycles: 150, 200,",
      "250 or 300 mg in J101 and 100, 150, 250, 500 or 700 mg in U102;",
      "56% of patients received 200 mg."
    ),
    regions = "Japan and USA (J101); USA (U102)",
    notes = paste(
      "Studies NCT02732275 (DS3201-A-J101) and NCT03110354 (DS3201-A-U102).",
      "769 triplicate-averaged ECGs with a PK sample within 0-15 min after the last ECG (568 from 72 J101 patients, 201 from 28 U102 patients), mostly Cycle 1 Days 1, 2, 8/15, 16 and 22 plus U102 Cycle 2 Day 1 and end of treatment.",
      "Baseline means (SD): heart rate 77.2 (15.4) bpm, QT 387 (38.2) ms, QTcF 417 (24.4) ms, calcium 2.29 (0.15), potassium 4.11 (0.42), magnesium 0.81 (0.11) mmol/L (Table 1).",
      "No placebo control. Fitted with lme4 / nlme in R 4.1.2; shrinkage 19.3% (intercept), 39.6% (slope), 5.35% (residual)."
    )
  )

  ini({
    # Tachibana 2025 Methods 'Data Analysis' equation, applied to DeltaQTcP:
    #   DeltaQTcP_i,t = (theta_int + eta_int,i) + (theta_slp + eta_slp,i) * Conc_i,t
    #                 + theta_b,l * (QTcP_i,0 - QTcP_0) + eps_i,t
    # All values from Table 2, DeltaQTcP columns. Every parameter is on the
    # linear ms scale and both random effects are additive.

    e0 <- -0.00395
    label("Intercept on DeltaQTcP at zero concentration and mean baseline QTcP (ms)")
    # Table 2 'Intercept, ms' = -0.00395 (SE 0.768), RSE 19,443%

    slope <- 0.000938
    label("Concentration-DeltaQTcP slope (ms per ng/mL)")
    # Table 2 'Slope in concentration' = 0.000938 (SE 0.000485), 95% CI 0.0000224 to 0.00191

    e_qtc_bl_e0 <- -0.287
    label("Effect of centered baseline QTcP on the intercept (ms per ms)")
    # Table 2 'Baseline QTcF/QTcP effect' = -0.287 (SE 0.0313), RSE 10.9%

    # Table 2 reports the random effects as SDs plus a correlation:
    #   SD(eta_int) = 5.35 ms, SD(eta_slp) = 0.00192 ms/(ng/mL), corr = 0.451.
    #   var(eta_int) = 5.35^2 = 28.6225
    #   cov          = 0.451 * 5.35 * 0.00192 = 0.004632672
    #   var(eta_slp) = 0.00192^2 = 3.6864e-06
    etae0 + etaslope ~ c(
      28.6225,
      0.004632672, 3.6864e-06
    )

    addSd <- 9.94
    label("Additive residual error SD on DeltaQTcP in study J101 (ms)")
    # Table 2 'SD (sigma) for study J101, ms' = 9.94 (SE 0.276)

    e_study_u102_addSd <- 1.20
    label("Multiplier on the residual SD for study U102 (fold)")
    # Table 2 'SD multiplier for study U102' = 1.20
  })

  model({
    # Baseline QTcP centered on 429 ms, derived from the Table 1 mean QTcF
    # and heart rate (the paper does not print the mean baseline QTcP).
    qtc_bl_ref <- 429
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    e0_i <- e0 + etae0 + e_qtc_bl_e0 * qtc_bl_centered
    slope_i <- slope + etaslope

    # Study-specific residual SD
    addSd_i <- addSd * (1 + (e_study_u102_addSd - 1) * STUDY_U102)

    # Change from baseline in QTcP (ms)
    QTcP <- e0_i + slope_i * CP_VALEMETOSTAT_NGML
    QTcP ~ add(addSd_i)
  })
}
