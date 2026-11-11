Tachibana_2025_valemetostat_QTcF <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-QTc model for the",
    "EZH1/EZH2 dual inhibitor valemetostat in 100 adults with",
    "relapsed/refractory non-Hodgkin lymphoma (phase 1 study DS3201-A-J101)",
    "or acute myeloid / acute lymphoblastic leukemia (phase 1 study",
    "DS3201-A-U102) receiving 100-700 mg once daily. The endpoint is the",
    "change from baseline in the Fridericia-corrected QT interval",
    "(DeltaQTcF, ms): DeltaQTcF = (e0 + etae0) + (slope + etaslope) *",
    "CP_VALEMETOSTAT_NGML + e_qtc_bl_e0 * (QTC_BL - 417) + eps, with",
    "intercept 1.05 ms, slope 0.000746 ms/(ng/mL) (95% CI includes 0) and",
    "baseline-QTcF coefficient -0.266 ms/ms. Additive, correlated random",
    "effects on the intercept (SD 5.52 ms) and slope (SD 0.00170",
    "ms/(ng/mL), correlation 0.0616). Additive residual SD 9.99 ms in",
    "J101, multiplied by 1.21 in U102 (STUDY_U102). No placebo arm and no",
    "nominal-time term, so the typical prediction at the cohort-mean",
    "baseline is e0 + slope * C (2.41 ms at the 1830 ng/mL steady-state",
    "Cmax of 200 mg once daily). PD-only model: the valemetostat plasma",
    "concentration is supplied as the time-varying covariate",
    "CP_VALEMETOSTAT_NGML (ng/mL). The companion population-corrected",
    "model is Tachibana_2025_valemetostat_QTcP.R.",
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
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF, in ms; driving covariate CP_VALEMETOSTAT_NGML is in ng/mL)"
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
      description = "Patient's pre-dose baseline Fridericia-corrected QT interval (QTcF_i,0), the triplicate-averaged ECG taken pre-dose on Cycle 1 Day 1. Time-fixed per patient. Enters the intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 417); set QTC_BL = 417 ms for the typical patient.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fridericia correction QTcF = QT / RR^(1/3), QT in ms and RR in s (Methods 'Study Design').",
        "Centering reference: the equation centers on QTcF_0, 'the mean of QTcF_i,0 over all patients' (Methods 'Data Analysis'); Table 1 reports that mean (SD) as 417 (24.4) ms for the 100 patients of the primary analysis set (417 (26.9) in J101, 416 (17.0) in U102), and Table S1 gives the same 417 ms for the baseline visit.",
        "Eligibility excluded QTcF > 450 ms at screening in both studies."
      ),
      source_name = "QTcF_i,0 (baseline QTcF)"
    ),
    STUDY_U102 = list(
      description = "Study indicator: 1 = record from DS3201-A-U102 (relapsed/refractory AML or ALL), 0 = record from DS3201-A-J101 (relapsed/refractory B- or T-cell NHL).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (DS3201-A-J101)",
      notes = paste(
        "Selects the residual SD only: addSd (9.99 ms) in J101 and addSd * e_study_u102_addSd (9.99 * 1.21 = 12.09 ms) in U102 (Table 2 'SD multiplier for study U102').",
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
      "No placebo control. Fitted with lme4 / nlme in R 4.1.2; shrinkage 21.6% (intercept), 57.5% (slope), 5.28% (residual)."
    )
  )

  ini({
    # Tachibana 2025 Methods 'Data Analysis' equation:
    #   DeltaQTcF_i,t = (theta_int + eta_int,i) + (theta_slp + eta_slp,i) * Conc_i,t
    #                 + theta_b,l * (QTcF_i,0 - QTcF_0) + eps_i,t
    # All values from Table 2, DeltaQTcF columns. Every parameter is on the
    # linear ms scale and both random effects are additive.

    e0 <- 1.05
    label("Intercept on DeltaQTcF at zero concentration and mean baseline QTcF (ms)")
    # Table 2 'Intercept, ms' = 1.05 (SE 0.779), RSE 74.2%

    slope <- 0.000746
    label("Concentration-DeltaQTcF slope (ms per ng/mL; 95% CI includes 0)")
    # Table 2 'Slope in concentration' = 0.000746 (SE 0.000460), 95% CI -0.000132 to 0.00164

    e_qtc_bl_e0 <- -0.266
    label("Effect of centered baseline QTcF on the intercept (ms per ms)")
    # Table 2 'Baseline QTcF/QTcP effect' = -0.266 (SE 0.0289), RSE 10.9%

    # Table 2 reports the random effects as SDs plus a correlation:
    #   SD(eta_int) = 5.52 ms, SD(eta_slp) = 0.00170 ms/(ng/mL), corr = 0.0616.
    # The slope SD is printed as 0.00710, but its SE 0.000518 with RSE 30.5%
    # gives 0.000518 / 0.305 = 0.00170, and its 95% CI 0.000865-0.00266 is
    # centered near 0.0017 and excludes 0.00710; 0.00170 is used.
    #   var(eta_int) = 5.52^2 = 30.4704
    #   cov          = 0.0616 * 5.52 * 0.00170 = 0.000578054
    #   var(eta_slp) = 0.00170^2 = 2.89e-06
    etae0 + etaslope ~ c(
      30.4704,
      0.000578054, 2.89e-06
    )

    addSd <- 9.99
    label("Additive residual error SD on DeltaQTcF in study J101 (ms)")
    # Table 2 'SD (sigma) for study J101, ms' = 9.99 (SE 0.278)

    e_study_u102_addSd <- 1.21
    label("Multiplier on the residual SD for study U102 (fold)")
    # Table 2 'SD multiplier for study U102' = 1.21
  })

  model({
    # Baseline QTcF centered on the published cohort mean (Table 1, 417 ms).
    qtc_bl_ref <- 417
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    e0_i <- e0 + etae0 + e_qtc_bl_e0 * qtc_bl_centered
    slope_i <- slope + etaslope

    # Study-specific residual SD
    addSd_i <- addSd * (1 + (e_study_u102_addSd - 1) * STUDY_U102)

    # Change from baseline in QTcF (ms)
    QTcF <- e0_i + slope_i * CP_VALEMETOSTAT_NGML
    QTcF ~ add(addSd_i)
  })
}
