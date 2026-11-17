Colin_2026_pediatricT2DM_hba1c_mbma <- function() {
  description <- paste(
    "MBMA (individual-participant-data longitudinal meta-analysis). HbA1c",
    "disease-progression model for paediatric patients aged 10 years and",
    "older with type 2 diabetes mellitus, fit to the placebo and active arms",
    "of six randomised trials submitted to the European Medicines Agency for",
    "seven antihyperglycaemic drugs (linagliptin, empagliflozin, alogliptin,",
    "exenatide, sitagliptin, saxagliptin, dapagliflozin). HbA1c is",
    "S0 * (1 + alpha * t) minus a symptomatic effect S0 * expit(theta + dTRT)",
    "* (1 - exp(-log(2) / t1/2 * t)) shared by placebo (dTRT = 0) and active",
    "treatment, where dTRT is a drug-specific logit-scale shift. Background",
    "insulin raises baseline HbA1c and the progression rate; time since",
    "diagnosis raises baseline HbA1c. Correlated between-subject",
    "variability on baseline, progression rate (Box-Cox transformed) and",
    "placebo logit, with between-study variability on baseline and",
    "progression rate nested above it. Log-additive residual error."
  )
  reference <- paste(
    "Colin PJ, Manolis E, Koomen JV.",
    "Model-Informed Approaches for Pediatric Type 2 Diabetes Evidence",
    "Generation: Insights from a Model-Based Longitudinal HbA1c Meta-Analysis",
    "from EMA Marketing Authorization Applications.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70468."
  )
  vignette <- "Colin_2026_pediatricT2DM_hba1c_mbma"
  units <- list(
    time = "day",
    dosing = "none (treatment enters through the TRT_<INN> arm indicators; the model has no dose events)",
    concentration = "HbA1c in % (NGSP)"
  )

  # The study-level etas are named eta_study_* (MBMA convention) and sit on
  # top of the subject-level etas, as in the NONMEM $LEVEL nesting.
  paper_specific_etas <- c("eta_study_lrbase", "eta_study_lslope")

  covariateData <- list(
    CONMED_INSULIN = list(
      description = "Background insulin therapy at study entry (1 = insulin or metformin + insulin, 0 = none or metformin only)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no background insulin)",
      notes = paste(
        "Time-fixed. Derived from the four-level background medication",
        "('none', 'insulin', 'metformin', 'metformin + insulin') as INS01 =",
        "1 when INS > 0 in the Text S1 control stream. Multiplies baseline",
        "HbA1c by exp(0.108) = 1.11 and the progression rate by",
        "exp(0.545) = 1.72."
      ),
      source_name = "INS (INS01 in Text S1)"
    ),
    T_DIAG_DIAB = list(
      description = "Time since type 2 diabetes diagnosis at study entry",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Enters baseline HbA1c as exp(0.0121 * (T_DIAG_DIAB - 2.5)),",
        "i.e. a 1.012-fold change per year from a 2.5-year reference. As in",
        "the Text S1 control stream (DIA = 2.5 unless DIAGNOSIS > 0), a value",
        "of 0 or below is treated as missing and replaced by the 2.5-year",
        "reference; code a missing duration as 0. Duration was unrecorded",
        "for the whole sitagliptin and dapagliflozin studies (Table 1)."
      ),
      source_name = "DIAGNOSIS (DIA in Text S1)"
    ),
    TRT_LINAGLIPTIN = list(
      description = "Linagliptin treatment-arm indicator (1 = randomised to linagliptin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 1 in Text S1 (THETA9). The seven TRT_ indicators are mutually exclusive and all 0 on placebo.",
      source_name = "TRT == 1"
    ),
    TRT_EMPAGLIFLOZIN = list(
      description = "Empagliflozin treatment-arm indicator (1 = randomised to empagliflozin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 2 in Text S1 (THETA10).",
      source_name = "TRT == 2"
    ),
    TRT_ALOGLIPTIN = list(
      description = "Alogliptin treatment-arm indicator (1 = randomised to alogliptin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 3 in Text S1 (THETA11).",
      source_name = "TRT == 3"
    ),
    TRT_EXENATIDE = list(
      description = "Exenatide (once-weekly 2 mg) treatment-arm indicator (1 = randomised to exenatide, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 4 in Text S1 (THETA12).",
      source_name = "TRT == 4"
    ),
    TRT_SITAGLIPTIN = list(
      description = "Sitagliptin treatment-arm indicator (1 = randomised to sitagliptin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 5 in Text S1 (THETA13).",
      source_name = "TRT == 5"
    ),
    TRT_SAXAGLIPTIN = list(
      description = "Saxagliptin treatment-arm indicator (1 = randomised to saxagliptin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 6 in Text S1 (THETA14).",
      source_name = "TRT == 6"
    ),
    TRT_DAPAGLIFLOZIN = list(
      description = "Dapagliflozin treatment-arm indicator (1 = randomised to dapagliflozin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or another drug)",
      notes = "TRT = 7 in Text S1 (THETA15).",
      source_name = "TRT == 7"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods, 'Model-based longitudinal meta-analysis') and not retained (Results: 'No other covariates were significant')."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained."
    ),
    CONMED_METFORMIN = list(
      description = "Background metformin therapy",
      units = "(binary)",
      type = "binary",
      notes = "Screened as part of the four-level background medication; only the insulin component was retained."
    )
  )

  compartmentData <- list()

  population <- list(
    species = "human",
    n_subjects = 802L,
    n_studies = 6L,
    age_range = "10-17 years in five studies, 10-24 years in the dapagliflozin study (Table S1); study means 14.0-16.2 years",
    age_mean = "14.6 years (SD 2.16)",
    weight_mean = "90.0 kg (SD 27.2)",
    sex_female_pct = 62,
    disease_state = "Type 2 diabetes mellitus with onset in youth; mean WHO BMI-for-age z-score 2.85 (SD 1.27); median time since diagnosis 2.2 years (IQR 1.0-3.7)",
    dose_range = paste(
      "Placebo or linagliptin 5 mg, empagliflozin 10/25 mg, alogliptin 25 mg,",
      "exenatide 2 mg once weekly, sitagliptin 100 mg, saxagliptin 2.5/5 mg,",
      "dapagliflozin 10 mg (Table S1); 20-26 weeks of treatment"
    ),
    co_medication = "Background metformin 41%, metformin + insulin 31%, insulin 10%, none 18%",
    regions = "Multinational trials submitted to the EMA in support of paediatric extensions of indication",
    notes = paste(
      "809 participants (456 active, 353 placebo) contributed 3295 HbA1c",
      "observations (Table 1, Results 'Data'); the final model excludes seven",
      "influential participants (IGNORE statements of Text S1), leaving 802.",
      "Study 1 randomised to both linagliptin and empagliflozin, so six",
      "studies carry seven drugs. Median follow-up 24 weeks (IQR 19-26)."
    )
  )

  ini({
    # Final estimates from the NONMEM output (Text S1, 'FINAL PARAMETER
    # ESTIMATE'); Table 2 prints the same values back-transformed.
    lrbase <- 2.01; label("Log baseline HbA1c without insulin and with a 2.5-year diabetes history (log %)") # Text S1 TH1 = 2.01; Table 2 beta0 = 7.46 %
    e_conmed_insulin_rbase <- 0.108; label("Log-scale effect of background insulin on baseline HbA1c (unitless)") # Text S1 TH6 = 0.108; Table 2 beta1 = 1.11
    e_t_diag_diab_rbase <- 0.0121; label("Log-scale effect per year of diabetes duration on baseline HbA1c (1/year)") # Text S1 TH8 = 0.0121; Table 2 beta2 = 1.012

    lslope <- -2.05; label("Log disease-progression rate alpha without insulin (log 1/year)") # Text S1 TH2 = -2.05; Table 2 beta3 = 0.13 /year
    e_conmed_insulin_slope <- 0.545; label("Log-scale effect of background insulin on the progression rate (unitless)") # Text S1 TH7 = 0.545; Table 2 beta4 = 1.72
    boxcox_lslope <- -0.315; label("Box-Cox shape of the progression-rate random effect (unitless)") # Text S1 TH16 = -0.315; Table 2 beta6 = -0.31

    logitpbo <- -2.61; label("Logit of the maximal placebo response as a fraction of baseline HbA1c (logit)") # Text S1 TH3 = -2.61; Table 2 theta = -2.61 (0.07 linear)
    lthalf_onset <- 3.87; label("Log half-life of onset of the placebo and treatment response (log day)") # Text S1 TH4 = 3.87 -> 47.9 day; Table 2 beta5 = 6.85 wks

    e_trt_linagliptin_logitpbo <- 0.265; label("Linagliptin shift of the response logit, placebo-corrected (logit)") # Text S1 TH9 = 0.265; Table 2 dTRT,1 = 0.26
    e_trt_empagliflozin_logitpbo <- 1.03; label("Empagliflozin shift of the response logit, placebo-corrected (logit)") # Text S1 TH10 = 1.03; Table 2 dTRT,2 = 1.03
    e_trt_alogliptin_logitpbo <- 0.532; label("Alogliptin shift of the response logit, placebo-corrected (logit)") # Text S1 TH11 = 0.532; Table 2 dTRT,3 = 0.53
    e_trt_exenatide_logitpbo <- 1.28; label("Exenatide shift of the response logit, placebo-corrected (logit)") # Text S1 TH12 = 1.28; Table 2 dTRT,4 = 1.28
    e_trt_sitagliptin_logitpbo <- 0.377; label("Sitagliptin shift of the response logit, placebo-corrected (logit)") # Text S1 TH13 = 0.377; Table 2 dTRT,5 = 0.38
    e_trt_saxagliptin_logitpbo <- 0.971; label("Saxagliptin shift of the response logit, placebo-corrected (logit)") # Text S1 TH14 = 0.971; Table 2 dTRT,6 = 0.97
    e_trt_dapagliflozin_logitpbo <- 0.954; label("Dapagliflozin shift of the response logit, placebo-corrected (logit)") # Text S1 TH15 = 0.954; Table 2 dTRT,7 = 0.95

    # Between-subject variability (Text S1 OMEGA covariance matrix, ETA1-ETA3).
    etalrbase + etalslope + etalogitpbo ~ c(
      0.0212,
      0.0277, 4.18,
      0.0333, 1.69, 1.91
    ) # Text S1 OMEGA ETA1-ETA3 block; Table 2 omega1 = 14.6 %CV, omega2 = 2.04 SD, omega3 = 1.38 SD

    # Between-study variability (Text S1 OMEGA covariance matrix, ETA4-ETA5).
    eta_study_lrbase + eta_study_lslope ~ c(
      4.77e-4,
      1.85e-3, 0.0120
    ) # Text S1 OMEGA ETA4-ETA5 block; Table 2 omega4 = 2.18 %CV

    expSd <- sqrt(3.42e-3); label("Additive residual SD of log HbA1c (log %)") # Text S1 SIGMA EPS1 = 3.42E-03 (SD 5.85E-02); Table 2 5.85 %CV
  })
  model({
    # Pre-randomisation records do not progress or respond (RTIME in Text S1).
    rtime <- t
    if (t < 0) rtime <- 0

    # Missing (coded 0) diabetes duration takes the 2.5-year reference.
    dia <- 2.5
    if (T_DIAG_DIAB > 0) dia <- T_DIAG_DIAB

    rbase_subj <- exp(lrbase + etalrbase)
    rbase <- rbase_subj * exp(eta_study_lrbase) *
      exp(e_conmed_insulin_rbase * CONMED_INSULIN) *
      exp(e_t_diag_diab_rbase * (dia - 2.5))

    # Box-Cox transform of the summed subject + study progression etas
    # (Text S1 ETA_ALPHA; paper Eq. 5).
    eta_slope_bc <- (exp(boxcox_lslope * (etalslope + eta_study_lslope)) - 1) / boxcox_lslope
    slope <- exp(lslope + eta_slope_bc) * exp(e_conmed_insulin_slope * CONMED_INSULIN)

    # Logit of the maximal symptomatic response; placebo when every TRT_ is 0.
    lgt_resp <- logitpbo + etalogitpbo +
      e_trt_linagliptin_logitpbo * TRT_LINAGLIPTIN +
      e_trt_empagliflozin_logitpbo * TRT_EMPAGLIFLOZIN +
      e_trt_alogliptin_logitpbo * TRT_ALOGLIPTIN +
      e_trt_exenatide_logitpbo * TRT_EXENATIDE +
      e_trt_sitagliptin_logitpbo * TRT_SITAGLIPTIN +
      e_trt_saxagliptin_logitpbo * TRT_SAXAGLIPTIN +
      e_trt_dapagliflozin_logitpbo * TRT_DAPAGLIFLOZIN
    resp_max <- expit(lgt_resp)

    thalf_onset <- exp(lthalf_onset)
    kon <- log(2) / thalf_onset
    onset <- 1 - exp(-kon * rtime)

    # The progression rate is per year; time is in days.
    Hba1c <- rbase * (1 + slope * rtime / 365) - rbase * resp_max * onset
    Hba1c ~ lnorm(expSd)
  })
}
