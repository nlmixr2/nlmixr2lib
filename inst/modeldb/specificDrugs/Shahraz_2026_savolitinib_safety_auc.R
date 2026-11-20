Shahraz_2026_savolitinib_safety_auc <- function() {
  description <- paste(
    "Exposure-safety logistic regressions for savolitinib in adults with",
    "locally advanced or metastatic solid tumours (Shahraz 2026, 998",
    "patients pooled from 10 studies). 23 independent binomial",
    "logistic regressions, one per adverse-event endpoint, each of the",
    "probability of the event on the individual steady-state savolitinib AUC over 24 h (AUC24,ss)",
    "supplied as data. There is no PK layer and no ODE; the outputs are",
    "algebraic probabilities and no observation endpoint is declared. The",
    "sibling model Shahraz_2026_savolitinib_safety_cmax holds the same",
    "endpoints regressed on the other exposure metric."
  )
  reference <- paste(
    "Shahraz A, Derippe T, Kode KK, Allmond J, Nelson S, Igwegbe I,",
    "Sawant-Basak A, Vishwanathan K, Miah K. Population Pharmacokinetics and",
    "Exposure-Response Analyses of Savolitinib Alone or in Combination With",
    "Osimertinib in Patients With Locally Advanced or Metastatic Solid",
    "Tumors. Clin Transl Sci. 2026;19:e70700. doi:10.1111/cts.70700.",
    "PMCID PMC13504015.",
    sep = " "
  )
  vignette <- "Shahraz_2026_savolitinib"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as a covariate column)",
    concentration = "probability of each adverse-event endpoint (0-1)"
  )

  covariateData <- list(
    AUC_SAVOLITINIB = list(
      description = "Individual steady-state savolitinib AUC over 24 h (AUC24,ss), supplied as data (this model has no PK layer).",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Shahraz 2026 Methods: AUC24,ss from the individual empirical Bayes estimates of the population PK model (modellib('Shahraz_2026_savolitinib')) at each patient's most prevalent dose and regimen; for twice-daily dosing it is the AUC over two dosing intervals. From that model, AUC_SAVOLITINIB = daily_dose_mg / cl * 1000. Enters each logit linearly and uncentred, except prob_edema_grade3 (natural log). Table S4 medians: 8829 (300 mg QD), 16220 (300 mg BID), 18354 (600 mg QD) ng*h/mL.",
      source_name = "AUC24ss"
    )
  )

  covariatesDataExcluded <- list(
    CMAX = list(
      description = "Individual steady-state savolitinib Cmax,ss",
      units = "ng/mL",
      type = "continuous",
      notes = "Each endpoint was also regressed on Cmax,ss separately; those fits are the sibling model Shahraz_2026_savolitinib_safety_cmax."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 998L,
    n_studies = 10L,
    age_range = "23-92 years",
    weight_range = "36.0-177 kg",
    sex_female_pct = 48.5,
    race_ethnicity = c(White = 40.5, `Black or African American` = 1.4, Asian = 55.7, Other = 2.4),
    disease_state = "Locally advanced or metastatic solid tumours (EGFR-mutated MET-altered NSCLC, papillary renal cell carcinoma, NSCLC with MET exon 14 alterations, other advanced solid tumours)",
    dose_range = "Savolitinib 100-1000 mg once daily or 100-600 mg twice daily (Table S4), alone or with osimertinib 80 mg once daily",
    regions = "Multinational; 55.7% Asian",
    notes = "The exposure-safety analysis set is the population PK analysis set (Shahraz 2026 Methods). Each endpoint was regressed separately on AUC24,ss and on Cmax,ss with a two-sided alpha of 0.001; the figure panels report the number of patients with and without each event (for example 762 / 236 for any related TEAE, Figure 5a)."
  )

  ini({
    logit_ref_teae_rel <- 0.9796
    label("Logit of the probability of any treatment-emergent adverse event causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Back-solved by least squares from the Results text: probabilities 0.76 [0.74-0.79] (300 mg BID), 0.75 [0.74-0.76] (300 mg QD), 0.76 [0.75-0.80] (600 mg QD) at the median [5th-95th percentile] AUC24,ss of Table S4; Figure 5a prints no coefficients
    e_auc_teae_rel <- 1.063e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Back-solved by least squares from the Results text: probabilities 0.76 [0.74-0.79] (300 mg BID), 0.75 [0.74-0.76] (300 mg QD), 0.76 [0.75-0.80] (600 mg QD) at the median [5th-95th percentile] AUC24,ss of Table S4; Figure 5a prints no coefficients
    logit_ref_teae_rel_g3 <- -1.4294
    label("Logit of the probability of any CTCAE grade >= 3 treatment-emergent adverse event causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Back-solved by least squares from the Results text: probabilities 0.24 [0.21-0.29], 0.22 [0.21-0.24], 0.25 [0.22-0.31] at the median [5th-95th percentile] AUC24,ss of Table S4; Figure 5c prints no coefficients
    e_auc_teae_rel_g3 <- 1.675e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Back-solved by least squares from the Results text: probabilities 0.24 [0.21-0.29], 0.22 [0.21-0.24], 0.25 [0.22-0.31] at the median [5th-95th percentile] AUC24,ss of Table S4; Figure 5c prints no coefficients
    logit_ref_dosemod <- -0.844
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose modification at AUC24,ss = 0 (unitless logit)") # Figure S11a footnote 'Predictor=-0.844+ AUC24ss*-1.24e-05' (TAEDMSRS)
    e_auc_dosemod <- -1.24e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S11a footnote 'Predictor=-0.844+ AUC24ss*-1.24e-05' (TAEDMSRS)
    logit_ref_dose_red <- -1.22
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose reduction at AUC24,ss = 0 (unitless logit)") # Figure S12a footnote 'Predictor=-1.22+ AUC24ss*-3.51e-05' (TAEDRSRS)
    e_auc_dose_red <- -3.51e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S12a footnote 'Predictor=-1.22+ AUC24ss*-3.51e-05' (TAEDRSRS)
    logit_ref_dose_int <- -1.5
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose interruption at AUC24,ss = 0 (unitless logit)") # Figure S13a footnote 'Predictor=-1.5+ AUC24ss*4.2e-07' (TAEDISRS)
    e_auc_dose_int <- 4.2e-07
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S13a footnote 'Predictor=-1.5+ AUC24ss*4.2e-07' (TAEDISRS)
    logit_ref_discont <- -2.19
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib discontinuation at AUC24,ss = 0 (unitless logit)") # Figure S14a footnote 'Predictor=-2.19+ AUC24ss*1.98e-05' (TAEDDSRS)
    e_auc_discont <- 1.98e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S14a footnote 'Predictor=-2.19+ AUC24ss*1.98e-05' (TAEDDSRS)
    logit_ref_heptox <- -1.17
    label("Logit of the probability of grade >= 1 hepatotoxicity (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S15a footnote 'Predictor=-1.17+ AUC24ss*2.57e-05' (HEPTOXG1)
    e_auc_heptox <- 2.57e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S15a footnote 'Predictor=-1.17+ AUC24ss*2.57e-05' (HEPTOXG1)
    logit_ref_heptox_g3 <- -2.44
    label("Logit of the probability of grade >= 3 hepatotoxicity (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S16a footnote 'Predictor=-2.44+ AUC24ss*2.49e-05' (HEPTOXG3)
    e_auc_heptox_g3 <- 2.49e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S16a footnote 'Predictor=-2.44+ AUC24ss*2.49e-05' (HEPTOXG3)
    logit_ref_hypsen <- -2.93
    label("Logit of the probability of grade >= 1 hypersensitivity (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S17a footnote 'Predictor=-2.93+ AUC24ss*-1.04e-05' (HYPSENG1)
    e_auc_hypsen <- -1.04e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S17a footnote 'Predictor=-2.93+ AUC24ss*-1.04e-05' (HYPSENG1)
    logit_ref_hypsen_g3 <- -3.53
    label("Logit of the probability of grade >= 3 hypersensitivity (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S18a footnote 'Predictor=-3.53+ AUC24ss*-2.83e-06' (HYPSENG3)
    e_auc_hypsen_g3 <- -2.83e-06
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S18a footnote 'Predictor=-3.53+ AUC24ss*-2.83e-06' (HYPSENG3)
    logit_ref_pyrexia <- -1.73
    label("Logit of the probability of grade >= 1 pyrexia (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S19a footnote 'Predictor=-1.73+ AUC24ss*6.07e-06' (PYRXG1)
    e_auc_pyrexia <- 6.07e-06
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S19a footnote 'Predictor=-1.73+ AUC24ss*6.07e-06' (PYRXG1)
    logit_ref_pyrexia_g3 <- -5
    label("Logit of the probability of grade >= 3 pyrexia (savolitinib adverse event of special interest) at AUC24,ss = 0 (unitless logit)") # Figure S20a footnote 'Predictor=-5+ AUC24ss*1.65e-05' (PYRXG3)
    e_auc_pyrexia_g3 <- 1.65e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S20a footnote 'Predictor=-5+ AUC24ss*1.65e-05' (PYRXG3)
    logit_ref_alt <- -2.47
    label("Logit of the probability of grade >= 1 ALT elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S21a footnote 'Predictor=-2.47+ AUC24ss*1.14e-05' (ALTAEG1S)
    e_auc_alt <- 1.14e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S21a footnote 'Predictor=-2.47+ AUC24ss*1.14e-05' (ALTAEG1S)
    logit_ref_alt_g3 <- -3.26
    label("Logit of the probability of grade >= 3 ALT elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S22a footnote 'Predictor=-3.26+ AUC24ss*-6.51e-08' (ALTAEG3S)
    e_auc_alt_g3 <- -6.51e-08
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S22a footnote 'Predictor=-3.26+ AUC24ss*-6.51e-08' (ALTAEG3S)
    logit_ref_ast <- -2.55
    label("Logit of the probability of grade >= 1 AST elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S23a footnote 'Predictor=-2.55+ AUC24ss*1.11e-05' (ASTAEG1S)
    e_auc_ast <- 1.11e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S23a footnote 'Predictor=-2.55+ AUC24ss*1.11e-05' (ASTAEG1S)
    logit_ref_ast_g3 <- -3.75
    label("Logit of the probability of grade >= 3 AST elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S24a footnote 'Predictor=-3.75+ AUC24ss*1.3e-05' (ASTAEG3S)
    e_auc_ast_g3 <- 1.3e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S24a footnote 'Predictor=-3.75+ AUC24ss*1.3e-05' (ASTAEG3S)
    logit_ref_alp <- -3.87
    label("Logit of the probability of grade >= 1 alkaline phosphatase elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S25a footnote 'Predictor=-3.87+ AUC24ss*1.15e-05' (ALPAEG1S)
    e_auc_alp <- 1.15e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S25a footnote 'Predictor=-3.87+ AUC24ss*1.15e-05' (ALPAEG1S)
    logit_ref_bili <- -4.58
    label("Logit of the probability of grade >= 1 total bilirubin elevation causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S26a footnote 'Predictor=-4.58+ AUC24ss*1.75e-05' (BILAEG1S)
    e_auc_bili <- 1.75e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S26a footnote 'Predictor=-4.58+ AUC24ss*1.75e-05' (BILAEG1S)
    logit_ref_edema <- -0.561
    label("Logit of the probability of grade >= 1 oedema causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S27a footnote 'Predictor=-0.561+ AUC24ss*-1.92e-05' (OEDAEG1S)
    e_auc_edema <- -1.92e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S27a footnote 'Predictor=-0.561+ AUC24ss*-1.92e-05' (OEDAEG1S)
    logit_ref_edema_g3 <- 1.4769
    label("Logit of the probability of grade >= 3 oedema causally related to savolitinib at AUC24,ss = 1 ng*h/mL (unitless logit)") # Digitised by the maintainers from the plotted curve of Figure S28a, whose footnote is missing; logistic in log(AUC24,ss) (RMS 0.0015 on the probability scale over 422 curve pixels), the form of every plotted curve in Figures S11-S31
    e_logauc_edema_g3 <- -0.52308
    label("Change in that logit per unit natural log of savolitinib AUC24,ss in ng*h/mL (unitless logit)") # Digitised by the maintainers from the plotted curve of Figure S28a, whose footnote is missing; logistic in log(AUC24,ss) (RMS 0.0015 on the probability scale over 422 curve pixels), the form of every plotted curve in Figures S11-S31
    logit_ref_hypoalb <- -3.02
    label("Logit of the probability of grade >= 1 hypoalbuminaemia causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S29a footnote 'Predictor=-3.02+ AUC24ss*5.54e-06' (ALBAEG1S)
    e_auc_hypoalb <- 5.54e-06
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S29a footnote 'Predictor=-3.02+ AUC24ss*5.54e-06' (ALBAEG1S)
    logit_ref_diarrhea <- -2.01
    label("Logit of the probability of grade >= 1 diarrhoea causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S30a footnote 'Predictor=-2.01+ AUC24ss*-1.5e-06' (DIAAEG1S)
    e_auc_diarrhea <- -1.5e-06
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S30a footnote 'Predictor=-2.01+ AUC24ss*-1.5e-06' (DIAAEG1S)
    logit_ref_vte <- -4.67
    label("Logit of the probability of grade >= 1 venous thromboembolism causally related to savolitinib at AUC24,ss = 0 (unitless logit)") # Figure S31a footnote 'Predictor=-4.67+ AUC24ss*2.21e-05' (VTEAEG1S)
    e_auc_vte <- 2.21e-05
    label("Change in that logit per 1 ng*h/mL of savolitinib AUC24,ss (logit per ng*h/mL)") # Figure S31a footnote 'Predictor=-4.67+ AUC24ss*2.21e-05' (VTEAEG1S)
  })

  model({
    # Each output is an independent binomial-logit regression (Figures 5
    # and S11-S31: 'Distribution of Y axis = Binomial; Link function =
    # Logit'), linear and uncentred in the exposure unless written with log().
    prob_teae_drug_related <- expit(logit_ref_teae_rel + e_auc_teae_rel * AUC_SAVOLITINIB)
    prob_teae_drug_related_grade3 <- expit(logit_ref_teae_rel_g3 + e_auc_teae_rel_g3 * AUC_SAVOLITINIB)
    prob_teae_dosemod <- expit(logit_ref_dosemod + e_auc_dosemod * AUC_SAVOLITINIB)
    prob_dose_reduction <- expit(logit_ref_dose_red + e_auc_dose_red * AUC_SAVOLITINIB)
    prob_dose_interruption <- expit(logit_ref_dose_int + e_auc_dose_int * AUC_SAVOLITINIB)
    prob_treatment_discontinuation <- expit(logit_ref_discont + e_auc_discont * AUC_SAVOLITINIB)
    prob_hepatotoxicity <- expit(logit_ref_heptox + e_auc_heptox * AUC_SAVOLITINIB)
    prob_hepatotoxicity_grade3 <- expit(logit_ref_heptox_g3 + e_auc_heptox_g3 * AUC_SAVOLITINIB)
    prob_hypersensitivity <- expit(logit_ref_hypsen + e_auc_hypsen * AUC_SAVOLITINIB)
    prob_hypersensitivity_grade3 <- expit(logit_ref_hypsen_g3 + e_auc_hypsen_g3 * AUC_SAVOLITINIB)
    prob_pyrexia <- expit(logit_ref_pyrexia + e_auc_pyrexia * AUC_SAVOLITINIB)
    prob_pyrexia_grade3 <- expit(logit_ref_pyrexia_g3 + e_auc_pyrexia_g3 * AUC_SAVOLITINIB)
    prob_alt_increase <- expit(logit_ref_alt + e_auc_alt * AUC_SAVOLITINIB)
    prob_alt_increase_grade3 <- expit(logit_ref_alt_g3 + e_auc_alt_g3 * AUC_SAVOLITINIB)
    prob_ast_increase <- expit(logit_ref_ast + e_auc_ast * AUC_SAVOLITINIB)
    prob_ast_increase_grade3 <- expit(logit_ref_ast_g3 + e_auc_ast_g3 * AUC_SAVOLITINIB)
    prob_alp_increase <- expit(logit_ref_alp + e_auc_alp * AUC_SAVOLITINIB)
    prob_bilirubin_increase <- expit(logit_ref_bili + e_auc_bili * AUC_SAVOLITINIB)
    prob_edema <- expit(logit_ref_edema + e_auc_edema * AUC_SAVOLITINIB)
    prob_edema_grade3 <- expit(logit_ref_edema_g3 + e_logauc_edema_g3 * log(AUC_SAVOLITINIB))
    prob_hypoalbuminemia <- expit(logit_ref_hypoalb + e_auc_hypoalb * AUC_SAVOLITINIB)
    prob_diarrhea <- expit(logit_ref_diarrhea + e_auc_diarrhea * AUC_SAVOLITINIB)
    prob_venous_thromboembolism <- expit(logit_ref_vte + e_auc_vte * AUC_SAVOLITINIB)
  })
}
