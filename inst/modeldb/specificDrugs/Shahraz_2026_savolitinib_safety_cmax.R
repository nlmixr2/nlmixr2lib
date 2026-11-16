Shahraz_2026_savolitinib_safety_cmax <- function() {
  description <- paste(
    "Exposure-safety logistic regressions for savolitinib in adults with",
    "locally advanced or metastatic solid tumours (Shahraz 2026, 998",
    "patients pooled from 10 studies). 22 independent binomial",
    "logistic regressions, one per adverse-event endpoint, each of the",
    "probability of the event on the individual steady-state savolitinib maximum concentration (Cmax,ss)",
    "supplied as data. There is no PK layer and no ODE; the outputs are",
    "algebraic probabilities and no observation endpoint is declared. The",
    "sibling model Shahraz_2026_savolitinib_safety_auc holds the same",
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
    CMAX = list(
      description = "Individual steady-state savolitinib maximum plasma concentration (Cmax,ss), supplied as data (this model has no PK layer).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL savolitinib plasma concentration at steady state, from the individual empirical Bayes estimates of the population PK model (modellib('Shahraz_2026_savolitinib')) at each patient's most prevalent dose and regimen (Shahraz 2026 Methods). Enters each logit linearly and uncentred, except prob_treatment_discontinuation (natural log). Table S4 medians: 1376 (300 mg QD), 1346 (300 mg BID), 2444 (600 mg QD) ng/mL.",
      source_name = "Cmax,ss"
    )
  )

  covariatesDataExcluded <- list(
    AUC_SAVOLITINIB = list(
      description = "Individual steady-state savolitinib AUC24,ss",
      units = "ng*h/mL",
      type = "continuous",
      notes = "Each endpoint was also regressed on AUC24,ss separately; those fits are the sibling model Shahraz_2026_savolitinib_safety_auc."
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
    notes = "The exposure-safety analysis set is the population PK analysis set (Shahraz 2026 Methods). Each endpoint was regressed separately on AUC24,ss and on Cmax,ss with a two-sided alpha of 0.001; the figure panels report the number of patients with and without each event (for example 762 / 236 for any related TEAE, Figure 5a). The grade >= 3 pyrexia regression on Cmax,ss (Figure S20b) is absent because the published supplement reprints the Figure S17b hypersensitivity panel in its place."
  )

  ini({
    logit_ref_teae_rel <- 0.9685
    label("Logit of the probability of any treatment-emergent adverse event causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Back-solved by least squares from the Results text: probabilities 0.75 [0.74-0.76], 0.75 [0.74-0.77], 0.77 [0.75-0.79] at the median [5th-95th percentile] Cmax,ss of Table S4; Figure 5b prints no coefficients
    e_cmax_teae_rel <- 9.529e-05
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Back-solved by least squares from the Results text: probabilities 0.75 [0.74-0.76], 0.75 [0.74-0.77], 0.77 [0.75-0.79] at the median [5th-95th percentile] Cmax,ss of Table S4; Figure 5b prints no coefficients
    logit_ref_teae_rel_g3 <- -1.6384
    label("Logit of the probability of any CTCAE grade >= 3 treatment-emergent adverse event causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Back-solved by least squares from the Results text: probabilities 0.21 [0.18-0.25], 0.22 [0.20-0.25], 0.26 [0.21-0.33] at the median [5th-95th percentile] Cmax,ss of Table S4; Figure 5d prints no coefficients
    e_cmax_teae_rel_g3 <- 2.428e-04
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Back-solved by least squares from the Results text: probabilities 0.21 [0.18-0.25], 0.22 [0.20-0.25], 0.26 [0.21-0.33] at the median [5th-95th percentile] Cmax,ss of Table S4; Figure 5d prints no coefficients
    logit_ref_dosemod <- -0.66
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose modification at Cmax,ss = 0 (unitless logit)") # Figure S11b footnote 'Predictor=-0.66+ Cmaxss*-0.000204' (TAEDMSRS)
    e_cmax_dosemod <- -0.000204
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S11b footnote 'Predictor=-0.66+ Cmaxss*-0.000204' (TAEDMSRS)
    logit_ref_dose_red <- -1.02
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose reduction at Cmax,ss = 0 (unitless logit)") # Figure S12b footnote 'Predictor=-1.02+ Cmaxss*-0.00041' (TAEDRSRS)
    e_cmax_dose_red <- -0.00041
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S12b footnote 'Predictor=-1.02+ Cmaxss*-0.00041' (TAEDRSRS)
    logit_ref_dose_int <- -1.39
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib dose interruption at Cmax,ss = 0 (unitless logit)") # Figure S13b footnote 'Predictor=-1.39+ Cmaxss*-5.07e-05' (TAEDISRS)
    e_cmax_dose_int <- -5.07e-05
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S13b footnote 'Predictor=-1.39+ Cmaxss*-5.07e-05' (TAEDISRS)
    logit_ref_discont <- -7.431
    label("Logit of the probability of any TEAE causally related to savolitinib leading to savolitinib discontinuation at Cmax,ss = 1 ng/mL (unitless logit)") # Digitised by the maintainers from the plotted curve of Figure S14b, whose footnote is missing; logistic in log(Cmax,ss) (RMS 0.001 on the probability scale over 421 curve pixels), the form of every plotted curve in Figures S11-S31
    e_logcmax_discont <- 0.7433
    label("Change in that logit per unit natural log of savolitinib Cmax,ss in ng/mL (unitless logit)") # Digitised by the maintainers from the plotted curve of Figure S14b, whose footnote is missing; logistic in log(Cmax,ss) (RMS 0.001 on the probability scale over 421 curve pixels), the form of every plotted curve in Figures S11-S31
    logit_ref_heptox <- -1.12
    label("Logit of the probability of grade >= 1 hepatotoxicity (savolitinib adverse event of special interest) at Cmax,ss = 0 (unitless logit)") # Figure S15b footnote 'Predictor=-1.12+ Cmaxss*0.000203' (HEPTOXG1)
    e_cmax_heptox <- 0.000203
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S15b footnote 'Predictor=-1.12+ Cmaxss*0.000203' (HEPTOXG1)
    logit_ref_heptox_g3 <- -2.39
    label("Logit of the probability of grade >= 3 hepatotoxicity (savolitinib adverse event of special interest) at Cmax,ss = 0 (unitless logit)") # Figure S16b footnote 'Predictor=-2.39+ Cmaxss*0.000204' (HEPTOXG3)
    e_cmax_heptox_g3 <- 0.000204
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S16b footnote 'Predictor=-2.39+ Cmaxss*0.000204' (HEPTOXG3)
    logit_ref_hypsen <- -3.65
    label("Logit of the probability of grade >= 1 hypersensitivity (savolitinib adverse event of special interest) at Cmax,ss = 0 (unitless logit)") # Figure S17b footnote 'Predictor=-3.65+ Cmaxss*0.000264' (HYPSENG1)
    e_cmax_hypsen <- 0.000264
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S17b footnote 'Predictor=-3.65+ Cmaxss*0.000264' (HYPSENG1)
    logit_ref_hypsen_g3 <- -4.37
    label("Logit of the probability of grade >= 3 hypersensitivity (savolitinib adverse event of special interest) at Cmax,ss = 0 (unitless logit)") # Figure S18b footnote 'Predictor=-4.37+ Cmaxss*0.000369' (HYPSENG3)
    e_cmax_hypsen_g3 <- 0.000369
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S18b footnote 'Predictor=-4.37+ Cmaxss*0.000369' (HYPSENG3)
    logit_ref_pyrexia <- -1.77
    label("Logit of the probability of grade >= 1 pyrexia (savolitinib adverse event of special interest) at Cmax,ss = 0 (unitless logit)") # Figure S19b footnote 'Predictor=-1.77+ Cmaxss*7.08e-05' (PYRXG1)
    e_cmax_pyrexia <- 7.08e-05
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S19b footnote 'Predictor=-1.77+ Cmaxss*7.08e-05' (PYRXG1)
    logit_ref_alt <- -2.6
    label("Logit of the probability of grade >= 1 ALT elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S21b footnote 'Predictor=-2.6+ Cmaxss*0.000164' (ALTAEG1S)
    e_cmax_alt <- 0.000164
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S21b footnote 'Predictor=-2.6+ Cmaxss*0.000164' (ALTAEG1S)
    logit_ref_alt_g3 <- -3.27
    label("Logit of the probability of grade >= 3 ALT elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S22b footnote 'Predictor=-3.27+ Cmaxss*8.23e-06' (ALTAEG3S)
    e_cmax_alt_g3 <- 8.23e-06
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S22b footnote 'Predictor=-3.27+ Cmaxss*8.23e-06' (ALTAEG3S)
    logit_ref_ast <- -2.76
    label("Logit of the probability of grade >= 1 AST elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S23b footnote 'Predictor=-2.76+ Cmaxss*0.000199' (ASTAEG1S)
    e_cmax_ast <- 0.000199
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S23b footnote 'Predictor=-2.76+ Cmaxss*0.000199' (ASTAEG1S)
    logit_ref_ast_g3 <- -4.27
    label("Logit of the probability of grade >= 3 AST elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S24b footnote 'Predictor=-4.27+ Cmaxss*0.000358' (ASTAEG3S)
    e_cmax_ast_g3 <- 0.000358
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S24b footnote 'Predictor=-4.27+ Cmaxss*0.000358' (ASTAEG3S)
    logit_ref_alp <- -3.87
    label("Logit of the probability of grade >= 1 alkaline phosphatase elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S25b footnote 'Predictor=-3.87+ Cmaxss*0.000101' (ALPAEG1S)
    e_cmax_alp <- 0.000101
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S25b footnote 'Predictor=-3.87+ Cmaxss*0.000101' (ALPAEG1S)
    logit_ref_bili <- -4.7
    label("Logit of the probability of grade >= 1 total bilirubin elevation causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S26b footnote 'Predictor=-4.7+ Cmaxss*0.000219' (BILAEG1S)
    e_cmax_bili <- 0.000219
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S26b footnote 'Predictor=-4.7+ Cmaxss*0.000219' (BILAEG1S)
    logit_ref_edema <- -0.267
    label("Logit of the probability of grade >= 1 oedema causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S27b footnote 'Predictor=-0.267+ Cmaxss*-0.000322' (OEDAEG1S)
    e_cmax_edema <- -0.000322
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S27b footnote 'Predictor=-0.267+ Cmaxss*-0.000322' (OEDAEG1S)
    logit_ref_edema_g3 <- -2.14
    label("Logit of the probability of grade >= 3 oedema causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S28b footnote 'Predictor=-2.14+ Cmaxss*-0.000768' (OEDAEG3S)
    e_cmax_edema_g3 <- -0.000768
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S28b footnote 'Predictor=-2.14+ Cmaxss*-0.000768' (OEDAEG3S)
    logit_ref_hypoalb <- -2.85
    label("Logit of the probability of grade >= 1 hypoalbuminaemia causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S29b footnote 'Predictor=-2.85+ Cmaxss*-3.77e-05' (ALBAEG1S)
    e_cmax_hypoalb <- -3.77e-05
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S29b footnote 'Predictor=-2.85+ Cmaxss*-3.77e-05' (ALBAEG1S)
    logit_ref_diarrhea <- -2
    label("Logit of the probability of grade >= 1 diarrhoea causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S30b footnote 'Predictor=-2+ Cmaxss*-2.08e-05' (DIAAEG1S)
    e_cmax_diarrhea <- -2.08e-05
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S30b footnote 'Predictor=-2+ Cmaxss*-2.08e-05' (DIAAEG1S)
    logit_ref_vte <- -4.98
    label("Logit of the probability of grade >= 1 venous thromboembolism causally related to savolitinib at Cmax,ss = 0 (unitless logit)") # Figure S31b footnote 'Predictor=-4.98+ Cmaxss*0.000345' (VTEAEG1S)
    e_cmax_vte <- 0.000345
    label("Change in that logit per 1 ng/mL of savolitinib Cmax,ss (logit per ng/mL)") # Figure S31b footnote 'Predictor=-4.98+ Cmaxss*0.000345' (VTEAEG1S)
  })

  model({
    # Each output is an independent binomial-logit regression (Figures 5
    # and S11-S31: 'Distribution of Y axis = Binomial; Link function =
    # Logit'), linear and uncentred in the exposure unless written with log().
    prob_teae_drug_related <- expit(logit_ref_teae_rel + e_cmax_teae_rel * CMAX)
    prob_teae_drug_related_grade3 <- expit(logit_ref_teae_rel_g3 + e_cmax_teae_rel_g3 * CMAX)
    prob_teae_dosemod <- expit(logit_ref_dosemod + e_cmax_dosemod * CMAX)
    prob_dose_reduction <- expit(logit_ref_dose_red + e_cmax_dose_red * CMAX)
    prob_dose_interruption <- expit(logit_ref_dose_int + e_cmax_dose_int * CMAX)
    prob_treatment_discontinuation <- expit(logit_ref_discont + e_logcmax_discont * log(CMAX))
    prob_hepatotoxicity <- expit(logit_ref_heptox + e_cmax_heptox * CMAX)
    prob_hepatotoxicity_grade3 <- expit(logit_ref_heptox_g3 + e_cmax_heptox_g3 * CMAX)
    prob_hypersensitivity <- expit(logit_ref_hypsen + e_cmax_hypsen * CMAX)
    prob_hypersensitivity_grade3 <- expit(logit_ref_hypsen_g3 + e_cmax_hypsen_g3 * CMAX)
    prob_pyrexia <- expit(logit_ref_pyrexia + e_cmax_pyrexia * CMAX)
    prob_alt_increase <- expit(logit_ref_alt + e_cmax_alt * CMAX)
    prob_alt_increase_grade3 <- expit(logit_ref_alt_g3 + e_cmax_alt_g3 * CMAX)
    prob_ast_increase <- expit(logit_ref_ast + e_cmax_ast * CMAX)
    prob_ast_increase_grade3 <- expit(logit_ref_ast_g3 + e_cmax_ast_g3 * CMAX)
    prob_alp_increase <- expit(logit_ref_alp + e_cmax_alp * CMAX)
    prob_bilirubin_increase <- expit(logit_ref_bili + e_cmax_bili * CMAX)
    prob_edema <- expit(logit_ref_edema + e_cmax_edema * CMAX)
    prob_edema_grade3 <- expit(logit_ref_edema_g3 + e_cmax_edema_g3 * CMAX)
    prob_hypoalbuminemia <- expit(logit_ref_hypoalb + e_cmax_hypoalb * CMAX)
    prob_diarrhea <- expit(logit_ref_diarrhea + e_cmax_diarrhea * CMAX)
    prob_venous_thromboembolism <- expit(logit_ref_vte + e_cmax_vte * CMAX)
  })
}
