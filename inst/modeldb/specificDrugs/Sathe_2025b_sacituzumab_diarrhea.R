Sathe_2025b_sacituzumab_diarrhea <- function() {
  description <- paste0(
    "Nonproportional-odds logistic exposure-safety model for the worst ",
    "grade of DIARRHEA in adults with metastatic triple-negative breast ",
    "cancer treated with sacituzumab govitecan 8 or 10 mg/kg on days 1 ",
    "and 8 of 21-day cycles (Sathe 2025b, n = 277, pooled IMMU-132-01 and ",
    "ASCENT). Two cumulative thresholds are returned: prob_diarrhea ",
    "(grade 1 or worse) and prob_diarrhea_grade3 (grade 3 or worse; no ",
    "grade 4 diarrhea was observed). Each threshold has its own ",
    "intercept and its own exposure slope on the logit, driven by CAV, ",
    "the individual average serum concentration of intact SG from the ",
    "first dose to the first report of the worst grade, in ug/mL. The ",
    "slopes are the printed grade-specific odds ratios (Table 3); the ",
    "intercepts are NOT printed: the grade >= 1 intercept was recovered ",
    "by digitising the fitted curve of Figure 3a with the slope held at ",
    "the printed value, and the grade >= 3 intercept was approximated ",
    "from the Table 3 predicted proportion at 10 mg/kg evaluated at the ",
    "Figure 3 mean exposure. Only the grade >= 1 slope is statistically ",
    "significant. No covariate was retained. There is no PK layer and no ",
    "ODE; exposure derives from the companion population PK model ",
    "Sathe_2024_sacituzumab. One of six exposure-response models in the ",
    "Sathe_2025b_sacituzumab_* family."
  )
  reference <- paste(
    "Sathe AG, Diderichsen PM, Fauchet F, Phan SC, Girish S, Othman AA.",
    "Exposure-Response Analyses of Sacituzumab Govitecan Efficacy and",
    "Safety in Patients With Metastatic Triple-Negative Breast Cancer.",
    "Clin Pharmacol Ther. 2025;117(2):570-578. doi:10.1002/cpt.3495.",
    "PMCID PMC11739744. Slopes from Table 3; intercepts are not printed",
    "and were digitised from Figure 3a (grade >= 1) or approximated from",
    "the Table 3 predicted proportion (grade >= 3).",
    sep = " "
  )
  vignette <- "Sathe_2025b_sacituzumab_exposure_response"

  units <- list(
    time = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_diarrhea and prob_diarrhea_grade3 (cumulative probabilities of a worst diarrhea grade of at least 1 and at least 3, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of INTACT sacituzumab govitecan (the antibody-drug conjugate, SG) from the first dose to the first report of the worst diarrhea grade. Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sathe 2025b Eq. 1 with t_event from Table S1: the time of the first reported worst AE grade, or 21 days after the last SG dose if no AE was reported within 21 days after the last dose. Because the window ends at an early event for most patients, these values are higher than the CR/ORR averages (Figure 3 caption: mean 19.6 ug/mL at 8 mg/kg and 29.7 ug/mL at 10 mg/kg for diarrhea, against 17.5 ug/mL at 10 mg/kg for ORR) and the distribution is strongly right-skewed (Figure 3a extends to about 160 ug/mL). Reproduce the column from the SG output (Cc) of modellib('Sathe_2024_sacituzumab'). Enters LINEARLY and UNCENTRED.",
      source_name = "CAVG_SG"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_ASCENT = list(
      description = "Study indicator (1 = ASCENT, 0 = IMMU-132-01)",
      units = "(binary)",
      type = "binary",
      notes = "Study was screened for the safety endpoints only (Sathe 2025b Methods, Covariate selection) but no covariate was retained: 'After the inclusion of CAVG_SG, none of the evaluated covariates met the selection criteria to be included in the final models for all AEs.' The same applies to baseline body weight, age, race, region, prior treatment, number of prior lines of therapy, ECOG performance status and UGT1A1 genotype. No point estimate is reported for any of them."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 277L,
    n_studies = 2L,
    age_range = "27-82 years (median 54)",
    weight_range = "37.2-132 kg (median 66.7)",
    sex_female_pct = 99.3,
    race_ethnicity = c(White = 80.5, Black = 11.2, Asian = 4.0, Other = 4.3),
    disease_state = "metastatic triple-negative breast cancer, previously treated",
    dose_range = "sacituzumab govitecan 8, 10 or 12 mg/kg starting dose IV on days 1 and 8 of 21-day cycles (10 mg/kg in about 92%)",
    regions = "North America 67.9%, Europe 32.1%",
    notes = "Observed worst diarrhea grade (Sathe 2025b Table S13): none 96 (34.7%), grade >= 1 181 (65.3%), grade >= 3 30 (10.8%), grade 4 none. Diarrhea grades follow CTCAE as reported in the source trials."
  )

  ini({
    # Sathe 2025b Table 3 prints only the grade-specific odds ratios per
    # 1 ug/mL CAVG_SG; Supplementary Eq. 4 defines the nonproportional-odds
    # form logit(P(grade >= k)) = beta_0k + beta_k * CAVG_SG. The
    # intercepts beta_0k are not printed anywhere in the paper or its
    # supplement.
    #
    # Grade >= 1: digitised from the fitted solid curve of Figure 3a.
    # The curve is steep (about 3 pixels per ug/mL), so it was scanned
    # both by row (41 points, probability 0.07-0.91) and by column. With
    # the slope held at log(1.45) the two scans give intercepts of -6.50
    # and -6.54; -6.52 is their mean, uncertain by about +/-0.2 (one
    # pixel). A free fit to the row scan returns an odds ratio of 1.452,
    # matching the printed 1.45, which confirms the digitisation.
    diarrhea_intercept <- -6.52; label("Logit of the probability of grade >= 1 diarrhea at CAVG_SG = 0 (unitless logit)") # digitised from Sathe 2025b Figure 3a with the slope held at the Table 3 value; no printed value exists
    e_cav_diarrhea <- log(1.45); label("Log-odds of grade >= 1 diarrhea per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table 3, Diarrhea Grade >= 1: odds ratio 1.45 (1.30-1.63), significant

    # Grade >= 3: neither printed nor plotted. Approximated as
    # logit(0.109) - log(1.01) * 29.7 = -2.40, from the Table 3 predicted
    # proportion at 10 mg/kg (0.109) and the Figure 3 caption mean
    # CAVG_SG for diarrhea at 10 mg/kg (29.7 ug/mL). This treats the
    # dataset-averaged prediction as the prediction at the mean exposure.
    diarrhea_grade3_intercept <- -2.40; label("Logit of the probability of grade >= 3 diarrhea at CAVG_SG = 0 (unitless logit)") # approximated from Sathe 2025b Table 3 (0.109 at 10 mg/kg) and the Figure 3 caption mean CAVG_SG 29.7 ug/mL; no printed value exists
    e_cav_diarrhea_grade3 <- log(1.01); label("Log-odds of grade >= 3 diarrhea per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table 3, Diarrhea Grade >= 3: odds ratio 1.01 (0.991-1.02), not significant

    addSd_prob_diarrhea <- fixed(0.001); label("Placeholder additive residual SD on prob_diarrhea; the source likelihood is multinomial") # not from source; see vignette Assumptions and deviations
    addSd_prob_diarrhea_grade3 <- fixed(0.001); label("Placeholder additive residual SD on prob_diarrhea_grade3; the source likelihood is multinomial") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Cumulative-threshold probabilities of the worst diarrhea grade.
    # With grade-specific slopes the two thresholds are not constrained
    # to stay ordered: below about 11 ug/mL the grade >= 3 curve lies
    # above the grade >= 1 curve, so per-grade probabilities must not be
    # obtained by differencing there.
    logit_diarrhea <- diarrhea_intercept + e_cav_diarrhea * CAV
    logit_diarrhea_grade3 <- diarrhea_grade3_intercept + e_cav_diarrhea_grade3 * CAV

    prob_diarrhea <- expit(logit_diarrhea)
    prob_diarrhea_grade3 <- expit(logit_diarrhea_grade3)

    prob_diarrhea ~ add(addSd_prob_diarrhea)
    prob_diarrhea_grade3 ~ add(addSd_prob_diarrhea_grade3)
  })
}
