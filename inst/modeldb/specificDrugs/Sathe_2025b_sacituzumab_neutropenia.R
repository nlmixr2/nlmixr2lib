Sathe_2025b_sacituzumab_neutropenia <- function() {
  description <- paste0(
    "Nonproportional-odds logistic exposure-safety model for the worst ",
    "grade of NEUTROPENIA (pooling the neutropenia, neutrophil count ",
    "decreased and febrile neutropenia adverse-event terms) in adults with metastatic triple-negative breast ",
    "cancer treated with sacituzumab govitecan 8 or 10 mg/kg on days 1 ",
    "and 8 of 21-day cycles (Sathe 2025b, n = 277, pooled IMMU-132-01 and ",
    "ASCENT). Three cumulative thresholds are returned: ",
    "prob_neutropenia (grade 1 or worse), prob_neutropenia_grade3 (grade ",
    "3 or worse) and prob_neutropenia_grade4 (grade 4). Each threshold has its own ",
    "intercept and its own exposure slope on the logit, driven by CAV, ",
    "the individual average serum concentration of intact SG from the ",
    "first dose to the first report of the worst grade, in ug/mL. The ",
    "slopes are the printed grade-specific odds ratios (Table 3); the ",
    "intercepts are NOT printed: the grade >= 1 and grade >= 3 intercepts ",
    "were recovered by digitising the fitted curves of Figures 3b and 3c ",
    "with the slopes held at the printed values, and the grade 4 ",
    "intercept was approximated from the Table 3 predicted proportion at ",
    "10 mg/kg evaluated at the Figure 3 mean exposure. The grade >= 1 and ",
    "grade >= 3 slopes are statistically significant; the grade 4 slope ",
    "is not. Because the slopes differ, the published grade >= 1 and ",
    "grade >= 3 curves cross near 17 ug/mL. No covariate was retained. There is no PK layer and no ",
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
    "and were digitised from Figures 3b and 3c (grade >= 1 and >= 3) or",
    "approximated from the Table 3 predicted proportion (grade 4).",
    sep = " "
  )
  vignette <- "Sathe_2025b_sacituzumab_exposure_response"

  units <- list(
    time = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_neutropenia, prob_neutropenia_grade3 and prob_neutropenia_grade4 (cumulative probabilities of a worst neutropenia grade of at least 1, at least 3 and 4, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of INTACT sacituzumab govitecan (the antibody-drug conjugate, SG) from the first dose to the first report of the worst neutropenia grade. Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sathe 2025b Eq. 1 with t_event from Table S1: the time of the first reported worst AE grade, or 21 days after the last SG dose if no AE was reported within 21 days after the last dose. Because the window ends at an early event for most patients, these values are higher than the CR/ORR averages (Figure 3 caption: mean 15.9 ug/mL at 8 mg/kg and 22.1 ug/mL at 10 mg/kg for neutropenia, against 29.7 ug/mL at 10 mg/kg for diarrhea and 17.5 ug/mL for ORR; the windows differ by endpoint, so the columns are not interchangeable). Figures 3b and 3c show values from about 7 to 52 ug/mL. Reproduce the column from the SG output (Cc) of modellib('Sathe_2024_sacituzumab'). Enters LINEARLY and UNCENTRED.",
      source_name = "CAVG_SG"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_ASCENT = list(
      description = "Study indicator (1 = ASCENT, 0 = IMMU-132-01)",
      units = "(binary)",
      type = "binary",
      notes = "Study was screened for the safety endpoints only (Sathe 2025b Methods, Covariate selection) but no covariate was retained: 'After the inclusion of CAVG_SG, none of the evaluated covariates met the selection criteria to be included in the final models for all AEs.' The same applies to baseline body weight, age, race, region, prior treatment, number of prior lines of therapy, ECOG performance status, UGT1A1 genotype, baseline neutrophil count and granulocyte colony-stimulating factor use (the last two screened for neutropenia only). No point estimate is reported for any of them."
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
    notes = "Observed worst neutropenia grade (Sathe 2025b Table S13): none 95 (34.3%), grade >= 1 182 (65.7%), grade >= 3 149 (53.8%), grade 4 51 (18.4%). The endpoint pools the neutropenia, neutrophil count decreased and febrile neutropenia adverse-event terms (Table S2)."
  )

  ini({
    # Sathe 2025b Table 3 prints only the grade-specific odds ratios per
    # 1 ug/mL CAVG_SG; Supplementary Eq. 4 defines the nonproportional-odds
    # form logit(P(grade >= k)) = beta_0k + beta_k * CAVG_SG. The
    # intercepts beta_0k are not printed anywhere in the paper or its
    # supplement.
    #
    # Grade >= 1: digitised from the fitted solid curve of Figure 3b
    # (column-wise scan, 146 points over CAVG_SG 7.9-28.1 ug/mL). With
    # the slope held at log(1.37) the least-squares intercept on the
    # logit is -5.63 (RMS residual 0.004 on the probability scale); a
    # free fit returns an odds ratio of 1.376 against the printed 1.37.
    neutropenia_intercept <- -5.63; label("Logit of the probability of grade >= 1 neutropenia at CAVG_SG = 0 (unitless logit)") # digitised from Sathe 2025b Figure 3b with the slope held at the Table 3 value; no printed value exists
    e_cav_neutropenia <- log(1.37); label("Log-odds of grade >= 1 neutropenia per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table 3, Neutropenia Grade >= 1: odds ratio 1.37 (1.26-1.50), significant

    # Grade >= 3: digitised from the fitted solid curve of Figure 3c
    # (235 points over CAVG_SG 7.0-51.3 ug/mL). With the slope held at
    # log(1.09) the intercept is -1.68 (RMS residual 0.006); a free fit
    # returns an odds ratio of 1.093 against the printed 1.09.
    neutropenia_grade3_intercept <- -1.68; label("Logit of the probability of grade >= 3 neutropenia at CAVG_SG = 0 (unitless logit)") # digitised from Sathe 2025b Figure 3c with the slope held at the Table 3 value; no printed value exists
    e_cav_neutropenia_grade3 <- log(1.09); label("Log-odds of grade >= 3 neutropenia per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table 3, Neutropenia Grade >= 3: odds ratio 1.09 (1.05-1.14), significant

    # Grade 4: neither printed nor plotted. Approximated as
    # logit(0.186) - log(1.03) * 22.1 = -2.13, from the Table 3 predicted
    # proportion at 10 mg/kg (0.186) and the Figure 3 caption mean
    # CAVG_SG for neutropenia at 10 mg/kg (22.1 ug/mL). The same shortcut
    # applied to grade >= 3 gives -1.73 against the digitised -1.68.
    neutropenia_grade4_intercept <- -2.13; label("Logit of the probability of grade 4 neutropenia at CAVG_SG = 0 (unitless logit)") # approximated from Sathe 2025b Table 3 (0.186 at 10 mg/kg) and the Figure 3 caption mean CAVG_SG 22.1 ug/mL; no printed value exists
    e_cav_neutropenia_grade4 <- log(1.03); label("Log-odds of grade 4 neutropenia per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table 3, Neutropenia Grade 4: odds ratio 1.03 (0.988-1.07), not significant

    addSd_prob_neutropenia <- fixed(0.001); label("Placeholder additive residual SD on prob_neutropenia; the source likelihood is multinomial") # not from source; see vignette Assumptions and deviations
    addSd_prob_neutropenia_grade3 <- fixed(0.001); label("Placeholder additive residual SD on prob_neutropenia_grade3; the source likelihood is multinomial") # not from source; see vignette Assumptions and deviations
    addSd_prob_neutropenia_grade4 <- fixed(0.001); label("Placeholder additive residual SD on prob_neutropenia_grade4; the source likelihood is multinomial") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Cumulative-threshold probabilities of the worst neutropenia grade.
    # With grade-specific slopes the thresholds are not constrained to
    # stay ordered, and the published curves are not: below about
    # 17 ug/mL the grade >= 3 curve (Figure 3c) lies above the grade >= 1
    # curve (Figure 3b), and below about 12 ug/mL the grade 4 curve does
    # too. Per-grade probabilities must not be obtained by differencing
    # in that range.
    logit_neutropenia <- neutropenia_intercept + e_cav_neutropenia * CAV
    logit_neutropenia_grade3 <- neutropenia_grade3_intercept + e_cav_neutropenia_grade3 * CAV
    logit_neutropenia_grade4 <- neutropenia_grade4_intercept + e_cav_neutropenia_grade4 * CAV

    prob_neutropenia <- expit(logit_neutropenia)
    prob_neutropenia_grade3 <- expit(logit_neutropenia_grade3)
    prob_neutropenia_grade4 <- expit(logit_neutropenia_grade4)

    prob_neutropenia ~ add(addSd_prob_neutropenia)
    prob_neutropenia_grade3 ~ add(addSd_prob_neutropenia_grade3)
    prob_neutropenia_grade4 ~ add(addSd_prob_neutropenia_grade4)
  })
}
