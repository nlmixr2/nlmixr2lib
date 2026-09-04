Ravva_2010_varenicline_car_w9_12 <- function() {
  description <- paste0(
    "Population exposure-response logistic-regression model for the ",
    "continuous abstinence rate (CAR) at weeks 9-12 in 1,892 adult ",
    "cigarette smokers treated with varenicline for smoking cessation ",
    "(Ravva 2010, final full model). The logit of the quit probability ",
    "is a baseline intercept multiplied by covariate factors for ",
    "nicotine dependence (Fagerstrom time to the first cigarette after ",
    "waking, 4-level score), age, sex and race, PLUS an additive linear ",
    "term in the individual steady-state varenicline exposure ",
    "AUC(0-24)ss (Ravva 2010 Equation 2). Exposure is a data covariate ",
    "computed as daily dose / (CL/F) from the companion population PK ",
    "model (see Ravva_2009_varenicline.R) and is set to 0 for placebo ",
    "subjects. Fitted in NONMEM V with the Laplacian method on the ",
    "dichotomous quit / no-quit outcome; the final model is naive ",
    "pooled, so no between-subject random effect and no residual-error ",
    "term are reported. Companion tolerability model in ",
    "Ravva_2010_varenicline_nausea.R; companion preliminary weeks 4-7 ",
    "models in Ravva_2010_varenicline_car_w4_7_study1.R and ",
    "Ravva_2010_varenicline_car_w4_7_study2.R."
  )
  reference <- paste(
    "Ravva P, Gastonguay MR, French JL, Tensfeldt TG, Faessel HM.",
    "Quantitative assessment of exposure-response relationships for the",
    "efficacy and tolerability of varenicline for smoking cessation.",
    "Clin Pharmacol Ther. 2010;87(3):336-344. doi:10.1038/clpt.2009.282.",
    "Companion population PK model supplying AUC(0-24)ss:",
    "Ravva P, Gastonguay MR, Tensfeldt TG, Faessel HM.",
    "Population pharmacokinetic analysis of varenicline in adult smokers.",
    "Br J Clin Pharmacol. 2009;68(5):669-681.",
    "doi:10.1111/j.1365-2125.2009.03520.x."
  )
  vignette <- "Ravva_2010_varenicline_exposure_response"
  units <- list(
    time          = "week",
    dosing        = "n/a (exposure-response model; varenicline exposure enters as the AUC_VAREN covariate rather than as a dosing event)",
    concentration = "p_car (probability of continuous abstinence at weeks 9-12, 0-1; also logit_car)"
  )

  covariateData <- list(
    AUC_VAREN = list(
      description        = "Individual varenicline steady-state daily exposure, AUC(0-24)ss. Ravva 2010 Methods (Pharmacokinetics): 'Individual 24-h daily exposure, measured as AUC(0-24)ss, was calculated as dose divided by CL/F_i (apparent clearance); the individual empirical Bayes estimate of apparent clearance was predicted from the final population PK model and parameters obtained from a pooled analysis in adult smokers.'",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the logit ADDITIVELY (not as a multiplicative factor on the intercept like the demographic covariates) -- see Ravva 2010 Equation 2. Set to 0 for placebo subjects (Ravva 2010 Figure 3a caption: 'exposure was set to zero for the placebo group'). Observed cohort range spans the three exposure bins reported in Results: 50-142, >142-184 and >184-408 ng*h/mL. The reference-population quit probability of 0.562 at 1 mg b.i.d. (Ravva 2010 Figure 3b) implies a typical AUC(0-24)ss of about 186 ng*h/mL at 1 mg b.i.d. and about 93 ng*h/mL at 0.5 mg b.i.d. Downstream users should compute the per-subject value from the companion Ravva 2009 varenicline population PK model as total daily dose / (CL/F).",
      source_name        = "AUC(0-24)ss"
    ),
    SMOKE_TTFC_SCORE = list(
      description        = "Fagerstrom Test for Nicotine Dependence item 1 ('How soon after you wake up do you smoke your first cigarette?') scored 0-3 as published by Ravva 2010 Methods: >60 min (0); 31-60 min (1); 6-30 min (2); within 5 min (3). Higher score = greater nicotine dependence.",
      units              = "(ordinal score 0-3)",
      type               = "categorical",
      reference_category = "0 = first cigarette more than 60 min after waking (the model reference level; all three derived indicators are 0)",
      notes              = "Decomposed inside model() into three 0/1 indicators, one per non-reference level, each carrying its own estimated exponent -- the source model does NOT impose a linear per-level effect. Note the deliberate distinction between the MODEL reference category (score 0, >60 min) and the REPORTING reference population (score 2, 6-30 min) used by Ravva 2010 to quote probabilities: Methods defines the reference population as 'Caucasian, 45-year-old, male smokers who smoke their first cigarette within 6-30 min after waking in the morning'. Cohort distribution for this endpoint (Ravva 2010 Table 1, CAR weeks 9-12 column): <=5 min 636; 6-30 min 832; 31-60 min 284; >60 min 140 (total 1,892). Preferred over the cigarettes-per-day item because the two are collinear (Results: 'Simultaneous inclusion of the correlated predictors FSQ1, FSQ4, and CO in the model was avoided by selecting the FSQ1 variable as the most relevant predictor of nicotine dependence').",
      source_name        = "FSQ1"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on the baseline logit with reference 45 years, (AGE/45)^e_age_base_logit (Ravva 2010 Equation 2). Cohort mean 43 years, range 18-75 years (Ravva 2010 Table 1, CAR weeks 9-12 column). The negative exponent means the baseline quit probability RISES with age: Results report a baseline quit probability increasing from 0.35 at age 18 to 0.64 at age 75.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male (the Ravva 2010 reference population is male)",
      notes              = "Ravva 2010 Equation 2 writes this term as theta7^(1-Sex) with the effect row labelled 'Female' in Table 3 and the reference population defined as male, so the paper's 'Sex' column is 1 = male and (1 - Sex) is identically the canonical SEXF (1 = female). Encoded directly as e_sexf_base_logit^SEXF with no value transformation. Cohort 889 female (47%) of 1,892 (Ravva 2010 Table 1). The estimate 1.02 is essentially null; Results conclude 'Gender had no influence on the efficacy of varenicline.'",
      source_name        = "Sex"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the Ravva 2010 reference population is Caucasian)",
      notes              = "Power-of-indicator form on the baseline logit: e_race_black_base_logit^RACE_BLACK (Ravva 2010 Equation 2). Cohort 210 (11%) of 1,892 (Ravva 2010 Table 1). Same canonical column as the companion Ravva_2009_varenicline.R population PK model. Because the intercept is negative, a factor above 1 LOWERS the quit probability; Results note 'a trend, albeit less precisely defined, toward a decreased baseline quit probability in blacks as compared with whites'.",
      source_name        = "Race (Black)"
    ),
    RACE_OTHER = list(
      description        = "Composite 'Other' race indicator pooling Asian, Hispanic and 'Other' races (the Ravva 2010 grouping)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the Ravva 2010 reference population is Caucasian)",
      notes              = "Power-of-indicator form on the baseline logit: e_race_other_base_logit^RACE_OTHER (Ravva 2010 Equation 2). Table 1 footnote a: 'Other includes asian, hispanic, and Other races.' Cohort 141 (8%) of 1,892 (Ravva 2010 Table 1). Same canonical column and same composite definition as the companion Ravva_2009_varenicline.R population PK model, which is why RACE_OTHER is used here in preference to the semantically narrower RACE_NONBLACK_NONWHITE.",
      source_name        = "Race (other)"
    )
  )

  # Screened by the source analysis but deliberately NOT retained in the final
  # model, so they are documentation only and are never referenced in model().
  covariatesDataExcluded <- list(
    SMOKE_CPD_SCORE = list(
      description        = "Fagerstrom Test for Nicotine Dependence item 4 ('How many cigarettes per day do you smoke?') scored 0-3 as published by Ravva 2010 Methods: 10 or less (0); 11-20 (1); 21-30 (2); >=31 (3).",
      units              = "(ordinal score 0-3)",
      type               = "categorical",
      reference_category = "0 = 10 or fewer cigarettes per day",
      notes              = "Screened as a candidate nicotine-dependence predictor but dropped for collinearity with SMOKE_TTFC_SCORE (Ravva 2010 Results: 'Simultaneous inclusion of the correlated predictors FSQ1, FSQ4, and CO in the model was avoided by selecting the FSQ1 variable as the most relevant predictor of nicotine dependence'). Discussion adds that time to the first cigarette 'was found to be more predictive of smoking cessation than the number of cigarettes smoked per day'. Cohort distribution (Ravva 2010 Table 1, CAR weeks 9-12 column): <=10 109; 11-20 1,061; 21-30 511; >30 211.",
      source_name        = "FSQ4"
    ),
    CO_EXHALED_PPM = list(
      description        = "Exhaled carbon monoxide concentration at baseline, a biochemical marker of recent smoking intensity",
      units              = "ppm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate baseline-smoking-status predictor but dropped for collinearity with SMOKE_TTFC_SCORE (same Results sentence as SMOKE_CPD_SCORE). Retained in the source analysis only for two non-modelling purposes: confirming self-reported abstinence (Methods: 'CAR was confirmed by an exhaled CO measurement of 10 p.p.m. or less') and imputing eight missing baseline values by multivariate regression on FSQ1 and gender. Cohort mean 22 ppm, range 1-81 ppm (Ravva 2010 Table 1, CAR weeks 9-12 column).",
      source_name        = "CO"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1892L,
    n_studies      = 3L,
    age_range      = "18-75 years",
    age_mean       = "43 years",
    weight_range   = NA_character_,
    sex_female_pct = 47.0,
    race_ethnicity = c(White = 81.0, Black = 11.0, Other = 8.0),
    disease_state  = "Adult cigarette smokers motivated to stop smoking. More than 75% reported smoking the first cigarette of the day within 30 min of waking; the majority smoked at least 11 cigarettes/day and more than 40% of those smoked at least 20 cigarettes/day.",
    dose_range     = "Varenicline 0.5 mg b.i.d. and 1 mg b.i.d. oral, plus placebo and (in studies 3 and 4) a bupropion SR 150 mg b.i.d. comparator arm, over a 12-week treatment period after an initial titration week.",
    regions        = "Multicenter, multinational.",
    smoking_marker = "Baseline exhaled carbon monoxide mean 22 ppm, range 1-81 ppm.",
    notes          = "Pooled from studies 2, 3 and 4 of the five randomized, double-blind, placebo-controlled trials summarized in Ravva 2010 Supplementary Table S1 (study 2 n = 609, study 3 n = 642, study 4 n = 641). Baseline demographics from Ravva 2010 Table 1, 'CAR at weeks 9-12' column. A total of 416 subjects (22% of this database) did not complete the treatment period and were assigned as quit failures; subjects with unknown smoking status or exhaled CO > 10 ppm were assumed to be smoking for the remainder of the study."
  )

  ini({
    # =====================================================================
    # All values are the FINAL full-model estimates from Ravva 2010 Table 3,
    # column 'CAR at weeks 9-12'. Percentages in the trailing comments are
    # the published %RSE (= SE / |estimate| * 100, per the Table 3 footnote).
    # Parameter numbering follows Ravva 2010 Equation 2.
    #
    # Model shape (Equation 2): demographic covariates enter as
    # MULTIPLICATIVE factors on the intercept, while the drug-exposure term
    # is ADDITIVE on the logit. That asymmetry is load-bearing -- it is what
    # reproduces the paper's reported quit probabilities exactly.
    # =====================================================================

    # ----- Baseline logit (intercept) at the reference covariates -----
    # Reference: White, male, 45 years old, FSQ1 score 0 (first cigarette
    # more than 60 min after waking), AUC(0-24)ss = 0.
    base_logit                    <- -0.657   ; label("Baseline logit of the weeks 9-12 quit probability at the reference covariates (unitless logit)")  # Ravva 2010 Table 3 'Intercept (theta1)' = -0.657 (25.7% RSE)

    # ----- Drug effect: linear in AUC(0-24)ss, ADDITIVE on the logit -----
    e_auc_varen_base_logit        <- 0.00813  ; label("Slope of the logit on varenicline AUC(0-24)ss (logit units per ng*h/mL)")  # Ravva 2010 Table 3 'Drug effect, AUC (theta2)' = 0.00813 (7.2% RSE)

    # ----- Nicotine dependence (FSQ1): one exponent per non-reference level -----
    # Multiplicative factors on the intercept. Because base_logit is negative,
    # a factor > 1 makes the logit more negative and LOWERS the quit
    # probability, which is the reported direction: quit probability falls
    # from 0.70 (>60 min) to 0.45 (<=5 min) at 1 mg b.i.d.
    e_smoke_ttfc_31_60_base_logit <- 1.54     ; label("Multiplicative factor on the baseline logit for first cigarette 31-60 min after waking, vs >60 min (unitless)")  # Ravva 2010 Table 3 'Effect of FSQ1: 31-60min' = 1.54 (26.9% RSE); theta3
    e_smoke_ttfc_6_30_base_logit  <- 1.92     ; label("Multiplicative factor on the baseline logit for first cigarette 6-30 min after waking, vs >60 min (unitless)")   # Ravva 2010 Table 3 'Effect of FSQ1: 6-30min' = 1.92 (24.6% RSE); theta4
    e_smoke_ttfc_le5_base_logit   <- 2.59     ; label("Multiplicative factor on the baseline logit for first cigarette within 5 min of waking, vs >60 min (unitless)")  # Ravva 2010 Table 3 'Effect of FSQ1: <=5min' = 2.59 (24.6% RSE); theta5

    # ----- Age: power function of AGE/45 on the intercept -----
    e_age_base_logit              <- -0.563   ; label("Power exponent of (AGE/45) on the baseline logit (unitless)")  # Ravva 2010 Table 3 'Effect of age, 18-75 years' = -0.563 (24.7% RSE); theta6

    # ----- Sex and race: power-of-indicator factors on the intercept -----
    e_sexf_base_logit             <- 1.02     ; label("Multiplicative factor on the baseline logit for female sex, vs male (unitless)")  # Ravva 2010 Table 3 'Effect of gender: Female' = 1.02 (7.1% RSE); theta7
    e_race_black_base_logit       <- 1.27     ; label("Multiplicative factor on the baseline logit for Black race, vs White (unitless)")  # Ravva 2010 Table 3 'Effect of race: Black' = 1.27 (10.8% RSE); theta8
    e_race_other_base_logit       <- 1.09     ; label("Multiplicative factor on the baseline logit for Other race, vs White (unitless)")  # Ravva 2010 Table 3 'Effect of race: Other' = 1.09 (13.9% RSE); theta9

    # ----- Residual error (placeholder; NOT from the source) -----
    # The source fit is a NONMEM Laplacian Bernoulli likelihood on the
    # dichotomous quit outcome and reports no $SIGMA. rxode2 / nlmixr2 expose
    # llikBinom / rxbinom but no Bernoulli ENDPOINT that parses in model(),
    # so the deterministic typical-value probability p_car is the declared
    # output and a tiny placeholder additive residual satisfies the
    # observation declaration. This does not change the typical-value
    # prediction. Mirrors Oniki_2018_nafld_risk.R and
    # Lindauer_2017_lacosamide_dropout.R. See vignette Assumptions and
    # deviations.
    addSd_p_car                   <- fixed(0.001) ; label("Placeholder additive residual SD on p_car; the source likelihood is Bernoulli, so no residual-error term is published")  # not paper-derived; see Assumptions and deviations
  })

  model({
    # ----- Nicotine-dependence multiplier (Ravva 2010 Equation 2) -----
    # FSQ1 is carried as the published 0-3 ordinal score and decomposed here
    # into the three non-reference indicators the source model uses. Score 0
    # (>60 min) makes every indicator 0, so the product collapses to 1.
    ttfc_mult <- e_smoke_ttfc_31_60_base_logit^(SMOKE_TTFC_SCORE == 1) *
                 e_smoke_ttfc_6_30_base_logit^(SMOKE_TTFC_SCORE == 2) *
                 e_smoke_ttfc_le5_base_logit^(SMOKE_TTFC_SCORE == 3)

    # ----- Individual baseline logit: covariates multiply the intercept -----
    base_logit_i <- base_logit *
                    ttfc_mult *
                    (AGE / 45)^e_age_base_logit *
                    e_sexf_base_logit^SEXF *
                    e_race_black_base_logit^RACE_BLACK *
                    e_race_other_base_logit^RACE_OTHER

    # ----- Drug effect adds linearly on the logit (Ravva 2010 Equation 2) -----
    logit_car <- base_logit_i + e_auc_varen_base_logit * AUC_VAREN

    # ----- Probability of continuous abstinence at weeks 9-12 -----
    p_car <- expit(logit_car)

    # Deterministic typical-value declaration with a placeholder residual
    # (see the ini() comment). Downstream callers can use p_car directly as
    # a probability, or draw Bernoulli outcomes externally with
    # rbinom(n, 1, p_car) on the rxSolve output.
    p_car ~ add(addSd_p_car)
  })
}
