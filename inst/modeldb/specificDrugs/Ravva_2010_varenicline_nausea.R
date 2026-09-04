Ravva_2010_varenicline_nausea <- function() {
  description <- paste0(
    "Population exposure-response logistic-regression model for nausea ",
    "incidence over a 12-week varenicline treatment period in 2,238 ",
    "adult cigarette smokers (Ravva 2010, final full model). The logit ",
    "of the nausea probability is a baseline intercept multiplied by ",
    "covariate factors for nicotine dependence (Fagerstrom time to the ",
    "first cigarette after waking, 4-level score), age, sex and race, ",
    "PLUS an additive linear term in the individual steady-state ",
    "varenicline exposure AUC(0-24)ss, PLUS the exponentially decaying ",
    "week-of-treatment term of Ravva 2010 Equation 2. Nausea is about ",
    "twofold more likely in women than in men regardless of treatment. ",
    "Exposure is a data covariate computed as daily dose / (CL/F) from ",
    "the companion population PK model (see Ravva_2009_varenicline.R) ",
    "and is set to 0 for placebo subjects. The published estimates are ",
    "the whole-treatment-period incidence model (one observation per ",
    "subject); the two parameters of the weekly time-course term and ",
    "the between-subject random effect are never reported and are ",
    "carried as fixed(0), leaving the printed equation structure ",
    "visible and user-perturbable. Companion efficacy model in ",
    "Ravva_2010_varenicline_car_w9_12.R."
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
    concentration = "p_nausea (probability of a nausea event, 0-1; also logit_nausea)"
  )

  covariateData <- list(
    AUC_VAREN = list(
      description        = "Individual varenicline steady-state daily exposure, AUC(0-24)ss. Ravva 2010 Methods (Pharmacokinetics): 'Individual 24-h daily exposure, measured as AUC(0-24)ss, was calculated as dose divided by CL/F_i (apparent clearance); the individual empirical Bayes estimate of apparent clearance was predicted from the final population PK model and parameters obtained from a pooled analysis in adult smokers.'",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the logit ADDITIVELY (not as a multiplicative factor on the intercept like the demographic covariates) -- see Ravva 2010 Equation 2. Set to 0 for placebo subjects (Ravva 2010 Figure 4a caption: 'exposure was set to zero for the placebo group'). The reference-population quit probability of 0.562 at 1 mg b.i.d. in the companion efficacy model (Ravva 2010 Figure 3b) implies a typical AUC(0-24)ss of about 186 ng*h/mL at 1 mg b.i.d. and about 93 ng*h/mL at 0.5 mg b.i.d.; those two values reproduce all four of the dose-specific nausea probabilities reported in Results. Downstream users should compute the per-subject value from the companion Ravva 2009 varenicline population PK model as total daily dose / (CL/F).",
      source_name        = "AUC(0-24)ss"
    ),
    SMOKE_TTFC_SCORE = list(
      description        = "Fagerstrom Test for Nicotine Dependence item 1 ('How soon after you wake up do you smoke your first cigarette?') scored 0-3 as published by Ravva 2010 Methods: >60 min (0); 31-60 min (1); 6-30 min (2); within 5 min (3). Higher score = greater nicotine dependence.",
      units              = "(ordinal score 0-3)",
      type               = "categorical",
      reference_category = "0 = first cigarette more than 60 min after waking (the model reference level; all three derived indicators are 0)",
      notes              = "Decomposed inside model() into three 0/1 indicators, one per non-reference level, each carrying its own estimated exponent -- the source model does NOT impose a linear per-level effect. Note the deliberate distinction between the MODEL reference category (score 0, >60 min) and the REPORTING reference population (score 2, 6-30 min) used by Ravva 2010 to quote probabilities. Cohort distribution for this endpoint (Ravva 2010 Table 1, 'Nausea incidence' column): <=5 min 759; 6-30 min 989; 31-60 min 324; >60 min 166 (total 2,238). Unlike the efficacy model, the three nausea FSQ1 factors (0.894, 0.867, 0.961) are all below 1 and close together, so nicotine dependence shifts the nausea probability far less than it shifts the quit probability.",
      source_name        = "FSQ1"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on the baseline logit with reference 45 years, (AGE/45)^e_age_base_logit (Ravva 2010 Equation 2). Cohort mean 44 years, range 18-75 years (Ravva 2010 Table 1, 'Nausea incidence' column). The positive exponent 0.374 is the opposite sign to the efficacy model's -0.563: because the intercept is negative, a positive exponent means the nausea probability DECREASES with increasing age.",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male (the Ravva 2010 reference population is male)",
      notes              = "Ravva 2010 Equation 2 writes this term as theta7^(1-Sex) with the effect row labelled 'Female' in Table 3 and the reference population defined as male, so the paper's 'Sex' column is 1 = male and (1 - Sex) is identically the canonical SEXF (1 = female). Encoded directly as e_sexf_base_logit^SEXF with no value transformation. Cohort 1,063 female (47%) of 2,238 (Ravva 2010 Table 1). This is the dominant tolerability covariate: the factor 0.704 on a negative intercept roughly doubles the odds of nausea in women, consistently in both arms (placebo 13.5% female vs 6.90% male, ratio 1.96; active 37.6% vs 19.0%, ratio 1.98), which Results interpret as 'women in general are more inclined that men to experience nausea as opposed to varenicline having a different effect in females'.",
      source_name        = "Sex"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the Ravva 2010 reference population is Caucasian)",
      notes              = "Power-of-indicator form on the baseline logit: e_race_black_base_logit^RACE_BLACK (Ravva 2010 Equation 2). Cohort 234 (10%) of 2,238 (Ravva 2010 Table 1). Same canonical column as the companion Ravva_2009_varenicline.R population PK model.",
      source_name        = "Race (Black)"
    ),
    RACE_OTHER = list(
      description        = "Composite 'Other' race indicator pooling Asian, Hispanic and 'Other' races (the Ravva 2010 grouping)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White (the Ravva 2010 reference population is Caucasian)",
      notes              = "Power-of-indicator form on the baseline logit: e_race_other_base_logit^RACE_OTHER (Ravva 2010 Equation 2). Table 1 footnote a: 'Other includes asian, hispanic, and Other races.' Cohort 158 (7%) of 2,238 (Ravva 2010 Table 1). Same canonical column and same composite definition as the companion Ravva_2009_varenicline.R population PK model.",
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
      notes              = "Screened as a candidate nicotine-dependence predictor but dropped for collinearity with SMOKE_TTFC_SCORE (Ravva 2010 Results: 'Simultaneous inclusion of the correlated predictors FSQ1, FSQ4, and CO in the model was avoided by selecting the FSQ1 variable as the most relevant predictor of nicotine dependence'). Cohort distribution (Ravva 2010 Table 1, 'Nausea incidence' column): <=10 122; 11-20 1,243; 21-30 606; >30 267.",
      source_name        = "FSQ4"
    ),
    CO_EXHALED_PPM = list(
      description        = "Exhaled carbon monoxide concentration at baseline, a biochemical marker of recent smoking intensity",
      units              = "ppm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a candidate baseline-smoking-status predictor but dropped for collinearity with SMOKE_TTFC_SCORE (same Results sentence as SMOKE_CPD_SCORE). Cohort mean 22 ppm, range 1-81 ppm (Ravva 2010 Table 1, 'Nausea incidence' column).",
      source_name        = "CO"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2238L,
    n_studies      = 4L,
    age_range      = "18-75 years",
    age_mean       = "44 years",
    weight_range   = NA_character_,
    sex_female_pct = 47.5,
    race_ethnicity = c(White = 83.0, Black = 10.0, Other = 7.0),
    disease_state  = "Adult cigarette smokers motivated to stop smoking. More than 75% reported smoking the first cigarette of the day within 30 min of waking; the majority smoked at least 11 cigarettes/day and more than 40% of those smoked at least 20 cigarettes/day.",
    dose_range     = "Varenicline 0.5 mg b.i.d. and 1 mg b.i.d. oral, plus placebo and (in studies 3 and 4) a bupropion SR 150 mg b.i.d. comparator arm. Treatment lasted 12 weeks in studies 2, 3 and 4 and 52 weeks in the long-term safety study 5.",
    regions        = "Multicenter, multinational.",
    smoking_marker = "Baseline exhaled carbon monoxide mean 22 ppm, range 1-81 ppm.",
    notes          = "Pooled from studies 2, 3, 4 and 5 of the five randomized, double-blind, placebo-controlled trials summarized in Ravva 2010 Supplementary Table S1 (study 2 n = 609, study 3 n = 642, study 4 n = 641, study 5 n = 346). Baseline demographics from Ravva 2010 Table 1, 'Nausea incidence' column. Nausea was derived from treatment-emergent adverse events, defined as events beginning or intensifying on or after the first day of study medication and within 7 days of the last dose. For subjects who discontinued early, nausea incidence was reported as observed up to the point of dropout."
  )

  ini({
    # =====================================================================
    # All estimated values are the FINAL full-model estimates from Ravva
    # 2010 Table 3, column 'Nausea incidence'. Percentages in the trailing
    # comments are the published %RSE (= SE / |estimate| * 100, per the
    # Table 3 footnote). Parameter numbering follows Ravva 2010 Equation 2.
    #
    # Model shape (Equation 2): demographic covariates enter as
    # MULTIPLICATIVE factors on the intercept, while the drug-exposure term
    # and the week-of-treatment term are ADDITIVE on the logit.
    # =====================================================================

    # ----- Baseline logit (intercept) at the reference covariates -----
    # Reference: White, male, 45 years old, FSQ1 score 0 (first cigarette
    # more than 60 min after waking), AUC(0-24)ss = 0.
    base_logit                    <- -2.93    ; label("Baseline logit of the nausea probability at the reference covariates (unitless logit)")  # Ravva 2010 Table 3 'Intercept (theta1)' = -2.93 (10.5% RSE)

    # ----- Drug effect: linear in AUC(0-24)ss, ADDITIVE on the logit -----
    e_auc_varen_base_logit        <- 0.00738  ; label("Slope of the logit on varenicline AUC(0-24)ss (logit units per ng*h/mL)")  # Ravva 2010 Table 3 'Drug effect, AUC (theta2)' = 0.00738 (8.5% RSE)

    # ----- Nicotine dependence (FSQ1): one exponent per non-reference level -----
    e_smoke_ttfc_31_60_base_logit <- 0.894    ; label("Multiplicative factor on the baseline logit for first cigarette 31-60 min after waking, vs >60 min (unitless)")  # Ravva 2010 Table 3 'Effect of FSQ1: 31-60min' = 0.894 (11.4% RSE); theta3
    e_smoke_ttfc_6_30_base_logit  <- 0.867    ; label("Multiplicative factor on the baseline logit for first cigarette 6-30 min after waking, vs >60 min (unitless)")   # Ravva 2010 Table 3 'Effect of FSQ1: 6-30min' = 0.867 (10.2% RSE); theta4
    e_smoke_ttfc_le5_base_logit   <- 0.961    ; label("Multiplicative factor on the baseline logit for first cigarette within 5 min of waking, vs >60 min (unitless)")  # Ravva 2010 Table 3 'Effect of FSQ1: <=5min' = 0.961 (10.3% RSE); theta5

    # ----- Age: power function of AGE/45 on the intercept -----
    e_age_base_logit              <- 0.374    ; label("Power exponent of (AGE/45) on the baseline logit (unitless)")  # Ravva 2010 Table 3 'Effect of age, 18-75 years' = 0.374 (26.0% RSE); theta6

    # ----- Sex and race: power-of-indicator factors on the intercept -----
    e_sexf_base_logit             <- 0.704    ; label("Multiplicative factor on the baseline logit for female sex, vs male (unitless)")  # Ravva 2010 Table 3 'Effect of gender: Female' = 0.704 (5.4% RSE); theta7
    e_race_black_base_logit       <- 1.20     ; label("Multiplicative factor on the baseline logit for Black race, vs White (unitless)")  # Ravva 2010 Table 3 'Effect of race: Black' = 1.20 (8.3% RSE); theta8
    e_race_other_base_logit       <- 0.887    ; label("Multiplicative factor on the baseline logit for Other race, vs White (unitless)")  # Ravva 2010 Table 3 'Effect of race: Other' = 0.887 (10.4% RSE); theta9

    # ----- Weekly time-course term: printed in Equation 2, values NEVER published -----
    # Ravva 2010 p.340: "This effect of time on the baseline probability of
    # nausea was modeled as an exponential function with two additional
    # parameters (Equation 2; theta_10 and theta_11)." Neither theta10 nor
    # theta11 appears in Table 3, in the supplements, or anywhere else in the
    # paper -- Table 3's nausea column is the WHOLE-TREATMENT-PERIOD
    # incidence model (one observation per subject), which is confirmed
    # arithmetically: all six reported nausea probabilities reproduce from
    # Table 3 with no time term active.
    #
    # Per the standing "unreported -> fixed(0) + erratum" policy the term is
    # encoded verbatim with both parameters fixed at 0, so the printed
    # equation structure stays visible at the parameter line and a user who
    # later obtains the values can set them without editing model(). With
    # logit_week_amp = 0 the whole term vanishes, so the model is numerically
    # identical to omitting it and every published answer key still
    # reproduces exactly. Note logit_week_kdes is unidentifiable once
    # logit_week_amp is 0 (it sits inside exp()); 0 is an arbitrary
    # placeholder, not a source value. See vignette Errata.
    logit_week_amp                <- fixed(0) ; label("Amplitude of the exponentially decaying week-of-treatment term on the nausea logit (unitless logit); value not published")  # Ravva 2010 Equation 2 theta11 -- never reported; see Errata
    logit_week_kdes               <- fixed(0) ; label("First-order decay rate constant of the week-of-treatment term on the nausea logit (1/week); value not published")  # Ravva 2010 Equation 2 theta10 -- never reported; see Errata

    # ----- Between-subject random effect on the logit -----
    # Equation 2 carries an additive eta on the logit for the repeated-measures
    # weekly analysis. Its variance is never reported, and the source authors
    # themselves rejected the mixed-effects fit: p.340, "the histogram plot of
    # random effects ... revealed a bimodal distribution, confirming the
    # inadequacy of mixed-effects modeling for dichotomous data", and Methods,
    # "a naive pooled analysis was also performed, fixing eta at 0". Fixed at 0
    # here, matching both the standing unreported-IIV policy and the authors'
    # own reported choice.
    etabase_logit ~ fixed(0)  # Ravva 2010 Equation 2 eta1; variance never reported, and the published results are the naive pooled fit in which the authors held eta at 0

    # ----- Residual error (placeholder; NOT from the source) -----
    # The source fit is a NONMEM Laplacian Bernoulli likelihood on the
    # dichotomous nausea outcome and reports no $SIGMA. rxode2 / nlmixr2
    # expose llikBinom / rxbinom but no Bernoulli ENDPOINT that parses in
    # model(), so the deterministic typical-value probability p_nausea is the
    # declared output and a tiny placeholder additive residual satisfies the
    # observation declaration. Mirrors Oniki_2018_nafld_risk.R.
    addSd_p_nausea                <- fixed(0.001) ; label("Placeholder additive residual SD on p_nausea; the source likelihood is Bernoulli, so no residual-error term is published")  # not paper-derived; see Assumptions and deviations
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

    # ----- Full Equation 2 logit -----
    # Model time is in weeks (see units$time), so `time` IS the "week" term of
    # Equation 2. With logit_week_amp fixed at 0 this contributes exactly 0
    # and the expression reduces to the published incidence model.
    logit_nausea <- base_logit_i +
                    e_auc_varen_base_logit * AUC_VAREN +
                    logit_week_amp * exp(-logit_week_kdes * time) +
                    etabase_logit

    # ----- Probability of a nausea event -----
    p_nausea <- expit(logit_nausea)

    # Deterministic typical-value declaration with a placeholder residual
    # (see the ini() comment). Downstream callers can use p_nausea directly
    # as a probability, or draw Bernoulli outcomes externally with
    # rbinom(n, 1, p_nausea) on the rxSolve output.
    p_nausea ~ add(addSd_p_nausea)
  })
}
