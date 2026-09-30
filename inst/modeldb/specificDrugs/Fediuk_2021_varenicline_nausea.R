Fediuk_2021_varenicline_nausea <- function() {
  description <- paste(
    "Logistic exposure-response model for the incidence of nausea or",
    "vomiting over the 12-week treatment period of the phase 4 varenicline",
    "study in adolescent smokers aged 12-20 years (Fediuk 2021, final",
    "model, one observation per subject). The logit is a baseline intercept",
    "multiplied by covariate factors for nicotine dependence (Fagerstrom",
    "time to first cigarette, 4-level score), age (power on AGE/16), female",
    "sex and race, PLUS an additive linear term in the individual",
    "steady-state varenicline AUC(0-24). Exposure is a data covariate",
    "computed as total daily dose / individual CL/F from the companion",
    "population PK model Fediuk_2021_varenicline.R, and is 0 for placebo.",
    "Nausea or vomiting is about 86% more likely in female adolescents."
  )
  reference <- paste(
    "Fediuk DJ, Sweeney K, Sahasrabudhe V, McRae T, Byon W.",
    "Population pharmacokinetics and exposure-response analyses of",
    "varenicline in adolescent smokers.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(7):769-781.",
    "doi:10.1002/psp4.12645"
  )
  vignette <- "Fediuk_2021_varenicline"
  units <- list(
    time = "week",
    dosing = "n/a (exposure-response model; varenicline exposure enters as the AUC_VAREN covariate, not as a dosing event)",
    concentration = "(probability, 0-1; prob_nausea_vomiting is the probability of at least one treatment-emergent nausea or vomiting event over the 12-week treatment period)"
  )

  covariateData <- list(
    AUC_VAREN = list(
      description = "Individual varenicline steady-state daily exposure, AUC(0-24)",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Fediuk 2021 Methods: 'The daily AUC24 was estimated from the empirical Bayes predictions of CL/F value and total daily dose for each subject. AUC24 was set to zero for subjects in the placebo group.' Compute as 1000 * (total daily dose, mg) / (individual CL/F, L/h) from Fediuk_2021_varenicline.R. Enters the logit ADDITIVELY (Text S2 '+THETA(2)*AUC'), unlike the demographic factors, which multiply the intercept. Figure 4 uses a representative AUC24 of 153 ng*h/mL.",
      source_name = "AUC"
    ),
    SMOKE_TTFC_SCORE = list(
      description = "Fagerstrom Test for Nicotine Dependence item 1 ('How soon after you wake up do you smoke your first cigarette?') scored 0-3: >60 min (0); 31-60 min (1); 6-30 min (2); within 5 min (3)",
      units = "(ordinal score 0-3)",
      type = "categorical",
      reference_category = "0 = first cigarette more than 60 min after waking (the model reference level; Methods 'zero was equal to >60 minutes')",
      notes = "Text S2 decomposes FSQ1 into three indicators, 'IF(FSQ1.EQ.1) FQ2=1', 'IF(FSQ1.EQ.2) FQ3=1', 'IF(FSQ1.EQ.3) FQ4=1', each with its own multiplicative factor on the intercept; same 0-3 orientation as the canonical, so no value transformation. The model reference level (score 0) held only 5 of 238 subjects (Table 2), and all three non-reference factors are about 0.075-0.078, so the intercept -30.9 is effectively a scale term and the covariate-adjusted baseline logit of a typical subject is about -2.3 to -2.4. Figure 4 reports ratios against a representative subject with score 1 (31-60 min), not the model reference level.",
      source_name = "FSQ1"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function (AGE/16)^e on the intercept (Text S2 '(AGE/16)**THETA(6)'). Reference 16 years. ER cohort median 16 years, range 12-20 (Table 2).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Text S2 'NSEX = SEX - 1' from a SEX column coded 1 = male, 2 = female, so NSEX is identically SEXF. The factor 0.554 on a negative intercept raises the probability: Figure 4 ratio about 1.86 at the representative subject (Results: 'a significant increase of ~86%').",
      source_name = "NSEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White",
      notes = "Text S2 'IF(RACE.EQ.2) RC2=1'. 16 of 238 ER subjects (Table 2).",
      source_name = "RC2 (RACE2)"
    ),
    RACE_OTHER = list(
      description = "Composite 'Other' race indicator pooling Asian, Hispanic, American Indian and mixed race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White",
      notes = "Text S2 'IF(RACE.GE.3) RC3=1'; Results: 'For consistency with the popPK analysis, subjects of Asian and other race were grouped together.' 46 of 238 ER subjects (Table 2: Asian 44 + Other 2).",
      source_name = "RC3 (RACE3)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 238L,
    n_studies = 1L,
    age_range = "12-20 years",
    age_median = "16 years",
    weight_range = NA_character_,
    sex_female_pct = 34.0,
    race_ethnicity = c(White = 73.9, Black = 6.72, Asian = 18.5, Other = 0.840),
    disease_state = "Nicotine-dependent adolescent smokers (Fagerstrom Test for Nicotine Dependence score >= 4, >= 5 cigarettes per day, >= 1 prior failed quit attempt) motivated to stop smoking.",
    dose_range = "Oral varenicline 0.5 mg q.d. or 0.5 mg b.i.d. (body weight <= 55 kg) and 0.5 mg b.i.d. or 1 mg b.i.d. (body weight > 55 kg), or placebo, for 12 weeks after a 1- or 2-week up-titration.",
    regions = "Multicenter phase 4 study NCT01312909.",
    smoking_dependence = "Time to first cigarette (FSQ1): within 5 min 98 (41.2%); 6-30 min 96 (40.3%); 31-60 min 39 (16.4%); >60 min 5 (2.10%) (Table 2).",
    notes = "Phase 4 study only (Fediuk 2021 Table 2, 'ER analyses' column): 238 subjects with one observation each, 99 placebo (AUC24 = 0) and 139 varenicline-treated subjects with at least one measurable concentration. Nausea and vomiting treatment-emergent adverse events from the first to the last dose of study medication, analysed as a naive-pooled incidence (one observation per subject). NONMEM Laplacian LIKELIHOOD estimation."
  )

  ini({
    # All values are the final estimates of Fediuk 2021 Table S2 and match the
    # Text S2 (run83.mod) $THETA block. Covariates multiply the intercept;
    # the exposure term is additive on the logit (Methods equation; Text S2
    # 'LGT = LGT2*THETA(9)**RC3+THETA(2)*AUC+ETA(1)').
    # Reference: 16-year-old White male, FSQ1 score 0 (>60 min), AUC24 = 0.
    base_logit <- -30.9; label("Baseline logit of the nausea-or-vomiting probability at the model reference covariates (unitless logit)") # Table S2 'Intercept (theta1)' = -30.9 (RSE 1.71%)
    e_auc_varen_base_logit <- 0.00911; label("Slope of the logit on varenicline AUC(0-24) (logit units per ng*h/mL)") # Table S2 'Effect of AUC (theta2)' = 0.00911 (RSE 27.4%)

    e_smoke_ttfc_31_60_base_logit <- 0.0750; label("Multiplicative factor on the baseline logit for first cigarette 31-60 min after waking, vs >60 min (unitless)") # Table S2 'Effect of FSQ1 31-60 min (theta3)' = 0.0750 (RSE 29.7%)
    e_smoke_ttfc_6_30_base_logit <- 0.0773; label("Multiplicative factor on the baseline logit for first cigarette 6-30 min after waking, vs >60 min (unitless)") # Table S2 'Effect of FSQ1 6-30 min (theta4)' = 0.0773 (RSE 18.0%)
    e_smoke_ttfc_le5_base_logit <- 0.0779; label("Multiplicative factor on the baseline logit for first cigarette within 5 min of waking, vs >60 min (unitless)") # Table S2 'Effect of FSQ1 within 5 min (theta5)' = 0.0779 (RSE 16.9%)

    e_age_base_logit <- -0.296; label("Power exponent of (AGE/16) on the baseline logit (unitless)") # Table S2 'Effect of age, years (theta6)' = -0.296 (RSE 245%)
    e_sexf_base_logit <- 0.554; label("Multiplicative factor on the baseline logit for female sex, vs male (unitless)") # Table S2 'Effect of sex (female) (theta7)' = 0.554 (RSE 21.5%)
    e_race_black_base_logit <- 0.870; label("Multiplicative factor on the baseline logit for Black race, vs White (unitless)") # Table S2 'Effect of race (black) (theta8)' = 0.870 (RSE 33.6%)
    e_race_other_base_logit <- 0.966; label("Multiplicative factor on the baseline logit for Other race, vs White (unitless)") # Table S2 'Effect of race (other) (theta9)' = 0.966 (RSE 22.0%)

    # Text S2 carries an additive ETA(1) on the logit with '$OMEGA 0 FIX'
    # (naive-pooled fit, one observation per subject).
    etabase_logit ~ fixed(0) # Text S2 $OMEGA 0 FIX

    # Placeholder residual (NOT from the source): the source likelihood is
    # Bernoulli (Text S2 'Y = P**DV*(1-P)**(1-DV)'), which has no residual-error
    # parameter; a tiny fixed additive SD lets the deterministic probability be
    # declared as the model output.
    addSd <- fixed(0.001); label("Placeholder additive residual SD on prob_nausea_vomiting; not from the source (Bernoulli likelihood)") # not paper-derived; see vignette Assumptions and deviations
  })

  model({
    # FSQ1 multiplier: three non-reference indicators (Text S2 FQ2 / FQ3 / FQ4)
    ttfc_mult <- e_smoke_ttfc_31_60_base_logit^(SMOKE_TTFC_SCORE == 1) *
      e_smoke_ttfc_6_30_base_logit^(SMOKE_TTFC_SCORE == 2) *
      e_smoke_ttfc_le5_base_logit^(SMOKE_TTFC_SCORE == 3)

    # Covariate-adjusted baseline logit (Text S2 LGT1, LGT2)
    base_logit_i <- base_logit *
      ttfc_mult *
      (AGE / 16)^e_age_base_logit *
      e_sexf_base_logit^SEXF *
      e_race_black_base_logit^RACE_BLACK *
      e_race_other_base_logit^RACE_OTHER

    # Full logit with the additive exposure term (Text S2 LGT)
    logit_nausea_vomiting <- base_logit_i + e_auc_varen_base_logit * AUC_VAREN + etabase_logit

    prob_nausea_vomiting <- expit(logit_nausea_vomiting)
    prob_nausea_vomiting ~ add(addSd)
  })
}
