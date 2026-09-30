Fediuk_2021_varenicline_car_w9_12 <- function() {
  description <- paste(
    "Logistic exposure-response model for the biochemically confirmed",
    "continuous abstinence rate at weeks 9-12 (CAR9-12) in the phase 4",
    "varenicline study in adolescent smokers aged 12-20 years (Fediuk 2021).",
    "Base model only: intercept plus a linear term in the individual",
    "steady-state varenicline AUC(0-24); the slope was not statistically",
    "significant (p = 0.303) and no covariates were added. The paper does",
    "not print the two estimates; they were recovered by digitizing the",
    "vector 'Model predicted' line of Figure 3a, so treat them as figure",
    "approximations. Exposure is a data covariate computed as total daily",
    "dose / individual CL/F from Fediuk_2021_varenicline.R, 0 for placebo."
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
    concentration = "(probability, 0-1; p_car is the probability of biochemically confirmed continuous abstinence over weeks 9-12)"
  )

  covariateData <- list(
    AUC_VAREN = list(
      description = "Individual varenicline steady-state daily exposure, AUC(0-24)",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Fediuk 2021 Methods: 'The daily AUC24 was estimated from the empirical Bayes predictions of CL/F value and total daily dose for each subject. AUC24 was set to zero for subjects in the placebo group.' Compute as 1000 * (total daily dose, mg) / (individual CL/F, L/h) from Fediuk_2021_varenicline.R. The Figure 3a line spans the observed exposure range 0 to about 314 ng*h/mL; do not extrapolate beyond it.",
      source_name = "AUC24"
    )
  )

  # Screened but never added: Fediuk 2021 Results, 'there was no statistically
  # significant trend (slope estimate; p = 0.303) and no further model
  # development was performed (Figure 3a)'. Documentation only.
  covariatesDataExcluded <- list(
    SMOKE_TTFC_SCORE = list(
      description = "Fagerstrom Test for Nicotine Dependence item 1 (time to first cigarette) scored 0-3",
      units = "(ordinal score 0-3)",
      type = "categorical",
      reference_category = "0 = first cigarette more than 60 min after waking",
      notes = "Part of the planned covariate set (FSQ1, age, sex, race) but not tested for CAR9-12 because the base exposure-response model showed no significant AUC24 effect.",
      source_name = "FSQ1"
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
    disease_state = "Nicotine-dependent adolescent smokers (Fagerstrom Test for Nicotine Dependence score >= 4, >= 5 cigarettes per day, >= 1 prior failed quit attempt) motivated to stop smoking, receiving brief age-appropriate cessation counselling at every visit.",
    dose_range = "Oral varenicline 0.5 mg q.d. or 0.5 mg b.i.d. (body weight <= 55 kg) and 0.5 mg b.i.d. or 1 mg b.i.d. (body weight > 55 kg), or placebo, for 12 weeks after a 1- or 2-week up-titration.",
    regions = "Multicenter phase 4 study NCT01312909.",
    notes = "Phase 4 study only (Fediuk 2021 Table 2, 'ER analyses' column): 238 subjects with one quit / non-quit observation each, 99 placebo (AUC24 = 0) and 139 varenicline-treated subjects with at least one measurable concentration. Self-reported abstinence was confirmed by urine cotinine; subjects who did not complete treatment were non-responders from discontinuation onward."
  )

  ini({
    # Fediuk 2021 reports only the slope p-value (0.303) for this model; Table
    # S2 and Text S2 cover the nausea/vomiting model only. Both values below
    # were recovered by digitizing the vector 'Model predicted' dotted line of
    # Figure 3a (PDF page 10, extracted with pdftocairo -svg, axes mapped on
    # the printed ticks). The drawn line is straight in probability from
    # p = 0.192 at AUC24 = 0 to p = 0.320 at AUC24 = 314 ng*h/mL; the logistic
    # below is the least-squares fit on the logit scale to the 99 digitized
    # dash centres and stays within 0.003 probability of the line over that
    # range. See the vignette Assumptions and deviations.
    base_logit <- -1.42; label("Logit of the CAR9-12 probability at AUC(0-24) = 0 (unitless logit)") # digitized from Figure 3a 'Model predicted' line; not printed in the paper
    e_auc_varen_base_logit <- 0.00211; label("Slope of the logit on varenicline AUC(0-24) (logit units per ng*h/mL)") # digitized from Figure 3a 'Model predicted' line; not printed (Results: slope p = 0.303)

    # Placeholder residual (NOT from the source): the source likelihood is
    # Bernoulli, which has no residual-error parameter.
    addSd <- fixed(0.001); label("Placeholder additive residual SD on p_car; not from the source (Bernoulli likelihood)") # not paper-derived; see vignette Assumptions and deviations
  })

  model({
    logit_car <- base_logit + e_auc_varen_base_logit * AUC_VAREN
    p_car <- expit(logit_car)
    p_car ~ add(addSd)
  })
}
