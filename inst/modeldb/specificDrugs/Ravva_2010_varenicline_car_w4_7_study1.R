Ravva_2010_varenicline_car_w4_7_study1 <- function() {
  description <- paste0(
    "Preliminary population exposure-response logistic-regression model ",
    "for the continuous abstinence rate (CAR) at weeks 4-7 in 490 adult ",
    "cigarette smokers from varenicline dose-ranging study 1 (Ravva ",
    "2010). The logit of the quit probability is an additive intercept ",
    "plus a linear term in the individual steady-state varenicline ",
    "exposure AUC(0-24)ss, with no covariates. This early-development ",
    "analysis supported dosage selection for the subsequent smoking-",
    "cessation trials; it was fitted separately in each of two ",
    "dose-ranging studies to demonstrate that the intercept and slope ",
    "reproduce across studies. See the companion study-2 fit in ",
    "Ravva_2010_varenicline_car_w4_7_study2.R and the final ",
    "covariate-bearing weeks 9-12 model in ",
    "Ravva_2010_varenicline_car_w9_12.R. Exposure is a data covariate ",
    "computed as daily dose / (CL/F) from the companion population PK ",
    "model (see Ravva_2009_varenicline.R) and is set to 0 for placebo ",
    "subjects. Fitted in NONMEM V with the Laplacian method on the ",
    "dichotomous quit / no-quit outcome; no between-subject random ",
    "effect and no residual-error term are reported."
  )
  reference <- paste(
    "Ravva P, Gastonguay MR, French JL, Tensfeldt TG, Faessel HM.",
    "Quantitative assessment of exposure-response relationships for the",
    "efficacy and tolerability of varenicline for smoking cessation.",
    "Clin Pharmacol Ther. 2010;87(3):336-344. doi:10.1038/clpt.2009.282.",
    "Study 1 is the 7-week dose-ranging trial reported in:",
    "Nides M, Oncken C, Gonzales D, Rennard S, Watsky EJ, Anziano R,",
    "Reeves KR. Smoking cessation with varenicline, a selective alpha-4",
    "beta-2 nicotinic receptor partial agonist.",
    "Arch Intern Med. 2006;166(15):1561-1568.",
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
    concentration = "p_car (probability of continuous abstinence at weeks 4-7, 0-1; also logit_car)"
  )

  covariateData <- list(
    AUC_VAREN = list(
      description        = "Individual varenicline steady-state daily exposure, AUC(0-24)ss. Ravva 2010 Methods (Pharmacokinetics): 'Individual 24-h daily exposure, measured as AUC(0-24)ss, was calculated as dose divided by CL/F_i (apparent clearance); the individual empirical Bayes estimate of apparent clearance was predicted from the final population PK model and parameters obtained from a pooled analysis in adult smokers.'",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only predictor in this preliminary model; it enters the logit additively alongside the intercept. Set to 0 for placebo subjects, consistent with the treatment of exposure throughout Ravva 2010. Study 1 randomized subjects to varenicline 0.3 mg q.d., 1 mg q.d. and 1 mg b.i.d., plus a bupropion SR comparator and placebo (Ravva 2010 Supplementary Table S1), so this fit spans a wider and lower daily-dose range than the 0.5 / 1 mg b.i.d. studies. Downstream users should compute the per-subject value from the companion Ravva 2009 varenicline population PK model as total daily dose / (CL/F).",
      source_name        = "AUC(0-24)ss"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 490L,
    n_studies      = 1L,
    age_range      = "19-66 years (pooled weeks 4-7 analysis database)",
    age_mean       = "43 years (pooled weeks 4-7 analysis database)",
    weight_range   = NA_character_,
    sex_female_pct = 51.0,
    race_ethnicity = c(White = 83.0, Black = 11.0, Other = 6.0),
    disease_state  = "Adult cigarette smokers motivated to stop smoking.",
    dose_range     = "Varenicline 0.3 mg q.d., 1 mg q.d. and 1 mg b.i.d. oral, plus bupropion SR 150 mg b.i.d. and placebo, over a 7-week treatment period.",
    regions        = "Multicenter.",
    smoking_marker = "Baseline exhaled carbon monoxide mean 24 ppm, range 1-98 ppm (pooled weeks 4-7 analysis database).",
    notes          = "Study 1 is the 7-week, randomized, double-blind, placebo-controlled, parallel-group, multicenter dose-ranging trial of Ravva 2010 Supplementary Table S1, which enrolled 626 male and female smokers; 490 contributed to the weeks 4-7 CAR analysis (Ravva 2010 Table 1). Baseline demographics in Ravva 2010 Table 1 are reported for the POOLED weeks 4-7 database (studies 1 and 2 together, n = 1,099), not per study, so the demographic entries above describe that pooled database and not study 1 alone. Pharmacokinetic sampling at weeks 1, 2 and 4 (or early termination)."
  )

  ini({
    # =====================================================================
    # Ravva 2010 Results, 'Population exposure-response relationships for
    # estimation of efficacy', p.337 verbatim: "comparable logit probability
    # response parameters (intercepts and slopes) were estimated with good
    # precision for the CAR end point for weeks 4-7 in each of the two
    # studies (typical value (% relative standard error (RSE))): study 1
    # (intercept, -1.44 (10.8%) and slope, 0.00623 (22.5%)); study 2
    # (intercept, -1.30 (11.6%) and slope, 0.00543 (19.9%))."
    #
    # Model shape: unlike the final weeks 9-12 model, this preliminary fit
    # has NO covariates, so the intercept enters ADDITIVELY on the logit
    # (the same shape as the Table 2 base-plus-AUC model,
    # logit(P) = theta1 + theta2 * AUC), not multiplicatively.
    # =====================================================================

    base_logit             <- -1.44    ; label("Baseline logit of the weeks 4-7 quit probability at zero varenicline exposure (unitless logit)")  # Ravva 2010 Results p.337: study 1 intercept = -1.44 (10.8% RSE)
    e_auc_varen_base_logit <- 0.00623  ; label("Slope of the logit on varenicline AUC(0-24)ss (logit units per ng*h/mL)")  # Ravva 2010 Results p.337: study 1 slope = 0.00623 (22.5% RSE)

    # ----- Residual error (placeholder; NOT from the source) -----
    # The source fit is a NONMEM Laplacian Bernoulli likelihood on the
    # dichotomous quit outcome and reports no $SIGMA. rxode2 / nlmixr2 expose
    # llikBinom / rxbinom but no Bernoulli ENDPOINT that parses in model(),
    # so the deterministic typical-value probability p_car is the declared
    # output and a tiny placeholder additive residual satisfies the
    # observation declaration. Mirrors Oniki_2018_nafld_risk.R.
    addSd_p_car            <- fixed(0.001) ; label("Placeholder additive residual SD on p_car; the source likelihood is Bernoulli, so no residual-error term is published")  # not paper-derived; see Assumptions and deviations
  })

  model({
    # Additive intercept plus linear exposure effect on the logit.
    logit_car <- base_logit + e_auc_varen_base_logit * AUC_VAREN

    # Probability of continuous abstinence at weeks 4-7.
    p_car <- expit(logit_car)

    # Deterministic typical-value declaration with a placeholder residual
    # (see the ini() comment).
    p_car ~ add(addSd_p_car)
  })
}
