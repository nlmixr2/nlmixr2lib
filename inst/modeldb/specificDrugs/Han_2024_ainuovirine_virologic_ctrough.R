Han_2024_ainuovirine_virologic_ctrough <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model relating the steady-state ",
    "trough concentration of the non-nucleoside reverse transcriptase ",
    "inhibitor ainuovirine (ANV) to virologic suppression -- HIV-RNA < 50 ",
    "copies/mL at week 48 -- in antiretroviral-therapy-naive adults living ",
    "with HIV-1 (Han 2024, the phase 3 trial ADYY-ACC007-301, ",
    "ChiCTR1800019041, 150 mg once daily at bedtime for 48 weeks). The ",
    "probability of suppression is expit(1.65 + 0.0038 * CTROUGH), where ",
    "CTROUGH is the individual steady-state trough ANV concentration in ",
    "ng/mL derived by Bayesian post-hoc estimation from the companion ",
    "population PK model (modellib('Han_2024_ainuovirine')). THE SLOPE IS ",
    "NOT SIGNIFICANT (P = 0.220): Han 2024's own conclusion is that the ",
    "exposure-efficacy relationship is FLAT over the studied exposure ",
    "range, adding nothing to a constant-probability model, which the ",
    "paper attributes to the antiviral effect having reached a plateau. ",
    "The model is packaged because a flat, non-significant slope with a ",
    "published intercept is itself a reusable quantitative result -- it ",
    "bounds how much virologic benefit further exposure could buy -- but ",
    "it must not be read as evidence that raising ANV exposure raises ",
    "efficacy. Contrast the companion safety models, in which the same ",
    "exposure metrics ARE significantly associated with adverse drug ",
    "reactions. There is no PK layer and no ODE: the exposure metric is ",
    "supplied as a data column. Three companion exposure-response models ",
    "in the Han_2024_ainuovirine_* family."
  )
  reference <- paste(
    "Han X, Sun J, Zhang Y, Jiang T, Zheng Q, Peng H, Wang Y, Xia W,",
    "Zhang T, Sun L, Yun X, Qin H, Wu H, Su B.",
    "Population pharmacokinetics of Ainuovirine and exposure-response",
    "analysis in human immunodeficiency virus-infected individuals.",
    "Chin Med J (Engl). 2024;137(20):2474-2482.",
    "doi:10.1097/CM9.0000000000002917.",
    "Individual exposure metrics derive from the companion population",
    "pharmacokinetic model reported in the same paper; see",
    "modellib('Han_2024_ainuovirine').",
    sep = " "
  )
  vignette <- "Han_2024_ainuovirine"
  units <- list(
    time          = "n/a (static landmark exposure-efficacy regression at week 48; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_hivrna_lt50 (probability of HIV-RNA < 50 copies/mL at week 48, 0-1; also logit_hivrna_lt50)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual steady-state trough plasma concentration of ainuovirine (Ctrough), per subject. Supplied as data: this model has no PK layer.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration at STEADY STATE. Han",
        "2024 Methods: 'Individual steady-state trough concentration",
        "(Ctrough) and area under the steady-state curve (AUCtau) were",
        "exposure factors evaluated by Bayesian analysis' -- i.e. these are",
        "empirical-Bayes post-hoc predictions from the companion population",
        "PK model in the same paper, not observed trough measurements.",
        "Reproduce them with modellib('Han_2024_ainuovirine') solved at",
        "150 mg once daily with MULTI_DOSE_PT = 1.",
        "Enters LINEARLY, not on a log scale: Han 2024 states 'Logistic",
        "regression indicated a linear effect of Ctrough, AUCtau, and",
        "curative effect'. The printed coefficient is therefore per ng/mL,",
        "not per e-fold.",
        "UNITS ARE LOAD-BEARING and Han 2024 never states them on Table 4;",
        "they are established from the paper's own unit system, in which",
        "every ANV concentration is reported in ng/mL (Figure 1A/1B y-axes,",
        "and the Table 3 additive residual 8.89 ng/mL). Supplying ug/mL",
        "instead would shrink the linear-predictor contribution by 1000-fold",
        "and collapse the model to its intercept. A cross-check confirms the",
        "scale: at the exposure the companion PK model predicts for the",
        "phase 3 regimen, this model and its AUCtau sibling agree closely",
        "on the ADR endpoint, which they could not do if either metric were",
        "off by a factor of 1000.",
        "Han 2024 does not tabulate the observed Ctrough distribution; the",
        "quartile boundaries used for Figure 4C are not printed."
      ),
      source_name        = "Ctrough (steady-state trough concentration)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 309L,
    n_studies      = 1L,
    n_observations = "One binary virologic-response record per subject at week 48. NOTE: Han 2024 never states the exposure-response analysis-set size. 309 is the number of phase 3 ADYY-ACC007-301 participants in the population PK dataset (Table 1, Table 2); the week-48 endpoint exists only in that trial, so the phase 1 subjects cannot contribute. The true analysis set may be smaller if any participant lacked a week-48 HIV-RNA result",
    age_range      = "mean 31.03 years (SD 9.95) in the phase 3 trial (Han 2024 Table 2)",
    weight_range   = "mean 67.06 kg (SD 11.42) in the phase 3 trial (Han 2024 Table 2)",
    sex_female_pct = 5.83,
    race_ethnicity = "Not reported; single Chinese centre (Beijing Youan Hospital, Capital Medical University)",
    disease_state  = "Antiretroviral-therapy-naive people living with HIV-1, randomised to an ainuovirine-based regimen (the comparator arm received an efavirenz-based regimen and is not part of this exposure-response analysis)",
    dose_range     = "Ainuovirine 150 mg orally once daily at bedtime on an empty stomach for 48 weeks",
    regions        = "Single centre, Beijing, China",
    notes          = paste0(
      "The phase 3 trial's primary objective was non-inferiority of the ",
      "ainuovirine regimen against the efavirenz regimen for the ",
      "proportion of participants reaching HIV-RNA < 50 copies/mL at ",
      "week 48. Because every participant in this analysis received the ",
      "same 150 mg dose, the exposure range driving the regression is ",
      "generated purely by between-subject pharmacokinetic variability ",
      "(omega(CL) = 30.9% in the companion PK model) rather than by dose ",
      "randomisation -- which is the usual reason a landmark ",
      "exposure-efficacy slope in a single-dose-level trial has little ",
      "power, and is worth keeping in view alongside the paper's ",
      "'effect has reached a plateau' reading."
    )
  )

  ini({
    # ==================================================================
    # Han 2024 Table 4, upper block "Exposure (Ctrough)-virological
    # response model". The model is a univariate binomial logistic
    # regression of the week-48 virologic-suppression indicator on the
    # individual steady-state trough:
    #
    #   logit(p) = beta0 + beta1 * Ctrough
    #
    # Table 4 prints an Estimate, an SD and a P-value per row and nothing
    # else -- no odds ratio, no confidence interval -- so unlike the
    # Chen 2021 lorlatinib family there is no redundant column against
    # which to cross-check each coefficient. Each value below is the
    # Table 4 Estimate column, and the SDs are recorded in the comments
    # for anyone who needs the uncertainty.
    #
    # The slope's P-value is quoted twice and consistently: Table 4 gives
    # 0.220, and the Results text gives "(slope only: P = 0.220; ...)"
    # for this model.
    #
    # The covariate is NOT centred and NOT scaled, so the intercept is
    # the logit at Ctrough = 0 ng/mL -- an extrapolated anchor outside
    # the observed exposure range, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- 1.65; label("Logit of the probability of HIV-RNA < 50 copies/mL at week 48 at CTROUGH = 0 ng/mL (unitless logit)")  # Han 2024 Table 4, beta 0 (Intercept) = 1.65, SD 0.55, P = 0.003

    # ----- Exposure effect on the logit -----
    # Printed as "0.38x10 -2" in Table 4, i.e. 0.38e-2 = 0.0038.
    e_ctrough_logit <- 0.0038; label("Log-odds of HIV-RNA < 50 copies/mL at week 48 per 1 ng/mL increase in steady-state ainuovirine trough concentration (unitless logit)")  # Han 2024 Table 4, beta 1 (Slope) = 0.38x10^-2, SD 0.31x10^-2, P = 0.220 -- NOT significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and no random effects are estimated. The
    # tiny fixed additive residual below exists only so rxode2 has an
    # error model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_hivrna_lt50 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_hivrna_lt50 <- logit_ref + e_ctrough_logit * CTROUGH

    prob_hivrna_lt50 <- expit(logit_hivrna_lt50)

    prob_hivrna_lt50 ~ add(addSd_prob_hivrna_lt50)
  })
}
