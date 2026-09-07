Han_2024_ainuovirine_virologic_auctau <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model relating the steady-state ",
    "dosing-interval exposure (AUCtau) of the non-nucleoside reverse ",
    "transcriptase inhibitor ainuovirine (ANV) to virologic suppression -- ",
    "HIV-RNA < 50 copies/mL at week 48 -- in antiretroviral-therapy-naive ",
    "adults living with HIV-1 (Han 2024, the phase 3 trial ",
    "ADYY-ACC007-301, ChiCTR1800019041, 150 mg once daily at bedtime for ",
    "48 weeks). The probability of suppression is expit(1.44 + 0.000012 * ",
    "AUC_ANV), where AUC_ANV is the individual steady-state AUC over the ",
    "24 h dosing interval in ng*h/mL derived by Bayesian post-hoc ",
    "estimation from the companion population PK model ",
    "(modellib('Han_2024_ainuovirine')). THE SLOPE IS NOT SIGNIFICANT ",
    "(P = 0.222): Han 2024's own conclusion is that the exposure-efficacy ",
    "relationship is FLAT over the studied exposure range, adding nothing ",
    "to a constant-probability model, which the paper attributes to the ",
    "antiviral effect having reached a plateau. The model is packaged ",
    "because a flat, non-significant slope with a published intercept is ",
    "itself a reusable quantitative result -- it bounds how much virologic ",
    "benefit further exposure could buy -- but it must not be read as ",
    "evidence that raising ANV exposure raises efficacy. This is the ",
    "AUCtau-parameterised sibling of ",
    "Han_2024_ainuovirine_virologic_ctrough: Han 2024 fitted the two ",
    "exposure metrics as separate univariate regressions and reported ",
    "both, so they are alternatives rather than a joint model, and must ",
    "not be applied together. There is no PK layer and no ODE: the ",
    "exposure metric is supplied as a data column. Three companion ",
    "exposure-response models in the Han_2024_ainuovirine_* family."
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
    dosing        = "n/a (no dose events; exposure enters as the AUC_ANV covariate column)",
    concentration = "prob_hivrna_lt50 (probability of HIV-RNA < 50 copies/mL at week 48, 0-1; also logit_hivrna_lt50)"
  )

  covariateData <- list(
    AUC_ANV = list(
      description        = "Individual steady-state area under the ainuovirine plasma concentration-time curve over the 24 h once-daily dosing interval (AUCtau), per subject. Supplied as data: this model has no PK layer.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma exposure at STEADY STATE, over the",
        "24 h dosing interval of a once-daily regimen (so AUCtau is",
        "AUC(0-24) at steady state, and equals AUCss numerically for this",
        "regimen). Han 2024 Methods: 'Individual steady-state trough",
        "concentration (Ctrough) and area under the steady-state curve",
        "(AUCtau) were exposure factors evaluated by Bayesian analysis' --",
        "i.e. these are empirical-Bayes post-hoc predictions from the",
        "companion population PK model in the same paper, not",
        "noncompartmental estimates from observed profiles. Reproduce them",
        "with modellib('Han_2024_ainuovirine') solved at 150 mg once daily",
        "with MULTI_DOSE_PT = 1; the closed form is F * dose / CL/F with",
        "the steady-state clearance, i.e. 0.716 * 150 mg / 15.96 L/h.",
        "Enters LINEARLY, not on a log scale: Han 2024 states 'Logistic",
        "regression indicated a linear effect of Ctrough, AUCtau, and",
        "curative effect'. The printed coefficient is therefore per",
        "ng*h/mL.",
        "UNITS ARE LOAD-BEARING and Han 2024 never states them on Table 4.",
        "ng*h/mL is established from the paper's own unit system: Figure 1D",
        "plots dose-normalised AUCinf on an axis labelled",
        "'h*ng*mL^-1*mg^-1', so the un-normalised AUC is in ng*h/mL, and",
        "every ANV concentration in the paper is in ng/mL. Two cross-checks",
        "agree. (1) Scale: the ratio of this model family's Ctrough slope",
        "to its AUCtau slope is 0.40e-2 / 0.14e-3 = 28.6 h on the ADR",
        "endpoint, the right order for an AUCtau/Ctrough ratio over a 24 h",
        "interval, which it would not be if the two metrics differed by a",
        "factor of 1000. (2) Agreement: at the exposure the companion PK",
        "model predicts for the phase 3 regimen, this model and its Ctrough",
        "sibling agree closely on the ADR endpoint.",
        "Han 2024 does not tabulate the observed AUCtau distribution; the",
        "quartile boundaries used for Figure 4D are not printed."
      ),
      source_name        = "AUCtau (area under the steady-state curve)"
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
      "Every participant in this analysis received the same 150 mg dose, ",
      "so the exposure range driving the regression is generated purely ",
      "by between-subject pharmacokinetic variability (omega(CL) = 30.9% ",
      "in the companion PK model) rather than by dose randomisation. ",
      "That is the usual reason a landmark exposure-efficacy slope in a ",
      "single-dose-level trial has little power, and is worth keeping in ",
      "view alongside the paper's 'effect has reached a plateau' reading."
    )
  )

  ini({
    # ==================================================================
    # Han 2024 Table 4, lower block "Exposure (AUCtau)-virological
    # response model". The model is a univariate binomial logistic
    # regression of the week-48 virologic-suppression indicator on the
    # individual steady-state dosing-interval AUC:
    #
    #   logit(p) = beta0 + beta1 * AUCtau
    #
    # Table 4 prints an Estimate, an SD and a P-value per row and nothing
    # else -- no odds ratio, no confidence interval -- so there is no
    # redundant column against which to cross-check each coefficient.
    # Each value below is the Table 4 Estimate column, and the SDs are
    # recorded in the comments for anyone who needs the uncertainty.
    #
    # The slope's P-value is quoted twice and consistently: Table 4 gives
    # 0.222, and the Results text gives "(slope only: P = 0.220;
    # P = 0.222, respectively)", the second of which is this model.
    #
    # The covariate is NOT centred and NOT scaled, so the intercept is
    # the logit at AUCtau = 0 -- an extrapolated anchor outside the
    # observed exposure range, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- 1.44; label("Logit of the probability of HIV-RNA < 50 copies/mL at week 48 at AUC_ANV = 0 ng*h/mL (unitless logit)")  # Han 2024 Table 4, beta 0 (Intercept) = 1.44, SD 0.71, P = 0.042

    # ----- Exposure effect on the logit -----
    # Printed as "0.12x10 -4" in Table 4, i.e. 0.12e-4 = 1.2e-5. Note the
    # SD (1.02e-4) is 8.5-fold the estimate, which is what a P-value of
    # 0.222 looks like: the slope is indistinguishable from zero.
    e_auc_anv_logit <- 0.000012; label("Log-odds of HIV-RNA < 50 copies/mL at week 48 per 1 ng*h/mL increase in steady-state ainuovirine AUCtau (unitless logit)")  # Han 2024 Table 4, beta 1 (Slope) = 0.12x10^-4, SD 1.02x10^-4, P = 0.222 -- NOT significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and no random effects are estimated. The
    # tiny fixed additive residual below exists only so rxode2 has an
    # error model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_hivrna_lt50 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_hivrna_lt50 <- logit_ref + e_auc_anv_logit * AUC_ANV

    prob_hivrna_lt50 <- expit(logit_hivrna_lt50)

    prob_hivrna_lt50 ~ add(addSd_prob_hivrna_lt50)
  })
}
