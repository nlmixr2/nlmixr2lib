Han_2024_ainuovirine_adr_auctau <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model relating the steady-state ",
    "dosing-interval exposure (AUCtau) of the non-nucleoside reverse ",
    "transcriptase inhibitor ainuovirine (ANV) to the incidence of adverse ",
    "drug reactions (ADR, composite across all reported reactions) in ",
    "antiretroviral-therapy-naive adults living with HIV-1 (Han 2024, the ",
    "phase 3 trial ADYY-ACC007-301, ChiCTR1800019041, 150 mg once daily at ",
    "bedtime for 48 weeks). The probability of an ADR is ",
    "expit(-0.27 + 0.00014 * AUC_ANV), where AUC_ANV is the individual ",
    "steady-state AUC over the 24 h dosing interval in ng*h/mL derived by ",
    "Bayesian post-hoc estimation from the companion population PK model ",
    "(modellib('Han_2024_ainuovirine')). This slope IS significant ",
    "(P = 0.015) and is one of the paper's two positive findings: the ",
    "probability of an adverse reaction rises with ANV exposure while the ",
    "companion virologic-response models are flat, which is the asymmetry ",
    "behind the paper's conclusion that 'optimization of ANV dose may be ",
    "warranted in clinical practice'. The intercept is not distinguishable ",
    "from zero (P = 0.500), so read the model as a slope result rather ",
    "than as a calibrated absolute-risk predictor. This is the ",
    "AUCtau-parameterised sibling of Han_2024_ainuovirine_adr_ctrough: Han ",
    "2024 fitted the two exposure metrics as separate univariate ",
    "regressions and reported both, so they are alternatives rather than a ",
    "joint model, and must not be applied together. There is no PK layer ",
    "and no ODE: the exposure metric is supplied as a data column. Three ",
    "companion exposure-response models in the Han_2024_ainuovirine_* ",
    "family."
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
    time          = "n/a (static landmark exposure-safety regression at week 48; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_ANV covariate column)",
    concentration = "prob_adr (probability of an adverse drug reaction, 0-1; also logit_adr)"
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
        "regression indicated a linear effect both in Ctrough with the",
        "incidence of ADR and AUCtau with the incidence of ADR'. The",
        "printed coefficient is therefore per ng*h/mL.",
        "UNITS ARE LOAD-BEARING and Han 2024 never states them on Table 5.",
        "ng*h/mL is established from the paper's own unit system: Figure 1D",
        "plots dose-normalised AUCinf on an axis labelled",
        "'h*ng*mL^-1*mg^-1', so the un-normalised AUC is in ng*h/mL, and",
        "every ANV concentration in the paper is in ng/mL. Two cross-checks",
        "agree. (1) Scale: the ratio of this endpoint's Ctrough slope to",
        "its AUCtau slope is 0.40e-2 / 0.14e-3 = 28.6 h, the right order",
        "for an AUCtau/Ctrough ratio over a 24 h interval, which it would",
        "not be if the two metrics differed by a factor of 1000.",
        "(2) Agreement: at the exposure the companion PK model predicts for",
        "the phase 3 regimen, this model and its Ctrough sibling agree",
        "closely on this endpoint.",
        "Han 2024 does not tabulate the observed AUCtau distribution; the",
        "quartile boundaries used for Figure 5D are not printed."
      ),
      source_name        = "AUCtau (area under the steady-state curve)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 309L,
    n_studies      = 1L,
    n_observations = "One binary adverse-drug-reaction record per subject. NOTE: Han 2024 never states the exposure-response analysis-set size. 309 is the number of phase 3 ADYY-ACC007-301 participants in the population PK dataset (Table 1, Table 2). Han 2024 also never reports the overall ADR incidence, the reaction types, or their severity grading, so this model cannot be calibrated against a published event rate",
    age_range      = "mean 31.03 years (SD 9.95) in the phase 3 trial (Han 2024 Table 2)",
    weight_range   = "mean 67.06 kg (SD 11.42) in the phase 3 trial (Han 2024 Table 2)",
    sex_female_pct = 5.83,
    race_ethnicity = "Not reported; single Chinese centre (Beijing Youan Hospital, Capital Medical University)",
    disease_state  = "Antiretroviral-therapy-naive people living with HIV-1, randomised to an ainuovirine-based regimen (the comparator arm received an efavirenz-based regimen and is not part of this exposure-response analysis)",
    dose_range     = "Ainuovirine 150 mg orally once daily at bedtime on an empty stomach for 48 weeks",
    regions        = "Single centre, Beijing, China",
    notes          = paste0(
      "The endpoint is a COMPOSITE 'incidence of adverse reactions' with ",
      "no preferred-term breakdown and no severity threshold reported, ",
      "so it is not comparable to a CTCAE-graded endpoint. Because every ",
      "participant received the same 150 mg dose, the exposure range ",
      "driving the regression comes purely from between-subject ",
      "pharmacokinetic variability (omega(CL) = 30.9% in the companion PK ",
      "model). Han 2024 separately reports that AUCtau differed between ",
      "participants with and without adverse events (P = 0.0141), a ",
      "two-group comparison distinct from -- and consistent with -- the ",
      "regression slope below."
    )
  )

  ini({
    # ==================================================================
    # Han 2024 Table 5, lower block "Exposure (AUCtau)-ADR model". The
    # model is a univariate binomial logistic regression of the
    # adverse-drug-reaction indicator on the individual steady-state
    # dosing-interval AUC:
    #
    #   logit(p) = beta0 + beta1 * AUCtau
    #
    # Table 5 prints an Estimate, an SD and a P-value per row and nothing
    # else -- no odds ratio, no confidence interval -- so there is no
    # redundant column against which to cross-check each coefficient.
    # Each value below is the Table 5 Estimate column, and the SDs are
    # recorded in the comments for anyone who needs the uncertainty.
    #
    # TWO DIFFERENT P-VALUES APPEAR IN THE PAPER FOR THIS ENDPOINT AND
    # THEY ARE NOT IN CONFLICT -- they test different things:
    #   * 0.0141 (Abstract, Results) is the P-value of the two-group
    #     comparison of AUCtau between participants who did and did not
    #     report an adverse event (Figure 5B).
    #   * 0.015 (Table 5) / 0.0153 (Results text, "P = 0.0192;
    #     P = 0.0153, respectively") is the P-value of the regression
    #     SLOPE below.
    # Only the second belongs to this model.
    #
    # The covariate is NOT centred and NOT scaled, so the intercept is
    # the logit at AUCtau = 0 -- an extrapolated anchor outside the
    # observed exposure range, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    # SD 0.40, P = 0.500: indistinguishable from zero.
    logit_ref <- -0.27; label("Logit of the probability of an adverse drug reaction at AUC_ANV = 0 ng*h/mL (unitless logit)")  # Han 2024 Table 5, beta 0 (Intercept) = -0.27, SD 0.40, P = 0.500

    # ----- Exposure effect on the logit -----
    # Printed as "0.14x10 -3" in Table 5, i.e. 0.14e-3 = 0.00014. The
    # odds of an ADR are multiplied by exp(0.00014) per ng*h/mL, i.e. by
    # exp(0.14) = 1.15 per 1000 ng*h/mL.
    e_auc_anv_logit <- 0.00014; label("Log-odds of an adverse drug reaction per 1 ng*h/mL increase in steady-state ainuovirine AUCtau (unitless logit)")  # Han 2024 Table 5, beta 1 (Slope) = 0.14x10^-3, SD 5.69x10^-5, P = 0.015 -- significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and no random effects are estimated. The
    # tiny fixed additive residual below exists only so rxode2 has an
    # error model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_adr <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_adr <- logit_ref + e_auc_anv_logit * AUC_ANV

    prob_adr <- expit(logit_adr)

    prob_adr ~ add(addSd_prob_adr)
  })
}
