Qin_2025_ropeginterferon_chr_week52 <- function() {
  description <- paste0(
    "Emax binomial logistic-regression exposure-efficacy model for ",
    "complete hematologic response (CHR) at WEEK 52 of subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 73 Chinese and Japanese ",
    "patients with polycythaemia vera (Qin 2025, phase II studies ",
    "A19-201 and A20-202). The probability of CHR is ",
    "expit(-7 + 8.397 * CAV / (1.98 + CAV)), where CAV is the ",
    "individual average total serum ropeg concentration over weeks ",
    "0-52 in ng/mL. THE BASELINE LOGIT IS FIXED AT -7 AND THE EC50 OF ",
    "1.98 ng/mL IS NOT SIGNIFICANT (p = 0.4231, standard error larger ",
    "than the estimate); the authors state plainly that 'confidence in ",
    "the Emax model was limited'. Because EC50 sits far below the ",
    "observed exposure range of roughly 10-70 ng/mL, the fitted curve ",
    "is effectively saturated everywhere the data live, which is how ",
    "this model encodes the FLAT week-52 exposure-response the ",
    "exploratory analysis found (p = 0.41). Treat it as a description ",
    "of an absent relationship, not as a usable potency estimate. ",
    "Contrast the week-24 companion ",
    "Qin_2025_ropeginterferon_chr_week24, which is linear and ",
    "significant. Companion models in the Qin_2025_ropeginterferon_* ",
    "family."
  )
  reference <- paste(
    "Qin A, Shimoda K, Suo S, Fu R, Kirito K, Wu D, Liao J, Chen H, Wu L,",
    "Su X, Gao Y, Sato T, Li Y, Zhang J, Shen W, Wang W, Zhang L, Jin J,",
    "Komatsu N.",
    "Population pharmacokinetics-pharmacodynamics and exposure-response of",
    "ropeginterferon alfa-2b in Chinese and Japanese patients with",
    "polycythemia vera.",
    "Pharmacol Res Perspect. 2025;13(3):e70109.",
    "doi:10.1002/prp2.70109.",
    sep = " "
  )
  vignette <- "Qin_2025_ropeginterferon"
  units <- list(
    time          = "n/a (static week-52 landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "prob_chr (probability of complete hematologic response at week 52, 0-1; also logit_chr)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over weeks 0 to 52.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg. Averaging window is",
        "weeks 0-52, i.e. from the FIRST DOSE to the landmark, so it",
        "includes the titration period and is NOT a steady-state",
        "interval average (Qin 2025 Methods 2.4.6.1). This is a",
        "different column from the week-24 companion model's CAV and the",
        "two must not be interchanged; the week-52 window dilutes the",
        "between-arm exposure difference, which Qin 2025 Results 3.3",
        "notes explicitly: 'At Week 52, the difference in exposure",
        "between the two titration regimens was less than that at Week",
        "24'. Derived by simulation from the actual dosing records and",
        "the individual post hoc (empirical Bayes) PK parameters of",
        "modellib('Qin_2025_ropeginterferon'). NOT centred and NOT",
        "scaled. Observed range about 10-70 ng/mL (Qin 2025 Figure 3D),",
        "with study medians near 28 ng/mL (A19-201) and 42 ng/mL",
        "(A20-202). Note that the whole observed range lies far ABOVE",
        "the fitted EC50 of 1.98 ng/mL, so the Emax term is saturated",
        "throughout and the model predicts a nearly constant",
        "probability."
      ),
      source_name        = "Cavg,0-52W (average concentration of participants from 0 to 52 weeks)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by forward inclusion at p < 0.05 and not retained (Qin 2025 Results 3.3: 'After covariate screening, no significant covariates were included in either model')."
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and not retained for this endpoint."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 73L,
    n_studies      = 2L,
    n_observations = "73 evaluable binary CHR records, one per patient (Qin 2025 Results 3.3: 'In two Phase II studies, 77 and 73 patients underwent CHR assessments at Weeks 24 and 52, respectively')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera; A20-202 enrolled patients resistant to or intolerant of hydroxyurea",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Complete hematologic response: hematocrit < 45% without phlebotomy in the previous 3 months, AND white blood cell count < 10 x 10^9/L, AND platelet count <= 400 x 10^9/L (Qin 2025 Methods 2.3)",
    notes          = paste0(
      "The Emax form was chosen because the exploratory analysis found ",
      "no monotone trend at this landmark: Qin 2025 Results 3.3 reports ",
      "'a flat trend between exposure and CHR rate was observed ",
      "(p = 0.41), and the probability plot showed identical CHR rates ",
      "in the second to fourth exposure quantiles and a lower CHR rate ",
      "in the first quantile.' The Discussion interprets this as ",
      "saturation rather than absence of effect: 'the CHR, which ",
      "consists of multiple blood parameter components, and ropeg ",
      "exposure had already reached a level close to maximum after 52 ",
      "weeks of treatment.' A user who wants a defensible ropeg ",
      "exposure-efficacy relationship should use the week-24 companion ",
      "model instead."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Table 4, row block "Emax logistic regression of CHR at
    # Week 52". Fitted in R 4.2.2 (Methods 2.4.8).
    #
    # Model form is Qin 2025 Equation (4):
    #
    #   logit(P_i,CHR) = E0 + Emax*Exposure_i/(EC50 + Exposure_i)
    #                       + beta^T * X_i
    #
    # with the beta^T * X_i covariate block empty because no covariate
    # survived forward inclusion at p < 0.05.
    #
    # HEALTH WARNING ON THIS PARAMETERISATION. The EC50 standard error
    # (2.472) EXCEEDS its point estimate (1.98) and p = 0.4231, so the
    # potency is not identified. Qin 2025 says so directly: "The
    # p-value for EC50 was not significant; therefore, confidence in
    # the Emax model was limited." The three parameters are also
    # strongly aliased -- with E0 held at -7 and EC50 far below the data
    # range, Emax simply absorbs whatever constant logit reproduces the
    # observed week-52 CHR rate. Reproducing the published curve
    # requires all three values together; none is interpretable alone.
    # ==================================================================

    # ----- Baseline logit (FIXED) -----
    # Table 4 prints "- 7 (FIX)" with no standard error and no p-value.
    # expit(-7) = 0.00091, i.e. the model asserts an essentially zero
    # CHR probability at zero exposure. This is a structural anchor
    # chosen to make the Emax term identifiable, not an estimate.
    logit_ref <- fixed(-7) ; label("Logit of the probability of complete hematologic response at week 52 at CAV = 0 ng/mL (unitless logit)")  # Qin 2025 Table 4: baseline logit E0 = -7 (FIX); no standard error or p-value is printed

    # ----- Saturable exposure effect on the logit -----
    lemax <- log(8.397) ; label("Maximum increment in the week-52 complete-hematologic-response logit attributable to ropeg exposure (unitless logit)")  # Qin 2025 Table 4: maximum effect Emax = 8.397, standard error 0.665, p < 0.0001
    lec50 <- log(1.98)  ; label("Average total serum ropeg concentration over weeks 0-52 producing half of Emax on the logit (EC50, ng/mL)")            # Qin 2025 Table 4: EC50 = 1.98, standard error 2.472, p = 0.4231 -- NOT significant, and the standard error exceeds the estimate

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and the analysis estimates no random
    # effects. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_chr <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # ----- Linear predictor, Qin 2025 Equation (4) -----
    logit_chr <- logit_ref + emax * CAV / (ec50 + CAV)

    prob_chr <- expit(logit_chr)

    # ----- Observation -----
    prob_chr ~ add(addSd_prob_chr)
  })
}
