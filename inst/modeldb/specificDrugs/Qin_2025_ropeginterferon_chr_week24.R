Qin_2025_ropeginterferon_chr_week24 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-efficacy model for complete ",
    "hematologic response (CHR) at WEEK 24 of subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 77 Chinese and Japanese ",
    "patients with polycythaemia vera (Qin 2025, phase II studies ",
    "A19-201 and A20-202). The probability of CHR is ",
    "expit(-1.2614 + 0.0411 * CAV), where CAV is the individual ",
    "average total serum ropeg concentration over weeks 0-24 in ng/mL, ",
    "derived from the empirical-Bayes parameters of the companion ",
    "population PK model Qin_2025_ropeginterferon. The exposure term ",
    "is SIGNIFICANT (p = 0.0284), which is the paper's headline ",
    "efficacy finding at this landmark and the quantitative basis for ",
    "preferring fast dose titration. No covariate survived screening. ",
    "There is no PK layer and no ODE. No between-subject random effect ",
    "and no residual error are estimated (Bernoulli likelihood). ",
    "Contrast the week-52 companion ",
    "Qin_2025_ropeginterferon_chr_week52, where the relationship has ",
    "flattened. Companion models in the Qin_2025_ropeginterferon_* ",
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
    time          = "n/a (static week-24 landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "prob_chr (probability of complete hematologic response at week 24, 0-1; also logit_chr)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over weeks 0 to 24.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg, matching the",
        "bioassay and the CS driver of the companion PK-PD models.",
        "Averaging window is weeks 0-24, i.e. from the FIRST DOSE to the",
        "landmark, not a single dosing interval at steady state (Qin",
        "2025 Methods 2.4.6.1: 'average concentration (Cavg,t) of",
        "patients from Weeks 0 to 24 and 0 to 52 was selected to",
        "approximate the actual average drug exposure of patients from",
        "the start of treatment to Weeks 24 and 52'). Because titration",
        "is still climbing over much of that window, this value is",
        "materially LOWER than the eventual steady-state Cav and the two",
        "must not be interchanged. Derived by simulation from the actual",
        "dosing records and the individual post hoc (empirical Bayes) PK",
        "parameters of modellib('Qin_2025_ropeginterferon'). NOT centred",
        "and NOT scaled, so the intercept is the logit at CAV = 0, an",
        "extrapolated anchor rather than a reference-patient",
        "probability. Observed range about 8-64 ng/mL (Qin 2025",
        "Figure 3A), with study medians near 20 ng/mL (A19-201, slow",
        "titration) and 37.5 ng/mL (A20-202, fast titration)."
      ),
      source_name        = "Cavg,0-24W (average concentration of participants from 0 to 24 weeks)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened by forward inclusion at p < 0.05 (Qin 2025 Methods",
        "2.4.6.2: 'Covariates tested were the same as those in the PopPK",
        "analysis') and NOT retained: 'After covariate screening, no",
        "significant covariates were included in either model.' Weight IS",
        "retained in the week-24 JAK2 V617F companion model",
        "(Qin_2025_ropeginterferon_jak2_week24), so this is an",
        "endpoint-specific null result rather than a reporting gap."
      )
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and not retained for this endpoint, despite being the single covariate retained on clearance in the companion population PK model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 2L,
    n_observations = "77 evaluable binary CHR records, one per patient (Qin 2025 Results 3.3: 'In two Phase II studies, 77 and 73 patients underwent CHR assessments at Weeks 24 and 52, respectively'); one of the 78 exposure-efficacy patients was not assessable at week 24",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline JAK2 V617F allele burden median 77.8% (A19-201) and 61.2% (A20-202)",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Complete hematologic response: hematocrit < 45% without phlebotomy in the previous 3 months, AND white blood cell count < 10 x 10^9/L, AND platelet count <= 400 x 10^9/L (Qin 2025 Methods 2.3)",
    notes          = paste0(
      "Model selection was evidence-led, not assumed. Qin 2025 ran a ",
      "graphical exploratory analysis first (Figure S5) and found a ",
      "significant exposure-CHR relationship at week 24 (p = 0.024 by ",
      "two-sided t-test on the CHR-versus-no-CHR exposure boxplots), ",
      "which is why a LINEAR logistic model was fitted here whereas an ",
      "Emax logistic was fitted at week 52. Observed CHR rates by ",
      "exposure quartile rise monotonically from about 0.35 in the ",
      "lowest quartile to about 0.68 in the highest (Figure 3A)."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Table 4, row block "Linear logistic regression of CHR at
    # Week 24". Fitted in R 4.2.2 (Methods 2.4.8), not in NONMEM: the
    # exposure-response layer is a landmark regression on the individual
    # post hoc exposures, not a mixed-effects model.
    #
    # Model form is Qin 2025 Equation (3):
    #
    #   logit(P_i,CHR) = log(P/(1-P)) = beta0 + beta1*Exposure_i
    #                                        + beta^T * X_i
    #
    # with the beta^T * X_i covariate block empty because no covariate
    # survived forward inclusion at p < 0.05.
    #
    # Table 4 prints an Estimate, a standard error and a p-value per row
    # and nothing else -- no odds ratio and no confidence interval -- so
    # there is no redundant printed column to cross-check the
    # coefficients against. Each value below is the Table 4 Estimate,
    # with the standard error recorded in the comment.
    #
    # The exposure regressor is NOT centred and NOT scaled.
    # ==================================================================

    # ----- Logit intercept -----
    # expit(-1.2614) = 0.221, the extrapolated CHR probability at zero
    # exposure. Marginally non-significant on its own (p = 0.0532),
    # which matters only for interpreting the anchor, not the slope.
    logit_ref <- -1.2614 ; label("Logit of the probability of complete hematologic response at week 24 at CAV = 0 ng/mL (unitless logit)")  # Qin 2025 Table 4: baseline logit beta0 = -1.2614, standard error 0.6525, p = 0.0532

    # ----- Exposure effect on the logit -----
    # exp(0.0411) = 1.042, so the odds of CHR at week 24 rise by 4.2%
    # per ng/mL and by exp(0.0411*20) = 2.27-fold across the roughly
    # 20 ng/mL gap between the slow- and fast-titration study medians.
    e_cav_logit <- 0.0411 ; label("Log-odds of complete hematologic response at week 24 per 1 ng/mL increase in average total serum ropeg concentration over weeks 0-24 (unitless logit)")  # Qin 2025 Table 4: exposure effect beta1 = 0.0411, standard error 0.0187, p = 0.0284 -- significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and the analysis estimates no random
    # effects. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity. See the vignette's
    # Assumptions and deviations.
    addSd_prob_chr <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (3) -----
    logit_chr <- logit_ref + e_cav_logit * CAV

    prob_chr <- expit(logit_chr)

    # ----- Observation -----
    prob_chr ~ add(addSd_prob_chr)
  })
}
