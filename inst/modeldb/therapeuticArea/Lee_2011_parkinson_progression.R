Lee_2011_parkinson_progression <- function() {
  description <- paste0(
    "Disease-progression model for the change-from-baseline in total ",
    "Unified Parkinson's Disease Rating Scale (UPDRS) score over study ",
    "time in early Parkinson's disease, fit by Lee and Gobburu (2011) ",
    "as the worked example of a Bayesian disease-drug-trial modeling ",
    "methodology paper. The predicted change from baseline is the ",
    "linear-disease-progression slope minus an asymptotic short-term ",
    "symptomatic-effect dip: mu = slope * time - symeff * (1 - exp(-ke0 * time)). ",
    "Both slope and symptomatic effect have separate placebo (beta0, ",
    "gamma0) and active-drug (beta1, gamma1) terms switched on by a ",
    "binary ON_TREATMENT indicator (0 = placebo, 1 = active drug arm). ",
    "Between-subject variability is additive on the linear slope ",
    "(etaslope) and on the linear symptomatic effect (etasymeff), with ",
    "the two etas assumed independent in the source paper. Residual ",
    "error is additive on the UPDRS scale. No drug PK input -- the ",
    "drug effect is encoded purely as a treatment-arm indicator. ",
    "Time is in 4-week months (the units of the model's parameters as ",
    "fit; the source paper's Table II column header 'Delta UPDRS score/week' ",
    "appears to be a typo for '/month' -- see the validation vignette's ",
    "Assumptions and deviations section for the math that resolves this ",
    "from the paper's own published placebo prediction at week 26). ",
    "Default parameter values in ini() are the TEMPO study rasagiline ",
    "(1 or 2 mg/day pooled) arm, Table II non-informative-prior ",
    "Bayesian posterior means. The validation vignette additionally ",
    "documents the ELLDOPA study (carbidopa/levodopa) Table II fit and ",
    "the TEMPO Table III power-prior sensitivity at weight parameter ",
    "alpha_0 in {0.1, 0.5, 1.0}."
  )
  reference <- paste(
    "Lee JY, Gobburu JVS.",
    "Bayesian Quantitative Disease-Drug-Trial Models for Parkinson's",
    "Disease to Guide Early Drug Development.",
    "The AAPS Journal 2011;13(4):508-518.",
    "doi:10.1208/s12248-011-9293-6.",
    sep = " "
  )
  vignette <- "Lee_2011_parkinson_progression"

  # Paper-specific etas: the source paper places the random effects
  # b_1i and b_2i additively on the linear slope and symeff
  # parameters (b_1i ~ N(0, w^2_1) added directly to slope_i;
  # b_2i ~ N(0, w^2_2) added directly to symeff_i). Because slope and
  # symeff are not log-transformed primary parameters (beta0/beta1
  # and gamma0/gamma1 are signed linear-scale fixed effects), the
  # canonical eta + l<param> pairing does not apply; declare the
  # eta names as paper-specific so checkModelConventions() accepts
  # them. See references/parameter-names.md "Paper-specific etas".
  paper_specific_etas <- c("etaslope", "etasymeff")

  units <- list(
    time          = "month (4 weeks; see description and vignette Errata for the source-paper week-vs-month unit discrepancy)",
    dosing        = "(none; the drug effect is encoded via the ON_TREATMENT covariate, not via drug input)",
    concentration = "(change from baseline in total UPDRS score, unitless; observation deltaUPDRS)"
  )

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Binary treatment-arm indicator. 0 = placebo arm; 1 = active drug arm. Time-fixed per subject in the source Lee 2011 parallel-group analyses (the indicator records the subject's randomization assignment, not a per-record drug-presence flag). In the source Lee 2011 paper this is the generic Trt indicator that toggles both the drug-effect-on-disease-progression slope (beta1) and the drug-effect-on-short-term-symptomatic-effect magnitude (gamma1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm).",
      notes              = "The default parameter values in this file's ini() are from the TEMPO study (rasagiline 1 or 2 mg/day pooled into a single active arm) Table II non-informative-prior Bayesian posterior means, so ON_TREATMENT = 1 in this file predicts the rasagiline-active disease trajectory. For levodopa-arm simulation, refer to the ELLDOPA Table II fit documented in the validation vignette's Errata / Assumptions and deviations section; for the TEMPO power-prior sensitivity (alpha_0 in {0.1, 0.5, 1.0}) refer to the same section.",
      source_name        = "Trt"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "Mean age 64 years (TEMPO study; Lee 2011 Data section). Range not reported in Lee 2011.",
    age_median     = "mean 64 years (TEMPO; Lee 2011 Data section)",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = "More than 90% Caucasian (Lee 2011 Data section). Exact percentages not reported.",
    disease_state  = "Early Parkinson's disease (TEMPO study: TVP-1012 in Early Monotherapy for Parkinson's disease Outpatients trial; double-blinded, randomized, fixed-dose parallel group). The full source Lee 2011 paper additionally analyses the ELLDOPA study (carbidopa/levodopa in de novo PD outpatients); the ELLDOPA-specific parameter estimates are documented in the validation vignette but not in this file's default ini().",
    dose_range     = "n/a (this model has no drug input; the drug effect is captured via the ON_TREATMENT covariate). In the underlying TEMPO trial the active arm was rasagiline 1 or 2 mg/day, pooled into a single 'active drug' arm in the Lee 2011 analysis because the two dose levels overlapped in disease-progression profile (Lee 2011 Data section). The companion ELLDOPA trial used carbidopa/levodopa 12.5/50, 25/100, or 50/200 mg three times daily.",
    regions        = "(not reported in Lee 2011)",
    notes          = "Exact subject counts, demographic breakdowns, and per-arm enrollment numbers are not reported in the Lee 2011 methodology paper (it does not retabulate the TEMPO or ELLDOPA baseline tables). Study durations: 26 weeks (= 6.5 months in the model's 4-week-month time axis) for TEMPO and 24 weeks (= 6 months) for ELLDOPA (Lee 2011 Data section). Both studies enrolled predominantly male patients (Lee 2011 Data section: 'more male patients were enrolled than female patients in both studies'); exact percentages are not given. The Lee 2011 analysis combined the two rasagiline dose levels (1 mg/day and 2 mg/day) into a single 'active drug' arm because the two doses overlapped in disease-progression time profile."
  )

  ini({
    # Population-mean structural parameters from Lee 2011 Table II,
    # TEMPO study, non-informative-prior Bayesian posterior means.
    # Slope-domain parameters (beta0, beta1) carry units of UPDRS units
    # per month (where 1 month = 4 weeks per the math that reproduces
    # the paper's published placebo prediction of 3.62 Delta UPDRS at week
    # 26: 0.73 * 6.5 - 1.12 * (1 - exp(-1.46 * 6.5)) = 3.625).
    # Symeff-domain parameters (gamma0, gamma1) carry units of UPDRS
    # units (the asymptotic value of the (1 - exp(-ke0 * t)) symptomatic
    # dip component). ke0 is the first-order rate constant (1/month)
    # for the approach to the asymptotic symptomatic effect.
    #
    # NOTE on source-paper unit typo: Table II column headers state
    # 'Delta UPDRS score/week' for the slope-domain parameters and
    # 'Delta UPDRS/week' for the symeff-domain parameters. The math of the
    # paper's own published placebo prediction (3.62 at week 26)
    # demonstrates that the parameters are actually per 4-week month
    # rather than per week (a week-based interpretation would give
    # 0.73 * 26 - 1.12 = 17.86, not 3.62, and ELLDOPA's natural-progression
    # slope of 0.99 per week would imply an implausibly fast progression
    # of about 50 UPDRS units per year vs the literature value of about
    # 10-12 per year for early PD). The symeff-domain '/week' header is
    # additionally dimensionally suspect because symeff appears as
    # symeff * (1 - exp(-ke0 * t)) in equation (1) and (1 - exp(...))
    # is unitless, so symeff must be in UPDRS units, not UPDRS/week. The
    # validation vignette's Assumptions and deviations section reproduces
    # this calculation and exposes the discrepancy explicitly.
    #
    # Sign convention: beta0 is positive (placebo arm exhibits disease
    # progression -- UPDRS increases over time), and beta1 is negative
    # (active drug slows progression). gamma0 and gamma1 are positive
    # because the model subtracts (gamma0 + gamma1 * ON_TREATMENT) * (1 - exp(-ke0 * t))
    # from the linear slope, so a positive gamma encodes a symptomatic
    # benefit (the predicted Delta UPDRS dips below the disease-progression
    # line early on). Because the source paper observes a between-subject
    # random effect on the linear scale rather than the log scale, the
    # population means are not log-transformed -- beta0/beta1/gamma0/gamma1
    # are signed linear parameters, not log-transformed positive
    # parameters.

    beta0  <- 0.73       ; label("Placebo natural disease-progression slope (UPDRS units / month)")
    # Lee 2011 Table II TEMPO 'beta_0 (placebo effect on slope)' = 0.73 (posterior SD 0.09). Source column header '/week' is a typo for '/month'; see header note.

    beta1  <- -0.30      ; label("Active-drug additional effect on disease-progression slope (UPDRS units / month)")
    # Lee 2011 Table II TEMPO 'beta_1 (drug effect on slope)' = -0.30 (posterior SD 0.11). Negative -> active drug slows progression.

    gamma0 <- 1.12       ; label("Placebo asymptotic short-term symptomatic dip (UPDRS units)")
    # Lee 2011 Table II TEMPO 'gamma_0 (placebo on symptomatic effect)' = 1.12 (posterior SD 0.48). Positive -> placebo arm exhibits early-on symptomatic benefit that decays toward the long-run disease-progression slope.

    gamma1 <- 1.61       ; label("Active-drug additional asymptotic short-term symptomatic dip (UPDRS units)")
    # Lee 2011 Table II TEMPO 'gamma_1 (drug on symptomatic effect)' = 1.61 (posterior SD 0.60). Positive -> active drug adds an extra symptomatic benefit on top of placebo.

    lke0   <- log(1.46)  ; label("Log first-order rate constant for approach to asymptotic symptomatic effect (1/month)")
    # Lee 2011 Table II TEMPO 'Ke_0 (speed to reach max symptomatic effect)' = 1.46 (posterior SD 0.61). Source paper text: 'The parameter ke0 measures the speed with which the initial symptomatic effect is achieved, which was assumed to be the same for both placebo and drug arms.'

    # Additive between-subject variability on the linear slope and on
    # the linear symeff. Source Lee 2011 reports the variances directly
    # (paper Table II 'w^2_1' and 'w^2_2' columns are labelled
    # 'between-subject variability'). The two etas are assumed
    # independent in the source paper (Methods text: 'For simplicity,
    # b_1i and b_2i were assumed to be independent'). nlmixr2 enters
    # them as additive normal random effects on the linear slope and
    # symeff; the variance is on the same linear scale and in the same
    # units as the parameter it modifies (UPDRS units/month for slope,
    # UPDRS units for symeff). The Lee 2011 paper's printed Sigma_2x2
    # matrix (page 510) shows the variances in the off-diagonal
    # positions and zeros on the diagonal -- this appears to be a
    # transcription error; the only reading consistent with the
    # 'assumed independent' prose has the variances on the diagonal.
    etaslope  ~ 0.47   # Lee 2011 Table II TEMPO 'w^2_1 (between-subject variability in slope, variance)' = 0.47 (posterior SD on the variance 0.07).
    etasymeff ~ 18.61  # Lee 2011 Table II TEMPO 'w^2_2 (between-subject variability in symptomatic effect, variance)' = 18.61 (posterior SD on the variance 2.31).

    # Additive residual error on the UPDRS scale. Source Lee 2011
    # reports sigma^2 (paper Table II 'sigma^2 (residual error)' = 8.53
    # for the TEMPO non-informative fit), which is a variance; the
    # additive SD for nlmixr2 is the square root.
    addSd <- sqrt(8.53) ; label("Additive residual error standard deviation on UPDRS scale (UPDRS units)")
    # Lee 2011 Table II TEMPO 'sigma^2 (residual error, variance)' = 8.53 (posterior SD on the variance 0.31). addSd = sqrt(8.53) ~ 2.921.
  })

  model({
    # Individual-subject linear parameters. Slope and symeff each carry
    # an additive placebo intercept (beta0 / gamma0), an additive
    # active-drug shift (beta1 / gamma1) toggled by ON_TREATMENT, and
    # an additive per-subject random effect (etaslope / etasymeff) on
    # the linear scale.
    slope  <- beta0  + beta1  * ON_TREATMENT + etaslope
    symeff <- gamma0 + gamma1 * ON_TREATMENT + etasymeff
    ke0    <- exp(lke0)

    # Predicted change from baseline in total UPDRS score at study
    # time t (months, where 1 month = 4 weeks per the source paper's
    # implicit time axis). Lee 2011 Eq. (1):
    #   mu_it = slope_i * t - symeff_i * (1 - exp(-ke0 * t))
    # The symptomatic-effect term is SUBTRACTED -- it produces an
    # initial dip below the linear disease-progression baseline that
    # rises asymptotically back toward the slope-only trajectory as the
    # symptomatic benefit saturates. This matches Lee 2011 Figure 1's
    # 'symptomatic effect' curve shape (initial dip, then parallel to
    # the natural-progression line).
    deltaUPDRS <- slope * time - symeff * (1 - exp(-ke0 * time))

    deltaUPDRS ~ add(addSd)
  })
}
