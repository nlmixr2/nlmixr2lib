Riggs_2012_albinterferon_svr_gt1 <- function() {
  description <- "Logistic-regression exposure-response model for sustained virologic response (SVR) in the HCV GENOTYPE 1 subpopulation treated with subcutaneous albinterferon alfa-2b plus daily oral ribavirin, from Riggs 2012. SVR (undetectable HCV RNA below 15 IU/mL at 24 weeks post-treatment) is modelled on the logit scale as a linear function of the individual average steady-state albIFN concentration (Cavg, a post-hoc Bayesian exposure metric from the companion population PK model) plus reduced treatment duration, baseline HCV-RNA category, Asian and Black race, weight at or above 75 kg, and whether the dose was reduced. Fitted with glm in R, not in NONMEM, and separately from the genotype 2/3 subpopulation. SVR is only weakly related to exposure across the studied range; baseline HCV RNA, treatment duration, and Black race dominate. IMPORTANT: the source tabulates fitted PROBABILITIES for one-covariate-at-a-time perturbations of a defined reference individual rather than printing the logit coefficients, so every coefficient here is recovered by exact algebraic inversion of Table II - see the model file's SOURCING NOTE and the vignette. Sister model files from the same paper: modellib('Riggs_2012_albinterferon') for the population PK model and modellib('Riggs_2012_albinterferon_svr_gt23') for the genotype 2/3 stratum."
  reference <- paste(
    "Riggs MM, Bergsma TT, Rogers JA, Gastonguay MR, Subramanian GM, Chen C,",
    "Devalaraja M, Corey AE, Sun H, Yu J, Stein DS.",
    "Population pharmacokinetics and exposure-response of albinterferon alfa-2b.",
    "J Clin Pharmacol. 2012;52(4):475-486.",
    "doi:10.1177/0091270011399576.",
    "SVR exposure-response methods are in Methods, 'Data Analysis' (p. 478);",
    "the fitted probabilities this model is recovered from are in Table II",
    "(p. 482) and are drawn in Figure 2 (p. 482).",
    "Sister model files from the same paper:",
    "modellib('Riggs_2012_albinterferon'),",
    "modellib('Riggs_2012_albinterferon_svr_gt23').",
    sep = " "
  )
  vignette <- "Riggs_2012_albinterferon"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # ------------------------------------------------------------------
  # SOURCING NOTE - how the coefficients in ini() were obtained.
  #
  # Riggs 2012 fits the SVR exposure-response with glm (logistic
  # regression) in R and reports the result as Table II: a fitted SVR
  # PROBABILITY, with a 95% confidence interval, for a reference
  # individual ("None") and for ten one-covariate-at-a-time perturbations
  # of it. The logit-scale coefficients themselves are never printed.
  #
  # Because every Table II row perturbs exactly ONE covariate away from a
  # single common reference (footnote a, reproduced verbatim in the
  # population$notes field below), each coefficient is recovered exactly
  # as the difference of logits:
  #
  #     b_<covariate> = logit(p_row) - logit(p_reference)
  #
  # and the intercept is logit(p_reference) itself. The exposure slope
  # comes from the two Cavg rows:
  #
  #     b_cavg = (logit(0.641) - logit(0.596)) / (90 - 60)  ng/mL^-1
  #
  # This inversion is SELF-VALIDATING. The intercept implied by the two
  # Cavg rows, extrapolated to the reference exposure of 75 ng/mL, is
  # 0.48427; the intercept read directly off the "None" row is 0.48531.
  # The two agree to 0.001 on the logit scale - i.e. to 0.6188 versus
  # 0.6190 in probability, which is the rounding granularity of the
  # three-significant-figure percentages Table II prints. That agreement
  # also confirms the two structural facts the inversion depends on:
  # that Cavg enters LINEARLY on the logit scale (rather than through an
  # Emax-type form, which Methods p. 478 allows a priori), and that all
  # Table II rows really are single-covariate perturbations of one
  # common reference. The "None"-row value is used for the intercept
  # because it is the direct reference estimate.
  #
  # Round-tripping the fitted model reproduces all three Cavg rows and
  # all seven covariate rows of the genotype 1 column of Table II to the
  # printed three significant figures; the vignette asserts this.
  #
  # Precision: Table II prints probabilities to three significant
  # figures, so each recovered coefficient carries roughly +/-0.005 on
  # the logit scale. The confidence intervals in Table II are NOT
  # invertible in the same way (they are intervals on a fitted
  # probability, not on a coefficient), so no standard errors are
  # carried here.
  # ------------------------------------------------------------------

  covariateData <- list(
    CSS_ALBIFN = list(
      description        = "Individual average steady-state serum albinterferon alfa-2b concentration (Cavg).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A post-hoc exposure metric, not an observation: Riggs 2012 obtained each patient's PK parameters from the companion population PK model by Bayesian maximum a posteriori estimation, then computed Cavg as the cumulative AUC over the patient's whole dosing history divided by the treatment duration - so it reflects dose adjustments, unlike Cmax, which reflects only the randomized first dose (Methods p. 478). Generate it for simulation with modellib('Riggs_2012_albinterferon'). The reference individual's value is 75 ng/mL. Observed medians of the individual estimates were 62.7 ng/mL at 900 ug every 2 weeks and 83.0 ng/mL at 1200 ug every 2 weeks; the safety analysis binned the whole cohort at 22.4 / 54.4 / 68.4 / 85.9 / 176 ng/mL (Figure 3 caption). Cmax was deliberately NOT used for the SVR analysis because albIFN is an indirect antiviral acting through receptor binding over time (Methods p. 478).",
      source_name        = "Cavg"
    ),
    T_TRT = list(
      description        = "Total duration of albinterferon treatment received.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only DICHOTOMIZED. The full genotype 1 course is 48 weeks and Table II footnote b defines the model's 'reduced treatment duration' level as 24 weeks, so model() forms trtReduced <- (T_TRT < 48): any course short of the full 48 weeks is treated as reduced. Riggs 2012 included this term specifically to absorb patients who stopped early for lack of early virologic response or for an adverse event (Discussion p. 485), and notes that their Cavg estimates were consistent with everyone else's and that Cavg had approached steady state before the week-12 / week-24 futility assessments, so duration and exposure are not confounded in this model.",
      source_name        = "treatment duration"
    ),
    HCV_VLOAD = list(
      description        = "Baseline hepatitis C virus RNA concentration in serum.",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a THREE-level categorization, not continuously: below 400,000 IU/mL; the 400,000-to-800,000 band, which is the reference (Table II footnote a); and 800,000 IU/mL or above. model() forms the two non-reference indicators inline. Note the popPK sister model dichotomizes the same column at a single 800,000 IU/mL cut-point instead. Baseline HCV RNA below 400,000 IU/mL is the strongest single predictor in this model, raising the reference SVR probability from 61.9% to 84.6%.",
      source_name        = "baseline HCV-RNA category"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator: 1 = Asian, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Other race, pooled; see notes)",
      notes              = "The SVR reference category is 'white or other race' POOLED (Table II footnote a and the Figure 2 caption), which is NOT the same reference as the companion popPK model, where White alone is the reference and RACE_OTHER carries its own coefficient. There is therefore deliberately no RACE_OTHER term in this model.",
      source_name        = "Asian race"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator: 1 = Black, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Other race, pooled)",
      notes              = "Lowers the reference SVR probability from 61.9% to 42.4%, the largest race effect in the model. Riggs 2012 emphasizes that this lower response rate in Black patients with genotype 1 appeared INDEPENDENT of any difference in albIFN exposure, and is consistent with the phase 3 reports and with prior literature on ethnic differences in interferon response (Results p. 482).",
      source_name        = "Black race"
    ),
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only DICHOTOMIZED at 75 kg: model() forms wtGe75 <- (WT >= 75), matching Table II's 'Weight >= 75 kg' row against a reference of below 75 kg (footnote a). The companion popPK model uses the same column continuously, as a power function of WT / 75 kg. The two roles are consistent - 75 kg is the reference weight throughout the paper - but note that the weight effect on SVR here (61.9% to 59.8%) is small even though weight is the one covariate with a clinically meaningful effect on exposure.",
      source_name        = "weight category"
    ),
    DOSE_REDUCED = list(
      description        = "Dose-reduction indicator: 1 = the patient's albinterferon dose was reduced at some point during treatment, 0 = no dose reduction.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no dose reduction; the source's reference individual)",
      notes              = "Time-fixed per subject (it summarizes the whole treatment course). Dose reductions arose both from the protocol amendments that moved every patient receiving 1200 ug every 2 weeks down to 900 ug after the data-monitoring committee flagged serious pulmonary adverse events, and from per-patient reductions mandated for neutropenia, thrombocytopenia, or anemia. The estimated effect is essentially null (61.9% to 62.9%), which Riggs 2012 presents as expected: without a marked exposure-response relationship, shifts in exposure caused by dose reduction should not move SVR (Discussion p. 485). Distinct from the covariate's effect on exposure, which is already absorbed by CSS_ALBIFN because Cavg is computed over the actual dosing history.",
      source_name        = "dose was reduced"
    )
  )

  # Covariates that Riggs 2012 names among the factors considered for the
  # SVR analysis but that do NOT appear as a row of Table II, as a series
  # in Figure 2, or in the reference-individual definition of Table II
  # footnote a. Documentation only; checkModelConventions() does not
  # require these to be referenced in model(). See the vignette's
  # Assumptions and deviations section.
  covariatesDataExcluded <- list(
    FIBROSIS_SCORE = list(
      description        = "Baseline liver fibrosis stage (METAVIR-style F0-F4 scoring; the paper contrasts F4 against below F4).",
      units              = "(stage)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "SCREENED, and repeatedly described as an important SVR predictor - the Abstract lists 'fibrosis score' among the important SVR predictors, and Results p. 482 states that a fibrosis score below 4 was a better predictor of SVR than albIFN Cavg - but NO coefficient is recoverable. Table II has no fibrosis row, Figure 2 draws no fibrosis series, and the reference-individual definition of Table II footnote a, which pins a level for every one of the six covariates that DO appear, is silent on fibrosis. Consequently the intercept recovered here is implicitly conditioned on whatever fibrosis level the reference individual carried, and this model cannot vary fibrosis. Every other coefficient is unaffected, because each is a difference of logits taken at a fixed fibrosis level. Fibrosis IS used elsewhere in the paper, as a stratifying variable (F4 versus below F4) for the tabular safety exposure-response summaries (Methods p. 478).",
      source_name        = "fibrosis score"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Named in the Discussion (p. 485) among the factors that, along with baseline viral load and genotype, predicted SVR, but there is no age row in Table II, no age series in Figure 2, and no age level in the reference-individual definition. Not recoverable; treated the same way as FIBROSIS_SCORE above.",
      source_name        = "age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "18-79 years (population PK analysis set)",
    weight_range   = "38-166 kg (population PK analysis set)",
    disease_state  = "Chronic hepatitis C virus GENOTYPE 1 infection, interferon-treatment-naive, with compensated liver disease.",
    dose_range     = "900 or 1200 ug albinterferon alfa-2b subcutaneously once every 2 weeks for 48 weeks, plus daily oral ribavirin. The 1200 ug arm was reduced to 900 ug during the phase 3 trial because of serious pulmonary adverse events.",
    endpoint       = "Sustained virologic response: undetectable HCV RNA (below 15 IU/mL) at 24 weeks after the end of treatment. Intent-to-treat SVR rates in the genotype 1 phase 3 study were 48.2% and 47.3% in the 900 and 1200 ug every-2-weeks groups.",
    notes          = "The exposure-response analysis set is the 1869 interferon-naive patients that remain after excluding phase 2 trial NCT00097435 (which enrolled prior interferon non-responders) from the 1984-patient population PK data set; the paper does not report how those 1869 split between the genotype 1 and genotype 2/3 analyses, so n_subjects is left NA rather than guessed. Reference individual for the Table II estimates, verbatim from footnote a and the Figure 2 caption: treatment duration 48 weeks; Cavg 75 ng/mL; baseline HCV RNA at or above 400,000 and below 800,000 IU/mL; white or other race; weight below 75 kg; and no dose reduction."
  )

  ini({
    # Every value below is recovered from the GENOTYPE 1 column of
    # Riggs 2012 Table II by the logit inversion described in the
    # SOURCING NOTE above. No coefficient is printed in the paper.
    logitsvr <- 0.485306 ; label("Logit of the SVR probability for the reference individual (unitless)")  # Table II, genotype 1, row "None": 61.9% (95% CI 51.7-71.1) -> qlogis(0.619)

    e_css_albifn_logitsvr    <- 0.0063627 ; label("Effect of average steady-state albIFN concentration on logit(SVR), per ng/mL above 75 ng/mL")  # Table II, genotype 1, rows "Cavg 60 ng/mL" 59.6% and "Cavg 90 ng/mL" 64.1%: (qlogis(0.641) - qlogis(0.596)) / 30
    e_trtdur_reduced_logitsvr <- -2.235680 ; label("Effect of reduced (24 week rather than 48 week) treatment duration on logit(SVR)")            # Table II, genotype 1, row "Reduced treatment duration" (footnote b: 24 weeks for genotype 1): 14.8% -> qlogis(0.148) - qlogis(0.619)
    e_hcv_vload_lt400k_logitsvr <- 1.218261 ; label("Effect of baseline HCV RNA below 400,000 IU/mL on logit(SVR)")                                # Table II, genotype 1, row "HCV RNA < 400 000 IU/mL": 84.6% -> qlogis(0.846) - qlogis(0.619)
    e_hcv_vload_ge800k_logitsvr <- 0.046911 ; label("Effect of baseline HCV RNA at or above 800,000 IU/mL on logit(SVR)")                          # Table II, genotype 1, row "HCV RNA >= 800 000 IU/mL": 63.0% -> qlogis(0.630) - qlogis(0.619)
    e_race_asian_logitsvr    <- 0.555472 ; label("Effect of Asian race on logit(SVR)")                                                             # Table II, genotype 1, row "Asian race": 73.9% -> qlogis(0.739) - qlogis(0.619)
    e_race_black_logitsvr    <- -0.791680 ; label("Effect of Black race on logit(SVR)")                                                            # Table II, genotype 1, row "Black race": 42.4% -> qlogis(0.424) - qlogis(0.619)
    e_wt_ge75_logitsvr       <- -0.088167 ; label("Effect of body weight at or above 75 kg on logit(SVR)")                                         # Table II, genotype 1, row "Weight >= 75 kg": 59.8% -> qlogis(0.598) - qlogis(0.619)
    e_dose_reduced_logitsvr  <- 0.042623 ; label("Effect of having had a dose reduction on logit(SVR)")                                            # Table II, genotype 1, row "Dose was reduced": 62.9% -> qlogis(0.629) - qlogis(0.619)

    # No random effects and no residual-error parameter. The source fits a
    # fixed-effects logistic regression with glm, one row per patient, so
    # all the stochasticity lives in the Bernoulli endpoint itself.
  })

  model({
    # ----------------------------------------------------------------
    # 1. Derived covariate terms. Every threshold is applied here rather
    #    than baked into the data, so the cut-points stay auditable.
    # ----------------------------------------------------------------
    # Full genotype 1 course is 48 weeks; Table II footnote b sets the
    # "reduced" level at 24 weeks.
    trtReduced <- (T_TRT < 48)

    # Three-level baseline HCV RNA; the 400,000-800,000 IU/mL band is the
    # reference and carries no term.
    hcvLo <- (HCV_VLOAD < 400000)
    hcvHi <- (HCV_VLOAD >= 800000)

    wtGe75 <- (WT >= 75)

    # ----------------------------------------------------------------
    # 2. Linear predictor on the logit scale. CSS_ALBIFN is centred at
    #    the reference individual's 75 ng/mL so that logitsvr is exactly
    #    the reference log-odds.
    # ----------------------------------------------------------------
    lpsvr <-
      logitsvr +
      e_css_albifn_logitsvr * (CSS_ALBIFN - 75) +
      e_trtdur_reduced_logitsvr * trtReduced +
      e_hcv_vload_lt400k_logitsvr * hcvLo +
      e_hcv_vload_ge800k_logitsvr * hcvHi +
      e_race_asian_logitsvr * RACE_ASIAN +
      e_race_black_logitsvr * RACE_BLACK +
      e_wt_ge75_logitsvr * wtGe75 +
      e_dose_reduced_logitsvr * DOSE_REDUCED

    # ----------------------------------------------------------------
    # 3. Probability of SVR, and the Bernoulli endpoint. psvr is exposed
    #    as an output column so the fitted probabilities can be compared
    #    against Table II directly, without simulating draws.
    # ----------------------------------------------------------------
    psvr <- expit(lpsvr)

    svr ~ dbinom(1, psvr)
  })
}
