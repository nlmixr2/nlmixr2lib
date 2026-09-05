Chen_2021_lorlatinib_weight_gain <- function() {
  description <- paste0(
    "Binomial logistic-regression safety model for CTCAE grade >= 2 ",
    "weight gain in adults with advanced ALK-positive or ROS1-positive ",
    "non-small cell lung cancer treated with the third-generation ",
    "ALK/ROS1 tyrosine kinase inhibitor lorlatinib (Chen 2021, ",
    "n = 328, the complete safety analysis set of the phase I/II study ",
    "B7461001, NCT01970865, 10-200 mg orally once daily or 35-100 mg ",
    "twice daily). The probability of the event is ",
    "expit(-4.757 + 0.029 * WT + 0.003 * TE) where WT is baseline body ",
    "weight in kg and TE is time on study from first dose up to the ",
    "event in days (the full on-study duration for patients in whom ",
    "the event never occurs) (Chen 2021 Table 3). NO LORLATINIB ",
    "EXPOSURE TERM IS PRESENT and that absence is the paper's result, ",
    "not a transcription gap: none of the nine candidate exposure ",
    "metrics reached significance in the univariate forward screen for ",
    "this endpoint, so the final model retains only the base-model ",
    "covariates. This is therefore a pure baseline-risk model -- ",
    "heavier patients and longer time on study are more likely to ",
    "record grade >= 2 weight gain, independent of lorlatinib plasma ",
    "concentration. There is no PK layer and no ODE. No between-subject ",
    "random effect and no residual error are estimated (Bernoulli ",
    "likelihood). Four companion exposure-response models in the ",
    "Chen_2021_lorlatinib_* family."
  )
  reference <- paste(
    "Chen J, Ruiz-Garcia A, James LP, Peltz G, Thurm H, Clancy J, Hibma J.",
    "Lorlatinib exposure-response analyses for safety and efficacy in a",
    "phase I/II trial to support benefit-risk assessment in non-small cell",
    "lung cancer.",
    "Clin Pharmacol Ther. 2021;110(5):1273-1281.",
    "doi:10.1002/cpt.2228.",
    sep = " "
  )
  vignette <- "Chen_2021_lorlatinib_exposure_response"
  units <- list(
    time          = "n/a (static landmark safety regression; no time dimension. The treatment-duration covariate T_FIRSTDOSE is carried in canonical hours and divided by 24 inside model() to recover the paper's days)",
    dosing        = "n/a (no dose events, and no exposure covariate: this endpoint retained no lorlatinib exposure metric)",
    concentration = "prob_weight_gain (probability of grade >= 2 weight gain, 0-1; also logit_weight_gain)"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight (BWT).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT centred and NOT allometrically scaled: Chen 2021 fits the raw",
        "kilogram value as a linear term on the logit, so the intercept",
        "-4.757 is the logit at WT = 0 kg and TE = 0 days rather than at",
        "any clinically meaningful patient. Safety analysis set (N = 328)",
        "median 66.79 kg, range 31.80-155.50 (Chen 2021 Table 1). The",
        "positive coefficient means heavier patients at baseline were MORE",
        "likely to record grade >= 2 weight gain, which is a plausible",
        "consequence of the CTCAE definition being a percentage change",
        "from baseline applied to a heterogeneous cohort. This is the only",
        "one of the four Chen 2021 safety endpoints that retains body",
        "weight."
      ),
      source_name        = "BWT (baseline body weight)"
    ),
    T_FIRSTDOSE = list(
      description        = "Time on study from the first lorlatinib dose up to the event (TE).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Canonical units are hours, so model() divides by 24 to recover",
        "the DAYS in which Chen 2021 estimates the coefficient (0.003 per",
        "day, odds ratio 1.003 per day) -- exactly the convention the",
        "T_FIRSTDOSE register entry prescribes. Right-censored at the",
        "event: Chen 2021 Methods states that 'for patients in whom the AE",
        "never occurs, the TE would be the entire length of on-study",
        "treatment', so this column is a per-subject scalar, not a",
        "time-varying clock, and the model is a landmark regression rather",
        "than a time-to-event model. Chen 2021 included TE a priori in all",
        "four safety base models on the rationale that 'the longer the",
        "patient is on therapy, the greater the risk of experiencing",
        "either of these safety end points, independent of the level of",
        "lorlatinib exposure'."
      ),
      source_name        = "TE (time from first dose up to the event, days)"
    )
  )

  covariatesDataExcluded <- list(
    CTROUGH = list(
      description = "Individual lorlatinib trough plasma concentration at steady state.",
      units       = "ng/mL",
      type        = "continuous",
      notes       = paste(
        "One of the nine lorlatinib exposure metrics screened for this",
        "endpoint (Chen 2021 Methods, 'Selection of lorlatinib exposure",
        "metrics'). Chen 2021 Results states plainly: 'For the safety end",
        "points weight gain grade >= 2 and hypertriglyceridemia grade",
        ">= 3, no significant E-R relationships were identified for",
        "lorlatinib.' No point estimate exists because no exposure term",
        "entered the final model -- this is a published null result, not a",
        "reporting gap, and the correct encoding is the omission of the",
        "term rather than a fixed(0) coefficient. The screened set was:",
        "Cmax cycle 1, CAUC cycle 1, AUC single dose, Ctrough cycle 1,",
        "Cmax event, CAUC prior, AUCss, Ctrough ss and CAUC complete, each",
        "in linear and log-transformed form."
      )
    ),
    TCHOL = list(
      description = "Baseline total serum cholesterol.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Chen 2021 Methods includes 'the baseline laboratory parameters",
        "associated with the respective safety end point' in each safety",
        "base model; baseline cholesterol is the paired laboratory value",
        "for the lipid endpoints, not for weight gain, and is not reported",
        "in this endpoint's Table 3 row."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator; 1 = Asian, 0 = other.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained for this endpoint; retained only in the",
        "companion hypertriglyceridemia model. Safety analysis set Asian",
        "110/328 (34%) (Chen 2021 Table 1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 328L,
    n_studies      = 1L,
    n_observations = "328 binary grade->=-2-weight-gain records (one per patient); the complete safety analysis set (Chen 2021 Table 3 n/N = 328/328), because body weight was recorded for every patient",
    age_range      = "median 53.00 years, range 19.00-85.00 (Chen 2021 Table 1, safety population)",
    weight_range   = "median 66.79 kg, range 31.80-155.50 (Chen 2021 Table 1, safety population)",
    sex_female_pct = 58.0,
    race_ethnicity = c(White = 51.0, Asian = 34.0, Other = 4.0, Black = 2.0, Missing = 9.0),
    disease_state  = "advanced anaplastic lymphoma kinase (ALK)-positive or c-ROS oncogene 1 (ROS1)-positive non-small cell lung cancer; ECOG performance status 0 (42%), 1 (54%), 2 (3%), 3 (0.3%); 85% had received at least one prior ALK inhibitor and 69% had CNS metastasis prior to or at any time on study",
    dose_range     = "lorlatinib 10, 25, 50, 75, 100, 150 or 200 mg orally once daily, or 35, 75 or 100 mg orally twice daily (phase I dose escalation); 100 mg orally once daily (phase II expansion, the labelled dose)",
    regions        = "Multicentre international phase I/II study B7461001 (NCT01970865); regional breakdown not reported in Chen 2021",
    notes          = paste0(
      "Weight gain was one of only four safety endpoints with an ",
      "incidence above 10% in all treated patients, the threshold Chen ",
      "2021 used to admit an endpoint to the exposure-response ",
      "analysis. The deviance difference versus the null model is ",
      "22.529 -- the smallest of the four safety models, consistent ",
      "with a model carrying no exposure term. Grade >= 2 (not >= 3) is ",
      "the severity threshold for this endpoint alone; the other three ",
      "safety endpoints use grade >= 3."
    )
  )

  ini({
    # ==================================================================
    # Chen 2021 Table 3, row block "Weight gain grade >= 2"
    # (n/N = 328/328, deviance difference vs null 22.529). This
    # endpoint has NO printed regression equation because it retained
    # no exposure term -- Eqs. 3, 4 and 5 cover hypercholesterolemia,
    # TEAE grade >= 3 and IC-ORR respectively. The linear predictor
    # below is assembled from the Table 3 coefficient rows, which are
    # the model's complete parameter set.
    #
    # Each value is the Table 3 "Estimate" column and is reported a
    # second time by its printed odds ratio, the Methods defining OR as
    # "the odds ratio determined by exponentiating the coefficient
    # estimate". BOTH columns are rounded to 3 decimals, so the correct
    # consistency test is an interval one -- does a real b exist with
    # round(b, 3) == Estimate AND round(exp(b), 3) == OR -- rather than
    # round(log(OR), 3) == Estimate, which wrongly treats the OR column
    # as exact. Both rows here pass the interval test; see the note on
    # the BWT row below.
    #
    # Covariates are NOT centred and NOT scaled; the intercept is an
    # extrapolated anchor, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -4.757     ; label("Logit of the grade >= 2 weight gain probability at WT = 0 kg and TE = 0 days (unitless logit)")  # Chen 2021 Table 3, Intercept -4.757 (95% CI -6.3244, -3.3327), P < 0.0001

    # ----- Covariate effects on the logit -----
    # The BWT row LOOKS inconsistent and is not. exp(0.029) = 1.02942
    # rounds to 1.029 while Table 3 prints OR 1.030, and
    # round(log(1.030), 3) = 0.030 does not return the printed 0.029 --
    # but both printed columns are rounded, and any true coefficient in
    # [0.02907, 0.02950) rounds to 0.029 AND exponentiates into
    # [1.0295, 1.0305), which rounds to 1.030. The two columns are
    # therefore mutually consistent; see vignette Errata item 3.
    #
    # Do NOT try to recover the unrounded estimate from the printed 95%
    # CI midpoint (0.0295): every interval in Chen 2021 is asymmetric
    # about its estimate, as profile-likelihood intervals from R's
    # confint.glm() are, so a CI midpoint is not an estimator here
    # (vignette Errata item 4).
    e_wt_logit          <- 0.029  ; label("Log-odds of grade >= 2 weight gain per 1 kg increase in baseline body weight (unitless logit)")  # Chen 2021 Table 3, Estimate 0.029 (0.0108, 0.0482), P = 0.0021; printed OR 1.030, consistent with the Estimate under the interval test above (exp(0.029) = 1.0294)
    e_t_firstdose_logit <- 0.003  ; label("Log-odds of grade >= 2 weight gain per additional day on study prior to the event (unitless logit)")  # Chen 2021 Table 3, Estimate 0.003 (0.0013, 0.0050), P = 0.0007; printed OR 1.003 = exp(0.003)

    # ----- No between-subject variability, no residual error -----
    # The source is a generalized binomial logistic regression fitted
    # with R's glm(family = "binomial"); a Bernoulli likelihood has no
    # sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_weight_gain <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # T_FIRSTDOSE is carried in canonical hours; Chen 2021 estimates the
    # TE coefficient per DAY (Table 3 footnote: "TE, time from first dose
    # up to the event (days)").
    te_days <- T_FIRSTDOSE / 24

    # ----- Linear predictor (Chen 2021 Table 3, weight-gain block) -----
    # No exposure term: none of the nine candidate lorlatinib exposure
    # metrics was significant for this endpoint.
    logit_weight_gain <- logit_ref +
      e_wt_logit          * WT +
      e_t_firstdose_logit * te_days

    prob_weight_gain <- expit(logit_weight_gain)

    # ----- Observation -----
    prob_weight_gain ~ add(addSd_prob_weight_gain)
  })
}
