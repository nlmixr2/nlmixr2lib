Chen_2021_lorlatinib_teae_grade3 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for ANY CTCAE ",
    "grade >= 3 treatment-emergent adverse event (composite across all ",
    "preferred terms) in adults with advanced ALK-positive or ",
    "ROS1-positive non-small cell lung cancer treated with the ",
    "third-generation ALK/ROS1 tyrosine kinase inhibitor lorlatinib ",
    "(Chen 2021, n = 328, the complete safety analysis set of the ",
    "phase I/II study B7461001, NCT01970865, 10-200 mg orally once ",
    "daily or 35-100 mg twice daily). The probability of the event is ",
    "expit(-7.995 + 0.012 * TCHOL + 0.012 * TE + 1.167 * ",
    "log(CTROUGH)) where TCHOL is baseline total cholesterol in mg/dL, ",
    "TE is time on study from first dose up to the event in days (the ",
    "full on-study duration for patients in whom the event never ",
    "occurs), and CTROUGH is the individual lorlatinib trough plasma ",
    "concentration at steady state in ng/mL (Chen 2021 Eq. 4, ",
    "Table 3). The log is a natural log, stated explicitly by the ",
    "paper: every unit increase in log(Ctrough,ss) multiplies the odds ",
    "of a grade >= 3 TEAE by 3.214. This composite endpoint overlaps ",
    "heavily with the companion hypercholesterolemia model -- Chen ",
    "2021 notes that because hypercholesterolemia was so common, many ",
    "of the grade >= 3 TEAEs WERE hypercholesterolemia events -- so ",
    "the two must not be treated as competing risks. There is no PK ",
    "layer and no ODE: the exposure metric is supplied as a data ",
    "column, derived in the source analysis from the companion ",
    "lorlatinib population PK model with time-varying clearance ",
    "(doi:10.1002/psp4.12585, packaged as Chen_2021_lorlatinib). No ",
    "between-subject random effect and no residual error are ",
    "estimated (Bernoulli likelihood). Four companion ",
    "exposure-response models in the Chen_2021_lorlatinib_* family."
  )
  reference <- paste(
    "Chen J, Ruiz-Garcia A, James LP, Peltz G, Thurm H, Clancy J, Hibma J.",
    "Lorlatinib exposure-response analyses for safety and efficacy in a",
    "phase I/II trial to support benefit-risk assessment in non-small cell",
    "lung cancer.",
    "Clin Pharmacol Ther. 2021;110(5):1273-1281.",
    "doi:10.1002/cpt.2228.",
    "Individual lorlatinib exposure metrics derive from the companion",
    "population pharmacokinetic model: Chen J, Houk B, Pithavala YK,",
    "Ruiz-Garcia A. Population pharmacokinetic model with time-varying",
    "clearance for lorlatinib using pooled data from patients with",
    "non-small cell lung cancer and healthy participants.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(2):148-160.",
    "doi:10.1002/psp4.12585.",
    sep = " "
  )
  vignette <- "Chen_2021_lorlatinib_exposure_response"
  units <- list(
    time          = "n/a (static landmark exposure-safety regression; no time dimension. The treatment-duration covariate T_FIRSTDOSE is carried in canonical hours and divided by 24 inside model() to recover the paper's days)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_teae_grade3 (probability of any grade >= 3 treatment-emergent adverse event, 0-1; also logit_teae_grade3)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual lorlatinib trough plasma concentration at steady state (Ctrough,ss), per subject. Supplied as data: this model has no PK layer, and the source analysis used individual predictions from the companion lorlatinib population PK model with time-varying clearance (Chen 2021, doi:10.1002/psp4.12585).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration, at STEADY STATE rather",
        "than at the end of cycle 1 -- Chen 2021 screened both Ctrough,ss",
        "and Ctrough,cycle 1 in the univariate exposure-metric forward",
        "screen and carried the steady-state metric into the final model.",
        "Because lorlatinib clearance auto-induces over roughly the first",
        "week (Chen 2021 popPK, doi:10.1002/psp4.12585), steady-state",
        "troughs are materially LOWER than day-1 troughs; supplying a",
        "single-dose trough would bias this model high. Enters on the",
        "NATURAL LOG scale, stated outright by the paper: 'the",
        "significantly related exposure metric, log(Ctrough ss) is natural",
        "log transformed; and thus, the increases are on the log scale'.",
        "The printed odds ratio 3.214 is therefore per e-fold, not per",
        "ng/mL. Units are load-bearing: the intercept absorbs the unit",
        "choice, so supplying ug/mL instead of ng/mL would shift the logit",
        "by 1.167 * log(1000) = 8.1. For calibration, the mean (SD)",
        "lorlatinib trough over cycle 1 was 114.97 (40.28) ng/mL in the 132",
        "ALK-inhibitor-pretreated patients with baseline CNS metastasis",
        "dosed at 100 mg q.d. (Chen 2021 Results, Exposure-efficacy)."
      ),
      source_name        = "Ctrough ss (trough concentration at steady-state)"
    ),
    TCHOL = list(
      description        = "Baseline total serum cholesterol (BCHOL).",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT centred: Chen 2021 fits the raw covariate value, so the",
        "intercept -7.995 is the logit at TCHOL = 0, TE = 0 and",
        "CTROUGH = 1 ng/mL rather than at any clinically meaningful",
        "patient. mg/dL, not mmol/L -- the coefficient 0.012 per unit is",
        "only sensible on the mg/dL scale, and the paper's own reference",
        "value of 193 mg/dL is quoted in mg/dL. Safety analysis set",
        "(N = 328) median 193.00 mg/dL, range 3.00-321.00, mean 192.95",
        "(SD 44.06), 28 patients missing (Chen 2021 Table 2). 193 mg/dL is",
        "the value the paper fixes this covariate to when drawing the",
        "Figure 2 exposure-response curve. Baseline cholesterol predicts a",
        "composite all-cause severe-TEAE endpoint because",
        "hypercholesterolemia itself is the dominant contributor to that",
        "composite in this trial."
      ),
      source_name        = "BCHOL (baseline cholesterol)"
    ),
    T_FIRSTDOSE = list(
      description        = "Time on study from the first lorlatinib dose up to the event (TE).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Canonical units are hours, so model() divides by 24 to recover",
        "the DAYS in which Chen 2021 estimates the coefficient (0.012 per",
        "day, odds ratio 1.012 per day) -- exactly the convention the",
        "T_FIRSTDOSE register entry prescribes. Right-censored at the",
        "event: Chen 2021 Methods states that 'for patients in whom the AE",
        "never occurs, the TE would be the entire length of on-study",
        "treatment', so this column is a per-subject scalar, not a",
        "time-varying clock, and the model is a landmark regression rather",
        "than a time-to-event model. Analysis-population median 38.75 days",
        "= 930 h, the value the paper fixes this covariate to for the",
        "Figure 2 curve. The TE effect here (0.012 per day) is three times",
        "the hypercholesterolemia model's (0.004 per day), i.e. cumulative",
        "on-study time matters far more for the composite endpoint."
      ),
      source_name        = "TE (time from first dose up to the event, days)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened across all four safety endpoints (Chen 2021 Table S4",
        "candidate-covariate list) but retained only in the companion",
        "weight-gain model. Safety analysis set median 66.79 kg, range",
        "31.80-155.50 (Chen 2021 Table 1)."
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
    n_observations = "328 binary any-grade->=-3-TEAE records (one per patient); this is the only one of the four Chen 2021 safety models fitted to the complete safety analysis set (Table 3 n/N = 328/328), because it needs no baseline laboratory covariate that could be missing",
    age_range      = "median 53.00 years, range 19.00-85.00 (Chen 2021 Table 1, safety population)",
    weight_range   = "median 66.79 kg, range 31.80-155.50 (Chen 2021 Table 1, safety population)",
    sex_female_pct = 58.0,
    race_ethnicity = c(White = 51.0, Asian = 34.0, Other = 4.0, Black = 2.0, Missing = 9.0),
    disease_state  = "advanced anaplastic lymphoma kinase (ALK)-positive or c-ROS oncogene 1 (ROS1)-positive non-small cell lung cancer; ECOG performance status 0 (42%), 1 (54%), 2 (3%), 3 (0.3%); 85% had received at least one prior ALK inhibitor and 69% had CNS metastasis prior to or at any time on study",
    dose_range     = "lorlatinib 10, 25, 50, 75, 100, 150 or 200 mg orally once daily, or 35, 75 or 100 mg orally twice daily (phase I dose escalation); 100 mg orally once daily (phase II expansion, the labelled dose)",
    regions        = "Multicentre international phase I/II study B7461001 (NCT01970865); regional breakdown not reported in Chen 2021",
    baseline_total_cholesterol = "median 193.00 mg/dL, range 3.00-321.00, mean 192.95 (SD 44.06), 28 missing (Chen 2021 Table 2, safety population)",
    concomitant_medication     = "concomitant statin therapy 266/328 (81%), steroid therapy 139 (42%), narcotics 164 (50%) (Chen 2021 Table 1)",
    notes          = paste0(
      "Chen 2021 assessed exposure-response only for safety endpoints ",
      "with an incidence above 10% in all treated patients. This ",
      "composite endpoint is a SUPERSET of the paper's three ",
      "specific-toxicity safety endpoints (hypercholesterolemia, ",
      "hypertriglyceridemia, weight gain) and must not be modelled ",
      "alongside them as a competing risk. Together with the ",
      "hypercholesterolemia model it is one of only two significant ",
      "exposure-safety relationships the paper identified; no ",
      "exposure-efficacy relationship was significant for either ",
      "endpoint."
    )
  )

  ini({
    # ==================================================================
    # Chen 2021 Table 3, row block "TEAE grade >= 3" (n/N = 328/328,
    # deviance difference vs null 55.680), and the printed regression
    # Eq. 4:
    #
    #   logit(p) = -7.995 + 0.012*BCHOL + 0.012*TE
    #              + 1.167*log(Ctrough ss)
    #
    # Every coefficient below is the Table 3 "Estimate" column. The
    # table is internally redundant -- its "OR" column is defined by
    # the Methods as "the odds ratio determined by exponentiating the
    # coefficient estimate" -- so each value is confirmed a second time
    # by its printed OR, and a third time by Eq. 4. Note that BOTH
    # printed columns are rounded to 3 decimals, so the consistency
    # test is an interval one -- does a real b exist with
    # round(b, 3) == Estimate AND round(exp(b), 3) == OR -- rather than
    # exp(printed Estimate) == printed OR, which does not hold for the
    # exposure row (exp(1.167) = 3.212 against a printed 3.214), nor
    # round(log(OR), 3) == Estimate, which is also too crude here
    # (log(3.214) rounds to 1.168, not the printed 1.167). All three
    # rows pass the interval test. The checks are re-run mechanically
    # in the validation vignette.
    #
    # Covariates are NOT centred and NOT scaled; the intercept is an
    # extrapolated anchor, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -7.995     ; label("Logit of the any grade >= 3 TEAE probability at TCHOL = 0 mg/dL, TE = 0 days and CTROUGH = 1 ng/mL (unitless logit)")  # Chen 2021 Table 3, Intercept -7.995 (95% CI -12.2153, -4.0263), P = 0.0001; Eq. 4

    # ----- Covariate effects on the logit -----
    e_tchol_logit       <- 0.012  ; label("Log-odds of any grade >= 3 TEAE per 1 mg/dL increase in baseline total cholesterol (unitless logit)")  # Chen 2021 Table 3, Estimate 0.012 (0.0058, 0.0191), P = 0.0003; printed OR 1.012 = exp(0.012)
    e_t_firstdose_logit <- 0.012  ; label("Log-odds of any grade >= 3 TEAE per additional day on study prior to the event (unitless logit)")  # Chen 2021 Table 3, Estimate 0.012 (0.0078, 0.0177), P < 0.0001; printed OR 1.012 = exp(0.012)
    e_ctrough_logit     <- 1.167  ; label("Log-odds of any grade >= 3 TEAE per e-fold increase in lorlatinib steady-state trough concentration (unitless logit)")  # Chen 2021 Table 3, Estimate 1.167 (0.4012, 1.9725), P = 0.0035; printed OR 3.214, consistent with the Estimate under the interval test above (exp(1.167) = 3.212; the OR column carries the unrounded estimate)

    # ----- No between-subject variability, no residual error -----
    # The source is a generalized binomial logistic regression fitted
    # with R's glm(family = "binomial"); a Bernoulli likelihood has no
    # sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_teae_grade3 <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # T_FIRSTDOSE is carried in canonical hours; Chen 2021 estimates the
    # TE coefficient per DAY (Table 3 footnote: "TE, time from first dose
    # up to the event (days)").
    te_days <- T_FIRSTDOSE / 24

    # ----- Linear predictor (Chen 2021 Eq. 4) -----
    # log() is rxode2's natural log, which the paper states explicitly
    # for this exposure metric.
    logit_teae_grade3 <- logit_ref +
      e_tchol_logit       * TCHOL +
      e_t_firstdose_logit * te_days +
      e_ctrough_logit     * log(CTROUGH)

    prob_teae_grade3 <- expit(logit_teae_grade3)

    # ----- Observation -----
    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
