Chen_2021_lorlatinib_hypertriglyceridemia <- function() {
  description <- paste0(
    "Binomial logistic-regression safety model for CTCAE grade >= 3 ",
    "hypertriglyceridemia in adults with advanced ALK-positive or ",
    "ROS1-positive non-small cell lung cancer treated with the ",
    "third-generation ALK/ROS1 tyrosine kinase inhibitor lorlatinib ",
    "(Chen 2021, n = 298 of the 328-patient safety analysis set of the ",
    "phase I/II study B7461001, NCT01970865, 10-200 mg orally once ",
    "daily or 35-100 mg twice daily). The probability of the event is ",
    "expit(-5.113 + 1.011 * RACE_ASIAN + 0.003 * TE + 0.018 * TRIG) ",
    "where RACE_ASIAN is an Asian-race indicator, TE is time on study ",
    "from first dose up to the event in days (the full on-study ",
    "duration for patients in whom the event never occurs), and TRIG ",
    "is baseline serum triglycerides in mg/dL (Chen 2021 Table 3). NO ",
    "LORLATINIB EXPOSURE TERM IS PRESENT and that absence is the ",
    "paper's result, not a transcription gap: none of the nine ",
    "candidate exposure metrics reached significance in the univariate ",
    "forward screen for this endpoint. This is therefore a pure ",
    "baseline-risk model, and it is the only Chen 2021 model to retain ",
    "a demographic covariate -- Asian patients had 2.749 times the ",
    "odds of grade >= 3 hypertriglyceridemia. There is no PK layer and ",
    "no ODE. No between-subject random effect and no residual error ",
    "are estimated (Bernoulli likelihood). Four companion ",
    "exposure-response models in the Chen_2021_lorlatinib_* family."
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
    concentration = "prob_hypertriglyceridemia (probability of grade >= 3 hypertriglyceridemia, 0-1; also logit_hypertriglyceridemia)"
  )

  covariateData <- list(
    TRIG = list(
      description        = "Baseline serum triglyceride concentration (BTG).",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT centred: Chen 2021 fits the raw covariate value, so the",
        "intercept -5.113 is the logit at TRIG = 0 mg/dL, TE = 0 days and",
        "a non-Asian patient rather than at any clinically meaningful",
        "patient. mg/dL, not mmol/L -- the coefficient 0.018 per unit is",
        "only sensible on the mg/dL scale, and Chen 2021 Table 2 reports",
        "triglycerides in mg/dL. Safety analysis set (N = 328) median",
        "107.50 mg/dL, range 28.00-451.40, mean 123.98 (SD 66.59), 30",
        "patients missing (Chen 2021 Table 2). This is the paired baseline",
        "laboratory value Chen 2021 Methods includes a priori in each",
        "safety base model ('the baseline laboratory parameters associated",
        "with the respective safety end point'), the triglyceride analogue",
        "of TCHOL in the companion hypercholesterolemia model, and it is",
        "the strongest predictor in this model (P < 0.0001)."
      ),
      source_name        = "BTG (baseline triglycerides)"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (ASN1); 1 = Asian, 0 = other.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; Chen 2021 Table 3 abbreviation list defines the covariate simply as 'ASN1, Asian', so the complement is the pooled non-Asian group)",
      notes              = paste(
        "The only demographic covariate retained in any of the five Chen",
        "2021 models. Enters as a plain log-odds shift, odds ratio 2.749",
        "for Asian versus non-Asian. Safety analysis set (N = 328) race:",
        "White 168 (51%), Asian 110 (34%), Other 13 (4%), Black 6 (2%),",
        "Missing 31 (9%) (Chen 2021 Table 1). Because 9% of the cohort has",
        "a MISSING race, the reference group of the fitted indicator is",
        "strictly 'recorded as non-Asian or missing', which a downstream",
        "user reproducing the analysis dataset should mirror rather than",
        "dropping missing-race patients. Note the direction is opposite to",
        "the companion lorlatinib population PK model",
        "(doi:10.1002/psp4.12585), which screened Asian race on clearance",
        "and did NOT retain it: the effect here is on baseline lipid",
        "susceptibility, not on drug exposure."
      ),
      source_name        = "ASN1 (Asian)"
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
        "than a time-to-event model."
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
        "term rather than a fixed(0) coefficient."
      )
    ),
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened across all four safety endpoints but retained only in",
        "the companion weight-gain model. Safety analysis set median",
        "66.79 kg, range 31.80-155.50 (Chen 2021 Table 1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 298L,
    n_studies      = 1L,
    n_observations = "298 evaluable binary grade->=-3-hypertriglyceridemia records (one per patient) out of the 328-patient safety analysis set; the 30-patient shortfall is the Chen 2021 Table 3 'n/N' column (298/328) and matches exactly the 30 patients with a missing baseline triglyceride value in Table 2",
    age_range      = "median 53.00 years, range 19.00-85.00 (Chen 2021 Table 1, safety population)",
    weight_range   = "median 66.79 kg, range 31.80-155.50 (Chen 2021 Table 1, safety population)",
    sex_female_pct = 58.0,
    race_ethnicity = c(White = 51.0, Asian = 34.0, Other = 4.0, Black = 2.0, Missing = 9.0),
    disease_state  = "advanced anaplastic lymphoma kinase (ALK)-positive or c-ROS oncogene 1 (ROS1)-positive non-small cell lung cancer; ECOG performance status 0 (42%), 1 (54%), 2 (3%), 3 (0.3%); 85% had received at least one prior ALK inhibitor and 69% had CNS metastasis prior to or at any time on study",
    dose_range     = "lorlatinib 10, 25, 50, 75, 100, 150 or 200 mg orally once daily, or 35, 75 or 100 mg orally twice daily (phase I dose escalation); 100 mg orally once daily (phase II expansion, the labelled dose)",
    regions        = "Multicentre international phase I/II study B7461001 (NCT01970865); regional breakdown not reported in Chen 2021",
    baseline_triglycerides = "median 107.50 mg/dL, range 28.00-451.40, mean 123.98 (SD 66.59), 30 missing (Chen 2021 Table 2, safety population)",
    concomitant_medication = "concomitant statin therapy 266/328 (81%), steroid therapy 139 (42%), narcotics 164 (50%) (Chen 2021 Table 1)",
    notes          = paste0(
      "Hypertriglyceridemia was the second most common ",
      "treatment-related adverse event in the phase II portion (60% ",
      "any grade, 16% grade 3-4) and so cleared Chen 2021's 10% ",
      "incidence threshold for admission to the exposure-response ",
      "analysis. Its deviance difference versus the null model is ",
      "62.215, the LARGEST of the four safety models, even though it ",
      "carries no exposure term -- the fit is driven almost entirely ",
      "by baseline triglycerides. Chen 2021 notes these lipid events ",
      "were well managed on study with statin / lipid-lowering therapy ",
      "and lorlatinib dose modifications."
    )
  )

  ini({
    # ==================================================================
    # Chen 2021 Table 3, row block "Hypertriglyceridemia grade >= 3"
    # (n/N = 298/328, deviance difference vs null 62.215). This
    # endpoint has NO printed regression equation because it retained
    # no exposure term -- Eqs. 3, 4 and 5 cover hypercholesterolemia,
    # TEAE grade >= 3 and IC-ORR respectively. The linear predictor
    # below is assembled from the Table 3 coefficient rows, which are
    # the model's complete parameter set.
    #
    # Each value is the Table 3 "Estimate" column and is confirmed a
    # second time by its printed odds ratio, the Methods defining OR as
    # "the odds ratio determined by exponentiating the coefficient
    # estimate". BOTH printed columns are rounded to 3 decimals, so the
    # consistency test is an interval one -- does a real b exist with
    # round(b, 3) == Estimate AND round(exp(b), 3) == OR -- rather than
    # exp(printed Estimate) == printed OR, which does not hold for the
    # race row (exp(1.011) = 2.748 against a printed 2.749). All three
    # rows pass the interval test.
    #
    # Covariates are NOT centred and NOT scaled; the intercept is an
    # extrapolated anchor, not a reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -5.113     ; label("Logit of the grade >= 3 hypertriglyceridemia probability for a non-Asian patient at TRIG = 0 mg/dL and TE = 0 days (unitless logit)")  # Chen 2021 Table 3, Intercept -5.113 (95% CI -6.4219, -3.9792), P < 0.0001

    # ----- Covariate effects on the logit -----
    e_race_asian_logit  <- 1.011  ; label("Log-odds of grade >= 3 hypertriglyceridemia for Asian race vs the non-Asian reference (unitless logit)")  # Chen 2021 Table 3, Estimate 1.011 (0.2592, 1.7841), P = 0.0089; printed OR 2.749, consistent with the Estimate under the interval test above (exp(1.011) = 2.748)
    e_t_firstdose_logit <- 0.003  ; label("Log-odds of grade >= 3 hypertriglyceridemia per additional day on study prior to the event (unitless logit)")  # Chen 2021 Table 3, Estimate 0.003 (0.0004, 0.0055), P = 0.0196; printed OR 1.003 = exp(0.003)
    e_trig_logit        <- 0.018  ; label("Log-odds of grade >= 3 hypertriglyceridemia per 1 mg/dL increase in baseline triglycerides (unitless logit)")  # Chen 2021 Table 3, Estimate 0.018 (0.0130, 0.0243), P < 0.0001; printed OR 1.018 = exp(0.018)

    # ----- No between-subject variability, no residual error -----
    # The source is a generalized binomial logistic regression fitted
    # with R's glm(family = "binomial"); a Bernoulli likelihood has no
    # sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_hypertriglyceridemia <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # T_FIRSTDOSE is carried in canonical hours; Chen 2021 estimates the
    # TE coefficient per DAY (Table 3 footnote: "TE, time from first dose
    # up to the event (days)").
    te_days <- T_FIRSTDOSE / 24

    # ----- Linear predictor (Chen 2021 Table 3, hypertriglyceridemia block) --
    # No exposure term: none of the nine candidate lorlatinib exposure
    # metrics was significant for this endpoint.
    logit_hypertriglyceridemia <- logit_ref +
      e_race_asian_logit  * RACE_ASIAN +
      e_t_firstdose_logit * te_days +
      e_trig_logit        * TRIG

    prob_hypertriglyceridemia <- expit(logit_hypertriglyceridemia)

    # ----- Observation -----
    prob_hypertriglyceridemia ~ add(addSd_prob_hypertriglyceridemia)
  })
}
