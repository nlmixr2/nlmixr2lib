Chen_2021_lorlatinib_hypercholesterolemia <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for CTCAE ",
    "grade >= 3 hypercholesterolemia in adults with advanced ",
    "ALK-positive or ROS1-positive non-small cell lung cancer treated ",
    "with the third-generation ALK/ROS1 tyrosine kinase inhibitor ",
    "lorlatinib (Chen 2021, n = 298 of the 328-patient safety analysis ",
    "set of the phase I/II study B7461001, NCT01970865, 10-200 mg ",
    "orally once daily or 35-100 mg twice daily). The probability of ",
    "the event is expit(-18.829 + 0.029 * TCHOL + 0.004 * TE + 1.659 * ",
    "log(CMAX)) where TCHOL is baseline total cholesterol in mg/dL, TE ",
    "is time on study from first dose up to the event in days (the ",
    "full on-study duration for patients in whom the event never ",
    "occurs), and CMAX is the maximum observed lorlatinib plasma ",
    "concentration prior to the adverse event in ng/mL (Chen 2021 ",
    "Eq. 3, Table 3). The log is a natural log. This is the stronger ",
    "of the paper's two significant exposure-safety relationships: ",
    "every e-fold rise in Cmax,event multiplies the odds of grade ",
    ">= 3 hypercholesterolemia by 5.256. There is no PK layer and no ",
    "ODE: the exposure metric is supplied as a data column, derived in ",
    "the source analysis from the companion lorlatinib population PK ",
    "model with time-varying clearance (doi:10.1002/psp4.12585, ",
    "packaged as Chen_2021_lorlatinib). No between-subject random ",
    "effect and no residual error are estimated (Bernoulli ",
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
    dosing        = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_hypercholesterolemia (probability of grade >= 3 hypercholesterolemia, 0-1; also logit_hypercholesterolemia)"
  )

  covariateData <- list(
    CMAX = list(
      description        = "Maximum observed lorlatinib plasma concentration prior to the adverse event (Cmax,event), per subject. Supplied as data: this model has no PK layer, and the source analysis used individual predictions from the companion lorlatinib population PK model with time-varying clearance (Chen 2021, doi:10.1002/psp4.12585).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration. NOT a fixed-interval",
        "Cmax: the averaging window ends at the adverse event, so the",
        "window length is subject-specific and equals the whole on-study",
        "duration for patients in whom the event never occurred (the same",
        "window convention as the paired T_FIRSTDOSE column). Enters on the",
        "NATURAL LOG scale (Chen 2021 Eq. 3; the Results narrative states",
        "for the companion TEAE model that the exposure metric 'is natural",
        "log transformed', and the same convention applies here). The",
        "printed odds ratio 5.256 is therefore per e-fold, not per ng/mL.",
        "Units are load-bearing: the intercept absorbs the unit choice, so",
        "supplying ug/mL instead of ng/mL would shift the logit by",
        "1.659 * log(1000) = 11.5 and destroy the fit. Chen 2021 reports",
        "every lorlatinib concentration in ng/mL; for calibration, the mean",
        "(SD) Cmax over cycle 1 was 687.11 (141.09) ng/mL in the 197",
        "ALK-inhibitor-pretreated patients dosed at 100 mg q.d.",
        "(Chen 2021 Results, Exposure-efficacy)."
      ),
      source_name        = "Cmax event (maximum observed concentration prior to the AE)"
    ),
    TCHOL = list(
      description        = "Baseline total serum cholesterol (BCHOL).",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "NOT centred: Chen 2021 fits the raw covariate value, so the",
        "intercept -18.829 is the logit at TCHOL = 0, TE = 0 and",
        "CMAX = 1 ng/mL rather than at any clinically meaningful patient.",
        "mg/dL, not mmol/L -- the coefficient 0.029 per unit is only",
        "sensible on the mg/dL scale (0.029 per mmol/L would be a",
        "negligible effect, and the paper's own reference value of",
        "193 mg/dL is quoted in mg/dL). Safety analysis set (N = 328)",
        "median 193.00 mg/dL, range 3.00-321.00, mean 192.95 (SD 44.06),",
        "28 patients missing (Chen 2021 Table 2). 193 mg/dL is the value",
        "the paper fixes this covariate to when drawing the Figure 1",
        "exposure-response curve."
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
        "the DAYS in which Chen 2021 estimates the coefficient (0.004 per",
        "day, odds ratio 1.004 per day) -- exactly the convention the",
        "T_FIRSTDOSE register entry prescribes. Right-censored at the",
        "event: Chen 2021 Methods states that 'for patients in whom the AE",
        "never occurs, the TE would be the entire length of on-study",
        "treatment', so this column is a per-subject scalar, not a",
        "time-varying clock, and the model is a landmark regression rather",
        "than a time-to-event model. Rationale for inclusion (Chen 2021",
        "Results): 'the longer the patient is on therapy, the more blood",
        "cholesterol will accumulate, until reaching hypercholesterolemia",
        "grade >= 3'. Analysis-population median 41 days = 984 h, the value",
        "the paper fixes this covariate to for the Figure 1 curve."
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
        "candidate-covariate list; Table 2 summarises it for every",
        "analysis population) but retained only in the companion weight-gain",
        "model, not here. Safety analysis set median 66.79 kg, range",
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
        "110/328 (34%), White 168 (51%), Black 6 (2%), Other 13 (4%),",
        "Missing 31 (9%) (Chen 2021 Table 1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 298L,
    n_studies      = 1L,
    n_observations = "298 evaluable binary grade->=-3-hypercholesterolemia records (one per patient) out of the 328-patient safety analysis set; the 30-patient shortfall is the Chen 2021 Table 3 'n/N' column (298/328) and reflects patients without a complete baseline-cholesterol record",
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
      "with an incidence above 10% in all treated patients. ",
      "Hypercholesterolemia was the most common treatment-related ",
      "adverse event in the phase II portion (81% any grade, 16% grade ",
      "3-4). The paper notes that because hypercholesterolemia was so ",
      "common, many of the grade >= 3 TEAEs captured by the companion ",
      "Chen_2021_lorlatinib_teae model WERE hypercholesterolemia events ",
      "-- the two endpoints overlap heavily and must not be treated as ",
      "competing risks."
    )
  )

  ini({
    # ==================================================================
    # Chen 2021 Table 3, row block "Hypercholesterolemia grade >= 3"
    # (n/N = 298/328, deviance difference vs null 53.361), and the
    # printed regression Eq. 3:
    #
    #   logit(p) = -18.829 + 0.029*BCHOL + 0.004*TE
    #              + 1.659*log(Cmax event)
    #
    # Every coefficient below is the Table 3 "Estimate" column. The
    # table is internally redundant -- its "OR" column is defined by
    # the Methods as "the odds ratio determined by exponentiating the
    # coefficient estimate" -- so each value is confirmed a second time
    # by its printed OR, and a third time by Eq. 3. Note that BOTH
    # printed columns are rounded to 3 decimals, so the consistency
    # test is an interval one -- does a real b exist with
    # round(b, 3) == Estimate AND round(exp(b), 3) == OR -- rather than
    # exp(printed Estimate) == printed OR, which does not hold for the
    # larger coefficients (exp(1.659) = 5.254 against a printed 5.256).
    # All three rows here pass. The checks are re-run mechanically in
    # the validation vignette.
    #
    # NOTE the covariates are NOT centred and NOT scaled: Chen 2021
    # fits raw BCHOL (mg/dL), raw TE (days) and log(Cmax) (ng/mL). The
    # intercept is therefore an extrapolated anchor, not a reference
    # patient probability, which is the structural difference from the
    # Fukae_2024_valemetostat_* logistic family.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -18.829     ; label("Logit of the grade >= 3 hypercholesterolemia probability at TCHOL = 0 mg/dL, TE = 0 days and CMAX = 1 ng/mL (unitless logit)")  # Chen 2021 Table 3, Intercept -18.829 (95% CI -30.4373, -8.1449), P = 0.0009; Eq. 3

    # ----- Covariate effects on the logit -----
    e_tchol_logit       <- 0.029  ; label("Log-odds of grade >= 3 hypercholesterolemia per 1 mg/dL increase in baseline total cholesterol (unitless logit)")  # Chen 2021 Table 3, Estimate 0.029 (0.0199, 0.0386), P < 0.0001; printed OR 1.029 = exp(0.029)
    e_t_firstdose_logit <- 0.004  ; label("Log-odds of grade >= 3 hypercholesterolemia per additional day on study prior to the event (unitless logit)")  # Chen 2021 Table 3, Estimate 0.004 (0.0001, 0.0069), P = 0.0413; printed OR 1.004 = exp(0.004)
    e_cmax_logit        <- 1.659  ; label("Log-odds of grade >= 3 hypercholesterolemia per e-fold increase in lorlatinib Cmax prior to the event (unitless logit)")  # Chen 2021 Table 3, Estimate 1.659 (0.0762, 3.3330), P = 0.0452; printed OR 5.256, consistent with the Estimate under the interval test above (exp(1.659) = 5.254; the OR column carries the unrounded estimate)

    # ----- No between-subject variability, no residual error -----
    # The source is a generalized binomial logistic regression fitted
    # with R's glm(family = "binomial"); a Bernoulli likelihood has no
    # sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_hypercholesterolemia <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # T_FIRSTDOSE is carried in canonical hours; Chen 2021 estimates the
    # TE coefficient per DAY (Table 3 footnote: "TE, time from first dose
    # up to the event (days)").
    te_days <- T_FIRSTDOSE / 24

    # ----- Linear predictor (Chen 2021 Eq. 3) -----
    # log() is rxode2's natural log, matching the paper's convention.
    logit_hypercholesterolemia <- logit_ref +
      e_tchol_logit       * TCHOL +
      e_t_firstdose_logit * te_days +
      e_cmax_logit        * log(CMAX)

    prob_hypercholesterolemia <- expit(logit_hypercholesterolemia)

    # ----- Observation -----
    prob_hypercholesterolemia ~ add(addSd_prob_hypercholesterolemia)
  })
}
