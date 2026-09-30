Chen_2021_luspatercept_hb_response <- function() {
  description <- paste0(
    "Landmark binomial logistic-regression exposure-response model relating the steady-state luspatercept AUC of the starting dose ",
    "to the probability of a >= 1 g/dL increase in trough hemoglobin on day 21 of the first dosing interval in adults with ",
    "beta-thalassemia and a low transfusion burden (< 12 RBC units/24 weeks) (Chen 2021, Figure 2B; n = 70, pooled phase 2 and 3). ",
    "prob = expit(logit_ref + e_auc_logit * AUCss), with the slope taken from the printed odds ratio of 3.73 per 50 ug*day/mL and ",
    "the intercept digitized from the fitted curve in Figure 2B. There is no PK layer and no ODE: exposure is supplied as the ",
    "AUC_LUSP data column (starting dose / individual CL/F from Chen_2021_luspatercept). No random effects and no residual error ",
    "are estimated (Bernoulli likelihood). One of four exposure-response models in the Chen_2021_luspatercept_* family."
  )
  reference <- paste(
    "Chen N, Kassir N, Laadem A, Giuseppi AC, Shetty J, Maxwell SE, Sriraman P, Ritland S, Linde PG, Budda B, Reynolds JG, Zhou S, Palmisano M.",
    "Population Pharmacokinetics and Exposure-Response Relationship of Luspatercept, an Erythroid Maturation Agent, in Anemic Patients With beta-Thalassemia.",
    "J Clin Pharmacol. 2021;61(1):52-63. doi:10.1002/jcph.1696. PMCID: PMC7754485.",
    "Logistic regression in Figure 2B (panel annotation 'OR (95% CI) = 3.73 (1.81, 9.16); P = 0.00015'; the caption defines OR per 50 units of AUCss);",
    "Results, 'Exposure-Response for Hb'.",
    sep = " "
  )
  vignette <- "Chen_2021_luspatercept"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_LUSP covariate column)",
    concentration = "prob_hb_increase_1g (probability of a >= 1 g/dL increase in trough hemoglobin on day 21 of cycle 1, 0-1; also logit_hb_increase_1g)"
  )

  covariateData <- list(
    AUC_LUSP = list(
      description = "Individual luspatercept steady-state area under the serum concentration-time curve over the 21-day dosing interval for the starting dose (AUCss). Supplied as data: this model has no PK layer.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Chen 2021 Methods: AUCss = starting dose / individual CL/F from the final population PK model (Chen_2021_luspatercept). Untransformed and uncentred; the intercept is an extrapolation to zero exposure (placebo patients were not part of the fit). Observed range in Figure 2B about 15-212 ug*day/mL.",
      source_name = "AUCss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_studies = 2L,
    n_observations = "70 binary responder records, one per patient",
    age_range = "Not reported for the exposure-response subset (PK population 18-66 years)",
    weight_range = "Not reported for the exposure-response subset (PK population 34.1-97.0 kg)",
    sex_female_pct = NA_real_,
    disease_state = "Adults with beta-thalassemia and a baseline transfusion burden < 12 RBC units/24 weeks with a trough (day 21) hemoglobin after the first dose; only hemoglobin values > 14 days after a transfusion were included.",
    dose_range = "First subcutaneous dose 0.2-1.25 mg/kg (q3w regimen).",
    regions = "Phase 2 study A536-04 (NCT01749540) and phase 3 study ACE-536-B-THAL-001 BELIEVE (NCT02604433).",
    notes = "Observed responders by AUCss quartile (Figure 2B): 1/18, 4/17, 5/17 and 10/18 (20/70 overall)."
  )

  ini({
    # Chen 2021 Figure 2B: univariate logistic regression
    #
    #   logit(p) = logit_ref + e_auc_logit * AUCss
    #
    # The panel prints OR = 3.73 and the caption defines 'OR, odds ratio for
    # 50 units of AUCss', so the slope is log(3.73)/50 = 0.02633 per
    # ug*day/mL. The intercept is not printed; it is read from the fitted
    # curve, which is drawn as a vector Bezier path in the PDF. At the
    # curve's middle node (AUCss 121.4, p = 0.4019) a curve of the printed
    # slope has intercept -3.593. Fitting both parameters to the three
    # on-curve nodes gives slope 0.02613 (0.8% from the printed OR) and
    # intercept -3.580. See the vignette source trace.
    logit_ref <- -3.593; label("Logit of the probability of a >= 1 g/dL day-21 hemoglobin increase extrapolated to zero luspatercept exposure (unitless logit)")  # digitized from Figure 2B fitted curve (not printed); see vignette
    e_auc_logit <- 0.02633; label("Log-odds of a >= 1 g/dL day-21 hemoglobin increase per 1 ug*day/mL of luspatercept AUCss (unitless logit)")  # Figure 2B: OR = 3.73 per 50 ug*day/mL; log(3.73)/50

    # No between-subject variability and no residual error in the source
    # (Bernoulli likelihood). The tiny fixed additive residual exists only so
    # rxode2 has an error model to attach to the typical-value probability.
    addSd_prob_hb_increase_1g <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_hb_increase_1g <- logit_ref + e_auc_logit * AUC_LUSP
    prob_hb_increase_1g <- expit(logit_hb_increase_1g)
    prob_hb_increase_1g ~ add(addSd_prob_hb_increase_1g)
  })
}
