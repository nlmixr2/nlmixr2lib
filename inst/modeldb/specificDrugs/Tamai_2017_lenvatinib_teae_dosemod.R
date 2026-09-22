Tamai_2017_lenvatinib_teae_dosemod <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-safety model for the ",
    "occurrence, during cycle 1, of a treatment-emergent adverse event ",
    "leading to lenvatinib withdrawal or dose reduction in adults with ",
    "advanced hepatocellular carcinoma Child-Pugh class A (Tamai 2017; ",
    "n = 45 of the phase 2 part of study 202, NCT00946153, lenvatinib ",
    "12 mg orally once daily in 4-week cycles). The probability of the ",
    "event is expit(-4.71 + 1.82 * AUC_LEN / 1000), where AUC_LEN is the ",
    "individual lenvatinib area under the plasma concentration-time ",
    "curve over the 24 h dosing interval at steady state based on the ",
    "STARTING dose, in ng*h/mL; the slope is estimated per 1000 ng*h/mL. ",
    "This is the paper's final exposure-response model: AUC entering ",
    "linearly beat a constant-effect, a log-linear and a saturable ",
    "model, and beat steady-state minimum concentration as the exposure ",
    "metric. No covariate was retained - demographics, liver-function ",
    "markers, baseline platelet count, ECOG status, Child-Pugh class, ",
    "hepatitis B or C aetiology, portal-vein involvement, prior systemic ",
    "chemotherapy, prior antihypertensive therapy and prior surgery were ",
    "all screened and none influenced the relationship. There is no PK ",
    "layer and no ODE: exposure enters as a per-subject data column that ",
    "the companion Tamai_2017_lenvatinib.R popPK model generates as ",
    "dose / (CL/F). No between-subject random effect and no residual ",
    "error are estimated (Bernoulli likelihood). Together with a ",
    "receiver-operating-characteristic analysis this model is what ",
    "produced the paper's recommendation of an 8 mg starting dose below ",
    "60 kg and 12 mg at or above 60 kg."
  )
  reference <- paste(
    "Tamai T, Hayato S, Hojo S, Suzuki T, Okusaka T, Ikeda K, Kumada H.",
    "Dose finding of lenvatinib in subjects with advanced hepatocellular",
    "carcinoma based on population pharmacokinetic and exposure-response",
    "analyses.",
    "J Clin Pharmacol. 2017;57(9):1138-1147.",
    "doi:10.1002/jcph.917.",
    sep = " "
  )
  vignette <- "Tamai_2017_lenvatinib"
  units <- list(
    time = "n/a (static landmark exposure-response regression evaluated once per subject over cycle 1; no time dimension)",
    dosing = "n/a (no dose events; the starting dose enters only through the AUC_LEN exposure column)",
    concentration = "prob_teae_dosemod (probability of a cycle-1 treatment-emergent adverse event leading to lenvatinib withdrawal or dose reduction, 0-1; also logit_teae_dosemod)"
  )

  covariateData <- list(
    AUC_LEN = list(
      description = "Individual lenvatinib area under the plasma concentration-time curve over the 24 h dosing interval at steady state, based on the starting dose.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "NOT centred and NOT log-transformed: Tamai 2017 fits AUC as a",
        "plain linear term, so the intercept -4.71 is the logit at",
        "AUC_LEN = 0 and is an extrapolated anchor rather than any real",
        "patient's risk. The slope 1.82 is estimated PER 1000 ng*h/mL, so",
        "model() divides the column by 1000. Derived by Tamai 2017 from the",
        "final population PK model as the starting dose divided by the",
        "model-predicted individual apparent clearance; the companion",
        "Tamai_2017_lenvatinib.R reproduces exactly this quantity. Both",
        "this paper and the Majid 2024 lenvatinib models that founded the",
        "AUC_LEN canonical use the once-daily steady-state 24 h AUC, so the",
        "canonical is reused rather than extended; the one convention",
        "difference to record is that Tamai 2017 freezes the column at the",
        "STARTING dose, whereas Majid 2024 updates it step-wise as the dose",
        "level changes. Freezing it is essential to the analysis here,",
        "because the modelled event IS the first dose reduction and a",
        "dose-tracking exposure column would be contaminated by the",
        "outcome. Observed distribution in the 45-subject",
        "exposure-response population (Tamai 2017 Results): median",
        "2950 ng*h/mL (range 1560-4250) in the 21 subjects who had the",
        "event, median 2050 ng*h/mL (range 1370-3270) in those who did not.",
        "The receiver-operating-characteristic analysis identified",
        "2430 ng*h/mL as the best cutoff for predicting the high-risk group",
        "(sensitivity 0.71, specificity 0.71, area under the ROC curve",
        "0.79)."
      ),
      source_name = "AUC (area under the plasma concentration-time curve at steady state based on the starting dose)"
    )
  )

  covariatesDataExcluded <- list(
    CTROUGH = list(
      description = "Individual lenvatinib minimum plasma concentration at steady state, taken as the cycle 1 day 15 trough of the individually predicted concentration-time profile.",
      units = "ng/mL",
      type = "continuous",
      notes = paste(
        "The rival exposure metric. Tamai 2017 Methods state that both AUC",
        "based on the starting dose and the steady-state minimum",
        "concentration were tested; Results report that 'lenvatinib AUC",
        "based on the starting dose as linear function was the best",
        "predictor'. No coefficient is published for the trough-based",
        "model because it was not selected, so the correct encoding is the",
        "absence of the term, not a fixed(0) coefficient."
      )
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened as a covariate ON the exposure-response relationship and",
        "NOT retained (Tamai 2017 Results: the effects of demographics",
        "'did not influence the exposure-response relationship'). This is",
        "not a null finding about weight and toxicity - weight matters a",
        "great deal in this paper, but it acts THROUGH exposure, via the",
        "CL/F body-weight term of the companion popPK model, and so adds",
        "nothing once AUC is in the logit. Observed: median 54.3 kg in the",
        "21 subjects with early dose modification (range 42.8-78.8, 1",
        "subject not evaluable) versus 67.6 kg in those without (range",
        "48.1-85.5). The receiver-operating-characteristic analysis put the",
        "best weight cutoff at 57.8 kg (sensitivity 0.77, specificity 0.67,",
        "area under the ROC curve 0.75), which the paper rounded to the",
        "60 kg dosing threshold it recommends."
      )
    ),
    AGE = list(
      description = "Subject age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened among the demographic covariates on the exposure-response relationship and not retained (Tamai 2017 Results)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened among the demographic covariates on the exposure-response relationship and not retained (Tamai 2017 Results). Tamai 2017 does not report the sex split of the 45-subject exposure-response population."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units = "g/L",
      type = "continuous",
      notes = "Screened among the liver-function markers on the exposure-response relationship and not retained (Tamai 2017 Methods and Results)."
    ),
    INR = list(
      description = "International normalized ratio of prothrombin time at baseline.",
      units = "(ratio)",
      type = "continuous",
      notes = "Screened among the liver-function markers on the exposure-response relationship and not retained (Tamai 2017 Methods and Results)."
    ),
    PLT = list(
      description = "Platelet count at baseline.",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened on the exposure-response relationship and not retained (Tamai 2017 Methods and Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    n_observations = paste(
      "45 binary cycle-1 records, one per subject; 21 subjects (46.7",
      "percent) experienced a treatment-emergent adverse event leading to",
      "lenvatinib dose reduction or discontinuation during cycle 1. The 45",
      "subjects are the exposure-response subset of the 65 study 202",
      "subjects who contributed to the population PK analysis, and of the",
      "46 subjects of the phase 2 expansion part."
    ),
    age_range = "Not reported for the 45-subject exposure-response population; the pooled 452-subject PK analysis set spans 18.0-85.0 years.",
    weight_range = "42.8-85.5 kg (union of the with-event range 42.8-78.8 and the without-event range 48.1-85.5 reported in Tamai 2017 Results)",
    weight_median = "54.3 kg in the 21 subjects with early dose modification; 67.6 kg in those without. Study 202 overall median 58.8 kg (Tamai 2017 Discussion).",
    disease_state = paste(
      "Advanced hepatocellular carcinoma, Child-Pugh class A. Study 202 was",
      "a multicentre open-label phase 1/2 study; its phase 1 part set the",
      "maximum tolerated dose at 12 mg once daily in 4-week continuous",
      "cycles, half the dose approved for radioiodine-refractory",
      "differentiated thyroid cancer."
    ),
    dose_range = "Lenvatinib 12 mg orally once daily in 4-week cycles (the phase 2 expansion dose of study 202).",
    regions = "Study 202 was conducted in Japan and Korea.",
    hepatic_function = "Child-Pugh class A; 71 percent of the study 202 cohort had alkaline phosphatase above the upper limit of normal (Tamai 2017 Discussion).",
    notes = paste(
      "The endpoint pools every treatment-emergent adverse event leading to",
      "withdrawal or dose reduction into a single event type; Tamai 2017",
      "Limitations state plainly that 'the same exposure-response",
      "relationship was assumed for all TEAEs', and that cumulative or",
      "chronic toxicity is out of scope because the window is cycle 1 only.",
      "The motivating observation is that 74 percent of the phase 2",
      "subjects (34/46) treated at 12 mg once daily required a dose",
      "reduction to 8 mg. A parallel Kaplan-Meier analysis of time to",
      "progression stratified by AUC tertiles found NO exposure-efficacy",
      "relationship, which is why the paper is willing to lower the",
      "starting dose in low-weight subjects. Fitted in NONMEM; the",
      "receiver-operating-characteristic work used R 3.1.0 with the pROC",
      "package."
    )
  )

  ini({
    # ==================================================================
    # Tamai 2017 Results, "Exposure-Response Relationship Analysis of
    # Occurrence of TEAEs Leading to Study Drug Withdrawal or Dose
    # Reduction", final paragraph. The model is printed there in words as
    #   Logit = intercept + slope * AUC
    # and the two estimates are quoted inline with their RSE and 95%
    # confidence interval. There is no parameter table for this model;
    # that paragraph is the complete parameter set.
    #
    # Both coefficients are reproduced a second time, independently, by
    # the paper's own reported group medians: at the with-event median
    # AUC of 2950 ng*h/mL the model gives p = 0.659, at the without-event
    # median of 2050 ng*h/mL it gives p = 0.273, and at the ROC cutoff of
    # 2430 ng*h/mL it gives p = 0.429 against an overall observed event
    # rate of 46.7 percent. The vignette runs these as assertions.
    # ==================================================================

    logit_ref <- -4.71; label("Logit of the cycle-1 dose-modification probability at zero lenvatinib exposure (unitless logit)") # Tamai 2017 Results: intercept = -4.71 (RSE 29.3%; 95% CI -7.41 to -2.01)
    e_auc_len_logit <- 1.82; label("Log-odds of a cycle-1 dose modification per 1000 ng*h/mL of lenvatinib steady-state AUC (unitless logit)") # Tamai 2017 Results: slope of lenvatinib AUC effect (per 1000 ng.h/mL) = 1.82 (RSE 28.8%; 95% CI 0.793 to 2.85)

    # No between-subject variability and no residual error. The source is
    # a logistic regression with a Bernoulli likelihood, which has no
    # sigma, and Tamai 2017 estimates no random effect on the logit. The
    # tiny fixed additive residual below exists only so that rxode2 has an
    # error model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette Assumptions and deviations.
    addSd_prob_teae_dosemod <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Tamai 2017 estimates the slope per 1000 ng*h/mL, while AUC_LEN is
    # carried in the canonical ng*h/mL of the register.
    auc_per1000 <- AUC_LEN / 1000

    logit_teae_dosemod <- logit_ref + e_auc_len_logit * auc_per1000

    prob_teae_dosemod <- expit(logit_teae_dosemod)

    prob_teae_dosemod ~ add(addSd_prob_teae_dosemod)
  })
}
