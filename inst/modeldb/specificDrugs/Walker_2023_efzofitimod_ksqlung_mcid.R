Walker_2023_efzofitimod_ksqlung_mcid <- function() {
  description <- paste(
    "Logistic-regression exposure-response model for the probability of a",
    "King's Sarcoidosis Questionnaire-Lung (KSQ-Lung) score increase of at",
    "least 4 points at Week 24 (the validated minimal clinically important",
    "difference) in adults with pulmonary sarcoidosis treated with",
    "efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks or placebo (Walker 2023,",
    "Phase 1b/2a, n = 32). The logit is linear in the time-averaged",
    "efzofitimod AUC from Day 1 to Week 24 (CAV, 0 for placebo) and in the",
    "baseline KSQ-Lung score. Static landmark model with no PK layer:",
    "exposure is supplied as data, which the source derived from the",
    "companion population PK model Walker_2023_efzofitimod."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Coefficients are from Supplementary Equation A.5."
  )
  vignette <- "Walker_2023_efzofitimod"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_ksqlung_mcid (probability of a KSQ-Lung increase of at least 4 points at Week 24, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Time-averaged efzofitimod serum exposure over the treatment period: individual AUC from Day 1 to Week 24 divided by the length of that window in days",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Walker 2023 Methods 2.3 ('Time averaged AUC = Total AUC in study duration / Study duration') and Results 3.3 (Day 1 to Week 24 window for every endpoint except OCS). The source prints the unit as ng*day/mL, but a cumulative AUC divided by a duration in days is numerically an average concentration in ng/mL. Individual values came from post hoc estimates of the companion popPK model Walker_2023_efzofitimod; the Table 3 tertile boundaries are 805, 3,638, 6,198 and 14,860. Set to 0 for placebo. Enters linearly and uncentred.",
      source_name = "AUC (Day 1-Week 24)"
    ),
    KSQ_LUNG_BL = list(
      description = "Baseline King's Sarcoidosis Questionnaire-Lung score",
      units = "points (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Last measurement on or before the first dose (Walker 2023 Table 4 footnote). Baseline mean 50.9 (SD 13) on efzofitimod and 42.4 (SD 10.6) on placebo. Higher is better health status. Enters linearly and uncentred with a NEGATIVE coefficient: a patient who starts higher has less room to improve by 4 points.",
      source_name = "Baseline KSQL Score (point)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "32 binary Week-24 responder records (one per efficacy-evaluable participant)",
    age_range = "mean 52.2 (SD 10.7) years in the efzofitimod arms; 53.2 (SD 10.4) on placebo",
    weight_range = "mean 98.5 (SD 26.7) kg in the efzofitimod arms; 89 (SD 14.15) kg on placebo",
    sex_female_pct = 53.1,
    disease_state = "Pulmonary sarcoidosis on chronic oral corticosteroids tapered to 5 mg/day prednisone-equivalent by Day 50",
    dose_range = "Efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks for six doses, or placebo; at least four doses received",
    regions = "United States (ATYR1923-C-002, NCT03824392)",
    notes = "Efficacy-evaluable population of Walker 2023 Table 1 (11 placebo, 21 efzofitimod). Seven of the 11 placebo participants met the 4-point MCID (Results 3.3), which drives the large intercept."
  )

  ini({
    logite0 <- 6.38; label("Logit intercept of the KSQ-Lung-responder probability (unitless logit)") # Supplementary Equation A.5 '6.38'
    e_cav_logit <- 0.000219; label("Log-odds of KSQ-Lung response per unit time-averaged AUC (logit per ng/mL)") # Supplementary Equation A.5 '0.000219 x AUC'
    e_ksq_lung_bl_logit <- -0.15; label("Log-odds of KSQ-Lung response per baseline KSQ-Lung point (logit per point)") # Supplementary Equation A.5 '-0.15 x Baseline KSQL Score (point)'

    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value probability.
    addSd_prob_ksqlung_mcid <- fixed(0.001); label("Placeholder additive residual SD on the typical-value probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplementary Equation A.5. Both predictors are uncentred, so
    # logite0 is the logit at zero exposure and a baseline score of 0.
    logit_ksqlung_mcid <- logite0 + e_cav_logit * CAV + e_ksq_lung_bl_logit * KSQ_LUNG_BL
    prob_ksqlung_mcid <- expit(logit_ksqlung_mcid)
    prob_ksqlung_mcid ~ add(addSd_prob_ksqlung_mcid)
  })
}
