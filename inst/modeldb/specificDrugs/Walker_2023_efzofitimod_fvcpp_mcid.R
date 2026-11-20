Walker_2023_efzofitimod_fvcpp_mcid <- function() {
  description <- paste(
    "Logistic-regression exposure-response model for the probability of a",
    "ppFVC increase of at least 2.5% at Week 24 (the minimal clinically",
    "important difference adopted from systemic-sclerosis-associated ILD) in",
    "adults with pulmonary sarcoidosis treated with efzofitimod 1, 3 or",
    "5 mg/kg IV every 4 weeks or placebo (Walker 2023, Phase 1b/2a, n = 32).",
    "The logit is linear in the time-averaged efzofitimod AUC from Day 1 to",
    "Week 24 (CAV, 0 for placebo) and in the baseline ppFVC. Static landmark",
    "model with no PK layer: exposure is supplied as data, which the source",
    "derived from the companion population PK model Walker_2023_efzofitimod."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Coefficients are from Supplementary Equation A.4."
  )
  vignette <- "Walker_2023_efzofitimod"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_fvcpp_incr2p5 (probability of a ppFVC increase of at least 2.5% at Week 24, 0-1)"
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
    FVC_PCTPRED = list(
      description = "Baseline forced vital capacity, percent of predicted",
      units = "% predicted",
      type = "continuous",
      reference_category = NULL,
      notes = "Last measurement on or before the first dose (Walker 2023 Table 4 footnote). Baseline mean 78.1 (SD 14.9) in the efzofitimod arms and 78.1 (SD 11.7) on placebo. Enters linearly and uncentred.",
      source_name = "Baseline ppFVC (%)"
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
    notes = "Efficacy-evaluable population of Walker 2023 Table 1 (11 placebo, 21 efzofitimod). The 2.5% threshold is a post hoc definition taken from systemic-sclerosis-associated ILD (Khanna 2009); Walker 2023 Methods 2.4 notes there is no validated ppFVC MCID for sarcoidosis."
  )

  ini({
    logite0 <- -10.25; label("Logit intercept of the ppFVC-responder probability (unitless logit)") # Supplementary Equation A.4 '-10.25'
    e_cav_logit <- 0.0000953; label("Log-odds of ppFVC response per unit time-averaged AUC (logit per ng/mL)") # Supplementary Equation A.4 '0.0000953 x AUC'
    e_fvc_pctpred_logit <- 0.11; label("Log-odds of ppFVC response per unit baseline ppFVC (logit per % predicted)") # Supplementary Equation A.4 '0.11 x Baseline ppFVC (%)'

    # The source is a binomial logistic regression (Bernoulli likelihood,
    # no random effects). The tiny additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value probability.
    addSd_prob_fvcpp_incr2p5 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplementary Equation A.4. Both predictors are uncentred, so
    # logite0 is the logit at zero exposure and zero baseline ppFVC, far
    # outside the data.
    logit_fvcpp_incr2p5 <- logite0 + e_cav_logit * CAV + e_fvc_pctpred_logit * FVC_PCTPRED
    prob_fvcpp_incr2p5 <- expit(logit_fvcpp_incr2p5)
    prob_fvcpp_incr2p5 ~ add(addSd_prob_fvcpp_incr2p5)
  })
}
