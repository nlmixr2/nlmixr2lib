Walker_2023_efzofitimod_ksqlung <- function() {
  description <- paste(
    "Linear-regression exposure-response model for the King's Sarcoidosis",
    "Questionnaire-Lung (KSQ-Lung) score at Week 24 in adults with pulmonary",
    "sarcoidosis treated with efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks",
    "or placebo (Walker 2023, Phase 1b/2a, n = 32). The Week-24 score is",
    "linear in the time-averaged efzofitimod AUC from Day 1 to Week 24 (CAV,",
    "0 for placebo). The source labels the equation as a percent change from",
    "baseline, but its values are the ABSOLUTE Week-24 score (0-100 points),",
    "so the output is the absolute ksqlung; the percent change is",
    "100 * (ksqlung - KSQ_LUNG_BL) / KSQ_LUNG_BL. Static landmark model with",
    "no PK layer: exposure is supplied as data, which the source derived from",
    "the companion population PK model Walker_2023_efzofitimod."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Coefficients are from Supplementary Equation A.3."
  )
  vignette <- "Walker_2023_efzofitimod"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "ksqlung (KSQ-Lung score at Week 24, points on a 0-100 scale)"
  )

  covariateData <- list(
    CAV = list(
      description = "Time-averaged efzofitimod serum exposure over the treatment period: individual AUC from Day 1 to Week 24 divided by the length of that window in days",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Walker 2023 Methods 2.3 ('Time averaged AUC = Total AUC in study duration / Study duration') and Results 3.3 (Day 1 to Week 24 window for every endpoint except OCS). The source prints the unit as ng*day/mL, but a cumulative AUC divided by a duration in days is numerically an average concentration in ng/mL. Individual values came from post hoc estimates of the companion popPK model Walker_2023_efzofitimod; the Table 3 tertile boundaries are 805, 3,638, 6,198 and 14,860. Set to 0 for placebo. Enters linearly and uncentred.",
      source_name = "AUC (Day 1-Week 24)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "32 subject-level Week-24 records (one per efficacy-evaluable participant)",
    age_range = "mean 52.2 (SD 10.7) years in the efzofitimod arms; 53.2 (SD 10.4) on placebo",
    weight_range = "mean 98.5 (SD 26.7) kg in the efzofitimod arms; 89 (SD 14.15) kg on placebo",
    sex_female_pct = 53.1,
    disease_state = "Pulmonary sarcoidosis on chronic oral corticosteroids tapered to 5 mg/day prednisone-equivalent by Day 50",
    dose_range = "Efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks for six doses, or placebo; at least four doses received",
    regions = "United States (ATYR1923-C-002, NCT03824392)",
    notes = "Efficacy-evaluable population of Walker 2023 Table 1 (11 placebo, 21 efzofitimod); baseline KSQ-Lung score mean 50.9 (SD 13) on efzofitimod and 42.4 (SD 10.6) on placebo (Table 4)."
  )

  ini({
    ksqlung_int <- 50.1; label("Intercept of the Week-24 KSQ-Lung regression: score at zero exposure (points)") # Supplementary Equation A.3 '+ 50.1'
    e_cav_ksqlung <- 0.001935; label("Slope of Week-24 KSQ-Lung score per unit time-averaged AUC (points per ng/mL)") # Supplementary Equation A.3 '0.001935 x AUC'

    # The source does not report the residual standard error of the
    # regression. The tiny fixed additive SD exists only so rxode2 has an
    # error model to attach to the typical-value prediction.
    addSd_ksqlung <- fixed(0.001); label("Placeholder additive residual SD on the typical-value prediction; the source residual is not reported (points)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplementary Equation A.3. Read as the absolute Week-24 score: at
    # zero exposure it gives 50.1 points, an 18% rise over the placebo
    # baseline of 42.4 and consistent with the Figure 6 intercept of about
    # +15%; read literally as a percent change it would give +50%.
    ksqlung <- ksqlung_int + e_cav_ksqlung * CAV
    ksqlung ~ add(addSd_ksqlung)
  })
}
