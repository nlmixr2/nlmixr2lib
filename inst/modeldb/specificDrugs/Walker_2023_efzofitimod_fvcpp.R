Walker_2023_efzofitimod_fvcpp <- function() {
  description <- paste(
    "Linear-regression exposure-response model for percent-predicted forced",
    "vital capacity (ppFVC) at Week 24 in adults with pulmonary sarcoidosis",
    "treated with efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks or placebo",
    "(Walker 2023, Phase 1b/2a, n = 32). The Week-24 ppFVC is linear in the",
    "time-averaged efzofitimod AUC from Day 1 to Week 24 (CAV, 0 for placebo)",
    "and in the baseline ppFVC. The source labels the equation as a percent",
    "change from baseline, but its values are the ABSOLUTE Week-24 ppFVC, so",
    "the output is the absolute fvcpp; the percent change is",
    "100 * (fvcpp - FVC_PCTPRED) / FVC_PCTPRED. Static landmark model with",
    "no PK layer: exposure is supplied as data, which the source derived from",
    "the companion population PK model Walker_2023_efzofitimod."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Coefficients are from Supplementary Equation A.2."
  )
  vignette <- "Walker_2023_efzofitimod"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "fvcpp (FVC percent predicted at Week 24, %)"
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
      source_name = "baseline ppFVC"
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
    notes = "Efficacy-evaluable population of Walker 2023 Table 1 (11 placebo, 21 efzofitimod); demographics and baseline disease characteristics in Tables 3 and 4."
  )

  ini({
    fvcpp_int <- -17.15; label("Intercept of the Week-24 ppFVC regression (% predicted)") # Supplementary Equation A.2 '-17.15'
    e_cav_fvcpp <- 0.0005061; label("Slope of Week-24 ppFVC per unit time-averaged AUC (% predicted per ng/mL)") # Supplementary Equation A.2 '0.0005061 x AUC'
    e_fvc_pctpred_fvcpp <- 1.182; label("Slope of Week-24 ppFVC per unit baseline ppFVC (unitless)") # Supplementary Equation A.2 '1.182 x baseline ppFVC'

    # The source does not report the residual standard error of the
    # regression. The tiny fixed additive SD exists only so rxode2 has an
    # error model to attach to the typical-value prediction.
    addSd_fvcpp <- fixed(0.001); label("Placeholder additive residual SD on the typical-value prediction; the source residual is not reported (% predicted)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplementary Equation A.2. Read as the absolute Week-24 ppFVC: at
    # the mean baseline of 78.1 and zero exposure it gives 75.2, a 3.7%
    # fall, consistent with the Figure 5 placebo intercept; read literally
    # as a percent change it would give +75%, off the Figure 5 scale.
    fvcpp <- fvcpp_int + e_cav_fvcpp * CAV + e_fvc_pctpred_fvcpp * FVC_PCTPRED
    fvcpp ~ add(addSd_fvcpp)
  })
}
