Walker_2023_efzofitimod_ocs <- function() {
  description <- paste(
    "Linear-regression exposure-response model for the percent change from",
    "baseline in mean daily oral corticosteroid (prednisone-equivalent) dose",
    "over the post-taper period (Day 51 to Week 24) in adults with pulmonary",
    "sarcoidosis treated with efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks or",
    "placebo (Walker 2023, Phase 1b/2a, n = 32). The response is linear in",
    "the time-averaged efzofitimod AUC over the same post-taper window (CAV,",
    "0 for placebo); no baseline covariate was retained. Static landmark",
    "model with no PK layer: exposure is supplied as data, which the source",
    "derived from the companion population PK model Walker_2023_efzofitimod."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Coefficients are from Supplementary Equation A.1."
  )
  vignette <- "Walker_2023_efzofitimod"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "ocsdosecfb (percent change from baseline in mean daily prednisone-equivalent dose, %)"
  )

  covariateData <- list(
    CAV = list(
      description = "Time-averaged efzofitimod serum exposure over the post-taper period: individual AUC from Day 51 to Week 24 divided by the length of that window in days",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Walker 2023 Methods 2.3 ('Time averaged AUC = Total AUC in study duration / Study duration') and 2.5 / Results 3.3 (Day 51 to Week 24 window for the OCS endpoint). The source prints the unit as ng*day/mL, but a cumulative AUC divided by a duration in days is numerically an average concentration in ng/mL. Individual values came from post hoc estimates of the companion popPK model Walker_2023_efzofitimod; Figure 4 shows an observed range of about 600-14,900. Set to 0 for placebo. Enters linearly and uncentred.",
      source_name = "AUC (Day 51-Week 24)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "32 subject-level post-taper records (one per efficacy-evaluable participant)",
    age_range = "mean 52.2 (SD 10.7) years in the efzofitimod arms; 53.2 (SD 10.4) on placebo",
    weight_range = "mean 98.5 (SD 26.7) kg in the efzofitimod arms; 89 (SD 14.15) kg on placebo",
    sex_female_pct = 53.1,
    disease_state = "Pulmonary sarcoidosis on chronic oral corticosteroids (10-25 mg/day prednisone or equivalent at entry), with a protocol-guided taper to 5 mg/day between Day 15 and Day 50",
    dose_range = "Efzofitimod 1, 3 or 5 mg/kg IV every 4 weeks for six doses, or placebo; at least four doses received",
    regions = "United States (ATYR1923-C-002, NCT03824392)",
    notes = "Efficacy-evaluable population of Walker 2023 Table 1 (11 placebo, 21 efzofitimod); demographics in Tables 3 and 4 (baseline prednisone-equivalent dose mean 13.1 mg/day on efzofitimod, 13.6 mg/day on placebo). Sex percentage is the efficacy-evaluable split (17 of 32 female)."
  )

  ini({
    ocsdosecfb_int <- -43.01; label("Intercept: percent change from baseline in post-taper mean daily OCS dose at zero exposure (%)") # Supplementary Equation A.1 '- 43.01'
    e_cav_ocsdosecfb <- -0.001731; label("Slope of percent change in post-taper mean daily OCS dose per unit time-averaged AUC (% per ng/mL)") # Supplementary Equation A.1 '-0.001731 x AUC'

    # The source does not report the residual standard error of the
    # regression. The tiny fixed additive SD exists only so rxode2 has an
    # error model to attach to the typical-value prediction.
    addSd_ocsdosecfb <- fixed(0.001); label("Placeholder additive residual SD on the typical-value prediction; the source residual is not reported (%)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplementary Equation A.1 (labelled 'drug-adjusted AUC percent change
    # from baseline'). Reproduces the fitted line of Figure 4.
    ocsdosecfb <- ocsdosecfb_int + e_cav_ocsdosecfb * CAV
    ocsdosecfb ~ add(addSd_ocsdosecfb)
  })
}
