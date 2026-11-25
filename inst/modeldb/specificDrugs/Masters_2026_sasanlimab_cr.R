Masters_2026_sasanlimab_cr <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-EFFICACY model for the probability",
    "of an investigator-assessed complete response (CR) as a function of the",
    "sasanlimab average serum concentration over the first 4-week dosing",
    "interval (Cavg,sd), in adults with BCG-naive high-risk non-muscle",
    "invasive bladder cancer with carcinoma in situ (CIS) at randomization",
    "treated with sasanlimab 300 mg SC every 4 weeks plus BCG (Masters 2026;",
    "CREST study B8011006 Cohort A, arms A and B, n = 179). The probability is",
    "expit(3.72 - 0.0670 * CAV), with Cavg,sd in ug/mL entering linearly and",
    "uncentred. The exposure effect was retained per the analysis plan but is",
    "not statistically significant (P = 0.15) and no covariate was retained.",
    "No PK layer: CAV is supplied as a data column (derived in the source",
    "analysis from the individual post hoc parameters of the companion",
    "population PK model, modellib('Masters_2026_sasanlimab')). No",
    "between-subject random effect and no residual error are estimated.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Coefficients are in Supporting Information Table S5.",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_cr_investigator (probability of investigator-assessed complete response, 0-1; also logit_cr_investigator)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual sasanlimab average serum concentration over the first dosing interval after the first 300 mg SC dose (Cavg,sd = AUC(0-28 days) / 28 days), predicted from the individual post hoc parameters of the companion population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "First-cycle (single-dose) metric, chosen in the source analysis following previous immuno-oncology exposure-efficacy practice (Methods, ER Efficacy). Enters LINEARLY and uncentred, so the intercept is the logit at zero exposure, far below the observed range. Across the 698 participants of Cohort A the median Cavg,sd was 21.8 ug/mL (range 5.2-38.0, mean 21.7, SD 5.63; Table S3); the CIS subset distribution is not reported separately.",
      source_name = "Cavg,sd"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 179L,
    n_studies = 1L,
    n_observations = "179 binary investigator-assessed complete-response records, one per participant",
    disease_state = "BCG-naive, high-risk non-muscle invasive bladder cancer with carcinoma in situ (CIS) at randomization; CREST trial (B8011006) Cohort A, arm A (sasanlimab + BCG induction + BCG maintenance) and arm B (sasanlimab + BCG induction only).",
    dose_range = "Sasanlimab 300 mg SC every 4 weeks, combined with BCG",
    notes = "Complete-response rate nearly 90% of 179 participants (Results, ER Efficacy). Demographics of the full Cohort A (n = 698 for EFS) are in Supporting Information Tables S1 and S2 (B8011006 column, n = 706); the CIS subset is not tabulated separately."
  )

  ini({
    # Supporting Information Table S5: logit(P) = intercept + slope * Cavg,sd
    # (Methods Equation 1 form). Odds ratios 41.4 and 0.935 per ug/mL.
    logit_ref <- 3.72; label("Logit of the probability of investigator-assessed complete response at zero Cavg,sd (unitless logit)") # Table S5: Intercept = 3.72 (95% CI 1.64, 6.09), P = 0.001
    e_cav_logit <- -0.0670; label("Log-odds of complete response per 1 ug/mL increase in Cavg,sd (unitless logit)") # Table S5: Cavg sd = -0.0670 (95% CI -0.162, 0.0225), P = 0.15

    # Bernoulli likelihood in the source; the additive residual below is a
    # placeholder so rxode2 has an error model to attach to the probability.
    addSd_prob_cr_investigator <- fixed(0.001); label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_cr_investigator <- logit_ref + e_cav_logit * CAV
    prob_cr_investigator <- expit(logit_cr_investigator)
    prob_cr_investigator ~ add(addSd_prob_cr_investigator)
  })
}
