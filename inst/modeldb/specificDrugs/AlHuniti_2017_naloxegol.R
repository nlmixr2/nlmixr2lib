AlHuniti_2017_naloxegol <- function() {
  description <- "Dose-response model of weekly spontaneous-bowel-movement (SBM) response and study dropout for oral naloxegol in adults with opioid-induced constipation and chronic noncancer pain (Al-Huniti 2017, phase II study 07-IN-NX003). A week counts as a response when it has at least 3 SBMs and at least 1 SBM more than baseline; its probability is a mixed-effects logistic regression that is linear in the daily naloxegol dose, logit(P) = E0 + alpha * DOSE + eta, with no time trend over the 4 treatment weeks. Time to study discontinuation follows a Weibull distribution parameterised by its median, which falls linearly with dose (median = B + beta * DOSE days). The two models were fitted jointly in NONMEM (Laplace). There is no pharmacokinetic layer: the nominal daily dose is the exposure metric. The linear median is positive only below about 92 mg/day; the phase II data cover 0 to 50 mg/day."
  reference <- paste(
    "Al-Huniti N, Xu H, Zhou D, Aksenov S, Fox R, Bui KH.",
    "Population Exposure-Response Modeling Supported Selection of Naloxegol",
    "Doses in Phase III Studies in Patients With Opioid-Induced Constipation.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(10):705-711.",
    "doi:10.1002/psp4.12229"
  )
  vignette <- "AlHuniti_2017_naloxegol"
  units <- list(
    time = "day",
    dosing = "mg (daily oral naloxegol dose, supplied as the DOSE covariate; no dosing events)",
    concentration = "(no concentration output; prob_sbm_response and sur are probabilities, 0-1)"
  )

  covariateData <- list(
    DOSE = list(
      description = "Randomized naloxegol daily dose (nominal arm dose)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a) of the DOSE canonical: the randomized once-daily dose, constant over the 4-week double-blind period (0 for placebo; 5, 25 or 50 mg/day in the phase II study). Enters linearly on the logit of the weekly response probability and on the median time to dropout. Set DOSE = 0 on any record before the first dose (the placebo run-in week), where the published visual predictive check shows the placebo-level response in every arm.",
      source_name = "Dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 185L,
    n_studies = 1L,
    age_range = "Adults; mean age 49.7 years",
    weight_range = "50.2-182.1 kg",
    sex_female_pct = 62.2,
    race_ethnicity = "86.5% White, 11.4% Black",
    disease_state = "Opioid-induced constipation in patients taking opioids for chronic noncancer pain.",
    dose_range = "Naloxegol 5, 25 or 50 mg orally once daily, or placebo, for 4 weeks after a 1-week single-blind placebo run-in (three sequential dose cohorts, each randomized 1:1 to naloxegol or placebo).",
    regions = "Not reported in this paper.",
    notes = "Modified intention-to-treat population of the randomized, double-blind, placebo-controlled, dose-escalation phase II study 07-IN-NX003 (Table 1: 95 placebo, 31 at 5 mg, 29 at 25 mg, 30 at 50 mg). The planned 100 mg cohort was not run after the 50 mg cohort's safety review. Demographics from the Results text."
  )

  ini({
    # Table 2 'Final model parameter estimates' (joint NONMEM 7.2 Laplace fit).
    logite0 <- -0.51
    label("Logit of the weekly SBM-response probability on placebo, E0 (unitless logit)") # Table 2: E0 = -0.51 (RSE 34%); Results: placebo weekly response probability 0.38
    e_dose_logit <- 0.061
    label("Increase in the logit of the weekly SBM-response probability per mg of daily naloxegol dose, alpha (1/mg)") # Table 2: alpha = 0.061 (RSE 15%); Results: odds ratio about 1.8 per 10 mg
    ltmed_drop <- log(110)
    label("Log of the median time to study discontinuation on placebo, B (log day)") # Table 2: B = 110 days (RSE 28%)
    e_dose_tmed_drop <- -1.2
    label("Change in the median time to study discontinuation per mg of daily naloxegol dose, beta (day/mg)") # Table 2: beta = -1.2 (RSE 57%); Results: median falls about 12 days per 10 mg
    lalfa_drop <- log(1.4)
    label("Log of the Weibull shape factor of the time to study discontinuation, gamma (unitless)") # Table 2: gamma = 1.4 (RSE 12%)

    etalogite0 ~ 2.6 # Table 2: interindividual variance of E0 = 2.6, additive on the logit

    # The source likelihoods are Bernoulli (weekly response) and the Weibull
    # event-time density (dropout); neither has a residual error. These fixed
    # placeholders only let nlmixr2 accept the two observation declarations.
    addSd_prob_sbm_response <- fixed(0.001)
    label("Placeholder additive residual SD on the weekly SBM-response probability (unitless)") # not from the source; see vignette Assumptions and deviations
    addSd_sur <- fixed(0.001)
    label("Placeholder additive residual SD on the probability of remaining in the study (unitless)") # not from the source; see vignette Assumptions and deviations
  })

  model({
    # Weekly SBM response: logit(P) = E0 + alpha * Dose + eta (Methods).
    prob_sbm_response <- expit(logite0 + e_dose_logit * DOSE + etalogite0)

    # Dropout: Weibull with median tmed_drop = B + beta * Dose (Methods),
    # p(t) = ln2 * gamma / tmed * (t / tmed)^(gamma - 1) * exp(-ln2 * (t / tmed)^gamma),
    # so the survivor function is exp(-ln2 * (t / tmed)^gamma). The cumulative
    # hazard is written in closed form (exact for a dose held constant from
    # t = 0) instead of integrating the hazard from t = 0.
    tmed_drop <- exp(ltmed_drop) + e_dose_tmed_drop * DOSE
    alfa_drop <- exp(lalfa_drop)
    hazard_drop <- log(2) * alfa_drop / tmed_drop * (t / tmed_drop)^(alfa_drop - 1)
    cumhaz_drop <- log(2) * (t / tmed_drop)^alfa_drop
    sur <- exp(-cumhaz_drop)

    prob_sbm_response ~ add(addSd_prob_sbm_response)
    sur ~ add(addSd_sur)
  })
}
