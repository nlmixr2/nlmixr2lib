Cosson_2018_basimglurant_dizziness <- function() {
  description <- paste(
    "Logistic-regression exposure-safety model for the probability of",
    "dizziness (first occurrence, any severity) in adults with major",
    "depressive disorder taking oral modified-release basimglurant",
    "(Cosson 2018, 310 phase II patients: placebo, 0.5 mg or 1.5 mg once",
    "daily). The linear predictor is -2.970 + 0.412 * Cmax,Day1, where",
    "Cmax,Day1 (ng/mL) is the individual predicted peak concentration over",
    "the first 24-hour dosing interval, set to 0 for placebo patients. No",
    "PK layer and no ODE: exposure enters as the static per-patient column",
    "CMAX, simulated with modellib('Cosson_2018_basimglurant'). No",
    "between-subject random effect and no residual error are estimated",
    "(Bernoulli likelihood, R glm)."
  )
  reference <- paste(
    "Cosson V, Schaedeli-Stark F, Arab-Alameddine M, Chavanne C, Guerini E,",
    "Derks M, Jaeschke G, Lindemann L, Umbricht D, Santarelli L.",
    "Population Pharmacokinetic and Exposure-dizziness Modeling for a",
    "Metabotropic Glutamate Receptor Subtype 5 Negative Allosteric Modulator",
    "in Major Depressive Disorder Patients.",
    "Clin Transl Sci. 2018;11(5):523-531. doi:10.1111/cts.12566.",
    "Exposure metric produced by modellib('Cosson_2018_basimglurant').",
    sep = " "
  )
  vignette <- "Cosson_2018_basimglurant_exposure_dizziness"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the covariate CMAX)",
    concentration = "prob_dizziness (probability of dizziness, 0-1; also logit_dizziness)"
  )

  covariateData <- list(
    CMAX = list(
      description = paste(
        "Individual predicted maximum basimglurant plasma concentration over",
        "the first dosing interval (0-24 h on Day 1), Cmax,Day1."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma concentration, SINGLE-DOSE (Day 1) peak, derived in",
        "Cosson 2018 from the individual predicted concentration-time",
        "profiles of the companion population PK model; simulate it with",
        "modellib('Cosson_2018_basimglurant') as the maximum of Cc over",
        "0-24 h after the first dose. Set to 0 for placebo patients, which",
        "leaves the intercept as the placebo logit. Enters linearly,",
        "neither centred nor scaled. Observed phase II medians",
        "(min-max): 1.21 (0.42-3.26) ng/mL at 0.5 mg and 3.47",
        "(1.04-9.62) ng/mL at 1.5 mg."
      ),
      source_name = "Cmax,Day1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 310L,
    n_studies = 1L,
    n_observations = "310 binary first-occurrence dizziness records (one per patient; landmark analysis)",
    disease_state = paste(
      "adults with major depressive disorder and inadequate response to",
      "ongoing antidepressant treatment (adjunctive basimglurant, phase II",
      "trial NCT01437657, 6-week treatment)"
    ),
    dose_range = "placebo (n = 110), 0.5 mg (n = 99) or 1.5 mg (n = 101) basimglurant modified-release once daily",
    notes = paste(
      "Cosson 2018 Results 'Exposure-safety results'. Dizziness was",
      "reported by 6/110 placebo, 3/99 (0.5 mg) and 26/101 (1.5 mg)",
      "analysed patients. Patients in active arms without PK observations",
      "were excluded. Only the first occurrence was analysed and severity",
      "was ignored because most events were mild and occurred early",
      "(>50% on Day 1). Besides Cmax,Day1 no covariate was statistically",
      "significant. Logit, probit and complementary log-log links were",
      "compared; the logit link is the reported final model."
    )
  )

  ini({
    # Cosson 2018 Table 2: logit(p) = intercept + slope * Cmax,Day1.
    logit_ref <- -2.970; label("Logit of the probability of dizziness at zero exposure (placebo) (unitless logit)") # Table 2 'Intercept' -2.970 (SE 0.301)
    e_cmax_logit <- 0.412; label("Log-odds of dizziness per 1 ng/mL increase in Cmax,Day1 (unitless logit per ng/mL)") # Table 2 'Slope' 0.412 (SE 0.087); odds ratio exp(0.412) = 1.51 matches the Table 2 / Discussion odds ratio 1.510

    # Bernoulli likelihood: no residual error in the source. This tiny fixed
    # additive residual exists only so rxode2 has an error model to attach to
    # the probability.
    addSd_prob_dizziness <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_dizziness <- logit_ref + e_cmax_logit * CMAX
    prob_dizziness <- expit(logit_dizziness)

    prob_dizziness ~ add(addSd_prob_dizziness)
  })
}
