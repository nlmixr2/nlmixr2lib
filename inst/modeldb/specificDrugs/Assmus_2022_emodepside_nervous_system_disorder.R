Assmus_2022_emodepside_nervous_system_disorder <- function() {
  description <- paste(
    "Binary logistic exposure-safety regression for the probability of",
    "a drug-related treatment-emergent nervous system disorder (MedDRA system organ class; mostly dizziness and headache)",
    "in healthy male volunteers given oral emodepside in three phase I",
    "studies (Assmus 2022, n = 142). The only predictor is the individual",
    "maximum plasma emodepside concentration (CMAX, ng/mL), entering",
    "linearly and uncentred: logit(p) = -3.2 + 0.0063 * CMAX",
    "(S6 Table). There is no PK layer and no ODE: CMAX is an individual",
    "prediction from the companion population PK model",
    "modellib('Assmus_2022_emodepside'). No random effect and no residual",
    "error are estimated (Bernoulli likelihood, fitted in R)."
  )
  reference <- paste(
    "Assmus F, Hoglund RM, Monnot F, Specht S, Scandale I, Tarning J.",
    "Drug development for the treatment of onchocerciasis: Population",
    "pharmacokinetic and adverse events modeling of emodepside.",
    "PLoS Negl Trop Dis. 2022;16(3):e0010219.",
    "doi:10.1371/journal.pntd.0010219. PMCID: PMC8912909.",
    "Coefficients from S6 Table (supporting information).",
    "Exposure metric produced by modellib('Assmus_2022_emodepside')."
  )
  vignette <- "Assmus_2022_emodepside"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the covariate CMAX)",
    concentration = "prob_nervous_system_disorder (probability of drug-related treatment-emergent nervous system disorder, 0-1; also logit_nervous_system_disorder)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual maximum plasma emodepside concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL venous plasma concentration, derived per subject from the",
        "final population PK model (Methods, 'Exposure-adverse events",
        "analysis'): the individual peak over the dosing regimen received",
        "(single dose in the single ascending dose and relative",
        "bioavailability studies; up to 10 days of once- or twice-daily",
        "dosing in the multiple ascending dose study). Simulate it with",
        "modellib('Assmus_2022_emodepside') and take the per-subject",
        "maximum of Cc. Enters linearly, neither centred nor scaled, so",
        "CMAX = 0 returns the intercept probability. Observed events",
        "occurred 20 min to 5.5 h after dosing, mostly shortly after Tmax."
      ),
      source_name = "Cmax"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 142L,
    n_studies = 3L,
    n_observations = "142 binary per-subject outcomes; 18 of 142 subjects (12.7%) had the event",
    age_range = "18-54 years",
    age_median = "32 years",
    weight_range = "53.2-105 kg",
    weight_median = "79.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "healthy male volunteers",
    dose_range = "oral emodepside 1-40 mg single doses and 5-10 mg once or twice daily for 10 days",
    regions = "United Kingdom",
    notes = paste(
      "Same 142 subjects as the population PK analysis (Assmus 2022",
      "Tables 1 and 2); placebo recipients were not included. All",
      "drug-related TEAEs were transient and mild, apart from two moderate",
      "headaches at 40 mg. Demographic covariates (age, BMI, smoking,",
      "drinking) did not improve the model. S6 Table diagnostics:",
      "AIC 90.0, ROC area 77.6%."
    )
  )

  ini({
    # S6 Table, column 'Nervous system disorder' of the Cmax block: binary logistic
    # regression of the event on the individual Cmax, logit(p) = b0 + b1 * Cmax.
    # Check: Results: an increase in Cmax from 300 to 400 ng/mL raises the probability from 21.1% to 33.3%, and 500 ng/mL gives 48%.
    logit_ref <- -3.2
    label("Logit of the probability of drug-related treatment-emergent nervous system disorder at a Cmax of 0 ng/mL (unitless logit)") # S6 Table 'Intercept (StError)' -3.20 (0.45)
    e_cmax_nervous_system_disorder <- 0.0063
    label("Log-odds of drug-related treatment-emergent nervous system disorder per 1 ng/mL increase in Cmax (unitless logit per ng/mL)") # S6 Table 'LogOdds (StError)' 0.0063 (0.0014); odds increase 0.63 (95% CI 0.35-0.91)% per ng/mL

    # No random effect and no residual error: Bernoulli likelihood. The tiny
    # fixed additive residual below exists only so rxode2 has an error model
    # to attach to the probability; it is not a published quantity.
    addSd_prob_nervous_system_disorder <- fixed(0.001)
    label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # CMAX must be supplied in ng/mL; see covariateData[['CMAX']]$notes.
    logit_nervous_system_disorder <- logit_ref + e_cmax_nervous_system_disorder * CMAX
    prob_nervous_system_disorder <- expit(logit_nervous_system_disorder)

    prob_nervous_system_disorder ~ add(addSd_prob_nervous_system_disorder)
  })
}
