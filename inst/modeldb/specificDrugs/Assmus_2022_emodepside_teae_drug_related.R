Assmus_2022_emodepside_teae_drug_related <- function() {
  description <- paste(
    "Binary logistic exposure-safety regression for the probability of",
    "any drug-related treatment-emergent adverse event (TEAE), of any system organ class",
    "in healthy male volunteers given oral emodepside in three phase I",
    "studies (Assmus 2022, n = 142). The only predictor is the individual",
    "maximum plasma emodepside concentration (CMAX, ng/mL), entering",
    "linearly and uncentred: logit(p) = -2.38 + 0.0064 * CMAX",
    "(S7 Table). There is no PK layer and no ODE: CMAX is an individual",
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
    "Coefficients from S7 Table (supporting information).",
    "Exposure metric produced by modellib('Assmus_2022_emodepside')."
  )
  vignette <- "Assmus_2022_emodepside"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the covariate CMAX)",
    concentration = "prob_teae_drug_related (probability of any drug-related TEAE, 0-1; also logit_teae_drug_related)"
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
    n_observations = "142 binary per-subject outcomes; 31 of 142 subjects (21.8%) had the event",
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
      "drinking) did not improve the model. S7 Table diagnostics:",
      "AIC 125.8, ROC area 75.2%."
    )
  )

  ini({
    # S7 Table, column 'Cmax' of the Cmax block: binary logistic
    # regression of the event on the individual Cmax, logit(p) = b0 + b1 * Cmax.
    # Check: S7 Table only; the paper prints no worked probabilities for this endpoint.
    logit_ref <- -2.38
    label("Logit of the probability of any drug-related TEAE at a Cmax of 0 ng/mL (unitless logit)") # S7 Table 'Intercept (StError)' -2.38 (0.34)
    e_cmax_teae_drug_related <- 0.0064
    label("Log-odds of any drug-related TEAE per 1 ng/mL increase in Cmax (unitless logit per ng/mL)") # S7 Table 'LogOdds (StError)' 0.0064 (0.0014); odds increase 0.64 (95% CI 0.37-0.92)% per ng/mL

    # No random effect and no residual error: Bernoulli likelihood. The tiny
    # fixed additive residual below exists only so rxode2 has an error model
    # to attach to the probability; it is not a published quantity.
    addSd_prob_teae_drug_related <- fixed(0.001)
    label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # CMAX must be supplied in ng/mL; see covariateData[['CMAX']]$notes.
    logit_teae_drug_related <- logit_ref + e_cmax_teae_drug_related * CMAX
    prob_teae_drug_related <- expit(logit_teae_drug_related)

    prob_teae_drug_related ~ add(addSd_prob_teae_drug_related)
  })
}
