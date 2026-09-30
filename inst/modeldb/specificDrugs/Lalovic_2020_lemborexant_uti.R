Lalovic_2020_lemborexant_uti <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-safety model for the probability of ",
    "reporting urinary tract infection as a treatment-emergent adverse event during lemborexant ",
    "treatment, as a linear function of the individual average steady-state ",
    "lemborexant concentration (Cav,ss) on the logit scale, in adults and elderly ",
    "subjects with insomnia disorder or irregular sleep-wake rhythm disorder from ",
    "studies 202, 303 (SUNRISE 2) and 304 (SUNRISE 1) (Lalovic 2020). The paper ",
    "prints no regression coefficients; it tabulates the fitted probability for ONE ",
    "reference subject, a 75-year-old White woman, at ten Cav,ss values and for ",
    "placebo (Table 4). The intercept and exposure slope here were back-solved by ",
    "the maintainers from those ten active-treatment rows, and the placebo logit is ",
    "the Table 4 placebo row, so the model reproduces Table 4 for that reference ",
    "subject only: the age, sex and race coefficients of the source model are not ",
    "reported and cannot be recovered. Cav,ss was not a statistically significant predictor of urinary tract infection (Results). Cav,ss enters as the CAV ",
    "covariate column (derive it from the companion population PK model ",
    "Lalovic_2020_lemborexant as daily dose / (24 h x CL/F)). No random effects ",
    "and no residual error are estimated (Bernoulli likelihood)."
  )
  reference <- paste(
    "Lalovic B, Majid O, Aluri J, Landry I, Moline M, Hussein Z. (2020).",
    "Population Pharmacokinetics and Exposure-Response Analyses for the Most",
    "Frequent Adverse Events Following Treatment With Lemborexant, an Orexin",
    "Receptor Antagonist, in Subjects With Insomnia Disorder.",
    "Journal of Clinical Pharmacology 60(12):1642-1654.",
    "doi:10.1002/jcph.1683.",
    sep = " "
  )
  vignette <- "Lalovic_2020_lemborexant"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_uti (probability of treatment-emergent urinary tract infection, 0-1; also logit_uti)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state lemborexant plasma concentration (Cav,ss)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Derived in the source from the individual population PK estimates",
        "(Lalovic 2020 Methods and Results); with the companion model",
        "Lalovic_2020_lemborexant, Cav,ss = daily dose / (24 h x CL/F). Table 4",
        "prints Cav,ss quantiles of 6.3-19 ng/mL (5th-95th) at 5 mg and 12-37",
        "ng/mL at 10 mg. The back-solved slope is only supported over that",
        "active-treatment range. For placebo subjects set PLACEBO = 1; CAV is then",
        "ignored."
      ),
      source_name = "Cav,ss"
    ),
    PLACEBO = list(
      description = "Placebo-arm indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (lemborexant-treated)",
      notes = paste(
        "Selects the Table 4 placebo-row probability. The placebo row does not lie",
        "on the straight logit line through the ten active-treatment rows",
        "extrapolated to Cav,ss = 0, so it is carried as its own logit rather than",
        "as CAV = 0. Its value (2.5%) equals the pooled observed incidence",
        "quoted in the Results -- which labels that figure 'active treatment',",
        "apparently swapped with the placebo value -- so it is most likely the",
        "observed placebo proportion rather than a covariate-adjusted prediction."
      ),
      source_name = "Placebo"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1664,
    n_studies = 3,
    age_range = "Adults and elderly; reference subject for this model 75 years",
    disease_state = paste(
      "Insomnia disorder (SUNRISE 1, study 304, n = 524; SUNRISE 2, study 303,",
      "n = 726) and irregular sleep-wake rhythm disorder (study 202, n = 62)"
    ),
    dose_range = "Placebo or lemborexant 5-15 mg nightly (Table 4 tabulates 5 and 10 mg)",
    notes = paste(
      "Exposure-response data set of 1664 subjects from studies 202, 303 and 304",
      "(Abstract). 44 subjects reported urinary tract infection (Results). The regression",
      "adjusted for age, sex and race (White, Black/African American, Japanese,",
      "others) with interaction terms considered, but those coefficients are not",
      "reported; this file encodes the 75-year-old White woman of Table 4."
    )
  )

  ini({
    # Lalovic 2020 Table 4, column 'Urinary Tract Infection'.
    #
    #   logit(p) = logit_active + e_cav_logit * CAV     (lemborexant-treated)
    #   logit(p) = logit_placebo                        (placebo)
    #
    # The source prints predicted probabilities, not coefficients. The two
    # active-treatment values below were back-solved by the maintainers as the
    # ordinary least-squares line of qlogis(Table 4 probability / 100) on the
    # printed Cav,ss over the ten 5 mg and 10 mg rows; the line reproduces every
    # row to within 0.4 percentage points, the rounding of the printed Cav,ss.
    # The validation vignette re-derives both values from the table.
    logit_active <- -4.0476; label("Logit of the urinary tract infection probability at Cav,ss = 0 on the active-treatment line, 75-year-old White woman (unitless logit)") # back-solved by the maintainers from Table 4 'Urinary Tract Infection', 5 mg and 10 mg rows (not a printed value)
    e_cav_logit <- 0.02590; label("Change in the logit of the urinary tract infection probability per ng/mL of Cav,ss (1/(ng/mL))") # back-solved by the maintainers from Table 4 'Urinary Tract Infection', 5 mg and 10 mg rows (not a printed value)
    logit_placebo <- -3.6636; label("Logit of the urinary tract infection probability on placebo, 75-year-old White woman (unitless logit)") # Table 4 'Urinary Tract Infection', Placebo row = 2.5%; qlogis(2.5 / 100)

    # The source is a logistic regression (Bernoulli likelihood, no sigma, no
    # random effects). This tiny fixed additive residual is NOT a published
    # quantity; it only gives rxode2 an error model to attach to the
    # probability output.
    addSd_prob_uti <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see the vignette Assumptions and deviations
  })

  model({
    logit_uti <- (1 - PLACEBO) * (logit_active + e_cav_logit * CAV) +
      PLACEBO * logit_placebo
    prob_uti <- expit(logit_uti)
    prob_uti ~ add(addSd_prob_uti)
  })
}
