Masters_2026_sasanlimab_irae_grade3 <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for the probability",
    "of a grade >= 3 immune-related adverse event (irAE) as a function of the",
    "sasanlimab maximum serum concentration after the first dose (Cmax,sd,",
    "linear and uncentred), in adults with non-muscle invasive bladder cancer",
    "and other solid tumors (Masters 2026; pooled studies B8011001, B8011006",
    "and B8011007, n = 1007). The exposure slope 0.00713 per ug/mL is the",
    "published estimate (not significant, P = 0.27). The intercept is NOT",
    "published: it was digitised from the fitted curve of Supporting Figure",
    "S2, so the model reproduces that plotted curve. The final published",
    "model also carries covariate terms that are not reported; the plotted",
    "curve (about 0.035 at the median exposure) sits near the observed",
    "incidence of the non-NMIBC monotherapy studies (4.7%, B8011001 +",
    "B8011007) and far UNDER the 19.7% incidence of patients with NMIBC",
    "treated with sasanlimab + BCG (study B8011006). Use it for non-NMIBC",
    "monotherapy only. No PK layer: CMAX is a data column. No random effect",
    "and no residual error are estimated.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Exposure slope from Table 4; intercept digitised from Supporting",
    "Information Figure S2 (panel 'irAE Grade>=3').",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_irae_grade3 (probability of a grade >= 3 immune-related adverse event, 0-1; also logit_irae_grade3)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual sasanlimab maximum serum concentration after the FIRST dose (Cmax,sd), predicted from the individual post hoc parameters of the companion population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL serum concentration; SINGLE-DOSE Cmax (first dose of the participant's regimen, IV 0.5-10 mg/kg or SC 300 / 600 mg); empirical-Bayes prediction from modellib('Masters_2026_sasanlimab'). Enters LINEARLY and uncentred in ug/mL (Table 4 row 'Cmax sd (ug/mL)'; Figure S2 axis 'Cmax sd (ug/mL)', 2-254). Observed 5th / 25th / 50th / 75th / 95th percentiles digitised from the Figure S2 reference lines: about 15.2 / 21.8 / 27.4 / 33.5 / 55.1 ug/mL.",
      source_name = "Cmax sd"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1007L,
    n_studies = 3L,
    n_events = 153L,
    n_observations = "1007 binary records, one per participant",
    age_range = "26-90 years",
    age_median = "66 years",
    weight_range = "38.2-173 kg",
    weight_median = "74 kg",
    sex_female_pct = 21.8,
    disease_state = "BCG-naive high-risk non-muscle invasive bladder cancer (70.1%, study B8011006 CREST Cohort A, sasanlimab + BCG), non-small cell lung cancer (19.7%), urothelial carcinoma (3.9%) and other advanced solid tumors (6.4%).",
    dose_range = "IV 0.5, 1, 3 or 10 mg/kg every 3 weeks (B8011001 dose escalation); SC 300 mg every 4 weeks (all three studies) and SC 600 mg every 6 weeks (B8011007).",
    notes = "Same 1007 participants as the population PK analysis (Results, ER Safety). Grade >= 3 irAE incidence 15.2% overall; 5.5% in B8011001 (n = 146), 19.7% in B8011006 (n = 706) and 3.9% in B8011007 (n = 155) (Supporting Table S6). A logistic model with an intercept averages to the observed incidence, so a curve at 0.035 cannot be the whole-population fit: the final model must hold a covariate (most likely NMIBC tumor type, study B8011006 or BCG combination) whose coefficient is not printed."
  )

  ini({
    # logit(P) = intercept + slope * Cmax,sd (Methods Equation 1 form).
    # Digitised by the maintainers from Figure S2 (panel 'irAE Grade>=3',
    # 500 curve points). The curve lies at probabilities 0.026-0.17, where a
    # single pixel is up to 0.1 on the logit scale; the free-slope fit gives
    # 0.0079 (printed 0.00713). With the slope held at 0.00713 the intercept
    # over the whole plotted curve is -3.51 (logit residual SD 0.06; -3.59
    # over the 5th-95th percentile exposure range alone).
    logit_ref <- -3.51; label("Logit of the probability of grade >= 3 irAE at Cmax,sd = 0 ug/mL, at the Figure S2 reference stratum (unitless logit)") # figure-derived: digitised from Supporting Figure S2 with the slope held at the Table 4 value; not printed in the paper
    e_cmax_logit <- 0.00713; label("Log-odds of grade >= 3 irAE per 1 ug/mL increase in Cmax,sd (unitless logit)") # Table 4: Cmax sd (ug/mL) = 0.00713 (95% CI -0.00762, 0.0185), P = 0.27, odds ratio 1.01

    # Bernoulli likelihood in the source; the additive residual below is a
    # placeholder so rxode2 has an error model to attach to the probability.
    addSd_prob_irae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_irae_grade3 <- logit_ref + e_cmax_logit * CMAX
    prob_irae_grade3 <- expit(logit_irae_grade3)
    prob_irae_grade3 ~ add(addSd_prob_irae_grade3)
  })
}
