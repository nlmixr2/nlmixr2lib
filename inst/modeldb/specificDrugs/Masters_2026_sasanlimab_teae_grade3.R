Masters_2026_sasanlimab_teae_grade3 <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for the probability",
    "of any grade >= 3 treatment-emergent adverse event (TEAE) as a function",
    "of the natural log of the sasanlimab maximum serum concentration at",
    "steady state (Cmax,ss), in adults with non-muscle invasive bladder",
    "cancer and other solid tumors (Masters 2026; pooled studies B8011001,",
    "B8011006 and B8011007, n = 1007). The exposure slope -0.344 per e-fold",
    "is the published estimate (P = 0.043, not significant at the analysis",
    "alpha of 0.01). The intercept is NOT published: it was digitised from",
    "the fitted curve of Supporting Figure S2, so the model reproduces that",
    "plotted curve. The final published model may also carry covariate terms",
    "that are not reported; the plotted curve sits near the observed",
    "incidence of the non-NMIBC monotherapy studies (B8011001 + B8011007).",
    "No PK layer: CMAX is a data column. No random effect and no residual",
    "error are estimated.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Exposure slope from Table 4; intercept digitised from Supporting",
    "Information Figure S2 (panel 'TEAE Grade >=3').",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_teae_grade3 (probability of any grade >= 3 treatment-emergent adverse event, 0-1; also logit_teae_grade3)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual sasanlimab maximum serum concentration at STEADY STATE (Cmax,ss), predicted from the individual post hoc parameters of the companion population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL serum concentration; STEADY-STATE Cmax of the participant's regimen (IV 0.5-10 mg/kg every 3 weeks, SC 300 mg every 4 weeks or SC 600 mg every 6 weeks); empirical-Bayes prediction from modellib('Masters_2026_sasanlimab'). Enters on the NATURAL LOG scale (Table 4 row 'Log(Cmax ss (ug/mL))'). Observed 5th / 25th / 50th / 75th / 95th percentiles digitised from the Figure S2 reference lines: about 28.5 / 42.1 / 53.8 / 66 / 89 ug/mL. Not interchangeable with the single-dose Cmax of the grade >= 2 TEAE model.",
      source_name = "Cmax ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1007L,
    n_studies = 3L,
    n_events = 437L,
    n_observations = "1007 binary records, one per participant",
    age_range = "26-90 years",
    age_median = "66 years",
    weight_range = "38.2-173 kg",
    weight_median = "74 kg",
    sex_female_pct = 21.8,
    disease_state = "BCG-naive high-risk non-muscle invasive bladder cancer (70.1%, study B8011006 CREST Cohort A, sasanlimab + BCG), non-small cell lung cancer (19.7%), urothelial carcinoma (3.9%) and other advanced solid tumors (6.4%).",
    dose_range = "IV 0.5, 1, 3 or 10 mg/kg every 3 weeks (B8011001 dose escalation); SC 300 mg every 4 weeks (all three studies) and SC 600 mg every 6 weeks (B8011007).",
    notes = "Same 1007 participants as the population PK analysis (Results, ER Safety). Grade >= 3 TEAE incidence 43.4% overall; 41.8% in B8011001 (n = 146), 44.8% in B8011006 (n = 706) and 38.7% in B8011007 (n = 155) (Supporting Table S6). The incidence differs little between studies, so the plotted curve lies within a few points of every stratum."
  )

  ini({
    # logit(P) = intercept + slope * log(Cmax,ss) (Methods Equation 1 form).
    # Digitised by the maintainers from Figure S2 (panel 'TEAE Grade >=3',
    # 500 curve points): with the slope free the curve gives -0.3444
    # (printed -0.344), and with the slope held at -0.344 the intercept is
    # 1.125 (logit residual SD 0.003).
    logit_ref <- 1.125; label("Logit of the probability of grade >= 3 TEAE at Cmax,ss = 1 ug/mL, at the Figure S2 reference stratum (unitless logit)") # figure-derived: digitised from Supporting Figure S2 with the slope held at the Table 4 value; not printed in the paper
    e_cmax_logit <- -0.344; label("Log-odds of grade >= 3 TEAE per e-fold increase in Cmax,ss (unitless logit)") # Table 4: Log(Cmax ss (ug/mL)) = -0.344 (95% CI -0.681, -0.0126), P = 0.043, odds ratio 0.709

    # Bernoulli likelihood in the source; the additive residual below is a
    # placeholder so rxode2 has an error model to attach to the probability.
    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_teae_grade3 <- logit_ref + e_cmax_logit * log(CMAX)
    prob_teae_grade3 <- expit(logit_teae_grade3)
    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
