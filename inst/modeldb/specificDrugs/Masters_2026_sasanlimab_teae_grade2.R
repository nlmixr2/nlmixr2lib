Masters_2026_sasanlimab_teae_grade2 <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for the probability",
    "of any grade >= 2 treatment-emergent adverse event (TEAE) as a function",
    "of the natural log of the sasanlimab maximum serum concentration after",
    "the first dose (Cmax,sd), in adults with non-muscle invasive bladder",
    "cancer and other solid tumors (Masters 2026; pooled studies B8011001,",
    "B8011006 and B8011007, n = 1007). The exposure slope 0.123 per e-fold is",
    "the published estimate (not significant, P = 0.62). The intercept is NOT",
    "published: it was digitised from the fitted curve of Supporting Figure",
    "S2, so the model reproduces that plotted curve. The final published",
    "model also carries covariate terms that are not reported; the plotted",
    "curve sits at the observed incidence of the non-NMIBC monotherapy",
    "studies (B8011001 + B8011007) and UNDER-predicts patients with NMIBC",
    "treated with sasanlimab + BCG (study B8011006). No PK layer: CMAX is a",
    "data column. No random effect and no residual error are estimated.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Exposure slope from Table 4; intercept digitised from Supporting",
    "Information Figure S2 (panel 'TEAE Grade >=2').",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_teae_grade2 (probability of any grade >= 2 treatment-emergent adverse event, 0-1; also logit_teae_grade2)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual sasanlimab maximum serum concentration after the FIRST dose (Cmax,sd), predicted from the individual post hoc parameters of the companion population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL serum concentration; SINGLE-DOSE Cmax (first dose of the participant's regimen, IV 0.5-10 mg/kg or SC 300 / 600 mg); empirical-Bayes prediction from modellib('Masters_2026_sasanlimab'). Enters on the NATURAL LOG scale (Table 4 row 'Log(Cmax sd (ug/mL))'; the Figure S2 axis 'Log(Cmax sd)' spans 1.1-5.5, i.e. 3-250 ug/mL, matching the linear 'Cmax sd' axis of the irAE grade >= 3 panel). Observed 5th / 25th / 50th / 75th / 95th percentiles digitised from the Figure S2 reference lines: about 14.9 / 21.7 / 27.1 / 33.5 / 55 ug/mL.",
      source_name = "Cmax sd"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1007L,
    n_studies = 3L,
    n_events = 852L,
    n_observations = "1007 binary records, one per participant",
    age_range = "26-90 years",
    age_median = "66 years",
    weight_range = "38.2-173 kg",
    weight_median = "74 kg",
    sex_female_pct = 21.8,
    disease_state = "BCG-naive high-risk non-muscle invasive bladder cancer (70.1%, study B8011006 CREST Cohort A, sasanlimab + BCG), non-small cell lung cancer (19.7%), urothelial carcinoma (3.9%) and other advanced solid tumors (6.4%).",
    dose_range = "IV 0.5, 1, 3 or 10 mg/kg every 3 weeks (B8011001 dose escalation); SC 300 mg every 4 weeks (all three studies) and SC 600 mg every 6 weeks (B8011007).",
    notes = "Same 1007 participants as the population PK analysis (Results, ER Safety). Grade >= 2 TEAE incidence 84.6% overall; 80.1% in B8011001 (n = 146), 86.7% in B8011006 (n = 706) and 79.4% in B8011007 (n = 155) (Supporting Table S6). The intercept reproduces the Figure S2 curve, which sits near the 79.7% pooled incidence of B8011001 + B8011007, not near the 84.6% overall incidence."
  )

  ini({
    # logit(P) = intercept + slope * log(Cmax,sd) (Methods Equation 1 form).
    # Digitised by the maintainers from Figure S2 (panel 'TEAE Grade >=2',
    # 500 curve points): with the slope free the curve gives 0.1233
    # (printed 0.123), and with the slope held at 0.123 the intercept is
    # 0.935 (logit residual SD 0.005). The intercept carries the unreported
    # covariate terms at the figure's (unstated) reference levels.
    logit_ref <- 0.935; label("Logit of the probability of grade >= 2 TEAE at Cmax,sd = 1 ug/mL, at the Figure S2 reference stratum (unitless logit)") # figure-derived: digitised from Supporting Figure S2 with the slope held at the Table 4 value; not printed in the paper
    e_cmax_logit <- 0.123; label("Log-odds of grade >= 2 TEAE per e-fold increase in Cmax,sd (unitless logit)") # Table 4: Log(Cmax sd (ug/mL)) = 0.123 (95% CI -0.363, 0.607), P = 0.62, odds ratio 1.13

    # Bernoulli likelihood in the source; the additive residual below is a
    # placeholder so rxode2 has an error model to attach to the probability.
    addSd_prob_teae_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_teae_grade2 <- logit_ref + e_cmax_logit * log(CMAX)
    prob_teae_grade2 <- expit(logit_teae_grade2)
    prob_teae_grade2 ~ add(addSd_prob_teae_grade2)
  })
}
