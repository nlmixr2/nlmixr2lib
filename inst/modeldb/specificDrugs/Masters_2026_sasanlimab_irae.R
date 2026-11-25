Masters_2026_sasanlimab_irae <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for the probability",
    "of an immune-related adverse event (irAE) of any grade as a function of",
    "the natural log of the sasanlimab trough serum concentration at steady",
    "state (Ctrough,ss), in adults with non-muscle invasive bladder cancer",
    "and other solid tumors (Masters 2026; pooled studies B8011001, B8011006",
    "and B8011007, n = 1007). The exposure slope -0.0788 per e-fold is the",
    "published estimate (not significant, P = 0.61). The intercept is NOT",
    "published: it was digitised from the fitted curve of Supporting Figure",
    "S2, so the model reproduces that plotted curve. The final published",
    "model also carries covariate terms that are not reported; the plotted",
    "curve (about 0.23 at the median exposure) sits near the observed",
    "incidence of the non-NMIBC monotherapy studies (25.6%, B8011001 +",
    "B8011007) and far UNDER the 67.0% incidence of patients with NMIBC",
    "treated with sasanlimab + BCG (study B8011006). Use it for non-NMIBC",
    "monotherapy only. No PK layer: CTROUGH is a data column. No random",
    "effect and no residual error are estimated.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Exposure slope from Table 4; intercept digitised from Supporting",
    "Information Figure S2 (panel 'irAE Any Grade').",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_irae (probability of an immune-related adverse event of any grade, 0-1; also logit_irae)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual sasanlimab trough serum concentration at STEADY STATE (Ctrough,ss), predicted from the individual post hoc parameters of the companion population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL serum concentration; STEADY-STATE trough of the participant's regimen (IV 0.5-10 mg/kg every 3 weeks, SC 300 mg every 4 weeks or SC 600 mg every 6 weeks); empirical-Bayes prediction from modellib('Masters_2026_sasanlimab'). Enters on the NATURAL LOG scale (Table 4 row 'Log(Ctrough ss (ug/mL))'). Observed 5th / 25th / 50th / 75th / 95th percentiles digitised from the Figure S2 reference lines: about 13.1 / 24.8 / 33.4 / 41.5 / 56.2 ug/mL.",
      source_name = "Ctrough ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1007L,
    n_studies = 3L,
    n_events = 550L,
    n_observations = "1007 binary records, one per participant",
    age_range = "26-90 years",
    age_median = "66 years",
    weight_range = "38.2-173 kg",
    weight_median = "74 kg",
    sex_female_pct = 21.8,
    disease_state = "BCG-naive high-risk non-muscle invasive bladder cancer (70.1%, study B8011006 CREST Cohort A, sasanlimab + BCG), non-small cell lung cancer (19.7%), urothelial carcinoma (3.9%) and other advanced solid tumors (6.4%).",
    dose_range = "IV 0.5, 1, 3 or 10 mg/kg every 3 weeks (B8011001 dose escalation); SC 300 mg every 4 weeks (all three studies) and SC 600 mg every 6 weeks (B8011007).",
    notes = "Same 1007 participants as the population PK analysis (Results, ER Safety). Any-grade irAE incidence 54.6% overall; 32.2% in B8011001 (n = 146), 67.0% in B8011006 (n = 706) and 19.4% in B8011007 (n = 155) (Supporting Table S6). A logistic model with an intercept averages to the observed incidence, so a curve at 0.23 cannot be the whole-population fit: the final model must hold a covariate (most likely NMIBC tumor type, study B8011006 or BCG combination, which coincide in 698 of 706 participants) whose coefficient is not printed."
  )

  ini({
    # logit(P) = intercept + slope * log(Ctrough,ss) (Methods Equation 1 form).
    # Digitised by the maintainers from Figure S2 (panel 'irAE Any Grade',
    # 500 curve points): with the slope free the curve gives -0.0784
    # (printed -0.0788), and with the slope held at -0.0788 the intercept is
    # -0.942 (logit residual SD 0.004). The intercept carries the unreported
    # covariate terms at the figure's (unstated) reference levels.
    logit_ref <- -0.942; label("Logit of the probability of any-grade irAE at Ctrough,ss = 1 ug/mL, at the Figure S2 reference stratum (unitless logit)") # figure-derived: digitised from Supporting Figure S2 with the slope held at the Table 4 value; not printed in the paper
    e_ctrough_logit <- -0.0788; label("Log-odds of any-grade irAE per e-fold increase in Ctrough,ss (unitless logit)") # Table 4: Log(Ctrough ss (ug/mL)) = -0.0788 (95% CI -0.382, 0.225), P = 0.61, odds ratio 0.924

    # Bernoulli likelihood in the source; the additive residual below is a
    # placeholder so rxode2 has an error model to attach to the probability.
    addSd_prob_irae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_irae <- logit_ref + e_ctrough_logit * log(CTROUGH)
    prob_irae <- expit(logit_irae)
    prob_irae ~ add(addSd_prob_irae)
  })
}
