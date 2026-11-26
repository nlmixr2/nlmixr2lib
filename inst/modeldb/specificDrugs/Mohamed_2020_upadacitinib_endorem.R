Mohamed_2020_upadacitinib_endorem <- function() {
  description <- paste(
    "Landmark exposure-response model for endoscopic remission (SES-CD <= 4",
    "and >= 2-point reduction from baseline, no subscore > 1) at week 12 or 16",
    "in adults with moderately-to-severely active Crohn's disease in the",
    "CELEST induction period on upadacitinib (Mohamed 2020). Hyperbolic Emax",
    "model on the probability scale, driven by the individual steady-state",
    "average plasma concentration (CAV); the placebo probability was",
    "estimated at zero and fixed. No covariates, no between-subject",
    "variability (Bernoulli likelihood)."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Lacerda AP, Othman AA. Exposure-Response",
    "Analyses for Upadacitinib Efficacy and Safety in the Crohn's Disease",
    "CELEST Study and Bridging to the Extended-Release Formulation. Clin",
    "Pharmacol Ther. 2020;107(3):639-649. doi:10.1002/cpt.1668"
  )
  vignette <- "Mohamed_2020_upadacitinib"
  units <- list(
    time = "n/a (static landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_endorem (probability of endoscopic remission at week 12 or 16, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual steady-state average upadacitinib plasma concentration over the dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL plasma concentration, entered untransformed in a hyperbolic Emax term, so it must be supplied in ng/mL. Table 2 footnote a: Cave = Dose / (CL/F x dosing interval) with the empirical-Bayes CL/F of each CELEST patient from the population PK model (Mohamed_2020_upadacitinib), and with simulated CL/F (and the 74.2% relative bioavailability) for extended-release regimens; 0 for placebo. Table 2 medians: 6.9, 14.9, 26.7 and 55.1 ng/mL for 3, 6, 12 and 24 mg IR twice daily, 25.1 ng/mL for 24 mg IR once daily, and 6.4, 12.9, 25.8, 38.7 and 51.6 ng/mL for 7.5, 15, 30, 45 and 60 mg ER once daily.",
      source_name = "Cave"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 216L,
    n_studies = 1L,
    age_range = "mean 40.8 years (SD 12.9)",
    weight_range = "mean 75.7 kg (SD 20.1)",
    sex_female_pct = 57,
    race_ethnicity = c(White = 86, Black = 6, Hispanic = 7, Asian = 1, Other = 0.5),
    disease_state = "Moderately-to-severely active Crohn's disease (CDAI 220-450, SES-CD >= 6 or >= 4 for ileal-only disease) with inadequate response or intolerance to an immunomodulator or anti-TNF therapy; 96% anti-TNF experienced. Baseline CDAI 300 (SD 59.9), SES-CD 14.9 (SD 8.4).",
    dose_range = "Placebo or immediate-release upadacitinib 3, 6, 12 or 24 mg twice daily or 24 mg once daily for 16 weeks (induction).",
    regions = "Multinational",
    notes = "CELEST (NCT02365649) induction period, Table 1. One ileocolonoscopy per patient after baseline, randomly assigned to week 12 or week 16, centrally read; the two visits are pooled. Nonlinear regression in SAS 9.4 on a Bernoulli outcome with nonresponder imputation."
  )

  ini({
    # Table S6 (endoscopic Endoscopic Remission row). Footnote b: each parameter was
    # estimated as exp(theta); values are the back-transformed estimates.
    e0 <- fixed(0)
    label("Probability of endoscopic remission without drug (fraction)") # Table S6 'Placebo' endoscopic remission = 0; footnote a: estimated to zero and subsequently fixed to zero
    lemax <- log(0.295)
    label("Maximum probability of endoscopic remission at saturating exposure (fraction)") # Table S6 'Emax' Endoscopic Remission = 0.295 (95% CI 0.0668-1.30); table note: the maximum observable added drug effect is Emax - Placebo
    lec50 <- log(26.9)
    label("Steady-state average concentration giving half of the maximum added effect (ng/mL)") # Table S6 'EC50 (ng/mL)' Endoscopic Remission = 26.9 (95% CI 1.44-502)

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_endorem <- fixed(0.001)
    label("Placeholder additive residual SD on prob_endorem (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Supplement Methods Equation 2:
    #   P(Y = 1) = Placebo + (Max - Placebo) * C / (EC50 + C),  C = Cave.
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    prob_endorem <- e0 + (emax - e0) * CAV / (ec50 + CAV)
    prob_endorem ~ add(addSd_prob_endorem)
  })
}
