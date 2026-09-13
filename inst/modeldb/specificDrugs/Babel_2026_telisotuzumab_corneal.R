Babel_2026_telisotuzumab_corneal <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for CTCAE grade 2 or ",
    "worse corneal epitheliopathy in adults with non-small cell lung ",
    "cancer treated with telisotuzumab vedotin (Babel 2026, n = 284 ",
    "pooled from a phase 1 study and the LUMINOSITY phase 2 study). ",
    "The probability of an event is expit(a + b * log(CAV)), where ",
    "CAV is the individual average Teliso-V CONJUGATE serum ",
    "concentration up to the time of event or up to the end of ",
    "treatment, in ug/mL. There is no PK layer and no ODE: exposure ",
    "is supplied as data and was obtained in the source from post hoc ",
    "estimates of the companion population PK model ",
    "Babel_2026_telisotuzumab using actual doses received. No ",
    "covariate was retained. This is the steepest of the three ",
    "conjugate-driven exposure-safety relationships in Babel 2026 ",
    "(nominal p = 3.69e-8), and the grade 3 or worse counterpart was ",
    "NOT modelled because only 2 such events occurred. NOTE: Babel ",
    "2026 does not tabulate the regression coefficients, so the two ",
    "values here were recovered by digitising the fitted curve of the ",
    "right panel of Figure 4; see the vignette for the digitisation ",
    "and its cross-checks."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "The regression coefficients are not tabulated by the source; they were",
    "digitised from the fitted line in the right panel of Figure 4 and validated",
    "against the simulated event probabilities in Table 1.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_corneal_epitheliopathy_grade2 (probability of a grade 2 or worse corneal epitheliopathy event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average serum concentration of the telisotuzumab vedotin CONJUGATE, computed up to the time of the event or up to the end of treatment if no event occurred. Supplied as data: this model has no PK layer.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Identical column and derivation to the companion model Babel_2026_telisotuzumab_neuropathy, on the same 284-patient exposure-safety analysis set. Babel 2026 Methods: post hoc estimates from the population PK model using actual doses received. Reproduce the column with modellib('Babel_2026_telisotuzumab'). Enters on the NATURAL-LOG scale and UNCENTRED. Babel 2026 Figure 4 shows binned quartile medians near 3.6, 5.1, 6.7 and 8.4 ug/mL, spanning to about 9.6 ug/mL.",
      source_name        = "CavgADC"
    )
  )

  covariatesDataExcluded <- list(
    LIVERMET = list(
      description = "Liver metastasis at baseline; 1 = present, 0 = absent.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed in Babel 2026 Table S2 as a safety covariate of interest but not retained: 'No covariates were found to have a significant effect on efficacy or safety'. No point estimate exists on disk. The remaining Table S2 safety covariates (age, sex, race, ethnicity, body weight, c-Met expression level, history of peripheral neuropathy, prior therapy, number of prior systemic therapies, treatment-emergent ADA status and nAb status) were screened and dropped on the same basis."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 284L,
    n_studies      = 2L,
    n_observations = "284 binary event records (one per patient; landmark analysis, no repeated measures)",
    disease_state  = "Advanced solid tumours, predominantly c-Met protein overexpressing non-small cell lung cancer; the phase 1 contribution is restricted to patients with NSCLC receiving monotherapy",
    dose_range     = "telisotuzumab vedotin 0.15-3.3 mg/kg every 3 weeks and 1.6-2.2 mg/kg every 2 weeks (phase 1) and 1.6 or 1.9 mg/kg every 2 weeks (LUMINOSITY)",
    notes          = paste0(
      "Baseline demographics of this analysis set are in Babel 2026 ",
      "Table S6. Babel 2026 states explicitly that 'grade >= 3 corneal ",
      "epitheliopathy events were not evaluated in the exposure-safety ",
      "analyses due to a very small number of events (N = 2)', so the ",
      "grade 2 threshold is the only one this endpoint supports. ",
      "Ocular toxicity is a recognised class effect of ",
      "microtubule-inhibitor antibody-drug conjugates."
    )
  )

  ini({
    # ==================================================================
    # Provenance: DIGITISED FROM THE RIGHT PANEL OF FIGURE 4, not read
    # from the text or a table - Babel 2026 tabulates no exposure-safety
    # regression coefficients. This is a reporting gap with no competing
    # printed value.
    #
    # Method: the fitted solid line was traced at 300 dpi over 206
    # usable pixel columns (CavgADC 0.4 to 9.6 ug/mL) and both candidate
    # forms from Babel 2026 Methods were fitted. The logarithmic form
    # fits with a root-mean-square residual of 0.08 percentage points
    # against 0.44 for the linear form.
    #
    # Cross-check against printed values: Babel 2026 Table 1 gives
    # simulated median probabilities of 15.0% at 1.9 mg/kg Q2W and
    # 7.80% at 1.6 mg/kg Q2W. Those imply an exposure ratio of
    # exp((logit(0.150) - logit(0.0780)) / b) = 1.190, against a
    # nominal dose ratio of 1.1875 - agreement to 0.2%, the closest of
    # the three endpoints, on a quantity that does not involve the
    # intercept. Solving for the intercept, the implied median CavgADC
    # is 6.78 ug/mL at 1.9 mg/kg; the companion neuropathy model, on
    # the same analysis set, implies 6.42 ug/mL. That 5.5% spread is
    # the practical accuracy of the digitised intercepts and is carried
    # forward as the tolerance of the vignette's absolute-probability
    # gate.
    # ==================================================================
    logit_ref <- -9.786; label("Logit of the event probability at CavgADC = 1 ug/mL (unitless logit)")                                     # digitised from Babel 2026 Figure 4 (right panel); no printed value exists
    e_cav_ce  <- 4.206;  label("Log-odds of grade 2 or worse corneal epitheliopathy per e-fold increase in conjugate CavgADC (unitless logit)") # digitised from Babel 2026 Figure 4 (right panel); no printed value exists

    # Binomial logistic regression with an exact Bernoulli likelihood:
    # no between-subject random effect and no residual error are
    # estimated. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_corneal_epitheliopathy_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # CAV is uncentred, so logit_ref is the logit at CavgADC = 1 ug/mL,
    # below the observed range; it is not a reference-patient
    # probability.
    logit_corneal_epitheliopathy_grade2 <- logit_ref + e_cav_ce * log(CAV)
    prob_corneal_epitheliopathy_grade2  <- expit(logit_corneal_epitheliopathy_grade2)

    prob_corneal_epitheliopathy_grade2 ~ add(addSd_prob_corneal_epitheliopathy_grade2)
  })
}
