Babel_2026_telisotuzumab_neuropathy <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for CTCAE grade 3 or ",
    "worse peripheral neuropathy in adults with non-small cell lung ",
    "cancer treated with telisotuzumab vedotin (Babel 2026, n = 284 ",
    "pooled from a phase 1 study and the LUMINOSITY phase 2 study). ",
    "The probability of an event is expit(a + b * log(CAV)), where ",
    "CAV is the individual average Teliso-V CONJUGATE serum ",
    "concentration up to the time of event or up to the end of ",
    "treatment, in ug/mL. There is no PK layer and no ODE: exposure ",
    "is supplied as data and was obtained in the source from post hoc ",
    "estimates of the companion population PK model ",
    "Babel_2026_telisotuzumab using actual doses received. No ",
    "covariate was retained. Peripheral sensory neuropathy is the ",
    "most common treatment-emergent adverse event with this ",
    "monomethyl auristatin E conjugate (30% all grades in ",
    "LUMINOSITY), and Babel 2026 finds it tracks CONJUGATE rather ",
    "than payload exposure. NOTE: Babel 2026 does not tabulate the ",
    "regression coefficients, so the two values here were recovered ",
    "by digitising the fitted curve of the left panel of Figure 4; ",
    "see the vignette for the digitisation and its cross-checks."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "The regression coefficients are not tabulated by the source; they were",
    "digitised from the fitted line in the left panel of Figure 4 and validated",
    "against the simulated event probabilities in Table 1.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_peripheral_neuropathy_grade3 (probability of a grade 3 or worse peripheral neuropathy event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average serum concentration of the telisotuzumab vedotin CONJUGATE, computed up to the time of the event or up to the end of treatment if no event occurred. Supplied as data: this model has no PK layer.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Same column and same derivation as in the companion efficacy model Babel_2026_telisotuzumab_orr, but computed on the larger exposure-SAFETY analysis set (284 patients from both studies rather than 193 LUMINOSITY patients). Babel 2026 Methods: post hoc estimates from the population PK model using actual doses received, averaged up to the time of the event. Reproduce the column with modellib('Babel_2026_telisotuzumab'). Enters on the NATURAL-LOG scale and UNCENTRED. Babel 2026 Figure 4 shows binned quartile medians near 3.6, 5.1, 6.6 and 8.2 ug/mL.",
      source_name        = "CavgADC"
    )
  )

  covariatesDataExcluded <- list(
    NEUROPATHY_HX = list(
      description = "History of peripheral neuropathy at baseline; 1 = yes, 0 = no.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Listed in Babel 2026 Table S2 as a safety covariate of interest and the most clinically plausible predictor for this particular endpoint, but not retained: 'No covariates were found to have a significant effect on efficacy or safety'. No point estimate exists on disk. The remaining Table S2 safety covariates (age, sex, race, ethnicity, body weight, c-Met expression level, prior therapy, number of prior systemic therapies, liver metastasis at baseline, treatment-emergent ADA status and nAb status) were screened and dropped on the same basis."
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
      "Table S6. Babel 2026 evaluated an endpoint further only when it ",
      "had more than 10 events and showed a trend in the quartile ",
      "plots. Grade 2 or worse peripheral neuropathy was also ",
      "significantly related to conjugate exposure and is plotted in ",
      "the same figure family, but only the grade 3 or worse panel is ",
      "shown in Figure 4, so only that one is recoverable. PK for one ",
      "patient was imputed from population estimates because the ",
      "patient discontinued before any PK sample was drawn."
    )
  )

  ini({
    # ==================================================================
    # Provenance: DIGITISED FROM THE LEFT PANEL OF FIGURE 4, not read
    # from the text or a table - Babel 2026 tabulates no exposure-safety
    # regression coefficients. This is a reporting gap with no competing
    # printed value.
    #
    # Method: the fitted solid line was traced at 300 dpi over 216
    # usable pixel columns (CavgADC 0.4 to 9.3 ug/mL) and both candidate
    # forms from Babel 2026 Methods ('linear and logarithmic logistic
    # regression analyses ... were evaluated') were fitted. The
    # logarithmic form fits with a root-mean-square residual of 0.08
    # percentage points against 0.25 for the linear form.
    #
    # Cross-check against printed values: Babel 2026 Table 1 gives
    # simulated median probabilities of 9.60% at 1.9 mg/kg Q2W and
    # 6.10% at 1.6 mg/kg Q2W. Those imply an exposure ratio of
    # exp((logit(0.0960) - logit(0.0610)) / b) = 1.194, against a
    # nominal dose ratio of 1.1875 - agreement to 0.6% on a quantity
    # that does not involve the intercept. Solving for the intercept,
    # the implied median CavgADC is 6.42 ug/mL at 1.9 mg/kg, consistent
    # with the 6.30 ug/mL implied independently by the ORR model.
    # ==================================================================
    logit_ref <- -7.437; label("Logit of the event probability at CavgADC = 1 ug/mL (unitless logit)")                                    # digitised from Babel 2026 Figure 4 (left panel); no printed value exists
    e_cav_pn  <- 2.793;  label("Log-odds of grade 3 or worse peripheral neuropathy per e-fold increase in conjugate CavgADC (unitless logit)") # digitised from Babel 2026 Figure 4 (left panel); no printed value exists

    # Binomial logistic regression with an exact Bernoulli likelihood:
    # no between-subject random effect and no residual error are
    # estimated. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_peripheral_neuropathy_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # CAV is uncentred, so logit_ref is the logit at CavgADC = 1 ug/mL,
    # below the observed range; it is not a reference-patient
    # probability.
    logit_peripheral_neuropathy_grade3 <- logit_ref + e_cav_pn * log(CAV)
    prob_peripheral_neuropathy_grade3  <- expit(logit_peripheral_neuropathy_grade3)

    prob_peripheral_neuropathy_grade3 ~ add(addSd_prob_peripheral_neuropathy_grade3)
  })
}
