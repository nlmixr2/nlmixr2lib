Babel_2026_telisotuzumab_teae <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for any CTCAE grade 3 ",
    "or worse treatment-emergent adverse event in adults with ",
    "non-small cell lung cancer treated with telisotuzumab vedotin ",
    "(Babel 2026, n = 284 pooled from a phase 1 study and the ",
    "LUMINOSITY phase 2 study). The probability of an event is ",
    "expit(a + b * CAV), where CAV is the individual average ",
    "unconjugated MMAE PAYLOAD plasma concentration up to the time of ",
    "event or up to the end of treatment, in ng/mL. This is the only ",
    "endpoint in Babel 2026 that tracks payload rather than conjugate ",
    "exposure: the Discussion states that 'grade >= 3 TEAEs were ",
    "correlated with unconjugated MMAE payload exposure metrics ",
    "only'. It is also the only one whose selected form is LINEAR in ",
    "exposure rather than logarithmic. There is no PK layer and no ",
    "ODE: exposure is supplied as data and was obtained in the source ",
    "from post hoc estimates of the companion payload population PK ",
    "model Babel_2026_telisotuzumab_mmae using actual doses received. ",
    "No covariate was retained. NOTE: Babel 2026 does not tabulate ",
    "the regression coefficients, so the two values here were ",
    "recovered by digitising the fitted curve of the CavgMMAE panel ",
    "of Figure 5; see the vignette for the digitisation, which is the ",
    "most precise of the four (root-mean-square residual 0.05 ",
    "percentage points)."
  )
  reference <- paste(
    "Babel H, Brunsdon P, Engelhardt B, Schmitt V, Ratajczak C, Mensing S,",
    "Menon RM, Parikh A. Population pharmacokinetics and exposure-response",
    "analyses for telisotuzumab vedotin in patients with c-Met protein",
    "overexpressing tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70219.",
    "doi:10.1002/psp4.70219. PMCID PMC12945708.",
    "The regression coefficients are not tabulated by the source; they were",
    "digitised from the fitted line in the upper-left (CavgMMAE) panel of",
    "Figure 5.",
    sep = " "
  )
  vignette <- "Babel_2026_telisotuzumab"

  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_teae_grade3 (probability of any grade 3 or worse treatment-emergent adverse event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average plasma concentration of the UNCONJUGATED MMAE payload, computed up to the time of the event or up to the end of treatment if no event occurred. Supplied as data: this model has no PK layer.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "PAYLOAD, not conjugate - the unit and the analyte both differ from the three companion conjugate-driven models (Babel_2026_telisotuzumab_orr, _neuropathy, _corneal), which carry CAV in ug/mL of conjugate. Supplying a conjugate Cavg here would overstate the predicted event rate enormously. Babel 2026 Methods: post hoc estimates from the payload population PK model using actual doses received. Reproduce the column with modellib('Babel_2026_telisotuzumab_mmae'), remembering that that model returns Cc in ug/mL, so multiply by 1000. Enters LINEARLY and UNCENTRED, so the intercept is interpretable directly as the logit at zero payload exposure. Babel 2026 Figure 5 shows the analysis-set range spanning roughly 0 to 4 ng/mL with binned quartile medians near 0.85, 1.35, 2.0 and 3.0 ng/mL.",
      source_name        = "CavgMMAE"
    )
  )

  covariatesDataExcluded <- list(
    ECOG_GE1 = list(
      description = "Baseline Eastern Cooperative Oncology Group performance status indicator; 1 = ECOG PS at least 1, 0 = ECOG PS 0.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Babel 2026 Table S2 lists baseline ECOG status as a covariate of interest for the exposure-EFFICACY regressions only, not for safety; it is recorded here because it is the covariate most often retained in oncology composite-tolerability models and its absence is a deliberate feature of this one. No point estimate exists on disk. Babel 2026 reports that 'no covariates were found to have a significant effect on efficacy or safety', so the Table S2 safety covariates (age, sex, race, ethnicity, body weight, c-Met expression level, history of peripheral neuropathy, liver metastasis at baseline, prior therapy, number of prior systemic therapies, treatment-emergent ADA status and nAb status) were screened and dropped."
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
      "Table S6. This composite endpoint is not mutually exclusive ",
      "with the companion grade-threshold endpoints: a grade 3 or ",
      "worse peripheral neuropathy event is also a grade 3 or worse ",
      "TEAE, so prob_teae_grade3 and ",
      "prob_peripheral_neuropathy_grade3 must not be treated as ",
      "competing risks. Babel 2026 reports the same endpoint against ",
      "three further payload exposure metrics in the remaining panels ",
      "of Figure 5 (CavgC1MMAE, CmingmMMAE and CmaxgmMMAE), all ",
      "significant; only the Cavg panel is encoded here so that the ",
      "exposure metric matches the companion conjugate-driven models."
    )
  )

  ini({
    # ==================================================================
    # Provenance: DIGITISED FROM THE UPPER-LEFT (CavgMMAE) PANEL OF
    # FIGURE 5, not read from the text or a table - Babel 2026
    # tabulates no exposure-safety regression coefficients. This is a
    # reporting gap with no competing printed value.
    #
    # Method: the fitted solid line was traced at 300 dpi over 288
    # usable pixel columns (CavgMMAE 0.02 to 4.6 ng/mL) and both
    # candidate forms from Babel 2026 Methods were fitted. Here the
    # LINEAR form is the selected one, and it is essentially exact: a
    # root-mean-square residual of 0.053 percentage points and a
    # maximum residual of 0.137 percentage points, against 7.6 and 19.9
    # for the logarithmic form. The reading is unambiguous in the
    # figure as well - the plotted curve meets the left axis at a
    # finite 25.3%, whereas any logarithmic form is pinned to zero
    # there.
    #
    # Because the exposure enters linearly and uncentred, logit_ref is
    # the genuine zero-exposure intercept and its digitised value is
    # directly checkable: expit(-1.0835) = 25.3%, matching the left
    # edge of the plotted curve.
    # ==================================================================
    logit_ref   <- -1.0835; label("Logit of the event probability at zero payload exposure (unitless logit)")                            # digitised from Babel 2026 Figure 5 (CavgMMAE panel); no printed value exists
    e_cav_teae  <- 0.7195;  label("Log-odds of any grade 3 or worse TEAE per 1 ng/mL increase in unconjugated MMAE CavgMMAE (unitless logit)") # digitised from Babel 2026 Figure 5 (CavgMMAE panel); no printed value exists

    # Binomial logistic regression with an exact Bernoulli likelihood:
    # no between-subject random effect and no residual error are
    # estimated. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_teae_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the natural exposure scale, uncentred, so
    # logit_ref is the logit at CavgMMAE = 0 ng/mL.
    logit_teae_grade3 <- logit_ref + e_cav_teae * CAV
    prob_teae_grade3  <- expit(logit_teae_grade3)

    prob_teae_grade3 ~ add(addSd_prob_teae_grade3)
  })
}
