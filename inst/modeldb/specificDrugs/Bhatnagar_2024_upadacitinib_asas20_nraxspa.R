Bhatnagar_2024_upadacitinib_asas20_nraxspa <- function() {
  description <- paste0(
    "Exposure-response logistic-regression model for the probability of an ",
    "ASAS20 response at week 14 in adults with non-radiographic axial ",
    "spondyloarthritis (nr-axSpA) treated with upadacitinib 15 mg once ",
    "daily or placebo in SELECT-AXIS 2 study 2. The paper's central finding ",
    "is that the log-odds of response carries NO upadacitinib exposure ",
    "term: a treatment-effect model, in which the active arm enters as a ",
    "pooled on/off switch, was selected by Akaike Information Criterion ",
    "over both linear and nonlinear functions of the model-estimated ",
    "average plasma concentration, and response rates did not rise with ",
    "increasing exposure across the range produced by the single 15 mg ",
    "once-daily regimen. The model is therefore a two-parameter logistic ",
    "regression in the treatment indicator, and exposure enters only ",
    "through the trial design that generated the exposure range. Fitted by ",
    "naive-pooled maximum likelihood with the R glm function; no ",
    "between-subject random effects and no covariates were retained. ",
    "Companion models from the same paper cover ASAS40 in nr-axSpA, ASAS20 ",
    "and ASAS40 in ankylosing spondylitis, and the population ",
    "pharmacokinetics that generated the exposure metric."
  )
  reference <- paste(
    "Bhatnagar S, Eckert D, Stodtmann S, Song I-H, Wung P, Liu W,",
    "Mohamed M-EF. Population pharmacokinetics and exposure-response",
    "analyses for efficacy and safety of upadacitinib in patients with",
    "axial spondyloarthritis. Clin Transl Sci. 2024;17(2):e13733.",
    "doi:10.1111/cts.13733.",
    "Parameter values are the 'ASAS20 Response in nr-axSpA' block of",
    "Table S4 of the Supporting Information; the model form is Bhatnagar",
    "2024 Equation 3. The companion population PK model that produced the",
    "exposure metric is modellib('Bhatnagar_2024_upadacitinib').",
    sep = " "
  )
  vignette <- "Bhatnagar_2024_upadacitinib"

  units <- list(
    time          = "week",
    dosing        = "(none; treatment enters as the binary ON_TREATMENT covariate, not as dosing events)",
    concentration = "(probability, 0-1; the observation prob_asas20 is the probability of an ASAS20 response at week 14, not a drug concentration)"
  )

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Randomized treatment-arm indicator (1 = upadacitinib 15 mg once daily, extended-release tablet; 0 = placebo).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm). The intercept is therefore the placebo log-odds of an ASAS20 response at week 14.",
      notes              = "Time-fixed per subject; patients were randomized 1:1 to upadacitinib 15 mg once daily or placebo, and the model uses only the double-blind placebo-controlled period through week 14. SELECT-AXIS 2 study 2 ran a 52-week placebo-controlled period, but the exposure-response analyses were deliberately restricted to the week-14 primary end point so that all three studies could be analysed on the same timeframe. This is the pooled on/off switch rather than an exposure-driven effect, which is exactly the case ON_TREATMENT is registered for: Bhatnagar 2024 evaluated treatment-effect, linear and nonlinear exposure models (Equations 3-5) and the treatment-effect model was selected. Model-estimated average plasma concentration was set to zero for placebo patients, so in the source the treatment indicator and a positive exposure are perfectly confounded and only one of them is identifiable.",
      source_name        = "Active treatment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 143L,
    n_studies      = 1L,
    age_range      = "19-79 years (Table S2, nr-axSpA column)",
    age_median     = "43 years (Table S2, nr-axSpA column; mean 44, SD 12)",
    weight_range   = "49.3-144 kg (Table S2, nr-axSpA column)",
    weight_median  = "79.8 kg (Table S2, nr-axSpA column; mean 81.1, SD 19.5)",
    sex_female_pct = 64,
    race_ethnicity = "Not tabulated in Bhatnagar 2024.",
    disease_state  = "Active non-radiographic axial spondyloarthritis (SELECT-AXIS 2 study 2). Unlike the ankylosing spondylitis cohorts, this population is predominantly female (64%).",
    dose_range     = "Upadacitinib 15 mg once daily, extended-release tablet, versus placebo (1:1 randomization).",
    regions        = "Multinational phase III study; not broken out in Bhatnagar 2024.",
    notes          = "ASAS20 is the Assessment of SpondyloArthritis International Society 20% response criterion, evaluated at week 14 alongside the primary ASAS40 end point; Bhatnagar 2024 Methods spells out the domain-level definition for ASAS40 only. The exposure-response analysis population is larger than the population PK analysis population (143 versus 71 patients with nr-axSpA) because PK samples were collected in only about 30% of SELECT-AXIS 2 patients. Fitted with the glm function in R 4.2.0; model selection by Akaike Information Criterion plus visual checks."
  )

  ini({
    # ======================================================================
    # Bhatnagar 2024 Table S4, "ASAS20 Response in nr-axSpA" block. The
    # model form is Equation 3 of the paper, the treatment-effect model:
    #   logit P(Y_i = 1) = alpha + beta_trt * TRT_i
    # No exposure term appears; Bhatnagar 2024 Results states that a
    # "treatment effect model best described the exposure-response
    # relationships for both the end points across all populations".
    #
    # SOURCE-TRACE CONFIRMATION. The point estimates give a placebo ASAS20
    # rate of expit(-0.511) = 37.5% and an upadacitinib rate of
    # expit(-0.511 + 1.18) = 66.1%, reproducing the clear separation from
    # placebo that Bhatnagar 2024 Figure 5 shows for the nr-axSpA panel.
    # Note the intercept is very close to the AS one (-0.456), so the two
    # populations differ mainly in the size of the treatment effect.
    #
    # No between-subject variability and no residual-error model are
    # reported: Table S4 lists fixed effects only and the endpoint is a
    # single binary record per subject. Encoded faithfully -- no invented
    # etas.
    # ======================================================================

    logitasas20 <- -0.511
    label("Logit of the probability of an ASAS20 response at week 14 on placebo (unitless)")
    # Table S4, ASAS20 Response in nr-axSpA, "Intercept": -0.511 (95% CI -0.988, -0.0337). expit(-0.511) = 0.375.

    e_on_treatment_logitasas20 <- 1.18
    label("Increase in the ASAS20 log-odds for upadacitinib 15 mg once daily versus placebo (unitless)")
    # Table S4, ASAS20 Response in nr-axSpA, "Active treatment": 1.18 (95% CI 0.498, 1.87), p = 0.000715. Odds ratio exp(1.18) = 3.25.

    # The source fits this model with a Bernoulli likelihood on a 0/1
    # response indicator, so there is no observation-error model to
    # translate. This placeholder additive residual is attached to the
    # probability output only so the nlmixr2 likelihood machinery accepts
    # the model for forward simulation. It is NOT from the source -- see the
    # vignette's Assumptions and deviations section. Same device as
    # Knebel_2012_istradefylline_dizziness.R.
    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on the probability output prob_asas20 (unitless); not from the source")
  })

  model({
    # Treatment-effect logistic regression, Bhatnagar 2024 Equation 3.
    lpasas20    <- logitasas20 + e_on_treatment_logitasas20 * ON_TREATMENT
    prob_asas20 <- expit(lpasas20)

    prob_asas20 ~ add(addSd)
  })
}
