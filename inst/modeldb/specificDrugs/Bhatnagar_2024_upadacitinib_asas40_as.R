Bhatnagar_2024_upadacitinib_asas40_as <- function() {
  description <- paste0(
    "Exposure-response logistic-regression model for the probability of an ",
    "ASAS40 response at week 14 in adults with ankylosing spondylitis (AS) ",
    "treated with upadacitinib 15 mg once daily or placebo in the ",
    "SELECT-AXIS 1 and SELECT-AXIS 2 trials. The paper's central finding is ",
    "that the log-odds of response carries NO upadacitinib exposure term: a ",
    "treatment-effect model, in which the active arm enters as a pooled ",
    "on/off switch, was selected by Akaike Information Criterion over both ",
    "linear and nonlinear functions of the model-estimated average plasma ",
    "concentration, and response rates did not rise with increasing ",
    "exposure across the range produced by the single 15 mg once-daily ",
    "regimen. The model is therefore a two-parameter logistic regression in ",
    "the treatment indicator, and exposure enters only through the trial ",
    "design that generated the exposure range. Fitted by naive-pooled ",
    "maximum likelihood with the R glm function; no between-subject random ",
    "effects and no covariates were retained. Companion models from the ",
    "same paper cover ASAS20 in AS, ASAS20 and ASAS40 in non-radiographic ",
    "axial spondyloarthritis, and the population pharmacokinetics that ",
    "generated the exposure metric."
  )
  reference <- paste(
    "Bhatnagar S, Eckert D, Stodtmann S, Song I-H, Wung P, Liu W,",
    "Mohamed M-EF. Population pharmacokinetics and exposure-response",
    "analyses for efficacy and safety of upadacitinib in patients with",
    "axial spondyloarthritis. Clin Transl Sci. 2024;17(2):e13733.",
    "doi:10.1111/cts.13733.",
    "Parameter values are the 'ASAS40 Response in AS' block of Table S4 of",
    "the Supporting Information; the model form is Bhatnagar 2024",
    "Equation 3. The companion population PK model that produced the",
    "exposure metric is modellib('Bhatnagar_2024_upadacitinib').",
    sep = " "
  )
  vignette <- "Bhatnagar_2024_upadacitinib"

  units <- list(
    time          = "week",
    dosing        = "(none; treatment enters as the binary ON_TREATMENT covariate, not as dosing events)",
    concentration = "(probability, 0-1; the observation prob_asas40 is the probability of an ASAS40 response at week 14, not a drug concentration)"
  )

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Randomized treatment-arm indicator (1 = upadacitinib 15 mg once daily, extended-release tablet; 0 = placebo).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm). The intercept is therefore the placebo log-odds of an ASAS40 response at week 14.",
      notes              = "Time-fixed per subject; patients were randomized 1:1 to upadacitinib 15 mg once daily or placebo in all three studies, and the model uses only the double-blind placebo-controlled period through week 14. This is the pooled on/off switch rather than an exposure-driven effect, which is exactly the case ON_TREATMENT is registered for: Bhatnagar 2024 evaluated treatment-effect, linear and nonlinear exposure models (Equations 3-5) and the treatment-effect model was selected. Model-estimated average plasma concentration was set to zero for placebo patients, so in the source the treatment indicator and a positive exposure are perfectly confounded and only one of them is identifiable.",
      source_name        = "Active treatment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 339L,
    n_studies      = 2L,
    age_range      = "21-82 years (Table S2, AS column)",
    age_median     = "45 years (Table S2, AS column; mean 45, SD 13)",
    weight_range   = "41.5-156 kg (Table S2, AS column)",
    weight_median  = "78.0 kg (Table S2, AS column; mean 80.4, SD 19.0)",
    sex_female_pct = 28,
    race_ethnicity = "Not tabulated in Bhatnagar 2024.",
    disease_state  = "Active ankylosing spondylitis fulfilling the modified New York criteria. Pooled from SELECT-AXIS 1 (bDMARD-naive, inadequate response or intolerance to at least two NSAIDs) and SELECT-AXIS 2 study 1 (inadequate response or intolerance to biologic DMARDs).",
    dose_range     = "Upadacitinib 15 mg once daily, extended-release tablet, versus placebo (1:1 randomization).",
    regions        = "Multinational phase II/III programme; not broken out in Bhatnagar 2024.",
    notes          = "ASAS40 response is defined in Bhatnagar 2024 Methods as at least 40% improvement and an absolute improvement of at least two units on a 0-10 numerical rating scale from baseline in at least three of four domains (patient global assessment of disease activity, patient assessment of back pain, Bath Ankylosing Spondylitis Functional Index, and inflammation taken as the mean of the two Bath Ankylosing Spondylitis Disease Activity Index morning-stiffness questions), with no worsening in the remaining domain. It was the primary end point of all three studies and was assessed at week 14. The exposure-response analysis population is larger than the population PK analysis population (339 versus 173 patients with AS) because PK samples were collected in only about 30% of SELECT-AXIS 2 patients. Fitted with the glm function in R 4.2.0; model selection by Akaike Information Criterion plus visual checks."
  )

  ini({
    # ======================================================================
    # Bhatnagar 2024 Table S4, "ASAS40 Response in AS" block. The model form
    # is Equation 3 of the paper, the treatment-effect model:
    #   logit P(Y_i = 1) = alpha + beta_trt * TRT_i
    # No exposure term appears; Bhatnagar 2024 Results states that a
    # "treatment effect model best described the exposure-response
    # relationships for both the end points across all populations".
    #
    # SOURCE-TRACE CONFIRMATION against the parent trial publications. The
    # point estimates give a placebo ASAS40 rate of expit(-1.31) = 21.2% and
    # an upadacitinib rate of expit(-1.31 + 1.06) = 43.8%. SELECT-AXIS 1
    # reported 26% versus 52% and SELECT-AXIS 2 study 1 reported 18% versus
    # 45%, so the pooled AS values sit between the two component studies as
    # they must.
    #
    # No between-subject variability and no residual-error model are
    # reported: Table S4 lists fixed effects only and the endpoint is a
    # single binary record per subject. Encoded faithfully -- no invented
    # etas.
    # ======================================================================

    logitasas40 <- -1.31
    label("Logit of the probability of an ASAS40 response at week 14 on placebo (unitless)")
    # Table S4, ASAS40 Response in AS, "Intercept": -1.31 (95% CI -1.69, -0.939). expit(-1.31) = 0.212.

    e_on_treatment_logitasas40 <- 1.06
    label("Increase in the ASAS40 log-odds for upadacitinib 15 mg once daily versus placebo (unitless)")
    # Table S4, ASAS40 Response in AS, "Active treatment": 1.06 (95% CI 0.579, 1.54), p = 0.0000147. Odds ratio exp(1.06) = 2.89.

    # The source fits this model with a Bernoulli likelihood on a 0/1
    # response indicator, so there is no observation-error model to
    # translate. This placeholder additive residual is attached to the
    # probability output only so the nlmixr2 likelihood machinery accepts
    # the model for forward simulation. It is NOT from the source -- see the
    # vignette's Assumptions and deviations section. Same device as
    # Knebel_2012_istradefylline_dizziness.R.
    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on the probability output prob_asas40 (unitless); not from the source")
  })

  model({
    # Treatment-effect logistic regression, Bhatnagar 2024 Equation 3.
    lpasas40    <- logitasas40 + e_on_treatment_logitasas40 * ON_TREATMENT
    prob_asas40 <- expit(lpasas40)

    prob_asas40 ~ add(addSd)
  })
}
