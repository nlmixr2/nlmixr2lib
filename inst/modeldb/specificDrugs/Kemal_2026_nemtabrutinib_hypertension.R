Kemal_2026_nemtabrutinib_hypertension <- function() {
  description <- "Logistic-regression exposure-safety model for any-grade hypertension events during nemtabrutinib monotherapy in patients with hematologic malignancies, from Kemal 2026. The probability of experiencing at least one any-grade hypertension event is modelled on the logit scale as a linear function of the individual average on-treatment concentration (Cavg, a post-hoc exposure metric from the companion population PK model). Time on treatment was tested and was NOT a significant covariate for this endpoint, so no time term appears; nor was any trend found for Grade 3+ hypertension, which is therefore not packaged. Fitted with glm in R, not in NONMEM, on all 578 treated patients regardless of primary diagnosis. IMPORTANT: the source prints no coefficients for either exposure-safety model, so both values here are recovered by digitizing the fitted curve in the right panel of Figure 4 - see the model file's SOURCING NOTE and the vignette. Sister model files from the same paper: modellib('Kemal_2026_nemtabrutinib') for the population PK model that generates Cavg, modellib('Kemal_2026_nemtabrutinib_ae') for the companion any-grade drug-related AE endpoint, and modellib('Kemal_2026_nemtabrutinib_bor') for the exposure-efficacy model."
  reference <- paste(
    "Kemal CC, Zweers TJ, Krekels EHJ, Chatterjee MS.",
    "Population Pharmacokinetic Modeling and Exposure-Response Analyses of",
    "Nemtabrutinib in Patients With Hematologic Malignancies.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5). doi:10.1002/psp4.70257.",
    "Exposure-safety methods are in Methods section 2.4;",
    "the fitted curve this model is recovered from is the RIGHT panel of",
    "Figure 4 (p. 10), and the narrative is in Results section 3.5.",
    "Sister model files from the same paper:",
    "modellib('Kemal_2026_nemtabrutinib'),",
    "modellib('Kemal_2026_nemtabrutinib_bor'),",
    "modellib('Kemal_2026_nemtabrutinib_ae').",
    sep = " "
  )
  vignette <- "Kemal_2026_nemtabrutinib"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # ------------------------------------------------------------------
  # SOURCING NOTE - how the coefficients in ini() were obtained.
  #
  # Identical situation, and identical method, to the sister model
  # modellib('Kemal_2026_nemtabrutinib_ae'); only the panel differs.
  #
  # Kemal 2026 reports the exposure-safety relationships ONLY as the two
  # fitted curves of Figure 4. Table S4 of the supplement tabulates the
  # exposure-EFFICACY model alone; no safety coefficient appears in the
  # main text, in any supplemental table, in the figure caption, or as an
  # annotation inside the Figure 4 panels (the only numbers printed
  # inside the right panel are the observed per-quartile event fractions
  # 13/145, 22/145, 28/144 and 25/144). The panels were inspected at
  # 600 dpi specifically to check for printed coefficients before falling
  # back to digitization.
  #
  # Both coefficients here are therefore DIGITIZED: the right panel of
  # Figure 4 was rendered at 600 dpi, the axes were calibrated on the
  # tick marks (x: 0/500/1000/1500/2000 ng/mL; y: 0.00/0.25/0.50/0.75/
  # 1.00), the black fitted line was traced at 1091 independent x
  # positions spanning Cavg 26-2164 ng/mL, and a straight line was fitted
  # to the traced points on the logit scale.
  #
  # As with the AE endpoint, two things make this more than a
  # curve-reading exercise.
  #
  # 1. THE FORM IS CONFIRMED, NOT ASSUMED. The traced points are straight
  #    on the logit scale to within 0.032 on that scale across the whole
  #    2100 ng/mL span. That independently confirms both that Cavg enters
  #    linearly on the logit scale and that no time-on-treatment term is
  #    present, matching the prose of Results 3.5.
  #
  # 2. THE METHOD IS CALIBRATED AGAINST AN ANSWER KEY IN THE SAME PAPER.
  #    Running this identical pipeline on Figure 3 - whose coefficients
  #    ARE printed, in Table S4 - recovers intercepts of -3.2399 and
  #    -1.7624 against printed values of -3.2466 and -1.7640, and slopes
  #    of 0.001922 and 0.001924 from the two independent curves against a
  #    printed 0.0019 (two significant figures). The pipeline is
  #    therefore accurate to better than 0.01 on an intercept and better
  #    than 1% on a slope, and the values below carry that uncertainty.
  #    The vignette reproduces this calibration as an assertion.
  #
  # Sanity check against the observed data drawn in the same panel: the
  # model gives 0.116 / 0.134 / 0.154 / 0.197 at the four exposure
  # quartile midpoints, against observed 13/145 = 0.090, 22/145 = 0.152,
  # 28/144 = 0.194 and 25/144 = 0.174.
  #
  # No standard errors are carried: the shaded band in Figure 4 is a
  # confidence interval on a fitted probability, which does not invert to
  # standard errors on the coefficients without the design matrix.
  # ------------------------------------------------------------------

  covariateData <- list(
    CSS_NEMTA = list(
      description        = "Individual average on-treatment plasma concentration of nemtabrutinib (Cavg).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A post-hoc exposure metric, not an observation. Kemal 2026 (Methods 2.3 and 2.4) simulated each participant's concentration-time profile from the companion population PK model using that participant's own dosing history - including dose interruptions and dose reductions - and their individual post hoc PK parameter estimates, then computed Cavg as the cumulative on-treatment AUC divided by the treatment duration. Generate it for simulation with modellib('Kemal_2026_nemtabrutinib'). Enters UNCENTRED, so the intercept is the log-odds at Cavg = 0. Observed distribution in the 578-patient safety cohort (Figure 4 reference boxplots): median about 600 ng/mL at 45 mg (n = 17), about 750 ng/mL at 65 mg (n = 434) and about 850 ng/mL at 80 mg (n = 106), with the pooled range extending to about 2160 ng/mL. The same column drives the sister AE model with a different slope; the two endpoints share the exposure metric, not the coefficients."
    )
  )

  # Covariate tested for this endpoint and explicitly rejected, so no
  # coefficient exists to carry. Documentation only; checkModelConventions()
  # does not require these to be referenced in model().
  covariatesDataExcluded <- list(
    T_TRT = list(
      description        = "Time on nemtabrutinib treatment (the paper's 'time on therapy').",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SCREENED AND REJECTED for this endpoint. Kemal 2026 Methods 2.4 states that time on treatment was investigated as a covariate for the safety endpoints, and Results 3.5 reports that 'Time on therapy was not found to be a significant covariate for these exposure-safety relationships'. It is therefore absent from this model, in deliberate contrast to the sister exposure-efficacy model modellib('Kemal_2026_nemtabrutinib_bor'), where the same covariate carries the single largest coefficient through a saturable term with ET50 = 200 days.",
      source_name        = "time on treatment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 578L,
    n_studies      = 2L,
    age_range      = "25-89 years",
    weight_range   = "41.2-147 kg",
    sex_female_pct = 34.1,
    disease_state  = "Hematologic malignancies, all primary diagnoses pooled: CLL/SLL (49.8%), other hematologic malignancy (30.6%), B-cell non-Hodgkin lymphoma (9.7%) and Waldenstrom's macroglobulinemia (9.7%). Unlike the sister exposure-efficacy model, which is restricted to the CLL/SLL subset, the exposure-safety analysis used every treated patient in the population PK analysis set regardless of indication.",
    dose_range     = "5-80 mg nemtabrutinib once daily orally; the doses drawn as reference boxplots in Figure 4 are 45, 65 and 80 mg (n = 17, 434 and 106 respectively).",
    endpoint       = "At least one any-grade hypertension event during nemtabrutinib monotherapy. Observed rates by exposure quartile (printed in the right panel of Figure 4): 13/145, 22/145, 28/144 and 25/144.",
    notes          = "Hypertension is a recognised class effect of BTK inhibitors and was one of the treatment-related events leading to discontinuation in BELLWAVE-001 (1.5%). Kemal 2026 found a significant exposure trend for any-grade hypertension but explicitly NO trend for Grade 3+ hypertension, so only the any-grade endpoint has a published curve and only it is packaged here."
  )

  ini({
    # BOTH values below are DIGITIZED from the fitted curve in the RIGHT
    # panel of Figure 4 (600 dpi trace, 1091 points, tick-mark axis
    # calibration); no coefficient for this model is printed anywhere in
    # the article or its supplement. Accuracy is bounded by running the
    # same pipeline on Figure 3, whose coefficients ARE printed in
    # Table S4: better than 0.01 on the intercept and better than 1% on
    # the slope. See the SOURCING NOTE above.
    #
    # No random effects and no residual-error parameter: the source fits
    # a fixed-effects logistic regression with glm, one row per patient.
    logithtn <- -2.3463; label("Logit of the probability of an any-grade hypertension event at Cavg = 0 (unitless)")  # digitized from Figure 4 right panel; intercept of the logit-scale trace

    e_css_nemta_logithtn <- 0.0007645; label("Effect of average on-treatment nemtabrutinib concentration on logit(any-grade hypertension), per ng/mL")  # digitized from Figure 4 right panel; slope of the logit-scale trace
  })

  model({
    # ----------------------------------------------------------------
    # Linear predictor on the logit scale. Cavg is UNCENTRED, and there
    # is deliberately no time-on-treatment term - Results 3.5 reports it
    # was not a significant covariate for this endpoint.
    # ----------------------------------------------------------------
    lphtn <- logithtn + e_css_nemta_logithtn * CSS_NEMTA

    # phtn is exposed as an output column so the fitted probabilities can
    # be compared against Figure 4 directly, without simulating draws.
    phtn <- expit(lphtn)

    htn ~ dbinom(1, phtn)
  })
}
