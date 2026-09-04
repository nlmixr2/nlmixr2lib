Kemal_2026_nemtabrutinib_ae <- function() {
  description <- "Logistic-regression exposure-safety model for any-grade investigator-assessed drug-related adverse events during nemtabrutinib monotherapy in patients with hematologic malignancies, from Kemal 2026. The probability of experiencing at least one any-grade drug-related AE is modelled on the logit scale as a linear function of the individual average on-treatment concentration (Cavg, a post-hoc exposure metric from the companion population PK model). Time on treatment was tested and was NOT a significant covariate for this endpoint, so no time term appears. Fitted with glm in R, not in NONMEM, on all 578 treated patients regardless of primary diagnosis. IMPORTANT: the source prints no coefficients for either exposure-safety model, so both values here are recovered by digitizing the fitted curve in the left panel of Figure 4 - see the model file's SOURCING NOTE and the vignette. Sister model files from the same paper: modellib('Kemal_2026_nemtabrutinib') for the population PK model that generates Cavg, modellib('Kemal_2026_nemtabrutinib_hypertension') for the companion hypertension endpoint, and modellib('Kemal_2026_nemtabrutinib_bor') for the exposure-efficacy model."
  reference <- paste(
    "Kemal CC, Zweers TJ, Krekels EHJ, Chatterjee MS.",
    "Population Pharmacokinetic Modeling and Exposure-Response Analyses of",
    "Nemtabrutinib in Patients With Hematologic Malignancies.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5). doi:10.1002/psp4.70257.",
    "Exposure-safety methods are in Methods section 2.4;",
    "the fitted curve this model is recovered from is the LEFT panel of",
    "Figure 4 (p. 10), and the narrative is in Results section 3.5.",
    "Sister model files from the same paper:",
    "modellib('Kemal_2026_nemtabrutinib'),",
    "modellib('Kemal_2026_nemtabrutinib_bor'),",
    "modellib('Kemal_2026_nemtabrutinib_hypertension').",
    sep = " "
  )
  vignette <- "Kemal_2026_nemtabrutinib"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # ------------------------------------------------------------------
  # SOURCING NOTE - how the coefficients in ini() were obtained.
  #
  # Kemal 2026 fits the exposure-safety relationships with glm in R and
  # reports them ONLY as the two fitted curves of Figure 4. Table S4 of
  # the supplement tabulates the exposure-EFFICACY model alone; there is
  # no safety equivalent, no coefficient appears in the main text, in any
  # supplemental table, in the figure caption, or as an annotation inside
  # the Figure 4 panels (the only numbers printed inside the panels are
  # the observed per-quartile event fractions 108/145, 109/145, 107/144
  # and 124/144). The panels were inspected at 600 dpi specifically to
  # check for printed coefficients before falling back to digitization.
  #
  # Both coefficients here are therefore DIGITIZED from the published
  # curve: the left panel of Figure 4 was rendered at 600 dpi, the axes
  # were calibrated on the tick marks (x: 0/500/1000/1500/2000 ng/mL;
  # y: 0.00/0.25/0.50/0.75/1.00), the black fitted line was traced at
  # 1042 independent x positions spanning Cavg 25-2163 ng/mL, and a
  # straight line was fitted to the traced points on the logit scale.
  #
  # Two things make this more than a curve-reading exercise.
  #
  # 1. THE FORM IS CONFIRMED, NOT ASSUMED. The traced points are
  #    straight on the logit scale to within 0.0097 on that scale across
  #    the whole 2100 ng/mL span - i.e. at the pixel-quantization limit.
  #    That independently confirms both that Cavg enters linearly on the
  #    logit scale and that no time-on-treatment term is present, which
  #    is what Results 3.5 states in prose ("Time on therapy was not
  #    found to be a significant covariate for these exposure-safety
  #    relationships").
  #
  # 2. THE METHOD IS CALIBRATED AGAINST AN ANSWER KEY IN THE SAME PAPER.
  #    Figure 3 of the same article draws the exposure-EFFICACY curve at
  #    two follow-up times, and that model's coefficients ARE printed, in
  #    Table S4. Running this identical digitization pipeline on Figure 3
  #    recovers intercepts of -3.2399 and -1.7624 against Table S4 values
  #    of -3.2466 and -1.7640 (errors +0.0066 and +0.0016 on the logit
  #    scale), and slopes of 0.001922 and 0.001924 from the two
  #    independent curves against a Table S4 value printed as 0.0019 to
  #    two significant figures. So the pipeline is accurate to better
  #    than 0.01 on an intercept and better than 1% on a slope, and the
  #    values below carry that uncertainty. (The Figure 3 exercise also
  #    shows the digitization is not merely as good as the paper here but
  #    slightly sharper on the slope than the two-significant-figure
  #    table.) The vignette reproduces this calibration as an assertion.
  #
  # Sanity check against the observed data drawn in the same panel: the
  # model gives 0.708 / 0.749 / 0.789 / 0.846 at the four exposure
  # quartile midpoints, against observed 108/145 = 0.745, 109/145 =
  # 0.752, 107/144 = 0.743 and 124/144 = 0.861.
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
      notes              = "A post-hoc exposure metric, not an observation. Kemal 2026 (Methods 2.3 and 2.4) simulated each participant's concentration-time profile from the companion population PK model using that participant's own dosing history - including dose interruptions and dose reductions - and their individual post hoc PK parameter estimates, then computed Cavg as the cumulative on-treatment AUC divided by the treatment duration. Generate it for simulation with modellib('Kemal_2026_nemtabrutinib'). Enters UNCENTRED, so the intercept is the log-odds at Cavg = 0. Observed distribution in the 578-patient safety cohort (Figure 4 reference boxplots): median about 600 ng/mL at 45 mg (n = 17), about 750 ng/mL at 65 mg (n = 434) and about 850 ng/mL at 80 mg (n = 106), with the pooled range extending to about 2160 ng/mL. Kemal 2026 also computed a per-patient Cmax and reports that Cavg and Cmax were significantly correlated, but used Cavg as the predictor for safety, tolerability and efficacy (Results 3.4); do not substitute a Cmax column here."
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
      notes              = "SCREENED AND REJECTED for this endpoint. Kemal 2026 Methods 2.4 states that time on treatment was investigated as a covariate for the safety endpoints, and Results 3.5 reports that 'Time on therapy was not found to be a significant covariate for these exposure-safety relationships'. It is therefore absent from this model, in deliberate contrast to the sister exposure-efficacy model modellib('Kemal_2026_nemtabrutinib_bor'), where the same covariate carries the single largest coefficient through a saturable term with ET50 = 200 days. The digitized Figure 4 curve is consistent with the exclusion: a single curve is drawn (rather than one per follow-up time as in Figure 3), and the traced points are straight on the logit scale.",
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
    endpoint       = "At least one any-grade drug-related adverse event during nemtabrutinib monotherapy, as assessed by the investigator. Observed rates by exposure quartile (printed in the left panel of Figure 4): 108/145, 109/145, 107/144 and 124/144.",
    notes          = "Kemal 2026 found a significant exposure trend for this endpoint and for any-grade hypertension only. No trend was found between exposure and Grade 3+ drug-related AEs, Grade 3+ hypertension, any-grade or Grade 3+ neutropenia, thrombocytopenia, anemia, infection, arrhythmia, diarrhea, rash or hemorrhage, nor between exposure and any tolerability endpoint (dose interruption, reduction, discontinuation, dose-limiting toxicity). Those negative endpoints have no fitted curve published and so are not packaged."
  )

  ini({
    # BOTH values below are DIGITIZED from the fitted curve in the LEFT
    # panel of Figure 4 (600 dpi trace, 1042 points, tick-mark axis
    # calibration); no coefficient for this model is printed anywhere in
    # the article or its supplement. Accuracy is bounded by running the
    # same pipeline on Figure 3, whose coefficients ARE printed in
    # Table S4: better than 0.01 on the intercept and better than 1% on
    # the slope. See the SOURCING NOTE above.
    #
    # No random effects and no residual-error parameter: the source fits
    # a fixed-effects logistic regression with glm, one row per patient.
    logitae <- 0.4733; label("Logit of the probability of an any-grade drug-related AE at Cavg = 0 (unitless)")  # digitized from Figure 4 left panel; intercept of the logit-scale trace

    e_css_nemta_logitae <- 0.001003; label("Effect of average on-treatment nemtabrutinib concentration on logit(any-grade drug-related AE), per ng/mL")  # digitized from Figure 4 left panel; slope of the logit-scale trace
  })

  model({
    # ----------------------------------------------------------------
    # Linear predictor on the logit scale. Cavg is UNCENTRED, and there
    # is deliberately no time-on-treatment term - Results 3.5 reports it
    # was not a significant covariate for this endpoint.
    # ----------------------------------------------------------------
    lpae <- logitae + e_css_nemta_logitae * CSS_NEMTA

    # pae is exposed as an output column so the fitted probabilities can
    # be compared against Figure 4 directly, without simulating draws.
    pae <- expit(lpae)

    ae ~ dbinom(1, pae)
  })
}
