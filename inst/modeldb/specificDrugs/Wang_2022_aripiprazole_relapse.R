Wang_2022_aripiprazole_relapse <- function() {
  description <- "Exponential (constant-hazard) time-to-event exposure-response model for exacerbation of psychotic symptoms / impending relapse in 615 adults with schizophrenia treated with aripiprazole once monthly (AOM), oral aripiprazole, or placebo in the phase 3 studies 31-07-246 and 31-07-247 (Wang 2022). The hazard of impending relapse is constant in time and takes one of two values according to whether the model-predicted aripiprazole minimum concentration proximate to the event, CMIN_ARI, reaches the 95 ng/mL efficacy threshold the paper identified: subjects at or above the threshold have a hazard 4.41-fold lower (equivalently a 4.41-fold longer expected time to relapse) than subjects below it. The source paper fits the model as an accelerated-failure-time exponential regression on log survival time; it is re-expressed here on the log-hazard scale, which is the exponential model's equivalent parameterisation and the one rxode2 forward simulation needs. The exposure driver CMIN_ARI is supplied as a data column and is generated in the source analysis from the companion population PK model modellib('Wang_2022_aripiprazole'). Forward simulation exposes `hazard` (instantaneous relapse rate per day) and `sur` (probability of no impending relapse since randomisation) as derived outputs."
  reference <- paste(
    "Wang X, Raoufinia A, Bihorel S, Passarell J, Mallikaarjun S, Phillips L.",
    "Population Pharmacokinetic Modeling and Exposure-Response Analysis for",
    "Aripiprazole Once Monthly in Subjects With Schizophrenia.",
    "Clin Pharmacol Drug Dev. 2022;11(2):150-164. doi:10.1002/cpdd.1022.",
    "Upstream popPK (driver of CMIN_ARI): same publication;",
    "see modellib('Wang_2022_aripiprazole').",
    sep = " "
  )
  vignette <- "Wang_2022_aripiprazole"

  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the exposure driver is supplied as the CMIN_ARI data covariate, in ng/mL)",
    concentration = "probability (the model output `sur` is the relapse-free survival probability, not a drug concentration)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    cumhaz = list(analyte = "cumulative hazard of impending relapse", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CMIN_ARI = list(
      description        = "Model-predicted minimum aripiprazole plasma concentration proximate to the impending-relapse event or censoring time (ng/mL).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Time-fixed per subject in the source analysis (one value per subject, taken proximate to the event or censoring time). Derived in Wang 2022 by applying the companion final population PK model -- modellib('Wang_2022_aripiprazole') -- to each subject's recorded dosing history together with that subject's empirical-Bayes PK parameters, by numerical integration in NONMEM (Wang 2022 Methods, 'Definition of PK Exposures'). Two prediction windows are used: if the event fell in the first AOM dose period the predicted concentration 24 h after the oral dose preceding the event is used; if it fell in the second or a later AOM dosing period the predicted concentration 672 h after the AOM dose preceding the event is used. Concentrations were also predicted for placebo subjects because washout of aripiprazole taken before randomisation was incomplete; the minimum model-predicted value over the whole analysis set was 0.000329 ng/mL, not zero (Wang 2022 Figure 5 note). Only the dichotomy at 95 ng/mL enters the model: Wang 2022 first fitted CMIN_ARI as a continuous predictor and found the fit biased above 95 ng/mL, because the observed time to relapse was similar across the three highest concentration quartiles, and therefore replaced the continuous term with the threshold indicator. Observed subject-level means (SD) in the analysis set were 47.8 (65.3) ng/mL for subjects assigned to placebo and 180.9 (83.1) ng/mL for subjects assigned to 400 mg AOM; 101.7 (90.0) ng/mL in subjects who relapsed and 161.5 (93.6) ng/mL in censored subjects. The observed range is 0 to 580 ng/mL (Wang 2022 Figure 6 panel headings), and 154 of the 615 subjects fell below the threshold.",
      source_name        = "predicted aripiprazole Cmin proximate to the event"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 615L,
    n_studies      = 2L,
    n_observations = "85 impending-relapse events and 530 censored subjects (Wang 2022 Results, 'Exposure-Response Analysis of Time to Relapse').",
    age_range      = "Adults; per-study demographics are in Wang 2022 Table S1 / S2 (supplement not on disk -- see the vignette Errata).",
    weight_range   = "Not reported in the main text.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in the main text for the exposure-response analysis set.",
    disease_state  = "Adults with schizophrenia who had completed the lead-in / stabilisation stages of studies 31-07-246 and 31-07-247 and had at least one documented pharmacodynamic end point plus a model-predicted aripiprazole Cmin proximate to the event.",
    dose_range     = "120 subjects on oral aripiprazole and 495 subjects on 400/300 mg AOM given every 28 days in the gluteal muscle, plus the placebo arms of the two studies. A further 121 subjects on 50/25 mg AOM were EXCLUDED because the companion population PK model did not adequately describe concentrations at that dose.",
    regions        = "Multicentre; not further specified in the main text.",
    notes          = "The end point is exacerbation of psychotic symptoms / impending relapse, defined during the randomised stage of each study by any of: a Clinical Global Impression of Improvement score >= 5 together with a specified worsening on individual Positive and Negative Syndrome Scale items; hospitalisation for worsening psychotic symptoms; a Clinical Global Impression of Severity of Suicide part 1 score of 4-5 or part 2 score of 6-7; or violent behaviour causing clinically significant injury or property damage (Wang 2022 Methods, 'Definition of Time to Relapse'). Time is measured from randomisation (stage 4 of study 31-07-246, stage 3 of study 31-07-247) and follow-up ran to about 365 days (Wang 2022 Figures 5 and 6). Of the 85 events, 48 (56.5%) were in placebo subjects and 37 (43.5%) in the 400/300 mg AOM group. The proportional-hazards assumption was tested and rejected, so a parametric model was used; among the exponential, Weibull, log-logistic and generalized gamma distributions the exponential fitted best. Fitted with the SAS 9.2 LIFEREG procedure."
  )

  ini({
    # ==================================================================
    # Wang 2022 Table 2 ("Parameter Estimates for Exponential Survival
    # Model of Impending Relapse as a Threshold Function of Model-
    # Predicted Aripiprazole Cmin").
    #
    # REPARAMETERISATION (the only non-transcription step in this file).
    # Table 2 reports an accelerated-failure-time (AFT) exponential
    # regression on log survival time, as fitted by SAS PROC LIFEREG:
    #
    #   log(T) = 6.256 + 1.484 * I(Cmin >= 95 ng/mL)
    #
    # For the exponential distribution the AFT and proportional-hazards
    # parameterisations coincide, with a constant hazard equal to the
    # reciprocal of the expected survival time:
    #
    #   h = 1 / E[T] = exp(-(6.256 + 1.484 * I))
    #
    # so the log-hazard intercept is the NEGATIVE of the Table 2
    # intercept and the log-hazard covariate coefficient is the NEGATIVE
    # of the Table 2 coefficient. Both signs are checked by the paper's
    # own derived quantity: exp(1.484) = 4.41, the Table 2 "calculated
    # hazard ratio of expected survival time", which must be a hazard
    # REDUCTION for subjects at or above the threshold.
    #
    # Time unit is days: the Wang 2022 Figure 5 and Figure 6 x-axes are
    # labelled "Time to Relapse (days)" and run to about 365 days.
    # ==================================================================

    lhaz_base  <- -6.256 ; label("Log baseline hazard of impending relapse when CMIN_ARI < 95 ng/mL (log(1/day))")       # Wang 2022 Table 2: intercept = 6.256 (SE 0.1474, 95% CI 5.97-6.55, P < .0001) on the AFT log-survival-time scale; negated to the log-hazard scale
    e_cmin_haz <- -1.484 ; label("Additive shift in the log hazard of impending relapse when CMIN_ARI >= 95 ng/mL (unitless)")  # Wang 2022 Table 2: predicted aripiprazole (Cmin >= 95 ng/mL) proximate to the event = 1.484 (SE 0.2177, 95% CI 1.06-1.91, P < .0001) on the AFT log-survival-time scale; negated to the log-hazard scale. exp(-e_cmin_haz) = exp(1.484) = 4.41 = the Table 2 hazard ratio (95% CI 2.89-6.75)

    # No interindividual variability and no residual error. Wang 2022
    # Table 2 reports only the two fixed effects and their standard
    # errors; a parametric survival regression of this kind carries its
    # likelihood in the event / censoring density itself rather than in
    # an observation-error model. Forward simulation in rxode2 exposes
    # `hazard` and `sur` as derived outputs.
  })

  model({
    # Efficacy threshold on the model-predicted minimum aripiprazole
    # concentration proximate to the event, in ng/mL. Wang 2022 Results,
    # "Exposure-Response Analysis of Time to Relapse": the continuous
    # exposure-response fit was biased above this value, so two groups
    # of predicted Cmin (< 95 ng/mL vs >= 95 ng/mL) were modelled.
    cminThreshold <- 95

    # The (CMIN_ARI >= cminThreshold) comparison yields 0 / 1 in the
    # generated C, reproducing the Table 2 threshold indicator.
    aboveThreshold <- (CMIN_ARI >= cminThreshold)

    # Constant (exponential) hazard, per day.
    hazard <- exp(lhaz_base + e_cmin_haz * aboveThreshold)

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
