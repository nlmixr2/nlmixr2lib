Yoshioka_2018_FXa_inhibitors_mbma <- function() {
  description <- "MBMA. PT-ratio-driven logistic event-rate model for direct oral factor Xa inhibitors (rivaroxaban, apixaban, edoxaban) in non-valvular atrial fibrillation. Inputs a population-mean prothrombin-time ratio (PTR) supplied per observation time; outputs per-arm probability of ischemic stroke/SE (p_isse) and of major bleeding (p_mb), plus a derived per-arm mortality probability. Fit by NONMEM 7.3 to per-arm event counts from 5 large RCTs (Yoshioka 2018; 57,655 patients). Suitable for simulating per-arm summary outcomes only; the upstream popPK -> PT-ratio layer for each FXa inhibitor is out of scope and PTR must be supplied externally."
  reference <- "Yoshioka H, Sato H, Hatakeyama H, Hisaka A. Model-based meta-analysis to evaluate optimal doses of direct oral factor Xa inhibitors in atrial fibrillation patients. Blood Adv. 2018;2(10):1066-1076. doi:10.1182/bloodadvances.2017013805."
  vignette <- "Yoshioka_2018_FXa_inhibitors_mbma"
  # Algebraic MBMA model: no rxode2 dose events are consumed (the input is the
  # PTR covariate column), and the model outputs p_isse / p_mb / p_death which
  # are dimensionless event probabilities in [0, 1] rather than drug
  # concentrations. The `units` entries below carry placeholder strings
  # ("probability" / "probability/probability") chosen so that
  # checkModelConventions() sees a dimensionally consistent dose-vs-concentration
  # pair (both "probability"); a more descriptive label would trip the
  # dimensional-compatibility check.
  units <- list(
    time          = "day",
    dosing        = "probability",
    concentration = "probability/probability"
  )

  covariateData <- list(
    PTR = list(
      description        = "Population-mean prothrombin time ratio (PT relative to pre-treatment baseline; unitless). Time-varying input that drives the FXa-inhibitor dose-response logistic curves.",
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "PTR = 1 reproduces placebo log-odds; PTR > 1 reflects FXa-inhibitor anticoagulation. Externally computed from an upstream popPK -> PT-ratio model for the FXa inhibitor of interest (Girgis 2014 ref 21 for rivaroxaban, Leil 2014 ref 22 and Chang 2016 ref 24 for apixaban, Krekels 2016 ref 23 and Koretsune 2015 ref 25 for edoxaban). PT values are corrected to RecombiplasTin reagent equivalence per Gosselin 2016 (ref 26) before computing the ratio. Must be > 0 because the major-bleeding equation evaluates log(PTR).",
      source_name        = "PT ratio (symbol x in Yoshioka 2018 Eq. 1 and Eq. 2)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 57655L,
    n_studies        = 5L,
    age_range        = "Adults with non-valvular atrial fibrillation; per-arm mean ages reported per trial.",
    weight_range     = "Per-arm baseline body-weight statistics reported per trial (Yoshioka 2018 Table 1); not used as a model covariate.",
    sex_female_pct   = NA_real_,
    disease_state    = "Non-valvular atrial fibrillation; primary efficacy endpoint is ischemic stroke / systemic embolism; primary safety endpoint is major bleeding.",
    dose_range       = "Pooled FXa-inhibitor RCT arms: rivaroxaban 15 mg OD (J-ROCKET AF) or 20 mg OD (ROCKET AF); apixaban 5 mg BID (ARISTOTLE, AVERROES); edoxaban 30 mg OD (low-dose) or 60 mg OD (high-dose) (ENGAGE AF-TIMI 48). Comparator arms: dose-adjusted warfarin in 4 trials and aspirin in AVERROES.",
    regions          = "Multinational; 5 RCTs -- ROCKET AF, J-ROCKET AF, ARISTOTLE, AVERROES, ENGAGE AF-TIMI 48 (Yoshioka 2018 Table 1).",
    trials_included  = "ROCKET AF (rivaroxaban vs dose-adjusted warfarin); J-ROCKET AF (rivaroxaban vs dose-adjusted warfarin in Japanese patients); ARISTOTLE (apixaban vs dose-adjusted warfarin); AVERROES (apixaban vs aspirin); ENGAGE AF-TIMI 48 (edoxaban high- and low-dose vs dose-adjusted warfarin).",
    notes            = "Summary-level MBMA: per-arm event counts and exposure-time were the modelled observations (no individual-patient data). The model predicts per-arm population-mean outcomes (probability of an event per patient-year of exposure), NOT individual concentrations or individual events. The only random effect is a between-STUDY variance on the placebo log-odds of ischemic stroke/SE (variance 0.00529; Table 2); no between-subject variability is estimated. Per-trial demographics and outcome rates are in Yoshioka 2018 Table 1."
  )

  ini({
    # Placebo (no-anticoagulation) log-odds intercepts on the logit scale
    # (Yoshioka 2018 Table 2; both estimated).
    e0_isse <- -2.93    ; label("Placebo log-odds of ischemic stroke / systemic embolism (E0_ISSE; logit scale)")   # Yoshioka 2018 Table 2
    e0_mb   <- -4.17    ; label("Placebo log-odds of major bleeding (E0_MB; logit scale)")                          # Yoshioka 2018 Table 2

    # FXa-inhibitor PT-ratio dose-response coefficients
    #   logit(P_ISSE) = E0_ISSE * exp((theta1 / theta2) * (exp(theta2 * (PTR - 1)) - 1))   -- Eq. 1
    #   logit(P_MB)   = E0_MB   + theta3 * log(PTR)                                         -- Eq. 2
    # At PTR = 1 both equations reduce to the placebo log-odds intercepts.
    theta1 <-   6.96    ; label("Slope exponent of PT-ratio dose-response for ischemic stroke/SE (unitless)")        # Yoshioka 2018 Table 2 (Eq. 1)
    theta2 <- -14.1     ; label("Steepness exponent of PT-ratio dose-response for ischemic stroke/SE (unitless)")    # Yoshioka 2018 Table 2 (Eq. 1)
    theta3 <-   1.92    ; label("Change in log-odds of major bleeding per 1-unit change in log(PT ratio)")           # Yoshioka 2018 Table 2 (Eq. 2)

    # Mortality-conversion weights (Eq. 3): the conditional probability of death
    # given the event. Derived in the source paper from the observed trial
    # mortality rates and held fixed when computing the optimal-dose mortality
    # outputs; encoded as fixed() because they are NOT estimated.
    w_isse <- fixed(0.23) ; label("Conditional mortality probability given an ischemic-stroke/SE event (fraction)")  # Yoshioka 2018 Methods, simulation-for-dose-optimization paragraph (Eq. 3 weight w1; "the weighting coefficient ... was set to 0.23 ... based on the mortality rates observed in the trials")
    w_mb   <- fixed(0.07) ; label("Conditional mortality probability given a major-bleeding event (fraction)")        # Yoshioka 2018 Methods, simulation-for-dose-optimization paragraph (Eq. 3 weight w2; "0.07 for major bleeding")

    # Between-STUDY (intertrial) random effect on the placebo ISSE log-odds.
    # MBMA convention: this is study-level variance (one draw per trial), NOT
    # individual between-subject variability. Named eta_study_* to mark the
    # scope difference from popPK etalcl / etalvc patterns. Per Yoshioka 2018
    # Results ("a statistically significant intertrial variability was found
    # only in the placebo response of ischemic stroke/SE"), this is the only
    # random effect in the published final model -- no other random effect
    # was retained; no covariate explained the remaining inter-trial variance.
    eta_study_e0_isse ~ 0.00529   # Yoshioka 2018 Table 2 (omega^2 = 0.00529; RSE 36.0%; 95% CI 0.00155 - 0.00902)
  })

  model({
    # Trial-arm placebo log-odds (apply the only random effect: study-level on ISSE intercept)
    plc_isse <- e0_isse + eta_study_e0_isse
    plc_mb   <- e0_mb

    # FXa-inhibitor dose-response on the logit scale
    fxa_mult_isse <- exp((theta1 / theta2) * (exp(theta2 * (PTR - 1)) - 1))
    logit_p_isse  <- plc_isse * fxa_mult_isse
    logit_p_mb    <- plc_mb + theta3 * log(PTR)

    # Event probabilities on the [0, 1] scale
    p_isse <- 1 / (1 + exp(-logit_p_isse))
    p_mb   <- 1 / (1 + exp(-logit_p_mb))

    # Derived mortality probabilities per Yoshioka 2018 Eq. 3:
    # P_death = w_isse * P_ISSE + w_mb * P_MB
    p_death_isse <- w_isse * p_isse
    p_death_mb   <- w_mb   * p_mb
    p_death      <- p_death_isse + p_death_mb
  })
}
