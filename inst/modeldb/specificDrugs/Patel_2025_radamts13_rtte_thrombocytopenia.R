Patel_2025_radamts13_rtte_thrombocytopenia <- function() {
  description <- "Parametric repeated time-to-event (RTTE) hazard model for thrombocytopenia in patients with congenital thrombotic thrombocytopenic purpura (cTTP) treated with recombinant ADAMTS13 (rADAMTS13; TAK-755) or plasma-based therapy. Constant baseline hazard with a between-subject random effect, multiplied by a sigmoid Emax inhibition acting on the log-hazard scale and driven by the average plasma ADAMTS13 activity (CAV, IU/mL). Fitted to longitudinal event data from all three prophylaxis periods of the pivotal phase III crossover study (NCT03393975), all ages. Maximal inhibition is nearly complete: Emax of -3.71 on the log hazard is a 97.5% reduction, reached with ECave50 = 0.0113 IU/mL (about 1.1% of normal ADAMTS13 activity) and a Hill coefficient fixed at 2.58 from the count model. Forward simulation exposes `hazard` (instantaneous event rate per day) and `sur` (probability of remaining thrombocytopenia-event-free since t = 0). Exposure comes from the companion population PK model modellib('Patel_2025_radamts13')."
  reference <- paste(
    "Patel M, Xu H, Barriere O, Diderichsen P, Patwari P, Zhu AZX,",
    "Marier JF, Peyret T, Wang LT, Mellgard B, Wang W, Bhattacharya I.",
    "Use of PopPK and E-R Analyses toward Explaining Causal Link Between",
    "ADAMTS13 in Recombinant vs. Plasma-Based Therapies and Clinical",
    "Effects in cTTP. Clin Pharmacol Ther. 2025;118(4):813-822.",
    "doi:10.1002/cpt.3720.",
    "Parameter values from Table S6 of the Supplementary Information;",
    "model form from Supplementary Methods S1, 'Longitudinal repeated",
    "time-to-event exposure-response modeling'.",
    "Exposure driver from the companion PopPK model in Table 1 of the same",
    "paper; see modellib('Patel_2025_radamts13').",
    sep = " "
  )
  vignette <- "Patel_2025_radamts13_exposure_response"

  # `etallambda0` carries the between-subject variability on a paper-
  # mechanistic log-hazard baseline rather than on a PK structural parameter.
  paper_specific_etas <- c("etallambda0")

  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the drug input is the CAV data covariate, in IU/mL)",
    concentration = "probability (the model output `sur` is the thrombocytopenia-event-free survival probability, not a drug concentration; the CAV covariate is in IU/mL)"
  )

  compartmentData <- list(
    cumhaz = list(analyte = "thrombocytopenia event hazard", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    CAV = list(
      description        = "Average plasma ADAMTS13 activity driving the instantaneous thrombocytopenia hazard (IU/mL).",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Per-record, and may be time-varying: the source model is driven by 'time-varying ADAMTS13 activity' and specifically by a daily Cave (paper Methods, 'Longitudinal repeated time-to-event exposure-response modeling'). The paper's own reported simulations, however, hold it constant at the treatment-arm median Cave - Table S8b's event-free probabilities are reproduced exactly by this model with a constant CAV, which is how the validation vignette gates it. 1 IU/mL = 100% of normal ADAMTS13 activity. The companion PopPK model modellib('Patel_2025_radamts13') returns Cc in IU/L, so a solve of that model must be divided by 1000 before it is used here. Set to 0 for an untreated period, which collapses the hazard to the baseline lambda0.",
      source_name        = "Cave"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 1L,
    age_range      = "All ages (paper Results, 'Using data from Periods 1 to 3 and all age groups').",
    weight_range   = "Not reported separately for the RTTE set; the parent PK analysis set spans 18.3-130.0 kg (Table 2).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported separately for the RTTE set.",
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP). Endpoint: thrombocytopenia, defined as a platelet count decreased by at least 25% from baseline or a platelet count <150,000/uL (paper Methods, 'Exposure-response analyses').",
    dose_range     = "rADAMTS13 40 IU/kg IV Q1W or Q2W; PBT approximately 10 IU/kg IV.",
    regions        = "Multinational phase III crossover study NCT03393975.",
    notes          = "The RTTE analysis uses longitudinal data from all three prophylaxis periods, unlike the count and Cox analyses which use Periods 1 and 2 only. n_subjects is given as 43 because Figure 1 footnote a states that the count exposure-response cohort of 41 excluded two patients who WERE included in the repeated time-to-event modeling; the paper does not print an explicit RTTE N, so 43 is inferred from that footnote and should be treated as the paper's implied cohort rather than a directly reported figure. The paper reports the RTTE estimates with 95% confidence intervals rather than RSEs (Table S6)."
  )

  ini({
    # ---- Constant baseline hazard ----
    # Supplementary Methods S1 specifies a constant baseline on the log-hazard
    # scale, f_baseline(t) = beta0. Table S6 reports the back-transformed
    # baseline rate lambda0 = 0.0349 with a 95% CI, so it is log-transformed
    # here per library convention. The rate is per DAY: with the published
    # LDH-endpoint parameters, exp(-hazard * 183) and exp(-hazard * 365)
    # reproduce every Table S8b event-free probability at Month 6 and Month 12
    # to within 0.0008, which fixes the time unit unambiguously.
    llambda0 <- log(0.0349); label("Log baseline thrombocytopenia hazard at zero ADAMTS13 activity (log events/day)")  # Table S6, repeated time-to-event, Thrombocytopenia: Lambda0 = 0.0349 (95% CI 0.0108, 0.116)

    # ---- Sigmoid Emax drug effect, on the log-hazard scale ----
    # Emax is negative because it reduces the hazard. Its magnitude is
    # confirmed by the Results text: 1 - exp(-3.71) = 97.55%, quoted as "a
    # maximum 97.5% reduction in the hazard of thrombocytopenia".
    emax_haz <- -3.71;  label("Maximum drug effect on the log baseline hazard (log-hazard units)")               # Table S6, repeated time-to-event, Thrombocytopenia: Emax = -3.71 (95% CI -5.13, -2.41)
    ec50     <- 0.0113; label("Average ADAMTS13 activity giving half the maximum hazard reduction (IU/mL)")      # Table S6, repeated time-to-event, Thrombocytopenia: EC50 = 0.0113 IU/mL (95% CI 0.00911, 0.0141). Quoted in Results as "an ECave50 of 0.0113 IU/mL"

    # Gamma was carried over from the count model rather than re-estimated;
    # Table S6 prints it as "2.58, Fixed".
    gamma <- fixed(2.58); label("Hill coefficient of the ADAMTS13 activity-hazard relationship (unitless)")      # Table S6, repeated time-to-event, Thrombocytopenia: Gamma = 2.58, Fixed

    # ---- Between-subject variability on the baseline hazard ----
    # Table S6 reports IIV as a percent coefficient of variation and gives the
    # back-transform in its own note: "IIV is presented as a % coefficient of
    # variation (CV%) derived from the standard deviation of the random effect
    # (eta_i) as (100 x (exp(omega^2)-1)^0.5)". Inverting that note:
    #   omega^2 = log(1 + 1.41^2) = 1.0946377
    etallambda0 ~ 1.0946377  # Table S6, repeated time-to-event, Thrombocytopenia: IIV Lambda0 = 141% (95% CI 82.5, 178); converted via the table note

    # No residual error. As a parametric hazard model the likelihood is the
    # survival / event density itself, not an observation-error model;
    # forward simulation exposes `hazard` and `sur` as derived outputs.
  })

  model({
    # 1. Subject-specific baseline hazard at zero ADAMTS13 activity.
    lambda0 <- exp(llambda0 + etallambda0)

    # 2. Sigmoid Emax inhibition. CAV = 0 gives drug = 0 exactly
    #    (0^gamma = 0), collapsing the hazard to lambda0.
    drug <- CAV^gamma / (CAV^gamma + ec50^gamma)

    # 3. Instantaneous hazard, events/day. The drug term acts on the log-
    #    hazard scale per Supplementary Methods S1:
    #      log h(t) = log(lambda0) + emax_haz * drug
    hazard <- lambda0 * exp(emax_haz * drug)

    # 4. Cumulative hazard and event-free survival.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
