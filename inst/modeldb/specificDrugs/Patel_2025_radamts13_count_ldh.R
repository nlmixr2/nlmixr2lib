Patel_2025_radamts13_count_ldh <- function() {
  description <- "Exposure-response count model for elevated lactate dehydrogenase (LDH) events, a marker of microangiopathic hemolytic anemia (MAHA), in patients with congenital thrombotic thrombocytopenic purpura (cTTP) treated with recombinant ADAMTS13 (rADAMTS13; TAK-755) or plasma-based therapy. Poisson model with a between-subject random effect on the baseline count parameter and a log-linear drug effect driven by the average plasma ADAMTS13 activity over a dosing interval (CAV, IU/mL). Fitted to event counts within prophylaxis Periods 1 and 2 of the pivotal phase III crossover study (NCT03393975) in 41 patients of all ages. Elevated LDH is defined as LDH above 1.5 times the baseline value or above 1.5 times the upper limit of normal. This is the LDH counterpart of modellib('Patel_2025_radamts13_count_thrombocytopenia'); the authors retained a linear rather than sigmoid Emax drug effect for this endpoint, so the two count models have different structures and are not a single multi-output fit. Exposure comes from the companion population PK model modellib('Patel_2025_radamts13')."
  reference <- paste(
    "Patel M, Xu H, Barriere O, Diderichsen P, Patwari P, Zhu AZX,",
    "Marier JF, Peyret T, Wang LT, Mellgard B, Wang W, Bhattacharya I.",
    "Use of PopPK and E-R Analyses toward Explaining Causal Link Between",
    "ADAMTS13 in Recombinant vs. Plasma-Based Therapies and Clinical",
    "Effects in cTTP. Clin Pharmacol Ther. 2025;118(4):813-822.",
    "doi:10.1002/cpt.3720.",
    "Parameter values from Table S5 of the Supplementary Information.",
    "Exposure driver from the companion PopPK model in Table 1 of the same",
    "paper; see modellib('Patel_2025_radamts13').",
    sep = " "
  )
  vignette <- "Patel_2025_radamts13_exposure_response"

  # `etab1` carries the between-subject variability on a paper-mechanistic
  # count-model baseline rather than on a standard PK structural parameter.
  paper_specific_etas <- c("etab1")

  units <- list(
    time          = "n/a (the count is accumulated over a whole prophylaxis period, roughly 6 months; the model has no time argument and no ODE states)",
    dosing        = "n/a (no drug-dosing events; the drug input is the CAV data covariate, in IU/mL)",
    concentration = "IU/mL (the CAV covariate; 1 IU/mL = 100% of normal plasma ADAMTS13 activity)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Mean average plasma ADAMTS13 activity over a dosing interval for the patient in the prophylaxis period being scored (IU/mL).",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. One value per patient per prophylaxis period, not a time-varying column: the paper derives a Cave after each dose in each patient from the companion PopPK model's individual post hoc parameters, then averages those within Period 1 and within Period 2 (paper Methods, 'Exposure-response count modeling'). 1 IU/mL = 100% of normal ADAMTS13 activity. Enters this model log-linearly on the untransformed IU/mL scale, so the fitted slope of -5.12 is per IU/mL: across the observed Cave range of roughly 0.015-0.43 IU/mL it multiplies the expected count by between 0.93 and 0.11. The companion PopPK model modellib('Patel_2025_radamts13') returns Cc in IU/L, so a solve of that model must be divided by 1000 before it is used here; the validation vignette does this explicitly.",
      source_name        = "Cave"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "All ages (the all-ages final model, as for the thrombocytopenia count model).",
    weight_range   = "Not reported separately for the 41-patient exposure-response set; the parent PK analysis set spans 18.3-130.0 kg (Table 2).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported separately for the 41-patient exposure-response set.",
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP). Endpoint: elevated LDH, defined as LDH raised to more than 1.5 times the baseline value or more than 1.5 times the upper limit of normal, used as a marker of microangiopathic hemolytic anemia (paper Methods, 'Exposure-response analyses').",
    dose_range     = "rADAMTS13 40 IU/kg IV Q1W or Q2W; PBT approximately 10 IU/kg IV.",
    regions        = "Multinational phase III crossover study NCT03393975.",
    notes          = "Same 41-patient cohort, periods and exposure metric as modellib('Patel_2025_radamts13_count_thrombocytopenia'); only the endpoint and the retained drug-effect form differ. Baseline LDH in the parent PK analysis set was a median of 178 U/L (range 106-1,030) (Table 2)."
  )

  ini({
    # ---- Baseline count ----
    # Table S5 reports B1 on the log scale (estimate -0.886) alongside its
    # back-transform in the "Untransformed values of estimate" column (0.412 =
    # exp(-0.886)), which fixes the scale unambiguously. Because the reported
    # value is already a log, the canonical `l<name> <- log(<value>)` form
    # collapses to a bare assignment here.
    b1 <- -0.886; label("Log baseline expected elevated-LDH count per period at zero ADAMTS13 activity (log counts)")  # Table S5, Elevated LDH counts: B1 = -0.886 (RSE 65.3%), untransformed 0.412

    # ---- Log-linear drug effect ----
    # The authors evaluated linear, exponential and sigmoid Emax exposure-
    # response forms (paper Methods) and retained a slope for this endpoint;
    # Table S5 reports "Slope" rather than the Emax / EC50 / Gamma triplet
    # that the thrombocytopenia model carries. The slope multiplies CAV on the
    # untransformed IU/mL scale and acts on the log-count scale, which is
    # confirmed by Table S7: exp(-exp(b1 + slope * CAV)) reproduces every
    # published probability-of-zero-events cell for this endpoint.
    slope_cav <- -5.12; label("Effect of average ADAMTS13 activity on the log expected elevated-LDH count (per IU/mL)")  # Table S5, Elevated LDH counts: Slope = -5.12 (RSE 49.8%)

    # ---- Between-subject variability ----
    # Read as the NONMEM OMEGA variance on the log-count scale, for the same
    # reason given in the sibling thrombocytopenia count model: Table S5
    # carries no percent-CV footnote (unlike Table S6), and 1.81 as a percent
    # CV would be an implausibly narrow 1.81% between-subject spread.
    etab1 ~ 1.81  # Table S5, Elevated LDH counts: BSV on B1 = 1.81 (RSE 55.6%)
  })

  model({
    # 1. Expected count over the prophylaxis period. The drug effect is
    #    log-linear, so CAV = 0 collapses to the baseline count exp(b1).
    lambda <- exp(b1 + etab1 + slope_cav * CAV)

    # 2. Poisson count observation.
    ldh_elevation_count ~ pois(lambda)
  })
}
