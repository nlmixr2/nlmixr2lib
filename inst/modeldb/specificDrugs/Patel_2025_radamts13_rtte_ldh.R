Patel_2025_radamts13_rtte_ldh <- function() {
  description <- "Parametric repeated time-to-event (RTTE) hazard model for elevated lactate dehydrogenase (LDH), a marker of microangiopathic hemolytic anemia (MAHA), in patients with congenital thrombotic thrombocytopenic purpura (cTTP) treated with recombinant ADAMTS13 (rADAMTS13; TAK-755) or plasma-based therapy. Constant baseline hazard with a between-subject random effect, multiplied by a sigmoid Emax inhibition acting on the log-hazard scale and driven by the average plasma ADAMTS13 activity (CAV, IU/mL). Fitted to longitudinal event data from the pivotal phase III crossover study (NCT03393975), all ages. Emax of -2.84 on the log hazard is a 94.2% maximum reduction, with ECave50 = 0.0133 IU/mL and a Hill coefficient fixed at 2.58. This is the LDH counterpart of modellib('Patel_2025_radamts13_rtte_thrombocytopenia'); the two share a structure but were fitted separately, with their own baseline hazards, Emax, EC50 and between-subject variability. Forward simulation exposes `hazard` (instantaneous event rate per day) and `sur` (probability of remaining elevated-LDH-event-free since t = 0). Exposure comes from the companion population PK model modellib('Patel_2025_radamts13')."
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
    concentration = "probability (the model output `sur` is the elevated-LDH-event-free survival probability, not a drug concentration; the CAV covariate is in IU/mL)"
  )

  compartmentData <- list(
    cumhaz = list(analyte = "elevated LDH event hazard", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    CAV = list(
      description        = "Average plasma ADAMTS13 activity driving the instantaneous elevated-LDH hazard (IU/mL).",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Per-record, and may be time-varying (paper Methods, 'Longitudinal repeated time-to-event exposure-response modeling'), but the paper's own reported simulations hold it constant at the treatment-arm median Cave: with a constant CAV this model reproduces all eight published Month-6 and Month-12 event-free probabilities in Table S8b to within 0.0008, which is how the validation vignette gates it. 1 IU/mL = 100% of normal ADAMTS13 activity. The companion PopPK model modellib('Patel_2025_radamts13') returns Cc in IU/L, so a solve of that model must be divided by 1000 before it is used here. Set to 0 for an untreated period, which collapses the hazard to the baseline lambda0. Table S8a median Cave values used for the published simulations: 0.0291 / 0.0582 IU/mL for PBT Q2W / Q1W and 0.192 / 0.384 IU/mL for rADAMTS13 Q2W / Q1W.",
      source_name        = "Cave"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 1L,
    age_range      = "All ages.",
    weight_range   = "Not reported separately for the RTTE set; the parent PK analysis set spans 18.3-130.0 kg (Table 2).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported separately for the RTTE set.",
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP). Endpoint: elevated LDH, defined as LDH raised to more than 1.5 times the baseline value or more than 1.5 times the upper limit of normal, used as a marker of microangiopathic hemolytic anemia (paper Methods, 'Exposure-response analyses').",
    dose_range     = "rADAMTS13 40 IU/kg IV Q1W or Q2W; PBT approximately 10 IU/kg IV.",
    regions        = "Multinational phase III crossover study NCT03393975.",
    notes          = "Same cohort and structure as modellib('Patel_2025_radamts13_rtte_thrombocytopenia'); only the endpoint and the fitted values differ. As for that model, n_subjects = 43 is inferred from Figure 1 footnote a (two patients in the RTTE analysis were not part of the 41-patient count cohort) rather than directly printed. Baseline LDH in the parent PK analysis set was a median of 178 U/L (range 106-1,030) (Table 2)."
  )

  ini({
    # ---- Constant baseline hazard ----
    # Per Supplementary Methods S1, f_baseline(t) = beta0. Table S6 reports
    # the back-transformed rate lambda0 = 0.00923 with a 95% CI, so it is
    # log-transformed here per library convention. The rate is per DAY:
    # exp(-hazard * 183) and exp(-hazard * 365) reproduce all eight Table S8b
    # Month-6 and Month-12 event-free probabilities to within 0.0008.
    llambda0 <- log(0.00923); label("Log baseline elevated-LDH hazard at zero ADAMTS13 activity (log events/day)")  # Table S6, repeated time-to-event, Elevated LDH: Lambda0 = 0.00923 (95% CI 0.00270, 0.0296)

    # ---- Sigmoid Emax drug effect, on the log-hazard scale ----
    # Emax is negative because it reduces the hazard; 1 - exp(-2.84) = 94.2%
    # maximum reduction. The paper quotes only the thrombocytopenia figure
    # (97.5%) in the Results text and refers to Table S8b for LDH.
    emax_haz <- -2.84;  label("Maximum drug effect on the log baseline hazard (log-hazard units)")           # Table S6, repeated time-to-event, Elevated LDH: Emax = -2.84 (95% CI -4.33, -1.34)
    ec50     <- 0.0133; label("Average ADAMTS13 activity giving half the maximum hazard reduction (IU/mL)")  # Table S6, repeated time-to-event, Elevated LDH: EC50 = 0.0133 IU/mL (95% CI 0.0100, 0.0180)

    gamma <- fixed(2.58); label("Hill coefficient of the ADAMTS13 activity-hazard relationship (unitless)")  # Table S6, repeated time-to-event, Elevated LDH: Gamma = 2.58, Fixed

    # ---- Between-subject variability on the baseline hazard ----
    # Converted from the percent CV using the Table S6 note,
    # CV% = 100 x (exp(omega^2) - 1)^0.5, so
    #   omega^2 = log(1 + 1.21^2) = 0.9018266
    etallambda0 ~ 0.9018266  # Table S6, repeated time-to-event, Elevated LDH: IIV Lambda0 = 121% (95% CI 0, 166); converted via the table note

    # No residual error; the likelihood is the survival / event density itself.
  })

  model({
    # 1. Subject-specific baseline hazard at zero ADAMTS13 activity.
    lambda0 <- exp(llambda0 + etallambda0)

    # 2. Sigmoid Emax inhibition. CAV = 0 gives drug = 0 exactly.
    drug <- CAV^gamma / (CAV^gamma + ec50^gamma)

    # 3. Instantaneous hazard, events/day, with the drug term on the
    #    log-hazard scale per Supplementary Methods S1.
    hazard <- lambda0 * exp(emax_haz * drug)

    # 4. Cumulative hazard and event-free survival.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
