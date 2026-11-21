Du_2025_repotrectinib_pfs_ros1 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-efficacy model for progression-free",
    "survival (progression or death) on repotrectinib, a ROS1 / TRK / ALK",
    "tyrosine kinase inhibitor, in adults with ROS1-positive non-small cell",
    "lung cancer from the TRIDENT-1 trial (Du 2025, n = 127, 56 events).",
    "Exposure enters as the natural log of Cavgc, the TIME-VARYING",
    "cumulative half-daily average repotrectinib concentration from day 1",
    "to the current time, so dose interruptions, reductions and the day-15",
    "switch to twice-daily dosing are all carried. Hazard ratio 0.389 per",
    "e-fold increase in Cavgc. Time-fixed covariates: age, body weight,",
    "baseline tumor size, sex, Asian race, ECOG performance status and prior",
    "ROS1 TKI treatment (hazard ratio 4.48). The model returns the RELATIVE",
    "hazard only: a Cox regression leaves the baseline hazard unspecified,",
    "so no survivor function is offered. Continuous covariates are centred",
    "on the Figure 3A reference values so hr is relative to a reference",
    "patient. Cavgc is supplied as a data column, derived in the source",
    "analysis from the companion population PK model, packaged as",
    "Du_2025_repotrectinib. One of nine exposure-response models in the",
    "Du_2025_repotrectinib_* family."
  )
  reference <- paste(
    "Du S, Hu Z, Shen J, Zhu L, Roy A, Lam J, Lu M, Kondic A, Hamuro L.",
    "Exposure-Response Analysis of Repotrectinib to Support the Dose",
    "Recommendation for Patients With ROS1-Positive NSCLC or NTRK-Positive",
    "Solid Tumors. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):1993-2005.",
    "doi:10.1002/psp4.70102. Individual repotrectinib exposures derive from",
    "the companion population PK model (doi:10.1002/psp4.70036); see",
    "modellib('Du_2025_repotrectinib')."
  )
  vignette <- "Du_2025_repotrectinib_exposure_response"
  units <- list(
    time = "day (the source analysis time scale; the model itself has no time dimension, Cavgc carries it)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "hr (relative hazard of progression or death, unitless; also lhr)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cavgc: time-varying cumulative half-daily average repotrectinib plasma concentration from day 1 to the current time",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Du 2025 Methods, Analysis Data: 'Time-varying cumulative half-daily",
        "average concentration (Cavgc) from Day 1 to the event or censor',",
        "simulated from the empirical Bayes estimates of the companion popPK",
        "model and the actual dosing history (Figure S1). It is the running",
        "mean AUC(0, t) / t, updated every half day, and enters the",
        "counting-process Cox model Surv(START, STOP, EVENT) as a",
        "time-varying covariate, so supply it per half-day interval. NATURAL",
        "log scale (Table S9 LCAVGC; Table S10 HR 0.3892 = exp(-0.9438))."
      ),
      source_name = "LCAVGC (log Cavgc)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 3A reference 57 years (P5-P95 35-75); exp(0.003978 * (35 - 57)) = 0.916 reproduces the printed P5 hazard ratio.",
      source_name = "AGE2"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 3A reference 64 kg (P5-P95 45.1-89); exp(-0.01772 * (45.1 - 64)) = 1.40 reproduces the printed P5 hazard ratio.",
      source_name = "WTB"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported in cm by Du 2025 (TSIZEB); model() divides the canonical mm column by 10 and centres on the Figure 3A reference 4.9 cm (P5-P95 1.44-14.7). exp(0.08498 * (14.7 - 4.9)) = 2.30 reproduces the printed P95 hazard ratio.",
      source_name = "TSIZEB"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table S10 'Sex [Female: Male]', hazard ratio 0.4626: female patients have a lower hazard of progression or death. 81 female of 127 (Figure 3A).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White or other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or other race)",
      notes = "Table S10 'Race [Asian:Caucasian/Others]'. 75 Asian of 127 (Figure 3A).",
      source_name = "RACE2"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes = "Table S10 'Baseline Performance Status [1:0]'. 85 of 127 with PS 1 (Figure 3A).",
      source_name = "PS"
    ),
    PRIOR_TKI = list(
      description = "Prior treatment with a tyrosine kinase inhibitor targeting the tumor's driver fusion (here a ROS1 TKI); 1 = TKI-pretreated, 0 = TKI-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (TKI-naive)",
      notes = "Table S10 'TKI pretreatment [TKI-pretreated:TKI-naive]', hazard ratio 4.476. 56 of 127 pretreated (Figure 3A).",
      source_name = "PRVTRTC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 127L,
    n_studies = 1L,
    n_observations = "127 progression-free survival records, 56 events and 71 censored (Table S1)",
    age_range = "28-80 years (median 57)",
    weight_range = "39.5-119.3 kg (median 63.7)",
    sex_female_pct = 63.8,
    race_ethnicity = c(Asian = 59.1, `White or other` = 40.9),
    disease_state = "Locally advanced or metastatic ROS1-positive non-small cell lung cancer; 71 ROS1 TKI-naive and 56 pretreated with one prior ROS1 TKI and no chemotherapy or immunotherapy",
    dose_range = "Repotrectinib orally; 117 of 127 on 160 mg once daily for 14 days then 160 mg twice daily, the rest on phase 1 regimens of 40-240 mg once daily or 160 mg twice daily",
    regions = "Multinational (TRIDENT-1, NCT03093116)",
    notes = "Du 2025 Table S1 (ER-efficacy analysis dataset, PFS columns, data cutoff March 2023) and Figure 3A."
  )

  ini({
    # ==================================================================
    # Du 2025 Table S10, 'Parameter Estimates of the E-R of PFS (Full
    # Models)', ROS1-positive NSCLC block; the full model is
    # Surv(START, STOP, EVENT) ~ LCAVGC + AGE2 + WTB + TSIZEB + SEX +
    # RACE2 + PS + PRVTRTC (Table S9, log.full.coxph, lowest BIC).
    # Coefficients are on the log-hazard scale as printed.
    #
    # NO BASELINE HAZARD IS ENCODED. A semiparametric Cox model never
    # estimates h0(t), so it is not an unreported parameter. The model
    # returns the relative hazard hr only. Survival software, no variance
    # components: no IIV, no residual error, no observation endpoint.
    # ==================================================================
    e_cav_haz        <- -0.9438   ; label("Log hazard ratio per natural-log unit of Cavgc (unitless)")                     # Du 2025 Table S10 ROS1 'Repo Log Cavgc [ng/mL]' -0.9438, SE 0.4524; HR 0.3892 (0.1603, 0.9445)
    e_age_haz        <- 0.003978  ; label("Log hazard ratio per year of age (1/year)")                                      # Du 2025 Table S10 ROS1 'Age [yr]' 0.003978, SE 0.01158; HR 1.004
    e_wt_haz         <- -0.01772  ; label("Log hazard ratio per kg of baseline body weight (1/kg)")                         # Du 2025 Table S10 ROS1 'Body Weight [kg]' -0.01772, SE 0.01561; HR 0.9824
    e_tumsz_haz      <- 0.08498   ; label("Log hazard ratio per cm of baseline tumor size (1/cm)")                          # Du 2025 Table S10 ROS1 'Baseline Tumor Size [cm]' 0.08498, SE 0.0328; HR 1.089 (1.021, 1.161)
    e_sexf_haz       <- -0.7708   ; label("Log hazard ratio, female vs male (unitless)")                                    # Du 2025 Table S10 ROS1 'Sex [Female: Male]' -0.7708, SE 0.3555; HR 0.4626 (0.2305, 0.9286)
    e_race_asian_haz <- -0.4399   ; label("Log hazard ratio, Asian vs Caucasian or other (unitless)")                       # Du 2025 Table S10 ROS1 'Race [Asian:Caucasian/Others]' -0.4399, SE 0.3655; HR 0.6441
    e_ecog_ge1_haz   <- 0.2098    ; label("Log hazard ratio, ECOG performance status 1 vs 0 (unitless)")                    # Du 2025 Table S10 ROS1 'Baseline Performance Status [1:0]' 0.2098, SE 0.307; HR 1.233
    e_prior_tki_haz  <- 1.499     ; label("Log hazard ratio, ROS1-TKI-pretreated vs TKI-naive (unitless)")                  # Du 2025 Table S10 ROS1 'TKI pretreatment [TKI-pretreated:TKI-naive]' 1.499, SE 0.3176; HR 4.476 (2.402, 8.342)
  })

  model({
    # Relative hazard against a reference patient: the Figure 3A
    # reference values for the continuous covariates (age 57 y, weight
    # 64 kg, tumor size 4.9 cm), the reference level of each categorical
    # covariate, and Cavgc = 1 ng/mL (exposure is left uncentred because
    # the paper gives no reference Cavgc). Under proportional hazards the
    # subject hazard is h0(t) * hr for any baseline h0(t); only ratios of
    # hr are meaningful, and centring changes nothing but h0.
    lhr <- e_cav_haz * log(CAV) +
      e_age_haz * (AGE - 57) +
      e_wt_haz * (WT - 64) +
      e_tumsz_haz * (TUMSZ / 10 - 4.9) +
      e_sexf_haz * SEXF +
      e_race_asian_haz * RACE_ASIAN +
      e_ecog_ge1_haz * ECOG_GE1 +
      e_prior_tki_haz * PRIOR_TKI

    hr <- exp(lhr)
  })
}
