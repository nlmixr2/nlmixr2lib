Du_2025_repotrectinib_pfs_ntrk <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-efficacy model for progression-free",
    "survival (progression or death) on repotrectinib, a ROS1 / TRK / ALK",
    "tyrosine kinase inhibitor, in adults with NTRK-fusion-positive advanced",
    "solid tumors from the TRIDENT-1 trial (Du 2025, n = 88, 53 events).",
    "Exposure enters as the natural log of Cavgc, the TIME-VARYING",
    "cumulative half-daily average repotrectinib concentration from day 1",
    "to the current time. Hazard ratio 0.607 per e-fold increase in Cavgc,",
    "not statistically significant (95% CI 0.29-1.29). Time-fixed",
    "covariates: age, body weight, baseline tumor size, sex, Asian race,",
    "ECOG performance status, prior TRK TKI treatment (hazard ratio 3.43),",
    "NTRK3 versus NTRK1 or NTRK2 fusion (hazard ratio 0.311) and non-NSCLC",
    "tumor type. The model returns the RELATIVE hazard only: a Cox",
    "regression leaves the baseline hazard unspecified. Continuous",
    "covariates are centred on the Figure 3B reference values. Cavgc is",
    "supplied as a data column, derived in the source analysis from the",
    "companion population PK model, packaged as Du_2025_repotrectinib. One",
    "of nine exposure-response models in the Du_2025_repotrectinib_* family."
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
      notes = "Du 2025 Methods, Analysis Data: running mean AUC(0, t) / t updated every half day, simulated from the empirical Bayes estimates of the companion popPK model and the actual dosing history (Figure S1); a time-varying covariate in the counting-process Cox model, so supply it per half-day interval. NATURAL log scale (Table S9 LCAVGC; Table S10 HR 0.6073 = exp(-0.4988)).",
      source_name = "LCAVGC (log Cavgc)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 3B reference 60 years (P5-P95 27-79); exp(0.01592 * (27 - 60)) = 0.591 reproduces the printed P5 hazard ratio.",
      source_name = "AGE2"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 3B reference 70 kg (P5-P95 48.3-96.1).",
      source_name = "WTB"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported in cm by Du 2025 (TSIZEB); model() divides the canonical mm column by 10 and centres on the Figure 3B reference 5.7 cm (P5-P95 1.51-17). exp(0.05532 * (1.51 - 5.7)) = 0.793 reproduces the printed P5 hazard ratio.",
      source_name = "TSIZEB"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table S10 'Sex [Female: Male]'. 47 female of 88 (Figure 3B).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White or other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or other race)",
      notes = "Table S10 'Race [Asian:Caucasian/Others]'. 33 Asian of 88 (Figure 3B).",
      source_name = "RACE2"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes = "Table S10 'Baseline Performance Status [1:0]'. 51 of 88 with PS 1 (Figure 3B).",
      source_name = "PS"
    ),
    PRIOR_TKI = list(
      description = "Prior treatment with a tyrosine kinase inhibitor targeting the tumor's driver fusion (here a TRK TKI); 1 = TKI-pretreated, 0 = TKI-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (TKI-naive)",
      notes = "Table S10 'TKI pretreatment [TKI-pretreated:TKI-naive]', hazard ratio 3.427. 48 of 88 pretreated (Figure 3B).",
      source_name = "PRVTRTC"
    ),
    TUM_NTRK3_MUT = list(
      description = "NTRK3 gene fusion indicator within NTRK-fusion-positive tumors; 1 = NTRK3 fusion, 0 = NTRK1 or NTRK2 fusion",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NTRK1 or NTRK2 fusion)",
      notes = "Table S10 'Genetic Mutation [NTRK3: NTRK1-2]', hazard ratio 0.3112. 51 NTRK3 of 88 (Figure 3B).",
      source_name = "GENMTC"
    ),
    TUMTP_OTHER = list(
      description = "Non-NSCLC tumor type indicator; 1 = any solid tumor other than non-small cell lung cancer, 0 = NSCLC",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-small cell lung cancer, the only named tumor type in this model)",
      notes = "Table S10 'Tumor Type [Others:NSCLC]'. 53 non-NSCLC of 88 (Figure 3B).",
      source_name = "TTYPEF"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 1L,
    n_observations = "88 progression-free survival records, 53 events and 35 censored (Table S1)",
    age_range = "20-84 years (median 60)",
    weight_range = "42-116.3 kg (median 70)",
    sex_female_pct = 53.4,
    race_ethnicity = c(Asian = 37.5, `White or other` = 62.5),
    disease_state = "Locally advanced or metastatic NTRK1-3-fusion-positive solid tumors (53 of 88 non-NSCLC); 40 TRK TKI-naive and 48 TRK TKI-pretreated; 51 NTRK3 and 37 NTRK1 or NTRK2 fusions",
    dose_range = "Repotrectinib orally; 80 of 88 on 160 mg once daily for 14 days then 160 mg twice daily, the rest on 40 or 160 mg once daily or 160 mg twice daily",
    regions = "Multinational (TRIDENT-1, NCT03093116)",
    notes = "Du 2025 Table S1 (ER-efficacy analysis dataset, PFS columns, data cutoff March 2023) and Figure 3B."
  )

  ini({
    # ==================================================================
    # Du 2025 Table S10, NTRK-positive solid tumors block; full model
    # Surv(START, STOP, EVENT) ~ LCAVGC + AGE2 + WTB + TSIZEB + SEX +
    # RACE2 + PS + PRVTRTC + GENMTC + TTYPEF (Table S9, lowest BIC).
    # Log-hazard coefficients as printed. No baseline hazard (Cox model;
    # relative hazard only), no IIV, no residual error.
    # ==================================================================
    e_cav_haz           <- -0.4988   ; label("Log hazard ratio per natural-log unit of Cavgc (unitless)")             # Du 2025 Table S10 NTRK 'Repo Log Cavgc [ng/mL]' -0.4988, SE 0.3829; HR 0.6073 (0.2867, 1.286)
    e_age_haz           <- 0.01592   ; label("Log hazard ratio per year of age (1/year)")                              # Du 2025 Table S10 NTRK 'Age [yr]' 0.01592, SE 0.01126; HR 1.016
    e_wt_haz            <- -0.001836 ; label("Log hazard ratio per kg of baseline body weight (1/kg)")                 # Du 2025 Table S10 NTRK 'Body Weight [kg]' -0.001836, SE 0.01077; HR 0.9982
    e_tumsz_haz         <- 0.05532   ; label("Log hazard ratio per cm of baseline tumor size (1/cm)")                  # Du 2025 Table S10 NTRK 'Baseline Tumor Size (cm)' 0.05532, SE 0.02353; HR 1.057 (1.009, 1.107)
    e_sexf_haz          <- -0.02179  ; label("Log hazard ratio, female vs male (unitless)")                            # Du 2025 Table S10 NTRK 'Sex [Female: Male]' -0.02179, SE 0.3367; HR 0.9784
    e_race_asian_haz    <- -0.1071   ; label("Log hazard ratio, Asian vs Caucasian or other (unitless)")               # Du 2025 Table S10 NTRK 'Race [Asian:Caucasian/Others]' -0.1071, SE 0.3953; HR 0.8984
    e_ecog_ge1_haz      <- 0.1702    ; label("Log hazard ratio, ECOG performance status 1 vs 0 (unitless)")            # Du 2025 Table S10 NTRK 'Baseline Performance Status [1:0]' 0.1702, SE 0.309; HR 1.186
    e_prior_tki_haz     <- 1.232     ; label("Log hazard ratio, TRK-TKI-pretreated vs TKI-naive (unitless)")           # Du 2025 Table S10 NTRK 'TKI pretreatment [TKI-pretreated:TKI-naive]' 1.232, SE 0.4189; HR 3.427 (1.508, 7.788)
    e_tum_ntrk3_mut_haz <- -1.167    ; label("Log hazard ratio, NTRK3 vs NTRK1 or NTRK2 fusion (unitless)")            # Du 2025 Table S10 NTRK 'Genetic Mutation [NTRK3: NTRK1-2]' -1.167, SE 0.3491; HR 0.3112 (0.157, 0.6169)
    e_tumtp_other_haz   <- 0.1612    ; label("Log hazard ratio, non-NSCLC vs NSCLC tumor type (unitless)")             # Du 2025 Table S10 NTRK 'Tumor Type [Others:NSCLC]' 0.1612, SE 0.3693; HR 1.175
  })

  model({
    # Relative hazard against a reference patient: Figure 3B reference
    # values for the continuous covariates (age 60 y, weight 70 kg, tumor
    # size 5.7 cm), the reference level of each categorical covariate,
    # and Cavgc = 1 ng/mL (no reference Cavgc is given in the paper).
    lhr <- e_cav_haz * log(CAV) +
      e_age_haz * (AGE - 60) +
      e_wt_haz * (WT - 70) +
      e_tumsz_haz * (TUMSZ / 10 - 5.7) +
      e_sexf_haz * SEXF +
      e_race_asian_haz * RACE_ASIAN +
      e_ecog_ge1_haz * ECOG_GE1 +
      e_prior_tki_haz * PRIOR_TKI +
      e_tum_ntrk3_mut_haz * TUM_NTRK3_MUT +
      e_tumtp_other_haz * TUMTP_OTHER

    hr <- exp(lhr)
  })
}
