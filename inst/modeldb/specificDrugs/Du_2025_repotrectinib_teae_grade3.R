Du_2025_repotrectinib_teae_grade3 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-safety model for time to the first grade 3 or higher treatment-emergent adverse event",
    "on repotrectinib, a ROS1 / TRK / ALK tyrosine kinase inhibitor, in 502",
    "adults with ALK-, ROS1- or NTRK1-3-positive advanced solid tumors from",
    "the TRIDENT-1 trial (Du 2025, pooled across tumor types). Exposure",
    "enters LINEARLY as Cavgc, the time-varying cumulative half-daily",
    "average repotrectinib concentration from day 1 to the current time.",
    "Exposure is NOT a statistically significant predictor (95% CI of the hazard ratio includes 1); ECOG performance status 0 (hazard ratio 0.568), Asian race (0.671 versus White) and male sex (1.34) are.",
    "Time-fixed covariates: age, body weight, sex, race (Asian, other versus",
    "White), ECOG performance status, prior TKI treatment, driver alteration",
    "(NTRK, ALK versus ROS1) and non-NSCLC tumor type, each coded with the",
    "source's reference level. The model returns the RELATIVE hazard only:",
    "a Cox regression leaves the baseline hazard unspecified. Age and body",
    "weight are centred on the Figure 5 reference values (56 years, 69 kg).",
    "Cavgc is supplied as a data column, derived in the source analysis from",
    "the companion population PK model, packaged as Du_2025_repotrectinib.",
    "One of nine exposure-response models in the Du_2025_repotrectinib_*",
    "family."
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
    concentration = "hr (relative hazard of a grade >= 3 treatment-emergent adverse event, unitless; also lhr)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cavgc: time-varying cumulative half-daily average repotrectinib plasma concentration from day 1 to the current time",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Du 2025 Methods, Analysis Data: running mean AUC(0, t) / t updated every half day, simulated from the empirical Bayes estimates of the companion popPK model and the actual dosing history (Figure S1); a time-varying covariate in the counting-process Cox model Surv(START, STOP, EVENT), so supply it per half-day interval. Enters LINEARLY (Table S13 model er.*.full.coxph with CAVGC, the lowest-BIC functional form for every safety endpoint).",
      source_name = "CAVGC"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 5 reference 56 years (P5-P95 32-75); the safety set median is 56, range 18-93 (Table S2).",
      source_name = "AGE2"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred on the Figure 5 reference 69 kg (P5-P95 47.6-105); the safety set median is 69, range 39.5-169.4 (Table S2).",
      source_name = "WTB"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female; the source codes Sex [Male:Female], so model() applies the printed coefficient to 1 - SEXF)",
      notes = "Du 2025 Table S14 'Sex [Male:Female]'. The printed coefficient is kept verbatim and multiplies the male indicator 1 - SEXF.",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = not Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_OTHER is also 0)",
      notes = "Du 2025 Table S14 'Race [Asian:White]'; paired with RACE_OTHER so both = 0 selects the White reference. Safety set 236 Asian of 502 (Table S2).",
      source_name = "RACEM"
    ),
    RACE_OTHER = list(
      description = "Other-race indicator (1 = neither White nor Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN is also 0)",
      notes = "Du 2025 Table S14 'Race [Others:White]'. Safety set 49 of 502 (Table S2).",
      source_name = "RACEM"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (ECOG performance status 1; the source codes Baseline Performance Status [0:1], so model() applies the printed coefficient to 1 - ECOG_GE1)",
      notes = "Du 2025 Table S14 'Baseline Performance Status [0:1]'; Figure 5 labels the levels 'Fully active' (PS 0) and 'Restricted' (PS 1). The printed coefficient is kept verbatim and multiplies the PS-0 indicator 1 - ECOG_GE1.",
      source_name = "PSM"
    ),
    PRIOR_TKI = list(
      description = "Prior treatment with a tyrosine kinase inhibitor targeting the tumor's driver alteration; 1 = TKI-pretreated, 0 = TKI-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (TKI-pretreated; the source codes TKI-pretreatment [TKI-Naive:TKI pretreated], so model() applies the printed coefficient to 1 - PRIOR_TKI)",
      notes = "Du 2025 Table S14 'TKI-pretreatment [TKI-Naive:TKI pretreated]'. Safety set 330 pretreated of 502 (Table S2). The printed coefficient is kept verbatim and multiplies the TKI-naive indicator 1 - PRIOR_TKI.",
      source_name = "PRVTRTM"
    ),
    TUM_NTRK_MUT = list(
      description = "NTRK1-3 gene fusion-positive tumor indicator; 1 = NTRK-fusion-positive, 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ROS1-positive when TUM_ALK_MUT is also 0)",
      notes = "Du 2025 Table S14 'Genetic mutation [NTRK:ROS1]'; paired with TUM_ALK_MUT so both = 0 selects the ROS1-positive reference. Safety set 113 NTRK of 502 (Table S2).",
      source_name = "GENMTM"
    ),
    TUM_ALK_MUT = list(
      description = "ALK rearrangement-positive tumor indicator; 1 = ALK-positive, 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ROS1-positive when TUM_NTRK_MUT is also 0)",
      notes = "Du 2025 Table S14 'Genetic mutation [ALK:ROS1]'. Safety set 30 ALK of 502 (Table S2).",
      source_name = "GENMTM"
    ),
    TUMTP_OTHER = list(
      description = "Non-NSCLC tumor type indicator; 1 = any solid tumor other than non-small cell lung cancer, 0 = NSCLC",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-small cell lung cancer, the only named tumor type in this model)",
      notes = "Du 2025 Table S14 'Tumor type [Others:NSCLC]'. Safety set 79 non-NSCLC of 502 (Table S2).",
      source_name = "CANTYPM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 502L,
    n_studies = 1L,
    n_observations = "502 time-to-first-event records in counting-process format (event counts not reported)",
    age_range = "18-93 years (median 56)",
    weight_range = "39.5-169.4 kg (median 69)",
    sex_female_pct = 56.5,
    race_ethnicity = c(White = 43.2, Asian = 47.0, Other = 9.8),
    disease_state = "Locally advanced or metastatic ROS1-positive (359), NTRK1-3-positive (113) or ALK-positive (30) solid tumors; 423 NSCLC and 79 other tumor types; 330 TKI-pretreated",
    dose_range = "Repotrectinib orally; 422 of 502 on 160 mg once daily for 14 days then 160 mg twice daily, the rest on phase 1 regimens of 40-240 mg once daily or 160-200 mg twice daily",
    regions = "Multinational (TRIDENT-1, NCT03093116)",
    notes = paste(
      "Du 2025 Table S2 (ER-safety analysis dataset) and Figure 5. Table S2",
      "prints 287 male / 215 female and 326 PS 0 / 176 PS 1, while Figure 5",
      "gives Male:Female 213:277 and Fully active:Restricted 170:320; both",
      "efficacy subsets are majority female and majority PS 1 (Table S1), so",
      "the Table S2 sex and performance-status rows are taken as transposed",
      "and the Figure 5 counts are used here (female 277/490 = 56.5%). See",
      "the vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Du 2025 Table S14, 'Parameter Estimates of the E-R of Safety
    # endpoints (Full models)', Gr3+ TEAEs block; full model
    # Surv(START, STOP, EVENT) ~ CAVGC + AGE2 + WTB + SEX + RACEM + PSM +
    # PRVTRTM + GENMTM + CANTYPM (Table S13, lowest-BIC parametric form).
    # Log-hazard coefficients as printed; each categorical effect keeps
    # the source's [comparator:reference] orientation, so model()
    # multiplies it by the comparator indicator.
    #
    # No baseline hazard (Cox model; relative hazard only), no IIV, no
    # residual error, no observation endpoint.
    # ==================================================================
    e_cav_haz          <- 0.0001849 ; label("Log hazard ratio per ng/mL of Cavgc (mL/ng)")  # Du 2025 Table S14 Gr3+ TEAEs 'Repo cumulative half- daily Cavg [ng/mL]' 0.0001849, SE 0.0005115; HR 1 (0.9992, 1.001)
    e_age_haz          <- 0.009629 ; label("Log hazard ratio per year of age (1/year)")  # Du 2025 Table S14 Gr3+ TEAEs 'Age [yr]' 0.009629, SE 0.005455; HR 1.01 (0.9989, 1.021)
    e_wt_haz           <- 0.006839 ; label("Log hazard ratio per kg of baseline body weight (1/kg)")  # Du 2025 Table S14 Gr3+ TEAEs 'Body Weight [kg]' 0.006839, SE 0.004668; HR 1.007 (0.9977, 1.016)
    e_male_haz         <- 0.293 ; label("Log hazard ratio, male vs female (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Sex [Male:Female]' 0.293, SE 0.1462; HR 1.34 (1.007, 1.785)
    e_race_asian_haz   <- -0.3985 ; label("Log hazard ratio, Asian vs White (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Race [Asian:White]' -0.3985, SE 0.171; HR 0.6713 (0.4801, 0.9387)
    e_race_other_haz   <- -0.2611 ; label("Log hazard ratio, other race vs White (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Race [Others:White]' -0.2611, SE 0.2379; HR 0.7702 (0.4831, 1.228)
    e_ecog0_haz        <- -0.5658 ; label("Log hazard ratio, ECOG performance status 0 vs 1 (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Baseline Performance Status [0:1]' -0.5658, SE 0.1502; HR 0.5679 (0.4231, 0.7623)
    e_tki_naive_haz    <- 0.08432 ; label("Log hazard ratio, TKI-naive vs TKI-pretreated (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'TKI-pretreatment [TKI-Naive:TKI pretreated]' 0.08432, SE 0.1471; HR 1.088 (0.8155, 1.451)
    e_tum_ntrk_mut_haz <- 0.1689 ; label("Log hazard ratio, NTRK-fusion-positive vs ROS1-positive (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Genetic mutation [NTRK:ROS1]' 0.1689, SE 0.2357; HR 1.184 (0.7459, 1.879)
    e_tum_alk_mut_haz  <- 0.4193 ; label("Log hazard ratio, ALK-positive vs ROS1-positive (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Genetic mutation [ALK:ROS1]' 0.4193, SE 0.2922; HR 1.521 (0.8577, 2.697)
    e_tumtp_other_haz  <- -0.01987 ; label("Log hazard ratio, non-NSCLC vs NSCLC tumor type (unitless)")  # Du 2025 Table S14 Gr3+ TEAEs 'Tumor type [Others:NSCLC]' -0.01987, SE 0.2645; HR 0.9803 (0.5837, 1.646)
  })

  model({
    # Relative hazard against the source's reference patient: Figure 5
    # reference age 56 y and weight 69 kg, female, White, ECOG PS 1,
    # TKI-pretreated, ROS1-positive NSCLC, at Cavgc = 0 ng/mL.
    isMale    <- 1 - SEXF
    isEcog0   <- 1 - ECOG_GE1
    isTkiNaive <- 1 - PRIOR_TKI

    lhr <- e_cav_haz * CAV +
      e_age_haz * (AGE - 56) +
      e_wt_haz * (WT - 69) +
      e_male_haz * isMale +
      e_race_asian_haz * RACE_ASIAN +
      e_race_other_haz * RACE_OTHER +
      e_ecog0_haz * isEcog0 +
      e_tki_naive_haz * isTkiNaive +
      e_tum_ntrk_mut_haz * TUM_NTRK_MUT +
      e_tum_alk_mut_haz * TUM_ALK_MUT +
      e_tumtp_other_haz * TUMTP_OTHER

    hr <- exp(lhr)
  })
}
