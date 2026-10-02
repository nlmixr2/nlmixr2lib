Yu_2022_igaNephropathy_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the time course of the change from ",
    "baseline in daily urinary protein excretion (g/day) in adults with IgA ",
    "nephropathy, fit to study-arm-level summary data from 40 clinical trials ",
    "(83 arms, 2288 participants) comparing placebo with six drug classes ",
    "grouped by pharmacological mechanism: corticosteroids, ",
    "immunosuppressants, renin-angiotensin system (RAS) blockers, antiplatelet ",
    "agents, N-3 fatty acids and 'other drugs' (agents outside the first five ",
    "classes and cross-class combinations). Every arm follows an Emax-in-time ",
    "model E(t) = Emax * t / (ET50 + t). Placebo arms have their own Emax ",
    "(-0.44 g/day) and a slow onset (ET50 = 27.2 months); the six drug classes ",
    "have class-specific Emax values and share a single ET50 of 5.59 months. ",
    "Arm-mean baseline urinary protein excretion is the only retained ",
    "covariate: it acts linearly on the drug-arm Emax (-0.63 g/day per 1 g/day ",
    "of baseline above the 1.82 g/day centring value) and not on placebo. ",
    "Between-STUDY-ARM (not between-subject) variability is carried on ET50 ",
    "only, as the scale form ET50_i = ET50 * (1 + eta); the between-arm eta on ",
    "Emax was estimated near zero and fixed to 0. The residual is additive ",
    "at unit study weight and the paper weights it by 1/sqrt(N) for an arm of ",
    "N participants. Suitable simulation scope is arm-mean proteinuria ",
    "time courses; the model is NOT suitable for individual-patient ",
    "simulation. Parameter values are Table 2 (NONMEM 7.4)."
  )

  reference <- paste(
    "Yu J, Luo J, Zhu H, Sui Z, Liu H, Li L, Zheng Q.",
    "Quantitative Comparison of the Clinical Efficacy of 6 Classes Drugs for",
    "IgA Nephropathy: A Model-Based Meta-Analysis of Drugs for Clinical",
    "Treatments. Front Immunol. 2022;13:825677.",
    "doi:10.3389/fimmu.2022.825677.",
    sep = " "
  )
  vignette <- "Yu_2022_igaNephropathy_mbma"

  # Yu 2022 Equation 2 places the between-arm random effect on ET50 in the
  # "scale" form P_i = P_typical * (1 + eta_i). The eta therefore has no
  # log-scale typical-value partner (the typical values are carried as
  # let50_placebo / let50_drug, one per arm type, and the same eta scales
  # both), which is the documented `paper_specific_etas` case.
  paper_specific_etas <- c("eta_study_et50")

  units <- list(
    time = "month (time since the start of treatment; ET50 is reported in months)",
    dosing = paste0(
      "(no dose events; each arm's treatment is identified by the drug-class ",
      "indicators TRT_*, and the model carries no dose or exposure term)"
    ),
    concentration = paste0(
      "g/day (change from baseline in daily urinary protein excretion; ",
      "negative values are a REDUCTION in proteinuria, i.e. benefit. The ",
      "observation is not a drug concentration, so the dosing string is ",
      "parenthesised to skip the dimensional check)"
    )
  )

  covariateData <- list(
    UPRO_BL = list(
      description = paste0(
        "Arm-level (study-arm median or mean) baseline daily urinary protein ",
        "excretion before treatment."
      ),
      units = "g/day",
      type = "continuous",
      reference_category = paste0(
        "n/a -- enters linearly as (UPRO_BL - 1.82) on the drug-arm Emax. ",
        "1.82 g/day is the centring value printed in Yu 2022 Equation 7."
      ),
      notes = paste0(
        "MBMA study-arm-level covariate: the published arm's baseline, not ",
        "an individual patient's. Range across the analysis 0.57-5.29 g/day ",
        "(Yu 2022 Table 1; overall median 1.9 g/day, placebo arms 1.86, ",
        "corticosteroids 1.6, immunosuppressants 2.77, RAS blockers 1.72, ",
        "antiplatelet agents 0.83, N-3 fatty acids 1.79, other drugs 2.1). ",
        "Acts on the drug-arm Emax only; Yu 2022 found no effect on the ",
        "placebo Emax (Results 'Model Building and Evaluation'). A higher ",
        "baseline gives a larger (more negative) drug effect: -0.63 g/day of ",
        "extra reduction per 1 g/day of baseline. The paper's simulations use ",
        "1.80 g/day (mild-to-moderate proteinuria) and 3.85 g/day (severe ",
        "proteinuria), stated as medians of the included studies."
      ),
      source_name = "Baseline daily urinary protein excretion (Yu 2022 Equation 7 'Baseline')"
    ),
    TRT_CORTICOSTEROID = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received a corticosteroid, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a different drug class).",
      notes = paste0(
        "MBMA study-arm-level indicator. Selects emax_corticosteroid = -1.47 ",
        "g/day. 8 trials / 9 arms / 288 participants (Yu 2022 Table 1). ",
        "Exactly one TRT_* indicator is 1 on a drug arm; all six are 0 on a ",
        "placebo arm."
      ),
      source_name = "Corticosteroids (Yu 2022 Table 2)"
    ),
    TRT_IMMUNOSUPPRESSANT = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received an immunosuppressant, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a different drug class).",
      notes = "MBMA study-arm-level indicator. Selects emax_immunosuppressant = -1.40 g/day. 11 trials / 13 arms / 415 participants (Yu 2022 Table 1).",
      source_name = "Immunosuppressant (Yu 2022 Table 2)"
    ),
    TRT_RAS_BLOCKER = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received a renin-angiotensin system blocker (ACE inhibitor or angiotensin receptor blocker), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a different drug class).",
      notes = "MBMA study-arm-level indicator. Selects emax_ras_blocker = -0.95 g/day. 18 trials / 27 arms / 634 participants (Yu 2022 Table 1).",
      source_name = "RAS blockers (Yu 2022 Table 2)"
    ),
    TRT_ANTIPLATELET = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received an antiplatelet agent, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a different drug class).",
      notes = "MBMA study-arm-level indicator. Selects emax_antiplatelet = -0.65 g/day. Only 2 trials / 2 arms / 28 participants (Yu 2022 Table 1), all with low baseline proteinuria (0.73-0.92 g/day), so this class's Emax is extrapolated well outside its own data at higher baselines.",
      source_name = "Antiplatelet agents (Yu 2022 Table 2)"
    ),
    TRT_OMEGA3_FA = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received N-3 (omega-3) polyunsaturated fatty acids (fish oil), 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a different drug class).",
      notes = "MBMA study-arm-level indicator. Selects emax_omega3_fa = -0.53 g/day. 4 trials / 5 arms / 157 participants (Yu 2022 Table 1).",
      source_name = "N-3 fatty acids (Yu 2022 Table 2)"
    ),
    TRT_OTHER_IGAN = list(
      description = "Binary study-arm drug-class indicator: 1 = the arm received an IgA-nephropathy treatment outside the five named classes, or a combination of drugs from different classes, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm or a named drug class).",
      notes = "MBMA study-arm-level indicator. Selects emax_other_igan = -1.31 g/day. 10 trials / 13 arms / 347 participants (Yu 2022 Table 1). A residual category defined by exclusion (Yu 2022 Results 'Characteristics of the Included Studies'), so it is heterogeneous and does not describe any single mechanism.",
      source_name = "Other drugs (Yu 2022 Table 2)"
    )
  )

  # Covariates Yu 2022 screened on the model parameters (Results 'Model
  # Building and Evaluation': "the influences of the daily urinary protein
  # excretion at baseline, age, and sex ratio on model parameters were
  # investigated") but did not retain. No coefficients are reported.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Arm-level median age.",
      units = "years",
      type = "continuous",
      notes = "Screened on the model parameters and not retained (Yu 2022 Results and Discussion). Arm-level ages 24.7-52 years, overall median 37 years (Yu 2022 Table 1). No coefficient reported."
    ),
    SEXF_PCT = list(
      description = "Arm-level percentage of participants who are female.",
      units = "%",
      type = "continuous",
      notes = "Yu 2022 screened the arm-level MALE percentage and did not retain it. Recorded here on the female-referenced scale, SEXF_PCT = 100 - male%. Source male% ranged 11.11-94.12 (median 58.33), i.e. female 5.88-88.89% (median 41.67%) (Yu 2022 Table 1). No coefficient reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2288L,
    n_studies = 40L,
    n_arms = 83L,
    age_range = "arm-level median age 24.7-52 years; overall median 37 years",
    sex_female_pct = NA_real_,
    sex_male_pct = "arm-level male percentage 11.11-94.12, overall median 58.33 (equivalently a female percentage of 5.88-88.89, median 41.67)",
    race_ethnicity = "Not reported at the arm level and not screened as a covariate.",
    disease_state = "Adults with IgA nephropathy (primary glomerular disease with mesangial IgA deposition) and persistent proteinuria.",
    baseline_proteinuria = "arm-level baseline daily urinary protein excretion 0.57-5.29 g/day; overall median 1.9 g/day",
    arms_by_class = paste0(
      "placebo 14 trials/14 arms/419 participants; corticosteroids 8/9/288; ",
      "immunosuppressants 11/13/415; RAS blockers 18/27/634; antiplatelet ",
      "agents 2/2/28; N-3 fatty acids 4/5/157; other drugs 10/13/347 ",
      "(Yu 2022 Table 1)"
    ),
    treatment_duration = "1-48 months across trials",
    dose_range = "not modelled; drug classes are pooled across agents and doses",
    regions = "International. PubMed and Embase searched up to 2019-11-18; English-language clinical trials published 1987-2017.",
    notes = paste0(
      "MBMA at the study-arm level: each data point is the arm-mean change ",
      "from baseline in daily urinary protein excretion (g/day) at one ",
      "follow-up time, digitised with Engauge Digitizer when only figures ",
      "were available. Drug classes follow the Japanese Society of ",
      "Nephrology 2014 guideline classification. The etas are BETWEEN-ARM, ",
      "not between-subject. The residual is weighted by the inverse square ",
      "root of the arm sample size (Yu 2022 Equation 3), so an arm of N ",
      "participants has residual SD addSd / sqrt(N); addSd in ini() is the ",
      "unit-weight value and the N weighting is applied downstream (same ",
      "convention as Guo_2025_glp1ReceptorAgonists_mbma and the other MBMA ",
      "models in this package). Covariates with more than 30 percent missing ",
      "values across studies (e.g. body weight and race) were not ",
      "investigated. Bootstrap success rate 98.1 percent over 1000 runs."
    )
  )

  ini({
    # ==================================================================
    # Placebo arm (Yu 2022 Equation 1 with the placebo parameters).
    # SIGN: negative = reduction in proteinuria, as printed in Table 2.
    # Emax is kept linear (signed) because it is negative.
    # ==================================================================
    emax_placebo <- -0.44
    label("Placebo maximum change from baseline in daily urinary protein excretion (g/day; negative = reduction)") # Yu 2022 Table 2 'Emax, placebo, g/day' = -0.44 (RSE 38.10%)

    let50_placebo <- log(27.2)
    label("Log time to half of the placebo maximum effect (log month)") # Yu 2022 Table 2 'ET50, placebo, month' = 27.20 (RSE 32.50%)

    # ==================================================================
    # Drug classes: class-specific Emax at the 1.82 g/day reference
    # baseline, one shared ET50 (Yu 2022 Results: "6 classes of drugs
    # with different pharmacological mechanisms shared the same ET50").
    # ==================================================================
    emax_corticosteroid <- -1.47
    label("Corticosteroid maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, Corticosteroids, g/day' = -1.47 (RSE 3.80%)

    emax_immunosuppressant <- -1.40
    label("Immunosuppressant maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, Immunosuppressant, g/day' = -1.40 (RSE 7.80%)

    emax_ras_blocker <- -0.95
    label("RAS blocker maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, RAS blockers, g/day' = -0.95 (RSE 12.00%)

    emax_antiplatelet <- -0.65
    label("Antiplatelet agent maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, Antiplatelet agents g/day' = -0.65 (RSE 18.20%)

    emax_omega3_fa <- -0.53
    label("N-3 fatty acid maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, N-3 fatty acids, g/day' = -0.53 (RSE 35.30%)

    emax_other_igan <- -1.31
    label("Other-drugs maximum change from baseline in urinary protein excretion at baseline 1.82 g/day (g/day)") # Yu 2022 Table 2 'Emax, Other drugs, g/day' = -1.31 (RSE 12.90%)

    let50_drug <- log(5.59)
    label("Log time to half of the drug maximum effect, shared by all six classes (log month)") # Yu 2022 Table 2 'ET50, Drug, typical, month' = 5.59 (RSE 28.40%)

    # ==================================================================
    # Baseline proteinuria on the drug-arm Emax (Yu 2022 Equation 7):
    #   Emax_drug_i = Emax_drug_typical - (Baseline - 1.82) * 0.63
    # written here as + e_upro_bl_emax * (UPRO_BL - 1.82) with the signed
    # Table 2 coefficient -0.63.
    # ==================================================================
    e_upro_bl_emax <- -0.63
    label("Linear effect of baseline urinary protein excretion on the drug-arm Emax (g/day per g/day)") # Yu 2022 Table 2 'Baseline on Emax, g/day' = -0.63 (RSE 18.30%); Equation 7

    # ==================================================================
    # Between-study-arm variability on ET50, scale form (Yu 2022
    # Equation 2: P_i = P_typical * (1 + eta_i), eta ~ N(0, omega^2)).
    # Table 2 prints 'eta (ET50)' = 0.65 with no stated scale. It is read
    # as the SD omega, so the ini() variance is 0.65^2, because its own
    # RSE rules out a variance: for 83 arms a variance estimate cannot
    # have an RSE below sqrt(2/83) = 15.5%, and the printed RSE is 7.4%.
    # The bootstrap 95% CI (0.55-0.90) implies an RSE of about 13%, also
    # below that floor. Read as a variance the value would be 0.65.
    #
    # 'eta (Emax)' = 0 FIXED in Table 2 (estimated near 0 and fixed "for
    # the stability of the model"). A zero-variance eta is omitted rather
    # than written as fixed(0), because a zero diagonal makes the omega
    # matrix singular for simulation; the model is otherwise identical.
    # ==================================================================
    eta_study_et50 ~ 0.4225 # Yu 2022 Table 2 'eta (ET50)' = 0.65 (RSE 7.40%; bootstrap median 0.68, 95% CI 0.55 to 0.90), read as the SD; 0.65^2 = 0.4225

    # ==================================================================
    # Residual error (Yu 2022 Equation 3):
    #   Y_obs_ij = Y_pred_ij + eps_ij / sqrt(N_ij), eps ~ N(0, sigma^2).
    # Table 2 'eps' = 1.54 is read as the SD sigma, on the same scale as
    # the eta row in the same table. Read as a variance, the unit-weight
    # SD would be sqrt(1.54) = 1.24 g/day.
    #
    # nlmixr2's add() takes a constant SD, so the 1/sqrt(N) arm-size
    # weighting is applied downstream: for an arm of N participants the
    # residual SD is addSd / sqrt(N).
    # ==================================================================
    addSd <- 1.54
    label("Additive residual SD on arm-mean change in urinary protein excretion (g/day) at UNIT study weight; the SD for an arm of N participants is addSd / sqrt(N)") # Yu 2022 Table 2 'eps' = 1.54 (RSE 8.00%; bootstrap median 1.48, 95% CI 1.16 to 1.75)
  })

  model({
    # Centring value of baseline urinary protein excretion (g/day), Yu 2022
    # Equation 7.
    ref_upro_bl <- 1.82

    # 1 on a drug arm, 0 on a placebo arm. The six TRT_* indicators are
    # mutually exclusive, so the sum is 0 or 1 in standard use.
    drug_arm <- TRT_CORTICOSTEROID + TRT_IMMUNOSUPPRESSANT + TRT_RAS_BLOCKER +
      TRT_ANTIPLATELET + TRT_OMEGA3_FA + TRT_OTHER_IGAN

    # Drug-arm Emax at the arm's baseline proteinuria (Yu 2022 Equation 7).
    emax_class <- emax_corticosteroid * TRT_CORTICOSTEROID +
      emax_immunosuppressant * TRT_IMMUNOSUPPRESSANT +
      emax_ras_blocker * TRT_RAS_BLOCKER +
      emax_antiplatelet * TRT_ANTIPLATELET +
      emax_omega3_fa * TRT_OMEGA3_FA +
      emax_other_igan * TRT_OTHER_IGAN
    emax_drug <- emax_class + e_upro_bl_emax * (UPRO_BL - ref_upro_bl)

    # Arm Emax and typical ET50: placebo parameters when drug_arm = 0, the
    # drug-class parameters when drug_arm = 1. Baseline proteinuria has no
    # effect on the placebo arm.
    emax_i <- drug_arm * emax_drug + (1 - drug_arm) * emax_placebo
    et50_typ <- drug_arm * exp(let50_drug) + (1 - drug_arm) * exp(let50_placebo)

    # Between-arm variability on ET50 in the scale form of Yu 2022
    # Equation 2. With omega = 0.65, about 6% of simulated arms draw
    # eta < -1 and hence a non-positive ET50, for which the Emax-in-time
    # curve is not defined; the vignette discusses this.
    et50_i <- et50_typ * (1 + eta_study_et50)

    # Arm-mean change from baseline in daily urinary protein excretion
    # (g/day), Yu 2022 Equation 1. Zero at t = 0, approaching emax_i.
    uprocfb <- emax_i * t / (et50_i + t)

    # Residual SD is the UNIT-WEIGHT value; divide by sqrt(arm N) downstream.
    uprocfb ~ add(addSd)
  })
}
