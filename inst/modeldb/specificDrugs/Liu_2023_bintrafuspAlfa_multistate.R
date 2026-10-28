Liu_2023_bintrafuspAlfa_multistate <- function() {
  description <- paste(
    "Six-state pharmacometric multistate model of tumour response, progression and survival in adults",
    "with second-line non-small cell lung cancer treated with bintrafusp alfa 500 or 1200 mg every 2 weeks",
    "(phase I trial NCT02517398, 2L NSCLC expansion cohort); final model with dose-independent transitions.",
    "Six ODE states carry the marginal state-occupancy probability of stable disease (s_stable, the",
    "initial state), response PR/CR (s_response), progression (s_progression), unknown tumour status after",
    "treatment stop without progression (s_unknown), dropout from survival follow-up (s_dropout, absorbing)",
    "and death (s_death, absorbing). Transitions are forward only. Constant hazards for stable -> response,",
    "stable -> unknown, response -> progression / unknown / death, unknown -> death and a shared dropout",
    "hazard from stable, progression and unknown (no response -> dropout); Weibull hazards in time since",
    "treatment initiation for stable -> progression and progression -> death (both decreasing, shape < 1);",
    "stable -> death fixed to an age-dependent gamma-Gompertz-Makeham background mortality. Absence of",
    "metastases at baseline lowers the stable -> progression hazard 4.18-fold. No between-subject",
    "variability and no observation-error model in the source."
  )
  reference <- paste(
    "Liu H, Milenkovic-Grisic AM, Krishnan SM, Jonsson S, Friberg LE, Girard P, Venkatakrishnan K,",
    "Vugmeyster Y, Khandelwal A, Karlsson MO.",
    "A multistate modeling and simulation framework to learn dose-response of oncology drugs:",
    "Application to bintrafusp alfa in non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1738-1750. doi:10.1002/psp4.12976.",
    "PMCID: PMC10681430. Equations 1-10 and Table 2 ('Final multistate model' column) of the main text;",
    "covariate coding from Supplementary Table S3; Gompertz-Makeham constants and the hazard code from the",
    "example NONMEM control stream in Appendix S2. Gompertz-Makeham form: Missov TI, Lenart A.",
    "Theor Popul Biol. 2013;90:29-35."
  )
  vignette <- "Liu_2023_bintrafuspAlfa"

  # Multistate state-occupancy probabilities, same family as the states of
  # Lin_2024_TB_multistate.R and Ibrahim_2023_ibrutinib_competing_risk.R.
  paper_specific_compartments <- c(
    "s_stable",
    "s_response",
    "s_progression",
    "s_unknown",
    "s_dropout",
    "s_death"
  )

  units <- list(
    time = "week",
    dosing = "n/a (no drug-dosing events; the model propagates state-occupancy probabilities)",
    concentration = "probability (all six states are occupancy probabilities, not drug concentrations)"
  )

  compartmentData <- list(
    s_stable = list(
      analyte = "probability of being in the stable-disease state (RECIST SD)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_response = list(
      analyte = "probability of being in the response state (RECIST PR or CR)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_progression = list(
      analyte = "probability of being in the progression state (RECIST PD)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_unknown = list(
      analyte = "probability of being in the unknown-tumour-status state (treatment and tumour evaluation stopped without documented progression)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_dropout = list(
      analyte = "probability of having dropped out of survival follow-up",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_death = list(
      analyte = "probability of having died",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    MET_GE1 = list(
      description = "Metastases at baseline: 1 = at least one metastatic site, 0 = no metastases.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (metastases present; 65 of 80 patients)",
      notes = paste(
        "Time-fixed per subject. The source column METAS is coded the OTHER way round: Supplementary",
        "Table S3 'Metastases (0:Yes or 1:No)'. The Appendix S2 control stream writes",
        "K13 = WB13 * EXP(METAS * THETA(12)) with THETA(12) = beta = -1.431 (Table 2), so the source",
        "reference (METAS = 0) is the metastatic majority, and patients WITHOUT metastases have a",
        "k13 lowered by exp(-1.431) = 0.239. Here METAS = 1 - MET_GE1, giving",
        "k13 = Weibull(t) * exp(e_metge1_k13 * (1 - MET_GE1)). Equivalently, metastases raise k13 by",
        "exp(1.431) = 4.18-fold relative to no metastases, as stated in the Results."
      ),
      source_name = "METAS (= 1 - MET_GE1)"
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only the fixed stable -> death background hazard k16 through the gamma-Gompertz-Makeham",
        "mortality curve (Appendix S2 control stream). Baseline age, as in the control stream $PK block;",
        "a time-varying age column is also accepted."
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_STEROID = list(
      description = "Concomitant corticosteroids (1 = yes).",
      units = "(binary)",
      type = "binary",
      notes = "Significant on k13 in the univariable analysis (dOFV = -4.93, coefficient -0.67) but driven by 2 influential individuals; not retained (Supplementary Table S3)."
    ),
    RACE_ASIAN = list(
      description = "Asian race (1 = Asian).",
      units = "(binary)",
      type = "binary",
      notes = "Significant on k23 and on the dropout hazard (all 13 dropouts were Asian patients from 4 Korean/Taiwanese sites, confounded by study site); not retained (Results; Supplementary Tables S3-S4)."
    ),
    CRP = list(
      description = "C-reactive protein at baseline.",
      units = "mg/L",
      type = "continuous",
      notes = "Significant on the dropout hazard (dOFV = -3.84), 1 influential individual; not retained (Supplementary Table S3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "38-85 years",
    age_mean = "63.9 years (500 mg arm), 62.1 years (1200 mg arm)",
    weight_range = "41.9-109.0 kg",
    sex_female_pct = 28.75,
    race_ethnicity = c(White = 48.75, Black = 1.25, Asian = 50),
    disease_state = "Advanced non-small cell lung cancer, second line (2L NSCLC); ECOG 0-1; 81% with metastases at baseline",
    dose_range = "500 or 1200 mg IV every 2 weeks (1:1 randomized, 40 patients per arm)",
    regions = "USA, United Kingdom, Belgium, Italy, Spain, South Korea, Taiwan (25 sites)",
    notes = paste(
      "Baseline characteristics from Liu 2023 Table 1. 406 multistate observations: 326 post-baseline",
      "RECIST 1.1 response categories per independent review (75 SD, 127 PR, 124 PD) and 80 absorbing",
      "events (56 deaths, 13 dropouts, 11 censored at data cutoff); 155 between-state transitions.",
      "Tumour assessments every 6 weeks in the first year and every 12 weeks thereafter. OS follow-up",
      "median 52 weeks (range 1-180)."
    )
  )

  ini({
    # All hazards are on a WEEK time axis. Table 2 reports the constant
    # hazards as rates in 1/week (footnote b) and the Weibull hazards as
    # k_ij(t) = scale * shape * (scale * t)^(shape - 1) with t the time since
    # treatment initiation in weeks (footnote c).
    lk12 <- log(0.032)
    label("Log hazard stable disease -> response k12 (1/week)") # Table 2 final model: k12 = 0.032 (95% CI 0.021, 0.050)
    lscale13 <- log(0.157)
    label("Log Weibull scale of the hazard stable disease -> progression (1/week)") # Table 2 final model: scale13 = 0.157 (0.086, 0.341)
    lshape13 <- log(0.427)
    label("Log Weibull shape of the hazard stable disease -> progression (unitless)") # Table 2 final model: shape13 = 0.427 (0.247, 0.631)
    lk14 <- log(0.012)
    label("Log hazard stable disease -> unknown k14 (1/week)") # Table 2 final model: k14 = 0.012 (0.006, 0.023)
    lk23 <- log(0.013)
    label("Log hazard response -> progression k23 (1/week)") # Table 2 final model: k23 = 0.013 (0.008, 0.022)
    lk24 <- log(0.003)
    label("Log hazard response -> unknown k24 (1/week)") # Table 2 final model: k24 = 0.003 (0.001, 0.008)
    lk26 <- log(0.002)
    label("Log hazard response -> death k26 (1/week)") # Table 2 final model: k26 = 0.002 (0.0008, 0.006)
    lscale36 <- log(0.023)
    label("Log Weibull scale of the hazard progression -> death (1/week)") # Table 2 final model: scale36 = 0.023 (0.015, 0.044)
    lshape36 <- log(0.631)
    label("Log Weibull shape of the hazard progression -> death (unitless)") # Table 2 final model: shape36 = 0.631 (0.439, 0.862)
    lk46 <- log(0.008)
    label("Log hazard unknown -> death k46 (1/week)") # Table 2 final model: k46 = 0.008 (0.004, 0.015)
    lkdrop <- log(0.003)
    label("Log dropout hazard k15 = k35 = k45 from stable, progression and unknown (1/week)") # Table 2 final model: k1,3,4-5 = 0.003 (0.002, 0.005)

    e_metge1_k13 <- -1.431
    label("Log-hazard coefficient on k13 for NO metastases at baseline, applied to (1 - MET_GE1)") # Table 2 final model: beta METAS on k13 = -1.431 (-2.570, -0.497); exp(1.431) = 4.18-fold (Results)

    # Background mortality for stable disease -> death (k16), fixed to the
    # gamma-Gompertz-Makeham distribution (Table 2 footnote d). Constants from
    # the Appendix S2 control stream $PK block; the curve gives an annual
    # hazard at age AGE.
    gm_alpha <- fixed(1.8645e-04)
    label("Gompertz-Makeham senescent-mortality level alpha at age 30 (1/year)") # Appendix S2: ALPH = 1.8645e-04
    gm_beta <- fixed(0.11425)
    label("Gompertz-Makeham rate of ageing beta (1/year of age)") # Appendix S2: BETA = 0.11425
    gm_c <- fixed(0.0007)
    label("Gompertz-Makeham age-independent (Makeham) hazard c (1/year)") # Appendix S2: CVALUE = 0.0007
    gm_gamma <- fixed(0.1)
    label("Gompertz-Makeham frailty variance gamma (unitless)") # Appendix S2: GAMA = 0.1

    # No IIV: the control stream declares BSV = ETA(1) as a place holder with
    # $OMEGA 0 FIX. The source maximises the multistate likelihood on the
    # observed categorical state (LAPLACE LIKE) and has no residual error. A
    # tiny additive residual on prob_death lets the nlmixr2 likelihood
    # machinery accept the model for forward simulation; not from the source.
    addSd <- 0.001
    label("Placeholder additive residual error on the death probability output (unitless); not from the source")
  })

  model({
    # Background mortality (Appendix S2 control stream):
    #   X1 = EXP(BETA*AGE-30); Y1 = ALPH*GAMA/BETA
    #   GGM = CVALUE + ALPH*X1/(1 + Y1*(X1-1)); MTT16 = 52.1429/GGM weeks
    # The printed EXP(BETA*AGE-30) is read as EXP(BETA*(AGE-30)), the
    # Missov-Lenart gamma-Gompertz-Makeham curve with age counted from 30;
    # see the vignette Assumptions section. 52.1429 weeks per year.
    gmX <- exp(gm_beta * (AGE - 30))
    gmY <- gm_alpha * gm_gamma / gm_beta
    gmHazYear <- gm_c + gm_alpha * gmX / (1 + gmY * (gmX - 1))
    k16 <- gmHazYear / 52.1429

    # Weibull hazards in time since treatment initiation (equation 9). The
    # 1e-16 offset is the control stream's DEL guard: with shape < 1 the
    # hazard is infinite at t = 0, although its integral is finite.
    scale13 <- exp(lscale13)
    shape13 <- exp(lshape13)
    scale36 <- exp(lscale36)
    shape36 <- exp(lshape36)
    k13 <- scale13 * shape13 * (scale13 * (t + 1e-16))^(shape13 - 1) *
      exp(e_metge1_k13 * (1 - MET_GE1))
    k36 <- scale36 * shape36 * (scale36 * (t + 1e-16))^(shape36 - 1)

    k12 <- exp(lk12)
    k14 <- exp(lk14)
    k23 <- exp(lk23)
    k24 <- exp(lk24)
    k26 <- exp(lk26)
    k46 <- exp(lk46)
    # Dropout completely at random from stable, progression and unknown
    # (k15 = k35 = k45); no response -> dropout transition was observed, so
    # k25 = 0 (control stream).
    k15 <- exp(lkdrop)
    k35 <- k15
    k45 <- k15

    # Equations 1-6; all patients start in stable disease at the first dose.
    d/dt(s_stable) <- -s_stable * (k12 + k13 + k14 + k15 + k16)
    s_stable(0) <- 1
    d/dt(s_response) <- s_stable * k12 - s_response * (k23 + k24 + k26)
    s_response(0) <- 0
    d/dt(s_progression) <- s_stable * k13 + s_response * k23 - s_progression * (k35 + k36)
    s_progression(0) <- 0
    d/dt(s_unknown) <- s_stable * k14 + s_response * k24 - s_unknown * (k45 + k46)
    s_unknown(0) <- 0
    d/dt(s_dropout) <- s_stable * k15 + s_progression * k35 + s_unknown * k45
    s_dropout(0) <- 0
    d/dt(s_death) <- s_stable * k16 + s_response * k26 + s_progression * k36 + s_unknown * k46
    s_death(0) <- 0

    # Marginal state-occupancy probabilities (Figure 2 VPC panels). The
    # source likelihood resets every state to 0 and the observed state to 1
    # after each observation; supply that conditioning as evid = 5 replace
    # events in the data when fitting.
    prob_stable_disease <- s_stable
    prob_tumor_response <- s_response
    prob_progressive_disease <- s_progression
    prob_unknown_status <- s_unknown
    prob_dropout <- s_dropout
    prob_death <- s_death
    # Mass-balance diagnostic, equation 7: must stay at 1.
    prob_state_sum <- s_stable + s_response + s_progression + s_unknown + s_dropout + s_death

    prob_death ~ add(addSd)
  })
}
