Liu_2026_sunitinib_multistate <- function() {
  description <- paste(
    "Six-state multistate model of tumour response, progression, treatment discontinuation and",
    "death in adults with cytokine-refractory metastatic renal cell carcinoma on sunitinib 50 mg",
    "once daily 4-weeks-on/2-weeks-off (Liu 2026; 165 patients, studies 1006 and 014). Six ODE",
    "states carry the occupancy probability of stable disease (s_stable, the initial state),",
    "response PR/CR (s_response), progression (s_progression), discontinuation due to adverse",
    "events (s_disc_ae), discontinuation for other reasons (s_disc_other) and death (s_death,",
    "absorbing); transitions are forward only. Constant hazards except stable -> response",
    "(Weibull, decreasing with time since the first dose) and progression -> death",
    "(log-logistic in the subject's time to progression T_PROGRESSION, so it is constant once the",
    "subject has progressed). Stable -> death and response -> death are fixed to a",
    "Gompertz-Makeham background mortality at baseline age; stable/response -> discontinuation",
    "for other reasons share one hazard, as do the two discontinuation -> death hazards.",
    "Stable -> progression falls with higher albumin, haemoglobin and age, rises with ECOG >= 1",
    "and is lower when the cancer was diagnosed more than 1 year earlier; progression -> death",
    "falls with higher albumin and with time since diagnosis > 1 year. No exposure term: drug",
    "exposure was not retained in the final model. No between-subject variability and no",
    "observation-error model in the source."
  )
  reference <- paste(
    "Liu H, Friberg LE. Time-Dependent Bias and Prognostic Confounding in Exposure-Response Analysis of",
    "Targeted Therapies: Lessons from Sunitinib in Metastatic Renal Cancer.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70447. PMID: 42587360; PMCID: PMC13468949.",
    "Final multistate model: Supplementary Material S3 (state equations, hazard functions and",
    "Table S4 with its footnotes). Gompertz-Makeham background mortality: Missov TI, Lenart A.",
    "Theor Popul Biol. 2013;90:29-35."
  )
  vignette <- "Liu_2026_sunitinib"

  # Multistate state-occupancy probabilities, same family as the states of
  # Liu_2023_bintrafuspAlfa_multistate.R and Lin_2024_TB_multistate.R.
  paper_specific_compartments <- c(
    "s_stable",
    "s_response",
    "s_progression",
    "s_disc_ae",
    "s_disc_other",
    "s_death"
  )

  units <- list(
    time = "week",
    dosing = "n/a (no drug-dosing events; the model propagates state-occupancy probabilities)",
    concentration = "probability (all six states are occupancy probabilities, not drug concentrations)"
  )

  compartmentData <- list(
    s_stable = list(
      analyte = "probability of being in the stable-disease state (S1, RECIST 1.0 SD)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_response = list(
      analyte = "probability of being in the response state (S2, RECIST 1.0 PR or CR)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_progression = list(
      analyte = "probability of being in the progression state (S3, RECIST 1.0 PD)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_disc_ae = list(
      analyte = "probability of having discontinued treatment because of adverse events (S4; entered 28 days after the last dose)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_disc_other = list(
      analyte = "probability of having discontinued treatment for other reasons (S5; sponsor decision, protocol violation, loss to follow-up or non-RECIST-confirmed progression)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_death = list(
      analyte = "probability of having died (S6)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log-hazard effect exp((ALB - 41) * theta) on stable -> progression (theta = -0.165) and",
        "progression -> death (theta = -0.0731), centred on the cohort median 41 g/L (Table S4 and",
        "footnote 2)."
      ),
      source_name = "ALB"
    ),
    HGB = list(
      description = "Haemoglobin at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Log-hazard effect exp((HB - 130) * theta) on stable -> progression, theta = -0.0255 per g/L (Table S4 and footnote 2).",
      source_name = "HB"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log-hazard effect exp((AGE - 57) * theta) on stable -> progression, theta = -0.0892 per",
        "year (Table S4 and footnote 2). Also drives the fixed Gompertz-Makeham background",
        "hazard of stable -> death and response -> death, evaluated at baseline age."
      ),
      source_name = "AGE"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG = 0)",
      notes = "Hazard ratio 1.48 on stable -> progression for ECOG >= 1 vs 0 (Table S4 'HR: ECOG >=1 vs 0'). Supplementary Table S1: ECOG 0 in 55%, 1 in 44%, > 1 in 1%.",
      source_name = "ECOG"
    ),
    T_DIAG_CANCER = list(
      description = "Time since the primary renal cell carcinoma diagnosis at the first sunitinib dose",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only as the dichotomy time since diagnosis > 1 year, derived in model() as",
        "T_DIAG_CANCER > 365.25: hazard ratio 0.950 on stable -> progression and 0.678 on",
        "progression -> death (Table S4 'HR: TDIAG >1y vs <=1y'; reference TDIAG <= 1 year,",
        "footnote 2). Supplementary Table S1 reports time since diagnosis in months, mean 38.0",
        "(SD 53.6), range 2.2-339."
      ),
      source_name = "TDIAG"
    ),
    T_PROGRESSION = list(
      description = "Subject's time of disease progression after the first dose",
      units = "week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time scale of the log-logistic progression -> death hazard (Table S4 'Time to Prog. on",
        "lambda36'). For a subject observed to progress, supply the observed time to progression",
        "on every record; for a subject who has not progressed its value only affects",
        "s_progression and s_death, which are then not meaningful. Population survival from",
        "baseline integrates over the time-to-progression distribution; the vignette shows how."
      ),
      source_name = "TTP"
    )
  )

  covariatesDataExcluded <- list(
    AUC_SUNITINIB = list(
      description = "Steady-state AUC of sunitinib plus SU12662 for 50 mg/day with cycle 1 apparent clearance (AUCss,DOSE50,CLcycle1)",
      units = "ug*h/mL",
      type = "continuous",
      notes = paste(
        "Tested as exp(beta * (AUC - 2.19)) on every transition (Table 2). Univariable HR per",
        "IQR (1.37 ug*h/mL): 0.71 on stable -> progression, 1.48 on stable -> disc. AEs, 0.78 on",
        "progression -> death; after adjustment for age and albumin the HRs were 1.13 (P = 0.27)",
        "and 0.94 (P = 0.41). Not part of the final model in Table S4."
      )
    ),
    SLD = list(
      description = "Sum of longest diameters of target lesions at baseline",
      units = "mm",
      type = "continuous",
      notes = "Listed in the Table S4 abbreviations and screened; not retained in the final multistate model."
    ),
    STUDY_014 = list(
      description = "Study 014 indicator (vs study 1006)",
      units = "(binary)",
      type = "binary",
      notes = "A study covariate was tested on all transitions (Supplementary S3); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 165L,
    n_studies = 2L,
    age_range = "24-87 years",
    age_mean = "57.1 years (SD 10.3)",
    weight_range = "45-150 kg",
    weight_mean = "85.0 kg (SD 18.6)",
    sex_female_pct = 35,
    race_ethnicity = c(White = 92, NonWhite = 8),
    disease_state = "Cytokine-refractory metastatic renal cell carcinoma",
    dose_range = "Sunitinib 50 mg orally once daily, 4 weeks on / 2 weeks off; 38.2% had a dose reduction to 37.5 or 25 mg",
    regions = "Multinational (Pfizer phase II trials NCT00077974 and NCT00054886)",
    albumin = "40.4 g/L (SD 4.8), range 20-49",
    haemoglobin = "124.3 g/L (SD 18.8), range 86-169",
    ecog = "0: 55%; 1: 44%; > 1: 1%",
    time_since_diagnosis = "38.0 months (SD 53.6), range 2.2-339",
    notes = paste(
      "Exposure-response analysis set of 165 of the 169 patients dosed in studies 1006 (n = 106",
      "enrolled) and 014 (n = 63); 4 excluded for incomplete dose-interruption records. Baseline",
      "characteristics from Supplementary Table S1. Figure 1 transition counts: 85 stable ->",
      "response, 56 stable -> progression, 62 response -> progression, 15 stable -> disc. AEs,",
      "3 response -> disc. AEs, 6 stable -> disc. other, 2 response -> disc. other, 2 stable ->",
      "death, 4 response -> death, 87 progression -> death, 15 disc. AEs -> death, 2 disc. other",
      "-> death; 110 deaths."
    )
  )

  ini({
    # Hazards on a WEEK time axis (Table S4 'Baseline transitions (week^-1)').
    # Baseline values correspond to HB 130 g/L, ALB 41 g/L, AGE 57 years,
    # ECOG 0 and time since diagnosis <= 1 year (footnote 2).
    lscale12 <- log(0.0471)
    label("Log Weibull scale of the hazard stable -> response (1/week)") # Table S4 lambda12 scale 0.0471 (RSE 13.4%)
    lshape12 <- log(0.844)
    label("Log Weibull shape of the hazard stable -> response, time since the first dose (unitless)") # Table S4 'Time since treatment initiation on lambda12' 0.844 (RSE 11.5%)
    lk13 <- log(0.00725)
    label("Log hazard stable -> progression lambda13 (1/week)") # Table S4 lambda13 0.00725 (RSE 30.5%)
    lk23 <- log(0.0171)
    label("Log hazard response -> progression lambda23 (1/week)") # Table S4 lambda23 0.0171 (RSE 13.5%)
    lk14 <- log(0.00758)
    label("Log hazard stable -> discontinuation due to AEs lambda14 (1/week)") # Table S4 lambda14 0.00758 (RSE 25.7%)
    lk24 <- log(0.000758)
    label("Log hazard response -> discontinuation due to AEs lambda24 (1/week)") # Table S4 lambda24 0.000758 (RSE 57.7%)
    lk15 <- log(0.00131)
    label("Log hazard stable or response -> discontinuation for other reasons, lambda15 = lambda25 (1/week)") # Table S4 lambda15 0.00131 (RSE 37.2%); lambda25 merged (footnote c)
    lscale36 <- log(0.0340)
    label("Log log-logistic scale of the hazard progression -> death (1/week)") # Table S4 lambda36 0.0340 (RSE 39.5%)
    lshape36 <- log(0.930)
    label("Log log-logistic shape of the hazard progression -> death, time to progression (unitless)") # Table S4 'Time to Prog. on lambda36' 0.930 (RSE 2.8%)
    lk46 <- log(0.0267)
    label("Log hazard discontinuation (AEs or other) -> death, lambda46 = lambda56 (1/week)") # Table S4 lambda46 0.0267 (RSE 22.4%); lambda56 merged (footnote e)

    # Covariate effects (Table S4 'Covariate effects').
    e_alb_k13 <- -0.165
    label("Log-hazard coefficient of (ALB - 41) on stable -> progression (1/(g/L))") # Table S4 ALB on lambda13 -0.165 (RSE 12.4%)
    e_hgb_k13 <- -0.0255
    label("Log-hazard coefficient of (HB - 130) on stable -> progression (1/(g/L))") # Table S4 HB on lambda13 -0.0255 (RSE 32.4%)
    e_age_k13 <- -0.0892
    label("Log-hazard coefficient of (AGE - 57) on stable -> progression (1/year)") # Table S4 AGE on lambda13 -0.0892 (RSE 13.9%)
    e_ecog_ge1_k13 <- log(1.48)
    label("Log hazard ratio of ECOG >= 1 vs 0 on stable -> progression (unitless)") # Table S4 ECOG on lambda13 HR 1.48 (RSE 29.3%)
    e_tdiag_gt1y_k13 <- log(0.950)
    label("Log hazard ratio of time since diagnosis > 1 year on stable -> progression (unitless)") # Table S4 TDIAG on lambda13 HR 0.950 (RSE 31.8%)
    e_alb_k36 <- -0.0731
    label("Log-hazard coefficient of (ALB - 41) on progression -> death (1/(g/L))") # Table S4 ALB on lambda36 -0.0731 (RSE 39.5%)
    e_tdiag_gt1y_k36 <- log(0.678)
    label("Log hazard ratio of time since diagnosis > 1 year on progression -> death (unitless)") # Table S4 TDIAG on lambda36 HR 0.678 (RSE 26.7%)

    # Stable -> death and response -> death fixed to the age-specific natural
    # death rate hGM(age) = 0.0007 + 0.00018645 * exp(0.11425 * (age - 30)) per
    # year (Table S4 footnote d; Missov and Lenart 2013).
    gm_c <- fixed(0.0007)
    label("Gompertz-Makeham age-independent (Makeham) hazard c (1/year)") # Table S4 footnote d: 0.0007
    gm_alpha <- fixed(0.00018645)
    label("Gompertz-Makeham senescent-mortality level alpha at age 30 (1/year)") # Table S4 footnote d: 0.00018645
    gm_beta <- fixed(0.11425)
    label("Gompertz-Makeham rate of ageing beta (1/year of age)") # Table S4 footnote d: 0.11425

    # The source maximises the multistate likelihood on the observed state and
    # has no residual error. A tiny additive residual on prob_death lets the
    # nlmixr2 likelihood machinery accept the model for forward simulation;
    # not from the source.
    addSd <- 0.001
    label("Placeholder additive residual error on the death probability output (unitless); not from the source")
  })

  model({
    # Covariate terms (Table S4 footnote 2 reference values).
    tdiag_gt1y <- T_DIAG_CANCER > 365.25
    cov13 <- exp(e_alb_k13 * (ALB - 41) + e_hgb_k13 * (HGB - 130) + e_age_k13 * (AGE - 57) +
      e_ecog_ge1_k13 * ECOG_GE1 + e_tdiag_gt1y_k13 * tdiag_gt1y)
    cov36 <- exp(e_alb_k36 * (ALB - 41) + e_tdiag_gt1y_k36 * tdiag_gt1y)

    # Weibull hazard in time since the first dose (Supplementary S3):
    # lambda(t) = scale * shape * (scale * t)^(shape - 1). With shape < 1 the
    # hazard is infinite at t = 0 although its integral is finite; the 1e-16
    # offset keeps the ODE right-hand side finite there.
    scale12 <- exp(lscale12)
    shape12 <- exp(lshape12)
    k12 <- scale12 * shape12 * (scale12 * (t + 1e-16))^(shape12 - 1)

    # Log-logistic hazard in the time to progression, in the same rate
    # parameterisation as the Weibull:
    # lambda36 = scale * shape * (scale * TTP)^(shape - 1) / (1 + (scale * TTP)^shape).
    scale36 <- exp(lscale36)
    shape36 <- exp(lshape36)
    ttp36 <- T_PROGRESSION + 1e-16
    k36 <- scale36 * shape36 * (scale36 * ttp36)^(shape36 - 1) /
      (1 + (scale36 * ttp36)^shape36) * cov36

    k13 <- exp(lk13) * cov13
    k23 <- exp(lk23)
    k14 <- exp(lk14)
    k24 <- exp(lk24)
    k15 <- exp(lk15)
    k25 <- k15
    k46 <- exp(lk46)
    k56 <- k46

    # Background mortality at baseline age, per year -> per week.
    gmHazYear <- gm_c + gm_alpha * exp(gm_beta * (AGE - 30))
    k16 <- gmHazYear / 52.1429
    k26 <- k16

    # Supplementary S3 state equations; all patients start in stable disease.
    d/dt(s_stable) <- -s_stable * (k12 + k13 + k14 + k15 + k16)
    s_stable(0) <- 1
    d/dt(s_response) <- s_stable * k12 - s_response * (k23 + k24 + k25 + k26)
    s_response(0) <- 0
    d/dt(s_progression) <- s_stable * k13 + s_response * k23 - s_progression * k36
    s_progression(0) <- 0
    d/dt(s_disc_ae) <- s_stable * k14 + s_response * k24 - s_disc_ae * k46
    s_disc_ae(0) <- 0
    d/dt(s_disc_other) <- s_stable * k15 + s_response * k25 - s_disc_other * k56
    s_disc_other(0) <- 0
    d/dt(s_death) <- s_stable * k16 + s_response * k26 + s_progression * k36 +
      s_disc_ae * k46 + s_disc_other * k56
    s_death(0) <- 0

    # State-occupancy probabilities. For fitting, the source likelihood resets
    # the states to the observed state at each transition; supply that
    # conditioning as evid = 5 replace events.
    prob_stable_disease <- s_stable
    prob_tumor_response <- s_response
    prob_progressive_disease <- s_progression
    prob_discontinuation_ae <- s_disc_ae
    prob_discontinuation_other <- s_disc_other
    prob_death <- s_death
    # Mass-balance diagnostic: must stay at 1.
    prob_state_sum <- s_stable + s_response + s_progression + s_disc_ae + s_disc_other + s_death

    prob_death ~ add(addSd)
  })
}
