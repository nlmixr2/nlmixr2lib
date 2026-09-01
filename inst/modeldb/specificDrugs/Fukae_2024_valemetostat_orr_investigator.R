Fukae_2024_valemetostat_orr_investigator <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-response model for overall ",
    "response rate (ORR) by INVESTIGATOR assessment in adults with ",
    "relapsed/refractory adult T-cell leukemia/lymphoma (ATLL) treated ",
    "with the oral EZH1/EZH2 dual inhibitor valemetostat (Fukae 2024, ",
    "n = 38 pooled across the J101 phase I and J201 phase II trials, ",
    "150-200 mg once daily). The probability of response is ",
    "expit(mu + bE*zE + x'b1 + x'b2*zE) where zE is the ",
    "centred-and-scaled unbound valemetostat steady-state AUC ",
    "(AUCU_VALE - 375)/250 ng*h/mL, x is the baseline covariate vector ",
    "(age, LDH, weight centred-and-scaled; sex and ECOG performance ",
    "status as 0/1 indicators), b1 are covariate effects on the logit ",
    "intercept and b2 are covariate effects on the exposure slope ",
    "(Fukae 2024 Methods equation, Table 2). Candidate covariate ",
    "effects were regularized via normal-mixture spike-and-slab priors ",
    "in a Bayesian framework (spike Normal(0, 0.1), slab ",
    "Normal(0, 2.5), mixing weight Beta(1, 1)); the intercept and the ",
    "exposure slope carried weakly informative priors and were not ",
    "regularized. Fitted with CmdStanR 0.4.0 (Hamiltonian Monte Carlo, ",
    "no-U-turn sampler); all parameters reached effective sample size ",
    ">= 1000 with Gelman-Rubin statistic < 1.05. This model estimates a ",
    "steeper exposure slope than its central-assessment companion ",
    "(odds ratio 1.22 vs 1.08 per 250 ng*h/mL) on a larger analysis ",
    "set, with shallower covariate effects. There is no PK layer and ",
    "no ODE: individual unbound AUCss is supplied as data from the ",
    "companion population PK model (doi:10.1002/psp4.13201). No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Companion central-assessment model in ",
    "Fukae_2024_valemetostat_orr_central; six companion ",
    "exposure-safety models in the Fukae_2024_valemetostat_* family."
  )
  reference <- paste(
    "Fukae M, Rogers J, Garcia R, Tachibana M, Shimizu T.",
    "Bayesian sparse regression for exposure-response analyses of efficacy",
    "and safety endpoints to justify the clinical dose of valemetostat for",
    "adult T-cell leukemia/lymphoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1655-1669.",
    "doi:10.1002/psp4.13203.",
    "Individual unbound AUCss inputs derive from the companion population PK",
    "model: Fukae M, Baron K, Tachibana M, et al. Population pharmacokinetics",
    "of total and unbound valemetostat and platelet dynamics in healthy",
    "volunteers and patients with non-Hodgkin lymphoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2024. doi:10.1002/psp4.13201.",
    sep = " "
  )
  vignette <- "Fukae_2024_valemetostat_exposure_response"
  units <- list(
    time          = "n/a (static landmark exposure-response model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUCU_VALE covariate column)",
    concentration = "prob_orr_investigator (probability of overall response, 0-1; also logit_orr_investigator)"
  )

  covariateData <- list(
    AUCU_VALE = list(
      description        = "Unbound (free) valemetostat plasma AUC over the once-daily 24 h dosing interval at steady state, per subject. Supplied as data: this model has no PK layer, and the source analysis used empirical-Bayes individual predictions from the companion population PK model (Fukae 2024, doi:10.1002/psp4.13201).",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 375 ng*h/mL and scaled by 250 ng*h/mL inside model(). 375 is the reference patient's exposure (Fukae 2024 Table 2 footnote), approximately the typical unbound AUCss at the approved 200 mg once-daily dose; 250 is the increment every published odds ratio is expressed per (Table 2 row 'Unbound valemetostat AUCSS: 250 ng*h/mL increase'). Observed 5th-95th percentile of the pooled exposure distribution 184-887 ng*h/mL. Unbound, not total: the coefficients would be badly mis-scaled if a total AUCss were supplied.",
      source_name        = "unbound valemetostat AUCSS"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 65 years and scaled by 10 years inside model() (Fukae 2024 Table 2 footnote reference patient; Table 2 row 'Age: 10 years increase'). ATLL efficacy analysis set median 69.0 years, range 37-84 (Table 1).",
      source_name        = "Age"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase concentration.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 250 U/L and scaled by 300 U/L inside model() (Fukae 2024 Table 2 footnote reference patient; Table 2 row 'LDH: 300 U/L increase'). ATLL efficacy analysis set median 315 U/L, range 143-2000 (Table 1).",
      source_name        = "LDH"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 63 kg and scaled by 20 kg inside model() (Fukae 2024 Table 2 footnote reference patient; Table 2 row 'Weight: 20 kg increase'). ATLL efficacy analysis set median 61.5 kg, range 34.5-111 (Table 1).",
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Binary covariates are NOT centred or scaled in this model -- Fukae 2024 Methods states that centring and scaling were applied to 'continuous covariates only'.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the reference patient is male per the Table 2 footnote)",
      notes              = "Enters as a log-odds shift on the logit intercept (e_sexf_logit) and as a modifier of the exposure slope (e_sexf_slope). ATLL efficacy analysis set 48.7% female (Table 1).",
      source_name        = "Sex: female"
    ),
    ECOG_GE1 = list(
      description        = "Baseline Eastern Cooperative Oncology Group performance status indicator; 1 = ECOG PS >= 1, 0 = ECOG PS 0. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG PS 0; the reference patient has ECOG PS 0 per the Table 2 footnote)",
      notes              = "The paper collapses ECOG PS to '0' versus '1+' (Fukae 2024 Methods and Table 1). The investigator-assessed model regularizes this effect much more strongly than the central-assessment model does (odds ratio 0.903 vs 0.378), which the authors attribute to the larger analysis set. ATLL efficacy analysis set 46.2% with ECOG PS >= 1 (Table 1).",
      source_name        = "ECOG PS score: >= 1"
    )
  )

  covariatesDataExcluded <- list(
    DIS_NHL_STAGE4 = list(
      description = "Non-Hodgkin lymphoma disease stage indicator; 1 = stage IV, 0 = stage I/II/III.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Named by Fukae 2024 Methods as a candidate categorical covariate for the exposure-efficacy analysis ('NHL disease stage (I-III/IV)'), but NOT reported in Table 2 and therefore not retained in the final model -- no point estimate exists anywhere on disk. Stage was heavily missing in the source cohort (14 of 39 ATLL patients, 35.9%; Table 1), the likely reason it was dropped. Documented here to preserve the covariate screen without carrying a convention warning."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38L,
    n_studies      = 2L,
    n_observations = "38 binary response records (one per patient; landmark analysis, no repeated measures)",
    age_range      = "median 69.0 years, range 37-84 across the 39-patient ATLL set (Fukae 2024 Table 1)",
    weight_range   = "median 61.5 kg, range 34.5-111 (Fukae 2024 Table 1)",
    sex_female_pct = 48.7,
    race_ethnicity = c(Asian = 87.2, Black = 12.8),
    disease_state  = "relapsed or refractory adult T-cell leukemia/lymphoma (ATLL); responses per modified 2009 ATLL response criteria, ORR = complete response (including uncertified CR) or partial response",
    dose_range     = "valemetostat 150-200 mg orally once daily until progressive disease or unacceptable toxicity (J101 dose escalation contributed 150 and 200 mg; J201 was 200 mg throughout)",
    regions        = "Japan and the United States (J101, phase I first-in-human) and Japan (J201, DS3201-A-J201, NCT04102150)",
    notes          = paste0(
      "Investigator assessment was available in both trials, so this ",
      "model's analysis set (38) is larger than the central-assessment ",
      "companion's (25). One of the 39 ATLL patients had a missing ORR ",
      "and was excluded (Fukae 2024 Table 1 footnote a). Observed ",
      "outcome: 22 responders (ORR 57.9%) and 16 non-responders (Fukae ",
      "2024 Results, Patient populations)."
    )
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Fukae 2024 Table 2, column
    # "ORR by investigator assessment". The table reports every effect
    # as an ODDS RATIO -- the footnote to the companion Table 3 states
    # the convention explicitly: "Effects are described as exp(beta)".
    # Each ini() value below is therefore log(printed odds ratio),
    # written in log() form so the published number is visible at the
    # trace site.
    #
    # The reference patient (Table 2 footnote) is a 65-year-old Asian
    # male weighing 63 kg, with normal hepatic function, ECOG PS 0,
    # baseline LDH 250 U/L and unbound valemetostat AUCss 375 ng*h/mL.
    # Those values are the centring constants applied in model(); the
    # per-row increments (250 ng*h/mL, 10 years, 300 U/L, 20 kg) are the
    # scaling constants. Both sets are corroborated by Table 1: the
    # centring values are the rounded cohort medians and the scaling
    # values are the rounded cohort standard deviations, exactly as
    # Fukae 2024 Methods describes ("centered and scaled at approximate
    # population median and standard deviation values").
    # ==================================================================

    # ----- Logit intercept (reference-patient response probability) -----
    # Cross-check: the Results narrative quotes "the reference patient
    # had a predicted 60% (95% CrI, 44%-77%) probability of achieving
    # clinical response", matching 0.602 (0.443, 0.770).
    logit_ref        <- log(0.602 / (1 - 0.602)) ; label("Logit of the response probability for the reference patient (unitless logit)")  # Fukae 2024 Table 2, "Population mean" 0.602 (0.443, 0.770); ESS 3881, Rhat 1.00

    # ----- Exposure slope on the logit (beta*_E) -----
    e_aucu_logit     <- log(1.22)    ; label("Log-odds of response per 250 ng*h/mL increase in unbound valemetostat AUCss (unitless logit)")   # Fukae 2024 Table 2, OR 1.22 (0.632, 2.87); ESS 3398, Rhat 1.00

    # ----- Covariate effects on the logit intercept (beta*_1) -----
    e_age_logit      <- log(1.01)    ; label("Log-odds shift on the response logit per 10-year increase in age (unitless logit)")               # Fukae 2024 Table 2, OR 1.01 (0.795, 1.28); ESS 9131, Rhat 1.00
    e_ldh_logit      <- log(0.992)   ; label("Log-odds shift on the response logit per 300 U/L increase in baseline LDH (unitless logit)")      # Fukae 2024 Table 2, OR 0.992 (0.866, 1.12); ESS 11283, Rhat 1.00
    e_wt_logit       <- log(0.984)   ; label("Log-odds shift on the response logit per 20 kg increase in body weight (unitless logit)")         # Fukae 2024 Table 2, OR 0.984 (0.765, 1.24); ESS 10884, Rhat 1.00
    e_sexf_logit     <- log(1.02)    ; label("Log-odds shift on the response logit for female sex vs male reference (unitless logit)")          # Fukae 2024 Table 2, OR 1.02 (0.788, 1.63); ESS 6066, Rhat 1.00
    e_ecog_ge1_logit <- log(0.903)   ; label("Log-odds shift on the response logit for ECOG PS >= 1 vs ECOG PS 0 reference (unitless logit)")    # Fukae 2024 Table 2, OR 0.903 (0.180, 1.18); ESS 2582, Rhat 1.00

    # ----- Covariate effects on the exposure slope (beta*_2) -----
    e_age_slope      <- log(0.948)   ; label("Shift in the unbound-AUCss exposure slope per 10-year increase in age (unitless logit)")          # Fukae 2024 Table 2, OR 0.948 (0.277, 1.12); ESS 1974, Rhat 1.00
    e_ldh_slope      <- log(1.02)    ; label("Shift in the unbound-AUCss exposure slope per 300 U/L increase in baseline LDH (unitless logit)") # Fukae 2024 Table 2, OR 1.02 (0.843, 1.48); ESS 5589, Rhat 1.00
    e_wt_slope       <- log(0.995)   ; label("Shift in the unbound-AUCss exposure slope per 20 kg increase in body weight (unitless logit)")    # Fukae 2024 Table 2, OR 0.995 (0.806, 1.21); ESS 10667, Rhat 1.00
    e_sexf_slope     <- log(1.01)    ; label("Shift in the unbound-AUCss exposure slope for female sex vs male reference (unitless logit)")     # Fukae 2024 Table 2, OR 1.01 (0.836, 1.33); ESS 7099, Rhat 1.00
    e_ecog_ge1_slope <- log(0.998)   ; label("Shift in the unbound-AUCss exposure slope for ECOG PS >= 1 vs ECOG PS 0 reference (unitless logit)")  # Fukae 2024 Table 2, OR 0.998 (0.797, 1.25); ESS 9591, Rhat 1.00

    # ----- No between-subject variability, no residual error -----
    # See Fukae_2024_valemetostat_orr_central.R for the full rationale:
    # the source model is a Bernoulli-likelihood logistic regression
    # with no omega and no sigma, so the deterministic probability is
    # emitted with a tiny placeholder additive residual.
    addSd_prob_orr_investigator <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    zauc  <- (AUCU_VALE - 375) / 250
    zage  <- (AGE - 65) / 10
    zldh  <- (LDH - 250) / 300
    zwt   <- (WT - 63) / 20

    # ----- Covariate effect on the logit intercept (x' beta*_1) -----
    cov_logit <- e_age_logit  * zage +
                 e_ldh_logit  * zldh +
                 e_wt_logit   * zwt  +
                 e_sexf_logit * SEXF +
                 e_ecog_ge1_logit * ECOG_GE1

    # ----- Covariate effect on the exposure slope (x' beta*_2) -----
    cov_slope <- e_age_slope  * zage +
                 e_ldh_slope  * zldh +
                 e_wt_slope   * zwt  +
                 e_sexf_slope * SEXF +
                 e_ecog_ge1_slope * ECOG_GE1

    # ----- Linear predictor (Fukae 2024 Methods equation) -----
    # logit P(Y = 1 | E, x) = mu + beta*_E E + x' beta*_1 + x' beta*_2 E
    logit_orr_investigator <- logit_ref +
                              e_aucu_logit * zauc +
                              cov_logit +
                              cov_slope * zauc

    prob_orr_investigator <- expit(logit_orr_investigator)

    # ----- Observation -----
    prob_orr_investigator ~ add(addSd_prob_orr_investigator)
  })
}
