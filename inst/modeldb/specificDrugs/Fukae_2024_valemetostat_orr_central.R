Fukae_2024_valemetostat_orr_central <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-response model for overall ",
    "response rate (ORR) by CENTRAL assessment in adults with ",
    "relapsed/refractory adult T-cell leukemia/lymphoma (ATLL) treated ",
    "with the oral EZH1/EZH2 dual inhibitor valemetostat (Fukae 2024, ",
    "n = 25, J201 phase II, 200 mg once daily). The probability of ",
    "response is expit(mu + bE*zE + x'b1 + x'b2*zE) where zE is the ",
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
    ">= 1000 with Gelman-Rubin statistic < 1.05. There is no PK layer ",
    "and no ODE: individual unbound AUCss is supplied as data from the ",
    "companion population PK model (doi:10.1002/psp4.13201). No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Companion investigator-assessed model in ",
    "Fukae_2024_valemetostat_orr_investigator; six companion ",
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
    concentration = "prob_orr_central (probability of overall response, 0-1; also logit_orr_central)"
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
      notes              = "The paper collapses ECOG PS to '0' versus '1+' (Fukae 2024 Methods and Table 1). This is the strongest covariate effect in the central-assessment model: odds ratio 0.378 on the logit intercept, i.e. a 62% lower odds of response for ECOG PS >= 1 at the reference exposure. ATLL efficacy analysis set 46.2% with ECOG PS >= 1 (Table 1).",
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
    n_subjects     = 25L,
    n_studies      = 1L,
    n_observations = "25 binary response records (one per patient; landmark analysis, no repeated measures)",
    age_range      = "ATLL efficacy analysis set (n = 39 including the one patient with missing ORR): median 69.0 years, range 37-84 (Fukae 2024 Table 1). J201-only central-assessment subset: median 69.0, range 59-84.",
    weight_range   = "J201 subset median 57.7 kg, range 34.5-82.9 (Fukae 2024 Table 1)",
    sex_female_pct = 52.0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "relapsed or refractory adult T-cell leukemia/lymphoma (ATLL); responses per modified 2009 ATLL response criteria, ORR = complete response (including uncertified CR) or partial response",
    dose_range     = "valemetostat 200 mg orally once daily until progressive disease or unacceptable toxicity",
    regions        = "Japan (J201, DS3201-A-J201, NCT04102150)",
    notes          = paste0(
      "Central assessment was performed only in the phase II J201 trial, ",
      "so this model's analysis set is the 25 J201 patients; the 14 J101 ",
      "ATLL patients have missing central ORR (Fukae 2024 Table 1). ",
      "Observed outcome: 12 responders (ORR 48.0%) and 13 non-responders ",
      "(Fukae 2024 Results, Patient populations). The companion ",
      "investigator-assessed model uses the larger 38-patient pooled set."
    )
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Fukae 2024 Table 2, column
    # "ORR by central assessment". The table reports every effect as an
    # ODDS RATIO -- the footnote to the companion Table 3 states the
    # convention explicitly: "Effects are described as exp(beta)". Each
    # ini() value below is therefore log(printed odds ratio), written in
    # log() form so the published number is visible at the trace site.
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
    # Table 2 reports the reference-patient PROBABILITY (0.619), not the
    # logit, so mu is recovered as log(p / (1 - p)). Cross-check: the
    # Results narrative quotes "a predicted 62% (95% CrI, 37%-87%)
    # probability of achieving clinical response in the reference
    # patient", matching 0.619 (0.373, 0.867).
    logit_ref        <- log(0.619 / (1 - 0.619)) ; label("Logit of the response probability for the reference patient (unitless logit)")  # Fukae 2024 Table 2, "Population mean" 0.619 (0.373, 0.867); ESS 1336, Rhat 1.00

    # ----- Exposure slope on the logit (beta*_E) -----
    e_aucu_logit     <- log(1.08)    ; label("Log-odds of response per 250 ng*h/mL increase in unbound valemetostat AUCss (unitless logit)")   # Fukae 2024 Table 2, OR 1.08 (0.396, 3.05); ESS 5543, Rhat 1.00

    # ----- Covariate effects on the logit intercept (beta*_1) -----
    e_age_logit      <- log(0.976)   ; label("Log-odds shift on the response logit per 10-year increase in age (unitless logit)")               # Fukae 2024 Table 2, OR 0.976 (0.259, 2.01); ESS 6800, Rhat 1.00
    e_ldh_logit      <- log(0.899)   ; label("Log-odds shift on the response logit per 300 U/L increase in baseline LDH (unitless logit)")      # Fukae 2024 Table 2, OR 0.899 (0.119, 1.11); ESS 1261, Rhat 1.00
    e_wt_logit       <- log(1.17)    ; label("Log-odds shift on the response logit per 20 kg increase in body weight (unitless logit)")         # Fukae 2024 Table 2, OR 1.17 (0.836, 9.23); ESS 1473, Rhat 1.00
    e_sexf_logit     <- log(0.969)   ; label("Log-odds shift on the response logit for female sex vs male reference (unitless logit)")          # Fukae 2024 Table 2, OR 0.969 (0.240, 2.06); ESS 5680, Rhat 1.00
    e_ecog_ge1_logit <- log(0.378)   ; label("Log-odds shift on the response logit for ECOG PS >= 1 vs ECOG PS 0 reference (unitless logit)")    # Fukae 2024 Table 2, OR 0.378 (0.0396, 1.13); ESS 1020, Rhat 1.00

    # ----- Covariate effects on the exposure slope (beta*_2) -----
    # These multiply the covariate by the centred-and-scaled exposure, so
    # they shift the steepness of the exposure-response relationship
    # rather than its level at the reference exposure.
    e_age_slope      <- log(1.05)    ; label("Shift in the unbound-AUCss exposure slope per 10-year increase in age (unitless logit)")          # Fukae 2024 Table 2, OR 1.05 (0.839, 5.02); ESS 1720, Rhat 1.00
    e_ldh_slope      <- log(1.00)    ; label("Shift in the unbound-AUCss exposure slope per 300 U/L increase in baseline LDH (unitless logit)") # Fukae 2024 Table 2, OR 1.00 (0.439, 1.93); ESS 8084, Rhat 1.00
    e_wt_slope       <- log(0.951)   ; label("Shift in the unbound-AUCss exposure slope per 20 kg increase in body weight (unitless logit)")    # Fukae 2024 Table 2, OR 0.951 (0.245, 1.20); ESS 2434, Rhat 1.00
    e_sexf_slope     <- log(0.991)   ; label("Shift in the unbound-AUCss exposure slope for female sex vs male reference (unitless logit)")     # Fukae 2024 Table 2, OR 0.991 (0.330, 2.28); ESS 5164, Rhat 1.00
    e_ecog_ge1_slope <- log(0.989)   ; label("Shift in the unbound-AUCss exposure slope for ECOG PS >= 1 vs ECOG PS 0 reference (unitless logit)")  # Fukae 2024 Table 2, OR 0.989 (0.343, 1.87); ESS 5630, Rhat 1.00

    # ----- No between-subject variability, no residual error -----
    # The source model is a Bernoulli-likelihood logistic regression: the
    # covariates and exposure fully determine each subject's response
    # probability, and no omega or sigma is estimated (Fukae 2024
    # Methods; Table 2 lists only the fixed effects). rxode2 requires an
    # observation declaration, so the deterministic probability is
    # emitted with a tiny placeholder additive residual, mirroring the
    # Oniki_2018_nafld_risk.R and Hansson_2013c_sunitinib.R pattern.
    # This does not perturb the predicted probability; see the vignette
    # Assumptions and deviations section.
    addSd_prob_orr_central <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    # Fukae 2024 Methods: "both E (exposure) and x (continuous
    # covariates only) were centered and scaled at approximate
    # population median and standard deviation values". Binary
    # covariates (SEXF, ECOG_GE1) are used untransformed.
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
    logit_orr_central <- logit_ref +
                         e_aucu_logit * zauc +
                         cov_logit +
                         cov_slope * zauc

    prob_orr_central <- expit(logit_orr_central)

    # ----- Observation -----
    # Deterministic probability of overall response by central
    # assessment. Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_orr_central) on the rxSolve output.
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
