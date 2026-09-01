Fukae_2024_valemetostat_dose_interruption <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-safety model for DOSE ",
    "INTERRUPTION due to a treatment-emergent adverse event in adults ",
    "with relapsed/refractory non-Hodgkin lymphoma, including adult ",
    "T-cell leukemia/lymphoma, treated with the oral EZH1/EZH2 dual ",
    "inhibitor valemetostat (Fukae 2024, n = 102 pooled across the ",
    "J101 phase I and J201 phase II trials, 150-300 mg once daily). ",
    "A tolerability rather than a toxicity endpoint: it measures ",
    "whether the assigned regimen could be maintained. The ",
    "probability of the event is expit(mu + bE*zE + x'b1 + x'b2*zE) ",
    "where zE is the centred-and-scaled unbound valemetostat ",
    "steady-state AUC (AUCU_VALE - 375)/250 ng*h/mL, x is the baseline ",
    "covariate vector (age, LDH and weight centred-and-scaled; sex, ",
    "mild hepatic impairment, White race, Black race and United States ",
    "enrollment as 0/1 indicators), b1 are covariate effects on the ",
    "logit intercept and b2 are covariate effects on the exposure ",
    "slope (Fukae 2024 Methods equation, Table 3). All candidate ",
    "covariate effects were regularized via normal-mixture ",
    "spike-and-slab priors (spike Normal(0, 0.1), slab Normal(0, 2.5), ",
    "mixing weight Beta(1, 1)); only the intercept and the exposure ",
    "slope carried weakly informative unregularized priors. Fitted ",
    "with CmdStanR 0.4.0 (Hamiltonian Monte Carlo, no-U-turn sampler). ",
    "There is no PK layer and no ODE: individual unbound AUCss is ",
    "supplied as data from the companion population PK model ",
    "(doi:10.1002/psp4.13201). No between-subject random effect and no ",
    "residual error are estimated (Bernoulli likelihood). Paired with ",
    "Fukae_2024_valemetostat_dose_reduction; five companion ",
    "exposure-safety models and two companion exposure-efficacy models ",
    "in the Fukae_2024_valemetostat_* family."
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
    time          = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUCU_VALE covariate column)",
    concentration = "prob_dose_interruption (probability of Dose interruption due to TEAE, 0-1; also logit_dose_interruption)"
  )

  covariateData <- list(
    AUCU_VALE = list(
      description        = "Unbound (free) valemetostat plasma AUC over the once-daily 24 h dosing interval at steady state, per subject. Supplied as data: this model has no PK layer, and the source analysis used empirical-Bayes individual predictions from the companion population PK model (Fukae 2024, doi:10.1002/psp4.13201).",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 375 ng*h/mL and scaled by 250 ng*h/mL inside model() (Fukae 2024 Table 3 footnote reference patient; Table 3 row 'Unbound valemetostat AUCSS: 250 ng*h/mL increase'). Observed 5th-95th percentile 184-887 ng*h/mL, which is the paper's modified region of practical equivalence. Unbound, not total.",
      source_name        = "unbound valemetostat AUCSS"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 65 years and scaled by 10 years inside model(). R/R NHL safety analysis set median 69.0 years, range 37-88 (Fukae 2024 Table 1).",
      source_name        = "Age"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase concentration.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 250 U/L and scaled by 300 U/L inside model(). R/R NHL safety analysis set median 244 U/L, range 119-2000 (Fukae 2024 Table 1).",
      source_name        = "LDH"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 63 kg and scaled by 20 kg inside model(). R/R NHL safety analysis set median 63.4 kg, range 34.5-114 (Fukae 2024 Table 1).",
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male. Binary covariates are NOT centred or scaled -- Fukae 2024 Methods applies centring and scaling to 'continuous covariates only'.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the reference patient is male per the Table 3 footnote)",
      notes              = "R/R NHL safety analysis set 44.1% female (Fukae 2024 Table 1).",
      source_name        = "Sex: female"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator; 1 = mild impairment, 0 = normal hepatic function. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; the reference patient has normal hepatic function per the Table 3 footnote)",
      notes              = "Fukae 2024 dichotomizes hepatic function as normal versus mild impairment; no moderate or severe stratum was enrolled. R/R NHL safety analysis set 13.7% mild (Table 1).",
      source_name        = "Hepatic function: mild impairment"
    ),
    RACE_WHITE = list(
      description        = "White race indicator; 1 = White, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian, when RACE_BLACK is also 0; the reference patient is Asian per the Table 3 footnote)",
      notes              = "Paired with RACE_BLACK to encode the paper's three-level race covariate (Asian / White / Black) with two binary indicators; both 0 gives the Asian reference. R/R NHL safety analysis set Asian 72.5%, White 21.6%, Black 5.9% (Fukae 2024 Table 1). Distinct from REGION_USA, which records enrollment country -- Fukae 2024 carries both simultaneously.",
      source_name        = "Race: White"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator; 1 = Black, 0 = otherwise. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Asian, when RACE_WHITE is also 0)",
      notes              = "Paired with RACE_WHITE; see that entry. Only 6 of 102 patients were Black (5.9%), which is why several Black-race credible intervals in Table 3 are extremely wide.",
      source_name        = "Race: Black"
    ),
    REGION_USA = list(
      description        = "United States enrollment indicator; 1 = enrolled and treated at a United States study site, 0 = enrolled in Japan. Not centred or scaled (binary).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Japan; the reference patient is 'an Asian male in Japan' per the Table 3 footnote)",
      notes              = "The pooled cohort has only two regions, so this is the exact complement of a Japan indicator. R/R NHL safety analysis set Japan 69.6%, United States 30.4% (Fukae 2024 Table 1). Distinct from the RACE_* indicators: 3 of the 102 patients were Asian enrolled in the United States, so region and race are not collinear.",
      source_name        = "Geographic region: United States"
    )
  )

  covariatesDataExcluded <- list(
    DIS_NHL_STAGE4 = list(
      description = "Non-Hodgkin lymphoma disease stage indicator; 1 = stage IV, 0 = stage I/II/III.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Named by Fukae 2024 Methods as a candidate categorical covariate for the exposure-safety analysis ('NHL disease stage (I-III/IV)'), but NOT reported in Table 3 and therefore not retained in the final model -- no point estimate exists anywhere on disk. Stage was heavily missing in the safety cohort (24 of 102 patients, 23.5%; Table 1), the likely reason it was dropped. Documented here to preserve the covariate screen without carrying a convention warning."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 102L,
    n_studies      = 2L,
    n_observations = "90 evaluable binary dose-interruption records (one per patient; 12 of the 102 patients had a missing outcome, Fukae 2024 Table 1)",
    age_range      = "median 69.0 years, range 37-88 (Fukae 2024 Table 1)",
    weight_range   = "median 63.4 kg, range 34.5-114 (Fukae 2024 Table 1)",
    sex_female_pct = 44.1,
    race_ethnicity = c(Asian = 72.5, White = 21.6, Black = 5.9),
    disease_state  = "relapsed or refractory non-Hodgkin lymphomas, including adult T-cell leukemia/lymphoma and other peripheral T-cell lymphomas",
    dose_range     = "valemetostat 150, 200, 250 or 300 mg orally once daily (J101 dose escalation and expansion; J201 was 200 mg throughout)",
    regions        = "Japan and the United States",
    notes          = "Observed outcome: 45 of 102 patients (44.1%) had a dose interruption due to a TEAE, 45 (44.1%) did not, and 12 (11.8%) were missing (Fukae 2024 Table 1). Paired with the dose-reduction endpoint; interruption is the more frequent and less severe of the two dose-modification endpoints."
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Fukae 2024 Table 3, column
    # "Dose interruption due to TEAE". The Table 3 footnote states the reporting
    # convention explicitly: "All estimates are expressed as odds
    # ratios, except the probability of a safety event for a reference
    # patient. Effects are described as exp(beta)." Each ini() value
    # below is therefore log(printed odds ratio), written in log() form
    # so the published number is visible at the trace site.
    #
    # Reference patient (Table 3 footnote): a 65-year-old Asian male in
    # Japan, weighing 63 kg, with normal hepatic function, baseline LDH
    # 250 U/L and unbound valemetostat AUCss 375 ng*h/mL. Centring
    # constants = those reference values; scaling constants = the
    # per-row increments, which match the rounded cohort standard
    # deviations in Table 1 as Fukae 2024 Methods describes.
    # ==================================================================

    # ----- Logit intercept (reference-patient event probability) -----
    logit_ref             <- log(0.514 / (1 - 0.514)) ; label("Logit of the dose interruption due to TEAE probability for the reference patient (unitless logit)")  # Fukae 2024 Table 3, "Population mean" 0.514 (0.404, 0.626); ESS 8021, Rhat 1.00

    # ----- Exposure slope on the logit (beta*_E) -----
    e_aucu_logit          <- log(1.32)     ; label("Log-odds of dose interruption due to TEAE per 250 ng*h/mL increase in unbound valemetostat AUCss (unitless logit)")  # Fukae 2024 Table 3, OR 1.32 (0.813, 2.15); ESS 5114, Rhat 1.00

    # ----- Covariate effects on the logit intercept (beta*_1) -----
    e_age_logit           <- log(0.967)    ; label("Log-odds shift on the dose interruption due to TEAE logit per 10-year increase in age (unitless logit)")  # Fukae 2024 Table 3, OR 0.967 (0.810, 1.13); ESS 13558, Rhat 1.00
    e_ldh_logit           <- log(0.991)    ; label("Log-odds shift on the dose interruption due to TEAE logit per 300 U/L increase in baseline LDH (unitless logit)")  # Fukae 2024 Table 3, OR 0.991 (0.827, 1.17); ESS 10754, Rhat 1.00
    e_wt_logit            <- log(1.01)     ; label("Log-odds shift on the dose interruption due to TEAE logit per 20 kg increase in body weight (unitless logit)")  # Fukae 2024 Table 3, OR 1.01 (0.852, 1.20); ESS 11344, Rhat 1.00
    e_sexf_logit          <- log(0.970)    ; label("Log-odds shift on the dose interruption due to TEAE logit for female sex vs male reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.970 (0.731, 1.25); ESS 9919, Rhat 1.00
    e_hepimp_mild_logit   <- log(0.992)    ; label("Log-odds shift on the dose interruption due to TEAE logit for mild hepatic impairment vs normal reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.992 (0.667, 1.50); ESS 11554, Rhat 1.00
    e_race_white_logit    <- log(1.03)     ; label("Log-odds shift on the dose interruption due to TEAE logit for White race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.03 (0.720, 1.60); ESS 6132, Rhat 1.00
    e_race_black_logit    <- log(0.836)    ; label("Log-odds shift on the dose interruption due to TEAE logit for Black race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.836 (0.160, 1.90); ESS 4660, Rhat 1.00
    e_region_usa_logit    <- log(1.04)     ; label("Log-odds shift on the dose interruption due to TEAE logit for United States vs Japan enrollment reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.04 (0.766, 1.85); ESS 3443, Rhat 1.00

    # ----- Covariate effects on the exposure slope (beta*_2) -----
    e_age_slope           <- log(0.964)    ; label("Shift in the unbound-AUCss exposure slope per 10-year increase in age (unitless logit)")  # Fukae 2024 Table 3, OR 0.964 (0.792, 1.12); ESS 9222, Rhat 1.00
    e_ldh_slope           <- log(0.972)    ; label("Shift in the unbound-AUCss exposure slope per 300 U/L increase in baseline LDH (unitless logit)")  # Fukae 2024 Table 3, OR 0.972 (0.750, 1.16); ESS 5358, Rhat 1.00
    e_wt_slope            <- log(1.03)     ; label("Shift in the unbound-AUCss exposure slope per 20 kg increase in body weight (unitless logit)")  # Fukae 2024 Table 3, OR 1.03 (0.855, 1.25); ESS 8317, Rhat 1.00
    e_sexf_slope          <- log(1.09)     ; label("Shift in the unbound-AUCss exposure slope for female sex vs male reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.09 (0.863, 3.31); ESS 2482, Rhat 1.00
    e_hepimp_mild_slope   <- log(0.907)    ; label("Shift in the unbound-AUCss exposure slope for mild hepatic impairment vs normal reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.907 (0.380, 1.22); ESS 4026, Rhat 1.00
    e_race_white_slope    <- log(1.20)     ; label("Shift in the unbound-AUCss exposure slope for White race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.20 (0.712, 10.5); ESS 3002, Rhat 1.00
    e_race_black_slope    <- log(1.76)     ; label("Shift in the unbound-AUCss exposure slope for Black race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.76 (0.609, 7.83e5); ESS 1561, Rhat 1.00
    e_region_usa_slope    <- log(1.21)     ; label("Shift in the unbound-AUCss exposure slope for United States vs Japan enrollment reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.21 (0.753, 15.2); ESS 2555, Rhat 1.00

    # ----- No between-subject variability, no residual error -----
    # See Fukae_2024_valemetostat_orr_central.R for the full rationale.
    addSd_prob_dose_interruption <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value dose interruption due to TEAE probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
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
                 e_hepimp_mild_logit * HEPIMP_MILD +
                 e_race_white_logit  * RACE_WHITE +
                 e_race_black_logit  * RACE_BLACK +
                 e_region_usa_logit  * REGION_USA

    # ----- Covariate effect on the exposure slope (x' beta*_2) -----
    cov_slope <- e_age_slope  * zage +
                 e_ldh_slope  * zldh +
                 e_wt_slope   * zwt  +
                 e_sexf_slope * SEXF +
                 e_hepimp_mild_slope * HEPIMP_MILD +
                 e_race_white_slope  * RACE_WHITE +
                 e_race_black_slope  * RACE_BLACK +
                 e_region_usa_slope  * REGION_USA

    # ----- Linear predictor (Fukae 2024 Methods equation) -----
    # logit P(Y = 1 | E, x) = mu + beta*_E E + x' beta*_1 + x' beta*_2 E
    logit_dose_interruption <- logit_ref +
                    e_aucu_logit * zauc +
                    cov_logit +
                    cov_slope * zauc

    prob_dose_interruption <- expit(logit_dose_interruption)

    # ----- Observation -----
    prob_dose_interruption ~ add(addSd_prob_dose_interruption)
  })
}

