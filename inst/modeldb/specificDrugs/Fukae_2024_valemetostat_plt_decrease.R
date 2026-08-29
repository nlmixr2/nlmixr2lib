Fukae_2024_valemetostat_plt_decrease <- function() {
  description <- paste0(
    "Bayesian logistic-regression exposure-safety model for CTCAE ",
    "grade >= 3 decrease in platelet count (assessed from laboratory ",
    "platelet values) in adults with relapsed/refractory non-Hodgkin ",
    "lymphoma, including adult T-cell leukemia/lymphoma, treated with ",
    "the oral EZH1/EZH2 dual inhibitor valemetostat (Fukae 2024, ",
    "n = 102 pooled across the J101 phase I and J201 phase II trials, ",
    "150-300 mg once daily). The probability of the event is ",
    "expit(mu + bE*zE + x'b1 + x'b2*zE) where zE is the ",
    "centred-and-scaled unbound valemetostat steady-state AUC ",
    "(AUCU_VALE - 375)/250 ng*h/mL, x is the baseline covariate vector ",
    "(age, LDH, weight and baseline platelet count ",
    "centred-and-scaled; sex, mild hepatic impairment, White race, ",
    "Black race and United States enrollment as 0/1 indicators), b1 ",
    "are covariate effects on the logit intercept and b2 are covariate ",
    "effects on the exposure slope (Fukae 2024 Methods equation, ",
    "Table 3). Candidate covariate effects were regularized via ",
    "normal-mixture spike-and-slab priors (spike Normal(0, 0.1), slab ",
    "Normal(0, 2.5), mixing weight Beta(1, 1)); the intercept, the ",
    "exposure slope and the BASELINE PLATELET effect on the intercept ",
    "were NOT regularized (slab-only prior). This is the steepest of ",
    "the paper's six exposure-safety relationships (odds ratio 2.03 ",
    "per 250 ng*h/mL of unbound AUCss, the only exposure effect whose ",
    "95% credible interval excludes 1) and it carries the largest ",
    "single covariate effect in the paper (baseline platelet odds ",
    "ratio 0.230 per 100 x 10^9/L increase). Fitted with CmdStanR ",
    "0.4.0 (Hamiltonian Monte Carlo, no-U-turn sampler). There is no ",
    "PK layer and no ODE: individual unbound AUCss is supplied as data ",
    "from the companion population PK model (doi:10.1002/psp4.13201). ",
    "No between-subject random effect and no residual error are ",
    "estimated (Bernoulli likelihood). Five companion exposure-safety ",
    "models and two companion exposure-efficacy models in the ",
    "Fukae_2024_valemetostat_* family."
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
    concentration = "prob_plt_decrease (probability of Grade >= 3 PLT decrease, 0-1; also logit_plt_decrease)"
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
    PLT = list(
      description        = "Baseline platelet count -- the laboratory value corresponding to this model's platelet-decrease endpoint.",
      units              = "10^9/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 200 x 10^9/L and scaled by 100 x 10^9/L inside model(). The 200 x 10^9/L reference is stated in the Fukae 2024 Results narrative ('Compared with the reference baseline value (200 x 10^9/L), a 100 x 10^9/L increase in baseline PLT') rather than in the Table 3 footnote, which omits the laboratory reference values. R/R NHL safety analysis set median 194 x 10^9/L, range 61.0-958 (Table 1). This effect was deliberately NOT regularized (slab-only prior) per Fukae 2024 Methods, and is the largest single covariate effect in the paper: odds ratio 0.230 per 100 x 10^9/L increase, a 77% reduction in the odds of a grade >= 3 platelet decrease.",
      source_name        = "Platelets"
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
    n_observations = "102 binary platelet-decrease records (one per patient; no missing outcomes for this endpoint, Fukae 2024 Table 1)",
    age_range      = "median 69.0 years, range 37-88 (Fukae 2024 Table 1)",
    weight_range   = "median 63.4 kg, range 34.5-114 (Fukae 2024 Table 1)",
    sex_female_pct = 44.1,
    race_ethnicity = c(Asian = 72.5, White = 21.6, Black = 5.9),
    disease_state  = "relapsed or refractory non-Hodgkin lymphomas, including adult T-cell leukemia/lymphoma and other peripheral T-cell lymphomas",
    dose_range     = "valemetostat 150, 200, 250 or 300 mg orally once daily (J101 dose escalation and expansion; J201 was 200 mg throughout)",
    regions        = "Japan and the United States",
    notes          = "Grade >= 3 PLT decreased was assessed from laboratory platelet values rather than from investigator-reported adverse-event terms (Fukae 2024 Table 1 note). Observed outcome: 25 of 102 patients (24.5%) had a grade >= 3 platelet decrease and 77 (75.5%) did not, with no missing values (Table 1)."
  )

  ini({
    # ==================================================================
    # All values are posterior medians from Fukae 2024 Table 3, column
    # "Grade >= 3 PLT decrease". The Table 3 footnote states the reporting
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
    logit_ref             <- log(0.165 / (1 - 0.165)) ; label("Logit of the grade >= 3 platelet decrease probability for the reference patient (unitless logit)")  # Fukae 2024 Table 3, "Population mean" 0.165 (0.0863, 0.257); ESS 2999, Rhat 1.00

    # ----- Exposure slope on the logit (beta*_E) -----
    e_aucu_logit          <- log(2.03)     ; label("Log-odds of grade >= 3 platelet decrease per 250 ng*h/mL increase in unbound valemetostat AUCss (unitless logit)")  # Fukae 2024 Table 3, OR 2.03 (1.23, 3.68); ESS 6209, Rhat 1.00

    # ----- Covariate effects on the logit intercept (beta*_1) -----
    e_age_logit           <- log(0.961)    ; label("Log-odds shift on the grade >= 3 platelet decrease logit per 10-year increase in age (unitless logit)")  # Fukae 2024 Table 3, OR 0.961 (0.764, 1.14); ESS 9788, Rhat 1.00
    e_ldh_logit           <- log(1.03)     ; label("Log-odds shift on the grade >= 3 platelet decrease logit per 300 U/L increase in baseline LDH (unitless logit)")  # Fukae 2024 Table 3, OR 1.03 (0.890, 1.22); ESS 12892, Rhat 1.00
    e_wt_logit            <- log(0.994)    ; label("Log-odds shift on the grade >= 3 platelet decrease logit per 20 kg increase in body weight (unitless logit)")  # Fukae 2024 Table 3, OR 0.994 (0.826, 1.20); ESS 13408, Rhat 1.00
    e_sexf_logit          <- log(1.06)     ; label("Log-odds shift on the grade >= 3 platelet decrease logit for female sex vs male reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.06 (0.828, 1.57); ESS 7790, Rhat 1.00
    e_plt_logit           <- log(0.230)    ; label("Log-odds shift on the grade >= 3 platelet decrease logit per 100 x 10^9/L increase in baseline platelet count (unitless logit)")  # Fukae 2024 Table 3, OR 0.230 (0.0877, 0.483); ESS 2001, Rhat 1.00; slab-only prior (not regularized)
    e_hepimp_mild_logit   <- log(1.12)     ; label("Log-odds shift on the grade >= 3 platelet decrease logit for mild hepatic impairment vs normal reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.12 (0.725, 1.96); ESS 9193, Rhat 1.00
    e_race_white_logit    <- log(1.01)     ; label("Log-odds shift on the grade >= 3 platelet decrease logit for White race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.01 (0.699, 1.58); ESS 9760, Rhat 1.00
    e_race_black_logit    <- log(0.945)    ; label("Log-odds shift on the grade >= 3 platelet decrease logit for Black race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.945 (0.370, 1.95); ESS 6004, Rhat 1.00
    e_region_usa_logit    <- log(0.954)    ; label("Log-odds shift on the grade >= 3 platelet decrease logit for United States vs Japan enrollment reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.954 (0.486, 1.30); ESS 5847, Rhat 1.00

    # ----- Covariate effects on the exposure slope (beta*_2) -----
    e_age_slope           <- log(0.983)    ; label("Shift in the unbound-AUCss exposure slope per 10-year increase in age (unitless logit)")  # Fukae 2024 Table 3, OR 0.983 (0.823, 1.16); ESS 12551, Rhat 1.00
    e_ldh_slope           <- log(1.02)     ; label("Shift in the unbound-AUCss exposure slope per 300 U/L increase in baseline LDH (unitless logit)")  # Fukae 2024 Table 3, OR 1.02 (0.852, 1.27); ESS 9304, Rhat 1.00
    e_wt_slope            <- log(1.03)     ; label("Shift in the unbound-AUCss exposure slope per 20 kg increase in body weight (unitless logit)")  # Fukae 2024 Table 3, OR 1.03 (0.849, 1.27); ESS 10392, Rhat 1.00
    e_sexf_slope          <- log(1.03)     ; label("Shift in the unbound-AUCss exposure slope for female sex vs male reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.03 (0.816, 1.38); ESS 7698, Rhat 1.00
    e_plt_slope           <- log(1.14)     ; label("Shift in the unbound-AUCss exposure slope per 100 x 10^9/L increase in baseline platelet count (unitless logit)")  # Fukae 2024 Table 3, OR 1.14 (0.908, 6.09); ESS 1060, Rhat 1.01; the only Rhat above 1.00 anywhere in Fukae 2024
    e_hepimp_mild_slope   <- log(0.872)    ; label("Shift in the unbound-AUCss exposure slope for mild hepatic impairment vs normal reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.872 (0.235, 1.20); ESS 2626, Rhat 1.00
    e_race_white_slope    <- log(0.992)    ; label("Shift in the unbound-AUCss exposure slope for White race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 0.992 (0.419, 1.70); ESS 4667, Rhat 1.00
    e_race_black_slope    <- log(1.47)     ; label("Shift in the unbound-AUCss exposure slope for Black race vs Asian reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.47 (0.558, 3.09e4); ESS 2482, Rhat 1.00
    e_region_usa_slope    <- log(1.21)     ; label("Shift in the unbound-AUCss exposure slope for United States vs Japan enrollment reference (unitless logit)")  # Fukae 2024 Table 3, OR 1.21 (0.752, 51.6); ESS 1741, Rhat 1.00

    # ----- No between-subject variability, no residual error -----
    # See Fukae_2024_valemetostat_orr_central.R for the full rationale.
    addSd_prob_plt_decrease <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value grade >= 3 platelet decrease probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Centre and scale the continuous predictors -----
    zauc  <- (AUCU_VALE - 375) / 250
    zage  <- (AGE - 65) / 10
    zldh  <- (LDH - 250) / 300
    zwt   <- (WT - 63) / 20
    zplt  <- (PLT - 200) / 100

    # ----- Covariate effect on the logit intercept (x' beta*_1) -----
    cov_logit <- e_age_logit  * zage +
                 e_ldh_logit  * zldh +
                 e_wt_logit   * zwt  +
                 e_plt_logit  * zplt +
                 e_sexf_logit * SEXF +
                 e_hepimp_mild_logit * HEPIMP_MILD +
                 e_race_white_logit  * RACE_WHITE +
                 e_race_black_logit  * RACE_BLACK +
                 e_region_usa_logit  * REGION_USA

    # ----- Covariate effect on the exposure slope (x' beta*_2) -----
    cov_slope <- e_age_slope  * zage +
                 e_ldh_slope  * zldh +
                 e_wt_slope   * zwt  +
                 e_plt_slope  * zplt +
                 e_sexf_slope * SEXF +
                 e_hepimp_mild_slope * HEPIMP_MILD +
                 e_race_white_slope  * RACE_WHITE +
                 e_race_black_slope  * RACE_BLACK +
                 e_region_usa_slope  * REGION_USA

    # ----- Linear predictor (Fukae 2024 Methods equation) -----
    # logit P(Y = 1 | E, x) = mu + beta*_E E + x' beta*_1 + x' beta*_2 E
    logit_plt_decrease <- logit_ref +
                    e_aucu_logit * zauc +
                    cov_logit +
                    cov_slope * zauc

    prob_plt_decrease <- expit(logit_plt_decrease)

    # ----- Observation -----
    prob_plt_decrease ~ add(addSd_prob_plt_decrease)
  })
}

