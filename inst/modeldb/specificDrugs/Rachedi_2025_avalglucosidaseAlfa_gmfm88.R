Rachedi_2025_avalglucosidaseAlfa_gmfm88 <- function() {
  description <- paste(
    "Nonlinear random-coefficient (NLRC) model of the change from baseline in",
    "the Gross Motor Function Measure-88 (GMFM-88) total percent score in",
    "infantile-onset (IOPD, Mini-COMET) and late-onset (LOPD, COMET) Pompe",
    "disease treated with enzyme replacement therapy (avalglucosidase alfa or",
    "alglucosidase alfa) (Rachedi 2025). No ODE and no dose: the treatment",
    "effect enters only through the observed or simulated percent change from",
    "baseline in urinary hexose tetrasaccharide (uHex4), supplied as a",
    "time-varying covariate, via an exponential asymptote",
    "max * (1 - exp(-rate * (-%dUHex4) * 0.01)) with a log-normal",
    "between-subject random rate. IOPD and LOPD were fitted jointly with",
    "population-specific intercepts, baseline-covariate effects and",
    "asymptotes and a shared random rate and residual error.",
    sep = " "
  )
  reference <- paste(
    "Rachedi F, Jreich R, Sparks S, Zaher A, An Haack K, Granados A,",
    "Meng Z. Clinical modeling of motor function to predict treatment",
    "efficacy and enable in silico treatment comparisons in infantile-onset",
    "Pompe disease. CPT Pharmacometrics Syst Pharmacol. 2025;14(11):1797.",
    "doi:10.1002/psp4.13287.",
    sep = " "
  )
  vignette <- "Rachedi_2025_avalglucosidaseAlfa"
  units <- list(
    time = "week",
    dosing = "(none; ERT effect enters through the uHex4 percent-change covariate)",
    concentration = "(observation gmfm88cfb is the GMFM-88 total percent score change from baseline, percentage points)"
  )

  covariateData <- list(
    DIS_INFANTILE_ONSET_POMPE = list(
      description = "Pompe disease phenotype: 1 = infantile-onset (IOPD, Mini-COMET), 0 = late-onset (LOPD, COMET).",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-fixed. Selects the population-specific intercept, baseline-covariate effects, medians and asymptote (Rachedi 2025 main-text IOPD equation vs Supplementary LOPD equation).",
      source_name = "population (IOPD / LOPD)"
    ),
    UHEX4_REL = list(
      description = "Relative change from baseline in urinary hexose tetrasaccharide (uHex4), as a fraction (e.g., -0.50 = 50% decrease).",
      units = "fraction",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying per visit. The paper's %dUHex4 divided by 100; the model term -%dUHex4 * 0.01 is therefore -UHEX4_REL. Supplied from observed uHex4 or from the Kaddi 2024 QSP model.",
      source_name = "%dUHex4"
    ),
    UHEX4_BASE = list(
      description = "Baseline urinary hexose tetrasaccharide (uHex4).",
      units = "mmol/mol creatinine",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Log-transformed and centred at the population median: 53.81 (IOPD) / 8.03 (LOPD) mmol/mol, Figure S2.",
      source_name = "uHex4 (mmol/mol) at baseline"
    ),
    SCORE_GMFM88 = list(
      description = "Baseline GMFM-88 total percent score (0-100).",
      units = "(percent score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Linear, centred at the population median: 67.74 (IOPD) / 84.67 (LOPD), Figure S2. IOPD scored on all five dimensions, LOPD on dimensions D and E (Table 1).",
      source_name = "GMFM-88 total percent score at baseline"
    ),
    AGE_DIAGNOSIS = list(
      description = "Age at diagnosis of Pompe disease.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Log-transformed and centred at the population median: 0.26 (IOPD) / 46.63 (LOPD) years, Figure S2.",
      source_name = "Age at diagnosis in years"
    ),
    T_PRIOR_ALGLUCOSIDASE_ALFA = list(
      description = "Duration of alglucosidase alfa treatment before study entry.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed; IOPD only (LOPD COMET patients were treatment-naive; supply 0). Linear, centred at 6.5 years (see vignette: the pooled IOPD median is not printed).",
      source_name = "Duration of previous ALGLU experience in years"
    ),
    SEXF = list(
      description = "Sex: 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "IOPD only: additive effect on the change score (Table 2 'Gender (female vs. male)').",
      source_name = "gender"
    ),
    CPK = list(
      description = "Baseline serum creatine kinase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "LOPD only; time-fixed baseline. Log-transformed and centred at the LOPD median 550 IU/L (Figure S2).",
      source_name = "Creatine kinase (IU/L) at baseline"
    ),
    CREAT = list(
      description = "Baseline serum creatinine.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "LOPD only; time-fixed baseline. Log-transformed and centred at the LOPD median 52 umol/L (Figure S2).",
      source_name = "Creatinine (umol/L) at baseline"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 121,
    n_studies = 2,
    age_range = "IOPD: 1-12 years at first study infusion (mean 4.7-7.7 by arm); LOPD: adults (mean 46.0-50.4 years by arm)",
    sex_female_pct = "IOPD 47.6 (10/21); LOPD 48.0 (48/100)",
    race_ethnicity = "IOPD: White 57%, Asian 33%, Black 10%; LOPD: White 94%, Asian 3%, Black 3% (Table S2)",
    disease_state = "Pompe disease (acid alpha-glucosidase deficiency): 21 infantile-onset patients previously treated with alglucosidase alfa (Mini-COMET, phase 2) and 100 treatment-naive late-onset patients (COMET, phase 3)",
    dose_range = "Avalglucosidase alfa 20 or 40 mg/kg IV Q2W; alglucosidase alfa 20 mg/kg Q2W (COMET) or current stable dose 20-40 mg/kg QW/Q2W (Mini-COMET)",
    regions = "Multinational",
    notes = "Visits at weeks 13, 25, 37, 49, (61,) 73 and 97 (Table 1). One Mini-COMET patient without post-baseline motor data was excluded."
  )

  ini({
    # Shared between IOPD and LOPD (Table 2 / Table S3 footnote a)
    lrate <- -0.97
    label("Mean of the log random rate on the uHex4 fractional decrease (meanr)") # Table 2 'Log rate (meanr)' = -0.97
    # Table 2 'Log error of the log rate logsigma_r' = 0.34; variance = exp(0.34)^2 = exp(0.68)
    etalrate ~ 1.97388
    addSd <- 3.32012
    label("Additive residual SD of the GMFM-88 change (percentage points)") # Table 2 'Log residual error logsigma_e' = 1.20; exp(1.20)

    # IOPD (Table 2)
    int_iopd <- 2.24
    label("IOPD intercept of the GMFM-88 change (percentage points)") # Table 2 'Intercept' = 2.24
    e_score_gmfm88_int_iopd <- -0.02
    label("IOPD effect of baseline GMFM-88 minus median (points per point)") # Table 2 'GMFM-88 total percent score at baseline' = -0.02
    e_uhex4_base_int_iopd <- -2.30
    label("IOPD effect of log baseline uHex4 minus log median (points)") # Table 2 'uHex4 (mmol/mol) at baseline (log)' = -2.30
    e_age_diagnosis_int_iopd <- 1.55
    label("IOPD effect of log age at diagnosis minus log median (points)") # Table 2 'Age at diagnosis in years (log)' = 1.55
    e_t_prior_alglucosidase_alfa_int_iopd <- -0.39
    label("IOPD effect of prior alglucosidase alfa duration minus median (points per year)") # Table 2 'Duration of previous ALGLU experience in years' = -0.39
    e_sexf_int_iopd <- 2.56
    label("IOPD effect of female sex (points)") # Table 2 'Gender (female vs. male)' = 2.56
    lmax_iopd <- 0.90
    label("IOPD log asymptote of the uHex4-driven GMFM-88 improvement (max0)") # Table 2 'Max0 (log)' = 0.90
    e_t_prior_alglucosidase_alfa_max_iopd <- -0.55
    label("IOPD effect of prior alglucosidase alfa duration on log asymptote (per year)") # Table 2 'Max [Duration of previous ALGLU experience in years]' = -0.55

    # LOPD (Table S3)
    int_lopd <- -0.78
    label("LOPD intercept of the GMFM-88 change (percentage points)") # Table S3 'Intercept' = -0.78
    e_score_gmfm88_int_lopd <- -0.00
    label("LOPD effect of baseline GMFM-88 minus median (points per point)") # Table S3 'GMFM total percent score at baseline' = -0.00
    e_uhex4_base_int_lopd <- 0.69
    label("LOPD effect of log baseline uHex4 minus log median (points)") # Table S3 'uHex4 (mmol/mol) at baseline (log)' = 0.69
    e_age_diagnosis_int_lopd <- -3.37
    label("LOPD effect of log age at diagnosis minus log median (points)") # Table S3 'Age at diagnosis in years (log)' = -3.37
    e_cpk_int_lopd <- -0.97
    label("LOPD effect of log baseline creatine kinase minus log median (points)") # Table S3 'Creatine kinase (IU/L) at baseline (log)' = -0.97
    e_creat_int_lopd <- -4.23
    label("LOPD effect of log baseline creatinine minus log median (points)") # Table S3 'Creatinine (umol/L) at baseline (log)' = -4.23
    lmax_lopd <- 2.40
    label("LOPD log asymptote of the uHex4-driven GMFM-88 improvement (max0)") # Table S3 'Max0 (log)' = 2.40
  })

  model({
    # Covariate centring values: population medians annotated on Figure S2,
    # except the prior-ALGLU duration (see covariateData notes).
    rate_uhex4 <- exp(lrate + etalrate)
    if (DIS_INFANTILE_ONSET_POMPE == 1) {
      shift <- int_iopd +
        e_score_gmfm88_int_iopd * (SCORE_GMFM88 - 67.74) +
        e_uhex4_base_int_iopd * (log(UHEX4_BASE) - log(53.81)) +
        e_age_diagnosis_int_iopd * (log(AGE_DIAGNOSIS) - log(0.26)) +
        e_t_prior_alglucosidase_alfa_int_iopd * (T_PRIOR_ALGLUCOSIDASE_ALFA - 6.5) +
        e_sexf_int_iopd * SEXF
      gmfmmax <- exp(lmax_iopd + e_t_prior_alglucosidase_alfa_max_iopd * (T_PRIOR_ALGLUCOSIDASE_ALFA - 6.5))
    } else {
      shift <- int_lopd +
        e_score_gmfm88_int_lopd * (SCORE_GMFM88 - 84.67) +
        e_uhex4_base_int_lopd * (log(UHEX4_BASE) - log(8.03)) +
        e_age_diagnosis_int_lopd * (log(AGE_DIAGNOSIS) - log(46.63)) +
        e_cpk_int_lopd * (log(CPK) - log(550)) +
        e_creat_int_lopd * (log(CREAT) - log(52))
      gmfmmax <- exp(lmax_lopd)
    }
    # -%dUHex4 * 0.01 = -UHEX4_REL: positive when uHex4 falls below baseline
    gmfm88cfb <- shift + gmfmmax * (1 - exp(-rate_uhex4 * (-UHEX4_REL)))
    gmfm88cfb ~ add(addSd)
  })
}
