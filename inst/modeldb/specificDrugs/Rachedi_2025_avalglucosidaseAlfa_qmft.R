Rachedi_2025_avalglucosidaseAlfa_qmft <- function() {
  description <- paste(
    "Nonlinear random-coefficient (NLRC) model of the change from baseline in",
    "the Quick Motor Function Test (QMFT) total raw score in infantile-onset",
    "(IOPD, Mini-COMET) and late-onset (LOPD, COMET) Pompe disease treated",
    "with enzyme replacement therapy (avalglucosidase alfa or alglucosidase",
    "alfa) (Rachedi 2025). No ODE and no dose: the treatment effect enters",
    "only through the observed or simulated percent change from baseline in",
    "urinary hexose tetrasaccharide (uHex4), supplied as a time-varying",
    "covariate, via an exponential asymptote",
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
    concentration = "(observation qmftcfb is the QMFT total raw score change from baseline, score points on the 0-64 scale)"
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
    URINARY_HEXOSE_TETRASACCHARIDE_REL = list(
      description = "Relative change from baseline in urinary hexose tetrasaccharide (uHex4), as a fraction (e.g., -0.50 = 50% decrease).",
      units = "fraction",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying per visit. The paper's %dUHex4 divided by 100; the model term -%dUHex4 * 0.01 is therefore -URINARY_HEXOSE_TETRASACCHARIDE_REL. Supplied from observed uHex4 or from the Kaddi 2024 QSP model.",
      source_name = "%dUHex4"
    ),
    URINARY_HEXOSE_TETRASACCHARIDE_BASE = list(
      description = "Baseline urinary hexose tetrasaccharide (uHex4).",
      units = "mmol/mol creatinine",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Log-transformed and centred at the population median: 53.81 (IOPD) / 8.03 (LOPD) mmol/mol, Figure S2.",
      source_name = "uHex4 (mmol/mol) at baseline"
    ),
    SCORE_QMFT = list(
      description = "Baseline QMFT total raw score (16 items scored 0-4; 0-64).",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Linear, centred at the population median: 24 (IOPD) / 41 (LOPD), Figure S2.",
      source_name = "QMFT score at baseline"
    ),
    AGE_DIAGNOSIS = list(
      description = "Age at diagnosis of Pompe disease.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "IOPD only; time-fixed. Log-transformed and centred at the IOPD median 0.26 years (Figure S2).",
      source_name = "Age at diagnosis in years"
    ),
    T_PRIOR_ALGLUCOSIDASE_ALFA = list(
      description = "Duration of alglucosidase alfa treatment before study entry.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed; IOPD only (LOPD COMET patients were treatment-naive; supply 0). Linear, centred at 7 years, a rounded value between the per-arm medians 2.8-8.3 years of Figure S2 (the pooled IOPD median is not printed).",
      source_name = "Duration of previous ALGLU experience in years"
    ),
    SEXF = list(
      description = "Sex: 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "IOPD only: additive effect on the change score (Table 3 'Gender (female vs. male)').",
      source_name = "gender"
    ),
    AGE = list(
      description = "Age at first study infusion.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "LOPD only; time-fixed. Log-transformed and centred at the LOPD median 47.84 years (Figure S2).",
      source_name = "Age at first alglucosidase alfa or avalglucosidase alfa infusion in years"
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "LOPD only; time-fixed. Log-transformed and centred at the LOPD median 26.3 kg/m^2 (Figure S2).",
      source_name = "BMI at baseline"
    ),
    CPK = list(
      description = "Baseline serum creatine kinase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "LOPD only; time-fixed baseline. Log-transformed and centred at the LOPD median 550 IU/L (Figure S2).",
      source_name = "Creatine kinase (IU/L) at baseline"
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
    # Shared between IOPD and LOPD (Table 3 / Table S4 footnote a)
    lrate <- -1.19
    label("Mean of the log random rate on the uHex4 fractional decrease (meanr)") # Table 3 'Log rate (meanr)' = -1.19
    # Table 3 'Log error of the log rate logsigma_r' = 0.52; variance = exp(0.52)^2 = exp(1.04)
    etalrate ~ 2.82922
    addSd <- 2.77319
    label("Additive residual SD of the QMFT change (score points)") # Table 3 'Log residual error logsigma_e' = 1.02; exp(1.02)

    # IOPD (Table 3)
    int_iopd <- 1.51
    label("IOPD intercept of the QMFT change (score points)") # Table 3 'Intercept' = 1.51
    e_score_qmft_int_iopd <- -0.03
    label("IOPD effect of baseline QMFT minus median (points per point)") # Table 3 'QMFT score at baseline' = -0.03
    e_urinary_hexose_tetrasaccharide_base_int_iopd <- -3.34
    label("IOPD effect of log baseline uHex4 minus log median (points)") # Table 3 'uHex4 (mmol/mol) at baseline (log)' = -3.34
    e_age_diagnosis_int_iopd <- 1.66
    label("IOPD effect of log age at diagnosis minus log median (points)") # Table 3 'Age at diagnosis in years (log)' = 1.66
    e_t_prior_alglucosidase_alfa_int_iopd <- -0.45
    label("IOPD effect of prior alglucosidase alfa duration minus median (points per year)") # Table 3 'Duration of previous ALGLU experience in years' = -0.45
    e_sexf_int_iopd <- 2.81
    label("IOPD effect of female sex (points)") # Table 3 'Gender (female vs. male)' = 2.81
    lmax_iopd <- 1.69
    label("IOPD log asymptote of the uHex4-driven QMFT improvement (max0)") # Table 3 'Max0 (log)' = 1.69

    # LOPD (Table S4)
    int_lopd <- 0.06
    label("LOPD intercept of the QMFT change (score points)") # Table S4 'Intercept' = 0.06
    e_score_qmft_int_lopd <- -0.10
    label("LOPD effect of baseline QMFT minus median (points per point)") # Table S4 'QMFT score at baseline' = -0.10
    e_urinary_hexose_tetrasaccharide_base_int_lopd <- 0.35
    label("LOPD effect of log baseline uHex4 minus log median (points)") # Table S4 'uHex4 (mmol/mol) at baseline (log)' = 0.35
    e_age_int_lopd <- -3.78
    label("LOPD effect of log age at first infusion minus log median (points)") # Table S4 'Age at first ... infusion in years (log)' = -3.78
    e_bmi_int_lopd <- -2.18
    label("LOPD effect of log baseline BMI minus log median (points)") # Table S4 'BMI at baseline (log)' = -2.18
    e_cpk_int_lopd <- -0.79
    label("LOPD effect of log baseline creatine kinase minus log median (points)") # Table S4 'Creatine kinase (IU/L) at baseline (log)' = -0.79
    lmax_lopd <- 2.35
    label("LOPD log asymptote of the uHex4-driven QMFT improvement (max0)") # Table S4 'Max0' = 2.35 (log scale per Methods)
  })

  model({
    # Covariate centring values: population medians annotated on Figure S2,
    # except the prior-ALGLU duration (see covariateData notes).
    rate_uhex4 <- exp(lrate + etalrate)
    if (DIS_INFANTILE_ONSET_POMPE == 1) {
      shift <- int_iopd +
        e_score_qmft_int_iopd * (SCORE_QMFT - 24) +
        e_urinary_hexose_tetrasaccharide_base_int_iopd * (log(URINARY_HEXOSE_TETRASACCHARIDE_BASE) - log(53.81)) +
        e_age_diagnosis_int_iopd * (log(AGE_DIAGNOSIS) - log(0.26)) +
        e_t_prior_alglucosidase_alfa_int_iopd * (T_PRIOR_ALGLUCOSIDASE_ALFA - 7) +
        e_sexf_int_iopd * SEXF
      qmftmax <- exp(lmax_iopd)
    } else {
      shift <- int_lopd +
        e_score_qmft_int_lopd * (SCORE_QMFT - 41) +
        e_urinary_hexose_tetrasaccharide_base_int_lopd * (log(URINARY_HEXOSE_TETRASACCHARIDE_BASE) - log(8.03)) +
        e_age_int_lopd * (log(AGE) - log(47.84)) +
        e_bmi_int_lopd * (log(BMI) - log(26.3)) +
        e_cpk_int_lopd * (log(CPK) - log(550))
      qmftmax <- exp(lmax_lopd)
    }
    # -%dUHex4 * 0.01 = -URINARY_HEXOSE_TETRASACCHARIDE_REL: positive when uHex4 falls below baseline.
    # The printed LOPD equation omits the leading '1 -'; the base-model form
    # stated in Methods is used for both populations.
    qmftcfb <- shift + qmftmax * (1 - exp(-rate_uhex4 * (-URINARY_HEXOSE_TETRASACCHARIDE_REL)))
    qmftcfb ~ add(addSd)
  })
}
