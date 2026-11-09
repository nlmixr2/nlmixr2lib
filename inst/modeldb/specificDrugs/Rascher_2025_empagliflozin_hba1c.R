Rascher_2025_empagliflozin_hba1c <- function() {
  description <- paste(
    "Exposure-response (indirect response) model for HbA1c in paediatric",
    "patients aged 10-17 years with type 2 diabetes mellitus treated with",
    "empagliflozin or placebo (DINAMO trial). Steady-state empagliflozin AUC",
    "(AUC_EMPA, from the companion Rascher_2025_empagliflozin popPK) inhibits",
    "HbA1c synthesis kin through an Imax function with AUC50 fixed at 703",
    "nmol*h/L; kin rises by a zero-order disease-progression rate in patients",
    "on background insulin. Covariates: insulin co-therapy on baseline HbA1c",
    "and on progression; eGFR (truncated at 120) and observed baseline HbA1c",
    "power effects on Imax. Correlated IIV on baseline and progression."
  )
  reference <- paste(
    "Rascher J, Cheng S, Johnston C, Haertter S, Marquard J, Tartakovsky I,",
    "Laffel L, Wojtyniak J-G.",
    "Pharmacokinetics and pharmacodynamics of empagliflozin in paediatric",
    "patients aged 10-17 years with type 2 diabetes mellitus.",
    "Br J Clin Pharmacol. 2025;91(8):2390-2400. doi:10.1002/bcp.70096."
  )
  vignette <- "Rascher_2025_empagliflozin"
  units <- list(
    time = "h",
    dosing = "none (empagliflozin exposure enters as the AUC_EMPA covariate)",
    concentration = "HbA1c in % (NGSP)"
  )

  compartmentData <- list(
    Hba1c = list(analyte = "HbA1c", units = "%", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    AUC_EMPA = list(
      description = "Individual steady-state empagliflozin AUC over the 24-h dosing interval",
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Derived by the authors as administered daily dose / individual CL/F",
        "(empirical Bayes CL/F from Rascher_2025_empagliflozin) and appended to",
        "the exposure-response data set (Methods 'PopPK', last paragraph).",
        "0 for placebo. Time-varying when the dose changes (e.g. 10 mg weeks",
        "1-14 then 25 mg weeks 15-26)."
      ),
      source_name = "AUCss"
    ),
    HBA1C = list(
      description = "Observed baseline HbA1c (per-subject, time-fixed)",
      units = "% (NGSP)",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters Imax as (HBA1C/7.9)^theta (Supplement Imax equations; reference 7.9 %, 'typical population reference value'). Distinct from the modelled Hba1c state, whose baseline is the estimated BASE.",
      source_name = "HbA1c (baseline)"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalised)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying eGFR (Zappitelli equation in DINAMO, Table S2 footnote), truncated at 120 mL/min/1.73 m^2 inside the model and entering Imax as (eGFR/90)^theta (Supplement Imax equations).",
      source_name = "eGFR"
    ),
    CONMED_INSULIN = list(
      description = "Insulin co-therapy at baseline (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no background insulin)",
      notes = paste(
        "Multiplies baseline HbA1c by exp(theta 6) = 1.15 and switches on the",
        "typical zero-order progression rate of kin (PROG_INS = 1; the typical",
        "progression is 0 without insulin) (Supplement PK/PD equations).",
        "49.4% of the PK/PD cohort (Results 'Data')."
      ),
      source_name = "INS (insulin co-therapy at baseline)"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_METFORMIN = list(
      description = "Concomitant metformin",
      units = "(binary)",
      type = "binary",
      notes = "90.3% of patients were on metformin at baseline, so metformin was not modelled as a covariate (Results 'Data')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 103L,
    n_studies = 1L,
    n_observations = 394L,
    studies = "DINAMO (NCT03429543), phase 3 randomised, double-blind trial in young people with type 2 diabetes.",
    age_range = "10-17 years",
    disease_state = "Type 2 diabetes mellitus (paediatric, 10-17 years)",
    dose_range = "Placebo (n = 51); empagliflozin 10 mg weeks 1-26 (n = 39) or 10 mg weeks 1-14 then 25 mg weeks 15-26 (n = 13), orally once daily (Results 'Data', Table 2).",
    co_medication = "Insulin at baseline 51 patients (49.4%); metformin at baseline 93 patients (90.3%).",
    regions = "Multinational (DINAMO).",
    notes = "Reference subject: no background insulin, eGFR 90 mL/min/1.73 m^2, baseline HbA1c 7.9 %."
  )

  ini({
    # Table 4 (DINAMO model posterior medians; log-domain parameters printed
    # back-transformed). kout is printed per day; the model runs in hours.
    lkout <- log(0.0489); label("HbA1c degradation rate constant kout (1/day; divided by 24 in model())") # Table 4 Kout = 0.0489 1/day (95% CDI 0.0346, 0.0717)
    lrbase <- log(7.35); label("Baseline HbA1c without background insulin (%)") # Table 4 BASE = 7.35 % (95% CDI 6.99, 7.72)
    ldp_slope <- log(5.64e-07); label("Zero-order disease-progression rate of kin with background insulin (%/h/h)") # Table 4 PROG = 5.64e-07 %/h/h (95% CDI 2.94e-07, 9.20e-07)
    logitimax <- log(0.101 / (1 - 0.101)); label("Logit of typical maximum fractional inhibition of kin (Imax = 0.101)") # Table 4 IMAX = exp(theta 4)/(1 + exp(theta 4)) = 10.1 % (95% CDI 7.39, 12.9)
    lauc50 <- fixed(log(703)); label("Steady-state AUC at 50% of Imax (nmol*h/L)") # Table 4 AUC50 = 703 nmol*h/L, fixed to the adult Baron 2016 value

    rbase_insulin <- 1.15; label("Multiplicative effect of insulin co-therapy on baseline HbA1c (ratio)") # Table 4 INS BASE = exp(theta 6) = 1.15
    e_crcl_imax <- 1.03; label("Power exponent of (min(eGFR, 120)/90) on Imax (unitless)") # Table 4 EGFR IMAX = 1.03 (95% CDI 0.917, 1.15)
    e_hba1c_imax <- 2.04; label("Power exponent of (baseline HbA1c/7.9) on Imax (unitless)") # Table 4 HbA1c IMAX = 2.04 (95% CDI 1.86, 2.21)

    # IIV block (Table 4). Diagonals printed as CV% = sqrt(exp(omega^2) - 1);
    # the off-diagonal 0.00112 is the printed covariance (correlation 0.095).
    # etadp_slope is ADDITIVE on the progression rate, scaled by 1/100000.
    etalrbase + etadp_slope ~ c(
      0.02559,
      0.00112, 0.005388
    ) # Table 4 IIV-BASE 16.1 %CV -> log(1 + 0.161^2) = 0.02559; IIV-PROG 7.35 -> log(1 + 0.0735^2) = 0.005388; cov 0.00112

    expSd <- 0.0632; label("Additive residual error on the log scale for HbA1c (SD)") # Table 4 Sigma 1,1 = 6.32 %CV (sqrt(sigma^2))
  })
  model({
    kout <- exp(lkout) / 24
    rbase <- exp(lrbase + etalrbase) * rbase_insulin^CONMED_INSULIN

    # Supplement: PROG_i = exp(theta 3) * PROG_INS + eta3 / 100000, with
    # PROG_INS = 1 for insulin co-therapy at baseline and 0 otherwise.
    dp_slope <- exp(ldp_slope) * CONMED_INSULIN + etadp_slope / 100000

    # Supplement Imax equations: time-varying eGFR truncated at 120.
    egfr_imax <- CRCL
    if (CRCL > 120) egfr_imax <- 120
    imax <- expit(logitimax) *
      (egfr_imax / 90)^e_crcl_imax *
      (HBA1C / 7.9)^e_hba1c_imax
    auc50 <- exp(lauc50)
    inh <- imax * AUC_EMPA / (auc50 + AUC_EMPA)

    # dkin/dt = PROG with kin(0) = BASE * kout, integrated in closed form;
    # t = 0 is the baseline (first) HbA1c observation.
    kin <- rbase * kout + dp_slope * t

    d / dt(Hba1c) <- kin * (1 - inh) - kout * Hba1c
    Hba1c(0) <- rbase
    Hba1c ~ lnorm(expSd)
  })
}
