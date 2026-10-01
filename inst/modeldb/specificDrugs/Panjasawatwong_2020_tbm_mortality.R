Panjasawatwong_2020_tbm_mortality <- function() {
  description <- "Weibull time-to-event model for death in Vietnamese children with tuberculous meningitis (TBM) treated with first-line antituberculosis drugs (Panjasawatwong 2020). The hazard h(t) = lambda * alpha * (lambda * t)^(alpha - 1) has one shape factor alpha (0.391, a hazard that falls with time) and a scale factor lambda estimated separately for each baseline TBM severity grade (I, II, III); grade was the only covariate retained, and first-day plasma and CSF drug exposures were tested but not retained. No drug input: the model gives the cumulative hazard and survival probability from the start of treatment as algebraic outputs. Companion PK models: Panjasawatwong_2020_isoniazid, _rifampicin, _pyrazinamide and _ethambutol."
  reference <- paste(
    "Panjasawatwong N, Wattanakul T, Hoglund RM, Bang ND, Pouplin T,",
    "Nosoongnoen W, Ngo VN, Day JN, Tarning J. (2020).",
    "Population pharmacokinetic properties of antituberculosis drugs in",
    "Vietnamese children with tuberculous meningitis.",
    "Antimicrob Agents Chemother 65(1):e00487-20.",
    "doi:10.1128/AAC.00487-20.",
    "Time-to-death model: supplemental material Table S3.",
    sep = " "
  )
  vignette <- "Panjasawatwong_2020_antituberculosis_tbm"
  units <- list(
    time = "h",
    dosing = "n/a (no drug-dosing events; the hazard depends only on baseline TBM grade)",
    concentration = "probability (the output `sur` is the probability of surviving to time t, not a drug concentration)"
  )

  covariateData <- list(
    TBM_GRADE_II = list(
      description = "Baseline TBM severity grade II indicator (1 = grade II, 0 = grade I or III)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-fixed, graded at enrolment by the Blantyre coma score (children < 5 years) or the GCS-based modified British MRC criterion (>= 5 years), Table 1 footnote b. Selects the grade II Weibull scale factor. Mutually exclusive with TBM_GRADE_III; grade I is the reference (both indicators 0).",
      source_name = "TBM severity"
    ),
    TBM_GRADE_III = list(
      description = "Baseline TBM severity grade III indicator (1 = grade III, 0 = grade I or II)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Time-fixed, graded at enrolment (Table 1 footnote b). Selects the grade III Weibull scale factor. Mutually exclusive with TBM_GRADE_II.",
      source_name = "TBM severity"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    n_observations = "15 deaths and 85 subjects censored at loss to follow-up or at the end of the study (240 days); 1 death in grade I, 4 in grade II and 10 in grade III (Table 1; Discussion).",
    age_range = "2 months to 15 years (0.167-15.0 years)",
    age_median = "3.0 years",
    weight_range = "4.0-43 kg",
    weight_median = "10.9 kg",
    sex_female_pct = 44,
    race_ethnicity = "Vietnamese",
    disease_state = "Suspected tuberculous meningitis (TBM); baseline severity grade I 58%, II 24%, III 18%; HIV positive 4%, negative 92%, unknown 4%.",
    dose_range = "Isoniazid 5, rifampicin 10, pyrazinamide 25 and ethambutol 15 mg/kg orally once daily (WHO 2006 paediatric regimen), with streptomycin for the first 2 months and adjunctive dexamethasone; treatment for 8 months.",
    regions = "Vietnam (Pham Ngoc Thach Hospital, Ho Chi Minh City), October 2009 to March 2011",
    notes = "Demographics from Panjasawatwong 2020 Table 1. Covariates screened on the baseline hazard: age, body weight, WAZ, HAZ, baseline TBM severity, HIV status, CRP, CSF protein, CSF lactate, CSF glucose and CSF/blood glucose ratio; first-day Cmax and AUC0-24 of each drug in plasma and CSF were also tested (Methods 'Population pharmacodynamic analysis'). Only TBM severity was retained (dOFV -19.0, Results). Most deaths occurred in the first week after enrolment."
  )

  ini({
    # Panjasawatwong 2020 supplemental Table S3 'Final parameter estimates
    # of the time-to-event (death) model'. Table S3 prints the Weibull
    # scale factor lambda itself for each grade, under 'Baseline (hr-1)',
    # with the equation h(t) = lambda*alpha*(lambda*t)^(alpha-1) and
    # lambda = (10^theta)(TBM Severity); the 10^theta form is how the
    # authors estimated lambda in NONMEM, and the printed values are
    # already on the lambda scale.
    llam_haz <- log(5.74e-9); label("Weibull scale factor lambda for TBM severity grade I, the reference grade (1/h)") # Table S3 'Baseline (hr-1) Grade I' = 5.74x10^-9 (RSE 73.8%; bootstrap 95% CI 6.38x10^-17 to 1.34x10^-8)
    llam_haz_tbm_grade_ii <- log(2.35e-6); label("Weibull scale factor lambda for TBM severity grade II (1/h)") # Table S3 'Baseline (hr-1) Grade II' = 2.35x10^-6 (RSE 87.0%; bootstrap 95% CI 1.73x10^-10 to 2.87x10^-6)
    llam_haz_tbm_grade_iii <- log(7.94e-5); label("Weibull scale factor lambda for TBM severity grade III (1/h)") # Table S3 'Baseline (hr-1) Grade III' = 7.94x10^-5 (RSE 11.7%; bootstrap 95% CI 5.70x10^-6 to 4.84x10^-4)
    lalfa_haz <- log(0.391); label("Weibull shape factor alpha (unitless)") # Table S3 'Slope' = 0.391 (RSE 15.3%; bootstrap 95% CI 0.304 to 0.545); the Table S3 legend calls alpha the shape factor
    # No interindividual variability and no residual error: Table S3
    # reports only the four fixed effects. A time-to-event likelihood is
    # carried by the event / censoring density, not by an observation-error
    # model; forward simulation uses the `sur` output.
  })
  model({
    # Scale factor for the subject's baseline grade (grade I when both
    # indicators are 0; Table S3 'lambda = (10^theta)(TBM Severity)').
    tbmGradeI <- 1 - TBM_GRADE_II - TBM_GRADE_III
    lam <- exp(llam_haz) *
      tbmGradeI +
      exp(llam_haz_tbm_grade_ii) * TBM_GRADE_II +
      exp(llam_haz_tbm_grade_iii) * TBM_GRADE_III
    alfa <- exp(lalfa_haz)

    # Weibull hazard (Table S3 equation). Because alpha < 1 the hazard is
    # infinite at t = 0 and falls with time, so the cumulative hazard is
    # written in its closed form (lambda * t)^alpha rather than integrated
    # as an ODE state from t = 0.
    hazard <- lam * alfa * (lam * t)^(alfa - 1)
    cumhaz <- (lam * t)^alfa
    sur <- exp(-cumhaz)
  })
}
