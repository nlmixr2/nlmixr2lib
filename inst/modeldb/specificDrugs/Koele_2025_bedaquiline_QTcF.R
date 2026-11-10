Koele_2025_bedaquiline_QTcF <- function() {
  description <- "Alternative population pharmacodynamic step-function model of the Fridericia-corrected QT interval (QTcF, ms) in adults with drug-susceptible pulmonary tuberculosis treated with bedaquiline, delamanid and moxifloxacin with or without sutezolid or delpazolid (SUDOCU and PanACEA DECODE phase IIb trials). Same structure as Koele_2025_bedaquiline_QTcTBT (the primary model of the paper) refitted to QTcF: QTcF is the individual baseline until day 10 of treatment and baseline plus a constant increase afterwards. Female sex and age shift the baseline; the day-14 AUC0-24 of the bedaquiline metabolite M2 shifts the on-treatment increase. Log-normal IIV on the baseline, normal proportional IIV on the increase (correlated), and log-normal IIV on each of two additive residual-error components (a per-timepoint error and a per-ECG-replicate error). No PK: the M2 exposure is a per-subject covariate."
  reference <- paste(
    "Koele SE, Stoycheva K, Mtweve C, Manyama C, Mpagama S, Mhimbira F, et al.; PanACEA consortium.",
    "Cardiac safety of bedaquiline, delamanid and moxifloxacin co-administered with or without",
    "varying doses of sutezolid or delpazolid for the treatment of drug-susceptible TB.",
    "J Antimicrob Chemother. 2025;80(8):2305-2313. doi:10.1093/jac/dkaf210"
  )
  vignette <- "Koele_2025_bedaquiline_QTc"
  paper_specific_etas <- c("eta_re_add", "eta_re_repl")
  paper_specific_residual_sds <- c("addSdRepl")
  units <- list(
    time = "day",
    dosing = "(none; PD-only model driven by time after the start of treatment and a per-subject exposure covariate)",
    concentration = "(observation QTcF is in ms; the exposure covariate AUC_BEDAQUILINE_M2 is in mg*h/L)"
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Male is the reference. Koele 2025 Text S1 codes SEX 0 = male, 1 = female, which is SEXF directly. Adds e_sexf_rbase = 11.8 ms to the baseline QTcF.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age at enrolment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Enters the baseline linearly as e_age_rbase * (AGE - 34), centred at the 34-year median (Table S4 footnote a 'parameterized as ms change per unit of deviation from median value'; the median is 34 in Text S1 of the primary QTcTBT model).",
      source_name = "AGE"
    ),
    AUC_BEDAQUILINE_M2 = list(
      description = "AUC0-24 of the bedaquiline N-desmethyl metabolite M2 at day 14 of treatment",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject, time-fixed, from non-compartmental analysis of the day-14 PK sampling (Koele 2025 Methods; Tables S5 and S6). Enters the on-treatment increase linearly as e_auc_bedaquiline_m2_emax * (AUC_BEDAQUILINE_M2 - 10.54), centred at the 10.54 mg*h/L median of Text S1. Impute a missing value to 10.54, as Text S1 does.",
      source_name = "AUCM2"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 149L,
    n_studies = 2L,
    n_observations = 2373L,
    age_range = "20-58 years",
    age_median = "33 years (SUDOCU), 34 years (DECODE)",
    weight_range = "40.0-84.0 kg",
    weight_median = "53.0 kg in both studies",
    sex_female_pct = 23.5,
    race_ethnicity = "Not tabulated; participants recruited in Tanzania and South Africa",
    disease_state = "Newly diagnosed smear-positive drug-susceptible pulmonary tuberculosis; 13 of 149 (8.7%) living with HIV. Pre-existing cardiac disease, QTcF > 450 ms and concomitant QT-prolonging medication were exclusion criteria.",
    dose_range = "Bedaquiline 400 mg daily for 2 weeks then 200 mg thrice weekly, delamanid 100 mg twice daily and moxifloxacin 400 mg daily, plus either sutezolid (0, 600 mg QD, 1200 mg QD, 600 mg BID or 800 mg BID; SUDOCU, 12 weeks) or delpazolid (0, 400 mg QD, 800 mg QD, 1200 mg QD or 800 mg BID; DECODE, 16 weeks)",
    regions = "Tanzania, South Africa",
    notes = "Same cohort and ECG data as Koele_2025_bedaquiline_QTcTBT (Koele 2025 Table 1; 2373 ECG observations). QTcF = QT / RR^(1/3). The paper reports this fit as a sensitivity analysis (Table S4) to show that Fridericia under-corrects the pre-treatment tachycardia of TB: the QTcF baseline (388 ms) is lower and the on-treatment increase (21.8 ms) larger than with QTcTBT (400 ms and 13.2 ms)."
  )

  ini({
    # Koele 2025 Table S4. The paper prints no code for this model; it states
    # the structure is that of the primary QTcTBT model (Text S1).
    lrbase <- log(388)
    label("Typical baseline QTcF for a 34-year-old man (ms)") # Table S4 'Baseline QTcF [ms]' = 388 (RSE 0.5%)
    emax <- 21.8
    label("Typical QTcF increase on treatment from day 10, at the median M2 AUC (ms)") # Table S4 'QTcF prolongation [ms]' = 21.8 (RSE 5.5%)
    e_sexf_rbase <- 11.8
    label("Additive effect of female sex on the baseline QTcF (ms)") # Table S4 'Effect of female sex [ms]' = 11.8 (RSE 25%)
    e_age_rbase <- 0.391
    label("Additive effect of age on the baseline QTcF, per year from 34 (ms/year)") # Table S4 'Effect of age [ms/year]' = 0.391 (RSE 39%)
    e_auc_bedaquiline_m2_emax <- 0.850
    label("Additive effect of M2 AUC0-24 on the QTcF increase, per mg*h/L from 10.54 (ms per mg*h/L)") # Table S4 'Effect of M2 exposure [ms/mg/L*h]' = 0.850 (RSE 28%)

    # Table S4 prints CV% = 100 * sqrt(exp(omega^2) - 1), the convention that
    # reproduces Table 2 from the Text S1 omegas of the QTcTBT model, so
    # omega^2 = log(1 + CV^2):
    #   baseline 4.30% -> 0.0018473; increase 43.5% -> 0.17330;
    #   covariance -0.481 * sqrt(0.0018473 * 0.17330) = -0.0086063.
    etalrbase + etaemax ~ c(0.0018473, -0.0086063, 0.17330) # Table S4 'IIV baseline QTcF (CV%)' 4.30, 'IIV QTcF prolongation (CV%)' 43.5, 'Correlation baseline QTcF and QTcF prolongation (%)' -48.1
    eta_re_add ~ 0.19728 # Table S4 'IIV additive error (CV%)' = 46.7 -> log(1 + 0.467^2)
    eta_re_repl ~ 0.18067 # Table S4 'IIV additive replicate error (CV%)' = 44.5 -> log(1 + 0.445^2)

    addSd <- 9.41
    label("Additive residual error SD, per ECG timepoint (ms)") # Table S4 'Additive error [ms]' = 9.41 (RSE 8.9%)
    addSdRepl <- 4.48
    label("Additive residual error SD, per ECG replicate (ms)") # Table S4 'Replicate additive error [ms]' = 4.48 (RSE 12%)
  })

  model({
    # Covariate reference values from Text S1: median age 34 years and
    # median M2 AUC0-24 10.54 mg*h/L.
    age_ref <- 34
    auc_m2_ref <- 10.54

    rbase <- exp(lrbase + etalrbase) + e_age_rbase * (AGE - age_ref) + e_sexf_rbase * SEXF
    qtc_inc <- emax * (1 + etaemax) + e_auc_bedaquiline_m2_emax * (AUC_BEDAQUILINE_M2 - auc_m2_ref)

    # Step function at day 10 of treatment, as in Text S1 of the QTcTBT model.
    QTcF <- rbase + qtc_inc * (t > 10)

    # Marginal SD of one reading: independent per-timepoint and per-replicate
    # errors, each scaled by its own exponential IIV. The within-timepoint
    # correlation NONMEM builds from the L2 data item is not represented.
    sd_qtc <- sqrt((addSd * exp(eta_re_add))^2 + (addSdRepl * exp(eta_re_repl))^2)
    QTcF ~ add(sd_qtc)
  })
}
