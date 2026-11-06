Koele_2025_bedaquiline_QTcTBT <- function() {
  description <- "Population pharmacodynamic step-function model of the TB-specific time-varying heart-rate-corrected QT interval (QTcTBT, ms) in adults with drug-susceptible pulmonary tuberculosis treated with bedaquiline, delamanid and moxifloxacin with or without sutezolid or delpazolid (SUDOCU and PanACEA DECODE phase IIb trials). QTcTBT is the individual baseline until day 10 of treatment and baseline plus a constant increase afterwards. Female sex and age shift the baseline; the day-14 AUC0-24 of the bedaquiline metabolite M2 shifts the on-treatment increase. Log-normal IIV on the baseline, normal proportional IIV on the increase (correlated), and log-normal IIV on each of two additive residual-error components (a per-timepoint error and a per-ECG-replicate error). No PK: the M2 exposure is a per-subject covariate."
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
    concentration = "(observation QTcTBT is in ms; the exposure covariate AUC_BEDAQUILINE_M2 is in mg*h/L)"
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Male is the reference. Koele 2025 Text S1 codes SEX 0 = male, 1 = female, which is SEXF directly. Adds e_sexf_rbase = 13.5 ms to the baseline QTcTBT.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age at enrolment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Enters the baseline linearly as e_age_rbase * (AGE - 34), centred at the 34-year median (Koele 2025 Text S1 'AGEeff = (THETA(5)*(AGE - 34))'; Figure S3 reference subject 'male, 34 years').",
      source_name = "AGE"
    ),
    AUC_BEDAQUILINE_M2 = list(
      description = "AUC0-24 of the bedaquiline N-desmethyl metabolite M2 at day 14 of treatment",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject, time-fixed. Koele 2025 derived it by non-compartmental analysis of the day-14 PK sampling (Methods, Covariate analysis; medians 10.1-12.4 mg*h/L by arm in Tables S5 and S6). Enters the on-treatment increase linearly as e_auc_bedaquiline_m2_emax * (AUC_BEDAQUILINE_M2 - 10.54), centred at the 10.54 mg*h/L median. Text S1 imputes a missing value to the median 10.54, which zeroes the term; do the same for subjects without an M2 AUC.",
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
    notes = "Koele 2025 Table 1: SUDOCU n = 73 (54 male), DECODE n = 76 (60 male); female share (19 + 16) / 149. 2373 ECG observations (1207 SUDOCU, 1166 DECODE), recorded weekly between 8 and 10 AM before the daily dose after 10 min supine rest, single readings except triplicates when a tracing showed > 50 ms change or > 480 ms. Baseline is the mean of three pre-treatment ECGs. Machine-read QT corrected by the TB-specific time-varying formula QTcTBT = QT / RR^(0.4081 - 0.0781 * (1 - exp(-log(2) * TAST / 7.74))), TAST = weeks after the start of treatment. NONMEM 7.5, FOCE-I; precision by SIR."
  )

  ini({
    # Koele 2025 Text S1 ($PRED) and Table 2. The $THETA / $OMEGA / $SIGMA
    # values in Text S1 reproduce every Table 2 estimate, so they are the
    # final estimates.
    lrbase <- log(400)
    label("Typical baseline QTcTBT for a 34-year-old man (ms)") # Table 2 'Baseline QTcTBT (ms)' = 400 (RSE 0.4%); Text S1 THETA(1)
    emax <- 13.2
    label("Typical QTcTBT increase on treatment from day 10, at the median M2 AUC (ms)") # Table 2 'QTcTBT prolongation (ms)' = 13.2 (RSE 8.8%); Text S1 THETA(2)
    e_sexf_rbase <- 13.5
    label("Additive effect of female sex on the baseline QTcTBT (ms)") # Table 2 'Effect of female sex (ms)' = 13.5 (RSE 18%); Text S1 THETA(3)
    e_age_rbase <- 0.367
    label("Additive effect of age on the baseline QTcTBT, per year from 34 (ms/year)") # Table 2 'Effect of age (ms/year)' = 0.367 (RSE 26%); Text S1 THETA(5)
    e_auc_bedaquiline_m2_emax <- 0.783
    label("Additive effect of M2 AUC0-24 on the QTcTBT increase, per mg*h/L from 10.54 (ms per mg*h/L)") # Table 2 'Effect of M2 exposure (ms/mg/L*h)' = 0.783 (RSE 30%); Text S1 THETA(4)

    # Text S1 $OMEGA BLOCK(2): variance 0.00162 (baseline, exponential) and
    # 0.481 (increase, TV * (1 + ETA)), covariance -0.0136. Table 2 prints
    # them as sqrt(exp(omega^2) - 1): 4.02% and 78.6% CV, correlation
    # -0.0136 / sqrt(0.00162 * 0.481) = -0.487.
    etalrbase + etaemax ~ c(0.00162, -0.0136, 0.481)
    # Text S1 $OMEGA 0.225 and 0.155, exponential on the two residual SDs;
    # Table 2 'IIV additive error (CV%)' 50.2 and 'IIV additive replicate error (CV%)' 40.9.
    eta_re_add ~ 0.225
    eta_re_repl ~ 0.155

    addSd <- 8.775
    label("Additive residual error SD, per ECG timepoint (ms)") # Text S1 $SIGMA 77 -> sqrt(77) = 8.775; Table 2 'Additive error (ms)' = 8.77 (RSE 9.5%)
    addSdRepl <- 4.919
    label("Additive residual error SD, per ECG replicate (ms)") # Text S1 $SIGMA BLOCK(1) 24.2 -> sqrt(24.2) = 4.919; Table 2 'Replicate additive error (ms)' = 4.92 (RSE 13%)
  })

  model({
    # Covariate reference values from Text S1: median age 34 years and
    # median M2 AUC0-24 10.54 mg*h/L.
    age_ref <- 34
    auc_m2_ref <- 10.54

    rbase <- exp(lrbase + etalrbase) + e_age_rbase * (AGE - age_ref) + e_sexf_rbase * SEXF
    qtc_inc <- emax * (1 + etaemax) + e_auc_bedaquiline_m2_emax * (AUC_BEDAQUILINE_M2 - auc_m2_ref)

    # Step function: Text S1 'IF(TIME.LE.10) IPRED = Intercept' and
    # 'IF(TIME.GT.10) IPRED = Intercept + QTincrease', TIME in days after the
    # start of treatment.
    QTcTBT <- rbase + qtc_inc * (t > 10)

    # Text S1: Y = IPRED + ERR(1)*EXP(ETA(3)) + repERR*EXP(ETA(4)). The two
    # epsilons are independent, so a single reading has this marginal SD.
    # NONMEM shares ERR(1) across the triplicate readings of one timepoint
    # (L2 data item); that within-timepoint correlation is not represented.
    sd_qtc <- sqrt((addSd * exp(eta_re_add))^2 + (addSdRepl * exp(eta_re_repl))^2)
    QTcTBT ~ add(sd_qtc)
  })
}
