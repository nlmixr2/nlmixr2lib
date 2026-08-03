Lindauer_2017_lacosamide_seizure <- function() {
  description <- "Repeated time-to-seizure model for adult patients newly diagnosed with focal or generalized tonic-clonic epilepsy on lacosamide (LCM) or carbamazepine controlled-release (CBZ-CR) monotherapy, from Lindauer 2017 (SP0993 trial; NCT01243177). Two Weibull sub-model hazards are exposed as separate outputs: hazard_1st for the time to the first seizure and hazard_2nd for the time to the second and subsequent seizures. Base Weibull hazard h(t) = lam * p * (lam * t)^(p - 1) with distinct scale (lam) and shape (p) for each sub-model. Baseline disease severity (NSP3M: number of seizures in the previous 3 months, categorised <2 / 2-6 reference / 7-50 / >50), daily AUC of the assigned drug (LCM or CBZ-CR; centred at 104 mg*h/L for LCM and 132 mg*h/L for CBZ), and age (LCM-only, on the first-seizure hazard) enter as log-linear covariate effects on lam. IIV on ln(lam_2nd) with SD 2.03 captures the substantial between-subject variability in the subsequent-seizure hazard. Sister dropout model for the same trial: modellib('Lindauer_2017_lacosamide_dropout')."
  reference <- paste(
    "Lindauer A, Laveille C, Stockis A.",
    "Time-to-seizure modeling of lacosamide used in monotherapy in patients",
    "with newly diagnosed epilepsy.",
    "Clin Pharmacokinet. 2017;56(11):1403-1413.",
    "doi:10.1007/s40262-017-0530-8.",
    "Open Access under CC BY-NC 4.0.",
    "Sister model file from the same paper:",
    "modellib('Lindauer_2017_lacosamide_dropout').",
    sep = " "
  )
  vignette <- "Lindauer_2017_lacosamide"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the daily-AUC covariate is a per-subject dose-step-varying summary of exposure)",
    concentration = "probability (the model outputs `sur_1st` and `sur_2nd` are survival probabilities, not drug concentrations)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    cumhaz_1st = list(analyte = "Hazard of first seizure", units = NA_character_, specimen = "not applicable", verified = FALSE),
    cumhaz_2nd = list(analyte = "Hazard of second and subsequent seizures", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age at enrolment (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the first-seizure log-hazard as e_age_1st * CONMED_LCM * (AGE - 41), i.e., the AGE effect is active only for patients on lacosamide (Lindauer 2017 Results Section 3.4: 'A significant effect of age was only identified for the LCM group ... indicating a treatment-age interaction'). Reference age 41 years is the cohort median (Lindauer 2017 Table 1).",
      source_name        = "AGE"
    ),
    NSP3M_LT2 = list(
      description        = "Baseline seizure-severity indicator: 1 = the patient reported fewer than 2 seizures in the 3 months before the trial; 0 = otherwise. One of three binary indicators derived from the paper's categorical NSP3M variable (baseline seizure count in prior 3 months split into <2, 2-6, 7-50, >50 with 2-6 as reference).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (2-6 seizures in the previous 3 months; the reference NSP3M category)",
      notes              = "Time-fixed per subject. Members of the reference category (NSP3M = 2-6) have NSP3M_LT2 = NSP3M_7_50 = NSP3M_GT50 = 0. Distribution in SP0993 (Lindauer 2017 Table 1): 25.7% of patients had <2 seizures. Enters both first-seizure and subsequent-seizure log-hazards.",
      source_name        = "NSP3M (categorical; the <2 bin)"
    ),
    NSP3M_7_50 = list(
      description        = "Baseline seizure-severity indicator: 1 = the patient reported 7 to 50 seizures in the 3 months before the trial; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (2-6 seizures in the previous 3 months; the reference NSP3M category)",
      notes              = "Time-fixed per subject. Distribution in SP0993 (Lindauer 2017 Table 1): 18.1% of patients had 7-50 seizures. Enters both first-seizure and subsequent-seizure log-hazards.",
      source_name        = "NSP3M (categorical; the 7-50 bin)"
    ),
    NSP3M_GT50 = list(
      description        = "Baseline seizure-severity indicator: 1 = the patient reported more than 50 seizures in the 3 months before the trial; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (2-6 seizures in the previous 3 months; the reference NSP3M category)",
      notes              = "Time-fixed per subject. Distribution in SP0993 (Lindauer 2017 Table 1): 4.8% of patients had >50 seizures (median 2, range 0-450 across the whole cohort). Enters both first-seizure and subsequent-seizure log-hazards.",
      source_name        = "NSP3M (categorical; the >50 bin)"
    ),
    CONMED_LCM = list(
      description        = "Treatment-arm indicator: 1 = subject assigned to the lacosamide (LCM) monotherapy arm; 0 = subject assigned to the carbamazepine controlled-release (CBZ-CR) monotherapy arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CBZ-CR arm)",
      notes              = "Time-fixed per subject in the SP0993 parallel-group monotherapy design. Enters the first-seizure log-hazard as a multiplier on the AGE effect (age effect is active only for LCM per Lindauer 2017 Section 3.4). Also used inside the vignette to decide which AUC column (AUC_LCM vs AUC_CBZ) is non-zero for the subject.",
      source_name        = "TYPE"
    ),
    AUC_LCM = list(
      description        = "Daily area under the plasma concentration-time curve of lacosamide at steady state (mg*h/L). Set to 0 for subjects randomised to CBZ-CR (CONMED_LCM = 0).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dose-step-varying: as a patient escalates from LCM 100 -> 200 -> 400 -> 600 mg/day, the AUC increases in proportion (see Lindauer 2017 Section 2.3 -- individual AUC derived from empirical-Bayes clearance from a previously published lacosamide popPK model). Enters both first-seizure and subsequent-seizure log-hazards as a centred-deviation form e_auc_lcm_1st * (AUC_LCM - 104) (analogous for the 2nd-event sub-model). Reference 104 mg*h/L is the typical daily AUC at the first target dose level for LCM (200 mg/day; Lindauer 2017 Table 3 note c). The paper reports typical daily AUC roughly 208 mg*h/L at 400 mg/day and 312 mg*h/L at 600 mg/day; document any per-simulation profile in `covariateData[[AUC_LCM]]$notes`.",
      source_name        = "AUC_LCM"
    ),
    AUC_CBZ = list(
      description        = "Daily area under the plasma concentration-time curve of carbamazepine at steady state (mg*h/L). Set to 0 for subjects randomised to LCM (CONMED_LCM = 1).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dose-step-varying: as a patient escalates from CBZ-CR 200 -> 400 -> 800 -> 1200 mg/day, the AUC increases in proportion (see Lindauer 2017 Section 2.3 -- individual AUC derived from empirical-Bayes clearance from a previously published carbamazepine popPK model). Enters both first-seizure and subsequent-seizure log-hazards as a centred-deviation form e_auc_cbz_1st * (AUC_CBZ - 132). Reference 132 mg*h/L is the typical daily AUC at the first target dose level for CBZ-CR (400 mg/day for a 70-kg patient; Lindauer 2017 Section 3.4 first paragraph).",
      source_name        = "AUC_CBZ"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 883L,
    n_studies      = 1L,
    age_range      = "16-87 years (median 40, IQR 26-55)",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = 46.3,
    race_ethnicity = NULL,
    disease_state  = "Adult patients (>=16 years) newly diagnosed with focal or generalized tonic-clonic seizures without signs of focal onset, provided they had no history or clinical or electroencephalographic findings suggestive of idiopathic generalized epilepsy (SP0993 inclusion criteria; ClinicalTrials.gov NCT01243177).",
    dose_range     = "LCM target dose levels 200, 400, or 600 mg/day BID (randomisation starting dose 100 mg/day); CBZ-CR target dose levels 400, 800, or 1200 mg/day BID (randomisation starting dose 200 mg/day). Dose escalation to the next level was triggered by seizure occurrence during the 26-week evaluation period at the current dose.",
    regions        = "multinational (SP0993)",
    biomarkers     = "Repeated time-to-seizure (seizure defined as one seizure-event day, so multiple seizures on the same day count as a single event). About 55% of patients had no seizures during the trial; 15% had one seizure; about 30% had two or more (Lindauer 2017 Section 3.1). Baseline covariate NSP3M (seizure count in the previous 3 months) had median 2 and range 0-450.",
    notes          = "Randomised 883 patients (LCM 443, CBZ-CR 440) analysed for seizures. Seizure model was originally developed on the historic N01061 dataset (comparing levetiracetam and CBZ-CR; NCT00150735) using the two-Weibull-sub-model approach of Abrantes et al. and re-estimated on SP0993 as described in Lindauer 2017 Section 2.2 'Modeling Strategy and Software'."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Lindauer 2017 Table 3 'Parameter estimates of
    # the final seizure model'. Scale parameters k1 and k2 are back-
    # transformed via note (b) 'Parameter back-transformed to normal
    # scale as exp(x)'; on the log scale they enter ini() as log(k_i).
    # Shape parameters p1 and p2 are on the linear scale in the paper --
    # log-transformed inside ini() so back-transformation gives the paper's
    # values directly. NSP3M / AUC / AGE covariate coefficients are on the
    # log-hazard scale (linear values, added to ln(lam)).
    # ----------------------------------------------------------------------

    # First-seizure Weibull baseline hazard (typical patient aged 41 years,
    # NSP3M = 2-6 reference category, typical CBZ AUC 132 mg*h/L).
    llam_1st <- log(0.000733); label("Log Weibull scale lambda for time-to-first-seizure hazard (1/day)")   # Lindauer 2017 Table 3 k1 = 0.000733 (90% CI 0.000534-0.000957)
    lp_1st   <- log(0.493);    label("Log Weibull shape p for time-to-first-seizure hazard (unitless)")     # Lindauer 2017 Table 3 p1 = 0.493 (90% CI 0.464-0.531)

    # Subsequent-seizure Weibull baseline hazard.
    llam_2nd <- log(0.00889);  label("Log Weibull scale lambda for time-to-second-and-subsequent-seizure hazard (1/day)") # Lindauer 2017 Table 3 k2 = 0.00889 (90% CI 0.00606-0.013)
    lp_2nd   <- log(0.713);    label("Log Weibull shape p for time-to-second-and-subsequent-seizure hazard (unitless)")   # Lindauer 2017 Table 3 p2 = 0.713 (90% CI 0.675-0.75)

    # NSP3M category effects on log(lam_1st) (Lindauer 2017 Table 3
    # 'NSP3M ( * ) * k1' rows; reference category 2-6 seizures).
    e_nsp3m_lt2_1st  <- -1.12; label("Log-hazard shift of NSP3M<2 (vs 2-6) on first-seizure lambda")   # Lindauer 2017 Table 3 NSP3M(<2)*k1 = -1.12 (90% CI -1.65 to -0.648)
    e_nsp3m_7_50_1st <-  1.94; label("Log-hazard shift of NSP3M 7-50 (vs 2-6) on first-seizure lambda") # Lindauer 2017 Table 3 NSP3M(7-50)*k1 = 1.94 (90% CI 1.43-2.43)
    e_nsp3m_gt50_1st <-  3.30; label("Log-hazard shift of NSP3M>50 (vs 2-6) on first-seizure lambda")   # Lindauer 2017 Table 3 NSP3M(>50)*k1 = 3.30 (90% CI 2.08-5.07)

    # NSP3M category effects on log(lam_2nd).
    e_nsp3m_lt2_2nd  <- -1.37; label("Log-hazard shift of NSP3M<2 (vs 2-6) on subsequent-seizure lambda")   # Lindauer 2017 Table 3 NSP3M(<2)*k2 = -1.37 (90% CI -2.06 to -0.659)
    e_nsp3m_7_50_2nd <-  1.36; label("Log-hazard shift of NSP3M 7-50 (vs 2-6) on subsequent-seizure lambda") # Lindauer 2017 Table 3 NSP3M(7-50)*k2 = 1.36 (90% CI 0.898-1.82)
    e_nsp3m_gt50_2nd <-  2.53; label("Log-hazard shift of NSP3M>50 (vs 2-6) on subsequent-seizure lambda")   # Lindauer 2017 Table 3 NSP3M(>50)*k2 = 2.53 (90% CI 1.98-3.1)

    # AUC covariate slopes (centred-deviation form; enter as
    # slope * (AUC - AUCref)). Paper Section 3.4 states 'AUC centered on
    # typical value at the first dose level for each drug'.
    e_auc_lcm_1st <- -0.00917; label("Slope of daily-AUC-LCM on first-seizure ln(lambda), centred at 104 mg*h/L (1/(mg*h/L))")     # Lindauer 2017 Table 3 AUC_LCM*k1 = -0.00917 (90% CI -0.0164 to -0.000334)
    e_auc_cbz_1st <- -0.00658; label("Slope of daily-AUC-CBZ on first-seizure ln(lambda), centred at 132 mg*h/L (1/(mg*h/L))")     # Lindauer 2017 Table 3 AUC_CBZ*k1 = -0.00658 (90% CI -0.015 to 0.00325)
    e_auc_lcm_2nd <- -0.00751; label("Slope of daily-AUC-LCM on subsequent-seizure ln(lambda), centred at 104 mg*h/L (1/(mg*h/L))") # Lindauer 2017 Table 3 AUC_LCM*k2 = -0.00751 (90% CI -0.00985 to -0.006)
    e_auc_cbz_2nd <- -0.0153;  label("Slope of daily-AUC-CBZ on subsequent-seizure ln(lambda), centred at 132 mg*h/L (1/(mg*h/L))") # Lindauer 2017 Table 3 AUC_CBZ*k2 = -0.0153 (90% CI -0.0183 to -0.0122)

    # AGE effect on log(lam_1st), active only for LCM patients (paper's
    # treatment-age interaction). Applied inside model() as
    # e_age_1st * CONMED_LCM * (AGE - 41).
    e_age_1st <- -0.0256; label("Slope of AGE on first-seizure ln(lambda) for LCM patients only, centred at 41 years (1/year)") # Lindauer 2017 Table 3 AGE*k1 = -0.0256 (90% CI -0.0449 to -0.0112)

    # IIV: SD of the random effect on ln(lam_2nd) only (paper: 'The
    # random-effect parameter eta is only included in the second hazard
    # equation'). Table 3 reports the SD 2.03 directly (not variance);
    # enter as var = SD^2 = 4.1209.
    etallam_2nd ~ 4.1209  # Lindauer 2017 Table 3 IIV ln(k2) SD = 2.03 (90% CI 1.85-2.2)

    # Residual error placeholder. The source NONMEM run uses the survival
    # / event-density likelihood (Laplace); no observation-error model.
    # This nlmixr2 translation is intended for forward simulation:
    # `hazard_1st`, `cumhaz_1st`, `sur_1st`, `hazard_2nd`, `cumhaz_2nd`,
    # and `sur_2nd` are exposed as derived outputs and a tiny additive
    # residual is attached to `sur_1st` so the nlmixr2 likelihood
    # machinery accepts the model.
    addSd <- 0.001; label("Placeholder additive residual error on the first-event survival-probability output `sur` (unitless); not from the source -- see vignette Assumptions")
  })

  model({
    # Back-transformed Weibull structural parameters.
    lam1_typical <- exp(llam_1st)
    p1           <- exp(lp_1st)
    lam2_typical <- exp(llam_2nd + etallam_2nd)
    p2           <- exp(lp_2nd)

    # NSP3M-category contribution to ln(lam).
    cov_nsp3m_1st <- e_nsp3m_lt2_1st * NSP3M_LT2 + e_nsp3m_7_50_1st * NSP3M_7_50 + e_nsp3m_gt50_1st * NSP3M_GT50
    cov_nsp3m_2nd <- e_nsp3m_lt2_2nd * NSP3M_LT2 + e_nsp3m_7_50_2nd * NSP3M_7_50 + e_nsp3m_gt50_2nd * NSP3M_GT50

    # Drug-specific AUC contributions (centred-deviation form). AUC_LCM is
    # 0 for CBZ-CR subjects and AUC_CBZ is 0 for LCM subjects, so only one
    # of the two centred terms is active at a time -- and the OTHER centred
    # term evaluates to slope * (0 - ref), which would spuriously modify the
    # reference-log-hazard. Gate each term by the treatment-arm indicator so
    # only the on-drug centred deviation contributes.
    cov_auc_1st <- e_auc_lcm_1st * CONMED_LCM * (AUC_LCM - 104) +
      e_auc_cbz_1st * (1 - CONMED_LCM) * (AUC_CBZ - 132)
    cov_auc_2nd <- e_auc_lcm_2nd * CONMED_LCM * (AUC_LCM - 104) +
      e_auc_cbz_2nd * (1 - CONMED_LCM) * (AUC_CBZ - 132)

    # AGE contribution (LCM-only, first-event hazard only; centred at 41 y).
    cov_age_1st <- e_age_1st * CONMED_LCM * (AGE - 41)

    # Effective log-scale Weibull scale parameters.
    lam1_eff <- lam1_typical * exp(cov_nsp3m_1st + cov_auc_1st + cov_age_1st)
    lam2_eff <- lam2_typical * exp(cov_nsp3m_2nd + cov_auc_2nd)

    # Weibull baseline hazards. The DEL = 1e-6 small-time offset keeps the
    # (lam * t)^(p - 1) term finite at t = 0 for p < 1 (both p1 = 0.493
    # and p2 = 0.713 are < 1) without affecting the integrated cumulative
    # hazard. Mirrors the Hansson_2013_sunitinib_os convention.
    del <- 1e-6
    hazard_1st <- lam1_eff * p1 * (lam1_eff * (t + del))^(p1 - 1)
    hazard_2nd <- lam2_eff * p2 * (lam2_eff * (t + del))^(p2 - 1)

    # Cumulative hazards and survivor functions. For simulation:
    #   sur(t) = P(no seizure in [0, t]) = probability of remaining
    #     seizure-free up to time t (the first-seizure Weibull survivor
    #     function; this is the observation variable);
    #   sur_2nd(t) = the analogous survivor function for the subsequent-
    #     seizure hazard, integrated from t = 0 (derived diagnostic
    #     output). In practice a subject enters the 2nd+ sub-model only
    #     after having their first seizure, so sur_2nd(t) applied
    #     conditionally on the first-seizure time T1 gives the
    #     probability of no further seizure between T1 and t.
    d/dt(cumhaz_1st) <- hazard_1st
    cumhaz_1st(0)    <- 0
    sur              <- exp(-cumhaz_1st)

    d/dt(cumhaz_2nd) <- hazard_2nd
    cumhaz_2nd(0)    <- 0
    sur_2nd          <- exp(-cumhaz_2nd)

    # Forward-simulation observation: placeholder additive residual on the
    # first-event survival-probability output so the nlmixr2 likelihood
    # machinery accepts the model.
    sur ~ add(addSd)
  })
}
