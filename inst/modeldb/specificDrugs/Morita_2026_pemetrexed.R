Morita_2026_pemetrexed <- function() {
  description <- paste(
    "Empirical piecewise-quadratic population model of the absolute",
    "neutrophil count (ANC) during the first 21-day cycle of",
    "pemetrexed-platinum (cisplatin or carboplatin) chemotherapy in Japanese",
    "adults with non-squamous non-small cell lung cancer (Morita 2026).",
    "ANC falls along one parabola from the patient's own pre-treatment",
    "count to a nadir at tnadir weeks and recovers along a second parabola",
    "of different curvature. The decline curvature depends on blood urea",
    "nitrogen, the recovery curvature on concomitant renin-angiotensin",
    "system inhibitor use, and the nadir time on whether the PBPK-predicted",
    "pemetrexed AUC0-24 exceeds 175 ug*h/mL. No ODEs and no dosing: the",
    "rxode2 time variable is the time since the cycle-1 pemetrexed dose in",
    "weeks, and the model is valid over that first cycle only.",
    sep = " "
  )
  reference <- paste(
    "Morita K, Shigetome K, Narise H, Kaneko T, Soejima N, Tanaka R, Jono H,",
    "Itoh H, Kadowaki D, Tokunaga K, Shibata A, Tanoue H, Ichikado K,",
    "Kajiwara-Morita A, Oniki K, Saruwatari J. Risk Assessment of",
    "Neutropenia Based on Exposure Information Without Plasma Concentration",
    "Measurement in Pemetrexed-Platinum-Based Chemotherapy: A Modeling",
    "Approach Using Real-World Clinical Data. Pharmaceutics.",
    "2026;18(8):946. doi:10.3390/pharmaceutics18080946.",
    sep = " "
  )
  vignette <- "Morita_2026_pemetrexed"
  units <- list(
    time = "week",
    dosing = "n/a (no dosing records; time is measured from the cycle-1 pemetrexed dose)",
    concentration = "10^3 cells/uL (absolute neutrophil count; equivalent to 10^9/L)"
  )

  covariateData <- list(
    NEUT = list(
      description = "Absolute neutrophil count before the cycle-1 pemetrexed dose (baseline)",
      units = "10^3 cells/uL (equivalent to 10^9/L)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed and REQUIRED: it is the model's value at time 0 and the",
        "level from which the decline parabola starts. Morita 2026 Methods",
        "2.4 defines Base as the neutrophil count before pemetrexed",
        "administration, typically within 1 week. Supplied in the",
        "paper's reporting unit (Figure 3 axis, x 10^3/uL), so divide a",
        "cells/uL value by 1000 before passing it."
      ),
      source_name = "BASE"
    ),
    BUN = list(
      description = "Blood urea nitrogen at the start of cycle 1",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value. Power effect on the decline curvature",
        "centred at 15.2 mg/dL, `(BUN / 15.2)^-0.975` (Morita 2026",
        "Equation 7; Text S1 `BUNP = (BUN/15.2) ** THETA(4)`). Development",
        "cohort mean 15.25 +/- 5.16 mg/dL (Table 1)."
      ),
      source_name = "BUN"
    ),
    CONMED_RAS_INHIBITOR = list(
      description = "Concomitant renin-angiotensin system inhibitor use (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no RAS inhibitor)",
      notes = paste(
        "Time-fixed, recorded at the start of cycle 1. Class indicator:",
        "the paper pools ACE inhibitors and angiotensin II receptor",
        "blockers without naming agents. Multiplies the recovery",
        "curvature by 0.0451 (Equation 8; Text S1 `rasP = ras1 + ras2 *",
        "THETA(5)`). 25 of 86 development-cohort patients (29.1%, Table 1)."
      ),
      source_name = "RAS"
    ),
    AUC_PEMETREXED = list(
      description = "Cycle-1 pemetrexed plasma AUC from 0 to 24 h predicted by the authors' PK-Sim PBPK model",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Enters ONLY through the indicator AUC_PEMETREXED > 175",
        "(Equation 9; Text S1 column AUC_2), which multiplies the nadir",
        "time by 1.44; the continuous value has no other effect. The paper",
        "derives it per patient from a PK-Sim PBPK model (Methods 2.2,",
        "Table S1) using race, age, sex, height, weight, BMI, BSA, eGFR",
        "and the pemetrexed dose, never from measured concentrations.",
        "Development-cohort values span about 108 to 285 ug*h/mL",
        "(Figure S2), about 80% of patients above 175."
      ),
      source_name = "AUC24"
    )
  )

  covariatesDataExcluded <- list(
    DOSE_PEMETREXED_MG = list(
      description = "Administered pemetrexed dose (mg)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Retained only in the Step-1 model (no PK covariates):",
        "A1 = 0.717 * (Dose/760)^1.27 (Equation 4). The final Step-2 model",
        "replaces it with the PBPK-derived AUC0-24 indicator and does not",
        "carry dose (Results 3.3), so it is documented here only."
      ),
      source_name = "Pem_DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 86L,
    n_studies = 1L,
    age_range = "mean 68.24 +/- 7.85 years (Table 1, development cohort)",
    weight_range = "mean 58.65 +/- 11.01 kg (Table 1)",
    sex_female_pct = 27.9,
    race_ethnicity = "Japanese (single centre)",
    disease_state = "Non-squamous non-small cell lung cancer, first cycle of pemetrexed-platinum chemotherapy, no prior pemetrexed or platinum",
    dose_range = paste(
      "Pemetrexed 500 mg/m2 standard (actual 469.54 +/- 51.35 mg/m2,",
      "760.64 +/- 127.05 mg) on day 1 of a 21-day cycle with cisplatin",
      "75 mg/m2 (38 patients) or carboplatin AUC 5 mg*min/mL (48 patients)"
    ),
    regions = "Japan (Saiseikai Kumamoto Hospital, January 2015 to December 2023)",
    renal_function = "CLcr 77.32 +/- 24.38 mL/min; eGFR (de-indexed) 71.09 +/- 19.14 mL/min; BUN 15.25 +/- 5.16 mg/dL (Table 1)",
    co_medication = "NSAIDs 26.7%, RAS inhibitors 29.1%, diuretics 9.3% (Table 1)",
    n_observations = 316L,
    notes = paste(
      "Retrospective single-centre real-world development cohort; 316",
      "neutrophil observations from the first cycle, censored after any",
      "G-CSF administration (Methods 2.1). An external validation cohort",
      "of 83 patients who also received pembrolizumab or bevacizumab",
      "(375 observations) was not used for estimation (Methods 2.6)."
    )
  )

  ini({
    # Typical values at BUN = 15.2 mg/dL, no RAS inhibitor, AUC0-24 <= 175 ug*h/mL.
    lcurv_decline <- log(0.334)
    label("Curvature of the ANC decline parabola, A1 (10^3 cells/uL/week^2)") # Table S7 Step 2 'A1' = 0.334; Equation 7
    lcurv_recovery <- log(0.0443)
    label("Curvature of the ANC recovery parabola, A2 (10^3 cells/uL/week^2)") # Table S7 Step 2 'A2' = 0.0443; Equation 8
    ltnadir <- log(1.39)
    label("Time of the ANC nadir and start of recovery, B (week)") # Table S7 Step 2 'B' = 1.39; Equation 9

    e_bun_curv_decline <- -0.975
    label("Power exponent of BUN/15.2 on the decline curvature (unitless)") # Table S7 Step 2 'BUN' = -0.975; Equation 7
    e_ras_curv_recovery <- 0.0451
    label("Multiplier on the recovery curvature with a concomitant RAS inhibitor (unitless)") # Table S7 Step 2 'RAS inhibitor' = 0.0451; Equation 8
    e_auc_tnadir <- 1.44
    label("Multiplier on the nadir time when pemetrexed AUC0-24 > 175 ug*h/mL (unitless)") # Table S7 Step 2 'AUC0-24' = 1.44; Equation 9

    # Table S7 reports omega as the SD of eta (Results 3.3: eta has variance omega^2);
    # the variances below are the squared SDs.
    etalcurv_decline ~ 0.697225 # Table S7 Step 2 'omega A1' = 0.835 (SD); 0.835^2
    etalcurv_recovery ~ 3.8416 # Table S7 Step 2 'omega A2' = 1.96 (SD); 1.96^2
    etaltnadir ~ 0.021904 # Table S7 Step 2 'omega B' = 0.148 (SD); 0.148^2

    propSd <- 0.271
    label("Proportional residual error (fraction)") # Table S7 Step 2 'epsilon' = 0.271; Text S1 Y = IPRED * exp(EPS(1))
  })

  model({
    # Individual parameters (Equations 7-9; Text S1 $PRED)
    auc_gt175 <- 0
    if (AUC_PEMETREXED > 175) auc_gt175 <- 1
    curv_decline <- exp(lcurv_decline + etalcurv_decline) *
      (BUN / 15.2)^e_bun_curv_decline
    curv_recovery <- exp(lcurv_recovery + etalcurv_recovery) *
      e_ras_curv_recovery^CONMED_RAS_INHIBITOR
    tnadir <- exp(ltnadir + etaltnadir) * e_auc_tnadir^auc_gt175

    # Both parabolas share the vertex (tnadir, anc_nadir), so the decline starts
    # at NEUT at time 0 and the two branches meet at the nadir (Text S1 N1, N2).
    # Nothing bounds the nadir below 0, as in the source.
    anc_nadir <- NEUT - curv_decline * tnadir^2
    if (time < tnadir) {
      ANC <- curv_decline * (time - tnadir)^2 + anc_nadir
    } else {
      ANC <- curv_recovery * (time - tnadir)^2 + anc_nadir
    }
    ANC ~ prop(propSd)
  })
}
