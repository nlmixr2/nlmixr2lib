Wattanakul_2026_paracetamol_fcta <- function() {
  description <- "Log-logistic time-to-event model for fever clearance time A (FCT-A, time from enrollment to the first temperature below 37.5 C) in Malaysian patients with Plasmodium knowlesi malaria and temperature >= 37.5 C at enrollment, randomised to oral paracetamol 1 g every 6 h for 72 h or no paracetamol (PACKNOW trial; Wattanakul 2026). The log-logistic hazard (scale 0.0940 1/h, shape 2.77) is multiplied by (1 + 0.158 * AUC_PARACETAMOL / 100), a 15.8% faster rate of fever clearance (FCT-A) per 100 mg*h/L of individual paracetamol AUC0-72h; untreated patients (AUC_PARACETAMOL = 0) take the baseline hazard. No IIV. Paracetamol exposure enters only as the per-subject AUC_PARACETAMOL covariate, computed from the companion PK model modellib('Wattanakul_2026_paracetamol'). Forward simulation exposes `hazard` and `sur` (probability that fever clearance (FCT-A) has not yet occurred)."
  reference <- paste(
    "Wattanakul T, Cooper DJ, Plewes K, Grigg MJ, Rajahram GS, William T,",
    "Dondorp AM, Edstein MD, Birrell GW, Anstey NM, Hoglund RM, Barber BE,",
    "Tarning J. Population Pharmacokinetics and Pharmacodynamics of",
    "Paracetamol in Malaysian Patients With Plasmodium knowlesi Malaria.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70283.",
    "doi:10.1002/psp4.70283.",
    "Exposure driver from the companion PK model; see",
    "modellib('Wattanakul_2026_paracetamol').",
    sep = " "
  )
  vignette <- "Wattanakul_2026_paracetamol"
  units <- list(
    time = "h",
    dosing = "n/a (no dosing events; paracetamol exposure is supplied as the AUC_PARACETAMOL covariate in mg*h/L)",
    concentration = "probability (sur is the probability that fever clearance (FCT-A) has not yet occurred)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of fever clearance (FCT-A)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_PARACETAMOL = list(
      description = "Individual model-predicted paracetamol AUC from 0 to 72 h after the first dose.",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject. Derived in the source from each patient's empirical-Bayes parameters in the final PK model (modellib('Wattanakul_2026_paracetamol')); 0 for patients who received no paracetamol. Enters linearly and uncentred as a multiplier on the hazard (supplement fever stream 'EFF = (1+BETA*AUC72MG)'). Paracetamol-arm median (IQR) 408 (331-504) mg*h/L (Table 2).",
      source_name = "AUC72MG"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 194L,
    n_studies = 1L,
    n_observations = "One fever clearance (FCT-A) time (event or censoring) per patient.",
    age_range = "5 years and older (inclusion criterion)",
    weight_range = "Not reported separately for the fever-clearance subset; see Table 1 for the full PK/PD cohort.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Malaysian (Sabah); not further broken down.",
    disease_state = "Hospitalised patients with PCR-confirmed Plasmodium knowlesi monoinfection and a temperature of at least 37.5 C at enrollment.",
    dose_range = "Paracetamol 1 g orally every 6 h for 72 h (n = 102) or no paracetamol (n = 92).",
    regions = "Sabah, Malaysia.",
    notes = "PACKNOW trial (NCT03056391). Of the 372 PK/PD patients, 195 had a temperature >= 37.5 C at enrollment; one control patient with missing data was excluded, leaving 194 (Results 3.2.2). FCT-A is the time from enrollment to the first recorded body temperature below 37.5 C (supplement S3). Both arms were fitted simultaneously in NONMEM with the Laplacian method."
  )

  ini({
    # Wattanakul 2026 supplement Table S1, 'Fever clearance time A'
    llambda <- log(0.0940)
    label("Scale parameter lambda of the log-logistic hazard (1/h)") # Table S1 Fever clearance time A: lambda = 0.0940 (RSE 6.9%)
    lgamma <- log(2.77)
    label("Shape parameter gamma of the log-logistic hazard (unitless)") # Table S1 Fever clearance time A: gamma = 2.77 (RSE 4.1%)
    e_auc_paracetamol_haz <- 0.00158
    label("Linear effect of paracetamol AUC0-72h on the hazard (fractional increase per mg*h/L)") # Table S1 Fever clearance time A: theta_AUC72H = 0.158 per 100 mg*h/L (RSE 45.7%), divided by 100
    # No IIV ('$OMEGA 0 FIX ; IIV_BASELINE HAZARD') and no residual error:
    # the likelihood is carried by the survival / hazard density.
  })

  model({
    lambda <- exp(llambda)
    gamma <- exp(lgamma)
    eff_auc <- 1 + e_auc_paracetamol_haz * AUC_PARACETAMOL

    # Supplement S3 log-logistic hazard; the 1e-6 offset is the stream's
    # 'DEL = 1E-6', which keeps (lambda * t)^(gamma - 1) finite at t = 0.
    hazard <- eff_auc * lambda * gamma * (lambda * t + 1e-6)^(gamma - 1) /
      (1 + (lambda * t + 1e-6)^gamma)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
