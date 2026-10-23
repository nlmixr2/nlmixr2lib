Wattanakul_2026_paracetamol_creatinine <- function() {
  description <- "Exposure-response model for the percentage change from baseline in serum creatinine over the first 72 h in Malaysian patients with Plasmodium knowlesi malaria randomised to oral paracetamol 1 g every 6 h for 72 h or no paracetamol (PACKNOW trial; Wattanakul 2026). Creatinine falls exponentially from baseline to a plateau of -9.9% at a rate constant of 0.0477 1/h. The rate constant rises by 25.4% per 10-fold increase in the individual paracetamol AUC0-72h (log10 scale, uncentred; zero in the control arm), and the plateau is scaled by the enrollment KDIGO acute-kidney-injury stage (multipliers 2.77, 4.93 and 7.52 for stages 1, 2 and 3). Additive IIV on the rate constant and plateau; additive residual error. Paracetamol exposure enters only as the per-subject AUC_PARACETAMOL covariate, computed from the companion PK model modellib('Wattanakul_2026_paracetamol')."
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
    concentration = "percent (crcfb is the percentage change from baseline in serum creatinine)"
  )

  covariateData <- list(
    AUC_PARACETAMOL = list(
      description = "Individual model-predicted paracetamol AUC from 0 to 72 h after the first dose.",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject. Derived in the source from each patient's empirical-Bayes parameters in the final PK model (modellib('Wattanakul_2026_paracetamol')). Enters as log10(AUC_PARACETAMOL), uncentred (supplement creatinine stream 'KAPPALAUC72 = (1 + THETA(4)*(LAUC72MG))'); set to 0 for patients who received no paracetamol, for whom the model uses log10-AUC = 0, i.e. no exposure effect. Paracetamol-arm median (IQR) 408 (331-504) mg*h/L (Table 2).",
      source_name = "AUC72MG (log10 form: LAUC72MG)"
    ),
    KDIGO_AKI_1 = list(
      description = "KDIGO acute-kidney-injury stage 1 at enrollment (1 = yes).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Time-fixed enrollment value, KDIGO serum-creatinine criteria only (urine output not collected). KDIGO_AKI_1 = (STAGE == 1). Prevalence 21.3% (paracetamol) and 27.5% (control), Table 1.",
      source_name = "STAGE"
    ),
    KDIGO_AKI_2 = list(
      description = "KDIGO acute-kidney-injury stage 2 at enrollment (1 = yes).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Time-fixed enrollment value. KDIGO_AKI_2 = (STAGE == 2). Prevalence 4.4% (paracetamol) and 2.6% (control), Table 1.",
      source_name = "STAGE"
    ),
    KDIGO_AKI_3 = list(
      description = "KDIGO acute-kidney-injury stage 3 at enrollment (1 = yes).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (stage 0 when KDIGO_AKI_1, KDIGO_AKI_2 and KDIGO_AKI_3 are all 0)",
      notes = "Time-fixed enrollment value. KDIGO_AKI_3 = (STAGE == 3). Two paracetamol-arm patients and no control patients (Table 1); patients who underwent renal replacement therapy were excluded from the analysis.",
      source_name = "STAGE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 372L,
    n_studies = 1L,
    n_observations = "2365 serum creatinine measurements, at enrollment and every 12 h up to 72 h, expressed as percentage change from each patient's enrollment value.",
    age_range = "5 years and older (inclusion criterion); IQR 26-49 (paracetamol) and 24-48 (control) years",
    age_median = "35 (paracetamol) and 36 (control) years",
    weight_range = "IQR 53.5-68.4 (paracetamol) and 52.0-66.5 (control) kg",
    weight_median = "61.0 (paracetamol) and 58.0 (control) kg",
    sex_female_pct = 15.6,
    race_ethnicity = "Malaysian (Sabah); not further broken down.",
    disease_state = "Hospitalised patients with PCR-confirmed Plasmodium knowlesi monoinfection; 7.5% WHO-defined severe malaria and 28.5% with acute kidney injury (KDIGO stage 1-3) at enrollment. Baseline creatinine median 83 (paracetamol) and 84 (control) umol/L.",
    dose_range = "Paracetamol 1 g orally every 6 h for 72 h (n = 183) or no paracetamol (n = 189).",
    regions = "Sabah, Malaysia.",
    notes = "PACKNOW trial (NCT03056391), both arms. Five patients who underwent renal replacement therapy and one with incomplete data were excluded (Results 3). Demographics from Wattanakul 2026 Table 1; the female percentage is (23 + 35) / 372."
  )

  ini({
    # Wattanakul 2026 Table 3. Random effects in Table 3 are on the
    # arithmetic scale (footnote c) and are reported as standard deviations:
    # the supplement creatinine stream's $OMEGA / $SIGMA values are the
    # squares of the Table 3 entries ($SIGMA 63 -> sqrt = 7.94 = Table 3
    # sigma; $OMEGA 0.00121 -> sqrt = 0.0348 vs Table 3 0.0377).
    kcr <- 0.0477
    label("Rate constant of the exponential change in creatinine towards its plateau, kappa (1/h)") # Table 3: kappa = 0.0477 1/h (RSE 5.6%)
    plateau_cr <- -9.9
    label("Plateau percentage change from baseline in creatinine for KDIGO stage 0 (%)") # Table 3: Plateau = -9.9% (RSE 8.0%)
    e_auc_paracetamol_kcr <- 0.254
    label("Fractional increase in kappa per unit log10 paracetamol AUC0-72h (per log10(mg*h/L))") # Table 3: theta_AUC72H on kappa = 0.254 (RSE 12.9%); Results 3.2.1: 25.4% per 10-fold increase in AUC0-72H
    e_kdigo_aki_1_plateau_cr <- 1.77
    label("Fractional increase in the plateau for KDIGO stage 1 at enrollment (unitless)") # Table 3: theta_KDIGO-1 = 1.77 (RSE 8.9%)
    e_kdigo_aki_2_plateau_cr <- 3.93
    label("Fractional increase in the plateau for KDIGO stage 2 at enrollment (unitless)") # Table 3: theta_KDIGO-2 = 3.93 (RSE 8.5%)
    e_kdigo_aki_3_plateau_cr <- 6.52
    label("Fractional increase in the plateau for KDIGO stage 3 at enrollment (unitless)") # Table 3: theta_KDIGO-3 = 6.52 (RSE 19.2%)

    # Additive IIV (variance = Table 3 SD squared)
    etakcr ~ 0.00142129 # Table 3: IIV kappa = 0.0377 (SD, arithmetic scale); 0.0377^2
    etaplateau_cr ~ 148.84 # Table 3: IIV Plateau = 12.2 (SD, arithmetic scale); 12.2^2

    addSd <- 7.94
    label("Additive residual SD of the percentage change from baseline in creatinine (%)") # Table 3: sigma = 7.94 (RSE 5.7%); = sqrt of stream '$SIGMA 63'
  })

  model({
    # log10 exposure, uncentred; 0 for untreated patients (no effect)
    lauc_paracetamol <- 0
    if (AUC_PARACETAMOL > 0) {
      lauc_paracetamol <- log10(AUC_PARACETAMOL)
    }

    # Stream 'PLATEAUSTAGE = (1 + THETA(.))' for each enrollment KDIGO stage
    kdigo_plateau <- 1 + e_kdigo_aki_1_plateau_cr * KDIGO_AKI_1 +
      e_kdigo_aki_2_plateau_cr * KDIGO_AKI_2 +
      e_kdigo_aki_3_plateau_cr * KDIGO_AKI_3

    kcr_i <- kcr * (1 + e_auc_paracetamol_kcr * lauc_paracetamol) + etakcr
    plateau_cr_i <- plateau_cr * kdigo_plateau + etaplateau_cr

    # Methods Eq. 7 with BaseCR fixed to 0 (stream
    # 'IPRED = (BASE-PLATEAU)*EXP(-KAPPA*TIME) + PLATEAU', BASE = 0)
    crcfb <- plateau_cr_i * (1 - exp(-kcr_i * t))
    crcfb ~ add(addSd)
  })
}
