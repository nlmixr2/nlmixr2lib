Perello_2024_hexasodiumFytate_hydroxyapatite <- function() {
  description <- paste0(
    "Simple Emax PK/PD model relating the individual average end-of-infusion ",
    "plasma concentration (Cmax) of hexasodium fytate (SNF472, measured as ",
    "the free acid myo-inositol hexaphosphate, IP6) to the individual average ",
    "ex vivo inhibition of hydroxyapatite crystallization in plasma, in adults ",
    "on maintenance hemodialysis with cardiovascular calcification (Perello ",
    "2024, phase 2b CaLIPSO PK/PD per-protocol subset, n = 55: placebo 15, ",
    "300 mg 20, 600 mg 20, infused thrice weekly during dialysis for 52 ",
    "weeks). Inhibition (%) = E0 + Emax * Cmax / (EC50 + Cmax) with the Hill ",
    "coefficient fixed at 1, E0 = 8.78%, Emax = 75.9% and EC50 = 5.48 uM ",
    "(Table 3); EC80 = 4 * EC50 = 21.9 uM. One record per subject (averages ",
    "over weeks 1, 10, 22 and 52), fitted by nonlinear least squares with an ",
    "additive residual; there is no between-subject random effect. There is ",
    "no PK layer and no ODE: the exposure metric is supplied as the CMAX ",
    "data column. The companion exposure-efficacy model is ",
    "modellib('Perello_2024_hexasodiumFytate_cac')."
  )
  reference <- paste(
    "Perello J, Alberti J, Torres JV, Ferrer MD, Perez MM, Bassissi F,",
    "Gold A, Raggi P, Chertow GM, Salcedo C.",
    "Hexasodium fytate exposure-response correlations in a randomized,",
    "placebo-controlled study of patients on dialysis with cardiovascular",
    "calcification. Front Pharmacol. 2024;15:1325186.",
    "doi:10.3389/fphar.2024.1325186.",
    sep = " "
  )
  vignette <- "Perello_2024_hexasodiumFytate"
  units <- list(
    time = "n/a (static regression of 52-week subject averages; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "hydroxyapatite_inhibition (% inhibition of hydroxyapatite crystallization, ex vivo plasma assay)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual average hexasodium fytate end-of-infusion plasma concentration (Cmax) across the week 1, 10, 22 and 52 visits. Supplied as data: this model has no PK layer.",
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma IP6 (hexasodium fytate free acid, molecular weight",
        "660 g/mol, so 1 uM = 0.660 mg/L) measured by LC-MS/MS in samples",
        "drawn within the last 10 min of a 2.5-h intradialytic infusion",
        "(Methods 2.1-2.2) -- an OBSERVED end-of-infusion concentration,",
        "not a model prediction, averaged arithmetically per subject over",
        "the evaluable visits. No accumulation over 52 weeks of thrice-weekly",
        "dosing (Results 3.2), so the average is effectively a steady-state",
        "Cmax. Values below the lower limit of quantification (0.758 uM)",
        "were set to LQ/2 = 0.379 uM, which is the value every placebo",
        "subject carries (Table 3 footnote). Observed group means 14.7 uM",
        "(SD 9.6) at 300 mg and 45.8 uM (SD 21.7) at 600 mg (Table 2)."
      ),
      source_name = "Cmax (average Cmax)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    n_observations = "55 subject-level records (one per subject: average PD inhibition vs average Cmax), from 116 paired PK/PD end-of-infusion samples (Table 1, Table 3)",
    disease_state = "adults with kidney failure on maintenance hemodialysis and coronary artery calcification (baseline Agatston CAC score 100-3,500), per-protocol PK/PD subset of the phase 2b CaLIPSO study (NCT02966028)",
    dose_range = "placebo (15), hexasodium fytate 300 mg (20) or 600 mg (20) infused over 2.5 h into the dialysis line thrice weekly for 52 weeks",
    notes = paste0(
      "Demographics of the PK/PD subset are not reported in Perello 2024 ",
      "(the full CaLIPSO cohort of 274 patients is described in Raggi 2020). ",
      "Only study sites with a -70 C freezer participated in the PK/PD ",
      "substudy (Methods 2.1). Fitted in SAS 9.4; Table 4 reports the ",
      "BLQ = 0 and BLQ = LQ sensitivity fits, which gave EC50 6.33 and ",
      "4.66 uM respectively; the packaged values are the primary BLQ = LQ/2 ",
      "fit."
    )
  )

  ini({
    # Perello 2024 Table 3 ('Final PK/PD model using average Cmax and
    # average inhibition of hydroxyapatite crystallization'), BLQ = LQ/2.
    # Equation: Methods 2.3, PD_i = E0 + Cmax_i^H * Emax / (Cmax_i^H + EC50^H)
    # + eps_i, with the final simple Emax model taking H = 1 (Results 3.4,
    # and EC80 = 21.9 = 4 * 5.48 in Table 3 confirms H = 1).

    # E0 can be negative (Table 4 BLQ = LQ fit: 2.59, 95% CI -8.43 to 13.6),
    # so it is kept on the linear scale.
    e0 <- 8.78
    label("Inhibition of hydroxyapatite crystallization at zero exposure, E0 (%)") # Table 3, E0 = 8.78 % inhibition, SE 3.82, 95% CI 1.13-16.4
    lemax <- log(75.9)
    label("Maximal drug-attributable inhibition of hydroxyapatite crystallization, Emax (%)") # Table 3, Emax = 75.9 % inhibition, SE 4.81, 95% CI 66.3-85.6
    lec50 <- log(5.48)
    label("Average Cmax producing half of Emax, EC50 (uM)") # Table 3, EC50 = 5.48 uM, SE 1.76, 95% CI 1.95-9.01

    # Additive residual on the % inhibition scale. The SD is the square
    # root of the Table 3 'Error' mean square, 126.4 on 52 DF:
    # sqrt(126.4) = 11.24 % inhibition.
    addSd <- 11.24
    label("Additive residual error SD on % inhibition of hydroxyapatite crystallization (%)") # Table 3, Error mean square 126.4 (SS 6,574 / DF 52); SD = sqrt(126.4)
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    hydroxyapatite_inhibition <- e0 + emax * CMAX / (ec50 + CMAX)

    hydroxyapatite_inhibition ~ add(addSd)
  })
}
