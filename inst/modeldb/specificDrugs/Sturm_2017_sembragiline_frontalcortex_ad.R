Sturm_2017_sembragiline_frontalcortex_ad <- function() {
  description <- paste(
    "Direct-effect Emax model of brain monoamine oxidase type B (MAO-B) enzyme inhibition in the frontal cortex as a function of sembragiline plasma concentration in patients with Alzheimer's disease (AD) (Sturm 2017).",
    "Inhibition was measured by [11C]-L-deprenyl-D2 PET as the percentage reduction of the regional specific-binding index lambda*k3 from the drug-free baseline scan, and was paired with the sembragiline plasma concentration drawn just before the second (steady-state, trough) scan.",
    "Inhibition = Emax * Cp / (Cp + EC50), with zero inhibition at zero concentration and no Hill coefficient.",
    "PD-only model: the sembragiline plasma concentration is supplied as the covariate CP_SEMBRAGILINE_NGML (ng/mL). The paper fit one Emax model per population and brain region of interest; this file is the AD / frontal cortex fit."
  )
  reference <- paste(
    "Sturm S, Forsberg A, Nave S, Stenkrona P, Seneca N, Varrone A, Comley RA, Fazio P,",
    "Jamois C, Nakao R, Ejduk Z, Al-Tawil N, Akenine U, Halldin C, Andreasen N, Ricci B (2017).",
    "Positron emission tomography measurement of brain MAO-B inhibition in patients with",
    "Alzheimer's disease and elderly controls after oral administration of sembragiline.",
    "Eur J Nucl Med Mol Imaging 44(3):382-391.",
    "doi:10.1007/s00259-016-3510-6.",
    sep = " "
  )
  vignette <- "Sturm_2017_sembragiline"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external sembragiline plasma-concentration covariate)",
    concentration = "(observation maob_inhibition is brain MAO-B enzyme inhibition in %; the driving covariate CP_SEMBRAGILINE_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_SEMBRAGILINE_NGML = list(
      description = "Sembragiline plasma concentration at the time of the MAO-B inhibition measurement, supplied as a (time-varying) covariate from observed samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives Inhibition = emax * CP_SEMBRAGILINE_NGML / (CP_SEMBRAGILINE_NGML + ec50).",
        "In Sturm 2017 this was the steady-state trough plasma concentration measured by validated LC-MS/MS just before the post-dose PET scan (Methods, 'Assessment of PK' and 'Sembragiline plasma concentration-MAO-B enzyme inhibition relationship'); Cp = 0 at the baseline scan.",
        "Pre-scan steady-state trough plasma concentration mean (SD) 2.55 (0.191) ng/mL at 0.1 mg, 7.06 (3.13) ng/mL at 0.2 mg, 41.9 (10.2) ng/mL at 1 mg and 96.7 (69.7) ng/mL at 5 mg (Supplementary Table 2).",
        "Set to 0 for drug-free records (inhibition is then 0)."
      ),
      source_name = "Cp"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "mean 65.0 (SD 8.4) years; protocol range 50-80 years",
    weight_range = "mean 72.3 (SD 14.8) kg",
    sex_female_pct = 80,
    race_ethnicity = "not reported",
    disease_state = "Probable Alzheimer's disease (NINCDS/ADRDA and DSM-IV), MMSE 17-26 (mean 22.1, SD 3.4), on stable background AD therapy (galantamine, donepezil, memantine, or galantamine + memantine)",
    dose_range = "Sembragiline orally once daily: 0.1 mg/day for 13-16 days after a 5 mg loading dose on day 1 (n = 2), 0.2 mg/day on days 2-7 after a 5 mg loading dose (n = 3), 1 mg/day for 14 days (n = 3; 2 with a post-dose scan), or 5 mg/day for 6-9 days (n = 2).",
    regions = "Sweden (Stockholm: enrolment at Karolinska Hospital, Huddinge; PET at the Karolinska Institutet PET Centre)",
    notes = paste(
      "Fit to 9 steady-state (inhibition, plasma concentration) pairs from 10 patients with probable Alzheimer's disease (AD, MMSE 17-26) in an open-label, parallel-group, multiple-dose PET study; demographics are Sturm 2017 Supplementary Table 1.",
      "The Emax model was fitted per population and per brain region of interest with Phoenix WinNonlin 6.2 (a naive-pooled nonlinear regression, not a population mixed-effects fit).",
      "All subjects were non-smokers and followed a low-tyramine diet."
    )
  )

  ini({
    # Simple Emax model fitted by Sturm 2017 (Methods, 'Sembragiline plasma
    # concentration-MAO-B enzyme inhibition relationship') with Phoenix
    # WinNonlin 6.2, one fit per population and region of interest. Both
    # typical values are Sturm 2017 Table 2, row 'AD / Frontal cortex'.
    # The parenthesised CV% in Table 2 is the precision of each estimate (the
    # WinNonlin coefficient of variation of the parameter estimate), not
    # between-subject variability.

    lemax <- log(81.1)
    label("Maximum MAO-B enzyme inhibition in the frontal cortex (Emax, %)")
    # Table 2, 'AD / Frontal cortex': Emax = 81.1% (CV% 3.2).

    lec50 <- log(1.06)
    label("Sembragiline plasma concentration giving half-maximal MAO-B inhibition in the frontal cortex (EC50, ng/mL)")
    # Table 2, 'AD / Frontal cortex': EC50 = 1.06 ng/mL (CV% 19.6).

    # Sturm 2017 reports no residual error model and no between-subject
    # variability for the WinNonlin fit, so the residual SD is fixed to 0 and
    # the model predicts the typical value only.
    addSd <- fixed(0)
    label("Additive residual error on MAO-B inhibition (percentage points; not reported by Sturm 2017)")
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Sturm 2017 Methods equation: Occ% = Emax * Cp / (Cp + EC50), Cp in
    # ng/mL; 0% at Cp = 0 and Emax at Cp = infinity.
    maob_inhibition <- emax * CP_SEMBRAGILINE_NGML / (CP_SEMBRAGILINE_NGML + ec50)
    maob_inhibition ~ add(addSd)
  })
}
