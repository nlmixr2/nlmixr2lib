Sturm_2017_sembragiline_hippocampus_ec <- function() {
  description <- paste(
    "Direct-effect Emax model of brain monoamine oxidase type B (MAO-B) enzyme inhibition in the hippocampus as a function of sembragiline plasma concentration in healthy elderly control (EC) subjects (Sturm 2017).",
    "Inhibition was measured by [11C]-L-deprenyl-D2 PET as the percentage reduction of the regional specific-binding index lambda*k3 from the drug-free baseline scan, and was paired with the sembragiline plasma concentration drawn just before the second (steady-state, trough) scan.",
    "Inhibition = Emax * Cp / (Cp + EC50), with zero inhibition at zero concentration and no Hill coefficient.",
    "PD-only model: the sembragiline plasma concentration is supplied as the covariate CP_SEMBRAGILINE_NGML (ng/mL). The paper fit one Emax model per population and brain region of interest; this file is the EC / hippocampus fit."
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
        "Pre-scan steady-state trough plasma concentration mean (SD) 3.29 (1.10) ng/mL at 0.2 mg and 33.9 (13.7) ng/mL at 1 mg (Supplementary Table 2).",
        "Set to 0 for drug-free records (inhibition is then 0)."
      ),
      source_name = "Cp"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "mean 70.5 (SD 10.2) years; protocol range 50-80 years",
    weight_range = "mean 84.8 (SD 8.0) kg",
    sex_female_pct = 16.7,
    race_ethnicity = "not reported",
    disease_state = "Healthy elderly controls with no evidence of Alzheimer's disease or other cognitive impairment (mean MMSE 28.3, SD 1.6)",
    dose_range = "Sembragiline orally once daily: 0.2 mg/day on days 2-7 after a 5 mg loading dose on day 1 (n = 3), or 1 mg/day for 14 days (n = 3).",
    regions = "Sweden (Stockholm: enrolment at Karolinska Hospital, Huddinge; PET at the Karolinska Institutet PET Centre)",
    notes = paste(
      "Fit to 6 steady-state (inhibition, plasma concentration) pairs from 6 healthy elderly control (EC) subjects without cognitive impairment in an open-label, parallel-group, multiple-dose PET study; demographics are Sturm 2017 Supplementary Table 1.",
      "The Emax model was fitted per population and per brain region of interest with Phoenix WinNonlin 6.2 (a naive-pooled nonlinear regression, not a population mixed-effects fit).",
      "All subjects were non-smokers and followed a low-tyramine diet."
    )
  )

  ini({
    # Simple Emax model fitted by Sturm 2017 (Methods, 'Sembragiline plasma
    # concentration-MAO-B enzyme inhibition relationship') with Phoenix
    # WinNonlin 6.2, one fit per population and region of interest. Both
    # typical values are Sturm 2017 Table 2, row 'EC / Hippocampus'.
    # The parenthesised CV% in Table 2 is the precision of each estimate (the
    # WinNonlin coefficient of variation of the parameter estimate), not
    # between-subject variability.

    lemax <- log(89.1)
    label("Maximum MAO-B enzyme inhibition in the hippocampus (Emax, %)")
    # Table 2, 'EC / Hippocampus': Emax = 89.1% (CV% 9.6).

    lec50 <- log(2.05)
    label("Sembragiline plasma concentration giving half-maximal MAO-B inhibition in the hippocampus (EC50, ng/mL)")
    # Table 2, 'EC / Hippocampus': EC50 = 2.05 ng/mL (CV% 33.4).

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
