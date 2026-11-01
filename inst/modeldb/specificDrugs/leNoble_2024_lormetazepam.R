leNoble_2024_lormetazepam <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for enteral (nasogastric) lormetazepam in mechanically ventilated ICU adults with COVID-19-associated ARDS (le Noble 2024), developed by an iterative two-stage Bayesian procedure in Edsim++/KinPop. Apparent volume scales linearly and apparent clearance allometrically (exponent 0.75) with body weight per 70 kg. Correlated inter-individual variability on Ka, V/F and CL/F; proportional residual error. Creatinine clearance and CRP were screened and not retained."
  reference <- "le Noble JLML, Shudofsky KN, Foudraine N, Punt N, Janssen PKJ. Pharmacokinetics of Enteral Lormetazepam in Mechanically Ventilated ICU Patients with COVID-19: An Adjunct Sedative Study. Clin Pharmacokinet. 2024;63(12):1769-1776. doi:10.1007/s40262-024-01455-3. PMCID PMC11649707. Parameter estimates from main-text Table 2; full-precision estimates and the population covariance matrix from Electronic Supplementary Material Part 1 (KinPop output of the final model '1-POP-All', 'Cycles = 118')."
  vignette <- "leNoble_2024_lormetazepam"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Most recent pre-ICU body weight (Methods, Clinical Parameters). Reference 70 kg, the unit denominator of the Table 2 parameters ('L/70 kg' for V/F, 'L/h/70 kg0.75' for CL/F): linear scaling of V/F and exponent 0.75 on CL/F. The paper does not report a weight summary; Table 1 gives median height 172.0 cm (IQR 163.5-180.5) and median ICU-admission BMI 29.1 kg/m^2 (IQR 27.4-31.1).",
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance on ICU admission",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested by linear regression against individual CL/F and not retained ('did not affect the total body clearance of lormetazepam and were not included in the final model', Results, Analysis of Covariates; Fig. 1). Median 79 mL/min (IQR 64-114), Table 1. Patients on continuous veno-venous hemofiltration were excluded.",
      source_name = "CLCr"
    ),
    CRP = list(
      description = "C-reactive protein on ICU admission (standard assay, baseline)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested by linear regression against individual CL/F and not retained (Results, Analysis of Covariates; Fig. 2). Median 164 mg/L (IQR 113.5-210.5), Table 1.",
      source_name = "CRP"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lormetazepam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lormetazepam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "mean 61.6 years (SD 13)",
    weight_range = "not reported; median height 172.0 cm (IQR 163.5-180.5), median ICU-admission BMI 29.1 kg/m^2 (IQR 27.4-31.1)",
    sex_female_pct = 26.7,
    disease_state = "Critically ill adults on invasive mechanical ventilation in the prone position with severe COVID-19-associated ARDS (PaO2/FiO2 < 200); median APACHE IV 57 (IQR 55-68)",
    dose_range = "Lormetazepam 2 mg twice daily (08:00 and 20:00) via nasogastric tube, as an adjunct to IV propofol, midazolam and dexmedetomidine",
    regions = "Netherlands (VieCuri Medical Center, Venlo)",
    renal_function = "Creatinine clearance median 79 mL/min (IQR 64-114); serum creatinine median 70 umol/L (IQR 58.5-83); patients on CVVH excluded",
    co_medication = "Hydroxychloroquine 80%, corticosteroids 60%; prokinetics (metoclopramide, erythromycin); fentanyl, acetaminophen, S-ketamine",
    notes = "Single-centre observational prospective cohort, 9-20 April 2020 (Table 1). 18 patients sampled, 3 on CVVH excluded. Residual arterial blood samples from routine sampling at 06:30, 16:00 and 22:00; median 10 samples per patient (range 2-21); 116 concentrations in the population fit (Supplement Part 1, N = 116). Assay range 1-20 ug/L by UHPLC-MS/MS; data below the LOQ excluded."
  )

  ini({
    # ------------------------------------------------------------------
    # Table 2 'Estimates' column (Mean, CV%, RSE%), carried at the full
    # precision of Supplement Part 1, 'Cycles = 118'. KinPop reports each
    # log-normally distributed parameter ('Log: Yes') as a population mean
    # +/- SD. The mean is carried as the typical value, the same convention
    # as the Edsim++/KinPop-built Meenks_2023_flucloxacillin model by the
    # same group. Fig. 3 confirms it: these means at 70 kg reproduce the
    # 6 ug/L steady-state troughs of all three published regimens.
    # ------------------------------------------------------------------
    lka <- log(6.3778)
    label("Absorption rate constant Ka (1/h)") # Table 2 Ka = 6.4 1/h (RSE 39.4%); Supplement Part 1 PO.k = 6.3778
    lvc <- log(207.3061)
    label("Apparent volume of distribution V/F for a 70 kg patient (L)") # Table 2 V/F = 207 L/70 kg (RSE 13.1%); Supplement Part 1 C01.V = 207.3061
    lcl <- log(14.5041)
    label("Apparent clearance CL/F for a 70 kg patient (L/h)") # Table 2 CL/F = 14.5 L/h/70 kg0.75 (RSE 12.9%); Supplement Part 1 TE.CL = 14.5041

    # Weight scaling is part of the parameter units, not estimated: V/F in
    # 'L/70 kg' and CL/F in 'L/h/70 kg0.75' (Table 2).
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V/F (unitless)") # Table 2 unit 'L/70 kg'
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)") # Table 2 unit 'L/h/70 kg0.75'

    # ------------------------------------------------------------------
    # IIV. Supplement Part 1 prints the population ('Pop') covariance of
    # the parameters relative to their means (diagonal = CV^2: V 0.256,
    # Ka 2.33, CL 0.250) and the correlations V-Ka 0.0681, V-CL 0.4279,
    # Ka-CL -0.0926. With CV = SD / mean from 'Cycles = 118'
    # (V 104.8195/207.3061 = 0.505627, Ka 9.7262/6.3778 = 1.525009,
    # CL 7.2573/14.5041 = 0.500362), each log-normal element is
    # omega_ij = log(1 + rho_ij * CV_i * CV_j):
    #   Ka   log(1 + 1.525009^2)                    = 1.201666
    #   Ka-V log(1 + 0.0681 * 1.525009 * 0.505627)  = 0.051179
    #   V    log(1 + 0.505627^2)                    = 0.227660
    #   Ka-CL log(1 - 0.0926 * 1.525009 * 0.500362) = -0.073280
    #   V-CL log(1 + 0.4279 * 0.505627 * 0.500362)  = 0.102789
    #   CL   log(1 + 0.500362^2)                    = 0.223433
    # ------------------------------------------------------------------
    etalka + etalvc + etalcl ~ c(
      1.201666,
      0.051179, 0.227660,
      -0.073280, 0.102789, 0.223433
    ) # Table 2 CV% Ka 153, V/F 51, CL/F 50; Supplement Part 1 'Pop' covariances and correlations

    propSd <- 0.4931
    label("Proportional residual error (fraction)") # Table 2 'Prop' = 49% (RSE 6.4%); Supplement Part 1 SD.Sd1 = 0.4931
  })

  model({
    # Individual parameters, scaled to a 70 kg reference (Table 2 units)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl

    kel <- cl / vc

    # One-compartment model with first-order absorption from the enteral
    # (nasogastric) dose; bioavailability is absorbed into V/F and CL/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and volume in L give mg/L; * 1000 reports ug/L as in the
    # paper's assay and figures.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
