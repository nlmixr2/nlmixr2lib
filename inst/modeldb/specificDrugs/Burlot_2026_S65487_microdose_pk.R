Burlot_2026_S65487_microdose_pk <- function() {
  description <- "Three-compartment population PK model for intravenous S65487 (Bcl2 inhibitor) with parallel linear and Michaelis-Menten elimination and saturable Michaelis-Menten distribution into the second peripheral compartment, fitted jointly to a 0.1 mg microdose in healthy volunteers and 25-1200 mg therapeutic doses in patients with hematological malignancies (microdose-optimized MO PK model; Burlot 2026)"
  reference <- "Burlot C, Riglet F, Romagnoli M, Leconte A, Delmas A, Schmitt A, Fouliard S. PK and PK/PD Modeling of Bcl2 Inhibitor S65487 in Patients With AML and Investigation of Nonlinearity With Microdosing. CPT Pharmacometrics Syst Pharmacol. 2026;15(7):e70288. doi:10.1002/psp4.70288"
  vignette <- "Burlot_2026_S65487"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    STUDY_S65487_EXT001 = list(
      description = "Indicator for the S65487 phase 0 microdose study EXT 001 (0.1 mg single IV dose in healthy volunteers); 0 = therapeutic-dose studies 002 and 003",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (therapeutic-dose studies 002 / 003)",
      notes = "Selects the residual-error model only: Burlot 2026 fitted separate combined error models to the microdose and therapeutic-dose data because the two used different LC-MS/MS assays (LLOQ 0.001 vs 5 ug/L). No covariate acts on a PK parameter ('No significant covariate was found').",
      source_name = "STUDY (EXT 001)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "S65487", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S65487", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "S65487", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 106,
    n_studies = 3,
    age_range = "19-87 years",
    age_median = "71.5 years",
    weight_range = "45-107 kg",
    weight_median = "73.9 kg",
    sex_female_pct = 35.8,
    disease_state = "Phase 0 (EXT 001, n = 5): healthy male volunteers, single 0.1 mg IV microdose. Phase I (002, n = 60): relapsed/refractory AML, non-Hodgkin lymphoma or multiple myeloma. Phase I/II (003, n = 41): previously untreated AML with azacitidine.",
    dose_range = "0.1 mg single IV microdose; 25-1200 mg IV infusion in the therapeutic-dose studies",
    regions = "Europe (EudraCT 2018-001117-33, 2018-004170-97 and 2020-003061-19)",
    notes = "Burlot 2026 Table 1 (overall column). 2031 S65487 concentrations (84 microdose, 1947 therapeutic dose); 12.7% BLQ handled by M3."
  )

  ini({
    # Table 3 (MO PK model, Monolix 2023R1, SAEM)
    lcl <- log(10.5)
    label("Linear clearance CL (L/h)") # Table 3 'Cl' 10.5 L/h (RSE 6.7%)
    lvc <- log(9.2)
    label("Central volume Vcentral (L)") # Table 3 'Vcentral' 9.2 L (RSE 9.3%)
    lvmax <- log(892.5)
    label("Maximum Michaelis-Menten elimination rate VmaxElimination (ug/h)") # Table 3 'VmaxElimination' 892.5 ug/h (RSE 30.8%)
    lkm <- log(17.7)
    label("Concentration at half-maximal elimination KmElimination (ug/L)") # Table 3 'KmElimination' 17.7 ug/L (RSE 32.9%)
    lq <- log(14.4)
    label("Intercompartmental clearance to the 1st peripheral compartment Q2 (L/h)") # Table 3 'Q2' 14.4 L/h (RSE 13.0%)
    lvp <- log(19.6)
    label("Volume of the 1st peripheral compartment V2 (L)") # Table 3 'V2' 19.6 L (RSE 8.7%)
    lk31 <- log(0.2)
    label("Transfer rate from the 2nd peripheral compartment to central k31 (1/h)") # Table 3 'k31' 0.2 1/h (RSE 5.3%)
    lvmax_peripheral2 <- log(5105.4)
    label("Maximum Michaelis-Menten distribution rate into the 2nd peripheral compartment VmaxDistribution (ug/h)") # Table 3 'VmaxDistribution' 5105.4 ug/h (RSE 16.8%)
    lkm_peripheral2 <- log(145.1)
    label("Concentration at half-maximal distribution KmDistribution (ug/L)") # Table 3 'KmDistribution' 145.1 ug/L (RSE 20.0%)

    # IIV -- Table 3 reports Monolix CV% = 100*sqrt(exp(omega^2) - 1);
    # variances below are omega^2 = log(1 + (CV/100)^2).
    # Correlation(Cl, Vcentral) = 0.5, so cov = 0.5*sqrt(0.34058*0.71215).
    etalcl + etalvc ~ c(0.34058, 0.24624, 0.71215) # Table 3 'IIV Cl' 63.7%, 'IIV V' 101.9%, 'Correlation Vcentral-Cl' 0.5
    etalq ~ 1.28322 # Table 3 'IIV Q2' 161.5%
    etalvp ~ 0.50936 # Table 3 'IIV V2' 81.5%

    # Residual error: one combined model per assay / dose group
    addSd <- 6.7
    label("Additive residual error, therapeutic-dose studies (ug/L)") # Table 3 'sigma Add2' 6.7 ug/L (RSE 19.7%)
    propSd <- 0.5
    label("Proportional residual error, therapeutic-dose studies (fraction)") # Table 3 'sigma Prop2' 0.5 (RSE 2.4%)
    addSd_microdose <- fixed(1.0e-3)
    label("Additive residual error, microdose study (ug/L)") # Table 3 'sigma Add - (fixed)' 1.0e-03 ug/L
    propSd_microdose <- 0.1
    label("Proportional residual error, microdose study (fraction)") # Table 3 'sigma Prop' 0.1 (RSE 11.3%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    k31 <- exp(lk31)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    vmax_peripheral2 <- exp(lvmax_peripheral2)
    km_peripheral2 <- exp(lkm_peripheral2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Dose in mg and volumes in L give mg/L; x 1000 gives ug/L (Km units)
    Cc <- central / vc * 1000

    # Vmax values are in ug/h; / 1000 converts the fluxes to mg/h
    elim_mm <- vmax / 1000 * Cc / (km + Cc)
    dist_mm <- vmax_peripheral2 / 1000 * Cc / (km_peripheral2 + Cc)

    # Equations 1-2: saturable elimination and saturable distribution into the
    # 2nd peripheral compartment, which returns drug first order (k31)
    d/dt(central) <- -kel * central - elim_mm - k12 * central + k21 * peripheral1 - dist_mm + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- dist_mm - k31 * peripheral2

    addSd_eff <- addSd * (1 - STUDY_S65487_EXT001) + addSd_microdose * STUDY_S65487_EXT001
    propSd_eff <- propSd * (1 - STUDY_S65487_EXT001) + propSd_microdose * STUDY_S65487_EXT001
    Cc ~ add(addSd_eff) + prop(propSd_eff) + combined1()
  })
}
