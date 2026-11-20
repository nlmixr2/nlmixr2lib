Burlot_2026_S65487_pkpd <- function() {
  description <- "Three-compartment population PK model for intravenous S65487 (Bcl2 inhibitor) linked to an indirect-response model of inhibited Bcl2/Bim complex formation in adults with hematological malignancies (PhI-II PK/PD model; Burlot 2026)"
  reference <- "Burlot C, Riglet F, Romagnoli M, Leconte A, Delmas A, Schmitt A, Fouliard S. PK and PK/PD Modeling of Bcl2 Inhibitor S65487 in Patients With AML and Investigation of Nonlinearity With Microdosing. CPT Pharmacometrics Syst Pharmacol. 2026;15(7):e70288. doi:10.1002/psp4.70288"
  vignette <- "Burlot_2026_S65487"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    STUDY_S65487_003 = list(
      description = "Indicator for the S65487 phase I/II study 003 (S65487 plus azacitidine in previously untreated AML); 0 = phase I study 002 (S65487 monotherapy in relapsed/refractory hematological malignancies)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I study 002)",
      notes = "Burlot 2026 Table 2 'beta-IC50 STUDY003' = -1.6 on log(IC50); the only covariate retained in the PhI-II PK/PD model. Subject-level (time-fixed).",
      source_name = "STUDY (003)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "S65487", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S65487", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "S65487", units = "mg", specimen = "tissue", verified = TRUE),
    bcl2_bim = list(
      analyte = "Bcl2/Bim complex (normalized by total protein)",
      units = "ug/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # Turnover state of the Bcl2/Bim complex measured in PBMC and normalized by
  # total protein; a drug-specific PD biomarker with no register entry.
  paper_specific_compartments <- c("bcl2_bim")

  population <- list(
    species = "human",
    n_subjects = 101,
    n_studies = 2,
    age_range = "28-87 years",
    age_median = "67 years (phase I), 75 years (phase I/II)",
    weight_range = "45.0-107.0 kg",
    weight_median = "71.7 kg (phase I), 74.4 kg (phase I/II)",
    sex_female_pct = 37.6,
    disease_state = "Phase I (study 002, n = 60): relapsed/refractory AML, non-Hodgkin lymphoma or multiple myeloma, S65487 monotherapy. Phase I/II (study 003, n = 41): previously untreated AML, S65487 in combination with azacitidine.",
    dose_range = "25-1200 mg IV infusion; phase I once weekly (with ramp-up) or days 1, 3, 5 then weekly over 21-day cycles; phase I/II weekly, days 1, 3, 5, 8 or days 1-5, 8, 9 over 28-day cycles",
    regions = "Europe (EudraCT 2018-004170-97 and 2020-003061-19)",
    notes = "Burlot 2026 Table 1 (therapeutic-dose columns). 1947 S65487 concentrations and 456 Bcl2/Bim complex concentrations (from 87 patients). The microdose phase 0 study in healthy volunteers is not part of this model; see Burlot_2026_S65487_microdose_pk."
  )

  ini({
    # Structural PK -- Table 2 'PhI-II PK/PD model' column (Monolix 2023R1, SAEM)
    lcl <- log(12.2)
    label("Clearance CL (L/h)") # Table 2 'Cl' 12.2 L/h (RSE 5.4%)
    lvc <- log(7.6)
    label("Central volume Vcentral (L)") # Table 2 'Vcentral' 7.6 L (RSE 10.7%)
    lq <- log(1.5)
    label("Intercompartmental clearance to the 1st peripheral compartment Q2 (L/h)") # Table 2 'Q2' 1.5 L/h (RSE 13.1%)
    lvp <- log(16.4)
    label("Volume of the 1st peripheral compartment V2 (L)") # Table 2 'V2' 16.4 L (RSE 14.0%)
    lq2 <- log(15.5)
    label("Intercompartmental clearance to the 2nd peripheral compartment Q3 (L/h)") # Table 2 'Q3' 15.5 L/h (RSE 18.6%)
    lvp2 <- log(14.8)
    label("Volume of the 2nd peripheral compartment V3 (L)") # Table 2 'V3' 14.8 L (RSE 7.5%)

    # Indirect-response PD (inhibition of Bcl2/Bim complex formation, kin)
    lrbase <- log(2.1e-4)
    label("Baseline normalized Bcl2/Bim complex R0 (ug/L)") # Table 2 'R0' 2.1e-04 ug/L (RSE 12.7%)
    lkout <- log(0.4)
    label("Bcl2/Bim complex degradation rate kout (1/h)") # Table 2 'kout' 0.4 1/h (RSE 8.8%)
    limax <- fixed(log(1))
    label("Maximum fractional inhibition of complex formation Imax (unitless)") # Results 'Imax parameter was set to 1'
    lic50 <- log(140.9)
    label("S65487 concentration for half-maximal inhibition IC50, phase I study 002 (ug/L)") # Table 2 'IC50' 140.9 ug/L (RSE 45.6%)
    e_study_s65487_003_ic50 <- -1.6
    label("Study 003 (phase I/II) effect on log(IC50) (unitless)") # Table 2 'beta-IC50 STUDY003' -1.6 (RSE 39.8%)

    # IIV -- Table 2 reports Monolix CV% = 100*sqrt(exp(omega^2) - 1);
    # variances below are omega^2 = log(1 + (CV/100)^2).
    # Correlation(Cl, Vcentral) = 0.6, so cov = 0.6*sqrt(0.29782*0.88806).
    etalcl + etalvc ~ c(0.29782, 0.30857, 0.88806) # Table 2 'IIV CL' 58.9%, 'IIV Vc' 119.6%, 'Correlation V-Cl' 0.6
    etalq ~ 2.03689 # Table 2 'IIV Q2' 258.2%
    etalvp ~ 2.66336 # Table 2 'IIV V2' 365.3%
    etalq2 ~ 2.66387 # Table 2 'IIV Q3' 365.4%
    etalvp2 ~ 0.25421 # Table 2 'IIV V3' 53.8%
    etalrbase ~ 1.14341 # Table 2 'IIV R0' 146.2%
    etalkout ~ 0.59532 # Table 2 'IIV kout' 90.2%
    etalic50 ~ 5.26392 # Table 2 'IIV IC50' 1386.5%

    # Residual error
    addSd <- 2.9
    label("Additive residual error for S65487 (ug/L)") # Table 2 'sigma Add' 2.9 ug/L (RSE 8.7%)
    propSd <- 0.5
    label("Proportional residual error for S65487 (fraction)") # Table 2 'sigma Prop' 0.5 (RSE 2.0%)
    propSd_bcl2_bim <- 0.7
    label("Proportional residual error for normalized Bcl2/Bim complex (fraction)") # Table 2 'sigma Prop2' 0.7 (RSE 4.6%)
  })

  model({
    # Individual PK parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Individual PD parameters; Monolix log-normal covariate model on IC50
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    imax <- exp(limax)
    ic50 <- exp(lic50 + e_study_s65487_003_ic50 * STUDY_S65487_003 + etalic50)
    kin <- rbase * kout

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg and volumes in L give mg/L; x 1000 gives ug/L
    Cc <- central / vc * 1000

    # Indirect response: S65487 inhibits complex formation (kin)
    bcl2_bim(0) <- rbase
    d/dt(bcl2_bim) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * bcl2_bim

    Cc ~ add(addSd) + prop(propSd) + combined1()
    bcl2_bim ~ prop(propSd_bcl2_bim)
  })
}
