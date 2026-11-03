Avachat_2025_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in critically ill, predominantly paediatric Asian patients (Taiwan; age 0.1-36.6 years, weight 0.7-93 kg) receiving therapeutic drug monitoring (Avachat 2025). Clearance scales allometrically with body weight (exponent fixed at 0.75) and volume linearly (exponent fixed at 1), both normalised to the 12.6 kg cohort median; clearance is further reduced by exp(-0.81) (about 56%) in patients on extracorporeal membrane oxygenation (ECMO) and by exp(-0.49) (about 39%) in patients whose vancomycin indication was a milder infection (cellulitis, surgical prophylaxis, neutropenic fever or other) rather than a severe one (pneumonia, bacteremia, osteomyelitis, meningitis or deep tissue infection). Exponential between-subject variability on CL and V; proportional residual error."
  reference <- "Avachat C, Hung CC, Birnbaum AK, Healy DP, Sherwin CMT, Lin YC. Pharmacokinetic effect of disease severity and use of extracorporeal membrane oxygenation in critically ill Asian patients receiving vancomycin. Front Pharmacol. 2025;16:1506793. doi:10.3389/fphar.2025.1506793. PMID 40078297; PMCID PMC11897554."
  vignette <- "Avachat_2025_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is given as an intravenous infusion, so the
  # dose enters `central` directly; Avachat 2025 Section 2.2 models "vancomycin
  # plasma concentrations".
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL as (WT/12.6)^0.75 and V as (WT/12.6)^1. The final-model equations in Avachat 2025 Section 3.3 print the ratio as WT / median(WT), and the text immediately below states 'a typical individual is someone weighing 12.6 kg, which is the median weight of the study population' (Table 1: median 12.6 kg, range 0.7-93 kg). The allometric exponents were fixed a priori (Section 2.2). The source does not say whether weight was time-varying across the sampling period.",
      source_name = "WT"
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator (1 = on ECMO, 0 = not on ECMO)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving ECMO)",
      notes = "Source column ECMO, same orientation ('ECMO = = 1 indicating patients on ECMO', Avachat 2025 Section 3.3). 12 of 90 patients (13.3%) were on ECMO (Table 1). Enters CL as the exponential multiplier exp(dCldECMO * ECMO) with dCldECMO = -0.81, a 55.5% reduction (the paper rounds to 56%). No ECMO effect on V was retained. The source does not resolve within-subject cannulation / decannulation timing, so the indicator is treated here as subject-level, as the paper reports it.",
      source_name = "ECMO"
    ),
    DIS_INFECT_MILD = list(
      description = "Mild-infection indication category for vancomycin therapy (1 = cellulitis, surgical prophylaxis, neutropenic fever or other; 0 = pneumonia, bacteremia, osteomyelitis, meningitis or deep tissue infection)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (severe infection: pneumonia, bacteremia, osteomyelitis, meningitis or deep tissue infection)",
      notes = "Source column DIAG, same orientation ('DIAG = = 1 indicating a patient belongs to mild infection category', Avachat 2025 Section 3.3). The nine admission diagnoses were grouped by the authors into two categories by the severity of the bacterial infection for which vancomycin was used (Section 2.1): severe = pneumonia (20), bacteremia (6), osteomyelitis (7), meningitis (4), deep tissue infection (1), 38 patients in total; mild = cellulitis (6), surgical prevention (35), neutropenic fever (4), others (7), 52 patients in total (Table 1). Enters CL as exp(dCldDiag * DIAG) with dCldDiag = -0.49, a 38.7% reduction (the paper rounds to 39%). Subject-level (diagnosis at admission).",
      source_name = "DIAG"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 90L,
    n_studies = 1L,
    n_samples = 263L,
    age_range = "0.1-36.6 years (median 3.7); 90% of the cohort were paediatric",
    weight_range = "0.7-93 kg (median 12.6)",
    sex_female_pct = 47.8,
    race_ethnicity = "Asian (Taiwanese; single centre)",
    disease_state = "Critically ill patients receiving IV vancomycin for treatment or prevention of MRSA or other bacterial infections with therapeutic drug monitoring. Indication: 4 MRSA, 8 severe MRSA, 27 non-MRSA, 50 empirical, 1 other. Severe-infection diagnoses (pneumonia 20, bacteremia 6, osteomyelitis 7, meningitis 4, deep tissue infection 1) versus mild (cellulitis 6, surgical prevention 35, neutropenic fever 4, others 7). 12 (13.3%) on ECMO; 33 (36.7%) underwent surgery (Avachat 2025 Table 1).",
    renal_function = "eGFR by the bedside Schwartz equation, median 70.4 (range 2.4-150.6); Table 1 labels the unit mL/min, the Limitations section mL/min/1.73 m^2. One patient on regular haemodialysis and one on CRRT. eGFR was not retained as a covariate.",
    dose_range = "IV infusion at 8-1500 mg/h (Section 3.1); individual doses and dosing intervals are not reported.",
    regions = "Taiwan (Kaohsiung Veterans General Hospital, single centre)",
    sampling = "Retrospective TDM data, 1 January 2014 - 31 December 2016: peak samples drawn > 1 h after the end of infusion and trough samples within 30 min before the next dose, predominantly troughs; 263 concentrations from 90 patients. Concentrations below the LLOD were excluded and those below the LLOQ were set to LLOQ/2.",
    notes = "Phoenix NLME 8.2, FOCE-ELS. Covariates screened by stepwise forward inclusion (p < 0.05) and backward elimination (p < 0.01): infection severity, ECMO status and surgery status; surgery on CL was removed at backward elimination (Table 2). eGFR was not significant. Model stability checked by a 1000-sample nonparametric bootstrap (Table 4)."
  )

  ini({
    # Structural parameters: Avachat 2025 Table 4, 'Final model' 'Estimate'
    # column. Typical subject: 12.6 kg, not on ECMO, severe infection.
    lcl <- log(1.03); label("Clearance for a 12.6 kg patient, not on ECMO, severe infection (L/h)")  # Table 4 'tvCL (L/h)' = 1.03 (SE 0.099, 95% CI 0.84-1.23)
    lvc <- log(11.61); label("Volume of distribution for a 12.6 kg patient (L)")                     # Table 4 'tvV (L)' = 11.61 (SE 1.22, 95% CI 9.2-14.03)

    # Allometric exponents fixed a priori (Section 2.2: 'fixed at 0.75 and 1,
    # respectively'); Table 4 rows 'dCldWT (fixed)' and 'dVdWT (fixed)'.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/12.6) for CL (unitless)")  # Table 4 'dCldWT (fixed)' = 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on (WT/12.6) for V (unitless)")      # Table 4 'dVdWT (fixed)' = 1

    # Categorical effects on CL, entering as exp(theta * indicator) per the
    # Section 3.3 final-model equation and the Table 3 worked example.
    e_ecmo_status_cl <- -0.81; label("Exponential coefficient of ECMO_STATUS on CL (unitless)")    # Table 4 'dCldECMO' = -0.81 (SE 0.17, 95% CI -1.14 to -0.47)
    e_infect_mild_cl <- -0.49; label("Exponential coefficient of DIS_INFECT_MILD on CL (unitless)") # Table 4 'dCldDiag' = -0.49 (SE 0.18, 95% CI -0.84 to -0.14)

    # Between-subject variability. Phoenix NLME reports the Omega matrix as
    # variances; read as variances, lognormal CV = 54% (CL) and 52% (V).
    etalcl ~ 0.26 # Table 4 'eta Cl' = 0.26 (SE 0.08, shrinkage 0.27)
    etalvc ~ 0.24 # Table 4 'eta V' = 0.24 (SE 0.07, shrinkage 0.52)

    # Proportional residual error (Section 3.2 'proportional error'; Section 5
    # 'A relatively high residual error of 44%').
    propSd <- 0.44; label("Proportional residual error (fraction)")  # Table 4 'Stdev0' = 0.44 (SE 0.03, 95% CI 0.39-0.5)
  })
  model({
    # Avachat 2025 Section 3.3:
    #   Cl = tvCl * (WT/median(WT))^dCldWT * exp(dCldECMO*(ECMO==1))
    #        * exp(dCldDIAG*(DIAG==1)) * exp(etaCl)
    #   V  = tvV * (WT/median(WT))^dVdWT * exp(etaV)
    # with median(WT) = 12.6 kg.
    cl <- exp(lcl + etalcl) * (WT / 12.6)^e_wt_cl * exp(e_ecmo_status_cl * ECMO_STATUS) * exp(e_infect_mild_cl * DIS_INFECT_MILD)
    vc <- exp(lvc + etalvc) * (WT / 12.6)^e_wt_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L, matching the Figure 3 axis.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
