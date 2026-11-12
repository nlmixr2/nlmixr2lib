Lanoiselee_2024_unfractionatedHeparin <- function() {
  description <- "One-compartment population PK model for intravenous unfractionated heparin (UFH) in adults on veno-arterial extracorporeal membrane oxygenation (VA-ECMO), with exposure measured as plasma anti-Xa activity (Lanoiselee 2024). Clearance rises with C-reactive protein and falls with serum creatinine; volume scales linearly with body weight, rises with C-reactive protein and is 1.91-fold larger when ECMO was instituted for post-cardiotomy shock. Combined additive and proportional residual error."
  reference <- "Lanoiselee J, Mourer J, Jungling M, Molliex S, Thellier L, Tabareau J, Jeanpierre E, Robin E, Susen S, Tavernier B, Vincentelli A, Ollier E, Moussa MD. Heparin Dosing Regimen Optimization in Veno-Arterial Extracorporeal Membrane Oxygenation: A Pharmacokinetic Analysis. Pharmaceutics. 2024;16(6):770. doi:10.3390/pharmaceutics16060770. PMC11207075."
  vignette <- "Lanoiselee_2024_unfractionatedHeparin"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  # UFH was given as intravenous bolus and continuous infusion; exposure was
  # measured as plasma anti-Xa activity by a chromogenic assay (Section 2.5),
  # so the central amount is UFH in anti-Xa international units.
  compartmentData <- list(
    central = list(analyte = "unfractionatedHeparin", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description = "Serum creatinine (time-varying)",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: collected regularly throughout UFH administration (Section 2.6). Enters CL as (CREAT/115)^-0.237 (Table 2); 115 umol/L is the development-cohort median (114.9 umol/L, Table 1). Patients on continuous renal replacement therapy (45% of the development cohort) are included, so the value reflects native plus CRRT clearance (Discussion).",
      source_name = "SCr"
    ),
    CRP = list(
      description = "C-reactive protein, standard assay (time-varying)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: collected regularly throughout UFH administration (Section 2.6). Enters CL as (CRP/100)^0.258 and V as (CRP/110)^0.191 (Table 2); the two different centring values are as printed in Table 2. Development-cohort median 113 (range 3-547) mg/L (Table 1).",
      source_name = "CRP"
    ),
    WT = list(
      description = "Baseline total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) total body weight (Section 2.6). Enters V as (WT/80)^1 with the exponent fixed to 1 (Table 2). Development-cohort median 75 (range 54-122) kg (Table 1).",
      source_name = "BW"
    ),
    ECMO_IND_POSTCARDIOTOMY = list(
      description = "ECMO indication: 1 = post-cardiotomy shock, 0 = medical indication",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (medical cardiogenic shock)",
      notes = "Table 2 footnote: 'Indic being equal to 1 for post-cardiotomy indication and 0 for medical indication'. Time-fixed. Post-cardiotomy low cardiac output syndrome in 26 of 64 development patients (Table 1). Enters V as exp(0.647 * ECMO_IND_POSTCARDIOTOMY) (Section 3.3 loading-dose equation).",
      source_name = "Indic"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as a time-fixed covariate in the stepwise search (Section 2.6) and not retained (Section 3.2)."
    ),
    FIB = list(
      description = "Plasma fibrinogen (time-varying)",
      units = "g/L",
      type = "continuous",
      notes = "Tested as a time-varying inflammation marker (Section 2.6) and not retained (Section 3.2)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy indicator (time-varying)",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a time-varying covariate (Section 2.6) and not retained (Section 3.2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 74L,
    n_studies = 1L,
    n_observations = 1703L,
    age_range = "52 +/- 13 years (mean +/- SD; Table 1)",
    weight_range = "54-122 kg",
    weight_median = "75 kg",
    sex_female_pct = 28.4,
    race_ethnicity = "Not reported (single-centre French cohort)",
    disease_state = "Adults on veno-arterial ECMO for refractory cardiogenic shock or cardiac arrest: post-cardiotomy low cardiac output syndrome 39%, myocardial infarction 32%, myocarditis 5%, acute-on-chronic heart disease 8%, other 15% (Table 1). Peripheral ECMO in 95%; CRRT during ECMO in 46%.",
    renal_function = "Serum creatinine median 114.9 (range 26.5-530.4) umol/L (Table 1).",
    dose_range = "100 IU/kg UFH loading dose before cannulation (unless coagulopathy, post-cardiotomy or prior UFH), then continuous IV infusion titrated to anti-Xa 0.2-0.7 IU/mL; mean dose 11.39 IU/kg/h (Section 2.4, Table 1).",
    regions = "France (Lille University Hospital)",
    notes = "Retrospective single-centre study, January 2020 - June 2021. 64 patients (1487 anti-Xa samples) in the development cohort and 10 (216 samples) in the validation cohort. Anti-Xa LLOQ 0.1 IU/mL and ULOQ 1.8 IU/mL; 25% of samples were left-censored and handled as censored data. Monolix 2023R1, SAEM."
  )

  ini({
    lcl <- log(3.41)
    label("Clearance at CREAT 115 umol/L and CRP 100 mg/L (L/h)") # Table 2: theta1 = 3.41 L/h (RSE 7.04%)
    lvc <- log(8.65)
    label("Volume at WT 80 kg, CRP 110 mg/L, medical indication (L)") # Table 2: theta4 = 8.65 L (RSE 15.7%)

    e_creat_cl <- -0.237
    label("Power exponent of serum creatinine on CL (unitless)") # Table 2: theta2 = -0.237 (RSE 0.18%)
    e_crp_cl <- 0.258
    label("Power exponent of CRP on CL (unitless)") # Table 2: theta3 = 0.258 (RSE 0.66%)
    e_ecmo_ind_postcardiotomy_vc <- 0.647
    label("Exponential effect of post-cardiotomy ECMO indication on V (unitless)") # Table 2: theta5 = 0.647 (no RSE printed)
    e_wt_vc <- fixed(1)
    label("Power exponent of body weight on V (unitless)") # Table 2: theta6 = 1 (fixed)
    e_crp_vc <- 0.191
    label("Power exponent of CRP on V (unitless)") # Table 2: theta7 = 0.191 (RSE 0.45%)

    # Table 2 reports the IIV as SDs of the log-normal random effects, so the
    # variances are their squares.
    etalcl ~ 0.2704 # Table 2: Omega CL (SD) = 0.52 (RSE 10.3%); variance 0.52^2
    etalvc ~ 0.5625 # Table 2: Omega V (SD) = 0.75 (RSE 13.2%); variance 0.75^2

    addSd <- 0.07
    label("Additive residual error (IU/mL)") # Table 2: additive residual variability (SD) = 0.07 (RSE 5.42%)
    propSd <- 0.41
    label("Proportional residual error (fraction)") # Table 2: proportional residual variability (SD) = 0.41 (RSE 3.62%)
  })
  model({
    # Covariate model, Table 2 and Section 3.3. The ECMO-indication term is
    # exp(theta5 * Indic), as in the Section 3.3 loading-dose equation; Table 2
    # prints it as theta4^(Indic x theta5), which would give a 1 L volume for
    # every medical patient.
    cl <- exp(lcl + etalcl) * (CREAT / 115)^e_creat_cl * (CRP / 100)^e_crp_cl
    vc <- exp(lvc + etalvc + e_ecmo_ind_postcardiotomy_vc * ECMO_IND_POSTCARDIOTOMY) *
      (WT / 80)^e_wt_vc * (CRP / 110)^e_crp_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in IU and volume in L give IU/L; anti-Xa is reported in IU/mL.
    Cc <- central / vc / 1000
    # Monolix combined error y = F + (a + b*F) * eps, i.e. nlmixr2 combined1.
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
