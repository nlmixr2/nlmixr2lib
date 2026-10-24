Wang_2022_venlafaxine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral venlafaxine (VEN)",
    "and its active metabolite O-desmethylvenlafaxine (ODV, desvenlafaxine)",
    "in Chinese healthy male volunteers and psychiatric patients (Wang 2022).",
    "One-compartment disposition for each analyte with first-order absorption",
    "(ka fixed at 0.63 1/h). The absorbed dose splits at absorption: a",
    "fraction fm_odv is converted pre-systemically (first pass) and enters",
    "the ODV compartment directly, and the remainder enters the VEN",
    "compartment. All systemic VEN clearance forms ODV (k23 = CL/V), and ODV",
    "is cleared with CL_odv. Clearances and volumes are apparent (relative to",
    "the unknown oral bioavailability). The model works on a MOLAR basis:",
    "doses are in umol of venlafaxine and both concentrations are in umol/L,",
    "as in the source figures. Patient status lowers VEN CL/F by 61.7 percent",
    "relative to healthy volunteers; concomitant amisulpride lowers VEN CL/F",
    "by 39.2 percent and raises ODV CL/F by 59.3 percent. Proportional",
    "residual error on each analyte."
  )
  reference <- paste(
    "Wang Z, Li L, Huang S, Wang X, Liu S, Li X, Kong W, Ni X, Zhang M,",
    "Huang S, Tan Y, Wen Y, Shang D. Joint population pharmacokinetic",
    "modeling of venlafaxine and O-desmethyl venlafaxine in healthy",
    "volunteers and patients to evaluate the impact of morbidity and",
    "concomitant medication. Front Pharmacol. 2022;13:978202.",
    "doi:10.3389/fphar.2022.978202"
  )
  vignette <- "Wang_2022_venlafaxine"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    depot = list(analyte = "venlafaxine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venlafaxine", units = "umol", specimen = "plasma", verified = TRUE),
    central_odv = list(analyte = "O-desmethylvenlafaxine", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator, 1 = healthy volunteer, 0 = psychiatric patient.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient) in the register; the source paper's reference is the healthy volunteer, so model() applies the 'morbid state' coefficient to the complement (1 - DIS_HEALTHY).",
      notes = "Wang 2022 Table 2 'theta morbid state on CL/F' = 0.617; Results 3.2: 'Morbidity state was associated with a decrease in the clearance of VEN by 61.7%', and the population-predicted VEN CL/F in patients is quoted as 31 L/h (= 80.9 * (1 - 0.617) = 31.0). Healthy volunteers are the 24 men of the bioequivalence Study 1; patients are the 127 psychiatric inpatients of the TDM Study 2 (diagnoses not itemised).",
      source_name = "morbid state (health status)"
    ),
    CONMED_AMISULPRIDE = list(
      description = "Concomitant amisulpride indicator, 1 = co-medicated with amisulpride, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no amisulpride).",
      notes = "Wang 2022 Table 2 'theta amisulpride on CL/F' = 0.392 and 'theta amisulpride on CLM/F' = 0.593; Results 3.2: VEN and ODV clearance with amisulpride 'was 39.2% lower ... and 59.3% higher' than without. 5.47 percent of Study 2 patients received amisulpride (Table 1); no healthy volunteer did. Combined multiplicatively with the patient effect (not stated in the text; the multiplicative form reproduces the Figure 3 median VEN CL of about 19 L/h in amisulpride-treated patients, 80.9 * 0.383 * 0.608 = 18.8).",
      source_name = "amisulpride"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Screened; VEN CL was higher in males in univariate analysis but the effect was not retained in the final model (Discussion)."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened; the age effect on VEN CL was small (p < 0.05) and not retained (Discussion)."
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 2.2); not retained."
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Effect on CL/F and CLM/F (p < 0.05) at forward inclusion; dropped at backward elimination (Discussion)."
    ),
    CONMED_QUETIAPINE = list(
      description = "Concomitant quetiapine indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Effect on V/F (p < 0.01) at forward inclusion; dropped at backward elimination (Discussion)."
    ),
    CONMED_CLOZAPINE = list(
      description = "Concomitant clozapine indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Effect on CLM/F (p < 0.01) at forward inclusion; dropped at backward elimination (Discussion)."
    ),
    CONMED_RISPERIDONE = list(
      description = "Concomitant risperidone indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Effect on V/F (p < 0.05) at forward inclusion; dropped at backward elimination (Discussion)."
    ),
    CONMED_OLANZAPINE = list(
      description = "Concomitant olanzapine indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 2.2); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 151L,
    n_studies = 2L,
    age_range = "Study 1: 18-27 years (mean 23.1); Study 2: 14-86 years (mean 37.9)",
    weight_range = "Study 1: 54-80 kg (mean 62.5); Study 2: 38-93 kg (mean 62.0)",
    sex_female_pct = 46.4,
    race_ethnicity = "Chinese (single-centre studies in Guangzhou, China)",
    disease_state = "Study 1: 24 healthy male volunteers (single-dose crossover bioequivalence study, VEN measured only). Study 2: 127 psychiatric inpatients on routine VEN therapy (retrospective TDM, VEN and ODV measured, mostly troughs).",
    dose_range = "Study 1: single oral 50 mg venlafaxine hydrochloride. Study 2: 25-300 mg/day venlafaxine (mean 132.9), once to three times daily, 83.6 percent sustained-release and 16.4 percent immediate-release.",
    regions = "China",
    n_observations = "664 VEN concentrations (Study 1) and 804 VEN + ODV concentrations (Study 2).",
    notes = "Wang 2022 Table 1. Female percentage is the pooled value: 55.12 percent of the 127 Study 2 patients (70 women) and none of the 24 Study 1 volunteers, 70 / 151 = 46.4 percent. Co-medications in Study 2: valproic acid 18.2 percent, quetiapine 23.4 percent, clozapine 15.7 percent, olanzapine 24.5 percent, risperidone 13.2 percent, amisulpride 5.5 percent."
  )

  ini({
    # Wang 2022 Table 2 (final joint PPK model). Healthy volunteer without
    # amisulpride is the reference for all typical values.
    lka <- fixed(log(0.63))
    label("Absorption rate constant (1/h)") # Table 2 'Ka' = 0.63 1/h, no RSE; Results 3.2 'fixed at 0.63 h-1'
    lcl <- log(80.9)
    label("Apparent VEN clearance CL/F, healthy volunteer (L/h)") # Table 2 'CL/F' = 80.9 L/h (RSE 9.7%)
    lvc <- log(628)
    label("Apparent VEN volume of distribution V/F (L)") # Table 2 'V/F' = 628 L (RSE 5.8%)
    lcl_odv <- log(22.1)
    label("Apparent ODV clearance CLM/F (L/h)") # Table 2 'CLM/F' = 22.1 L/h (RSE 6.6%)
    lvc_odv <- log(238)
    label("Apparent ODV volume of distribution VM/F (L)") # Table 2 'VM/F' = 238 L (RSE 33.1%)
    fm_odv <- 0.048
    label("Fraction of the absorbed VEN dose converted to ODV pre-systemically, FP (fraction)") # Table 2 'FP' = 0.048 (RSE 18.3%)

    e_patient_cl <- 0.617
    label("Fractional decrease in VEN CL/F in patients vs healthy volunteers (fraction)") # Table 2 'theta morbid state on CL/F' = 0.617 (RSE 5.9%)
    e_conmed_amisulpride_cl <- 0.392
    label("Fractional decrease in VEN CL/F with concomitant amisulpride (fraction)") # Table 2 'theta amisulpride on CL/F' = 0.392 (RSE 17.8%)
    e_conmed_amisulpride_cl_odv <- 0.593
    label("Fractional increase in ODV CL/F with concomitant amisulpride (fraction)") # Table 2 'theta amisulpride on CLM/F' = 0.593 (RSE 28.7%)

    # Table 2 heads the IIV block 'Interindividual variability, %CV' but prints
    # fractions. They are the NONMEM OMEGA variances: the spread of the
    # individual clearances in Figure 3 (healthy VEN CL IQR about 57-105 L/h,
    # log-SD about 0.45) matches sqrt(0.219) = 0.47, not 0.219. See vignette.
    etalcl ~ 0.219 # Table 2 IIV 'CL/F' = 0.219 (omega^2; RSE 13.4%)
    etalvc ~ 0.106 # Table 2 IIV 'V/F' = 0.106 (omega^2; RSE 35.2%)
    etalcl_odv ~ 0.156 # Table 2 IIV 'CLM/F' = 0.156 (omega^2; RSE 19.4%)
    etalvc_odv ~ 1.38 # Table 2 IIV 'VM/F' = 1.38 (omega^2; RSE 39.3%)

    # Residual rows of the same column are read on the same (variance) scale:
    # SD = sqrt(SIGMA). See vignette.
    propSd <- sqrt(0.123)
    label("Proportional residual error SD, VEN (fraction)") # Table 2 'Proportional error on VEN' = 0.123 (sigma^2; RSE 8.5%)
    propSd_odv <- sqrt(0.101)
    label("Proportional residual error SD, ODV (fraction)") # Table 2 'Proportional error on ODV' = 0.101 (sigma^2; RSE 13.6%)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 - e_patient_cl * (1 - DIS_HEALTHY)) *
      (1 - e_conmed_amisulpride_cl * CONMED_AMISULPRIDE)
    vc <- exp(lvc + etalvc)
    cl_odv <- exp(lcl_odv + etalcl_odv) *
      (1 + e_conmed_amisulpride_cl_odv * CONMED_AMISULPRIDE)
    vc_odv <- exp(lvc_odv + etalvc_odv)

    # Figure 1: K23 = CL/V carries all systemic VEN elimination into ODV;
    # K30 = CLM/VM eliminates ODV.
    kel <- cl / vc
    kel_odv <- cl_odv / vc_odv

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * (1 - fm_odv) * depot - kel * central
    d/dt(central_odv) <- ka * fm_odv * depot + kel * central - kel_odv * central_odv

    Cc <- central / vc
    Cc_odv <- central_odv / vc_odv

    Cc ~ prop(propSd)
    Cc_odv ~ prop(propSd_odv)
  })
}
