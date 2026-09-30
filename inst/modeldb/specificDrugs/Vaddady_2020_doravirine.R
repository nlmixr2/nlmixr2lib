Vaddady_2020_doravirine <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral doravirine in healthy participants, treatment-naive adults with HIV-1, and virologically suppressed adults with HIV-1 switching to doravirine/lamivudine/tenofovir disoproxil fumarate (DRIVE-SHIFT immediate-switch group); linear age effect on CL/F, linear weight and healthy-versus-HIV-1 effects on V/F, and dose-band relative bioavailability"
  reference <- paste(
    "Vaddady P, Kandala B, Yee KL. Population Pharmacokinetic and Pharmacodynamic",
    "Analysis To Evaluate a Switch to Doravirine/Lamivudine/Tenofovir Disoproxil",
    "Fumarate in People Living with HIV-1. Antimicrob Agents Chemother.",
    "2020;64(11):e00590-20. doi:10.1128/AAC.00590-20.",
    "Covariate-equation forms and the residual-error structure are those of the",
    "predecessor model it re-estimates: Yee KL, Ouerdani A, Claussen A, de Greef R,",
    "Wenning L. Population Pharmacokinetics of Doravirine and Exposure-Response",
    "Analysis in Individuals with HIV-1. Antimicrob Agents Chemother.",
    "2019;63(4):e02502-18. doi:10.1128/AAC.02502-18."
  )
  vignette <- "Vaddady_2020_doravirine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL/F centred at 34 years (the median age of the Yee 2019 analysis data set): CL = TVCL * (1 + theta * (AGE - 34)).",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on V/F centred at 75 kg: V = TVV * (1 + theta * (WT - 75)).",
      source_name = "Weight"
    ),
    DIS_HEALTHY = list(
      description = "Healthy participant indicator (1 = healthy volunteer, 0 = person living with HIV-1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (people living with HIV-1, treatment-naive or virologically suppressed switch participants)",
      notes = "Multiplicative (1 + theta * DIS_HEALTHY) factor on V/F. Yee 2019 Table 2 footnote: 'flag is equal to 1 for healthy volunteers and 0 for HIV-infected individuals (the most common population)'; same orientation as the canonical, no transformation.",
      source_name = "flag (Subject status)"
    ),
    DOSE_DORAVIRINE_MG = list(
      description = "Administered doravirine dose amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Selects the relative-bioavailability dose band: < 30 mg, 30-120 mg (reference, F1 = 1), > 120 mg. The analysis data covered 6-200 mg; the clinical dose is 100 mg (reference band).",
      source_name = "dose (F1 dose band)"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2b study record indicator (1 = phase 2b, 0 = other phase)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 when STUDY_PHASE3 is also 0)",
      notes = "Changes the residual error only: phase 2b and phase 3 records share one log-scale SD, phase 1 records use the two time-after-dose-split SDs. For simulating a phase 2b/3-like sparse-sampled cohort set STUDY_PHASE2 = 0 and STUDY_PHASE3 = 1.",
      source_name = "Phase 2b/3"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study record indicator (1 = phase 3, 0 = other phase)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 when STUDY_PHASE2 is also 0)",
      notes = "Residual-error stratum only (see STUDY_PHASE2). P018 (DRIVE-FORWARD), P021 (DRIVE-AHEAD) and P024 (DRIVE-SHIFT) are phase 3.",
      source_name = "Phase 2b/3"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "doravirine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "doravirine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1743L,
    n_studies = 24L,
    age_median = "34 years (centring value; Yee 2019 analysis data set)",
    weight_median = "75 kg (centring value)",
    sex_female_pct = 19.2,
    race_ethnicity = c(White = 65.2, Black = 21.6, Asian = 6.7, Multiracial = 5.1, Other = 1.4),
    disease_state = "341 healthy participants (phase 1), 959 treatment-naive adults with HIV-1 (phase 1b/2b/3) and 443 virologically suppressed adults with HIV-1 switching to doravirine/lamivudine/tenofovir disoproxil fumarate (DRIVE-SHIFT immediate-switch group)",
    dose_range = "6-200 mg oral doravirine (single and multiple dose) in phase 1; 25-200 mg once daily in phase 2b; 100 mg once daily in phase 3",
    regions = "Multinational",
    notes = paste(
      "Vaddady 2020 main text: 341 healthy + 959 treatment-naive HIV-1 participants of the",
      "original Yee 2019 data set plus 443 DRIVE-SHIFT immediate-switch-group participants",
      "(1,402 with HIV-1). Studies: 20 phase 1 trials, phase 2b P007, and phase 3 P018, P021,",
      "P024. Sex and race percentages are from Yee 2019 Table 1 (the 1,300 original",
      "participants); Vaddady 2020 does not tabulate demographics for the combined set."
    )
  )

  # Three residual-error magnitudes (phase 1 <= 0.5 h postdose, phase 1 > 0.5 h
  # postdose, phase 2b/3), combined into expSdCc inside model(). Declared here so
  # checkModelConventions() does not read the stratum suffixes as deviant
  # residual-error names.
  paper_specific_residual_sds <- c("expSdP1Early", "expSdP1Late", "expSdP23")

  ini({
    lka <- log(1.42); label("Absorption rate constant Ka (1/h)") # Vaddady 2020 Table S1: Ka = 1.42 1/h (RSE 4.45%)
    lcl <- log(6.21); label("Apparent clearance CL/F for a 34-year-old (L/h)") # Vaddady 2020 Table S1: CL/F = 6.21 L/h (RSE 1.19%)
    lvc <- log(159); label("Apparent volume V/F for a 75-kg person living with HIV-1 (L)") # Vaddady 2020 Table S1: V/F = 159 L (RSE 2.72%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 for 30-120 mg doses (reference)") # Vaddady 2020 Table S1: 'F1 30-120 mg (reference)' = 1, no RSE

    e_dose_lt30_fdepot <- 1.20; label("Relative bioavailability F1 for doses < 30 mg vs 30-120 mg (ratio)") # Vaddady 2020 Table S1: F1 <30 mg = 1.20 (RSE 5.72%)
    e_dose_gt120_fdepot <- 0.882; label("Relative bioavailability F1 for doses > 120 mg vs 30-120 mg (ratio)") # Vaddady 2020 Table S1: F1 >120 mg = 0.882 (RSE 9.10%)
    e_age_cl <- -0.00540; label("Linear effect of age on CL/F, per year from 34 years (1/year)") # Vaddady 2020 Table S1: Age on CL = -0.00540 (RSE 13.0%)
    e_wt_vc <- 0.00788; label("Linear effect of body weight on V/F, per kg from 75 kg (1/kg)") # Vaddady 2020 Table S1: Weight on V = 0.00788 (RSE 12.9%)
    e_dis_healthy_vc <- -0.205; label("Fractional change in V/F for healthy participants vs people living with HIV-1 (fraction)") # Vaddady 2020 Table S1: Subject status on V = -0.205 (RSE 12.1%)

    etalcl ~ 0.104 # Vaddady 2020 Table S1: omega2 CL/F = 0.104 (33.1% CV)
    etalvc ~ 0.098 # Vaddady 2020 Table S1: omega2 V/F = 0.098 (32.1% CV)

    expSdP1Early <- 0.224; label("Log-scale additive residual SD, phase 1 records <= 0.5 h postdose") # Vaddady 2020 Table S1: SD Phase 1 <=0.5 h postdose = 0.224
    expSdP1Late <- 1.25; label("Log-scale additive residual SD, phase 1 records > 0.5 h postdose") # Vaddady 2020 Table S1: SD Phase 1 >0.5 h postdose = 1.25
    expSdP23 <- 0.504; label("Log-scale additive residual SD, phase 2b/3 records") # Vaddady 2020 Table S1: SD Phase 2b/3 = 0.504
  })

  model({
    # Dose-band relative bioavailability (Yee 2019 / Vaddady 2020 Table S1:
    # < 30 mg, 30-120 mg reference, > 120 mg)
    doseLt30 <- DOSE_DORAVIRINE_MG < 30
    doseGt120 <- DOSE_DORAVIRINE_MG > 120
    fdose <- 1 + doseLt30 * (e_dose_lt30_fdepot - 1) + doseGt120 * (e_dose_gt120_fdepot - 1)

    # Individual parameters (Vaddady 2020 Table S1 footnote; Yee 2019 Table 2
    # footnote for the healthy-status factor)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 + e_age_cl * (AGE - 34))
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - 75)) * (1 + e_dis_healthy_vc * DIS_HEALTHY)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot) * fdose

    Cc <- central / vc

    # Residual error: additive on log-transformed concentrations (Yee 2019),
    # with separate SDs for phase 1 <= 0.5 h postdose, phase 1 > 0.5 h
    # postdose, and phase 2b/3.
    tsld <- tad()
    isEarly <- tsld <= 0.5
    isPh23 <- STUDY_PHASE2 + STUDY_PHASE3
    expSdCc <- (1 - isPh23) * (isEarly * expSdP1Early + (1 - isEarly) * expSdP1Late) +
      isPh23 * expSdP23

    Cc ~ lnorm(expSdCc)
  })
}
