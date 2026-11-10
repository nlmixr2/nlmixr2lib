Huestis_2025_midomafetamine <- function() {
  description <- "One-compartment oral population PK model for midomafetamine (MDMA) with an absorption lag followed by sequential zero-order then first-order absorption, linked to a two-compartment model of its metabolite tenamfetamine (MDA, 3,4-methylenedioxyamphetamine) formed from a fixed 10% of MDMA clearance, with body-weight effects on MDMA CL/F and V/F and on MDA CL and Q, a sex effect on MDA peripheral volume, high-fat-meal and fasted effects on MDMA ka, and an empirical study effect on MDMA CL/F and V/F, in healthy adult volunteers (Huestis 2025)"
  reference <- paste(
    "Huestis MA, Smith WB, Leonowens C, Blanchard R, Viaccoz A, Spargo E,",
    "Miner NB, Yazar-Klosinski B. MDMA pharmacokinetics: A population and",
    "physiologically based pharmacokinetics model-informed analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):376-388.",
    "doi:10.1002/psp4.13282.",
    "Only the population PK model (Table 2) is encoded; the paper's Simcyp",
    "minimal-PBPK model with CYP2D6 mechanism-based inhibition is not.",
    sep = " "
  )
  vignette <- "Huestis_2025_midomafetamine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Estimated power effects on MDMA CL/F (exponent 1.03) and V2/F (0.885) and on MDA CLM (0.580) and QM (1.00), all normalised to 70 kg (Huestis 2025 Table 2; Methods Equation 3 names 70 kg as the standard body-weight normaliser and the Figure 2 reference individual weighs 70 kg). Overall analysis-population median 73.8 kg, range 43.1-113.1 kg (Supplementary Table S4).",
      source_name = "body weight"
    ),
    SEXF = list(
      description = "Biological sex, female indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Proportional effect on the MDA peripheral volume V4: V4 * (1 - 0.181 * SEXF), i.e. 18% lower in females (Huestis 2025 Table 2 'Sex effect on V4' = -0.181 and Results 'PopPK model analysis'). Not retained on any MDMA parameter.",
      source_name = "sex"
    ),
    STUDY_MPKF = list(
      description = "MPKF food-effect study indicator (1 = record from study MPKF, 0 = record from the NIDA study)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NIDA study)",
      notes = "Empirical study effect on MDMA CL/F (-36%) and V2/F (-26.3%), applied proportionally as (1 + theta * STUDY_MPKF), and the selector between the two study-specific additive residual errors of MDMA (Huestis 2025 Table 2 and Results 'PopPK model analysis'). The paper attributes the effect to an idiosyncrasy of the pooled dataset (washout, drug-use history or other between-study differences) and does not consider it meaningful for patients; set STUDY_MPKF = 0 to simulate the NIDA reference. No study effect was retained in the MDA model.",
      source_name = "study"
    ),
    FED = list(
      description = "Any meal before dosing (1 = light breakfast or high-fat meal, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "The paper's three-level prandial factor (Results 'PopPK model analysis') is decomposed as: NIDA light breakfast (the paper's reference level) FED = 1, FED_HIGHFAT = 0; MPKF high-fat meal FED = 1, FED_HIGHFAT = 1; MPKF fasted FED = 0, FED_HIGHFAT = 0. The fasted effect on ka is therefore applied as (1 + 1.93 * (1 - FED)), so the light-breakfast reference keeps the tabulated typical ka of 1.17 1/h. MPKF fasting: at least 10 h before dosing; NIDA: light meal before dosing (Supplementary Materials).",
      source_name = "fed status (fasted)"
    ),
    FED_HIGHFAT = list(
      description = "High-fat/high-calorie meal before dosing (1 = high-fat meal, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (light breakfast or fasted)",
      notes = "Proportional effect on MDMA ka: ka * (1 - 0.484 * FED_HIGHFAT), i.e. 48% slower absorption than the NIDA light-breakfast reference (Huestis 2025 Table 2 'Fed status effect on ka'). The MPKF meal was high fat (~50% of calories) and high calorie (~800-1000 kcal), eaten within 30 min before dosing (Supplementary Materials). Records with FED_HIGHFAT = 1 must also carry FED = 1.",
      source_name = "fed status (fed)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and central volume of MDMA and MDA in the stepwise covariate search but not retained (Huestis 2025 Table S1 and Results).",
      source_name = "age"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL and central volume but not retained (Huestis 2025 Table S1). Creatinine clearance (mL/min, not BSA-normalised) was likewise tested on clearance and not retained; the Figure 2 reference individual has 113 mL/min.",
      source_name = "albumin"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on clearance but not retained (Huestis 2025 Table S1).",
      source_name = "ALT"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White)",
      notes = "Race was tested on CL and central volume but not retained (Huestis 2025 Table S1).",
      source_name = "race"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "midomafetamine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midomafetamine", units = "mg", specimen = "plasma", verified = TRUE),
    central_tenamfetamine = list(analyte = "tenamfetamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_tenamfetamine = list(analyte = "tenamfetamine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 2L,
    age_range = "18.0-60.1 years",
    age_median = "25.0 years",
    weight_range = "43.1-113.1 kg",
    weight_median = "73.8 kg",
    sex_female_pct = 38.5,
    race_ethnicity = c(White = 40.0, Black = 58.5, Missing = 1.5),
    disease_state = "Healthy adult volunteers",
    dose_range = "Single oral MDMA HCl doses: 120 mg fasted and after a high-fat meal (MPKF crossover); 1.0 and 1.6 mg/kg after a light meal (NIDA crossover; 43-150 mg)",
    regions = "United States",
    notes = "Pooled from two phase I crossover studies: MPKF (N = 16; open-label randomised-sequence food-effect study, 120 mg MDMA HCl fasted vs high-fat/high-calorie meal, sampling to 72 h) and the NIDA study (N = 49; double-blind placebo-controlled within-subject study, 1.0 and 1.6 mg/kg MDMA HCl after a light meal, sampling to 167 h). 2012 MDMA and 1871 MDA plasma concentrations (Huestis 2025 Results; demographics Supplementary Table S4). NONMEM 7.4, FOCE-I."
  )

  ini({
    # MDMA (midomafetamine) structural parameters. Reference individual: 70 kg,
    # NIDA study, light breakfast before dosing. Huestis 2025 Table 2.
    lcl <- log(41.5); label("Apparent clearance of MDMA CL/F at 70 kg, NIDA study (L/h)") # Table 2 MDMA 'CL/F (L/h)' = 41.5 (RSE 4.70%)
    lvc <- log(462); label("Apparent central volume of MDMA V2/F at 70 kg, NIDA study (L)") # Table 2 MDMA 'V2/F (L)' = 462 (RSE 3.46%)
    lka <- log(1.17); label("First-order absorption rate constant of MDMA after a light breakfast (1/h)") # Table 2 MDMA 'ka (1/h)' = 1.17 (RSE 8.26%)
    ltlag <- log(0.322); label("Absorption lag time ALAG1 (h)") # Table 2 MDMA 'ALAG1 (h)' = 0.322 (RSE 5.15%)
    ld1 <- log(0.304); label("Duration of the zero-order input into the depot D1 (h)") # Table 2 MDMA 'D1 (h)' = 0.304 (RSE 18.2%)

    # MDMA covariate effects. Continuous effects are power functions of WT/70
    # (Methods Equation 3); categorical effects are proportional, (1 + theta * X)
    # (Methods Equation 4).
    e_wt_cl <- 1.03; label("Power exponent of body weight on MDMA CL/F (unitless)") # Table 2 'Body weight effect on CL/F' = 1.03 (RSE 17.2%)
    e_wt_vc <- 0.885; label("Power exponent of body weight on MDMA V2/F (unitless)") # Table 2 'Body weight effect on V2/F' = 0.885 (RSE 13.1%)
    e_study_mpkf_cl <- -0.36; label("Proportional change in MDMA CL/F for study MPKF vs NIDA (fraction)") # Table 2 'Study effect on CL/F' = -0.36 (RSE 16.40%)
    e_study_mpkf_vc <- -0.263; label("Proportional change in MDMA V2/F for study MPKF vs NIDA (fraction)") # Table 2 'Study effect on V2/F' = -0.263 (RSE 14.40%)
    e_fed_highfat_ka <- -0.484; label("Proportional change in MDMA ka after a high-fat meal vs light breakfast (fraction)") # Table 2 'Fed status effect on ka' = -0.484 (RSE 18.6%)
    e_fasted_ka <- 1.93; label("Proportional change in MDMA ka when fasted vs light breakfast (fraction)") # Table 2 'Fasted status effect on ka' = 1.93 (RSE 72.3%)

    # Fraction of MDMA clearance forming MDA, held fixed.
    fm <- fixed(0.1); label("Fraction of MDMA clearance converted to MDA (fraction)") # Table 2 'fmet' = 0.1 FIXED

    # MDA (tenamfetamine) structural parameters, apparent relative to fm. Reference
    # individual: 70 kg, male. Huestis 2025 Table 2.
    lcl_tenamfetamine <- log(34.9); label("Apparent clearance of MDA CLM at 70 kg (L/h)") # Table 2 MDA 'CLM (L/h)' = 34.9 (RSE 3.7%)
    lvc_tenamfetamine <- log(20.8); label("Apparent central volume of MDA V3 (L)") # Table 2 MDA 'V3 (L)' = 20.8 (RSE 22.6%)
    lq_tenamfetamine <- log(79.6); label("Apparent intercompartmental clearance of MDA QM at 70 kg (L/h)") # Table 2 MDA 'QM (L/h)' = 79.6 (RSE 10.3%)
    lvp_tenamfetamine <- log(213); label("Apparent peripheral volume of MDA V4 in males (L)") # Table 2 MDA 'V4 (L)' = 213 (RSE 7.66%)

    # MDA covariate effects.
    e_wt_cl_tenamfetamine <- 0.580; label("Power exponent of body weight on MDA CLM (unitless)") # Table 2 'Body weight effect on CLM' = 0.580 (RSE 26.2%)
    e_wt_q_tenamfetamine <- 1.00; label("Power exponent of body weight on MDA QM (unitless)") # Table 2 'Body weight effect on QM' = 1.00 (RSE 29.8%)
    e_sexf_vp_tenamfetamine <- -0.181; label("Proportional change in MDA V4 for females vs males (fraction)") # Table 2 'Sex effect on V4' = -0.181 (RSE 38%)

    # IIV. Table 2 prints %CV = sqrt(omega^2) * 100 (Methods, below Equation 1),
    # so omega^2 = (CV/100)^2; covariances = r * sd_i * sd_j with the r values
    # printed beneath Table 2.
    etalcl + etalvc ~ c(
      0.093025,
      0.04066565, 0.039601
    ) # Table 2: IIV CL/F 30.5%, IIV V2/F 19.9%; 'Correlation of CL/F, V2/F: r = 0.67' -> 0.67 * 0.305 * 0.199
    etalka ~ 0.362404 # Table 2: IIV on KA 60.2% CV -> 0.602^2
    etald1 ~ 2.3409 # Table 2: IIV on D1 153% CV -> 1.53^2
    etalcl_tenamfetamine + etalq_tenamfetamine + etalvp_tenamfetamine ~ c(
      0.1024,
      0.1228032, 0.219024,
      0.0968704, 0.1153152, 0.123904
    ) # Table 2: IIV CLM 32.0%, QM 46.8%, V4 35.2%; r(CLM,QM) = 0.82, r(CLM,V4) = 0.86, r(QM,V4) = 0.70

    # Residual error, Equation 2: Y = C * (1 + eps1) + eps2.
    propSd <- 0.207; label("Proportional residual error of MDMA (fraction)") # Table 2 MDMA 'Prop. Error' = 0.207 (RSE 9.24%)
    addSd <- 8.11; label("Additive residual error of MDMA, NIDA study (ng/mL)") # Table 2 MDMA 'Add. Error NIDA' = 8.11 (RSE 20.8%)
    addSd_mpkf <- fixed(2.5); label("Additive residual error of MDMA, MPKF study (ng/mL)") # Table 2 MDMA 'Add. Error MPKF' = 2.5 FIX
    propSd_tenamfetamine <- 0.122; label("Proportional residual error of MDA (fraction)") # Table 2 MDA 'Prop. Error' = 0.122 (RSE 32.8%)
    addSd_tenamfetamine <- 0.776; label("Additive residual error of MDA (ng/mL)") # Table 2 MDA 'Add. Error' = 0.776 (RSE 19.6%)
  })

  model({
    # MDMA individual parameters
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_study_mpkf_cl * STUDY_MPKF)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (1 + e_study_mpkf_vc * STUDY_MPKF)
    ka <- exp(lka + etalka) * (1 + e_fed_highfat_ka * FED_HIGHFAT) * (1 + e_fasted_ka * (1 - FED))
    tlag <- exp(ltlag)
    d1 <- exp(ld1 + etald1)

    # MDA individual parameters
    cl_tenamfetamine <- exp(lcl_tenamfetamine + etalcl_tenamfetamine) * (WT / 70)^e_wt_cl_tenamfetamine
    vc_tenamfetamine <- exp(lvc_tenamfetamine)
    q_tenamfetamine <- exp(lq_tenamfetamine + etalq_tenamfetamine) * (WT / 70)^e_wt_q_tenamfetamine
    vp_tenamfetamine <- exp(lvp_tenamfetamine + etalvp_tenamfetamine) * (1 + e_sexf_vp_tenamfetamine * SEXF)

    kel <- cl / vc
    kel_tenamfetamine <- cl_tenamfetamine / vc_tenamfetamine
    k12_tenamfetamine <- q_tenamfetamine / vc_tenamfetamine
    k21_tenamfetamine <- q_tenamfetamine / vp_tenamfetamine

    # The dose enters the depot after the lag as a zero-order input of duration
    # D1 (dose records carry rate = -2), then is absorbed first-order.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_tenamfetamine) <- fm * kel * central - kel_tenamfetamine * central_tenamfetamine -
      k12_tenamfetamine * central_tenamfetamine + k21_tenamfetamine * peripheral1_tenamfetamine
    d/dt(peripheral1_tenamfetamine) <- k12_tenamfetamine * central_tenamfetamine - k21_tenamfetamine * peripheral1_tenamfetamine

    alag(depot) <- tlag
    dur(depot) <- d1

    # Study-specific additive residual error for MDMA
    addSd_i <- addSd * (1 - STUDY_MPKF) + addSd_mpkf * STUDY_MPKF

    # Amounts in mg, volumes in L -> mg/L; x 1000 gives ng/mL
    Cc <- 1000 * central / vc
    Cc_tenamfetamine <- 1000 * central_tenamfetamine / vc_tenamfetamine
    Cc ~ add(addSd_i) + prop(propSd)
    Cc_tenamfetamine ~ add(addSd_tenamfetamine) + prop(propSd_tenamfetamine)
  })
}
