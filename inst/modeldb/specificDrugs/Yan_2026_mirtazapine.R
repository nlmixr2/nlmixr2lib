Yan_2026_mirtazapine <- function() {
  description <- "One-compartment population PK model for oral mirtazapine with first-order absorption in Chinese inpatients with depression, built from routine therapeutic-drug-monitoring trough samples (Yan 2026). Apparent oral clearance is reduced by 29.4% in obese patients (BMI >= 28 kg/m^2), by 26.9% with concomitant paroxetine and by 51.1% with concomitant fluvoxamine (multiplicative); the absorption rate constant is fixed to a literature value. Between-subject variability on CL/F and V/F, proportional residual error."
  reference <- paste(
    "Yan H, Huang W, Xia H, Luo Y, Li Y, Shang D.",
    "Optimizing Mirtazapine Initial Dosing: A Population Analysis of the Effects of",
    "BMI, Paroxetine and Fluvoxamine.",
    "Drug Design, Development and Therapy. 2026;20:601238.",
    "doi:10.2147/DDDT.S601238.",
    "Final model equations in Results ('Development of Population Pharmacokinetic Model');",
    "parameter estimates Table 2.",
    "The fixed absorption rate constant is attributed to the paper's reference 13:",
    "Grasmader K, Verwohlt PL, Kuhn KU, Dragicevic A, et al. Population pharmacokinetic",
    "analysis of mirtazapine. Eur J Clin Pharmacol. 2004;60(7):473-480.",
    "doi:10.1007/s00228-004-0737-0.",
    sep = " "
  )
  vignette <- "Yan_2026_mirtazapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "mirtazapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mirtazapine", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    DIS_OBESE = list(
      description = "Obesity indicator (1 = BMI >= 28 kg/m^2, 0 = BMI < 28 kg/m^2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BMI < 28 kg/m^2)",
      notes = paste(
        "The paper's BMI covariate is a dichotomised body-mass index, not the continuous",
        "value (Methods, Covariate Model: 'For BMI covariates, a COV of 1 represents",
        "BMI >= 28 kg/m2 and 0 represents BMI < 28 kg/m2'). The 28 kg/m^2 cut-off is the",
        "Chinese adult obesity threshold of the Working Group on Obesity in China",
        "(Zhou 2002, the paper's reference 25). 20 of 105 patients (19.05%) were",
        "classed as obese (Table 1). Time-fixed per patient in this source."
      ),
      source_name = "BMI"
    ),
    CONMED_PAROXETINE = list(
      description = "Concomitant paroxetine indicator (1 = co-administered, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant paroxetine)",
      notes = paste(
        "The paper's `PRX` term (Results, CL/F equation; Methods: 'For concomitant",
        "medication, a COV of 1 represents co-administration and 0 represents no",
        "co-administration'). 13 of 105 patients (12.38%; Table 1). The paper attributes",
        "the effect to paroxetine's CYP2D6 inhibition (Discussion). Time-fixed per patient",
        "in this source."
      ),
      source_name = "PRX"
    ),
    CONMED_FLUVOXAMINE = list(
      description = "Concomitant fluvoxamine indicator (1 = co-administered, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant fluvoxamine)",
      notes = paste(
        "The paper's `FLV` term (Results, CL/F equation). Only 6 of 105 patients",
        "(5.71%; Table 1), which the paper flags as a limitation. The paper attributes the",
        "effect to fluvoxamine's strong CYP1A2 inhibition plus its CYP3A4 and CYP2C19",
        "inhibition, CYP1A2 and CYP3A4 being main mirtazapine pathways (Discussion).",
        "Time-fixed per patient in this source."
      ),
      source_name = "FLV"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F but not retained (Supplementary Table 1, first forward round, dOFV -3.600, P > 0.01). 62 female / 43 male (Table 1); the paper codes 1 = female, the same orientation as SEXF."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F (linear, median-centred) but not retained (Supplementary Table 1, dOFV -0.896). Cohort mean 48.52 years (SD 20.80), range 13-79 (Table 1)."
    ),
    CONMED_VENLAFAXINE = list(
      description = "Concomitant venlafaxine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Supplementary Table 1, dOFV -3.592). 26 of 105 patients (Table 1)."
    ),
    CONMED_FLUOXETINE = list(
      description = "Concomitant fluoxetine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Supplementary Table 1, dOFV -0.488). 5 of 105 patients (Table 1)."
    ),
    CONMED_DULOXETINE = list(
      description = "Concomitant duloxetine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Supplementary Table 1, dOFV -0.723). 12 of 105 patients (Table 1)."
    ),
    CONMED_RISPERIDONE = list(
      description = "Concomitant risperidone indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Supplementary Table 1, dOFV -4.033). 10 of 105 patients (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (Supplementary Table 1, dOFV -0.109). Cohort mean 26.84 U/L (SD 21.56; Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (Supplementary Table 1, dOFV +0.28). Cohort mean 23.18 U/L (SD 12.13; Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened but not retained (Supplementary Table 1, dOFV +0.131). Cohort mean 67.11 umol/L (SD 13.63; Table 1)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalised)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened but not retained (Supplementary Table 1, 'GFR', dOFV -0.054). Cohort mean 99.88 (SD 22.60; Table 1); the estimating equation is not stated."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 105,
    n_studies = 1,
    n_observations = 210,
    age_range = "13-79 years",
    age_mean = "48.52 years (SD 20.80)",
    weight_range = "25.4-104.0 kg",
    weight_mean = "62.19 kg (SD 13.27)",
    height_mean = "1.63 m (SD 0.09)",
    sex_female_pct = 59.05,
    race_ethnicity = "Chinese (inclusion criterion)",
    disease_state = "Inpatients with depressive disorders (major depressive disorder, recurrent depressive disorder, post-schizophrenic depression, mixed anxiety and depressive disorder, depressive disorder with anxious distress or severe depressive episode), ICD-11 / DSM-5.",
    dose_range = "7.5-45 mg mirtazapine tablets once nightly (mean daily dose 31.92 mg, SD 11.71)",
    regions = "China (single centre: The Affiliated Brain Hospital of Guangzhou Medical University)",
    notes = paste(
      "Retrospective therapeutic-drug-monitoring data collected 1 May 2024 - 15 September",
      "2025; steady-state serum samples drawn 08:00-09:00 the morning after the nightly",
      "dose, i.e. roughly 12 h post-dose (Methods). Patients with severe hepatic disease or",
      "eGFR < 30 mL/min/1.73 m^2 were excluded; no smokers were included (Limitations).",
      "20 of 105 patients had BMI >= 28 kg/m^2; 13 took paroxetine and 6 fluvoxamine",
      "(Table 1). Estimation: NONMEM 7.3, ADVAN2, FOCE-I. Evaluated by GOF, a 1000-run",
      "bootstrap (92.4% convergence), NPDE, an independent 49-patient TDM cohort and",
      "literature mean concentrations (Figure 4)."
    )
  )

  ini({
    # Structural parameters -- Results, final-model equations and Table 2.
    # Reference subject: BMI < 28 kg/m^2, no paroxetine, no fluvoxamine.
    lcl <- log(28.9); label("Apparent oral clearance CL/F, non-obese without paroxetine or fluvoxamine (L/h)") # Table 2: CL/F = 28.9 L/h (RSE 4%; bootstrap median 28.515 [25.423, 31.343])
    lvc <- log(310); label("Apparent volume of distribution V/F (L)") # Table 2: V/F = 310 L (RSE 22%; bootstrap median 302.315 [191.154, 462.256])

    # Methods, Base Model: Ka 'was fixed at a literature-based value of 1.2 h-1'
    # (Discussion cites reference 13, Grasmader 2004). Table 2 prints the unit
    # as 'L/h', a typo for 1/h.
    lka <- fixed(log(1.2)); label("First-order absorption rate constant ka (1/h)") # Table 2: Ka = 1.2 FIX

    # Covariate effects on CL/F, proportional (linear-shift) form for
    # categorical covariates (Methods, Covariate Model), combined
    # multiplicatively in the Results equation:
    #   CL/F = 28.9 * (1 - 0.294*BMI) * (1 - 0.269*PRX) * (1 - 0.511*FLV)
    e_obese_cl <- -0.294; label("Proportional change in CL/F when BMI >= 28 kg/m^2 (fraction)") # Table 2: Theta CL-BMI = -0.294 (RSE 17%)
    e_par_cl <- -0.269; label("Proportional change in CL/F with concomitant paroxetine (fraction)") # Table 2: Theta CL-Paroxetine = -0.269 (RSE 24%)
    e_fluvox_cl <- -0.511; label("Proportional change in CL/F with concomitant fluvoxamine (fraction)") # Table 2: Theta CL-Fluvoxamine = -0.511 (RSE 13%)

    # IIV, exponential model (Methods Eq. P_i = P_pop * exp(eta_i), eta ~ N(0,
    # omega^2)). The Table 2 'IIV' column is read as omega^2 (variance); see
    # the vignette for the check against the paper's Figure 5 PTA values.
    etalcl ~ 0.0613 # Table 2: IIV CL/F = 0.0613 (omega^2; 25.1% CV)
    etalvc ~ 0.131 # Table 2: IIV V/F = 0.131 (omega^2; 37.4% CV)

    # Proportional residual error Y = F * (1 + eps), eps ~ N(0, sigma^2)
    # (Methods). Table 2 'PRO' = 0.0747 is read as sigma^2, so the SD is
    # sqrt(0.0747) = 0.2733.
    propSd <- sqrt(0.0747); label("Proportional residual error (fraction)") # Table 2: PRO = 0.0747 (RSE 16%), read as a variance
  })

  model({
    # 1. Individual parameters (Results, final-model equations)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (1 + e_obese_cl * DIS_OBESE) *
      (1 + e_par_cl * CONMED_PAROXETINE) *
      (1 + e_fluvox_cl * CONMED_FLUVOXAMINE)
    vc <- exp(lvc + etalvc)

    # 2. Micro-constants
    kel <- cl / vc

    # 3. One-compartment model with first-order oral absorption (NONMEM ADVAN2).
    #    Oral-only data: cl and vc are the apparent CL/F and V/F, so no
    #    bioavailability term is applied.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Observation. Dose in mg and vc in L give mg/L; x1000 converts to the
    #    ng/mL in which the paper reports serum mirtazapine.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
