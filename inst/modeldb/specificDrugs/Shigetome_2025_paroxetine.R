Shigetome_2025_paroxetine <- function() {
  description <- "One-compartment population PK model for paroxetine in 179 Japanese adults with major depressive disorder treated at Hirosaki University Hospital or Dokkyo Medical University, fitted to 329 steady-state plasma concentrations (therapeutic drug monitoring). Neither an absorption process nor saturable metabolism improved the fit, so oral doses enter the central compartment directly and CL/F and Vd/F are apparent oral parameters. Apparent clearance carries power covariate effects of age (exponent -0.536, higher age = lower CL) and total body weight (exponent 0.563), both normalized to the study means (41.5 years, 58.6 kg); CYP2D6 genotype was tested in three classifications and had no significant effect. Vd/F was fixed at 1010 L (mean body weight 58.6 kg times a literature Vd/F of 17.2 L/kg) because sparse trough sampling could not support its estimation. Exponential inter-individual variability on CL and a proportional residual error."

  reference <- paste(
    "Shigetome K, Egashira T, Tomita T, Higa N, Iwashita K, Morita K, Nishimura M,",
    "Kaneko T, Maeda H, Yamada KD, Kajiwara-Morita A, Oniki K, Yasui-Furukori N,",
    "Saruwatari J.",
    "Effect of Cumulative Exposure on the Efficacy of Paroxetine:",
    "A Population Pharmacokinetic-Pharmacodynamic and Machine Learning Analyses.",
    "CPT Pharmacometrics Syst Pharmacol. 2025 Jun;14(6):1119-1127.",
    "doi:10.1002/psp4.70032.",
    sep = " "
  )
  vignette <- "Shigetome_2025_paroxetine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Age at the time of paroxetine treatment",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters CL/F as the power term (AGE / 41.5)^e_age_cl, normalized to the mean age of the 179-patient PK population (Shigetome 2025 Table 1). The exponent is negative, so older patients have lower apparent clearance. Age reduced the objective function by 17.838 in the univariate screen and was retained through backward elimination (Table S1).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters CL/F as the power term (WT / 58.6)^e_wt_cl, normalized to the mean body weight of the 179-patient PK population (Shigetome 2025 Table 1). Sex was also a significant univariate covariate on CL but was dropped because it correlates with body weight (mean weight greater in males); the final model describes CL by age and body weight only. Vd/F was fixed and carries no covariates.",
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    central = list(analyte = "paroxetine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 179L,
    n_studies      = 1L,
    n_observations = 329L,
    age_range      = "16-80 years",
    age_mean       = "41.5 years (SD 14.2)",
    weight_range   = "40-85 kg",
    weight_mean    = "58.6 kg (SD 10.5)",
    sex_female_pct = 59.2,
    race_ethnicity = c(Asian = 100),
    disease_state  = "major depressive disorder (DSM-IV), antidepressant-naive at paroxetine initiation, normal renal and hepatic function, no interacting co-medication",
    dose_range     = "10-40 mg/day orally (mean 22.3 mg/day, SD 11.1); initial dose 10-20 mg/day then weekly increases of 10 mg/day to the maximum tolerated dose",
    regions        = "Japan (Hirosaki University Hospital and Dokkyo Medical University School of Medicine)",
    genotype       = "CYP2D6 intermediate metabolizers 121 (67.6%), extensive metabolizers 58 (32.4%), poor metabolizers 0 (0%); *1/*2 functional, *10/*41 reduced function, *5 non-functional",
    notes          = "Retrospective therapeutic-drug-monitoring cohort; 1.8 samples per patient on average (range 1-5), all drawn at steady state after at least 9 days on a stable dose. Observed plasma paroxetine concentrations 51.9 ng/mL (SD 49.0, range 2.6-286.5). Baseline demographics are Shigetome 2025 Table 1. A separate 50-patient subset with serial MADRS assessments supported the PK-PD model in Shigetome_2025_paroxetine_madrs."
  )

  ini({
    # ================================================================
    # Final popPK model (Shigetome 2025 Results section 3.2 printed
    # equation, parameter estimates in Table S2 of the Supporting
    # Information; NONMEM control stream in Text S2):
    #
    #   CL (L/h) = 23.6 * (AGE / 41.5)^-0.536 * (BW / 58.6)^0.563
    #   Vd/F (L) = 1010   (fixed)
    #
    # NONMEM ADVAN1 TRANS2 (one compartment, no absorption process),
    # S1 = V / 1000 so that concentrations are ng/mL for doses in mg.
    # ================================================================

    lcl <- log(23.6)
    label("Apparent oral clearance CL/F at the reference age (41.5 y) and weight (58.6 kg) (L/h)")  # Shigetome 2025 Table S2, CL = 23.6 (RSE 4.75%, bootstrap 95% CI 21.4-25.9)

    lvc <- fixed(log(1010))
    label("Apparent central volume of distribution Vd/F (L); mean body weight 58.6 kg times the literature Vd/F 17.2 L/kg, not estimated here")  # Shigetome 2025 Methods 2.3 and Table S2, Vd/F = 1010 L, NONMEM $THETA "1010 FIX"

    e_age_cl <- -0.536
    label("Power exponent on (AGE / 41.5) for CL/F (unitless)")  # Shigetome 2025 Results 3.2 equation and Table S2 (RSE 25.4%, bootstrap 95% CI -0.790 to -0.270)

    e_wt_cl <- 0.563
    label("Power exponent on (WT / 58.6) for CL/F (unitless)")  # Shigetome 2025 Results 3.2 equation and Table S2 (RSE 48.9%, bootstrap 95% CI -0.01 to 1.19)

    # ================================================================
    # Inter-individual variability: exponential on CL (Results 3.2;
    # NONMEM $PK "CL = ... * EXP(ETA(1))").
    #
    # Table S2 reports 0.307 in a column headed "%CV". The value is the
    # raw NONMEM OMEGA (a variance), not a CV: the same tables report
    # the popPK/PD additive residual as 190 with the endpoint measured
    # in percent, which can only be a variance (SD 13.8 percentage
    # points), and only the variance reading reproduces the observed
    # spread of concentrations. Var(ln C) = Var(ln dose) + omega^2 +
    # covariate variance + sigma^2 = 0.222 + 0.307 + 0.042 + 0.155 =
    # 0.725, i.e. CV(C) = 103%, against the observed 49.0/51.9 = 94%
    # (Table 1). Reading 0.307 and 0.155 as SDs instead gives CV(C) =
    # 68%, which the data exclude.
    # ================================================================
    etalcl ~ 0.307  # Shigetome 2025 Table S2, inter-individual variability on CL (RSE 15.7%, shrinkage 7.83%, bootstrap 95% CI 0.197-0.392)

    # ================================================================
    # Residual error: proportional (NONMEM $ERROR "W = IPRED;
    # Y = IPRED + W * EPS(1)"). SIGMA = 0.155 is a variance, so the
    # proportional SD is sqrt(0.155) = 0.394.
    # ================================================================
    propSd <- 0.394
    label("Proportional residual error SD (fraction)")  # Shigetome 2025 Table S2, proportional error 0.155 (variance; RSE 10.8%, shrinkage 19.3%, bootstrap 95% CI 0.123-0.192); sqrt(0.155) = 0.3937
  })

  model({
    # Apparent oral clearance with power effects of age and body weight,
    # normalized to the PK-population means (Shigetome 2025 Table 1).
    cl <- exp(lcl + etalcl) * (AGE / 41.5)^e_age_cl * (WT / 58.6)^e_wt_cl

    # Apparent central volume, fixed at 1010 L for every subject.
    vc <- exp(lvc)

    kel <- cl / vc

    # One compartment, no absorption process: the oral dose is placed
    # directly into central, so CL and Vd are apparent (CL/F, Vd/F).
    d/dt(central) <- -kel * central

    # central is in mg and vc in L, so central/vc is mg/L; the factor
    # 1000 converts to ng/mL and reproduces the NONMEM scaling
    # S1 = V / 1000.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
