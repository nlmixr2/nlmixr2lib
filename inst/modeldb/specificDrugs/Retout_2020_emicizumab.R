Retout_2020_emicizumab <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and first-order elimination (no lag time) for emicizumab, a bispecific anti-FIXa/FX humanized monoclonal antibody, in adult, adolescent and pediatric (1 year and older) persons with hemophilia A with or without factor VIII inhibitors (Retout 2020; phase I/II + HAVEN 1-4). Body weight (power, 70 kg) and albumin (linear, 45 g/L) on CL/F, body weight (power) and Black race (fractional) on V/F, and a linear decline of the apparent bioavailability F above 30 years of age."
  reference <- "Retout S, Schmitt C, Petry C, Mercier F, Frey N. Population Pharmacokinetic Analysis and Exploratory Exposure-Bleeding Rate Relationship of Emicizumab in Adult and Pediatric Persons with Hemophilia A. Clin Pharmacokinet. 2020;59(12):1611-1625. doi:10.1007/s40262-020-00904-z"
  vignette <- "Retout_2020_emicizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F (exponent 0.911) and V/F (exponent 1.00) normalized to 70 kg (Retout 2020 Table 4 header 'Fixed effects (BW 70 kg; ...)' and the covariate equations in Section 3.1). Analysis-population range 9.50-156 kg, median 69.1 kg (Table 2).",
      source_name = "BW"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL/F: CL/F multiplied by (1 - 0.0157 * (ALB - 45)), i.e. LOWER clearance at HIGHER albumin (Retout 2020 Section 3.1 equation). Median 45.0 g/L, range 16.8-56.6 g/L (Table 2); two PwHA with abnormally low albumin (16.8 and 27.0 g/L) were excluded from covariate model development.",
      source_name = "ALB"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Hockey-stick effect on the apparent bioavailability F: F = 1 for AGE <= 30 years, else F = 1 - 0.00651 * (AGE - 30) (Retout 2020 Section 3.1 equation). Range 1.22-77.0 years, median 30.0 years (Table 2).",
      source_name = "AGE"
    ),
    RACE_BLACK = list(
      description = "Black race indicator: 1 = Black, 0 = any other race (White, Asian including Japanese, other or unknown).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Fractional effect on V/F: V/F multiplied by (1 - 0.215 * BLK) (Retout 2020 Section 3.1: 'where BLK is 1 if persons are black, and 0 otherwise'). 31/389 (8.0%) of the analysis population were Black (Table 3).",
      source_name = "BLK"
    )
  )

  covariatesDataExcluded <- list(
    FVIII_INH = list(
      description = "Factor VIII inhibitor status (1 = with FVIII inhibitors)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (without FVIII inhibitors)",
      notes = "Tested in the covariate search and not retained (Retout 2020 Section 4: 'None of the other covariates tested either statistically (i.e., FVIII inhibitor or non-inhibitor status, body mass index, body surface area) ... was found to further explain the PK variability'). 194/389 (49.9%) had FVIII inhibitors (Table 3).",
      source_name = "status"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested and not retained (Retout 2020 Section 4).",
      source_name = "BMI"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested and not retained (Retout 2020 Section 4).",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "emicizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "emicizumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 389L,
    n_studies = 5L,
    age_range = "1.22-77.0 years",
    age_median = "30.0 years",
    weight_range = "9.50-156 kg",
    weight_median = "69.1 kg",
    race_ethnicity = c(White = 62.7, Black = 8.0, Asian = 22.9, Other_or_unknown = 6.4),
    disease_state = "Hemophilia A with (49.9%) or without (50.1%) factor VIII inhibitors; pediatric (1 to < 12 years), adolescent and adult persons with hemophilia A (PwHA), all emicizumab-naive at study entry.",
    dose_range = "Subcutaneous emicizumab: 3 mg/kg QW for 4 weeks followed by 1.5 mg/kg QW, 3 mg/kg Q2W or 6 mg/kg Q4W (HAVEN 1-4); 6 mg/kg Q4W without loading (HAVEN 4 run-in); 1 mg/kg f/b 0.3 mg/kg QW, 3 mg/kg f/b 1 mg/kg QW, or 3 mg/kg QW (Japanese phase I/II). Eleven HAVEN patients were up-titrated to 3 mg/kg QW.",
    regions = "Multinational (HAVEN 1-4) and Japan (phase I/II).",
    studies = "Japanese phase I/II (JapicCTI-121934 / JapicCTI-132195, n = 18), HAVEN 1 (NCT02622321, n = 112), HAVEN 2 (NCT02795767, pediatric, n = 63), HAVEN 3 (NCT02847637, n = 148), HAVEN 4 (NCT03020160, n = 48).",
    regimen_n = "QW 292 (75.1%), Q2W 49 (12.6%), Q4W run-in 7 (1.8%), Q4W expansion 41 (10.5%) (Table 3).",
    albumin = "45.0 g/L median (16.8-56.6)",
    notes = "Demographics from Retout 2020 Tables 1-3. 4966 plasma concentrations from 383 evaluable PK profiles were used for base-model development; the covariate model used 381 PwHA after excluding two with abnormally low albumin. One PwHA with an anti-drug-antibody-associated decline in exposure was excluded from model development. Sex is not tabulated in the paper (hemophilia A is X-linked; the cohort is expected to be male)."
  )

  ini({
    # Structural parameters - Retout 2020 Table 4 final-model fixed effects for
    # the reference PwHA (BW 70 kg; ALB 45 g/L; age <= 30 years; non-Black).
    lka <- log(0.536); label("First-order SC absorption rate constant ka (1/day)") # Retout 2020 Table 4: KA = 0.536 1/day (RSE 7.1%)
    lcl <- log(0.272); label("Apparent clearance CL/F (L/day)") # Retout 2020 Table 4: CL/F = 0.272 L/day (RSE 1.9%)
    lvc <- log(10.4); label("Apparent volume of distribution V/F (L)") # Retout 2020 Table 4: V/F = 10.4 L (RSE 1.9%)

    # Covariate effects - Retout 2020 Table 4 and the Section 3.1 equations
    #   CL/F = 0.272 * (BW/70)^0.911 * (1 - 0.0157 * (ALB - 45))
    #   V/F  = 10.4  * (BW/70)^1.00  * (1 - 0.215 * BLK)
    #   F    = 1 if AGE <= 30 else 1 - 0.00651 * (AGE - 30)
    # Table 4 prints the ALB and age coefficients as positive magnitudes and
    # the Black-race coefficient as -0.215; the printed equations carry the
    # minus sign for all three. The model() block follows the equations
    # (confirmed by the paper's Cav,SS predictions: 46.4 ug/mL at ALB 33 g/L,
    # 68.0 ug/mL at ALB 57 g/L, -31% exposure at age 77 years).
    e_wt_cl <- 0.911; label("Power exponent of (WT/70) on CL/F (unitless)") # Retout 2020 Table 4: Effect of BW on CL/F = 0.911 (RSE 3.2%)
    e_alb_cl <- 0.0157; label("Fractional decrease in CL/F per g/L of albumin above 45 g/L (1/(g/L))") # Retout 2020 Table 4: Effect of ALB on CL/F = 1.57e-2 (RSE 28.4%); sign from the Section 3.1 equation '(1 - 0.0157 x (ALB - 45))'
    e_wt_vc <- 1.00; label("Power exponent of (WT/70) on V/F (unitless)") # Retout 2020 Table 4: Effect of BW on V/F = 1.00 (RSE 3.0%)
    e_race_black_vc <- -0.215; label("Fractional change in V/F for Black race (unitless)") # Retout 2020 Table 4: Effect of Black on V/F = -0.215 (RSE 19.7%)
    e_age_fdepot <- 0.00651; label("Fractional decrease in bioavailability F per year of age above 30 years (1/year)") # Retout 2020 Table 4: Effect of AGE>30 years on F = 6.51e-3 (RSE 16.3%); sign from the Section 3.1 equation 'F = 1 - 0.00651 x (AGE - 30)'

    # Between-person variability - Retout 2020 Table 4 (exponential model,
    # reported as CV%); omega^2 = log(1 + CV^2). Covariances from the reported
    # correlations: cov = r * sqrt(omega_i^2 * omega_j^2). The V/F-KA
    # correlation is not reported in Table 4 and is set to 0.
    #   CL/F: CV 28.7%  -> omega^2 = log(1 + 0.287^2) = 0.079152
    #   V/F:  CV 25.9%  -> omega^2 = log(1 + 0.259^2) = 0.064927
    #   KA:   CV 72.5%  -> omega^2 = log(1 + 0.725^2) = 0.422400
    #   r(CL/F, V/F) =  0.217  -> cov =  0.015556
    #   r(CL/F, KA)  = -0.341  -> cov = -0.062352
    etalcl + etalvc + etalka ~ c(
      0.079152,
      0.015556, 0.064927,
      -0.062352, 0, 0.422400
    ) # Retout 2020 Table 4: BPV CL/F 28.7%, V/F 25.9%, KA 72.5% CV; correlations CL/F-V/F 0.217, CL/F-KA -0.341

    # Residual error - Retout 2020 Table 4 (combined additive + proportional).
    addSd <- fixed(0.025); label("Additive residual error (ug/mL)") # Retout 2020 Table 4: sigma1 (additive) = 0.025 ug/mL Fix (half the phase I/II LLOQ of 50 ng/mL, Section 3.1)
    propSd <- 0.146; label("Proportional residual error (fraction)") # Retout 2020 Table 4: sigma2 (proportional) = 14.6% (RSE 2.0%)
  })

  model({
    # Individual PK parameters (Retout 2020 Section 3.1 equations).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 - e_alb_cl * (ALB - 45))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (1 + e_race_black_vc * RACE_BLACK)

    # Apparent bioavailability: 1 up to 30 years, then a linear decline.
    fdepot <- 1 - e_age_fdepot * (AGE > 30) * (AGE - 30)

    kel <- cl / vc

    # One-compartment model with first-order SC absorption.
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central
    f(depot) <- fdepot

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
