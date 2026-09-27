Chen_2021_luspatercept <- function() {
  description <- "One-compartment population PK model for subcutaneous luspatercept (activin receptor type IIB / IgG1 Fc-fusion protein) in adults with transfusion-dependent beta-thalassemia (Chen 2021), with first-order absorption and first-order linear elimination parameterised in CL/F and V1/F; body weight and baseline albumin power covariates and an exponential baseline red-blood-cell transfusion burden covariate on CL/F, and body weight power and exponential baseline transfusion burden covariates on V1/F."
  reference <- "Chen N, Kassir N, Laadem A, Giuseppi AC, Shetty J, Maxwell SE, Sriraman P, Ritland S, Linde PG, Budda B, Reynolds JG, Zhou S, Palmisano M. Population Pharmacokinetics and Exposure-Response Relationship of Luspatercept, an Erythroid Maturation Agent, in Anemic Patients With beta-Thalassemia. J Clin Pharmacol. 2021;61(1):52-63. doi:10.1002/jcph.1696. PMCID: PMC7754485."
  vignette <- "Chen_2021_luspatercept"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "luspatercept", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "luspatercept", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 70 kg per the Chen 2021 final-model covariate equations for CL/F and V1/F (Results). Power exponents 0.806 (CL/F) and 0.705 (V1/F) per Table 2. Observed median 57.1 kg (range 34.1-97.0) per Table 1.",
      source_name = "Weight"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 46 g/L (the dataset median, Table 1) per the Chen 2021 final-model CL/F equation. Power exponent -0.881 per Table 2. Observed range 30.0-56.0 g/L.",
      source_name = "Albumin"
    ),
    RBCT_BL = list(
      description = "Baseline red-blood-cell transfusion burden",
      units = "RBC units/24 weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred at 14 RBC units/24 weeks in the Chen 2021 final-model equations (Results: 'e^(-0.0118 x [RBCT - 14])' on CL/F and 'e^(-0.0141 x [RBCT - 14])' on V1/F); 14 is the Table 1 median of 14.1 rounded. Observed range 0-34.0 RBC units/24 weeks.",
      source_name = "RBCT (RBCT burden, units/24 weeks)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the full-model covariate analysis and dropped as insignificant or of low clinical relevance (Results; Discussion: 'In this young patient population, age was not a significant covariate of luspatercept PK'). Median 32 years (18-66)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Tested and dropped as insignificant or of low clinical relevance (Results). 56.8% female."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "White",
      notes = "Asian versus White tested and dropped as insignificant or of low clinical relevance (Results). 28.8% Asian."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Mild to moderate renal impairment tested and dropped as insignificant or of low clinical relevance (Results). Median 120.0 (53.7-314.0)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 285L,
    n_studies = 2L,
    n_observations = 3680L,
    age_range = "18-66 years",
    age_median = "32 years",
    weight_range = "34.1-97.0 kg",
    weight_median = "57.1 kg",
    sex_female_pct = 56.8,
    race_ethnicity = c(White = 63.5, Asian = 28.8, Other = 7.7),
    disease_state = "Adults with beta-thalassemia requiring regular red-blood-cell transfusions (baseline transfusion burden median 14.1, range 0-34.0 RBC units/24 weeks). 59.6% splenectomised; genotype beta0/beta0 23.5%, non-beta0/beta0 53.7%, missing 22.8%. 84.9% on concurrent iron chelation therapy.",
    dose_range = "Subcutaneous luspatercept 0.2-1.25 mg/kg once every 3 weeks (q3w). Dose-escalation cohorts received a single dose level (0.2-1.25 mg/kg); expansion and phase 3 patients started at 0.8 or 1 mg/kg with stepwise escalation to 1 and 1.25 mg/kg. 79.8% started at 1 mg/kg; 34% escalated to 1.25 mg/kg in the first year.",
    regions = "Multinational: phase 2 dose-finding/expansion study A536-04 (NCT01749540, n = 64) with its extension A536-06 (NCT02268409), and pivotal phase 3 study ACE-536-B-THAL-001 BELIEVE (NCT02604433, n = 221).",
    notes = "Baseline demographics from Chen 2021 Table 1. Renal function: 86.0% normal (eGFR >= 90), 13.0% mild, 1.1% moderate impairment. Albumin median 46.0 g/L (30.0-56.0). 3680 quantifiable serum luspatercept concentrations collected on days 5-610 after the first dose; ELISA range 50-600 ng/mL; 0.6% of postdose samples below the limit of quantitation were excluded."
  )

  ini({
    # Structural PK parameters - Chen 2021 Table 2 final-model NONMEM
    # estimates. Reference subject: 70 kg, 46 g/L baseline albumin,
    # baseline transfusion burden 14 RBC units/24 weeks.
    lcl <- log(0.532); label("Apparent clearance CL/F (L/day) at reference covariates")        # Chen 2021 Table 2: CL/F = 0.532 L/day
    lvc <- log(8.39);  label("Apparent central volume V1/F (L) at reference covariates")       # Chen 2021 Table 2: V1/F = 8.39 L
    lka <- log(0.409); label("Absorption rate Ka (1/day)")                                     # Chen 2021 Table 2: Ka = 0.409 1/day

    # Covariate effects - Chen 2021 Table 2 and the final-model covariate
    # equations in Results: CL/F = 0.532 * (Weight/70)^0.806 *
    # (Albumin/46)^-0.881 * e^(-0.0118 * [RBCT - 14]) and
    # V1/F = 8.39 * (Weight/70)^0.705 * e^(-0.0141 * [RBCT - 14]).
    e_wt_cl     <-  0.806;  label("Power exponent of (WT/70 kg) on CL/F (unitless)")                          # Chen 2021 Table 2: Weight on CL/F = 0.806
    e_alb_cl    <- -0.881;  label("Power exponent of (ALB/46 g/L) on CL/F (unitless)")                        # Chen 2021 Table 2: Albumin on CL/F = -0.881
    e_rbct_cl   <- -0.0118; label("Exponential coefficient of (RBCT_BL - 14) on CL/F (per RBC unit/24 weeks)") # Chen 2021 Table 2: RBCT burden on CL/F = -0.0118
    e_wt_vc     <-  0.705;  label("Power exponent of (WT/70 kg) on V1/F (unitless)")                          # Chen 2021 Table 2: Weight on V1/F = 0.705
    e_rbct_vc   <- -0.0141; label("Exponential coefficient of (RBCT_BL - 14) on V1/F (per RBC unit/24 weeks)") # Chen 2021 Table 2: RBCT burden on V1/F = -0.0141

    # Inter-individual variability - Chen 2021 Table 2 reports IIV of CL/F
    # 34.7% and of V1/F 27.6% for an exponential IIV model (Methods). These
    # are read as sqrt(omega^2) x 100, the convention of the sibling Chen
    # 2020 luspatercept (MDS) analysis by the same group. Cross-check: with
    # sqrt(omega^2_CL) = 0.347, 100*sqrt(exp(0.347^2)-1) = 35.8%, matching
    # the 36% descriptive CV of AUCss reported in Results. No IIV on Ka
    # ('Inclusion of IIV for Ka led to large shrinkage'); no correlation
    # reported.
    etalcl ~ 0.1204  # Chen 2021 Table 2: IIV CL/F 34.7% read as sqrt(omega^2); 0.347^2
    etalvc ~ 0.0762  # Chen 2021 Table 2: IIV V1/F 27.6% read as sqrt(omega^2); 0.276^2

    # Residual error - Chen 2021 Table 2 'Residual variability' 20.8%.
    # Methods: concentrations were natural-log-transformed and residual
    # variability modelled with an additive error model, i.e. a
    # proportional error in linear space.
    propSd <- 0.208; label("Proportional residual error (fraction)")                          # Chen 2021 Table 2: residual variability = 20.8% (log-additive)
  })

  model({
    # Individual PK parameters with the Chen 2021 final-model covariate
    # equations (Results). Reference subject: 70 kg, 46 g/L baseline
    # albumin, baseline transfusion burden 14 RBC units/24 weeks.
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      (ALB / 46)^e_alb_cl *
      exp(e_rbct_cl * (RBCT_BL - 14))
    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc *
      exp(e_rbct_vc * (RBCT_BL - 14))
    ka <- exp(lka)

    # One-compartment SC model with first-order absorption and elimination
    # (Results: 'A 1-compartment model with first-order absorption and
    # elimination best described the concentration-time profiles of
    # luspatercept after subcutaneous injection'). Bioavailability is
    # absorbed into the apparent CL/F and V1/F.
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, volume in L: concentration in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
