Li_2025_levamlodipine <- function() {
  description <- "One-compartment population PK model with first-order absorption (ka fixed) and a direct-effect inhibitory Imax PD model with baseline for systolic blood pressure after oral levamlodipine besylate once daily in 43 Chinese adults with hypertension (9 smokers, 34 non-smokers) sampled at steady-state trough (Li 2025). Apparent clearance scales as a power of serum uric acid (exponent -0.69, reference 311.0 umol/L) and of MDRD eGFR (exponent 0.95, reference 90.5 mL/min/1.73 m^2); smoking, plasma nicotine and the CYP1A1, CYP1A2, CYP2A6, CYP3A4 and CYP3A5 polymorphisms were screened and not retained. SBP = E0 - Imax * Cc / (IC50 + Cc). The paper reports exponential inter-individual variability on V, CL, E0, Imax and IC50 but none of the variances, so all five etas are fixed to zero."
  reference <- paste(
    "Li G, Guan Y, Yang Y, Xiang Q, Chen S, Shao J, Chen Y, Yu X.",
    "Model-Informed Precision Dosing of Levamlodipine Besylate in Smoking Patients.",
    "Drug Des Devel Ther. 2025;19:3193-3207.",
    "doi:10.2147/DDDT.S501762"
  )
  vignette <- "Li_2025_levamlodipine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "levamlodipine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levamlodipine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    URIC_ACID = list(
      description = "Serum uric acid concentration",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value (blood biochemistry sampled on D1, Li 2025 Methods 'Study Design'). Power effect on CL/F: (URIC_ACID / 311.0)^-0.69 (Li 2025 final-model CL equation; 311.0 umol/L is the stated cohort mean and equals the Table 1 non-smoker mean). Cohort means (SD): non-smokers 311.0 (101.20), smokers 357.8 (116.70) umol/L (Table 1). A negative exponent means CL/F falls as uric acid rises; the abstract's 'lower eGFR and UA levels exhibited lower CL' reads the sign the other way, but the printed equation and Table 5 both give -0.69 and are followed here.",
      source_name = "UA"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD, BSA-normalised)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value. Power effect on CL/F: (CRCL / 90.5)^0.95 (Li 2025 final-model CL equation; 90.5 is the stated cohort mean, the Table 1 non-smoker mean is 90.46). Estimated with the Chinese MDRD form printed under Li 2025 Table 1: 186 x (Scr/88.4)^-1.154 x Age^-0.203 x 0.742 (if female), Scr in umol/L. Cohort means (SD): non-smokers 90.46 (17.60), smokers 106.45 (20.15) mL/min/1.73 m^2 (Table 1).",
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    SMOKE = list(
      description = "Current smoker indicator (1 = smoker, 0 = non-smoker)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on PK and PD parameters together with the plasma nicotine concentration (smokers only, three troughs per subject) and not retained (Li 2025 Results 'Covariates and Final Model'). 9 of 43 subjects smoked; all were male (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract and Results). Means 62.33 (non-smokers) and 64 (smokers) years (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained; coded male = 0, female = 1 in the source (Li 2025 Methods)."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Initially significant on CL in the forward step (Li 2025 Results) but absent from the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened and not retained (Li 2025 abstract)."
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 rs776746 genotype",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened as rs776746 CC = 0, TC = 1, TT = 2 and not retained (Li 2025 Methods and Results). CYP3A4 rs2740574, CYP1A1 rs4646903 and rs1048943, CYP1A2 rs762551 and CYP2A6 rs12459249 were also screened and not retained (Table 3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 43L,
    n_studies = 1L,
    n_observations = "125 steady-state trough levamlodipine concentrations (2 below the limit of detection excluded) and 129 systolic blood pressure values",
    age_range = "18-70 years by inclusion; means 62.33 (non-smokers) and 64 (smokers) years",
    weight_range = "means 61.08 (non-smokers) and 66.05 (smokers) kg",
    sex_female_pct = 62.8,
    race_ethnicity = "Chinese",
    disease_state = "Essential hypertension on long-term levamlodipine besylate monotherapy; 9 current smokers and 34 non-smokers",
    dose_range = "Levamlodipine besylate tablets orally once daily at the patient's fixed regimen; the Discussion simulates 2.5 and 5.0 mg daily",
    regions = "Chongqing, China (NCT05126381)",
    notes = "Prospective study; three steady-state trough concentrations and three SBP readings per patient on days 1-3, plus a baseline SBP. Baseline demographics from Li 2025 Table 1; genotypes from Table 3."
  )

  ini({
    lka <- fixed(log(0.79)); label("Absorption rate constant ka (1/h)") # Table 5 'tvKa' = 0.79; fixed to the literature value (Methods 'Pop PK / PD Model Construction')
    lvc <- log(3383); label("Apparent volume of distribution V/F (L)") # Table 5 'tvV' = 3383 L; final-model V equation
    lcl <- log(58.76); label("Apparent clearance CL/F at the reference uric acid and eGFR (L/h)") # Table 5 'tvCl' = 58.76 L/h; final-model CL equation

    e_uric_acid_cl <- -0.69; label("Power exponent of serum uric acid on CL/F (unitless)") # Table 5 'dCLdUA' = -0.69; final-model CL equation
    e_crcl_cl <- 0.95; label("Power exponent of eGFR on CL/F (unitless)") # Table 5 'dCLdeGFR' = 0.95; final-model CL equation

    le0 <- log(168); label("Baseline systolic blood pressure E0 (mmHg)") # Table 5 'tvE0' = 168; final-model E0 equation
    limax <- log(31.1); label("Maximal SBP reduction Imax (mmHg)") # Table 5 'tvImax' = 31.1; final-model Imax equation
    lic50 <- log(1.71); label("Levamlodipine concentration at half-maximal SBP reduction IC50 (ng/mL)") # Table 5 'tvIC50' = 1.71 ng/mL; final-model IC50 equation

    # Methods equations 1-5: exponential IIV on V, CL, E0, Imax and IC50; no
    # omega estimate is reported anywhere in the paper.
    etalvc ~ fixed(0)
    etalcl ~ fixed(0)
    etale0 ~ fixed(0)
    etalimax ~ fixed(0)
    etalic50 ~ fixed(0)

    propSd <- 0.06; label("Proportional residual error on levamlodipine concentration (fraction)") # Table 5 'stdev0' = 0.06
    propSd_sbp <- 0.06; label("Proportional residual error on SBP (fraction)") # Table 5 'stdev0' = 0.06, the single residual reported for the model
  })

  model({
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    # CL = 58.76 x (UA/311.0)^-0.69 x (eGFR/90.5)^0.95 x exp(nCL)
    cl <- exp(lcl + etalcl) * (URIC_ACID / 311.0)^e_uric_acid_cl * (CRCL / 90.5)^e_crcl_cl

    e0 <- exp(le0 + etale0)
    imax <- exp(limax + etalimax)
    ic50 <- exp(lic50 + etalic50)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # central in mg, vc in L -> mg/L; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    # Direct-effect inhibitory Imax model with baseline (Results 'Base Model')
    sbp <- e0 - imax * Cc / (ic50 + Cc)

    Cc ~ prop(propSd)
    sbp ~ prop(propSd_sbp)
  })
}
