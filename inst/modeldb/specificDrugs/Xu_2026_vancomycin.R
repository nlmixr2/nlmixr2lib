Xu_2026_vancomycin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for vancomycin in children aged 3 months to",
    "18 years (Australian multi-centre data). All four disposition parameters are allometrically",
    "scaled to a 70 kg adult (exponent 0.75 on CL and Q, 1 on V1 and V2); clearance additionally",
    "carries a Rhodin sigmoid renal-maturation function of postmenstrual age (TM50 47.7 weeks,",
    "Hill 3.4, both fixed) and an exponential serum-creatinine effect on the deviation of serum",
    "creatinine from the Ceriotti age-expected value. Exponential inter-individual variability on",
    "CL and V1, inter-occasion variability on CL, and a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Xu H, Yang W, Lai T, McMullan B, Yeoh D, Wilkins A, Goff Z, Zhu X, Gwee A.",
    "Model-informed precision dosing of vancomycin in children 3 months to 18 years of age using",
    "Australia-wide data. Antimicrob Agents Chemother. 2026. doi:10.1128/aac.01840-25.",
    "Final estimates from Table 2; structural equations from Methods Equations 1-6 and the Table 2",
    "footnote equations; model-development history from supplemental Table S2.",
    sep = " "
  )
  vignette <- "Xu_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling to a 70 kg adult (Methods 'Population pharmacokinetic modeling';",
        "Table 2 footnote equations): (WT/70)^0.75 on CL and Q, (WT/70)^1 on V1 and V2.",
        "Development-cohort median 22.9 kg (range 2.70-115; Table 1).",
        sep = " "
      ),
      source_name = "weight"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "WEEKS, not the register default of months, because the Rhodin 2009 maturation constants",
        "(TM50 47.7 weeks, Hill 3.4) are on that scale (Table 2 'TM50 (weeks)'). Drives the sigmoid",
        "renal-maturation factor on CL, PMA^3.4 / (PMA^3.4 + 47.7^3.4) (Methods Equation 3; Table 2",
        "footnote). Development-cohort median 397 weeks (range 39.4-998; Table 1). The paper does",
        "not state how PMA was derived when gestational age was unknown; for a term-born child",
        "PMA = 40 weeks + postnatal age in weeks.",
        sep = " "
      ),
      source_name = "PMA"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only the Ceriotti age-expected serum creatinine,",
        "SCRstd = -2.37330 - 12.91367 * ln(AGE) + 23.93581 * AGE^0.5 (Methods Equation 6), in",
        "years. Development-cohort median 6.80 years (range 0.20-18.4; Table 1).",
        sep = " "
      ),
      source_name = "age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Exponential effect on CL of the deviation from the age-expected value,",
        "F_renal = exp(-0.0188 * (CREAT - SCRstd)) (Methods Equation 5; Table 2). Unit umol/L per",
        "Table 1 and the Table S2 footnote ('SCR: serum creatinine (umol/L)'). Measured by",
        "enzymatic assays at all four centres (Table S1). Development-cohort median 31.0 umol/L",
        "(range 10.0-328; Table 1).",
        sep = " "
      ),
      source_name = "SCR"
    ),
    OCC = list(
      description = "Occasion index (treatment episode) for inter-occasion variability on CL",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Methods: IOV in CL 'to account for changes in drug clearance across multiple treatment",
        "episodes'. The paper does not report how many occasions were defined; four are provided",
        "(OCC = 1-4). Any other OCC value contributes no IOV.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  # Covariates Xu 2026 screened by stepwise covariate modelling and did not retain
  # ('After SCM, there were no significant covariates to be included in the model').
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 'Population pharmacokinetic modeling'); not retained after SCM (Results)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Methods); not retained after SCM (Results). Table 1 median 31.0 g/L."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened as 'Urea' (Methods); not retained after SCM (Results). Table 1 median 3.80 mmol/L."
    ),
    DIS_CRITILL = list(
      description = "Intensive-care-unit admission indicator (1 = admitted to ICU)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened as 'ICU admission' (Methods); not retained after SCM (Results). 27.8% of the",
        "development cohort. The paper also screened a central-line indicator, which has no",
        "register canonical and is recorded here only: not retained.",
        sep = " "
      )
    ),
    PNA = list(
      description = "Postnatal age for the O'Hanlon birth-transition maturation factor (days in the source)",
      units = "months",
      type = "continuous",
      notes = paste(
        "The 'updated' maturation function with a PNA birth-transition factor (Methods Equation 4;",
        "Table S2 PNAMAX 0.588, PNAT50 6.94 days, fixed from O'Hanlon 2023) was evaluated and",
        "rejected in favour of the Rhodin PMA function (Table S2: AIC 12569.38 vs 12682.54).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 626L,
    n_studies = 1L,
    n_concentrations = 3687L,
    age_range = "0.20-18.4 years (3 months to 18 years by inclusion criteria)",
    age_median = "6.80 years",
    weight_range = "2.70-115 kg",
    weight_median = "22.9 kg",
    pma_range = "39.4-998 weeks",
    sex_female_pct = 45.0,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Hospitalised children receiving intravenous vancomycin for suspected or proven gram-positive",
      "infection; 27.8% admitted to intensive care. Patients on ECMO, peritoneal dialysis or renal",
      "replacement therapy were excluded.",
      sep = " "
    ),
    renal_function = "Serum creatinine median 31.0 umol/L (range 10.0-328)",
    dose_range = paste(
      "Standard 15 mg/kg every 6 h as 1-h intermittent infusions, or continuous infusion",
      "(Table S1); routine therapeutic drug monitoring.",
      sep = " "
    ),
    regions = paste(
      "Australia. Model development: The Royal Children's Hospital, Melbourne (March 2020 to May",
      "2023). External validation: Perth Children's Hospital (35 children), Sydney Children's",
      "Hospital (17) and The Children's Hospital at Westmead (170).",
      sep = " "
    ),
    notes = paste(
      "Table 1 (development cohort, RCH). The Results text gives 627 patients; Table 1 lists 626.",
      "Median 3 concentrations per child (range 1-57).",
      sep = " "
    )
  )

  ini({
    # Structural parameters for a 70 kg, fully mature subject whose serum
    # creatinine equals the Ceriotti age-expected value.
    lcl <- log(7.56); label("Clearance for a 70 kg mature subject (L/h)") # Table 2 'CL (L/h/70 kg)' = 7.56 (RSE 2.8%)
    lvc <- log(85.6); label("Central volume for a 70 kg subject (L)") # Table 2 'V1 (L/70 kg)' = 85.6 (RSE 0.9%)
    lq <- log(1.30); label("Intercompartmental clearance for a 70 kg subject (L/h)") # Table 2 'Q (L/h/70 kg)' = 1.30 (RSE 14.7%)
    lvp <- log(286); label("Peripheral volume for a 70 kg subject (L)") # Table 2 'V2 (L/70 kg)' = 286 (RSE 31.2%)

    # Allometric exponents, fixed (Results 'CL, Q: 0.75; V1, V2: 1.0'; Table 2
    # footnote equations).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/70) for CL and Q (unitless)") # Results; Table 2 footnote
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/70) for V1 and V2 (unitless)") # Results; Table 2 footnote

    # Renal maturation (Rhodin 2009, paper reference 26), fixed.
    hill_mat <- fixed(3.4); label("Hill coefficient of renal maturation (unitless)") # Table 2 'Hill' = 3.4 FIX
    tmat50 <- fixed(47.7); label("Postmenstrual age at 50 percent renal maturation (weeks)") # Table 2 'TM50 (weeks)' = 47.7 FIX

    # Serum creatinine effect on CL: exp(-e_creat_cl * (CREAT - SCRstd)).
    e_creat_cl <- 0.0188; label("Exponential coefficient of serum creatinine excess on CL (per umol/L)") # Table 2 'Effect of creatinine on CL' = 0.0188 (RSE 3.9%)

    # IIV. Table S2 footnote: CV% = sqrt(omega^2) x 100%, so variance = (CV/100)^2.
    etalcl ~ 0.048841 # Table 2 IIV CL 22.1% -> 0.221^2
    etalvc ~ 0.108241 # Table 2 IIV V1 32.9% -> 0.329^2

    # IOV on CL, 18.7% -> 0.187^2 = 0.034969. Occasions after the first share
    # the variance (NONMEM BLOCK SAME analogue) and are fixed().
    etaiov_cl_1 ~ 0.034969 # Table 2 IOV CL 18.7%
    etaiov_cl_2 ~ fixed(0.034969)
    etaiov_cl_3 ~ fixed(0.034969)
    etaiov_cl_4 ~ fixed(0.034969)

    propSd <- 0.233; label("Proportional residual error (fraction)") # Table 2 'Proportional error (%)' = 23.3%; Table S2 'Prop.error' = 0.233
  })

  model({
    # Ceriotti age-expected serum creatinine (umol/L), age in years (Methods
    # Equation 6). The Table 2 footnote prints the constant as +2.37330; Equation
    # 6 and the Ceriotti equation as used elsewhere carry -2.37330.
    creat_std <- -2.37330 - 12.91367 * log(AGE) + 23.93581 * AGE^0.5
    f_renal <- exp(-e_creat_cl * (CREAT - creat_std))

    # Sigmoid renal maturation on postmenstrual age in weeks (Equation 3).
    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    iov_cl <- (OCC == 1) * etaiov_cl_1 + (OCC == 2) * etaiov_cl_2 +
      (OCC == 3) * etaiov_cl_3 + (OCC == 4) * etaiov_cl_4

    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl_q * fmat * f_renal
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
