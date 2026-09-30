Xiang_2020_bemarituzumab <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten elimination for bemarituzumab (anti-FGFR2b antibody) in adults with advanced solid tumours including gastric and gastroesophageal junction adenocarcinoma (Xiang 2020)"
  reference <- "Xiang H, Liu L, Gao Y, Ahene A, Macal M, Hsu AW, Dreiling L, Collins H. Population pharmacokinetic analysis of phase 1 bemarituzumab data to support phase 2 gastroesophageal adenocarcinoma FIGHT trial. Cancer Chemother Pharmacol. 2020;86(5):595-606. doi:10.1007/s00280-020-04139-4"
  vignette <- "Xiang_2020_bemarituzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Bemarituzumab was measured in SERUM by a validated ELISA (Methods,
  # 'Determination of serum concentration of bemarituzumab in humans').
  compartmentData <- list(
    central = list(analyte = "bemarituzumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "bemarituzumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL and Vc centred on 61 kg (Equations 3 and 4; the simulation 'typical",
        "population' is a 61 kg male with albumin 3.7 g/dL). Cohort median 61.4 kg, range 35.5-148 kg",
        "(Table 1). Also sets the mg/kg dose.",
        sep = " "
      ),
      source_name = "Weight"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper reports albumin in g/dL and centres the power effect on CL at 3.7 g/dL",
        "(Equation 3). The canonical column is g/L, converted inside model() via alb_gdL = ALB * 0.1.",
        "Cohort median 3.7 g/dL, range 1.9-4.6 g/dL (Table 1).",
        sep = " "
      ),
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Equation 4 enters 'female' as an exponential shift on Vc: Vc = exp(theta2 + ... +",
        "theta10 * female); theta10 = -0.191, so females have a 17.4% lower Vc. 42/75 (56%) female",
        "(Table 1).",
        sep = " "
      ),
      source_name = "female"
    )
  )

  # Covariates listed in Table 1 and screened on the CL and Vc etas but not
  # retained (Methods, 'Human population pharmacokinetic analysis'; Discussion:
  # 'Among 13 covariates evaluated (Table 1), only 3 had a significant impact').
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age; screened, not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 58 (range 25-86) years; Table 1."
    ),
    CRCL = list(
      description = "Creatinine clearance; screened, not retained.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 75.3 (range 26.4-200) mL/min; Table 1."
    ),
    FGFR2B_HIGH = list(
      description = paste(
        "FGFR2b overexpression (>= 10% of tumour cells with 3+ membranous IHC staining) in the",
        "gastric/GEJ subgroup; not a covariate for PK (Discussion; Supplementary Fig. 2).",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FGFR2b other)",
      notes = "26/53 (49%) FGFR2b high; Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 75L,
    n_studies = 1L,
    n_observations = 814L,
    age_range = "25-86 years",
    age_median = "58 years",
    weight_range = "35.5-148 kg",
    weight_median = "61.4 kg",
    sex_female_pct = 56,
    race_ethnicity = c(
      Asian = 58.7,
      White = 38.7,
      `American Indian or Alaska Native` = 1.33,
      `Black or African American` = 1.33
    ),
    disease_state = paste(
      "Advanced solid tumours: gastric and gastroesophageal junction adenocarcinoma (53/75, 70.7%)",
      "and other solid tumours (22/75, 29.3%).",
      sep = " "
    ),
    dose_range = paste(
      "0.3, 1, 3, 6, 10 and 15 mg/kg IV Q2W as a 30-min infusion (dose escalation, n = 27);",
      "15 mg/kg Q2W (expansion, n = 48).",
      sep = " "
    ),
    regions = "Phase 1 study FPA144-001 (NCT02318329)",
    albumin_range = "1.9-4.6 g/dL (median 3.7)",
    notes = paste(
      "Baseline covariates in Table 1. No patient developed post-dose anti-drug antibodies.",
      "Assay LLOQ 0.125 ug/mL.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: typical values for a 61 kg male with albumin 3.7 g/dL.
    lcl <- log(0.331); label("Linear clearance (L/day)") # Table 2, row 'Linear clearance, CL (L/day)' final = 0.331 (3.55%)
    lvc <- log(3.70); label("Central volume of distribution (L)") # Table 2, row 'Volume of central compartment, Vc (L)' final = 3.70 (2.98%)
    lq <- log(0.788); label("Intercompartmental clearance (L/day)") # Table 2, row 'Distribution clearance, Q (L/day)' final = 0.788 (8.99%)
    lvp <- log(2.05); label("Peripheral volume of distribution (L)") # Table 2, row 'Volume of peripheral compartment, Vp (L)' final = 2.05 (7.73%)
    # Vmax is printed as 1.70 ug/day (Table 2 and Results text) = 0.0017 mg/day in
    # this file's mg dosing units. Taken at its printed unit, the typical
    # Ctrough,ss of 123.8 ug/mL and the -6.2% female effect in Fig. 2 are both
    # reproduced; reading it as mg/day gives 118.8 ug/mL and -6.4%. See the
    # vignette section 'Units of Vmax'.
    lvmax <- log(0.0017); label("Maximum Michaelis-Menten elimination rate (mg/day)") # Table 2, row 'V max (ug/day)' final = 1.70 (14.3%)
    lkm <- log(4.58); label("Michaelis-Menten constant (ug/mL)") # Table 2, row 'K M (ug/mL)' final = 4.58 (15.1%)

    # Covariate effects (Equations 3 and 4 with ln() of the normalised covariate;
    # see the vignette section 'Form of the covariate equations').
    e_wt_cl <- 0.601; label("Power exponent of WT/61 on CL (unitless)") # Table 2, row 'Influence of body weight on CL' = 0.601 (19.5%)
    e_alb_cl <- -0.776; label("Power exponent of ALB/3.7 g/dL on CL (unitless)") # Table 2, row 'Influence of albumin on CL' = -0.776 (19.9%)
    e_wt_vc <- 0.303; label("Power exponent of WT/61 on Vc (unitless)") # Table 2, row 'Influence of body weight on V c' = 0.303 (24.8%)
    e_sexf_vc <- -0.191; label("Log-scale shift in Vc for females (unitless)") # Table 2, row 'Influence of sex on V c' = -0.191 (23.4%)

    # IIV. Table 2 prints the diagonals as percentages and the CL-Vc
    # off-diagonal as a bare covariance (0.0141); a covariance only exists on
    # the raw OMEGA scale, so the diagonals are read as 100 * sqrt(OMEGA):
    # variance = (percent / 100)^2.
    etalcl + etalvc ~ c(
      0.073441, # 0.271^2; Table 2, row 'Interindividual variability of CL' = 27.1 (19.7%)
      0.0141, # Table 2, row 'Covariance between CL and V c' = 0.0141 (46.1%), as printed
      0.029929 # 0.173^2; Table 2, row 'Interindividual variability of V c' = 17.3 (19.1%)
    )
    etalvp ~ 0.36 # 0.600^2; Table 2, row 'Interindividual variability of V p' = 60.0 (24.1%)
    etalvmax ~ 1.6384 # 1.28^2; Table 2, row 'Interindividual variability of V max' = 128 (23.3%)

    propSd <- 0.145; label("Proportional residual error (fraction)") # Table 2, row 'Residual variability (%CV)' = 14.5 (5.89%)
  })
  model({
    # Albumin: canonical column is g/L; the model was calibrated on g/dL.
    alb_gdL <- ALB * 0.1

    cl <- exp(lcl + e_wt_cl * log(WT / 61) + e_alb_cl * log(alb_gdL / 3.7) + etalcl)
    vc <- exp(lvc + e_wt_vc * log(WT / 61) + e_sexf_vc * SEXF + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Equations 1 and 2: parallel linear and Michaelis-Menten elimination from
    # the central compartment; Km is a concentration (ug/mL = mg/L).
    d/dt(central) <- -(vmax / (km + central / vc)) / vc * central - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
