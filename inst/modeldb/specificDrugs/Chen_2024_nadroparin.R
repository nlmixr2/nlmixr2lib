Chen_2024_nadroparin <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and first-order elimination for nadroparin (a low-molecular-weight heparin) in 40 preterm and term neonates and infants under 8 months of age treated for arterial or venous thromboembolic disease (Chen 2024). The model is fitted to plasma anti-Xa activity (IU/mL) rather than to a drug concentration, so the disposition parameters are apparent anti-Xa clearance (CL/F) and apparent anti-Xa volume of distribution (Vd/F). The single retained covariate is Schwartz-formula creatinine clearance on CL/F, entered as a median-normalised power model CL/F = 0.211 * (CRCL / 51.1)^0.238, consistent with nadroparin being cleared predominantly by renal excretion. Inter-individual variability was retained on CL/F only; the variances on Vd/F and ka approached zero during model building and were fixed to zero. Body-weight allometric scaling with a fixed 0.75 exponent and a sigmoid Emax postmenstrual-age maturation function were both tested but were not retained in the final model once creatinine clearance was included."
  reference <- paste(
    "Chen Y, Lan J, Zhu L, Dong M, Wang Y, Li Z.",
    "Is the current therapeutic dosage of nadroparin adequate for neonates and infants",
    "under 8 months with thromboembolic disease? a population pharmacokinetic study",
    "from a national children's medical center.",
    "Front Pharmacol. 2024 Jan 31;15:1331673.",
    "doi:10.3389/fphar.2024.1331673",
    sep = " "
  )
  vignette <- "Chen_2024_nadroparin"

  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/mL"
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the pediatric Schwartz formula, CLcr (mL/min/1.73 m^2) = k * HT (cm) / SCr (mg/dL), with k = 0.33 for preterm infants and k = 0.45 for term infants throughout the first year of life (Schwartz 1984)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Chen 2024 Methods 'Model development'",
        "equation (5) specifies the median-normalised power (exponential) covariate model",
        "P_i = P_p * (Cov_i / Cov_median)^theta; Table 2 gives the exponent theta = 0.238",
        "(RSE 29.2%) for CLCR on CL/F, and Table 1 gives the cohort CLCR median of",
        "51.1 mL/min/1.73 m^2 that supplies Cov_median. Adding CLCR to CL/F produced the",
        "largest reduction in the objective function value of all 17 screened covariates and",
        "the largest reduction in the CL/F inter-individual variance (Chen 2024 Results",
        "'Population PK modeling'), which the authors attribute to nadroparin being cleared",
        "predominantly by renal excretion (Discussion). Table 1 reports the cohort CLCR range",
        "as 5.3-314.7 mL/min/1.73 m^2 (mean 73.8, SD 62.8); over that range the covariate",
        "factor (CRCL / 51.1)^0.238 spans 0.58 to 1.54. Values are the Schwartz creatinine-based",
        "estimate, already normalised to 1.73 m^2 by construction of the formula, matching the",
        "canonical CRCL units without further conversion. The same Schwartz k * HT / SCr form",
        "is the registered CRCL alias used by MedellinGaribay_2015_gentamicin.R.",
        sep = " "
      ),
      source_name        = "CLCR (Chen 2024 Table 1 and Table 2 row 'CL/F_CLCR')"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight at the time of the anti-Xa sample",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chen 2024 Methods 'Model development' states that weight was selected a priori as a",
        "size descriptor and that clearance was scaled allometrically to a 70 kg adult with a",
        "fixed exponent of 0.75 (equation (3), CL/F = CLstd * (BW/70)^0.75). The Discussion",
        "reports that the allometric weight term did not survive backward elimination (the",
        "increase in the objective function value on its removal was not significant at",
        "p > 0.01), and Table 2 lists CL/F_CLCR as the only entry under 'Covariate model'.",
        "The final CL/F of 0.211 L/h is therefore an absolute typical value for this cohort,",
        "not a 70 kg-standardised CLstd: the Discussion's independently reported weight-",
        "normalised clearance of 0.068 L/h/kg equals 0.211 L/h divided by the Table 1 median",
        "body weight of 3.1 kg, which confirms that no allometric term remains in the final",
        "model. Table 1 cohort weight: mean (SD) 3.7 (2.0) kg, median 3.1 kg, range 1.2-9.3 kg.",
        sep = " "
      ),
      source_name        = "BW (Chen 2024 Table 1)"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "months (canonical); Chen 2024 reports postmenstrual age in weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chen 2024 Methods 'Model development' equation (4) tested a sigmoid Emax maturation",
        "function of postmenstrual age on clearance, parameterised by TM50 (the postmenstrual",
        "age at which clearance is 50% of the mature value) and the Hill coefficient gamma.",
        "The Discussion reports that postmenstrual age was statistically significant when",
        "tested on its own but lost significance once creatinine clearance was in the model,",
        "which the authors attribute to the high correlation between creatinine clearance and",
        "maturation-related covariates in pediatric patients. Neither TM50 nor gamma is",
        "reported anywhere in the paper, so the maturation function is not reproducible even",
        "as a documented alternative. Table 1 cohort PMA: mean (SD) 43.5 (9.5) weeks,",
        "median 40.4 weeks, range 30.7-69.0 weeks.",
        sep = " "
      ),
      source_name        = "PMA (Chen 2024 Table 1)"
    ),
    PNA = list(
      description        = "Postnatal (chronological) age since birth",
      units              = "months (canonical); Chen 2024 reports postnatal age in days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a candidate covariate on clearance. Chen 2024 Discussion reports that",
        "postnatal age was statistically significant when tested individually but showed no",
        "significant effect once creatinine clearance was included. No point estimate is",
        "reported. Table 1 cohort PNA: mean (SD) 58.4 (60.4) days, median 40.0 days,",
        "range 3.0-224.0 days.",
        sep = " "
      ),
      source_name        = "PNA (Chen 2024 Table 1)"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the covariate analysis (Chen 2024 Methods 'Model development') and not",
        "retained; no point estimate is reported. Table 1 cohort GA: mean (SD) 35.2 (4.3)",
        "weeks, median 36.8 weeks, range 25.0-41.3 weeks, with 1 subject (2.5%) below 28",
        "weeks, 9 (22.5%) at 28-32 weeks, 10 (25%) at 32-37 weeks and 20 (50%) at 37 weeks",
        "or above.",
        sep = " "
      ),
      source_name        = "GA (Chen 2024 Table 1)"
    ),
    WT_BIRTH = list(
      description        = "Birth body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the covariate analysis and not retained; no point estimate is reported. Table 1 cohort BBW: mean (SD) 2.6 (1.1) kg, median 2.9 kg, range 0.6-4.6 kg.",
      source_name        = "BBW (Chen 2024 Table 1)"
    ),
    HT = list(
      description        = "Body length / height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the covariate analysis and not retained as a covariate in its own right; height nevertheless enters the model indirectly because the Schwartz creatinine clearance CRCL is computed as k * HT / SCr. Table 1 cohort HT: mean (SD) 49.0 (9.9) cm, median 50.0 cm, range 28.0-70.0 cm.",
      source_name        = "HT (Chen 2024 Table 1)"
    ),
    BSA = list(
      description        = "Body surface area, computed as sqrt(HT (cm) * BW (kg) / 3600)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the covariate analysis and not retained; no point estimate is reported. Table 1 cohort BSA: mean (SD) 0.2 (0.1) m^2, median 0.2 m^2, range 0.1-0.4 m^2.",
      source_name        = "BSA (Chen 2024 Table 1)"
    ),
    SEXF = list(
      description        = "Sex; 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened as 'gender' in the covariate analysis via the dichotomous-covariate model of Chen 2024 equation (6) and not retained; no point estimate is reported. Table 1 cohort: 23 boys and 17 girls (42.5% female).",
      source_name        = "Gender (Chen 2024 Table 1)"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the covariate analysis and not retained as a covariate in its own right; serum creatinine nevertheless enters the model indirectly because the Schwartz creatinine clearance CRCL is computed as k * HT / SCr. Table 1 cohort SCR: mean (SD) 42.4 (58.1) umol/L, median 29.7 umol/L, range 6.7-374.8 umol/L. Note that the Schwartz formula as printed in Chen 2024 Methods requires SCr in mg/dL, whereas Table 1 tabulates SCR in umol/L (1 mg/dL = 88.4 umol/L).",
      source_name        = "SCR (Chen 2024 Table 1)"
    ),
    CYSC = list(
      description        = "Serum cystatin C",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as an alternative renal-function marker and not retained; creatinine clearance produced the larger objective-function reduction. Table 1 cohort CysC: mean (SD) 1.2 (0.4) mg/L, median 1.1 mg/L, range 0.7-3.0 mg/L.",
      source_name        = "CysC (Chen 2024 Table 1)"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the covariate analysis and not retained; no point estimate is reported. Table 1 cohort ALB: mean (SD) 32.3 (5.0) g/L, median 33.4 g/L, range 21.8-44.1 g/L.",
      source_name        = "ALB (Chen 2024 Table 1)"
    ),
    ALT = list(
      description        = "Alanine transaminase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a hepatic-function marker and not retained; no point estimate is reported. Table 1 cohort ALT: mean (SD) 61.9 (120.5) U/L, median 17.7 U/L, range 3.2-429.3 U/L.",
      source_name        = "ALT (Chen 2024 Table 1)"
    ),
    AST = list(
      description        = "Aspartate transaminase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a hepatic-function marker and not retained; no point estimate is reported. Table 1 cohort AST: mean (SD) 71.8 (68.4) U/L, median 41.7 U/L, range 16.1-277.1 U/L.",
      source_name        = "AST (Chen 2024 Table 1)"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a hepatic-function marker and not retained; no point estimate is reported. Table 1 cohort TBIL: mean (SD) 53.2 (65.3) umol/L, median 24.6 umol/L, range 1.8-299.6 umol/L.",
      source_name        = "TBIL (Chen 2024 Table 1)"
    ),
    DBIL = list(
      description        = "Direct (conjugated) serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a cholestasis marker alongside total bilirubin and not retained; no point estimate is reported. Table 1 cohort DBIL: mean (SD) 14.6 (39.8) umol/L, median 7.0 umol/L, range 0.5-250.2 umol/L.",
      source_name        = "DBIL (Chen 2024 Table 1)"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a renal-function marker and not retained; creatinine clearance produced the larger objective-function reduction. Table 1 cohort BUN: mean (SD) 3.2 (2.7) mmol/L, median 2.2 mmol/L, range 0.4-12.0 mmol/L.",
      source_name        = "BUN (Chen 2024 Table 1)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "nadroparin",
      units    = "IU (anti-Xa activity)",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "nadroparin",
      units    = "IU (anti-Xa activity)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    n_observations = "56 anti-Xa samples from 40 patients (Chen 2024 Table 1). 51 patients were enrolled; 11 were excluded because every one of their samples was below the limit of quantification.",
    age_range      = "postnatal age 3-224 days (median 40 days); postmenstrual age 30.7-69.0 weeks (median 40.4 weeks); gestational age at birth 25.0-41.3 weeks (median 36.8 weeks)",
    age_median     = "postnatal age 40 days (postmenstrual age 40.4 weeks)",
    weight_range   = "1.2-9.3 kg",
    weight_median  = "3.1 kg",
    sex_female_pct = 42.5,
    race_ethnicity = "Chinese (single centre in Shanghai, China; race and ethnicity are not tabulated in the paper)",
    disease_state  = "Preterm or term neonates and infants under 8 months of age with suspected or diagnosed arterial or venous thrombosis. Table 1: 23 patients had venous thrombosis (lower limb 7, external iliac vein 5, inferior vena cava 4, femoral vein 2, portal vein 2, jugular vein 2, umbilical vein 1) and 17 had arterial thrombosis (neonatal cerebral infarction 10, arterial thromboembolism 6, thrombosis of the abdominal aorta 1).",
    dose_range     = "Nadroparin (Fraxiparine) 150-200 IU/kg subcutaneously every 12 h per the local treatment protocol",
    regions        = "China (Children's Hospital of Fudan University, National Children's Medical Center, Shanghai; single centre)",
    renal_function = "Schwartz creatinine clearance 5.3-314.7 mL/min/1.73 m^2 (median 51.1, mean 73.8, SD 62.8)",
    notes          = paste(
      "Retrospective single-centre chart review of patients treated between July 2021 and",
      "December 2023 (Chen 2024 Methods 'Patients and data collection'). Anti-Xa activity was",
      "measured with an anti-Xa clotting assay (STA-liquid ANTI-Xa, Diagnostica Stago) with a",
      "lower limit of quantification of 0.1 IU/mL and a calibration range of 0.1-2.00 IU/mL.",
      "Sampling was sparse and, importantly, essentially all samples were peak samples drawn",
      "4 h after a dose (Methods 'Blood sampling and anti-Xa determination'; the Discussion",
      "states explicitly that all samples obtained were peak concentrations), so the absorption",
      "and terminal phases of the profile are only weakly informed by the data. Samples below",
      "the limit of quantification were discarded rather than modelled; a sensitivity analysis",
      "using the M3 method is summarised in the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # Structural model: one compartment with first-order absorption from a
    # subcutaneous depot and first-order elimination, fitted to plasma anti-Xa
    # activity. Chen 2024 Results 'Population PK modeling'; estimates and the
    # bootstrap validation are in Table 2.
    #
    # Because nadroparin is dosed subcutaneously and no absolute bioavailability
    # was estimated, clearance and volume are apparent (CL/F, Vd/F) and no
    # separate bioavailability parameter appears in the model.
    # ========================================================================
    lka <- log(0.495); label("Absorption rate constant ka (1/h)")                                                  # Chen 2024 Table 2: ka = 0.495 1/h (RSE 16.4%; bootstrap median 0.495, 95% CI 0.073-1.09)
    lcl <- log(0.211); label("Apparent anti-Xa clearance CL/F (L/h) at the cohort median CRCL of 51.1 mL/min/1.73 m^2")  # Chen 2024 Table 2: CL/F = 0.211 L/h (RSE 9.4%; bootstrap median 0.207, 95% CI 0.063-0.259)
    lvc <- log(1.55);  label("Apparent anti-Xa volume of distribution Vd/F (L)")                                   # Chen 2024 Table 2: Vd/F = 1.55 L (RSE 13.7%; bootstrap median 1.51, 95% CI 0.02-2.79)

    # ---- Retained covariate effect on CL/F ---------------------------------
    # Chen 2024 equation (5) is the median-normalised power model
    #     P_i = P_p * (Cov_i / Cov_median)^theta
    # and Table 2 reports theta = 0.238 for CLCR on CL/F. The normalising
    # constant Cov_median is the Table 1 cohort CLCR median of
    # 51.1 mL/min/1.73 m^2, so the retained relationship is
    #     CL/F = 0.211 * (CRCL / 51.1)^0.238.
    e_crcl_cl <- 0.238; label("Exponent for creatinine clearance on apparent clearance (unitless)")                # Chen 2024 Table 2 row 'CL/F_CLCR': 0.238 (RSE 29.2%; bootstrap median 0.235, 95% CI 0.056-0.381)

    # ========================================================================
    # Inter-individual variability
    # Chen 2024 equation (1) specifies an exponential (log-normal) IIV model,
    # P_i = P_p * exp(eta_i). Table 2 reports the IIV as a percent coefficient
    # of variation and its footnote defines CV as the coefficient of variation
    # of the parameter values, i.e. the CV of the log-normally distributed
    # individual parameter, so the internal variance is omega^2 = log(1 + CV^2).
    #
    # Chen 2024 Results 'Population PK modeling' states that the variance for
    # volume distribution and the absorption rate constant were fixed to zero,
    # as they approached zero during the model-building process. A variance of
    # exactly zero is mathematically identical to the absence of the random
    # effect, and encoding it as etalvc ~ fixed(0) would make the omega matrix
    # singular, so no eta is declared on lvc or lka. Table 2 accordingly reports
    # an inter-individual variability entry for CL/F only. This provenance is
    # recorded in the vignette Assumptions and deviations section.
    # ========================================================================
    etalcl ~ 0.0678689     # Chen 2024 Table 2: omega_CL/F CV 26.5% (RSE 29.6%; bootstrap median 24, 95% CI 0.7-41.6); log(1 + 0.265^2) = 0.0678689

    # ========================================================================
    # Residual unexplained variability
    # Chen 2024 equation (2) writes the general combined form
    # OBS_i = IPRED * exp(eps_1) + eps_2, but Results 'Population PK modeling'
    # states that residual variability was best described by a proportional
    # error model, and Table 2 reports a single proportional residual error
    # with no additive term. For a proportional model the tabulated %CV is the
    # standard deviation of the proportional term, which is what prop() takes.
    # ========================================================================
    propSd <- 0.355; label("Proportional residual error on anti-Xa activity (fraction)")                           # Chen 2024 Table 2: proportional residual error 35.5% CV (RSE 10.0%; bootstrap median 34.6, 95% CI 27.2-41.5)
  })

  model({
    # ---- Individual PK parameters ------------------------------------------
    # Chen 2024 equation (5) with Cov_median = 51.1 mL/min/1.73 m^2, the
    # Table 1 cohort median creatinine clearance:
    #     CL/F = 0.211 * (CRCL / 51.1)^0.238
    # No allometric weight term and no postmenstrual-age maturation term appear
    # here: both were tested (equations (3) and (4)) but neither was retained in
    # the final model, and Table 2 lists CL/F_CLCR as the only covariate effect.
    # No eta on ka or vc: Chen 2024 fixed both variances to zero (see ini()).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 51.1)^e_crcl_cl
    vc <- exp(lvc)

    kel <- cl / vc

    # ---- ODE system --------------------------------------------------------
    # Subcutaneous dose enters the depot; first-order absorption into the
    # central compartment and first-order elimination from it.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and residual error ------------------------------------
    # The depot and central amounts carry anti-Xa activity in IU and vc is in
    # L, so central / vc is IU/L. Divide by 1000 to report IU/mL, the unit in
    # which the anti-Xa assay was calibrated and in which the 0.5-1 IU/mL
    # therapeutic target is expressed (Chen 2024 Methods 'Blood sampling and
    # anti-Xa determination').
    Cc <- central / vc / 1000
    Cc ~ prop(propSd)
  })
}
