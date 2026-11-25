Kim_2018_doripenem <- function() {
  description <- "Two-compartment IV population PK model for doripenem in 11 healthy Korean adult volunteers after a single 250 mg 1-hour infusion (Kim 2018). Body weight enters every structural parameter through fixed allometric exponents (0.75 on CL and Q, 1 on Vc and Vp) normalised to the cohort median of 57 kg; no other covariate was retained. Correlated log-normal IIV on CL and Vc and a proportional residual error."
  reference <- paste(
    "Kim SW, Choe S, Lee DH, Jang HJ, Lee SJ, Hong JH, Yang HJ, Jang S,",
    "Park JY, Kim A, Kim YK, Lee D.",
    "Pharmacokinetics of doripenem in healthy Koreans and Monte Carlo",
    "simulations to explore optimal dosage regimens in patients with",
    "normal and enhanced renal function.",
    "Ther Drug Monit. 2018;40(4):425-433.",
    "doi:10.1097/FTD.0000000000000528.",
    sep = " "
  )
  vignette <- "Kim_2018_doripenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "doripenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doripenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling theta = theta_TV * (WT / WT_MED)^k with WT_MED",
        "= 57 kg, the cohort median (Kim 2018 Results, allometric equation;",
        "Table 1 median 57.0 kg, range 50-80). Exponents fixed at 0.75 on",
        "CL and Q and 1 on Vc and Vp because the cohort was too small to",
        "estimate them (Table 3 footnote). Individual values in",
        "Supplementary Table 1."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL, Q, Vc and Vp in stepwise covariate modelling; not retained (Kim 2018 Methods / Results)."
    ),
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Tested on CL, Q, Vc and Vp; not retained (Kim 2018 Methods / Results)."
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL; not retained (Kim 2018 Methods / Results). Cohort range 101-139 mL/min (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL; not retained (Kim 2018 Methods / Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    age_range = "22-30 years",
    age_median = "24 years",
    weight_range = "50-80 kg",
    weight_median = "57 kg",
    sex_female_pct = 54.5,
    race_ethnicity = "Korean (100%)",
    disease_state = "Healthy adult volunteers (Cockcroft-Gault CLCR 101-139 mL/min)",
    dose_range = "250 mg doripenem single IV infusion over 1 hour",
    regions = "Republic of Korea (Pusan National University Hospital)",
    n_concentrations = 132L,
    crcl_range = "101-139 mL/min (Cockcroft-Gault; median 123)",
    notes = paste(
      "Kim 2018 Table 1 and Supplementary Table 1: 12 subjects enrolled",
      "(6 men, 6 women); one man did not attend, leaving 5 men and 6 women.",
      "Twelve plasma samples per subject at 0 (pre-dose), 0.33, 0.67, 1,",
      "1.5, 2, 2.5, 3, 4, 5, 6 and 8 h after the start of the infusion.",
      "UHPLC-MS/MS assay, LLOQ 0.01 mg/L. NONMEM 7.3 FOCE-I."
    )
  )

  ini({
    # Structural parameters for a 57 kg subject (Kim 2018 Table 3 'Estimates').
    lcl <- log(22.0); label("Clearance at WT = 57 kg (L/h)")                     # Table 3: CL = 22.0 L/h (RSE 10.2%)
    lvc <- log(16.0); label("Central volume of distribution at WT = 57 kg (L)")  # Table 3: VC = 16.0 L (RSE 13.0%)
    lq  <- log(1.83); label("Intercompartmental clearance at WT = 57 kg (L/h)")  # Table 3: Q = 1.83 L/h (RSE 22.3%)
    lvp <- log(3.13); label("Peripheral volume of distribution at WT = 57 kg (L)") # Table 3: VP = 3.13 L (RSE 14.7%)

    # Allometric exponents, fixed (Results; Table 3 footnote final-model equations).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL and Q (unitless)") # Table 3 footnote: CL, Q x (WT/57)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent of WT on Vc and Vp (unitless)") # Table 3 footnote: VC, VP x (WT/57)

    # IIV reported as CV% (Table 3); omega^2 = log(1 + CV^2).
    # CL: CV 31.6% -> 0.09518; Vc: CV 35.3% -> 0.11796;
    # covariance = 0.871 * sqrt(0.09518 * 0.11796) = 0.09218.
    etalcl + etalvc ~ c(0.09518, 0.09218, 0.11796) # Table 3: omega_CL 31.6%, omega_V 35.3%, correlation 0.871

    # Residual error: Y = PRED + PRED * eps (Results); Table 3 sigma_Proportional = 0.301 read as an SD.
    propSd <- 0.301; label("Proportional residual error (fraction)") # Table 3: sigma_Proportional = 0.301 (RSE 7.63%)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 57)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 57)^e_wt_vc
    q  <- exp(lq) * (WT / 57)^e_wt_cl
    vp <- exp(lvp) * (WT / 57)^e_wt_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
