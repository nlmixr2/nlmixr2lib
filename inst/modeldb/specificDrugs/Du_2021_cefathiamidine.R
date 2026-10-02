Du_2021_cefathiamidine <- function() {
  description <- "One-compartment IV population PK model for cefathiamidine in infants 0.35-1.86 years with augmented renal clearance (Du 2021), with fixed allometric body-weight scaling on CL and V (reference 10.25 kg) and a power-form age effect on CL (reference 1.25 years)."
  reference <- paste(
    "Du B, Zhou Y, Tang BH, Wu YE, Yang XM, Shi HY, Yao BF, Hao GX, You DP,",
    "van den Anker J, Zheng Y, Zhao W. Population Pharmacokinetic Study of",
    "Cefathiamidine in Infants With Augmented Renal Clearance.",
    "Front Pharmacol. 2021;12:630047.",
    "doi:10.3389/fphar.2021.630047.",
    sep = " "
  )
  vignette <- "Du_2021_cefathiamidine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefathiamidine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight on the day of first sampling",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL (fixed exponent 0.75) and V (fixed exponent 1) with reference 10.25 kg, the cohort median current weight on the day of first sampling (Du 2021 Table 2 formulae and footnote).",
      source_name = "CW"
    ),
    AGE = list(
      description = "Postnatal age in years on the day of first sampling",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL, F_age = (AGE / 1.25)^theta_3, with reference 1.25 years, the cohort median age on the day of first sampling (Du 2021 Table 2 formula and footnote). Supported range 0.35-1.86 years.",
      source_name = "AGE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "0.35-1.86 years",
    age_median = "1.25 years",
    weight_range = "8.0-13.0 kg",
    weight_median = "10.25 kg",
    sex_female_pct = 50,
    race_ethnicity = "Chinese (Du 2021 Table 1)",
    disease_state = "Infants (<= 2 years) with augmented renal clearance (Schwartz eGFR >= 130 mL/min/1.73 m^2) and hematologic disease (immune thrombocytopenia 6, leukemia 3, anemia 3, infectious mononucleosis syndrome 2, agranulocytosis 2, other 4) receiving cefathiamidine for suspected or confirmed bacterial infection",
    dose_range = "100 mg/kg/day cefathiamidine IV q12h as a 30-min infusion (median 50 mg/kg/dose, range 40-100; median 500 mg/dose, range 400-1000)",
    regions = "China (Children's Hospital of Hebei Province, Shijiazhuang; single centre)",
    renal_function = "eGFR (Schwartz) median 197 mL/min/1.73 m^2 (range 132-413); serum creatinine median 20 umol/L (range 10-26)",
    notes = "Baseline demographics per Du 2021 Table 1. 36 scavenged plasma samples (0.15-222 ug/mL, all above the 30 ng/mL LLOQ of the UPLC-MS/MS assay). NONMEM 7.4, FOCE-I. The Results text quotes a weight range of 8.0-12.5 kg while Table 1 prints 8.00-13.00 kg; the table value is recorded here."
  )

  ini({
    # Structural parameters at the reference infant (WT = 10.25 kg, AGE = 1.25 years);
    # Du 2021 Table 2 'Full dataset' final estimates.
    lcl <- log(2.20); label("Clearance at 10.25 kg and 1.25 years (L/h)") # Du 2021 Table 2: theta_1 = 2.20 (RSE 8.30%)
    lvc <- log(3.36); label("Volume of distribution at 10.25 kg (L)") # Du 2021 Table 2: theta_2 = 3.36 (RSE 8.2%)

    # Allometric exponents fixed a priori (Du 2021 Methods and Results, Covariate
    # Analysis: 'fixed allometric exponents of 0.75 and 1 for CL and V').
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Du 2021 Results, Covariate Analysis
    e_wt_vc <- fixed(1); label("Allometric exponent on V (unitless)") # Du 2021 Results, Covariate Analysis; Table 2 V formula

    # Power-form age effect on CL: F_age = (AGE / 1.25)^theta_3.
    e_age_cl <- 0.662; label("Power exponent of (AGE/1.25) on CL (unitless)") # Du 2021 Table 2: theta_3 = 0.662 (RSE 21.6%)

    # Inter-individual variability, exponential (Du 2021 Methods: theta_i = theta * exp(eta_i)).
    # Table 2 reports CV%; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.063478 # log(0.256^2 + 1); Du 2021 Table 2 IIV CL = 25.6% (shrinkage 15.1%)
    etalvc ~ 0.048958 # log(0.224^2 + 1); Du 2021 Table 2 IIV V = 22.4% (shrinkage 28.5%)

    # Exponential residual error (Du 2021 Results: 'An exponential model best
    # described residual variability'), Y = IPRED * exp(eps) -> lnorm(). Table 2
    # ERR(1) = 22.6% as CV; SD on the log scale = sqrt(log(0.226^2 + 1)).
    expSd <- 0.223191; label("Lognormal residual error (SD on the log scale)") # Du 2021 Table 2: ERR(1) = 22.6% (shrinkage 35.4%)
  })
  model({
    # Du 2021 Table 2: CL = theta_1 x (CW/10.25)^0.75 x F_age; F_age = (AGE/1.25)^theta_3;
    # V = theta_2 x (CW/10.25).
    cl <- exp(lcl + etalcl) * (WT / 10.25)^e_wt_cl * (AGE / 1.25)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 10.25)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
