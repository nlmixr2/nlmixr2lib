Li_2021_vancomycin_normal <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Chinese infants (1-24 months) with normal renal function, eGFR 30-86 mL/min/1.73 m^2 (Li 2021, Model 1). CL scales as a power of body weight (reference 2.25 kg, estimated exponent 1.24) and exponentially with serum creatinine (exp(-0.533 * SCR / 27.1)); V scales as a power of body weight (reference 2.25 kg, estimated exponent 1.28). IIV on CL only; proportional residual error."
  reference <- "Li DY, Li L, Li GZ, Hu YH, Guo HL, Jing X, Chen F, Ji X, Xu J, Dai HR. Population Pharmacokinetics Modeling of Vancomycin Among Chinese Infants With Normal and Augmented Renal Function. Front Pediatr. 2021;9:713588. doi:10.3389/fped.2021.713588"
  vignette <- "Li_2021_vancomycin_renal_function"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 2 (Model 1 group): median 2.25 kg (range 1.15-13). Power-scaled on CL and V with reference 2.25 kg (Table 5 header equation and Results 'Model 1' equation).",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 2 (Model 1 group): median 27.1 umol/L (range 14-315). Enters CL as exp(e_creat_cl * CREAT / 27.1) (Table 5 header equation 'e**[theta5*SCR/27.1]'). This form is NOT centred: at CREAT = 27.1 umol/L the factor is exp(-0.533) = 0.587, not 1.",
      source_name = "SCR"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 61L,
    n_studies = 1L,
    age_range = "1-24 months postnatal age",
    age_median = "1 month (mean 2.31 months)",
    weight_range = "1.15-13 kg",
    weight_median = "2.25 kg",
    sex_female_pct = 39.3,
    race_ethnicity = "Chinese (single centre)",
    disease_state = "Infants receiving IV vancomycin for at least 3 days with normal renal function, defined as modified-Schwartz eGFR between 30 and 86 mL/min/1.73 m^2 (eGFR median 57.56, range 30-85.56); 29 of 61 born preterm; 47 of 61 co-administered meropenem or imipenem.",
    dose_range = "Intermittent IV infusion over >= 60 min, two to four times daily, TDM-adjusted; daily dose median 75 mg/day (range 24-320).",
    regions = "China (Children's Hospital of Nanjing Medical University, Nanjing)",
    renal_function = "Serum creatinine median 27.1 umol/L (range 14-315); eGFR 30-86 mL/min/1.73 m^2 by study design.",
    n_concentrations = 135L,
    notes = "Li 2021 Table 2, 'Model 1' column (the table header labels this group 'Reduced renal function'; the abstract, Methods and Results call it 'normal renal function'). Patients enrolled January 2017 - July 2021. One trough (30 min before the fifth dose) and one peak (30 min after the fifth dose) per enrolment: 69 troughs and 66 peaks. Whole-blood samples assayed by EMIT (calibration 2-50 ug/mL). NONMEM 7.4, ADVAN1 TRANS2. Screened but not retained on CL: sex, age, height, ALT, AST, BUN, cystatin C, albumin, total protein, eGFR, and meropenem/imipenem co-administration (Table 4)."
  )

  ini({
    # Structural parameters (Li 2021 Table 5 'NONMEM estimate'; equation
    # CL = theta1 * (WT/2.25)^theta3 * exp(theta5 * SCR/27.1),
    # V = theta2 * (WT/2.25)^theta4).
    lcl <- log(0.407); label("CL coefficient theta1 at WT = 2.25 kg before the creatinine factor (L/h)") # Li 2021 Table 5: theta1 = 0.407 (RSE 9.3%)
    lvc <- log(1.86); label("Typical V at WT = 2.25 kg (L)") # Li 2021 Table 5: theta2 = 1.86 (RSE 8.4%)

    e_wt_cl <- 1.24; label("Power exponent of WT/2.25 on CL (unitless)") # Li 2021 Table 5: theta3 = 1.24 (RSE 7.6%)
    e_wt_vc <- 1.28; label("Power exponent of WT/2.25 on V (unitless)") # Li 2021 Table 5: theta4 = 1.28 (RSE 13.3%)
    e_creat_cl <- -0.533; label("Exponential coefficient of CREAT/27.1 on CL (unitless)") # Li 2021 Table 5: theta5 = -0.533 (RSE 15.0%)

    # IIV on CL. Table 5 'BSV_CL' = 0.315 is read as the SD of eta (omega), so
    # the variance is 0.315^2 = 0.099225. The variance reading is ruled out
    # by the PROP_RV row of the same table; see the vignette.
    etalcl ~ 0.099225 # Li 2021 Table 5: BSV_CL = 0.315 (RSE 36.5%), squared

    # Proportional residual error (Methods 'Base Model': Y = F x (1 + eps)).
    propSd <- 0.319; label("Proportional residual error (fraction)") # Li 2021 Table 5: PROP_RV = 0.319 (RSE 17.9%)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 2.25)^e_wt_cl * exp(e_creat_cl * CREAT / 27.1)
    vc <- exp(lvc) * (WT / 2.25)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
