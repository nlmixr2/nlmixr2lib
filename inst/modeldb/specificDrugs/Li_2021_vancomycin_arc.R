Li_2021_vancomycin_arc <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Chinese infants (1-24 months) with augmented renal clearance, eGFR >= 86 mL/min/1.73 m^2 (Li 2021, Model 2). Body weight is the only covariate: CL and V scale as powers of body weight (reference 4.6 kg, estimated exponents 1.03 and 0.918). IIV on CL only; proportional residual error."
  reference <- "Li DY, Li L, Li GZ, Hu YH, Guo HL, Jing X, Chen F, Ji X, Xu J, Dai HR. Population Pharmacokinetics Modeling of Vancomycin Among Chinese Infants With Normal and Augmented Renal Function. Front Pediatr. 2021;9:713588. doi:10.3389/fped.2021.713588"
  vignette <- "Li_2021_vancomycin_renal_function"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 2 (Model 2 group): median 4.60 kg (range 2.2-14). Power-scaled on CL and V with reference 4.6 kg (Table 6 header equation and Results 'Model 2' equation).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 1L,
    age_range = "1-24 months postnatal age",
    age_median = "1 month (mean 3.36 months)",
    weight_range = "2.2-14 kg",
    weight_median = "4.60 kg",
    sex_female_pct = 35.9,
    race_ethnicity = "Chinese (single centre)",
    disease_state = "Infants receiving IV vancomycin for at least 3 days with augmented renal clearance (ARC), defined as modified-Schwartz eGFR >= 86 mL/min/1.73 m^2 (eGFR median 128, range 90.8-280); 17 of 64 born preterm; 54 of 64 co-administered meropenem or imipenem.",
    dose_range = "Intermittent IV infusion over >= 60 min, two to four times daily, TDM-adjusted; daily dose median 216 mg/day (range 50-640).",
    regions = "China (Children's Hospital of Nanjing Medical University, Nanjing)",
    renal_function = "Serum creatinine median 16.8 umol/L (range 8-29); eGFR >= 86 mL/min/1.73 m^2 by study design.",
    n_concentrations = 139L,
    notes = "Li 2021 Table 2, 'Model 2' column (the table header labels this group 'Normal renal function'; the abstract, Methods and Results call it the ARC group). Patients enrolled January 2017 - July 2021. One trough (30 min before the fifth dose) and one peak (30 min after the fifth dose) per enrolment: 88 troughs and 51 peaks. Whole-blood samples assayed by EMIT (calibration 2-50 ug/mL). NONMEM 7.4, ADVAN1 TRANS2. Serum creatinine and the other screened covariates were not retained (Results 'Model Building')."
  )

  ini({
    # Structural parameters (Li 2021 Table 6 'NONMEM estimate'; equation
    # CL = theta1 * (WT/4.6)^theta3, V = theta2 * (WT/4.6)^theta4).
    lcl <- log(0.756); label("Typical CL at WT = 4.6 kg (L/h)") # Li 2021 Table 6: theta1 = 0.756 (RSE 5.2%)
    lvc <- log(4.89); label("Typical V at WT = 4.6 kg (L)") # Li 2021 Table 6: theta2 = 4.89 (RSE 8.8%)

    e_wt_cl <- 1.03; label("Power exponent of WT/4.6 on CL (unitless)") # Li 2021 Table 6: theta3 = 1.03 (RSE 13.1%)
    e_wt_vc <- 0.918; label("Power exponent of WT/4.6 on V (unitless)") # Li 2021 Table 6: theta4 = 0.918 (RSE 17.8%)

    # IIV on CL. Table 6 'BSV_CL' = 0.312 is read as the SD of eta (omega), so
    # the variance is 0.312^2 = 0.097344. See the vignette for why the
    # variance reading is ruled out.
    etalcl ~ 0.097344 # Li 2021 Table 6: BSV_CL = 0.312 (RSE 27.5%), squared

    # Proportional residual error (Methods 'Base Model': Y = F x (1 + eps)).
    propSd <- 0.319; label("Proportional residual error (fraction)") # Li 2021 Table 6: PROP_RV = 0.319 (RSE 16.5%)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 4.6)^e_wt_cl
    vc <- exp(lvc) * (WT / 4.6)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
