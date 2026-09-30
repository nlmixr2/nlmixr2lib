Li_2021_vancomycin_all <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Chinese infants (1-24 months) across all levels of renal function, eGFR >= 30 mL/min/1.73 m^2 (Li 2021, Model 3, pooled). CL scales as a power of body weight (reference 3.45 kg, estimated exponent 1.23) and exponentially with serum creatinine (exp(-0.377 * SCR / 19)); V scales as a power of body weight (reference 3.45 kg, estimated exponent 1.29). IIV on CL only; proportional residual error."
  reference <- "Li DY, Li L, Li GZ, Hu YH, Guo HL, Jing X, Chen F, Ji X, Xu J, Dai HR. Population Pharmacokinetics Modeling of Vancomycin Among Chinese Infants With Normal and Augmented Renal Function. Front Pediatr. 2021;9:713588. doi:10.3389/fped.2021.713588"
  vignette <- "Li_2021_vancomycin_renal_function"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 2 (Model 3 group): median 3.3 kg (range 1.15-14). Power-scaled on CL and V with reference 3.45 kg as printed in both the Table 7 header equation and the Results 'Model 3' equation (not the Table 2 median of 3.3 kg).",
      source_name = "WT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Table 2 (Model 3 group): median 18.95 umol/L (range 8-147). Enters CL as exp(e_creat_cl * CREAT / 19) (Table 7 header equation 'e**[theta5*SCR/19]'). This form is NOT centred: at CREAT = 19 umol/L the factor is exp(-0.377) = 0.686, not 1.",
      source_name = "SCR"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 115L,
    n_studies = 1L,
    age_range = "1-24 months postnatal age",
    age_median = "1 month (mean 2.87 months)",
    weight_range = "1.15-14 kg",
    weight_median = "3.3 kg",
    sex_female_pct = 36.5,
    race_ethnicity = "Chinese (single centre)",
    disease_state = "Infants receiving IV vancomycin for at least 3 days with modified-Schwartz eGFR >= 30 mL/min/1.73 m^2 (eGFR median 98.7, range 30-280), pooling the normal-renal-function and ARC groups; 46 of 115 born preterm; 91 of 115 co-administered meropenem or imipenem.",
    dose_range = "Intermittent IV infusion over >= 60 min, two to four times daily, TDM-adjusted; daily dose median 145 mg/day (range 24-640).",
    regions = "China (Children's Hospital of Nanjing Medical University, Nanjing)",
    renal_function = "Serum creatinine median 18.95 umol/L (range 8-147); eGFR median 98.7 mL/min/1.73 m^2 (range 30-280).",
    n_concentrations = 276L,
    notes = "Li 2021 Table 2, 'Model 3' column. Model 3 pools the Model 1 (n = 61) and Model 2 (n = 64) data; the group sizes sum to more than 115 because patients whose renal function changed during treatment contributed to both groups (Results 'Patients'). Patients enrolled January 2017 - July 2021. 158 troughs (30 min before the fifth dose) and 118 peaks (30 min after the fifth dose). Whole-blood samples assayed by EMIT (calibration 2-50 ug/mL). NONMEM 7.4, ADVAN1 TRANS2."
  )

  ini({
    # Structural parameters (Li 2021 Table 7 'NONMEM estimate'; equation
    # CL = theta1 * (WT/3.45)^theta3 * exp(theta5 * SCR/19),
    # V = theta2 * (WT/3.45)^theta4).
    lcl <- log(0.707); label("CL coefficient theta1 at WT = 3.45 kg before the creatinine factor (L/h)") # Li 2021 Table 7: theta1 = 0.707 (RSE 7.1%)
    lvc <- log(3.39); label("Typical V at WT = 3.45 kg (L)") # Li 2021 Table 7: theta2 = 3.39 (RSE 7%)

    e_wt_cl <- 1.23; label("Power exponent of WT/3.45 on CL (unitless)") # Li 2021 Table 7: theta3 = 1.23 (RSE 5.3%)
    e_wt_vc <- 1.29; label("Power exponent of WT/3.45 on V (unitless)") # Li 2021 Table 7: theta4 = 1.29 (RSE 8.1%)
    e_creat_cl <- -0.377; label("Exponential coefficient of CREAT/19 on CL (unitless)") # Li 2021 Table 7: theta5 = -0.377 (RSE 13%)

    # IIV on CL. Table 7 'BSV_CL' = 0.311 is read as the SD of eta (omega), so
    # the variance is 0.311^2 = 0.096721. See the vignette for why the
    # variance reading is ruled out.
    etalcl ~ 0.096721 # Li 2021 Table 7: BSV_CL = 0.311 (RSE 22.9%), squared

    # Proportional residual error (Methods 'Base Model': Y = F x (1 + eps)).
    propSd <- 0.335; label("Proportional residual error (fraction)") # Li 2021 Table 7: PROP_RV = 0.335 (RSE 11.1%)
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 3.45)^e_wt_cl * exp(e_creat_cl * CREAT / 19)
    vc <- exp(lvc) * (WT / 3.45)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
