Fu_2021_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in adolescent and adult Chinese patients with haematological diseases and neutropenia (Fu 2021). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 70 kg) and clearance additionally by a power of CKD-EPI creatinine clearance (reference 116 mL/min/1.73 m^2); both volumes scale linearly with body weight. Between-subject variability on CL and V1 only; the residual-error magnitude is not reported."
  reference <- "Fu X, Lin L, Huang L, Guo L. Clinical application of vancomycin population pharmacokinetics model in patients with hematological diseases and neutropenia. Biopharm Drug Dispos. 2021;42(9):427-434. doi:10.1002/bdd.2303"
  vignette <- "Fu_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column BW. Reference 70 kg, printed as the literal denominator of Fu 2021 equations (6)-(9); this is a rounded standard, not the cohort mean of 56.84 kg (SD 11.09, range 33.00-83.50; Table 1, modelling patients). Exponent 0.75 on CL and Q and 1 on V1 and V2, carried in the equations without any estimate, RSE or confidence interval in Table 3, so both are treated as theory-based fixed exponents.",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column CLCR, computed with the CKD-EPI equation (Fu 2021 Methods section 2.2). Reference 116 mL/min/1.73 m^2, printed in equation (6), close to the cohort mean of 118.78 (SD 22.69, range 45.60-163.80; Table 1). Enters clearance as the power term (CRCL / 116)^0.895 with the exponent estimated (Table 3: 0.895, RSE 14.3%). CRCL was the only covariate retained of the age, sex, body weight, serum creatinine, CRCL, WBC, ANC, haemoglobin, platelet, total protein, albumin, ALT and AST screen. The cohort skews towards augmented renal clearance, so values below about 45 mL/min/1.73 m^2 are extrapolation.",
      source_name = "CLCR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    age_range = "17-83 years",
    age_median = "43.28 years (mean, SD 15.88)",
    weight_range = "33.00-83.50 kg",
    weight_median = "56.84 kg (mean, SD 11.09)",
    sex_female_pct = 45.5,
    race_ethnicity = "Chinese (single-centre study in Haikou, Hainan, China; race not otherwise reported)",
    disease_state = "Patients aged >= 14 years with haematological diseases and neutropenia (ANC < 0.5 x 10^9/L) receiving intermittent intravenous vancomycin: acute myeloid leukaemia 34, acute lymphoblastic leukaemia 16, lymphoma 10, aplastic anaemia 5, myelodysplastic syndrome 5, chronic myeloid leukaemia 4, multiple myeloma 2, chronic monocytic leukaemia 1. Patients on any blood purification treatment were excluded.",
    dose_range = "Intermittent IV vancomycin; 80.74% of regimens 1 g q12h, 7.34% 1 g q8h, 8.26% 0.5 g q6h, 1.83% 0.5 g q8h and 1.83% 0.5 g q24h. Daily dose 2.05 +/- 0.32 g/day (range 0.50-3.00); infusion rate 894.83 +/- 283.72 mg/h (range 250-1000) (Fu 2021 Table 2).",
    regions = "China (Department of Hematology, Hainan General Hospital, Haikou; 1 January 2018 - 1 January 2020)",
    renal_function = "CKD-EPI creatinine clearance 118.78 +/- 22.69 mL/min/1.73 m^2 (range 45.60-163.80); serum creatinine 53.66 +/- 18.36 umol/L (range 23-125)",
    n_concentrations = 152L,
    notes = "42 male and 35 female patients; 109 trough and 43 peak serum vancomycin concentrations (Siemens Viva-E immunoassay, quantitative range 2.0-50.0 ug/mL). Troughs drawn 0.5 h before an infusion, peaks 0.5-1 h after its end. Fit in NONMEM 7.3 with PsN 3.4.2; 1000-replicate bootstrap (90% convergence, Table 4) and external validation on 26 further patients (MDPE -4.68%, MAPE 18.74%, F20 52.27%, F30 68.18%). The paper's clinical-application arm (74 further patients with CRCL >= 90 mL/min/1.73 m^2, Table 5) used the model to choose 1 g q8h starting doses. The model itself was first published in Chinese by Lin et al. (Chinese Journal of Clinical Pharmacology and Therapeutics, 2021); all values here are those printed in Fu 2021."
  )

  ini({
    # Structural parameters, typical values at WT = 70 kg and
    # CRCL = 116 mL/min/1.73 m^2 (Fu 2021 Table 3; equations (6)-(9)).
    lcl <- log(6.84); label("Clearance at WT = 70 kg and CRCL = 116 mL/min/1.73 m^2 (L/h)")  # Fu 2021 Table 3: CL = 6.84 L/h (RSE 4.5%); equation (6)
    lvc <- log(20.5); label("Central volume of distribution at WT = 70 kg (L)")               # Fu 2021 Table 3: V1 = 20.5 L (RSE 17.3%); equation (7)
    lq <- log(15.2); label("Intercompartmental clearance at WT = 70 kg (L/h)")               # Fu 2021 Table 3: Q = 15.2 L/h (RSE 23.4%); equation (8)
    lvp <- log(50); label("Peripheral volume of distribution at WT = 70 kg (L)")             # Fu 2021 Table 3: V2 = 50 L (RSE 27.2%); equation (9)

    # Body-weight exponents: printed in equations (6)-(9) with no estimate,
    # RSE or CI in Table 3, i.e. theory-based allometric constants.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL (unitless)")  # Fu 2021 equation (6): (BW/70)^0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent of (WT/70) on Q (unitless)")    # Fu 2021 equation (8): (BW/70)^0.75
    e_wt_vc <- fixed(1); label("Exponent of (WT/70) on V1 (unitless)")                # Fu 2021 equation (7): (BW/70), power 1
    e_wt_vp <- fixed(1); label("Exponent of (WT/70) on V2 (unitless)")                # Fu 2021 equation (9): (BW/70), power 1

    # Renal-function effect on clearance (estimated).
    e_crcl_cl <- 0.895; label("Power exponent of (CRCL/116) on CL (unitless)")  # Fu 2021 Table 3: theta CLCR_CL = 0.895 (RSE 14.3%); equation (6)

    # Between-subject variability, exponential model (Fu 2021 Methods 2.2).
    # Table 3 reports eta1 and eta2 as percentages; they are read as
    # omega x 100 (variance = (P/100)^2), the convention of the base-model
    # report in Results 3.2.1, which prints omega^2 directly.
    etalcl ~ 0.031684 # Fu 2021 Table 3: eta1 = 17.8% (RSE 15%); 0.178^2
    etalvc ~ 0.1089 # Fu 2021 Table 3: eta2 = 33% (RSE 35.8%); 0.33^2
    # Equations (8) and (9) carry eta3 on Q and eta4 on V2, but neither the
    # final model (Table 3) nor the base model (Results 3.2.1) reports a
    # magnitude for them; held at zero rather than invented.
    etalq ~ fixed(0) # Fu 2021 equation (8) eta3; magnitude not reported
    etalvp ~ fixed(0) # Fu 2021 equation (9) eta4; magnitude not reported

    # Residual error: additive, proportional and combined ('mixed') forms
    # were compared (Methods 2.2) and the combined form carried into the
    # base model (Results 3.2.1), but no magnitude is reported anywhere in
    # the paper; both components are held at zero (typical/individual
    # predictions only).
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported in the source)") # Fu 2021 Methods 2.2 / Results 3.2.1: combined error, magnitude not reported
    addSd <- fixed(0); label("Additive residual error (ug/mL; not reported in the source)")         # Fu 2021 Methods 2.2 / Results 3.2.1: combined error, magnitude not reported
  })
  model({
    # Individual PK parameters (Fu 2021 equations (6)-(9)).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 116)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc is in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
