Zhao_2021_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in Chinese adult inpatients with renal function ranging from impaired to augmented renal clearance (Zhao 2021). Clearance is a sigmoid (Hill) function of Cockcroft-Gault creatinine clearance that saturates at a maximum of 5.58 L/h (CL = 5.58 x CrCl^1.5 / (93.8^1.5 + CrCl^1.5)); the central volume is 8.02 L in non-ICU and 35.7 L in ICU patients. Exponential inter-individual variability on CL and Vc; proportional residual error."
  reference <- "Zhao S, He N, Zhang Y, Wang C, Zhai S, Zhang C. Population Pharmacokinetic Modeling and Dose Optimization of Vancomycin in Chinese Patients with Augmented Renal Clearance. Antibiotics (Basel). 2021;10(10):1238. doi:10.3390/antibiotics10101238"
  vignette <- "Zhao_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Zhao 2021 Equation 1 symbol CG ('the creatinine clearance estimated by Cockcroft-Gault equation'); Table 1 footnote: 'CLcr, creatinine clearance calculated by the Cockcroft-Gault equation'. Stored under the canonical CRCL column in raw mL/min (Cockcroft-Gault is not a BSA-normalized estimator), per the CRCL register entry's provision for raw mL/min, matching the Dorajoo_2019_vancomycin.R and Goti_2018_vancomycin.R precedents. The paper does not state which body weight entered the Cockcroft-Gault calculation; only total body weight (TBW) was recorded (Methods 4.1). Table 1 summarises 'the mean value of CLcr for each patient during admission' (median 86.7, range 18.4-390.7 mL/min); the paper does not say whether the NONMEM dataset carried CrCl as a time-varying or a per-patient value, and the model accepts either. Enters CL through a Hill function with CrCl50 = 93.8 mL/min and Hill coefficient 1.5 (no centering). Applicability: patients with CrCl < 15 mL/min or on renal replacement therapy were excluded.",
      source_name = "CG"
    ),
    DIS_CRITILL = list(
      description = "Admission to the intensive care unit (1 = ICU inpatient, 0 = inpatient on another ward)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Zhao 2021 Equation 2: 'Vc = 8.02 (non-ICU patients) or 35.7 (ICU patients)'. The paper's covariate is 'the admission at the ICU' (Section 2.2). 82 of 209 patients (39.2%) were ICU inpatients (Table 1). Patients of the Surgical Department of the ICU were excluded (Methods 4.1). Time-fixed per patient.",
      source_name = "ICU"
    )
  )

  # Covariates screened during forward inclusion / backward elimination
  # (Methods 4.3) that were NOT retained in the final model. The paper
  # publishes no point estimates for them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 4.3); not retained. Cohort mean 66.0 years, SD 16.4 (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Source paper reports 'gender'; screened (Methods 4.3), not retained. 126 of 209 (60.3%) male, so 39.7% female (Table 1). Sex enters the model only indirectly through the Cockcroft-Gault CrCl calculation."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Source symbol TBW; screened (Methods 4.3), not retained as a direct covariate. Cohort mean 63.4 kg, SD 12.9 (Table 1). Weight enters the model only indirectly through the Cockcroft-Gault CrCl calculation."
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 209L,
    n_studies = 1L,
    age_range = "adults >= 18 years; mean 66.0, SD 16.4",
    weight_range = "mean 63.4 kg, SD 12.9 (total body weight)",
    sex_female_pct = 39.7,
    race_ethnicity = "Chinese (single-centre cohort at Peking University Third Hospital, Beijing)",
    disease_state = "Adult inpatients receiving intermittent intravenous vancomycin, with renal function from impaired (CrCl >= 15 mL/min) to augmented (CrCl >= 130 mL/min in 24.4%). 39.2% were ICU inpatients; 18.7% had shock and 3.3% multiple organ failure. End-stage renal disease, renal replacement therapy, acute kidney injury, and Hematology or Surgical-ICU inpatients were excluded.",
    dose_range = "Intermittent IV vancomycin at varied clinical doses; median daily dose 1875 mg (IQR 1461.9-2352.0).",
    regions = "China (Beijing)",
    renal_function = "Cockcroft-Gault CrCl (per-patient mean during admission) median 86.7 mL/min, range 18.4-390.7 (Table 1); 51 patients (24.4%) had CrCl >= 130 mL/min.",
    n_concentrations = 424L,
    notes = "Demographics from Zhao 2021 Table 1. Retrospective therapeutic-drug-monitoring data collected January 2010 - June 2018; 424 serum concentrations (CMIA immunoassay, Abbott ARCHITECT) from 209 patients, 49.3% of whom contributed a single sample; 69.3% of samples were drawn 5-12 h after the start of infusion. NONMEM 7.4.4 FOCE-I; two-compartment preferred over one-compartment on AIC (2089.5 vs 2162.9). Bootstrap success rate 84.6%."
  )

  ini({
    # Clearance (Zhao 2021 Equation 1 and Table 2 'Final Model Estimate'):
    #   CL = 5.58 * CG^1.5 / (93.8^1.5 + CG^1.5) * exp(eta1)   (L/h)
    # Check against the paper's own prose: at CrCl = 130 mL/min (the ARC
    # threshold) CL = 5.58 * 0.620 = 3.46 L/h, the lower end of the Abstract's
    # '3.46 and 5.58 L/h in patients with ARC'.
    lclmax <- log(5.58)
    label("Maximum (asymptotic) clearance as CrCl -> infinity (L/h)") # Table 2, 'CL max (L/h)' = 5.58 (RSE 17%)
    lcrcl50 <- log(93.8)
    label("Creatinine clearance giving half-maximal clearance (mL/min)") # Table 2, 'CG CLmax50' = 93.8 (RSE 24%); unit printed as L/h in Table 2 but mL/min in the Section 2.2 prose
    lhill <- log(1.5)
    label("Hill coefficient of the clearance / creatinine clearance relationship (unitless)") # Table 2, 's' = 1.5 (RSE 14%)

    # Central volume (Equation 2): 8.02 L non-ICU, 35.7 L ICU. The ICU value
    # is carried as a log-ratio effect on the non-ICU reference so that the
    # two printed typical values are recovered exactly.
    lvc <- log(8.02)
    label("Central volume of distribution in non-ICU patients (L)") # Table 2, 'V c non-ICU (L)' = 8.02 (RSE 12%)
    e_dis_critill_vc <- log(35.7 / 8.02)
    label("Log ratio of ICU to non-ICU central volume (unitless)") # Table 2, 'V c ICU (L)' = 35.7 (RSE 13%) vs 'V c non-ICU (L)' = 8.02

    lq <- log(2.66)
    label("Intercompartmental clearance (L/h)") # Table 2, 'Q(L/h)' = 2.66 (RSE 12%); Equation 3
    lvp <- log(36.8)
    label("Peripheral volume of distribution (L)") # Table 2, 'V p (L)' = 36.8 (RSE 15%); Equation 4

    # Inter-individual variability (Equation 5, exponential). Table 2 prints
    # NONMEM omega^2 variances: estimate +/- 1.96 * RSE * estimate reproduces
    # the printed bootstrap 95% CI (CL: 0.053-0.101 vs 0.05-0.10). The eta on
    # CL multiplies the whole Hill expression, which is identical to an eta
    # on the maximum clearance.
    etalclmax ~ 0.0771 # Table 2, 'IIV CL' = 0.0771 (RSE 16%; bootstrap 0.075, 95% CI 0.05-0.10); 28.3% CV
    etalvc ~ 0.223 # Table 2, 'IIV V c' = 0.223 (RSE 56%; bootstrap 0.20, 95% CI 0.0039-0.53); 51.1% CV

    # Residual error. Table 2 labels the row 'Additive residual error', but
    # Methods Equation 6 is Cobs = Cpred + Cpred * eps with eps described as
    # 'the proportional error', so the error model is proportional. The
    # printed 0.0466 is the NONMEM sigma^2 variance (0.0466 +/- 1.96 * 0.14 *
    # 0.0466 = 0.034-0.059 reproduces the bootstrap CI 0.032-0.060), so
    # propSd = sqrt(0.0466) = 0.2159.
    propSd <- 0.2159
    label("Proportional residual error (fraction)") # Table 2, 'Additive residual error' = 0.0466 (RSE 14%) as sigma^2; Equation 6 proportional form
  })
  model({
    # Individual parameters (Equations 1-4).
    crcl50 <- exp(lcrcl50)
    hill <- exp(lhill)
    cl <- exp(lclmax + etalclmax) * CRCL^hill / (crcl50^hill + CRCL^hill)
    vc <- exp(lvc + e_dis_critill_vc * DIS_CRITILL + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central / vc in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
