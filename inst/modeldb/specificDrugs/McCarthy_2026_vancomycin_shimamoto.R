McCarthy_2026_vancomycin_shimamoto <- function() {
  description <- paste(
    "Two-compartment IV population PK model for vancomycin in paediatric cardiovascular intensive",
    "care unit (CVICU) patients: the Shimamoto 2024 congenital-heart-disease structure re-estimated",
    "by McCarthy 2026 on 219 paediatric CVICU patients, with the allometric weight exponent on",
    "clearance estimated instead of fixed at 0.75. Clearance scales with (WT/70)^theta, with a",
    "power of age-banded Schwartz creatinine clearance (CRCL/120) computed inside the model from",
    "age, weight, height and serum creatinine, and with a sigmoidal postmenstrual-age maturation",
    "(Hill 1.73, PMA50 36.9 weeks, both held at the published values); both volumes scale linearly",
    "with WT/70 and Q with (WT/70)^0.75 (Q held at the published value). Best-performing model for",
    "paediatric CVICU patients in McCarthy 2026.",
    sep = " "
  )
  reference <- paste(
    "McCarthy MG, Keizer RJ, Hughes JH. Evaluation and Improvement of Specialized Vancomycin",
    "Pharmacokinetic Models for Pediatric Cardiovascular Intensive Care Unit and Pediatric Oncology",
    "Patients. Clin Pharmacokinet. 2026;65(4):527. doi:10.1007/s40262-026-01623-7.",
    "(Table 6 'Shimamoto 2024 (refit)'; Supplementary file 2 'Shimamoto 2023 (refit)' NONMEM",
    "control stream; Table S2 structure.) Structural model from Shimamoto Y, Fukushima K, Mizuno T,",
    "et al. Model-informed vancomycin dosing optimization to address delayed renal maturation in",
    "infants and young children with critical congenital heart disease. Clin Pharmacol Ther.",
    "2024;115:239-247. doi:10.1002/cpt.3095.",
    sep = " "
  )
  vignette <- "McCarthy_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

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
      notes = paste(
        "(WT/70)^THETA(10) on CL (estimated 0.560), WT/70 on V1 and V2 and (WT/70)^0.75 on Q in",
        "the control stream. Also selects the Schwartz coefficient for infants under 1 year",
        "(0.244 below 2.5 kg, 0.277 at or above). McCarthy 2026 Table 1, Study 2 CVICU: median",
        "9 kg (range 3-91).",
        sep = " "
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the Schwartz coefficient k of the in-model creatinine clearance (control stream",
        "'eGFR calculated by revised Schwartz'): 0.244 for AGE < 1 year and WT < 2.5 kg, 0.277",
        "for AGE < 1 year and WT >= 2.5 kg, 0.296 for 1 <= AGE < 12 years, 0.413 for AGE >= 12",
        "years. Table 1, Study 2 CVICU: median 1.1 years (range 0.1-17.8).",
        sep = " "
      ),
      source_name = "AGE"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only in the in-model creatinine clearance CRCL = k * HT / CR (mL/min/1.73 m^2 with",
        "HT in cm and CR in mg/dL). Table 1, Study 2 CVICU: median 73 cm (range 51-187).",
        sep = " "
      ),
      source_name = "HT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Denominator of the in-model creatinine clearance CRCL = k * HT / CR. Table 1, Study 2",
        "CVICU: median 0.31 mg/dL (range 0.06-2.44).",
        sep = " "
      ),
      source_name = "CR"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maturation FMAT = PMA^HILL / (PMA^HILL + T50^HILL) with PMA in weeks (Table 6 footnote",
        "'T50 postmenstrual age at which 50% of adult clearance is reached (weeks)'). Where",
        "gestational age was unknown McCarthy 2026 imputed PMA as chronological age plus 40",
        "weeks. Table 1, Study 2 CVICU: median 99 weeks (range 46-970).",
        sep = " "
      ),
      source_name = "PMA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 219L,
    n_studies = 1L,
    age_range = "0.1-17.8 years",
    age_median = "1.1 years",
    weight_range = "3-91 kg",
    weight_median = "9 kg",
    sex_female_pct = 46.1,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Paediatric patients (< 20 years) tagged as cardiovascular intensive care unit (CVICU)",
      "patients and treated with intermittent intravenous vancomycin. Patients on haemodialysis,",
      "CRRT, ECMO or a ventricular assist device, or also tagged NICU or PICU, were excluded.",
      sep = " "
    ),
    dose_range = "Routine-care intermittent IV vancomycin; median dose 125 mg (range 25-2000), median interval 7.9 h, median infusion 1 h",
    regions = "United States (11 sites)",
    renal_function = "Serum creatinine median 0.31 mg/dL (range 0.06-2.44)",
    pma_range = "46-970 weeks",
    n_concentrations = 1136L,
    notes = paste(
      "McCarthy 2026 Table 1, Study 2 (CVICU only; the matched non-CVICU controls were not used",
      "in the refit). 391 treatment courses, each treated as a separate individual in NONMEM",
      "(Methods, Pharmacokinetic Analysis). 101/219 female. Routine therapeutic drug monitoring",
      "records from the InsightRX Nova dosing software, January 2018 to April 2025.",
      sep = " "
    )
  )

  ini({
    lcl <- log(3.56332); label("Clearance at WT 70 kg and CRCL 120 mL/min/1.73 m^2, before maturation (L/h)")
    # Table 6 refit thetaCL 3.56 (SE 0.10); control stream $THETA(1) 3.56332
    lvc <- log(12.7238); label("Central volume at WT 70 kg (L)")
    # Table 6 refit thetaVD1 12.72 (SE 0.31); control stream $THETA(2) 12.7238
    lq <- fixed(log(2.95)); label("Intercompartmental clearance at WT 70 kg (L/h)")
    # Table 6 thetaQ 2.95 (refit blank = fixed to published); control stream $THETA(3) 2.95 FIX
    lvp <- log(26.7204); label("Peripheral volume at WT 70 kg (L)")
    # Table 6 refit thetaVD2 26.72 (SE 0.08); control stream $THETA(4) 26.7204

    e_wt_cl <- 0.559713; label("Allometric exponent on WT/70 for CL (unitless)")
    # Table 6 refit thetaWT 0.56 (SE 0.05); control stream $THETA(10) 0.559713
    e_wt_q <- fixed(0.75); label("Allometric exponent on WT/70 for Q (unitless)")
    # control stream 'Q = TVQ * (WT/70)**0.75'; Table S2
    e_crcl_cl <- 0.630524; label("Exponent on CRCL/120 for CL (unitless)")
    # Table 6 refit thetaCRCL 0.63 (SE 0.07); control stream $THETA(5) 0.630524
    hill_mat <- fixed(1.73); label("Hill coefficient of clearance maturation (unitless)")
    # Table 6 thetaHILL 1.73 (refit blank = fixed to published); control stream $THETA(6) 1.73 FIX
    tmat50 <- fixed(36.9); label("Postmenstrual age at 50% clearance maturation (weeks)")
    # Table 6 thetaT50 36.90 (refit blank = fixed to published); control stream $THETA(7) 36.9 FIX

    # IIV: $OMEGA variances; Table 6 prints their square roots (omega). No IIV
    # on V1 or Q (V1 = TVV1, Q = TVQ * ... in the control stream).
    # sqrt(0.068973) = 0.263 -> Table 6 refit omegaCL 0.26
    etalcl ~ 0.068973
    # sqrt(0.217882) = 0.467 -> Table 6 refit omegaVD2 0.47
    etalvp ~ 0.217882

    # Residual error: W = SQRT(IPRED**2 * PROP**2 + ADD**2), Y = IPRED + W*ERR(1),
    # $SIGMA 1 FIX, so PROP and ADD are standard deviations.
    propSd <- 0.207734; label("Proportional residual error (fraction)")
    # Table 6 refit thetaPROP 0.21 (SE 0.06); control stream $THETA(8) 0.207734
    addSd <- 0.890438; label("Additive residual error (mg/L)")
    # Table 6 refit thetaADD 0.89 (SE 0.23); control stream $THETA(9) 0.890438
  })
  model({
    # Age-banded Schwartz coefficient, as in the control stream
    # ('eGFR calculated by revised Schwartz'); the later IF wins at AGE = 12.
    k_schwartz <- 0.244
    if (AGE < 1 & WT >= 2.5) k_schwartz <- 0.277
    if (AGE >= 1 & AGE <= 12) k_schwartz <- 0.296
    if (AGE >= 12) k_schwartz <- 0.413
    crcl <- k_schwartz * HT / CREAT

    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (crcl / 120)^e_crcl_cl * fmat
    vc <- exp(lvc) * (WT / 70)
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The control stream's third state, DADT(3) = A(1)/V1, is an AUC
    # accumulator with no effect on the PK; it is omitted here.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
