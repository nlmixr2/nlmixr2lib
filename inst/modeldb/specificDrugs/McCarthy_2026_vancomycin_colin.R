McCarthy_2026_vancomycin_colin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for vancomycin in children with and without cancer:",
    "the Colin 2019 pooled 'throughout life' structure with its cancer-status effect on clearance,",
    "re-estimated by McCarthy 2026 on 371 paediatric oncology patients and 371 propensity-score",
    "matched paediatric non-oncology patients. Volumes scale linearly with WT/70; clearance scales",
    "with the individual central volume to the 0.75 power and with sigmoidal postmenstrual-age",
    "maturation and decline functions and an exponential serum-creatinine term centred on an",
    "age-standardised creatinine; Q scales with the individual peripheral volume to the 0.75 power.",
    "CL, V1, V2, Q, the cancer-status effect, IIV and residual error were re-estimated; the",
    "maturation, decline and creatinine shape parameters were held at the Colin 2019 values.",
    sep = " "
  )
  reference <- paste(
    "McCarthy MG, Keizer RJ, Hughes JH. Evaluation and Improvement of Specialized Vancomycin",
    "Pharmacokinetic Models for Pediatric Cardiovascular Intensive Care Unit and Pediatric Oncology",
    "Patients. Clin Pharmacokinet. 2026;65(4):527. doi:10.1007/s40262-026-01623-7.",
    "(Table 4 'Colin 2019 (Onc) refit'; Supplementary file 2 'Colin 2019 (Onc) refit' NONMEM",
    "control stream; Table S2 structure.) Structural model from Colin PJ, Allegaert K, Thomson AH,",
    "et al. Vancomycin pharmacokinetics throughout life: results from a pooled population analysis",
    "and evaluation of current dosing recommendations. Clin Pharmacokinet. 2019;58:767-780.",
    "doi:10.1007/s40262-018-0727-5.",
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
        "Fsize = WT/70 in the control stream, applied linearly to V1 and V2. CL and Q are",
        "not scaled by WT directly but by the individual volume ratios (V1/V1std)^0.75 and",
        "(V2/V2std)^0.75. McCarthy 2026 Table 1, Study 1: oncology median 26 kg (range 5-112),",
        "matched non-oncology median 28 kg (range 5-124).",
        sep = " "
      ),
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied in weeks (data item PMA) and converted to years inside the model as",
        "PMA_YRS = PMA/52, the scale on which the control stream evaluates the maturation,",
        "decline and standardised-creatinine functions (comment in the stream: 'paper",
        "calculates Fmat with units of weeks but model files use units of years. Kept years').",
        "Reference PMA 35 + 40/52 years (a 35-year-old born at 40 weeks). Where gestational age",
        "was unknown McCarthy 2026 imputed PMA as chronological age plus 40 weeks (Methods, Data",
        "Collection). Table 1, Study 1: oncology median 413 weeks (range 51-1081), non-oncology",
        "median 454 weeks (range 56-998).",
        sep = " "
      ),
      source_name = "PMA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as Fscr = exp(-EXPO2 * (CR - SCRstd)), where SCRstd is the PMA-standardised",
        "creatinine of the control stream (0.84 mg/dL for a 35-year-old). Units fixed by",
        "Table 2 footnote 'CR serum creatinine (mg/dL)'. Table 1, Study 1: oncology median",
        "0.35 mg/dL (range 0.06-1.15), non-oncology median 0.34 mg/dL (range 0.06-1.08).",
        sep = " "
      ),
      source_name = "CR"
    ),
    DIS_CANCER_PED = list(
      description = "Paediatric oncology patient indicator (1 = oncology patient, 0 = non-oncology patient)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (paediatric non-oncology patient)",
      notes = paste(
        "CANC in the control stream, entering as TVCL * (1 + CANC * THETA(12)). 1 = patient tagged",
        "as an oncology patient in the dosing software (82% haematological malignancy, 3.2% solid",
        "tumour, 15% unspecified; Table 1); 0 = propensity-score matched patient tagged as neither",
        "oncology nor CVICU (Methods, Propensity Score Matching). The refit estimate is",
        "0.0563 (Table 4: 0.06, SE 0.40), against 0.294 in the original Colin 2019 adult",
        "haematological-malignancy data.",
        sep = " "
      ),
      source_name = "CANC"
    )
  )

  covariatesDataExcluded <- list(
    SAMPLE_CAPILLARY = list(
      description = "Heel-prick (capillary) sampling indicator of the original Colin 2019 model",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "HEEL in the control stream, entering as (1 - HEEL * 0.312) on V1 and (1 - HEEL * 0.597)",
        "on Q, both FIXED at the Colin 2019 values (THETA(10), THETA(11)). McCarthy 2026 set",
        "HEEL = 0 for every patient because heel-prick status was not collected (Table 2",
        "footnote a), so both factors are identically 1 in the refit and the terms are omitted",
        "here. Use the original Colin 2019 model for neonatal heel-prick data.",
        sep = " "
      ),
      source_name = "HEEL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 742L,
    n_studies = 1L,
    age_range = "0.2-20.0 years",
    age_median = "7.2 years (oncology), 7.9 years (non-oncology)",
    weight_range = "5-124 kg",
    weight_median = "26 kg (oncology), 28 kg (non-oncology)",
    sex_female_pct = 46.0,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Paediatric patients (< 20 years) treated with intermittent intravenous vancomycin: 371",
      "oncology patients (82% haematological malignancy, 3.2% solid tumour, 15% unspecified; 13%",
      "neutropenic) and 371 propensity-score matched non-oncology patients. Patients on",
      "haemodialysis, CRRT, ECMO or a ventricular assist device, or tagged NICU or PICU, were",
      "excluded.",
      sep = " "
    ),
    dose_range = "Routine-care intermittent IV vancomycin; median dose 428 mg (oncology) and 500 mg (non-oncology), median interval about 6.5 h, median infusion 1 h",
    regions = "United States (20 sites for the oncology patients, 58 for the matched controls)",
    renal_function = "Serum creatinine median 0.35 mg/dL (oncology) and 0.34 mg/dL (non-oncology)",
    n_concentrations = 2649L,
    notes = paste(
      "McCarthy 2026 Table 1, Study 1 (oncology + matched non-oncology). The refit used 1094",
      "treatment courses (592 + 502), each treated as a separate individual in NONMEM",
      "(Methods, Pharmacokinetic Analysis). Female percentage pooled from 166/371 and 175/371.",
      "Data are routine therapeutic drug monitoring records from the InsightRX Nova dosing",
      "software, January 2018 to April 2025.",
      sep = " "
    )
  )

  ini({
    # Reference subject: WT 70 kg, PMA 35 + 40/52 years, CREAT equal to the
    # age-standardised value, DIS_CANCER_PED = 0. The control stream
    # parameterises CL, V1, V2 and Q as EXP(THETA), so the $THETA values are
    # already on the log scale.
    lcl <- 1.51077; label("Clearance at the reference subject (L/h)")
    # Table 4 refit thetaCL 1.51 (SE 0.02); control stream $THETA(1) 1.51077 (exp = 4.53 L/h)
    lvc <- 3.57148; label("Central volume at 70 kg (L)")
    # Table 4 refit thetaVD1 3.57 (SE 0.03); control stream $THETA(2) 3.57148 (exp = 35.6 L)
    lvp <- 3.9105; label("Peripheral volume at 70 kg (L)")
    # Table 4 refit thetaVD2 3.91 (SE 0.06); control stream $THETA(3) 3.9105 (exp = 49.9 L)
    lq <- 0.807181; label("Intercompartmental clearance at the reference peripheral volume (L/h)")
    # Table 4 refit thetaQ2 0.81 (SE 0.13); control stream $THETA(4) 0.807181 (exp = 2.24 L/h)

    # Size: exponents written as constants in the control stream
    # (TVCL = ... * (V1/RV1)**0.75, TVQ = ... * (V2/RV2)**0.75).
    e_wt_cl <- fixed(0.75); label("Exponent on the individual V1/V1std ratio for CL (unitless)")
    # control stream TVCL '(V1/RV1) ** 0.75'; Table S2 '(V1/thetaV1)^0.75'
    e_wt_q <- fixed(0.75); label("Exponent on the individual V2/V2std ratio for Q (unitless)")
    # control stream TVQ '(V2/RV2)**0.75'; Table S2 '(V2/thetaV2)^0.75'

    # Maturation, decline and creatinine shape: held at the Colin 2019 values
    # (control stream $THETA(5)-(9) FIX), log scale as in the stream.
    ltdec50 <- fixed(4.12); label("Postmenstrual age at 50% clearance decline (log years)")
    # control stream $THETA(5) DEG50 4.12 FIX (exp = 61.6 years)
    ltmat50 <- fixed(-0.113); label("Postmenstrual age at 50% clearance maturation (log years)")
    # control stream $THETA(6) MAT50 -0.113 FIX (exp = 0.893 years)
    lhill_mat <- fixed(1.06); label("Hill coefficient of clearance maturation (log unitless)")
    # control stream $THETA(7) GAM 1.06 FIX (exp = 2.89)
    lhill_dec <- fixed(0.807); label("Hill coefficient of clearance decline (log unitless)")
    # control stream $THETA(8) GAM2 0.807 FIX, used as GAM2 = -EXP(THETA(8)) (exp = 2.24)
    le_creat_cl <- fixed(-0.433); label("Exponential serum-creatinine coefficient on CL (log dL/mg)")
    # control stream $THETA(9) EXPO2 -0.433 FIX (exp = 0.649 dL/mg)

    e_cancer_ped_cl <- 0.0563011; label("Fractional increase in CL for paediatric oncology patients (unitless)")
    # Table 4 refit thetaCANC 0.06 (SE 0.40); control stream $THETA(12) 0.0563011

    # IIV: $OMEGA variances; Table 4 prints their square roots (omega).
    # sqrt(0.0652814) = 0.256 -> Table 4 refit omegaCL 0.26
    etalcl ~ 0.0652814
    # sqrt(0.171886) = 0.415 -> Table 4 refit omegaVD1 0.41
    etalvc ~ 0.171886
    # sqrt(0.862485) = 0.929 -> Table 4 refit omegaVD2 0.93
    etalvp ~ 0.862485

    # Residual error: Y = IPRED*(1+EPS(1)) + EPS(2); $SIGMA are variances.
    propSd <- 0.222634; label("Proportional residual error (fraction)")
    # $SIGMA(1) 0.0495662, sqrt = 0.2226 -> Table 4 refit thetaPROP 0.22
    addSd <- 0.805160; label("Additive residual error (mg/L)")
    # $SIGMA(2) 0.648283, sqrt = 0.8052 -> Table 4 refit thetaADD 0.81
  })
  model({
    # Postmenstrual age in years and the reference age (control stream
    # PMA_YRS = PMA/52, PMA_R = 35+(40/52)).
    pma_yrs <- PAGE / 52
    pma_ref <- 35 + 40 / 52

    tdec50 <- exp(ltdec50)
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    hill_dec <- exp(lhill_dec)
    e_creat_cl <- exp(le_creat_cl)

    # Standardised serum creatinine (control stream LnSCRstd, which uses
    # LOG10 of PMA in years; 0.84 mg/dL at the reference age).
    creat_std <- exp(-1.22839 + log10(pma_yrs) * 0.67190 + 6.27017 * exp(-3.10940 * pma_yrs))

    fsize <- WT / 70
    fmat <- pma_yrs^hill_mat / (pma_yrs^hill_mat + tmat50^hill_mat)
    fmat_ref <- pma_ref^hill_mat / (pma_ref^hill_mat + tmat50^hill_mat)
    # Fage = PMA^GAM2 / (PMA^GAM2 + DEG50^GAM2) with GAM2 = -hill_dec
    fdec <- pma_yrs^(-hill_dec) / (pma_yrs^(-hill_dec) + tdec50^(-hill_dec))
    fdec_ref <- pma_ref^(-hill_dec) / (pma_ref^(-hill_dec) + tdec50^(-hill_dec))
    fscr <- exp(-e_creat_cl * (CREAT - creat_std))

    vc <- exp(lvc + etalvc) * fsize
    vp <- exp(lvp + etalvp) * fsize
    # CL and Q scale with the INDIVIDUAL volumes (V1 and V2 include their
    # etas), as in the control stream.
    cl <- exp(lcl + etalcl) * (vc / exp(lvc))^e_wt_cl * (fmat / fmat_ref) * (fdec / fdec_ref) *
      fscr * (1 + DIS_CANCER_PED * e_cancer_ped_cl)
    q <- exp(lq) * (vp / exp(lvp))^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # S1 = V1: dose in mg over volume in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
