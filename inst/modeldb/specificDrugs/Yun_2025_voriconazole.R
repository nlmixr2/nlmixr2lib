Yun_2025_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with fixed first-order absorption and linear elimination for voriconazole in Chinese adults with liver cirrhosis (Yun 2025); clearance is reduced in Child-Pugh B and C relative to Child-Pugh A and scales as a power function of serum retinol-binding protein (RBP4) normalised to the cohort median of 23.9 mg/L."
  reference <- "Yun J, Bai H, Wang Z, Ma Y, Liu W. Serum retinol binding protein as a novel marker for clearance and dosage optimization: pharmacokinetics study of voriconazole in a cirrhosis population. Front Pharmacol. 2025;16:1543323. doi:10.3389/fphar.2025.1543323"
  vignette <- "Yun_2025_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    RBP4 = list(
      description = "Serum retinol-binding protein (RBP4) concentration",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on clearance, (RBP4 / 23.9)^0.428, where 23.9 mg/L is the cohort",
        "median (Yun 2025 Supplementary Table S1; Discussion 'the median concentration",
        "was 23.9 mg/L'). Cohort mean (SD) 27.1 (24.9) mg/L, range 2.50-137 mg/L;",
        "7 of 78 patients (9.0%) had no RBP value and the paper does not state how",
        "they were imputed. RBP is synthesised by hepatocytes and falls with",
        "cirrhosis severity (median 31.8 / 26.4 / 6.4 mg/L in Child-Pugh A / B / C,",
        "Supplementary Table S6), so it acts as a quantitative liver-function marker",
        "alongside the categorical Child-Pugh grade."
      ),
      source_name = "RBP"
    ),
    HEPIMP_MOD = list(
      description = "Child-Pugh class B (moderate hepatic impairment) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Child-Pugh A or C; Child-Pugh A is the model reference when HEPIMP_MOD = HEPIMP_SEV = 0)",
      notes = paste(
        "Child-Pugh classification scheme (Yun 2025 Methods 2.3: score from hepatic",
        "encephalopathy, ascites, serum bilirubin, serum albumin and prothrombin",
        "time). Every subject had cirrhosis, so the reference category is Child-Pugh",
        "A, not normal hepatic function. Multiplicative factor 0.673 on CL",
        "(Table 2 'Child-Pugh B on CL'). 20 of 78 patients (25.6%) were Child-Pugh B."
      ),
      source_name = "Child-Pugh classification = B"
    ),
    HEPIMP_SEV = list(
      description = "Child-Pugh class C (severe hepatic impairment) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Child-Pugh A or B; Child-Pugh A is the model reference when HEPIMP_MOD = HEPIMP_SEV = 0)",
      notes = paste(
        "Child-Pugh classification scheme. Multiplicative factor 0.581 on CL",
        "(Table 2 'Child-Pugh C on CL'). 28 of 78 patients (35.9%) were Child-Pugh C.",
        "HEPIMP_MOD and HEPIMP_SEV are mutually exclusive."
      ),
      source_name = "Child-Pugh classification = C"
    )
  )

  # Covariates screened in the forward-inclusion step (Yun 2025 Supplementary
  # Table S4) but not retained in the final model.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Table S4 model13 (CLWT), dOFV -1.657, not significant. Cohort median 63.0 kg (41.0-120), Table S1."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table S4 model11 (CLAGE), dOFV -4.073, not significant. Cohort median 55.5 years (27.0-91.0), Table S1."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Table S4 model12 (CLHT), dOFV -1.41, not significant."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Table S4 model08 (CLBMI), dOFV -1.349, not significant."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Table S4 model09 (CLBSA), dOFV -1.488, not significant."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Table S4 model14 (CLCRP), dOFV -0.001, not significant."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9 cells/L",
      type = "continuous",
      notes = "Table S4 model15 (CLWBC), dOFV -0.847, not significant."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Table S4 model16 (CLHb), dOFV -3.271, not significant."
    ),
    PROCALCITONIN = list(
      description = "Procalcitonin",
      units = "ug/L",
      type = "continuous",
      notes = "Table S4 model17 (CLPCT), dOFV -1.279, not significant."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Table S4 model18 (CLALB), dOFV -0.106, not significant."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Table S4 model19 (CLALT), dOFV -0.206, not significant."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Table S4 model20 (CLAST), dOFV -6.008, not significant."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Table S4 model21 (CLTBIL), significant in the forward step (dOFV -8.935) but not retained once RBP was in the model (model33, dOFV -0.481)."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Table S4 model22 (CLDBIL), significant in the forward step (dOFV -11.246) but not retained once RBP was in the model (model34, dOFV -1.187)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Table S4 model24 (CLALP), dOFV -0.289, not significant."
    ),
    TBA = list(
      description = "Total serum bile acids",
      units = "umol/L",
      type = "continuous",
      notes = "Table S4 model25 (CLTBA), significant in the forward step (dOFV -28.24) but not retained once RBP was in the model (model35, dOFV -5.695)."
    ),
    GGT = list(
      description = "Gamma-glutamyl transpeptidase",
      units = "U/L",
      type = "continuous",
      notes = "Table S4 model23 (CLGGT), dOFV -2.556, not significant."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Table S4 model26 (CLCREAT), dOFV -0.116, not significant."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Table S4 model10 (CLCRCL), dOFV -1.164, not significant; eGFR (model27) was likewise not significant."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Table S4 model29 (CLPPI), dOFV +0.305, not significant. 48 of 78 patients (61.5%) received a PPI (Table S1)."
    ),
    TX_LIVER = list(
      description = "Liver transplant recipient",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'After liver transplantation' (Methods 2.3) and not retained. 11 of 78 patients (14.1%) were post-transplant (Table S1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    n_observations = 151L,
    age_range = "27-91 years",
    age_median = "55.5 years",
    weight_range = "41-120 kg",
    weight_median = "63.0 kg",
    sex_female_pct = 26.9,
    race_ethnicity = c(Chinese = 100),
    disease_state = paste(
      "Adults with liver cirrhosis treated with voriconazole for fungal infection:",
      "Child-Pugh A 30 (38.5%), B 20 (25.6%), C 28 (35.9%); 11 (14.1%) after liver",
      "transplantation; 8 (10.3%) on concomitant antituberculosis drugs."
    ),
    dose_range = paste(
      "Voriconazole 100-400 mg q12h in the first 24 h and 100-200 mg q12h",
      "maintenance (Yun 2025 Table 1); route not tabulated."
    ),
    regions = "Single centre: Beijing You'An Hospital, Capital Medical University, Beijing, China.",
    notes = paste(
      "Retrospective therapeutic-drug-monitoring study, September 2022 to August 2024.",
      "151 steady-state trough concentrations (0.1-15 mg/L UPLC-MS/MS assay; BLQ",
      "excluded) from 78 patients. External validation in 22 further patients",
      "(48 troughs). Model built in NONMEM 7.5.0 with PsN 5.2.6; covariates entered",
      "at p < 0.01 and retained at p < 0.001."
    )
  )

  ini({
    # Reference subject: Child-Pugh A (HEPIMP_MOD = HEPIMP_SEV = 0) with serum
    # RBP4 = 23.9 mg/L, the cohort median.

    # Trough-only data; ka fixed from Shi et al. 2019 (Yun 2025 Section 3.3,
    # and Supplementary Table S3 shows CL and V move < 0.6% across 0.5-2x ka).
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h)") # Yun 2025 Table 2: Ka (Fixed, 1/h) = 1.1
    lcl <- log(6.96); label("Clearance for Child-Pugh A at RBP4 = 23.9 mg/L (L/h)") # Yun 2025 Table 2: CL = 6.96 L/h (95% CI 5.7-8.5); Section 3.3 equation
    lvc <- log(745); label("Central volume of distribution (L)") # Yun 2025 Table 2: Vc = 745 L (95% CI 532-1040); Section 3.3 equation

    # Covariate effects on CL (Section 3.3 equation: CL = 6.96 * exp(eta) *
    # {1, 0.673, 0.581} * (RBP/23.9)^0.428). Child-Pugh factors are encoded as
    # log-scale shifts so the published fractions stay visible inside log().
    e_hepimp_mod_cl <- log(0.673); label("Log-scale CL shift for Child-Pugh B vs A (unitless)") # Yun 2025 Table 2: Child-Pugh B on CL = 0.673 (95% CI 0.514-0.881)
    e_hepimp_sev_cl <- log(0.581); label("Log-scale CL shift for Child-Pugh C vs A (unitless)") # Yun 2025 Table 2: Child-Pugh C on CL = 0.581 (95% CI 0.457-0.739)
    e_rbp4_cl <- 0.428; label("Power exponent of RBP4 / 23.9 on CL (unitless)") # Yun 2025 Table 2: RBP on CL = 0.428 (95% CI 0.288-0.568)

    # IIV on CL only: Supplementary Table S4 model02 = model01 - IIVVC and the
    # final model28 builds on it, so the printed eta_VC is not in the final
    # model. Table 2 'IIV of CL' 45.9 (95% CI 31.1-58) is read as a CV and
    # back-transformed: log(1 + 0.459^2) = 0.1912. The CI midpoint on the
    # variance scale is 0.1912 under the CV reading versus 0.2166 against
    # 0.459^2 = 0.2107 under the SD reading.
    etalcl ~ 0.1912 # Yun 2025 Table 2: IIV of CL = 45.9% (95% CI 31.1-58); log(1 + 0.459^2)

    # Combined residual error (Methods 2.5 'additive, proportional, and mixed').
    propSd <- 0.361; label("Proportional residual error (fraction)") # Yun 2025 Table 2: RSV_CV = 36.1% (95% CI 23.4-45.3)
    addSd <- 0.634; label("Additive residual error (mg/L)") # Yun 2025 Table 2: RSV_SD = 0.634 mg/L (95% CI 0.00-0.994)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + e_hepimp_mod_cl * HEPIMP_MOD + e_hepimp_sev_cl * HEPIMP_SEV + etalcl) *
      (RBP4 / 23.9)^e_rbp4_cl
    vc <- exp(lvc)

    # Oral doses go to depot (F = 1; Table S4 model04 '+F' was not significant);
    # intravenous doses go to central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
