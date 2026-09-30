Zhou_2020_paclitaxel_liposomal <- function() {
  description <- paste(
    "Three-compartment population PK model for total plasma paclitaxel",
    "(liposome-encapsulated plus released drug) after a 3-h intravenous",
    "infusion of paclitaxel liposome (Lipusu) 175 mg/m^2 in adults with",
    "squamous non-small cell lung cancer (Zhou 2020, n = 45). Linear",
    "first-order elimination from the central compartment, a deep and a",
    "shallow peripheral compartment, and between-subject variability on",
    "clearance only. No covariate was retained: age, sex, body weight,",
    "total bilirubin, albumin, serum creatinine, creatinine clearance and",
    "co-administration of aidi injection were all screened and rejected.",
    "The companion exposure-safety model is",
    "Zhou_2020_paclitaxel_liposomal_neutropenia."
  )
  reference <- paste(
    "Zhou H, Yan J, Chen W, Yang J, Liu M, Zhang Y, Shen X, Ma Y, Hu X,",
    "Wang Y, Du K, Li G. Population Pharmacokinetics and Exposure-Safety",
    "Relationship of Paclitaxel Liposome in Patients With Non-small Cell",
    "Lung Cancer. Front Oncol. 2020;10:1731 (issue dated 5 February 2021).",
    "doi:10.3389/fonc.2020.01731.",
    "Final estimates from Table 2; model structure confirmed against the",
    "NONMEM control stream in Supplementary Data Sheet 1.",
    sep = " "
  )
  vignette <- "Zhou_2020_paclitaxel_liposomal"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age. Screened on the PK parameters (Methods 'Model Development and Evaluation'; Results 'Covariate Analysis') and not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 59 (range 36-75) years, Table 1. NONMEM $INPUT column AGE.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex. Screened and not retained; 91% of the cohort was male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "41 of 45 patients male, Table 1. NONMEM $INPUT column GEND (coding not printed).",
      source_name = "GEND"
    ),
    WT = list(
      description = "Body weight. Screened and not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 71.5 (range 45-100) kg, Table 1. NONMEM $INPUT column WT.",
      source_name = "WT"
    ),
    TBILI = list(
      description = "Total bilirubin. Screened and not retained.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 8.4 (range 2.3-35.5) umol/L, Table 1. NONMEM $INPUT column TB.",
      source_name = "TB"
    ),
    ALB = list(
      description = "Serum albumin. Screened and not retained.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 42.6 (range 30.7-50.5) g/L, Table 1. NONMEM $INPUT column ALB.",
      source_name = "ALB"
    ),
    CREAT = list(
      description = "Serum creatinine. Screened and not retained.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 69 (range 30.1-152) umol/L, Table 1. NONMEM $INPUT column SCR.",
      source_name = "SCR"
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (not BSA-normalised). Screened and not retained.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Median 101.1 (range 53.4-159.1) mL/min, Table 1. Estimating equation not stated. NONMEM $INPUT column CLCR.",
      source_name = "CLCR"
    ),
    CONMED_AIDI = list(
      description = "Co-administration of aidi injection (a traditional Chinese medicine compound preparation, 100 mL) before the cycle-2 paclitaxel liposome infusion; 1 = given, 0 = not given. Screened and not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no aidi injection (cycle 1)",
      notes = "Aidi was given in cycle 2 only (43 of 45 patients), so it is confounded with cycle. Figure 4 shows overlapping dose-normalised profiles with and without it. NONMEM $INPUT column ADAD.",
      source_name = "ADAD"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "paclitaxel (total: liposome-encapsulated plus released)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "paclitaxel (total: liposome-encapsulated plus released)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "paclitaxel (total: liposome-encapsulated plus released)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    n_observations = 349L,
    age_range = "36-75 years",
    age_median = "59 years",
    weight_range = "45-100 kg",
    weight_median = "71.5 kg",
    sex_female_pct = 8.9,
    race_ethnicity = "Chinese (single centre in Beijing); race not tabulated",
    disease_state = "Cytologically or histologically confirmed squamous non-small cell lung cancer, Karnofsky score > 70, normal hepatic and renal function",
    dose_range = "Paclitaxel liposome 175 mg/m^2 as a 3-h IV infusion on day 1 of a 3-week cycle, rounded to vial size (210, 240, 270 or 300 mg), followed on day 2 by cisplatin 75 mg/m^2 or carboplatin AUC 4-5; two cycles sampled",
    regions = "China",
    renal_function = "Creatinine clearance median 101.1 (53.4-159.1) mL/min",
    hepatic_function = "Total bilirubin median 8.4 (2.3-35.5) umol/L",
    co_medication = "Premedication with dexamethasone, diphenhydramine and cimetidine; aidi injection before chemotherapy in cycle 2",
    notes = "Baseline demographics from Table 1. Sparse sampling at 1.5 h (mid-infusion), 3 h (end of infusion), 4, 6 and 21 h after the start of the infusion in each of two cycles; 43 of 45 patients completed cycle 2. Total (encapsulated plus unencapsulated) paclitaxel measured by LC-MS/MS, LLOQ 10 ng/mL; no sample was below the LLOQ."
  )

  ini({
    # Table 2 of Zhou 2020 prints the exponentiated typical values. The
    # supplementary NONMEM control stream (Data Sheet 1) parameterises
    # every structural parameter as EXP(THETA(n) + ETA(n)) and confirms the
    # compartment mapping: its V2/CL2 pair is the DEEP peripheral
    # compartment (Table 2 Vp1/Q1) and its V3/CL3 pair is the SHALLOW one
    # (Table 2 Vp2/Q2). The $THETA block of that stream holds initial
    # estimates (exp(3.09) = 21.98 L/h vs the final 21.55 L/h), so the final
    # values below come from Table 2, not from the stream.
    lcl <- log(21.55); label("Clearance from the central compartment, CL (L/h)") # Table 2, CL = 21.55 (95% CI 18.59-25) L/h
    lvc <- log(0.9248); label("Central volume of distribution, Vc (L)") # Table 2, Vc = 0.9248 (95% CI 0.7733-1.106) L
    lq <- log(4.62); label("Intercompartmental clearance between central and deep peripheral compartments, Q1 (L/h)") # Table 2, Q1 = 4.62 (95% CI 3.915-5.453) L/h
    lvp <- log(44.15); label("Deep peripheral volume of distribution, Vp1 (L)") # Table 2, Vp1 = 44.15 (95% CI 36.5-53.39) L
    lq2 <- log(15.85); label("Intercompartmental clearance between central and shallow peripheral compartments, Q2 (L/h)") # Table 2, Q2 = 15.85 (95% CI 10.95-22.96) L/h
    lvp2 <- log(5.577); label("Shallow peripheral volume of distribution, Vp2 (L)") # Table 2, Vp2 = 5.577 (95% CI 4.838-6.429) L

    # Table 2 row 'Inter-individual variable (%) CL' = 20.65%, read as the
    # standard deviation of eta on the log scale: omega^2 = 0.2065^2. The
    # control stream's $OMEGA for CL1 (0.0427, sqrt = 20.66%) supports the
    # SD reading; the other five $OMEGA diagonals are 0 FIX.
    etalcl ~ 0.04264 # Table 2, IIV CL 20.65% -> 0.2065^2

    # Table 2 row 'Residual error(%) s1' = 44.55%; control stream
    # Y = F*(1+EPS(1)) is purely proportional.
    propSd <- 0.4455; label("Proportional residual error (fraction)") # Table 2, sigma1 = 44.55%
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg, volume in L -> total paclitaxel concentration in mg/L
    # (the units of Figures 3 and 4).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
