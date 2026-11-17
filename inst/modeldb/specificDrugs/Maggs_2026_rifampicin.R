Maggs_2026_rifampicin <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral rifampicin (300 mg",
    "twice daily) and 25-desacetyl rifampicin given with oral fusidic acid",
    "(500 mg three times daily) to adults with staphylococcal bone, joint or",
    "endovascular infection (Maggs 2026). One-compartment disposition for",
    "each analyte with first-order absorption after an absorption lag. The",
    "absorbed dose splits at absorption: a fraction (1 - FP = 1.9 percent) is",
    "converted to the metabolite pre-systemically and the rest enters the",
    "rifampicin compartment; all systemic rifampicin clearance forms",
    "25-desacetyl rifampicin. Clearances and volumes are apparent and",
    "allometrically scaled to 70 kg (exponents fixed at 0.75 and 1).",
    "Co-administered fusidic acid lowers the metabolite's clearance by 39.2",
    "percent. Inter-individual variability on ka and metabolite clearance,",
    "inter-occasion variability (four occasions) on metabolite clearance,",
    "and proportional residual error on log-transformed concentrations of",
    "each analyte."
  )
  reference <- paste(
    "Maggs C, Manning L, Yoo O, Salman S, Page-Sharp M, Aboltins C. (2026).",
    "A prospective pharmacokinetic interaction study between rifampicin and",
    "fusidic acid for the treatment of staphylococcal infections.",
    "J Antimicrob Chemother 81(10):959-967.",
    "doi:10.1093/jac/dkag293."
  )
  vignette <- "Maggs_2026_rifampicin_fusidicAcid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    central_desacetylrif = list(
      analyte = "25-desacetyl rifampicin",
      units = "mg (rifampicin mass equivalents)",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight, used for a priori allometric scaling of clearances and volumes.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 70 kg (Table 3 units 'L h-1.70 kg-1' and footnote 'for an",
        "individual with body weight of 70 kg'); exponents 0.75 on clearance",
        "and 1 on volume fixed a priori (Methods; Table S2 run 1302 prints",
        "'CLM = TVCLM * (WT/70)**0.75 * ...'). The cohort's body weights are",
        "not reported (median BMI 31 kg/m2)."
      ),
      source_name = "WT"
    ),
    CONMED_FUSIDIC = list(
      description = "Concomitant fusidic acid indicator, 1 = fusidic acid co-administered, 0 = rifampicin alone.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rifampicin alone)",
      notes = paste(
        "Time-varying. Table 3 'Interaction (1 - INT) x CLM/FDesRiF in the",
        "presence of FA' = 0.392; Results: 'In the presence of fusidic acid,",
        "clearance of 25-desacetyl rifampicin decreased by nearly 40%'. In the",
        "study, rifampicin was given alone only during the first 24 h of the",
        "six participants randomised to rifampicin as index drug; every other",
        "rifampicin record is under combination therapy. Table S2 run 1302",
        "prints the effect as '(1 - THETA(8) * (RF_only - 1))' without",
        "defining the coding of RF_only; the direction is taken from Table 3",
        "and the Results text."
      ),
      source_name = "FA (RF_only in Table S2)"
    ),
    OCC = list(
      description = "Study occasion for inter-occasion variability on 25-desacetyl rifampicin clearance (1 = Day 1, 2 = Day 2, 3 = Days 7-21, 4 = Day 28).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Table 3 'IOV in CLM/FDesRiF' = 71.2%. The occasions are the four",
        "sampling blocks the paper defines in Methods ('Day 1, Day 2, Days",
        "7-21 and Day 28') and Table S1 ('Occ 1: Day 0; Occ 2: Day 1; Occ 3:",
        "Days 7-21; Occ 4: Day 28', with the study's first dosing day counted",
        "as Day 0). Each occasion has its own eta; occasions 2-4 are fixed to",
        "the occasion-1 variance (the NONMEM BLOCK SAME convention). Pass OCC",
        "= 0 to switch inter-occasion variability off, e.g. for a",
        "typical-value simulation."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance.",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on CLM (Table S2 'Joint Screen'); not retained."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened on ka (Table S2 'Joint Screen'); not retained."
    ),
    HGB = list(
      description = "Haemoglobin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened on ka (Table S2 'Joint Screen'); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "median 59 years (IQR 48-68)",
    weight_range = "not reported; median BMI 31 kg/m2 (IQR 29-35)",
    sex_female_pct = 30,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Presumed or confirmed staphylococcal infection: osteoarticular infection",
      "of the axial skeleton (4), large-joint septic arthritis (4, including 2",
      "prosthetic joint infections), endovascular infection (1) and small-joint",
      "infection (1). Baseline median creatinine 78 umol/L (IQR 63-89), albumin",
      "27 g/L (21-35), CRP 44 (22-112)."
    ),
    dose_range = paste(
      "Rifampicin 300 mg orally twice daily with fusidic acid 500 mg orally",
      "three times daily. Randomised index drug for the first 24 h (rifampicin",
      "n = 6, fusidic acid n = 4), then the combination for the treatment",
      "course."
    ),
    regions = "Melbourne, Victoria, Australia (St Vincent's Hospital Melbourne and Northern Health), 2019-2022.",
    n_observations = "101 rifampicin and 96 25-desacetyl rifampicin concentrations.",
    notes = paste(
      "Maggs 2026 Results paragraph 1. Concentrations were measured by",
      "LC-MS/MS in capillary dried blood spots (LLOQ 0.015 mg/L rifampicin,",
      "0.014 mg/L 25-desacetyl rifampicin) on Day 1, Day 2 and Day 28 +/- 7",
      "(pre-dose, 30-40 min, 1.5-2 h and 3.5-4.5 h, plus 7-8 h on Days 1 and",
      "2; Table 1). NONMEM 7.5.1 (ADVAN6), FOCE-I."
    )
  )

  ini({
    # Maggs 2026 Table 3, final joint rifampicin / 25-desacetyl rifampicin
    # model (Table S2 run 1523, OFV 7.017). Typical values are for a 70 kg
    # individual taking rifampicin alone.
    ltlag <- log(0.512)
    label("Absorption lag time (h)") # Table 3 'ALAG' = 0.512 (RSE 0%; bootstrap 95% CI 0.487-0.582)
    lka <- log(3.11)
    label("Absorption rate constant ka (1/h)") # Table 3 'Ka, h-1' = 3.11 (RSE 62%)
    lcl <- log(14)
    label("Apparent rifampicin clearance CL/F at 70 kg (L/h)") # Table 3 'CL/FRif, L h-1.70 kg-1' = 14 (RSE 13%)
    lvc <- log(69.1)
    label("Apparent rifampicin volume V/F at 70 kg (L)") # Table 3 'V/FRif, l.70 kg-1' = 69.1 (RSE 26%)
    lcl_desacetylrif <- log(109)
    label("Apparent 25-desacetyl rifampicin clearance CLM/F at 70 kg without fusidic acid (L/h)") # Table 3 'CLM/FDesRiF, L h-1.70 kg-1' = 109 (RSE 13%)
    lvc_desacetylrif <- log(46.2)
    label("Apparent 25-desacetyl rifampicin volume VM/F at 70 kg (L)") # Table 3 'V/FDesRif, l.70 kg-1' = 46.2 (RSE 53%)
    fm_desacetylrif <- 0.019
    label("Fraction of the absorbed rifampicin dose converted to 25-desacetyl rifampicin pre-systemically, 1 - FP (fraction)") # Table 3 'FP' (fraction absorbed without first-pass metabolism) = 0.981 (RSE 2%); 1 - 0.981 = 0.019

    e_conmed_fusidic_cl_desacetylrif <- -0.392
    label("Fractional change in 25-desacetyl rifampicin CLM/F with co-administered fusidic acid (fraction)") # Table 3 'Interaction (1 - INT) x CLM/FDesRiF in the presence of FA' = 0.392 (RSE 25%)

    # IIV, IOV and residual error are printed as 100% x sqrt(variance
    # estimate) (Table 3 footnote), so each variance is the squared fraction.
    etalka ~ 1.21 # Table 3 'IIV in Ka' = 110% (shrinkage 12%); 1.10^2
    etalcl_desacetylrif ~ 0.09 # Table 3 'IIV in CLM/FDesRiF' = 30% (shrinkage 46%); 0.30^2
    etaiov_cl_desacetylrif_1 ~ 0.506944 # Table 3 'IOV in CLM/FDesRiF' = 71.2% (shrinkage 47%); 0.712^2
    etaiov_cl_desacetylrif_2 ~ fixed(0.506944) # shares the occasion-1 IOV variance
    etaiov_cl_desacetylrif_3 ~ fixed(0.506944) # shares the occasion-1 IOV variance
    etaiov_cl_desacetylrif_4 ~ fixed(0.506944) # shares the occasion-1 IOV variance

    propSd <- 0.66
    label("Proportional residual error SD, rifampicin (fraction)") # Table 3 'Rif, proportional error (%)' = 66.0 (shrinkage 7%)
    propSd_desacetylrif <- 0.351
    label("Proportional residual error SD, 25-desacetyl rifampicin (fraction)") # Table 3 'DesRif, proportional error (%)' = 35.1 (shrinkage 7%)
  })
  model({
    # Occasion indicators multiplex the IOV etas on metabolite clearance;
    # OCC = 0 switches IOV off.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl_desacetylrif <- oc1 * etaiov_cl_desacetylrif_1 +
      oc2 * etaiov_cl_desacetylrif_2 +
      oc3 * etaiov_cl_desacetylrif_3 +
      oc4 * etaiov_cl_desacetylrif_4

    tlag <- exp(ltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl) * (WT / 70)^0.75
    vc <- exp(lvc) * (WT / 70)
    cl_desacetylrif <- exp(lcl_desacetylrif + etalcl_desacetylrif + iov_cl_desacetylrif) *
      (WT / 70)^0.75 *
      (1 + e_conmed_fusidic_cl_desacetylrif * CONMED_FUSIDIC)
    vc_desacetylrif <- exp(lvc_desacetylrif) * (WT / 70)

    # All systemic rifampicin clearance forms 25-desacetyl rifampicin
    # (Results: 'under the assumption that all rifampicin is converted');
    # amounts are carried in rifampicin mass units with no molecular-weight
    # conversion.
    kel <- cl / vc
    kel_desacetylrif <- cl_desacetylrif / vc_desacetylrif

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * (1 - fm_desacetylrif) * depot - kel * central
    d/dt(central_desacetylrif) <- ka * fm_desacetylrif * depot + kel * central -
      kel_desacetylrif * central_desacetylrif
    alag(depot) <- tlag

    Cc <- central / vc
    Cc_desacetylrif <- central_desacetylrif / vc_desacetylrif

    Cc ~ prop(propSd)
    Cc_desacetylrif ~ prop(propSd_desacetylrif)
  })
}
