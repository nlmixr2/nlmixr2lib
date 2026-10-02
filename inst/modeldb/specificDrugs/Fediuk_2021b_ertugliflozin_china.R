Fediuk_2021b_ertugliflozin_china <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ertugliflozin in healthy",
    "adults and adults with type 2 diabetes mellitus, fitted to 17 phase 1-3",
    "studies (2620 subjects, 16,018 concentrations) with ethnicity covariates",
    "for Asian subjects from mainland China and Asian subjects from the rest",
    "of the world versus non-Asian subjects (Fediuk 2021b analysis 2, data",
    "set 2). First-order absorption with a lag time and first-order",
    "elimination. Allometric body-weight scaling fixed at 0.75 on CL/F and Q/F",
    "and 1 on Vc/F and Vp/F (85 kg reference). Covariate effects on CL/F (eGFR",
    "power term referenced to 90 mL/min/1.73 m^2 and capped at 120; T2DM,",
    "female sex, the two Asian groups), on Vc/F (age power term referenced to",
    "65 years; female sex, the two Asian groups), and on absorption (fed and",
    "without-regard-to-food multipliers on ka and fractional decreases in",
    "relative bioavailability). IIV on CL/F only; log-scale additive residual",
    "error estimated separately for the phase 1 and the phase 2/3 studies."
  )
  reference <- paste(
    "Fediuk DJ, Sahasrabudhe V, Dawra VK, Zhou S, Sweeney K (2021).",
    "Population Pharmacokinetic Analyses of Ertugliflozin in Select Ethnic",
    "Populations. Clinical Pharmacology in Drug Development",
    "10(11):1297-1306. doi:10.1002/cpdd.970. Absorption covariate",
    "parameterization from Fediuk DJ et al. (2021) Clin Pharmacol Drug Dev",
    "10(7):696-706, doi:10.1002/cpdd.885 (Supplementary Appendix Equations",
    "7-8)."
  )
  vignette <- "Fediuk_2021b_ertugliflozin_ethnicity"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "ertugliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ertugliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ertugliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value. Power scaling (WT/85)^0.75 on CL/F and",
        "Q/F and (WT/85)^1 on Vc/F and Vp/F, exponents fixed (Table 2",
        "'0.750 (fixed)' and '1.00 (fixed)' rows; Table S3). Reference 85 kg",
        "(Methods, Data Set 2). Cohort range 42.6-197 kg (Table 1)."
      ),
      source_name = "BWT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (AGE/65)^-0.603 on Vc/F only (Table S3). Reference",
        "65 years (Methods, Data Set 2). Cohort range 18-87 years (Table 1)."
      ),
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Baseline eGFR, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (eGFR/90)^0.449 on CL/F (Table S3). Reference",
        "90 mL/min/1.73 m^2; values above 120 mL/min/1.73 m^2 were fixed to",
        "120 (Methods, Data Set 2), applied in model() as min(CRCL, 120) so an",
        "uncapped eGFR may be supplied. The eGFR equation (4-variable MDRD) is",
        "inherited from the parent analysis (Fediuk 2021, doi:10.1002/cpdd.885).",
        "Cohort median 88.0, range 6.8-196 mL/min/1.73 m^2 (Table 1)."
      ),
      source_name = "eGFR"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus patient status (1 = T2DM, 0 = healthy subject)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject)",
      notes = paste(
        "Source PTST (healthy subjects = 0, T2DM patients = 1; Table S3), same",
        "orientation as the canonical. Multiplicative 0.856^DIS_DIAB on CL/F",
        "(Table 2). 2412 of 2620 subjects (92.1%) had T2DM (Table 1)."
      ),
      source_name = "PTST"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Source NSEX (male = 0, female = 1; Table S3), same orientation as the",
        "canonical. Multiplicative 0.970^SEXF on CL/F and 1.66^SEXF on Vc/F",
        "(Table 2)."
      ),
      source_name = "NSEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race/ethnicity indicator (1 = Asian, 0 = not Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = paste(
        "Self-reported Asian ethnicity. The paper's two indicators are",
        "RACE2 = RACE_ASIAN * REGION_CHINA (Asian subjects from mainland",
        "China, N = 277) and RACE3 = RACE_ASIAN * (1 - REGION_CHINA) (Asian",
        "subjects from the rest of the world, N = 382); non-Asian subjects",
        "(N = 1961) are the reference whatever their enrollment region",
        "(Methods, Data Set 2; Table S3)."
      ),
      source_name = "RACE2 / RACE3 (with REGION_CHINA)"
    ),
    REGION_CHINA = list(
      description = "From mainland China (1 = enrolled in mainland China, 0 = elsewhere)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rest of the world)",
      notes = paste(
        "Mainland China only (Hong Kong and Taiwan are rest of the world).",
        "Contributed by the phase 1 study in healthy Chinese subjects (study",
        "I-10) and the Chinese sites of the phase 3 study in Asian patients",
        "with T2DM (study III-5; Table S1). Acts only together with",
        "RACE_ASIAN: multiplicative 1.04 on CL/F and 1.44 on Vc/F for",
        "RACE_ASIAN * REGION_CHINA, against 1.08 and 2.15 for Asian subjects",
        "from the rest of the world (Table 2)."
      ),
      source_name = "RACE2 (with RACE_ASIAN)"
    ),
    FED = list(
      description = "Fed-state dose-record indicator (1 = administered with food, 0 = fasted or food status not documented)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, when FED_MISSING is also 0)",
      notes = paste(
        "Source FOODEFF1 of the parent analysis (Fediuk 2021, Supplementary",
        "Appendix Equation 8). Per dose record. Multiplicative 0.639^FED on ka",
        "and a fractional decrease F1 = 1 - 0.157 * FED on relative",
        "bioavailability (Table 2). Mutually exclusive with FED_MISSING."
      ),
      source_name = "FOODEFF1"
    ),
    FED_MISSING = list(
      description = "Food-status-not-documented dose-record indicator ('without regard to food')",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (food status recorded)",
      notes = paste(
        "Source FOODEFF2 of the parent analysis. The phase 3 studies (including",
        "study III-5, Table S1) dosed 'without regard to food' and did not",
        "record food status, so this is the canonical FED_MISSING level rather",
        "than a meal state. Multiplicative 0.796^FED_MISSING on ka and a",
        "fractional decrease F1 = 1 - 0.148 * FED_MISSING (Table 2). Set to 0",
        "when simulating a defined fasted or fed state."
      ),
      source_name = "FOODEFF2"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2 study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 when STUDY_PHASE3 is also 0)",
      notes = paste(
        "Selects the residual-error magnitude only: the phase 2 and phase 3",
        "studies share one log-scale residual SD (0.831) and the phase 1",
        "studies another (0.486) (Table 2). Use 0 (with STUDY_PHASE3 = 0) to",
        "simulate a richly sampled phase 1 profile."
      ),
      source_name = "study phase (not a published column)"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study-stratum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 when STUDY_PHASE2 is also 0)",
      notes = paste(
        "Selects the residual-error magnitude only, pooled with phase 2",
        "(Table 2 row 'Phase 2/3 residual error'). Mutually exclusive with",
        "STUDY_PHASE2."
      ),
      source_name = "study phase (not a published column)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2620L,
    n_studies = 17L,
    n_observations = 16018L,
    age_range = "18-87 years",
    age_median = "57.0 years (mean 55.6, SD 11.5)",
    weight_range = "42.6-197 kg",
    weight_median = "82.4 kg (mean 84.7, SD 19.6)",
    sex_female_pct = 43.4,
    race_ethnicity = c(
      `Asian from mainland China` = 10.6,
      `Asian from rest of world` = 14.6,
      `non-Asian` = 74.8
    ),
    disease_state = paste(
      "2412 adults with type 2 diabetes mellitus (92.1%) and 208 healthy",
      "subjects (7.94%)."
    ),
    renal_function = "eGFR median 88.0 (range 6.8-196) mL/min/1.73 m^2 (Table 1).",
    subgroups = paste(
      "Asian from mainland China (N = 277): body weight mean 70.5 kg (SD",
      "10.7), age 54.5 years, eGFR median 101, 43.0% female, 94.2% T2DM.",
      "Asian from the rest of the world (N = 382): body weight mean 72.8 kg",
      "(SD 15.4), age 54.9 years, eGFR median 87.2, 40.8% female, 96.1% T2DM.",
      "Non-Asian (N = 1961): body weight mean 89.1 kg (SD 19.4), age 55.8",
      "years, eGFR median 86.2, 44.0% female, 91.0% T2DM (Table 1)."
    ),
    dose_range = paste(
      "Single and multiple oral doses across 10 phase 1, 2 phase 2 and 5",
      "phase 3 studies; the two studies added to the parent data set dosed",
      "5 and 15 mg tablets (Table S1)."
    ),
    regions = "Multinational, including mainland China (studies I-10 and III-5).",
    notes = paste(
      "The 15-study data set of the parent analysis (Fediuk 2021,",
      "doi:10.1002/cpdd.885) plus a phase 1 single- and multiple-dose study in",
      "healthy Chinese subjects (I-10, N = 16) and a phase 3 study in Asian",
      "patients with T2DM (III-5, NCT02630706, N = 330) (Table S1). LLOQ",
      "0.5 ng/mL (0.1 ng/mL in one study); BLQ records removed. Baseline",
      "demographics from Table 1; final estimates from Table 2; covariate",
      "equations from Table S3."
    )
  )

  # The phase 1 and phase 2/3 studies carry separately estimated residual
  # SDs; nlmixr2 takes one residual symbol per output, so they are combined
  # into expSdCc inside model() by the study-phase indicators (the
  # Fediuk_2021_ertugliflozin.R pattern).
  paper_specific_residual_sds <- c("expSdPh1", "expSdPh23")

  ini({
    # Structural parameters for the reference subject: 65-year-old healthy
    # non-Asian man, 85 kg, eGFR 90 mL/min/1.73 m^2, fasted (Methods, Data
    # Set 2; Table S3).
    lcl <- log(11.8); label("Apparent clearance CL/F (L/h)") # Table 2 data set 2 'CL/F, L/h' 11.8
    lvc <- log(5.23); label("Apparent central volume Vc/F (L)") # Table 2 data set 2 'Vc/F, L' 5.23
    lvp <- log(113); label("Apparent peripheral volume Vp/F (L)") # Table 2 data set 2 'Vp/F, L' 113
    lq <- log(7.17); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 data set 2 'Q/F, L/h' 7.17
    lka <- log(0.323); label("Absorption rate constant, fasted (1/h)") # Table 2 data set 2 'ka, h-1' 0.323
    ltlag <- log(0.227); label("Absorption lag time (h)") # Table 2 data set 2 'Lag time (ALAG1), h' 0.227
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, fasted (unitless)") # Table 2 data set 2 'Relative bioavailability (F1)' 1.00 (fixed)

    # Allometric exponents (fixed)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of (WT/85) on CL/F and Q/F (unitless)") # Table 2 data set 2 'Effect of body weight' 0.750 (fixed), CL/F and Q/F rows
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of (WT/85) on Vc/F and Vp/F (unitless)") # Table 2 data set 2 'Effect of body weight' 1.00 (fixed), Vc/F and Vp/F rows

    # Covariate effects on CL/F (Table S3, data set 2)
    e_crcl_cl <- 0.449; label("Power exponent of (eGFR/90) on CL/F (unitless)") # Table 2 data set 2 CL/F 'Effect of eGFR' 0.449
    e_dis_diab_cl <- 0.856; label("Multiplicative effect of T2DM on CL/F (ratio)") # Table 2 data set 2 CL/F 'Effect of T2DM patient status' 0.856
    e_sexf_cl <- 0.970; label("Multiplicative effect of female sex on CL/F (ratio)") # Table 2 data set 2 CL/F 'Effect of female sex' 0.970
    e_asian_china_cl <- 1.04; label("Multiplicative effect of Asian from mainland China on CL/F (ratio)") # Table 2 data set 2 CL/F 'Effect of Asian from mainland China' 1.04
    e_asian_row_cl <- 1.08; label("Multiplicative effect of Asian from the rest of the world on CL/F (ratio)") # Table 2 data set 2 CL/F 'Effect of Asian from ROW' 1.08

    # Covariate effects on Vc/F (Table S3, data set 2)
    e_age_vc <- -0.603; label("Power exponent of (AGE/65) on Vc/F (unitless)") # Table 2 data set 2 Vc/F 'Effect of age' -0.603
    e_sexf_vc <- 1.66; label("Multiplicative effect of female sex on Vc/F (ratio)") # Table 2 data set 2 Vc/F 'Effect of female sex' 1.66
    e_asian_china_vc <- 1.44; label("Multiplicative effect of Asian from mainland China on Vc/F (ratio)") # Table 2 data set 2 Vc/F 'Effect of Asian from mainland China' 1.44
    e_asian_row_vc <- 2.15; label("Multiplicative effect of Asian from the rest of the world on Vc/F (ratio)") # Table 2 data set 2 Vc/F 'Effect of Asian from ROW' 2.15

    # Food effects on absorption (parent analysis Supplementary Appendix
    # Equations 7-8, same parameterization)
    e_fed_ka <- 0.639; label("Multiplicative effect of fed state on ka (ratio)") # Table 2 data set 2 ka 'Effect of food' 0.639
    e_fed_missing_ka <- 0.796; label("Multiplicative effect of without-regard-to-food dosing on ka (ratio)") # Table 2 data set 2 ka 'Effect of without regard to food' 0.796
    e_fed_fdepot <- 0.157; label("Fractional decrease in F1 with food (fraction)") # Table 2 data set 2 F1 'Effect of food' 0.157
    e_fed_missing_fdepot <- 0.148; label("Fractional decrease in F1 with without-regard-to-food dosing (fraction)") # Table 2 data set 2 F1 'Effect of without regard to food' 0.148

    # IIV on CL/F only
    etalcl ~ 0.0999 # Table 2 data set 2 'omega2 CL/F' 0.0999 (31.6% CV in Results)

    # Log-scale additive residual error, reported on the SD scale (Results:
    # 'Residual error estimates expressed as a coefficient of variation were
    # 48.6% ... and 83.1%').
    expSdPh1 <- 0.486; label("Log-scale residual SD, phase 1 studies (unitless)") # Table 2 data set 2 'Phase 1 residual error' 0.486
    expSdPh23 <- 0.831; label("Log-scale residual SD, phase 2/3 studies (unitless)") # Table 2 data set 2 'Phase 2/3 residual error' 0.831
  })
  model({
    # eGFR above 120 mL/min/1.73 m^2 was fixed to 120 (Methods, Data Set 2).
    egfr_capped <- min(CRCL, 120)

    # Table S3 RACE2 (Asian from mainland China) and RACE3 (Asian from the
    # rest of the world); non-Asian subjects are the reference.
    asian_china <- RACE_ASIAN * REGION_CHINA
    asian_row <- RACE_ASIAN * (1 - REGION_CHINA)

    cl <- exp(lcl + etalcl) * (WT / 85)^e_wt_cl_q * (egfr_capped / 90)^e_crcl_cl *
      e_dis_diab_cl^DIS_DIAB * e_sexf_cl^SEXF *
      e_asian_china_cl^asian_china * e_asian_row_cl^asian_row
    vc <- exp(lvc) * (WT / 85)^e_wt_vc_vp * (AGE / 65)^e_age_vc * e_sexf_vc^SEXF *
      e_asian_china_vc^asian_china * e_asian_row_vc^asian_row
    vp <- exp(lvp) * (WT / 85)^e_wt_vc_vp
    q <- exp(lq) * (WT / 85)^e_wt_cl_q
    ka <- exp(lka) * e_fed_ka^FED * e_fed_missing_ka^FED_MISSING
    tlag <- exp(ltlag)
    # Parent analysis Supplementary Appendix Equation 8:
    # F1 = 1 - FOODEFF1*theta - FOODEFF2*theta
    fdepot <- exp(lfdepot) * (1 - e_fed_fdepot * FED - e_fed_missing_fdepot * FED_MISSING)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc

    ph23 <- STUDY_PHASE2 + STUDY_PHASE3
    expSdCc <- expSdPh1 * (1 - ph23) + expSdPh23 * ph23
    Cc ~ lnorm(expSdCc)
  })
}
