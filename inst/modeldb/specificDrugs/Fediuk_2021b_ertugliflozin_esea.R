Fediuk_2021b_ertugliflozin_esea <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ertugliflozin in healthy",
    "adults and adults with type 2 diabetes mellitus, refitted to the 15-study",
    "data set of Fediuk 2021 (2276 subjects, 13,692 concentrations) with an",
    "East/Southeast (E/SE) Asian versus non-E/SE Asian ethnicity covariate",
    "(Fediuk 2021b analysis 1, data set 1). First-order absorption with a lag",
    "time and first-order elimination. Allometric body-weight scaling fixed at",
    "0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F (85 kg reference). Covariate",
    "effects on CL/F (eGFR power term referenced to 90 mL/min/1.73 m^2 and",
    "capped at 120; T2DM, female sex, E/SE Asian ethnicity), on Vc/F (age",
    "power term referenced to 65 years; female sex, E/SE Asian ethnicity), and",
    "on absorption (fed and without-regard-to-food multipliers on ka and",
    "fractional decreases in relative bioavailability). IIV on CL/F only;",
    "log-scale additive residual error estimated separately for the phase 1",
    "and the phase 2/3 studies."
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
        "(Methods, Data Set 1). Cohort range 42.6-197 kg (Table 1)."
      ),
      source_name = "BWT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (AGE/65)^-0.192 on Vc/F only (Table S3). Reference",
        "65 years (Methods, Data Set 1). Cohort range 18-87 years (Table 1)."
      ),
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Baseline eGFR, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (eGFR/90)^0.458 on CL/F (Table S3). Reference",
        "90 mL/min/1.73 m^2; values above 120 mL/min/1.73 m^2 were fixed to",
        "120 (Methods, Data Set 1), applied in model() as min(CRCL, 120) so an",
        "uncapped eGFR may be supplied. The eGFR equation (4-variable MDRD) is",
        "inherited from the parent analysis (Fediuk 2021, doi:10.1002/cpdd.885).",
        "Cohort median 86.6, range 6.8-196 mL/min/1.73 m^2 (Table 1)."
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
        "orientation as the canonical. Multiplicative 0.921^DIS_DIAB on CL/F",
        "(Table 2). 2084 of 2276 subjects (91.6%) had T2DM (Table 1)."
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
        "canonical. Multiplicative 0.960^SEXF on CL/F and 1.46^SEXF on Vc/F",
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
        "Self-reported Asian ethnicity. Enters the model only through the",
        "paper's EASIA indicator, RACE_ASIAN * (REGION_EASTASIA +",
        "REGION_SOUTHEASTASIA): Asian subjects enrolled at US or European",
        "sites were categorized as non-E/SE Asian (Methods, Data Set 1), so",
        "RACE_ASIAN = 1 alone does not select the E/SE Asian effect."
      ),
      source_name = "EASIA (with the region columns)"
    ),
    REGION_EASTASIA = list(
      description = "Enrolled at an East Asian study site (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (enrolled outside East and Southeast Asia when REGION_SOUTHEASTASIA is also 0)",
      notes = paste(
        "East Asian sites of the phase 2/3 studies: Hong Kong, Republic of",
        "Korea and Taiwan (Methods, Data Set 1). The paper also counts every",
        "Japanese subject of the Japanese/Western phase 1 PK/PD study (Li",
        "2021, doi:10.1002/cpdd.908) as E/SE Asian; code those subjects",
        "REGION_EASTASIA = 1 with RACE_ASIAN = 1. Mutually exclusive with",
        "REGION_SOUTHEASTASIA."
      ),
      source_name = "EASIA (with RACE_ASIAN)"
    ),
    REGION_SOUTHEASTASIA = list(
      description = "Enrolled at a Southeast Asian study site (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (enrolled outside East and Southeast Asia when REGION_EASTASIA is also 0)",
      notes = paste(
        "Southeast Asian sites of the phase 2/3 studies: Malaysia, the",
        "Philippines and Thailand (Methods, Data Set 1). Mutually exclusive",
        "with REGION_EASTASIA; the two are summed into the E/SE Asian region",
        "term of model()."
      ),
      source_name = "EASIA (with RACE_ASIAN)"
    ),
    FED = list(
      description = "Fed-state dose-record indicator (1 = administered with food, 0 = fasted or food status not documented)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, when FED_MISSING is also 0)",
      notes = paste(
        "Source FOODEFF1 of the parent analysis (Fediuk 2021, Supplementary",
        "Appendix Equation 8). Per dose record. Multiplicative 0.725^FED on ka",
        "and a fractional decrease F1 = 1 - 0.0702 * FED on relative",
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
        "Source FOODEFF2 of the parent analysis. The phase 3 studies dosed",
        "'without regard to food' and did not record food status, so this is",
        "the canonical FED_MISSING level rather than a meal state.",
        "Multiplicative 0.654^FED_MISSING on ka and a fractional decrease",
        "F1 = 1 - 0.0681 * FED_MISSING (Table 2). Set to 0 when simulating a",
        "defined fasted or fed state."
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
        "studies share one log-scale residual SD (0.837) and the phase 1",
        "studies another (0.389) (Table 2). Use 0 (with STUDY_PHASE3 = 0) to",
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
    n_subjects = 2276L,
    n_studies = 15L,
    n_observations = 13692L,
    age_range = "18-87 years",
    age_median = "57.0 years (mean 55.7, SD 11.6)",
    weight_range = "42.6-197 kg",
    weight_median = "84.8 kg (mean 86.9, SD 19.7)",
    sex_female_pct = 43.5,
    race_ethnicity = c(`E/SE Asian` = 6.8, `non-E/SE Asian` = 93.2),
    disease_state = paste(
      "2084 adults with type 2 diabetes mellitus (91.6%) and 192 healthy",
      "subjects (8.44%)."
    ),
    renal_function = "eGFR median 86.6 (range 6.8-196) mL/min/1.73 m^2 (Table 1).",
    subgroups = paste(
      "E/SE Asian (N = 154): body weight mean 69.7 kg (SD 14.6), age 54.1",
      "years, 48.7% female, 92.2% T2DM. Non-E/SE Asian (N = 2122): body",
      "weight mean 88.2 kg (SD 19.4), age 55.8 years, 43.1% female, 91.5%",
      "T2DM (Table 1)."
    ),
    dose_range = "Single and multiple oral doses across 10 phase 1, 2 phase 2 and 3 phase 3 studies of the parent analysis; phase 3 regimens 5 and 15 mg once daily.",
    regions = "Multinational; E/SE Asian sites in Hong Kong, Republic of Korea, Malaysia, the Philippines, Taiwan and Thailand.",
    notes = paste(
      "Same analysis data set as the parent popPK analysis (Fediuk 2021,",
      "doi:10.1002/cpdd.885); this refit replaces the White/Black/Asian/other",
      "race covariates by the E/SE Asian indicator. LLOQ 0.5 ng/mL (0.1 ng/mL",
      "in one study); BLQ records removed. Baseline demographics from Table 1;",
      "final estimates from Table 2; covariate equations from Table S3."
    )
  )

  # The phase 1 and phase 2/3 studies carry separately estimated residual
  # SDs; nlmixr2 takes one residual symbol per output, so they are combined
  # into expSdCc inside model() by the study-phase indicators (the
  # Fediuk_2021_ertugliflozin.R pattern).
  paper_specific_residual_sds <- c("expSdPh1", "expSdPh23")

  ini({
    # Structural parameters for the reference subject: 65-year-old healthy
    # non-E/SE Asian man, 85 kg, eGFR 90 mL/min/1.73 m^2, fasted (Methods,
    # Data Set 1; Table S3).
    lcl <- log(11.9); label("Apparent clearance CL/F (L/h)") # Table 2 data set 1 'CL/F, L/h' 11.9
    lvc <- log(6.51); label("Apparent central volume Vc/F (L)") # Table 2 data set 1 'Vc/F, L' 6.51
    lvp <- log(107); label("Apparent peripheral volume Vp/F (L)") # Table 2 data set 1 'Vp/F, L' 107
    lq <- log(7.76); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 data set 1 'Q/F, L/h' 7.76
    lka <- log(0.329); label("Absorption rate constant, fasted (1/h)") # Table 2 data set 1 'ka, h-1' 0.329
    ltlag <- log(0.226); label("Absorption lag time (h)") # Table 2 data set 1 'Lag time (ALAG1), h' 0.226
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, fasted (unitless)") # Table 2 data set 1 'Relative bioavailability (F1)' 1.00 (fixed)

    # Allometric exponents (fixed)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of (WT/85) on CL/F and Q/F (unitless)") # Table 2 data set 1 'Effect of body weight' 0.750 (fixed), CL/F and Q/F rows
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of (WT/85) on Vc/F and Vp/F (unitless)") # Table 2 data set 1 'Effect of body weight' 1.00 (fixed), Vc/F and Vp/F rows

    # Covariate effects on CL/F (Table S3, data set 1)
    e_crcl_cl <- 0.458; label("Power exponent of (eGFR/90) on CL/F (unitless)") # Table 2 data set 1 CL/F 'Effect of eGFR' 0.458
    e_dis_diab_cl <- 0.921; label("Multiplicative effect of T2DM on CL/F (ratio)") # Table 2 data set 1 CL/F 'Effect of T2DM patient status' 0.921
    e_sexf_cl <- 0.960; label("Multiplicative effect of female sex on CL/F (ratio)") # Table 2 data set 1 CL/F 'Effect of female sex' 0.960
    e_asian_esea_cl <- 1.17; label("Multiplicative effect of E/SE Asian ethnicity on CL/F (ratio)") # Table 2 data set 1 CL/F 'Effect of E/SE Asian ethnicity' 1.17

    # Covariate effects on Vc/F (Table S3, data set 1)
    e_age_vc <- -0.192; label("Power exponent of (AGE/65) on Vc/F (unitless)") # Table 2 data set 1 Vc/F 'Effect of age' -0.192
    e_sexf_vc <- 1.46; label("Multiplicative effect of female sex on Vc/F (ratio)") # Table 2 data set 1 Vc/F 'Effect of female sex' 1.46
    e_asian_esea_vc <- 2.48; label("Multiplicative effect of E/SE Asian ethnicity on Vc/F (ratio)") # Table 2 data set 1 Vc/F 'Effect of E/SE Asian ethnicity' 2.48

    # Food effects on absorption (parent analysis Supplementary Appendix
    # Equations 7-8, same parameterization)
    e_fed_ka <- 0.725; label("Multiplicative effect of fed state on ka (ratio)") # Table 2 data set 1 ka 'Effect of food' 0.725
    e_fed_missing_ka <- 0.654; label("Multiplicative effect of without-regard-to-food dosing on ka (ratio)") # Table 2 data set 1 ka 'Effect of without regard to food' 0.654
    e_fed_fdepot <- 0.0702; label("Fractional decrease in F1 with food (fraction)") # Table 2 data set 1 F1 'Effect of food' 0.0702
    e_fed_missing_fdepot <- 0.0681; label("Fractional decrease in F1 with without-regard-to-food dosing (fraction)") # Table 2 data set 1 F1 'Effect of without regard to food' 0.0681

    # IIV on CL/F only
    etalcl ~ 0.101 # Table 2 data set 1 'omega2 CL/F' 0.101 (31.8% CV in Results)

    # Log-scale additive residual error, reported on the SD scale (Results:
    # 'Residual error estimates were 38.9% ... and 83.7%').
    expSdPh1 <- 0.389; label("Log-scale residual SD, phase 1 studies (unitless)") # Table 2 data set 1 'Phase 1 residual error' 0.389
    expSdPh23 <- 0.837; label("Log-scale residual SD, phase 2/3 studies (unitless)") # Table 2 data set 1 'Phase 2/3 residual error' 0.837
  })
  model({
    # eGFR above 120 mL/min/1.73 m^2 was fixed to 120 (Methods, Data Set 1).
    egfr_capped <- min(CRCL, 120)

    # Table S3 EASIA: Asian subjects enrolled at E/SE Asian sites (or in the
    # Japanese phase 1 study); Asian subjects at US/European sites are 0.
    easia <- RACE_ASIAN * (REGION_EASTASIA + REGION_SOUTHEASTASIA)

    cl <- exp(lcl + etalcl) * (WT / 85)^e_wt_cl_q * (egfr_capped / 90)^e_crcl_cl *
      e_dis_diab_cl^DIS_DIAB * e_sexf_cl^SEXF * e_asian_esea_cl^easia
    vc <- exp(lvc) * (WT / 85)^e_wt_vc_vp * (AGE / 65)^e_age_vc * e_sexf_vc^SEXF *
      e_asian_esea_vc^easia
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
