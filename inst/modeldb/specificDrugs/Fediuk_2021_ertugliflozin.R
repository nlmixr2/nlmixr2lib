Fediuk_2021_ertugliflozin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ertugliflozin in healthy",
    "adults and adults with type 2 diabetes mellitus, pooled across 15",
    "phase 1-3 studies (Fediuk 2021; 2276 subjects, 13,691 concentrations).",
    "First-order absorption with a lag time and first-order elimination.",
    "Allometric body-weight scaling fixed at 0.75 on CL/F and Q/F and 1 on",
    "Vc/F and Vp/F (85 kg reference). Full-model covariate effects on CL/F",
    "(eGFR power term referenced to 90 mL/min/1.73 m^2 and capped at 120;",
    "T2DM, female sex, Black, Asian and other race), on Vc/F (age power",
    "term referenced to 65 years; female sex, Black, Asian and other race),",
    "and on absorption (fed and without-regard-to-food multipliers on ka",
    "and fractional decreases in relative bioavailability). IIV on CL/F",
    "only; log-scale additive residual error estimated separately for the",
    "phase 1 and the phase 2/3 studies."
  )
  reference <- paste(
    "Fediuk DJ, Zhou S, Dawra VK, Sahasrabudhe V, Sweeney K (2021).",
    "Population Pharmacokinetic Model for Ertugliflozin in Healthy Subjects",
    "and Patients With Type 2 Diabetes Mellitus. Clinical Pharmacology in",
    "Drug Development 10(7):696-706. doi:10.1002/cpdd.885."
  )
  vignette <- "Fediuk_2021_ertugliflozin"
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
        "Q/F and (WT/85)^1 on Vc/F and Vp/F, exponents fixed (Table 2 'FIX'",
        "rows; Supplementary Appendix Equation 8). Reference 85 kg 'based",
        "upon the population median of 84.8 kg' (Supplementary Appendix,",
        "Covariate Evaluation). Cohort range 42.6-197 kg (Table 1)."
      ),
      source_name = "BWT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (AGE/65)^-0.243 on Vc/F only (Supplementary Appendix",
        "Equation 8). Reference 65 years 'the minimum age considered as",
        "elderly'. Age was deliberately not tested on CL/F because it was",
        "collinear with eGFR (Methods, Covariate Evaluation). Cohort range",
        "18-87 years (Table 1)."
      ),
      source_name = "AGE"
    ),
    CRCL = list(
      description = paste(
        "Baseline eGFR by the 4-variable MDRD equation, BSA-normalized"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (eGFR/90)^0.455 on CL/F (Supplementary Appendix",
        "Equation 8). Reference 90 mL/min/1.73 m^2 'the minimum value",
        "considered to be normal renal function'. The authors set eGFR",
        "values above 120 mL/min/1.73 m^2 to 120 in the analysis data set",
        "(Methods, Covariate Evaluation, citing hyperfiltration in early",
        "T2DM); model() applies that cap as min(CRCL, 120), so an",
        "uncapped eGFR may be supplied. Cohort median 86.6, range",
        "6.8-196 mL/min/1.73 m^2 (Table 1)."
      ),
      source_name = "eGFR"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus patient status (1 = T2DM, 0 = healthy subject)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject)",
      notes = paste(
        "Source PTST (healthy subjects = 0, T2DM = 1; Supplementary Appendix",
        "Equation 8), same orientation as the canonical. Type 2 diabetes",
        "specifically. Multiplicative 0.904^DIS_DIAB on CL/F (Table 2).",
        "2084 of 2276 subjects (91.6%) had T2DM (Table 1)."
      ),
      source_name = "PTST"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Source NSEX (male = 0, female = 1; Supplementary Appendix",
        "Equation 8), same orientation as the canonical. Multiplicative",
        "0.962^SEXF on CL/F and 1.36^SEXF on Vc/F (Table 2)."
      ),
      source_name = "NSEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = not Black)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference race is White)",
      notes = paste(
        "Source RACE2 (Supplementary Appendix Equation 8). Multiplicative",
        "0.985^RACE_BLACK on CL/F and 0.917^RACE_BLACK on Vc/F (Table 2).",
        "Mutually exclusive with RACE_ASIAN and RACE_OTHER; all three 0 for",
        "White subjects."
      ),
      source_name = "RACE2"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = not Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference race is White)",
      notes = paste(
        "Source RACE3 (Supplementary Appendix Equation 8). Multiplicative",
        "1.08^RACE_ASIAN on CL/F and 2.12^RACE_ASIAN on Vc/F (Table 2)."
      ),
      source_name = "RACE3"
    ),
    RACE_OTHER = list(
      description = "Other race indicator (1 = other race, 0 = White, Black or Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference race is White)",
      notes = paste(
        "Source RACE4 (Supplementary Appendix Equation 8). Multiplicative",
        "0.992^RACE_OTHER on CL/F and 1.15^RACE_OTHER on Vc/F (Table 2)."
      ),
      source_name = "RACE4"
    ),
    FED = list(
      description = "Fed-state dose-record indicator (1 = administered with food, 0 = fasted or food status not documented)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, when FED_MISSING is also 0)",
      notes = paste(
        "Source FOODEFF1 (Supplementary Appendix Equation 8). Per dose",
        "record: some phase 1 studies dosed the same subjects fasted and fed",
        "in different periods (Table 1 footnote a). The phase 2 studies",
        "dosed with the morning meal (Methods, Covariate Evaluation).",
        "Multiplicative 0.726^FED on ka and a fractional decrease",
        "F1 = 1 - 0.0683 * FED on relative bioavailability (Table 2,",
        "Supplementary Appendix Equation 7). Mutually exclusive with",
        "FED_MISSING."
      ),
      source_name = "FOODEFF1"
    ),
    FED_MISSING = list(
      description = "Food-status-not-documented dose-record indicator ('without regard to food')",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (food status recorded)",
      notes = paste(
        "Source FOODEFF2. In the phase 3 studies 'food status was not",
        "documented and per protocol ertugliflozin could be administered",
        "without regard to food' (Methods, Covariate Evaluation), so this",
        "stratum is the canonical FED_MISSING level rather than a",
        "physiological meal state. Multiplicative 0.663^FED_MISSING on ka",
        "and a fractional decrease F1 = 1 - 0.0809 * FED_MISSING (Table 2).",
        "Set to 0 when simulating a defined fasted or fed state; set to 1",
        "only to reproduce the phase 3 outpatient setting. FED must be 0",
        "on any record with FED_MISSING = 1."
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
        "studies share one log-scale residual SD (0.836) and the phase 1",
        "studies another (0.387) (Table 2; Methods, PopPK Model). It",
        "touches no structural or covariate parameter. Use 0 (with",
        "STUDY_PHASE3 = 0) to simulate a richly sampled phase 1 profile."
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
    age_range = "18-87 years",
    age_median = "57.0 years (mean 55.7, SD 11.6)",
    weight_range = "42.6-197 kg",
    weight_median = "84.8 kg (mean 86.9, SD 19.7)",
    sex_female_pct = 43.5,
    race_ethnicity = c(White = 71.8, Black = 8.74, Asian = 13.8, Other = 5.62),
    disease_state = paste(
      "2084 adults with type 2 diabetes mellitus (91.6%) and 192 healthy",
      "subjects (8.4%)."
    ),
    renal_function = paste(
      "MDRD eGFR median 86.6 (range 6.8-196) mL/min/1.73 m^2; about 44%",
      "normal (>= 90), 41% mild (60 to < 90), 14% moderate (30 to < 60)",
      "and 1% severe (< 30) renal impairment."
    ),
    dose_range = paste(
      "Oral solution / suspension or tablets; single doses 0.5-300 mg and",
      "once- or twice-daily multiple doses; phase 3 regimens 5 and 15 mg",
      "once daily (Table S1)."
    ),
    food_status = "Fed 473 (19.3%), without regard to food 1697 (69.4%), remainder fasted (Table 1).",
    regions = "Multinational (Pfizer / Merck ertugliflozin development program).",
    notes = paste(
      "9 phase 1 (rich sampling), 2 phase 2 and 4 phase 3 (sparse",
      "sampling) studies (Tables S1, S2). LLOQ 0.5 ng/mL (0.1 ng/mL in one",
      "study); 848 BLQ records (5%) removed. NONMEM 7.3 FOCEI on",
      "log-transformed concentrations. Baseline demographics from Table 1;",
      "final estimates from Table 2; covariate equations from the",
      "Supplementary Appendix (Equations 5-8)."
    )
  )

  # The phase 1 and phase 2/3 studies carry separately estimated residual
  # SDs; nlmixr2 takes one residual symbol per output, so they are combined
  # into expSdCc inside model() by the study-phase indicators (the
  # Rich_2026_momelotinib.R / Ooi_2026_elafibranor.R pattern).
  paper_specific_residual_sds <- c("expSdPh1", "expSdPh23")

  ini({
    # Structural parameters for the reference subject: healthy White male,
    # 85 kg, 65 years, eGFR 90 mL/min/1.73 m^2, fasted (Results, Parameter
    # Estimate Results).
    lcl <- log(12.0); label("Apparent clearance CL/F (L/h)") # Table 2 final model 'CL/F (L/h)' 12.0
    lvc <- log(6.54); label("Apparent central volume Vc/F (L)") # Table 2 final model 'Vc/F (L)' 6.54
    lvp <- log(107); label("Apparent peripheral volume Vp/F (L)") # Table 2 final model 'Vp/F (L)' 107
    lq <- log(7.77); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 final model 'Q/F (L/h)' 7.77
    lka <- log(0.329); label("Absorption rate constant, fasted (1/h)") # Table 2 final model 'ka (h-1)' 0.329
    ltlag <- log(0.228); label("Absorption lag time (h)") # Table 2 final model 'Lag time (h)' 0.228
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, fasted (unitless)") # Table 2 'Relative bioavailability (F1)' 1.00 FIX

    # Allometric exponents (fixed a priori)
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of (WT/85) on CL/F and Q/F (unitless)") # Table 2 'Effect of body weight' 0.750 FIX (CL/F and Q/F rows)
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of (WT/85) on Vc/F and Vp/F (unitless)") # Table 2 'Effect of body weight' 1.00 FIX (Vc/F and Vp/F rows)

    # Covariate effects on CL/F (Supplementary Appendix Equation 8)
    e_crcl_cl <- 0.455; label("Power exponent of (eGFR/90) on CL/F (unitless)") # Table 2 CL/F 'Effect of eGFR' 0.455
    e_dis_diab_cl <- 0.904; label("Multiplicative effect of T2DM on CL/F (ratio)") # Table 2 CL/F 'Effect of T2DM patient status' 0.904
    e_sexf_cl <- 0.962; label("Multiplicative effect of female sex on CL/F (ratio)") # Table 2 CL/F 'Effect of female sex' 0.962
    e_race_black_cl <- 0.985; label("Multiplicative effect of Black race on CL/F (ratio)") # Table 2 CL/F 'Effect of Black race' 0.985
    e_race_asian_cl <- 1.08; label("Multiplicative effect of Asian race on CL/F (ratio)") # Table 2 CL/F 'Effect of Asian race' 1.08
    e_race_other_cl <- 0.992; label("Multiplicative effect of other race on CL/F (ratio)") # Table 2 CL/F 'Effect of other race' 0.992

    # Covariate effects on Vc/F (Supplementary Appendix Equation 8)
    e_age_vc <- -0.243; label("Power exponent of (AGE/65) on Vc/F (unitless)") # Table 2 Vc/F 'Effect of age' -0.243
    e_sexf_vc <- 1.36; label("Multiplicative effect of female sex on Vc/F (ratio)") # Table 2 Vc/F 'Effect of female sex' 1.36
    e_race_black_vc <- 0.917; label("Multiplicative effect of Black race on Vc/F (ratio)") # Table 2 Vc/F 'Effect of Black race' 0.917
    e_race_asian_vc <- 2.12; label("Multiplicative effect of Asian race on Vc/F (ratio)") # Table 2 Vc/F 'Effect of Asian race' 2.12
    e_race_other_vc <- 1.15; label("Multiplicative effect of other race on Vc/F (ratio)") # Table 2 Vc/F 'Effect of other race' 1.15

    # Food effects on absorption (Supplementary Appendix Equations 6-8)
    e_fed_ka <- 0.726; label("Multiplicative effect of fed state on ka (ratio)") # Table 2 ka 'Effect of food' 0.726
    e_fed_missing_ka <- 0.663; label("Multiplicative effect of without-regard-to-food dosing on ka (ratio)") # Table 2 ka 'Effect of without regard to food' 0.663
    e_fed_fdepot <- 0.0683; label("Fractional decrease in F1 with food (fraction)") # Table 2 F1 'Effect of food' 0.0683
    e_fed_missing_fdepot <- 0.0809; label("Fractional decrease in F1 with without-regard-to-food dosing (fraction)") # Table 2 F1 'Effect of without regard to food' 0.0809

    # IIV on CL/F only (Discussion: reduced omega structure)
    etalcl ~ 0.102 # Table 2 final model 'omega2 (CL/F)' 0.102 (32% CV)

    # Log-scale additive residual error (Supplementary Appendix Equation 4),
    # reported on the SD scale: Results give 38.7% and 83.6%.
    expSdPh1 <- 0.387; label("Log-scale residual SD, phase 1 studies (unitless)") # Table 2 'Phase 1 residual error' 0.387
    expSdPh23 <- 0.836; label("Log-scale residual SD, phase 2/3 studies (unitless)") # Table 2 'Phase 2/3 residual error' 0.836
  })
  model({
    # eGFR above 120 mL/min/1.73 m^2 was set to 120 in the analysis data set
    # (Methods, Covariate Evaluation).
    egfr_capped <- min(CRCL, 120)

    cl <- exp(lcl + etalcl) * (WT / 85)^e_wt_cl_q * (egfr_capped / 90)^e_crcl_cl *
      e_dis_diab_cl^DIS_DIAB * e_sexf_cl^SEXF *
      e_race_black_cl^RACE_BLACK * e_race_asian_cl^RACE_ASIAN * e_race_other_cl^RACE_OTHER
    vc <- exp(lvc) * (WT / 85)^e_wt_vc_vp * (AGE / 65)^e_age_vc * e_sexf_vc^SEXF *
      e_race_black_vc^RACE_BLACK * e_race_asian_vc^RACE_ASIAN * e_race_other_vc^RACE_OTHER
    vp <- exp(lvp) * (WT / 85)^e_wt_vc_vp
    q <- exp(lq) * (WT / 85)^e_wt_cl_q
    ka <- exp(lka) * e_fed_ka^FED * e_fed_missing_ka^FED_MISSING
    tlag <- exp(ltlag)
    # Supplementary Appendix Equation 8: F1 = 1 - FOODEFF1*theta - FOODEFF2*theta
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
