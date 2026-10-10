Hayato_2022_lecanemab <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenous lecanemab in subjects with early Alzheimer's disease (and",
    "mild-to-moderate AD in the phase 1 studies), pooled from studies 101,",
    "104 and the phase 2 study 201 Core and open-label extension. Covariate",
    "effects: body weight (power, reference 73.7 kg) and serum albumin",
    "(power, reference 42.9 g/L) on CL; female sex and time-varying",
    "ADA-positive status as ratios on CL; body weight (power) and female",
    "sex on V1; Japanese race as a ratio on V2. Manufacturing Process B",
    "carries a relative bioavailability of 0.998 with between-subject",
    "variability that applies to Process B doses only. Proportional",
    "residual error differs by study. The individual PK parameters of this",
    "model drive the companion PK/PD models",
    "Hayato_2022_lecanemab_suvr, Hayato_2022_lecanemab_abeta4240 and",
    "Hayato_2022_lecanemab_ptau181.",
    sep = " "
  )
  reference <- paste(
    "Hayato S, Takenaka O, Sreerama Reddy SH, Landry I, Reyderman L,",
    "Koyama A, Swanson C, Yasuda S, Hussein Z (2022).",
    "Population pharmacokinetic-pharmacodynamic analyses of amyloid",
    "positron emission tomography and plasma biomarkers for lecanemab in",
    "subjects with early Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol 11(12):1578-1591.",
    "doi:10.1002/psp4.12862.",
    sep = " "
  )
  vignette <- "Hayato_2022_lecanemab"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    central = list(
      analyte = "lecanemab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lecanemab",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on CL (exponent 0.403) and V1 (exponent 0.606),",
        "normalised to 73.7 kg (supplement Model S1 $PK: '(WGT/73.7)'; the",
        "PK analysis set mean weight in Table S2). The Figure 2 caption's",
        "'reference 73.4 kg' subject is a typographical variant; the",
        "control stream and the equation printed under Table 1 both use",
        "73.7 kg.",
        sep = " "
      ),
      source_name = "WGT (Model S1 $INPUT); BW (Table 1 equations)"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL with a NEGATIVE exponent (-0.243): CL falls as",
        "albumin rises. Normalised to 42.9 g/L (Model S1 $PK and the Table",
        "1 equation; the PK analysis set mean in Table S2).",
        sep = " "
      ),
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Source column SEXN is coded 0 = male, 1 = female (Table 1",
        "footnote), so SEXF = SEXN with no inversion. Ratio 0.792 on CL and",
        "0.893 on V1, applied as ratio^SEXF.",
        sep = " "
      ),
      source_name = "SEXN"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = paste(
        "Time-varying categorical covariate (Methods: 'ADA status as a",
        "categorical time-variant covariate'). ADA-negative conclusive and",
        "ADA-negative inconclusive were pooled as negative, and PK",
        "observations with missing ADA status were assumed ADA-negative.",
        "Ratio 1.09 on CL, applied as ratio^ADA_POS.",
        sep = " "
      ),
      source_name = "ADA"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator (1 = Japanese, 0 = non-Japanese)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese)",
      notes = paste(
        "Model S1 $PK: 'RACE=0; IF (RACEN.EQ.3.1) RACE=1'. Ratio 0.455 on",
        "V2, applied as ratio^RACE_JAPANESE.",
        sep = " "
      ),
      source_name = "RACEN == 3.1 (raw); RACE (Model S1 $PK); JPN (Table 1 equations)"
    ),
    FORM_LEC_PROCESSB = list(
      description = "Lecanemab manufacturing Process B drug-product indicator (1 = Process B, 0 = Process A)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Process A; F fixed at 1)",
      notes = paste(
        "Per-dose indicator: Process B material was used in the study 201",
        "open-label extension (Table S1). Model S1 $PK codes 'F1=1; IF",
        "(FORM.EQ.1) F1=THETA(11)*EXP(ETA(4))', so the 0.998 relative",
        "bioavailability and its 34.2% between-subject variability apply to",
        "Process B doses only. PK analysis set: 8568 Process A and 459",
        "Process B observations (Table S2).",
        sep = " "
      ),
      source_name = "FORM"
    ),
    STUDY_LEC101 = list(
      description = "Lecanemab phase 1 study 101 cohort indicator (1 = study 101)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 201 when STUDY_LEC104 is also 0)",
      notes = paste(
        "Selects the study-specific proportional residual error (14.0%);",
        "Model S1 $ERROR uses ERR(1) for study 101. Set to 0 for",
        "simulation of the study 201 / clinical setting.",
        sep = " "
      ),
      source_name = "STUD"
    ),
    STUDY_LEC104 = list(
      description = "Lecanemab phase 1 study 104 cohort indicator (1 = study 104)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 201 when STUDY_LEC101 is also 0)",
      notes = paste(
        "Selects the study-specific proportional residual error (19.7%);",
        "Model S1 $ERROR: 'IF (STUD.EQ.104) Y = W + W*ERR(2)'.",
        sep = " "
      ),
      source_name = "STUD"
    )
  )

  covariatesDataExcluded <- list(
    ADA_TITER = list(
      description = "Anti-drug antibody titer (time-varying)",
      units = "(titer)",
      type = "continuous",
      notes = paste(
        "Tested on CL as a continuous time-varying covariate and not",
        "retained (Table S3 model 3, dOFV 2.99, p > 0.01); the binary",
        "ADA_POS status was retained instead.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 725L,
    n_studies = 3L,
    n_observations = paste(
      "9027 serum lecanemab concentrations: 661 (7.3%) from study 101,",
      "371 (4.1%) from study 104, 7995 (88.6%) from study 201 Core and OLE",
      sep = " "
    ),
    age_range = "50-93 years (mean 71.0, SD 8.2; Table S2)",
    weight_range = "41.2-124.7 kg (mean 73.7, SD 14.7; Table S2)",
    albumin_range = "35.0-53.0 g/L (mean 42.9, SD 2.9; Table S2)",
    sex_female_pct = 46.9,
    race_ethnicity = c(
      White = 86.3,
      `Black/African American` = 3.6,
      `Asian/Other Asian (excluding Chinese and Japanese)` = 2.1,
      Japanese = 6.9,
      `American Indian/Alaskan/other/missing` = 1.1
    ),
    disease_state = paste(
      "Early Alzheimer's disease (mild cognitive impairment due to AD or",
      "mild AD dementia) in studies 104 and 201; mild-to-moderate AD in",
      "study 101",
      sep = " "
    ),
    dose_range = paste(
      "Intravenous infusion over 60 +/- 10 min: single doses 0.1-15 mg/kg",
      "(study 101 SAD); 0.3-3 mg/kg every 4 weeks and 10 mg/kg every 2",
      "weeks (study 101 MAD); 2.5, 5 and 10 mg/kg every 4 weeks (study",
      "104); 2.5, 5 and 10 mg/kg biweekly or 5 and 10 mg/kg monthly for 18",
      "months (study 201 Core); 10 mg/kg biweekly (study 201 OLE)",
      sep = " "
    ),
    regions = "Multinational (North America, Europe, Japan, South Korea)",
    ada_status = "observations: 7208 (79.9%) ADA-negative, 1247 (13.8%) ADA-positive, 572 (6.3%) missing (assumed negative) (Table S2)",
    notes = paste(
      "Baseline demographics from supplement Table S2 (PK column). The",
      "structural model and its sex, weight, albumin and Japanese-race",
      "effects were carried from the previously published lecanemab PK",
      "model and re-estimated on this dataset with the study 201 OLE data",
      "added; ADA status and the Process B relative bioavailability were",
      "added in this analysis (Table S3). Eta shrinkage: CL 9.96%, V1",
      "30.5%, V2 31.7%, F 63.2% (Table 1 footnote).",
      sep = " "
    )
  )

  ini({
    # Final estimates from Hayato 2022 Table 1; functional forms from the
    # equations printed under Table 1 and supplement Model S1 $PK.
    lcl <- log(0.0181); label("Clearance for the reference subject (L/h)") # Table 1 CL = 0.0181 L/h (%RSE 2.55)
    lvc <- log(3.22); label("Central volume of distribution for the reference subject (L)") # Table 1 V1 = 3.22 L (%RSE 1.18)
    lq <- log(0.0349); label("Intercompartmental clearance (L/h)") # Table 1 Q = 0.0349 L/h (%RSE 8.02)
    lvp <- log(2.19); label("Peripheral volume of distribution for the reference subject (L)") # Table 1 V2 = 2.19 L (%RSE 7.21)

    # Process A bioavailability is the structural anchor: Model S1 $PK
    # sets 'F1=1' with no THETA.
    lfcentral <- fixed(log(1)); label("Relative bioavailability of the intravenous dose for Process A (unitless)") # Model S1 $PK 'F1=1'
    e_processb_f <- 0.998; label("Relative bioavailability for Process B vs Process A (unitless)") # Table 1 'F for process B' = 0.998 (%RSE 4.07)

    e_wt_cl <- 0.403; label("Power exponent on (WT/73.7 kg) for CL (unitless)") # Table 1 'Weight ~ CL (exponent)' = 0.403
    e_alb_cl <- -0.243; label("Power exponent on (ALB/42.9 g/L) for CL (unitless)") # Table 1 'Albumin ~ CL (exponent)' = -0.243
    e_female_cl <- 0.792; label("CL ratio for females vs males, applied as ratio^SEXF (unitless)") # Table 1 'Females ~ CL (ratio)' = 0.792
    e_ada_cl <- 1.09; label("CL ratio for ADA-positive vs ADA-negative, applied as ratio^ADA_POS (unitless)") # Table 1 'ADA positive ~ CL (ratio to ADA negative)' = 1.09
    e_wt_vc <- 0.606; label("Power exponent on (WT/73.7 kg) for V1 (unitless)") # Table 1 'Weight ~ V1 (exponent)' = 0.606
    e_female_vc <- 0.893; label("V1 ratio for females vs males, applied as ratio^SEXF (unitless)") # Table 1 'Females ~ V1 (ratio)' = 0.893
    e_japanese_vp <- 0.455; label("V2 ratio for Japanese vs non-Japanese, applied as ratio^RACE_JAPANESE (unitless)") # Table 1 'Japanese race ~ V2 (ratio)' = 0.455

    # Table 1 footnote defines CV% as 'square root of variance x 100', so
    # omega^2 = (CV%/100)^2. Model S1 $OMEGA is diagonal.
    etalcl ~ 0.151321 # Table 1 IIV CL 38.9% -> 0.389^2
    etalvc ~ 0.0196 # Table 1 IIV V1 14.0% -> 0.140^2
    etalvp ~ 0.990025 # Table 1 IIV V2 99.5% -> 0.995^2
    etalfcentral ~ 0.116964 # Table 1 IIV F 34.2% -> 0.342^2; applies to Process B doses only (Model S1 $PK)

    # Model S1 $ERROR: W = F + 0.01; Y = W + W*ERR(k), with k chosen by
    # study. Table 1 'Residual variability (CV%)' is sqrt(sigma^2) x 100.
    propSd <- 0.303; label("Proportional residual error, study 201 (fraction)") # Table 1 'Proportional: study 201' = 30.3%
    propSdStudy101 <- 0.140; label("Proportional residual error, study 101 (fraction)") # Table 1 'Proportional: study 101' = 14.0%
    propSdStudy104 <- 0.197; label("Proportional residual error, study 104 (fraction)") # Table 1 'Proportional: study 104' = 19.7%
  })

  model({
    # Individual parameters (Table 1 equations; Model S1 $PK)
    cl <- exp(lcl + etalcl) * (WT / 73.7)^e_wt_cl * (ALB / 42.9)^e_alb_cl *
      e_female_cl^SEXF * e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 73.7)^e_wt_vc * e_female_vc^SEXF
    vp <- exp(lvp + etalvp) * e_japanese_vp^RACE_JAPANESE
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Model S1 $PK: 'F1=1; IF (FORM.EQ.1) F1=THETA(11)*EXP(ETA(4))' --
    # Process A has F = 1 exactly with no variability; the ratio and its
    # eta apply to Process B doses only.
    fprocessb <- exp(lfcentral + etalfcentral) * e_processb_f
    f(central) <- (1 - FORM_LEC_PROCESSB) + FORM_LEC_PROCESSB * fprocessb

    # Dose (mg) / volume (L) = mg/L = ug/mL
    Cc <- central / vc

    # Model S1 $ERROR: W = F + 0.01; Y = W + W*ERR(k). The residual SD is
    # sigma_k * (Cc + 0.01), i.e. a combined1 error whose additive part is
    # 0.01 ug/mL times the study's proportional SD.
    propSdCc <- propSd * (1 - STUDY_LEC101 - STUDY_LEC104) +
      propSdStudy101 * STUDY_LEC101 + propSdStudy104 * STUDY_LEC104
    addSdCc <- 0.01 * propSdCc
    Cc ~ add(addSdCc) + prop(propSdCc) + combined1()
  })
}
