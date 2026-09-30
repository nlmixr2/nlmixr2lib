Kuchimanchi_2022_diroximelFumarate <- function() {
  description <- paste(
    "Joint population PK model for the two diroximel fumarate (DRF)",
    "metabolites, monomethyl fumarate (MMF, the active moiety) and",
    "2-hydroxyethyl succinimide (HES, inactive), after oral DRF in 341",
    "healthy volunteers and 48 patients with relapsing-remitting multiple",
    "sclerosis across 11 phase I and III studies. DRF itself is not",
    "measurable in plasma, so each dose enters two parallel absorption",
    "chains as its molar equivalent of MMF and of HES; each chain is a",
    "dose compartment followed by eight transit compartments (nine",
    "first-order transfers at the metabolite's absorption rate constant)",
    "into a one-compartment disposition with first-order elimination. One",
    "central volume is shared by both metabolites; HES bioavailability is",
    "fixed at 0.6 from a mass-balance study and MMF bioavailability is",
    "estimated relative to it. Body weight scales both clearances and the",
    "volume, baseline eGFR scales HES clearance, and patients with MS have",
    "lower clearance of both metabolites. Meal fat content and evening",
    "dosing slow absorption, meal fat content lowers MMF bioavailability,",
    "and HES absorption carries a lag after an evening dose or a low-fat",
    "meal. The log-scale residual error is stratified by metabolite, meal",
    "state and dose time."
  )
  reference <- paste(
    "Kuchimanchi M, Bockbrader H, Dolphin N, Epling D, Quinlan L, Chapel S,",
    "Penner N. Development of a Population Pharmacokinetic Model for the",
    "Diroximel Fumarate Metabolites Monomethyl Fumarate and 2-Hydroxyethyl",
    "Succinimide Following Oral Administration of Diroximel Fumarate in",
    "Healthy Participants and Patients with Multiple Sclerosis. Neurol Ther.",
    "2022;11(1):353-371. doi:10.1007/s40120-021-00316-6"
  )
  vignette <- "Kuchimanchi_2022_diroximelFumarate"
  # Doses are milligrams of DIROXIMEL FUMARATE; model() converts each dose to
  # its molar-equivalent mass of MMF and of HES. Both outputs (Cc = MMF,
  # Cc_hes = HES) are in ug/mL.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # The residual SD is stratified by meal state and dose time within each
  # metabolite (eight thetas in the ESM control stream $ERROR block); the
  # canonical matcher recognises only the bare expSd / expSd_hes names.
  paper_specific_residual_sds <- c(
    "expSdFed",
    "expSdEvening",
    "expSdFedMissing",
    "expSdFed_hes",
    "expSdEvening_hes",
    "expSdFedMissing_hes"
  )

  compartmentData <- list(
    depot = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "monomethyl fumarate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "monomethyl fumarate", units = "mg", specimen = "plasma", verified = TRUE),
    depot_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit6_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit7_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit8_hes = list(
      analyte = "2-hydroxyethyl succinimide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_hes = list(analyte = "2-hydroxyethyl succinimide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed) body weight, source column BWT (ESM control",
        "stream $INPUT; the stream imputes the 78 kg median when missing).",
        "Power scaling about the 78 kg cohort median on MMF clearance",
        "(exponent 0.831), HES clearance (0.335) and the shared central",
        "volume (0.878)."
      ),
      source_name = "BWT"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (MDRD, BSA-denormalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline eGFR from the Modification of Diet in Renal Disease",
        "equation, 'expressed in absolute units (mL/min) following",
        "denormalization using individual participant body [surface] area'",
        "(ESM Table S2 footnote b). This is NOT the register's default",
        "BSA-normalized mL/min/1.73 m^2 unit: multiply an MDRD eGFR by",
        "BSA / 1.73 before supplying it. Source column BEGFR; the stream",
        "imputes the 111.9 mL/min median when missing. Enters HES clearance",
        "only, as (CRCL / 111.9)^0.547."
      ),
      source_name = "BEGFR"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy volunteer, 0 = patient with relapsing-remitting MS)",
      units = "(binary)",
      type = "binary",
      reference_category = 1,
      notes = paste(
        "The source column PTST is the complement (PTST = 1 for a patient",
        "with MS from EVOLVE-MS-1 / EVOLVE-MS-2, 0 for a healthy or",
        "renally impaired phase I participant), so the model uses",
        "PTST = 1 - DIS_HEALTHY. The healthy volunteer is therefore the",
        "reference: patients have 28.4% lower MMF clearance and 12.2% lower",
        "HES clearance. Participants with renal impairment in study A108",
        "were not patients with MS and take DIS_HEALTHY = 1. The ESM notes",
        "that the patient-status covariate is confounded with the",
        "unknown-meal-state stratum (FED_MISSING), which applies only to",
        "patients."
      ),
      source_name = "PTST"
    ),
    FED_LOWFAT = list(
      description = "Dose taken with a low-fat meal (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column FAT = 1 (ESM stream BFAT1). Study A109: low-fat,",
        "low-calorie meal consumed 30 min before the dose (Table 1).",
        "Mutually exclusive with FED_MEDFAT, FED_HIGHFAT and FED_MISSING;",
        "all four are 0 for a fasted dose, which is the reference. Slows",
        "MMF and HES absorption, lowers MMF bioavailability and adds a",
        "0.421 h HES absorption lag."
      ),
      source_name = "FAT"
    ),
    FED_MEDFAT = list(
      description = "Dose taken with a medium-fat meal (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column FAT = 2 (ESM stream BFAT2; 'Moderate' in ESM Table",
        "S2). Study A109: medium-fat, medium-calorie meal consumed 30 min",
        "before the dose (Table 1). Mutually exclusive with the other meal",
        "indicators."
      ),
      source_name = "FAT"
    ),
    FED_HIGHFAT = list(
      description = "Dose taken with a high-fat meal (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column FAT = 3 (ESM stream BFAT3). Studies A102 part A and",
        "A104: high-fat, high-calorie meal consumed 30 min before the dose",
        "(Table 1). Mutually exclusive with the other meal indicators."
      ),
      source_name = "FAT"
    ),
    FED_MISSING = list(
      description = "Meal state at dosing unknown (1 = unknown, 0 = recorded)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column FAT = 4 (ESM stream BFAT4; Table 2 'UNK'",
        "administration with or without meal of unknown fat content, only in",
        "patients). The phase III patients were told to take DRF with or",
        "without food while avoiding a high-fat meal, so the stratum is a",
        "nuisance average over unrecorded meal states. Its estimated effects",
        "(faster absorption of both metabolites and its own residual SDs)",
        "belong to the published fit; set it to 0 when simulating a defined",
        "meal state."
      ),
      source_name = "FAT"
    ),
    DOSETIME_EVENING = list(
      description = "Evening (PM) dose indicator (1 = evening dose, 0 = morning dose)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column PM, used in the stream as PM1 after zeroing it for",
        "studies coded 110 and above ('PM IS ONLY FOR STUD 102'): only the",
        "serial evening-dose profiles of study A102 part B informed the",
        "evening effects. Slows MMF and HES absorption, adds a 1.96 h HES",
        "absorption lag and selects the evening residual SDs. The model",
        "reads the covariate at every record (as NONMEM does), so hold it",
        "constant across a dosing interval when simulating."
      ),
      source_name = "PM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 389L,
    n_studies = 11L,
    age_range = "18-75 years",
    age_median = "35 years",
    weight_range = "47.4-126.3 kg",
    weight_median = "78 kg",
    sex_female_pct = 49.4,
    race_ethnicity = c(White = 66.3, Black = 30.8, Asian = 1.5, Other = 1.2),
    hispanic_pct = 30.5,
    disease_state = paste(
      "341 healthy volunteers (including 32 participants with normal to",
      "severely impaired renal function in study A108) and 48 patients",
      "with relapsing-remitting multiple sclerosis (EVOLVE-MS-1 and",
      "EVOLVE-MS-2)"
    ),
    renal_function = paste(
      "eGFR category at baseline: normal 75.5%, mild impairment 20.0%,",
      "moderate 2.3%, severe 2.0% (MDRD eGFR, ESM Table S2); median",
      "eGFR 111.4 mL/min (range 15.2-185.8, ESM Table S1)"
    ),
    dose_range = paste(
      "DRF 49-980 mg single doses and 210-924 mg twice daily; 69% received",
      "the approved 462 mg dose. Fasted (n = 252) or with low-fat (n = 47),",
      "medium-fat (n = 47) or high-fat (n = 58) meals; patients with an",
      "unknown meal state (n = 48)"
    ),
    n_observations = "4694 MMF and 8088 HES concentrations",
    notes = paste(
      "Demographics from the Results 'Study Population' paragraph and ESM",
      "Tables S1-S2 (192 of 389 female, 49.4%). Studies 001, A102-A106,",
      "A108-A110 (phase I) and EVOLVE-MS-1 / EVOLVE-MS-2 (phase III), Table 1."
    )
  )

  ini({
    # Structural parameters (Table 2; ESM control stream $THETA order). The
    # typical values refer to a 78 kg healthy volunteer with an eGFR of
    # 111.9 mL/min, dosed fasted in the morning.
    lcl <- log(13.5); label("MMF clearance CLMMF (L/h)") # Table 2 theta 1 CLMMF = 13.5 L/h
    lvc <- log(30.4); label("Central volume shared by MMF and HES, Vc (L)") # Table 2 theta 2 Vc = 30.4 L (stream V3 = V2)
    lka <- log(5.04); label("MMF absorption / transit rate constant KaMMF (1/h)") # Table 2 theta 3 KaMMF = 5.04 1/h
    lka_hes <- log(3.24); label("HES absorption / transit rate constant KaHES (1/h)") # Table 2 theta 6 KaHES = 3.24 1/h
    lcl_hes <- log(1.49); label("HES clearance CLHES (L/h)") # Table 2 theta 7 CLHES = 1.49 L/h
    lfdepot_hes <- fixed(log(0.6)); label("HES bioavailability F4, from the mass-balance study (fraction)") # Table 2 theta 8 F4 = 0.6 FIXED
    lfdepot <- log(0.162); label("MMF bioavailability F1 relative to the HES-scaled volume (fraction)") # Table 2 theta 9 F1 = 0.162

    # Continuous covariate effects (power functions)
    e_wt_vc <- 0.878; label("Power exponent of (WT/78) on Vc (unitless)") # Table 2 theta 10 'WT on Vc' = 0.878
    e_crcl_cl_hes <- 0.547; label("Power exponent of (eGFR/111.9) on CLHES (unitless)") # Table 2 theta 26 'eGFR on CLHES' = 0.547
    e_wt_cl <- 0.831; label("Power exponent of (WT/78) on CLMMF (unitless)") # Table 2 theta 27 = 0.831 (row mislabelled 'eGFR on CLHES'; the stream comment and model equation give WT on CLM)
    e_wt_cl_hes <- 0.335; label("Power exponent of (WT/78) on CLHES (unitless)") # Table 2 theta 28 'WT on CLHES' = 0.335

    # Patient status (PTST = 1 - DIS_HEALTHY): fractional change, (1 + theta * PTST)
    e_patient_cl <- -0.284; label("Fractional change in CLMMF for patients with MS vs healthy volunteers (unitless)") # Table 2 theta 35 'PTST on CLMMF' = -0.284
    e_patient_cl_hes <- -0.122; label("Fractional change in CLHES for patients with MS vs healthy volunteers (unitless)") # Table 2 theta 36 'PTST on CLHES' = -0.122

    # Meal and dose-time effects on KaMMF: fractional change, (1 + theta * indicator).
    # All except the unknown-meal effect were fixed at the phase I base-model
    # estimates (Table 2 footnote; ESM Table S5).
    e_dosetime_evening_ka <- fixed(-0.592); label("Fractional change in KaMMF for an evening dose (unitless)") # Table 2 theta 11 'PM dosing on KaMMF' = -0.592 FIXED
    e_fed_lowfat_ka <- fixed(-0.368); label("Fractional change in KaMMF with a low-fat meal (unitless)") # Table 2 theta 12 'LOW on KaMMF' = -0.368 FIXED
    e_fed_medfat_ka <- fixed(-0.512); label("Fractional change in KaMMF with a medium-fat meal (unitless)") # Table 2 theta 13 'MED on KaMMF' = -0.512 FIXED
    e_fed_highfat_ka <- fixed(-0.666); label("Fractional change in KaMMF with a high-fat meal (unitless)") # Table 2 theta 14 'HI on KaMMF' = -0.666 FIXED
    e_fed_missing_ka <- 0.843; label("Fractional change in KaMMF with an unknown meal state (unitless)") # Table 2 theta 15 'UNK on KaMMF' = 0.843

    # Meal effects on MMF bioavailability F1
    e_fed_lowfat_fdepot <- fixed(-0.296); label("Fractional change in F1 with a low-fat meal (unitless)") # Table 2 theta 16 'LOW on F1' = -0.296 FIXED
    e_fed_medfat_fdepot <- fixed(-0.301); label("Fractional change in F1 with a medium-fat meal (unitless)") # Table 2 theta 17 = -0.301 FIXED (row mislabelled 'LOW on F1'; stream comment 'MED FAT ON F1')
    e_fed_highfat_fdepot <- fixed(-0.131); label("Fractional change in F1 with a high-fat meal (unitless)") # Table 2 theta 18 'HI on F1' = -0.131 FIXED

    # Meal and dose-time effects on KaHES
    e_dosetime_evening_ka_hes <- fixed(-0.267); label("Fractional change in KaHES for an evening dose (unitless)") # Table 2 theta 19 'PM dosing on KaHES' = -0.267 FIXED
    e_fed_lowfat_ka_hes <- fixed(-0.335); label("Fractional change in KaHES with a low-fat meal (unitless)") # Table 2 theta 20 'LOW on KaHES' = -0.335 FIXED
    e_fed_medfat_ka_hes <- fixed(-0.492); label("Fractional change in KaHES with a medium-fat meal (unitless)") # Table 2 theta 21 'MED on KaHES' = -0.492 FIXED
    e_fed_highfat_ka_hes <- fixed(-0.621); label("Fractional change in KaHES with a high-fat meal (unitless)") # Table 2 theta 22 'HI on KaHES' = -0.621 FIXED
    e_fed_missing_ka_hes <- 0.399; label("Fractional change in KaHES with an unknown meal state (unitless)") # Table 2 theta 23 'UNK on KaHES' = 0.399

    # HES absorption lag: additive contributions on a zero baseline,
    # ALAG4 = 0 + PM1 * theta24 + BFAT1 * theta25 (ESM control stream $PK)
    e_dosetime_evening_tlag_hes <- fixed(1.96); label("HES absorption lag time for an evening dose (h)") # Table 2 theta 24 'HES ALAG4 PM dosing' = 1.96 h FIXED
    e_fed_lowfat_tlag_hes <- fixed(0.421); label("HES absorption lag time with a low-fat meal (h)") # Table 2 theta 25 'HES ALAG4 LOW' = 0.421 h FIXED

    # IIV: log-normal, diagonal OMEGA. The ESM (Random Effects Model
    # Development) reports %CV = sqrt(omega^2) * 100, so omega^2 = (CV/100)^2.
    etalcl ~ 0.056169 # Table 2 ETA1 'CLMMF' 23.7 %CV -> 0.237^2
    etalvc ~ 0.039204 # Table 2 ETA2 'Vc' 19.8 %CV -> 0.198^2
    etalcl_hes ~ 0.0324 # Table 2 ETA4 'CLHES' 18.0 %CV -> 0.180^2
    etalka ~ 0.1369 # Table 2 IIV row 5 = 37.0 %CV -> 0.370^2 (row mislabelled 'ETA4-CLHES'; stream ETA5 is KaMMF, Results text '37% for KaMMF')
    etalka_hes ~ 0.179776 # Table 2 ETA8 'KaHES' 42.4 %CV -> 0.424^2

    # Residual error: log-transform-both-sides, ln(Y) = ln(C) + W * EPS with
    # $SIGMA 1 FIX, so each theta is the SD on the log scale. Table 2 prints
    # the thetas x 100 as '%'.
    expSd <- 0.895; label("MMF log-scale residual SD, morning dose fasted") # Table 2 theta 4 'RE MMF (AM dose, fasted)' = 89.5%
    expSd_hes <- 0.252; label("HES log-scale residual SD, morning dose fasted") # Table 2 theta 5 'RE HES (AM dose, fasted)' = 25.2%
    expSdFed <- 1.03; label("MMF log-scale residual SD, morning dose with a low-, medium- or high-fat meal") # Table 2 theta 29 'RE MMF (AM dose, fed)' = 103%
    expSdEvening <- 1.12; label("MMF log-scale residual SD, evening dose fasted") # Table 2 theta 30 'RE MMF (PM dose, fasted)' = 112%
    expSdFedMissing <- 1.02; label("MMF log-scale residual SD, unknown meal state") # Table 2 theta 31 'RE MMF (UNK)' = 102%
    expSdFed_hes <- 0.468; label("HES log-scale residual SD, morning dose with a low-, medium- or high-fat meal") # Table 2 theta 32 'RE HES (AM dose, fed)' = 46.8%
    expSdEvening_hes <- 0.184; label("HES log-scale residual SD, evening dose fasted") # Table 2 theta 33 'RE HES (PM dose, fasted)' = 18.4%
    expSdFedMissing_hes <- 0.372; label("HES log-scale residual SD, unknown meal state") # Table 2 theta 34 'RE HES (UNK)' = 37.2%
  })

  model({
    # Molecular weights (g/mol) from the molecular formulae and IUPAC standard
    # atomic weights: DRF C11H13NO6, MMF C5H6O4, HES C6H9NO3. The paper enters
    # each DRF dose as its molar equivalent of MMF and of HES ('dose of MMF or
    # HES (molar) = dose of DRF (g)/molecular weight (DRF)', Model
    # Development), with concentrations in mass units (stream S2 = V2/1000).
    mw_drf <- 255.22
    mw_mmf <- 130.10
    mw_hes <- 143.14

    # Derived covariate terms
    patient <- 1 - DIS_HEALTHY
    fed_any <- FED_LOWFAT + FED_MEDFAT + FED_HIGHFAT

    # Individual parameters (ESM control stream $PK)
    cl <- exp(lcl + etalcl) * (WT / 78)^e_wt_cl * (1 + e_patient_cl * patient)
    cl_hes <- exp(lcl_hes + etalcl_hes) * (CRCL / 111.9)^e_crcl_cl_hes *
      (WT / 78)^e_wt_cl_hes * (1 + e_patient_cl_hes * patient)
    vc <- exp(lvc + etalvc) * (WT / 78)^e_wt_vc

    ka <- exp(lka + etalka) *
      (1 + e_dosetime_evening_ka * DOSETIME_EVENING) *
      (1 + e_fed_lowfat_ka * FED_LOWFAT) *
      (1 + e_fed_medfat_ka * FED_MEDFAT) *
      (1 + e_fed_highfat_ka * FED_HIGHFAT) *
      (1 + e_fed_missing_ka * FED_MISSING)
    ka_hes <- exp(lka_hes + etalka_hes) *
      (1 + e_dosetime_evening_ka_hes * DOSETIME_EVENING) *
      (1 + e_fed_lowfat_ka_hes * FED_LOWFAT) *
      (1 + e_fed_medfat_ka_hes * FED_MEDFAT) *
      (1 + e_fed_highfat_ka_hes * FED_HIGHFAT) *
      (1 + e_fed_missing_ka_hes * FED_MISSING)

    fdepot <- exp(lfdepot) *
      (1 + e_fed_lowfat_fdepot * FED_LOWFAT) *
      (1 + e_fed_medfat_fdepot * FED_MEDFAT) *
      (1 + e_fed_highfat_fdepot * FED_HIGHFAT)
    fdepot_hes <- exp(lfdepot_hes)
    tlag_hes <- e_dosetime_evening_tlag_hes * DOSETIME_EVENING +
      e_fed_lowfat_tlag_hes * FED_LOWFAT

    kel <- cl / vc
    kel_hes <- cl_hes / vc

    # MMF: dose compartment -> 8 transit compartments -> plasma (stream
    # compartments 1 -> 5 -> 7 -> ... -> 19 -> 2, nine transfers all at KAM)
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(transit4) <- ka * transit3 - ka * transit4
    d/dt(transit5) <- ka * transit4 - ka * transit5
    d/dt(transit6) <- ka * transit5 - ka * transit6
    d/dt(transit7) <- ka * transit6 - ka * transit7
    d/dt(transit8) <- ka * transit7 - ka * transit8
    d/dt(central) <- ka * transit8 - kel * central

    # HES: dose compartment -> 8 transit compartments -> plasma (stream
    # compartments 4 -> 6 -> 8 -> ... -> 20 -> 3, nine transfers all at KAH)
    d/dt(depot_hes) <- -ka_hes * depot_hes
    d/dt(transit1_hes) <- ka_hes * depot_hes - ka_hes * transit1_hes
    d/dt(transit2_hes) <- ka_hes * transit1_hes - ka_hes * transit2_hes
    d/dt(transit3_hes) <- ka_hes * transit2_hes - ka_hes * transit3_hes
    d/dt(transit4_hes) <- ka_hes * transit3_hes - ka_hes * transit4_hes
    d/dt(transit5_hes) <- ka_hes * transit4_hes - ka_hes * transit5_hes
    d/dt(transit6_hes) <- ka_hes * transit5_hes - ka_hes * transit6_hes
    d/dt(transit7_hes) <- ka_hes * transit6_hes - ka_hes * transit7_hes
    d/dt(transit8_hes) <- ka_hes * transit7_hes - ka_hes * transit8_hes
    d/dt(central_hes) <- ka_hes * transit8_hes - kel_hes * central_hes

    # Each DRF dose (mg) must be given twice, once to depot and once to
    # depot_hes; bioavailability carries the molar conversion to metabolite mass.
    f(depot) <- fdepot * mw_mmf / mw_drf
    f(depot_hes) <- fdepot_hes * mw_hes / mw_drf
    alag(depot_hes) <- tlag_hes

    # Observations: mg / L = ug/mL
    Cc <- central / vc
    Cc_hes <- central_hes / vc

    # Residual SD by condition, in the order of the stream's $ERROR IF-chain:
    # unknown meal state overrides dose time, dose time overrides meal. The
    # stream leaves an evening dose WITH a known meal undefined (no such
    # records); it takes the evening SD here.
    if (FED_MISSING == 1) {
      sdMmf <- expSdFedMissing
      sdHes <- expSdFedMissing_hes
    } else if (DOSETIME_EVENING == 1) {
      sdMmf <- expSdEvening
      sdHes <- expSdEvening_hes
    } else if (fed_any == 1) {
      sdMmf <- expSdFed
      sdHes <- expSdFed_hes
    } else {
      sdMmf <- expSd
      sdHes <- expSd_hes
    }

    Cc ~ lnorm(sdMmf)
    Cc_hes ~ lnorm(sdHes)
  })
}
