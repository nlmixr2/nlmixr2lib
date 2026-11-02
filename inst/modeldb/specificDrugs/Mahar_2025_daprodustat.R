Mahar_2025_daprodustat <- function() {
  description <- paste(
    "Longitudinal dose-hemoglobin model of daprodustat in adults with anemia",
    "of chronic kidney disease, fitted to five global phase III studies",
    "(Mahar 2025, Dose-Hgb Model 3). The current daprodustat dose level,",
    "allometrically scaled by baseline body weight, stimulates the production",
    "rate of a precursor-cell compartment through a power function",
    "(1 + SLOPE * dose^GAMMA, with GAMMA proportional to SLOPE). Precursors",
    "(mean transit t 3.19 days, fixed) mature into twelve red-cell age",
    "compartments whose transit times sum to the red-cell lifespan (82.8",
    "days, fixed, with covariate effects), and total hemoglobin is the sum of",
    "the twelve age compartments. Hemoglobin production declines",
    "exponentially toward a (possibly negative) floor with a disease",
    "progression rate that depends on dialysis status and prior ESA use and",
    "has year-specific random effects for the first, second and later years.",
    "Prior ESA users start at an ESA-elevated apparent baseline and wash out",
    "to a separately estimated true baseline. Patient-group-specific",
    "baselines (non-dialysis, hemodialysis, peritoneal dialysis, starting",
    "dialysis) with CKD stage, body weight, baseline hepcidin, race, prior",
    "ESA dose and ESA hyporesponder effects; study-specific additive",
    "residual error on hemoglobin (g/L)."
  )
  reference <- paste(
    "Mahar KM, Yang S, Visser S, van Noort M, van den Berg P, Post TM.",
    "Integrated Longitudinal Population Dose-Hemoglobin Response of",
    "Daprodustat Following Dose Titration in Patients With Anemia in Chronic",
    "Kidney Disease. Clin Pharmacol Ther. 2025;117(3):846-855.",
    "doi:10.1002/cpt.3544"
  )
  vignette <- "Mahar_2025_daprodustat"
  units <- list(time = "day", dosing = "mg", concentration = "g/L")

  covariateData <- list(
    DOSE_DAPRODUSTAT_MG = list(
      description = "Current daprodustat dose level per administration (0 when no daprodustat is being given)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. The model is driven by the dose level, not by",
        "concentrations: the control stream data item CDOSE (`DDR = CDOSE /",
        "WTSL * CLOP_DOSE` in $DES) holds the dose currently assigned by the",
        "titration algorithm and is carried forward between titration steps.",
        "For once-daily regimens this is the daily dose (1-24 mg in the",
        "studies); for the three-times-weekly regimen of study 204837 it is",
        "the per-administration dose (2-48 mg) and the effect is multiplied by",
        "the study-204837 factor e_study_204837_slope (Mahar 2025 Table 2,",
        "0.422, close to 3/7). No dosing events are needed."
      ),
      source_name = "CDOSE"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value (WTBL; the control stream imputes the",
        "median 73.3 kg when missing). Reference 73.3 kg for all three uses:",
        "power effect on both hemoglobin baselines, power effect on red-cell",
        "lifespan, and allometric scaling of the dose with exponent 0.75",
        "(fixed)."
      ),
      source_name = "WTBL"
    ),
    CKD_STAGE = list(
      description = "Chronic kidney disease stage at baseline (1-5 by CKD-EPI eGFR for non-dialysis subjects; 6 for dialysis-dependent subjects)",
      units = "(stage)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Mahar 2025 Table S2 defines stages 1-5 by CKD-EPI eGFR and adds",
        "stage 6 for dialysis-dependent subjects (HD or PD) so that stage-5",
        "non-dialysis and dialysis subjects can be told apart. Enters as a",
        "power effect `(CKD_STAGE / 5)^e_ckd_stage_rbase` on both the",
        "apparent baseline of ESA non-users and the true baseline of ESA users",
        "(control stream `BLSTAG = (STAGE/medStage)**THETA(11)`, medStage =",
        "5)."
      ),
      source_name = "STAGE"
    ),
    HEPCIDIN_BL = list(
      description = "Baseline serum hepcidin concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Power effect `(HEPCIDIN_BL / 121.356)^e_hepcidin_bl_rbase`",
        "on the same baselines as CKD_STAGE (control stream `BLHEP`, medBhep =",
        "121.356 ng/mL; the median is imputed when missing)."
      ),
      source_name = "BHEP"
    ),
    ESAD = list(
      description = "Prior erythropoiesis-stimulating-agent (ESA) dose at baseline; ESAD > 0 identifies a prior ESA user",
      units = "IU/week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Computed as PRIEPDB (IU/kg/week) x WTBL (kg) (control stream",
        "`PEPOD2 = PRIEPDB * WTBL`). Enters the true baseline of ESA users as",
        "the uncentred exponential factor `exp(e_esad_rbase_true / 1000 *",
        "ESAD)`; Table S8 reports -1.2% to -14.0% over the 5th-95th",
        "percentiles 1521-18994 IU/week. ESAD > 0 also serves as the ESA-user",
        "indicator (control stream EPO = 1 historic or 2 recent; both were",
        "modelled identically), selecting the apparent-to-true baseline",
        "washout and the ESA-user effects on red-cell lifespan and disease",
        "progression. For an ESA user with unrecorded prior dose the control",
        "stream imputes PRIEPDB = 71.299 IU/kg/week, i.e. ESAD = 71.299 x WT.",
        "Set ESAD = 0 for ESA non-users; the model also forces the ESA-user",
        "indicator to 0 in studies 205270 and 201410, as the control stream",
        "does."
      ),
      source_name = "PRIEPDB x WTBL (EPO > 0 as the user flag)"
    ),
    ESA_HYPORESPONDER = list(
      description = "ESA hyporesponder indicator (1 = prior ESA user classified as hyporesponder)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not a hyporesponder, or ESA non-user / unknown)",
      notes = paste(
        "Hyporesponder per Mahar 2025 Table S4 footnote c: ESA resistance",
        "index >= 1.5 IU/kg/week/(g/L) or prior ESA dose >= 450 IU/kg/week.",
        "Lowers the apparent (ESA-elevated) baseline of ESA users by",
        "e_esa_hyporesponder_rbase = -3.07%. The control stream codes unknown",
        "status as 222, which carries no effect; encode it as 0 here."
      ),
      source_name = "HYPORESP"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Hemodialysis indicator (1 = on hemodialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-dialysis when PERIT_DIAL is also 0)",
      notes = paste(
        "Three-level dialysis population with non-dialysis (ND) as the",
        "reference, together with PERIT_DIAL (control stream POP 1 = ND,",
        "2 = HD, 3 = PD). Chronic HD (STUDY_201410 = 0) selects the HD",
        "baselines, +13.8% red-cell lifespan and -57.3% disease progression.",
        "In study 201410 (starting dialysis) it only selects the",
        "starting-HD apparent baseline; the starting-dialysis effects then",
        "replace the chronic-dialysis ones."
      ),
      source_name = "POP == 2"
    ),
    PERIT_DIAL = list(
      description = "Peritoneal-dialysis indicator (1 = on peritoneal dialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-dialysis when RRT_HEMODIAL_STATUS is also 0)",
      notes = paste(
        "Paired with RRT_HEMODIAL_STATUS against a non-dialysis reference.",
        "Chronic PD selects the PD baselines and shares the -57.3%",
        "disease-progression effect with HD (control stream `PopHDF *",
        "THETA(31) + PopPDF * THETA(31)`); PD has no red-cell-lifespan",
        "effect."
      ),
      source_name = "POP == 3"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, the largest group, and other races)",
      notes = paste(
        "Control stream RaceBlackF (IARACEN = 11). Lowers the apparent",
        "baseline of ESA non-users and the true baseline of ESA users by",
        "3.76%, and the drug-effect slope by 12.8%."
      ),
      source_name = "IARACEN == 11"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, the largest group, and other races)",
      notes = paste(
        "Control stream RaceAsianF (IARACEN 13-16 = Central/South, East,",
        "Japanese, South-East Asian, or 98 = mixed Asian). Lowers the true",
        "baseline of ESA users by 4.67%; no effect on ESA non-users."
      ),
      source_name = "IARACEN in 13-16, 98"
    ),
    CONMED_CLOPIDOGREL = list(
      description = "Concomitant clopidogrel indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no clopidogrel)",
      notes = paste(
        "Time-varying (Table S6). Multiplies the effective dose by",
        "(1 + e_conmed_clopidogrel_dose) = 0.898 (control stream",
        "`CLOP_DOSE`)."
      ),
      source_name = "CLPDGREL"
    ),
    STUDY_205270 = list(
      description = "Study 205270 indicator (ESA non-users, not on dialysis, once daily; placebo-controlled)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 200807, the reference for residual error)",
      notes = paste(
        "Scales the apparent-baseline random effect by e_study_205270_etarbase",
        "= 0.729 and selects the study's additive residual SD (5.72 g/L).",
        "The ESA-user indicator is forced to 0 in this study."
      ),
      source_name = "STUDYID == 205270"
    ),
    STUDY_204837 = list(
      description = "Study 204837 indicator (ESA users on hemodialysis, daprodustat three times weekly)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (once-daily dosing; study 200807 for residual error)",
      notes = paste(
        "The only three-times-weekly study. Multiplies the drug effect by",
        "e_study_204837_slope = 0.422 (Table 2 'Factor on CDOSE for t.i.w.",
        "dosing in 204,837'; the control stream applies it to the effect term,",
        "`EFF = 1 + SLOPE * (DDR**GAMMA) * TIW`), and selects the study's",
        "additive residual SD (6.28 g/L). Set to 1 to simulate t.i.w. dosing",
        "with DOSE_DAPRODUSTAT_MG as the per-administration dose."
      ),
      source_name = "STUDYID == 204837"
    ),
    STUDY_201410 = list(
      description = "Study 201410 indicator (ESA non-users just starting or about to start hemodialysis or peritoneal dialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (chronic dialysis or non-dialysis; study 200807 for residual error)",
      notes = paste(
        "The control stream assigns every subject of this study to the",
        "'starting dialysis' population (`IF(STUDYID.EQ.201410) POP = 4`),",
        "so this indicator carries the starting-dialysis effects: +30.9%",
        "red-cell lifespan, +96.8% disease progression, and the",
        "starting-HD / starting-PD apparent baselines (modality from",
        "RRT_HEMODIAL_STATUS / PERIT_DIAL). It also scales the",
        "apparent-baseline random effect by 1.13 and selects the study's",
        "additive residual SD (7.24 g/L). The ESA-user indicator is forced to",
        "0 in this study."
      ),
      source_name = "STUDYID == 201410"
    ),
    STUDY_200808 = list(
      description = "Study 200808 indicator (ESA non-users and users not on dialysis, once daily)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 200807)",
      notes = "Selects the study's additive residual SD (6.39 g/L); no structural effect.",
      source_name = "STUDYID == 200808"
    )
  )

  compartmentData <- list(
    precursor1 = list(analyte = "hemoglobin", units = "g/L", specimen = "not applicable", verified = TRUE),
    hb1 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb2 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb3 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb4 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb5 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb6 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb7 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb8 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb9 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb10 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb11 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE),
    hb12 = list(analyte = "hemoglobin", units = "g/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2770L,
    n_studies = 5L,
    n_observations = 53535L,
    age_range = "17-96 years",
    age_median = "63 years",
    weight_range = "36-194 kg",
    weight_median = "73 kg",
    sex_female_pct = 50,
    race_ethnicity = "White 62%, Asian 18%, Black 13%, American Indian / Alaska Native 4% (Table S7 percentages as printed)",
    disease_state = paste(
      "Adults with anemia of chronic kidney disease: non-dialysis (ND) ESA",
      "non-users (studies 205270, 200808) and ESA users (200808); prevalent",
      "hemodialysis or peritoneal dialysis ESA users (204837, 200807); and",
      "ESA non-users just starting or about to start dialysis (201410).",
      "Mean baseline hemoglobin 100 g/L (SD 10)."
    ),
    renal_function = "Non-dialysis 1497 (54%; CKD stages 2-5), hemodialysis 1140 (41%), peritoneal dialysis 133 (5%)",
    dose_range = paste(
      "Daprodustat titrated to a hemoglobin target (10-11 g/dL; 11-12 g/dL",
      "in study 205270): once daily 1-24 mg, or three times weekly 2-48 mg",
      "(study 204837). Starting doses 1-4 mg QD (ESA non-users), 4-12 mg QD",
      "or 8-24 mg t.i.w. (ESA users), Table S4."
    ),
    regions = "Global (five pivotal phase III studies)",
    notes = paste(
      "Table S7 sex and race counts do not add up to N; the printed percentages are used. Analysis data set A: a random two-thirds of the evaluable daprodustat",
      "subjects of studies 205270 (n = 205), 204837 (180), 201410 (105),",
      "200807 (988) and 200808 (1292); the remaining 1384 subjects formed",
      "the internal validation data set B (Table S4). Treatment duration 28",
      "weeks to about 4 years. Control (ESA, placebo) arms were not",
      "analysed. Demographics from the Results and Table S7; Table S7 lists",
      "1701 non-dialysis subjects (57%), which does not match its own",
      "per-study counts (205 + 1292 = 1497) and is not used here."
    )
  )

  ini({
    # Apparent hemoglobin baseline at study start (control stream tvBL,
    # THETA(1)-THETA(6)). For ESA users this is the ESA-elevated value before
    # washout; for ESA non-users it is also the true baseline.
    rbase_nd_esanaive <- 93.3; label("Apparent Hgb baseline, non-dialysis ESA non-users (g/L)") # Table 2 final: 93.3 (RSE 0.4%)
    rbase_nd_esa <- 103; label("Apparent Hgb baseline, non-dialysis ESA users (g/L)") # Table 2 final: 103 (RSE 0.3%)
    rbase_hd_esa <- 105; label("Apparent Hgb baseline, hemodialysis ESA users (g/L)") # Table 2 final: 105 (RSE 0.2%)
    rbase_pd_esa <- 103; label("Apparent Hgb baseline, peritoneal-dialysis ESA users (g/L)") # Table 2 final: 103 (RSE 0.7%)
    rbase_hd_start <- 95.5; label("Hgb baseline, just started or starting HD in study 201410 (g/L)") # Table 2 final: 95.5 (RSE 1.1%)
    rbase_pd_start <- 101; label("Hgb baseline, just started or starting PD in study 201410 (g/L)") # Table 2 final: 101 (RSE 2.6%)

    # True baseline of ESA users after ESA washout (control stream tvBHB2,
    # THETA(7)-THETA(8); HD and PD share THETA(8)).
    rbase_true_nd_esa <- 95.3; label("True Hgb baseline, non-dialysis ESA users (g/L)") # Table 2 final: 95.3 (RSE 0.8%)
    rbase_true_dial_esa <- 94.9; label("True Hgb baseline, HD/PD ESA users (g/L)") # Table 2 final: 94.9 (RSE 0.7%)

    # Baseline covariate effects
    e_esad_rbase_true <- -0.00797; label("Prior ESA dose on true baseline, exp(theta/1000 * ESAD) (per IU/week)") # Table 2 final: -0.00797 (RSE 6.8%)
    e_esa_hyporesponder_rbase <- -0.0307; label("Relative effect of ESA hyporesponder on apparent baseline of ESA users (fraction)") # Table 2 final: -0.0307 (RSE 13.3%)
    e_ckd_stage_rbase <- -0.107; label("Power exponent of (CKD stage / 5) on true and apparent baseline (unitless)") # Table 2 final: -0.107 (RSE 12.9%)
    e_wt_rbase <- 0.0520; label("Power exponent of (WT / 73.3 kg) on true and apparent baseline (unitless)") # Table 2 final: 0.0520 (RSE 16.8%)
    e_hepcidin_bl_rbase <- -0.0170; label("Power exponent of (baseline hepcidin / 121.356 ng/mL) on true and apparent baseline (unitless)") # Table 2 final: -0.0170 (RSE 16.0%)
    e_race_asian_rbase_true <- -0.0467; label("Relative effect of Asian race on true baseline of ESA users (fraction)") # Table 2 final: -0.0467 (RSE 21.1%)
    e_race_black_rbase <- -0.0376; label("Relative effect of Black race on true and apparent baseline (fraction)") # Table 2 final: -0.0376 (RSE 17.5%)

    # Precursor and red-cell lifespans (both fixed to the Dose-Hgb Model 1/2
    # values; control stream THETA(16) = 3.19 FIX, THETA(17) = 82.819 FIX).
    lmtt <- fixed(log(3.19)); label("Precursor-cell mean transit t (day)") # Table 2: 3.19 (Fixed)
    lmtt_rbc <- fixed(log(82.819)); label("Red-blood-cell lifespan over the 12 age compartments (day)") # Table 2: 82.8 (Fixed); control stream THETA(17) = 82.819 FIX
    e_esa_mtt_rbc <- 0.100; label("Relative effect of prior ESA use on red-cell lifespan (fraction)") # Table 2 final: 0.100 (RSE 33.0%)
    e_wt_mtt_rbc <- 0.271; label("Power exponent of (WT / 73.3 kg) on red-cell lifespan (unitless)") # Table 2 final: 0.271 (RSE 20.7%)
    e_hd_mtt_rbc <- 0.138; label("Relative effect of chronic hemodialysis on red-cell lifespan (fraction)") # Table 2 final: 0.138 (RSE 30.3%)
    e_study_201410_mtt_rbc <- 0.309; label("Relative effect of just starting/started HD or PD (study 201410) on red-cell lifespan (fraction)") # Table 2 final: 0.309 (RSE 41.3%)

    # Drug effect: EFF = 1 + SLOPE * DDR^GAMMA, GAMMA = ratio * SLOPE
    lslope <- log(0.0580); label("Drug-effect slope SLOPE (1/mg)") # Table 2 final: 0.0580 (RSE 1.9%)
    gamma_slope_ratio <- 14.9; label("Ratio GAMMA / SLOPE of the dose power exponent to the slope (mg)") # Table 2 final: 'Correlation Slope * Gamma' 14.9 (RSE 2.8%)
    e_race_black_slope <- -0.128; label("Relative effect of Black race on SLOPE (fraction)") # Table 2 final: -0.128 (RSE 11.1%)
    e_wt_dose <- fixed(0.75); label("Allometric exponent of (WT / 73.3 kg) dividing the dose (unitless)") # Table 2: 0.750 (Fixed)
    e_study_204837_slope <- 0.422; label("Multiplicative factor on the drug effect for t.i.w. dosing in study 204837 (unitless)") # Table 2 final: 0.422 (RSE 9.5%)
    e_conmed_clopidogrel_dose <- -0.102; label("Relative effect of clopidogrel use on the effective dose (fraction)") # Table 2 final: -0.102 (RSE 16.7%)

    # Disease progression of the precursor production rate
    kprog <- -0.0851; label("Disease progression rate (1/(1000 day))") # Table 2 final: -0.0851 (RSE 8.3%)
    kin_min <- -0.678; label("Minimum (time-infinity) precursor production rate (g/L/day)") # Table 2 final: -0.678 (RSE 11.1%)
    e_esa_kprog <- -0.318; label("Relative effect of prior ESA use on disease progression rate (fraction)") # Table 2 final: -0.318 (RSE 32.7%)
    e_dial_kprog <- -0.573; label("Relative effect of chronic HD or PD on disease progression rate (fraction)") # Table 2 final: -0.573 (RSE 17.7%)
    e_study_201410_kprog <- 0.968; label("Relative effect of just starting/started HD or PD (study 201410) on disease progression rate (fraction)") # Table 2 final: 0.968 (RSE 34.2%)

    # Study-specific scale of the apparent-baseline random effect (OM1F)
    e_study_205270_etarbase <- 0.729; label("Scale factor on the apparent-baseline random effect in study 205270 (unitless)") # Table 2 final: 'Study 205270 on omega 1' 0.729 (RSE 7.2%)
    e_study_201410_etarbase <- 1.13; label("Scale factor on the apparent-baseline random effect in study 201410 (unitless)") # Table 2 final: 'Study 201410 on omega 1' 1.13 (RSE 9.3%)

    # IIV. Baselines and disease-progression rates are additive (g/L and
    # 1/(1000 day)); SLOPE and red-cell lifespan are exponential. Table 2
    # reports variances and covariances.
    etarbase + etarbase_true ~ c(
      50.1,
      40.3, 131
    ) # Table 2 final: var apparent baseline 50.1, cov 40.3 (r = 0.497), var true baseline 131
    etalslope ~ 0.0280 # Table 2 final: variance slope (exponential) 0.0280
    etalmtt_rbc ~ 0.206 # Table 2 final: variance TRBC (exponential) 0.206
    etabvv_kprog_1 + etabvv_kprog_2 + etabvv_kprog_3 ~ c(
      0.0504,
      0.0212, 0.0178,
      0.0121, 0.0108, 0.0104
    ) # Table 2 final: year-1, year-2 and year-3+ disease progression variances 0.0504, 0.0178, 0.0104; covariances 0.0212 (1-2), 0.0121 (1-3), 0.0108 (2-3)

    # Additive residual SD by study; study 200807 is the reference.
    addSd <- 7.27; label("Additive residual SD, study 200807 (g/L)") # Table 2 final: 7.27 (RSE 0.2%)
    addSd_205270 <- 5.72; label("Additive residual SD, study 205270 (g/L)") # Table 2 final: 5.72 (RSE 1.0%)
    addSd_204837 <- 6.28; label("Additive residual SD, study 204837 (g/L)") # Table 2 final: 6.28 (RSE 0.9%)
    addSd_201410 <- 7.24; label("Additive residual SD, study 201410 (g/L)") # Table 2 final: 7.24 (RSE 1.2%)
    addSd_200808 <- 6.39; label("Additive residual SD, study 200808 (g/L)") # Table 2 final: 6.39 (RSE 0.2%)
  })

  model({
    # Patient group (control stream POP / EPO recoding). Study 201410 is the
    # starting-dialysis population (POP = 4); its modality comes from
    # RRT_HEMODIAL_STATUS / PERIT_DIAL. Studies 205270 and 201410 force the
    # ESA-user flag to 0.
    dial_start <- STUDY_201410
    hd_chronic <- RRT_HEMODIAL_STATUS * (1 - dial_start)
    pd_chronic <- PERIT_DIAL * (1 - dial_start)
    nd <- (1 - RRT_HEMODIAL_STATUS) * (1 - PERIT_DIAL)
    esa <- (ESAD > 0) * (1 - STUDY_205270) * (1 - STUDY_201410)

    # Typical apparent baseline. Combinations absent from the studies
    # (chronic dialysis ESA non-users, starting-dialysis ESA users) have no
    # estimated baseline and give 0.
    rbase_app_tv <- nd * (1 - esa) * rbase_nd_esanaive +
      nd * esa * rbase_nd_esa +
      hd_chronic * esa * rbase_hd_esa +
      pd_chronic * esa * rbase_pd_esa +
      dial_start * RRT_HEMODIAL_STATUS * rbase_hd_start +
      dial_start * PERIT_DIAL * rbase_pd_start
    rbase_true_tv <- nd * rbase_true_nd_esa + (hd_chronic + pd_chronic) * rbase_true_dial_esa

    # Baseline covariate factors (control stream BLSTAG, BLWT, BLHEP, PEPL,
    # BLHYP, RACE_TBSL, RACE_BBSL)
    bl_stage <- (CKD_STAGE / 5)^e_ckd_stage_rbase
    bl_wt <- (WT / 73.3)^e_wt_rbase
    bl_hep <- (HEPCIDIN_BL / 121.356)^e_hepcidin_bl_rbase
    bl_esad <- exp(e_esad_rbase_true / 1000 * ESAD)
    bl_hyp <- 1 + e_esa_hyporesponder_rbase * ESA_HYPORESPONDER
    bl_asian <- 1 + e_race_asian_rbase_true * RACE_ASIAN
    bl_black <- 1 + e_race_black_rbase * RACE_BLACK
    om1f <- 1 + (e_study_205270_etarbase - 1) * STUDY_205270 +
      (e_study_201410_etarbase - 1) * STUDY_201410

    # Apparent baseline (BHB) and true baseline after ESA washout (BHB2);
    # for ESA non-users the two coincide.
    rbase_app <- (1 - esa) * rbase_app_tv * bl_stage * bl_wt * bl_hep * bl_black +
      esa * rbase_app_tv * bl_hyp + om1f * etarbase
    rbase_true <- (1 - esa) * rbase_app +
      esa * (rbase_true_tv * bl_stage * bl_wt * bl_hep * bl_esad * bl_asian * bl_black + etarbase_true)

    # Lifespans and transit rates (NRBC = 12 age compartments)
    mtt <- exp(lmtt)
    mtt_rbc <- exp(lmtt_rbc + etalmtt_rbc) * (1 + e_esa_mtt_rbc * esa) *
      (WT / 73.3)^e_wt_mtt_rbc *
      (1 + e_hd_mtt_rbc * hd_chronic + e_study_201410_mtt_rbc * dial_start)
    kprec <- 1 / mtt
    krbc <- 12 / mtt_rbc

    # Production rate: the apparent baseline holds before the first dose
    # (t < 0) and the true baseline from the first dose on, so ESA users
    # wash out to their true baseline (control stream POSTD1 switch).
    kin <- ((t >= 0) * rbase_true + (t < 0) * rbase_app) / mtt_rbc

    # Drug effect on precursor production
    slope <- exp(lslope + etalslope) * (1 + e_race_black_slope * RACE_BLACK)
    gamma <- gamma_slope_ratio * slope
    tiw <- 1 + (e_study_204837_slope - 1) * STUDY_204837
    dose_eff <- DOSE_DAPRODUSTAT_MG / (WT / 73.3)^e_wt_dose *
      (1 + e_conmed_clopidogrel_dose * CONMED_CLOPIDOGREL)
    eff <- 1 + slope * dose_eff^gamma * tiw

    # Disease progression (Mahar 2025 Eq. 1; t in days since the first
    # dose). Year-specific additive random effects switch at 367 and 732
    # days after the first dose (control stream RTFD 8808 h and 17568 h).
    yr1 <- (t < 367)
    yr2 <- (t >= 367) * (t < 732)
    yr3 <- (t >= 732)
    dispr <- (kprog * (1 + e_esa_kprog * esa) *
      (1 + e_dial_kprog * (hd_chronic + pd_chronic) + e_study_201410_kprog * dial_start) +
      yr1 * etabvv_kprog_1 + yr2 * etabvv_kprog_2 + yr3 * etabvv_kprog_3) / 1000
    prod <- max((kin - kin_min) * exp(dispr * t) + kin_min, 0)

    d/dt(precursor1) <- prod * eff - kprec * precursor1
    d/dt(hb1) <- kprec * precursor1 - krbc * hb1
    d/dt(hb2) <- krbc * (hb1 - hb2)
    d/dt(hb3) <- krbc * (hb2 - hb3)
    d/dt(hb4) <- krbc * (hb3 - hb4)
    d/dt(hb5) <- krbc * (hb4 - hb5)
    d/dt(hb6) <- krbc * (hb5 - hb6)
    d/dt(hb7) <- krbc * (hb6 - hb7)
    d/dt(hb8) <- krbc * (hb7 - hb8)
    d/dt(hb9) <- krbc * (hb8 - hb9)
    d/dt(hb10) <- krbc * (hb9 - hb10)
    d/dt(hb11) <- krbc * (hb10 - hb11)
    d/dt(hb12) <- krbc * (hb11 - hb12)

    # Steady state at the apparent baseline (control stream A_0)
    precursor1(0) <- rbase_app * mtt / mtt_rbc
    hb1(0) <- rbase_app / 12
    hb2(0) <- rbase_app / 12
    hb3(0) <- rbase_app / 12
    hb4(0) <- rbase_app / 12
    hb5(0) <- rbase_app / 12
    hb6(0) <- rbase_app / 12
    hb7(0) <- rbase_app / 12
    hb8(0) <- rbase_app / 12
    hb9(0) <- rbase_app / 12
    hb10(0) <- rbase_app / 12
    hb11(0) <- rbase_app / 12
    hb12(0) <- rbase_app / 12

    hb <- hb1 + hb2 + hb3 + hb4 + hb5 + hb6 + hb7 + hb8 + hb9 + hb10 + hb11 + hb12

    addSd_i <- addSd * (1 - STUDY_205270 - STUDY_204837 - STUDY_201410 - STUDY_200808) +
      addSd_205270 * STUDY_205270 + addSd_204837 * STUDY_204837 +
      addSd_201410 * STUDY_201410 + addSd_200808 * STUDY_200808
    hb ~ add(addSd_i)
  })
}
