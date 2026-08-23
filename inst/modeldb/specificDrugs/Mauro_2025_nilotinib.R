Mauro_2025_nilotinib <- function() {
  description <- "Two-compartment population PK model for nilotinib (Danziten reduced-dose tablets and Tasigna capsules) in healthy adults, with lagged zero-order absorption into the central compartment, a delayed secondary depot for late re-absorption, and formulation-, dose- and prandial-state-dependent relative bioavailability"
  reference <- paste(
    "Mauro M, Radich J, Jain P, Sequeira D, Bellanti F, Douer D.",
    "Pharmacokinetic profile of novel reduced-dose Danziten (nilotinib tablets)",
    "versus Tasigna (nilotinib capsules): in vivo bioequivalence and population",
    "pharmacokinetic analysis. Cancer Chemother Pharmacol. 2025;95(1):56.",
    "doi:10.1007/s00280-025-04777-6.",
    "Parameters are the 'Final Model' of Supplemental Table 2, fitted to 13 of the",
    "14 pooled studies and used for the published fasted steady-state bioequivalence",
    "simulations (main text Table 1). The companion 'Updated Model' that adds the",
    "modified-fasted low-fat-meal study is",
    "modellib('Mauro_2025_nilotinib_lowfat').",
    sep = " "
  )
  vignette <- "Mauro_2025_nilotinib"

  # Every nilotinib administration is entered as TWO parallel dose records with the
  # same amt (Supplemental Information, 'Initial Model Development'):
  #   (a) cmt = "central", rate = -2 -- the primary lagged zero-order input, with
  #       modelled duration D1 and lag ALAG1 supplied by dur()/alag() below;
  #   (b) cmt = "depot2"            -- the empirical secondary depot that carries the
  #       late (TSLD ~ 12-24 h) re-absorption, first-order with its own lag.
  # Declared explicitly because the registry's default detection only recognises
  # doses into depot and central.
  dosing <- c("central", "depot2")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units             = "kg",
      type              = "continuous",
      reference_category = NULL,
      notes             = paste(
        "Power (allometric-style) effect on CL. Mauro 2025 does not state the",
        "reference weight; 70 kg is assumed here per the rounded-standard",
        "convention. NOTE the published exponent is NEGATIVE (-0.727, 95% CI",
        "-0.934 to -0.532), i.e. heavier subjects have LOWER nilotinib clearance;",
        "this is not a sign transcription error -- both Supplemental Table 2 and",
        "Supplemental Table 3 report it below zero with confidence intervals",
        "entirely below zero.",
        sep = " "
      ),
      source_name       = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (male)",
      notes             = paste(
        "Two separate effects. (1) On CL for every subject (-22.7%: female CL is",
        "lower). (2) On relative bioavailability, but ONLY for female subjects",
        "receiving nilotinib TABLETS under FASTED conditions (+11.6%); the",
        "Supplemental Table 2 footnote reads 'Only female receiving nilotinib",
        "Tablets under fasted conditions', and the supplement text states no",
        "male/female absorption difference was seen under fed or modified fasted",
        "conditions. The gating is implemented as SEXF * FORM_TABLET * (1 - FED).",
        "Only 39 of the 502 pooled completers were female (studies 006-23, 007-23",
        "and 128-22).",
        sep = " "
      ),
      source_name       = "SEX"
    ),
    FORM_TABLET = list(
      description        = "Formulation indicator, 1 = nilotinib tablet, 0 = nilotinib capsule",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (Tasigna nilotinib hydrochloride monohydrate capsule -- the F1 = 1 anchor)",
      notes             = paste(
        "The comparator here is a CAPSULE, not the non-tablet oral liquid named in",
        "the FORM_TABLET register default (same situation as Wada 2023 sparsentan).",
        "1 = the novel Danziten film-coated tablet (nilotinib tartrate, modified",
        "composition and process for bioavailability enhancement); 0 = the approved",
        "Tasigna capsule (nilotinib hydrochloride monohydrate). Per dose record, so",
        "a crossover subject contributes both levels. Carries effects on relative",
        "bioavailability (+98.8% at the tablet reference dose) and on the zero-order",
        "absorption duration D1 (-42.3%), and selects between the formulation-",
        "specific food effects.",
        sep = " "
      ),
      source_name       = "FORM"
    ),
    DOSE = list(
      description        = "Administered nilotinib dose for the current dose record",
      units             = "mg",
      type              = "continuous",
      reference_category = NULL,
      notes             = paste(
        "The TOTAL dose administered at the event, not the per-unit strength: a",
        "190 mg tablet dose is 2 x 95 mg tablets and enters as DOSE = 190.",
        "Two uses. (1) A shared power effect on relative bioavailability",
        "normalised to 400 mg for capsules and 200 mg for tablets, with no further",
        "dose effect for tablets above 200 mg (implemented with min(DOSE, 200)),",
        "reproducing the less-than-dose-proportional capsule exposure and the",
        "dose-proportional tablet exposure above 190 mg. (2) Selection among the",
        "dose-specific tablet food effects on relative bioavailability.",
        sep = " "
      ),
      source_name       = "DOSE"
    ),
    FED = list(
      description        = "Any-food-vs-fasted indicator for the dose record, 1 = a meal preceded the dose",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (fasted: overnight fast of >= 10 h with no meal before dosing)",
      notes             = paste(
        "Carries the food effects that Mauro 2025 estimated WITHOUT distinguishing",
        "the prandial sub-state: on the zero-order absorption duration D1 and the",
        "absorption lag ALAG1 (both formulation-specific) and on the central volume",
        "V1 and inter-compartmental clearance Q (both shared across formulations).",
        "Also gates the female tablet bioavailability effect to fasted records.",
        "In this Final Model every fed record is a high-fat meal, so FED and",
        "FED_HIGHFAT coincide; both are carried because the register prescribes an",
        "any-food indicator alongside the composition indicators, and because the",
        "companion low-fat model separates them.",
        sep = " "
      ),
      source_name       = "FOOD"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal indicator for the dose record, 1 = a high-fat meal preceded the dose",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (fasted)",
      notes             = paste(
        "Mauro 2025 Methods, 'Prandial state definitions': the high-fat meal is a",
        "high-fat high-calorie breakfast of 800-1000 kcal, non-vegetarian, about 50%",
        "fat content, consumed within 30 min. Paired with MEAL_PREDOSE_2H, which",
        "says whether the dose followed the START of that meal by about 0.5 h",
        "(the paper's 'fed' state) or by about 2 h (the paper's 'modified fasting",
        "with a high-fat meal' state). Both timings carry DIFFERENT bioavailability",
        "effects for tablets, which is why the timing indicator is required.",
        sep = " "
      ),
      source_name       = "FOOD"
    ),
    MEAL_PREDOSE_2H = list(
      description        = "Dose-to-meal interval indicator, 1 = the dose was taken about 2 h after the start of the meal rather than about 0.5 h after it",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (dose taken about 0.5 h after the start of the meal, or fasted)",
      notes             = paste(
        "Mauro 2025 calls the 2 h interval 'modified fasting' (permitted by the",
        "approved nilotinib capsule label) and the 0.5 h interval 'fed'. The meal",
        "composition is identical between the two high-fat states, so this",
        "indicator is the ONLY thing separating them, and their estimated effects",
        "on tablet bioavailability differ (at 142 mg: +48.4% at 0.5 h vs +58.5% at",
        "2 h; at 190 mg: +61.7% vs +60.3%). NOTE that the phrase 'modified fasting'",
        "is not portable between papers -- read each paper's own definition.",
        sep = " "
      ),
      source_name       = "FOOD"
    ),
    OCC = list(
      description        = "Occasion (crossover period) index for the dose record",
      units             = "(integer)",
      type              = "categorical",
      reference_category = NULL,
      notes             = paste(
        "Drives the inter-occasion variability on the central volume V1",
        "(Supplemental Table 2, 'Inter-Occasion Variability'). Mauro 2025 does not",
        "state how many occasions the IOV spans; three are implemented here because",
        "the largest crossover design among the pooled studies is a three-way,",
        "three-period crossover (Supplemental Table 1). Set OCC = 1, 2 or 3; for a",
        "single-occasion simulation set OCC = 1 throughout.",
        sep = " "
      ),
      source_name       = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Mauro 2025 Methods, 'Analyses') but not retained in the final model; no point estimate is published."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis but not retained; no point estimate is published."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis but not retained; no point estimate is published."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis but not retained; no point estimate is published."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis but not retained; no point estimate is published."
    )
  )

  compartmentData <- list(
    depot2 = list(analyte = "nilotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nilotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nilotinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 502,
    n_studies = 13,
    age_range = "adults (not reported)",
    weight_range = "not reported",
    sex_female_pct = 7.8,
    disease_state = "healthy",
    dose_range = "50-300 mg single oral dose (nilotinib tablets); 50-400 mg single oral dose (nilotinib capsules)",
    regions = "India (Phase 1 units, ZenRise, Hyderabad); analysis by GK Analytics, Hyderabad and Certara",
    notes = paste(
      "Healthy adult men and women in 14 single-dose relative-bioavailability and",
      "food-effect crossover studies (Supplemental Table 1). The paper reports the",
      "pooled analysis dataset as 30,089 nilotinib observations from 502 subjects",
      "with 1,311 below-limit-of-quantification records excluded, and does not break",
      "the counts down per model; this Final Model was fitted to 13 of the 14 studies",
      "(study 052-23, the modified-fasted low-fat-meal study, is the one added in the",
      "companion Updated Model). 39 female participants contributed, from studies",
      "006-23, 007-23 and 128-22; every other study enrolled healthy men only, which",
      "is why sex is identifiable only through those three studies. Baseline",
      "demographic distributions (age, weight, laboratory values) are not tabulated",
      "in the paper or its Supplemental Information. The published steady-state",
      "simulations resampled 50 subjects (50% male / 50% female) from the analysis",
      "dataset.",
      sep = " "
    )
  )

  ini({
    # ---- Disposition -- Supplemental Table 2, 'Typical Values' ----
    lcl <- log(34.5); label("Clearance (CL, L/h)") # Suppl. Table 2: CL = 34.5 L/h (RSE 2.85%, 95% CI 32.4-36.2; bootstrap 34.5)
    lvc <- log(557); label("Central volume of distribution (V1, L)") # Suppl. Table 2: V1 = 557 L (RSE 2.66%, 95% CI 525-586; bootstrap 556)
    lq <- log(4.54); label("Inter-compartmental clearance (Q, L/h)") # Suppl. Table 2: Q = 4.54 L/h (RSE 16.6%, 95% CI 3.26-6.73; bootstrap 4.58)
    lvp <- log(223); label("Peripheral volume of distribution (V2, L)") # Suppl. Table 2: V2 = 223 L (RSE 8.78%, 95% CI 185-273; bootstrap 224)

    # ---- Primary absorption: lagged zero-order input into central ----
    # Supplemental Information, 'Initial Model Development': "The oral absorption of
    # nilotinib was well described by a zero-order process with a lag-time"; first-
    # order, mixed zero-then-first-order and transit-compartment alternatives were
    # tested and rejected. D1/ALAG1/F1 are NONMEM compartment-1 attributes and
    # compartment 1 is the central compartment in a zero-order-input ADVAN.
    ld1 <- log(2.32); label("Duration of the zero-order absorption input into central (D1, h)") # Suppl. Table 2: D1 = 2.32 h (RSE 1.63%, 95% CI 2.07-2.40; bootstrap 2.32)
    ltlag <- log(0.0988); label("Absorption lag time of the zero-order input (ALAG1, h)") # Suppl. Table 2: ALAG1 = 0.0988 h (RSE 14.8%, 95% CI 0.0719-0.136; bootstrap 0.102)
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the capsule reference at 400 mg fasted (F1, fraction)") # Suppl. Info F1 equations: F1_capsule = 1 * (DOSE/400)^theta10, i.e. F1 is anchored at 1 for the capsule reference dose

    # ---- Secondary depot: delayed late re-absorption ----
    # Supplemental Information, 'Initial Model Development': an increase in
    # concentrations at TSLD of about 12 and 24 h indicated delayed absorption or
    # re-absorption, described empirically by routing an additional dose record of
    # the same amount through a secondary depot with first-order absorption and its
    # own lag (delta OFV -575.6).
    lka2 <- log(0.140); label("First-order absorption rate constant from the secondary depot (1/h)") # Suppl. Table 2: absorption from secondary depot = 0.140 1/h (RSE 5.46%, 95% CI 0.127-0.159; bootstrap 0.141)
    ltlag2 <- log(5.35); label("Absorption lag time of the secondary depot (h)") # Suppl. Table 2: absorption lag re-absorbed from secondary depot = 5.35 h (RSE 6.24%, 95% CI 4.65-6.08; bootstrap 5.34)
    lfdepot2 <- log(1 - 0.888); label("Fraction of F1 that is re-absorbed through the secondary depot (fraction)") # Suppl. Table 2: fraction re-absorbed from secondary depot (proportional) = -0.888 (RSE 0.656%, 95% CI -0.900 to -0.877), i.e. 1 - 0.888 = 0.112 of F1

    # ---- Covariate effects on relative bioavailability -- Suppl. Table 2 ----
    e_form_tablet_fdepot <- 0.988; label("Proportional effect of the tablet formulation on F1 (unitless)") # Suppl. Table 2: formulation (nilotinib tablets) effect on F1 (proportional) = 0.988 (RSE 3.84%, 95% CI 0.913-1.06); matches the main text "98.8% increase at 190 mg as compared to 400 mg capsules"
    e_dose_fdepot <- -0.701; label("Power exponent of dose on F1, shared by both formulations (unitless)") # Suppl. Table 2: dose effect on F1 (power) = -0.701 (RSE 6.79%, 95% CI -0.811 to -0.609); Suppl. Info: same exponent for tablets and capsules, reference doses 200 mg (tablets) and 400 mg (capsules)
    e_sexf_fdepot <- 0.116; label("Proportional effect of female sex on F1 for tablets under fasted conditions (unitless)") # Suppl. Table 2: sex (female*) covariate effect on F1 (proportional) = 0.116 (RSE 36.8%, 95% CI 0.0287-0.202); footnote * = "Only female receiving nilotinib Tablets under fasted conditions"

    # High-fat meal taken about 0.5 h before the dose (the paper's 'fed' state).
    # Tablet effects are dose-specific; the single capsule estimate is shared with
    # the 2 h ('modified fasting') timing -- see e_hfmeal_cap_fdepot.
    e_hfmeal05h_tab50_fdepot <- 0.166; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 50 mg tablets (unitless)") # Suppl. Table 2: food (fed) effect on F1 at 50 mg for nilotinib tablets = 0.166 (RSE 33.9%, 95% CI 0.0363-0.268)
    e_hfmeal05h_tab142_fdepot <- 0.484; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 142 mg tablets (unitless)") # Suppl. Table 2: food (fed) effect on F1 at 142 mg for nilotinib tablets = 0.484 (RSE 11.4%, 95% CI 0.384-0.603); Suppl. Info text quotes 48.4% at 142 mg
    e_hfmeal05h_tab190_fdepot <- 0.617; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 190 mg tablets (unitless)") # Suppl. Table 2: food (fed) effect on F1 at 190 mg for nilotinib tablets = 0.617 (RSE 13.1%, 95% CI 0.471-0.774); Suppl. Info text quotes 61.7% at 190 mg
    e_hfmeal05h_tab240_fdepot <- 0.748; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 240 mg tablets (unitless)") # Suppl. Table 2: food (fed) effect on F1 at 240 mg for nilotinib tablets = 0.748 (RSE 10.0%, 95% CI 0.589-0.909)
    e_hfmeal05h_tab300_fdepot <- 0.614; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 300 mg tablets (unitless)") # Suppl. Table 2: food (fed) effect on F1 at 300 mg for nilotinib tablets = 0.614 (RSE 12.0%, 95% CI 0.485-0.783)

    # High-fat meal taken about 2 h before the dose (the paper's 'modified fasting
    # with a high-fat meal' state), tablets only; capsules share one estimate.
    e_hfmeal2h_tab142_fdepot <- 0.585; label("Proportional effect of a high-fat meal 2 h pre-dose on F1, 142 mg tablets (unitless)") # Suppl. Table 2: food (modified fasted) effect on F1 at 142 mg for nilotinib tablets = 0.585 (RSE 10.1%, 95% CI 0.481-0.700)
    e_hfmeal2h_tab190_fdepot <- 0.603; label("Proportional effect of a high-fat meal 2 h pre-dose on F1, 190 mg tablets (unitless)") # Suppl. Table 2 labels this row "at 190 mg for nilotinib capsules" = 0.603 (RSE 12.0%, 95% CI 0.461-0.748); the label is a typo -- Suppl. Table 3 gives the identically-placed row as "for nilotinib tablets", and capsules carry the single non-dose-specific 1.81 effect instead
    e_hfmeal_cap_fdepot <- 1.81; label("Proportional effect of a high-fat meal on F1, capsules, any dose-to-meal interval (unitless)") # Suppl. Table 2: food (modified fasted) effect on F1 for nilotinib capsules = 1.81 (RSE 7.19%, 95% CI 1.55-2.07); capsules were dosed under BOTH the 0.5 h and 2 h high-fat conditions (Suppl. Table 1) yet only one capsule F1 food effect was estimated, so it is shared across the two timings

    # ---- Covariate effects on absorption and distribution -- Suppl. Table 2 ----
    # These rows are labelled "Food effect ..." with no prandial sub-state, so they
    # apply to every non-fasted record.
    e_form_tablet_d1 <- -0.423; label("Proportional effect of the tablet formulation on D1 (unitless)") # Suppl. Table 2: formulation (nilotinib tablets) effect on D1 (proportional) = -0.423 (RSE 3.21%, 95% CI -0.450 to -0.371)
    e_fed_tab_d1 <- 2.91; label("Proportional effect of food on D1, tablets (unitless)") # Suppl. Table 2: food effect on D1 for nilotinib tablets (proportional) = 2.91 (RSE 3.58%, 95% CI 2.68-3.54)
    e_fed_cap_d1 <- 0.746; label("Proportional effect of food on D1, capsules (unitless)") # Suppl. Table 2: food effect on D1 for nilotinib capsules (proportional) = 0.746 (RSE 9.71%, 95% CI 0.634-0.978)
    e_fed_tab_tlag <- 2.98; label("Proportional effect of food on ALAG1, tablets (unitless)") # Suppl. Table 2: food effect on ALAG1 for nilotinib tablets (proportional) = 2.98 (RSE 20.5%, 95% CI 2.01-4.57)
    e_fed_cap_tlag <- 7.79; label("Proportional effect of food on ALAG1, capsules (unitless)") # Suppl. Table 2: food effect on ALAG1 for nilotinib capsules (proportional) = 7.79 (RSE 17.1%, 95% CI 5.43-11.4)
    e_fed_vc <- -0.248; label("Proportional effect of food on V1 (unitless)") # Suppl. Table 2: food effect on V1 (proportional) = -0.248 (RSE 11.8%, 95% CI -0.381 to -0.195); Suppl. Info: same effect for both formulations
    e_fed_q <- 2.07; label("Proportional effect of food on Q (unitless)") # Suppl. Table 2: food effect on Q (proportional) = 2.07 (RSE 22.8%, 95% CI 1.29-3.51); Suppl. Info: same effect for both formulations

    # ---- Covariate effects on clearance -- Suppl. Table 2 ----
    e_wt_cl <- -0.727; label("Power exponent of body weight on CL (unitless)") # Suppl. Table 2: bodyweight covariate effect on CL (power) = -0.727 (RSE 13.4%, 95% CI -0.934 to -0.532; bootstrap -0.722)
    e_sexf_cl <- -0.227; label("Proportional effect of female sex on CL (unitless)") # Suppl. Table 2: sex (female) covariate effect on CL (proportional) = -0.227 (RSE 16.2%, 95% CI -0.298 to -0.156)

    # ---- Between-subject variability -- Suppl. Table 2, 'Between Subject Variability' ----
    # Reported on the log-scale variance scale; correlation = 0.0628 / sqrt(0.108 *
    # 0.144) = 0.504. Suppl. Info: "a correlation was visible between CL and V1 and
    # the inclusion of an omega block between the two parameters resulted in a
    # significant reduction on OFV (-92.8 points)".
    etalcl + etalvc ~ c(0.108, 0.0628, 0.144) # Suppl. Table 2: BSV CL = 0.108 (RSE 8.89%, shrinkage 2.40%), correlation CL-V1 = 0.0628 (RSE 12.7%), BSV V1 = 0.144 (RSE 9.30%, shrinkage 1.60%)

    # ---- Inter-occasion variability on V1 -- Suppl. Table 2 ----
    # One variance is reported; the per-occasion etas share it, the NONMEM
    # $OMEGA BLOCK(1) SAME idiom. The number of occasions is not stated (see
    # covariateData$OCC notes).
    etaiov_vc_1 ~ 0.157 # Suppl. Table 2: IOV V1 = 0.157 (RSE 8.52%, shrinkage 4.60%; bootstrap 0.158)
    etaiov_vc_2 ~ fixed(0.157) # equal to occasion 1 per $OMEGA BLOCK(1) SAME
    etaiov_vc_3 ~ fixed(0.157) # equal to occasion 1 per $OMEGA BLOCK(1) SAME

    # ---- Residual error -- Suppl. Table 2, 'Residual Error' ----
    # Suppl. Info: "A combined (additive and proportional) error model was used to
    # describe the residual variability." Supplemental Table 2 labels the
    # proportional term "Exponential Error"; the narrative wording is followed here
    # and the term is encoded as proportional.
    propSd <- 0.230; label("Proportional residual error (fraction)") # Suppl. Table 2: exponential error = 0.230 (RSE 2.87%, 95% CI 0.216-0.244; bootstrap 0.230)
    addSd <- 36.6; label("Additive residual error (ng/mL)") # Suppl. Table 2: additive error = 36.6 ng/mL (RSE 6.75%, 95% CI 31.6-41.4; bootstrap 36.8)
  })

  model({
    # ---- 1. Prandial-state and formulation strata ----
    # Mauro 2025 Methods, 'Prandial state definitions', defines four mutually
    # exclusive states. FED_HIGHFAT carries the meal composition and
    # MEAL_PREDOSE_2H the dose-to-meal interval, so the two high-fat states are
    # separable:
    #   fasted                        FED = 0
    #   high-fat meal 0.5 h pre-dose  FED = 1, FED_HIGHFAT = 1, MEAL_PREDOSE_2H = 0
    #   high-fat meal 2 h pre-dose    FED = 1, FED_HIGHFAT = 1, MEAL_PREDOSE_2H = 1
    tab <- FORM_TABLET
    cap <- 1 - FORM_TABLET
    hf05h <- FED_HIGHFAT * (1 - MEAL_PREDOSE_2H)
    hf2h <- FED_HIGHFAT * MEAL_PREDOSE_2H

    # ---- 2. Relative bioavailability ----
    # Supplemental Information, 'Initial Model Development':
    #   F1_capsule = 1 * (DOSE/400)^theta10
    #   F1_tablet  = 1 * (1 + theta6) * (DOSE/200)^theta10   for DOSE <= 200 mg
    #   F1_tablet  = 1 * (1 + theta6)                        for DOSE >  200 mg
    # The tablet branch is the min(DOSE, 200) clamp below: above 200 mg the dose
    # ratio is pinned at 1 so tablet exposure stays dose-proportional.
    doseratio <- tab * min(DOSE, 200) / 200 + cap * DOSE / 400
    f1_fasted <- exp(lfdepot) * (1 + e_form_tablet_fdepot * tab) * doseratio^e_dose_fdepot

    # Dose-specific tablet food effects. Mauro 2025 estimated them only at 50, 142,
    # 190, 240 and 300 mg (0.5 h pre-dose) and at 142 and 190 mg (2 h pre-dose), so
    # the step boundaries below are midpoints between the published anchors and are
    # this encoding's choice, not the paper's. The effects are non-monotonic in dose
    # (0.617 at 190 mg vs 0.614 at 300 mg), so interpolating between anchors would
    # not be justified.
    fed05h_tab <- e_hfmeal05h_tab50_fdepot * (DOSE < 96) +
      e_hfmeal05h_tab142_fdepot * (DOSE >= 96) * (DOSE < 166) +
      e_hfmeal05h_tab190_fdepot * (DOSE >= 166) * (DOSE < 215) +
      e_hfmeal05h_tab240_fdepot * (DOSE >= 215) * (DOSE < 270) +
      e_hfmeal05h_tab300_fdepot * (DOSE >= 270)
    fed2h_tab <- e_hfmeal2h_tab142_fdepot * (DOSE < 166) +
      e_hfmeal2h_tab190_fdepot * (DOSE >= 166)

    foodf1 <- hf05h * (tab * fed05h_tab + cap * e_hfmeal_cap_fdepot) +
      hf2h * (tab * fed2h_tab + cap * e_hfmeal_cap_fdepot)

    # The female bioavailability effect applies only to tablets under fasted
    # conditions (Suppl. Table 2 footnote).
    f1 <- f1_fasted * (1 + foodf1) * (1 + e_sexf_fdepot * SEXF * tab * (1 - FED))

    # ---- 3. Individual parameters ----
    iov_vc <- (OCC == 1) * etaiov_vc_1 + (OCC == 2) * etaiov_vc_2 +
      (OCC == 3) * etaiov_vc_3

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc + iov_vc) * (1 + e_fed_vc * FED)
    q <- exp(lq) * (1 + e_fed_q * FED)
    vp <- exp(lvp)
    d1 <- exp(ld1) * (1 + e_form_tablet_d1 * tab) *
      (1 + FED * (tab * e_fed_tab_d1 + cap * e_fed_cap_d1))
    tlag <- exp(ltlag) * (1 + FED * (tab * e_fed_tab_tlag + cap * e_fed_cap_tlag))
    ka2 <- exp(lka2)
    tlag2 <- exp(ltlag2)
    fdepot2 <- exp(lfdepot2)

    # ---- 4. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 5. ODE system ----
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 6. Dose routing ----
    # Primary lagged zero-order input into central (requires rate = -2 on the
    # central dose record) and the secondary depot carrying fdepot2 * F1.
    f(central) <- f1
    dur(central) <- d1
    alag(central) <- tlag
    f(depot2) <- f1 * fdepot2
    alag(depot2) <- tlag2

    # ---- 7. Observation and residual error ----
    # Dose is in mg and vc in L, so central/vc is mg/L; the factor 1000 converts to
    # the ng/mL scale of Supplemental Table 2's additive error and of main-text
    # Tables 1 and 2.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
