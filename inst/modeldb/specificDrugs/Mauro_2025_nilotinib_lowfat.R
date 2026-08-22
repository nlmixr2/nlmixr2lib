Mauro_2025_nilotinib_lowfat <- function() {
  description <- "Two-compartment population PK model for nilotinib (Danziten reduced-dose tablets and Tasigna capsules) in healthy adults, refitted with the modified-fasted low-fat-meal study added, so all four published prandial states (fasted, high-fat meal 0.5 h pre-dose, high-fat meal 2 h pre-dose, low-fat meal 2 h pre-dose) are covered"
  reference <- paste(
    "Mauro M, Radich J, Jain P, Sequeira D, Bellanti F, Douer D.",
    "Pharmacokinetic profile of novel reduced-dose Danziten (nilotinib tablets)",
    "versus Tasigna (nilotinib capsules): in vivo bioequivalence and population",
    "pharmacokinetic analysis. Cancer Chemother Pharmacol. 2025;95(1):56.",
    "doi:10.1007/s00280-025-04777-6.",
    "Parameters are the 'Updated Model' of Supplemental Table 3, which adds study",
    "052-23 (modified fasting with a low-fat meal) to the 13 studies of the Final",
    "Model and was used for the published low-fat-meal steady-state simulations",
    "(main text Table 2). The companion Final Model that drives the fasted",
    "bioequivalence simulations is modellib('Mauro_2025_nilotinib').",
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
        "convention. NOTE the published exponent is NEGATIVE (-0.670, 95% CI",
        "-0.854 to -0.486), i.e. heavier subjects have LOWER nilotinib clearance;",
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
        "Two separate effects. (1) On CL for every subject (-26.6%: female CL is",
        "lower). (2) On relative bioavailability, but ONLY for female subjects",
        "receiving nilotinib TABLETS under FASTED conditions (+11.1%); the",
        "Supplemental Table 3 footnote reads 'Only female receiving nilotinib",
        "Tablets under fasted conditions', and the supplement text states no",
        "male/female absorption difference was seen under fed or modified fasted",
        "conditions. The gating is implemented as SEXF * FORM_TABLET * (1 - FED).",
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
        "bioavailability (+99.2% at the tablet reference dose) and on the zero-order",
        "absorption duration D1 (-43.2%), and selects between the formulation-",
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
        "dose effect for tablets above 200 mg (implemented with min(DOSE, 200)).",
        "(2) Selection among the dose-specific tablet high-fat-meal effects on",
        "relative bioavailability. The LOW-fat-meal effects are NOT dose-specific:",
        "study 052-23 tested one dose per formulation (190 mg tablets, 400 mg",
        "capsules), so a single effect per formulation is applied at every dose.",
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
        "Under the low-fat 2 h condition the absorption parameters D1 and ALAG1 are",
        "governed by the study-052-23 estimates instead (tablet D1, capsule ALAG1)",
        "or revert to their fasted values (capsule D1, tablet ALAG1); the V1 and Q",
        "food effects, carried over unchanged from the Final Model, still apply.",
        sep = " "
      ),
      source_name       = "FOOD"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal indicator for the dose record, 1 = a high-fat meal preceded the dose",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (fasted or low-fat meal)",
      notes             = paste(
        "Mauro 2025 Methods, 'Prandial state definitions': the high-fat meal is a",
        "high-fat high-calorie breakfast of 800-1000 kcal, non-vegetarian, about 50%",
        "fat content, consumed within 30 min. Paired with MEAL_PREDOSE_2H, which",
        "says whether the dose followed the START of that meal by about 0.5 h (the",
        "paper's 'fed' state) or by about 2 h (the paper's 'modified fasting with a",
        "high-fat meal' state).",
        sep = " "
      ),
      source_name       = "FOOD"
    ),
    FED_LOWFAT = list(
      description        = "Low-fat-meal indicator for the dose record, 1 = a low-fat meal preceded the dose",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (fasted or high-fat meal)",
      notes             = paste(
        "Mauro 2025 Methods, 'Prandial state definitions': the low-fat meal is a",
        "light-fat high-calorie breakfast of 400-500 kcal with about 25% fat",
        "content. In this paper the low-fat meal was only ever given 2 h before the",
        "dose (study 052-23, 'modified fasting with a low-fat meal'), so the",
        "low-fat effects are applied as FED_LOWFAT * MEAL_PREDOSE_2H; a low-fat",
        "meal 0.5 h before the dose was not studied and carries no estimate.",
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
        "on tablet bioavailability differ (at 142 mg: +46.7% at 0.5 h vs +57.1% at",
        "2 h; at 190 mg: +60.4% vs +59.0%). NOTE that the phrase 'modified fasting'",
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
        "(Supplemental Table 3, 'Inter-Occasion Variability'). Mauro 2025 does not",
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
    n_studies = 14,
    age_range = "adults (not reported)",
    weight_range = "not reported",
    sex_female_pct = 7.8,
    disease_state = "healthy",
    dose_range = "50-300 mg single oral dose (nilotinib tablets); 50-400 mg single oral dose (nilotinib capsules)",
    regions = "India (Phase 1 units, ZenRise, Hyderabad); analysis by GK Analytics, Hyderabad and Certara",
    notes = paste(
      "Healthy adult men and women in all 14 single-dose relative-bioavailability",
      "and food-effect crossover studies of Supplemental Table 1, including study",
      "052-23 (190 mg tablets vs 400 mg capsules under modified fasting with a",
      "low-fat meal, 90 enrolled / 86 completed) which the companion Final Model",
      "does not use. The paper reports the pooled analysis dataset as 30,089",
      "nilotinib observations from 502 subjects with 1,311 below-limit-of-",
      "quantification records excluded. 39 female participants contributed, from",
      "studies 006-23, 007-23 and 128-22. Baseline demographic distributions (age,",
      "weight, laboratory values) are not tabulated in the paper or its",
      "Supplemental Information. The published steady-state simulations resampled",
      "50 subjects (50% male / 50% female) from the analysis dataset.",
      sep = " "
    )
  )

  ini({
    # ---- Disposition -- Supplemental Table 3, 'Typical Values' ----
    lcl <- log(34.3); label("Clearance (CL, L/h)") # Suppl. Table 3: CL = 34.3 L/h (RSE 3.18%, 95% CI 32.1-36.4)
    lvc <- log(551); label("Central volume of distribution (V1, L)") # Suppl. Table 3: V1 = 551 L (RSE 2.71%, 95% CI 522-581)
    lq <- log(5.56); label("Inter-compartmental clearance (Q, L/h)") # Suppl. Table 3: Q = 5.56 L/h (RSE 19.1%, 95% CI 3.48-7.63)
    lvp <- log(217); label("Peripheral volume of distribution (V2, L)") # Suppl. Table 3: V2 = 217 L (RSE 6.74%, 95% CI 188-246)

    # ---- Primary absorption: lagged zero-order input into central ----
    ld1 <- log(2.41); label("Duration of the zero-order absorption input into central (D1, h)") # Suppl. Table 3: D1 = 2.41 h (RSE 1.61%, 95% CI 2.33-2.48)
    ltlag <- log(0.0866); label("Absorption lag time of the zero-order input (ALAG1, h)") # Suppl. Table 3: ALAG1 = 0.0866 h (RSE 21.0%, 95% CI 0.0510-0.122)
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the capsule reference at 400 mg fasted (F1, fraction)") # Suppl. Info F1 equations: F1_capsule = 1 * (DOSE/400)^theta10, i.e. F1 is anchored at 1 for the capsule reference dose

    # ---- Secondary depot: delayed late re-absorption ----
    lka2 <- log(0.141); label("First-order absorption rate constant from the secondary depot (1/h)") # Suppl. Table 3: absorption from secondary depot = 0.141 1/h (RSE 6.16%, 95% CI 0.124-0.158)
    ltlag2 <- log(4.72); label("Absorption lag time of the secondary depot (h)") # Suppl. Table 3: absorption lag re-absorbed from secondary depot = 4.72 h (RSE 5.58%, 95% CI 4.21-5.24)
    lfdepot2 <- log(1 - 0.890); label("Fraction of F1 that is re-absorbed through the secondary depot (fraction)") # Suppl. Table 3: fraction re-absorbed from secondary depot (proportional) = -0.890 (RSE 0.681%, 95% CI -0.901 to -0.878), i.e. 1 - 0.890 = 0.110 of F1

    # ---- Covariate effects on relative bioavailability -- Suppl. Table 3 ----
    e_form_tablet_fdepot <- 0.992; label("Proportional effect of the tablet formulation on F1 (unitless)") # Suppl. Table 3: formulation (nilotinib tablets) effect on F1 (proportional) = 0.992 (RSE 3.84%, 95% CI 0.917-1.07)
    e_dose_fdepot <- -0.698; label("Power exponent of dose on F1, shared by both formulations (unitless)") # Suppl. Table 3: dose effect on F1 (power) = -0.698 (RSE 6.79%, 95% CI -0.791 to -0.605)
    e_sexf_fdepot <- 0.111; label("Proportional effect of female sex on F1 for tablets under fasted conditions (unitless)") # Suppl. Table 3: sex (female*) covariate effect on F1 (proportional) = 0.111 (RSE 38.3%, 95% CI 0.0277-0.195); footnote * = "Only female receiving nilotinib Tablets under fasted conditions"

    # High-fat meal taken about 0.5 h before the dose (the paper's 'fed' state).
    e_hfmeal05h_tab50_fdepot <- 0.155; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 50 mg tablets (unitless)") # Suppl. Table 3: food (fed) effect on F1 at 50 mg for nilotinib tablets = 0.155 (RSE 36.1%, 95% CI 0.0451-0.264)
    e_hfmeal05h_tab142_fdepot <- 0.467; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 142 mg tablets (unitless)") # Suppl. Table 3: food (fed) effect on F1 at 142 mg for nilotinib tablets = 0.467 (RSE 11.6%, 95% CI 0.361-0.573)
    e_hfmeal05h_tab190_fdepot <- 0.604; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 190 mg tablets (unitless)") # Suppl. Table 3: food (fed) effect on F1 at 190 mg for nilotinib tablets = 0.604 (RSE 13.2%, 95% CI 0.448-0.760)
    e_hfmeal05h_tab240_fdepot <- 0.733; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 240 mg tablets (unitless)") # Suppl. Table 3: food (fed) effect on F1 at 240 mg for nilotinib tablets = 0.733 (RSE 10.0%, 95% CI 0.589-0.877)
    e_hfmeal05h_tab300_fdepot <- 0.610; label("Proportional effect of a high-fat meal 0.5 h pre-dose on F1, 300 mg tablets (unitless)") # Suppl. Table 3: food (fed) effect on F1 at 300 mg for nilotinib tablets = 0.610 (RSE 12.0%, 95% CI 0.467-0.753)

    # High-fat meal taken about 2 h before the dose (the paper's 'modified fasting
    # with a high-fat meal' state).
    e_hfmeal2h_tab142_fdepot <- 0.571; label("Proportional effect of a high-fat meal 2 h pre-dose on F1, 142 mg tablets (unitless)") # Suppl. Table 3: food (modified fasted) effect on F1 at 142 mg for nilotinib tablets = 0.571 (RSE 10.1%, 95% CI 0.457-0.684)
    e_hfmeal2h_tab190_fdepot <- 0.590; label("Proportional effect of a high-fat meal 2 h pre-dose on F1, 190 mg tablets (unitless)") # Suppl. Table 3: food (modified fasted) effect on F1 at 190 mg for nilotinib tablets = 0.590 (RSE 12.0%, 95% CI 0.452-0.729); this is the row that Suppl. Table 2 mislabels as "for nilotinib capsules"
    e_hfmeal_cap_fdepot <- 1.79; label("Proportional effect of a high-fat meal on F1, capsules, any dose-to-meal interval (unitless)") # Suppl. Table 3: food (modified fasted) effect on F1 for nilotinib capsules = 1.79 (RSE 7.15%, 95% CI 1.54-2.04); capsules were dosed under BOTH the 0.5 h and 2 h high-fat conditions (Suppl. Table 1) yet only one capsule high-fat F1 effect was estimated, so it is shared across the two timings

    # Low-fat meal taken about 2 h before the dose (study 052-23). Not dose-specific:
    # the study gave 190 mg tablets and 400 mg capsules only.
    e_lfmeal2h_tab_fdepot <- 0.294; label("Proportional effect of a low-fat meal 2 h pre-dose on F1, tablets (unitless)") # Suppl. Table 3: food (modified fasted low-fat meal) effect on F1 for nilotinib (proportional) = 0.294 (RSE 20.6%, 95% CI 0.175-0.413); Suppl. Info text: "The resulting increase in relative bioavailability as compared to fasted state was 29.4% and 51.2% for nilotinib tablets and capsules respectively", so the unqualified row is the TABLET estimate
    e_lfmeal2h_cap_fdepot <- 0.512; label("Proportional effect of a low-fat meal 2 h pre-dose on F1, capsules (unitless)") # Suppl. Table 3: food (modified fasted low-fat meal) effect on F1 for nilotinib capsules (proportional) = 0.512 (RSE 16.4%, 95% CI 0.347-0.677); matches the 51.2% in the Suppl. Info text

    # ---- Covariate effects on absorption and distribution -- Suppl. Table 3 ----
    e_form_tablet_d1 <- -0.432; label("Proportional effect of the tablet formulation on D1 (unitless)") # Suppl. Table 3: formulation (nilotinib tablets) effect on D1 (proportional) = -0.432 (RSE 3.03%, 95% CI -0.458 to -0.406)
    e_fed_tab_d1 <- 2.95; label("Proportional effect of food on D1, tablets (unitless)") # Suppl. Table 3: food effect on D1 for nilotinib tablets (proportional) = 2.95 (RSE 3.80%, 95% CI 2.73-3.17)
    e_fed_cap_d1 <- 0.706; label("Proportional effect of food on D1, capsules (unitless)") # Suppl. Table 3: food effect on D1 for nilotinib capsules (proportional) = 0.706 (RSE 9.48%, 95% CI 0.575-0.837)
    e_lfmeal2h_tab_d1 <- 0.948; label("Proportional effect of a low-fat meal 2 h pre-dose on D1, tablets (unitless)") # Suppl. Table 3: food (modified fasted low-fat meal) effect on D1 for nilotinib (proportional) = 0.948 (RSE 7.20%, 95% CI 0.814-1.08); Suppl. Info text: "the inclusion of an effect on D1 for nilotinib tablets and on lag-time for capsules allowed improving the overall description of the data", so the D1 estimate is the TABLET one
    e_fed_tab_tlag <- 3.62; label("Proportional effect of food on ALAG1, tablets (unitless)") # Suppl. Table 3: food effect on ALAG1 for nilotinib tablets (proportional) = 3.62 (RSE 26.7%, 95% CI 1.72-5.52)
    e_fed_cap_tlag <- 9.10; label("Proportional effect of food on ALAG1, capsules (unitless)") # Suppl. Table 3: food effect on ALAG1 for nilotinib capsules (proportional) = 9.10 (RSE 23.4%, 95% CI 4.93-13.3)
    e_lfmeal2h_cap_tlag <- 3.25; label("Proportional effect of a low-fat meal 2 h pre-dose on ALAG1, capsules (unitless)") # Suppl. Table 3: food (modified fasted low-fat meal) effect on ALAG1 for nilotinib capsules (proportional) = 3.25 (RSE 27.5%, 95% CI 1.50-5.00)
    e_fed_vc <- -0.286; label("Proportional effect of food on V1 (unitless)") # Suppl. Table 3: food effect on V1 (proportional) = -0.286 (RSE 10.3%, 95% CI -0.343 to -0.228)
    e_fed_q <- 2.42; label("Proportional effect of food on Q (unitless)") # Suppl. Table 3: food effect on Q (proportional) = 2.42 (RSE 25.5%, 95% CI 1.21-3.64)

    # ---- Covariate effects on clearance -- Suppl. Table 3 ----
    e_wt_cl <- -0.670; label("Power exponent of body weight on CL (unitless)") # Suppl. Table 3: bodyweight covariate effect on CL (power) = -0.670 (RSE 14.0%, 95% CI -0.854 to -0.486)
    e_sexf_cl <- -0.266; label("Proportional effect of female sex on CL (unitless)") # Suppl. Table 3: sex (female) covariate effect on CL (proportional) = -0.266 (RSE 11.4%, 95% CI -0.325 to -0.206)

    # ---- Between-subject variability -- Suppl. Table 3 ----
    # Log-scale variances; correlation = 0.0621 / sqrt(0.108 * 0.146) = 0.494.
    etalcl + etalvc ~ c(0.108, 0.0621, 0.146) # Suppl. Table 3: BSV CL = 0.108 (RSE 7.76%, shrinkage 2.60%), correlation CL-V1 = 0.0621 (RSE 11.4%), BSV V1 = 0.146 (RSE 8.22%, shrinkage 1.80%)

    # ---- Inter-occasion variability on V1 -- Suppl. Table 3 ----
    etaiov_vc_1 ~ 0.157 # Suppl. Table 3: IOV V1 = 0.157 (RSE 7.96%, shrinkage 5.80%)
    etaiov_vc_2 ~ fixed(0.157) # equal to occasion 1 per $OMEGA BLOCK(1) SAME
    etaiov_vc_3 ~ fixed(0.157) # equal to occasion 1 per $OMEGA BLOCK(1) SAME

    # ---- Residual error -- Suppl. Table 3 ----
    # Suppl. Info: "A combined (additive and proportional) error model was used";
    # Supplemental Table 3 labels the proportional term "Exponential Error".
    propSd <- 0.232; label("Proportional residual error (fraction)") # Suppl. Table 3: exponential error = 0.232 (RSE 2.45%, 95% CI 0.221-0.243)
    addSd <- 37.1; label("Additive residual error (ng/mL)") # Suppl. Table 3: additive error = 37.1 ng/mL (RSE 6.22%, 95% CI 32.6-41.6)
  })

  model({
    # ---- 1. Prandial-state and formulation strata ----
    # Mauro 2025 Methods, 'Prandial state definitions', defines four mutually
    # exclusive states. FED_HIGHFAT / FED_LOWFAT carry the meal composition and
    # MEAL_PREDOSE_2H the dose-to-meal interval:
    #   fasted                        FED = 0
    #   high-fat meal 0.5 h pre-dose  FED = 1, FED_HIGHFAT = 1, MEAL_PREDOSE_2H = 0
    #   high-fat meal 2 h pre-dose    FED = 1, FED_HIGHFAT = 1, MEAL_PREDOSE_2H = 1
    #   low-fat meal 2 h pre-dose     FED = 1, FED_LOWFAT  = 1, MEAL_PREDOSE_2H = 1
    tab <- FORM_TABLET
    cap <- 1 - FORM_TABLET
    hf05h <- FED_HIGHFAT * (1 - MEAL_PREDOSE_2H)
    hf2h <- FED_HIGHFAT * MEAL_PREDOSE_2H
    lf2h <- FED_LOWFAT * MEAL_PREDOSE_2H

    # ---- 2. Relative bioavailability ----
    # Supplemental Information, 'Initial Model Development':
    #   F1_capsule = 1 * (DOSE/400)^theta10
    #   F1_tablet  = 1 * (1 + theta6) * (DOSE/200)^theta10   for DOSE <= 200 mg
    #   F1_tablet  = 1 * (1 + theta6)                        for DOSE >  200 mg
    doseratio <- tab * min(DOSE, 200) / 200 + cap * DOSE / 400
    f1_fasted <- exp(lfdepot) * (1 + e_form_tablet_fdepot * tab) * doseratio^e_dose_fdepot

    # Dose-specific tablet high-fat effects. Mauro 2025 estimated them only at 50,
    # 142, 190, 240 and 300 mg (0.5 h pre-dose) and at 142 and 190 mg (2 h
    # pre-dose), so the step boundaries below are midpoints between the published
    # anchors and are this encoding's choice, not the paper's. The effects are
    # non-monotonic in dose (0.604 at 190 mg vs 0.610 at 300 mg), so interpolating
    # between anchors would not be justified.
    fed05h_tab <- e_hfmeal05h_tab50_fdepot * (DOSE < 96) +
      e_hfmeal05h_tab142_fdepot * (DOSE >= 96) * (DOSE < 166) +
      e_hfmeal05h_tab190_fdepot * (DOSE >= 166) * (DOSE < 215) +
      e_hfmeal05h_tab240_fdepot * (DOSE >= 215) * (DOSE < 270) +
      e_hfmeal05h_tab300_fdepot * (DOSE >= 270)
    fed2h_tab <- e_hfmeal2h_tab142_fdepot * (DOSE < 166) +
      e_hfmeal2h_tab190_fdepot * (DOSE >= 166)

    foodf1 <- hf05h * (tab * fed05h_tab + cap * e_hfmeal_cap_fdepot) +
      hf2h * (tab * fed2h_tab + cap * e_hfmeal_cap_fdepot) +
      lf2h * (tab * e_lfmeal2h_tab_fdepot + cap * e_lfmeal2h_cap_fdepot)

    # The female bioavailability effect applies only to tablets under fasted
    # conditions (Suppl. Table 3 footnote).
    f1 <- f1_fasted * (1 + foodf1) * (1 + e_sexf_fdepot * SEXF * tab * (1 - FED))

    # ---- 3. Individual parameters ----
    iov_vc <- (OCC == 1) * etaiov_vc_1 + (OCC == 2) * etaiov_vc_2 +
      (OCC == 3) * etaiov_vc_3

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_sexf_cl * SEXF)
    vc <- exp(lvc + etalvc + iov_vc) * (1 + e_fed_vc * FED)
    q <- exp(lq) * (1 + e_fed_q * FED)
    vp <- exp(lvp)

    # The low-fat 2 h condition got its absorption re-described from scratch:
    # Supplemental Information, 'Updated model with modified fasted low-fat meal
    # data' -- "First, effects of the new food conditions were evaluated on the
    # bioavailability of the two products relative to their reference bioavailability
    # under fasted conditions. Secondly, differences on other absorption parameters
    # such as D1 and lag-time were evaluated", of which only a tablet D1 effect and
    # a capsule lag-time effect were retained. So under low-fat each absorption
    # parameter takes EITHER its own low-fat estimate (tablet D1, capsule ALAG1) OR
    # its fasted value (capsule D1, tablet ALAG1) -- a non-retained effect means no
    # change from the fasted reference, not the general fed value. The food effects
    # on V1 and Q are structural effects carried over unchanged from the Final Model
    # ("while maintaining the previously identified structural and covariate
    # effects") and so still apply. Checked against main-text Table 2: this reading
    # reproduces the published low-fat Cmax closely, whereas letting the general
    # D1 / ALAG1 food effects persist under the low-fat condition roughly triples
    # the capsule Cmax shortfall. The vignette carries the side-by-side numbers.
    d1_tab_food <- lf2h * e_lfmeal2h_tab_d1 + (1 - lf2h) * FED * e_fed_tab_d1
    d1_cap_food <- (1 - lf2h) * FED * e_fed_cap_d1
    tlag_tab_food <- (1 - lf2h) * FED * e_fed_tab_tlag
    tlag_cap_food <- lf2h * e_lfmeal2h_cap_tlag + (1 - lf2h) * FED * e_fed_cap_tlag

    d1 <- exp(ld1) * (1 + e_form_tablet_d1 * tab) *
      (1 + tab * d1_tab_food + cap * d1_cap_food)
    tlag <- exp(ltlag) * (1 + tab * tlag_tab_food + cap * tlag_cap_food)
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
    f(central) <- f1
    dur(central) <- d1
    alag(central) <- tlag
    f(depot2) <- f1 * fdepot2
    alag(depot2) <- tlag2

    # ---- 7. Observation and residual error ----
    # Dose is in mg and vc in L, so central/vc is mg/L; the factor 1000 converts to
    # the ng/mL scale of Supplemental Table 3's additive error and of main-text
    # Tables 1 and 2.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
