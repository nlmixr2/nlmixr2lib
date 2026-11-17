Lohmer_2026_buprenorphine_dropout <- function() {
  description <- paste(
    "Gompertz time-to-event model of treatment discontinuation (dropout)",
    "in adults with moderate-to-severe opioid use disorder receiving",
    "monthly subcutaneous extended-release buprenorphine (BUP-XR,",
    "SUBLOCADE) or placebo, integrated across a pre-fentanyl phase 3 study",
    "(Study 1, NCT02357901) and a high-risk fentanyl-era study (Study 2,",
    "NCT04995029). hazard(t) = beta0 * beta1[craving category] *",
    "(AGE/40)^beta2 * beta3^RACE_BLACK * exp(-ke * t), with separate",
    "beta0, ke and age exponent for the Study 1 placebo arm. The craving",
    "category is the last observed Opioid Craving VAS category carried",
    "forward (0, 1-5, 6-20, >20 mm), supplied as a time-varying covariate;",
    "craving >20 mm doubles the hazard. The same parameters describe both",
    "studies once craving is accounted for. Survival = exp(-cumhaz)."
  )
  reference <- paste(
    "Lohmer LRL, Greenwald MK, Kakar R, Laffont CM. Optimal",
    "buprenorphine exposure to treat patients with opioid use disorder in",
    "the era of fentanyl and polysubstance use. Front Pharmacol.",
    "2026;17:1839067. doi:10.3389/fphar.2026.1839067.",
    "Parameter estimates from Table 4; model equations from the",
    "Supplementary Material ('Treatment Retention Model Equations').",
    "Sister models from the same paper:",
    "modellib('Lohmer_2026_buprenorphine_opioidUse'),",
    "modellib('Lohmer_2026_buprenorphine_craving')."
  )
  vignette <- "Lohmer_2026_buprenorphine"

  units <- list(
    time = "day",
    dosing = "N/A (time-to-event model; no dosing compartment)",
    concentration = "N/A (survival probability output)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of treatment discontinuation",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    OCVAS_CAT = list(
      description = "Last observed current-craving category on the 0-100 mm Opioid Craving VAS, carried forward: 1 = 0 mm, 2 = 1-5 mm, 3 = 6-20 mm, 4 = >20 mm.",
      units = "(category 1-4)",
      type = "categorical",
      reference_category = "1 (0 mm, no craving)",
      notes = "Time-varying, last observation carried forward from the start of the trial (Methods, 'Modeling efficacy outcomes'). Same coding as the `ocvas_cat` endpoint of Lohmer_2026_buprenorphine_craving, so a simulated craving category can be fed in directly. Multiplies the hazard by 0.897, 1.04 and 2.28 for categories 2, 3 and 4 (Table 4).",
      source_name = "Craving LOCF"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model (AGE/40)^beta2 with a separate exponent for placebo (Supplementary Material; Table 4).",
      source_name = "Age"
    ),
    RACE_BLACK = list(
      description = "Black / African American race (1 = yes, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or other race)",
      notes = "Multiplies the hazard by 0.663 ('34% lower in Black/African American participants', Results).",
      source_name = "RACE"
    ),
    PLACEBO = list(
      description = "Placebo-treated subject (Study 1 placebo arm) = 1; BUP-XR-treated subject (either study) = 0.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BUP-XR)",
      notes = "Selects the placebo baseline hazard, hazard-decline rate and age exponent (Table 4).",
      source_name = "Placebo / BUP-XR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1071L,
    n_studies = 2L,
    age_range = "18-72 years (Study 1 median 38, range 19-64; Study 2 median 40, range 18-72)",
    age_mean = "Study 1 39.7 years (SD 11.1); Study 2 41.2 years (SD 10.6)",
    disease_state = "Treatment-seeking adults with moderate-to-severe opioid use disorder (DSM-5); Study 2 additionally required high-risk opioid use.",
    dose_range = "Study 1: BUP-XR 300/300 mg, 300/100 mg, or placebo every 4 weeks for 24 weeks. Study 2: 300 mg at weeks 1 and 2, then 100 or 300 mg every 4 weeks to week 38.",
    regions = "USA (Study 1); USA and Canada (Study 2)",
    notes = "Dropout dataset of 1,071 participants (489 Study 1, 582 Study 2 including participants who discontinued before the week-6 randomisation), Results 'Dropout model'. Demographics in Table 1 are for the 989-participant efficacy population."
  )

  ini({
    # Baseline hazard and its exponential decline (Gompertz), by treatment
    lhaz0_drug <- log(0.00332); label("Baseline dropout hazard beta0 with BUP-XR (1/day)") # Table 4 'beta0 BUP-XR' = 0.00332 (RSE 12%)
    lhaz0_placebo <- log(0.00860); label("Baseline dropout hazard beta0 with placebo, Study 1 (1/day)") # Table 4 'beta0 Placebo' = 0.00860 (RSE 23%)
    lkhaz_drug <- log(0.00271); label("Rate of exponential decrease of the hazard with BUP-XR (1/day)") # Table 4 'ke BUP-XR' = 0.00271 1/day (RSE 27%)
    lkhaz_placebo <- log(0.00906); label("Rate of exponential decrease of the hazard with placebo (1/day)") # Table 4 'ke Placebo' = 0.00906 1/day (RSE 36%)

    # Covariate effects on the hazard
    e_ocvas_cat2_haz <- 0.897; label("Hazard multiplier for last-observed craving 1-5 mm vs 0 mm (unitless)") # Table 4 'beta1 craving=1-5' = 0.897 (RSE 13%)
    e_ocvas_cat3_haz <- 1.04; label("Hazard multiplier for last-observed craving 6-20 mm vs 0 mm (unitless)") # Table 4 'beta1 craving=6-20' = 1.04 (RSE 15%)
    e_ocvas_cat4_haz <- 2.28; label("Hazard multiplier for last-observed craving >20 mm vs 0 mm (unitless)") # Table 4 'beta1 craving >20' = 2.28 (RSE 12%)
    e_age_haz_drug <- -0.581; label("Power exponent of AGE/40 on the hazard with BUP-XR (unitless)") # Table 4 'beta2 age, BUP-XR' = -0.581 (RSE 32%)
    e_age_haz_placebo <- -1.70; label("Power exponent of AGE/40 on the hazard with placebo (unitless)") # Table 4 'beta2 age, Placebo' = -1.70 (RSE 27%)
    e_race_black_haz <- 0.663; label("Hazard multiplier for Black/African American participants (unitless)") # Table 4 'beta3 race' = 0.663 (RSE 13%)
  })

  model({
    haz0 <- exp(lhaz0_drug) * (1 - PLACEBO) + exp(lhaz0_placebo) * PLACEBO
    khaz <- exp(lkhaz_drug) * (1 - PLACEBO) + exp(lkhaz_placebo) * PLACEBO
    e_age <- e_age_haz_drug * (1 - PLACEBO) + e_age_haz_placebo * PLACEBO

    # beta1^(Craving LOCF): one multiplier per category, category 1 = 1
    hr_craving <- 1
    if (OCVAS_CAT == 2) hr_craving <- e_ocvas_cat2_haz
    if (OCVAS_CAT == 3) hr_craving <- e_ocvas_cat3_haz
    if (OCVAS_CAT == 4) hr_craving <- e_ocvas_cat4_haz

    hazard <- haz0 * hr_craving * (AGE / 40)^e_age *
      e_race_black_haz^RACE_BLACK * exp(-khaz * t)
    d/dt(cumhaz) <- hazard
    sur <- exp(-cumhaz)
  })
}
