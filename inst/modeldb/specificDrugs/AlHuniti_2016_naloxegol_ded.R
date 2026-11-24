AlHuniti_2016_naloxegol_ded <- function() {
  description <- "Time-to-event model of diary entry discontinuation (DED, the day after the last non-missing eDiary entry) over the 12-week double-blind treatment period in adults with noncancer pain and opioid-induced constipation randomized to placebo or oral naloxegol 12.5 or 25 mg once daily (Al-Huniti 2016). Log-linear hazard in time with a separate baseline log-hazard and time slope for each randomized arm, plus an additional steep slope from day 80 that reproduces final study visits held up to 4 days before day 84. No between-subject variability; no covariate effects were retained."
  reference <- paste(
    "Al-Huniti N, Nielsen JC, Hutmacher MM, Lappalainen J, Cantagallo K, Sostek M.",
    "Population Exposure-Response Modeling of Naloxegol in Patients With",
    "Noncancer-Related Pain and Opioid-Induced Constipation.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(7):359-366.",
    "doi:10.1002/psp4.12099"
  )
  vignette <- "AlHuniti_2016_naloxegol"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "(no concentration output; sur is the probability of continued diary entry, 0-1)"
  )

  covariateData <- list(
    DOSE = list(
      description = "Randomized naloxegol daily dose (nominal arm dose)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a) of the DOSE canonical: the randomized nominal dose, constant over the treatment period (source column NDOSE). Only the three levels 0 (placebo), 12.5 and 25 mg/day are defined; each selects its own baseline log-hazard and time slope. Any other value falls back to the placebo parameters.",
      source_name = "NDOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on the baseline log-hazard and the time slope (Supplementary Table S2); the WAM procedure selected no covariate.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on the baseline log-hazard and the time slope; not retained.",
      source_name = "SEXF"
    ),
    OPIOID_WEAK_ONLY = list(
      description = "Baseline maintenance opioid regimen consisted of weak opioids only",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (strong opioid, or weak plus strong opioids, at baseline)",
      notes = "Baseline opioid potency screened on the baseline log-hazard and the time slope; not retained.",
      source_name = "BOTP"
    )
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "diary-entry-discontinuation hazard",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1331L,
    n_studies = 2L,
    age_range = "Adults; mean age 51.9-52.9 years by study arm (Table 1)",
    sex_female_pct = 62,
    race_ethnicity = "About 79% White, 19% Black, 1% Asian, 2% other (Table 1)",
    disease_state = "Chronic noncancer pain treated with stable maintenance opioids and opioid-induced constipation.",
    dose_range = "Naloxegol 12.5 or 25 mg orally once daily, or placebo, for 12 weeks.",
    regions = "Multinational (phase III KODIAC 4 and KODIAC 5).",
    notes = "Day 0 to day 84 of the double-blind treatment period. 543 patients (40.8%) made their last diary entry before day 84."
  )

  ini({
    # Supplementary Table S6 (final DED model, control stream run 3).
    lhaz_placebo <- -6.28
    label("Baseline log-hazard of diary entry discontinuation, placebo arm (log 1/day)") # Table S6: theta base-pbo = -6.28
    lhaz_nal12p5 <- -6.18
    label("Baseline log-hazard of diary entry discontinuation, naloxegol 12.5 mg arm (log 1/day)") # Table S6: theta base-12.5 = -6.18
    lhaz_nal25 <- -5.54
    label("Baseline log-hazard of diary entry discontinuation, naloxegol 25 mg arm (log 1/day)") # Table S6: theta base-25 = -5.54
    slope_haz_placebo <- 0.00970
    label("Change in log-hazard per day of treatment, placebo arm (1/day)") # Table S6: theta slp-pbo = 0.00970
    slope_haz_nal12p5 <- 0.00818
    label("Change in log-hazard per day of treatment, naloxegol 12.5 mg arm (1/day)") # Table S6: theta slp-12.5 = 0.00818
    slope_haz_nal25 <- -0.00114
    label("Change in log-hazard per day of treatment, naloxegol 25 mg arm (1/day)") # Table S6: theta slp-25 = -0.00114 (Table S5 run 3 prints -0.00143)
    slope_haz_end <- 1.08
    label("Additional change in log-hazard per day from day 80 onward, all arms (1/day)") # Table S6: theta end = 1.08

    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on the survival probability (unitless)") # not from the source; see vignette Assumptions and deviations
  })

  model({
    nal12p5 <- 1 * (DOSE == 12.5)
    nal25 <- 1 * (DOSE == 25)
    pbo <- 1 - nal12p5 - nal25

    lhaz0 <- pbo * lhaz_placebo + nal12p5 * lhaz_nal12p5 + nal25 * lhaz_nal25
    slope <- pbo * slope_haz_placebo + nal12p5 * slope_haz_nal12p5 + nal25 * slope_haz_nal25

    # The end-of-study slope switches on at day 80 (I_T80 = 1 for t >= 80).
    end80 <- 1 * (time >= 80)
    hazard <- exp(lhaz0 + slope * time + end80 * slope_haz_end * (time - 80))

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
