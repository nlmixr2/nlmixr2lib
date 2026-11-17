Masters_2026_sasanlimab_efs <- function() {
  description <- paste(
    "Parametric log-normal time-to-event model for investigator-assessed",
    "event-free survival (EFS) in adults with BCG-naive high-risk non-muscle",
    "invasive bladder cancer treated with sasanlimab 300 mg SC every 4 weeks",
    "plus BCG (Masters 2026; CREST study B8011006 Cohort A, arms A and B,",
    "n = 698, 158 events). log(EFS time) is normal with location mu and scale",
    "sigma; mu is the typical location multiplied by covariate factors for age",
    "(linear), baseline body weight (power) and treatment arm (BCG induction",
    "only versus induction plus maintenance). Sasanlimab exposure was tested",
    "and was NOT retained, so the model carries no exposure term. No ODE",
    "state, no random effects and no residual error; it exposes the",
    "event-free probability `sur`, the cumulative hazard `cumhaz` and the",
    "median EFS time `tmed`. Time is in months.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Parameter estimates are in Supporting Information Table S4; the visual",
    "predictive check is Supporting Information Figure S4.",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(
    time = "month",
    dosing = "n/a (no dosing events; sasanlimab exposure was not retained in the final EFS model)",
    concentration = "n/a (outputs are an event-free probability, a cumulative hazard and a median EFS time in months)"
  )

  covariateData <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear multiplicative effect on the log-normal location, centred at 67 years: mu * (1 - 0.0108 * (AGE - 67)). The centring value is the median age of the B8011006 popPK participants (Table S1, n = 706); the EFS analysis set (n = 698) median is not reported. Older age shortens EFS.",
      source_name = "Age"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power multiplicative effect on the log-normal location, normalised to 76 kg: mu * (WT / 76)^-0.336. The normalising value is the median baseline weight of the B8011006 popPK participants (Table S1, n = 706). Higher weight shortens EFS.",
      source_name = "Body Weight"
    ),
    BCG_MAINT = list(
      description = "BCG maintenance indicator for the CREST Cohort A treatment arm: 1 = arm A (sasanlimab + BCG induction + BCG maintenance), 0 = arm B (sasanlimab + BCG induction only, no maintenance).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (arm A, induction plus maintenance)",
      notes = "Masters 2026 reports a 'Treatment Arm on mu' coefficient of -0.144 and states that receipt of BCG induction + maintenance is beneficial to EFS, so the negative coefficient applies to the arm WITHOUT maintenance: mu * (1 - 0.144 * (1 - BCG_MAINT)). Arm A n = 350, arm B n = 348 (Table S3).",
      source_name = "Treatment Arm"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 698L,
    n_studies = 1L,
    n_events = 158L,
    age_range = "31-90 years (B8011006 popPK participants, Table S1)",
    age_median = "67 years",
    weight_range = "40.5-173 kg",
    weight_median = "76 kg",
    sex_female_pct = 17.7,
    race_ethnicity = c(`Non-Asian` = 62.6, Asian = 35.1, Missing = 2.3),
    disease_state = "BCG-naive, high-risk non-muscle invasive bladder cancer; CREST trial (B8011006) Cohort A arm A (sasanlimab + BCG induction + maintenance, n = 350) and arm B (sasanlimab + BCG induction only, n = 348).",
    dose_range = "Sasanlimab 300 mg SC every 4 weeks, combined with BCG",
    notes = "158 investigator-assessed EFS events in 698 participants; median EFS not reached at data cutoff. The log-normal distribution fitted best among exponential, Weibull, Gompertz, log-normal and log-logistic. Fitted in NONMEM 7.5 with stepwise covariate modelling (forward alpha 0.05, backward alpha 0.001). Demographics (Tables S1, S2) are tabulated for the 706 B8011006 popPK participants, not for the 698 EFS participants."
  )

  ini({
    # Supporting Information Table S4. Time scale: months (Figure 4 and Figure
    # S4 axes); with these values the population-averaged event-free
    # probability reproduces the Figure S4 VPC band (see the vignette).
    ltmed0 <- 5.56; label("Typical location mu of the log-normal EFS distribution (log-months) at age 67 years, 76 kg, arm A") # Table S4: theta mu = 5.56 (RSE 4.0%)
    e_age_tmed <- -0.0108; label("Linear fractional effect of age on mu (per year)") # Table S4: theta Age on mu = -0.0108 (RSE 21.3%)
    e_wt_tmed <- -0.336; label("Power exponent of baseline body weight on mu (unitless)") # Table S4: theta Body Weight on mu = -0.336 (RSE 34.9%)
    e_bcgnomaint_tmed <- -0.144; label("Fractional change in mu for arm B (BCG induction only) vs arm A (unitless)") # Table S4: theta Treatment Arm on mu = -0.144 (RSE 28.0%)
    sdlogt <- 2.28; label("Scale sigma of the log-normal EFS distribution (SD of log EFS time)") # Table S4: theta sigma = 2.28 (RSE 4.6%)
  })

  model({
    # Covariates act multiplicatively on the location (stepwise covariate
    # modelling parameterisation): linear 1 + theta * (cov - median), power
    # (cov / median)^theta, categorical 1 + theta * indicator.
    mu_logt <- ltmed0 * (1 + e_age_tmed * (AGE - 67)) * (WT / 76)^e_wt_tmed *
      (1 + e_bcgnomaint_tmed * (1 - BCG_MAINT))

    # Median EFS, months
    tmed <- exp(mu_logt)

    # Event-free probability; the 1e-6 offset keeps log(t) finite at t = 0.
    sur <- 1 - pnorm((log(t + 1e-6) - mu_logt) / sdlogt)
    cumhaz <- -log(sur)
  })
}
