dePillis_2025_covid19Vaccine_lme <- function() {
  description <- paste(
    "Linear mixed-effects model of the decay of SARS-CoV-2 neutralizing",
    "antibody (NAb) activity from day 37 after the first dose of a two-dose",
    "mRNA COVID-19 vaccine primary series (BNT162b2 or mRNA-1273) in",
    "COVID-19-naive and convalesced adults. logit(NAb) is a cubic polynomial",
    "in centred time since the first dose, tc = (t - 178)/365 years, with a",
    "correlated random intercept and random linear slope; prior COVID-19",
    "infection (convalesced) shifts the intercept and the quadratic term.",
    "Algebraic (no ODE); time is days since the first vaccine dose, and the",
    "model describes only the decay region (t >= 37 days) up to the first",
    "booster or breakthrough infection.",
    sep = " "
  )
  reference <- paste(
    "de Pillis L, Caffrey R, Chen G, Dela M, Eldevik L, Liu H, McConnell JP,",
    "Shabahang S, Varvel SA. (2025). Comparing neutralizing antibody activity",
    "over time between naive and convalesced COVID-19 vaccinated individuals.",
    "Sci Rep 15:34800. doi:10.1038/s41598-025-18673-9.",
    "Correction (author-name only): doi:10.1038/s41598-025-26404-3.",
    sep = " "
  )
  vignette <- "dePillis_2025_covid19Vaccine"
  units <- list(
    time = "day",
    dosing = "none",
    concentration = "fraction neutralization (0-1)"
  )

  covariateData <- list(
    SARS_SEROPOS = list(
      description = "Convalesced status: 1 = history of COVID-19 infection before the first vaccine dose (convalesced), 0 = COVID-19-naive",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper's indicator C (data column isConvalesced = Yes / No).",
        "Subjects without a known infection but with nucleocapsid antibodies",
        "in their pre-vaccine sample were excluded, so in this vaccine-naive",
        "population 1 means seropositive from prior infection and 0 means",
        "nucleocapsid-seronegative at baseline."
      ),
      source_name = "isConvalesced (C in Eq. 4)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_studies = 1L,
    age_range = "adults; mean (SD) age 53.6 (13.3) years naive, 53.7 (10.5) years convalesced",
    weight_range = "mean (SD) BMI 27.4 (5.8) kg/m^2 naive, 29.5 (4.2) kg/m^2 convalesced",
    sex_female_pct = 47.8,
    race_ethnicity = NA_character_,
    disease_state = "Adults vaccinated against COVID-19: 68 COVID-19-naive and 24 convalesced (prior COVID-19 infection)",
    dose_range = "Two-dose mRNA primary series: mRNA-1273 (Moderna) or BNT162b2 (Pfizer-BioNTech)",
    regions = "United States",
    notes = paste(
      "Aditxt prospective cohort NCT05379478. Only observations from day 37",
      "after the first vaccine dose until the first booster or breakthrough",
      "infection enter this model. Table 1 of de Pillis 2025."
    )
  )

  ini({
    # Fixed effects of the optimal model (Eqs. 3-4), on the logit NAb scale.
    # Time enters as tc = (t - tbar)/365 with t in days since the first dose.
    tbar <- fixed(178)
    label("Centring time, mean days since the first dose over the decay-region observations (day)") # Results: 'tbar approx 178 days'
    b00 <- -0.082
    label("Intercept of logit(NAb) at t = tbar, naive (unitless)") # Table 2 beta00 = -0.082 (SE 0.140)
    e_sars_seropos_b00 <- 1.535
    label("Convalesced shift in the intercept of logit(NAb), paper beta01 (unitless)") # Table 2 beta01 = 1.535 (SE 0.276)
    b10 <- -4.015
    label("Linear slope of logit(NAb) in centred time (1/year)") # Table 2 beta10 = -4.015 (SE 0.361)
    b20 <- 5.711
    label("Quadratic coefficient of logit(NAb) in centred time, naive (1/year^2)") # Table 2 beta20 = 5.711 (SE 0.757)
    e_sars_seropos_b20 <- 4.890
    label("Convalesced shift in the quadratic coefficient, paper beta21 (1/year^2)") # Table 2 beta21 = 4.890 (SE 0.929)
    b30 <- -5.715
    label("Cubic coefficient of logit(NAb) in centred time (1/year^3)") # Table 2 beta30 = -5.715 (SE 1.215)

    # Correlated random intercept b0i and random slope b1i (Eq. 4). Table 2
    # reports SDs and their correlation; covariance = rho * sd0 * sd1.
    etab00 + etab10 ~ c(
      1.032^2,
      -0.420 * 1.032 * 2.103,
      2.103^2
    ) # Table 2 sd(b0) = 1.032, sd(b1) = 2.103, rho(b1,b0) = -0.420

    addSd <- 0.327
    label("Residual error SD on the logit NAb scale (unitless)") # Table 2 'Model Error Estimate' sigma_epsilon = 0.327
  })

  model({
    tc <- (t - tbar) / 365

    beta0 <- b00 + e_sars_seropos_b00 * SARS_SEROPOS + etab00
    beta1 <- b10 + etab10
    beta2 <- b20 + e_sars_seropos_b20 * SARS_SEROPOS
    beta3 <- b30

    # Eq. 3: logit(y) = beta0 + beta1*tc + beta2*tc^2 + beta3*tc^3 + eps
    nab_activity <- expit(beta0 + beta1 * tc + beta2 * tc^2 + beta3 * tc^3)
    nab_activity ~ logitNorm(addSd, 0, 1)
  })
}
