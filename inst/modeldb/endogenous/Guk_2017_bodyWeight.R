Guk_2017_bodyWeight <- function() {
  description <- paste(
    "Body-weight versus postmenstrual age (PMA) growth model for Korean",
    "children under 12 years of age (Guk 2017 Transl Clin Pharmacol). This",
    "is a covariate (anthropometry) model, not a drug model: it has no drug,",
    "no dosing and no ODEs, and the body-weight output is an algebraic",
    "function of the rxode2 time variable, which the model interprets as",
    "PMA in WEEKS (fitted range 24.6 to 718.3 weeks). Weight is a sigmoid",
    "(Hill) function of PMA multiplied by a logistic function of PMA,",
    "BW = WMAX * PMA^hill / (PMA50^hill + PMA^hill) /",
    "(exp(-alpha * PMA) + beta / alpha), with sex-specific WMAX, PMA50,",
    "Hill coefficient and logistic maximum alpha/beta and a shared logistic",
    "rate alpha. One weight per child, so no between-subject variability",
    "was estimated; residual error is proportional. SEXF is required.",
    "Intended for imputing missing pediatric weights from PMA and for",
    "simulating a time-varying weight covariate.",
    sep = " "
  )
  reference <- paste(
    "Guk J, Chae D, Park K.",
    "Relationship between body weight and postmenstrual age in a Korean",
    "pediatric population.",
    "Transl Clin Pharmacol. 2017;25(2):101-105.",
    "doi:10.12793/tcp.2017.25.2.101.",
    sep = " "
  )
  vignette <- "Guk_2017_bodyWeight"
  units <- list(
    time = "week (postmenstrual age, PMA = gestational age + postnatal age; fitted range 24.6 to 718.3)",
    dosing = "n/a (no exogenous dosing; anthropometric covariate model)",
    concentration = "n/a (no drug concentration). Output is body weight (kg)."
  )

  covariateData <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Time-fixed and REQUIRED. Guk 2017 reports a separate estimate for",
        "boys and for girls of WMAX, PMA50, the Hill coefficient and the",
        "logistic maximum alpha/beta; the logistic rate alpha is shared",
        "(Table 2; Results: 'gender affected all parameters except for the",
        "rate of change in body weight alpha'). The model selects the girls'",
        "value when SEXF is 1 and the boys' value when SEXF is 0. Cohort:",
        "47,030 boys (56.7%) and 35,984 girls (43.3%) (Table 1)."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83014L,
    n_studies = 1L,
    age_range = "Under 12 years of age; PMA 24.6 to 718.3 weeks (median 193.9, mean 251.0, SD 192.3; Table 1)",
    weight_range = "0.4 to 97.0 kg (median 14.0, mean 17.3, SD 11.9; Table 1)",
    sex_female_pct = 43.3,
    race_ethnicity = "Korean",
    disease_state = paste(
      "Pediatric patients treated at Severance Hospital (Seoul, Korea),",
      "January 2012 to March 2017; any diagnosis. 34,851 infants (42.0%) and",
      "48,163 children (58.0%) (Table 1)."
    ),
    dose_range = "n/a (no drug)",
    regions = "Korea (Severance Hospital, Seoul)",
    notes = paste(
      "Electronic-medical-record data with a single weight per patient",
      "(Methods, Data). PMA = gestational age + chronological age; where the",
      "gestational age was missing the child was assumed to be born at term",
      "(40 weeks). Patients with missing weight, PMA or sex were excluded.",
      "Fitted in NONMEM 7.3.0 with FOCE."
    )
  )

  ini({
    # Eq. (2): BW = WMAX * PMA^gamma / (PMA50^gamma + PMA^gamma)
    #               * 1 / (exp(-alpha * PMA) + beta / alpha)
    # Table 2 reports each sex-specific parameter as a separate row.
    lbwmax_male <- log(8.42)
    label("Maximum of the sigmoid component WMAX, boys (kg)") # Table 2, WMAX (kg), Boys = 8.42 (RSE 0.31%)
    lbwmax_female <- log(8.22)
    label("Maximum of the sigmoid component WMAX, girls (kg)") # Table 2, WMAX (kg), Girls = 8.22 (RSE 37.2%)

    # The Results text gives PMA50 to two decimals (44.13 and 44.88 weeks);
    # Table 2 prints 44.1 and 44.8, the same values truncated to one decimal.
    ltmat50_male <- log(44.13)
    label("PMA at half of the sigmoid maximum PMA50, boys (week)") # Results paragraph 2: 'PMA50 were estimated as 44.13 weeks for boys'; Table 2 row PMA50 (weeks), Boys = 44.1 (RSE 2.26%)
    ltmat50_female <- log(44.88)
    label("PMA at half of the sigmoid maximum PMA50, girls (week)") # Results paragraph 2: '44.88 weeks for girls'; Table 2 row PMA50 (weeks), Girls = 44.8 (RSE 22.7%)

    lhill_male <- log(4.53)
    label("Hill coefficient gamma of the sigmoid component, boys (unitless)") # Table 2, gamma, Boys = 4.53 (RSE 0.84%); Results paragraph 2 = 4.53
    lhill_female <- log(4.11)
    label("Hill coefficient gamma of the sigmoid component, girls (unitless)") # Table 2, gamma, Girls = 4.11 (RSE 11.9%); Results paragraph 2 = 4.11

    lkgrow <- log(0.003299)
    label("Rate constant alpha of the logistic component, both sexes (1/week)") # Table 2, alpha (/week) = 0.003299 (RSE 1.27%)

    # The logistic factor tends to alpha/beta as PMA grows; the paper
    # estimated beta/alpha and reported alpha/beta by the delta method.
    llogismax_male <- log(14.3)
    label("Maximum of the logistic component alpha/beta, boys (unitless)") # Table 2, alpha/beta, Boys = 14.3 (RSE 4.45%)
    llogismax_female <- log(12.7)
    label("Maximum of the logistic component alpha/beta, girls (unitless)") # Table 2, alpha/beta, Girls = 12.7 (RSE 32.2%)

    # No IIV: one weight per patient (Results paragraph 1).
    propSd <- 0.179
    label("Proportional residual error on body weight (fraction)") # Table 2, Residual variability, sigma^2 (CV%) = 17.9 (RSE 1.23%)
  })

  model({
    # rxode2 time is postmenstrual age in weeks
    pma <- time

    bwmax <- SEXF * exp(lbwmax_female) + (1 - SEXF) * exp(lbwmax_male)
    tmat50 <- SEXF * exp(ltmat50_female) + (1 - SEXF) * exp(ltmat50_male)
    hill <- SEXF * exp(lhill_female) + (1 - SEXF) * exp(lhill_male)
    kgrow <- exp(lkgrow)
    logismax <- SEXF * exp(llogismax_female) + (1 - SEXF) * exp(llogismax_male)

    # Eq. (2): sigmoid component times logistic component
    fsig <- pma^hill / (tmat50^hill + pma^hill)
    flogis <- 1 / (exp(-kgrow * pma) + 1 / logismax)
    bw <- bwmax * fsig * flogis

    bw ~ prop(propSd)
  })
}
