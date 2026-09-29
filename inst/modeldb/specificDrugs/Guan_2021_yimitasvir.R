Guan_2021_yimitasvir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral yimitasvir (an HCV NS5A",
    "inhibitor) in Chinese healthy volunteers and patients with chronic HCV",
    "genotype 1 infection (Guan 2021; N = 219 across six studies), with",
    "sequential zero-order (duration Td) then first-order absorption,",
    "first-order elimination, a linear decrease in relative bioavailability",
    "of 12.9% per 100 mg above 100 mg, food effects on ka and F, sex and",
    "baseline ALT effects on CL/F, and a patient effect on Td."
  )
  reference <- paste(
    "Guan X, Tang X, Zhang Y, Xie H, Luo L, Wu D, Chen R, Hu P. (2021).",
    "Population Pharmacokinetic Analysis of Yimitasvir in Chinese Healthy",
    "Volunteers and Patients With Chronic Hepatitis C Virus Infection.",
    "Front Pharmacol 11:617122. doi:10.3389/fphar.2020.617122"
  )
  vignette <- "Guan_2021_yimitasvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "yimitasvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "yimitasvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "yimitasvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Guan 2021 Table 3 'CL/F ~ Female' = -0.251 in the exponential",
        "categorical form of Eq. 5, CL/F x exp(-0.251 * SEXF). 93 of 219",
        "subjects (42.5%) were female (Table 2)."
      ),
      source_name = "Gender"
    ),
    ALT = list(
      description = "Baseline serum alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline value (the paper did not test time-varying ALT;",
        "Discussion). Power effect normalised by the population median",
        "31.6 IU/L (Eq. 4; Table 2 total median; Results 'typical values",
        "(for a healthy male volunteer with ALT value of 31.6 IU/L...)').",
        "IU/L and U/L are numerically identical."
      ),
      source_name = "ALT"
    ),
    HCV_POS = list(
      description = "Chronic hepatitis C patient (vs healthy volunteer) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = paste(
        "Guan 2021 'disease status' covariate: 1 = patient with chronic HCV",
        "genotype 1 infection (phase 1b and phase 2 studies, n = 147),",
        "0 = healthy volunteer (phase 1 studies, n = 72). Acts on the",
        "zero-order absorption duration only: Td x exp(-0.416 * HCV_POS),",
        "2.17 h in healthy volunteers vs 1.43 h in patients (Table 3 'Td ~",
        "Patient'; Results)."
      ),
      source_name = "disease status (healthy volunteers vs patients)"
    ),
    DOSE = list(
      description = "Administered yimitasvir dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-dose amount (30-600 mg in the analysis dataset; once-daily",
        "dosing, so per-dose and daily dose coincide). Enters the linear",
        "bioavailability model of Eq. 3: F = 1 for DOSE <= 100 mg and",
        "F = 1 - Alpha * (DOSE - 100) / 100 for DOSE > 100 mg, with Alpha",
        "fixed to 0.129 (Table 3). Supply the same value as the dosing",
        "record's amt, and place the DOSE column after amt in the event",
        "table (rxode2 does not pass a DOSE column that precedes amt to",
        "the model)."
      ),
      source_name = "Dose"
    ),
    FED_HIGHFAT = list(
      description = "High-fat meal at dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not high-fat fed; reference is fasted at least 10 h)",
      notes = paste(
        "Guan 2021 food status level 2 (high-fat meal; food-effect",
        "crossover study CTR20150123, Table 1 / Table 2). Effects",
        "exp(-2.40) on ka and exp(-0.485) on F (Table 3 'Ka ~ Food2',",
        "'F ~ Food2'; ka -90.9%, F -38.5% per the abstract). Food status",
        "0 (the reference) = fasted at least 10 h before a morning dose;",
        "the paper also assigned the phase 2 sparse-sampling patients",
        "(dosed at least 2 h before or after a meal) to level 0 for",
        "convenience (Discussion)."
      ),
      source_name = "Food (= 2)"
    ),
    MEAL_PREDOSE_4H = list(
      description = "Dose taken at least 4 h after the start of a meal indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted at least 10 h before a morning dose, or high-fat fed)",
      notes = paste(
        "Guan 2021 food status level 1: phase 1b patients (CTR20150549)",
        "took yimitasvir in the evening, at least 4 h after dinner",
        "(Table 1 'Fasted at least 4 h'; Discussion 'administered",
        "yimitasvir in the evening after dinner more than 4 h'). Effects",
        "exp(-1.71) on ka and exp(-0.341) on F (Table 3 'Ka ~ Food1',",
        "'F ~ Food1'). Mutually exclusive with FED_HIGHFAT."
      ),
      source_name = "Food (= 1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 219L,
    n_studies = 6L,
    n_observations = 3540L,
    age_range = "18-75 years",
    age_median = "38 years",
    weight_range = "44-100 kg",
    weight_median = "62 kg",
    sex_female_pct = 42.5,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "72 healthy volunteers and 147 patients with chronic HCV genotype 1",
      "infection (DAA-naive)"
    ),
    dose_range = paste(
      "30-600 mg single dose; 30-400 mg once daily for 7 days; 100 or",
      "200 mg once daily for 12 weeks with sofosbuvir 400 mg"
    ),
    regions = "China",
    co_medication = "sofosbuvir 400 mg in 129 phase 2 patients (not a significant covariate)",
    notes = paste(
      "Four phase 1 studies (SAD, MAD, high-dose supplement, food-effect",
      "crossover), one phase 1b and one phase 2 study (Guan 2021 Table 1);",
      "baseline demographics in Table 2. Median ALT 31.6 IU/L (range",
      "5.0-666). Phoenix NLME 8.1, FOCE-ELS."
    )
  )

  ini({
    lka <- log(1.31); label("First-order absorption rate constant ka (1/h)") # Table 3 'Ka' = 1.31 h^-1
    lcl <- log(13.8); label("Apparent clearance CL/F (L/h)") # Table 3 'CL/F' = 13.8 L/h
    lvc <- log(188); label("Apparent central volume V1/F (L)") # Table 3 'V1/F' = 188 L
    lq <- log(3.96); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F' = 3.96 L/h
    lvp <- log(58.6); label("Apparent peripheral volume V2/F (L)") # Table 3 'V2/F' = 58.6 L
    ld1 <- log(2.17); label("Zero-order absorption duration Td in healthy volunteers (h)") # Table 3 'Td' = 2.17 h
    lfdepot <- fixed(log(1)); label("Relative bioavailability at doses <= 100 mg (unitless)") # Table 3 'F' = 1 (fixed); Eq. 3 theta_F

    e_dose_fdepot <- fixed(0.129); label("Fractional decrease in F per 100 mg above 100 mg (Alpha, unitless)") # Table 3 'Alpha' = 0.129 (fixed); Eq. 3

    e_meal_predose_4h_ka <- -1.71; label("Food status 1 (dosed >= 4 h after dinner) effect on ka, exponential (unitless)") # Table 3 'Ka ~ Food1' = -1.71
    e_fed_highfat_ka <- -2.40; label("Food status 2 (high-fat meal) effect on ka, exponential (unitless)") # Table 3 'Ka ~ Food2' = -2.40
    e_meal_predose_4h_fdepot <- -0.341; label("Food status 1 (dosed >= 4 h after dinner) effect on F, exponential (unitless)") # Table 3 'F ~ Food1' = -0.341
    e_fed_highfat_fdepot <- -0.485; label("Food status 2 (high-fat meal) effect on F, exponential (unitless)") # Table 3 'F ~ Food2' = -0.485
    e_sexf_cl <- -0.251; label("Female effect on CL/F, exponential (unitless)") # Table 3 'CL/F ~ Female' = -0.251
    e_alt_cl <- -0.0950; label("Power exponent of ALT/31.6 on CL/F (unitless)") # Table 3 'CL/F ~ ALT' = -0.0950
    e_hcv_pos_d1 <- -0.416; label("HCV patient effect on Td, exponential (unitless)") # Table 3 'Td ~ Patient' = -0.416

    # Table 3 Random effect rows 'CL/F' = 0.485, 'V1/F' = 0.736 read as SDs
    # (omega) and 'CL/F ~ V1/F' = 0.565 read as their correlation; the
    # covariance is 0.565 * 0.485 * 0.736 = 0.20168. See the vignette
    # 'IIV scale' section: the Figure 5 AUCss 5th percentile (3.26) fixes
    # this reading.
    etalcl + etalvc ~ c(0.235225, 0.201682, 0.541696) # Table 3; variances 0.485^2 and 0.736^2
    propSd <- 0.305; label("Proportional residual error (fraction)") # Table 3 'sigma' = 0.305
  })

  model({
    # Eq. 5 exponential categorical effects; Eq. 4 median-normalised power
    ka <- exp(lka + e_meal_predose_4h_ka * MEAL_PREDOSE_4H + e_fed_highfat_ka * FED_HIGHFAT)
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF) * (ALT / 31.6)^e_alt_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    d1 <- exp(ld1 + e_hcv_pos_d1 * HCV_POS)

    # Eq. 3: F = theta_F for DOSE <= 100 mg, theta_F - Alpha * (DOSE - 100) / 100 above
    fdose <- exp(lfdepot) - e_dose_fdepot * (DOSE - 100) / 100
    if (DOSE <= 100) fdose <- exp(lfdepot)
    fdepot <- fdose * exp(e_meal_predose_4h_fdepot * MEAL_PREDOSE_4H + e_fed_highfat_fdepot * FED_HIGHFAT)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Sequential zero-order input into the depot over d1, then first-order ka
    dur(depot) <- d1
    f(depot) <- fdepot

    # mg / L = ug/mL; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
