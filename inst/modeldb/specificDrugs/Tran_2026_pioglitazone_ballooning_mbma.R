Tran_2026_pioglitazone_ballooning_mbma <- function() {
  description <- "MBMA. Proportional-odds (cumulative-logit) ordered-categorical model of the NASH CRN hepatocyte ballooning score (0-2) under pioglitazone in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to patient-level virtual data that the authors resampled from the arm-level histology of five published trials. The cumulative logit of P(score <= m) is the running sum of the thresholds plus a drug term alpha * dose * t / (ET50 + t), linear in the daily dose (which fitted better than log(dose), dOFV = 56) with a slow onset (ET50 = 17.4 months). The drug term is zero for placebo and pre-treatment records (DOSE_PIO_MGD = 0). Random effect: an additive eta on the logit, shared by both thresholds. Because the random effect was estimated on resampled virtual patients, simulations reproduce arm-level score distributions, not real between-patient spread. Companion endpoint models from the same paper: modellib('Tran_2026_pioglitazone_steatosis_mbma'), modellib('Tran_2026_pioglitazone_inflammation_mbma'), modellib('Tran_2026_pioglitazone_fibrosis_mbma'), modellib('Tran_2026_pioglitazone_alt_mbma') and modellib('Tran_2026_pioglitazone_ast_mbma')."

  reference <- paste(
    "Tran QT, Bui TT, Ngo LT, Yang BR, Baek IH, Nguyen VH, Lee KA, Yun HY,",
    "Chae JW, Lee S, Kim JH, Jung W. Model-Based Meta-Analysis of the",
    "Relationship Between Pioglitazone and Histological Outcomes in Metabolic",
    "Dysfunction-Associated Steatohepatitis Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70034.",
    "doi:10.1002/psp4.70034.",
    sep = " "
  )
  vignette <- "Tran_2026_pioglitazone_mash"

  units <- list(
    time = "month",
    dosing = "(none; the daily pioglitazone dose enters as the DOSE_PIO_MGD covariate in mg/day, not as dosing events)",
    concentration = "(ordinal hepatocyte ballooning score 0-2; the observation ballooning is a histology grade and the outputs p0-p2 are the probabilities of each grade, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_PIO_MGD = list(
      description = "Daily pioglitazone dose of the treatment arm.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: 0 on placebo records AND on the pre-treatment (time 0) biopsy record of a pioglitazone arm, and the arm's daily dose (30 or 45 mg/day in the source trials) on post-treatment records. Enters linearly (alpha * dose), so a dose of 0 gives no drug effect.",
      source_name = "dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 540L,
    n_studies = 8L,
    age_range = "arm mean ages 40.4-55 years (Tran 2026 Table S2)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "biopsy-proven metabolic dysfunction-associated steatohepatitis (MASH, formerly NASH), with or without type 2 diabetes or prediabetes",
    dose_range = "pioglitazone 30 or 45 mg/day orally (Belfort 2006 and Cusi 2016 started at 30 mg/day and increased to 45 mg/day after 2 months), or placebo; treatment duration 6-24 months",
    regions = "not reported at study level",
    notes = "Ballooning model: five of the eight trials contributed ballooning data on the 0-2 scale (Belfort 2006, Cusi 2016, Huang 2021, Sanyal 2010, Sharma 2012; Tran 2026 Table S3). Pre- and post-treatment arm-level data (score counts or mean and SD) were converted to 500 virtual patients per arm and time by bootstrap resampling or a beta distribution (Tran 2026 Methods 2.2.1 and the supplementary code). The 540 subjects are the total across all eight trials in the meta-analysis."
  )

  ini({
    # Tran 2026 Methods 2.2.1, Results 3.2.1 and Table 1 (Ballooning):
    #   logit[P(Y <= m)] = sum_{L=0}^{m} beta_L
    #                      + alpha * Dose * Time / (ET50 + Time) + eta
    # beta_0 is the cumulative logit of the lowest grade; beta_1 steps it to
    # the cumulative logit of P(score <= 1).
    logite0 <- -1.12
    label("Cumulative logit of P(ballooning score <= 0) without drug (beta_0; unitless)") # Tran 2026 Table 1, Ballooning beta0 = -1.12 (RSE 2.7%)
    logitstep1 <- 1.7
    label("Increment from the logit of P(score <= 0) to the logit of P(score <= 1) (beta_1; unitless)") # Tran 2026 Table 1, Ballooning beta1 = 1.7 (RSE 1.8%)

    e_dose_pio_mgd_ballooning <- 0.0835
    label("Slope of the cumulative logit on the daily dose (alpha; per mg/day)") # Tran 2026 Table 1, Ballooning alpha = 0.0835 (RSE 15%); linear dose form per Results 3.2.1
    let50 <- log(17.4)
    label("Log time to half of the maximal drug effect (ET50; month)") # Tran 2026 Table 1, Ballooning 'EC50 (month)' = 17.4 (RSE 27.1%); the cell's equation and the Discussion both call it ET50

    # Table 1 OMEGA row. IIV_beta is additive on the logit (the Methods 2.2.1
    # equation adds eta to the linear predictor). No IIV on alpha was
    # reported for ballooning.
    etalogite0 ~ 0.232 # Tran 2026 Table 1, Ballooning 'IIV_beta' = 0.232 (RSE 16.7%)
  })

  model({
    et50 <- exp(let50)
    drug <- e_dose_pio_mgd_ballooning * DOSE_PIO_MGD * t / (et50 + t)

    # Cumulative probabilities P(score <= m), then per-grade probabilities
    # (Tran 2026 Methods 2.2.1).
    lp <- logite0 + etalogite0 + drug
    pcum0 <- expit(lp)
    pcum1 <- expit(lp + logitstep1)
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- 1 - pcum1

    # rxode2 ordinal endpoint: categories 1-3 carry probabilities p0, p1 and
    # 1 - p0 - p1, i.e. simulated category k is score k - 1.
    ballooning ~ c(p0, p1)
  })
}
