Tran_2026_pioglitazone_steatosis_mbma <- function() {
  description <- "MBMA. Proportional-odds (cumulative-logit) ordered-categorical model of the NASH CRN hepatic steatosis score (0-3) under pioglitazone in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to patient-level virtual data that the authors resampled from the arm-level histology of five published trials. The cumulative logit of P(score <= m) is the running sum of the thresholds plus a drug term alpha * log(dose) * t / (ET50 + t), so the improvement is dose- and time-dependent with an ET50 of 0.584 months. The drug term is zero for placebo and pre-treatment records (DOSE_PIO_MGD = 0). Random effects: an additive eta on the logit (shared by every threshold) and an exponential eta on alpha. Because the random effects were estimated on resampled virtual patients, simulations reproduce arm-level score distributions, not real between-patient spread. Companion endpoint models from the same paper: modellib('Tran_2026_pioglitazone_inflammation_mbma'), modellib('Tran_2026_pioglitazone_ballooning_mbma'), modellib('Tran_2026_pioglitazone_fibrosis_mbma'), modellib('Tran_2026_pioglitazone_alt_mbma') and modellib('Tran_2026_pioglitazone_ast_mbma')."

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
    concentration = "(ordinal steatosis score 0-3; the observation steatosis is a histology grade and the outputs p0-p3 are the probabilities of each grade, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_PIO_MGD = list(
      description = "Daily pioglitazone dose of the treatment arm.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: 0 on placebo records AND on the pre-treatment (time 0) biopsy record of a pioglitazone arm, and the arm's daily dose (30 or 45 mg/day in the source trials) on post-treatment records. The drug term uses the natural log of the dose and is defined as zero when the dose is 0 (log(0) is undefined; the placebo arms of the source data show no drug effect).",
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
    notes = "Steatosis model: five of the eight trials contributed steatosis data (Belfort 2006, Cusi 2016, Huang 2021, Sanyal 2010, Sharma 2012; Tran 2026 Table S3). Pre- and post-treatment arm-level data (score counts or mean and SD) were converted to 500 virtual patients per arm and time by bootstrap resampling or a beta distribution (Tran 2026 Methods 2.2.1 and the supplementary code), and the model was fitted to those virtual patients. Scores reported on a 0-4 scale were excluded so that all steatosis data are on the 0-3 NASH CRN scale. The 540 subjects are the total across all eight trials in the meta-analysis."
  )

  ini({
    # Tran 2026 Methods 2.2.1 and Table 1 (Steatosis):
    #   logit[P(Y <= m)] = sum_{L=0}^{m} beta_L
    #                      + alpha * log(dose) * Time / (ET50 + Time) + eta
    # beta_0 is the cumulative logit of the lowest grade; beta_1 and beta_2
    # are the increments that step it to the next cumulative threshold.
    # Table 1 lists one beta per threshold below the top grade (three for
    # the 0-3 scale), which fixes the 0-indexed reading of the sum.
    logite0 <- -3.07
    label("Cumulative logit of P(steatosis score <= 0) without drug (beta_0; unitless)") # Tran 2026 Table 1, Steatosis beta0 = -3.07 (RSE 1.6%)
    logitstep1 <- 2.12
    label("Increment from the logit of P(score <= 0) to the logit of P(score <= 1) (beta_1; unitless)") # Tran 2026 Table 1, Steatosis beta1 = 2.12 (RSE 1.9%)
    logitstep2 <- 1.31
    label("Increment from the logit of P(score <= 1) to the logit of P(score <= 2) (beta_2; unitless)") # Tran 2026 Table 1, Steatosis beta2 = 1.31 (RSE 1.8%)

    e_dose_pio_mgd_steatosis <- 0.49
    label("Slope of the cumulative logit on the natural log of the daily dose (alpha; per log(mg/day))") # Tran 2026 Table 1, Steatosis alpha = 0.49 (RSE 4.4%)
    let50 <- log(0.584)
    label("Log time to half of the maximal drug effect (ET50; month)") # Tran 2026 Table 1, Steatosis 'EC50 (month)' = 0.584 (RSE 11.8%); the cell's equation and the Discussion both call it ET50

    # Table 1 OMEGA rows. IIV_beta is the eta added to the linear predictor in
    # the Methods 2.2.1 equation, so it is additive on the logit and shared by
    # every threshold; IIV_alpha is taken as exponential (see vignette).
    etalogite0 ~ 0.0298 # Tran 2026 Table 1, Steatosis 'IIV_beta' = 0.0298 (RSE 14.9%)
    etae_dose_pio_mgd_steatosis ~ 0.244 # Tran 2026 Table 1, Steatosis 'IIV_alpha' = 0.244 (RSE 7.7%)
  })

  model({
    et50 <- exp(let50)

    # Natural log of the daily dose; the drug term is zero without drug.
    ldose <- 0
    if (DOSE_PIO_MGD > 0) {
      ldose <- log(DOSE_PIO_MGD)
    }
    drug <- e_dose_pio_mgd_steatosis * exp(etae_dose_pio_mgd_steatosis) * ldose * t / (et50 + t)

    # Cumulative probabilities P(score <= m), then per-grade probabilities
    # (Tran 2026 Methods 2.2.1).
    lp <- logite0 + etalogite0 + drug
    pcum0 <- expit(lp)
    pcum1 <- expit(lp + logitstep1)
    pcum2 <- expit(lp + logitstep1 + logitstep2)
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- pcum2 - pcum1
    p3 <- 1 - pcum2

    # rxode2 ordinal endpoint: categories 1-4 carry probabilities p0, p1,
    # p2 and 1 - p0 - p1 - p2, i.e. simulated category k is score k - 1.
    # (The named form c(p0 = 0, ...) would add the category values to
    # ini() as pseudo-parameters.)
    steatosis ~ c(p0, p1, p2)
  })
}
