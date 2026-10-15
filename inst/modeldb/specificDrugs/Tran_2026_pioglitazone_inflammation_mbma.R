Tran_2026_pioglitazone_inflammation_mbma <- function() {
  description <- "MBMA. Proportional-odds (cumulative-logit) ordered-categorical model of the NASH CRN lobular inflammation score (0-3) under pioglitazone in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to patient-level virtual data that the authors resampled from the arm-level histology of four published trials. The cumulative logit of P(score <= m) is the running sum of the thresholds plus a time-independent drug term alpha * log(dose), so any post-treatment record on drug carries the full effect. The drug term is zero for placebo and pre-treatment records (DOSE_PIO_MGD = 0). Random effects: an additive eta on the logit (shared by every threshold) and an exponential eta on alpha. Because the random effects were estimated on resampled virtual patients, simulations reproduce arm-level score distributions, not real between-patient spread. Companion endpoint models from the same paper: modellib('Tran_2026_pioglitazone_steatosis_mbma'), modellib('Tran_2026_pioglitazone_ballooning_mbma'), modellib('Tran_2026_pioglitazone_fibrosis_mbma'), modellib('Tran_2026_pioglitazone_alt_mbma') and modellib('Tran_2026_pioglitazone_ast_mbma')."

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
    concentration = "(ordinal lobular inflammation score 0-3; the observation inflammation is a histology grade and the outputs p0-p3 are the probabilities of each grade, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_PIO_MGD = list(
      description = "Daily pioglitazone dose of the treatment arm.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: 0 on placebo records AND on the pre-treatment (time 0) biopsy record of a pioglitazone arm, and the arm's daily dose (30 or 45 mg/day in the source trials) on post-treatment records. The drug term uses the natural log of the dose and is defined as zero when the dose is 0 (log(0) is undefined; the placebo arms of the source data show no drug effect). Because this endpoint has no time term, the pre-treatment record must carry 0 or it would already show the full drug effect.",
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
    notes = "Inflammation model: four of the eight trials contributed inflammation data on the 0-3 scale (Cusi 2016, Huang 2021, Sanyal 2010, Sharma 2012; Tran 2026 Table S3). Only Cusi 2016 contributed a 45 mg/day arm, so the 45 mg/day prediction rests on that single trial (Tran 2026 Discussion). Pre- and post-treatment arm-level data (score counts or mean and SD) were converted to 500 virtual patients per arm and time by bootstrap resampling or a beta distribution (Tran 2026 Methods 2.2.1 and the supplementary code). The 540 subjects are the total across all eight trials in the meta-analysis."
  )

  ini({
    # Tran 2026 Methods 2.2.1 and Table 1 (Inflammation):
    #   logit[P(Y <= m)] = sum_{L=0}^{m} beta_L + alpha * log(dose) + eta
    # beta_0 is the cumulative logit of the lowest grade; beta_1 and beta_2
    # are the increments that step it to the next cumulative threshold.
    logite0 <- -2.44
    label("Cumulative logit of P(inflammation score <= 0) without drug (beta_0; unitless)") # Tran 2026 Table 1, Inflammation beta0 = -2.44 (RSE 3%)
    logitstep1 <- 2.58
    label("Increment from the logit of P(score <= 0) to the logit of P(score <= 1) (beta_1; unitless)") # Tran 2026 Table 1, Inflammation beta1 = 2.58 (RSE 2.1%)
    logitstep2 <- 3.06
    label("Increment from the logit of P(score <= 1) to the logit of P(score <= 2) (beta_2; unitless)") # Tran 2026 Table 1, Inflammation beta2 = 3.06 (RSE 2.6%)

    e_dose_pio_mgd_inflammation <- 0.413
    label("Slope of the cumulative logit on the natural log of the daily dose (alpha; per log(mg/day))") # Tran 2026 Table 1, Inflammation alpha = 0.413 (RSE 10%)

    # Table 1 OMEGA rows. IIV_beta is additive on the logit (the Methods
    # 2.2.1 equation adds eta to the linear predictor); IIV_alpha is taken as
    # exponential (see vignette).
    etalogite0 ~ 0.35 # Tran 2026 Table 1, Inflammation 'IIV_beta' = 0.35 (RSE 18.2%)
    etae_dose_pio_mgd_inflammation ~ 0.329 # Tran 2026 Table 1, Inflammation 'IIV_alpha' = 0.329 (RSE 17.1%)
  })

  model({
    # Natural log of the daily dose; the drug term is zero without drug.
    ldose <- 0
    if (DOSE_PIO_MGD > 0) {
      ldose <- log(DOSE_PIO_MGD)
    }
    drug <- e_dose_pio_mgd_inflammation * exp(etae_dose_pio_mgd_inflammation) * ldose

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
    inflammation ~ c(p0, p1, p2)
  })
}
