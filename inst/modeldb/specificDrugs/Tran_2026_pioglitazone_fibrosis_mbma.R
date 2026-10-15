Tran_2026_pioglitazone_fibrosis_mbma <- function() {
  description <- "MBMA. Proportional-odds (cumulative-logit) ordered-categorical model of the NASH CRN fibrosis stage (0-4) under pioglitazone in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to patient-level virtual data that the authors resampled from the arm-level histology of all eight published trials. The cumulative logit of P(stage <= m) is the running sum of the thresholds plus a time-independent drug term alpha * log(dose), so any post-treatment record on drug carries the full effect. The drug term is zero for placebo and pre-treatment records (DOSE_PIO_MGD = 0). Random effect: an additive eta on the logit, shared by every threshold. Because the random effect was estimated on resampled virtual patients, simulations reproduce arm-level stage distributions, not real between-patient spread. Companion endpoint models from the same paper: modellib('Tran_2026_pioglitazone_steatosis_mbma'), modellib('Tran_2026_pioglitazone_inflammation_mbma'), modellib('Tran_2026_pioglitazone_ballooning_mbma'), modellib('Tran_2026_pioglitazone_alt_mbma') and modellib('Tran_2026_pioglitazone_ast_mbma')."

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
    concentration = "(ordinal fibrosis stage 0-4; the observation fibrosis is a histology stage and the outputs p0-p4 are the probabilities of each stage, not a drug concentration)"
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
    notes = "Fibrosis model: all eight trials contributed fibrosis data on the 0-4 scale (Aithal 2008, Belfort 2006, Cusi 2016, Huang 2021, Lutchman 2007, Promrat 2006, Sanyal 2010, Sharma 2012; Tran 2026 Table S3). Cusi 2016 reported stages 3 and 4 pooled. Pre- and post-treatment arm-level data (stage counts or mean and SD) were converted to 500 virtual patients per arm and time by bootstrap resampling or a beta distribution (Tran 2026 Methods 2.2.1 and the supplementary code). The 540 subjects are the total across all eight trials in the meta-analysis."
  )

  ini({
    # Tran 2026 Methods 2.2.1 and Table 1 (Fibrosis):
    #   logit[P(Y <= m)] = sum_{L=0}^{m} beta_L + alpha * log(dose) + eta
    # beta_0 is the cumulative logit of the lowest stage; beta_1 to beta_3
    # are the increments that step it to the next cumulative threshold.
    logite0 <- -0.946
    label("Cumulative logit of P(fibrosis stage <= 0) without drug (beta_0; unitless)") # Tran 2026 Table 1, Fibrosis beta0 = -0.946 (RSE 2.5%)
    logitstep1 <- 1.36
    label("Increment from the logit of P(stage <= 0) to the logit of P(stage <= 1) (beta_1; unitless)") # Tran 2026 Table 1, Fibrosis beta1 = 1.36 (RSE 1.6%)
    logitstep2 <- 0.946
    label("Increment from the logit of P(stage <= 1) to the logit of P(stage <= 2) (beta_2; unitless)") # Tran 2026 Table 1, Fibrosis beta2 = 0.946 (RSE 2%)
    logitstep3 <- 1.28
    label("Increment from the logit of P(stage <= 2) to the logit of P(stage <= 3) (beta_3; unitless)") # Tran 2026 Table 1, Fibrosis beta3 = 1.28 (RSE 2.5%)

    e_dose_pio_mgd_fibrosis <- 0.168
    label("Slope of the cumulative logit on the natural log of the daily dose (alpha; per log(mg/day))") # Tran 2026 Table 1, Fibrosis alpha = 0.168 (RSE 5.8%)

    # Table 1 OMEGA row. IIV_beta is additive on the logit (the Methods 2.2.1
    # equation adds eta to the linear predictor). No IIV on alpha was
    # reported for fibrosis.
    etalogite0 ~ 0.086 # Tran 2026 Table 1, Fibrosis 'IIV_beta' = 0.086 (RSE 37.6%)
  })

  model({
    # Natural log of the daily dose; the drug term is zero without drug.
    ldose <- 0
    if (DOSE_PIO_MGD > 0) {
      ldose <- log(DOSE_PIO_MGD)
    }
    drug <- e_dose_pio_mgd_fibrosis * ldose

    # Cumulative probabilities P(stage <= m), then per-stage probabilities
    # (Tran 2026 Methods 2.2.1).
    lp <- logite0 + etalogite0 + drug
    pcum0 <- expit(lp)
    pcum1 <- expit(lp + logitstep1)
    pcum2 <- expit(lp + logitstep1 + logitstep2)
    pcum3 <- expit(lp + logitstep1 + logitstep2 + logitstep3)
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- pcum2 - pcum1
    p3 <- pcum3 - pcum2
    p4 <- 1 - pcum3

    # rxode2 ordinal endpoint: categories 1-5 carry probabilities p0 to p3
    # and 1 - p0 - p1 - p2 - p3, i.e. simulated category k is stage k - 1.
    fibrosis ~ c(p0, p1, p2, p3)
  })
}
