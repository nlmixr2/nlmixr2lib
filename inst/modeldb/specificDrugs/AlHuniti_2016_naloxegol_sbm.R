AlHuniti_2016_naloxegol_sbm <- function() {
  description <- "Longitudinal ordered-categorical (cumulative-logit) Markov dose-response model of the daily number of spontaneous bowel movements (SBMs, 0 to >= 6 per day) on oral naloxegol in adults with noncancer pain and opioid-induced constipation (Al-Huniti 2016). Separate threshold logits for the 2-week baseline period and the 12-week treatment period, each conditional on whether the patient had at least one SBM on the previous day; an exponential-plateau placebo (non-drug) effect from the first day of treatment; and an Emax function of the daily naloxegol dose, reduced in patients whose baseline maintenance opioids were weak opioids only. Four correlated additive random effects on the post-baseline and baseline intercepts, the placebo magnitude and the Emax."
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
    concentration = "(no concentration output; the observation sbm is the daily spontaneous-bowel-movement count, ordinal 0 to >= 6, and p0-p6 are the probabilities of each count)"
  )

  covariateData <- list(
    DOSE = list(
      description = "Naloxegol daily dose administered on the day of the record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (b) of the DOSE canonical: the current daily naloxegol dose (0 for placebo and for every day of the 2-week baseline period; 12.5 or 25 mg/day on treatment). Enters the logit through the Emax function DOSE * emax / (ed50 + DOSE). The paper selected dose over individual predicted Cavg (Supplementary Table S3, runs 8 vs 11b).",
      source_name = "DOSE"
    ),
    PDV = list(
      description = "Number of spontaneous bowel movements recorded on the previous day (Markov state)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = "Source column PDVORD (previous day's SBM count, truncated at 6). Only PDV > 0 versus PDV = 0 matters: it selects which of the two sets of threshold logits applies. In simulation the column must be updated from the previously simulated count; see the vignette.",
      source_name = "PDVORD"
    ),
    OPIOID_WEAK_ONLY = list(
      description = "Baseline maintenance opioid regimen consisted of weak opioids only (1) versus a strong opioid alone or weak plus strong opioids (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (strong opioid, or weak plus strong opioids, at baseline)",
      notes = "Source BOPWK = 1 when BOTP = 1 (weak opioid only); BOTP = 2 (strong) and BOTP = 3 (weak + strong) both give 0. Lowers Emax by 0.978 logit units (67% of the typical 1.46).",
      source_name = "BOPWK"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "In the working full model on the placebo magnitude and on Emax; removed by backward elimination (Results).",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "In the working full model on Emax; removed by backward elimination.",
      source_name = "BMI"
    ),
    RACE_BLACK = list(
      description = "Black or African American race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not Black)",
      notes = "Black or Asian race on Emax was in the working full model; removed by backward elimination.",
      source_name = "RACEID"
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened in the full model (Supplementary Table S1); not retained.",
      source_name = "SEXF"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1331L,
    n_studies = 2L,
    age_range = "Adults; mean age 51.9-52.9 years by study arm (Table 1)",
    sex_female_pct = 62,
    race_ethnicity = "About 79% White, 19% Black, 1% Asian, 2% other (Table 1)",
    disease_state = "Chronic noncancer pain treated with stable maintenance opioids (30-1000 mg/day oral morphine equivalents for >= 4 weeks) and opioid-induced constipation (< 3 SBMs per week).",
    dose_range = "Naloxegol 12.5 or 25 mg orally once daily, or placebo, for 12 weeks.",
    regions = "Multinational (phase III KODIAC 4 and KODIAC 5).",
    notes = "Intention-to-treat population of the two identical phase III studies KODIAC 4 (NCT01309841) and KODIAC 5 (NCT01323790). 108,051 daily diary entries from day -13 to day 84; 55.1% recorded 0 and 40.1% recorded 1 SBM. About two thirds of patients took strong opioids at baseline."
  )

  ini({
    # Threshold logits, Supplementary Table S4 and the final control stream.
    # The model is written for P(SBM >= m): b1 is the logit of P(SBM >= 1) and
    # each later threshold subtracts exp(lb<m>). 'post' = treatment period
    # (time > 0), 'base' = baseline period (time <= 0); 'pdv0' / 'pdv1' =
    # previous-day SBM count zero / non-zero.
    b1_post_pdv0 <- -0.344
    label("Logit of P(SBM >= 1), treatment period, no SBM on the previous day (logit)") # Table S4: beta1 post baseline [prior SBM = 0] = -0.344
    b1_post_pdv1 <- -0.676
    label("Logit of P(SBM >= 1), treatment period, >= 1 SBM on the previous day (logit)") # Table S4: beta1 post baseline [prior SBM != 0] = -0.676
    b1_base_pdv0 <- -1.42
    label("Logit of P(SBM >= 1), baseline period, no SBM on the previous day (logit)") # Table S4: beta1 baseline [prior SBM = 0] = -1.42
    b1_base_pdv1 <- -2.15
    label("Logit of P(SBM >= 1), baseline period, >= 1 SBM on the previous day (logit)") # Table S4: beta1 baseline [prior SBM != 0] = -2.15

    lb2_post_pdv0 <- 1.21
    label("Log of the logit decrement from threshold 1 to 2, treatment period, previous-day SBM = 0 (log logit)") # Table S4: LN beta2 post baseline [prior SBM = 0] = 1.21
    lb2_post_pdv1 <- 1.29
    label("Log of the logit decrement from threshold 1 to 2, treatment period, previous-day SBM > 0 (log logit)") # Table S4: LN beta2 post baseline [prior SBM != 0] = 1.29
    lb2_base_pdv0 <- 1.3
    label("Log of the logit decrement from threshold 1 to 2, baseline period, previous-day SBM = 0 (log logit)") # Table S4: LN beta2 baseline [prior SBM = 0] = 1.3
    lb2_base_pdv1 <- 1.17
    label("Log of the logit decrement from threshold 1 to 2, baseline period, previous-day SBM > 0 (log logit)") # Table S4: LN beta2 baseline [prior SBM != 0] = 1.17
    lb3_post_pdv0 <- 0.634
    label("Log of the logit decrement from threshold 2 to 3, treatment period, previous-day SBM = 0 (log logit)") # Table S4: LN beta3 post baseline [prior SBM = 0] = 0.634
    lb3_post_pdv1 <- 0.704
    label("Log of the logit decrement from threshold 2 to 3, treatment period, previous-day SBM > 0 (log logit)") # Table S4: LN beta3 post baseline [prior SBM != 0] = 0.704
    lb3_base_pdv0 <- 0.995
    label("Log of the logit decrement from threshold 2 to 3, baseline period, previous-day SBM = 0 (log logit)") # Table S4: LN beta3 baseline [prior SBM = 0] = 0.995
    lb3_base_pdv1 <- 0.845
    label("Log of the logit decrement from threshold 2 to 3, baseline period, previous-day SBM > 0 (log logit)") # Table S4: LN beta3 baseline [prior SBM != 0] = 0.845
    lb4_post_pdv0 <- 0.381
    label("Log of the logit decrement from threshold 3 to 4, treatment period, previous-day SBM = 0 (log logit)") # Table S4: LN beta4 post baseline [prior SBM = 0] = 0.381
    lb4_post_pdv1 <- 0.293
    label("Log of the logit decrement from threshold 3 to 4, treatment period, previous-day SBM > 0 (log logit)") # Table S4: LN beta4 post baseline [prior SBM != 0] = 0.293
    lb5_post_pdv0 <- 0.018
    label("Log of the logit decrement from threshold 4 to 5, treatment period, previous-day SBM = 0 (log logit)") # Table S4: LN beta5 post baseline [prior SBM = 0] = 0.018
    lb5_post_pdv1 <- 0.331
    label("Log of the logit decrement from threshold 4 to 5, treatment period, previous-day SBM > 0 (log logit)") # Table S4: LN beta5 post baseline [prior SBM != 0] = 0.331
    lb6_post_pdv0 <- 0.0122
    label("Log of the logit decrement from threshold 5 to 6, treatment period, previous-day SBM = 0 (log logit)") # Table S4: LN beta6 post baseline [prior SBM = 0] = 0.0122
    lb6_post_pdv1 <- 0.0558
    label("Log of the logit decrement from threshold 5 to 6, treatment period, previous-day SBM > 0 (log logit)") # Table S4: LN beta6 post baseline [prior SBM != 0] = 0.0558

    # Non-drug (placebo) effect, exponential plateau from the first treatment day.
    emax_placebo <- 0.189
    label("Maximum non-drug (placebo) change in the SBM logit, NDM (logit)") # Table S4: NDM = 0.189
    lthalf_placebo <- 2.92
    label("Log of the half-life of the non-drug effect (log day)") # Table S4: LN ND t1/2 = 2.92 (18.5 days; plateau in about 93 days, Results)

    # Naloxegol dose-response (Emax on dose).
    emax <- 1.46
    label("Maximum naloxegol effect on the SBM logit, Emax (logit)") # Table S4: Emax = 1.46
    e_opioid_weak_only_emax <- -0.978
    label("Additive shift in Emax when baseline opioids were weak opioids only (logit)") # Table S4: Weak baseline opioid type on Emax = -0.978 (67% lower Emax, Results)
    led50 <- 2.94
    label("Log of the naloxegol daily dose giving half the maximum effect, ED50 (log mg)") # Table S4: LN ED50 (mg) = 2.94 (ED50 = 18.9 mg, Results)

    # Random effects: Table S4 variances (diagonal) and covariances (lower
    # triangle), in the order post-baseline intercept, baseline intercept,
    # NDM, Emax. The control stream estimates them through a Cholesky
    # factor (THETA 23-32) on unit-variance ETAs; these are its derived values.
    etab1_post + etab1_base + etaemax_placebo + etaemax ~ c(
      1.91,
      0.347, 0.288,
      -0.566, -0.0225, 1.59,
      -1.56, -0.207, -0.472, 7.14
    )
  })

  model({
    post <- 1 * (time > 0)
    prev <- 1 * (PDV > 0)

    # Threshold logits for the four period x previous-day strata. In the
    # baseline period the control stream sets the decrements for thresholds
    # 4-6 to the constant 20, so P(SBM >= 4) is effectively zero there.
    b1_post <- b1_post_pdv0 + etab1_post + prev * (b1_post_pdv1 - b1_post_pdv0)
    b1_base <- b1_base_pdv0 + etab1_base + prev * (b1_base_pdv1 - b1_base_pdv0)
    b1 <- post * b1_post + (1 - post) * b1_base
    b2 <- post * ((1 - prev) * exp(lb2_post_pdv0) + prev * exp(lb2_post_pdv1)) +
      (1 - post) * ((1 - prev) * exp(lb2_base_pdv0) + prev * exp(lb2_base_pdv1))
    b3 <- post * ((1 - prev) * exp(lb3_post_pdv0) + prev * exp(lb3_post_pdv1)) +
      (1 - post) * ((1 - prev) * exp(lb3_base_pdv0) + prev * exp(lb3_base_pdv1))
    b4 <- post * ((1 - prev) * exp(lb4_post_pdv0) + prev * exp(lb4_post_pdv1)) + (1 - post) * 20
    b5 <- post * ((1 - prev) * exp(lb5_post_pdv0) + prev * exp(lb5_post_pdv1)) + (1 - post) * 20
    b6 <- post * ((1 - prev) * exp(lb6_post_pdv0) + prev * exp(lb6_post_pdv1)) + (1 - post) * 20

    # Non-drug effect: zero in the baseline period, rising from day 1 (the
    # first treatment day, TSFD = 1) toward emax_placebo with half-life thalf.
    ndm <- emax_placebo + etaemax_placebo
    ndk <- log(2) / exp(lthalf_placebo)
    fnd <- post * ndm * (1 - exp(-ndk * (time - 1)))

    # Drug effect: Emax function of the daily naloxegol dose.
    emax_i <- emax + e_opioid_weak_only_emax * OPIOID_WEAK_ONLY + etaemax
    ed50 <- exp(led50)
    fd <- emax_i * DOSE / (ed50 + DOSE)

    lge1 <- b1 + fnd + fd
    lge2 <- lge1 - b2
    lge3 <- lge2 - b3
    lge4 <- lge3 - b4
    lge5 <- lge4 - b5
    lge6 <- lge5 - b6

    pge1 <- expit(lge1)
    pge2 <- expit(lge2)
    pge3 <- expit(lge3)
    pge4 <- expit(lge4)
    pge5 <- expit(lge5)
    pge6 <- expit(lge6)

    # Probability of exactly m SBMs on the day (m = 6 means >= 6).
    p0 <- 1 - pge1
    p1 <- pge1 - pge2
    p2 <- pge2 - pge3
    p3 <- pge3 - pge4
    p4 <- pge4 - pge5
    p5 <- pge5 - pge6
    p6 <- pge6

    # Expected daily SBM count (counts above 6 truncated at 6).
    sbm_mean <- p1 + 2 * p2 + 3 * p3 + 4 * p4 + 5 * p5 + 6 * p6

    # rxode2 ordinal endpoint: simulated category k is an SBM count of k - 1.
    sbm ~ c(p0, p1, p2, p3, p4, p5)
  })
}
