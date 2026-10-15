Okubo_2021_apremilast_spga <- function() {
  description <- paste(
    "Longitudinal logistic exposure-response model for the probability of a",
    "static Physician Global Assessment response (sPGA 0 'clear' or 1",
    "'almost clear') in Japanese and non-Japanese adults with moderate to",
    "severe plaque psoriasis treated with oral apremilast or placebo",
    "(Okubo 2021). The logit is a baseline intercept plus a placebo effect",
    "with first-order onset in time plus an Emax function of the individual",
    "steady-state 12 h AUC of apremilast with its own first-order onset;",
    "Japanese race multiplies the intercept, the placebo effect and Emax. No",
    "between-subject variability was retained. The AUC is a covariate column",
    "produced by the companion popPK model Okubo_2021_apremilast."
  )
  reference <- paste(
    "Okubo Y, Ohtsuki M, Komine M, Imafuku S, Kassir N, Petric R, Nemoto O.",
    "Population pharmacokinetic and exposure-response analysis of apremilast",
    "in Japanese subjects with moderate to severe psoriasis.",
    "J Dermatol. 2021;48(11):1652-1664. doi:10.1111/1346-8138.16068"
  )
  vignette <- "Okubo_2021_apremilast"
  units <- list(
    time = "week",
    dosing = "mg",
    concentration = "prob_spga01 (probability of sPGA 0 or 1 response, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    AUC_APREMILAST = list(
      description = "Individual steady-state AUC of apremilast over the 12 h twice-daily dosing interval",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual AUCtau,ss predicted from the final PPK model",
        "(Okubo_2021_apremilast; Methods 2.5). 0 for placebo. Time-fixed per",
        "subject; in this model the drug term is ramped in by its own",
        "first-order onset. Median week-16 values 1065.4 (10 mg BID), 2030.6",
        "(20 mg BID) and 3169.2 (30 mg BID) ng*h/mL (Results 3.3)."
      ),
      source_name = "AUCtau,ss"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator, 1 = Japanese, 0 = non-Japanese",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese; the paper's comparator is predominantly Caucasian)",
      notes = paste(
        "Okubo 2021 Results 3.3: the sPGA model was the only E-R model improved",
        "by an effect of race (Japanese vs Caucasian), acting on both the",
        "placebo and the drug response (Table S2 'If Japanese' rows)."
      ),
      source_name = "race (Japanese vs Caucasian)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1442L,
    n_studies = 3L,
    age_range = "Adults; study means 46.5-51.4 years in the PPK subset (Okubo 2021 Table 3)",
    weight_range = "Study means 70.2-96.3 kg in the PPK subset (Okubo 2021 Table 3)",
    disease_state = "Moderate to severe chronic plaque psoriasis (PASI >= 12, BSA >= 10%)",
    dose_range = "Placebo or oral apremilast 10, 20 or 30 mg twice daily",
    regions = "Japan (PSOR-011); North America, Europe and Australia (PSOR-005, PSOR-008/ESTEEM 1)",
    notes = paste(
      "9094 sPGA observations from 1442 subjects (Okubo 2021 Results 3.3).",
      "Analysis windows: weeks 2-24 in PSOR-005, 2-16 in PSOR-008/ESTEEM 1 and",
      "2-40 in PSOR-011 (Methods 2.5). Fitted in NONMEM 7.3 (Appendix S4)."
    )
  )

  ini({
    # Okubo 2021 Table S2, sPGA (0 or 1) column ('typical value, in the
    # absence of covariate effects'). Logit-scale parameters are printed
    # untransformed.
    bsl_pbo <- -8.81; label("Baseline logit of sPGA 0/1 response at time 0, non-Japanese (logit units)") # Table S2: Intercept -8.81 (RSE 16.9%)
    asym_pbo <- 5.87; label("Asymptotic placebo effect on the sPGA logit, non-Japanese (logit units)") # Table S2: Placebo effect, logit 5.87 (RSE 25.5%)
    lkpbo <- log(0.579); label("First-order onset rate constant of the placebo effect (1/week)") # Table S2: K placebo 0.579 1/week (RSE 14.2%)
    emax <- 4.58; label("Maximum apremilast effect on the sPGA logit, non-Japanese (logit units)") # Table S2: Emax 4.58 (RSE 22.1%)
    lauc50 <- log(1110); label("Steady-state AUCtau giving half the maximum drug effect (ng*h/mL)") # Table S2: E50 1110 ng*h/mL (RSE 47.4%)

    # Onset rate of the drug effect. Appendix S5 states that the sPGA drug
    # effect has an exponential delay component, but Table S2 does not print
    # its rate. Back-solved by the maintainers from the ten week-16 sPGA
    # probabilities printed in Table S1 and Results 3.3: every one requires
    # the drug term to be 0.5784 of its full size at week 16 (rounding
    # interval 0.5783-0.5786), so kdrug = -log(1 - 0.5784) / 16 = 0.0540/week.
    # Figure S3C (population predictions over weeks 2-24) confirms the gradual
    # onset; a constant 0.58 scaling does not reproduce it.
    lkdrug <- log(0.0540); label("First-order onset rate constant of the drug effect (1/week); back-solved, not printed") # back-solved from Table S1 and Results 3.3; see vignette Assumptions and deviations

    # Japanese race multipliers (Table S2 'If Japanese' rows)
    e_race_japanese_bsl_pbo <- 0.682; label("Multiplicative factor on the intercept for Japanese subjects (unitless)") # Table S2: If Japanese x 0.682 (RSE 20.2%)
    e_race_japanese_asym_pbo <- 0.716; label("Multiplicative factor on the placebo effect for Japanese subjects (unitless)") # Table S2: If Japanese x 0.716 (RSE 31.7%)
    e_race_japanese_emax <- 0.661; label("Multiplicative factor on Emax for Japanese subjects (unitless)") # Table S2: If Japanese x 0.661 (RSE 16.2%)

    # Table S2 reports a between-subject variability fixed to 0 for every
    # parameter. Placeholder residual: the source likelihood is Bernoulli.
    addSd_prob_spga01 <- fixed(0.001); label("Placeholder additive residual SD on prob_spga01 (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    kpbo <- exp(lkpbo)
    kdrug <- exp(lkdrug)
    auc50 <- exp(lauc50)

    bsl <- bsl_pbo * e_race_japanese_bsl_pbo^RACE_JAPANESE
    asym <- asym_pbo * e_race_japanese_asym_pbo^RACE_JAPANESE
    emax_i <- emax * e_race_japanese_emax^RACE_JAPANESE

    # Appendix S5: placebo model with an exponential delay component and an
    # Emax drug-effect model with an exponential delay component; time in weeks.
    pbo <- asym * (1 - exp(-kpbo * time))
    drug <- emax_i * AUC_APREMILAST / (auc50 + AUC_APREMILAST) * (1 - exp(-kdrug * time))
    logit_spga01 <- bsl + pbo + drug

    prob_spga01 <- expit(logit_spga01)
    prob_spga01 ~ add(addSd_prob_spga01)
  })
}
