Okubo_2021_apremilast_pasi50 <- function() {
  description <- paste(
    "Longitudinal logistic exposure-response model for the probability of a",
    "PASI-50 response (>= 50% reduction from baseline in Psoriasis Area and",
    "Severity Index) in Japanese and non-Japanese adults with moderate to",
    "severe plaque psoriasis treated with oral apremilast or placebo",
    "(Okubo 2021). The logit is a baseline intercept plus a placebo effect",
    "with first-order onset in time plus an Emax function of the individual",
    "steady-state 12 h AUC of apremilast. No race effect and no",
    "between-subject variability were retained. The AUC is a covariate column",
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
    concentration = "prob_pasi50 (probability of PASI-50 response, 0-1; not a drug concentration)"
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
        "subject: the drug term enters at full size from the first visit.",
        "Median week-16 values 1065.4 (10 mg BID), 2030.6 (20 mg BID) and",
        "3169.2 (30 mg BID) ng*h/mL (Results 3.3)."
      ),
      source_name = "AUCtau,ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1433L,
    n_studies = 3L,
    age_range = "Adults; study means 46.5-51.4 years in the PPK subset (Okubo 2021 Table 3)",
    weight_range = "Study means 70.2-96.3 kg in the PPK subset (Okubo 2021 Table 3)",
    disease_state = "Moderate to severe chronic plaque psoriasis (PASI >= 12, BSA >= 10%)",
    dose_range = "Placebo or oral apremilast 10, 20 or 30 mg twice daily",
    regions = "Japan (PSOR-011); North America, Europe and Australia (PSOR-005, PSOR-008/ESTEEM 1)",
    notes = paste(
      "9087 PASI-50 observations from 1433 subjects (Okubo 2021 Results 3.3).",
      "Analysis windows: weeks 2-24 in PSOR-005, 2-16 in PSOR-008/ESTEEM 1 and",
      "2-40 in PSOR-011 (Methods 2.5). Fitted in NONMEM 7.3 (Appendix S4)."
    )
  )

  ini({
    # Okubo 2021 Table S2, PASI-50 column ('typical value, in the absence of
    # covariate effects'). Logit-scale parameters are printed untransformed.
    bsl_pbo <- -5.84; label("Baseline logit of PASI-50 response at time 0 (logit units)") # Table S2: Intercept -5.84 (RSE 9.4%)
    asym_pbo <- 4.50; label("Asymptotic placebo effect on the PASI-50 logit (logit units)") # Table S2: Placebo effect, logit 4.50 (RSE 10.2%)
    lkpbo <- log(0.282); label("First-order onset rate constant of the placebo effect (1/week)") # Table S2: K placebo 0.282 1/week (RSE 15.6%)
    emax <- 2.91; label("Maximum apremilast effect on the PASI-50 logit (logit units)") # Table S2: Emax 2.91 (RSE 16.3%)
    lauc50 <- log(1656); label("Steady-state AUCtau giving half the maximum drug effect (ng*h/mL)") # Table S2: E50 1656 ng*h/mL (RSE 47.5%)

    # Table S2 reports 'no effect observed' for Japanese race and a
    # between-subject variability fixed to 0 for every parameter.
    # Placeholder residual: the source likelihood is Bernoulli.
    addSd_prob_pasi50 <- fixed(0.001); label("Placeholder additive residual SD on prob_pasi50 (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    kpbo <- exp(lkpbo)
    auc50 <- exp(lauc50)

    # Appendix S5: the PASI-50 model has the same form as the PASI-75 model --
    # placebo effect with an exponential onset in time (weeks) plus an Emax
    # function of the individual AUCtau,ss.
    pbo <- asym_pbo * (1 - exp(-kpbo * time))
    drug <- emax * AUC_APREMILAST / (auc50 + AUC_APREMILAST)
    logit_pasi50 <- bsl_pbo + pbo + drug

    prob_pasi50 <- expit(logit_pasi50)
    prob_pasi50 ~ add(addSd_prob_pasi50)
  })
}
