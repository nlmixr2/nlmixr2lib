Boucher_2016_topiramate_mbma <- function() {
  description <- paste0(
    "MBMA. Landmark logistic dose-response model-based meta-analysis of the ",
    "paresthesia incidence rate with topiramate in episodic migraine ",
    "prophylaxis, fitted to 17 study arms from six randomized ",
    "placebo-controlled trials (1650 subjects; Edwards 2003, Silberstein ",
    "2004, Brandes 2004, Diener 2004, Storey 2001 and Silberstein 2006) at ",
    "daily doses of 0, 50, 100 and 200 mg. On the logit scale the arm's ",
    "log-odds of paresthesia is a study-specific placebo intercept plus a ",
    "three-parameter Emax function of the assigned daily dose: ",
    "logit(p) = E0 + eta_study + Emax * DOSE / (ED50 + DOSE). The source ",
    "likelihood is binomial on the arm's event count, so the model carries ",
    "NO residual error of its own and the shipped additive residual is a ",
    "placeholder. Variability is BETWEEN-STUDY only: one normally ",
    "distributed random effect on the placebo intercept (SD 0.17 on the ",
    "logit scale), so the simulation scope is the study-arm paresthesia ",
    "rate, NOT an individual patient's risk. ED50 is 17.5 mg/day, well ",
    "below the lowest active dose studied, which is why the incidence ",
    "rises steeply between placebo and 50 mg and then nearly plateaus: the ",
    "typical rate is 7.2% at placebo, 40.0% at 50 mg, 47.9% at 100 mg and ",
    "52.9% at 200 mg. Parameter values are the frequentist (NONMEM) column ",
    "of Table 1; the paper's Bayesian (OpenBUGS) fit of the same model gave ",
    "closely agreeing fixed effects (E0 -2.51, Emax 2.95, ED50 18.18) but a ",
    "materially larger between-study SD (0.38). This model is one of two ",
    "landmark examples in the same tutorial; the companion is ",
    "modellib('Boucher_2016_naproxen_mbma')."
  )

  reference <- paste(
    "Boucher M, Bennetts M.",
    "The Many Flavors of Model-Based Meta-Analysis:",
    "Part I-Introduction and Landmark Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2016 Feb;5(2):54-64.",
    "doi:10.1002/psp4.12041.",
    "Structural model: 'Models for binary data' section (binomial likelihood",
    "and the displayed logit equation), and the Supplementary Materials",
    "NONMEM control stream ($PRED block of PSP4-5-54-s006.txt) and OpenBUGS",
    "model (PSP4-5-54-s005.txt), which agree exactly.",
    "Parameter values: Table 1, 'Frequentist approach' column.",
    "Included trials, arm sizes and event counts: Supplementary Materials",
    "Table 2, supplied as the dataset PSP4-5-54-s008.csv.",
    sep = " "
  )

  vignette <- "Boucher_2016_landmark_mbma"

  units <- list(
    time = paste0(
      "not applicable (LANDMARK meta-analysis: every modelled observation is ",
      "the paresthesia incidence accumulated over one whole trial, so the ",
      "model has no time term and its prediction is constant in time. The ",
      "source does not report the included trials' durations)"
    ),
    dosing = paste0(
      "mg/day (assigned daily topiramate dose for the study arm, supplied ",
      "through the DOSE_TPM_MGD covariate column and NOT as rxode2 dose ",
      "events; the arms studied were 0, 50, 100 and 200 mg/day)"
    ),
    concentration = paste0(
      "fraction/arm (probability that a subject in the study arm reports ",
      "paresthesia during the trial; the output prob_paresthesia is NOT a ",
      "drug concentration - the slash only satisfies checkModelConventions ",
      "unit parsing)"
    )
  )

  covariateData <- list(
    DOSE_TPM_MGD = list(
      description = paste0(
        "Assigned total daily topiramate dose for the study arm, in mg/day; ",
        "0 in the placebo arms."
      ),
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate here (a property of the trial arm), ",
        "in contrast to the per-patient per-record usage of the same ",
        "canonical in Lee_2024_topiramate.R. It is the dose-response driver ",
        "rather than a covariate on a structural parameter: it enters the ",
        "Emax term directly as Emax * DOSE / (ED50 + DOSE), so a value of 0 ",
        "makes that term exactly 0 and leaves the placebo intercept. ",
        "Supplied as a covariate column and not as an rxode2 amt / EVID = 1 ",
        "event because the model is purely algebraic with no PK compartment ",
        "(the DOSE_EE_UG / DOSE_AGT_UG dose-as-covariate precedent). The ",
        "meta-analysis observed only 0, 50, 100 and 200 mg/day, and only two ",
        "of the six trials contributed a 50 mg arm and two a 100 mg arm ",
        "(Supplementary Materials Table 2), so dose-response above 200 ",
        "mg/day is extrapolation. The source calls the column 'dose' in both ",
        "the NONMEM control stream ($INPUT ID DOSE NTOT DV) and the ",
        "OpenBUGS data list; the canonical name is used here because ",
        "rxode2::etTrans consumes a column named bare 'DOSE' before model() ",
        "sees it."
      ),
      source_name = "dose (Boucher 2016 Supplementary Materials Table 2 / PSP4-5-54-s008.csv)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1650L,
    n_studies = 6L,
    n_arms = 17L,
    disease_state = paste0(
      "Adults with episodic migraine receiving topiramate for migraine ",
      "prophylaxis. The modelled endpoint is a SAFETY endpoint - the ",
      "incidence of paresthesia, topiramate's commonest dose-limiting ",
      "adverse event - not an efficacy endpoint."
    ),
    dose_range = "0 (placebo), 50, 100 and 200 mg/day topiramate",
    design = paste0(
      "Six randomized placebo-controlled episodic migraine prophylaxis ",
      "trials, pooled at the study-arm level: Edwards 2003 (0/200 mg, ",
      "n = 15/15), Silberstein 2004 (0/50/100/200 mg, n = 116/118/126/113), ",
      "Brandes 2004 (0/50/100/200 mg, n = 113/117/119/117), Diener 2004 ",
      "(0/100/200 mg, n = 143/141/144), Storey 2001 (0/200 mg, n = 21/19) ",
      "and Silberstein 2006 (0/200 mg, n = 73/140). Every trial contributed ",
      "a placebo arm and a 200 mg arm; only Silberstein 2004 and Brandes ",
      "2004 contributed 50 mg arms."
    ),
    notes = paste0(
      "MBMA AT THE STUDY-ARM LEVEL: each modelled data point is one trial ",
      "arm's paresthesia event count out of that arm's subject total, and ",
      "the model is intended for simulating study-arm incidence rates - it ",
      "is NOT suitable for individual-subject simulation. n_subjects is the ",
      "sum of the 17 arm sizes in Supplementary Materials Table 2 ",
      "(15+15+116+118+126+113+113+117+119+117+143+141+144+21+19+73+140). ",
      "The observed placebo paresthesia rate ranged from 4.4% (Brandes ",
      "2004) to 26.7% (Edwards 2003); the two extreme placebo arms are the ",
      "two smallest trials (Edwards 2003 n = 15 and Storey 2001 n = 21), ",
      "which the source notes 'are less influential than they would be with ",
      "large numbers of subjects' and which are the reason the Bayesian ",
      "between-study SD (0.38) is more than twice the frequentist one ",
      "(0.17). No covariates other than dose were tested; the source ",
      "states only that there 'may be covariates that explain some of this ",
      "variability'. Demographics (age, sex, weight) are not reported for ",
      "the pooled cohort because the meta-analysis was built from published ",
      "summary incidence tables."
    )
  )

  ini({
    # ========================================================================
    # STRUCTURE. The NONMEM $PRED block and the OpenBUGS model in the
    # Supplementary Materials specify the identical model:
    #
    #   NONMEM   E0  = THETA(1) + ETA(1)
    #            LGT = E0 + EMAX*DOSE/(ED50 + DOSE)
    #            PR1 = 1/(1+EXP(-LGT))
    #   OpenBUGS logit(p.col[i]) <- mu + ((Emax*dose[i])/(ED50 + dose[i]))
    #                               + eta1[ID[i]]
    #            col[i] ~ dbin(p.col[i], n[i])
    #
    # i.e. a study-level random intercept on the placebo log-odds plus a
    # three-parameter Emax function of the assigned daily dose. There is no
    # Hill term: the main text says the sigmoidal Emax model "adds a fourth,
    # Hill, parameter" and this example did not use it.
    #
    # VALUES are Table 1, 'Frequentist approach (95% confidence interval)'
    # column, i.e. the NONMEM fit. The Bayesian column of the same table is
    # quoted in each comment for contrast but is NOT what is encoded.
    # ========================================================================

    e0 <- -2.56
    label("Placebo log-odds of paresthesia in a typical study (logit units); expit(-2.56) = 7.2% paresthesia at dose 0")
    # Table 1, row E0, frequentist: -2.56 (-3.03, -2.09). Bayesian: -2.51
    # (-3.06, -1.90). The CI half-width 0.47 implies SE 0.24.

    emax <- 2.91
    label("Maximal topiramate effect on the log-odds of paresthesia (logit units); the asymptotic rate is expit(-2.56 + 2.91) = 58.7%")
    # Table 1, row E max, frequentist: 2.91 (2.56, 3.26). Bayesian: 2.95
    # (2.50, 3.44). Cross-checked against the NONMEM variance-covariance
    # matrix quoted in the Supplementary Materials R script
    # (PSP4-5-54-s005.txt): Var(Emax) = 0.0302, so SE = 0.174 and the Wald
    # half-width is 1.96 * 0.174 = 0.34, matching the printed 0.35.

    ed50 <- 17.5
    label("Daily topiramate dose giving half the maximal effect on the log-odds of paresthesia (mg/day)")
    # Table 1, row ED50, mg, frequentist: 17.5 (11.64, 23.36). Bayesian:
    # 18.18 (6.08, 36.46), a much wider credible interval. Cross-checked
    # against the NONMEM variance-covariance matrix quoted in
    # PSP4-5-54-s005.txt: Var(ED50) = 8.58, so SE = 2.93 and the Wald
    # half-width is 1.96 * 2.93 = 5.74, matching the printed 5.86.
    #
    # NATURAL SCALE, not log. Both source implementations estimated ED50 on
    # the natural scale under a positivity constraint (NONMEM
    # '$THETA (0, 20)'; OpenBUGS 'ED50 ~ dunif(0.0001, 1000)'), and Table 1
    # prints a symmetric natural-scale Wald interval, so an led50 <-
    # log(17.5) parameterisation would misrepresent the published
    # uncertainty. The published covariance cov(Emax, ED50) = 0.0453 is also
    # on this scale and is used by the vignette's delta-method check.

    # ---- Between-STUDY variability ------------------------------------------
    # 'Where E0i is the log odds of paresthesia on placebo for study i and
    # assumed to be normally distributed with mean E0 and variance s^2'
    # ('Models for binary data'). ONE random effect, on the placebo intercept
    # only - Emax and ED50 are common to every study. Named eta_study_* per
    # the package MBMA convention so it cannot be mistaken for a
    # between-SUBJECT eta; this model has no between-subject variability at
    # all, and the source explicitly contrasts its placement 'on the overall
    # function, rather than the on-treatment effect ... in contrast to the
    # first WOMAC pain example'.
    eta_study_e0 ~ 0.0289
    # Table 1, row s, frequentist: 0.17 (0.04, 0.30). ini() takes the
    # VARIANCE, so 0.17^2 = 0.0289.
    #
    # SD, NOT VARIANCE - three independent confirmations:
    #  (1) The OpenBUGS code defines the same quantity as a standard
    #      deviation: 'tau.e0 ~ dunif(0,10)' with
    #      'prec.e0 <- 1/(tau.e0*tau.e0)'. Table 1 reports 's' for both
    #      columns as one quantity, and the main text compares them directly
    #      ('s, which is estimated as 0.17 in the frequentist results and
    #      0.38 using the Bayesian approach'), so the frequentist entry is on
    #      the same SD scale.
    #  (2) Data check. The four LARGE trials' observed placebo logits are
    #      -2.603, -3.073, -2.702 and -2.847 (Supplementary Materials Table
    #      2), whose sample SD is 0.204 - close to 0.17 once binomial
    #      sampling error is removed, and nowhere near the 0.41 that reading
    #      0.17 as a variance would imply.
    #  (3) The printed interval (0.04, 0.30) is symmetric about 0.17 with
    #      half-width 0.13, the shape of a NONMEM standard error reported on
    #      the estimated scale.

    # ---- Residual error ------------------------------------------------------
    addSd_prob_paresthesia <- fixed(0.001)
    label("Placeholder additive residual SD on the study-arm paresthesia probability; the source likelihood is binomial on the arm event count, so there is no source residual")
    # NOT FROM SOURCE. Both source implementations use a binomial likelihood
    # on the arm's event count out of its subject total ('Yij ~
    # Binomial(Nij, pij)'; NONMEM '-2LL' with an explicit binomial Y;
    # OpenBUGS 'col[i] ~ dbin(p.col[i], n[i])'), so arm-level sampling error
    # is a property of the arm size N rather than an estimated parameter.
    # A downstream simulation that needs the arm's observed rate should draw
    # rbinom(1, N_arm, prob_paresthesia) / N_arm rather than use this term.
    # See the vignette Assumptions and deviations section.
  })

  model({
    # ---- Logit-scale dose response ------------------------------------------
    # The between-study eta is added on the logit scale, before the
    # transform, and sits on the placebo intercept only.
    logit_paresthesia <-
      e0 + eta_study_e0 +
      emax * DOSE_TPM_MGD / (ed50 + DOSE_TPM_MGD)

    # ---- Probability of paresthesia -----------------------------------------
    # The study-arm paresthesia incidence rate, in [0, 1]. NOT an individual
    # patient's risk - see description and population$notes.
    prob_paresthesia <- expit(logit_paresthesia)

    # ---- Observation ---------------------------------------------------------
    prob_paresthesia ~ add(addSd_prob_paresthesia)
  })
}
