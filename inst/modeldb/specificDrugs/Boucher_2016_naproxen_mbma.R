Boucher_2016_naproxen_mbma <- function() {
  description <- paste0(
    "MBMA. Landmark random-effects meta-regression of the week-2 treatment ",
    "difference in Western Ontario and McMaster Universities (WOMAC) pain ",
    "score (0-20 scale) between naproxen 500 mg twice daily and placebo in ",
    "osteoarthritis of the knee or hip, fitted to the single week-2 ",
    "summary estimate from each of 13 randomized double-blind ",
    "placebo-controlled parallel-group trials (10 internal unpublished, 3 ",
    "published; 9 flare designs, 4 non-flare). The modelled quantity is the ",
    "already-differenced naproxen-minus-placebo mean change from baseline, ",
    "so the model is a study-level linear predictor with no dose, no time ",
    "course and no placebo arm of its own: difference = E0 + dE0_flare * ",
    "FLARE + eta_study. Flare design is carried as a study-level covariate ",
    "but was NOT statistically significant (-0.14, 95% CI -0.47 to 0.18), ",
    "which was the source's answer to its own research question about ",
    "whether to run a flare design. Variability is BETWEEN-STUDY only ",
    "(one random effect, SD 0.1371 WOMAC units), and the residual is the ",
    "study's own reported standard error supplied with the data rather ",
    "than an estimated parameter, so sigma is fixed to 1 and reweighted ",
    "per study downstream. The simulation scope is a future trial's ",
    "study-level treatment difference, NOT an individual patient's pain ",
    "score. LANDMARK, single timepoint: for the longitudinal time-course ",
    "Emax model of the same endpoint by the same authors (18 trials, WOMAC ",
    "on a 0-10 scale, study-ARM means rather than differences) see ",
    "modellib('Boucher_2018_naproxen_mbma'). The companion landmark model ",
    "from this same tutorial is modellib('Boucher_2016_topiramate_mbma')."
  )

  reference <- paste(
    "Boucher M, Bennetts M.",
    "The Many Flavors of Model-Based Meta-Analysis:",
    "Part I-Introduction and Landmark Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2016 Feb;5(2):54-64.",
    "doi:10.1002/psp4.12041.",
    "Structural model: 'Model descriptions' section, 'Meta-regression",
    "model' paragraph, and the Supplementary Materials OpenBUGS",
    "meta-regression model (linebugs_MR in PSP4-5-54-s002.txt).",
    "Parameter values: the frequentist REML metafor fit printed by the",
    "authors in the Supplementary Materials R script PSP4-5-54-s002.txt",
    "('# -0.92 (-1.20, -0.63) / # FLARE: -0.14 (-0.47, 0.18) /",
    "# tau: 0.1371'); the flare estimate is also quoted in the main text",
    "'WOMAC PAIN RESULTS' section.",
    "Included trials, treatment differences, standard errors and flare",
    "flags: Supplementary Materials Table 1, supplied as the dataset",
    "PSP4-5-54-s007.csv.",
    sep = " "
  )

  vignette <- "Boucher_2016_landmark_mbma"

  units <- list(
    time = paste0(
      "not applicable (LANDMARK meta-analysis: the single modelled ",
      "observation per trial is the treatment difference at week 2, so the ",
      "model has no time term and its prediction is constant in time)"
    ),
    dosing = paste0(
      "not applicable (every active arm in every included trial used the ",
      "same regimen, naproxen 500 mg twice daily, so no dose-response is ",
      "identifiable and treatment is implicit in the endpoint - the ",
      "modelled quantity is already a naproxen-minus-placebo difference. ",
      "There are no rxode2 dose events)"
    ),
    concentration = paste0(
      "WOMAC pain units/study (the output Cc is one trial's ",
      "naproxen-minus-placebo difference in mean change from baseline in ",
      "WOMAC pain at week 2, on the 0-20 total scale, so a NEGATIVE value ",
      "is a benefit; it is NOT a drug concentration - the slash only ",
      "satisfies checkModelConventions unit parsing)"
    )
  )

  covariateData <- list(
    FLARE = list(
      description = paste0(
        "Study-level flare-design indicator: 1 if the trial used a flare ",
        "design (subjects were washed out of their pain medications and ",
        "were required to have a predefined increase in pain before being ",
        "eligible for randomization), 0 for a non-flare design."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-flare design)",
      notes = paste0(
        "MBMA study-LEVEL covariate (a property of the trial design, not of ",
        "an individual patient). Same canonical and same role as in the ",
        "sibling Boucher_2018_naproxen_mbma.R, but note the different ",
        "target: there flare shifts the study-ARM baseline E0 and Emax ",
        "separately, whereas here the modelled quantity is already a ",
        "within-trial difference, so the single coefficient is flare's ",
        "effect on the naproxen-placebo CONTRAST. Encoded as 'Fi' in the ",
        "source's meta-regression equation (theta_i = theta_i + theta_FL * ",
        "Fi + e_i, 'where Fi is 0 for non-flare and 1 for flare') and as ",
        "'Flare' in the supplied dataset. Of the 13 trials, 9 were flare ",
        "designs and 4 were non-flare. THE EFFECT WAS NOT SIGNIFICANT: the ",
        "source reports -0.14 (95% CI -0.47, 0.18) frequentist and -0.13 ",
        "(95% CrI -0.47, 0.24) Bayesian, and concludes 'flare did not seem ",
        "to be a significant covariate'. It is retained in the model file ",
        "because answering that question was one of the paper's two stated ",
        "research questions and because the coefficient is a published ",
        "point estimate, not because the source recommends using it."
      ),
      source_name = "Flare (Boucher 2016 Supplementary Materials Table 1 / PSP4-5-54-s007.csv); Fi in the 'Meta-regression model' equation"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 13L,
    disease_state = paste0(
      "Adults with osteoarthritis pain of the knee or hip. The endpoint is ",
      "the WOMAC pain subscale, five pain-related questions each scored 0 ",
      "(no pain) to 4 (maximum pain) and summed to a total between 0 and ",
      "20, analysed as a continuous variable."
    ),
    dose_range = "naproxen 500 mg twice daily vs placebo (a single regimen; no dose-ranging)",
    design = paste0(
      "13 double-blind, placebo-controlled, randomized parallel-group ",
      "trials in which both a naproxen 500 mg twice-daily arm and a placebo ",
      "arm were present. 9 used a flare design and 4 did not. 10 of the 13 ",
      "were internal unpublished Pfizer trials; the 3 published ones are ",
      "Schnitzer et al. 2005, Baerwald et al. 2010 and Schnitzer et al. ",
      "2010. The week-2 arithmetic means were the ones analysed."
    ),
    timepoints = "week 2 only (landmark analysis of a single timepoint)",
    notes = paste0(
      "MBMA AT THE STUDY LEVEL: each modelled data point is one trial's ",
      "naproxen-minus-placebo difference in mean change from baseline in ",
      "WOMAC pain at week 2, together with that difference's reported ",
      "standard error. The model is intended for simulating a future ",
      "trial's treatment difference - for example to set a proof-of-concept ",
      "target value, the source's first research question - and is NOT ",
      "suitable for individual-subject simulation. n_subjects is NA because ",
      "the source reports only per-trial treatment differences and their ",
      "standard errors, never the arm sizes. The 13 observed differences ",
      "range from -1.4776 to -0.49 WOMAC units with standard errors from ",
      "0.153 to 0.431 (Supplementary Materials Table 1). ",
      "OTHER MODELS FITTED TO THE SAME 13 TRIALS, reported by the source ",
      "for comparison but not encoded here because the meta-regression ",
      "nests them: the fixed-effects pooled estimate is -1.025 (Bayesian ",
      "posterior SD 0.057), and the random-effects pooled estimate without ",
      "the flare covariate is -1.0232 (95% CI -1.1652 to -0.8811, 95% ",
      "prediction interval -1.3494 to -0.6970) frequentist REML and ",
      "-1.0245 (posterior SD 0.0792, tau 0.1661) Bayesian. The source's ",
      "closing recommendation is that a random-effects model be the ",
      "starting point and that the PREDICTION interval, not the confidence ",
      "interval, be used when designing a future trial. Publication bias ",
      "was assessed and not detected (Begg and Mazumdar rank correlation ",
      "P = 0.95, Egger regression P = 0.99), though the source cautions ",
      "this is hard to judge when 10 of 13 trials are unpublished internal ",
      "studies."
    )
  )

  ini({
    # ========================================================================
    # STRUCTURE ('Model descriptions' / 'Meta-regression model'):
    #
    #   theta_hat_i = theta_i + theta_FL * Fi + e_i
    #   theta_i ~ N(delta, s^2)          between-study random effect
    #   e_i     ~ N(0, nu_i^2)           within-study sampling error, where
    #                                    nu_i is the study's OWN reported SE
    #                                    and is treated as KNOWN
    #
    # The equivalent OpenBUGS code in the Supplementary Materials is
    #   mu[i]    <- delta[i] + FLCOV*Flare[i]
    #   delta[i] ~ dnorm(d, prec.bsv)
    #   y[i]     ~ dnorm(mu[i], prec.y[i]);  prec.y[i] <- 1/(se_wp[i]^2)
    #
    # so 'd' is the mean treatment difference in NON-FLARE trials and FLCOV
    # is the additional effect of a flare design.
    #
    # VALUES are the frequentist REML fit (metafor rma with mods = Flare)
    # whose output the authors print verbatim in the Supplementary Materials
    # R script PSP4-5-54-s002.txt. Bayesian counterparts are quoted per line.
    # ========================================================================

    e0 <- -0.92
    label("Mean naproxen-minus-placebo difference in week-2 WOMAC pain in a typical NON-FLARE trial (WOMAC units, 0-20 scale; negative = pain reduction)")
    # PSP4-5-54-s002.txt, printed metafor meta-regression output:
    # '# -0.92 (-1.20, -0.63)'. This is the meta-regression INTERCEPT, i.e.
    # the non-flare mean - not the overall pooled estimate, which is -1.0232
    # (see population$notes). Data check: the 4 non-flare trials' differences
    # are -1, -0.49, -1.4776 and -1.1, whose inverse-variance weighted mean
    # is -0.88; REML shrinks toward the unweighted mean (-1.02), landing at
    # -0.92.

    e_flare_e0 <- -0.14
    label("Additional naproxen-minus-placebo difference in a FLARE-design trial (WOMAC units); the flare-design mean difference is -0.92 + -0.14 = -1.06")
    # PSP4-5-54-s002.txt, printed metafor meta-regression output:
    # '# FLARE: -0.14 (-0.47, 0.18)'; also quoted in the main text 'WOMAC
    # PAIN RESULTS' section, alongside the Bayesian estimate -0.13 (95% CrI
    # -0.47, 0.24). NOT SIGNIFICANT - the interval spans 0 - and the source
    # concludes 'flare did not seem to be a significant covariate'. Retained
    # as a published point estimate; see covariateData$FLARE$notes. Data
    # check: the 9 flare trials' differences average -1.03 against the 4
    # non-flare trials' -1.02, so the modest fitted separation comes from
    # the inverse-variance weighting, not from a visible group difference.

    # ---- Between-STUDY variability ------------------------------------------
    # 'Note that s^2 is the between-study variance of theta_i'
    # ('Random effects model'). ONE random effect on the treatment
    # DIFFERENCE. Named eta_study_* per the package MBMA convention so it
    # cannot be mistaken for a between-SUBJECT eta; this model has no
    # between-subject variability at all. The source is explicit that
    # because a within-trial contrast is being pooled, 'the main effect of
    # trial has been eliminated and so this between-study variability
    # reflects the treatment-by-study interaction'.
    eta_study_e0 ~ 0.01879641
    # PSP4-5-54-s002.txt, printed metafor meta-regression output:
    # '# tau: 0.1371'. metafor's 'tau' is the between-study standard
    # deviation (its square, tau2, is the variance), and ini() takes the
    # VARIANCE, so 0.1371^2 = 0.01879641. This is the RESIDUAL heterogeneity
    # after the flare covariate; the flare-free random-effects fit of the
    # same data gave a Bayesian posterior tau of 0.1661.

    # ---- Residual error ------------------------------------------------------
    addSd <- fixed(1)
    label("Residual standard deviation multiplier; the per-study residual SD is addSd * SE_i, the study's own reported standard error of its treatment difference")
    # NOT ESTIMATED. The source treats each study's within-study variance as
    # KNOWN and equal to the reported squared standard error ('Usually, we
    # treat nu_i^2 as known and equal to the estimated variance of
    # theta_hat_i'), which the OpenBUGS code implements as
    # 'prec.y[i] <- 1/(se_wp[i] * se_wp[i])'. SE_i is therefore data, not a
    # parameter, and sigma is fixed to 1 so that the operative residual SD
    # is exactly SE_i.
    #
    # PER-STUDY REWEIGHTING. SE_i is a property of the trial being
    # simulated rather than of the model, so the 1 -> SE_i rescaling is left
    # to the user, exactly as in the sibling Boucher_2018_naproxen_mbma.
    # The model as shipped emits the SE = 1 residual; MULTIPLY THE SIMULATED
    # RESIDUAL BY SE_i, or equivalently set addSd to SE_i, to reproduce a
    # given trial. The 13 observed standard errors (0.153 to 0.431) are in
    # Supplementary Materials Table 1. The vignette demonstrates the
    # rescaling.
  })

  model({
    # Study-level covariate (binary): FLARE (1 = flare design, 0 = non-flare).
    #
    # Linear predictor for one trial's naproxen-minus-placebo difference in
    # week-2 WOMAC pain. There is no dose term (every trial used naproxen
    # 500 mg twice daily) and no time term (single landmark timepoint), so
    # the prediction is constant in time by construction.
    #
    # Cc is the package's canonical single-output observation variable
    # (R/conventions.R observationVar); it is a WOMAC pain DIFFERENCE on the
    # 0-20 scale, NOT a concentration - same naming convention as the
    # sibling Boucher_2018_naproxen_mbma. Negative values are a benefit.
    Cc <- e0 + e_flare_e0 * FLARE + eta_study_e0

    Cc ~ add(addSd)
  })
}
