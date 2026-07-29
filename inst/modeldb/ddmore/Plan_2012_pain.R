Plan_2012_pain <- function() {
  description <- "Markov Integer Model for placebo time-course of Likert (0-10) pain scores in adults; pooled placebo arm of three Phase III neuropathic-pain trials (Plan 2012; DDMODEL00000194)"
  reference <- paste(
    "Plan EL, Elshoff JP, Stockis A, Sargentini-Maier ML, Karlsson MO. (2012).",
    "Likert pain score modeling: a Markov integer model and an autoregressive",
    "continuous model. Clin Pharmacol Ther 91(4):820-828.",
    "doi:10.1038/clpt.2011.301.",
    "DDMORE Foundation Model Repository: DDMODEL00000194."
  )
  vignette <- "Plan_2012_pain"
  units <- list(
    time = "day",
    dosing = "(none; placebo arm)",
    concentration = "(11-point Likert pain score, 0-10, unitless)"
  )
  ddmore_id <- "DDMODEL00000194"
  replicate_of <- NULL

  covariateData <- list(
    CONMED_PARA = list(
      description        = "Concomitant paracetamol (acetaminophen) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant paracetamol)",
      notes              = "Source column PCM. Acts as additive shift on logit(lambda/10): PHL = logit(TVLAM) + e_conmed_para * CONMED_PARA. New canonical entry registered in inst/references/covariate-columns.md alongside this model.",
      source_name        = "PCM"
    )
  )

  population <- list(
    n_subjects     = 231,
    n_studies      = 3,
    age_range      = "(not extracted; Plan 2012 publication not on disk for cross-check)",
    weight_range   = "(not extracted)",
    sex_female_pct = "(not extracted)",
    disease_state  = "Neuropathic pain (placebo arm)",
    dose_range     = "(placebo arm only; no drug exposure)",
    notes          = "Pooled placebo arm of three Phase III neuropathic-pain trials. Daily 11-point Likert pain measurements over 18 weeks; 22,492 measurements total. Demographic detail (age, weight, sex split) not derivable from the DDMORE bundle; the linked Plan 2012 publication (doi:10.1038/clpt.2011.301) was not on disk in /home/bill/github/mab_human_consensus/literature for cross-check at extraction time. n_subjects (231) is taken from the DDMORE RDF model-has-description-long field."
  )

  ini({
    # All values are taken from the .mod $THETA / $OMEGA blocks, which equal
    # the FINAL PARAMETER ESTIMATE in Output_real_likert_pain_count.lst because
    # the .mod was run with $ESTIMATION MAXEVAL=0 — i.e. NONMEM evaluates the
    # objective at the supplied THETA/OMEGA without estimating. The .mod
    # therefore carries the publication's final estimates as its initials.
    # See ddmore-source.md for the convention.

    logitbas <- log(0.620667 / (1 - 0.620667))
    label("Logit of typical baseline pain score / 10 (PHI in source; back-transforms via 10*expit() to typical baseline 6.21 on 0-10 scale)")
    # .mod $THETA(1) = 6.20667 (BASELINE, bounded 0-10); .lst final TH 1 = 6.21E+00.

    logitpef <- log(0.18986 / (1 - 0.18986))
    label("Logit of maximum placebo effect (typical 0.190 fraction of baseline)")
    # .mod $THETA(2) = 0.18986 (PLACEBO_EFFECT, 0-1); .lst final TH 2 = 1.90E-01.

    lpha <- log(27.7045)
    label("Log of placebo half-time (days)")
    # .mod $THETA(3) = 27.7045 (PLACEBO_HALF_TIME); .lst final TH 3 = 2.77E+01.

    logitpi00 <- log(0.554617 / (1 - 0.554617))
    label("Logit of Markov stay-probability baseline (unitless logit), conditional on previous score = 0")
    # .mod $THETA(4) = 0.554617 (PROBABILITY_OF_INFLATION_0/0); .lst final TH 4.

    logitpi09 <- log(0.119517 / (1 - 0.119517))
    label("Logit of Markov stay-probability baseline (unitless logit), conditional on previous score = 9")
    # .mod $THETA(5) = 0.119517 (PROBABILITY_OF_INFLATION_0/9); .lst final TH 5.

    logitpi10 <- log(0.443759 / (1 - 0.443759))
    label("Logit of Markov stay-probability baseline (unitless logit), conditional on previous score = 10")
    # .mod $THETA(6) = 0.443759 (PROBABILITY_OF_INFLATION_0/10); .lst final TH 6.

    logitpi1 <- log(0.359302 / (1 - 0.359302))
    label("Logit of Markov inflation probability (unitless logit) for |score - previous| = 1")
    # .mod $THETA(7) = 0.359302 (PROBABILITY_OF_INFLATION_1); .lst final TH 7.

    logitpi2 <- log(0.00472972 / (1 - 0.00472972))
    label("Logit of Markov inflation probability (unitless logit) for |score - previous| = 2")
    # .mod $THETA(8) = 0.00472972 (PROBABILITY_OF_INFLATION_2); .lst final TH 8.

    logitpi3 <- log(0.000403033 / (1 - 0.000403033))
    label("Logit of Markov inflation probability (unitless logit) for |score - previous| = 3")
    # .mod $THETA(9) = 0.000403033 (PROBABILITY_OF_INFLATION_3); .lst final TH 9.

    logitdis <- log(0.99286 / (1 - 0.99286))
    label("Logit of underdispersion parameter DIS (0-1)")
    # .mod $THETA(10) = 0.99286 (DIS); .lst final TH 10.

    lte0 <- log(0.00643534)
    label("Log of Markov transition-element drift TE0 (1/day)")
    # .mod $THETA(11) = 0.00643534 (TE0); .lst final TH 11.

    e_conmed_para <- 0.36374
    label("Additive shift on logit(lambda/10) when CONMED_PARA = 1")
    # .mod $THETA(12) = 0.36374 (COV); .lst final TH 12.

    # IIV variances (parameter scale = transformed scale).
    etalogitbas ~ 0.568985
    # .mod $OMEGA(1) on ETA(1); .lst diag OMEGA ETA1 = 5.69E-01.
    etalogitpef ~ 3.77567
    # .mod $OMEGA(2) on ETA(2); .lst diag OMEGA ETA2 = 3.78E+00.
    etalpha ~ 0.352913
    # .mod $OMEGA(3) on ETA(3); .lst diag OMEGA ETA3 = 3.53E-01.

    # NONMEM $OMEGA BLOCK(3) on ETA(4), ETA(5), ETA(6) — correlated IIV across
    # the three Markov-inflation logits {pi00 / pi1 / pi2}. NONMEM ETA(4) is
    # SHARED across pi00, pi09, pi10 in $PRED (same eta term added to each
    # baseline logit); we attach the shared eta to logitpi00 here and replicate
    # it onto logitpi09 / logitpi10 inside model().
    etalogitpi00 + etalogitpi1 + etalogitpi2 ~ c(2.69723,
                                                 2.45376, 3.56524,
                                                 -0.755447, 0.805685, 1.92485)
    # .mod $OMEGA BLOCK(3) on ETA(4..6); .lst OMEGA OM44/OM55/OM66 diags and
    # OM45 / OM46 / OM56 covariances.

    etalogitpi3 ~ fix(0)
    # .mod $OMEGA(7) = 0 FIX on ETA(7) (no IIV on logitpi3).

    etalogitdis ~ 24.3896
    # .mod $OMEGA(8) on ETA(8); .lst diag OMEGA ETA8 = 2.44E+01.

    etalte0 ~ 1.33918
    # .mod $OMEGA(9) on ETA(9); .lst diag OMEGA ETA9 = 1.34E+00.
  })

  model({
    # Time variable: source uses TIME (h) and a derived DAYS = TIME/24 column;
    # this model declares units$time = "day", so `time` here is days directly.
    days <- time

    # Individual baseline pain score (0-10), placebo effect (0-1), placebo
    # half-time (d). `expit(x) = 1/(1 + exp(-x))` is rxode2's logit-inverse.
    bas <- 10 * expit(logitbas + etalogitbas)
    pef <-       expit(logitpef + etalogitpef)
    pha <- exp(lpha + etalpha)

    # Placebo decay: fraction of baseline at time t (1 at t=0; (1 - PEF) as t -> Inf).
    plc <- 1 - pef * (1 - exp(-log(2) / pha * days))

    # Mean count lambda on 0-10 scale, with the paracetamol additive logit shift.
    tvlam <- bas * plc / 10
    phl   <- log(tvlam / (1 - tvlam)) + e_conmed_para * CONMED_PARA
    lam   <- 10 * exp(phl) / (1 + exp(phl))

    # Time-varying Markov drift (TE * days; TEF in source).
    te  <- exp(lte0 + etalte0)
    tef <- te * days

    # Markov / inflated-truncated-Poisson likelihood structure (documented for
    # source-trace fidelity; not exercised by the simplified observation model
    # below). The full publication likelihood is a truncated Poisson with
    # underdispersion `dis` and Markov inflations `pi0` / `pi1` / `pi2` / `pi3`
    # conditional on the previous score. nlmixr2 / rxode2 do not natively
    # express that joint distribution, so the formal observation here is a
    # plain Poisson on `lam`. See vignette "Assumptions and deviations" for
    # the deviation rationale.
    pi00 <- expit(logitpi00 + etalogitpi00 + tef)
    pi09 <- expit(logitpi09 + etalogitpi00 + tef)
    pi10 <- expit(logitpi10 + etalogitpi00 + tef)
    pi1  <- 0.5 * expit(logitpi1 + etalogitpi1)
    pi2  <- 0.5 * expit(logitpi2 + etalogitpi2)
    pi3  <- 0.5 * expit(logitpi3 + etalogitpi3)
    dis  <- expit(logitdis + etalogitdis)

    # Observation: typical-value mean pain score (0-10), modelled as Poisson(lam)
    # for nlmixr2 fitting compatibility. F.3 mechanistic-sanity validation in
    # the vignette compares the typical-value `score` trajectory against the
    # placebo-decay analytic form Plan 2012 reports.
    score <- lam
    score ~ pois(lam)
  })
}
