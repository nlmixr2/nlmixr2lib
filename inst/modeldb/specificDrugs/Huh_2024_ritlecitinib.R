Huh_2024_ritlecitinib <- function() {
  description <- paste0(
    "Longitudinal exposure-response model of the Severity of Alopecia ",
    "Tool (SALT) score during ritlecitinib treatment in adolescents and ",
    "adults with alopecia areata, pooled across studies B7931005 (phase ",
    "2a), B7981015 (phase 2b/3) and B7981032 (phase 3 long-term). SALT ",
    "is a continuous bounded outcome (CBO) on 0-100 (100 = complete ",
    "scalp hair loss); Huh 2024 handles the boundary mass with the ",
    "Hutmacher censoring approach after an Aranda-Ordaz transformation ",
    "of the non-boundary data, x = log(((100 - SALT)/100)^(-alpha) - 1) ",
    "- log(alpha)). The transformed score is described by a latent-",
    "variable indirect-response model: the conditional mean is ",
    "mu = BASE - (PBO(t) - 1) - (E(t) - 1), where PBO(t) is a placebo ",
    "latent variable driven by an on-treatment indicator through a ",
    "chain of three transit compartments, and E(t) is a drug latent ",
    "variable driven by an Emax/EC50 function of the average ",
    "ritlecitinib concentration between adjacent SALT records through a ",
    "chain of two transit compartments. Each chain has kin = kout, so ",
    "all latent states start at kin/kout = 1 and the paper's mean ",
    "transit time is (number of transits + 1)/kout. Disease severity ",
    "(alopecia totalis / universalis) is a structural covariate on ",
    "BASE, Pmax and the drug-chain half-life; study B7981032 shifts ",
    "BASE in the non-AT/AU group only. No covariate survived the ",
    "stepwise screen. There is no PK component -- exposure enters as ",
    "the CAV column, the empirical-Bayes average concentration from the ",
    "companion ritlecitinib population PK model of Wojciechowski 2023 ",
    "(modellib('Wojciechowski_2023_ritlecitinib_final'))."
  )

  reference <- paste(
    "Huh Y, Wojciechowski J, Purohit VS.",
    "Moving Beyond Boundaries: Utilization of Longitudinal",
    "Exposure-Response Model for Bounded Outcome Score to Inform",
    "Decision Making in the Accelerated Drug Development Paradigm.",
    "Clin Pharmacokinet 2024;63(3):381-394.",
    "doi:10.1007/s40262-024-01347-6.",
    sep = " "
  )

  vignette <- "Huh_2024_ritlecitinib"

  # There is no PK component: ritlecitinib exposure must be supplied from
  # outside the model as the per-interval average concentration CAV. Huh 2024
  # took it from the empirical Bayes estimates of the companion population PK
  # model, packaged here as `Wojciechowski_2023_ritlecitinib_final`.
  depends <- "CAV"

  # Paper-mechanistic latent states and outputs.
  #
  # `placebo1` .. `placebo4` are the placebo latent variable PBO(t) of
  # Huh 2024 Eq 4 preceded by the three transit compartments the final
  # placebo model added (Section 3.2). `placebo4` IS PBO(t); the mean
  # transit time the paper quotes, (number of transits + 1)/kout, counts
  # all four states. There is no canonical compartment for a placebo-
  # response delay chain, so the chain is declared paper-specific rather
  # than overloading `lat<n>` (a physiologic latency) or `transit<n>` (an
  # absorption chain). The drug-side chain uses the blessed `effect<n>`
  # prefix, with `effect3` = E(t).
  #
  # `salt_transformed` is the Aranda-Ordaz-transformed SALT score on
  # which the model was fit and on which residual error acts (Eq 2);
  # `salt` is its back-transform onto the published 0-100 SALT scale.
  # The canonical PK output `Cc` does not apply -- there is no drug
  # concentration state in this model.
  paper_specific_compartments <- c(
    "placebo1", "placebo2", "placebo3", "placebo4",
    "salt_transformed", "salt"
  )

  # `etapmax` is additive on the linear-scale Pmax (Section 2.3: "Pmax
  # with an additive model (Pi = Ppop + eta_i) to allow both disease
  # worsening and improving"), so the canonical `etal<name>` + `l<name>`
  # log-normal pairing does not apply. Same pattern as
  # `Wojciechowski_2015_rheumatoidArthritis`'s `etabase` / `etaex1`.
  paper_specific_etas <- "etapmax"

  # Every state is a dimensionless latent variable, not a measurable
  # analyte: the two chains carry the placebo and drug arms of the Eq 4
  # indirect-response system on the Aranda-Ordaz-transformed SALT
  # domain. Each starts at its untreated steady state of kin/kout = 1.
  compartmentData <- list(
    placebo1 = list(analyte = "placebo latent variable (transit 1 of 3)", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    placebo2 = list(analyte = "placebo latent variable (transit 2 of 3)", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    placebo3 = list(analyte = "placebo latent variable (transit 3 of 3)", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    placebo4 = list(analyte = "placebo latent variable PBO(t) of Huh 2024 Eq 4", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    effect1  = list(analyte = "drug-effect latent variable (transit 1 of 2)", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    effect2  = list(analyte = "drug-effect latent variable (transit 2 of 2)", units = "(dimensionless)", specimen = "not applicable", verified = FALSE),
    effect3  = list(analyte = "drug-effect latent variable E(t) of Huh 2024 Eq 4", units = "(dimensionless)", specimen = "not applicable", verified = FALSE)
  )

  units <- list(
    time          = "week (time since the first study treatment record)",
    dosing        = "(none; no PK component -- ritlecitinib exposure enters as the CAV covariate column, not as dosing events)",
    concentration = "(Severity of Alopecia Tool score, 0-100 percent scalp hair loss; the fitted output salt_transformed is on the Aranda-Ordaz-transformed scale of Huh 2024 Eq 1, with the back-transformed salt also emitted)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Average ritlecitinib plasma concentration over the interval between the previous SALT record and the current SALT record.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per record. Huh 2024 Section 2.2: derived from the empirical Bayes estimates of the final ritlecitinib population PK model of Wojciechowski 2023 (this package's `Wojciechowski_2023_ritlecitinib_final`) combined with the patient's own dosing diary, then averaged over the inter-SALT-record interval. The averaging window is therefore the observation interval, NOT a dosing interval and NOT a steady-state Cavg computed from dose and clearance -- Huh 2024 Section 4.1 argues this choice specifically so that treatment interruptions are reflected in the exposure metric. Set to 0 for placebo periods and for any interval in which no ritlecitinib was taken; the Emax term then collapses to 0. For orientation, Huh 2024 Section 3.2 reports the Cavg of the 50 mg QD regimen as 52 ng/mL, essentially at the estimated EC50 of 53.6 ng/mL.",
      source_name        = "Cavg"
    ),
    TRT_PHASE = list(
      description        = "Indicator that the record falls within a period in which study treatment (ritlecitinib OR matching placebo) is being administered.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (off study treatment: pre-treatment, drug holiday, or after treatment withdrawal).",
      notes              = "Time-varying per record. This is the I_PBO indicator of Huh 2024 Eq 4, 'an indicator variable that equals 1 if treatment was given and equals 0 otherwise'. In this model it gates the PLACEBO latent chain only; the drug latent chain is gated separately by CAV = 0, so a placebo-arm record carries TRT_PHASE = 1 and CAV = 0. Setting TRT_PHASE = 0 (together with CAV = 0) is how the paper's treatment-interruption simulation (Section 2.6 / Table 3) switches both effects off. The canonical entry's description says the indicator switches drug AND placebo effects on; here only the placebo effect is gated by the column itself, which is a narrowing consistent with the canonical rather than a different concept.",
      source_name        = "I_PBO"
    ),
    DIS_ALOPECIA_AT_AU = list(
      description        = "Indicator that the participant's alopecia areata has progressed to alopecia totalis (complete scalp hair loss) or alopecia universalis (complete scalp and body hair loss), the severe end of the alopecia areata spectrum.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-AT/AU alopecia areata: patchy disease not meeting AT or AU criteria).",
      notes              = "Time-fixed per subject (disease classification at study entry). Huh 2024 Table 1: 532/1268 (42.0%) of the analysis population were AT/AU. AT/AU participants have a baseline SALT of 100 by definition of the classification, versus a mean of 74.0 (median 77.9) in the non-AT/AU group. Huh 2024 retained this as a STRUCTURAL covariate -- it is the only covariate in the final model, entering on three parameters: BASE (separate typical values 11.6 vs 1.92 on the transformed scale), Pmax (2.75 vs 0 fix) and the drug-chain half-life (a -60.1% fractional shift, 3.11 vs 7.80 weeks). Note that the B7981032 study effect on BASE applies to the non-AT/AU group only, because the AT/AU baseline is pinned at SALT = 100 regardless of study entry criteria.",
      source_name        = "AT/AU status"
    ),
    STUDY_B7981032 = list(
      description        = "Indicator that the SALT record comes from study B7981032 (NCT04006457), the phase 3 long-term ritlecitinib study.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (studies B7931005 and B7981015, whose baseline-SALT inclusion criterion was >= 50%).",
      notes              = "Per-subject time-fixed within this analysis. Huh 2024 Section 3.2: 'Study B7981032 effect on BASE was incorporated to address differences in inclusion criteria between studies (baseline SALT >= 50 for B7931005 and B7981015 vs baseline SALT >= 25 for B7981032)'. Table 1: 458/1268 (36.1%) of participants. The effect is a MULTIPLICATIVE fractional shift on the non-AT/AU BASE, confirmed against the paper's own worked value: 1.92 * (1 + (-0.645)) = 0.68, the figure Section 3.2 quotes for the non-AT/AU group in B7981032 (an additive reading would give 1.28 and is falsified).",
      source_name        = "B7981032"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female-sex indicator at baseline.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on BASE and Emax (Huh 2024 Section 2.3) but not retained: 'the forward addition step of covariate analysis did not identify any important covariates' (Section 3.2). Table 1: 805/1268 (63.5%) female."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on BASE, Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1: mean 70.7 kg (SD 17.6), median 68.4 kg (range 29.6-200.0)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on BASE, Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1: mean 33.8 years (SD 14.2), median 32 (range 12-72); 170/1268 (13.4%) were adolescents 12 to <18 years."
    ),
    RACE_ASIAN = list(
      description = "Asian-race indicator at baseline.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Race was screened on BASE and Emax but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1 reports White 889 (70.1%), Asian 286 (22.6%), African American 46 (3.6%), Other 47 (3.7%); no per-level coefficient is reported because no level entered the final model."
    ),
    PRIORTRT = list(
      description = "Indicator of prior pharmacological treatment for alopecia areata.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on BASE, Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1: Yes 499 (39.4%), No 22 (1.7%), Unknown 747 (58.9%) -- the large unknown stratum limits what this covariate could have resolved."
    ),
    DURF = list(
      description = "Duration of alopecia areata since first diagnosis.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1: mean 10.0 years (SD 10.5), median 6.67 (range 0.04-60.1)."
    ),
    DURC = list(
      description = "Duration of the current alopecia areata episode.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). Table 1: mean 3.22 years (SD 2.82), median 2.25 (range 0.02-29.5). The paper singles this one out in the Discussion: a >10-year current episode has been hypothesised to predict non-response, but B7981015 and B7981032 restricted the current episode to <= 10 years, which the authors suggest is why no trend was visible."
    ),
    BL_SALT = list(
      description = "Baseline Severity of Alopecia Tool score.",
      units       = "percent scalp hair loss (0-100)",
      type        = "continuous",
      notes       = "Screened on Emax but not retained (Huh 2024 Sections 2.3 and 3.2). Its structural content is already carried by DIS_ALOPECIA_AT_AU and STUDY_B7981032 acting on BASE. Table 1: mean 84.9 (SD 21.4) overall; 74.0 (SD 22.5) in non-AT/AU and exactly 100 (SD 0) in AT/AU."
    ),
    REGION = list(
      description = "Geographic region of the enrolling site.",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on BASE, Emax and the drug-chain half-life but not retained (Huh 2024 Sections 2.3 and 3.2). The paper does not tabulate the region levels or their sizes."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1268L,
    n_studies      = 3L,
    n_observations = 11857L,
    age_range      = "12-72 years (Table 1 median 32, min 12, max 72); 170 (13.4%) adolescents 12 to <18 years and 1098 (86.6%) adults >= 18 years",
    weight_range   = "29.6-200.0 kg (Table 1 median 68.4, mean 70.7, SD 17.6)",
    sex_female_pct = 63.5,
    race_ethnicity = "White 889 (70.1%), Asian 286 (22.6%), Other 47 (3.7%), African American 46 (3.6%) (Table 1)",
    disease_state  = "Alopecia areata. 736 (58.0%) non-AT/AU and 532 (42.0%) alopecia totalis or alopecia universalis (Table 1). Baseline SALT mean 84.9 (SD 21.4) overall; 74.0 (SD 22.5) non-AT/AU and 100 (SD 0) AT/AU. Duration since first diagnosis median 6.67 years; duration of the current episode median 2.25 years. Prior pharmacological treatment for alopecia areata: 499 (39.4%) yes, 22 (1.7%) no, 747 (58.9%) unknown.",
    dose_range     = "Ritlecitinib 10, 30 or 50 mg QD, with or without a 200 mg QD 4-week loading dose, plus matching placebo, across the three pooled studies. B7931005 (phase 2a, NCT02974868, n = 95): 200 mg QD x 4 weeks then 50 mg QD x 20 weeks, or placebo, with a single-blind extension after a 4-week drug holiday. B7981015 (phase 2b/3, NCT03732807, n = 715): 200/50, 200/30, 50, 30 and 10 mg QD arms plus two placebo-to-active switch arms. B7981032 (phase 3 long-term, NCT04006457, n = 458): 50 mg QD for rollover participants, 200 mg QD x 4 weeks then 50 mg QD for de novo participants.",
    regions        = "Multiregional (region was screened as a covariate but the source does not tabulate the levels).",
    notes          = "Inclusion required >= 50% scalp hair loss in B7931005 and B7981015 and >= 25% in B7981032 -- the difference the STUDY_B7981032 effect on BASE absorbs. B7981015 and B7981032 enrolled participants aged >= 12 years; B7931005 enrolled adults only. B7981032 was ongoing at the time of the analysis, so only the available data cut was included. Boundary SALT records were common: 5.39% of observations at SALT = 0 and 24.9% at SALT = 100, and they are handled by the paper's censored likelihood rather than by the structural model encoded here."
  )

  ini({
    # ================================================================
    # Structural model -- Huh 2024 Equations 1-4 and Section 3.2.
    #
    # Bounded-outcome transformation (Eq 1), with y the SALT score on
    # 0-100 and alpha the estimated Aranda-Ordaz transformation factor:
    #
    #   z = (100 - y) / 100
    #   x = h(y, alpha) = log((z^(-alpha) - 1) / alpha)
    #
    # so x -> +Inf as SALT -> 100 (complete hair loss) and x -> -Inf as
    # SALT -> 0 (complete regrowth). The back-transform used below is
    #
    #   y = 100 * (1 - (alpha * exp(x) + 1)^(-1 / alpha)).
    #
    # ORIENTATION NOTE. As typeset, Eq 1 pairs z = (100 - y)/100 with a
    # (1 - z) inside the transformation, i.e. it composes to y/100. That
    # composition is falsified by the paper's own numbers, twice over:
    # (i) BASE for the AT/AU group is +11.6 and Section 3.2 states the
    # AT/AU baseline is SALT = 100, but y/100 -> 1 sends the printed
    # transform to -Inf, not +Inf; and (ii) Section 3.2 reports that the
    # Emax of 15.8 "is translated into a complete recovery in SALT
    # score", which requires the drug term (which SUBTRACTS from mu per
    # Eq 3) to move the score toward 0. Both are satisfied only by the
    # (100 - y)/100 composition encoded above, so one of the two halves
    # of Eq 1 as typeset carries a sign/definition slip. See the
    # vignette Errata for the full adjudication.
    #
    # Conditional mean on the transformed scale (Eq 2 / Eq 3):
    #
    #   x | eta = mu(eta) + sigma * eps
    #   mu(eta) = f_b(eta) - f_placebo(t) - f_drug(t)
    #
    # Latent-variable indirect response with transit chains (Eq 4 plus
    # the three placebo / two drug transit compartments added in
    # Section 3.2). kin1 = kout1 and kin2 = kout2 (Section 3.2: separate
    # estimation "was not supported"), so every latent state has an
    # initial condition of kin/kout = 1 and
    #
    #   f_placebo(t) = placebo4 - 1,   f_drug(t) = effect3 - 1.
    #
    # All point estimates below are Huh 2024 Table 2 (final model,
    # NONMEM 7.5.0 SAEM + IMP with Laplace/interaction, ADVAN13); the
    # SIR medians in the same table agree with the asymptotic estimates
    # to the printed precision.
    # ================================================================

    lbase_nonatau <- log(1.92)
    label("Log of the typical baseline transformed SALT score in the non-AT/AU group in studies B7931005 and B7981015 (transformed scale; bare-name base_nonatau = 1.92)")
    # Table 2 "BASE for non-AT/AU" = 1.92 (RSE 4.04%; 95% CI 1.77, 2.07; SIR median 1.92, SIR 95% CI 1.77, 2.07).

    lbase_atau <- log(11.6)
    label("Log of the typical baseline transformed SALT score in the AT/AU group (transformed scale; bare-name base_atau = 11.6, corresponding to SALT = 100 by definition of the AT/AU classification)")
    # Table 2 "BASE for AT/AU" = 11.6 (RSE 2.83%; 95% CI 10.9, 12.2; SIR median 11.6, SIR 95% CI 10.9, 12.1).

    e_study_b7981032_base <- -0.645
    label("Multiplicative fractional shift on the non-AT/AU baseline for records from study B7981032, relative to B7931005 / B7981015; enters as base_nonatau * (1 + e_study_b7981032_base * STUDY_B7981032)")
    # Table 2 "B7981032 effect on BASE for non-AT/AU" = -0.645 (RSE 5.76%; 95% CI -0.718, -0.572; SIR median -0.644, SIR 95% CI -0.717, -0.574). Multiplicative, not additive: Section 3.2 quotes 0.68 for the non-AT/AU group in B7981032, and 1.92 * (1 - 0.645) = 0.68 while 1.92 - 0.645 = 1.28.

    pmax_nonatau <- fixed(0)
    label("Maximum placebo effect on the transformed SALT score in the non-AT/AU group (transformed scale)")
    # Table 2 "Pmax for non-AT/AU" = 0 fix. Section 3.2: "The Pmax parameter was not precisely estimated for the non-AT/AU group and therefore was fixed to zero". Kept as an explicit fixed(0) theta so the estimated-then-fixed provenance is not lost.

    pmax_atau <- 2.75
    label("Maximum placebo effect on the transformed SALT score in the AT/AU group (transformed scale)")
    # Table 2 "Pmax for AT/AU" = 2.75 (RSE 8.52%; 95% CI 2.29, 3.21; SIR median 2.77, SIR 95% CI 2.28, 3.25). Linear scale, not log: the additive IIV in Section 2.3 is explicitly intended "to allow both disease worsening and improving", so an individual Pmax may be negative.

    lthalfrec_pbo <- fixed(log(1.93))
    label("Log of the placebo-chain turnover half-life (weeks; bare-name thalfrec_pbo = 1.93, giving kout1 = kin1 = log(2)/1.93 per week)")
    # Table 2 "kout1 half-life (wk)" = 1.93 fix; Table 2 footnote b: "kout1 half-life was fixed from the final placebo model". With the three transit compartments the placebo mean transit time is (3 + 1)/kout1 = 4 * 1.93/log(2) = 11.1 weeks.

    lemax <- log(15.8)
    label("Log of the maximum ritlecitinib effect on the transformed SALT score (transformed scale; bare-name emax = 15.8, shared by the non-AT/AU and AT/AU groups)")
    # Table 2 "Emax" = 15.8 (RSE 6.23%; 95% CI 13.9, 17.8; SIR median 15.8, SIR 95% CI 13.8, 17.6). Section 3.2: the same Emax applies to both disease-severity groups.

    lec50 <- log(53.6)
    label("Log of the average ritlecitinib concentration producing half of Emax (ng/mL; bare-name ec50 = 53.6)")
    # Table 2 "EC50 (ng/mL)" = 53.6 (RSE 9.83%; 95% CI 43.3, 64.0; SIR median 53.8, SIR 95% CI 43.6, 65.4). Section 3.2 notes this is close to the Cavg of the 50 mg QD regimen (52 ng/mL).

    lthalfrec_drug <- log(7.80)
    label("Log of the drug-chain turnover half-life in the non-AT/AU group (weeks; bare-name thalfrec_drug = 7.80, giving kout2 = kin2 = log(2)/7.80 per week)")
    # Table 2 "kout2 half-life (wk)" = 7.80 (RSE 7.61%; 95% CI 6.64, 8.96; SIR median 7.71, SIR 95% CI 6.51, 9.01). With the two transit compartments the non-AT/AU drug mean transit time is (2 + 1)/kout2 = 3 * 7.80/log(2) = 33.8 weeks, the figure quoted in Section 3.2.

    e_atau_thalfrec_drug <- -0.601
    label("Multiplicative fractional shift on the drug-chain turnover half-life for AT/AU participants; enters as thalfrec_drug * (1 + e_atau_thalfrec_drug * DIS_ALOPECIA_AT_AU)")
    # Table 2 "AT/AU effect on kout2 half-life (wk)" = -0.601 (RSE 5.10%; 95% CI -0.661, -0.541; SIR median -0.598, SIR 95% CI -0.660, -0.542). Multiplicative: 7.80 * (1 - 0.601) = 3.11 weeks, the AT/AU half-life quoted in Section 3.2, and 3 * 3.11/log(2) = 13.5 weeks, the AT/AU mean transit time quoted there.

    alpha <- 1.19
    label("Aranda-Ordaz transformation factor applied to the scaled non-boundary SALT data (dimensionless)")
    # Table 2 "Transformation factor" = 1.19 (RSE 2.73%; 95% CI 1.12, 1.25; SIR median 1.19, SIR 95% CI 1.13, 1.25).

    # ================================================================
    # Inter-individual variability. Section 2.3: BASE, Emax and
    # kin2/kout2 carry multiplicative exponential IIV (Pi = Ppop *
    # exp(eta_i)); Pmax carries additive IIV (Pi = Ppop + eta_i).
    # Variances and covariances are Table 2's omega^2 and
    # covariance rows, which are already on the variance scale.
    #
    # Table 2 reports covariances only for the (Emax, Pmax) and
    # (Emax, kout2) pairs, so the (Pmax, kout2) element is zero and BASE
    # is diagonal. The resulting 3x3 block is positive definite (leading
    # minors 7.62, 5.550 and 3.104).
    #
    # The kout2 eta is attached to the HALF-LIFE, matching the parameter
    # the paper actually estimated (Section 2.3: "t1/2 was estimated
    # instead of kin or kout") and matching the level the AT/AU
    # covariate acts on ("AT/AU effect on kout2 half-life"). Placing it
    # on the rate constant instead would flip the sign of its covariance
    # with Emax; see the vignette Errata.
    # ================================================================

    etalbase ~ 0.225
    # Table 2 "omega^2 BASE" = 0.225 (RSE 6.54%; 95% CI 0.196, 0.254; SIR median 0.224). Shrinkage 16% (Table 2 footnote b).

    etapmax + etalemax + etalthalfrec_drug ~ c(
      7.62,
      -0.442, 0.754,
      0.000,  0.549, 0.973
    )
    # Diagonals: Table 2 "omega^2 Pmax" = 7.62 (RSE 8.27%; 95% CI 6.39, 8.86), "omega^2 Emax" = 0.754 (RSE 12.2%; 95% CI 0.573, 0.934), "omega^2 Kout2" = 0.973 (RSE 9.14%; 95% CI 0.799, 1.15).
    # Off-diagonals: Table 2 "Covariance-Emax and Pmax" = -0.442 (RSE 19.5%; 95% CI -0.611, -0.273) and "Covariance-Emax and kout2" = 0.549 (RSE 12.2%; 95% CI 0.418, 0.681). No Pmax-kout2 covariance is reported, hence 0.
    # Shrinkages (Table 2 footnote b): 33% for Pmax, 37% for Emax, 37% for kout2.

    # ================================================================
    # Residual error -- additive on the transformed scale (Eq 2:
    # x|eta = mu(eta) + sigma * eps with eps ~ N(0, 1)).
    # ================================================================

    addSd <- 1.18
    label("Additive residual standard deviation on the transformed SALT scale (transformed-scale units; the sigma of Huh 2024 Eq 2)")
    # Table 2 "Residual error" = 1.18 (RSE 1.51%; 95% CI 1.14, 1.21; SIR median 1.18, SIR 95% CI 1.14, 1.21). Epsilon shrinkage 25% (Table 2 footnote b).
  })

  model({
    # --- Individual parameters ----------------------------------------
    # BASE: separate typical values by disease severity, with the
    # B7981032 inclusion-criterion shift applied to the non-AT/AU group
    # only (the AT/AU baseline is pinned at SALT = 100 by the disease
    # classification, so study entry criteria cannot move it).
    base_i <-
      ((1 - DIS_ALOPECIA_AT_AU) * exp(lbase_nonatau) *
         (1 + e_study_b7981032_base * STUDY_B7981032) +
         DIS_ALOPECIA_AT_AU * exp(lbase_atau)) *
      exp(etalbase)

    # Pmax: linear scale with additive IIV so an individual placebo
    # effect may be negative (disease worsening).
    pmax_i <-
      (1 - DIS_ALOPECIA_AT_AU) * pmax_nonatau +
      DIS_ALOPECIA_AT_AU * pmax_atau +
      etapmax

    emax_i <- exp(lemax + etalemax)
    ec50 <- exp(lec50)

    # Turnover half-lives, and the rate constants they imply. kin equals
    # kout in both chains, so only kout appears below.
    thalfrec_pbo <- exp(lthalfrec_pbo)
    thalfrec_drug_i <- exp(lthalfrec_drug + etalthalfrec_drug) *
      (1 + e_atau_thalfrec_drug * DIS_ALOPECIA_AT_AU)

    kout1 <- log(2) / thalfrec_pbo
    kout2 <- log(2) / thalfrec_drug_i

    # --- Placebo latent chain (Eq 4 + three transit compartments) -----
    # The on-treatment indicator stimulates production into the head of
    # the chain; the delay is carried by the three transits.
    d/dt(placebo1) <- kout1 * (1 + TRT_PHASE * pmax_i) - kout1 * placebo1
    d/dt(placebo2) <- kout1 * placebo1 - kout1 * placebo2
    d/dt(placebo3) <- kout1 * placebo2 - kout1 * placebo3
    d/dt(placebo4) <- kout1 * placebo3 - kout1 * placebo4

    # --- Drug latent chain (Eq 4 + two transit compartments) ----------
    # CAV is the average ritlecitinib concentration over the current
    # inter-SALT-record interval; CAV = 0 collapses the Emax term, which
    # is how placebo periods and treatment interruptions are encoded.
    drugStim <- emax_i * CAV / (ec50 + CAV)
    d/dt(effect1) <- kout2 * (1 + drugStim) - kout2 * effect1
    d/dt(effect2) <- kout2 * effect1 - kout2 * effect2
    d/dt(effect3) <- kout2 * effect2 - kout2 * effect3

    # All latent states start at their untreated steady state kin/kout,
    # which equals 1 because kin1 = kout1 and kin2 = kout2 (Eq 4
    # initial conditions PBO(0) = kin1/kout1 and E(0) = kin2/kout2).
    placebo1(0) <- 1
    placebo2(0) <- 1
    placebo3(0) <- 1
    placebo4(0) <- 1
    effect1(0) <- 1
    effect2(0) <- 1
    effect3(0) <- 1

    # --- Conditional mean and observation -----------------------------
    # Eq 3: mu = f_b - f_placebo - f_drug, with f_placebo = PBO(t) -
    # kin1/kout1 and f_drug = E(t) - kin2/kout2 (both baselines = 1).
    salt_transformed <- base_i - (placebo4 - 1) - (effect3 - 1)

    # Back-transform onto the published 0-100 SALT scale (inverse of
    # Eq 1). The 0 and 100 boundaries are approached asymptotically;
    # Huh 2024 handles the observed boundary mass with a censored
    # likelihood, which is a fitting-time construct and is not part of
    # the structural model encoded here.
    salt <- 100 * (1 - (alpha * exp(salt_transformed) + 1)^(-1 / alpha))

    # Eq 2: additive residual error on the transformed scale. To fit
    # observed SALT with this model, transform the observations first:
    #   salt_transformed_obs <- log((((100 - SALT) / 100)^(-alpha) - 1) / alpha)
    salt_transformed ~ add(addSd)
  })
}
