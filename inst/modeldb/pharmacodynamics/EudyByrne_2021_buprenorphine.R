EudyByrne_2021_buprenorphine <- function() {
  description <- "Indirect-response pharmacodynamic model of sublingual buprenorphine effect on MOTHER NAS severity scores in neonates with neonatal opioid withdrawal syndrome (NOWS). Natural postnatal-age-driven withdrawal decay (NOWSMAX, NOWSM) sets a time-varying baseline NAS score via Kin * (1 + NOWST); buprenorphine central-compartment concentration Cbuprenorphine stimulates NAS-score elimination (Kout) through a Hill Emax term (EMAX * C / (EC50 + C) + 1). The model has NO PK ODE; the buprenorphine concentration Cbuprenorphine (ng/mL) is a required time-varying input covariate carried in the event table. Upstream buprenorphine PK from Moore et al. 2018 (Clin Pharmacol Ther 103:1029-1037; DOI 10.1002/cpt.1064) is not packaged in nlmixr2lib at extraction time; users needing to simulate the PD from a dose regimen must supply concentrations from an external PK source. The initial NAS score is set to its drug-free quasi-steady state, NOWS0 = Kin * (1 + NOWST) / Kout, at the baseline PNA."
  reference <- paste(
    "Eudy-Byrne R, Zane N, Adeniyi-Jones SC, Kaushal G, Ruiz-Garcia A,",
    "Gastonguay MR, Kraft WK.",
    "Pharmacometric dose optimization of buprenorphine in neonatal opioid",
    "withdrawal syndrome. Clin Transl Sci. 2021;14(6):2171-2183.",
    "doi:10.1111/cts.13074.",
    "PD model builds on the upstream buprenorphine PK model from",
    "Moore JN, Gastonguay MR, Ng CM, et al. Clin Pharmacol Ther",
    "2018;103(6):1029-1037 (doi:10.1002/cpt.1064), which is not packaged",
    "in nlmixr2lib at extraction time; the paper's NAS-score simulations",
    "reused BBORN observed buprenorphine exposures.",
    sep = " "
  )
  vignette <- "EudyByrne_2021_buprenorphine"
  units <- list(
    time          = "h",
    dosing        = "N/A (PD-only; buprenorphine concentration is a required input covariate)",
    concentration = "MOTHER NAS score (observation, unitless integer scale); ng/mL (buprenorphine input covariate)"
  )

  depends <- c("Cbuprenorphine", "PNA")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    nows = list(analyte = "MOTHER NAS severity scores", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    Cbuprenorphine = list(
      description        = "Buprenorphine plasma concentration (central compartment) driving the drug effect in the indirect-response NAS-score PD model. Time-varying; carried in the event-table `covariates` block.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Paper equation (Results / Model development): EFFECTdrug = EMAX * C2 / (EC50 + C2) + 1, with `C2` the buprenorphine concentration and EC50 = 0.942 ng/mL. Observed BBORN concentrations ranged from < 0.1 ng/mL (LLQ) to 0.6 ng/mL, mean 0.249 (SD 0.101) ng/mL; BPHORE observed mean 0.275 (SD 0.243) ng/mL, LLQ 0.05 ng/mL (paper Methods). This is an EXTERNAL input to the PD model -- the upstream Moore 2018 PK model (2-compartment, PNA covariate on CL) is not packaged in nlmixr2lib. Not in inst/references/covariate-columns.md because it is a paper-specific PD driver rather than a demographic column.",
      source_name        = "C2 (Eudy-Byrne 2021 Results / Model development, equation for EFFECTdrug)"
    ),
    PNA = list(
      description        = "Postnatal age (chronological age since birth). Time-varying; increments with the simulation clock and drives the natural NAS decay term NOWST = NOWSMAX * exp(-NOWSM * PNA_days).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Canonical PNA in inst/references/covariate-columns.md is in MONTHS. The Eudy-Byrne 2021 paper reports NOWSM in units of 1/DAY (Table S3), so inside model() PNA (months) is converted to days via `pna_days <- PNA * 30.4375` before use in the exponential decay term. Observed BBORN age at last dose 21.1 (SD 11.6) days = 0.693 (SD 0.381) months; BPHORE 25.5 (SD 9.2) days = 0.838 (SD 0.302) months (Table S2).",
      source_name        = "PNA (Eudy-Byrne 2021 Results / Model development, NOWST equation)"
    )
  )

  population <- list(
    species        = "human (neonates)",
    n_subjects     = 28L,
    n_studies      = 1L,
    age_range      = "postnatal age at last dose: BBORN 21.1 (SD 11.6) days; BPHORE 25.5 (SD 9.2) days",
    weight_range   = "birth weight: BBORN 3.10 (SD 0.430) kg; BPHORE 3.10 (SD 0.292) kg",
    sex_female_pct = 39,
    race_ethnicity = NULL,
    disease_state  = "neonatal opioid withdrawal syndrome (NOWS) / neonatal abstinence syndrome (NAS) in infants greater than 36 weeks gestation exposed to opioids in utero (predominantly maternal methadone).",
    dose_range     = "sublingual buprenorphine; BBORN starting 5.3 microgram/kg q8h (max 20 microgram/kg q8h = 60 microgram/kg/day) with 25% up-titration and 10% down-titration; BPHORE starting 8 microgram/kg q8h (max 25 microgram/kg q8h = 75 microgram/kg/day) with 33% up-titration and 15% down-titration (Table S1). PD model uses buprenorphine concentrations directly, not doses.",
    regions        = "USA (Thomas Jefferson University, Philadelphia; BBORN NCT01452789 and BPHORE NCT03608696)",
    n_pk_obs       = 117L,
    n_pd_obs       = 3609L,
    notes          = "PD model was estimated on N = 28 BBORN infants who received buprenorphine (paper 'Clinical data used in model-based analysis'). 117 buprenorphine concentrations and 3609 MOTHER NAS score observations were used; observations coincident with phenobarbital or clonidine adjunct therapy, or below LLQ, were excluded (paper Methods). BPHORE (N = 10) was used to validate simulation-based dosing recommendations, not to re-estimate the model. Reference categories for covariates are not applicable (no covariates on structural parameters were reported)."
  )

  ini({
    # =====================================================================
    # Structural PD parameters (Eudy-Byrne 2021 Table S3, Supplementary
    # Information cts13074-sup-0007-tables3.docx). All fixed effects were
    # estimated with good precision; 95% CI shown in comments does not
    # encompass zero for any structural parameter.
    # =====================================================================

    # Maximum baseline (birth) NAS score in the natural-history term:
    # NOWST = NOWSMAX * exp(-NOWSM * PNA_days). Paper reports the estimate
    # on the natural (untransformed) scale; encoded here as log for a
    # positively-constrained parameter with the standard log-normal IIV
    # (Table S3 %CV column back-transforms as CV = sqrt(exp(omega^2)-1)).
    lnowsmax <- log(1.92)
    label("NOWSMAX: maximum baseline NAS-natural-history term at PNA = 0 (score, unitless)")  # Eudy-Byrne 2021 Table S3: NOWSMAX = 1.92 (95% CI 1.76-2.08)

    # First-order natural decay rate of the withdrawal-history term with
    # postnatal age. Units 1/day.
    lnowsm <- log(0.107)
    label("NOWSM: natural NAS-withdrawal decay rate with postnatal age (1/day)")  # Eudy-Byrne 2021 Table S3: NOWSM = 0.107 (95% CI 0.102-0.112)

    # Maximum stimulation of NAS-score elimination (Kout) by buprenorphine.
    lemax <- log(1.85)
    label("EMAX: maximum stimulation of NAS-score elimination by buprenorphine (unitless multiplier on Kout above baseline)")  # Eudy-Byrne 2021 Table S3: EMAX = 1.85 (95% CI 1.83-1.87)

    # Concentration eliciting half-maximum drug effect.
    lec50 <- log(0.942)
    label("EC50: buprenorphine concentration for half-maximum drug effect (ng/mL)")  # Eudy-Byrne 2021 Table S3: EC50 = 0.942 (95% CI 0.870-1.01)

    # Zero-order NAS production rate.
    lkin <- log(0.139)
    label("Kin: NAS score production rate (score/h)")  # Eudy-Byrne 2021 Table S3: Kin = 0.139 (95% CI 0.128-0.151)

    # First-order NAS elimination rate constant.
    lkout <- log(0.0301)
    label("Kout: first-order NAS score elimination rate constant (1/h)")  # Eudy-Byrne 2021 Table S3: Kout = 0.0301 (95% CI 0.0300-0.0302)

    # =====================================================================
    # Inter-individual variability (Table S3 omega-block). Table S3 lists
    # the variances/covariances on the log-normal scale. NOWSMAX and NOWSM
    # are correlated (2x2 block), Kout and EMAX are independent (diagonal).
    # Back-checks:
    #   omega^2 = 1.14  -> CV = sqrt(exp(1.14) - 1) * 100 = 146%   (matches Table S3 %CV)
    #   omega^2 = 1.42  -> CV = 177%                              (matches)
    #   omega^2 = 0.108 -> CV = 33.8%                             (matches)
    #   omega^2 = 0.726 -> CV = 103%                              (matches)
    #   cov(NOWSMAX,NOWSM) = 0.990; corr = 0.990/sqrt(1.14*1.42) = 0.778
    #   (matches Table S3 "0.778 (corr)")
    # Table S3 flags NOWSMAX/NOWSM as poorly estimated (95% CI encompasses
    # zero); values are retained as reported per the "carry basic-model IIV
    # forward when the final model is silent" convention. See vignette
    # Assumptions and deviations.
    # =====================================================================
    etalnowsmax + etalnowsm ~ c(1.14,
                                0.990, 1.42)  # Eudy-Byrne 2021 Table S3: omega1,1 = 1.14; omega2,1 = 0.990; omega2,2 = 1.42
    etalkout ~ 0.108   # Eudy-Byrne 2021 Table S3: omega3,3 = 0.108 (%CV 33.8, shrinkage 9.26%)
    etalemax ~ 0.726   # Eudy-Byrne 2021 Table S3: omega4,4 = 0.726 (%CV 103, shrinkage 15.6%)

    # =====================================================================
    # Residual error. Table S3 reports SIGMAadd = 2.30 in units of "score"
    # (linear, not squared), consistent with an additive SD on the NAS-score
    # scale. The paper Methods (Adaptive dose simulations) describes
    # simulations drawing residuals from NID(0, sigma^2); the tight 95% CI
    # (2.29, 2.30) is consistent with a well-estimated SD parameter on the
    # score scale.
    # =====================================================================
    addSd_nows <- 2.30
    label("Additive residual SD on MOTHER NAS score (score units)")  # Eudy-Byrne 2021 Table S3: SIGMAadd = 2.30 (95% CI 2.29-2.30)
  })

  model({
    # -------------------------------------------------------------------
    # Individual parameters (typical value * exp(eta)).
    # -------------------------------------------------------------------
    nowsmax <- exp(lnowsmax + etalnowsmax)
    nowsm   <- exp(lnowsm   + etalnowsm)     # 1/day
    emax    <- exp(lemax    + etalemax)
    ec50    <- exp(lec50)
    kin     <- exp(lkin)                     # score/hr
    kout    <- exp(lkout    + etalkout)      # 1/hr

    # -------------------------------------------------------------------
    # Natural NAS-withdrawal course with postnatal age. The canonical PNA
    # covariate is in MONTHS (inst/references/covariate-columns.md), but
    # NOWSM is estimated in 1/DAY per Table S3. Convert PNA (months) to
    # days using the standard 30.4375 days/month factor. Same pattern used
    # by Zhao 2018 (see PNA canonical notes).
    # -------------------------------------------------------------------
    pna_days <- PNA * 30.4375
    nowst    <- nowsmax * exp(-nowsm * pna_days)

    # -------------------------------------------------------------------
    # Drug effect on Kout: buprenorphine stimulates NAS-score elimination.
    # Paper equation (Results / Model development):
    #   EFFECTdrug = EMAX * C2 / (EC50 + C2) + 1
    # The "+ 1" makes EFFECTdrug = 1 at C2 = 0 so d/dt(nows) reduces to
    # kin*(1+nowst) - kout*nows (drug-free indirect-response baseline).
    # -------------------------------------------------------------------
    effect_drug <- 1 + emax * Cbuprenorphine / (ec50 + Cbuprenorphine)

    # -------------------------------------------------------------------
    # ODE for the MOTHER NAS score. Time unit = hours (matches Kin/Kout).
    # -------------------------------------------------------------------
    d/dt(nows) <- kin * (1 + nowst) - kout * nows * effect_drug

    # -------------------------------------------------------------------
    # Initial condition (paper Results / Model development):
    #   NOWS0 = Kin * (1 + NOWST) / Kout
    # This is the drug-free quasi-steady-state NAS score at the baseline
    # PNA. As PNA -> infinity, NOWST -> 0 and NOWS0 -> Kin/Kout = 4.62
    # score (untreated long-term NAS floor).
    # -------------------------------------------------------------------
    nows(0) <- kin * (1 + nowst) / kout

    # -------------------------------------------------------------------
    # Additive residual error on the NAS-score scale.
    # -------------------------------------------------------------------
    nows ~ add(addSd_nows)
  })
}
