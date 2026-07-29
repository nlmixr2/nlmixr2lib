Conrado_2014_alzheimer <- function() {
  description <- "Updated Alzheimer's disease progression model for ADAS-Cog total score (0-70) over time, fit to ~25,000 ADAS-Cog observations from 4,494 subjects across 15 studies in the Coalition Against Major Diseases (CAMD) Alzheimer's database. The expected score is described by a Richards three-parameter logistic growth function whose asymptote is the upper boundary score 70; the residual distribution on the bounded ADAS-Cog/70 scale is a beta distribution with precision parameter TAU. Subject-level baseline and slope inter-individual variability are correlated through a 2x2 BLOCK; covariates of sex on baseline, APOE4-allele count on baseline and slope, age on slope, and concomitant Alzheimer's-symptomatic medication on slope are reproduced from the source. The original publication adds a third-level (study) random effect on baseline and slope plus study-1131-specific scalers; nlmixr2 does not natively support multi-level random effects, so those layers are dropped here and documented in the validation vignette's Errata section."
  reference <- paste(
    "Conrado DJ, Denney WS, Chen S, Ito K. (2014).",
    "An updated Alzheimer's disease progression model:",
    "incorporating non-linearity, beta regression, and a third-level random effect in NONMEM.",
    "J Pharmacokinet Pharmacodyn 41(6):581-598.",
    "doi:10.1007/s10928-014-9375-z. PMID: 25168488.",
    "DDMORE Foundation Model Repository: DDMODEL00000290.",
    sep = " "
  )
  vignette <- "Conrado_2014_alzheimer"
  units <- list(
    time = "day",
    dosing = "(none; disease-progression model, no drug input)",
    concentration = "(ADAS-Cog total score, 0-70, unitless)"
  )
  ddmore_id    <- "DDMODEL00000290"
  replicate_of <- NULL

  covariateData <- list(
    AGE = list(
      description        = "Subject age at baseline (years)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in Conrado 2014. Centred on 75 years inside model() (the source-paper centring); the slope-coefficient e_sl_age is the fractional change in disease-progression slope per year above 75. Source column AGE in the bundle's NONMEM input dataset.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator with 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) -- the source paper centres the baseline-coefficient on the female cohort, which it labels the 'most common' category",
      notes              = "Source column SEX in the bundle is coded 1 = female, 2 = male. The canonical SEXF (1 = female) is derived as `SEXF = as.integer(SEX == 1)`. The source-paper coefficient 0.953 multiplies the typical baseline score for males (SEXF = 0); the female reference cohort takes a multiplier of 1.",
      source_name        = "SEX"
    ),
    APOE4_COUNT = list(
      description        = "APOE-epsilon-4 allele count: 0 = non-carrier, 1 = heterozygous, 2 = homozygous (continuous representation, with the population mean used as the centring value)",
      units              = "(count, 0 / 1 / 2 alleles per subject)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred on 0.72 inside model() -- the population-mean APOE-epsilon-4 count in the CAMD cohort, taken directly from the source `.mod`'s `APOE4C - 0.72` re-centring expressions. Source column APOE4C in the bundle is the cleaned continuous version of the upstream APOE4 column (0 = non-carrier, 1 = heterozygous, 2 = homozygous, 3 = unknown), with the unknown subjects recoded to the population mean.",
      source_name        = "APOE4C"
    ),
    CONMED_AD = list(
      description        = "Concomitant Alzheimer's-symptomatic medication indicator at baseline (typically a cholinesterase inhibitor and / or memantine), 1 = on treatment, 0 = not on treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (on concomitant Alzheimer's-symptomatic medication) -- the source paper labels this as the 'most common' category and centres the slope-effect coefficient on the on-treatment cohort, so the multiplicative slope factor is 1 for CONMED_AD = 1 and 1 + e_sl_conmed_ad_off for CONMED_AD = 0. The non-standard most-common-as-reference convention is preserved from the source for traceability.",
      notes              = "Source column COMED2 in the bundle. Time-fixed per subject in the CAMD-derived dataset. The source `.mod` does not name the symptomatic-medication class beyond the binary flag; the publication context (CAMD ADAS-Cog disease-progression dataset, 2014) makes cholinesterase-inhibitor / memantine the dominant interpretation. See the canonical `CONMED_AD` entry in inst/references/covariate-columns.md for the full notes on the reference-category inversion.",
      source_name        = "COMED2"
    )
  )

  population <- list(
    n_subjects     = 4494L,
    n_studies      = 15L,
    age_range      = "Adults with Alzheimer's disease, mild cognitive impairment, or healthy elderly comparators; specific age-range and median not extracted from the source bundle and the linked publication PDF was not on disk in /home/bill/github/mab_human_consensus/literature for cross-check at extraction time.",
    weight_range   = "(not extracted; not a covariate in the disease-progression model)",
    sex_female_pct = NA_real_,
    disease_state  = "Alzheimer's disease (mild-to-moderate AD predominates), pooled across the 15 randomised-controlled-trial arms contributing to the CAMD ADAS-Cog disease-progression dataset (2014 release).",
    dose_range     = "(not applicable; disease-progression model on the placebo and active-treatment arms pooled, with concomitant Alzheimer's-symptomatic medication as a covariate rather than a dosed input)",
    regions        = "(not extracted; CAMD pools randomised-controlled-trial arms across multiple international sponsors)",
    notes          = "Subject and study counts taken directly from the bundle's Output_real_CPathAD.lst header (TOT. NO. OF INDIVIDUALS = 4494; ETABAR N = 15 for the study-level etas, indicating 15 randomised-controlled-trial arms in the source dataset). Demographic detail (age range / median, weight range, sex split, regional breakdown) is described in the Conrado 2014 publication's Methods / Results tables but the publication PDF was not on disk at extraction time, so finer-grained population descriptors are recorded as NA. The dataset is the Coalition Against Major Diseases (CAMD) ADAS-Cog database aggregated by the Critical Path Institute for the published model-fitting exercise."
  )

  ini({
    # All values come from Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE
    # block (after MINIMIZATION SUCCESSFUL at line 872). The .mod $THETA /
    # $OMEGA blocks are the initial estimates only -- the .lst final estimates
    # are the canonical source per the DDMORE-source convention. The
    # publication's Tables 4-5 were not available on disk for cross-check;
    # this is documented in the vignette's Errata section. The TVBL parameter
    # is log-transformed because the typical baseline score must remain
    # strictly positive; TVSL is left on its natural scale because the
    # disease-progression slope can in principle be negative (the source-
    # paper IIV is additive on slope, consistent with this).

    # Structural parameters -- Richards three-parameter logistic for the
    # population mean ADAS-Cog/70.
    lrbase <- log(22.2) ; label("Log of typical baseline ADAS-Cog total score (0-70 scale)")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 1 = 2.22E+01.

    sl <- 0.153 ; label("Typical disease-progression slope (1/year on the bounded 0-1 ADAS-Cog/70 scale, before back-transformation through the Richards growth function)")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 2 = 1.53E-01.

    lshape <- log(6.91) ; label("Log of Richards three-parameter logistic shape factor (unitless; controls the curvature of the bounded growth trajectory)")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 3 = 6.91E+00.

    ltau <- log(87.5) ; label("Log of beta-distribution precision parameter TAU = alpha + beta (unitless; larger values mean a tighter residual distribution on the bounded 0-1 scale)")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 4 = 8.75E+01.

    # Covariate effects.
    e_bl_male <- 0.953 ; label("Multiplicative factor on typical baseline ADAS-Cog for males (SEXF = 0); the female reference cohort takes a factor of 1")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 5 = 9.53E-01.

    e_sl_age <- -0.024 ; label("Fractional change in typical disease-progression slope per year of age above 75 (multiplicative on slope: factor = 1 + e_sl_age * (AGE - 75))")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 6 = -2.40E-02.

    e_sl_apoe4 <- 0.195 ; label("Fractional change in typical disease-progression slope per APOE-epsilon-4 allele count unit above 0.72 (multiplicative on slope: factor = 1 + e_sl_apoe4 * (APOE4_COUNT - 0.72))")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 7 = 1.95E-01.

    e_sl_conmed_ad_off <- -0.302 ; label("Additive shift to the unit slope-multiplier when CONMED_AD = 0 (off concomitant AD-symptomatic medication); slope_factor = 1 + (1 - CONMED_AD) * e_sl_conmed_ad_off, so off-treatment subjects take a slope multiplier of 1 + (-0.302) = 0.698")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 8 = -3.02E-01.

    e_bl_apoe4 <- 0.0372 ; label("Fractional change in typical baseline ADAS-Cog per APOE-epsilon-4 allele count unit above 0.72 (multiplicative on baseline: factor = 1 + e_bl_apoe4 * (APOE4_COUNT - 0.72))")
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE TH 9 = 3.72E-02.

    # Subject-level inter-individual variability -- BLOCK(2) on (etalrbase, etasl).
    # etalrbase is on the log-baseline scale (multiplicative log-normal IIV on BL);
    # etasl is on the natural slope scale (additive IIV on SL). The off-diagonal
    # is the cross-scale covariance reported by NONMEM directly.
    etalrbase + etasl ~ c(0.156,
                       0.0224, 0.0413)
    # Output_real_CPathAD.lst FINAL PARAMETER ESTIMATE OMEGA BLOCK(2) #1:
    # ETA1,ETA1 = 1.56E-01, ETA1,ETA2 = 2.24E-02, ETA2,ETA2 = 4.13E-02.

    # No residual-error sigma parameter: the source uses a -2*log-likelihood
    # of the beta distribution as Y, with the precision parameter TAU acting
    # as the dispersion control. nlmixr2 / rxode2 expose dbeta as a built-in
    # observation likelihood (`obs ~ dbeta(alpha, beta)`), which reproduces
    # the source likelihood exactly with no need for a separate residual-error
    # parameter.
  })

  model({
    # The three-level (study x subject) random-effects structure of the source
    # paper is collapsed to the subject level here; nlmixr2 does not natively
    # support multi-level random effects, so the study-level (level-4) etas
    # ETA(3), ETA(4) and the study-1131-specific scalers THETA(10), THETA(11)
    # are dropped. The fixed-effect THETAs (TH1-TH9) and the subject-level
    # BLOCK(2) BSV (TH/ETA1, ETA2) are reproduced exactly from the .lst final
    # estimates. See the validation vignette's Errata section for the deviation
    # rationale.

    # 1. Time-in-years scaling. The source dataset uses TIME in days (max ~5 *
    # 365.25 = 1826 days across the disease-progression follow-up window); the
    # source `.mod` rescales to YTIME = TIME/365.25 to avoid sub-percent slope
    # values. Reproduce the same rescaling here.
    yt <- time / 365.25

    # 2. Covariate factors on baseline and slope. Each factor is exactly the
    # source's BLAPOE4C / SLAPOE4C / SLAGE / SLCOMED2 / BLSEX construction,
    # rewritten with the canonical covariate-column names.
    bl_apoe4_factor   <- 1 + e_bl_apoe4 * (APOE4_COUNT - 0.72)
    sl_apoe4_factor   <- 1 + e_sl_apoe4 * (APOE4_COUNT - 0.72)
    sl_age_factor     <- 1 + e_sl_age   * (AGE        - 75.00)
    sl_conmed_factor  <- 1 + (1 - CONMED_AD) * e_sl_conmed_ad_off
    bl_sex_factor     <- SEXF * 1 + (1 - SEXF) * e_bl_male

    bl_cov <- bl_sex_factor * bl_apoe4_factor
    sl_cov <- sl_age_factor * sl_apoe4_factor * sl_conmed_factor

    # 3. Individual baseline and slope. The source applies a multiplicative
    # log-normal IIV on baseline (etalrbase on the log scale) and an additive IIV
    # on slope (etasl on the natural slope scale). The covariate factors enter
    # multiplicatively in both cases; on the log-baseline scale, the
    # multiplicative bl_cov is added as log(bl_cov).
    rbase <- exp(lrbase + log(bl_cov) + etalrbase)
    sl_indiv <- (sl + etasl) * sl_cov

    # 4. Richards three-parameter logistic mean on the bounded 0-1 ADAS-Cog/70
    # scale. The source equations live in the .mod $PRED block:
    #   DEN1 = BL**SHAPE
    #   DEN2 = (70**SHAPE) - (BL**SHAPE)
    #   DEN3 = EXP(-SHAPE * SL * YTIME)
    #   MUR  = BL / (DEN1 + DEN2*DEN3)**(1/SHAPE)
    # Note: the source defines BL on the 0-70 ADAS-Cog scale (typical 22.2),
    # so MUR is BL / (...)**(1/SHAPE) which yields a value in (0, 1) only when
    # the right-hand side denominator is in (BL, 70). Algebraically MUR is
    # equivalent to the bounded Richards growth from initial value BL/70 toward
    # asymptote 1; we keep the literal source form for traceability.
    shape <- exp(lshape)
    den1  <- rbase^shape
    den2  <- 70^shape - rbase^shape
    den3  <- exp(-shape * sl_indiv * yt)
    denn  <- den1 + den2 * den3
    mur   <- rbase / denn^(1/shape)

    # 5. Beta-distribution shape parameters from the source $PRED block.
    tau   <- exp(ltau)
    alpha <- mur * tau
    beta  <- (1 - mur) * tau

    # 6. Observation. The source observation Y is a -2*log-likelihood of the
    # beta distribution evaluated at WDV = bounded ADAS-Cog/70 (the
    # `ADASTRANS = (ADAS / 70)` column in the input dataset). The natural
    # nlmixr2 / rxode2 expression of this likelihood is `dbeta(alpha, beta)`
    # on the bounded outcome ADAS_NORM = ADAS / 70. Predicted values are on
    # the bounded 0-1 scale; multiply by 70 outside the model to get back to
    # the 0-70 ADAS-Cog total score for plotting and validation.
    ADAS_NORM <- mur
    ADAS_NORM ~ dbeta(alpha, beta)
  })
}
