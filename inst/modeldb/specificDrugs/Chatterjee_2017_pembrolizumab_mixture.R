Chatterjee_2017_pembrolizumab_mixture <- function() {
  description <- paste(
    "Four-subpopulation MIXTURE tumor-size (sum of longest diameters, SLD)",
    "model for pembrolizumab in advanced melanoma, developed by Chatterjee",
    "et al. (Merck) on the KEYNOTE-001 melanoma cohorts (n = 364, October",
    "2013 cutoff). Every subject carries the same bi-component structure --",
    "a 'shallow' portion of the tumor that grows at KL and shrinks at KD",
    "plus a static 'deep' portion -- and the four latent classes differ only",
    "in which of those terms are switched on:",
    "escape (KD = 0, no deep portion, shallow baseline x 2.00),",
    "monophasic slow (KD, no deep portion),",
    "biphasic (KD x 4.17, deep portion present) and",
    "monophasic fast (KD x 4.17, no deep portion), giving",
    "SLD(t) = BASEL1 * exp((KL - KD) * t) + BASEL2.",
    "Number of target lesions (power) and ECOG performance status were",
    "retained on the shallow baseline. The paper's second-stage multinomial",
    "logistic regression on the class logits is carried in the same file:",
    "it returns prob_escape, prob_monophasic_slow, prob_biphasic and",
    "prob_monophasic_fast as derived outputs from AUCss-6weeks, the number",
    "of affected lymph nodes and baseline tumor size, so a simulated subject",
    "can be assigned a class before the tumor-size model is solved.",
    "Exposure was NOT a significant predictor of tumor shrinkage and enters",
    "only the class probabilities; the paper's conclusion is that response",
    "is flat over 2-10 mg/kg. There is no PK input -- exposure is supplied",
    "per subject as the covariate AUC_PEMBRO, which the source analysis",
    "obtained as dose/CL from the companion pembrolizumab population-PK",
    "model (Ahamadi 2017; available in this library as",
    "Ahamadi_2017_pembrolizumab)."
  )
  reference <- paste(
    "Chatterjee MS, Elassaiss-Schaap J, Lindauer A, Turner DC, Sostelly A,",
    "Freshwater T, Mayawala K, Ahamadi M, Stone JA, de Greef R, Kondic AG,",
    "de Alwis DP.",
    "Population pharmacokinetic/pharmacodynamic modeling of tumor size",
    "dynamics in pembrolizumab-treated advanced melanoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(1):29-39.",
    "doi:10.1002/psp4.12140. PMID: 27896901. PMCID: PMC5270297.",
    "Structural equation from the main-article Methods ('Initial",
    "exposure-response tumor size (mixture) model'); per-class",
    "parameterization from supplementary Table 1; final tumor-size parameter",
    "values from main-article Table 2; multinomial-regression parameter",
    "values from supplementary Table 2B; covariate parameterization and",
    "centering constants from the supplementary NONMEM control streams",
    "('Final Tumor Size Reduction Model' and 'Multinomial Logistic",
    "Regression Model').",
    sep = " "
  )
  vignette <- "Chatterjee_2017_pembrolizumab_melanoma"

  units <- list(
    time = "day",
    dosing = "n/a (no PK input; pembrolizumab exposure enters as the per-subject covariate AUC_PEMBRO in mg*day/L)",
    concentration = "mm (the observable `TS` is the RECIST 1.1 sum of the longest diameters of target lesions)"
  )

  covariateData <- list(
    MIX_MONO_SLOW = list(
      description = "Latent mixture-class indicator: 1 = the subject belongs to the 'monophasic slow' subpopulation (steady first-order tumor shrinkage at the base kill rate), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference class is 'escape' (all three MIX_ indicators 0)",
      source_name = "MIXNUM = 2",
      notes = paste(
        "Supplementary Table 1 and the supplementary NONMEM control stream IF(MIXNUM.EQ.2) block: growth rate KL, kill rate KD, shallow compartment BASEL1, deep compartment fixed to 0.",
        "Exactly one of MIX_MONO_SLOW, MIX_BIPHASIC, MIX_MONO_FAST may be 1; all three 0 means the escape class.",
        "Class prevalence in the source cohort (main-article Table 2, P2): 0.389.",
        "For simulation, draw the class from the multinomial probabilities this model returns (prob_escape, prob_monophasic_slow, prob_biphasic, prob_monophasic_fast).",
        sep = " "
      )
    ),
    MIX_BIPHASIC = list(
      description = "Latent mixture-class indicator: 1 = the subject belongs to the 'biphasic' subpopulation (fast initial shrinkage of the accessible tumor onto a static residual plateau), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference class is 'escape' (all three MIX_ indicators 0)",
      source_name = "MIXNUM = 3",
      notes = paste(
        "Supplementary Table 1 and the supplementary NONMEM control stream IF(MIXNUM.EQ.3) block: growth rate KL, kill rate KD x KD_Factor, shallow compartment BASEL1, deep compartment BASEL2.",
        "This is the ONLY class in which the static 'deep' tumor portion BASEL2 is non-zero; it is what produces the post-nadir plateau the authors describe as uncharacteristic of conventional chemotherapy.",
        "Class prevalence in the source cohort (main-article Table 2, P3): 0.255.",
        sep = " "
      )
    ),
    MIX_MONO_FAST = list(
      description = "Latent mixture-class indicator: 1 = the subject belongs to the 'monophasic fast' subpopulation (fast first-order tumor shrinkage, no static residual), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference class is 'escape' (all three MIX_ indicators 0)",
      source_name = "MIXNUM = 4",
      notes = paste(
        "Supplementary Table 1 and the supplementary NONMEM control stream IF(MIXNUM.EQ.4) block: growth rate KL, kill rate KD x KD_Factor, shallow compartment BASEL1, deep compartment fixed to 0.",
        "Shares the kill-rate multiplier KD_Factor with MIX_BIPHASIC and differs from it only by the absence of the static deep portion.",
        "Class prevalence in the source cohort (main-article Table 2, P4): 0.0628.",
        sep = " "
      )
    ),
    NTARGET = list(
      description = "Baseline number of target lesions per RECIST 1.1.",
      units = "(count)",
      type = "continuous",
      reference_category = NULL,
      source_name = "NTARGET",
      notes = paste(
        "Power effect on the shallow baseline: TVBASEL1 * (NTARGET / 3.00)^0.654 (supplementary NONMEM control stream BAS1NTARGET block; exponent from main-article Table 2).",
        "The reference count 3.00 is the denominator in the source control stream, not a reported cohort median.",
        "The source stream maps a missing count (NTARGET = -99) to a covariate multiplier of exactly 1. Reproduce that here by setting NTARGET = 3 for a subject with no lesion count, which gives (3/3)^0.654 = 1 identically; no separate missingness column is needed.",
        "Strongest single covariate relationship in the mixture analysis (dOFV -102, P < 0.0001; supplementary Table 3).",
        sep = " "
      )
    ),
    WHO_PS = list(
      description = "Baseline ECOG (Eastern Cooperative Oncology Group) performance status.",
      units = "(integer score)",
      type = "categorical",
      reference_category = "0 (normal activity; the most frequent category)",
      source_name = "BECOGN",
      notes = paste(
        "Two-category fractional-deviation effect on the shallow baseline: BASEL1 * (1 + 0.344) for ECOG 1 versus ECOG 0 (supplementary NONMEM control stream BAS1BECOGN block; coefficient from main-article Table 2, 'Fractional change in BASEL1 for ECOG status 1 compared to ECOG status 0').",
        "Trial eligibility restricted enrolment to ECOG 0-1, so only those two levels occur; model() gates on WHO_PS == 1 and any other value falls to the reference.",
        "Selected in the second forward step of the mixture covariate search (dOFV -12.3, P = 0.0004; supplementary Table 3).",
        sep = " "
      )
    ),
    TUM_SLD = list(
      description = "Observed baseline sum of the longest diameters of target lesions per RECIST 1.1, measured at screening.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      source_name = "BASE",
      notes = paste(
        "Used ONLY by the second-stage multinomial regression on the class logits, as a linear term centered at 98.15 mm (supplementary NONMEM control stream COEFFBASE block: THETA(6)*(BASE - 98.15); coefficient from supplementary Table 2B, betaBASE = -6.72e-3).",
        "The negative coefficient means a LARGER baseline tumor lowers the odds of every responder class relative to the escape class.",
        "The tumor-size model itself does NOT take the observed baseline as an input: BASEL1 and BASEL2 are estimated population parameters with their own between-subject variability, which is the main structural difference from the paper's consolidated model (Chatterjee_2017_pembrolizumab_consolidated).",
        "The centering constant 98.15 mm is the source stream's value for the KEYNOTE-001 mixture dataset and is not printed in the article.",
        sep = " "
      )
    ),
    NNODAL = list(
      description = "Baseline number of affected (involved) lymph-node lesions.",
      units = "(count)",
      type = "continuous",
      reference_category = NULL,
      source_name = "NNODAL",
      notes = paste(
        "Used ONLY by the second-stage multinomial regression on the class logits, as an uncentered linear term (supplementary NONMEM control stream COEFFNNODAL block: THETA(5)*(NNODAL - 0.00); coefficient from supplementary Table 2B, betaNNODAL = 0.112).",
        "The positive coefficient means more involved nodes raise the odds of every responder class relative to the escape class.",
        "Selected in the second forward step of the multinomial covariate search (dOFV -15.3, P < 0.0001; supplementary Table 3, lower panel).",
        sep = " "
      )
    ),
    AUC_PEMBRO = list(
      description = "Per-subject pembrolizumab area under the serum concentration-time curve at steady state over a 6-week interval (AUCss-6weeks).",
      units = "mg*day/L (equivalently ug*day/mL; numerically identical)",
      type = "continuous",
      reference_category = NULL,
      source_name = "AUC2",
      notes = paste(
        "Linear term on the class logits, centered at 6000 mg*day/L (supplementary NONMEM control stream COEFFAUC2 block: THETA(4)/10000*(AUC2 - 6000.00); coefficient from supplementary Table 2B, betaAUC2 = 5.61e-5 per mg*day/L).",
        "In the tumor-size model exposure was tested on KD and was NOT statistically significant; the authors retained an exposure relationship only so the magnitude of a potential effect could be simulated (main-article Methods and Results). In the FINAL mixture parameterization reported in Table 2 no AUC term appears on KD at all, so exposure enters this packaged model only through the class probabilities.",
        "The source analysis did not model pembrolizumab PK here: AUC estimates came from a separate population-PK analysis (Ahamadi 2017, packaged as Ahamadi_2017_pembrolizumab). The supplementary simulation methodology computes it as AUCss-6weeks = (dose x weight / CL) x 6 / REGIMEN, with a median CL of 0.2065 L/day and a median weight of 78.9 kg in this cohort.",
        sep = " "
      )
    )
  )

  # Screened during the mixture covariate search but not retained in the final
  # model (main-article Methods, "Covariate exploration", and supplementary
  # Table 3, which records that no relationship was dropped in backward
  # elimination -- i.e. everything below failed the forward-inclusion step).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "year",
      type = "continuous",
      notes = "Tested and not selected (main-article Methods, 'Covariate exploration')."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Tested and not selected (main-article Methods). Median 78.9 kg in this cohort per the supplementary multinomial control stream. Weight still enters indirectly through the mg/kg dose that determines AUC_PEMBRO."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not selected (main-article Methods, listed as 'gender')."
    ),
    TUM_BRAF_MUT = list(
      description = "Tumor BRAF (B-type Raf) mutation indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not selected in the mixture analysis (main-article Methods). It WAS retained on kgrowth in the paper's later consolidated model; see Chatterjee_2017_pembrolizumab_consolidated."
    ),
    DIS_STAGE = list(
      description = "Melanoma disease stage.",
      units = "(categorical)",
      type = "categorical",
      notes = "Tested and not selected (main-article Methods, 'disease stage'; source column STAGEN)."
    ),
    PRIOR_IPI = list(
      description = "Prior ipilimumab treatment indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not selected in the mixture analysis (main-article Methods, 'IPI pretreatment status'; source column IPIN). It WAS retained on f in the paper's later consolidated model."
    ),
    STUDY_PART = list(
      description = "KEYNOTE-001 study part / cohort (B1, B2, B3, D).",
      units = "(categorical)",
      type = "categorical",
      notes = "Tested and not selected (main-article Methods, 'study part'; source column PRT)."
    ),
    RANDOMIZED = list(
      description = "Randomization status (randomized versus sequentially assigned).",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not selected (main-article Methods, 'randomization status'; source column RAND)."
    ),
    REGIMEN = list(
      description = "Dosing interval in weeks (2 for Q2W, 3 for Q3W).",
      units = "week",
      type = "categorical",
      notes = "Tested and not selected (main-article Methods, 'regimen'; source column REGIMEN). It survives in the source data only as the divisor that converts a steady-state interval AUC into AUCss-6weeks."
    )
  )

  population <- list(
    species = "human (adults with unresectable or metastatic melanoma)",
    n_subjects = 364L,
    n_studies = 1L,
    disease_state = "advanced (unresectable stage III or stage IV) melanoma; 168 of 364 ipilimumab-naive, the remainder previously treated with ipilimumab",
    dose_range = "pembrolizumab 10 mg/kg IV Q2W (n = 51), 10 mg/kg IV Q3W (n = 167), or 2 mg/kg IV Q3W (n = 146); not a model input, and enters only through AUC_PEMBRO",
    regions = "KEYNOTE-001 (NCT01295827), phase Ib, multicenter open-label",
    notes = paste(
      "Mixture-model dataset (main-article Methods, 'Patients included in the tumor-size model datasets'): advanced melanoma tumor-size data from KEYNOTE-001 with an October 2013 cutoff.",
      "Only patients with at least one measurable lesion at baseline who were evaluable for pharmacokinetics were included.",
      "Estimated class prevalences (main-article Table 2): escape 0.294, monophasic slow 0.389, biphasic 0.255, monophasic fast 0.0628. Table 2 footnote c states these are 'derived from estimates of the logits and corrected for the frequency of patients with missing post-baseline scans'; the raw multinomial of the Table 2 logits gives 0.182 / 0.450 / 0.295 / 0.0727, and mixing 13.6% of subjects a priori into the escape class reproduces all four reported values to three significant figures. The vignette runs that arithmetic as an assertion.",
      "Patients who dropped out before any post-baseline scan were assigned to the escape class a priori (supplementary NONMEM control stream, MSPSBSFL flag). That assignment is what motivated the two-stage covariate search: those subjects cannot contribute to the estimated class probabilities.",
      "Estimation: NONMEM 7.2.0, FOCE with interaction; bootstrap statistics in Table 2 from 695 successfully minimized replicates of 1000 (multinomial regression: 904 of 1000). Covariate search by PsN 3.5.3 stepwise covariate modeling (forward P < 0.01, backward P < 0.001).",
      "The paper reports a BASE mixture model in supplementary Table 2A; per library policy only the FINAL model (Table 2) is packaged here.",
      sep = "\n"
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters of the tumor-size model. All values from
    # main-article Table 2 ('Final mixture model'), 'Estimate' column.
    # Table 2 footnote b: KL and KD were estimated log-transformed and are
    # reported back-transformed, in units of 10^-3/day.
    # ---------------------------------------------------------------------
    lkgrowth <- log(2.76e-3)
    label("log first-order tumor growth rate constant KL (1/day)") # Table 2: KL = 2.76 x 10^-3/day, bootstrap median 2.77, 90% CI 2.26-3.9

    lkdeath <- log(3.57e-3)
    label("log first-order tumor kill rate constant KD of the slow monophasic class (1/day)") # Table 2: KD = 3.57 x 10^-3/day, bootstrap median 3.61, 90% CI 2.96-4.74

    lrbase_shallow <- log(56.3)
    label("log baseline size of the shallow (growing/shrinking) tumor portion BASEL1 (mm)") # Table 2: BASEL1 = 56.3 mm, bootstrap median 56.4, 90% CI 49.8-63

    lrbase_deep <- log(25.1)
    label("log baseline size of the deep (static) tumor portion BASEL2 (mm), biphasic class only") # Table 2: BASEL2 = 25.1 mm, bootstrap median 25.8, 90% CI 18-37.9

    # ---------------------------------------------------------------------
    # Covariate effects on the shallow baseline. Functional forms from the
    # supplementary NONMEM control stream; values from Table 2.
    # ---------------------------------------------------------------------
    e_ntarget_lrbase_shallow <- 0.654
    label("Power exponent of (NTARGET / 3) on the shallow baseline BASEL1 (unitless)") # Table 2: 'Exponent of relationship between number of target lesions and BASEL1' = 0.654, 90% CI 0.555-0.765

    e_ecog1_lrbase_shallow <- 0.344
    label("Fractional change in the shallow baseline BASEL1 for ECOG 1 versus ECOG 0 (unitless)") # Table 2: 'Fractional change in BASEL1 for ECOG status 1 compared to ECOG status 0' = 0.344, 90% CI 0.151-0.577

    # ---------------------------------------------------------------------
    # Between-class multipliers (supplementary Table 1). These are
    # MULTIPLICATIVE factors, not fractional deviations: the escape class
    # carries 2.00 x BASEL1 and the biphasic / monophasic-fast classes carry
    # 4.17 x KD. The escape class additionally has KD fixed to 0 and the
    # deep portion fixed to 0, and every class except biphasic has the deep
    # portion fixed to 0 -- those are structural zeros, not parameters.
    # ---------------------------------------------------------------------
    e_mix_escape_lrbase_shallow <- 2
    label("Multiplicative factor on the shallow baseline BASEL1 in the escape class (unitless)") # Table 2: 'Difference of BASEL1 in escape group relative to all other groups' = 2, bootstrap median 2, 90% CI 1.66-2.35

    e_mix_fast_kdeath <- 4.17
    label("Multiplicative factor on the kill rate KD in the biphasic and monophasic-fast classes (unitless)") # Table 2: 'Difference in KD in biphasic and fast monophasic groups relative to slow monophasic group' = 4.17, 90% CI 3.42-4.85

    # ---------------------------------------------------------------------
    # Second-stage multinomial logistic regression on the class logits.
    # All values from supplementary Table 2B ('Multinomial regression on the
    # subgroup probabilities'), which is the covariate-containing stage; the
    # escape class is the reference. The control stream adds ONE shared
    # covariate term to all three logits (a proportional-odds form), which
    # is reproduced exactly below.
    #
    # Cross-check of this reading: the raw multinomial of the three baseline
    # logits below gives 0.2859 / 0.3957 / 0.2630 / 0.05545, matching the
    # P1-P4 row of supplementary Table 2B (0.286 / 0.396 / 0.263 / 0.0553)
    # to three significant figures. The vignette asserts it.
    # ---------------------------------------------------------------------
    lgt_mono_slow <- 0.325
    label("Baseline logit of the monophasic-slow class versus the escape class (unitless)") # Supplementary Table 2B: BL2 = 0.325, bootstrap median 0.323, 90% CI 0.085-0.591

    lgt_biphasic <- -0.0834
    label("Baseline logit of the biphasic class versus the escape class (unitless)") # Supplementary Table 2B: BL3 = -0.0834, bootstrap median -0.0847, 90% CI -0.362 to 0.202

    lgt_mono_fast <- -1.64
    label("Baseline logit of the monophasic-fast class versus the escape class (unitless)") # Supplementary Table 2B: BL4 = -1.64, bootstrap median -1.64, 90% CI -2.13 to -1.24

    e_auc_lgt <- 5.61e-5
    label("Shared additive effect on each class logit per mg*day/L of AUCss-6weeks above 6000 (1/(mg*day/L))") # Supplementary Table 2B: beta_AUC2 = 5.61 x 10^-5, 90% CI 0.785 x 10^-5 to 11.1 x 10^-5

    e_nnodal_lgt <- 0.112
    label("Shared additive effect on each class logit per affected lymph node (unitless)") # Supplementary Table 2B: beta_NNODAL = 0.112, bootstrap median 0.116, 90% CI 0.0723-0.173

    e_tumsld_lgt <- -6.72e-3
    label("Shared additive effect on each class logit per mm of baseline tumor size above 98.15 mm (1/mm)") # Supplementary Table 2B: beta_BASE = -6.72 x 10^-3, bootstrap median -6.92 x 10^-3, 90% CI -9.11 x 10^-3 to -4.89 x 10^-3

    # ---------------------------------------------------------------------
    # Between-subject variability. Table 2 reports %CV; Table 2 footnote d
    # gives the transform CV% = 100 * sqrt(exp(omega^2) - 1), so the
    # variances below are omega^2 = log(1 + CV^2). Reading confirmed against
    # the supplementary control stream's $OMEGA initial estimates, which
    # back-transform to CV values bracketing the base and final tables.
    #
    # There is deliberately NO eta on lkgrowth: the supplementary control
    # stream declares `$OMEGA 0 FIX ; IIV_KL`, i.e. the source FIXED the
    # between-subject variance of the growth rate to zero, and Table 2
    # reports no IIV row for KL. A zero-variance eta is omitted rather than
    # written as fixed(0) because a singular OMEGA block breaks rxode2's
    # Cholesky decomposition at solve time.
    # ---------------------------------------------------------------------
    etalkdeath ~ 0.109392 # Table 2: 'IIV KD, %CV' = 34 (90% CI 25.4-42.2, shrinkage 26.7%) -> omega^2 = log(1 + 0.34^2)

    # Full 2x2 block on the two baselines. Table 2 reports the two %CV values
    # and their correlation coefficient separately; the covariance below is
    # corr * sqrt(var1 * var2) on the eta scale, which is the scale NONMEM
    # prints a correlation on.
    etalrbase_shallow + etalrbase_deep ~ c(
      0.490796,
      0.811621, 1.802122
    ) # Table 2: 'IIV BL1, %CV' = 79.6 (shrinkage 2.8%), 'IIV BL2, %CV' = 225 (shrinkage 8.1%), 'Corr, BL1~BL2' = 0.863 -> vars log(1+0.796^2), log(1+2.25^2), cov 0.863*sqrt(product)

    # Between-subject scaling of the residual-error magnitude (NONMEM's
    # eta-on-epsilon). Anchor fixed at 1 so the typical subject carries the
    # Table 2 residual exactly; the eta scales it log-normally.
    lrv <- fixed(log(1))
    label("Residual-variability scaling anchor (unitless)") # structural anchor so etalrv scales the Table 2 residual; not a source value

    etalrv ~ 0.065901 # Table 2: 'ETA_EPS' = 26.1 (90% CI 15.9-33.1, shrinkage 28.9%) -> omega^2 = log(1 + 0.261^2). Supplementary Table 2A heads this row '%CV'; the Table 2 heading '(variance)' is inconsistent with both that and with the control stream's $OMEGA 0.0625 initial, which back-transforms to 25.4% CV.

    # ---------------------------------------------------------------------
    # Residual error. The source fits log-transformed tumor size with a
    # combined error whose standard deviation on the LOG scale is
    #   W = sqrt(propSd^2 + addSd^2 / SLD^2)
    # (supplementary NONMEM control stream: W = SQRT(THETA(10)**2 +
    # THETA(11)**2 / IPRED1**2), $SIGMA 1 FIX), i.e. a proportional and an
    # additive component in linear SLD space.
    # ---------------------------------------------------------------------
    propSd <- 0.103
    label("Proportional component of the residual error on tumor size (unitless)") # Table 2: 'Residual error, Proportional, %CV' = 10.3 (bootstrap median 10.2, 90% CI 8.84-11.8)

    addSd <- 3.29
    label("Additive component of the residual error on tumor size (mm)") # Table 2: 'Residual error, Additive, mm' = 3.29 (bootstrap median 3.24, 90% CI 2.65-4.1)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Latent class bookkeeping. Exactly one of the three indicators is
    #    1, or all three are 0 for the reference 'escape' class.
    # -------------------------------------------------------------------
    mix_escape <- 1 - MIX_MONO_SLOW - MIX_BIPHASIC - MIX_MONO_FAST

    # -------------------------------------------------------------------
    # 2. Individual tumor-dynamics parameters
    # -------------------------------------------------------------------
    # Growth rate: common to all four classes, no between-subject variance
    # (source fixed IIV_KL to 0).
    kgrowth <- exp(lkgrowth)

    # Kill rate: 0 in the escape class, the base KD in the monophasic-slow
    # class, and KD x 4.17 in the biphasic and monophasic-fast classes
    # (supplementary Table 1).
    kdeath_factor <- MIX_MONO_SLOW + e_mix_fast_kdeath * (MIX_BIPHASIC + MIX_MONO_FAST)
    kdeath <- exp(lkdeath + etalkdeath) * kdeath_factor

    # Shallow baseline: covariate-adjusted, then scaled by 2.00 in the
    # escape class. A missing lesion count should be carried as NTARGET = 3,
    # which makes the power term exactly 1 (see covariateData notes).
    base_cov <- (NTARGET / 3)^e_ntarget_lrbase_shallow *
      (1 + e_ecog1_lrbase_shallow * (WHO_PS == 1))
    base_class <- 1 + (e_mix_escape_lrbase_shallow - 1) * mix_escape
    basel_shallow <- exp(lrbase_shallow + etalrbase_shallow) * base_cov * base_class

    # Deep (static) baseline: present in the biphasic class only.
    basel_deep <- exp(lrbase_deep + etalrbase_deep) * MIX_BIPHASIC

    # -------------------------------------------------------------------
    # 3. Second-stage multinomial logistic regression on the class logits
    #    (supplementary Table 2B). One shared covariate term is added to
    #    all three non-reference logits, exactly as the source control
    #    stream does. These are derived outputs, not observations: use them
    #    to assign MIX_MONO_SLOW / MIX_BIPHASIC / MIX_MONO_FAST before
    #    solving the tumor-size model.
    # -------------------------------------------------------------------
    cov_lgt <- e_auc_lgt * (AUC_PEMBRO - 6000) +
      e_nnodal_lgt * NNODAL +
      e_tumsld_lgt * (TUM_SLD - 98.15)
    lgt_slow <- lgt_mono_slow + cov_lgt
    lgt_bi <- lgt_biphasic + cov_lgt
    lgt_fast <- lgt_mono_fast + cov_lgt
    lgt_denom <- 1 + exp(lgt_slow) + exp(lgt_bi) + exp(lgt_fast)

    prob_escape <- 1 / lgt_denom
    prob_monophasic_slow <- exp(lgt_slow) / lgt_denom
    prob_biphasic <- exp(lgt_bi) / lgt_denom
    prob_monophasic_fast <- exp(lgt_fast) / lgt_denom

    # -------------------------------------------------------------------
    # 4. Tumor size and residual error
    # -------------------------------------------------------------------
    # Main-article Methods:
    #   y_i(t) = BASEL1_i * exp(KL * t - KD_i * t) + BASEL2_i
    # No ODE state is needed; the source itself is a NONMEM $PRED model.
    TS <- basel_shallow * exp((kgrowth - kdeath) * time) + basel_deep

    # Log-scale residual SD, scaled per subject by exp(etalrv). The source
    # additionally floors W at 0.01 and the prediction at 1e-5; neither
    # floor can bind here because propSd = 0.103 > 0.01 and the prediction
    # is a sum of positive exponentials.
    sdlog <- sqrt(propSd^2 + addSd^2 / TS^2) * exp(lrv + etalrv)

    TS ~ lnorm(sdlog)
  })
}
attr(Chatterjee_2017_pembrolizumab_mixture, "message") <-
  "Four-class mixture tumor-size model for pembrolizumab in advanced melanoma (Chatterjee 2017; KEYNOTE-001, n = 364). Observable `TS` is the RECIST 1.1 sum of longest diameters in mm, TS(t) = BASEL1 * exp((KL - KD) * t) + BASEL2. Assign each subject a latent class with the binary covariates MIX_MONO_SLOW / MIX_BIPHASIC / MIX_MONO_FAST (all three 0 = the reference 'escape' class); the file also returns the paper's second-stage multinomial class probabilities prob_escape / prob_monophasic_slow / prob_biphasic / prob_monophasic_fast for that assignment. Required covariates: MIX_MONO_SLOW, MIX_BIPHASIC, MIX_MONO_FAST, NTARGET, WHO_PS, TUM_SLD, NNODAL, AUC_PEMBRO. No PK input: supply exposure per subject as AUC_PEMBRO (mg*day/L), which the source computed as dose/CL from the companion popPK model packaged here as Ahamadi_2017_pembrolizumab. Exposure affects only the class probabilities, and the paper's conclusion is that response is flat across 2-10 mg/kg."
Chatterjee_2017_pembrolizumab_mixture
