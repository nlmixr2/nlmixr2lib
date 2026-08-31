Hanan_2026_peginterferon_alfa_24wk_mbma <- function() {
  description <- paste0(
    "MBMA. Logistic model-based meta-analysis of hepatitis B surface ",
    "antigen (HBsAg) loss 24 WEEKS AFTER cessation of pegylated-",
    "interferon-alpha (Peg-IFNalpha) - the paper's surrogate for ",
    "functional cure - in chronic hepatitis B virus infection, fit to ",
    "58 study-strata-arms (28 studies, 4267 participants) identified ",
    "by a PRISMA systematic literature review of Embase / MEDLINE / ",
    "Cochrane, January 2000 - July 2022 (Hanan 2026 Table 1). HBsAg ",
    "loss is a binary endpoint (HBsAg below the assay limit of ",
    "detection); the study-strata-arm probability is logit-linear in a ",
    "regimen-specific intercept plus two retained covariates - ",
    "Peg-IFNalpha treatment duration (centered at 48 weeks) and ",
    "baseline HBsAg (centered at 3 log10 IU/mL) - with a ",
    "between-study-strata-arm random effect on the logit scale (Hanan ",
    "2026 Section 2.3.3 final covariate model). Baseline HBeAg status ",
    "was screened here and NOT retained, which is the one structural ",
    "difference from the paper's end-of-treatment model. Simulation ",
    "scope is the AGGREGATE study-arm HBsAg-loss proportion for a new ",
    "trial, NOT an individual patient outcome and NOT a ",
    "concentration; the model contains no PK, no ODE state and no dose ",
    "events, and the 'time' axis is not used. The companion model for ",
    "the same paper's first, independently fit endpoint (HBsAg loss at ",
    "end of treatment) is Hanan_2026_peginterferon_alfa_eot_mbma."
  )

  reference <- paste(
    "Hanan NJ, Zierhut ML, Nader A, Mahajan A, Kaur A, Kumar K,",
    "Dixon SA, Das J, Magee M, Theodore D, Gielen V. A Systematic",
    "Review and Model-Based Meta-Analysis of",
    "Pegylated-Interferon-alpha-Induced HBsAg Loss in Chronic",
    "Hepatitis B Virus Infection. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15:e70164. doi:10.1002/psp4.70164.",
    sep = " "
  )

  vignette <- "Hanan_2026_peginterferon_alfa_hbsag_loss"

  # Between-STUDY-STRATA-ARM effect on the logit-scale response (Hanan
  # 2026 Section 2.3.1, eta_isa ~ N(0, sigma^2_eta)), not a
  # between-SUBJECT effect on a structural PK parameter. Declared
  # paper-specific per the Goteti_2024_SLE_mbma precedent.
  paper_specific_etas <- c("eta_study")

  units <- list(
    time          = "week",
    dosing        = paste(
      "(no rxode2 dose events; the Peg-IFNalpha regimen enters only as",
      "the binary CONMED_NUC combination indicator selecting the",
      "regimen intercept, and as the T_PEGIFN planned treatment",
      "duration covariate. Non-standard Peg-IFNalpha dosing was an",
      "exclusion criterion, so dose amount is not a model input.)"
    ),
    concentration = paste(
      "(dimensionless proportion, 0-1: p_hbsag_loss, the study-arm",
      "probability of HBsAg loss 24 weeks after Peg-IFNalpha",
      "cessation. Aggregate study-arm summary, not an",
      "individual-subject concentration.)"
    )
  )

  covariateData <- list(
    HBSAG_BL_LOG10 = list(
      description        = paste(
        "Study-strata-arm baseline (pre-treatment) hepatitis B surface",
        "antigen concentration on the log10 IU/mL scale. Enters the",
        "logit-scale response centered at 3 log10 IU/mL."
      ),
      units              = "log10 IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-strata-arm-level covariate: the aggregate baseline",
        "HBsAg summary for the arm, not an individual measurement.",
        "Hanan 2026 Section 2.2 derivation: the reported MEDIAN log10",
        "IU/mL was used when available, the mean when only a mean was",
        "reported, and the stratum midpoint when baseline HBsAg was",
        "reported as tightly defined strata (e.g. 0-500, 500-1000,",
        "1000-1500 IU/mL). Studies reporting only broad non-continuous",
        "categories (e.g. <3000 vs >3000 IU/mL) were EXCLUDED as",
        "unreliable to impute.",
        "",
        "Centering value 3 log10 IU/mL is the across-arm population",
        "median (Hanan 2026 Table 2 row 'Baseline HBsAg (per log10",
        "IU/mL, ref.=3)'). Each 1 log10 IU/mL increase multiplies the",
        "odds of HBsAg loss by exp(-1.58) = 0.206, i.e. reduces the",
        "odds by 79% (Hanan 2026 Results Section 3.3) - within",
        "rounding of the -1.57 estimated in the EOT model, which the",
        "paper notes as 'the magnitudes of covariate effects were",
        "similar across models'.",
        "",
        "The log10 scale is explicit in the name so that a raw IU/mL",
        "column is never silently substituted; the paper itself quotes",
        "raw IU/mL values elsewhere (e.g. 2485 IU/mL and 2081 IU/mL",
        "virtual-cohort baselines in Section 3.5), which are",
        "log10 = 3.395 and 3.318. Same discipline as the",
        "BACT_PTT_LOG10CFU register entry."
      ),
      source_name        = "Baseline HBsAg (log10 IU/mL)"
    ),
    T_PEGIFN = list(
      description        = paste(
        "Planned duration of the Peg-IFNalpha treatment course for the",
        "study-strata-arm, in weeks. Enters the logit-scale response",
        "centered at 48 weeks."
      ),
      units              = "week",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MBMA study-strata-arm-level covariate. Centering value 48",
        "weeks is the standard Peg-IFNalpha regimen and the reference",
        "used in Hanan 2026 Table 2 ('Peg-IFNalpha treatment duration",
        "(per week, ref.=48 weeks)'). Each additional week of",
        "Peg-IFNalpha increases the odds of HBsAg loss 24 weeks",
        "post-treatment by 3.1% (exp(0.0303) = 1.0308; Hanan 2026",
        "Results Section 3.3) - more than twice the 1.3% per week seen",
        "at end of treatment.",
        "",
        "This is the planned/administered treatment DURATION of the",
        "arm, a design characteristic - it is not elapsed time on the",
        "model integration axis and it is not a per-subject",
        "time-varying covariate. In particular it is NOT the 24-week",
        "post-treatment follow-up interval, which is fixed by the",
        "endpoint definition and is not a covariate. The model has no",
        "time axis at all (Hanan 2026 Section 2.2: two cross-sectional",
        "models were developed rather than a longitudinal trajectory",
        "model).",
        "",
        "Member of the T_<event> duration family alongside T_CPB",
        "(cardiopulmonary bypass duration) and T_PUMP (breast-milk",
        "expression session duration)."
      ),
      source_name        = "Peg-IFNalpha treatment duration (weeks)"
    ),
    CONMED_NUC = list(
      description        = paste(
        "1 = the arm received Peg-IFNalpha PLUS a nucleos(t)ide",
        "analogue (NA) as combination therapy, 0 = Peg-IFNalpha",
        "monotherapy. Selects which of the two regimen-specific",
        "logit-scale intercepts applies."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Peg-IFNalpha monotherapy)",
      notes              = paste(
        "MBMA study-strata-arm-level covariate. Hanan 2026 Section",
        "2.3.2 extends the base model to a regimen-specific fixed",
        "effect drug.eff_d distinguishing Peg-IFNalpha monotherapy from",
        "Peg-IFNalpha + NA combination therapy; this indicator selects",
        "d. It is a treatment-arm definition, NOT an additive",
        "covariate effect, which is why the model file carries two",
        "intercepts (logitp_mono, logitp_comb) rather than one",
        "intercept plus an e_conmed_nuc_* coefficient - that is the",
        "form the paper estimated and reports in Table 2.",
        "",
        "The class token NUC is used rather than the paper's",
        "abbreviation NA because NA is R's missing-value literal.",
        "Class-level member of the CONMED_ family alongside",
        "CONMED_AZOLE, CONMED_ABX, CONMED_AED and CONMED_DIURETIC,",
        "which are likewise drug classes rather than single INNs.",
        "",
        "This indicator marks NA co-administration DURING the",
        "Peg-IFNalpha course. Whether the study protocol permitted NA",
        "to be CONTINUED after Peg-IFNalpha cessation is a separate",
        "covariate that Hanan 2026 screened in this 24-week model only,",
        "and did not retain; see covariatesDataExcluded below."
      ),
      source_name        = "Treatment group (Peg-IFNalpha vs Peg-IFNalpha + NA)"
    )
  )

  # Covariates that Hanan 2026 SCREENED for association with 24-week
  # post-treatment HBsAg loss but did NOT retain in the final model
  # (Results Section 3.3: "In the 24-week post-treatment model, two
  # covariates were significant: baseline HBsAg (p < 0.0001) and
  # Peg-IFNalpha treatment duration (p = 0.0002). No significant
  # associations were observed for age, gender, race, study design, or
  # continuation of NA."). The paper reports no point estimate for any
  # of them, so there is nothing to encode - they are recorded here to
  # preserve the provenance of the covariate screen without carrying a
  # "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    HBEAG_POS = list(
      description        = paste(
        "1 = baseline hepatitis B e-antigen (HBeAg) positive,",
        "0 = HBeAg negative."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HBeAg negative)",
      notes              = paste(
        "THE ONE STRUCTURAL DIFFERENCE between this model and the",
        "companion EOT model. Baseline HBeAg status was among the",
        "eight covariates screened in both models; it was retained at",
        "end of treatment (OR 2.37, logit 0.863) but was NOT retained",
        "here, and Hanan 2026 reports no point estimate for it in the",
        "24-week model. Consequently the 24-week reference population",
        "is defined WITHOUT an HBeAg qualifier - Table 2 gives it as",
        "'baseline HBsAg 3 log10 IU/mL, 48 weeks Peg-IFNalpha' - and",
        "this model predicts the same probability for HBeAg-positive",
        "and HBeAg-negative arms at equal baseline HBsAg and duration.",
        "",
        "Do NOT carry the EOT model's 0.863 coefficient over to this",
        "endpoint; the paper did not estimate it here.",
        "",
        "See covariateData$HBEAG_POS in",
        "Hanan_2026_peginterferon_alfa_eot_mbma for the retained form."
      ),
      source_name        = "Baseline HBeAg status"
    ),
    CONMED_NUC_CONTINUED = list(
      description        = paste(
        "1 = the study protocol permitted nucleos(t)ide analogue (NA)",
        "therapy to be CONTINUED after Peg-IFNalpha cessation,",
        "0 = NA was stopped at Peg-IFNalpha cessation."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NA stopped at Peg-IFNalpha cessation)",
      notes              = paste(
        "Screened in the 24-WEEK MODEL ONLY (it is undefined at end of",
        "treatment) and not retained; no point estimate reported.",
        "Hanan 2026 Section 2.2 explains the rationale for testing it:",
        "some patients achieve HBsAg loss following NA cessation, which",
        "could influence sustained HBsAg loss post-Peg-IFNalpha. Coded",
        "from whether study protocols permitted NA continuation, not",
        "from observed per-patient behaviour.",
        "",
        "CONMED_NUC_CONTINUED is a descriptive name used for",
        "documentation only - it is not a ratified entry in",
        "inst/references/covariate-columns.md, because the covariate is",
        "never referenced in model()."
      ),
      source_name        = "Continuation of NA after Peg-IFNalpha cessation"
    ),
    AGE = list(
      description        = "Study-strata-arm mean or median age.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as one of the eight candidate covariates (Hanan 2026",
        "Results Section 3.3) and not retained; no point estimate",
        "reported. Centered to the across-arm population median per",
        "Section 2.2, but that median is not published."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Study-strata-arm proportion female.",
      units              = "(proportion, 0-1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The paper's 'gender' covariate, screened and not retained;",
        "no point estimate reported. Hanan 2026 Section 2.2 lists",
        "gender among the binary covariates, which at the",
        "study-strata-arm level is the arm's sex composition."
      ),
      source_name        = "Gender"
    ),
    RACE_ASIAN_PCT = list(
      description        = paste(
        "Study-strata-arm racial composition. Hanan 2026 Section 2.2",
        "defines its 'race' covariate as the proportion of NON-Asian",
        "participants, ranging from 0 (100% Asian) to 1 (100%",
        "non-Asian) - the complement of this canonical."
      ),
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained; no point estimate reported. NOTE",
        "THE ORIENTATION: the source column is the proportion",
        "NON-Asian, so RACE_ASIAN_PCT = 100 * (1 - source value). The",
        "Discussion attributes the null result to the predominance of",
        "Asian studies in the analysis set, where HBV genotypes B and",
        "C are most common, leaving little contrast to estimate."
      ),
      source_name        = "Race (proportion non-Asian)"
    ),
    DESIGN_RCT = list(
      description        = paste(
        "1 = randomized controlled trial, 0 = non-RCT (prospective",
        "cohort, retrospective cohort, or single-arm trial)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-RCT)",
      notes              = paste(
        "Trial-design covariate screened and not retained; no point",
        "estimate reported. 43 of the 58 24-week study-strata-arms",
        "were from RCTs, 8 from prospective cohorts, 3 from",
        "retrospective cohorts and 3 from single-arm trials (Hanan",
        "2026 Results Section 3.1). DESIGN_RCT is a descriptive name",
        "used for documentation only - it is not a ratified entry in",
        "inst/references/covariate-columns.md, because the covariate is",
        "never referenced in model()."
      ),
      source_name        = "Study design (RCT vs non-RCT)"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 4267L,
    n_studies            = 28L,
    n_study_strata       = 35L,
    n_study_strata_arms  = 58L,
    age_range            = "adults (Hanan 2026 Table S1 PICOS: age group 'Adults')",
    disease_state        = "chronic hepatitis B virus (HBV) infection",
    dose_range           = paste(
      "standard Peg-IFNalpha regimens only (studies investigating",
      "non-standard Peg-IFNalpha dosing were excluded); Peg-IFNalpha",
      "monotherapy or Peg-IFNalpha + nucleos(t)ide analogue",
      "combination therapy"
    ),
    regions              = paste(
      "not tabulated per arm; the Discussion states the analysis set is",
      "predominantly Asian studies, where HBV genotypes B and C",
      "predominate"
    ),
    notes                = paste(
      "Aggregate published data, not individual patient data. The",
      "observation unit is the STUDY-STRATA-ARM: a treatment arm within",
      "a study, optionally stratified by HBeAg status, baseline HBsAg",
      "or baseline HBV DNA (Hanan 2026 Section 2.2). Of the 58",
      "24-week study-strata-arms, 43 were from RCTs, 8 from",
      "prospective cohorts, 3 from retrospective cohorts and 3 from",
      "single-arm trials; three studies provided results stratified by",
      "HBeAg status, baseline HBV DNA and baseline HBsAg (Section 3.1,",
      "Table 1).",
      "",
      "The endpoint is HBsAg loss 24 weeks after cessation of",
      "Peg-IFNalpha, used as a surrogate for functional cure (sustained",
      "HBsAg and HBV DNA loss for 6 months after cessation of all",
      "antiviral therapy). Results at other post-treatment time points",
      "(e.g. 3 months) were excluded.",
      "",
      "This analysis set OVERLAPS the paper's end-of-treatment set -",
      "some studies contributed to both subsets if they reported both",
      "endpoints - but the two models were fit independently (Section",
      "2.3), which is why they are separate model files. Too few trials",
      "reported both endpoints, or reported intermediate or",
      "longer-term follow-up, to support a longitudinal trajectory",
      "model.",
      "",
      "Analysis-set exclusions (Section 2.2 and Table S1): per-protocol",
      "only results (ITT / mITT prioritized), non-standard Peg-IFNalpha",
      "dosing, non-standard HBsAg-loss definitions, missing baseline",
      "HBeAg status, and baseline HBsAg reported only as broad",
      "non-continuous categories. Table S1 additionally excludes",
      "pregnant or post-partum women, HBV/HIV coinfection, post",
      "liver-transplant and immunosuppressed participants, so the model",
      "does not describe those groups.",
      "",
      "Missingness for HBV genotype, baseline HBV DNA and baseline ALT",
      "prevented those covariates from being evaluated at all",
      "(Discussion), which the authors name as the principal limitation."
    )
  )

  ini({
    # ================================================================
    # Hanan 2026 Section 2.3.3 final covariate model, 24-week
    # post-treatment endpoint:
    #
    #   p_isa = logit^-1( drug.eff_d
    #                     + sum_m beta_m * theta_m,isa
    #                     + eta_isa )
    #
    # for study i, stratum s, treatment arm a. All values below are
    # from Hanan 2026 Table 2, "24-week post-treatment model" block,
    # "Estimate (logit)" column. The paper's "Transformed estimate"
    # column is the back-transform reproduced in model() and asserted
    # in the vignette.
    # ================================================================

    logitp_mono <- -2.08
    label(paste(
      "Peg-IFNalpha monotherapy logit-scale 24-week post-treatment",
      "HBsAg-loss intercept in the reference population (baseline",
      "HBsAg 3 log10 IU/mL, 48 weeks of Peg-IFNalpha; no HBeAg",
      "qualifier - HBeAg was not retained in this model);",
      "back-transform logit^-1(-2.08) = 11.1%"
    ))
    # Table 2, 24-week post-treatment model, "Drug effect
    # (Peg-IFNalpha)": estimate (logit) -2.08, transformed estimate
    # 11.1%, %RSE 9.1, 95% PI 2.7-36.6, p < 0.001. Interpretation
    # column: "11.1% HBsAg loss in reference population (baseline HBsAg
    # 3 log10 IU/mL, 48 weeks Peg-IFNalpha)."

    logitp_comb <- -1.81
    label(paste(
      "Peg-IFNalpha + nucleos(t)ide analogue logit-scale 24-week",
      "post-treatment HBsAg-loss intercept in the reference",
      "population; back-transform logit^-1(-1.81) = 14.1%"
    ))
    # Table 2, 24-week post-treatment model, "Drug effect
    # (Peg-IFNalpha +NA)": estimate (logit) -1.81, transformed estimate
    # 14.1%, %RSE 11.9, 95% PI 3.4-43.4, p < 0.001. (Two typographical
    # slips in this row's Interpretation cell - a doubled percent sign
    # "14.1% %" and a regimen label reading "48 weeks Peg-IFNalpha"
    # where the row is the Peg-IFNalpha + NA arm - are cosmetic and do
    # not affect any estimate; see the vignette Errata.)

    e_t_pegifn_logitp <- 0.0303
    label(paste(
      "Peg-IFNalpha treatment-duration effect on the logit-scale",
      "24-week post-treatment HBsAg-loss probability (per week above",
      "the 48-week reference); OR 1.03 per week"
    ))
    # Table 2, 24-week post-treatment model, "Peg-IFNalpha treatment
    # duration (per week, ref.=48 weeks)": estimate (logit) 0.0303,
    # transformed estimate OR 1.03, %RSE 24.6, 95% PI 1.02-1.05,
    # p < 0.001 (Results Section 3.3 gives the same term as
    # p = 0.0002).

    e_hbsag_bl_log10_logitp <- -1.58
    label(paste(
      "Baseline HBsAg effect on the logit-scale 24-week post-treatment",
      "HBsAg-loss probability (per log10 IU/mL above the 3 log10 IU/mL",
      "reference); OR 0.21 per log10 IU/mL"
    ))
    # Table 2, 24-week post-treatment model, "Baseline HBsAg (per log10
    # IU/mL, ref.=3)": estimate (logit) -1.58, transformed estimate
    # OR 0.21, %RSE 13.0, 95% PI 0.14-0.30, p < 0.001. Results Section
    # 3.3: "For each log10 unit increase in baseline HBsAg, the odds of
    # achieving HBsAg loss decreased by 79%" (1 - exp(-1.58) = 1 -
    # 0.206 = 0.794).

    # ================================================================
    # Between-study-strata-arm random effect, on the logit scale.
    # NOTE: this is BETWEEN-TRIAL variance, not popPK between-subject
    # variance. It describes how much a NEW STUDY-ARM's underlying
    # HBsAg-loss log-odds differs from the typical arm; it says nothing
    # about individual patients within an arm. It is larger here than
    # in the EOT model (0.572 vs 0.484) despite a much larger
    # proportional reduction from its base model, because the 24-week
    # base-model heterogeneity was itself far larger (2.544 vs 1.112).
    # ================================================================
    eta_study ~ 0.572
    # Table 2, 24-week post-treatment model, "Between-trial variance":
    # 0.572 (variance on the logit scale). Reduced from the base
    # model's 2.544 by the treatment and covariate effects (Results
    # Section 3.3: "from 2.544 to 0.572 in the 24-week model, a 77.5%
    # reduction"). The Abstract and Discussion instead quote 77.6% for
    # the same pair of variances; (2.544 - 0.572) / 2.544 = 77.5%, so
    # the Results figure is the arithmetically correct one - see the
    # vignette Errata.

    # ================================================================
    # Residual error. Hanan 2026 Section 2.3.1 fixes the residual
    # variance to sigma^2 = 1 and applies a binomial variance function
    # w_isa = sqrt(p_isa * (1 - p_isa) / N_isa), so the residual acts
    # as a standardized binomial sampling error whose magnitude is set
    # by the arm size N_isa rather than by an estimated parameter.
    # There is therefore no estimated residual-error parameter to
    # encode. addSd is fixed at a negligible value purely for
    # nlmixr2 endpoint compatibility (same device as
    # Goteti_2024_SLE_mbma). model() exposes w_binom_unit =
    # sqrt(p*(1-p)) so a user can recover the paper's weight as
    # w_binom_unit / sqrt(N_isa). The paper's own Table 3 simulations
    # do NOT include this term - they resample the fixed effects and
    # eta_study only (Section 2.5).
    # ================================================================
    addSd <- fixed(0.001)
    label("Placeholder additive residual error (see the note above; not a paper-estimated parameter)")
  })

  model({
    # -------- Centered covariates (Hanan 2026 Section 2.3.3) --------
    # Continuous covariates enter centered at the across-arm median so
    # that the drug.eff_d intercepts describe a "typical"
    # study-strata-arm.
    dur_dev   <- T_PEGIFN - 48
    hbsag_dev <- HBSAG_BL_LOG10 - 3

    # -------- Regimen-specific intercept (Section 2.3.2) ------------
    # drug.eff_d: one fixed effect per regimen d, selected by the
    # combination-therapy indicator rather than estimated as a
    # contrast.
    drug_eff <- logitp_mono * (1 - CONMED_NUC) +
                logitp_comb * CONMED_NUC

    # -------- Final covariate model on the logit scale --------------
    # No HBeAg term: baseline HBeAg status was screened in this model
    # and not retained (Results Section 3.3).
    logit_p <- drug_eff +
               e_t_pegifn_logitp       * dur_dev   +
               e_hbsag_bl_log10_logitp * hbsag_dev +
               eta_study

    # -------- Study-arm HBsAg-loss probability ----------------------
    p_hbsag_loss <- exp(logit_p) / (1 + exp(logit_p))

    # -------- Binomial sampling-error weight (Section 2.3.1) --------
    # The paper's variance function is w_isa = sqrt(p*(1-p)/N_isa).
    # N_isa is a property of the arm being simulated rather than of the
    # model, so the N-free factor is exposed here and the user divides
    # by sqrt(N_isa).
    w_binom_unit <- sqrt(p_hbsag_loss * (1 - p_hbsag_loss))

    # -------- Observation -------------------------------------------
    # Cc is the aggregate study-arm HBsAg-loss proportion (0-1), NOT a
    # concentration. p_hbsag_loss and w_binom_unit are also emitted as
    # derived variables in the rxSolve output.
    Cc <- p_hbsag_loss
    Cc ~ add(addSd)
  })
}
