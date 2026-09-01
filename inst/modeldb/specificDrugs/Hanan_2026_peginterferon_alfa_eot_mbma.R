Hanan_2026_peginterferon_alfa_eot_mbma <- function() {
  description <- paste0(
    "MBMA. Logistic model-based meta-analysis of hepatitis B surface ",
    "antigen (HBsAg) loss at END OF TREATMENT (EOT) with pegylated-",
    "interferon-alpha (Peg-IFNalpha)-based regimens in chronic ",
    "hepatitis B virus infection, fit to 83 study-strata-arms (45 ",
    "studies, 11,493 participants) identified by a PRISMA systematic ",
    "literature review of Embase / MEDLINE / Cochrane, January 2000 - ",
    "July 2022 (Hanan 2026 Table 1). HBsAg loss is a binary endpoint ",
    "(HBsAg below the assay limit of detection); the study-strata-arm ",
    "probability is logit-linear in a regimen-specific intercept plus ",
    "three retained covariates - Peg-IFNalpha treatment duration ",
    "(centered at 48 weeks), baseline HBsAg (centered at 3 log10 ",
    "IU/mL) and baseline HBeAg status (reference negative) - with a ",
    "between-study-strata-arm random effect on the logit scale (Hanan ",
    "2026 Section 2.3.3 final covariate model). Simulation scope is ",
    "the AGGREGATE study-arm HBsAg-loss proportion for a new trial, ",
    "NOT an individual patient outcome and NOT a concentration; the ",
    "model contains no PK, no ODE state and no dose events, and the ",
    "'time' axis is not used. The companion model for the same paper's ",
    "second, independently fit endpoint (HBsAg loss 24 weeks after ",
    "Peg-IFNalpha cessation) is ",
    "Hanan_2026_peginterferon_alfa_24wk_mbma."
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

  # The single random effect is a between-STUDY-STRATA-ARM effect on the
  # logit-scale response (Hanan 2026 Eq. in Section 2.3.1, eta_isa ~
  # N(0, sigma^2_eta)), not a between-SUBJECT effect on a structural PK
  # parameter, so it does not follow the popPK eta<transformed-param>
  # convention. Declared paper-specific per the Goteti_2024_SLE_mbma
  # precedent.
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
      "probability of HBsAg loss at end of treatment. Aggregate",
      "study-arm summary, not an individual-subject concentration.)"
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
        "odds of HBsAg loss by exp(-1.57) = 0.208, i.e. reduces the",
        "odds by 79% (Hanan 2026 Results Section 3.3).",
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
    HBEAG_POS = list(
      description        = paste(
        "1 = baseline hepatitis B e-antigen (HBeAg) positive,",
        "0 = HBeAg negative. For a study-strata-arm this is the",
        "stratum's HBeAg status (studies not reporting baseline HBeAg",
        "status were excluded from the analysis set)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HBeAg negative)",
      notes              = paste(
        "MBMA study-strata-arm-level covariate. Hanan 2026 Table 2",
        "reports the effect as 'HBeAg status (ref=negative)' with a",
        "logit-scale coefficient of 0.863, i.e. OR 2.37 (95% PI",
        "1.34-4.17, p = 0.008): HBeAg-positive arms show a more than",
        "twofold higher odds of HBsAg loss at EOT. The Table 2",
        "Interpretation column states the resulting reference-population",
        "rates directly - 18.1% vs 8.5% for Peg-IFNalpha monotherapy",
        "and 23.8% vs 11.7% for Peg-IFNalpha + NA - and the model file",
        "reproduces both exactly.",
        "",
        "Retained ONLY in the EOT model. The companion 24-week",
        "post-treatment model",
        "(Hanan_2026_peginterferon_alfa_24wk_mbma) screened the same",
        "covariate and did NOT retain it, so it appears there in",
        "covariatesDataExcluded rather than covariateData.",
        "",
        "This is a disease-PHASE / serologic marker within the primary",
        "indication (chronic HBV), not a coinfection flag - contrast",
        "HCV_POS / HIV_POS / TB_POS, which mark comorbid infection on a",
        "non-hepatitis-B primary indication.",
        "",
        "When simulating a MIXED arm (a cohort containing both",
        "HBeAg-positive and HBeAg-negative participants) note that the",
        "model is logit-linear in this indicator, so the arm-level",
        "proportion is the HBeAg-weighted average of the two subgroup",
        "probabilities, NOT the probability evaluated at the mean",
        "indicator value. See the vignette."
      ),
      source_name        = "Baseline HBeAg status"
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
        "Peg-IFNalpha increases the odds of HBsAg loss at EOT by 1.3%",
        "(exp(0.0133) = 1.0134; Hanan 2026 Results Section 3.3).",
        "",
        "This is the planned/administered treatment DURATION of the",
        "arm, a design characteristic - it is not elapsed time on the",
        "model integration axis and it is not a per-subject",
        "time-varying covariate. The model has no time axis at all",
        "(both Hanan 2026 models are cross-sectional; Section 2.2",
        "states that two cross-sectional models were developed rather",
        "than a longitudinal trajectory model because too few trials",
        "reported both endpoints or intermediate follow-up).",
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
        "Participants in the combination arms may have stopped OR",
        "continued NA after Peg-IFNalpha cessation (Hanan 2026 Section",
        "2.2). That distinction is a SEPARATE covariate that was",
        "screened only in the 24-week model and was not retained; see",
        "covariatesDataExcluded in",
        "Hanan_2026_peginterferon_alfa_24wk_mbma."
      ),
      source_name        = "Treatment group (Peg-IFNalpha vs Peg-IFNalpha + NA)"
    )
  )

  # Covariates that Hanan 2026 SCREENED for association with EOT HBsAg
  # loss but did NOT retain in the final model (Results Section 3.3:
  # "No significant associations were observed for age, gender, race,
  # study design, or continuation of NA after Peg-IFNalpha"). The paper
  # reports no point estimate for any of them, so there is nothing to
  # encode - they are recorded here to preserve the provenance of the
  # covariate screen without carrying a "declared but not referenced"
  # convention warning.
  covariatesDataExcluded <- list(
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
        "estimate reported. 63 of the 83 EOT study-strata-arms were",
        "from RCTs (Hanan 2026 Results Section 3.1). DESIGN_RCT is a",
        "descriptive name used for documentation only - it is not a",
        "ratified entry in inst/references/covariate-columns.md,",
        "because the covariate is never referenced in model()."
      ),
      source_name        = "Study design (RCT vs non-RCT)"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 11493L,
    n_studies            = 45L,
    n_study_strata       = 52L,
    n_study_strata_arms  = 83L,
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
      "or baseline HBV DNA (Hanan 2026 Section 2.2). Of the 83 EOT",
      "study-strata-arms, 63 were from RCTs and the remainder from",
      "prospective cohorts, retrospective cohorts or single-arm trials;",
      "four studies provided results stratified by HBeAg status and",
      "baseline HBsAg (Section 3.1, Table 1).",
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
    # Hanan 2026 Section 2.3.3 final covariate model, EOT endpoint:
    #
    #   p_isa = logit^-1( drug.eff_d
    #                     + sum_m beta_m * theta_m,isa
    #                     + eta_isa )
    #
    # for study i, stratum s, treatment arm a. All values below are
    # from Hanan 2026 Table 2, "EOT model" block, "Estimate (logit)"
    # column. The paper's "Transformed estimate" column is the
    # back-transform reproduced in model() and asserted in the
    # vignette.
    # ================================================================

    logitp_mono <- -2.37
    label(paste(
      "Peg-IFNalpha monotherapy logit-scale EOT HBsAg-loss intercept",
      "in the reference population (HBeAg-negative, baseline HBsAg 3",
      "log10 IU/mL, 48 weeks of Peg-IFNalpha); back-transform",
      "logit^-1(-2.37) = 8.5%"
    ))
    # Table 2, EOT model, "Drug effect (Peg-IFNalpha)": estimate
    # (logit) -2.37, transformed estimate 8.54%, %RSE 7.7, 95% PI
    # 2.3-27.7, p < 0.001. Interpretation column: "8.5% HBsAg loss in
    # reference population (HBeAg-, baseline HBsAg 3 log10 IU/mL,
    # 48 weeks Peg-IFNalpha)."

    logitp_comb <- -2.02
    label(paste(
      "Peg-IFNalpha + nucleos(t)ide analogue logit-scale EOT",
      "HBsAg-loss intercept in the reference population;",
      "back-transform logit^-1(-2.02) = 11.7%"
    ))
    # Table 2, EOT model, "Drug effect (Peg-IFNalpha +NA)": estimate
    # (logit) -2.02, transformed estimate 11.7%, %RSE 10.3, 95% PI
    # 3.1-35.5, p < 0.001.

    e_t_pegifn_logitp <- 0.0133
    label(paste(
      "Peg-IFNalpha treatment-duration effect on the logit-scale EOT",
      "HBsAg-loss probability (per week above the 48-week reference);",
      "OR 1.013 per week"
    ))
    # Table 2, EOT model, "Peg-IFNalpha treatment duration (per week,
    # ref.=48 weeks)": estimate (logit) 0.0133, transformed estimate
    # OR 1.013, %RSE 32.9, 95% PI 1.00-1.02, p = 0.005.

    e_hbsag_bl_log10_logitp <- -1.57
    label(paste(
      "Baseline HBsAg effect on the logit-scale EOT HBsAg-loss",
      "probability (per log10 IU/mL above the 3 log10 IU/mL",
      "reference); OR 0.21 per log10 IU/mL"
    ))
    # Table 2, EOT model, "Baseline HBsAg (per log10 IU/mL, ref.=3)":
    # estimate (logit) -1.57, transformed estimate OR 0.21, %RSE 15.5,
    # 95% PI 0.13-0.33, p < 0.001. Results Section 3.3: "For each log10
    # unit increase in baseline HBsAg, the odds of achieving HBsAg loss
    # decreased by 79%" (1 - exp(-1.57) = 1 - 0.208 = 0.792).

    e_hbeag_pos_logitp <- 0.863
    label(paste(
      "Baseline HBeAg-positive effect on the logit-scale EOT",
      "HBsAg-loss probability (reference HBeAg-negative); OR 2.37"
    ))
    # Table 2, EOT model, "HBeAg status (ref=negative)": estimate
    # (logit) 0.863, transformed estimate OR 2.37, %RSE 34.7, 95% PI
    # 1.34-4.17, p = 0.008. (The Abstract and Results Section 3.3 quote
    # p = 0.007 for the same term; Table 2 prints 0.008. The
    # discrepancy is in the p-value only and does not affect any
    # estimate - see the vignette Errata.)

    # ================================================================
    # Between-study-strata-arm random effect, on the logit scale.
    # NOTE: this is BETWEEN-TRIAL variance, not popPK between-subject
    # variance. It describes how much a NEW STUDY-ARM's underlying
    # HBsAg-loss log-odds differs from the typical arm; it says nothing
    # about individual patients within an arm.
    # ================================================================
    eta_study ~ 0.484
    # Table 2, EOT model, "Between-trial variance": 0.484 (variance on
    # the logit scale). Reduced from the base model's 1.112 by the
    # treatment and covariate effects (Results Section 3.3: "from 1.112
    # to 0.484 in the EOT model, a 56.5% reduction"). The Abstract and
    # Discussion instead quote 58.1% for the same pair of variances;
    # (1.112 - 0.484) / 1.112 = 56.5%, so the Results figure is the
    # arithmetically correct one - see the vignette Errata.

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
    logit_p <- drug_eff +
               e_t_pegifn_logitp       * dur_dev   +
               e_hbsag_bl_log10_logitp * hbsag_dev +
               e_hbeag_pos_logitp      * HBEAG_POS +
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
