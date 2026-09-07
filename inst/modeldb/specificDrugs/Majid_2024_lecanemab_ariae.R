Majid_2024_lecanemab_ariae <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-safety model for the incidence",
    "of amyloid-related imaging abnormalities with edema/effusion (ARIA-E)",
    "in 2641 subjects with early Alzheimer's disease pooled from the",
    "lecanemab phase II Study 201 Core (n = 852) and phase III Clarity AD",
    "Study 301 Core (n = 1789); 177 subjects experienced ARIA-E. The",
    "probability of a subject experiencing ARIA-E at any time during the",
    "study is expit(-4.89 + 0.00666 * CMAX + 0.640 * APOE4_HET + 1.91 *",
    "APOE4_HOM), where CMAX is the individual model-predicted maximum",
    "serum lecanemab concentration at steady state in ug/mL, taken as an",
    "empirical-Bayes prediction from the companion population PK model",
    "Majid_2024_lecanemab. Exposure enters LINEARLY and UNTRANSFORMED, so",
    "the odds ratio is 1.95 per 100 ug/mL. APOE4 genotype is the only",
    "covariate retained from the univariate screen, entered as two binary",
    "indicators against the non-carrier reference (odds ratios 1.90",
    "heterozygous, 6.75 homozygous). This is a static per-subject landmark",
    "regression, not a time-to-event model: there is no ODE, no PK layer",
    "and no time dimension. NONMEM fits it by LAPLACE LIKE with the single",
    "eta fixed to zero, so no between-subject random effect and no",
    "residual error are estimated. The companion isolated-ARIA-H endpoint",
    "was analysed GRAPHICALLY ONLY and yielded no model -- the published",
    "finding is that isolated ARIA-H incidence is independent of lecanemab",
    "exposure and similar between placebo and treated subjects -- so it is",
    "deliberately absent here rather than omitted by transcription.",
    sep = " "
  )
  reference <- paste(
    "Majid O, Cao Y, Willis BA, Hayato S, Takenaka O, Lalovic B,",
    "Sreerama Reddy SH, Penner N, Reyderman L, Yasuda S, Hussein Z (2024).",
    "Population pharmacokinetics and exposure-response analyses of safety",
    "(ARIA-E and isolated ARIA-H) of lecanemab in subjects with early",
    "Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2111-2123.",
    "doi:10.1002/psp4.13224.",
    sep = " "
  )
  vignette <- "Majid_2024_lecanemab"

  units <- list(
    time          = "n/a (static per-subject landmark logistic regression; the outcome is 'ARIA-E occurred at any time during the study', so the model carries no time dimension and no dosing events)",
    dosing        = "n/a (no dose events; lecanemab exposure enters through the CMAX covariate, which an upstream population PK model supplies)",
    concentration = "prob_ariae (probability that a subject experiences ARIA-E during the study, 0-1; also logit_ariae, the untransformed linear predictor)"
  )

  covariateData <- list(
    CMAX = list(
      description        = "Individual model-predicted maximum serum lecanemab concentration at steady state (Css,max).",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "(a) TOTAL serum concentration, not unbound -- the IP/LC-MS/MS",
        "assay measures total lecanemab. (b) STEADY-STATE maximum over a",
        "dosing interval, not a single-dose or cycle-1 peak. (c) Derived as",
        "an individual EMPIRICAL BAYES prediction from the companion",
        "population PK model modellib('Majid_2024_lecanemab'): the Methods",
        "state that 'Individual empirical Bayes estimates of the PK",
        "parameters from the final model were used to derive Css,max',",
        "which is exactly the derivation the CMAX register entry requires",
        "to be documented per model. Enters LINEARLY and UNCENTRED, so the",
        "intercept is the logit at Css,max = 0 in an APOE4 non-carrier and",
        "is an extrapolated anchor rather than any real patient's risk.",
        "UNITS ARE LOAD-BEARING: the slope 0.00666 is per 1 ug/mL, and",
        "Table 2 gives the same effect a second time as an odds ratio of",
        "1.95 per ug/DECILITRE, i.e. per 100 ug/mL (exp(0.00666 * 100) =",
        "1.946). Feeding this model a Css,max in ng/mL would inflate the",
        "linear predictor a thousandfold. Study 301 model-predicted",
        "Css,max spans roughly 58-544 ug/mL with a mean of 305 ug/mL at 10",
        "mg/kg bi-weekly (Discussion); placebo subjects contribute",
        "Css,max = 0 and so sit at the intercept."
      ),
      source_name        = "CMAXSS (supplement Text S2 $INPUT and $PRED); Css,max (paper narrative and Table 2)"
    ),
    APOE4_HET = list(
      description        = "APOE-epsilon4 heterozygote indicator; 1 = exactly one epsilon4 allele, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (APOE-epsilon4 non-carrier; APOE4_HET and APOE4_HOM are mutually exclusive and a non-carrier has both set to 0)",
      notes              = paste(
        "Supplement Text S2 $PRED derives the pair as 'AP1=0; AP2=0; IF",
        "(APOEGEN.EQ.1) AP1=1; IF (APOEGEN.EQ.2) AP2=1', confirming the",
        "APOEGEN source alias named in the register entry and confirming",
        "that the omitted third level (non-carrier) is the reference. The",
        "two indicators are estimated separately rather than as an allele",
        "count because the effect is strongly NON-additive: the homozygote",
        "log-odds shift (1.91) is about three times the heterozygote shift",
        "(0.640), not twice it, so APOE4_COUNT could not represent this",
        "model. ARIA-E analysis set: 1423 heterozygous carriers of 2641",
        "(non-carrier 803, homozygous 415; Table S2)."
      ),
      source_name        = "APOEGEN == 1 (raw); AP1 (derived indicator in $PRED)"
    ),
    APOE4_HOM = list(
      description        = "APOE-epsilon4 homozygote indicator; 1 = two epsilon4 alleles, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (APOE-epsilon4 non-carrier)",
      notes              = paste(
        "Companion to APOE4_HET; see that entry for the derivation and the",
        "non-additivity argument. This is the largest effect in the model",
        "by a wide margin (odds ratio 6.75) and is the quantitative basis",
        "for the labelled recommendation that APOE-epsilon4 homozygotes",
        "carry the highest ARIA risk. ARIA-E analysis set: 415 homozygous",
        "carriers of 2641 (Table S2)."
      ),
      source_name        = "APOEGEN == 2 (raw); AP2 (derived indicator in $PRED)"
    )
  )

  covariatesDataExcluded <- list(
    CTROUGH = list(
      description = "Individual model-predicted minimum serum lecanemab concentration at steady state (Css,min).",
      units       = "ug/mL",
      type        = "continuous",
      notes       = paste(
        "Screened as an alternative exposure metric and NOT retained. This",
        "is a discrimination result, not a null result: the Results state",
        "that 'all exposure metrics were statistically significant and",
        "correlated with increasing incidence of ARIA-E' but that Css,max",
        "'was the statistically best predictor of ARIA-E in this analysis'.",
        "Css,max, Css,min and Css,av are mutually correlated, so only one",
        "could enter. No point estimate is published for a Css,min slope,",
        "and substituting Css,min into this model's 0.00666 slope would be",
        "a units-and-scale error rather than an approximation."
      )
    ),
    CAV = list(
      description = "Individual model-predicted average serum lecanemab concentration at steady state (Css,av).",
      units       = "ug/mL",
      type        = "continuous",
      notes       = paste(
        "Screened as an alternative exposure metric and not retained; see",
        "the CTROUGH entry for the reasoning. Note that the sibling model",
        "Cao_2026_lecanemab drives an Alzheimer's-disease QSP model from",
        "this same quantity under the canonical name CSS_LEC, so a user",
        "holding a Css,av column has a use for it -- just not in this",
        "model, whose slope was estimated against the peak."
      )
    ),
    APOE4_CARRIER = list(
      description = "APOE-epsilon4 carrier indicator, collapsing heterozygotes and homozygotes.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Carrier-vs-non-carrier status was screened alongside the",
        "three-level genotype and the genotype won: because the homozygote",
        "effect is roughly three times the heterozygote effect, collapsing",
        "the two carrier levels loses most of the signal. ARIA-E analysis",
        "set: 1838 carriers (70%) vs 803 non-carriers (30%) (Table S2)."
      )
    ),
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a univariate predictor of ARIA-E incidence and not retained. ARIA-E analysis set median 71.0 kg, range 29.2-130 (Table S2). Carried as column BWGT in supplement Text S2 $INPUT."
    ),
    AGE = list(
      description = "Age at study entry.",
      units       = "year",
      type        = "continuous",
      notes       = "Screened as a univariate predictor of ARIA-E incidence and not retained. ARIA-E analysis set median 72 years, range 50-90 (Table S2)."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a univariate predictor of ARIA-E incidence and not retained. ARIA-E analysis set: 1357 females (51.4%), 1284 males (48.6%) (Table S2)."
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status at the SUBJECT level.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a univariate predictor of ARIA-E incidence and not",
        "retained. Note the level differs from the companion population PK",
        "model, which uses a time-varying SAMPLE-level ADA status: here it",
        "is a single per-subject flag. ARIA-E analysis set: 423",
        "ADA-positive (16%), 2218 ADA-negative (84%) (Table S2)."
      )
    ),
    MMSE = list(
      description = "Baseline Mini-Mental State Examination total score.",
      units       = "(points, 0-30)",
      type        = "continuous",
      notes       = "Screened as a univariate predictor of ARIA-E incidence and not retained. Carried as column BMMSE in supplement Text S2 $INPUT."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2641L,
    n_studies      = 2L,
    n_observations = "2641 per-subject binary ARIA-E records (one per subject): 1789 from Study 301 Core and 852 from Study 201 Core. 177 subjects were ARIA-E positive, of whom 160 received lecanemab (129 at the approved 10 mg/kg bi-weekly regimen) and 17 received placebo",
    age_range      = "median 72 years, range 50-90 (Table S2)",
    weight_range   = "median 71.0 kg, range 29.2-130 (Table S2)",
    sex_female_pct = 51.4,
    race_ethnicity = c(
      White = 81.3, Japanese = 7.0, Korean = 4.9,
      `Black/African American` = 2.5,
      `Asian excluding Chinese/Japanese/Korean` = 1.0,
      Chinese = 0.5,
      `American Indian/Alaskan/Other/Missing` = 2.7
    ),
    disease_state  = "early Alzheimer's disease with confirmed amyloid pathology: mild cognitive impairment due to AD 1654 (62.6%), mild AD dementia 987 (37.4%) (Table S2)",
    dose_range     = "placebo 1142 subjects; lecanemab 1499 subjects -- bi-weekly 2.5 mg/kg (52), 5 mg/kg (92), 10 mg/kg (1053); monthly 5 mg/kg (51), 10 mg/kg (251) (Table S2)",
    regions        = "multicentre international; Study 201 Core and Study 301 Core (Clarity AD, NCT03887455)",
    apoe4_genotype = "non-carrier 803, heterozygous carrier 1423, homozygous carrier 415 (Table S2)",
    notes          = paste(
      "Precision was high for the core parameters (%RSE < 9.9%) and",
      "acceptable for the covariate effects (< 36%), and bootstrap medians",
      "agree with the point estimates throughout (Table 2). The model was",
      "evaluated by non-parametric bootstrap and by overlaying the",
      "model-predicted proportion experiencing ARIA-E on the observed",
      "proportion within each Css,max quartile, by APOE4 genotype (Figure",
      "3). Placebo subjects are included in the fit and enter at",
      "Css,max = 0, which is what identifies the intercept."
    )
  )

  ini({
    # ==================================================================
    # All values are FINAL estimates from Majid 2024 Table 2 ("Population
    # parameters and bootstrap CIs for the final logistic regression model
    # for incidence of ARIA-E"), whose header row prints the model form
    # itself:
    #   Logit = INT + SLP*Css,max + Cov_APOE4Hetero + Cov_APOE4Homo
    # Supplement Text S2 $PRED confirms it line for line:
    #   LOGIT = THETA(1) + THETA(2)*CMAXSS + THETA(3)*AP1 + THETA(4)*AP2
    #   A = EXP(LOGIT); PROB = A/(1+A)
    #
    # CLOSED-FORM CHECKS. Every published summary of this model is
    # reproduced by these four numbers alone:
    #   (1) Odds ratios, Table 2. exp(0.00666*100) = 1.946 -> printed 1.95
    #       per ug/dL; exp(0.640) = 1.896 -> printed 1.90; exp(1.91) =
    #       6.754 -> printed 6.75.
    #   (2) Discussion exposure-range odds ratios at the Study 301
    #       Css,max limits: exp(0.00666*58) = 1.471 -> printed 1.47, and
    #       exp(0.00666*544) = 37.45 -> printed 37.5.
    #   (3) Discussion incidence at the mean Css,max of 305 ug/mL, where
    #       the linear predictor is -4.89 + 0.00666*305 = -2.8587:
    #       non-carrier expit(-2.8587) = 5.42% -> printed 5.45%;
    #       heterozygous expit(-2.2187) = 9.81% -> printed 9.85%;
    #       homozygous expit(-0.9487) = 27.9% -> printed 28.0%.
    #   All three families agree to within rounding, which pins both the
    #   ug/mL unit of the slope and the non-carrier reference category.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -4.89; label("Logit of the ARIA-E probability for an APOE-epsilon4 non-carrier at a steady-state Cmax of 0 ug/mL (unitless logit)")  # Majid 2024 Table 2 "Intercept (INT)" -4.89, %RSE 5.38, bootstrap median -4.91 (95% CI -5.55 to -4.39)

    # ----- Exposure effect -----
    e_cmax_ariae <- 0.00666; label("Log-odds of ARIA-E per 1 ug/mL increase in steady-state maximum serum lecanemab concentration (unitless logit per ug/mL)")  # Majid 2024 Table 2 "Slope of lecanemab exposure effect (SLP; per Css,max unit in ug/mL)" 0.00666, %RSE 9.82, bootstrap median 0.00670 (95% CI 0.00540-0.00790); confirmed by the printed odds ratio 1.95 per ug/deciliter = exp(0.00666*100)

    # ----- APOE4 genotype effects (reference: non-carrier) -----
    e_apoe4_het_ariae <- 0.640; label("Log-odds of ARIA-E for APOE-epsilon4 heterozygotes vs the non-carrier reference (unitless logit)")  # Majid 2024 Table 2 "Cov APOE4 Hetero" 0.640, %RSE 35.8, bootstrap median 0.630 (95% CI 0.236-1.15); printed odds ratio 1.90 = exp(0.640)
    e_apoe4_hom_ariae <- 1.91 ; label("Log-odds of ARIA-E for APOE-epsilon4 homozygotes vs the non-carrier reference (unitless logit)")    # Majid 2024 Table 2 "Cov APOE4 Homo" 1.91, %RSE 12.7, bootstrap median 1.92 (95% CI 1.52-2.51); printed odds ratio 6.75 = exp(1.91)

    # ----- No between-subject variability, no residual error -----
    # Supplement Text S2 fits this with "$OMEGA 0 FIX" under
    # "$ESTIMATION ... LAPLACE LIKE": the single eta is fixed at zero and
    # the likelihood is Bernoulli, so there is no omega and no sigma to
    # transcribe. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity. Same convention as the
    # Chen_2021_lorlatinib_* exposure-response family.
    addSd_prob_ariae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value ARIA-E probability; the source likelihood is Bernoulli and reports no residual error")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor (Majid 2024 Table 2 header row) ----------
    # CMAX must be supplied in ug/mL; see covariateData[[CMAX]]$notes.
    # Set CMAX = 0 to recover the placebo / no-exposure baseline risk.
    logit_ariae <- logit_ref +
      e_cmax_ariae      * CMAX +
      e_apoe4_het_ariae * APOE4_HET +
      e_apoe4_hom_ariae * APOE4_HOM

    prob_ariae <- expit(logit_ariae)

    # ----- Observation ----------------------------------------------
    prob_ariae ~ add(addSd_prob_ariae)
  })
}
