Moein_2025_etrolizumab_induction_endoimp <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for ENDOSCOPIC ",
    "IMPROVEMENT at the END OF INDUCTION (week 14) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 384, the induction exposure-response ",
    "analysis set of the phase 3 BERGAMOT study, NCT02394028). ",
    "Endoscopic improvement is a >= 50% reduction from the baseline ",
    "Simple Endoscopic Score for Crohn's Disease (SES-CD). The ",
    "etrolizumab exposure slope is the largest of the three induction ",
    "endpoints but still does not reach significance (0.0682 log-odds ",
    "per ug/mL, RSE 53.4%, P = 0.061); TNF-naive status and ",
    "ileum-only disease location are the retained covariates and both ",
    "are significant. There is no PK layer and no ODE -- the exposure ",
    "metric arrives as the CTROUGH data column, the individual ",
    "predicted week-4 trough after a SINGLE 105 or 210 mg dose from ",
    "the companion model modellib('Moein_2025_etrolizumab'). One of ",
    "six landmark models in the Moein_2025_etrolizumab_* family."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S8, 'Endoscopic Improvement /",
    "Covariate-adjusted final model' column; the exposure slope is also",
    "printed in main-text Table 2. The reference-patient covariate",
    "values are Table S8 footnote b.",
    "Exposure metric supplied by modellib('Moein_2025_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(
    time          = "n/a (static landmark logistic regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_endoimp (probability of endoscopic improvement at end of induction, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual predicted etrolizumab serum trough concentration at week 4 following a SINGLE dose (Ctrough,W4,adjusted)",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "TOTAL etrolizumab, entering the logit LINEARLY and ",
        "UNCENTERED, so logite0 is the logit at zero exposure -- the ",
        "placebo arm (the paper fits a single intercept shared by ",
        "placebo and active arms). Deliberately a SINGLE-DOSE ",
        "prediction so that the exposure metric is not confounded by ",
        "the outcome-correlated decline in clearance. For INDUCTION the ",
        "dose is the patient's assigned 105 or 210 mg; the maintenance ",
        "companion models always use 105 mg, so the two columns are NOT ",
        "interchangeable. Obtain it by solving a single dose to day 28 ",
        "with modellib('Moein_2025_etrolizumab'). Placebo subjects take ",
        "CTROUGH = 0. Observed distribution (Table S6): mean 3.85, ",
        "median 3.23, range 0-16.6 ug/mL."
      ),
      source_name        = "Ctrough,W4,adjusted"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF biologic therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (TNF-experienced) -- the PAPER's reference, the complement of the canonical column's 0 level",
      notes              = paste0(
        "The published coefficient is on the TNF-NAIVE side of the ",
        "contrast (Table S8 'TNF-naive' = 1.14, RSE 23.0%, P < 0.01), ",
        "with TNF-experienced absorbed into the intercept (footnote b). ",
        "The canonical PRIOR_TNF is 1 for TNF-experienced, so the naive ",
        "indicator is formed in model() as (1 - PRIOR_TNF) with the ",
        "published coefficient and intercept carried UNCHANGED -- the ",
        "register's documented SMOKE_NEVER idiom, which preserves the ",
        "published parameter covariance. This is the largest covariate ",
        "effect anywhere in the induction models: a TNF-naive patient's ",
        "odds of endoscopic improvement are exp(1.14) = 3.1-fold those ",
        "of an otherwise identical TNF-experienced patient. Analysis ",
        "set (Table S7): 205 experienced (53%), 179 naive (47%)."
      ),
      source_name        = "TNF status"
    ),
    DISLOC_ILEUM = list(
      description        = "Crohn's disease located in the ileum only (Montreal L1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3, the most common group and the paper's reference; implies DISLOC_COLON = 0 as well)",
      notes              = paste0(
        "Log-odds shift versus ileocolonic disease. Paired with ",
        "DISLOC_COLON to encode the three-level Montreal disease-location ",
        "categorical: DISLOC_ILEUM = 1 is ileum only (L1), ",
        "DISLOC_COLON = 1 is colon only (L2), and both zero is ileum and ",
        "colon (L3), which Table S8 footnote b names as the reference. ",
        "The two indicators are mutually exclusive. Ileal disease ",
        "predicts markedly LOWER efficacy (-1.17 log-odds here), which ",
        "the Discussion notes is consistent with published findings for ",
        "certolizumab pegol and infliximab. Colon-only was not retained ",
        "for this endpoint, so it is predicted identically to the ",
        "ileocolonic reference (stated explicitly in the Figure 4g ",
        "caption). Analysis set (Table S7): 230 ileum and colon (60%), ",
        "72 ileum only (19%), 82 colon only (21%)."
      ),
      source_name        = "Disease location"
    ),
    DISLOC_COLON = list(
      description        = "Crohn's disease located in the colon only (Montreal L2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3)",
      notes              = paste0(
        "Declared for completeness of the three-level disease-location ",
        "encoding and referenced in model() with a coefficient of ",
        "exactly zero, because Table S8 retained NO colon-only effect ",
        "for endoscopic improvement at induction: the Figure 4g caption ",
        "states 'Colon only was not significantly different from ileum ",
        "and colon, and therefore the two are predicted to be ",
        "identical'. Carrying the term explicitly with a zero ",
        "coefficient keeps the data-column contract identical across ",
        "all six models of the family, so one cohort data frame drives ",
        "every model. Mutually exclusive with DISLOC_ILEUM."
      ),
      source_name        = "Disease location"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 384L,
    n_studies      = 1L,
    n_observations = "384 binary outcome records, one per patient (landmark analysis at end of induction)",
    age_range      = "18.0-79.0 years",
    age_median     = "35.0 years",
    weight_range   = "40.4-160 kg",
    weight_median  = "71.8 kg",
    sex_female_pct = 46,
    disease_state  = "Moderately-to-severely active Crohn's disease; baseline CDAI median 322, SES-CD median 12.0",
    dose_range     = "Placebo, etrolizumab 105 mg SC Q4W, or etrolizumab 210 mg SC Q4W with an additional 210 mg loading dose at week 2; 14-week induction phase",
    regions        = "Multinational (BERGAMOT, NCT02394028)",
    notes          = paste0(
      "Same induction analysis set as the two companion induction ",
      "models. Baseline characteristics from Moein 2025 Tables S6 and ",
      "S7, induction column. Table S8 footnote b reference patient: ",
      "CDAI 322, SES-CD 12.0, fecal calprotectin 805 ug/g, CRP ",
      "8.34 mg/L, albumin 42.0 g/L, white blood cells 8.14 x10^9/L, ",
      "neutrophils 5.65 x10^9/L, MAdCAM-1 17.7 U, disease location ",
      "ileum and colon, non-smoker, male, TNF-experienced. Only the ",
      "TNF-status and disease-location entries of that list are ",
      "retained covariates for this endpoint; the rest were screened by ",
      "the full covariate model and dropped."
    )
  )

  ini({
    # ==================================================================
    # Moein 2025 Supplementary Table S8, "Endoscopic Improvement ->
    # Covariate-adjusted final model" column.
    #
    # The intercept is printed on the PROBABILITY scale -- Table S8's
    # Units line reads "Intercept (Probability)" and footnote b says it
    # "reflects the probability for the outcome of a patient treated
    # with placebo" at the listed reference covariates -- so it is
    # logit-transformed here. All other coefficients are already
    # log-odds ("LO") per the same Units line.
    # ==================================================================
    logite0 <- logit(0.144); label("Logit of the probability of endoscopic improvement at end of induction for the reference placebo patient (unitless logit)")  # Table S8: Intercept = 0.144 (RSE 13.5%, P < 0.01), on the PROBABILITY scale per the Units line and footnote b

    # ----- Exposure slope (log-odds per ug/mL) -----
    e_ctrough_endoimp <- 0.0682; label("Log-odds of endoscopic improvement at end of induction per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S8 'ER effect' = 0.0682 (RSE 53.4%); main-text Table 2 gives the same value with P = 0.061

    # ----- Covariate effects on the logit -----
    e_tnfnaive_endoimp <- 1.14;  label("Log-odds shift for TNF-naive vs TNF-experienced patients (unitless logit)")            # Table S8 'TNF-naive' = 1.14 (RSE 23.0%, P < 0.01)
    e_ileum_endoimp    <- -1.17; label("Log-odds shift for ileum-only vs ileocolonic disease location (unitless logit)")       # Table S8 'Disease location Ileum only' = -1.17 (RSE 33.3%, P < 0.05)

    # Colon-only was NOT retained for this endpoint; Table S8 leaves the
    # cell blank and Figure 4g states colon-only is predicted identically
    # to the ileocolonic reference. Held at exactly zero (fixed, so it is
    # unmistakably a structural zero rather than an estimate) to keep the
    # data-column contract uniform across the six-model family.
    e_colon_endoimp <- fixed(0); label("Log-odds shift for colon-only vs ileocolonic disease location, not retained by the source model and held at a structural zero (unitless logit)")  # Table S8: no colon-only row for endoscopic improvement; Figure 4g caption confirms colon only is predicted identically to ileum and colon

    # ----- No between-subject variability, no residual error -----
    # Bernoulli likelihood; the placeholder SD exists only because
    # rxode2 requires an observation declaration. Not a source value.
    addSd_prob_endoimp <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the log-odds scale. Covariate terms are
    # differences from the reference category (TNF-experienced,
    # ileocolonic), so they vanish for a reference patient.
    logit_endoimp <-
      logite0 +
      e_ctrough_endoimp * CTROUGH +
      e_tnfnaive_endoimp * (1 - PRIOR_TNF) +
      e_ileum_endoimp * DISLOC_ILEUM +
      e_colon_endoimp * DISLOC_COLON

    prob_endoimp <- expit(logit_endoimp)

    prob_endoimp ~ add(addSd_prob_endoimp)
  })
}
