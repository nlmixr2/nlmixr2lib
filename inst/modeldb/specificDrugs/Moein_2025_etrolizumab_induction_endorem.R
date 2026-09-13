Moein_2025_etrolizumab_induction_endorem <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for ENDOSCOPIC ",
    "REMISSION at the END OF INDUCTION (week 14) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 384, the induction exposure-response ",
    "analysis set of the phase 3 BERGAMOT study, NCT02394028). ",
    "Endoscopic remission is an SES-CD <= 4 (<= 2 for ileal patients) ",
    "with no segment scoring above 1. This is the most heavily ",
    "covariate-adjusted of the six models in the family, retaining ",
    "five covariates: TNF-naive status, both disease-location ",
    "indicators, baseline SES-CD and baseline CDAI. The etrolizumab ",
    "exposure slope is not significant (0.0495 log-odds per ug/mL, RSE ",
    "95.8%, P = 0.297). There is no PK layer and no ODE -- the ",
    "exposure metric arrives as the CTROUGH data column, the ",
    "individual predicted week-4 trough after a SINGLE 105 or 210 mg ",
    "dose from modellib('Moein_2025_etrolizumab'). One of six landmark ",
    "models in the Moein_2025_etrolizumab_* family."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S8, 'Endoscopic Remission /",
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
    concentration = "prob_endorem (probability of endoscopic remission at end of induction, 0-1)"
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
        "placebo arm (a single intercept is shared by placebo and ",
        "active arms). Deliberately a SINGLE-DOSE prediction so the ",
        "metric is not confounded by the outcome-correlated decline in ",
        "clearance. For INDUCTION the dose is the patient's assigned ",
        "105 or 210 mg; the maintenance companions always use 105 mg, ",
        "so the columns are NOT interchangeable. Obtain it by solving a ",
        "single dose to day 28 with ",
        "modellib('Moein_2025_etrolizumab'). Placebo subjects take ",
        "CTROUGH = 0. Observed distribution (Table S6): mean 3.85, ",
        "median 3.23, range 0-16.6 ug/mL. The Discussion notes a ",
        "positive exposure trend for endoscopic remission at BOTH ",
        "induction and maintenance, even though only the maintenance ",
        "slope reaches significance."
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
        "contrast (Table S8 'TNF-naive' = 0.785, RSE 48.4%, P < 0.05), ",
        "with TNF-experienced absorbed into the intercept (footnote b). ",
        "The canonical PRIOR_TNF is 1 for TNF-experienced, so the naive ",
        "indicator is formed in model() as (1 - PRIOR_TNF) with the ",
        "published coefficient and intercept UNCHANGED -- the ",
        "register's SMOKE_NEVER idiom, preserving the published ",
        "parameter covariance. Analysis set (Table S7): 205 ",
        "TNF-experienced (53%), 179 TNF-naive (47%)."
      ),
      source_name        = "TNF status"
    ),
    DISLOC_ILEUM = list(
      description        = "Crohn's disease located in the ileum only (Montreal L1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3, the paper's reference; implies DISLOC_COLON = 0 as well)",
      notes              = paste0(
        "Log-odds shift versus ileocolonic disease (-1.03 here). Paired ",
        "with DISLOC_COLON to encode the three-level Montreal ",
        "disease-location categorical, with both indicators zero ",
        "meaning ileum and colon (L3) -- the reference named in ",
        "Table S8 footnote b. Mutually exclusive with DISLOC_COLON. ",
        "This endpoint is the only induction model that retains BOTH ",
        "location indicators, and they point in opposite directions: ",
        "ileal disease lowers the probability while colon-only raises ",
        "it. Analysis set (Table S7): 230 ileum and colon (60%), 72 ",
        "ileum only (19%), 82 colon only (21%)."
      ),
      source_name        = "Disease location"
    ),
    DISLOC_COLON = list(
      description        = "Crohn's disease located in the colon only (Montreal L2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3)",
      notes              = paste0(
        "Log-odds shift versus ileocolonic disease (+0.881 here). ",
        "Mutually exclusive with DISLOC_ILEUM; both zero means ",
        "ileocolonic. Colon-only disease predicts HIGHER endoscopic ",
        "remission than ileocolonic, the mirror image of the ileum-only ",
        "penalty, consistent with the Discussion's observation of lower ",
        "efficacy in ileal Crohn's disease."
      ),
      source_name        = "Disease location"
    ),
    SCORE_SESCD = list(
      description        = "Baseline Simple Endoscopic Score for Crohn's Disease (SES-CD)",
      units              = "(score, 0-56)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Linear effect on the logit, CENTERED at the reference score of ",
        "12.0: e_sescd_endorem * (SCORE_SESCD - 12.0). The centering is ",
        "load-bearing and is what makes the printed intercept a ",
        "probability AT the reference covariates -- Table S8 footnote b ",
        "defines the intercept patient as having 'SES-CD score: 12.0', ",
        "and the paper's Figure 4c shows the percentile-stratified ",
        "curves passing through the reference intercept at the median. ",
        "A higher baseline endoscopic burden lowers the probability of ",
        "reaching endoscopic remission, which is mechanically sensible ",
        "given the endpoint's absolute SES-CD threshold. Analysis-set ",
        "distribution (Table S6): median 12.0, range 4.00-43.0."
      ),
      source_name        = "SES-CD score"
    ),
    SCORE_CDAI = list(
      description        = "Baseline Crohn's Disease Activity Index (CDAI) score",
      units              = "(score, 0-600)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Linear effect on the logit, CENTERED at the INDUCTION ",
        "reference score of 322: e_cdai_endorem * (SCORE_CDAI - 322). ",
        "Note that the reference differs between phases -- Table S8 ",
        "(induction) uses 322 while Table S9 (maintenance) uses 320 -- ",
        "so the maintenance companion models center at 320. Higher ",
        "baseline disease activity lowers the probability of a ",
        "favourable outcome. Analysis-set distribution (Table S6): ",
        "median 322, range 215-481."
      ),
      source_name        = "CDAI score"
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
      "ileum and colon, non-smoker, male, TNF-experienced."
    )
  )

  ini({
    # ==================================================================
    # Moein 2025 Supplementary Table S8, "Endoscopic Remission ->
    # Covariate-adjusted final model" column.
    #
    # The intercept is printed on the PROBABILITY scale (Table S8 Units
    # line: "Intercept (Probability)"; footnote b: "reflects the
    # probability ... of a patient treated with placebo" at the listed
    # reference covariates) and is therefore logit-transformed here. All
    # other coefficients are already log-odds.
    #
    # Both continuous covariates are CENTERED at their footnote-b
    # reference values, which is what makes that statement true.
    # ==================================================================
    logite0 <- logit(0.0577); label("Logit of the probability of endoscopic remission at end of induction for the reference placebo patient (unitless logit)")  # Table S8: Intercept = 0.0577 (RSE 14.0%, P < 0.01), on the PROBABILITY scale per the Units line and footnote b

    # ----- Exposure slope (log-odds per ug/mL) -----
    e_ctrough_endorem <- 0.0495; label("Log-odds of endoscopic remission at end of induction per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S8 'ER effect' = 0.0495 (RSE 95.8%); main-text Table 2 gives the same value with P = 0.297

    # ----- Covariate effects on the logit -----
    e_tnfnaive_endorem <- 0.785;    label("Log-odds shift for TNF-naive vs TNF-experienced patients (unitless logit)")             # Table S8 'TNF-naive' = 0.785 (RSE 48.4%, P < 0.05)
    e_ileum_endorem    <- -1.03;    label("Log-odds shift for ileum-only vs ileocolonic disease location (unitless logit)")        # Table S8 'Disease location Ileum only' = -1.03 (RSE 50.8%, P < 0.05)
    e_colon_endorem    <- 0.881;    label("Log-odds shift for colon-only vs ileocolonic disease location (unitless logit)")        # Table S8 'Disease location Colon only' = 0.881 (RSE 46.9%, P < 0.05)
    e_sescd_endorem    <- -0.113;   label("Log-odds shift per SES-CD point above the reference score of 12.0 (unitless logit)")    # Table S8 'SES-CD effect' = -0.113 (RSE 34.6%, P < 0.05); units 'LO per score'
    e_cdai_endorem     <- -0.00578; label("Log-odds shift per CDAI point above the induction reference score of 322 (unitless logit)")  # Table S8 'CDAI effect' = -0.00578 (RSE 51.8%); units 'LO per score'

    # ----- No between-subject variability, no residual error -----
    addSd_prob_endorem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the log-odds scale. Continuous covariates are
    # centered at their reference values and the categorical ones are
    # differences from the reference category, so every covariate term
    # vanishes for a reference patient and logite0 alone remains.
    logit_endorem <-
      logite0 +
      e_ctrough_endorem * CTROUGH +
      e_tnfnaive_endorem * (1 - PRIOR_TNF) +
      e_ileum_endorem * DISLOC_ILEUM +
      e_colon_endorem * DISLOC_COLON +
      e_sescd_endorem * (SCORE_SESCD - 12.0) +
      e_cdai_endorem * (SCORE_CDAI - 322)

    prob_endorem <- expit(logit_endorem)

    prob_endorem ~ add(addSd_prob_endorem)
  })
}
