Moein_2025_etrolizumab_maintenance_endorem <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for ENDOSCOPIC ",
    "REMISSION at the END OF MAINTENANCE (week 66) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 434, the maintenance ",
    "exposure-response analysis set of the phase 3 BERGAMOT study, ",
    "NCT02394028, placebo or 105 mg SC Q4W). Endoscopic remission is ",
    "an SES-CD <= 4 (<= 2 for ileal patients) with no segment scoring ",
    "above 1. This is the steepest exposure slope of the six models ",
    "(0.257 log-odds per ug/mL, RSE 29%, P = 0.000555) and it acts on ",
    "the lowest baseline probability (4.18% for the reference placebo ",
    "patient), so it produces the largest relative exposure effect in ",
    "the paper. Two covariates are retained: baseline SES-CD and ",
    "colon-only disease location. There is no PK layer and no ODE -- ",
    "the exposure metric arrives as the CTROUGH data column, the ",
    "individual predicted week-4 trough after a SINGLE 105 mg dose ",
    "from modellib('Moein_2025_etrolizumab'). Reproduces panels (c) ",
    "and (h) of Moein 2025 Figure 4. One of six landmark models in the ",
    "Moein_2025_etrolizumab_* family."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S9, 'Endoscopic Remission /",
    "Covariate-adjusted final model' column; the exposure slope is also",
    "printed in main-text Table 2. The reference-patient covariate",
    "values are Table S9 footnote b. Graphical check: Figure 4c and 4h.",
    "Exposure metric supplied by modellib('Moein_2025_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(
    time          = "n/a (static landmark logistic regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_endorem (probability of endoscopic remission at end of maintenance, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual predicted etrolizumab serum trough concentration at week 4 following a SINGLE 105 mg dose (Ctrough,W4,adjusted)",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "TOTAL etrolizumab, entering the logit LINEARLY and ",
        "UNCENTERED, so logite0 is the logit at zero exposure -- the ",
        "placebo arm (a single intercept is shared by placebo and ",
        "active arms). For MAINTENANCE the metric is ALWAYS predicted ",
        "from a single 105 mg dose, even for patients who received ",
        "210 mg during induction; the induction companion models use ",
        "the patient's assigned 105 or 210 mg, so the columns are NOT ",
        "interchangeable. The single-dose basis removes the ",
        "confounding an on-treatment trough would carry, because ",
        "etrolizumab clearance declines over time in a way that ",
        "correlates with clinical improvement. Obtain it by solving a ",
        "single 105 mg dose to day 28 with ",
        "modellib('Moein_2025_etrolizumab'). Placebo subjects take ",
        "CTROUGH = 0. Observed distribution (Table S6, maintenance): ",
        "mean 1.66, median 0.164, range 0-9.14 ug/mL."
      ),
      source_name        = "Ctrough,W4,adjusted"
    ),
    DISLOC_ILEUM = list(
      description        = "Crohn's disease located in the ileum only (Montreal L1)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3)",
      notes              = paste0(
        "Declared for completeness of the three-level disease-location ",
        "encoding and referenced in model() with a coefficient of ",
        "exactly zero, because Table S9 retained no ileum-only effect ",
        "for endoscopic remission at maintenance. The Figure 4h caption ",
        "states it directly: 'Ileum only was not significantly ",
        "different from ileum and colon, and therefore the two are ",
        "predicted to be identical', with the ileocolonic curve hidden ",
        "beneath the ileum-only one. Note this is the MIRROR of the ",
        "endoscopic-improvement model, where ileum-only was retained ",
        "and colon-only was the structural zero. Carrying both terms ",
        "explicitly keeps the data-column contract identical across all ",
        "six models of the family. Mutually exclusive with ",
        "DISLOC_COLON. Analysis set (Table S7): 267 ileum and colon ",
        "(62%), 78 ileum only (18%), 89 colon only (21%)."
      ),
      source_name        = "Disease location"
    ),
    DISLOC_COLON = list(
      description        = "Crohn's disease located in the colon only (Montreal L2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ileocolonic disease, Montreal L3, the paper's reference; implies DISLOC_ILEUM = 0 as well)",
      notes              = paste0(
        "Log-odds shift versus ileocolonic disease (+1.03 here) -- the ",
        "largest retained covariate effect in the maintenance models. ",
        "Colon-only disease predicts markedly HIGHER endoscopic ",
        "remission, consistent with the Discussion's observation of ",
        "lower efficacy in ileal Crohn's disease. Paired with ",
        "DISLOC_ILEUM to encode the three-level Montreal ",
        "disease-location categorical, with both indicators zero ",
        "meaning ileum and colon (L3) -- the reference named in ",
        "Table S9 footnote b. Reproduced in Figure 4h."
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
        "12.0: e_sescd_endorem * (SCORE_SESCD - 12.0). Table S9 ",
        "footnote b defines the intercept patient as having 'SES-CD ",
        "score: 12.0', and the centering is what makes the printed ",
        "intercept a probability at the reference. A higher baseline ",
        "endoscopic burden lowers the probability of reaching ",
        "endoscopic remission (-0.0765 log-odds per point), which is ",
        "mechanically sensible given the endpoint's absolute SES-CD ",
        "threshold. Note the induction counterpart estimates a steeper ",
        "-0.113 per point on the same covariate. Reproduced in ",
        "Figure 4c. Analysis-set distribution (Table S6): median 12.0, ",
        "range 3.00-38.0."
      ),
      source_name        = "SES-CD score"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 434L,
    n_studies      = 1L,
    n_observations = "434 binary outcome records, one per patient (landmark analysis at end of maintenance)",
    age_range      = "18.0-76.0 years",
    age_median     = "37.0 years",
    weight_range   = "35.3-154 kg",
    weight_median  = "70.4 kg",
    sex_female_pct = 50,
    disease_state  = "Moderately-to-severely active Crohn's disease with a CDAI-70 response at the end of induction; baseline CDAI median 320, SES-CD median 12.0",
    dose_range     = "Placebo or etrolizumab 105 mg SC Q4W over a 52-week maintenance phase, following a 14-week induction phase on 105 or 210 mg SC",
    regions        = "Multinational (BERGAMOT, NCT02394028)",
    notes          = paste0(
      "Same maintenance analysis set as the two companion maintenance ",
      "models. Entry required a CDAI-70 response at the end of ",
      "induction, so the set is RESPONDER-ENRICHED; induction-placebo ",
      "patients were not randomised into maintenance and are excluded. ",
      "Baseline characteristics from Moein 2025 Tables S6 and S7, ",
      "maintenance column. Table S9 footnote b reference patient: CDAI ",
      "320, SES-CD 12.0, fecal calprotectin 836 ug/g, CRP 8.35 mg/L, ",
      "albumin 42.0 g/L, white blood cells 7.88 x10^9/L, neutrophils ",
      "5.41 x10^9/L, MAdCAM-1 17.8 U, disease location ileum and colon, ",
      "non-smoker, male, TNF-experienced."
    )
  )

  ini({
    # ==================================================================
    # Moein 2025 Supplementary Table S9, "Endoscopic Remission ->
    # Covariate-adjusted final model" column.
    #
    # The intercept is printed on the PROBABILITY scale (Table S9 Units
    # line: "Intercept (Probability)"; footnote b) and is
    # logit-transformed here; all other coefficients are already
    # log-odds.
    #
    # Figure 4h confirms the reading at the placebo end of the x-axis:
    # the panel reads about 0.04 for the ileocolonic / ileum-only
    # reference and about 0.10 for colon-only, against
    # expit(logit(0.0418)) = 0.0418 and
    # expit(logit(0.0418) + 1.03) = 0.109 from this parameterisation.
    # ==================================================================
    logite0 <- logit(0.0418); label("Logit of the probability of endoscopic remission at end of maintenance for the reference placebo patient (unitless logit)")  # Table S9: Intercept = 0.0418 (RSE 9.53%, P < 0.01), on the PROBABILITY scale per the Units line and footnote b

    # ----- Exposure slope (log-odds per ug/mL) -----
    # The steepest slope of the six models, acting on the lowest
    # baseline probability.
    e_ctrough_endorem <- 0.257; label("Log-odds of endoscopic remission at end of maintenance per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S9 'ER effect' = 0.257 (RSE 29.0%, P < 0.01); main-text Table 2 gives the same value with P = 0.000555

    # ----- Covariate effects on the logit -----
    e_colon_endorem <- 1.03;    label("Log-odds shift for colon-only vs ileocolonic disease location (unitless logit)")             # Table S9 'Disease location Colon only' = 1.03 (RSE 39.7%, P < 0.05)
    e_sescd_endorem <- -0.0765; label("Log-odds shift per SES-CD point above the reference score of 12.0 (unitless logit)")         # Table S9 'SES-CD effect' = -0.0765 (RSE 43.3%, P < 0.05); units 'LO per score'

    # Ileum-only was NOT retained for this endpoint (Table S9 leaves the
    # cell blank; Figure 4h states ileum-only is predicted identically to
    # the ileocolonic reference). Held at exactly zero, and marked fixed
    # so it reads unmistakably as a structural zero rather than an
    # estimate, to keep the data-column contract uniform across the
    # six-model family.
    e_ileum_endorem <- fixed(0); label("Log-odds shift for ileum-only vs ileocolonic disease location, not retained by the source model and held at a structural zero (unitless logit)")  # Table S9: no ileum-only row for endoscopic remission; Figure 4h caption confirms ileum only is predicted identically to ileum and colon

    # ----- No between-subject variability, no residual error -----
    addSd_prob_endorem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the log-odds scale. The continuous covariate
    # is centered at its reference value and the categorical ones are
    # differences from the reference category, so every covariate term
    # vanishes for a reference patient and logite0 alone remains.
    logit_endorem <-
      logite0 +
      e_ctrough_endorem * CTROUGH +
      e_ileum_endorem * DISLOC_ILEUM +
      e_colon_endorem * DISLOC_COLON +
      e_sescd_endorem * (SCORE_SESCD - 12.0)

    prob_endorem <- expit(logit_endorem)

    prob_endorem ~ add(addSd_prob_endorem)
  })
}
