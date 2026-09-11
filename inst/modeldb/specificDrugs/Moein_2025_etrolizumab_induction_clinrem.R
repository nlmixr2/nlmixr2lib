Moein_2025_etrolizumab_induction_clinrem <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for CLINICAL ",
    "REMISSION at the END OF INDUCTION (week 14) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 384, the induction exposure-response ",
    "analysis set of the phase 3 BERGAMOT study, NCT02394028, ",
    "placebo / 105 mg SC Q4W / 210 mg SC Q4W with a week-2 loading ",
    "dose). Clinical remission is a CDAI score < 150 with stool ",
    "frequency mean daily score <= 3 and abdominal pain mean daily ",
    "score <= 1 and no worsening in either subscore. The etrolizumab ",
    "exposure slope is NOT statistically significant for this endpoint ",
    "(0.0259 log-odds per ug/mL, RSE 126%, P = 0.426): at induction ",
    "prior anti-TNF status, not exposure, drives the outcome. Only ",
    "TNF-naive status was retained as a covariate. There is no PK ",
    "layer and no ODE -- the exposure metric arrives as the CTROUGH ",
    "data column, which in the source is the individual predicted ",
    "week-4 trough after a SINGLE 105 or 210 mg dose from the ",
    "companion population PK model modellib('Moein_2025_etrolizumab'). ",
    "One of six landmark models in the Moein_2025_etrolizumab_* family ",
    "(three endpoints x induction/maintenance)."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S8, 'Clinical Remission /",
    "Covariate-adjusted final model' column; the exposure slope is also",
    "printed in main-text Table 2. The reference-patient covariate",
    "values are Table S8 footnote a.",
    "Exposure metric supplied by modellib('Moein_2025_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(
    time          = "n/a (static landmark logistic regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_clinrem (probability of clinical remission at end of induction, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual predicted etrolizumab serum trough concentration at week 4 following a SINGLE dose (Ctrough,W4,adjusted)",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "TOTAL (not unbound) etrolizumab, entering the logit LINEARLY ",
        "and UNCENTERED, so logite0 is the logit at zero exposure -- ",
        "which is exactly the placebo arm. Moein 2025 found that 'linear ",
        "models with a single intercept for placebo and active ",
        "treatments provided an adequate fit', so no separate placebo ",
        "intercept exists. The metric is deliberately a SINGLE-DOSE ",
        "prediction: because etrolizumab clearance falls with time in a ",
        "way that correlates with clinical improvement, an on-treatment ",
        "trough would be confounded by outcome, so the authors ",
        "back-predicted the week-4 trough that each patient would have ",
        "had from one dose alone. For INDUCTION the dose used is the ",
        "patient's actual assigned dose, 105 or 210 mg (for MAINTENANCE ",
        "it is always 105 mg, so the induction and maintenance columns ",
        "are NOT interchangeable). Obtain it from ",
        "modellib('Moein_2025_etrolizumab') by solving a single dose to ",
        "day 28; that model's clearance is constant over a single dose, ",
        "so the time-dependent term is inactive by construction. ",
        "Placebo subjects take CTROUGH = 0. Units are load-bearing: the ",
        "intercept absorbs the unit choice. Observed distribution in ",
        "this analysis set (Table S6): mean 3.85, median 3.23, range ",
        "0-16.6 ug/mL."
      ),
      source_name        = "Ctrough,W4,adjusted"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF biologic therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (TNF-experienced) -- see notes; this is the PAPER's reference category, which is the complement of the canonical column's 0 level",
      notes              = paste0(
        "The published coefficient is stated on the TNF-NAIVE side of ",
        "the contrast (Table S8 row 'TNF-naive' = 0.580, RSE 39.2%, ",
        "P < 0.05), with TNF-EXPERIENCED as the reference absorbed into ",
        "the intercept (Table S8 footnote a: 'TNF status: ",
        "TNF-experienced'). The canonical column PRIOR_TNF is 1 for ",
        "TNF-experienced, so the naive indicator is formed inside ",
        "model() as (1 - PRIOR_TNF) and the published coefficient and ",
        "intercept are carried UNCHANGED. This follows the register's ",
        "documented idiom for SMOKE_NEVER, and it preserves the ",
        "correspondence with the published parameter covariance. Do NOT ",
        "instead flip the sign onto PRIOR_TNF and shift the intercept. ",
        "TNF-naive patients had a significantly higher probability of a ",
        "favourable outcome at end of induction for ALL THREE endpoints ",
        "-- the paper's headline induction finding. Analysis set split ",
        "(Table S7): 205 TNF-experienced (53%), 179 TNF-naive (47%)."
      ),
      source_name        = "TNF status"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 384L,
    n_studies      = 1L,
    n_observations = "384 binary outcome records, one per patient (landmark analysis at end of induction; no repeated measures)",
    age_range      = "18.0-79.0 years",
    age_median     = "35.0 years",
    weight_range   = "40.4-160 kg",
    weight_median  = "71.8 kg",
    sex_female_pct = 46,
    disease_state  = "Moderately-to-severely active Crohn's disease; baseline CDAI median 322 (range 215-481), SES-CD median 12.0 (4.00-43.0)",
    dose_range     = "Placebo, etrolizumab 105 mg SC Q4W, or etrolizumab 210 mg SC Q4W with an additional 210 mg loading dose at week 2; 14-week induction phase",
    regions        = "Multinational (BERGAMOT, NCT02394028)",
    notes          = paste0(
      "Baseline characteristics from Moein 2025 Tables S6 (continuous) ",
      "and S7 (categorical), induction column. Smoking status: 53% ",
      "non-smoker, 25% previously smoked, 22% current smoker. Disease ",
      "location: 60% ileum and colon, 19% ileum only, 21% colon only. ",
      "Other baseline medians: albumin 42.0 g/L, CRP 8.32 mg/L, fecal ",
      "calprotectin 800 ug/g, MAdCAM-1 17.6 U, white blood cells ",
      "8.10 x10^9/L, neutrophils 5.62 x10^9/L, ADA titer 0. This is a ",
      "SUBSET of the 864 CD patients in the companion population PK ",
      "analysis: only BERGAMOT induction patients with an evaluable ",
      "week-14 outcome are included."
    )
  )

  ini({
    # ==================================================================
    # Moein 2025 Supplementary Table S8, "Clinical Remission ->
    # Covariate-adjusted final model" column.
    #
    # SCALE OF THE INTERCEPT -- this is the one thing that is easy to get
    # wrong. Table S8's "Units" line reads "Intercept (Probability)",
    # and footnote a says the intercept "reflects the PROBABILITY for
    # the outcome of a patient treated with placebo" at the listed
    # reference covariates. So the printed 0.233 is a PROBABILITY, not a
    # log-odds, and it must be logit-transformed to become the linear
    # predictor's intercept. Everything else in the table is already on
    # the log-odds ("LO") scale.
    #
    # The probability reading is confirmed independently by the paper's
    # own Figure 4: reading the placebo (x = 0) intercepts off all eight
    # panels of that figure reproduces expit(logite0 + covariate shifts)
    # from Table S9 to within the line width. See the vignette.
    # ==================================================================
    logite0 <- logit(0.233); label("Logit of the probability of clinical remission at end of induction for the reference placebo patient (unitless logit)")  # Table S8: Intercept = 0.233 (RSE 16.4%, P < 0.01), reported on the PROBABILITY scale per the table's Units line and footnote a

    # ----- Exposure slope (log-odds per ug/mL) -----
    # NOT significant for this endpoint; retained because it is the
    # model's structural exposure term and the paper reports it as the
    # final-model estimate.
    e_ctrough_clinrem <- 0.0259; label("Log-odds of clinical remission at end of induction per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S8 'ER effect' = 0.0259 (RSE 126%); main-text Table 2 gives the same value with P = 0.426

    # ----- Covariate effect on the logit -----
    e_tnfnaive_clinrem <- 0.580; label("Log-odds shift for TNF-naive vs TNF-experienced patients (unitless logit)")  # Table S8 'TNF-naive' = 0.580 (RSE 39.2%, P < 0.05)

    # ----- No between-subject variability, no residual error -----
    # The source likelihood is Bernoulli: no omega and no sigma is
    # estimated. rxode2 requires an observation declaration, so a tiny
    # placeholder additive SD is fixed here purely to satisfy the
    # parser. It is NOT a source value; see the vignette Assumptions
    # and deviations section.
    addSd_prob_clinrem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor on the log-odds scale. The exposure term is
    # uncentered; the covariate term is a difference from the reference
    # category (TNF-experienced), so it vanishes for a reference
    # patient and logite0 alone remains.
    #
    # (1 - PRIOR_TNF) is the TNF-naive indicator -- see the PRIOR_TNF
    # covariateData note for why the contrast is formed this way rather
    # than by flipping the coefficient's sign.
    # ------------------------------------------------------------------
    logit_clinrem <-
      logite0 +
      e_ctrough_clinrem * CTROUGH +
      e_tnfnaive_clinrem * (1 - PRIOR_TNF)

    prob_clinrem <- expit(logit_clinrem)

    # Deterministic probability of clinical remission. Downstream
    # callers can draw binary outcomes with
    # rbinom(n, 1, prob_clinrem) on the rxSolve output.
    prob_clinrem ~ add(addSd_prob_clinrem)
  })
}
