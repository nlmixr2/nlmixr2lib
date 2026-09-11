Moein_2025_etrolizumab_maintenance_endoimp <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for ENDOSCOPIC ",
    "IMPROVEMENT at the END OF MAINTENANCE (week 66) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 434, the maintenance ",
    "exposure-response analysis set of the phase 3 BERGAMOT study, ",
    "NCT02394028, placebo or 105 mg SC Q4W). Endoscopic improvement is ",
    "a >= 50% reduction from the baseline Simple Endoscopic Score for ",
    "Crohn's Disease (SES-CD). This endpoint carries the most ",
    "significant exposure slope of the six models (0.234 log-odds per ",
    "ug/mL, RSE 26%, P = 0.00012). Three covariates are retained: ",
    "ileum-only disease location and the two smoking-status ",
    "indicators, whose effects run in OPPOSITE directions -- former ",
    "smokers do worse and current smokers do better than never ",
    "smokers. There is no PK layer and no ODE -- the exposure metric ",
    "arrives as the CTROUGH data column, the individual predicted ",
    "week-4 trough after a SINGLE 105 mg dose from ",
    "modellib('Moein_2025_etrolizumab'). Reproduces panels (f) and (g) ",
    "of Moein 2025 Figure 4. One of six landmark models in the ",
    "Moein_2025_etrolizumab_* family."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S9, 'Endoscopic Improvement /",
    "Covariate-adjusted final model' column; the exposure slope is also",
    "printed in main-text Table 2. The reference-patient covariate",
    "values are Table S9 footnote b. Graphical check: Figure 4f and 4g.",
    "Exposure metric supplied by modellib('Moein_2025_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(
    time          = "n/a (static landmark logistic regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_endoimp (probability of endoscopic improvement at end of maintenance, 0-1)"
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
      reference_category = "0 (ileocolonic disease, Montreal L3, the paper's reference; implies DISLOC_COLON = 0 as well)",
      notes              = paste0(
        "Log-odds shift versus ileocolonic disease (-0.937 here). ",
        "Paired with DISLOC_COLON to encode the three-level Montreal ",
        "disease-location categorical, with both indicators zero ",
        "meaning ileum and colon (L3) -- the reference named in Table S9 ",
        "footnote b. Mutually exclusive with DISLOC_COLON. Colon-only ",
        "was not retained for this endpoint, so it is predicted ",
        "identically to the ileocolonic reference; the Figure 4g caption ",
        "says so explicitly and notes the colon-only curve is hidden ",
        "beneath the ileocolonic one. Analysis set (Table S7): 267 ",
        "ileum and colon (62%), 78 ileum only (18%), 89 colon only ",
        "(21%)."
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
        "exactly zero, because Table S9 retained no colon-only effect ",
        "for endoscopic improvement at maintenance (Figure 4g caption: ",
        "'Colon only was not significantly different from ileum and ",
        "colon, and therefore the two are predicted to be identical'). ",
        "Carrying the term explicitly with a zero coefficient keeps the ",
        "data-column contract identical across all six models of the ",
        "family. Mutually exclusive with DISLOC_ILEUM."
      ),
      source_name        = "Disease location"
    ),
    SMOKE_NEVER = list(
      description        = "Never-smoker indicator at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (never smoker) -- see notes; the PAPER's reference is the never-smoker group, which is the canonical pair's implicit THIRD level",
      notes              = paste0(
        "Used together with SMOKE_CURRENT to encode the three-level ",
        "smoking-status categorical. NOTE the reference-category ",
        "difference from the register's documented pairing: the ",
        "canonical SMOKE_NEVER / SMOKE_CURRENT pair leaves FORMER ",
        "smoker as the implicit reference, whereas Moein 2025 uses ",
        "NEVER smoker as the reference (Table S9 footnote b: 'Smoking ",
        "status: Non-smoker') and reports coefficients for the other ",
        "two levels. Rather than re-parameterise onto the register's ",
        "reference -- which would require shifting the intercept and ",
        "would break correspondence with the published parameter ",
        "covariance -- the former-smoker indicator is formed inside ",
        "model() from the two canonical columns as ",
        "(1 - SMOKE_NEVER - SMOKE_CURRENT), and both published ",
        "coefficients and the intercept are carried UNCHANGED. This is ",
        "the same construction the register documents for deriving an ",
        "ever-smoker indicator as (1 - SMOKE_NEVER). The three ",
        "indicators are mutually exclusive and exhaustive, so exactly ",
        "one of never / former / current holds per subject. Analysis ",
        "set (Table S7): 223 non-smoker (51%), 120 previously smoked ",
        "(28%), 91 current smoker (21%)."
      ),
      source_name        = "Smoking status"
    ),
    SMOKE_CURRENT = list(
      description        = "Current-smoker indicator at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a current smoker); the effect is stated versus the NEVER-smoker reference, so SMOKE_NEVER must be supplied alongside it",
      notes              = paste0(
        "Log-odds shift for current versus never smoker (+0.636). ",
        "Paired with SMOKE_NEVER as described in that entry's notes; ",
        "both columns are required because the former-smoker indicator ",
        "is derived from the pair. Current smokers had a HIGHER ",
        "predicted probability of endoscopic improvement while former ",
        "smokers had a LOWER one (-0.903), a non-monotone pattern the ",
        "Discussion flags as incompletely understood and warranting ",
        "further investigation -- so treat the direction with caution ",
        "rather than as a mechanistic claim. Reproduced in Figure 4f."
      ),
      source_name        = "Smoking status"
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
    # Moein 2025 Supplementary Table S9, "Endoscopic Improvement ->
    # Covariate-adjusted final model" column.
    #
    # The intercept is printed on the PROBABILITY scale (Table S9 Units
    # line: "Intercept (Probability)"; footnote b) and is
    # logit-transformed here; all other coefficients are already
    # log-odds.
    #
    # Figure 4f confirms the reading at the placebo end of the x-axis:
    # the panel reads about 0.14 for non-smokers, 0.23 for current
    # smokers and 0.06 for former smokers, against
    # expit(logit(0.138)) = 0.138,
    # expit(logit(0.138) + 0.636) = 0.233 and
    # expit(logit(0.138) - 0.903) = 0.061 from this parameterisation.
    # ==================================================================
    logite0 <- logit(0.138); label("Logit of the probability of endoscopic improvement at end of maintenance for the reference placebo patient (unitless logit)")  # Table S9: Intercept = 0.138 (RSE 12.8%, P < 0.01), on the PROBABILITY scale per the Units line and footnote b

    # ----- Exposure slope (log-odds per ug/mL) -----
    # The strongest exposure-response signal in the paper.
    e_ctrough_endoimp <- 0.234; label("Log-odds of endoscopic improvement at end of maintenance per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S9 'ER effect' = 0.234 (RSE 26.0%, P < 0.01); main-text Table 2 gives the same value with P = 0.00012

    # ----- Covariate effects on the logit -----
    e_ileum_endoimp       <- -0.937; label("Log-odds shift for ileum-only vs ileocolonic disease location (unitless logit)")  # Table S9 'Disease location Ileum only' = -0.937 (RSE 49.3%, P < 0.05)
    e_smokeformer_endoimp <- -0.903; label("Log-odds shift for former vs never smokers (unitless logit)")                     # Table S9 'Previously smoked' = -0.903 (RSE 43.8%, P < 0.05)
    e_smokecurrent_endoimp <- 0.636; label("Log-odds shift for current vs never smokers (unitless logit)")                    # Table S9 'Current smoker' = 0.636 (RSE 47.7%, P < 0.05)

    # Colon-only was NOT retained for this endpoint (Table S9 leaves the
    # cell blank; Figure 4g states colon-only is predicted identically to
    # the ileocolonic reference). Held at exactly zero, and marked fixed
    # so it reads unmistakably as a structural zero rather than an
    # estimate, to keep the data-column contract uniform across the
    # six-model family.
    e_colon_endoimp <- fixed(0); label("Log-odds shift for colon-only vs ileocolonic disease location, not retained by the source model and held at a structural zero (unitless logit)")  # Table S9: no colon-only row for endoscopic improvement; Figure 4g caption confirms colon only is predicted identically to ileum and colon

    # ----- No between-subject variability, no residual error -----
    addSd_prob_endoimp <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the log-odds scale. Every categorical term is
    # a difference from the reference category (ileocolonic, never
    # smoker), so all covariate terms vanish for a reference patient
    # and logite0 alone remains.
    #
    # The former-smoker indicator is derived from the canonical pair --
    # see the SMOKE_NEVER covariateData note for why the contrast is
    # built this way rather than by re-referencing the coefficients.
    logit_endoimp <-
      logite0 +
      e_ctrough_endoimp * CTROUGH +
      e_ileum_endoimp * DISLOC_ILEUM +
      e_colon_endoimp * DISLOC_COLON +
      e_smokeformer_endoimp * (1 - SMOKE_NEVER - SMOKE_CURRENT) +
      e_smokecurrent_endoimp * SMOKE_CURRENT

    prob_endoimp <- expit(logit_endoimp)

    prob_endoimp ~ add(addSd_prob_endoimp)
  })
}
