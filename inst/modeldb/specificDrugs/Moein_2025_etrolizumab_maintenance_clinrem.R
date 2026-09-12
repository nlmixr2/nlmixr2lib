Moein_2025_etrolizumab_maintenance_clinrem <- function() {
  description <- paste0(
    "Landmark binomial logistic exposure-response model for CLINICAL ",
    "REMISSION at the END OF MAINTENANCE (week 66) in adults with ",
    "moderately-to-severely active Crohn's disease treated with ",
    "etrolizumab (Moein 2025, n = 434, the maintenance ",
    "exposure-response analysis set of the phase 3 BERGAMOT study, ",
    "NCT02394028, placebo or 105 mg SC Q4W). Clinical remission is a ",
    "CDAI score < 150 with stool frequency mean daily score <= 3 and ",
    "abdominal pain mean daily score <= 1 and no worsening in either ",
    "subscore. UNLIKE the induction counterpart, the etrolizumab ",
    "exposure slope here IS significant (0.147 log-odds per ug/mL, RSE ",
    "36%, P = 0.00544) -- the paper's central finding is that ",
    "exposure-response is evident at the end of maintenance but not at ",
    "the end of induction. Four covariates are retained: baseline ",
    "albumin, female sex, TNF-naive status and baseline CDAI. There is ",
    "no PK layer and no ODE -- the exposure metric arrives as the ",
    "CTROUGH data column, the individual predicted week-4 trough after ",
    "a SINGLE 105 mg dose from ",
    "modellib('Moein_2025_etrolizumab'). Reproduces panels (a), (b), ",
    "(d) and (e) of Moein 2025 Figure 4. One of six landmark models in ",
    "the Moein_2025_etrolizumab_* family."
  )
  reference <- paste(
    "Moein A, Ribbing J, Ibrahim MMA, Zhang W, Kassir N.",
    "Population pharmacokinetics and exposure-response relationships of",
    "etrolizumab in patients with moderately-to-severely active Crohn's",
    "disease.",
    "J Clin Pharmacol. 2025;65(10):1208-1219.",
    "doi:10.1002/jcph.70043.",
    "Coefficients from Supplementary Table S9, 'Clinical Remission /",
    "Covariate-adjusted final model' column; the exposure slope is also",
    "printed in main-text Table 2. The reference-patient covariate",
    "values are Table S9 footnote a. Graphical check: Figure 4a, 4b, 4d",
    "and 4e.",
    "Exposure metric supplied by modellib('Moein_2025_etrolizumab').",
    sep = " "
  )
  vignette <- "Moein_2025_etrolizumab"
  units <- list(
    time          = "n/a (static landmark logistic regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_clinrem (probability of clinical remission at end of maintenance, 0-1)"
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
        "210 mg during induction, because every maintenance patient was ",
        "given 105 mg in that phase; the induction companion models ",
        "instead use the patient's assigned 105 or 210 mg, so the two ",
        "columns are NOT interchangeable. The single-dose basis removes ",
        "the confounding that an on-treatment trough would carry, since ",
        "etrolizumab clearance declines over time in a way that ",
        "correlates with clinical improvement. Obtain it by solving a ",
        "single 105 mg dose to day 28 with ",
        "modellib('Moein_2025_etrolizumab'). Placebo subjects take ",
        "CTROUGH = 0. Observed distribution (Table S6, maintenance): ",
        "mean 1.66, median 0.164, range 0-9.14 ug/mL -- the median sits ",
        "near zero because roughly half of this analysis set was ",
        "randomised to placebo for maintenance. The Figure 4 panels ",
        "span 0-9 ug/mL."
      ),
      source_name        = "Ctrough,W4,adjusted"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Linear effect on the logit, CENTERED at the reference of ",
        "42.0 g/L: e_alb_clinrem * (ALB - 42.0). Table S9 footnote a ",
        "defines the intercept patient as having 'Albumin (g/L): 42.0', ",
        "and the centering is what makes the printed intercept a ",
        "probability at the reference. Higher baseline albumin predicts ",
        "a HIGHER probability of clinical remission (+0.0822 log-odds ",
        "per g/L). Albumin is one of the paper's pre-specified ",
        "CONFOUNDING covariates -- it is associated with both ",
        "etrolizumab clearance (lower albumin raises CL and so lowers ",
        "exposure) and with the clinical outcome, which is precisely ",
        "why the full covariate model adjusts for it. Reproduced in ",
        "Figure 4a. Analysis-set distribution (Table S6): median 42.0, ",
        "range 27.0-53.0."
      ),
      source_name        = "Albumin"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the paper's reference; Table S9 footnote a states 'Sex: Male')",
      notes              = paste0(
        "Log-odds shift for female versus male (-0.563). The canonical ",
        "SEXF orientation (1 = female) matches the published ",
        "coefficient's orientation exactly, so no sign flip or ",
        "intercept shift is needed. Female patients had a LOWER ",
        "predicted probability of clinical remission; the Discussion ",
        "notes that the role of sex in IBD progression remains ",
        "inconclusive. Reproduced in Figure 4d. Analysis set ",
        "(Table S7): 218 male (50%), 216 female (50%)."
      ),
      source_name        = "Sex"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF biologic therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (TNF-experienced) -- the PAPER's reference, the complement of the canonical column's 0 level",
      notes              = paste0(
        "The published coefficient is on the TNF-NAIVE side of the ",
        "contrast (Table S9 'TNF-naive' = 0.564, RSE 40.0%, P < 0.05), ",
        "with TNF-experienced absorbed into the intercept (footnote a). ",
        "The canonical PRIOR_TNF is 1 for TNF-experienced, so the naive ",
        "indicator is formed in model() as (1 - PRIOR_TNF) with the ",
        "published coefficient and intercept UNCHANGED -- the ",
        "register's SMOKE_NEVER idiom, preserving the published ",
        "parameter covariance. Together with albumin, prior TNF status ",
        "is called out in the Conclusions as a significant PROGNOSTIC ",
        "factor for clinical remission during maintenance. Reproduced ",
        "in Figure 4e. Analysis set (Table S7): 254 TNF-experienced ",
        "(59%), 180 TNF-naive (41%)."
      ),
      source_name        = "TNF status"
    ),
    SCORE_CDAI = list(
      description        = "Baseline Crohn's Disease Activity Index (CDAI) score",
      units              = "(score, 0-600)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Linear effect on the logit, CENTERED at the MAINTENANCE ",
        "reference score of 320: e_cdai_clinrem * (SCORE_CDAI - 320). ",
        "The reference differs by phase -- Table S9 (maintenance) uses ",
        "320 while Table S8 (induction) uses 322 -- so the induction ",
        "companion models center at 322. Higher baseline disease ",
        "activity lowers the probability of remission, consistent with ",
        "the Discussion's comparison to certolizumab pegol. Reproduced ",
        "in Figure 4b. Analysis-set distribution (Table S6): median ",
        "320, range 215-481."
      ),
      source_name        = "CDAI score"
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
      "Entry into the maintenance phase required a CDAI-70 response at ",
      "the end of induction, so this is a RESPONDER-ENRICHED set and is ",
      "not a random subset of the induction population. Patients who ",
      "received placebo during induction were not randomised into the ",
      "maintenance phase and are excluded from this analysis. Baseline ",
      "characteristics from Moein 2025 Tables S6 and S7, maintenance ",
      "column. Smoking status: 51% non-smoker, 28% previously smoked, ",
      "21% current smoker. Disease location: 62% ileum and colon, 18% ",
      "ileum only, 21% colon only. Table S9 footnote a reference ",
      "patient: CDAI 320, SES-CD 12.0, fecal calprotectin 832 ug/g, CRP ",
      "8.30 mg/L, albumin 42.0 g/L, white blood cells 7.85 x10^9/L, ",
      "neutrophils 5.40 x10^9/L, MAdCAM-1 18.0 U, disease location ",
      "ileum and colon, non-smoker, male, TNF-experienced."
    )
  )

  ini({
    # ==================================================================
    # Moein 2025 Supplementary Table S9, "Clinical Remission ->
    # Covariate-adjusted final model" column.
    #
    # The intercept is printed on the PROBABILITY scale (Table S9 Units
    # line: "Intercept (Probability)"; footnote a: "reflects the
    # probability for the outcome of a patient treated with placebo" at
    # the listed reference covariates) and is logit-transformed here.
    # All other coefficients are already log-odds.
    #
    # This reading is confirmed by the paper's own Figure 4: at the
    # placebo end of the x-axis, panel (e) reads about 0.245 for
    # TNF-experienced and 0.355 for TNF-naive, against
    # expit(logit(0.244)) = 0.244 and expit(logit(0.244) + 0.564) =
    # 0.362 from this parameterisation. Panels (a), (b) and (d) agree
    # to the same tolerance.
    # ==================================================================
    logite0 <- logit(0.244); label("Logit of the probability of clinical remission at end of maintenance for the reference placebo patient (unitless logit)")  # Table S9: Intercept = 0.244 (RSE 18.1%, P < 0.01), on the PROBABILITY scale per the Units line and footnote a

    # ----- Exposure slope (log-odds per ug/mL) -----
    # SIGNIFICANT for this endpoint, unlike the induction counterpart.
    e_ctrough_clinrem <- 0.147; label("Log-odds of clinical remission at end of maintenance per ug/mL of single-dose week-4 trough (unitless logit)")  # Table S9 'ER effect' = 0.147 (RSE 36%, P < 0.01); main-text Table 2 gives the same value with P = 0.00544

    # ----- Covariate effects on the logit -----
    e_alb_clinrem      <- 0.0822;   label("Log-odds shift per g/L of albumin above the reference of 42.0 g/L (unitless logit)")     # Table S9 'Albumin Effect' = 0.0822 (RSE 35.3%, P < 0.01)
    e_sexf_clinrem     <- -0.563;   label("Log-odds shift for female vs male sex (unitless logit)")                                 # Table S9 'Female' = -0.563 (RSE 39.9%, P < 0.05)
    e_tnfnaive_clinrem <- 0.564;    label("Log-odds shift for TNF-naive vs TNF-experienced patients (unitless logit)")               # Table S9 'TNF-naive' = 0.564 (RSE 40.0%, P < 0.05)
    e_cdai_clinrem     <- -0.00379; label("Log-odds shift per CDAI point above the maintenance reference score of 320 (unitless logit)")  # Table S9 'CDAI effect' = -0.00379 (RSE 50.1%, P < 0.05); units 'LO per score'

    # ----- No between-subject variability, no residual error -----
    addSd_prob_clinrem <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor on the log-odds scale. Continuous covariates are
    # centered at their reference values and the categorical ones are
    # differences from the reference category, so every covariate term
    # vanishes for a reference patient and logite0 alone remains.
    logit_clinrem <-
      logite0 +
      e_ctrough_clinrem * CTROUGH +
      e_alb_clinrem * (ALB - 42.0) +
      e_sexf_clinrem * SEXF +
      e_tnfnaive_clinrem * (1 - PRIOR_TNF) +
      e_cdai_clinrem * (SCORE_CDAI - 320)

    prob_clinrem <- expit(logit_clinrem)

    prob_clinrem ~ add(addSd_prob_clinrem)
  })
}
