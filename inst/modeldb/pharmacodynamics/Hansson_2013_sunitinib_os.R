Hansson_2013_sunitinib_os <- function() {
  description <- "Weibull time-to-event model for overall survival (OS) in adults with imatinib-resistant gastrointestinal stromal tumours (GIST) on sunitinib. The hazard function is h(t) = lambda * alpha * (lambda * t)^(alpha - 1) * exp(beta_anc * ANC + beta_dbprel * DBP_REL + beta_tumor * TUMSZ), where the three log-linear hazard modulators are the time-varying absolute neutrophil count ANC(t) (from the upstream Hansson 2013 myelosuppression model), the relative change in diastolic blood pressure from baseline DBP_REL(t) (from the upstream Hansson 2013 dBP indirect-response model), and the time-fixed baseline tumor size TUMSZ in mm (paper Table 2 'beta3 Tumor base'). All three modulators are consumed as data covariates. Time units inside the model are hours; the source paper reports the Weibull lambda in per-week units, converted to per-hour inside ini() so the parameter values match Table 2 row 'lambda (/week)'. A separate Weibull censoring distribution (lambdacens, alphacens) is exposed as a derived output for use in simulation-based dropout."
  reference <- paste(
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO.",
    "PKPD modeling of predictors for adverse effects and overall survival",
    "in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "Sister model files from the same paper:",
    "modellib('Hansson_2013_sunitinib_myelosuppression'),",
    "modellib('Hansson_2013_sunitinib_dbp'),",
    "modellib('Hansson_2013c_sunitinib') [fatigue],",
    "modellib('Hansson_2013_sunitinib_hfs').",
    sep = " "
  )
  vignette <- "Hansson_2013_sunitinib_os"
  units <- list(
    time          = "hour",
    dosing        = "n/a (no drug-dosing events; ANC(t) and DBP_REL(t) enter as time-varying covariates from the upstream Hansson 2013 myelosuppression and dBP models)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list(
    ANC = list(
      description        = "Time-varying absolute neutrophil count (10^9/L), typically simulated from the upstream Hansson_2013_sunitinib_myelosuppression model. Enters the hazard via beta_anc * ANC. Lower ANC -> lower hazard (paper: 'A more pronounced decrease in ANC over time ... decreased the hazard risk of death').",
      units              = "10^9/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Required input. The paper extrapolated ANC(t) from the individual myelosuppression-model predictions using AUC as the predictor, assuming dosing and schedule per protocol until time of censoring/death. For forward simulations either (a) simulate ANC(t) from `Hansson_2013_sunitinib_myelosuppression` and pass it in as a time-varying column, or (b) hold ANC at a typical-value level for sensitivity analyses.",
      source_name        = "ANC"
    ),
    DBP_REL = list(
      description        = "Time-varying relative change in diastolic blood pressure from baseline (unitless fraction; e.g., 0.10 = +10% above baseline). Typically simulated from the upstream Hansson_2013_sunitinib_dbp model as DBP_REL(t) = (dbp(t) - dbp0) / dbp0. Enters the hazard via beta_dbprel * DBP_REL. Larger relative increase -> lower hazard (paper: 'patients with a greater relative change in blood pressure ... displayed the longest OS').",
      units              = "fraction (unitless; positive when dBP elevated above baseline)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Required input. The paper extrapolated dBP_REL(t) from the individual dBP indirect-response model predictions. For forward simulations either (a) simulate dBP(t) from `Hansson_2013_sunitinib_dbp` and compute DBP_REL = (dbp - dbp0) / dbp0, or (b) hold DBP_REL at a typical-value level for sensitivity analyses.",
      source_name        = "DBPREL"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size at start of treatment (sum of longest diameters of target lesions, mm). Time-fixed per subject. Enters the hazard via beta_tumor * TUMSZ. Larger baseline tumor size -> higher hazard (paper: 'smaller tumor size at baseline, displayed the longest OS').",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Required input. Median baseline SLD by study (from Hansson 2013 e84 companion paper Table 1): 194 mm in study 1004, 108 mm in study 1047, 166 mm in study 1045, 255 mm in study 013. For typical-cohort simulations a midpoint of ~150 mm is a reasonable single-value choice.",
      source_name        = "TUMSZ"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies (Demetri 2006, George 2009, Shirao 2010, Maki 2005).",
    dose_range     = "sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 or continuous schedule (Table 1).",
    regions        = "multinational (study 1004) and Japanese (study 1045).",
    biomarkers     = "Overall survival (weeks from treatment start to death or censoring). Median (range) survival per study (Hansson 2013 Table 1): 61 (4-226) weeks in study 1004, 31 (15-81) in study 1047, 37 (27-48) in study 1045, 39 (4-96) in study 013. Total events: 163 deaths observed across the 303-patient pooled cohort.",
    notes          = "Hansson 2013 e85 Methods: 'The underlying distribution of the observed survival data was evaluated by exponential, Weibull, log-logistic, extreme value, and Gompertz probability density functions.' The Weibull was selected as the baseline-hazard form. Censoring was described by a separate Weibull (lambdacens = 0.0019/week, alphacens = 1.27 per Table 2) applied in simulations. n_events = 163 reported in the paper text ('The median number of simulated events (n = 151, range: 126-176) was in accordance with the number of observed events (n = 163)')."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Hansson 2013 e85 Table 2 'Survival model' block.
    # Paper Table 2 reports the Weibull scale lambda in per-week units.
    # Conversion to per-hour: lambda_h = lambda_w / (7 * 24) = lambda_w / 168.
    # ----------------------------------------------------------------------

    # Weibull scale (lambda).
    llam_haz  <- log(0.0079 / 168); label("Weibull baseline-hazard scale lambda (1/h; Hansson 2013 Table 2 row 'lambda (/week) = 0.0079')") # Table 2 lambda = 0.0079/week (RSE 55%); converted /168 -> 4.7024e-5 /h
    # Weibull shape (alpha).
    lalfa_haz <- log(1.15);         label("Weibull baseline-hazard shape alpha (unitless; Hansson 2013 Table 2 row 'alpha')")              # Table 2 alpha = 1.15 (RSE 9.1%)

    # Log-linear hazard modulators.
    e_anc_haz     <- 4.76;     label("Hazard log-linear coefficient on ANC (L/10^9 cells; Hansson 2013 Table 2 row 'beta1 ANC (l/.10^9)')") # Table 2 beta1 ANC = 4.76 (RSE 31%)
    e_dbprel_haz  <- -1.29;    label("Hazard log-linear coefficient on relative change in dBP (unitless per unit fractional dBP increase; Hansson 2013 Table 2 row 'beta2 dBPREL')") # Table 2 beta2 dBPREL = -1.29 (RSE 27%)
    e_tumor_haz   <- -0.00172; label("Hazard log-linear coefficient on baseline tumor size (1/mm; Hansson 2013 Table 2 row 'beta3 Tumor base (/mm)')") # Table 2 beta3 Tumor base = -0.00172 (RSE 46%)

    # Censoring Weibull (used for forward dropout simulation). Fixed in
    # this model file (the paper does not report estimated uncertainty in
    # a way that would parameterise the dropout sub-model separately;
    # values are point-only). Wrapped in fixed() per parameter-names.md
    # 'Fixed parameters' guidance.
    llamcens_haz  <- fixed(log(0.0019 / 168)); label("Censoring Weibull scale lambda_cens (1/h; Hansson 2013 Table 2 row 'lambda_cens (/week) = 0.0019')") # Table 2 lambda_cens = 0.0019/week (RSE 6.6%); converted /168 -> 1.131e-5 /h
    lalfacens_haz <- fixed(log(1.27));         label("Censoring Weibull shape alpha_cens (unitless; Hansson 2013 Table 2 row 'alpha_cens')")               # Table 2 alpha_cens = 1.27 (RSE 44%)

    # ----------------------------------------------------------------------
    # IIV. The paper does not report estimated IIV for the OS hazard
    # parameters (Table 2 'IIV CV%' column is dash for all survival-model
    # rows). No eta* parameters are added.
    # ----------------------------------------------------------------------

    # Residual error placeholder. The source NONMEM run uses LIKE (the
    # survival / event-density likelihood); no observation-error model.
    # This nlmixr2 translation is intended for forward simulation:
    # `hazard`, `cumhaz`, `sur`, `hazard_cens`, and `sur_cens` are exposed
    # as derived outputs and a tiny additive residual is attached to `sur`
    # so the nlmixr2 likelihood machinery accepts the model.
    addSd_sur <- 0.001; label("Placeholder additive residual error on the survival-probability output (unitless); not from the source -- see vignette Assumptions")
  })

  model({
    # Back-transformed structural parameters
    lam_haz       <- exp(llam_haz)
    alfa_haz      <- exp(lalfa_haz)
    lamcens_haz   <- exp(llamcens_haz)
    alfacens_haz  <- exp(lalfacens_haz)

    # Weibull baseline hazard with log-linear covariate modulators.
    # The DEL = 1e-6 small-time offset keeps the (lam * t)^(alpha - 1) term
    # finite at t = 0 without affecting the integrated cumulative hazard.
    # Mirrors the Zecchin_2016_survival convention.
    del <- 1e-6
    hazard <- lam_haz * alfa_haz * (lam_haz * (t + del))^(alfa_haz - 1) *
      exp(e_anc_haz * ANC + e_dbprel_haz * DBP_REL + e_tumor_haz * TUMSZ)

    # Cumulative hazard and survival probability (death).
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)

    # Censoring Weibull (no covariate modulators).
    hazard_cens <- lamcens_haz * alfacens_haz * (lamcens_haz * (t + del))^(alfacens_haz - 1)
    d/dt(cumhaz_cens) <- hazard_cens
    cumhaz_cens(0)    <- 0
    sur_cens          <- exp(-cumhaz_cens)

    # Forward-simulation observation: survival probability with a tiny
    # placeholder additive residual so the nlmixr2 likelihood machinery
    # accepts the model.
    sur ~ add(addSd_sur)
  })
}
