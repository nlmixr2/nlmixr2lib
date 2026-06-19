Sherer_2012_AAA <- function() {
  description <- "Hierarchical Bayesian disease-progression model of abdominal aortic aneurysm (AAA) diameter in men with small (30-49 mm) AAA identified during ultrasound screening (Sherer 2012). The expected AAA diameter follows the ODE dY/dt = beta1 + beta2 * (Y - beta0) with Y(0) = beta0, so the growth rate changes at a constant rate with change in size; the three individual-level parameters (baseline size beta0, baseline growth rate beta1, and constant first derivative of growth rate with size beta2) are drawn from a multivariate normal distribution with a full 3x3 covariance. Covariate effects: baseline AAA diameter on all three parameters, log10 plasma D-dimer on beta1 and beta2, and a diabetes-mellitus binary on beta2. Disease-progression model with no drug dosing."

  reference <- paste(
    "Sherer EA, Bies RR, Clancy P, Norman PE, Golledge J.",
    "Growth of screen-detected abdominal aortic aneurysms in men: a Bayesian analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2012 Oct 24;1(10):e12.",
    "doi:10.1038/psp.2012.13.",
    sep = " "
  )
  vignette <- "Sherer_2012_AAA"
  units <- list(
    time          = "year",
    dosing        = "n/a (disease-progression model with no drug dosing)",
    concentration = "mm (abdominal aortic aneurysm diameter, observation aaaSize)"
  )

  covariateData <- list(
    AAA_DIAM = list(
      description        = "Baseline (screening) abdominal aortic aneurysm diameter; the per-subject time-fixed measurement that anchors the ODE initial condition and enters the regression equations for all three individual-level parameters (beta0, beta1, beta2).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; the value is the screening ultrasound diameter (men's HIMS cohort: 30-49 mm small AAA inclusion criterion). Enters the typical-value regressions as the proportional term AAA_DIAM / 32.7, where 32.7 mm is the cohort median baseline diameter (Sherer 2012 Table 1). The same column also seeds the ODE state at t = 0 indirectly via beta0 = e_aaadiam_b0 * (AAA_DIAM / 32.7) + etabeta0.",
      source_name        = "Y(0)"
    ),
    DDIMER = list(
      description        = "Plasma D-dimer concentration measured at the HIMS follow-up blood sampling (2001-2004), used as a baseline covariate on the AAA growth rate parameters via log10-transformation.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (single measurement at follow-up rather than at AAA-screening baseline; Sherer 2012 Discussion explicitly flags this as a limitation). Enters the typical-value regressions as the proportional term log10(DDIMER) / log10(326), where 326 ng/mL is the cohort median D-dimer (Sherer 2012 Table 1) and log10(326) approx 2.513. The cohort interquartile range is 142-785 ng/mL; the model was developed against this range and extrapolation outside it is not validated.",
      source_name        = "C^(D-dimer)"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator (Type 1 or Type 2 not distinguished). 1 = patient has diabetes; 0 = no diabetes.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Time-fixed at HIMS screening per subject (self-reported on the lifestyle questionnaire alongside hypertension, coronary heart disease, stroke). Affects only beta2 (first derivative of growth rate with size); the Sherer 2012 forward-selection / backward-elimination did not find a significant effect on beta1. Cohort prevalence 14% (42 / 299 men, Table 1).",
      source_name        = "Diabetes"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 299L,
    n_studies      = 1L,
    age_range      = "65-83 years (HIMS screening 1996-1999 enrolment criterion)",
    age_median     = "72 years",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = "predominantly white (Sherer 2012 Discussion)",
    disease_state  = "Men with a small (30-49 mm) abdominal aortic aneurysm detected on ultrasound screening, with both serial AAA-diameter measurements and a plasma D-dimer measurement available. Cohort drawn from the Health in Men Study (HIMS) screening study in Perth, Western Australia.",
    dose_range     = "Not applicable -- disease-progression model with no drug dosing.",
    regions        = "Perth, Western Australia",
    aaa_bl_range   = "30-49 mm (small AAA inclusion criterion); cohort median 32.7 mm (q1 30.8, q3 36.0)",
    ddimer_range   = "cohort median 326 ng/mL (q1 142, q3 785)",
    diab_prev_pct  = 14.0,
    htn_prev_pct   = 49.0,
    chd_prev_pct   = 38.0,
    smoke_prev_pct = 84.0,
    followup_dur   = "median 5.5 years (q1 5, q3 6)",
    n_observations = "1,732 AAA size measurements (median 6 per patient; q1 6, q3 7)",
    notes          = "Subset of 875 men diagnosed with small AAA during the Western Australia HIMS screening study who had both serial diameter measurements and a D-dimer measurement (Sherer 2012 Results). Demographics, medical conditions, and blood biochemistry summarised in Sherer 2012 Table 1; the model file's covariate set keeps the three covariates retained in the final model (AAA_DIAM, DDIMER, DIS_DIAB)."
  )

  ini({
    # ============================================================
    # Final model: 'First derivative of AAA growth rate with size is
    # constant' (Sherer 2012 Methods, page 7). Individual-level
    # parameters beta0 (baseline size, mm), beta1 (baseline growth
    # rate, mm/year), beta2 (constant first derivative of growth
    # rate with size, 1/year) are drawn from MVN(beta_bar, Sigma);
    # beta_bar values are reproduced from Sherer 2012 Table 3
    # (fixed-parameter median posterior estimates) and the
    # supporting regression equations in the Results section
    # (page 2).
    # ============================================================

    # ---------------- Beta0 (baseline AAA size, mm) ----------------
    # Equation (Sherer 2012 page 2): beta_bar_0 = e_aaadiam_b0 * (Y(0) / median(Y(0)))
    # The slope is reported in Results as "10.0 mm per 10 mm
    # increase" (95% CrI 9.95-10.01 mm), i.e. a proportional
    # coefficient of approximately 1.0 between observed baseline
    # diameter and beta0. The Table 3 value of 32.6 mm equals the
    # value of beta_bar_0 evaluated at Y(0) = median(Y(0)) = 32.7
    # mm.
    e_aaadiam_b0 <- 32.6  ; label("Baseline AAA size parameter at median Y(0) = 32.7 mm (mm)")  # Table 3 row 1 (median 32.6, 95% CrI 32.5-32.7)

    # ---------------- Beta1 (baseline growth rate, mm/year) --------
    # Equation: beta_bar_1 = b1_int
    #                       + e_aaadiam_b1 * (Y(0) / median(Y(0)))
    #                       + e_ddimer_b1  * (log10(C) / median(log10(C)))
    # Population baseline growth rate at median covariates is
    # b1_int + e_aaadiam_b1 + e_ddimer_b1 = -1.61 + 2.03 + 0.90
    # = 1.32 mm/year (matches Sherer 2012 Discussion p. 4: "a
    # relatively low baseline growth rate to the measurement
    # variability ratio (1.32 mm/year vs. 0.97 mm)").
    b1_int        <- -1.61 ; label("Beta1 offset from covariate effects (mm/year)")                       # Table 3 row 2 (median -1.61, 95% CrI -3.08, -0.30)
    e_aaadiam_b1  <-  2.03 ; label("Beta1 effect of (Y(0)/median Y(0)) (mm/year)")                        # Table 3 row 3 (median 2.03, 95% CrI 0.87, 3.40)
    e_ddimer_b1   <-  0.90 ; label("Beta1 effect of (log10 D-dimer / median log10 D-dimer) (mm/year)")    # Table 3 row 4 (median 0.90, 95% CrI 0.11, 1.64)

    # ---------------- Beta2 (first derivative of growth rate
    # with size, 1/year) -------------------------------------------
    # Equation: beta_bar_2 = b2_int
    #                       + e_aaadiam_b2 * (Y(0) / median(Y(0)))
    #                       + e_ddimer_b2  * (log10(C) / median(log10(C)))
    #                       + e_diab_b2    * Diabetes
    b2_int        <- -1.05 ; label("Beta2 offset from covariate effects (1/year)")                        # Table 3 row 5 (median -1.05, 95% CrI -1.52, -0.53)
    e_aaadiam_b2  <-  0.59 ; label("Beta2 effect of (Y(0)/median Y(0)) (1/year)")                         # Table 3 row 6 (median 0.59, 95% CrI 0.11, 1.03)
    e_ddimer_b2   <-  0.37 ; label("Beta2 effect of (log10 D-dimer / median log10 D-dimer) (1/year)")     # Table 3 row 7 (median 0.37, 95% CrI 0.13, 0.62)
    e_diab_b2     <- -0.32 ; label("Beta2 effect for diabetes (DIS_DIAB = 1) (1/year)")                       # Table 3 row 8 (median -0.32, 95% CrI -0.45, -0.18)

    # ---------------- Random effects ------------------------------
    # MVN(0, Sigma) on the additive deviations of (beta0, beta1,
    # beta2) from their typical values (Sherer 2012 Methods page 7,
    # "multivariate normal distribution ... covariance matrix").
    # The etas are named to pair with the ini() anchor for each
    # individual-level parameter so checkModelConventions() can
    # match them: etae_aaadiam_b0 sits next to e_aaadiam_b0 (the
    # purely-proportional regression slope for beta0), etab1_int
    # next to b1_int, etab2_int next to b2_int. Inside model() each
    # eta is added at the end of the typical-value regression so
    # the random effect attaches to the full beta_i, matching the
    # paper's beta_i = beta_bar_i + eta_i, eta_i ~ MVN(0, Sigma).
    # Lower-triangular Sigma in row-major order: var(beta0);
    # cov(beta0,beta1), var(beta1); cov(beta0,beta2),
    # cov(beta1,beta2), var(beta2). Values from Sherer 2012
    # Table 3 'Parameter covariance matrix' rows.
    etae_aaadiam_b0 + etab1_int + etab2_int ~ c(
      0.19,
      0.30,  1.11,
     -0.06, -0.15,  0.10
    )  # Table 3 Sigma_11..Sigma_33 (mm^2, mm^2/year, 1/year^2 and the three covariances)

    # ---------------- Residual error ------------------------------
    addSd <- 0.97  ; label("Additive residual SD on AAA diameter (mm)")  # Table 3 sigma_epsilon (median 0.97, 95% CrI 0.93-1.00)
  })

  model({
    # ----- Reference values used in the typical-value regressions
    # for (beta0, beta1, beta2). The Sherer 2012 regression
    # equations divide the raw covariate by its cohort median so
    # that each effect coefficient represents the value of beta_bar
    # contributed by a subject at the median covariate level.
    ref_aaadiam        <- 32.7              # median(Y(0))  -- Sherer 2012 Table 1 ('AAA diameter at baseline, median')
    ref_log10_ddimer   <- log10(326)        # median(log10 D-dimer) -- Sherer 2012 Table 1 ('D-dimer protein, median 326 ng/mL'); log10(326) approx 2.5132

    # ----- Normalised covariates (paper Eq. on page 2)
    aaadiam_norm <- AAA_DIAM / ref_aaadiam
    ddimer_norm  <- log10(DDIMER) / ref_log10_ddimer

    # ----- Subject-specific structural parameters. The random
    # deviation (etae_aaadiam_b0, etab1_int, etab2_int) is added at
    # the end of each typical-value regression so the eta attaches
    # to the full beta_i, not to the single ini() anchor parameter
    # whose name it shares for convention-pairing purposes.
    beta0 <- e_aaadiam_b0 * aaadiam_norm + etae_aaadiam_b0
    beta1 <- b1_int + e_aaadiam_b1 * aaadiam_norm + e_ddimer_b1 * ddimer_norm + etab1_int
    beta2 <- b2_int + e_aaadiam_b2 * aaadiam_norm + e_ddimer_b2 * ddimer_norm + e_diab_b2 * DIS_DIAB + etab2_int

    # ----- AAA growth ODE (Sherer 2012 Methods page 7,
    # 'First derivative of AAA growth rate with size is constant').
    # dY/dt = r(Y), with dr/dY = beta2 (a constant) and r(Y(0)) =
    # beta1, hence r(Y) = beta1 + beta2 * (Y - beta0) and dY/dt =
    # beta1 + beta2 * (Y - beta0). Initial condition is the
    # subject-specific beta0. The closed-form is
    # Y(t) = beta0 + (beta1 / beta2) * (exp(beta2 * t) - 1), which
    # reduces to beta0 + beta1 * t in the limit beta2 -> 0; the ODE
    # form is used here because beta2's population typical value is
    # near zero (Sherer 2012 Discussion: median beta2 approx
    # 0.06/year) and individual draws can be of either sign.
    d/dt(aaa) <- beta1 + beta2 * (aaa - beta0)
    aaa(0)    <- beta0

    aaaSize <- aaa
    aaaSize ~ add(addSd)
  })
}
