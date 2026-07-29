Zhu_2023_omalizumab_pediatric <- function() {
  description <- "Population pharmacodynamic indirect-response (IDR Type IV) model for forced expiratory volume in 1 second (FEV1, percent predicted on a 0-1 fractional scale) driven by serum free IgE in pediatric patients (6-11 years) with moderate to severe persistent inadequately controlled allergic asthma treated with omalizumab (Zhu 2023). PD-only -- free IgE concentration enters as the exogenous time-varying covariate IGE_FREE (per-FEV1-observation interpolated value, ng/mL). Adapted from the adult/adolescent IgE-FEV1 model (Lowe et al. 2009) with simplifications for the sparser pediatric data: observed baseline FEV1 per subject, common IIV magnitude on Imax and FEV1max, and Hill coefficient gamma fixed at 9. The estimated IC50 in pediatrics (39.4 ng/mL) is higher than in adults/adolescents (19.8 ng/mL); the model still supports the 25 ng/mL free IgE target underpinning the Xolair dosing table."

  reference <- paste(
    "Zhu R, Wang X, Anderson E, Deng M, Pivirotto S, Jin J, Kassir N, Owen R.",
    "Population-Based Pharmacodynamic Modeling of Omalizumab in Pediatric Patients",
    "with Moderate to Severe Persistent Inadequately Controlled Allergic Asthma.",
    "AAPS J. 2023;25(4):56. doi:10.1208/s12248-023-00823-4 (PMID 37266853).",
    sep = " "
  )

  vignette <- "Zhu_2023_omalizumab_pediatric"

  units <- list(
    time          = "week",
    dosing        = "n/a (PD-only model; omalizumab and free IgE PK provided externally via the IGE_FREE covariate)",
    concentration = "fraction predicted (FEV1 percent predicted on the 0-1 fractional scale; multiply by 100 to display on the conventional 0-100 percent scale)"
  )

  covariateData <- list(
    IGE_FREE = list(
      description        = "Serum free IgE concentration at the FEV1 observation time. The exogenous time-varying PD driver: it enters the inhibition Hill term imax * IGE_FREE^hill / (ec50^hill + IGE_FREE^hill) to suppress steady-state FEV1pp.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per FEV1 observation. The Zhu 2023 source paper preprocessed the IgE input as follows (Methods 'Data Used in the Population Pediatric IgE-FEV1 Model'). For omalizumab-treated subjects: total IgE at week 0 was used pretreatment (free IgE = total IgE before anti-IgE drug), assayed free IgE was used post-treatment, and the free IgE value at week 0.1 was imputed by back-extrapolation from the first post-treatment free IgE observation to capture the rapid suppression of free IgE within one day of subcutaneous omalizumab dosing. Linear (and log-linear, tested as a sensitivity analysis) interpolation between observed IgE values produced the per-FEV1-observation input value. For placebo subjects: the subject's AVERAGE total IgE over the 24-week steroid-stable period was used (constant per subject across all observation rows), because free IgE equals total IgE in the absence of anti-IgE drug. Downstream nlmixr2lib users combining this PD model with an upstream omalizumab popPK / IgE-binding model can populate IGE_FREE from the freeIgE observable of Hayashi_2007_omalizumab.R or any equivalent sequential PK source.",
      source_name        = "C_IgE,i(t) (paper symbol in the structural-model equation; NONMEM column name not disclosed)"
    )
  )

  # Documented-but-unused covariates -- the paper screened these but did
  # not retain any in the final model. They are kept here so the
  # provenance of the paper's covariate screen is preserved without
  # triggering an unused-covariate convention warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age at study entry",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the pediatric IgE-FEV1 model with a power model (age normalized by reference value) on Imax, IC50, and FEV1max (Zhu 2023 Methods 'Model Development and Covariate Analysis'). The covariate models 'either did not converge or appeared to reach local minimums during estimation, indicating that the data did not support identification of covariate effects'. Dropped from the final model.",
      source_name        = "Age"
    ),
    WT = list(
      description        = "Body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed among covariates considered for evaluation but not formally tested because of collinearity with age in the pediatric cohort. Paper Methods: 'Given the different response pattern (in both magnitude and onset) observed in pediatrics and the collinearity between covariates (e.g., high correlation between age and body weight in pediatric subjects), only the age effect was examined during the pediatric model development'.",
      source_name        = "Body weight"
    ),
    FEV1PP_BASELINE = list(
      description        = "Baseline FEV1 (percent predicted) at study entry",
      units              = "fraction predicted (0-1 scale)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed as a covariate of interest. Observed baseline FEV1 was used directly per-subject as the initial condition for the FEV1 ODE in the fitting process (Zhu 2023 Results: 'observed baseline FEV1 was used rather than an estimated value for each individual'), so the baseline does not enter the final model as a separate covariate effect on a structural parameter.",
      source_name        = "Baseline FEV1 (% predicted)"
    ),
    IGE_BASELINE = list(
      description        = "Baseline serum total IgE at screening",
      units              = "ng/mL (paper also reports IU/mL; conversion 1 IU/mL = 2.42 ng/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed as a covariate of interest. Visual inspection of IIV versus baseline IgE 'did not reveal any clear trends' (Supplementary Figure S1) and the covariate was not retained in the final model.",
      source_name        = "Baseline IgE (screening)"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 535L,
    n_studies         = 1L,
    study_names       = c("IA05 (omalizumab in pediatric moderate-to-severe asthma; 24-week steroid-stable phase used for the IgE-FEV1 model)"),
    n_observations    = 1515L,
    n_omalizumab_obs  = 960L,
    n_placebo_obs     = 555L,
    n_omalizumab_subj = 351L,
    n_placebo_subj    = 184L,
    age_range         = "6-11 years",
    age_median        = "9 years",
    weight_range      = "19.3-81.3 kg",
    weight_median     = "31 kg",
    height_range      = "104-168 cm",
    height_median     = "134 cm",
    sex_female_pct    = NA_real_,
    race_ethnicity    = NA_character_,
    disease_state     = "Children aged 6-11 with moderate-to-severe, persistent, inadequately controlled allergic asthma. Baseline FEV1 (% predicted) mean (range) 87% (25-148%); baseline serum IgE mean (range) 950 (36-4500) ng/mL.",
    dose_range        = "Omalizumab 75-375 mg SC every 2 or 4 weeks (dose and dosing frequency determined per the US package insert table from body weight + baseline IgE at screening, designed to achieve individualized free IgE suppression).",
    regions           = "Multinational pediatric trial (Study IA05).",
    notes             = "n=535 pediatric subjects from Study IA05 (Phase III randomized, double-blind, placebo-controlled; full enrolment 627 patients, of which 535 contributed data over the 24-week steroid-stable period used for the IgE-FEV1 model). Demographics from Table I of Zhu 2023. The 24-week steroid-stable phase was used to control for the confounding effect of steroid use on FEV1."
  )

  ini({
    # -----------------------------------------------------------
    # Structural model (Zhu 2023 PDF page 2; same form as the
    # adult/adolescent model of Lowe et al. 2009):
    #
    #   dFEV1_i/dt = Kdet_i * FEV1max_i
    #                  * (1 - Imax_i * C_IgE_i(t)^gamma / (IC50_i^gamma + C_IgE_i(t)^gamma))
    #                - Kdet_i * FEV1_i(t)
    #
    # This is an IDR Type IV form: the drug signal (here, free IgE)
    # inhibits the FEV1pp production rate kin = Kdet * FEV1max while the
    # loss rate kout = Kdet is unchanged. Mapping to nlmixr2lib
    # canonicals: Kdet -> lkout (canonical IDR rate constant), FEV1max
    # -> lrbase (canonical baseline value of a turnover state), Imax ->
    # logitimax (extending the logitemax pattern from
    # Kim_2017_fimasartan.R; paper's NONMEM parameterisation is Imax =
    # exp(theta2) / (1 + exp(theta2))), IC50 -> lec50 (canonical sigmoid
    # half-effect concentration), gamma -> lhill (canonical Hill
    # exponent, fixed at log(9) per Table II).
    # -----------------------------------------------------------
    lkout     <- log(0.0198)     ; label("FEV1pp turnover rate constant Kdet (1/week)")                                  # Zhu 2023 Table II: Kdet = 0.0198 1/week (95% CI 0.0163, 0.0240; %RSE 9.95)
    lrbase    <- log(0.879)      ; label("Maximum FEV1pp at free IgE = 0 (FEV1max, fraction predicted)")                 # Zhu 2023 Table II: FEV1max = 0.879 (95% CI 0.844, 0.915; %RSE 2.08; fraction predicted -- multiply by 100 for the 0-100 percent scale)
    logitimax <- qlogis(0.0717)  ; label("Logit-scale maximum IgE inhibitory effect (Imax, unitless; constrained to (0,1))")  # Zhu 2023 Table II: Imax = 0.0717 (95% CI 0.0288, 0.168); logit transform per paper Methods constrains Imax to (0,1); theta2 = qlogis(0.0717) approx -2.560
    lec50     <- log(39.4)       ; label("Free IgE concentration causing 50% maximum inhibition (IC50, ng/mL)")           # Zhu 2023 Table II: IC50 = 39.4 ng/mL (95% CI 24.3, 63.9; %RSE 25.0)
    lhill     <- fixed(log(9))   ; label("Hill coefficient gamma (unitless) -- FIXED at 9")                                # Zhu 2023 Table II: gamma = 9 (FIXED); sensitivity analysis examined values 1-9 and found 'values greater than 5 gave similar fits' (Results paragraph); the adult/adolescent estimate was 15.9 (%RSE 38.0) but a value of 9 was used in the pediatric model for stability with the sparser data

    # -----------------------------------------------------------
    # IIV: the same variance is estimated on logit-scale Imax and
    # log-scale FEV1max (Zhu 2023 Results 'Model Development and
    # Covariate Analysis': 'the IIV of model parameters was reduced to a
    # common parameter, i.e., the same IIV magnitude was estimated on
    # Imax and FEV1max'). Table II reports IIV variance 0.0952
    # (SD 0.309) on both Omega(2,2) and Omega(3,3) with %RSE 16.3 on the
    # variance estimate. Shrinkage was high for Imax IIV (96.5%; sparse
    # data) and moderate for FEV1max IIV (20.3%). The two etas are
    # declared with the same starting value here so that a re-fit honors
    # the published shared-variance constraint; users wanting to relax
    # the constraint can simply re-fit with the etas free.
    # -----------------------------------------------------------
    etalogitimax ~ 0.0952           # Zhu 2023 Table II Omega(2,2); shared with Omega(3,3) per paper-stated equality constraint
    etalrbase    ~ 0.0952           # Zhu 2023 Table II Omega(3,3); same starting value as Omega(2,2)

    # -----------------------------------------------------------
    # Residual error: additive on the 0-1 fractional FEV1pp scale (Zhu
    # 2023 Table II additive SD = 0.0704; 95% CI 0.0674, 0.0734;
    # %RSE 2.18; shrinkage 12.0%). Equivalent to 7.04 percentage points
    # on the conventional 0-100 percent-predicted scale.
    # -----------------------------------------------------------
    addSd     <- 0.0704          ; label("Additive residual SD on FEV1pp (fraction predicted)")                          # Zhu 2023 Table II: additive SD = 0.0704 (95% CI 0.0674, 0.0734; %RSE 2.18); on the 0-1 fractional FEV1pp scale
  })

  model({
    # ----- Individual parameters -----
    # Imax via the inverse-logit (expit) so the individual Imax is
    # bounded to (0, 1). Paper Methods (footnote to Table II):
    # 'parameters modeled using logit and lognormal transformations'.
    imax  <- expit(logitimax + etalogitimax)
    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout)
    ec50  <- exp(lec50)
    hill  <- exp(lhill)

    # ----- Sigmoid Imax inhibition driven by serum free IgE -----
    # Paper Eq. PDF page 2 (inside the parentheses):
    #   inhibition fraction = Imax * C_IgE(t)^gamma / (IC50^gamma + C_IgE(t)^gamma)
    # IGE_FREE is the time-varying free-IgE covariate column (ng/mL)
    # supplied per-observation by the dataset; canonical name registered
    # in inst/references/covariate-columns.md (founding example: this
    # model). For sequential PK-PD use, IGE_FREE is populated from an
    # upstream omalizumab popPK / binding model (Hayashi_2007_omalizumab
    # emits freeIgE in ng/mL); for standalone use, the column carries
    # the interpolated assayed free-IgE values per the paper's
    # preprocessing (see covariateData$IGE_FREE$notes).
    inhibition <- imax * IGE_FREE^hill / (ec50^hill + IGE_FREE^hill)

    # ----- IDR Type IV ODE -----
    # Paper Eq. PDF page 2 (verbatim form):
    #   d/dt(fev1pp) = Kdet * FEV1max * (1 - inhibition) - Kdet * fev1pp
    d/dt(fev1pp) <- kout * rbase * (1 - inhibition) - kout * fev1pp

    # Initial condition: the IgE-free maximum FEV1pp = rbase. The Zhu
    # 2023 fitting workflow used each subject's observed baseline
    # FEV1pp instead (Results: 'observed baseline FEV1 was used rather
    # than an estimated value for each individual'); for nlmixr2lib
    # simulations without observed baselines, rbase is a self-consistent
    # default (the model's structural maximum FEV1pp). To replicate the
    # paper's per-subject baseline-as-initial-condition behaviour,
    # downstream callers can override fev1pp(0) via rxode2 event-table
    # mechanisms or via a Bayesian posterior of the observed baseline.
    fev1pp(0) <- rbase

    # ----- Observation and residual error -----
    fev1pp ~ add(addSd)
  })
}
