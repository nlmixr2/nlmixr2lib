Zecchin_2016_survival <- function() {
  description <- "Time-to-event model for overall survival (OS) in advanced epithelial ovarian cancer (Zecchin 2016 / DDMODEL00000218): Weibull baseline hazard with covariate effects of normalised baseline SLD (TUM_SLD / 70 mm), tumour-size-ratio TSR(t) capped at week 12, time-varying new-lesion indicator (NWLS), and binary ECOG performance status, with the underlying SLD trajectory (subject-specific tumour-growth and drug-cytotoxicity rate constants from the upstream Zecchin 2016 SLD model) integrated inline."
  reference <- paste(
    "Zecchin C, Gueorguieva I, Enas NH, Friberg LE.",
    "Models for change in tumour size, appearance of new lesions",
    "and survival probability in patients with advanced epithelial",
    "ovarian cancer.",
    "Br J Clin Pharmacol. 2016;82(3):717-727.",
    "doi:10.1111/bcp.12994.",
    "DDMORE Foundation Model Repository: DDMODEL00000218.",
    "Subject-specific tumour-dynamics covariates (KG, KD0, KD1, IBASE)",
    "are empirical-Bayes outputs from the upstream SLD model (DDMODEL00000217;",
    "see modellib('Zecchin_2016_tumorovarian')).",
    sep = " "
  )
  vignette <- "Zecchin_2016_survival"
  paper_specific_compartments <- c("wts")

  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; drug exposure enters as time-varying per-cycle AUC_CARBO and AUC_GEM)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  ddmore_id    <- "DDMODEL00000218"
  replicate_of <- NULL

  covariateData <- list(
    KG = list(
      description        = "Subject-specific tumour-size first-order growth rate constant carried over (as an empirical-Bayes posterior) from the upstream Zecchin 2016 SLD model (DDMODEL00000217). Used inside the inline SLD ODE term `KG / 1000 * tumor_size`.",
      units              = "(1/day) * 1000 (source NONMEM convention; the `/1000` rescaling is applied inside the model)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Supplied directly in the dataset (column `KG` in DDMODEL00000218 NONMEM `$INPUT`; identical column shipped in the bundle's Simulated_OS.csv). When this OS model is used standalone, populate `KG` per subject by first fitting `Zecchin_2016_tumorovarian` and extracting the per-subject empirical-Bayes posterior of `lkg + etalkg`. Bundle simulated-dataset typical values: ~0.85 (bundle subject 1).",
      source_name        = "KG"
    ),
    KD0 = list(
      description        = "Subject-specific carboplatin-driven tumour-size death rate constant carried over (as an empirical-Bayes posterior) from the upstream Zecchin 2016 SLD model. Pairs with `AUC_CARBO` inside the SLD ODE term `KD0 / 1000 * AUC_CARBO * tumor_size`.",
      units              = "(1/day per AUC_CARBO unit) * 1000 (source NONMEM convention; the `/1000` rescaling is applied inside the model)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Supplied directly in the dataset (column `KD0` in DDMODEL00000218 NONMEM `$INPUT`; identical column shipped in the bundle's Simulated_OS.csv). The source SLD model (DDMODEL00000217) wires a single shared `ETA(2)` onto both `KD0` and `KD1`, so per-subject `KD0` and `KD1` are correlated 1:1 in the empirical-Bayes posterior. Bundle simulated-dataset typical values: ~0.06 (bundle subject 1).",
      source_name        = "KD0"
    ),
    KD1 = list(
      description        = "Subject-specific gemcitabine-driven tumour-size death rate constant carried over (as an empirical-Bayes posterior) from the upstream Zecchin 2016 SLD model. Pairs with `AUC_GEM` inside the SLD ODE term `KD1 / 100 * AUC_GEM * tumor_size`.",
      units              = "(1/day per AUC_GEM unit) * 100 (source NONMEM convention; the `/100` rescaling is applied inside the model)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Supplied directly in the dataset (column `KD1` in DDMODEL00000218 NONMEM `$INPUT`; identical column shipped in the bundle's Simulated_OS.csv). Shares ETA(2) with `KD0` in the source SLD fit (see KD0 notes). Bundle simulated-dataset typical values: ~0.02 (bundle subject 1).",
      source_name        = "KD1"
    ),
    IBASE = list(
      description        = "Subject-specific baseline tumour-size estimate carried over (as an empirical-Bayes posterior) from the upstream Zecchin 2016 SLD model. Sets the SLD ODE initial state (`tumor_size(0) = IBASE * 1000`, in mm) and the denominator of the time-varying tumour-size ratio (`mmbas = IBASE * 1000`).",
      units              = "metres (source NONMEM convention; the `*1000` conversion to mm is applied inside the model)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Distinct from the canonical `TUM_SLD` covariate, which carries the *measured* baseline SLD in mm. `IBASE` is the *fitted* baseline from the upstream IPP run and is the value the SLD ODE integrates from. Bundle simulated-dataset typical values: ~0.04-0.50 m (bundle subject 1: 0.285 m = 285 mm).",
      source_name        = "IBASE"
    ),
    AUC_CARBO = list(
      description        = "Per-cycle average AUC of carboplatin (time-varying drug-exposure covariate driving the carboplatin cytotoxic-death term inside the inline SLD ODE).",
      units              = "carboplatin AUC units (mg*min/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Held step-wise constant within each chemotherapy cycle and reset at the start of the next cycle. Set to 0 in cycles where carboplatin is not administered. The DDMORE bundle's Simulated_OS.csv encodes this column as `AUC0`; the source `$INPUT` NM-TRAN column is `CB`.",
      source_name        = "CB"
    ),
    AUC_GEM = list(
      description        = "Per-cycle average AUC of gemcitabine (parent plus active intracellular metabolite per Zecchin 2016 Methods); time-varying drug-exposure covariate driving the gemcitabine cytotoxic-death term inside the inline SLD ODE.",
      units              = "gemcitabine AUC units (paper composite mol*day / 10^6 cells)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Held step-wise constant within each chemotherapy cycle and reset at the start of the next cycle. Set to 0 in carboplatin-monotherapy cycles. The DDMORE bundle's Simulated_OS.csv encodes this column as `AUC1`; the source `$INPUT` NM-TRAN column is `G`.",
      source_name        = "G"
    ),
    NWLS = list(
      description        = "Time-varying binary indicator of whether a new (non-target) RECIST lesion has appeared since enrolment. 1 = new lesion present at the current observation time; 0 = no new lesion as of the current time. Once `NWLS` flips to 1 it stays 1 for subsequent times in that subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no new lesion appeared as of the current time)",
      notes              = "Step-function flag, not a transient pulse. The bundle's Simulated_OS.csv encodes this column as `NWLSCOV`; the source `$INPUT` column is `NWLS`. Hazard effect: multiplicative `exp(e_nwls_haz * NWLS)` with `e_nwls_haz = 1.23` (Output_real_OS.lst FINAL TH5 / Zecchin 2016 Table 2 gamma_NewLes(t) = 1.23).",
      source_name        = "NWLS"
    ),
    TUM_SLD = list(
      description        = "Measured baseline sum of longest diameters (SLD) at enrolment per RECIST 1.1.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Used in the OS hazard via the time-fixed normalised covariate `NSLD0 = TUM_SLD / 70`. The reference value 70 mm is the average baseline SLD in the Zecchin 2016 cohort (Methods / Table 2). Distinct from `IBASE`, which is the empirical-Bayes posterior of the *fitted* SLD baseline from the upstream IPP run (the two values differ at the per-subject level). Source dataset column: `SLD0`.",
      source_name        = "SLD0"
    ),
    ECOG_GE1 = list(
      description        = "Binary baseline Eastern Cooperative Oncology Group (ECOG) performance-status indicator: 1 if ECOG >= 1 at enrolment, 0 if ECOG = 0.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG = 0; fully active / asymptomatic)",
      notes              = "Time-fixed per subject. The Zecchin 2016 paper dichotomizes ECOG to 0 vs >=1 at enrolment because of the small number of patients with ECOG > 1 in the trial cohort (Zecchin 2016 Methods). The bundle's Simulated_OS.csv ships an already-binarized `ECOG` column with values in {0, 1}, so the model's `ECOG_GE1` covariate maps directly from the bundle's `ECOG` column. Hazard effect: multiplicative `exp(e_ecog_haz * ECOG_GE1)` with `e_ecog_haz = 0.516` (Output_real_OS.lst FINAL TH6 / Zecchin 2016 Table 2 gamma_ECOG = 0.518).",
      source_name        = "ECOG"
    )
  )

  population <- list(
    n_subjects     = 336L,
    n_studies      = 1L,
    age_range      = "median ~59 years (advanced epithelial ovarian cancer cohort; Zecchin 2016 Table S2 / paper text)",
    weight_range   = "not transcribed in this extraction (Zecchin 2016 Table S2 captures the demographic distributions; the WebFetch summary did not report the weight quantiles)",
    sex_female_pct = 100,
    disease_state  = "advanced (FIGO stage III/IV) epithelial ovarian cancer (recurrent / platinum-sensitive cohort, randomised Phase III chemotherapy trial)",
    dose_range     = "Phase III chemotherapy: carboplatin monotherapy (target AUC 5.0 mg*min/mL Q3W) or carboplatin (target AUC 4.0 mg*min/mL Q3W) plus gemcitabine, per the trial protocol referenced by Zecchin 2016",
    notes          = "336 patients pooled from a randomised Phase III trial in advanced epithelial ovarian cancer (Zecchin 2016, BJCP 82(3):717-727; PMID 27136318). The current OS model was fit using the IPP (Iterative Population PK) approach: the upstream SLD model (DDMODEL00000217) supplies subject-level empirical-Bayes posteriors of KG, KD0, KD1, IBASE, which feed into the OS model via the dataset. Median baseline SLD ~70 mm (used as the reference TVSLD0). The publication PDF was not on disk for this extraction; the Methods / Table 2 cross-check was performed via PMC HTML (PMC5338128)."
  )

  ini({
    # Final estimates from Output_real_OS.lst FINAL PARAMETER ESTIMATE block
    # (LAPLACIAN CONDITIONAL ESTIMATION, MINIMIZATION SUCCESSFUL).
    # Cross-checked against Zecchin 2016 Table 2 (PMC5338128). The .lst
    # carries lambda in units of 1/day; Table 2 reports lambda = 0.036 in
    # 1/month (paper convention), which converts to ~0.001183/day via the
    # 30.4375-day average month — matching the .lst FINAL TH1.
    llam_haz  <- log(0.001183); label("Weibull baseline-hazard scale parameter, lambda_OS (1/day)")  # Output_real_OS.lst FINAL TH1 = 1.18E-03; Zecchin 2016 Table 2 lambda_OS = 0.036 (1/month) = 0.001183 (1/day)
    lalfa_haz <- log(2.004);    label("Weibull baseline-hazard shape parameter, alpha_OS (unitless)") # Output_real_OS.lst FINAL TH2 = 2.00E+00; Zecchin 2016 Table 2 alpha_OS = 1.99
    e_sld0_haz <- 0.2084;       label("Hazard log-linear coefficient on baseline SLD/70 (unitless)")  # Output_real_OS.lst FINAL TH3 = 2.08E-01; Zecchin 2016 Table 2 gamma_SLD0 = 0.189 (small rounding)
    e_tsr_haz  <- 0.9036;       label("Hazard log-linear coefficient on tumour-size ratio TSR(t) capped at week 12 (unitless)") # Output_real_OS.lst FINAL TH4 = 9.04E-01; Zecchin 2016 Table 2 gamma_TSR = 0.893
    e_nwls_haz <- 1.230;        label("Hazard log-linear coefficient on new-lesion indicator NWLS (unitless)") # Output_real_OS.lst FINAL TH5 = 1.23E+00; Zecchin 2016 Table 2 gamma_NewLes = 1.23
    e_ecog_haz <- 0.5162;       label("Hazard log-linear coefficient on binary ECOG_GE1 indicator (unitless)") # Output_real_OS.lst FINAL TH6 = 5.16E-01; Zecchin 2016 Table 2 gamma_ECOG = 0.518

    # IIV. Source NONMEM run wires ETA(1) on LAM with $OMEGA 0 FIX
    # (Output_real_OS.lst FINAL OMEGA(1,1) = 0.00E+00, ETAshrink 100%);
    # the OMEGA was a placeholder and the OS model has no estimated IIV.
    # No eta* parameters are added.

    # No residual error. The source NONMEM `$ERROR` block sets
    # `Y = SUR` (censored) or `Y = SUR * HAZN` (event) and uses the
    # F_FLAG/LIKE branch, i.e., the likelihood is the survival /
    # event-density itself rather than an observation-error model. This
    # nlmixr2 translation is intended for forward simulation: `sur` and
    # `hazard` are exposed as derived outputs.
  })
  model({
    # Back-transformed Weibull parameters
    lam_haz  <- exp(llam_haz)
    alfa_haz <- exp(lalfa_haz)

    # Fitted-baseline SLD in mm and the time-fixed normalised baseline-SLD
    # covariate term. TVSLD0 = 70 mm is the cohort-average reference
    # (Zecchin 2016 Methods / source $PK block).
    mmbas <- IBASE * 1000
    nsld0 <- TUM_SLD / 70

    # Inline SLD ODE — direct port of $DES from Executable_OS.mod with the
    # source `/1000` and `/100` numerical scalings preserved verbatim.
    # KG, KD0, KD1 are subject-level empirical-Bayes posteriors from the
    # upstream Zecchin 2016 SLD model (DDMODEL00000217); AUC_CARBO and
    # AUC_GEM are time-varying drug-exposure covariates.
    dtum_dt <- KG / 1000 * tumor_size -
      (KD0 / 1000 * AUC_CARBO + KD1 / 100 * AUC_GEM) * tumor_size
    d/dt(tumor_size) <- dtum_dt
    tumor_size(0) <- IBASE * 1000  # mm

    # Tumour-size ratio and the time-varying covariate fed into the
    # hazard. Source $DES freezes WTS at the week-12 (84-day) value, so
    # the wts state integrates d(TSR)/dt for t <= 84 and 0 thereafter.
    # Since TSR(0) = (mmbas - mmbas)/mmbas = 0, wts(0) = 0 reproduces the
    # source `IF(T.EQ.0) WTS = 0` initialisation exactly.
    tsr <- (tumor_size - mmbas) / mmbas
    d/dt(wts) <- ifelse(t <= 84, dtum_dt / mmbas, 0)
    wts(0) <- 0

    # Weibull baseline hazard with covariate-effect log-linear modulation.
    # The DEL = 1e-6 small-time offset matches the source $DES / $ERROR
    # blocks and keeps the (lam_haz * t)^(alfa_haz - 1) term finite at
    # t = 0 without affecting the integrated cumulative hazard.
    del <- 1e-6
    hazard <- lam_haz * alfa_haz * (lam_haz * (t + del))^(alfa_haz - 1) *
      exp(e_sld0_haz * nsld0 + e_tsr_haz * wts +
          e_nwls_haz * NWLS + e_ecog_haz * ECOG_GE1)

    # Cumulative hazard and survival probability
    d/dt(cumHazard) <- hazard
    cumHazard(0) <- 0
    sur <- exp(-cumHazard)
  })
}
