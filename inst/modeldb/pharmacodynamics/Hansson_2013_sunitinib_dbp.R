Hansson_2013_sunitinib_dbp <- function() {
  description <- "Indirect-response model of sunitinib-induced increase in diastolic blood pressure (dBP) in adults with imatinib-resistant gastrointestinal stromal tumours (GIST). The dBP state turns over via a stimulated zero-order production rate Kin and a first-order removal rate Kout (= 1 / MRT), where the per-cycle drug-exposure summary AUC = DOSE / CLI linearly stimulates Kin via a slope factor dBP_slope. Kin is parameterised as dBP0 * Kout so the dBP steady state without drug equals dBP0. A separate higher baseline dBP0_placebo is recorded as a placebo-arm typical value (Results: 'this group had a significantly higher baseline dBP (dBP0) when estimated separately'). The PD model has no PK ODE; sunitinib exposure enters as the AUC summary computed from time-varying DOSE and per-subject CLI."
  reference <- paste(
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO.",
    "PKPD modeling of predictors for adverse effects and overall survival",
    "in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "Sister model files from the same paper:",
    "modellib('Hansson_2013_sunitinib_myelosuppression'),",
    "modellib('Hansson_2013c_sunitinib') [fatigue],",
    "modellib('Hansson_2013_sunitinib_hfs'),",
    "modellib('Hansson_2013_sunitinib_os').",
    "Indirect-response form adapted from Keizer RJ et al.",
    "J Pharmacokinet Pharmacodyn 2010;37(4):347-363,",
    "doi:10.1007/s10928-010-9163-3.",
    sep = " "
  )
  vignette <- "Hansson_2013_sunitinib_dbp"
  units <- list(time = "hour", dosing = "mg", concentration = "mmHg (diastolic blood pressure)")

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper Methods describes sunitinib administered at 25-75 mg PO QD on 4/2, 2/2, 2/1, or continuous schedules across studies 1004, 1047, 1045, and 013 (Table 1). For typical-cohort vignette simulations the value is held at 50 mg during on-cycles of a 4/2 schedule and 0 mg during off-cycles, matching the largest cohort (Demetri 2006 / Study 1004).",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e85 Methods describes the upstream popPK as a previously developed 2-compartment model (Houk et al. 2009 Clin Cancer Res 15:2497-2506; that popPK is not packaged in nlmixr2lib at extraction time). The companion Hansson_2013a / Hansson_2013c sunitinib model files use a typical-value reference of 32.819 L/h.",
      source_name        = "CL"
    ),
    PLACEBO = list(
      description        = "Binary indicator: 1 = placebo-arm subject (Study 1004 placebo run-in n=47), 0 = active sunitinib arm. Switches the typical baseline dBP0 from 71.8 mmHg (active) to 77.6 mmHg (placebo) per Hansson 2013 Table 2 row 'dBP0 placebo'.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (active sunitinib arm)",
      notes              = "The paper Results: 'No increase in dBP could be identified for placebo patients. However, this group had a significantly higher baseline dBP (dBP0) when estimated separately.' Placebo-arm subjects also have DOSE = 0 (and therefore AUC = 0 and no drug-effect term). For typical-cohort simulations set every subject to 0 (active).",
      source_name        = "(derived from treatment arm; Study 1004 placebo run-in subjects only)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST (paper text reports n=303 pooled across phases I-III; per-cohort baseline-demographics table not in the trimmed PDF section that includes Methods + Results + Tables)",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies (Demetri 2006 study 1004, George 2009 study 1047, Shirao 2010 study 1045, Maki 2005 study 013).",
    dose_range     = "sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 or continuous schedule (Table 1). Placebo arm: no sunitinib.",
    regions        = "phase III multinational (study 1004); Japanese phase I/II (study 1045); other studies' regions not stated in the trimmed paper text.",
    biomarkers     = "Diastolic blood pressure (dBP) measured serially during treatment cycles. Median (range) observed dBP during treatment: 80 (20-120) mmHg in study 1004, 78 (40-120) in study 1047, 79 (40-120) in study 1045, 80 (50-130) in study 013 (Hansson 2013 Table 1).",
    notes          = "n_subjects = 303 reported in Hansson 2013 e85 Methods. The dBP model was fit with Methods-specified blood pressure observations assumed to occur in the morning ('The actual times of the day for blood pressure measurements were not available and were therefore assumed to occur in the morning for all observations')."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PD parameters from Hansson 2013 e85 Table 2 'Blood pressure
    # model' block. The paper Table 2 reports point estimates with relative
    # standard error (RSE %) for the typical population; IIV is reported as
    # CV% which back-transforms to omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------

    # Baseline dBP for actively treated subjects.
    ldbp0          <- log(71.8); label("Baseline dBP for sunitinib-treated subjects (mmHg; Hansson 2013 Table 2 row 'dBP0')") # Table 2 dBP0 = 71.8 mmHg (RSE 1.0%)

    # Multiplicative shift for placebo subjects: dBP0_placebo / dBP0 = 77.6 / 71.8 = 1.0808.
    # Encoded as an exponential effect on log(dBP0): log(77.6/71.8) = 0.0777.
    e_placebo_dbp0 <- log(77.6 / 71.8); label("Multiplicative shift on log(dBP0) for placebo-arm subjects yielding dBP0_placebo = 77.6 mmHg (Hansson 2013 Table 2 row 'dBP0 placebo')") # Table 2 dBP0 placebo = 77.6 mmHg (RSE 1.6%)

    # Mean residence time (= 1 / kout) in hours.
    lmrt           <- log(361);  label("Mean residence time MRT = 1/kout for dBP (h; Hansson 2013 Table 2 row 'MRT (= 1/kout)')") # Table 2 MRT = 361 h (RSE 17%)

    # Linear slope from AUC to the production-rate stimulation factor.
    ldbp_slope     <- log(0.119); label("Linear slope factor dBP_slope linking sunitinib AUC to the Kin stimulation (L/(mg*h); Hansson 2013 Table 2 row 'dBP slope')") # Table 2 dBP slope = 0.119 L/mg.h (RSE 9.4%)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports IIV as CV%; conversion
    # to log-normal omega^2 = log((CV/100)^2 + 1).
    # dBP0       IIV = 12% -> omega2 = log(0.12^2 + 1) = 0.01435
    # MRT        IIV = 83% -> omega2 = log(0.83^2 + 1) = 0.5328
    # dBP_slope  IIV = 65% -> omega2 = log(0.65^2 + 1) = 0.3578
    # dBP0 and dBP_slope are reported as correlated 65% (Results section).
    # cov(dBP0, dBP_slope) = 0.65 * sqrt(0.01435 * 0.3578) = 0.04657
    # ------------------------------------------------------------------
    etaldbp0 + etaldbp_slope ~ c(0.01435,
                                 0.04657, 0.3578) # Table 2 dBP0 CV=12%, dBP_slope CV=65%, corr=0.65
    etalmrt ~ 0.5328  # Table 2 MRT CV=83%

    # ------------------------------------------------------------------
    # Residual error from Hansson 2013 Table 2: combined additive
    # (6.24 mmHg) and proportional (6.97%) error.
    # ------------------------------------------------------------------
    addSd_dbp  <- 6.24;   label("Additive residual error on dBP (mmHg; Hansson 2013 Table 2 row 'Residual error (mmHg)')") # Table 2 = 6.24 mmHg (RSE 16%)
    propSd_dbp <- 0.0697; label("Proportional residual error on dBP (fraction; Hansson 2013 Table 2 row 'Residual error (%)')") # Table 2 = 6.97% (RSE 24%)
  })

  model({
    # ---- 1. Per-cycle drug-exposure summary (mg*h/L AUC) ----
    # mg / (L/h) = mg*h/L; matches the Hansson_2013_sunitinib_myelosuppression
    # and Hansson_2013c_sunitinib convention.
    auc <- DOSE / CLI

    # ---- 2. Individual structural parameters ----
    # dBP0 covaries between treated (71.8 mmHg) and placebo (77.6 mmHg) arms
    # via the binary PLACEBO indicator. Methods/Results: the paper estimated
    # the placebo dBP0 as a separately fitted typical value because placebo
    # patients had higher baseline dBP. PLACEBO = 1 -> dBP0 = 77.6;
    # PLACEBO = 0 -> dBP0 = 71.8.
    dbp0      <- exp(ldbp0 + e_placebo_dbp0 * PLACEBO + etaldbp0)
    mrt       <- exp(lmrt + etalmrt)
    dbp_slope <- exp(ldbp_slope + etaldbp_slope)

    # ---- 3. Indirect-response rate constants ----
    # kout = 1 / MRT; kin is derived as dBP0 * kout so that the drug-free
    # steady state of dBP equals dBP0 (Methods: 'dBP0 and kout ... were
    # estimated, and Kin was derived as dBP0 x kout').
    kout <- 1 / mrt
    kin  <- dbp0 * kout

    # ---- 4. Drug effect: linear stimulation of Kin ----
    # Methods: 'Sunitinib AUC was linked to the production rate (Kin) by a
    # linear slope factor (dBP_Drug effect)'. The standard indirect-response
    # encoding is Kin * (1 + dBP_slope * AUC).
    drug_effect <- 1 + dbp_slope * auc

    # ---- 5. ODE ----
    d/dt(dbp) <- kin * drug_effect - kout * dbp
    dbp(0)    <- dbp0

    # ---- 6. Observation ----
    # Combined additive + proportional residual error (Table 2).
    dbp ~ add(addSd_dbp) + prop(propSd_dbp)
  })
}
