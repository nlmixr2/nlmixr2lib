Lin_2020_glasdegib_exposure <- function() {
  description <- "Parametric exponential time-to-event model for overall survival in the exposure-response analysis subset of BRIGHT AML 1003 Phase 2 (Lin 2020): older adults with newly diagnosed AML ineligible for intensive chemotherapy who were randomised to glasdegib 100 mg QD + LDAC and had glasdegib PK data available (n = 75). Hazard h(t) = lambda. Seven glasdegib exposure metrics (first-dose Cmax, end-of-cycle-1 Cmax, end-of-cycle-1 Cmin, cycle-1 cumulative AUC, cycle-1 Cavg, average AUC over dosing interval, overall Cavg; raw scale and log-transformed) and baseline ECOG performance status and cytogenetic risk were tested in SCM forward inclusion (alpha < 0.05); ECOG and cytogenetic risk entered the full model but no covariates survived backward elimination (alpha < 0.001). The final exposure-response model is therefore an intercept-only exponential survival on the glasdegib + LDAC subset."
  reference <- paste(
    "Lin S, Shaik N, Chan G, Cortes JE, Ruiz-Garcia A.",
    "An evaluation of overall survival in patients with newly diagnosed",
    "acute myeloid leukemia and the relationship with glasdegib treatment",
    "and exposure.",
    "Cancer Chemother Pharmacol. 2020;86(4):451-459.",
    "doi:10.1007/s00280-020-04132-x.",
    "BRIGHT AML 1003: ClinicalTrials.gov NCT01546038.",
    sep = " "
  )
  vignette <- "Lin_2020_glasdegib_AML_overall_survival"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; all subjects are in the glasdegib + LDAC arm and exposure metrics were not retained as covariates)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    CMAX_FIRST_DOSE = list(
      description        = "First-dose maximum concentration (ng/mL) of glasdegib at the 100 mg QD clinical dose.",
      units              = "ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion on both the raw scale and after natural-log transformation; did not reach significance (alpha < 0.05). Lin 2020 Table 2 cohort summary: mean 592.0, median 602.2 (range 0.8-1437.8), n = 75."
    ),
    CMAX_EOC1 = list(
      description        = "End-of-cycle-1 maximum concentration (ng/mL) of glasdegib.",
      units              = "ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Cycle definition = 28 days. Lin 2020 Table 2 cohort summary: mean 1308.3, median 1139.9 (range 275.3-3612.5)."
    ),
    CMIN_EOC1 = list(
      description        = "End-of-cycle-1 minimum concentration (ng/mL) of glasdegib.",
      units              = "ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Lin 2020 Table 2 cohort summary: mean 750.7, median 577.6 (range 4.9-2762.0)."
    ),
    AUC_CYCLE1 = list(
      description        = "Cycle-1 cumulative AUC (h*ng/mL) of glasdegib.",
      units              = "h*ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Lin 2020 Table 2 cohort summary: mean 631571.2, median 575220.0 (range 54368.0-1797600.0)."
    ),
    CAVG_CYCLE1 = list(
      description        = "Cycle-1 average concentration (ng/mL) of glasdegib (= cycle-1 AUC / time).",
      units              = "ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Lin 2020 Table 2 cohort summary: mean 1009.7, median 927.1 (range 161.9-3505.4)."
    ),
    AUC_TAU_AVG = list(
      description        = "Average AUC over the dosing interval (h*ng/mL) of glasdegib.",
      units              = "h*ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Lin 2020 Table 2 cohort summary: mean 18580.2, median 16571.3 (range 2927.9-61527.3)."
    ),
    CAVG_OVERALL = list(
      description        = "Overall average concentration (ng/mL) of glasdegib across the treatment duration.",
      units              = "ng/mL",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (raw and log scale); not significant. Lin 2020 Table 2 cohort summary: mean 910.1, median 783.8 (range 218.5-3505.4)."
    ),
    ECOG_GE1 = list(
      description        = "Binary baseline ECOG performance status indicator (1 = ECOG >= 1, 0 = ECOG = 0).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Entered the SCM full model on the exposure-response cohort during forward inclusion (alpha < 0.05) but did not survive backward elimination (alpha < 0.001)."
    ),
    CYTO_POOR = list(
      description        = "Poor cytogenetic risk indicator (1 = poor, 0 = good/intermediate).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Entered the SCM full model on the exposure-response cohort during forward inclusion (alpha < 0.05) but did not survive backward elimination (alpha < 0.001)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 75L,
    n_studies      = 1L,
    age_range      = "Subset of BRIGHT AML 1003 Phase 2 glasdegib + LDAC arm (n = 78) restricted to patients with glasdegib PK information available (n = 75); per-subset age distribution not reported separately. Cohort-wide median 76 years, range 58-92 (Lin 2020 Table 1)",
    weight_range   = "Per-subset distribution not reported; cohort-wide median 78.2 kg, range 47.5-118.0",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "newly diagnosed, previously untreated acute myeloid leukemia in adults ineligible for intensive chemotherapy",
    dose_range     = "Glasdegib 100 mg orally once daily in 28-day cycles + LDAC 20 mg SC BID for 10 days per 28-day cycle; dose reductions allowed for adverse-event management (14/75 subjects had >=1 dose reduction; 13/75 due to treatment-related AEs)",
    regions        = "Multicenter (NCT01546038)",
    notes          = "75 of the 78 Phase 2 glasdegib + LDAC subjects with available glasdegib PK. Per-subject glasdegib exposure metrics derived from empirical Bayes estimates of an upstream population PK model (Lin 2020 reference [16]); the exposure metrics did not impact OS (no metric reached alpha < 0.05 in SCM forward inclusion). Final TTE model = base model = exponential with no covariates."
  )

  ini({
    # Exponential time-to-event hazard for overall survival in the
    # exposure-response analysis subset (BRIGHT AML 1003 Phase 2, n = 75).
    # Final estimates from Lin 2020 page 5: "The estimated base hazard
    # (RSE) was 0.00246 (14.19%) and the survival probability was
    # described by: S(t) = exp(-0.00246 * t)." Fig. 3 annotation
    # cross-checks the same expression.
    llam_haz <- log(0.00246); label("Exponential baseline hazard for overall survival in the glasdegib + LDAC exposure-response cohort (lambda, log(1/day))")  # Lin 2020 Results page 5: lambda = 0.00246, RSE 14.19%; matches Fig. 3 annotation S(t) = exp(-0.00246 t)

    # No estimated IIV (parametric population TTE).
    # No residual error (NONMEM LIKE / F_FLAG branch).
    # No covariates were retained after SCM backward elimination.
  })
  model({
    # Constant exponential baseline hazard (intercept-only model).
    lam_haz <- exp(llam_haz)
    hazard  <- lam_haz

    # Cumulative hazard and survival probability.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
