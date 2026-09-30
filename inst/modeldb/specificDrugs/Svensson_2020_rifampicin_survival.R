Svensson_2020_rifampicin_survival <- function() {
  description <- paste(
    "Parametric time-to-event model for 6-month all-cause mortality in",
    "adults with tuberculous meningitis treated with standard- or high-dose",
    "rifampicin in three Indonesian phase 2 trials (Svensson 2020). The",
    "hazard declines exponentially from a base hazard of 0.0286 1/day at a",
    "rate of 0.0333 1/day, rises linearly as baseline Glasgow Coma Scale",
    "falls below 13 and as a power (1.04) of age relative to 30 years, and",
    "is reduced by an Emax function of the individual day-2 rifampicin",
    "plasma AUC0-24 with maximal effect fixed to a 100% reduction and AUC50",
    "171 mg*h/L. The AUC covariate is the individual day 2 +/- 1 plasma",
    "AUC0-24 predicted by the companion PK model Svensson_2020_rifampicin."
  )
  reference <- paste(
    "Svensson EM, Dian S, te Brake L, Ganiem AR, Yunivita V, van Laarhoven A,",
    "van Crevel R, Ruslami R, Aarnoutse RE (2020).",
    "Model-Based Meta-analysis of Rifampicin Exposure and Mortality in",
    "Indonesian Tuberculous Meningitis Trials.",
    "Clin Infect Dis 71(8):1817-1823. doi:10.1093/cid/ciz1071.",
    "Parameter estimates and hazard equation from Table 2 and its footnote;",
    "cross-checked against the NONMEM control stream in the Online Data",
    "Supplement ('NONMEM code survival model').",
    sep = " "
  )
  vignette <- "Svensson_2020_rifampicin"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    cumhaz = list(
      analyte = "Cumulative hazard of death",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age at study entry",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the hazard, (AGE/30)^1.04, centred on the cohort median of 30 years (Table 1; range 16-81).",
      source_name = "AGE"
    ),
    SCORE_GCS = list(
      description = "Glasgow Coma Scale total score at baseline (3-15; lower = more severe impairment of consciousness)",
      units = "(score points, 3-15)",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on the hazard, (1 - 0.256 * (SCORE_GCS - 13)), centred on the cohort median of 13 (Table 1; range 3-15). Time-fixed baseline value. Source column GCSB.",
      source_name = "GCSB"
    ),
    AUC_RIF = list(
      description = "Individual rifampicin plasma AUC0-24 on day 2 +/- 1 of study treatment",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical Bayes AUC0-24 at the first PK occasion predicted by the",
        "final PK model (Svensson_2020_rifampicin, OCC = 1). Time-fixed per",
        "patient. Enters as EFF = 1 - AUC_RIF / (171 + AUC_RIF). For the 15",
        "patients without PK data the control stream imputes typical values",
        "by dose: 47.61 (450 mg oral), 119.1 (600 mg intravenous), 123.1",
        "(750 mg oral), 183 (900 mg oral) and 233 mg*h/L (1350 mg oral).",
        "Source column AUC."
      ),
      source_name = "AUC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Not significant in univariate analysis (P >= .05; Results, Survival Model). 45% female (Table 1)."
    ),
    HIV_POS = list(
      description = "HIV infection (1 = positive)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (HIV-negative)",
      notes = "Not significant in univariate analysis; 12% of patients were HIV-infected (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not significant in univariate analysis. Median 46 kg, range 34-78 (Table 1)."
    ),
    CSF_TPRO = list(
      description = "Cerebrospinal-fluid total protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Not significant on the hazard (Results, Survival Model). Median 165 mg/dL = 1.65 g/L (Table 1). CSF neutrophil count was also screened and not retained; the length of intensified treatment (14 or 30 days) was tested and not significant (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 148L,
    n_studies = 3L,
    n_events = "58 deaths and 15 dropouts within 6 months",
    age_range = "16-81 years (median 30)",
    age_median = "30 years",
    weight_range = "34-78 kg (median 46)",
    weight_median = "46 kg",
    sex_female_pct = 45,
    race_ethnicity = "Indonesian (all three trials were conducted in Bandung, Indonesia)",
    disease_state = paste(
      "Adults with definite (56%), probable or possible tuberculous",
      "meningitis; 95% grade 2 or 3; baseline Glasgow Coma Scale median 13",
      "(3-15); 12% HIV-infected. All received adjunctive dexamethasone."
    ),
    dose_range = "Rifampicin once daily orally 450, 750, 900 or 1350 mg, or 600 mg as a 1.5 h intravenous infusion, for the first 14 or 30 days",
    regions = "Indonesia (Bandung)",
    notes = "Individual-patient pooled analysis of three randomized phase 2 trials (Ruslami 2013, Yunivita 2016, Dian 2018) with 6-month follow-up of survival. Demographics are Table 1."
  )

  ini({
    # Table 2 final estimates; control stream $THETA carries more digits.
    lbase <- log(0.0285801)
    label("Base hazard of death at time zero for the typical patient without rifampicin effect (1/day)")  # Table 2 BASE 0.0286 1/day (RSE 41%); $THETA (0,0.0285801) ; 1 Base hazard [per day]
    lkdec <- log(0.0333392)
    label("Rate constant of the exponential decline in hazard (1/day)")                                    # Table 2 k 0.0333 1/day (17%); $THETA (0,0.0333392,5) ; 2 coefficient exponential decline
    e_gcs_haz <- -0.255676
    label("Linear effect of baseline Glasgow Coma Scale on the hazard, per point above 13 (1/point)")      # Table 2 theta GCS -0.256 (28%); $THETA (-0.5,-0.255676,5) ; 3 GCSB effect linear
    e_age_haz <- 1.0438
    label("Power exponent of age / 30 years on the hazard (unitless)")                                     # Table 2 theta age 1.04 (39%); $THETA (-1,1.0438,5) ; 4 AGE effect power
    lec50_auc <- log(171.358)
    label("Day-2 rifampicin plasma AUC0-24 giving half of the maximal (100%) hazard reduction (mg*h/L)")   # Table 2 theta RIF EC50 171 mg/L*h (86%); $THETA (0,171.358) ; 5 RIF EC50 effect

    # No random effects: the control stream carries $OMEGA 0 FIX ('the ETA is
    # a placeholder here'), so none is added.
    # The source likelihood is the survival / event density (LIKE); there
    # is no residual error. A tiny placeholder additive error on `sur` lets
    # the nlmixr2 machinery accept the forward-simulation model.
    addSd_sur <- 0.001
    label("Placeholder additive residual error on the survival probability (unitless); not from the source")
  })

  model({
    base <- exp(lbase)
    kdec <- exp(lkdec)
    ec50_auc <- exp(lec50_auc)

    # Table 2 footnote a and control stream $DES:
    # h(t) = BASE*exp(-k*t)*(1 + thGCS*(GCS - 13))*(AGE/30)^thAGE*
    #        (1 - AUC/(EC50 + AUC))
    eff_rif <- 1 - AUC_RIF / (ec50_auc + AUC_RIF)
    hazard <- base * exp(-kdec * t) * (1 + e_gcs_haz * (SCORE_GCS - 13)) *
      (AGE / 30)^e_age_haz * eff_rif

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    sur ~ add(addSd_sur)
  })
}
