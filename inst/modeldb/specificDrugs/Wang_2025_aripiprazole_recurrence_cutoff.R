Wang_2025_aripiprazole_recurrence_cutoff <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for the time to recurrence",
    "of any mood episode in adults with bipolar I disorder maintained on",
    "aripiprazole once monthly (AOM) 400 mg or withdrawn to placebo in the",
    "randomised double-blind withdrawal phase of trial NCT01567527 (Wang 2025,",
    "categorical model). This is the separately fitted model in which the",
    "exposure driver CTROUGH -- the model-predicted aripiprazole plasma",
    "concentration 672 h after the first double-blind injection ('Ctau'), from",
    "the companion population PK model modellib('Wang_2022_aripiprazole') --",
    "enters only through the 95 ng/mL efficacy threshold previously established",
    "in schizophrenia: subjects at or above the threshold have a hazard of",
    "recurrence 0.6450-fold that of subjects below it (a 36 percent, or",
    "1.55-fold, reduction). The authors' preferred final model is the",
    "continuous one, modellib('Wang_2025_aripiprazole_recurrence'); this",
    "categorical model was fitted to allow direct comparison with the",
    "schizophrenia threshold model modellib('Wang_2022_aripiprazole_relapse').",
    "Because the fit is a semiparametric Cox regression, the baseline hazard is",
    "unspecified by construction and is not reported; the model therefore",
    "returns the RELATIVE hazard hr against the below-threshold reference group,",
    "not an absolute survivor function. The model is algebraic and",
    "deterministic (no ODE state, no drug input, no IIV, no residual error)."
  )
  reference <- paste(
    "Wang X, Phillips L, Harlin M, Bell Lynum KS, Larsen F, Such P,",
    "Madera-McDonough J, Yildirim M, Procyshyn RM, Chepke C, Passarell J.",
    "Exposure-Response Analysis for Aripiprazole Once-Monthly in Patients",
    "Diagnosed With Bipolar I Disorder.",
    "Clin Pharmacol Drug Dev. 2025;14(11):869-878. doi:10.1002/cpdd.1580.",
    "The hazard ratio is from the Results section 'Exposure-Response Model'.",
    "The upstream population PK model that generates the CTROUGH covariate is",
    "Wang X, Raoufinia A, Bihorel S, et al. Clin Pharmacol Drug Dev.",
    "2022;11(2):150-164, doi:10.1002/cpdd.1022; see",
    "modellib('Wang_2022_aripiprazole').",
    sep = " "
  )
  vignette <- "Wang_2025_aripiprazole_bipolar"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method; the source analysis measured time in days from randomisation)",
    dosing = "n/a (exposure-response model; aripiprazole exposure enters as the CTROUGH covariate, not as a dose record)",
    concentration = "hr (hazard of mood-episode recurrence relative to CTROUGH < 95 ng/mL, unitless); not a drug concentration"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted aripiprazole plasma concentration 672 h (28 days) after the first injection of the double-blind withdrawal phase (Wang 2025 'Ctau'), i.e. the end-of-interval trough of the first double-blind dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "(a) TOTAL aripiprazole plasma concentration (parent drug only).",
        "(b) Landmark: 672 h after the first double-blind (Phase D) injection;",
        "one time-fixed value per subject over the whole follow-up.",
        "(c) Empirical-Bayes prediction from the Wang 2022 population PK model",
        "(modellib('Wang_2022_aripiprazole'), all parameters fixed) applied to",
        "each subject's recorded dosing history (Wang 2025 Methods, Appendix S3).",
        "Placebo records are NOT zero: they carry the residual concentration",
        "washing out from the last open-label AOM dose, about 1344 h earlier.",
        "The model consumes this column ONLY through the indicator",
        "(CTROUGH >= ctrough_cut), with ctrough_cut = 95 ng/mL."
      ),
      source_name = "Ctau"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 265L,
    n_studies = 1L,
    n_events = 103L,
    age_range = "18-65 years by inclusion criterion; mean (SD) 40.1 (11.0) years in the exposure-response set (Wang 2025 Table S5)",
    weight_range = "Body mass index mean (SD) 30.2 (7.2) kg/m^2 in the exposure-response set (Table S5)",
    sex_female_pct = 57.4,
    race_ethnicity = c(
      White = 54.3,
      `Black or African American` = 27.9,
      Asian = 14.0,
      `American Indian or Alaska Native` = 0.8,
      Other = 3.0
    ),
    disease_state = "Adults with bipolar I disorder, stabilised on open-label AOM after a manic index episode, randomised to continue AOM or switch to placebo. See modellib('Wang_2025_aripiprazole_recurrence') for the full description of the cohort and the event definition.",
    dose_range = "AOM 400 mg (or 300 mg after a tolerability reduction) gluteal IM every 4 weeks (n = 132) or placebo injection (n = 133) in the 52-week double-blind phase.",
    regions = "103 sites in the USA, Canada, Poland, Romania, Taiwan, South Korea and Japan (Appendix S1, Table S1)",
    notes = "Same 265-subject analysis set, event definition and follow-up as modellib('Wang_2025_aripiprazole_recurrence'). The threshold was tested as a new categorical model (P = .0270); several other cut points were assessed (data not shown) and none described the data better than the continuous model."
  )

  ini({
    # ==================================================================
    # Wang 2025 categorical exposure-response model (Results, 'Exposure-
    # Response Model'):
    #
    #   h_i(t) = h0(t) * exp(beta * I(Ctau_i >= 95 ng/mL))
    #
    # Only the hazard ratio exp(beta) = 0.6450 (95% CI 0.4372-0.9514) is
    # printed, so beta is encoded as log(0.6450). The midpoint of the log
    # CI, (log(0.4372) + log(0.9514)) / 2 = -0.4386, agrees with
    # log(0.6450) = -0.4385, confirming the CI is symmetric on the log
    # scale as a Wald interval for a Cox coefficient should be.
    #
    # NO BASELINE HAZARD IS ENCODED (semiparametric Cox regression; see
    # Wang_2025_aripiprazole_recurrence.R). No IIV and no residual error.
    # ==================================================================

    # Threshold carried over from the schizophrenia exposure-response
    # analysis (Wang 2022) and tested, not estimated, here; hence fixed().
    ctrough_cut <- fixed(95); label("Aripiprazole Ctau threshold separating the reference (below) and above-threshold groups (ng/mL)")  # Wang 2025 Results, 'Exposure-Response Model': 'An aripiprazole Ctau of 95 ng/mL was tested as a new model (categorical)'

    e_ctrough_cut_haz <- log(0.6450); label("Log hazard ratio for recurrence of any mood episode, Ctau >= 95 ng/mL versus < 95 ng/mL (log scale; HR 0.6450)")  # Wang 2025 Results, 'Exposure-Response Model': HR 0.6450 (95% CI 0.4372-0.9514), P = .0270
  })

  model({
    # Boundary convention follows the paper's wording exactly: the
    # above-threshold group is Ctau >= 95 ng/mL; the reference group is
    # Ctau < 95 ng/mL.
    aboveCut <- (CTROUGH >= ctrough_cut)

    lhr <- e_ctrough_cut_haz * aboveCut
    hr <- exp(lhr)
  })
}
