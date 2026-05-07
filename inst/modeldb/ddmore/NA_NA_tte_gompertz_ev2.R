NA_NA_tte_gompertz_ev2 <- function() {
  description <- "Parametric time-to-event Gompertz hazard model for Event 2 in the BAST PTTE 2017 four-event teaching dataset (DDMODEL00000243). Hazard h(t) = (lam/1000) * exp((alpha/1000)*t) * exp((coef_auc/1000)*(AUC_BAST_FW - 3065.5)). Event 2 in the bundle's simulated dataset is interval-censored (CENSORING = 2), assessed at scheduled visits rather than observed exactly; the BAST guiding-document Section  2.4.1 (Figure 2-2) selected Gompertz as the base distribution by AIC, then Section  2.4.2 (Table 2-3) retained first-week AUC as the only covariate."
  reference <- paste(
    "BAST Inc Limited.",
    "BAST approach to parametric time-to-event (PTTE) modelling.",
    "Loughborough, UK; 12 July 2017.",
    "Internal guiding document (BAST_PTTE_modelling.pdf) shipped",
    "with DDMORE bundle DDMODEL00000243; no peer-reviewed publication.",
    "Run prepared by Jon Moss (Command.txt; runEV2_105).",
    "DDMORE Foundation Model Repository: DDMODEL00000243.",
    sep = " "
  )
  vignette <- "NA_NA_tte_gompertz_ev2"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the AUC_BAST_FW covariate is a per-subject baseline first-week-AUC summary)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  ddmore_id    <- "DDMODEL00000243"
  replicate_of <- NULL

  covariateData <- list(
    AUC_BAST_FW = list(
      description        = "First-week AUC of the unspecified hypothetical drug used in the BAST PTTE 2017 teaching dataset.",
      units              = "ug*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Centred at 3065.5 ug*h/L inside the hazard via `exp((coef_auc/1000) * (AUC_BAST_FW - 3065.5))`. The bundle's Simulated_event_data.csv carries the source-named `AUC` column for 200 hypothetical patients with range 858.6-7673.3 ug*h/L (mean 3189); 3065.5 is the cohort median used by the BAST PTTE guiding-document. Renamed from the source-data column `AUC` to the canonical `AUC_BAST_FW` in the register so that a future model using a generic `AUC` column with different drug semantics will not silently collide.",
      source_name        = "AUC"
    )
  )

  population <- list(
    n_subjects     = 200L,
    n_studies      = 1L,
    age_range      = "24-84 years (mean 58.7) in the BAST PTTE 2017 simulated cohort",
    weight_range   = "not reported (the BAST PTTE 2017 simulated cohort does not include body weight)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Hypothetical / unspecified clinical population (the BAST PTTE 2017 guiding document is a methodological teaching example with simulated event data; no real drug, indication, or patient cohort).",
    dose_range     = "Not applicable (no drug administration is modelled; the AUC_BAST_FW covariate is a per-subject baseline summary of the unspecified drug's first-week exposure).",
    regions        = "Not applicable (simulated data).",
    notes          = "200 simulated patients; 104 (52%) had Event 2. Event 2 is interval-censored: exact event times are unknown, only that the event occurred between two scheduled assessment visits (BAST guiding document Section  2.2.1). For the other three events in the same bundle see NA_NA_tte_gompertz.R, NA_NA_tte_lognormal.R, and NA_NA_tte_loglogistic.R. Source: BAST Inc Limited, 'BAST approach to parametric time-to-event (PTTE) modelling', 12 July 2017 (BAST_PTTE_modelling.pdf); base-distribution selection Section  2.4.1 / Figure 2-2; covariate selection Section  2.4.2 / Table 2-3; final-fit listing Output_simulated_runEV2_105.res."
  )

  ini({
    # Final estimates from Output_simulated_runEV2_105.res FINAL PARAMETER
    # ESTIMATE block (LAPLACIAN CONDITIONAL ESTIMATION, MIN SUCCESSFUL,
    # OBJV = 1040443.295, NSIG = 3.7). The .mod uses internal numerical
    # rescalings (lam/1000, alpha/1000, coef_auc/1000); the rescalings
    # are kept inside model() so that ini() values match the .lst FINAL
    # THETA values one-to-one.
    llam_ev2    <- log(3.28); label("Event 2 Gompertz lambda scale parameter (.mod TH1; rescaled to lam/1000 inside model())")            # Output_simulated_runEV2_105.res FINAL TH1 = 3.28E+00
    lalpha_ev2  <- log(4.05); label("Event 2 Gompertz alpha shape parameter, time-decay rate (.mod TH2; rescaled to alpha/1000 inside model())")  # Output_simulated_runEV2_105.res FINAL TH2 = 4.05E+00
    e_auc_ev2   <- 0.309;     label("Event 2 hazard log-linear coefficient on (AUC_BAST_FW - 3065.5) ug*h/L (.mod TH3; rescaled to coef/1000 inside model())")  # Output_simulated_runEV2_105.res FINAL TH3 = 3.09E-01

    # No IIV. Source NONMEM run wires ETA(1) on lam with `$OMEGA 0 FIX`
    # (Output_simulated_runEV2_105.res FINAL OMEGA(1,1) = 0 FIXED). No
    # eta* parameters are added.

    # No residual error. Source $ERROR block uses `$EST ... LIKE` (the
    # likelihood is the survival/event density). Forward-simulation
    # translation: `hazard` and `sur` are exposed as derived outputs.
  })
  model({
    # Numerical rescalings carried verbatim from Executable_runEV2_105.mod
    # $PK / $DES blocks (`Lam1 = Lam/1000`, `alpha1 = alpha/1000`,
    # `EXP((TH3/1000)*(AUC-3065.5))`). Preserved here so ini() values
    # match the .lst FINAL THETA values directly.
    lam     <- exp(llam_ev2)
    alpha   <- exp(lalpha_ev2)
    lam_1   <- lam   / 1000
    alpha_1 <- alpha / 1000

    # Covariate effect on the hazard (centred-deviation form, multiplicative).
    val <- exp((e_auc_ev2 / 1000) * (AUC_BAST_FW - 3065.5))

    # Gompertz hazard: h(t) = lam_1 * exp(alpha_1 * t) * val.
    hazard <- val * lam_1 * exp(alpha_1 * t)

    # Cumulative hazard and survival. The Gompertz integrates to
    # H(t) = lam_1 * (exp(alpha_1 * t) - 1) / alpha_1, but we use the ODE
    # form to keep the source-trace one-to-one with the .mod $DES block.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur             <- exp(-cumhaz)
  })
}
