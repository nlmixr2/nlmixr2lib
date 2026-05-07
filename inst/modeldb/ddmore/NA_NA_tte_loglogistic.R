NA_NA_tte_loglogistic <- function() {
  description <- "Parametric time-to-event log-logistic hazard model for Competing Event 2 in the BAST PTTE 2017 four-event teaching dataset (DDMODEL00000243). Hazard h(t) = alpha * (lam^alpha) * t^(alpha-1) / (1 + (lam*t)^alpha), where lam = lambda/1000 and alpha is unscaled. No covariates were retained. Note: the BAST guiding-document text (Section  2.4.1, Figure 2-4) selected log-normal as the base distribution for Competing Event 2, but the executable supplied with the bundle (Executable_runCOMPEV2_005.mod, $PROBLEM 'Log_logistic model') and the final-fit listing in the bundle are log-logistic, not log-normal. This file follows the executable; the publication-vs-bundle distribution discrepancy is documented in the validation vignette's Errata."
  reference <- paste(
    "BAST Inc Limited.",
    "BAST approach to parametric time-to-event (PTTE) modelling.",
    "Loughborough, UK; 12 July 2017.",
    "Internal guiding document (BAST_PTTE_modelling.pdf) shipped",
    "with DDMORE bundle DDMODEL00000243; no peer-reviewed publication.",
    "Run prepared by Jon Moss (Command.txt; runCOMPEV2_005).",
    "DDMORE Foundation Model Repository: DDMODEL00000243.",
    sep = " "
  )
  vignette <- "NA_NA_tte_loglogistic"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; no covariates in this final model)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  ddmore_id    <- "DDMODEL00000243"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    n_subjects     = 200L,
    n_studies      = 1L,
    age_range      = "24-84 years (mean 58.7) in the BAST PTTE 2017 simulated cohort",
    weight_range   = "not reported (the BAST PTTE 2017 simulated cohort does not include body weight)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Hypothetical / unspecified clinical population (the BAST PTTE 2017 guiding document is a methodological teaching example with simulated event data; no real drug, indication, or patient cohort).",
    dose_range     = "Not applicable (no drug administration is modelled).",
    regions        = "Not applicable (simulated data).",
    notes          = "200 simulated patients; 142 (71%) had Competing Event 2. Competing Event 2 is right-censored only (CENSORING = 1; exact event times are observed). For the other three events in the same bundle see NA_NA_tte_gompertz.R, NA_NA_tte_gompertz_ev2.R, and NA_NA_tte_lognormal.R. Source: BAST Inc Limited, 'BAST approach to parametric time-to-event (PTTE) modelling', 12 July 2017 (BAST_PTTE_modelling.pdf); base-distribution discussion Section  2.4.1 / Figure 2-4 (selected log-normal in text); covariate-test result Table 2-5 (no covariate retained); final-fit executable Executable_runCOMPEV2_005.mod (log-logistic) and listing Output_simulated_runCOMPEV2_005.res."
  )

  ini({
    # Final estimates from Output_simulated_runCOMPEV2_005.res FINAL
    # PARAMETER ESTIMATE block (LAPLACIAN CONDITIONAL ESTIMATION,
    # MIN SUCCESSFUL, OBJV = 1846.946, NSIG = 3.2). The .mod applies a
    # /1000 rescaling on lambda only (`Lam1 = Lam/1000`); alpha is used
    # without rescaling (`alpha1 = alpha`). ini() values match the .lst
    # FINAL THETA values one-to-one.
    llam_compev2   <- log(5.93); label("Competing Event 2 log-logistic lambda scale parameter (.mod TH1; rescaled to lam/1000 inside model())")  # Output_simulated_runCOMPEV2_005.res FINAL TH1 = 5.93E+00
    lalpha_compev2 <- log(1.60); label("Competing Event 2 log-logistic alpha shape parameter (.mod TH2; unscaled)")                              # Output_simulated_runCOMPEV2_005.res FINAL TH2 = 1.60E+00

    # No IIV. Source NONMEM run wires ETA(1) on lambda with `$OMEGA 0 FIX`
    # (Output_simulated_runCOMPEV2_005.res FINAL OMEGA(1,1) = 0 FIXED).
    # No eta* parameters are added.

    # No residual error. Source $ERROR block uses `$EST ... LIKE`
    # (the likelihood is the survival / event density). Forward-simulation
    # translation: `hazard` and `sur` are exposed as derived outputs.
  })
  model({
    # Back-transformed shape and scale.
    lam_compev2   <- exp(llam_compev2)
    alpha_compev2 <- exp(lalpha_compev2)

    # Numerical rescaling carried verbatim from
    # Executable_runCOMPEV2_005.mod $PK (`Lam1 = Lam/1000`,
    # `alpha1 = alpha`).
    lam_1   <- lam_compev2 / 1000
    alpha_1 <- alpha_compev2

    # Log-logistic hazard:
    #   h(t) = alpha * lam^alpha * t^(alpha - 1) / (1 + (lam * t)^alpha)
    # Source .mod $DES does NOT add the BAST guiding-document `DEL = 1e-8`
    # offset to t inside `t^(alpha - 1)`, so the hazard is computed without
    # an offset here too. At t = 0 this evaluates to 0 (alpha > 1), which
    # is finite and correct.
    val_1 <- alpha_1 * (lam_1^alpha_1) * t^(alpha_1 - 1)
    val_2 <- (lam_1^alpha_1) * (t^alpha_1)
    hazard <- val_1 / (1 + val_2)

    # Cumulative hazard and survival.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur             <- exp(-cumhaz)
  })
}
