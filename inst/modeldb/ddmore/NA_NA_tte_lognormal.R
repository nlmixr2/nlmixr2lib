NA_NA_tte_lognormal <- function() {
  description <- "Parametric time-to-event log-normal hazard model for Competing Event 1 in the BAST PTTE 2017 four-event teaching dataset (DDMODEL00000243). Hazard h(t) = val * pdf(t) / (1 - Phi(log((t+DEL)/alpha) / lambda)), where pdf(t) is the log-normal probability density with shape lambda (sigma) and time-scale alpha, val = exp((coef_age/1000) * (AGE - 55)) is the multiplicative AGE effect, and Phi is the standard normal CDF. Competing Event 1 is interval-censored; the BAST guiding-document Section  2.4.1 (Figure 2-3) selected log-normal as the base distribution, then Section  2.4.2 (Table 2-4) retained AGE as the only covariate."
  reference <- paste(
    "BAST Inc Limited.",
    "BAST approach to parametric time-to-event (PTTE) modelling.",
    "Loughborough, UK; 12 July 2017.",
    "Internal guiding document (BAST_PTTE_modelling.pdf) shipped",
    "with DDMORE bundle DDMODEL00000243; no peer-reviewed publication.",
    "Run prepared by Jon Moss (Command.txt; runCOMPEV1_101).",
    "DDMORE Foundation Model Repository: DDMODEL00000243.",
    sep = " "
  )
  vignette <- "NA_NA_tte_lognormal"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the AGE covariate is a per-subject baseline value)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  ddmore_id    <- "DDMODEL00000243"
  replicate_of <- NULL

  covariateData <- list(
    AGE = list(
      description        = "Subject age at enrolment (years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Centred at 55 years inside the hazard via `exp((coef_age/1000) * (AGE - 55))`. The bundle's Simulated_event_data.csv carries AGE for 200 hypothetical patients with range 24-84 years (mean 58.7); 55 years is approximately the cohort median rounded to the nearest 5.",
      source_name        = "AGE"
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
    dose_range     = "Not applicable (no drug administration is modelled).",
    regions        = "Not applicable (simulated data).",
    notes          = "200 simulated patients; 36 (18%) had Competing Event 1. Competing Event 1 is interval-censored: exact event times are unknown, only that the event occurred between two scheduled assessment visits (BAST guiding document Section  2.2.1). For the other three events in the same bundle see NA_NA_tte_gompertz.R, NA_NA_tte_gompertz_ev2.R, and NA_NA_tte_loglogistic.R. Source: BAST Inc Limited, 'BAST approach to parametric time-to-event (PTTE) modelling', 12 July 2017 (BAST_PTTE_modelling.pdf); base-distribution selection Section  2.4.1 / Figure 2-3; covariate selection Section  2.4.2 / Table 2-4; final-fit listing Output_simulated_runCOMPEV1_101.res."
  )

  ini({
    # Final estimates from Output_simulated_runCOMPEV1_101.res FINAL
    # PARAMETER ESTIMATE block (LAPLACIAN CONDITIONAL ESTIMATION,
    # MIN SUCCESSFUL, OBJV = 360253.678, NSIG = 4.2). The .mod applies
    # internal numerical rescaling only on the AGE covariate
    # (coef_age/1000); lambda and alpha enter the hazard at their raw
    # THETA values. ini() values therefore match the .lst FINAL THETA
    # values one-to-one.
    llambda_compev1 <- log(1.42);   label("Competing Event 1 log-normal lambda parameter (sigma; .mod TH1)")        # Output_simulated_runCOMPEV1_101.res FINAL TH1 = 1.42E+00
    lalpha_compev1  <- log(984);    label("Competing Event 1 log-normal alpha parameter (median time scale, days; .mod TH2)")  # Output_simulated_runCOMPEV1_101.res FINAL TH2 = 9.84E+02
    e_age_compev1   <- 50.9;        label("Competing Event 1 hazard log-linear coefficient on (AGE - 55) years (.mod TH3; rescaled to coef/1000 inside model())")  # Output_simulated_runCOMPEV1_101.res FINAL TH3 = 5.09E+01

    # No IIV. Source NONMEM run wires ETA(1) on lambda with `$OMEGA 0 FIX`
    # (Output_simulated_runCOMPEV1_101.res FINAL OMEGA(1,1) = 0 FIXED).
    # No eta* parameters are added.

    # No residual error. Source $ERROR block uses `$EST ... LIKE`
    # (the likelihood is the survival / event density). Forward-simulation
    # translation: `hazard` and `sur` are exposed as derived outputs.
  })
  model({
    # Back-transformed shape and scale.
    lambda_compev1 <- exp(llambda_compev1)
    alpha_compev1  <- exp(lalpha_compev1)

    # Covariate effect on the hazard (centred-deviation form, multiplicative).
    val <- exp((e_age_compev1 / 1000) * (AGE - 55))

    # Log-normal hazard. Source .mod $DES uses a small-time offset
    # `DEL = 1e-8` to keep `log(t/alpha)` finite at t = 0; we preserve
    # that here. `pnorm(.)` is the standard normal CDF (NONMEM `phi(.)`).
    pi_const <- 3.1415927
    fac      <- sqrt(2 * pi_const)
    del      <- 1e-8
    num      <- log((t + del) / alpha_compev1)
    pdf_t    <- (1 / (fac * lambda_compev1 * (t + del))) *
      exp(-(num^2) / (2 * lambda_compev1 * lambda_compev1))
    surv_ln  <- 1 - pnorm(num / lambda_compev1)

    hazard <- val * pdf_t / surv_ln

    # Cumulative hazard and survival.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur             <- exp(-cumhaz)
  })
}
