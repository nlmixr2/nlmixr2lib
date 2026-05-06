NA_NA_tte_gompertz <- function() {
  description <- "Parametric time-to-event base hazard model for Event 1 in the BAST PTTE 2017 four-event teaching dataset (DDMODEL00000243). The .mod $PROBLEM line names this a 'Gompertz hazard model' but the equation has no time-varying alpha*t term, so the realised hazard is constant: h(t) = (lam/1000) * exp((coef_neut/10000)*(NEUT-4133)) * exp((coef_age/100)*(AGE-55)). The BAST guiding-document text (Figure 2-1, page 13) confirms an exponential distribution was selected for Event 1; the .mod / file name retain the 'Gompertz' label per the source $PROBLEM line and the operator's selected option NA_NA_tte_gompertz.R."
  reference <- paste(
    "BAST Inc Limited.",
    "BAST approach to parametric time-to-event (PTTE) modelling.",
    "Loughborough, UK; 12 July 2017.",
    "Internal guiding document (BAST_PTTE_modelling.pdf) shipped",
    "with DDMORE bundle DDMODEL00000243; no peer-reviewed publication.",
    "Run prepared by Jon Moss (Command.txt; runEV1_201).",
    "DDMORE Foundation Model Repository: DDMODEL00000243.",
    sep = " "
  )
  vignette <- "NA_NA_tte_gompertz"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; covariates AGE and NEUT enter as time-fixed columns)",
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
      notes              = "Time-fixed per subject. Centred at 55 years inside the hazard via `exp((coef_age/100) * (AGE - 55))`. The bundle's Simulated_event_data.csv carries AGE for 200 hypothetical patients with range 24-84 years (mean 58.7); 55 years is approximately the cohort median rounded to the nearest 5.",
      source_name        = "AGE"
    ),
    NEUT = list(
      description        = "Baseline absolute neutrophil count (cells/mm³).",
      units              = "cells/mm^3 (equivalent to cells/uL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Centred at 4133 cells/mm³ inside the hazard via `exp((coef_neut/10000) * (NEUT - 4133))`. The bundle's Simulated_event_data.csv carries NEUT for 200 hypothetical patients with range 1030-14,888 cells/mm³ (mean 4424); 4133 is the cohort median used by the BAST PTTE guiding-document.",
      source_name        = "NEUT"
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
    dose_range     = "Not applicable (no drug administration is modelled; covariates AGE and NEUT enter the hazard at baseline values).",
    regions        = "Not applicable (simulated data).",
    notes          = "200 simulated patients with four timed event types (Event 1, Event 2, Competing Event 1, Competing Event 2) and six baseline covariates (AGE, NEUT, PRE_TRE, MAX_LEG, AUC, CMAX). This file extracts only the Event 1 final hazard model (runEV1_201). For the other three events see NA_NA_tte_gompertz_ev2.R, NA_NA_tte_lognormal.R, and NA_NA_tte_loglogistic.R. Source: BAST Inc Limited, 'BAST approach to parametric time-to-event (PTTE) modelling', 12 July 2017 (BAST_PTTE_modelling.pdf in the DDMORE bundle); guiding document § 2.2.1, Table 2-1; covariate selection result Table 2-2; final-fit listing Output_simulated_runEV1_201.res."
  )

  ini({
    # Final estimates from Output_simulated_runEV1_201.res FINAL PARAMETER
    # ESTIMATE block (LAPLACIAN CONDITIONAL ESTIMATION, MIN SUCCESSFUL,
    # OBJV = 1001.926, NSIG = 4.3). The .mod uses internal numerical
    # rescalings (lam/1000, coef_neut/10000, coef_age/100); the rescalings
    # are kept inside model() so that ini() values match the .lst FINAL
    # THETA values one-to-one.
    llam_ev1   <- log(2.80); label("Event 1 baseline-rate scale parameter, lambda (.mod TH1; rescaled to lam/1000 inside model() to give a 1/day hazard)")  # Output_simulated_runEV1_201.res FINAL TH1 = 2.80E+00
    e_neut_ev1 <- -1.56;     label("Event 1 hazard log-linear coefficient on (NEUT - 4133)/mm^3 (.mod TH2; rescaled to coef/10000 inside model())")        # Output_simulated_runEV1_201.res FINAL TH2 = -1.56E+00
    e_age_ev1  <- 3.20;      label("Event 1 hazard log-linear coefficient on (AGE - 55) years (.mod TH3; rescaled to coef/100 inside model())")           # Output_simulated_runEV1_201.res FINAL TH3 = 3.20E+00

    # No IIV. Source NONMEM run wires ETA(1) on lam with `$OMEGA 0 FIX`
    # (Output_simulated_runEV1_201.res FINAL OMEGA(1,1) = 0 FIXED); the
    # OMEGA was a placeholder and the EV1 hazard model has no estimated
    # IIV. No eta* parameters are added.

    # No residual error. The source $ERROR block computes
    # Y = (1 - CS) * HAZNOW * SURV + CS * SURV with `$EST ... LIKE` so the
    # likelihood is the survival/event density itself, not an
    # observation-error model. This nlmixr2 translation is intended for
    # forward simulation; `hazard` and `sur` are exposed as derived
    # outputs.
  })
  model({
    # Numerical rescalings carried verbatim from Executable_runEV1_201.mod
    # $PK / $DES blocks (`LamC = Lam/1000`, `EXP((TH2/10000)*(NEUT-4133))`,
    # `EXP((TH3/100)*(AGE-55))`). Preserved here so ini() values match the
    # .lst FINAL THETA values directly.
    lam   <- exp(llam_ev1)
    lam_c <- lam / 1000

    # Covariate effect on the hazard (centred-deviation form, multiplicative).
    val <- exp((e_neut_ev1 / 10000) * (NEUT - 4133)) *
      exp((e_age_ev1 / 100) * (AGE - 55))

    # Effective hazard. The .mod $PROBLEM names this "Gompertz" but the
    # equation has no `exp(alpha*t)` factor, so the realised hazard is
    # constant in time (i.e., exponential survival). The BAST guiding
    # document confirms an exponential distribution was selected for
    # Event 1 (§ 2.4.1, Figure 2-1).
    hazard <- val * lam_c

    # Cumulative hazard and survival.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur             <- exp(-cumhaz)
  })
}
