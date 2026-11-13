# Fixture models for the solve probe in helper-solveProbe.R. Each one carries
# the single feature a probe rule exists for, so the rule can be tested without
# a registry model.

# Between-occasion variability on clearance: rxSolve needs an OCC column with
# at least two levels (.probeGroupingVars()).
.pbOccModel <- function() {
  description <- "One-compartment model with between-occasion variability on clearance"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    etalcl ~ 0.1
    etaiovlcl ~ 0.05 | OCC
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl + etalcl + etaiovlcl)
    vc <- exp(lvc)
    d/dt(central) <- -cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# An algebraic concentration endpoint alongside an ODE-state PD endpoint, for
# the `first` observation form of probeEvents().
.pbTwoEndpointModel <- function() {
  description <- "One-compartment PK with an indirect-response PD state"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    lkout <- -1; label("Response turnover rate (kout, 1/h)")
    etalcl ~ 0.1
    propSd <- 0.1; label("Proportional residual error (fraction)")
    addSd_effect <- 0.1; label("Additive residual error on the response (unitless)")
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    kout <- exp(lkout)
    d/dt(central) <- -cl / vc * central
    Cc <- central / vc
    effect(0) <- 1
    d/dt(effect) <- kout * (1 + Cc) - kout * effect
    Cc ~ prop(propSd)
    effect ~ add(addSd_effect)
  })
}

# A 100 h absorption lag, longer than the probe's default 48 h grid.
.pbLagModel <- function() {
  description <- "One-compartment model with a long absorption lag"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    ltlag <- log(100); label("Absorption lag (h)")
    etalcl ~ 0.1
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    alag(depot) <- exp(ltlag)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# Neonatal covariates with no reference anywhere in the model, so each falls
# through to .probeNameDefaults: no reference category, no normalising
# constant, no numeric metadata.
.pbNeonateModel <- function() {
  description <- "Neonatal model whose covariates have no reference value"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  covariateData <- list(
    PNA = list(description = "Postnatal age", units = "months", type = "continuous"),
    GA = list(description = "Gestational age at birth", units = "weeks", type = "continuous"),
    HCT = list(description = "Hematocrit", units = "(fraction)", type = "continuous"),
    T_INFUSION_LAG = list(description = "Infusion-line lag", units = "h", type = "continuous")
  )
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    e_pna_cl <- 0.1; label("Effect of postnatal age on CL (1/month)")
    e_ga_vc <- 0.01; label("Effect of gestational age on Vc (1/week)")
    etalcl ~ 0.1
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl + etalcl + e_pna_cl * PNA) * (1 - HCT)
    vc <- exp(lvc + e_ga_vc * GA)
    d/dt(central) <- -cl / vc * central
    alag(central) <- T_INFUSION_LAG
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# FPG appears in `FPG - 3.5`, a floor rather than a reference, while its
# covariateData states the reference 8.2; WT has no metadata and is read from
# `WT / 70`.
.pbMetadataReferenceModel <- function() {
  description <- "Model whose covariate metadata states a reference the equations do not"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  covariateData <- list(
    FPG = list(description = "Fasting plasma glucose", units = "mmol/L", type = "continuous", reference_value = 8.2)
  )
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    etalcl ~ 0.1
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl + etalcl) * (FPG - 3.5) / 4.7 * (WT / 70)
    vc <- exp(lvc)
    d/dt(central) <- -cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
