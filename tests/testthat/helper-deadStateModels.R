# Fixture models for .checkDeadStates(). Each differs from .dsSplitDoseModel
# only in how its second depot can leave zero.

# A dose split between two depots by f(), with only the first declared: depot2
# has no input, no initial condition and no `dosing` entry.
.dsSplitDoseModel <- function() {
  description <- "Dose split between two depots, second one undeclared"
  reference <- "R"
  dosing <- "depot"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lka2 <- -1; label("Second absorption rate (ka2, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    ka2 <- exp(lka2)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    f(depot) <- 0.6
    f(depot2) <- 0.4
    d/dt(central) <- ka * depot + ka2 * depot2 - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# The same model with depot2 declared.
.dsSplitDoseDeclaredModel <- function() {
  description <- "Dose split between two depots, both declared"
  reference <- "R"
  dosing <- c("depot", "depot2")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lka2 <- -1; label("Second absorption rate (ka2, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    ka2 <- exp(lka2)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    f(depot) <- 0.6
    f(depot2) <- 0.4
    d/dt(central) <- ka * depot + ka2 * depot2 - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# A second depot given an initial amount, a transit state fed by the depot,
# and an effect state with a production term: each can leave zero.
.dsAliveModel <- function() {
  description <- "States that leave zero by initial condition, inflow or production"
  reference <- "R"
  dosing <- "depot"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    lkout <- -1; label("Response turnover rate (kout, 1/h)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kout <- exp(lkout)
    depot2(0) <- 5
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka * depot2
    d/dt(transit1) <- ka * (depot + depot2) - ka * transit1
    d/dt(central) <- ka * transit1 - cl / vc * central
    Cc <- central / vc
    d/dt(effect) <- kout * (1 + Cc) - kout * effect
    Cc ~ prop(propSd)
  })
}

# The input to the peripheral compartment sits on a continuation line that
# starts with `+`. rxode2 joins it to the d/dt() above it, R does not: the
# check must read the equations rxode2 compiles.
.dsContinuationModel <- function() {
  description <- "Two-compartment model with an operator-led continuation line"
  reference <- "R"
  dosing <- "central"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    lq <- 1; label("Intercompartmental clearance (Q, L/h)")
    lvp <- 4; label("Peripheral volume (Vp, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    d/dt(central) <- -cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- -q / vp * peripheral1
    + q / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# A drug concentration held constant (`d/dt(c2) <- 0`) and never dosed.
.dsConstantZeroModel <- function() {
  description <- "Static concentration that is never set"
  reference <- "R"
  dosing <- "c1"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")
  ini({
    lkg <- -1; label("Bacterial growth rate (kg, 1/h)")
    addSd <- 0.1; label("Additive residual error (log10 CFU/mL)")
  })
  model({
    kg <- exp(lkg)
    d/dt(c1) <- 0
    d/dt(c2) <- 0
    bact(0) <- 1e6
    d/dt(bact) <- kg * bact - (c1 + c2) / (1 + c1 + c2) * bact
    log_cfu <- log10(bact)
    log_cfu ~ add(addSd)
  })
}

# States whose only d/dt() term vanishes at 0 in less literal ways: a power, a
# function that is 0 at 0, a variable defined as a multiple of the state, and
# an explicit zero initial condition. aux_exp's term, exp(aux_exp), is 1 at 0,
# so that state is fed.
.dsFactorFormsModel <- function() {
  description <- "States that vanish through powers, functions and aliases"
  reference <- "R"
  dosing <- "central"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    lk <- -1; label("Rate constant (k, 1/h)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    k <- exp(lk)
    d/dt(central) <- -cl / vc * central
    d/dt(aux_pow) <- -k * aux_pow^2
    d/dt(aux_sqrt) <- -k * sqrt(aux_sqrt)
    d/dt(aux_exp) <- -k * exp(aux_exp)
    loss <- k * aux_alias
    d/dt(aux_alias) <- -loss
    aux_zero(0) <- 0
    d/dt(aux_zero) <- -k * aux_zero
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
