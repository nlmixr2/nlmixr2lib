# Fixture models for .checkBareStatements(): the same two-compartment model
# with the peripheral input written on an operator-led continuation line, and
# with the operator moved to the end of the line above.

.bsContinuationModel <- function() {
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

.bsJoinedModel <- function() {
  description <- "Two-compartment model with the operator at the line end"
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
    d/dt(peripheral1) <- -q / vp * peripheral1 +
      q / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
