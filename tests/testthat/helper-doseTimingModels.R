# Fixture models for the dose-history compartment rule
# (.checkDoseTimingCompartments() and .modelDosingCmts()). They differ only in
# the `dosing` metadata and in the compartment a timing function names, so each
# test isolates one branch of the rule.

.dtStomachModel <- function() {
  description <- "Stomach release model with no dosing metadata"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    release <- podo(stomach) * ka * exp(-ka * tad(stomach))
    d/dt(stomach) <- -release
    d/dt(central) <- release - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

.dtStomachDosedModel <- function() {
  description <- "Stomach release model dosed into the stomach"
  reference <- "R"
  dosing <- "stomach"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    release <- podo(stomach) * ka * exp(-ka * tad(stomach))
    d/dt(stomach) <- -release
    d/dt(central) <- release - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# No `dosing` metadata, so the fallback is depot and central, which is where
# tad() looks.
.dtFallbackDepotModel <- function() {
  description <- "Lagged depot model relying on the depot/central fallback"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka) * (tad(depot) > 0.5)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# Bare tad() and podo() read the most recent dose into any compartment, so
# they name nothing to check even though `stomach` is not a dosing compartment.
.dtBareTimingModel <- function() {
  description <- "Stomach release model using compartment-free timing"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    release <- podo() * ka * exp(-ka * tad())
    d/dt(stomach) <- -release
    d/dt(central) <- release - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# Morning and evening doses in two depots, both declared; a third depot read
# by tafd() is not.
.dtTwoDepotModel <- function() {
  description <- "Two declared depots and one undeclared"
  reference <- "R"
  dosing <- c("depot", "depot2")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka) * (tad(depot) > 0.5) * (tad0(depot2) > 0.5)
    late <- tafd(depot3) > 24
    cl <- exp(lcl) * (1 + late)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka * depot2
    d/dt(depot3) <- -ka * depot3
    d/dt(central) <- ka * (depot + depot2 + depot3) - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# `dosing <- character(0)` declares that nothing is dosed (a model driven by
# a covariate), so every compartment-named timing read is undosed.
.dtEmptyDosingModel <- function() {
  description <- "Stomach release model declaring no dosing compartment"
  reference <- "R"
  dosing <- character(0)
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    release <- podo(stomach) * ka * exp(-ka * tad(stomach))
    d/dt(stomach) <- -release
    d/dt(central) <- release - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# Neither depot nor central, and no `dosing`: the registry records NA.
.dtNoFallbackModel <- function() {
  description <- "Stomach release model with no depot or central"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    release <- podo(stomach) * ka * exp(-ka * tad(stomach))
    d/dt(stomach) <- -release
    d/dt(plasma) <- release - cl / vc * plasma
    Cc <- plasma / vc
    Cc ~ prop(propSd)
  })
}

# Fallback fixtures: depot and central, no `dosing`, differing in how central
# is reached and whether it is a route of its own.

# depot feeds central, which also takes an infusion of its own (`dur()`).
.dtInfusedCentralModel <- function() {
  description <- "Oral and infused model relying on the depot/central fallback"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    dur(central) <- 0.5
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# depot feeds central, and the model reads central's dose history.
.dtCentralTimingModel <- function() {
  description <- "Oral model reading the time since a central dose"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl) * (1 + 0.1 * (tad(central) < 1))
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# depot reaches central only through a transit compartment.
.dtTransitFallbackModel <- function() {
  description <- "Oral model with a transit compartment relying on the fallback"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lktr <- 0; label("Transit rate (ktr, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    ktr <- exp(lktr)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# depot feeds a metabolite pool, not central, so central is a route of its
# own.
.dtSeparateDepotModel <- function() {
  description <- "Depot feeding a metabolite, central dosed on its own"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(metabolite) <- ka * depot - cl / vc * metabolite
    d/dt(central) <- -cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# A linCmt() solution with an absorption rate: its depot feeds its central.
.dtLinCmtDepotModel <- function() {
  description <- "One-compartment linCmt() model with first-order absorption"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}

# depot feeds a metabolite pool; in d/dt(central) it appears only in a term
# that removes from central, which is not a feed.
.dtLossOnlyModel <- function() {
  description <- "Depot removing from central without feeding it"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(metabolite) <- ka * depot - cl / vc * metabolite
    d/dt(central) <- -ka * depot - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# depot feeds central through a variable assigned once, which is followed to
# its definition.
.dtAliasFeedModel <- function() {
  description <- "Oral model whose absorption rate is a named variable"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    absorb <- ka * depot
    d/dt(depot) <- -absorb
    d/dt(central) <- absorb - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# depot feeds central through a variable assigned in each branch of an
# if/else. A variable assigned more than once is not followed, so the feed is
# not seen.
.dtIfElseAliasFeedModel <- function() {
  description <- "Oral model whose absorption rate is assigned in an if/else"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    if (t > 1) {
      absorb <- ka * depot
    } else {
      absorb <- 0.5 * ka * depot
    }
    d/dt(depot) <- -absorb
    d/dt(central) <- absorb - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# d/dt(central) is assigned in each branch of an if/else, and depot feeds it
# in both.
.dtIfElseFeedModel <- function() {
  description <- "Oral model whose central equation is assigned in an if/else"
  reference <- "R"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    if (t > 1) {
      d/dt(central) <- ka * depot - cl / vc * central
    } else {
      d/dt(central) <- 0.5 * ka * depot - cl / vc * central
    }
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# .dtFallbackDepotModel, in which depot feeds central, with one statement
# appended to its model() block, given as text ("dur(central) <- 0.5").
.dtFallbackVariant <- function(statement) {
  ui <- suppressMessages(rxode2::rxode2(.dtFallbackDepotModel))
  suppressMessages(eval(bquote(rxode2::model(ui, .(str2lang(statement)), append = TRUE))))
}

# Reading this as an active binding fails, as reading an unparseable
# model() block does.
.dtUnreadableModelBlock <- function() {
  stop("the model() block could not be read")
}

# The registry's dosing compartments for a fixture model.
.dtDosingCmts <- function(model) {
  nlmixr2lib:::.modelDosingCmts(suppressMessages(nlmixr2est::nlmixr(model)))
}

# A read of an undeclared compartment inside one branch of an if/else.
.dtIfElseModel <- function() {
  description <- "Depot model with a stomach timing read inside an if/else"
  reference <- "R"
  dosing <- "depot"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    if (t > 1) {
      kgut <- ka * (tad(stomach) > 0.5)
    } else {
      kgut <- ka
    }
    d/dt(depot) <- -ka * depot
    d/dt(stomach) <- -kgut * stomach
    d/dt(central) <- ka * depot + kgut * stomach - cl / vc * central
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}

# Variables named after dose-history functions: `tad <- tad(depot)` reads a
# declared compartment, and `dose <- 5` is an assignment, not a call. Neither
# is flagged.
.dtShadowNamesModel <- function() {
  description <- "Depot model with variables named tad and dose"
  reference <- "R"
  dosing <- "depot"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  ini({
    lka <- 0; label("Absorption rate (ka, 1/h)")
    lcl <- 1; label("Clearance (CL, L/h)")
    lvc <- 3; label("Central volume (Vc, L)")
    propSd <- 0.1; label("Proportional residual error (fraction)")
  })
  model({
    tad <- tad(depot)
    dose <- 5
    ka <- exp(lka) * (tad > 0.5)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    Cc <- central / vc * dose / 5
    Cc ~ prop(propSd)
  })
}
