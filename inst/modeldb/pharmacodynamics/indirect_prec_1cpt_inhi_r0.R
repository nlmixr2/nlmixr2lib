indirect_prec_1cpt_inhi_r0 <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response (effect). PK parameterized with clearance and volume; response parameterized with its drug-free baseline R0, from which the precursor production rate is derived."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "drug", units = NA_character_, specimen = "administration site", verified = FALSE),
    central = list(analyte = "drug", units = NA_character_, specimen = "plasma", verified = FALSE),
    precursor1 = list(analyte = "precursor", units = NA_character_, specimen = "not applicable", verified = FALSE),
    effect = list(analyte = "effect", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl   <- 0.04 ; label("Clearance (CL)")
    lrbase     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
    lkin <- 0.48; label("First-order rate constant for transfer of precursor to the response (kin)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response (kout)")
    lkpout <- 0.45 ; label("First order rate constant for loss of precursor (kpout)")
    limax <- log(0.8) ; label("Maximum inhibitory factor attributed to drug (Imax)")
    lic50 <- 0.67; label("Drug concentration producing 50% of maximum inhibition (IC50)")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    cl <- exp(lcl)
    rbase <- exp(lrbase)
    kin <- exp(lkin)
    kout <- exp(lkout)
    kpout <- exp(lkpout)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    fdepot   <- exp(lfdepot)

    kel <- cl/vc
    # Drug-free, the precursor rests at P0 = kpin/(kin + kpout) and the
    # response at R0 = kin*P0/kout. Solved for kpin, that is the production
    # rate below, positive for any positive kin, kout and kpout. At any
    # constant concentration the response returns to R0, because the
    # precursor pool adapts, so the drug moves the response only transiently.
    kpin <- rbase*kout*(kin + kpout)/kin
    Cc <-  central/vc

    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    precursor1(0) <- kpin/(kin + kpout)
    d/dt(precursor1)  <- kpin -(kin + kpout)*(1-imax*Cc/(Cc + ic50))*precursor1
    effect(0) <- kin*kpin/(kout*(kin + kpout))
    d/dt(effect)   <- kin*(1-imax*Cc/(Cc + ic50))*precursor1-kout*effect

    Cc ~ prop(propSd)
  })
}
