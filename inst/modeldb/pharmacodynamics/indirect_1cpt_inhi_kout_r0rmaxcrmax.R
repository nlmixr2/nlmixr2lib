indirect_1cpt_inhi_kout_r0rmaxcrmax <- function() {
  description <- "One compartment indirect response model with inhibition of kout."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "drug", units = NA_character_, specimen = "administration site", verified = FALSE),
    central = list(analyte = "drug", units = NA_character_, specimen = "plasma", verified = FALSE),
    effect = list(analyte = "response", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  ini({
    lka     <- 0.45 ; label("Absorption rate (Ka)")
    lvc     <- 3.45 ; label("Central volume of distribution (Vc)")
    lcl     <- 0.85 ; label("Clearance (Cl)")
    lrbase     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
    lrmax   <- 0.9  ; label("Highest response, reached when the concentration is CRmax (Rmax)")
    ls1     <- 1.0  ; label("Magnitude of the initial slope of the response versus time curve at high concentration (|S1|)")
    limax   <- log(0.8) ; label("Maximum inhibitory factor attributed to drug (Imax)")
    lcrmax  <- 0.67 ; label("Plasma concentration of drug at the time of maximal response (CRmax)")
    lfdepot <- 0.4  ; label("Bioavailability (F)")
    propSd  <- 0.5  ; label("Proportional residual error (fraction)")
  })
  model({
    ka     <- exp(lka)
    vc     <- exp(lvc)
    cl     <- exp(lcl)
    rbase     <- exp(lrbase)
    rmax   <- exp(lrmax)
    s1     <- exp(ls1)
    imax   <- exp(limax)
    crmax  <- exp(lcrmax)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    # Drug-free, the response rests at R0 = kin/kout. At high concentration
    # it first rises at kin*Imax, so |S1| = kin*Imax. It is highest, Rmax,
    # when dR/dt = 0 at concentration CRmax:
    #   kin = kout*(1 - Imax*CRmax/(CRmax + IC50))*Rmax,
    # which gives IC50 below and needs R0 < Rmax < R0/(1 - Imax).
    kin <- s1/imax
    kout <- kin/rbase
    IC50 <- crmax*(rbase-(1-imax)*rmax)/(rmax-rbase)
    Cc <-  central/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(kel)*central
    d/dt(effect) <- kin - kout*(1-imax*Cc/(Cc + IC50))*effect
    
    Cc ~ prop(propSd)
  })
}
