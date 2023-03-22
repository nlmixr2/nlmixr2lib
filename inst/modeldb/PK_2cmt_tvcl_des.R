PK_2cmt_tvcl_des <- function() {
  description <- "Two compartment PK model with time-dependent clearance using differential equations"
  reference <- "Justin to add"
  ini({
    lcl <- log(0.2) ; label("Time-stationary clearance (CLTS)")
    lcltmax <- log(0.22) ; label("Typical value of the maximal change of clearance relative to baseline (Tmax)")
    lclgamma <- log(1) ; label("Hill coefficient for time-dependent clearance")
    lclt50 <- log(30) ; label("Time for 50% of maximal CL change")
    lvc  <- log(20) ; label("Central volume of distribution (V)")
    lvp  <- log(150) ; label("Peripheral volume of distribution (Vp)")
    lq  <- log(0.75) ; label("Intercompartmental clearance (Q)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    clts <- exp(lcl)
    cltmax <- exp(lcltmax)
    clgamma <- exp(lclgamma)
    clt50 <- exp(lclt50)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    cl <- clts*exp(cltmax*time^clgamma/(clt50^clgamma+time^clgamma))

    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(central)     <- - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    cp <- central / vc

    cp ~ prop(propSd)
  })
}
