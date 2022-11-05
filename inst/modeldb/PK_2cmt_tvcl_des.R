PK_2cmt_tvcl_des <- function() {
  description <- "Two compartment PK model with time-dependent clearance using differential equations"
  ini({
    lclts <- -1.6 ; label("Time-stationary clearance (CLTS)")
    ltmax <- -1.5 ; label("Typical value of the maximal change of clearance relative to baseline (Tmax)")
    lgamma <- 0 ; label("Hill coefficient")
    lt50 <- 3.4 ; label("Time for 50% of maximal CL change")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lq  <- -0.3 ; label("Intercompartmental clearance (Q)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    clts <- exp(lclts)
    tmax <- exp(ltmax)
    gamma <- exp(lgamma)
    t50 <- exp(lt50)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    cl <- clts*exp(tmax*time^gamma/(t50^gamma+time^gamma))
    
    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(central)     <- - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    cp <- central / vc

    cp ~ prop(propSd)
  })
}
