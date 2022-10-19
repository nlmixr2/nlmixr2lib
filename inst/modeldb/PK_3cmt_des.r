PK_3cmt_des <- function() {
  description <- "Three compartment PK model with linear clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lvp2  <- 8 ; label("Second peripheral volume of distribution (Vp2)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    lq2  <- 0.5 ; label("Second intercompartmental clearance (Q2)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q  <- exp(lq)
    q2  <- exp(lq2)

    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp
    k13 <- q2/vc
    k31 <- q2/vp2

    d/dt(depot) <- -ka*depot
    d/dt(center) <-  ka*depot - kel*center - k12*center + k21*peripheral1 - k13*center + k31*peripheral2
    d/dt(peripheral1) <- k12*center - k21*peripheral1
    d/dt(peripheral2) <- k13*center - k31*peripheral2
    cp <- center / vc

    cp ~ prop(prop.err)
  })
}
