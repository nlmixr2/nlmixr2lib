# Description: Two compartment PK model with linear clearance using differential equations
PK_2cmt_des <- function() {
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lv  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    v <- exp(lv)
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl/v
    k12 <- q/v
    k21 <- q/vp

    d/dt(depot) <- -ka*depot
    d/dt(center) <-  ka*depot - kel*center - k12*center + k21*peripheral1
    d/dt(peripheral1) <- k12*center - k21*peripheral1
    cp <- center / v

    cp ~ prop(prop.err)
  })
}
