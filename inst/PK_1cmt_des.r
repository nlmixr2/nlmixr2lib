# Description: One compartment PK model with linear clearance using differential equations
PK_1cmt_des <- function() {
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lv  <- 3.45 ; label("Central volume of distribution (V)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    v  <- exp(lv)
    kel <- cl / v;
    d/dt(depot)  = -ka*depot;
    d/dt(centr)  =  ka*depot-kel*centr;
    cp = centr / v;
    cp ~ prop(prop.err)
  })
}
