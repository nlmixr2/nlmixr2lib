PK_1cmt_des <- function() {
  description <- "One compartment PK model with linear clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    prop.err <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka*depot
    d/dt(center) <- ka*depot-kel*center

    cp = center / vc
    cp ~ prop(prop.err)
  })
}
