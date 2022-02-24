PK_1cmt <- function() {
  ini({
    ka <- 0.45 # Log Ka
    cl <- 1    # Log Cl
    v  <- 3.45  # log V
    prop.err <- 0.5
  })
  model({
    ka <- exp(ka)
    cl <- exp(cl)
    v  <- exp(v)
    kel <- cl / v;
    d/dt(depot)  = -ka*depot;
    d/dt(centr)  =  ka*depot-kel*centr;
    cp = centr / v;
    cp ~ prop(prop.err)
  })
}
