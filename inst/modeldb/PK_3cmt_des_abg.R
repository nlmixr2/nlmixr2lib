PK_3cmt_des_abg <- function() {
  description <- "Three compartment PK model with linear clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lk21 <- 0.24; label("rate constant from peripheral1 to central")
    lk31 <- 0.39; label ("rate constant from peripheral2 to central")
    lalpha <- 0.4; label("first rate constant")
    lbeta  <- 0.5;  label("second rate constant")
    lgamma <- 0.45; label("third rate constant")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    k21 <- exp(lk21)
    k31 <- exp(lk31)
    alpha <- exp(lalpha)
    beta <- exp (lbeta)
    gamma <- exp (lgamma)
    
    kel <- alpha*beta*gamma/k31
    vc <- alpha+beta+gamma
    vp <- alpha*beta+alpha*gamma+beta*gamma
    k13 <- (vp+k31*k31-k31*vc-kel*k21)/(k21-k31)
    k12 <- vc-kel-k13-k21-k31
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - kel*central - k12*central + k21*peripheral1 - k13*central + k31*peripheral2
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    d/dt(peripheral2) <- k13*central - k31*peripheral2
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
