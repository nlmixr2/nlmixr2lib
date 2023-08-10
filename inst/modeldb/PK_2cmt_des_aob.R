PK_2cmt_des_aob <- function() {
  description <- "Two compartment PK model with linear clearance using differential equations"
  ini({
    lka    <- 0.45 ;label("Absorption rate (Ka)")
    laob   <- 0.3  ;    label("First macro-constant (A) over second macro-constant (B)")
    lalpha <- 0.4;     label("first rate constant")
    lbeta  <- 0.5;  label("second rate constant")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka   <- exp(lka)
    aob  <-exp(laob)
    alpha<-exp(lalpha)
    beta <-exp(lbeta)
    
    k21 <- (aob*beta+alpha)/(aob+1)
    kel <- alpha*beta/k21
    k12 <- alpha+beta-k21-kel
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}