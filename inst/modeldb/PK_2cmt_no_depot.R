PK_2cmt_no_depot <- function() {
  description <- "Two compartment PK model with linear clearance using differential equations"
  ini({
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    
    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp
    
  
    d/dt(central) <-  - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
