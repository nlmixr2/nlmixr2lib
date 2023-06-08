PK_2cmt_vss   <- function() {
  description <- "Two compartment PK model with linear clearance"
  ini({
    lka    <- 0.45 ; label("Absorption rate (Ka)")
    lcl    <- 1 ; label("Clearance (CL)")
    lvc    <- 3 ; label("Central volume of distribution (V)")
    lvss   <- 5 ; label("Volume at steady-state (Vss)")
    lq     <- 0.1 ; label("Intercompartmental clearance (Q)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    vss<- exp(lvss)
    q  <- exp(lq)
    
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
