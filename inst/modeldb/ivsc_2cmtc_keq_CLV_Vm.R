ivsc_2cmtc_keq_CLV_Vm <- function() {
  description <- "Two compartment TMDD model with QSS approximation. Parameterized with clearances and volumes."
  ini({
    lka    <- 0.45 ; label("Absorption rate (Ka)")
    lvc    <- 0.45 ; label("Central volume of distribution (Vc)")
    lcl    <- 0.04 ; label("Clearance (CL)")
    lq     <- 0.1  ; label("Intercompartmental clearance (Q)")
    lvp    <- 0.38 ; label("Peripheral volume of distribution (Vp)")
    lvm    <- 0.1  ; label("maximum target-mediated rate of elimination (mg/L/d)")
    lkeT  <- 0.1 ; label("Target elimination rate (1/d)")
    lkeDT <- 0.4 ; label("Internalisation rate (1/d)")
    lkeq <- 0.35 ; label("")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl)
    q <- exp(lq)
    vp <- exp(lvp)
    vm  <- exp(lvm)
    keT <- exp(lkeT)
    keDT <- exp(lkeDT)
    keq <- exp(lkeq)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    k12 <- q /vc
    k21 <- q/vp    
    ksyn<- vm/vc
    
    B <- Cc  - Ttot - keq
    D <- 0.5*( B + (B^2 + 4*keq*Cc)^0.5 ) 
    Ac<- D*vc
    DT<- Ttot*D/(keq+D) 
    T <- Ttot - DT    
    
    d/dt(depot)     <- -ka *depot
    f(depot)        <- fdepot
    d/dt(central)  <- ka *Ad - k12*Ac + k21*Ap - kel*Ac - keDT*DT*vc
    d/dt(Ap)     <- k12*Ac - k21*Ap
    d/dt(Ttot)   <- ksyn - keT *T  - keDT*DT
    
    Cc <- central/vc
    Cc ~ prop(propSd)
  })
}