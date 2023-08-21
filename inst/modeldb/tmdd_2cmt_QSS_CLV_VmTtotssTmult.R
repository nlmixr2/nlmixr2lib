tmdd_2cmt_QSS_CLV_VmTtotssTmult <- function() {
  description <- "Two compartment TMDD model with QSS approximation. Parameterized with clearances and volumes. "
  ini({
    lka    <- 0.45 ; label("Absorption rate (Ka)")
    lvc    <- 0.45 ; label("Central volume of distribution (Vc)")
    lcl    <- 0.04 ; label("Clearance (CL)")
    lq     <- 0.1  ; label("Intercompartmental clearance (Q)")
    lvp    <- 0.38 ; label("Peripheral volume of distribution (Vp)")
    lvm    <- 0.1  ; label("maximum target-mediated rate of elimination (mg/L/d)")
    lT0    <- 0.1  ; label("")
    lTmult <- 0.4  ; label("")
    lTtotss<- 0.5  ; label("Total target concentration")
    lkeq   <- 0.35 ; label("")
    lfdepot<- 0.4  ; label("Bioavailability (F)")
    propSd <- 0.5  ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl)
    q <- exp(lq)
    vp <- exp(lvp)
    vm <- exp(lvm)
    T0 <- exp(lT0)
    Ttmult <- exp(lTmult)
    Ttotss <- exp(lTtotss)
    keq <- exp(lkeq)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    k12 <- q /vc
    k21 <- q /Vp    
    ksyn<- vm/vc
    T0  <- Ttotss*Tmult
    keT <- ksyn/T0
    keDT<- ksyn/Ttotss
    
    
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
