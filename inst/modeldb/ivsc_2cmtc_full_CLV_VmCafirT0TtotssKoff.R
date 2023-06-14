ivsc_2cmtc_full_CLV_VmCafirT0TtotssKoff <- function() {
  description <- "Two compartment TMDD model with full approximation. Parameterized with clearances and volumes. "
  ini({
    lka    <- 0.45 ; label("Absorption rate (Ka)")
    lvc    <- 0.45 ; label("Central volume of distribution (Vc)")
    lcl    <- 0.04 ; label("Clearance (CL)")
    lq     <- 0.1  ; label("Intercompartmental clearance (Q)")
    lvp    <- 0.38 ; label("Peripheral volume of distribution (Vp)")
    lvm    <- 0.1  ; label("maximum target-mediated rate of elimination (mg/L/d)")
    lT0    <- 0.36 ; label("Intial target concentration (T0)")
    lTtotss<- 0.5  ; label("Total target concentration(Ttotss)")
    lkoff <- 0.35 ; label("Drug-target dissociation rate (Koff)")
    lcafir  <- 0.1; label ("Averaged Free target concentration to Initial target concentration Ratio")
    lfdepot   <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl)
    q <- exp(lq)
    vp <- exp(lvp)
    vm <- exp(lvm)
    T0 <- exp(lT0)
    Ttotss <- exp(lTtotss)
    koff  <- exp(lkoff)
    fdepot<- exp(lfdepot)
    cafir <- exp(lcafir)
    
    kel <- cl/vc
    k12 <- q /vc
    k21 <- q /Vp    
    ksyn<- vm/vc
    keT <- ksyn/T0
    keDT<- ksyn/Ttotss
    keq <- cafir/(Ttotss/T0)
    kon <- (koff+keDT)/keq
    
    
    d/dt(depot)      <- -ka *depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka *depot - k12*central + k21*peripheral1- keD *central  + (- kon*D*T + koff*DT)*vc
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    d/dt(T)          <- ksyn- keT *T - kon*D*T + koff*DT
    d/dt(DT)         <- -keDT*DT     + kon*D*T - koff*DT
    
    Cc <-  central/Vc
    Cc ~ prop(propSd)
  })
}