tmdd_2cmt_full_CLV <- function() {
  description <- "Two compartment TMDD model with full approximation. Parameterized with clearances and volumes."
  ini({
    lka   <- 0.45; label("Absorption rate (Ka)")
    lvc   <- 0.45; label("Central volume of distribution (Vc)")
    lcl    <- 0.04 ; label("Clearance (CL)")
    lq     <- 0.1  ; label("Intercompartmental clearance (Q)")
    lvp    <- 0.38 ; label("Peripheral volume of distribution (Vp)")
    lksyn <- 0.1 ; label("Target Synthesis rate (1/d)") 
    lkeDT <- 0.4 ; label("Internalisation rate (1/d)")
    lkd   <- 0.65; label("Dissociation constant (Kd) ")
    lkon  <- 0.37; label("Drug-target association rate (Kon)")
    lT0   <- 0.36; label("Intial target concentration (T0)")
    lfdepot<- 0.4; label("Bioavailability (F)")
    propSd <- 0.5; label("Proportional residual error (frcentraltion)")
  })
  model({
    ka   <- exp(lka)
    vc   <- exp(lvc)
    cl   <- exp(lcl)
    q    <- exp(lq)
    vp   <- exp(lvp)
    ksyn <- exp(lksyn)
    keDT <- exp(lkeDT)
    kon  <- exp(lkon)
    kd   <- exp(lkd)
    T0   <- exp(lT0)
    fdepot <- exp(lfdepot)
    
    keD <- cl/vc
    k12 <- q/vc
    k21 <- q/vp
    koff <- kd*kon
    keT  <- ksyn/T0
    
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
