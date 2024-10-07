PK_2cmt_mAb_Davda_2014 <- function() {
  description <- "Two compartment PK model with linear clearance for average monoclonal antibodies (Davda 2014)"
  reference <- "Davda JP, Dodds MG, Gibbs MA, Wisdom W, Gibbs JP. A model-based meta-analysis of monoclonal antibody pharmacokinetics to guide optimal first-in-human study design. MAbs. 2014;6(4):1094-1102. doi:10.4161/mabs.29095"
  units = list(time = "day", dosing = "mg")
  ini({
    lfdepot <- log(0.744) ; label("Subcutaneous bioavailability (fraction)")
    lka <- log(0.282) ; label("Absorption rate (Ka, 1/day)")
    lcl <- log(0.200) ; label("Clearance (CL, L/day)")
    lv  <- log(3.61) ; label("Central volume of distribution (V, L)")
    lvp  <- log(2.75) ; label("Peripheral volume of distribution (Vp, L)")
    lq  <- log(0.747) ; label("Intercompartmental clearance (Q, L/day)")

    allocl <- 0.865 ; label("Allometric exponent on clearance and intercompartmental clearance (unitless)")
    allov <- 0.957 ; label("Allometric exponent on volumes of distribution (unitless)")

    etafdepot ~ 0
    etaka ~ 0.416
    etacl + etav + etavp ~ c(0.0987,
                             0.0786, 0.116,
                             0.0377, 0.0619, 0.0789)
    etaq ~ 0.699

    prop.err <- 0.144 ; label("Proportional residual error (fraction)")
  })
  model({
    # WT is body weight in kg
    fdepot <- exp(lfdepot + etafdepot)
    ka <- exp(lka + etaka)
    wtnorm <- log(WT/70)
    cl <- exp(lcl + allocl*wtnorm + etacl)
    q  <- exp(lq + allocl*wtnorm + etaq)
    v <- exp(lv + allov*wtnorm + etav)
    vp <- exp(lvp + allov*wtnorm + etavp)

    Cc <- linCmt()
    
    f(depot) <- fdepot  # Units are dosing units/L (typically mg/L = ug/mL)
    Cc ~ prop(prop.err)
  })
}
