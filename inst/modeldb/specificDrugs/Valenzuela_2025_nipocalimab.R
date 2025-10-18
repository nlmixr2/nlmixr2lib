Valenzuela_2025_nipocalimab <- function() {
  description <- "Nipocalimab PK and PK/PD models in myasthenia gravis (Valenzuela 2025)"
  reference <- "Valenzuela B, Neyens M, Zhu Y, et al. Nipocalimab Dose Selection in Generalized Myasthenia Gravis. CPT: Pharmacometrics & Systems Pharmacology. n/a(n/a). doi:10.1002/psp4.70109"
  units <- list(time="hr", dosing="mg")
  covariates <-
    list(
      ELISA = "1 for ELISA PK assay; 0 for ECLIA assay",
      PHASE1 = "1 for Phase 1 study; 0 for Phase 2 study"
    )

  # All parameters are from Table 3
  ini({
    # PK model
    lcl <- log(0.655); label("Linear serum clearance (L/d)")
    bsv_cl ~ cvToBsv(24.1)
    allo_cl <- fixed(0.75); label("Allometric exponent for clearance (unitless)")
    lvc <- log(3.23); label("Volume of the central compartment (L)")
    bsv_vc ~ cvToBsv(15.0)
    allo_v <- fixed(1); label("Allometric exponent for volume (unitless)")
    lq <- log(0.250); label("Intercompartmental clearance (L/d)")
    lvp <- log(0.622); label("Volume of the peripheral compartment (L)")
    lFcRn0 <- log(143); label("Total (free) FcRn at baseline (nmol/L)")
    bsv_FcRn0 ~ cvToBsv(25.0)
    FRmax <- 0.947; label("Maximal fraction of FcRn bound to nipocalimab (fraction)")
    lKss <- log(6.05); label("Quasi-steady-state dissociation equilibrium constant (ug/mL)")
    lkint <- log(62.4); label("Internalization rate constant for FcRn-nipocalimab complex (1/d)")
    lkdeg <- log(1.3); label("Degradation rate constant of the free FcRn (1/d)")
    addSdPKELISA <- 0.445; label("additive RUV for PK model, ELISA assay (nmol/L)")
    addSdPKECLIA <- 0.0342; label("additive RUV for PK model, ECLIA assay (nmol/L)")
    propSdPKPhase1 <- 0.0834; label("proportional RUV for Phase 1 studies (fraction)")
    propSdPKPhase2 <- 0.367; label("proportional RUV for Phase 2 studies (fraction)")

    # PK-IgG Model
    lIgG0 <- log(11.4); label("IgG at baseline (g/L)")
    bsv_IgG0 ~ cvToBsv(21.9)
    FRIgG0_M281_004 <- 0.777; label("Fraction of IgG at baseline for study MOM-M281-004 (fraction)")
    lkdeg_IgG <- log(0.217); label("IgG degradation rate constant (without FcRn recycling) (1/d)")
    bsv_kdeg_IgG ~ cvToBsv(16.9)
    IgK <- 5.08; label("FcRn recycling (unitless)")
    bsv_IgK ~ cvToBsv(26.3)
    # lKtotal <- log(0.036); label("IgG degradation rate constant (with full FcRn recycling) (1/d)")
    # IgGmaxRed <- 83.6; label("Maximum IgG reduction (at full inhibition of FcRn recycling) (%)")
    propSdIgG <- 0.0858; label("Proportional RUV for IgG (fraction)")

    addSdRO <- 2.98; label("Additive RUV for RO (%)")
    propSdRO <- 0.227; label("Proportional RUV for RO (fraction)")

    # IgG/MG-ADL Model
    lke0 <- log(0.414); label("IgG effect compartment rate constant (1/d)")
    SIgG <- -0.216; label("Slope between MG-ADL and IgG reduction (points/10% IgG reduction)")
    EIgG <- 0.871; label("Exponent of MG-ADL baseline effect on SIgG (unitless)")
    IDecplacebo <- -1.08; label("Initial placebo decrease in MGADL after start of treatment (points)")
    bsv_SIgG + bsv_IDecplacebo ~ c(cvToBsv(3.76), -0.733, cvToBsv(1.89))
    EADL <- 1.23; label("Exponent of MG-ADL baseline effect on IDecplacebo (unitless)")
    Splacebo <- -0.0594; label("Slope between MG-ADL and time (placebo effect) (points/week)")
    bsv_Splacebo ~ cvToBsv(0.125)
    addSdMGADL <- 1.50; label("Additive RUV for MG-ADL (points)")
  })
  model({
    # From supplement S4

    # PK model
    cl <- exp(lcl + bsv_cl) * (WT/75)^allo_cl
    vc <- exp(lvc + bsv_vc) * (WT/75)^allo_vc
    q <- exp(lq) * (WT/75)^allo_cl
    vp <- exp(lvp) * (WT/75)^allo_vc

    FcRn0 <- exp(lFcRn0 + bsv_FcRn0)
    Kss <- exp(lKss)
    kint <- exp(lkint)
    kdeg <- exp(lkdeg)

    addSdPK = addSdPKELISA*ELISA + addSdPKECLIA*(1 - ELISA)
    propSdPK = propSdPKPhase1*PHASE1 + propSdPKPhase2*(1 - PHASE1)

    d/dt(central) <- CL*Cfree - kint * vc * (Rtotal*Cfree)/(Kss + Cfree) - Q*Cfree + Q/vp*periph1
    d/dt(Rtotal) <- Ksyn - kdeg*Rtotal - (kint - kdeg)*(Rtotal*Cfree)/(Kss + Cfree)
    d/dt(periph1) <- Q*Cfree - Q/vp*periph1

    # PK/PD model
    IgG0 <- exp(lIgG0 + bsv_IgG0)
  })
}
