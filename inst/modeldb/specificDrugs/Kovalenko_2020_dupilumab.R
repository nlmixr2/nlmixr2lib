Kovalenko_2020_dupilumab <- function() {
  description <- "Dupilumab model (Kovalenko 2020)"
  reference <- "Kovalenko P, Davis JD, Li M, et al. Base and Covariate Population Pharmacokinetic Analyses of Dupilumab Using Phase 3 Data. Clinical Pharmacology in Drug Development. 2020;9(6):756-767. doi:10.1002/cpdd.780"
  # Model 1 from table 1 and supplementary Table 2 in the publication and its
  # supplement.
  covariateData <-
    list(
      WT = "Body weight in kg"
    )
  ini({
    lvc <- log(2.48); label("central volume (L)")
    lke <- log(0.0534); label("elimination rate (1/d)")
    lkcp <- log(0.213); label("central-to-peripheral rate (1/d)")
    Mpc <- 0.686; label("ratio of kcp and kpc (kpc is peripheral to central rate with units of 1/d)")
    lka <- log(0.256); label("absorption rate (1/d)")
    lMTT <- log(0.105); label("mean transit time (d)")
    lVm <- log(1.07); label("maximum target-mediated rate of elimination (mg/L/d)")
    Km <- fixed(0.01); label("Michaelis-Menten constant (mg/L)")
    lfdepot <- log(0.643); label("Bioavailability (fraction)")
    e_wt_vc <- 0.711; label("Exponent of weight on central volume (unitless)")

    etalvc ~ 0.192
    etalke ~ 0.285
    etalka ~ 0.474
    etalvm ~ 0.236
    etamtt ~ 0.525 # etamtt is assumed to be on log-scale MTT to prevent negative values; this is a difference relative to Supplementary Table 2

    cppropSd <- 0.15; label("Proportional residual error (fraction)")
    cpaddSd <- fixed(0.03); label("Additive residual error (mg/L)")
  })
  model({
    # Weight normalization to 75 kg is assumed based on prior publication.  It
    # is not specified in the current publication:
    # Kovalenko P, DiCioccio AT, Davis JD, et al. Exploratory Population PK
    # Analysis of Dupilumab, a Fully Human Monoclonal Antibody Against
    # IL-4Ralpha, in Atopic Dermatitis Patients and Normal Volunteers. CPT
    # Pharmacometrics Syst Pharmacol. 2016;5(11):617-624. doi:10.1002/psp4.12136
    vc <- exp(lvc + etalvc)*(WT/75)^e_wt_vc
    ke <- exp(lke + etalke)
    kcp <- exp(lkcp)
    ka <- exp(lka + etalka)
    MTT <- exp(lMTT + etamtt)
    Vm <- exp(lVm + etalvm)

    # Derived parameters
    kpc <- kcp/Mpc
    ktr <- (3 + 1)/MTT

    d/dt(depot) <- -ktr*depot
    d/dt(transit1) <- ktr*(depot - transit1)
    d/dt(transit2) <- ktr*(transit1 - transit2)
    d/dt(transit3) <- ktr*transit2 - ka*transit3
    # Linear and Michaelis-Menten clearance
    d/dt(central) <-                 ka*transit3 - ke*central - kcp*central + kpc*periph - central*(Vm/(Km + central/vc))
    d/dt(periph) <-                                             kcp*central - kpc*periph

    f(depot) <- exp(lfdepot)
    # No unit conversion is required to change mg/L (dosing amount/central
    # volume unit) to mg/L (measurement unit)
    cp <- central/vc
    cp ~ add(cpaddSd) + prop(cppropSd)
  })
}
