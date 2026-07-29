Dao_2020_sultiame <- function() {
  description <- "Population PK model for sultiame in healthy adult volunteers with non-linear distribution into erythrocytes (saturable binding to a putative red-blood-cell carrier). Four-compartment structure: depot (oral absorption, KA fixed at 1/h), central (plasma), erythrocytes (drug bound to a saturable carrier in red blood cells parameterised by KON, KOFF, BTOT), and urine (cumulative urinary excretion as a fraction QREN of total elimination). Drug binding to erythrocytes is written in mass-action form on amounts (KON in 1/(h*mg)). DDMORE Foundation Model Repository entry DDMODEL00000298, fit on 4 healthy volunteers (433 observations) by NONMEM ADVAN13 FOCEI."
  reference <- paste(
    "Dao K, Thoueille P, Decosterd LA, Mercier T, Guidi M, Bardinet C,",
    "Lebon S, Choong E, Castang A, Guittet C, Granier LA, Buclin T (2020).",
    "Pharmacokinetic profile of sultiame in healthy volunteers with in vitro",
    "characterization of its uptake by red blood cells.",
    "Pharmacology Research & Perspectives 8(2):e00558.",
    "doi:10.1002/prp2.558.",
    "DDMORE Foundation Model Repository: DDMODEL00000298.",
    sep = " "
  )
  vignette <- "Dao_2020_sultiame"
  ddmore_id <- "DDMODEL00000298"
  replicate_of <- NULL
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list()

  population <- list(
    n_subjects     = 4L,
    n_studies      = 1L,
    age_range      = "Healthy adult volunteers",
    weight_range   = "Healthy adult volunteers",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy adult volunteers receiving single oral doses of sultiame.",
    dose_range     = "Single oral doses of 50, 100, and 200 mg sultiame in the DDMORE bundle's simulated event table; the linked publication (Dao 2020, Pharmacology Research & Perspectives 8(2):e00558) describes a Phase 1 dose-ranging PK study in healthy volunteers.",
    regions        = "Switzerland (University Hospital of Lausanne, per .lst NONMEM license)",
    notes          = "n_subjects = 4 is the count of distinct individuals in the NONMEM run (Output_real_sultiame_nonlinear_PK.lst, 'TOT. NO. OF INDIVIDUALS: 4', 433 observations). The full Dao 2020 publication is not on disk for this extraction; the population characteristics above summarise what the DDMORE bundle and the publication's title disclose. Individual demographics (age, weight, sex) are not reproduced in the bundle's simulated dataset (Simulated_data_PK_sultiame.csv); the underlying Export_PK_Nonmem_urine.csv referenced by the .mod has WT and AGE columns but they are not used by the structural model."
  )

  ini({
    # Structural parameters - final estimates from Output_real_sultiame_nonlinear_PK.lst
    # FINAL PARAMETER ESTIMATE block (after MINIMIZATION SUCCESSFUL on line 858).
    # Each TH<i> maps to the labelled $THETA slot in Executable_sultiame_nonlinear_PK.mod.
    lcl   <- log(11.0);   label("Plasma clearance CL (L/h)")                                              # Output_real_*.lst FINAL THETA TH 1 = 1.10E+01 (CL)
    lvc   <- log(56.3);   label("Plasma volume of distribution V2 (L)")                                   # Output_real_*.lst FINAL THETA TH 2 = 5.63E+01 (V2 in source; nlmixr2 canonical lvc for the central / plasma volume)
    lvp   <- log(2.93);   label("Erythrocyte distribution volume V3 (L)")                                 # Output_real_*.lst FINAL THETA TH 3 = 2.93E+00 (V3 in source; nlmixr2 canonical lvp for the first peripheral-distribution volume - the saturable erythrocyte compartment); .mod fixes V3 with no IIV
    lkon  <- log(0.949);  label("Plasma->erythrocyte association rate KON (1/(h*mg))")                    # Output_real_*.lst FINAL THETA TH 4 = 9.49E-01 (KON)
    lbtot <- log(97.1);   label("Maximal erythrocyte binding capacity BTOT (mg)")                         # Output_real_*.lst FINAL THETA TH 5 = 9.71E+01 (BTOT)
    lkoff <- log(0.796);  label("Erythrocyte->plasma dissociation rate KOFF (1/h)")                       # Output_real_*.lst FINAL THETA TH 6 = 7.96E-01 (KOFF)
    lka   <- fixed(log(1)); label("First-order absorption rate KA (1/h, fixed)")                          # Output_real_*.lst FINAL THETA TH 7 = 1.00E+00 (KA, '1 FIX' in $THETA)
    lqren <- log(0.247);  label("Renal extraction fraction QREN (unitless, bounded 0-1 in source)")       # Output_real_*.lst FINAL THETA TH 8 = 2.47E-01 (QREN)

    # Inter-individual variability - diagonal $OMEGA on the linear parameter scale
    # (NONMEM source: 'CL = THETA(1) * EXP(ETA(1))', etc.). The variance applies on the
    # log-transformed scale used here (lcl, lv2, lqren, lbtot).
    etalcl   ~ 0.0761    # Output_real_*.lst FINAL OMEGA(1,1) = 7.61E-02 (IIV CL)
    etalvc   ~ 0.00858   # Output_real_*.lst FINAL OMEGA(2,2) = 8.58E-03 (IIV V2 in source = lvc here)
    etalqren ~ 0.0995    # Output_real_*.lst FINAL OMEGA(3,3) = 9.95E-02 (IIV QREN); ETA3 attached to QREN per .mod $PK ordering
    etalbtot ~ 0.0145    # Output_real_*.lst FINAL OMEGA(4,4) = 1.45E-02 (IIV BTOT); ETA4 attached to BTOT per .mod $PK ordering

    # Residual error - in the source $ERROR block EPS(1) is fixed to 1 (SIGMA = 1 FIX)
    # and the per-output proportional / additive standard deviations are encoded as
    # THETA(9..14). W = SQRT(THETApro^2 * IPRED^2 + THETAadd^2) and Y = IPRED + W*EPS(1)
    # i.e. additive plus proportional residual error, with the THETA values being SDs.
    propSd        <- 0.567   ; label("Plasma proportional residual SD (fraction)")                         # Output_real_*.lst FINAL THETA TH 9  = 5.67E-01 (Prop.RE plasma)
    # Plasma additive SD (TH10) is fixed to 0 in the source and therefore omitted here.
    propSd_Crbc   <- 0.263   ; label("Erythrocyte proportional residual SD (fraction)")                    # Output_real_*.lst FINAL THETA TH 11 = 2.63E-01 (Prop.RE eryth)
    addSd_Crbc    <- 0.012   ; label("Erythrocyte additive residual SD (mg/L)")                            # Output_real_*.lst FINAL THETA TH 12 = 1.20E-02 (Add.RE eryth)
    propSd_Aurine <- 0.433   ; label("Cumulative urine proportional residual SD (fraction)")               # Output_real_*.lst FINAL THETA TH 13 = 4.33E-01 (Prop.RE urine)
    # Urine additive SD (TH14) is fixed to 0 in the source and therefore omitted here.
  })

  model({
    # Individual parameters - back-transformed from the log-scale typicals plus etas.
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp)
    kon  <- exp(lkon)
    btot <- exp(lbtot + etalbtot)
    koff <- exp(lkoff)
    ka   <- exp(lka)
    qren <- exp(lqren + etalqren)

    # Plasma elimination rate (1/h) - 'KE = CL/V2' in source $PK; here ke = cl/vc.
    ke <- cl / vc

    # Mass-action ODEs - amounts in mg. Mirrors $DES of Executable_sultiame_nonlinear_PK.mod
    # (DADT(1)..DADT(4)). depot, central, erythrocytes, urine map to A(1)..A(4) respectively.
    # The non-linear plasma <-> erythrocyte binding is written on amounts: KON has units
    # 1/(h*mg), so KON * central * (btot - erythrocytes) = mg/h. Drug leaves plasma at
    # rate KE * central; only fraction QREN of that is captured in the urine compartment
    # (the (1 - QREN) fraction represents non-renal elimination).
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - ke * central - kon * central * (btot - erythrocytes) + koff * erythrocytes
    d/dt(erythrocytes) <-  kon * central * (btot - erythrocytes) - koff * erythrocytes
    d/dt(urine)        <-  ke * central * qren

    # Observation outputs.
    #   Cc      - plasma concentration (mg/L), canonical central-compartment output
    #   Crbc    - erythrocyte-bound drug concentration (mg/L), paper-named multi-output
    #   Aurine  - cumulative urinary amount (mg). The source uses urinary concentration
    #             (IPRED4 = A(4) / S4 with S4 = UVOL/1000), but UVOL is a per-record
    #             collection volume that must be supplied alongside each urine
    #             observation. To keep the packaged model self-contained and free of
    #             a per-observation external column, the urine output is exposed as
    #             cumulative amount; consumers can divide by their own urine volume to
    #             recover concentration. A purely proportional residual error is the
    #             same in fractional terms whether applied to amount or concentration,
    #             so AurinepropSd carries over unchanged.
    Cc     <- central       / vc
    Crbc   <- erythrocytes  / vp
    Aurine <- urine

    Cc     ~ prop(propSd)
    Crbc   ~ add(addSd_Crbc) + prop(propSd_Crbc)
    Aurine ~ prop(propSd_Aurine)
  })
}
