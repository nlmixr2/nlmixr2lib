Jung_2023_donepezil_singledose <- function() {
  description <- paste(
    "Two-compartment population PK model for single-dose donepezil given orally or as a",
    "transdermal patch, with a fractal (Kopelman) release rate on the patch arm. This is",
    "Model Case 1 of Jung 2023, whose contribution is to replace a constant first-order rate",
    "by the time-dependent instantaneous rate coefficient rate = q / time^h (Eq. 1), where h",
    "is a heterogeneity exponent bounded in (0, 1). Oral dose enters a gut depot absorbed at",
    "Ka; patch dose enters a formulation reservoir that empties into a two-compartment transit",
    "chain and thence into the same central compartment. Only the reservoir-to-transit1 rate",
    "is fractal (kf = ktr / time^h); the two downstream transit rates stay constant at ktr.",
    "Both routes share one central and one peripheral compartment and one linear clearance."
  )
  reference <- paste(
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304. Model Case 1 (fractal model):",
    "Supplementary Table S1 (estimates) and Code S3 (NONMEM control stream).",
    "The underlying base model and the clinical study are reported in",
    "Jung W, Jung H, Vu N-AT, Kim G-Y, Kim G-W, Chae J-w, Kim T, Yun H-y.",
    "Model-Based Equivalent Dose Optimization to Develop New Donepezil Patch Formulation.",
    "Pharmaceutics. 2022;14(2):244. doi:10.3390/pharmaceutics14020244."
  )
  vignette <- "Jung_2023_fractal_kinetics"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Dose enters either the oral depot or the transdermal reservoir; neither is
  # named `depot`, so the dosing targets are declared explicitly.
  dosing <- c("depot_oral", "depot_td")

  compartmentData <- list(
    depot_oral  = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    depot_td    = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 18,
    n_studies      = 1,
    age_range      = "24-33 years",
    age_median     = "30.0 years",
    weight_range   = "55.7-80.9 kg",
    weight_median  = "63.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy adult Korean male volunteers",
    dose_range     = "Single dose: donepezil 10 mg oral tablet (Aricept) or a 108 mg / 96 cm2 transdermal patch",
    regions        = "Republic of Korea",
    n_observations = 383,
    notes          = paste(
      "Jung 2023 Table 2 reports Model Case 1 as 18 subjects and 383 observations, fitted over",
      "312 h. The underlying study (Jung 2022, ref [3] of Jung 2023) is a randomised, open-label,",
      "two-treatment, two-sequence, two-period crossover bioequivalence study in healthy Korean",
      "men: 12 enrolled, 9 evaluable in both periods (Jung 2022 Table 1), so the 18 modelled",
      "subjects correspond to 9 subjects x 2 crossover periods entered as separate records.",
      "Demographics above are the 9 evaluable subjects of Jung 2022 Table 1."
    )
  )

  ini({
    # Structural parameters. Values are the FINAL fractal-model estimates: Code S3
    # $THETA carries the converged values (they reproduce Supplementary Table S1 to
    # three significant figures), so these are final estimates and not initial ones.
    lka   <- log(0.0834578) ; label("Oral absorption rate constant Ka (1/h)")               # Code S3 $THETA 1 (0,0.0834578); Table S1 fractal Ka = 0.0835 (RSE 54.4%)
    lcl   <- log(9.68116)   ; label("Clearance CL (L/h)")                                   # Code S3 $THETA 2 (0,9.68116); Table S1 fractal CL = 9.68 (RSE 8.89%)
    lvc   <- log(51.0518)   ; label("Central volume of distribution Vc (L)")                # Code S3 $THETA 3 (0,51.0518); Table S1 fractal Vc = 51.1 (RSE 62.3%)
    lvp   <- log(564.22)    ; label("Peripheral volume of distribution Vp (L)")              # Code S3 $THETA 4 (0,564.22); Table S1 fractal Vp = 564 (RSE 9.27%)
    lq    <- log(28.4)      ; label("Inter-compartmental clearance Q (L/h)")                 # Code S3 $THETA 5 (0,28.4); Table S1 fractal Q = 28.4 (RSE 52.3%)
    lktr  <- log(0.043948)  ; label("Transdermal transit rate constant Kt (1/h)")            # Code S3 $THETA 7 (0,0.043948,1); Table S1 fractal Kt = 0.0439 (RSE 21.3%)

    # Fractal-kinetics heterogeneity exponent. Bounded (0, 1) in the control
    # stream; h = 0 recovers the constant-rate (Fick) base model.
    h_abs <- 0.320373       ; label("Heterogeneity exponent h of the fractal transdermal release rate (unitless)")  # Code S3 $THETA 6 (0,0.320373,1); Table S1 fractal h = 0.32 (RSE 24.3%)

    # IIV. Code S3 $OMEGA holds variances of log-normal etas; the CV% printed in
    # Table S1 is recovered as sqrt(exp(omega) - 1).
    etalka  ~ 0.0241974     # Code S3 $OMEGA 1; Table S1 fractal Ka IIV 15.7% CV [Shr 44.85%]
    etalcl  ~ 0.114149      # Code S3 $OMEGA 2; Table S1 fractal CL IIV 34.7% CV [Shr 0.00%]
    etalvc  ~ 0.161231      # Code S3 $OMEGA 3; Table S1 fractal Vc IIV 41.8% CV [Shr 42.73%]
    etalktr ~ 0.0309662     # Code S3 $OMEGA 4; Table S1 fractal Kt IIV 17.7% CV [Shr 30.57%]

    # Residual error. Code S3 $ERROR: W = SQRT(THETA(8)^2 + THETA(9)^2 * IPRED^2),
    # i.e. combined additive + proportional on the ng/mL scale, with $SIGMA 1 FIX.
    addSd  <- 2.53394       ; label("Additive residual error (ng/mL)")                       # Code S3 $THETA 8 (0,2.53394); Table S1 fractal Add-error = 2.53 (RSE 16.5%)
    propSd <- 0.0909088     ; label("Proportional residual error (fraction)")                # Code S3 $THETA 9 (0,0.0909088); Table S1 fractal Prop-error = 0.0909 (RSE 25.1%)
  })

  model({
    # 1. Individual parameters
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ktr <- exp(lktr + etalktr)

    # 2. Micro-constants
    kel <- cl / vc
    kcp <- q  / vc
    kpc <- q  / vp

    # 3. Fractal release rate from the patch reservoir. Jung 2023 Eq. (1):
    #    rate = q_rate / time^h, an instantaneous rate coefficient rather than a
    #    rate constant. Code S3 $PK: KF = EXP(LOG(KT) - H*LOG(TIME)), i.e. the
    #    transit rate constant Kt divided by time^h. `time` is the modelled time
    #    from the point of dosing (Jung 2023 Section 2.1); Case 1 is single-dose,
    #    so absolute solver time equals time after dose. h < 1 makes the
    #    singularity at time = 0 integrable.
    kf <- ktr / time^h_abs

    # 4. ODE system. Code S3 $MODEL / $DES, compartments in the published order:
    #    1 GUT, 2 CENT, 3 PERI, 4 SKIN, 5 TRAN1, 6 TRAN2.
    d/dt(depot_oral)  <- -ka * depot_oral
    d/dt(central)     <-  ka * depot_oral + kpc * peripheral1 - kcp * central -
      kel * central + ktr * transit2
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1
    d/dt(depot_td)    <- -kf * depot_td
    d/dt(transit1)    <-  kf * depot_td - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2

    # 5. Observation. Amounts are mg and vc is L, so central/vc is mg/L;
    #    1 mg/L = 1000 ng/mL. Code S3 $ERROR: IPRED = A(2)/(VC/1000).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
