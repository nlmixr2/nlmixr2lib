Jung_2023_im_depot <- function() {
  description <- paste(
    "Two-compartment population PK model for a controlled-release intramuscular injection with",
    "parallel fast and slow release, and a fractal (Kopelman) rate on the fast-release arm.",
    "This is Model Case 3 of Jung 2023, whose contribution is to replace a constant first-order",
    "rate by the time-dependent instantaneous rate coefficient rate = q / time^h (Eq. 1).",
    "A fraction frac of the dose stays in a fast-release muscle depot absorbed at the fractal",
    "rate kf = ka_fast / time^h; the remaining (1 - frac) is delivered into a slow-release",
    "muscle depot through a Savic (2007) transit-compartment chain evaluated in closed form",
    "with Stirling's approximation, and is absorbed at the constant rate ka_slow. NOTE: the",
    "drug is not named in Jung 2023 -- Model Case 3 is described only as in-house",
    "controlled-release intramuscular injection data."
  )
  reference <- paste(
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304. Model Case 3 (fractal model):",
    "Supplementary Table S3 (estimates) and Code S7 (NONMEM control stream).",
    "The transit-compartment closed form follows Savic RM, Jonker DM, Kerbusch T,",
    "Karlsson MO. Implementation of a transit compartment model for describing drug",
    "absorption in pharmacokinetic studies. J Pharmacokinet Pharmacodyn. 2007;34(5):711-726."
  )
  vignette <- "Jung_2023_fractal_kinetics"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # The dose is given to the fast-release depot only (Code S7 F1 = FRAC, F2 = 0);
  # the slow arm is fed by the transit closed form rather than by a dose record.
  dosing <- c("depot")

  compartmentData <- list(
    depot       = list(analyte = "undisclosed", units = "ug", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "undisclosed", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "undisclosed", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "undisclosed", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    disease_state  = "not reported",
    dose_range     = "Single intramuscular injection at four different dose amounts (amounts not reported)",
    n_observations = 339,
    notes          = paste(
      "Jung 2023 Table 2 reports Model Case 3 as 20 subjects and 339 observations; Jung 2023",
      "Section 2.3 states the intramuscular injection was administered once, the model was",
      "fitted over 672 h, and four different drug amounts were dosed. The data are described",
      "only as in-house, so the drug, the dose amounts and the demographics are all",
      "unreported. Species is recorded as human because Jung 2023 presents Case 3 alongside",
      "four other human clinical popPK models and the estimated volumes (Vc 0.356 L,",
      "Vp 2.62 L) and clearance (0.175 L/h) are consistent with a plasma-confined peptide in",
      "an adult; Jung 2023 does not state the species explicitly."
    )
  )

  ini({
    # Structural parameters. Code S7 $THETA carries the converged fractal-model
    # values; they reproduce Supplementary Table S3 to three significant figures.
    lvc      <- log(0.356114)   ; label("Central volume of distribution V1 (L)")                # Code S7 $THETA 1 (0,0.356114); Table S3 fractal V1 = 0.356
    lvp      <- log(2.61518)    ; label("Peripheral volume of distribution V2 (L)")             # Code S7 $THETA 2 (0,2.61518); Table S3 fractal V2 = 2.62
    lcl      <- log(0.17515)    ; label("Clearance CL (L/h)")                                   # Code S7 $THETA 3 (0,0.17515); Table S3 fractal CL = 0.175
    lq       <- log(0.659805)   ; label("Inter-compartmental clearance Q (L/h)")                # Code S7 $THETA 4 (0,0.659805); Table S3 fractal Q = 0.66
    lka_fast <- log(0.103596)   ; label("Fast-release absorption rate coefficient Ka1 (1/h)")   # Code S7 $THETA 6 (0,0.103596); Table S3 fractal Ka1 = 0.104
    lka_slow <- log(0.00179072) ; label("Slow-release absorption rate constant Ka2 (1/h)")      # Code S7 $THETA 7 (0,0.00179072); Table S3 fractal Ka2 = 0.00179
    lmtt     <- log(136.31)     ; label("Mean transit time MTT of the slow-release chain (h)")  # Code S7 $THETA 8 (100,136.31,400); Table S3 fractal MTT = 136
    lntr     <- log(83.7549)    ; label("Number of transit compartments N (unitless)")          # Code S7 $THETA 9 (0,83.7549); Table S3 fractal N = 83.8
    lfrac    <- log(0.172219)   ; label("Fraction of dose to the fast-release depot FRAC (unitless)") # Code S7 $THETA 10 (0,0.172219,1); Table S3 fractal FRAC = 0.172

    # Fractal-kinetics heterogeneity exponent, applied to the fast-release rate.
    h_abs    <- 0.026835        ; label("Heterogeneity exponent h of the fractal fast-release rate (unitless)") # Code S7 $THETA 5 (0,0.026835,1); Table S3 fractal h = 0.0268

    # IIV. Code S7 $OMEGA holds variances; Table S3 CV% = sqrt(exp(omega) - 1) for
    # the log-normal etas. h carries an ADDITIVE eta: Code S7 H = THETA(5) + ETA(6).
    etalvc      ~ 0.0682854     # Code S7 $OMEGA 1; Table S3 fractal V1 IIV 26.6% CV [Shr 45.59%]
    etalq       ~ 0.201086      # Code S7 $OMEGA 2; Table S3 fractal Q IIV 47.2% CV [Shr 32.69%]
    etalka_fast ~ 0.0364168     # Code S7 $OMEGA 3; Table S3 fractal Ka1 IIV 19.3% CV [Shr 23.54%]
    etalntr     ~ 0.00649553    # Code S7 $OMEGA 4; Table S3 fractal N IIV 8.08% CV [Shr 28.27%]
    etalfrac    ~ 0.163404      # Code S7 $OMEGA 5; Table S3 fractal FRAC IIV 42.1% CV [Shr 0.00%]
    etah_abs    ~ 0.926984      # Code S7 $OMEGA 6; Table S3 fractal h IIV 124% [Shr 59.25%]

    # Residual error. Code S7 $ERROR is W = SQRT(THETA(11)^2 * IPRED^2), i.e.
    # PROPORTIONAL only, with $SIGMA 1 FIX -- see the vignette Errata for the
    # disagreement with Supplementary Table S3, which additionally prints an
    # "Add-error" cell that this single-theta error block cannot produce.
    propSd <- 0.112738          ; label("Proportional residual error (fraction)")               # Code S7 $THETA 11 (0,0.112738), labelled "Add.err,centr" but used proportionally
  })

  model({
    # 1. Individual parameters
    vc      <- exp(lvc      + etalvc)
    vp      <- exp(lvp)
    cl      <- exp(lcl)
    q       <- exp(lq       + etalq)
    ka_fast <- exp(lka_fast + etalka_fast)
    ka_slow <- exp(lka_slow)
    mtt     <- exp(lmtt)
    ntr     <- exp(lntr     + etalntr)
    frac    <- exp(lfrac    + etalfrac)

    # Code S7: H = THETA(5) + ETA(6) -- an ADDITIVE eta on a bare exponent, so h
    # is reproduced by addition rather than by the usual exp(l... + eta) form.
    h <- h_abs + etah_abs

    # 2. Micro-constants
    kel <- cl / vc
    kcp <- q  / vc
    kpc <- q  / vp

    # Transit-chain rate constant. Code S7: KT = (N+1)/MTT.
    ktr <- (ntr + 1) / mtt

    # log(N!) by Stirling's approximation. Code S7:
    # LNFAC = LOG(2.5066) + (N+0.5)*LOG(N) - N, where 2.5066 is sqrt(2*pi).
    lnfac <- log(2.5066) + (ntr + 0.5) * log(ntr) - ntr

    # 3. Fractal fast-release rate. Jung 2023 Eq. (1): rate = q_rate / time^h.
    #    Code S7 $PK: KA1 = EXP(LOG(THETA(6)) - H*LOG(TIME) + ETA(3)). Case 3 is
    #    single-dose, so absolute solver time equals time after dose.
    kf <- ka_fast / time^h

    # 4. Savic (2007) transit chain in closed form, delivering the slow fraction
    #    (1 - frac) of the dose into the slow-release depot. Code S7 $DES:
    #    IPT1 = EXP(LOG((1-FRAC)*DOSE1) + LOG(KT) + N*LOG(KT*(TIME-T1))
    #               - KT*(TIME-T1) - LNFAC)
    #    with DOSE1 the dose amount and T1 its time, i.e. podo(depot) and tad()
    #    here. podo() returns the amount before bioavailability, as NONMEM's
    #    PODO() does. The compartment argument is REQUIRED: a bare podo() silently
    #    evaluates to NA through the rxUi / nlmixr2 function-model path (it works
    #    only inside a plain rxode2({}) block), which would make every downstream
    #    state NA.
    tdose <- tad()
    inpt  <- exp(log((1 - frac) * podo(depot)) + log(ktr) + ntr * log(ktr * tdose) -
                   ktr * tdose - lnfac)

    # 5. ODE system. Code S7 $MODEL / $DES, compartments in the published order:
    #    1 D_FAST, 2 D_SLOW, 3 CENT, 4 PERI. Both depots are intramuscular; they
    #    are parallel release arms of one injection, not two dosing routes.
    d/dt(depot)       <- -kf * depot
    d/dt(depot2)      <-  inpt - ka_slow * depot2
    d/dt(central)     <-  kf * depot + ka_slow * depot2 + kpc * peripheral1 -
      kcp * central - kel * central
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1

    # 6. Dose partitioning. Code S7: F1 = FRAC, F2 = 0.
    f(depot) <- frac

    # 7. Observation. Code S7 $ERROR: C_P = A(3)/V1, with the comment
    #    "DV=concentration [ug/L] AMT =ug", so no unit rescaling is needed.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
