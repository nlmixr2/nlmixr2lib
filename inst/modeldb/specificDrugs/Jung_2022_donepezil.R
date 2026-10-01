Jung_2022_donepezil <- function() {
  description <- paste(
    "Two-compartment population PK model for donepezil given orally or as a once-weekly",
    "transdermal patch in healthy adult Korean men (Jung 2022). Oral dose enters a gut",
    "depot absorbed at first-order Ka; patch dose enters a formulation reservoir that",
    "empties at first-order Kt through two transit compartments (same Kt) into the same",
    "central compartment. Both routes share one central and one peripheral compartment",
    "and one linear clearance. The fraction of the nominal patch strength that reaches",
    "the skin is fixed from the in-vitro dissolution equation (Jung 2022 Eq. 3) at the",
    "one-week (168 h) wear time used throughout the study and its simulations."
  )
  reference <- paste(
    "Jung W, Jung H, Vu N-AT, Kim G-Y, Kim G-W, Chae J-w, Kim T, Yun H-y.",
    "Model-Based Equivalent Dose Optimization to Develop New Donepezil Patch Formulation.",
    "Pharmaceutics. 2022;14(2):244. doi:10.3390/pharmaceutics14020244.",
    "Unrounded final estimates and the NONMEM control stream are published as",
    "Code S2 ('case 1, base model', OFV 1443.703) in the supplement of",
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304."
  )
  vignette <- "Jung_2022_donepezil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Dose enters either the oral depot or the transdermal reservoir; neither is
  # named `depot`, so the dosing targets are declared explicitly.
  dosing <- c("depot_oral", "depot_td")

  compartmentData <- list(
    depot_oral = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    depot_td = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 9,
    n_studies = 1,
    age_range = "24-33 years",
    age_median = "30.0 years",
    weight_range = "55.7-80.9 kg",
    weight_median = "63.1 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy adult male volunteers",
    dose_range = paste(
      "Donepezil 10 mg oral tablet (Aricept) once daily for 7 days, or one",
      "108 mg / 96 cm2 transdermal patch worn for 1 week, in a two-period crossover"
    ),
    regions = "Not stated; ethics approval by Raptim Research Ltd (Mumbai, India), sponsor and analysts in the Republic of Korea",
    n_observations = 383,
    notes = paste(
      "Randomised, open-label, two-treatment, two-sequence, two-period crossover",
      "bioequivalence study (TL/WZ/19/001141) with a washout of at least 21 days.",
      "Twelve subjects enrolled and nine completed both periods; only the nine were",
      "modelled (Jung 2022 Section 3.1, Table 1). The 383 observations are from Jung 2023",
      "Table 2, which refits the same dataset (Case 1) and reproduces this model's OFV",
      "(1443.70). Sampling to 312 h after the first dose."
    )
  )

  ini({
    # Structural parameters. Values are the unrounded final estimates of Code S2
    # (Jung 2023 supplement), which is this model: its OFV 1443.703 equals Jung 2022
    # Table 3, and every $THETA rounds to the Table 3 estimate.
    lka  <- log(0.0496711) ; label("Oral absorption rate constant Ka (1/h)")          # Jung 2022 Table 3 Ka = 0.0497 (RSE 25%); Code S2 $THETA 1 = 0.0496711
    lcl  <- log(10.0268)   ; label("Clearance CL (L/h)")                              # Jung 2022 Table 3 CL = 10 (RSE 9%); Code S2 $THETA 2 = 10.0268
    lvc  <- log(26.2189)   ; label("Central volume of distribution Vc (L)")           # Jung 2022 Table 3 Vc = 26.2 (RSE 35%); Code S2 $THETA 3 = 26.2189
    lvp  <- log(562.037)   ; label("Peripheral volume of distribution Vp (L)")        # Jung 2022 Table 3 Vp = 562 (RSE 11%); Code S2 $THETA 4 = 562.037
    lq   <- log(15.6292)   ; label("Inter-compartmental clearance Q (L/h)")           # Jung 2022 Table 3 Q = 15.6 (RSE 33%); Code S2 $THETA 5 = 15.6292
    lktr <- log(0.0270191) ; label("Transdermal release and transit rate constant Kt (1/h)") # Jung 2022 Table 3 Kt = 0.027 (RSE 9%); Code S2 $THETA 6 = 0.0270191

    # Fraction of the nominal patch strength that enters the skin reservoir.
    # Jung 2022 Eq. 3: Drug dissolution = 78.257 * Duration / (Duration + 8.481) %
    # of patch dose (fitted to the in-vitro dissolution data, Table S3), and
    # Disposed amount in skin = 0.74 * Drug dissolution. Duration is the wear time
    # in hours; the study and all simulations used a one-week patch, Duration = 168,
    # giving 0.74 * 0.78257 * 168 / 176.481 = 0.5513. Not estimated: in Code S2 this
    # fraction is pre-applied to the patch AMT and no bioavailability term appears.
    lfdepot_td <- fixed(log(0.74 * 0.78257 * 168 / (168 + 8.481))) ; label("Fraction of nominal patch dose delivered to skin over a 168 h wear (fraction)") # Jung 2022 Eq. 3 with Duration = 168 h (Section 2.1 'patches ... for 1 week')

    # IIV. Code S2 $OMEGA holds variances of exponential (log-normal) etas; the
    # CV% column of Table 3 is recovered as sqrt(exp(omega) - 1).
    etalka  ~ 0.00968106   # Jung 2022 Table 3 IIV Ka 0.00968 (9.9% CV) [Shr 51%]; Code S2 $OMEGA 1
    etalcl  ~ 0.130076     # Jung 2022 Table 3 IIV CL 0.13 (37.3% CV) [Shr 0%]; Code S2 $OMEGA 2
    etalvc  ~ 0.197923     # Jung 2022 Table 3 IIV Vc 0.198 (46.8% CV) [Shr 42%]; Code S2 $OMEGA 3
    etalktr ~ 0.020045     # Jung 2022 Table 3 IIV Kt 0.02 (14.2% CV) [Shr 31%]; Code S2 $OMEGA 4

    # Residual error. Code S2 $ERROR: W = SQRT(THETA(7)^2 + THETA(8)^2 * IPRED^2)
    # with $SIGMA 1 FIX, i.e. combined additive + proportional on the ng/mL scale.
    addSd  <- 2.89074   ; label("Additive residual error (ng/mL)")        # Jung 2022 Table 3 additive error 2.89 (RSE 13%); Code S2 $THETA 7 = 2.89074
    propSd <- 0.0795173 ; label("Proportional residual error (fraction)") # Jung 2022 Table 3 proportional error 0.0795 (RSE 29%); Code S2 $THETA 8 = 0.0795173
  })

  model({
    # 1. Individual parameters (Code S2 $PK; no IIV on Vp or Q)
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ktr <- exp(lktr + etalktr)

    # 2. Micro-constants (Jung 2022 Eq. 4)
    kel <- cl / vc
    kcp <- q / vc
    kpc <- q / vp

    # 3. ODE system. Jung 2022 Eqs. 1, 2 and 4; Code S2 $DES, compartments in the
    #    published order 1 GUT, 2 CENT, 3 PERI, 4 SKIN, 5 TRAN1, 6 TRAN2. Eq. 4 prints
    #    the oral input to CENT as KA*SKIN; Code S2 and Eq. 1 use KA*GUT.
    d/dt(depot_oral)  <- -ka * depot_oral
    d/dt(central)     <-  ka * depot_oral + kpc * peripheral1 - kcp * central -
      kel * central + ktr * transit2
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1
    d/dt(depot_td)    <- -ktr * depot_td
    d/dt(transit1)    <-  ktr * depot_td - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2

    # 4. Patch delivered fraction (Jung 2022 Eq. 3); dose the NOMINAL patch strength.
    f(depot_td) <- exp(lfdepot_td)

    # 5. Observation. Amounts are mg and vc is L, so central/vc is mg/L;
    #    1 mg/L = 1000 ng/mL. Code S2 $ERROR: IPRED = A(2)/(VC/1000).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
