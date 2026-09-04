Jung_2023_donepezil_multipledose <- function() {
  description <- paste(
    "Two-compartment population PK model for multiple-dose donepezil given as an oral",
    "titration followed by three transdermal patches, with a fractal (Kopelman) release rate",
    "on the patch arm that is RESET at every new dose. This is Model Case 2 of Jung 2023,",
    "whose contribution is to replace a constant first-order rate by the time-dependent",
    "instantaneous rate coefficient rate = q / time^h (Eq. 1), where h is a heterogeneity",
    "exponent bounded in (0, 1). Oral dose enters a gut depot absorbed at Ka; patch dose",
    "enters a formulation reservoir that empties at the fractal release rate into a single",
    "transit compartment and thence into the same central compartment. Both routes share one",
    "central and one peripheral compartment and one linear clearance. Case 2 is the only one",
    "of the five cases whose fractal clock restarts at each dose, and the only one for which",
    "the authors specify an onset-time regulariser (tau = 1e-7 h). Jung 2023 does not name the",
    "drug for this case; donepezil is inferred - see the note in `population`."
  )
  reference <- paste(
    "Jung W, Ryu H-j, Chae J-w, Yun H-y. Fractal Kinetic Implementation in Population",
    "Pharmacokinetic Modeling. Pharmaceutics. 2023;15(1):304.",
    "doi:10.3390/pharmaceutics15010304. Model Case 2 (fractal model):",
    "Supplementary Table S2 (estimates) and Code S5 (NONMEM control stream).",
    "The single-dose companion model and the donepezil patch programme are reported in",
    "Jung W, Jung H, Vu N-AT, Kim G-Y, Kim G-W, Chae J-w, Kim T, Yun H-y.",
    "Model-Based Equivalent Dose Optimization to Develop New Donepezil Patch Formulation.",
    "Pharmaceutics. 2022;14(2):244. doi:10.3390/pharmaceutics14020244;",
    "see modellib('Jung_2023_donepezil_singledose') for Model Case 1."
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
    transit1    = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 44,
    n_studies      = 1,
    disease_state  = "Healthy adult volunteers",
    dose_range     = paste(
      "Multiple dose: 7 days of oral titration followed by three transdermal-patch doses,",
      "observed to 2496 h. Jung 2023 Section 2.3 states only that two different amounts were",
      "dosed for the oral and transdermal-patch routes; the amounts themselves are not reported."
    ),
    n_observations = 3024,
    notes          = paste(
      "Jung 2023 Table 2 reports Model Case 2 as 44 subjects and 3024 observations.",
      "DRUG IDENTITY. Jung 2023 never names the drug for Case 2, and unlike Case 1 the",
      "description carries no citation. Donepezil is an inference, adopted per operator sidecar",
      "oare_PMC9867137 request-001 q2 (option A), on this evidence: the case comes from the same",
      "first author and group as Case 1; it has the same integrated oral-plus-patch structure;",
      "Codes S3 and S5 share the identical concentration scaling IPRED = A(2)/(VC/1000), i.e. mg",
      "dosing with ng/mL observations; and the disposition estimates closely match Case 1's",
      "donepezil fit (CL 8.97 vs 9.68 L/h, Vss 657 vs 615 L). Ref [3] of Jung 2023 also",
      "anticipates this exact design ('The model can handle complicated dosing plans such as",
      "giving an oral titration period in patch study'). Against the inference: the paper never",
      "states it, and three patches over 2496 h is about 35 days per patch, which does not match",
      "the weekly patch of ref [3]. Treat the drug label as provisional; the parameter values and",
      "structure are transcribed exactly as published and do not depend on the identification.",
      "Species is likewise not stated by Jung 2023; human is inferred from the sibling cases,",
      "all of which are human clinical studies."
    )
  )

  ini({
    # Structural parameters. Values are the FINAL fractal-model estimates: Code S5
    # $THETA carries the converged values (they reproduce Supplementary Table S2 to
    # three significant figures), so these are final estimates and not initial ones.
    lka   <- log(0.914658)  ; label("Oral absorption rate constant Ka (1/h)")                  # Code S5 $THETA 1 (0,0.914658); Table S2 fractal Ka = 0.915 (RSE 11.9%)
    lcl   <- log(8.97291)   ; label("Clearance CL (L/h)")                                      # Code S5 $THETA 2 (0,8.97291); Table S2 fractal CL = 8.97 (RSE 3.44%)
    lvc   <- log(243.65)    ; label("Central volume of distribution Vc (L)")                   # Code S5 $THETA 3 (0,243.65); Table S2 fractal Vc = 244 (RSE 5.6%)
    lvp   <- log(413.263)   ; label("Peripheral volume of distribution Vp (L)")                # Code S5 $THETA 4 (0,413.263); Table S2 fractal Vp = 413 (RSE 3.69%)
    lq    <- log(56.6586)   ; label("Inter-compartmental clearance Q (L/h)")                   # Code S5 $THETA 5 (0,56.6586); Table S2 fractal Q = 56.7 (RSE 2.51%)
    lktr  <- log(0.0186256) ; label("Transdermal transit rate constant Kt (1/h)")              # Code S5 $THETA 7 (0,0.0186256); Table S2 fractal Kt = 0.0186 (RSE 5.15%)

    # Base rate of the fractal patch-release step. Unlike Case 1 -- where the
    # fractal rate is built from the transit constant Kt -- Case 2 estimates a
    # separate release rate Kf for the reservoir-emptying step.
    lkrel <- log(0.398403)  ; label("Release rate constant Kf from the patch formulation reservoir (1/h)") # Code S5 $THETA 8 (0,0.398403); Table S2 fractal Kf = 0.398 (RSE 20.9%)

    # Fractal-kinetics heterogeneity exponent. Bounded (0, 1) in the control
    # stream; h = 0 recovers the constant-rate (Fick) base model.
    h_abs <- 0.894192       ; label("Heterogeneity exponent h of the fractal transdermal release rate (unitless)") # Code S5 $THETA 6 (0,0.894192,1); Table S2 fractal h = 0.894 (RSE 7.01%)

    # Onset-time regulariser, fixed by the authors rather than estimated. Case 2
    # is the only one of the five cases that carries one; it keeps the release
    # rate finite at the instant of each dose, where the elapsed-time clock is 0.
    tau_abs <- fixed(1e-7)  ; label("Onset-time regulariser tau of the fractal release rate (h)")          # Code S5 $PK "TAU = 0.0000001" (a hardcoded constant, not a $THETA)

    # IIV. Code S5 $OMEGA holds variances of log-normal etas; the CV% printed in
    # Table S2 is recovered as sqrt(exp(omega) - 1). Note the $OMEGA order (Ka,
    # CL, Vc, Q, Kt, h) differs from the ETA() numbering used in $PK, where Q
    # takes ETA(4) and Kt takes ETA(5); the mapping below follows $PK.
    etalka   ~ 0.236701     # Code S5 $OMEGA 1 (ETA(1) on Ka); Table S2 fractal Ka IIV 51.7% CV [Shr 34.24%]
    etalcl   ~ 0.043698     # Code S5 $OMEGA 2 (ETA(2) on CL); Table S2 fractal CL IIV 21.1% CV [Shr 0.49%]
    etalvc   ~ 0.162774     # Code S5 $OMEGA 3 (ETA(3) on Vc); Table S2 fractal Vc IIV 42.1% CV [Shr 12.17%]
    etalq    ~ 0.00557957   # Code S5 $OMEGA 4 (ETA(4) on Q);  Table S2 fractal Q IIV 7.48% CV [Shr 83.55%]
    etalktr  ~ 0.105941     # Code S5 $OMEGA 5 (ETA(5) on Kt); Table S2 fractal Kt IIV 33.4% CV [Shr 11.67%]
    etah_abs ~ 0.0269217    # Code S5 $OMEGA 6 (ETA(6) on h);  Table S2 fractal h IIV 16.5% CV [Shr 13%]

    # Residual error. Code S5 $ERROR: W = SQRT(THETA(9)^2 + THETA(10)^2 * IPRED^2),
    # i.e. combined additive + proportional on the ng/mL scale, with $SIGMA 1 FIX.
    addSd  <- 0.678153      ; label("Additive residual error (ng/mL)")                          # Code S5 $THETA 9 (0,0.678153); Table S2 fractal Add-error = 0.678 (RSE 31.9%)
    propSd <- 0.16358       ; label("Proportional residual error (fraction)")                   # Code S5 $THETA 10 (0,0.16358); Table S2 fractal Prop-error = 0.164 (RSE 7.7%)
  })

  model({
    # 1. Individual parameters. Code S5 $PK puts etas on Ka, CL, Vc, Q, Kt and h;
    #    Vp and the release rate Kf are typical-value only.
    ka   <- exp(lka   + etalka)
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp)
    q    <- exp(lq    + etalq)
    ktr  <- exp(lktr  + etalktr)
    krel <- exp(lkrel)

    # Code S5: H = THETA(6) * EXP(ETA(6)) -- a multiplicative log-normal eta on a
    # bare (untransformed) exponent, so h is reproduced as h_abs * exp(eta) rather
    # than through the usual exp(l... + eta) idiom.
    h <- h_abs * exp(etah_abs)

    # 2. Micro-constants
    kel <- cl / vc
    kcp <- q  / vc
    kpc <- q  / vp

    # 3. Fractal release rate from the patch reservoir, with the clock RESET at
    #    each dose. Code S5 $PK carries TD, the time of the most recent dosing
    #    record, and forms KF = EXP(LOG(THETA(8)) - H*LOG(TIME - TD + TAU)),
    #    i.e. Kf / (time since the last dose + tau)^h.
    #
    #    `tad(depot_td)` is the time since the last dose into the patch reservoir.
    #    In the control stream TD is reset by ANY dosing record (IF DOSENO.GT.0),
    #    including the oral ones, which is a broader trigger. The two are
    #    equivalent for this study's design: Jung 2023 Section 2.3 describes the
    #    oral titration and the patch doses as SEQUENTIAL (7 days of oral, then
    #    three patches), so no oral dose is given while the reservoir holds drug.
    #    tau keeps the coefficient finite at the instant of each dose.
    #
    #    The release FLUX is formed inside a guard rather than as a bare rate
    #    constant. Before the first patch dose rxode2's tad() is NA, and although
    #    the reservoir is empty then, NA * 0 is NA in C and would silently poison
    #    every downstream state for the whole solve. Guarding on depot_td > 0
    #    means tad() is evaluated only once the reservoir actually holds drug.
    #    (rxode2 5.1.7 has no tad0()/tafd0() zero-returning variant.)
    rel_td <- 0
    if (depot_td > 0) {
      rel_td <- krel * depot_td / (tad(depot_td) + tau_abs)^h
    }

    # 4. ODE system. Code S5 $MODEL / $DES, compartments in the published order:
    #    1 GUT, 2 CENT, 3 PERI, 4 SKIN, 5 DEPOT (transit).
    #
    #    NOTE on the paper's own labels: $MODEL names compartment 4 "SKIN" and
    #    compartment 5 "DEPOT", while the $DES comments call 4 "Skin
    #    (formulation)" and 5 "Depot, transit", and the Table S2 legend glosses Kf
    #    as the "rate from formulation to skin" -- which reverses the direction.
    #    The ODEs are unambiguous and are what is implemented here: Kf drains
    #    compartment 4 into compartment 5, and Kt drains 5 into central. Compartment
    #    4 is therefore the patch formulation reservoir (`depot_td`) and 5 the single
    #    transit compartment, exactly as in Case 1, which uses the identical $MODEL
    #    naming with two transits.
    d/dt(depot_oral)  <- -ka * depot_oral
    d/dt(central)     <-  ka * depot_oral + kpc * peripheral1 - kcp * central -
      kel * central + ktr * transit1
    d/dt(peripheral1) <-  kcp * central - kpc * peripheral1
    d/dt(depot_td)    <- -rel_td
    d/dt(transit1)    <-  rel_td - ktr * transit1

    # 5. Observation. Amounts are mg and vc is L, so central/vc is mg/L;
    #    1 mg/L = 1000 ng/mL. Code S5 $ERROR: IPRED = A(2)/(VC/1000).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
