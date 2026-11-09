Oggianu_2023_trazodone_human <- function() {
  description <- paste(
    "Two-compartment population PK model of oral trazodone (immediate-release",
    "tablet and 6% oral drops) with first-order absorption and a lag time in",
    "healthy adult volunteers, carrying the authors' translational writhing",
    "layer: the human plasma concentration, rescaled to the mouse unbound",
    "fraction, drives the mouse acetic-acid writhing PD model (stimulus,",
    "nine transit compartments, writhing compartment) with the mouse",
    "trazodone potency, as used to predict the percent inhibition of a pain",
    "stimulus at low human doses. The writhing layer is a cross-species",
    "projection, not a fit to human PD data."
  )
  reference <- paste(
    "Oggianu L, Ottaviani M, Maddalena F, Leoni M, Pellizzoni G,",
    "Gaudesi E, Piccinini L, Comerio M, Furlan R, Salvatore C, Tongiani S,",
    "Calisti F. PK/PD analysis of trazodone and gabapentin in neuropathic",
    "pain rodent models: Translational PK-PD modeling from nonclinical to",
    "clinical development. Clin Transl Sci. 2023;16(4):606-617.",
    "doi:10.1111/cts.13472. Human PK: Supplementary Table S2; translation",
    "method: Methods 'Prediction of human pharmacological effect'; writhing",
    "PD: Table 2, Table 4B and Supplementary Table S1."
  )
  vignette <- "Oggianu_2023_trazodone_gabapentin"
  units <- list(
    time = "h",
    dosing = "mg (trazodone); unit stimulus (acetic acid)",
    concentration = "ng/mL"
  )

  # The acetic-acid stimulus site and the writhing count are specific to
  # the writhing-test model; transit1-transit9 follow the registered
  # numbered-chain pattern.
  paper_specific_compartments <- c("stimulus", "writhing")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "trazodone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "trazodone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "trazodone", units = "mg", specimen = "tissue", verified = TRUE),
    stimulus = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit8 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit9 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    writhing = list(analyte = "writhing count", units = "writhings", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "22-54 years",
    disease_state = "healthy volunteers, fasted",
    dose_range = paste(
      "Single oral doses: 50 mg immediate-release tablet and 30, 60 and",
      "90 mg as 6% oral drops, in a randomized four-way crossover with",
      "7-day washouts."
    ),
    notes = paste(
      "Phase I study (reference 16 of the paper). Sex, weight and race of",
      "the volunteers are not reported in this paper. Formulation (tablet",
      "vs. drops) was investigated as a covariate; Supplementary Table S2",
      "reports no formulation effect, so the final model has none. The",
      "writhing parameters are the mouse estimates (CD-1 male mice) and",
      "were not fitted to human data."
    )
  )

  ini({
    # --- Trazodone PK in humans (Oggianu 2023 Supplementary Table S2) ----
    lcl <- log(9.86); label("Apparent clearance CL/F (L/h)") # Table S2 CL/F = 9.86 L/h
    lvc <- log(50.8); label("Apparent central volume V2/F (L)") # Table S2 V2/F = 50.8 L
    lq <- log(31.6); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S2 Q/F = 31.6 L/h
    lvp <- log(76.9); label("Apparent peripheral volume V3/F (L)") # Table S2 V3/F = 76.9 L
    lka <- log(3.30); label("First-order absorption rate constant ka (1/h)") # Table S2 ka = 3.30 1/h
    ltlag <- log(0.0787); label("Absorption lag time ALAG1 (h)") # Table S2 ALAG1 = 0.0787 h

    # --- Mouse writhing PD carried into the human projection -------------
    # The paper assumes the mouse PD parameters, including the trazodone
    # potency k2, hold in humans ('a similar k2 across species').
    lktr <- fixed(log(1.36)); label("Stimulus transit rate constant ktr (1/min)") # Table 2 ktr = 1.36 1/min (mouse)
    lstim <- fixed(log(26.10)); label("Acetic-acid stimulus amplitude S at injection (writhings)") # Table 2 S = 26.10 (mouse)
    lkout <- fixed(log(0.104)); label("First-order loss rate constant k of the writhing compartment (1/min)") # Table 2 k = 0.104 1/min (mouse)
    lslope <- fixed(log(10.9)); label("Trazodone potency k2 on stimulus loss (L/umol/min)") # Table 4B trazodone k2 = 10.9 L/umol/min (mouse)

    # --- IIV (variances; %CV in Table S2 is sqrt(omega^2)) ---------------
    etalcl ~ 0.0909 # Table S2 omega2 CL/F = 0.0909 (30.1 %CV)
    etalvc ~ 0.0813 # Table S2 omega2 V2/F = 0.0813 (28.5 %CV)
    etalq ~ 0.0854 # Table S2 omega2 Q/F = 0.0854 (29.2 %CV)
    etalvp ~ 0.0577 # Table S2 omega2 V3/F = 0.0577 (24.0 %CV)
    etalka ~ 0.454 # Table S2 omega2 ka = 0.454 (67.4 %CV)
    etaltlag ~ 0.00148 # Table S2 omega2 alag = 0.00148 (3.85 %CV)
    etalstim ~ fixed(0.0231) # Table 2 omega2 S = 0.0231 (mouse)
    etalkout ~ fixed(0.0349) # Table 2 omega2 k = 0.0349 (mouse)
    etalslope ~ fixed(1.42) # Table 4B trazodone omega2 k2 = 1.42 (mouse)

    # --- Residual error ---------------------------------------------------
    propSd <- sqrt(0.0634); label("Proportional residual error, trazodone plasma concentration (fraction)") # Table S2 sigma2 prop = 0.0634 (25.2 %CV)
  })

  model({
    # The writhing rate constants are per minute; the model clock is hours.
    minPerHour <- 60
    # Free-base molecular weight (g/mol): ng/mL divided by MW is umol/L.
    mwTrazodone <- 371.86
    # Plasma protein binding 97% in humans and 89% in mice (Methods,
    # 'Prediction of human pharmacological effect'); the human
    # concentration is rescaled to the same unbound concentration in mouse
    # plasma before it drives the mouse-fitted potency.
    fuHuman <- 1 - 0.97
    fuMouse <- 1 - 0.89

    # --- Trazodone PK -----------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose mg and volume L give mg/L; x1000 reports ng/mL.
    Cc <- 1000 * central / vc
    # Mouse-equivalent molar concentration (Table S3 'AUCeq').
    cumolEq <- Cc * (fuHuman / fuMouse) / mwTrazodone

    # --- Writhing PD (Supplementary Table S1B) ----------------------------
    ktr <- exp(lktr) * minPerHour
    stim <- exp(lstim + etalstim)
    kout <- exp(lkout + etalkout) * minPerHour
    slope <- exp(lslope + etalslope) * minPerHour

    d/dt(stimulus) <- -slope * cumolEq * stimulus - ktr * stimulus
    d/dt(transit1) <- ktr * stimulus - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ktr * transit8
    d/dt(transit9) <- ktr * transit8 - ktr * transit9
    d/dt(writhing) <- ktr * transit9 - kout * writhing

    # The stimulus is dosed as amt = 1 into `stimulus`, scaled to the
    # subject's amplitude S (the paper's a0(0) = S).
    f(stimulus) <- stim

    Cc ~ prop(propSd)
  })
}
