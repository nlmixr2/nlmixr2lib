Valenzuela_2025_loperamide_QTcF_parent <- function() {
  description <- paste(
    "Population PK / concentration-QTc model for loperamide (parent) in",
    "healthy adults from a randomized, placebo- and moxifloxacin-controlled",
    "four-way crossover thorough-QT study (NCT04225078; single oral",
    "loperamide 8 mg and 48 mg). Fit-for-purpose PK: two-compartment",
    "disposition with linear elimination and sequential zero- then",
    "first-order absorption after a lag time; the 48 mg dose has its own",
    "relative bioavailability, absorption rate constant and lag time. IIV on",
    "Vc/F, F and ka; inter-occasion variability on the lag time and the",
    "zero-order duration. An effect compartment (ke0) drives the",
    "prespecified linear model for the placebo-adjusted change from baseline",
    "in Fridericia-corrected QT (ddQTcF = (b0 + eta) + (b1 + eta) * Ce +",
    "b2 * (QTcF baseline - mean)). This loperamide-driven exposure-response",
    "model is the paper's sensitivity analysis; the primary model is driven",
    "by the N-desmethyl metabolite M1 (Valenzuela_2025_loperamide_QTcF_M1).",
    sep = " "
  )

  reference <- paste(
    "Valenzuela B, Gisleskog PO, Cirillo I, Coenen E, Ariyawansa J, Ali SR,",
    "Takhtoukh S, Perez-Ruixo JJ, Ackaert O. Evaluation of the Effect of",
    "Loperamide on the Cardiac Repolarization Interval Using",
    "Exposure-Response Analysis. Clin Transl Sci. 2025;18(2):e70114.",
    "doi:10.1111/cts.70114.",
    sep = " "
  )

  vignette <- "Valenzuela_2025_loperamide"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (loperamide); observation QTcF is the placebo-adjusted change from baseline in the Fridericia-corrected QT interval (ddQTcF), in ms"
  )

  compartmentData <- list(
    depot = list(analyte = "loperamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "loperamide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "loperamide", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "loperamide", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    DOSE_HIGH = list(
      description = "Supratherapeutic-dose indicator: 1 = the 48 mg loperamide dose, 0 = the 8 mg therapeutic dose.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (8 mg)",
      notes = paste(
        "Selects the 48 mg typical values of relative bioavailability (1.20 vs 1.00 fixed), absorption rate constant (4.21 vs 1.19 1/h) and lag time (0.271 vs 0.149 h) (Valenzuela 2025 Table 1; Figure S3 footnote 'Covariate of dose on F, ka/F and Alag').",
        "Only the 8 and 48 mg single doses were studied, so the indicator names a dose level rather than binning a continuous dose; doses between or outside them are extrapolation.",
        "Must be set on the dose record (f() and alag() are evaluated there) and on every observation row of the same period."
      ),
      source_name = "dose (8 mg vs 48 mg)"
    ),
    OCC = list(
      description = "Occasion index for the inter-occasion variability on the absorption lag time and the zero-order absorption duration: 1 and 2 = the two loperamide treatment periods of the crossover.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Valenzuela 2025 Table 1 footnote a marks the Alag and D1 variabilities as inter-occasion. Each participant received loperamide in two periods (8 mg and 48 mg); placebo and moxifloxacin periods carry no loperamide.",
        "Two occasions are encoded (OCC = 1, 2), sharing one variance per parameter. Any other value (for example 0) switches the inter-occasion random effects off, which is the typical-value prediction.",
        "Set on the dose record, where alag() and dur() are evaluated."
      ),
      source_name = "occasion (treatment period)"
    ),
    QTC_BL = list(
      description = "Participant's pre-dose baseline Fridericia-corrected QT interval for the treatment period (QTcBas). Enters the ddQTcF intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 400).",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per participant and per period: 'For each participant and treatment, the average of three sets of triplicate pre-dose ECGs was used as the baseline for that treatment period' (Methods, ECG data collection). Fridericia correction.",
        "Centering value: Equation 1 centers on the population mean baseline QTcF, which the paper does not print. 400 ms is used as a rounded standard for healthy adults (the same value as Wang_2022_abrocitinib_QTcF). Set QTC_BL = 400 for a typical participant, which removes the term."
      ),
      source_name = "QTcBas"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 53L,
    n_studies = 1L,
    n_observations = 1408L,
    age_range = "18-55 years (eligibility criterion; observed range not reported)",
    weight_range = ">= 50 kg; BMI 18-30 kg/m^2 (eligibility criteria)",
    sex_female_pct = 66.2,
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "Single oral loperamide 8 mg (therapeutic; 51 participants) and 48 mg",
      "(supratherapeutic; 50 participants), placebo and moxifloxacin 400 mg",
      "in a randomized double-blind four-way crossover with >= 7-day washouts."
    ),
    regions = "Belgium",
    notes = paste(
      "Thorough-QT study NCT04225078, conducted January 2020 to January 2022.",
      "65 participants were enrolled (22 male, 43 female; sex_female_pct is for the enrolled cohort); the exposure-response dataset holds 1408 time-matched PK and ddQTcF records from 53 participants (Results, Dataset).",
      "ECGs were extracted from continuous 12-lead Holter recordings in triplicate at 0.5, 1, 2, 3.5, 4, 5, 5.5, 6, 6.5, 7, 8, 12, 24 and 36 h post dose, time-matched to PK samples.",
      "34 loperamide (2.4%) concentrations were below the quantification limit; their effect-compartment concentrations were predicted from the individual PK parameters."
    )
  )

  ini({
    # ---- Loperamide PK (Valenzuela 2025 Table 1, 'Loperamide' columns) ----
    lcl <- log(293)
    label("Apparent clearance CL/F (L/h)") # Table 1 'CL/F (L/h)' = 293 (SE 5.3%)
    lvc <- log(3380)
    label("Apparent central volume Vc/F (L)") # Table 1 'Vc/F (L)' = 3380 (SE 6.8%)
    lvp <- log(1770)
    label("Apparent peripheral volume Vp/F (L)") # Table 1 'Vp/F (L)' = 1770 (SE 6.2%)
    lq <- log(219)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Q/F (L/h)' = 219 (SE 10.9%)

    lfdepot <- fixed(log(1))
    label("Relative bioavailability, 8 mg dose (fraction)") # Table 1 'F 8 mg' = 1.00 Fixed
    lfdepot_highdose <- log(1.20)
    label("Relative bioavailability, 48 mg dose (fraction)") # Table 1 'F 48 mg' = 1.20 (SE 3.26%)
    lka <- log(1.19)
    label("First-order absorption rate constant, 8 mg dose (1/h)") # Table 1 'ka 8 mg (h-1)' = 1.19 (SE 10.2%)
    lka_highdose <- log(4.21)
    label("First-order absorption rate constant, 48 mg dose (1/h)") # Table 1 'ka 48 mg (h-1)' = 4.21 (SE 12.1%)
    ltlag <- log(0.149)
    label("Absorption lag time, 8 mg dose (h)") # Table 1 'Alag 8 mg (h)' = 0.149 (SE 14.7%)
    ltlag_highdose <- log(0.271)
    label("Absorption lag time, 48 mg dose (h)") # Table 1 'Alag 48 mg (h)' = 0.271 (SE 10.2%)
    ld1 <- log(0.551)
    label("Duration of the zero-order absorption input D1 (h)") # Table 1 'D1 (h)' = 0.551 (SE 4.34%)

    # IIV and IOV are printed as %CV; variance = log(1 + CV^2).
    etalvc ~ 0.0376965 # Table 1 'Vc/F' IIV = 19.6%
    etalfdepot ~ 0.173302 # Table 1 'F 8 mg' row IIV = 43.5% (applies to F at both doses)
    etalka ~ 0.679148 # Table 1 'ka 8 mg' row IIV = 98.6% (applies to ka at both doses)
    etaiov_ltlag_1 ~ 0.457832 # Table 1 'Alag 8 mg' row IOV (footnote a) = 76.2%
    etaiov_ltlag_2 ~ fixed(0.457832) # same variance as occasion 1
    etaiov_ld1_1 ~ 0.426209 # Table 1 'D1' row IOV (footnote a) = 72.9%
    etaiov_ld1_2 ~ fixed(0.426209) # same variance as occasion 1

    propSd <- 0.177
    label("Residual error on log-concentration, loperamide (SD)") # Table 1 'Residual error' = 0.177, footnote b 'Additive in log scale'

    # ---- Loperamide concentration-ddQTcF model (Table 2, 'Loperamide Model') ----
    lke0 <- log(0.024)
    label("Effect-compartment equilibration rate constant ke0 (1/h)") # Table 2 'ke0 (h-1)' = 0.024 (RSE 6.8%), estimated on the log scale (footnote a)
    e0 <- -1.76
    label("Intercept of the linear concentration-ddQTcF model (ms)") # Table 2 'Intercept (msec)' = -1.76 (RSE 35.2%)
    slope <- 2.91
    label("Slope of ddQTcF on loperamide effect-compartment concentration (ms per ng/mL)") # Table 2 'Slope (msec/ng.mL)' = 2.91 (RSE 13.7%)
    e_qtc_bl_e0 <- -0.173
    label("Effect of centered baseline QTcF on the intercept (ms per ms)") # Table 2 'Effect centered baseline (msec)' = -0.173 (RSE 31.8%)

    etae0 ~ 17.7 # Table 2 'omega2 INT (add)' = 17.7 (variance, ms^2)
    etaslope ~ fixed(0.00001) # Table 2 'omega2 slope (add)' = 0.00001, no RSE reported (held at a negligible value)
    etalke0 ~ 0.947 # Table 2 'omega2 ke0 (exp)' = 0.947

    addSd <- 6.54217
    label("Additive residual error on ddQTcF (ms)") # Table 2 'sigma2 (add)' = 42.8 ms^2; SD = sqrt(42.8)
  })

  model({
    # Occasion indicators for the inter-occasion variability (two
    # loperamide periods per participant).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ltlag <- oc1 * etaiov_ltlag_1 + oc2 * etaiov_ltlag_2
    iov_ld1 <- oc1 * etaiov_ld1_1 + oc2 * etaiov_ld1_2

    # The 48 mg dose carries its own typical F, ka and lag time (Table 1):
    # lx + (lx_highdose - lx) * DOSE_HIGH returns exp(lx) at 8 mg and
    # exp(lx_highdose) at 48 mg.
    fdepot <- exp(lfdepot + (lfdepot_highdose - lfdepot) * DOSE_HIGH + etalfdepot)
    ka <- exp(lka + (lka_highdose - lka) * DOSE_HIGH + etalka)
    tlag <- exp(ltlag + (ltlag_highdose - ltlag) * DOSE_HIGH + iov_ltlag)
    d1 <- exp(ld1 + iov_ld1)

    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    ke0 <- exp(lke0 + etalke0)

    # Sequential zero- then first-order absorption (Figure S3): after the lag,
    # the dose enters the depot at a constant rate over d1, and leaves it by
    # first-order absorption. Dose records need rate = -2.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1

    # Dose in mg, volumes in L: mg/L * 1000 = ng/mL.
    Cc <- central / vc * 1000

    # Effect compartment in concentration units (ng/mL).
    d/dt(effect) <- ke0 * (Cc - effect)

    # Prespecified linear model (Equation 1); placebo-adjusted change from
    # baseline in QTcF. Baseline QTcF centered on 400 ms (mean not printed).
    qtc_bl_centered <- QTC_BL - 400
    QTcF <- (e0 + etae0) + (slope + etaslope) * effect + e_qtc_bl_e0 * qtc_bl_centered

    Cc ~ prop(propSd)
    QTcF ~ add(addSd)
  })
}
