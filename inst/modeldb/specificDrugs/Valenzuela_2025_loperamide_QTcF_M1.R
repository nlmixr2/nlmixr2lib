Valenzuela_2025_loperamide_QTcF_M1 <- function() {
  description <- paste(
    "Population PK / concentration-QTc model for N-desmethyl loperamide (M1),",
    "the major metabolite of loperamide, in healthy adults from a randomized,",
    "placebo- and moxifloxacin-controlled four-way crossover thorough-QT",
    "study (NCT04225078; single oral loperamide 8 mg and 48 mg). This is the",
    "paper's primary exposure-response model. Fit-for-purpose M1 PK with the",
    "metabolite treated as if dosed independently: the loperamide dose enters",
    "a depot by a zero-order input after a lag time and reaches a",
    "two-compartment M1 disposition model by first-order absorption; the 48",
    "mg dose has its own relative bioavailability, absorption rate constant",
    "and lag time. IIV on Vc/F, F and ka; inter-occasion variability on the",
    "lag time and the zero-order duration. An effect compartment (ke0) on the",
    "M1 concentration drives the prespecified linear model for the",
    "placebo-adjusted change from baseline in Fridericia-corrected QT",
    "(ddQTcF = (b0 + eta) + (b1 + eta) * Ce + b2 * (QTcF baseline - mean)).",
    "The loperamide-driven sensitivity model is",
    "Valenzuela_2025_loperamide_QTcF_parent.",
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
    dosing = "mg (loperamide dose)",
    concentration = "ng/mL (M1); observation QTcF is the placebo-adjusted change from baseline in the Fridericia-corrected QT interval (ddQTcF), in ms"
  )

  compartmentData <- list(
    depot = list(
      analyte = "loperamide (dose driving the apparent M1 input)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_m1 = list(analyte = "N-desmethyl loperamide (M1)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m1 = list(analyte = "N-desmethyl loperamide (M1)", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(
      analyte = "N-desmethyl loperamide (M1)",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    DOSE_HIGH = list(
      description = "Supratherapeutic-dose indicator: 1 = the 48 mg loperamide dose, 0 = the 8 mg therapeutic dose.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (8 mg)",
      notes = paste(
        "Selects the 48 mg typical values of the M1 relative bioavailability (1.10 vs 1.00 fixed), absorption rate constant (0.991 vs 0.258 1/h) and lag time (0.376 vs 0.162 h) (Valenzuela 2025 Table 1; Figure S3 footnote 'Covariate of dose on F, ka/F and Alag').",
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
      "51 M1 (3.5%) concentrations were below the quantification limit; their effect-compartment concentrations were predicted from the individual PK parameters."
    )
  )

  ini({
    # ---- M1 PK (Valenzuela 2025 Table 1, 'M1' columns) ----
    lcl_m1 <- log(52.4)
    label("Apparent M1 clearance CL/F (L/h)") # Table 1 'CL/F (L/h)' M1 = 52.4 (SE 5.12%)
    lvc_m1 <- log(1650)
    label("Apparent M1 central volume Vc/F (L)") # Table 1 'Vc/F (L)' M1 = 1650 (SE 5.98%)
    lvp_m1 <- log(805)
    label("Apparent M1 peripheral volume Vp/F (L)") # Table 1 'Vp/F (L)' M1 = 805 (SE 11.2%)
    lq_m1 <- log(96.4)
    label("Apparent M1 intercompartmental clearance Q/F (L/h)") # Table 1 'Q/F (L/h)' M1 = 96.4 (SE 14.7%)

    # Input of M1 from the loperamide dose (metabolite treated as if dosed).
    lfdepot <- fixed(log(1))
    label("Relative bioavailability of the M1 input, 8 mg dose (fraction)") # Table 1 'F 8 mg' M1 = 1.00 Fixed
    lfdepot_highdose <- log(1.10)
    label("Relative bioavailability of the M1 input, 48 mg dose (fraction)") # Table 1 'F 48 mg' M1 = 1.10 (SE 3.13%)
    lka <- log(0.258)
    label("First-order rate constant of the M1 input, 8 mg dose (1/h)") # Table 1 'ka 8 mg (h-1)' M1 = 0.258 (SE 7.98%)
    lka_highdose <- log(0.991)
    label("First-order rate constant of the M1 input, 48 mg dose (1/h)") # Table 1 'ka 48 mg (h-1)' M1 = 0.991 (SE 12.9%)
    ltlag <- log(0.162)
    label("Lag time of the M1 input, 8 mg dose (h)") # Table 1 'Alag 8 mg (h)' M1 = 0.162 (SE 12.0%)
    ltlag_highdose <- log(0.376)
    label("Lag time of the M1 input, 48 mg dose (h)") # Table 1 'Alag 48 mg (h)' M1 = 0.376 (SE 6.86%)
    ld1 <- log(0.714)
    label("Duration of the zero-order M1 input D1 (h)") # Table 1 'D1 (h)' M1 = 0.714 (SE 10.7%)

    # IIV and IOV are printed as %CV; variance = log(1 + CV^2).
    etalvc_m1 ~ 0.0992379 # Table 1 'Vc/F' M1 IIV = 32.3%
    etalfdepot ~ 0.0423543 # Table 1 'F 8 mg' row M1 IIV = 20.8% (applies to F at both doses)
    etalka ~ 0.15258 # Table 1 'ka 8 mg' row M1 IIV = 40.6% (applies to ka at both doses)
    etaiov_ltlag_1 ~ 0.301325 # Table 1 'Alag 8 mg' row M1 IOV (footnote a) = 59.3%
    etaiov_ltlag_2 ~ fixed(0.301325) # same variance as occasion 1
    etaiov_ld1_1 ~ 0.664151 # Table 1 'D1' row M1 IOV (footnote a) = 97.1%
    etaiov_ld1_2 ~ fixed(0.664151) # same variance as occasion 1

    propSd_m1 <- 0.171
    label("Residual error on log-concentration, M1 (SD)") # Table 1 'Residual error' M1 = 0.171, footnote b 'Additive in log scale'

    # ---- M1 concentration-ddQTcF model (Table 2, 'M1 Model') ----
    lke0 <- log(0.094)
    label("Effect-compartment equilibration rate constant ke0 (1/h)") # Table 2 'ke0 (h-1)' = 0.094 (RSE 10.3%), estimated on the log scale (footnote a)
    e0 <- -1.66
    label("Intercept of the linear concentration-ddQTcF model (ms)") # Table 2 'Intercept (msec)' = -1.66 (RSE 36.1%)
    slope <- 0.544
    label("Slope of ddQTcF on M1 effect-compartment concentration (ms per ng/mL)") # Table 2 'Slope (msec/ng.mL)' = 0.544 (RSE 16.7%)
    e_qtc_bl_e0 <- -0.155
    label("Effect of centered baseline QTcF on the intercept (ms per ms)") # Table 2 'Effect centered baseline (msec)' = -0.155 (RSE 34.1%)

    etae0 ~ 15.7 # Table 2 'omega2 INT (add)' = 15.7 (variance, ms^2)
    etaslope ~ 0.143 # Table 2 'omega2 slope (add)' = 0.143 (variance)
    etalke0 ~ 0.748 # Table 2 'omega2 ke0 (exp)' = 0.748

    addSd <- 6.45755
    label("Additive residual error on ddQTcF (ms)") # Table 2 'sigma2 (add)' = 41.7 ms^2; SD = sqrt(41.7)
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

    cl_m1 <- exp(lcl_m1)
    vc_m1 <- exp(lvc_m1 + etalvc_m1)
    vp_m1 <- exp(lvp_m1)
    q_m1 <- exp(lq_m1)

    kel_m1 <- cl_m1 / vc_m1
    k12_m1 <- q_m1 / vc_m1
    k21_m1 <- q_m1 / vp_m1

    ke0 <- exp(lke0 + etalke0)

    # Sequential zero- then first-order input (Figure S3), with M1 treated as
    # if dosed independently (Results): the loperamide dose (mg) enters the
    # depot at a constant rate over d1 after the lag, and reaches the M1
    # central compartment by first-order absorption. The apparent parameters
    # absorb the formation fraction and any molecular-weight difference.
    # Dose records need rate = -2.
    d/dt(depot) <- -ka * depot
    d/dt(central_m1) <- ka * depot - kel_m1 * central_m1 - k12_m1 * central_m1 + k21_m1 * peripheral1_m1
    d/dt(peripheral1_m1) <- k12_m1 * central_m1 - k21_m1 * peripheral1_m1

    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1

    # Dose in mg, volumes in L: mg/L * 1000 = ng/mL.
    Cc_m1 <- central_m1 / vc_m1 * 1000

    # Effect compartment in concentration units (ng/mL M1).
    d/dt(effect) <- ke0 * (Cc_m1 - effect)

    # Prespecified linear model (Equation 1); placebo-adjusted change from
    # baseline in QTcF. Baseline QTcF centered on 400 ms (mean not printed).
    qtc_bl_centered <- QTC_BL - 400
    QTcF <- (e0 + etae0) + (slope + etaslope) * effect + e_qtc_bl_e0 * qtc_bl_centered

    Cc_m1 ~ prop(propSd_m1)
    QTcF ~ add(addSd)
  })
}
