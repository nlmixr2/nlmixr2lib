Kratochwil_2021_petesicatib <- function() {
  description <- "One-compartment population PK model for oral petesicatib (RO5459072, a cathepsin-S inhibitor) in healthy adults, with dose-dependent proximal-intestine bioavailability, a lagged first-order distal-intestine (colonic) absorption route in the fasted state, and a longer transit chain with dose-independent bioavailability in the fed state"
  reference <- paste(
    "Kratochwil NA, Stillhart C, Diack C, Nagel S, Al Kotbi N, Frey N.",
    "Population pharmacokinetic analysis of RO5459072, a low water-soluble",
    "drug exhibiting complex food-drug interactions.",
    "Br J Clin Pharmacol. 2021;87(9):3550-3560.",
    "doi:10.1111/bcp.14771",
    sep = " "
  )
  vignette <- "Kratochwil_2021_petesicatib"
  paper_specific_residual_sds <- c("expSd", "expSdEarlyFed")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FED = list(
      description = "Fed-state-at-dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, overnight fast of at least 8 h)",
      notes = paste(
        "1 = dose administered with food, 0 = fasted. Source column FOOD in",
        "the supplementary NONMEM control stream. Switches the whole",
        "absorption structure: fasted doses use a depot + 3-transit chain",
        "(plus the lagged distal route at doses >= 10 mg); fed doses use a",
        "depot + 9-transit chain, a dose-independent bioavailability of",
        "1.18 and no distal route. The fed arms were the SAD 100 mg dose",
        "after a high-fat high-calorie breakfast and the MAD 50/100/200 mg",
        "first doses given with food (meal composition not specified), so",
        "the general FED indicator is used rather than FED_HIGHFAT. FED is",
        "treated as constant within a subject: in the analysis each SAD",
        "cross-over period was a separate subject, and the stream switches",
        "its $DES on FOOD, so a within-subject change of FED is outside the",
        "fitted model. FED also selects the larger early residual error for",
        "the first 2 h after the first dose."
      ),
      source_name = "FOOD"
    ),
    DOSE = list(
      description = "Administered petesicatib dose level for the current dose record",
      units = "mg",
      type = "continuous",
      reference_category = "10 mg (fasted proximal bioavailability fixed to 1 at 10 mg)",
      notes = paste(
        "Per-record dose in mg (source column DOSE). Drives the fasted",
        "bioavailability: 1 mg and 3 mg each have an estimated proximal",
        "bioavailability and no distal absorption; doses >= 10 mg use the",
        "Emax form FFaPI = 1 - (D-10)^h / ((D-10)^h + (D50-10)^h) with",
        "h = 1 and distal bioavailability FFaDI = FFaPI * (1 - FFaPI)",
        "(Kratochwil 2021 Equations 1-2; supplementary control stream). It",
        "also selects the fasted absorption rate (2.05 1/h at >= 10 mg,",
        "2.95 1/h below 10 mg). Unused when FED = 1. Only 1, 3 and",
        "10-600 mg were studied fasted; the stream defines the low-dose",
        "branch only at exactly 1 and 3 mg, so doses between 3 and 10 mg",
        "are outside the model (here DOSE < 2 takes the 1 mg value and",
        "2 <= DOSE < 10 the 3 mg value). Supply DOSE on the dose rows."
      ),
      source_name = "DOSE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit9 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "petesicatib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "petesicatib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 39L,
    n_studies = 2L,
    n_observations = 816L,
    age_range = "adults (not reported)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "single oral doses of 1, 3, 10, 30, 100, 300 and 600 mg fasted and",
      "100 mg fed (SAD, NCT02295332); first doses of 50, 100 and 200 mg",
      "with food (MAD, NCT02521610)"
    ),
    regions = "Netherlands (PRA Health Sciences, Groningen)",
    notes = paste(
      "17 SAD and 22 MAD healthy male and female volunteers; the dataset",
      "was the first-dose PK up to 48-50 h. Each SAD cross-over period was",
      "treated as a separate subject. About 10% of samples were below the",
      "1 ng/mL LLOQ and were excluded. Formulation: micronized drug",
      "substance in hard-gelatine capsules without excipients. Demographic",
      "summaries (age, weight, sex split) are not reported in the paper or",
      "its supplement (Kratochwil 2021 Methods 2.1-2.3, Results 3.1)."
    )
  )

  ini({
    # Disposition (fasted and fed)
    lcl <- log(9.38); label("Apparent clearance CL/F (L/h)") # Table 1 'CL/F' = 9.38 (3.65% RSE); stream THETA(2)
    lvc <- log(109); label("Apparent central volume Vc/F (L)") # Table 1 'Vc/F' = 109 (4.76% RSE); stream THETA(1)

    # Proximal-intestine absorption (depot -> transit chain -> central)
    lka_fastge10 <- log(2.05); label("Proximal absorption / transit rate constant, fasted doses >= 10 mg (1/h)") # Table 1 'kaFaPI for doses > 3 mg' = 2.05 (6.29% RSE); stream THETA(3), applied at DOSE >= 10
    lka_fastlt10fed <- log(2.95); label("Proximal absorption / transit rate constant, fasted doses < 10 mg and all fed doses (1/h)") # Table 1 'kaFa for doses of 1 and 3 mg while fasted and kaFe while fed' = 2.95 (3.12% RSE); stream THETA(4)

    # Distal-intestine (colonic) absorption, fasted doses >= 10 mg only
    lka2 <- log(0.065); label("Distal-intestine first-order absorption rate constant kaFaDI (1/h)") # Table 1 'kaFaDI for doses > 3 mg' = 0.065 (9.31% RSE); stream THETA(5) KA4
    ltlag2 <- log(11.2); label("Lag time of the distal-intestine absorption LagFaDI (h)") # Table 1 'LagFaDI for doses > 3 mg' = 11.2 (6.92% RSE); stream THETA(6) ALAG4

    # Apparent bioavailability (relative to fasted 10 mg = 1)
    led50 <- log(216); label("D50: dose giving half-maximal reduction of fasted proximal bioavailability (mg)") # Table 1 'D50' = 216 mg (13.4% RSE); stream THETA(7)
    hill_fdepot <- fixed(1); label("Hill coefficient of the dose-Emax function on fasted proximal bioavailability (unitless)") # stream THETA(12) HILL = 1 FIX; Equation 1 printed without an exponent
    lfdepot_fast1mg <- log(0.486); label("Proximal bioavailability FFaPI, fasted 1 mg dose (fraction of 10 mg)") # Table 1 'FFaPI for a dose of 1 mg' = 0.486 (8.77% RSE); stream THETA(9)
    lfdepot_fast3mg <- log(0.735); label("Proximal bioavailability FFaPI, fasted 3 mg dose (fraction of 10 mg)") # Table 1 'FFaPI for a dose of 3 mg' = 0.735 (7.85% RSE); stream THETA(10)
    lfdepot_fed <- log(1.18); label("Apparent bioavailability FFe, fed state, dose independent (fraction of fasted 10 mg)") # Table 1 'FFe' = 1.18 (4.70% RSE); stream THETA(11)

    # IIV (variances from the supplementary control stream $OMEGA; Table 1
    # reports the same values as CV% = sqrt(exp(omega^2) - 1)). The CL-Vc
    # covariance is in the stream's $OMEGA BLOCK(2) but not in Table 1.
    etalcl + etalvc ~ c(0.0171, 0.0201, 0.0561) # stream $OMEGA BLOCK(2): Vc 5.61E-2, cov 2.01E-2, CL 1.71E-2; Table 1 CV 13.1% (CL/F), 24.0% (Vc/F)
    etalka_fastge10 ~ 0.108 # stream $OMEGA ETA(3) = 0.108; Table 1 'kaFaPI for doses > 3 mg' CV 33.8%
    etalka_fastlt10fed ~ 0.143 # stream $OMEGA ETA(5) = 0.143; Table 1 'kaFe' CV 39.2%; applied to fed doses only
    etalfdepot2 ~ 0.0979 # stream $OMEGA ETA(6) = 9.79E-2 on F4; Table 1 'FFaDI' CV 32.1%
    etaled50 ~ 0.138 # stream $OMEGA ETA(7) = 0.138; Table 1 'D50' CV 38.5%

    # Residual error: log-transform-both-sides, Y = log(F) + W * EPS, SIGMA 1 FIX
    expSd <- 0.244; label("Log-scale residual SD, fasted data and fed data > 2 h after the first dose") # Table 1 'Proportional (fasted and time after first dose > 2 hr in fed condition)' = 24.4%; stream THETA(8) W1
    expSdEarlyFed <- 1.99; label("Log-scale residual SD, fed data <= 2 h after the first dose") # Table 1 'Proportional (time after first dose <= 2 hr in fed condition)' = 199%; stream THETA(13) W2
  })

  model({
    # ---- 1. Covariate-derived switches ----
    fasted_high <- (1 - FED) * (DOSE >= 10)
    fasted_low <- (1 - FED) * (DOSE < 10)

    # ---- 2. Individual parameters ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # Proximal chain rate: one rate constant for every step of the chain
    # (the stream's DADT uses KA1 / KA3 / KAF for depot, transit and
    # transit-to-central steps alike). The kaFe eta applies to fed doses
    # only; the fasted < 10 mg rate is the typical value (stream KA3 = TVKA3).
    ka_high <- exp(lka_fastge10 + etalka_fastge10)
    ka_fed <- exp(lka_fastlt10fed + etalka_fastlt10fed)
    ka_low <- exp(lka_fastlt10fed)
    ka <- fasted_high * ka_high + fasted_low * ka_low + FED * ka_fed

    ka2 <- exp(lka2)
    tlag2 <- exp(ltlag2)

    # Fasted dose >= 10 mg: Emax reduction of proximal bioavailability,
    # referenced to 10 mg (Equation 1), and distal bioavailability equal to
    # FFaPI times the fraction not absorbed proximally (Equation 2).
    d50 <- exp(led50 + etaled50)
    dose_above <- max(DOSE - 10, 0)
    frac_dose <- dose_above^hill_fdepot / (dose_above^hill_fdepot + (d50 - 10)^hill_fdepot)
    f_prox_high <- 1 - frac_dose
    f_dist_high <- f_prox_high * (1 - f_prox_high) * exp(etalfdepot2)

    # Fasted 1 mg and 3 mg: separately estimated proximal bioavailability,
    # no distal absorption.
    f_prox_low <- exp(lfdepot_fast1mg) * (DOSE < 2) + exp(lfdepot_fast3mg) * (DOSE >= 2)

    f_prox <- fasted_high * f_prox_high + fasted_low * f_prox_low + FED * exp(lfdepot_fed)
    f_dist <- fasted_high * f_dist_high

    # ---- 3. ODEs ----
    # Proximal intestine: depot + 3 transit compartments when fasted; the
    # fed chain continues through transit4-transit9 (depot + 9 transits).
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(transit4) <- FED * ka * transit3 - ka * transit4
    d/dt(transit5) <- ka * transit4 - ka * transit5
    d/dt(transit6) <- ka * transit5 - ka * transit6
    d/dt(transit7) <- ka * transit6 - ka * transit7
    d/dt(transit8) <- ka * transit7 - ka * transit8
    d/dt(transit9) <- ka * transit8 - ka * transit9
    # Distal intestine (colon), fasted doses >= 10 mg; dosed in parallel.
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- (1 - FED) * ka * transit3 + ka * transit9 + ka2 * depot2 - kel * central

    # ---- 4. Bioavailability and lag ----
    # Every dose is given twice, as in the source dataset: once to depot and
    # once to depot2 (NONMEM compartments 1 and 4, F1 and F4).
    f(depot) <- f_prox
    f(depot2) <- f_dist
    alag(depot2) <- tlag2

    # ---- 5. Observation and time-dependent residual error ----
    # Central amount in mg, volume in L: mg/L * 1000 = ng/mL.
    Cc <- central / vc * 1000

    t_first_dose <- tafd()
    early_fed <- FED * (t_first_dose <= 2)
    expSdObs <- expSd * (1 - early_fed) + expSdEarlyFed * early_fed
    Cc ~ lnorm(expSdObs)
  })
}
