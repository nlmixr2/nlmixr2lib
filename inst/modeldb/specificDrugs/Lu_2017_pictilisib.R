Lu_2017_pictilisib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral pictilisib (GDC-0941, a",
    "pan-class I PI3K inhibitor) in 31 healthy volunteers given a single",
    "40 mg dose in a four-period crossover of fasted, high-fat-fed, and",
    "rabeprazole (proton-pump inhibitor) pretreatment conditions. First-order",
    "absorption with an absorption lag and linear elimination from the",
    "central compartment. Food and PPI both lower the absorption rate",
    "constant (log-additive effects); relative bioavailability is additive in",
    "PPI (lower) and food (higher), and its between-subject variability is",
    "estimated separately for each of the four food/PPI states."
  )
  reference <- paste(
    "Lu T, Fraczkiewicz G, Salphati L, Budha N, Dalziel G, Smelick GS,",
    "Morrissey KM, Davis JD, Jin JY, Ware JA. Combining 'Bottom-up' and",
    "'Top-down' Approaches to Assess the Impact of Food and Gastric pH on",
    "Pictilisib (GDC-0941) Pharmacokinetics. CPT Pharmacometrics Syst",
    "Pharmacol. 2017;6(11):747-755. doi:10.1002/psp4.12228"
  )
  vignette <- "Lu_2017_pictilisib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "pictilisib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pictilisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pictilisib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    FED = list(
      description = "Fed-versus-fasted state at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Per-period indicator. Fed means a standard high-fat meal eaten",
        "within 30 min of the pictilisib dose (Methods, 'Clinical",
        "pharmacology study'); the paper calls this state 'food' throughout",
        "and tests it as a general fasted/fed covariate, so FED is used",
        "rather than FED_HIGHFAT. Effects (supplement 7 control stream,",
        "Table 1): log-additive on KA (KA = exp(THETA(5) + THETA(7) * FOOD",
        "+ THETA(9) * PPI + ETA(6))), additive on relative bioavailability",
        "(F1 = THETA(8) + THETA(10) * PPI + THETA(11) * FOOD), and selects",
        "which of the four state-specific F1 etas applies. Crossover design,",
        "so each subject contributes both levels."
      ),
      source_name = "FOOD"
    ),
    CONMED_PPI = list(
      description = "Proton-pump inhibitor (rabeprazole) pretreatment",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI)",
      notes = paste(
        "Per-period indicator: 1 = 20 mg rabeprazole daily for 4 days before",
        "and coadministered with the single pictilisib dose (Methods,",
        "'Clinical pharmacology study'). Effects (supplement 7 control",
        "stream, Table 1): log-additive on KA and additive on relative",
        "bioavailability, and selects which of the four state-specific F1",
        "etas applies. Crossover design, so each subject contributes both",
        "levels."
      ),
      source_name = "PPI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    n_observations = 1202L,
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Single 40 mg oral pictilisib tablet dose in each of four periods",
      "(fasted, fed, fasted with PPI, fed with PPI)"
    ),
    study_design = paste(
      "Single-center, four-period, two-sequence, open-label, randomized",
      "crossover (NCT00999128). Sequence 1: fasted -> fed -> fasted/PPI ->",
      "fed/PPI; sequence 2: fed -> fasted -> fed/PPI -> fasted/PPI. 32",
      "subjects randomized, 31 dosed."
    ),
    co_medication = "Rabeprazole 20 mg once daily for 4 days before the pictilisib dose in the PPI periods",
    notes = paste(
      "Lu 2017 does not tabulate age, weight, sex or race for the 31",
      "volunteers (the clinical study is reported in Ware 2013, Mol Pharm",
      "10:4074-4081); the GastroPlus typical subject is described as a",
      "32-year-old, 71 kg male. No body-size or demographic covariate enters",
      "the final popPK model."
    )
  )

  ini({
    # Disposition (Lu 2017 Table 1). All are apparent values (/F) for the
    # fasted, no-PPI reference state, whose relative bioavailability is 1.
    lcl <- log(26.0); label("Apparent clearance CL/F, fasted without PPI (L/h)")                        # Table 1 'Clearance (CL/F)' = 26.0 L/h; footnote 'apparent clearance for fasted state without PPI'
    lvc <- log(304.9); label("Apparent central volume Vc/F (L)")                                       # Table 1 'Central volume (Vc/F)' = 304.9 L
    lq <- log(13.9); label("Apparent intercompartmental clearance Q/F (L/h)")                           # Table 1 'Distribution clearance (Q/F)' = 13.9 L/h
    lvp <- log(175.9); label("Apparent peripheral volume Vp/F (L)")                                     # Table 1 'Peripheral volume (Vp/F)' = 175.9 L

    # Absorption (Table 1; control stream supplement 7:
    # KA = EXP(THETA(5) + THETA(7)*FOOD + THETA(9)*PPI + ETA(6))).
    lka <- log(2.63); label("Absorption rate constant, fasted without PPI (1/h)")                        # Table 1 'Ka, PPI=0, Fasted' = 2.63 1/h
    e_fed_ka <- log(0.297 / 2.63); label("Log-additive food effect on ka (unitless)")                   # Table 1 'Ka, PPI=0, Fed' = 0.297 vs 'PPI=0, Fasted' = 2.63; exp(e) = 0.113; check: 1.361 * exp(e) = 0.154 = 'PPI=1, Fed'
    e_conmed_ppi_ka <- log(1.361 / 2.63); label("Log-additive PPI effect on ka (unitless)")             # Table 1 'Ka, PPI=1, Fasted' = 1.361 vs 'PPI=0, Fasted' = 2.63; exp(e) = 0.518
    ltlag <- log(0.47); label("Absorption lag time (h)")                                                 # Table 1 'Absorption lag time (Tlag)' = 0.47 h

    # Relative bioavailability, additive in the control stream:
    # F1 = (THETA(8) + THETA(10)*PPI + THETA(11)*FOOD) * EXP(FETA), THETA(8) = 1 FIX.
    lfdepot <- fixed(log(1)); label("Relative bioavailability, fasted without PPI (fraction)")          # Table 1 'F, PPI=0, Fasted' = 100%; supplement 7 '$THETA ... 1 FIX' (THETA(8))
    e_conmed_ppi_f <- -0.571; label("Additive PPI effect on relative bioavailability (fraction)")       # Table 1 'F, PPI=1, Fasted' = 42.9% -> 0.429 - 1
    e_fed_f <- 0.184; label("Additive food effect on relative bioavailability (fraction)")              # Table 1 'F, PPI=0, Fed' = 118.4% -> 1.184 - 1; check: 1 - 0.571 + 0.184 = 0.613 = 'PPI=1, Fed' 61.3%

    # Between-subject variability (Table 1 'BSV' column, read as
    # omega = BSV/100 on the log scale, i.e. omega^2 = (BSV/100)^2). Q/F and
    # Vp/F carry no eta (supplement 7 '$OMEGA 0 FIX 0 FIX').
    etalcl ~ 0.0529                   # Table 1 'CL/F' BSV 23.0% -> 0.230^2
    etalvc ~ 0.142129                 # Table 1 'Vc/F' BSV 37.7% -> 0.377^2
    etalka ~ 0.269361                 # Table 1 'Ka' BSV 51.9% (one eta shared by all four states) -> 0.519^2
    etaltlag ~ 0.0009                 # Table 1 'Tlag' BSV 3.00% -> 0.030^2
    # State-specific etas on F1 (supplement 7: FETA = ETA(7)*OCC1 + ETA(8)*OCC2 + ETA(9)*OCC3 + ETA(10)*OCC4)
    etalfdepot_fasted ~ 0.0729        # Table 1 'F, PPI=0, Fasted' BSV 27% -> 0.27^2 (ETA(7), OCC1)
    etalfdepot_fastedppi ~ 0.4096     # Table 1 'F, PPI=1, Fasted' BSV 64% -> 0.64^2 (ETA(8), OCC2)
    etalfdepot_fed ~ 0.007921         # Table 1 'F, PPI=0, Fed' BSV 8.9% -> 0.089^2 (ETA(9), OCC3)
    etalfdepot_fedppi ~ 0.1024        # Table 1 'F, PPI=1, Fed' BSV 32% -> 0.32^2 (ETA(10), OCC4)

    # Residual error: additive on log-transformed concentrations
    # (supplement 7 '$ERROR Y = LOG(F) + EPS(1)'), i.e. log-normal.
    expSd <- 0.561; label("Residual error SD on the log scale (unitless)")                               # Table 1 'Residual error (proportional)' = 56.1%
  })
  model({
    # Indicators of the four food/PPI states (supplement 7 OCC1-OCC4).
    st_fasted <- (1 - FED) * (1 - CONMED_PPI)
    st_fastedppi <- (1 - FED) * CONMED_PPI
    st_fed <- FED * (1 - CONMED_PPI)
    st_fedppi <- FED * CONMED_PPI
    eta_f <- etalfdepot_fasted * st_fasted + etalfdepot_fastedppi * st_fastedppi +
      etalfdepot_fed * st_fed + etalfdepot_fedppi * st_fedppi

    ka <- exp(lka + e_fed_ka * FED + e_conmed_ppi_ka * CONMED_PPI + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- (exp(lfdepot) + e_conmed_ppi_f * CONMED_PPI + e_fed_f * FED) * exp(eta_f)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and vc in L give mg/L; the control stream scales with
    # S2 = V2/1000, reporting ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
