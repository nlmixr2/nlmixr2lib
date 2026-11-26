Ahn_2020_PF06648671 <- function() {
  description <- paste(
    "Population PK/PD model for the oral gamma-secretase modulator",
    "PF-06648671 and four cerebrospinal-fluid amyloid-beta species in",
    "healthy young and elderly adults. Two-compartment PK with first-order",
    "absorption and linear elimination drives four parallel indirect-response",
    "models on the zero-order production rate of each CSF amyloid-beta",
    "species: a hyperbolic Imax inhibition for Abeta42 and Abeta40 and a",
    "sigmoid Emax stimulation for Abeta38 and Abeta37. An empirical",
    "power-of-time placebo drift, (Pmax + eta) * time^Kp, is added to each",
    "observed species in the serial-lumbar-catheter study B7991003 only.",
    sep = " "
  )
  reference <- paste(
    "Ahn JE, Carrieri C, Dela Cruz F, Fullerton T, Hajos-Korcsok E, He P,",
    "Kantaridis C, Leurent C, Liu R, Mancuso J, Mendes da Costa L, Qiu R",
    "(2020). Pharmacokinetic and Pharmacodynamic Effects of a",
    "gamma-Secretase Modulator, PF-06648671, on CSF Amyloid-beta Peptides in",
    "Randomized Phase I Studies. Clin Pharmacol Ther 107(1):211-220.",
    "doi:10.1002/cpt.1570.",
    sep = " "
  )
  vignette <- "Ahn_2020_PF06648671"

  # CSF amyloid-beta turnover states, one per species, in the same
  # <species>_csf form as abeta_csf in Lin_2022_aducanumab_qsp.
  paper_specific_compartments <- c(
    "abeta42_csf",
    "abeta40_csf",
    "abeta38_csf",
    "abeta37_csf"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "PF-06648671", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "PF-06648671", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "PF-06648671", units = "mg", specimen = "tissue", verified = TRUE),
    abeta42_csf = list(analyte = "amyloid-beta 1-42", units = "pg/mL", specimen = "CSF", verified = TRUE),
    abeta40_csf = list(analyte = "amyloid-beta 1-40", units = "pg/mL", specimen = "CSF", verified = TRUE),
    abeta38_csf = list(analyte = "amyloid-beta 1-38", units = "pg/mL", specimen = "CSF", verified = TRUE),
    abeta37_csf = list(analyte = "amyloid-beta 1-37", units = "pg/mL", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    STUDY_B7991003 = list(
      description = "Indicator for the single-dose serial-CSF-sampling study B7991003 (1 = B7991003, 0 = any other study)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single lumbar puncture, e.g. B7991002)",
      notes = paste(
        "Switches on the empirical placebo drift of Eq. (1), which Ahn 2020",
        "applied to the single-dose serial-CSF data only ('Overall PD response",
        "was the sum of baseline, placebo (single-dose PD data), and treatment",
        "effects'). The drift is attributed to frequent CSF sampling through a",
        "lumbar catheter; the paper reports a negligible placebo effect with the",
        "single lumbar puncture of B7991002."
      ),
      source_name = "study B7991003"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 120,
    n_studies = 3,
    age_range = "18-75 years (part 2 of B7991002 enrolled 65-75-year-olds)",
    weight_range = "52.4-107.2 kg",
    sex_female_pct = 4.2,
    race_ethnicity = c(White = 84.2, Black = 10.0, Asian = 4.2, Other = 1.7),
    disease_state = "Healthy volunteers (young and elderly)",
    dose_range = "2-360 mg oral single dose (B7991001, B7991003); 4-360 mg oral once daily for 14 days (B7991002)",
    regions = "Belgium (B7991001, B7991002) and United States (B7991003)",
    notes = paste(
      "Three phase I studies (Ahn 2020 Table 1 and Results 'Study",
      "populations'): B7991001 single ascending dose (n = 18, all white men),",
      "B7991003 single-dose 150/300 mg with serial CSF sampling via lumbar",
      "catheter over 36 h (n = 22), and B7991002 14-day multiple ascending",
      "dose with a single lumbar puncture at baseline and day 14 (n = 70",
      "young, n = 10 elderly). Five subjects were women (one placebo subject",
      "in B7991003, four elderly subjects in B7991002). The PD model was fit",
      "to the CSF data of B7991003 and B7991002 (Figure S2)."
    )
  )

  ini({
    # Population PK (Table S4); PK parameters were held fixed while the PD
    # parameters were estimated (PPPD approach, Methods 'PK/PD modeling').
    lka <- log(1.49); label("First-order absorption rate constant ka (1/h)") # Table S4 'ka' = 1.49 /hr
    lcl <- log(3.56); label("Apparent clearance CL/F (L/h)") # Table S4 'CL/F' = 3.56 L/hr
    lvc <- log(46.3); label("Apparent central volume Vc/F (L)") # Table S4 'Vc/F' = 46.3 L
    lq <- log(9.17); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 'Q/F' = 9.17 L/hr
    lvp <- log(32.8); label("Apparent peripheral volume Vp/F (L)") # Table S4 'Vp/F' = 32.8 L

    # CSF Abeta42: production inhibited (Table 2, Eq. 2)
    lrbase_abeta42 <- log(1001); label("Baseline CSF Abeta42 (pg/mL)") # Table 2 Abeta42 'Baseline' = 1,001 pg/mL
    lkout_abeta42 <- log(0.121); label("First-order turnover rate constant of CSF Abeta42 kout (1/h)") # Table 2 Abeta42 'kout' = 0.121 /hour
    limax_abeta42 <- log(0.744); label("Maximum fractional inhibition of Abeta42 production Imax (unitless)") # Table 2 Abeta42 'Imax' = 0.744
    lic50_abeta42 <- log(421); label("Plasma PF-06648671 concentration at half Imax for Abeta42 IC50 (ng/mL)") # Table 2 Abeta42 'IC50' = 421 ng/mL

    # CSF Abeta40: production inhibited (Table 2, Eq. 2)
    lrbase_abeta40 <- log(7584); label("Baseline CSF Abeta40 (pg/mL)") # Table 2 Abeta40 'Baseline' = 7,584 pg/mL
    lkout_abeta40 <- log(0.0835); label("First-order turnover rate constant of CSF Abeta40 kout (1/h)") # Table 2 Abeta40 'kout' = 0.0835 /hour
    limax_abeta40 <- log(0.725); label("Maximum fractional inhibition of Abeta40 production Imax (unitless)") # Table 2 Abeta40 'Imax' = 0.725
    lic50_abeta40 <- log(1228); label("Plasma PF-06648671 concentration at half Imax for Abeta40 IC50 (ng/mL)") # Table 2 Abeta40 'IC50' = 1,228 ng/mL

    # CSF Abeta38: production stimulated (Table 2, Eq. 3)
    lrbase_abeta38 <- log(2404); label("Baseline CSF Abeta38 (pg/mL)") # Table 2 Abeta38 'Baseline' = 2,404 pg/mL
    lkout_abeta38 <- log(0.144); label("First-order turnover rate constant of CSF Abeta38 kout (1/h)") # Table 2 Abeta38 'kout' = 0.144 /hour
    lemax_abeta38 <- log(0.869); label("Maximum fractional stimulation of Abeta38 production Emax (unitless)") # Table 2 Abeta38 'Emax' = 0.869
    lec50_abeta38 <- log(1765); label("Plasma PF-06648671 concentration at half Emax for Abeta38 EC50 (ng/mL)") # Table 2 Abeta38 'EC50' = 1,765 ng/mL
    lhill_abeta38 <- log(1.4); label("Sigmoidicity of the Abeta38 stimulation gamma (unitless)") # Table 2 Abeta38 'Gamma' = 1.4

    # CSF Abeta37: production stimulated (Table 2, Eq. 3)
    lrbase_abeta37 <- log(702); label("Baseline CSF Abeta37 (pg/mL)") # Table 2 Abeta37 'Baseline' = 702 pg/mL
    lkout_abeta37 <- log(0.0964); label("First-order turnover rate constant of CSF Abeta37 kout (1/h)") # Table 2 Abeta37 'kout' = 0.0964 /hour
    lemax_abeta37 <- log(5.63); label("Maximum fractional stimulation of Abeta37 production Emax (unitless)") # Table 2 Abeta37 'Emax' = 5.63
    lec50_abeta37 <- log(1285); label("Plasma PF-06648671 concentration at half Emax for Abeta37 EC50 (ng/mL)") # Table 2 Abeta37 'EC50' = 1,285 ng/mL
    lhill_abeta37 <- log(1.42); label("Sigmoidicity of the Abeta37 stimulation gamma (unitless)") # Table 2 Abeta37 'Gamma' = 1.42

    # Empirical placebo drift, Eq. (1): (Pmax + eta) * time^Kp. Fit to the
    # B7991003 placebo data first and held fixed while the drug effects
    # were estimated (Methods). Pmax is additive on the pg/mL scale of the
    # observed species; time is hours after the dose.
    pmax_abeta42 <- fixed(87.4); label("Placebo drift magnitude Pmax for Abeta42 (pg/mL per h^Kp)") # Table S3 Abeta42 'Pmax' = 87.4
    pow_abeta42 <- fixed(0.517); label("Placebo drift time exponent Kp for Abeta42 (unitless)") # Table S3 Abeta42 'Kp' = 0.517
    pmax_abeta40 <- fixed(176); label("Placebo drift magnitude Pmax for Abeta40 (pg/mL per h^Kp)") # Table S3 Abeta40 'Pmax' = 176
    pow_abeta40 <- fixed(0.799); label("Placebo drift time exponent Kp for Abeta40 (unitless)") # Table S3 Abeta40 'Kp' = 0.799
    pmax_abeta38 <- fixed(70.7); label("Placebo drift magnitude Pmax for Abeta38 (pg/mL per h^Kp)") # Table S3 Abeta38 'Pmax' = 70.7
    pow_abeta38 <- fixed(0.707); label("Placebo drift time exponent Kp for Abeta38 (unitless)") # Table S3 Abeta38 'Kp' = 0.707
    pmax_abeta37 <- fixed(35.0); label("Placebo drift magnitude Pmax for Abeta37 (pg/mL per h^Kp)") # Table S3 Abeta37 'Pmax' = 35.0
    pow_abeta37 <- fixed(0.598); label("Placebo drift time exponent Kp for Abeta37 (unitless)") # Table S3 Abeta37 'Kp' = 0.598

    # PK IIV (Table S4 reports inter-individual VARIANCES)
    etalka ~ 0.189 # Table S4 'IIV on ka' = 0.189
    etalcl + etalvc ~ c(0.180, 0.0679, 0.0671) # Table S4 'IIV on CL/F' = 0.180, 'Covariance CL and Vc' = 0.0679, 'IIV on Vc/F' = 0.0671
    etalvp ~ 0.0867 # Table S4 'IIV on Vp/F' = 0.0867

    # Additive IIV on Pmax (Eq. 1); the variances are not reported
    # anywhere in the paper or supplement, so they are held at zero.
    etapmax_abeta42 ~ fixed(0) # Eq. 1 eta on Pmax; variance not reported
    etapmax_abeta40 ~ fixed(0) # Eq. 1 eta on Pmax; variance not reported
    etapmax_abeta38 ~ fixed(0) # Eq. 1 eta on Pmax; variance not reported
    etapmax_abeta37 ~ fixed(0) # Eq. 1 eta on Pmax; variance not reported

    # PK residual error (Table S4)
    propSd <- 0.164; label("Proportional residual error on plasma PF-06648671 (fraction)") # Table S4 'Proportional error, CV (%)' = 16.4
    addSd <- 0.462; label("Additive residual error on plasma PF-06648671 (ng/mL)") # Table S4 'Additive error, SD' = 0.462

    # PD residual error: magnitude and form are not reported; held at zero
    # so simulations of the CSF endpoints are prediction-only.
    propSd_abeta42 <- fixed(0); label("Proportional residual error on CSF Abeta42 (fraction; not reported)") # not reported
    propSd_abeta40 <- fixed(0); label("Proportional residual error on CSF Abeta40 (fraction; not reported)") # not reported
    propSd_abeta38 <- fixed(0); label("Proportional residual error on CSF Abeta38 (fraction; not reported)") # not reported
    propSd_abeta37 <- fixed(0); label("Proportional residual error on CSF Abeta37 (fraction; not reported)") # not reported
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma concentration: mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc

    # Indirect-response parameters; baseline = kin / kout (Methods)
    rbase_abeta42 <- exp(lrbase_abeta42)
    kout_abeta42 <- exp(lkout_abeta42)
    kin_abeta42 <- kout_abeta42 * rbase_abeta42
    imax_abeta42 <- exp(limax_abeta42)
    ic50_abeta42 <- exp(lic50_abeta42)

    rbase_abeta40 <- exp(lrbase_abeta40)
    kout_abeta40 <- exp(lkout_abeta40)
    kin_abeta40 <- kout_abeta40 * rbase_abeta40
    imax_abeta40 <- exp(limax_abeta40)
    ic50_abeta40 <- exp(lic50_abeta40)

    rbase_abeta38 <- exp(lrbase_abeta38)
    kout_abeta38 <- exp(lkout_abeta38)
    kin_abeta38 <- kout_abeta38 * rbase_abeta38
    emax_abeta38 <- exp(lemax_abeta38)
    ec50_abeta38 <- exp(lec50_abeta38)
    hill_abeta38 <- exp(lhill_abeta38)

    rbase_abeta37 <- exp(lrbase_abeta37)
    kout_abeta37 <- exp(lkout_abeta37)
    kin_abeta37 <- kout_abeta37 * rbase_abeta37
    emax_abeta37 <- exp(lemax_abeta37)
    ec50_abeta37 <- exp(lec50_abeta37)
    hill_abeta37 <- exp(lhill_abeta37)

    # Eq. (2): hyperbolic inhibition of production (Abeta42, Abeta40)
    inh_abeta42 <- imax_abeta42 * Cc / (ic50_abeta42 + Cc)
    inh_abeta40 <- imax_abeta40 * Cc / (ic50_abeta40 + Cc)
    # Eq. (3): sigmoid stimulation of production (Abeta38, Abeta37)
    stim_abeta38 <- emax_abeta38 * Cc^hill_abeta38 / (ec50_abeta38^hill_abeta38 + Cc^hill_abeta38)
    stim_abeta37 <- emax_abeta37 * Cc^hill_abeta37 / (ec50_abeta37^hill_abeta37 + Cc^hill_abeta37)

    d/dt(abeta42_csf) <- kin_abeta42 * (1 - inh_abeta42) - kout_abeta42 * abeta42_csf
    d/dt(abeta40_csf) <- kin_abeta40 * (1 - inh_abeta40) - kout_abeta40 * abeta40_csf
    d/dt(abeta38_csf) <- kin_abeta38 * (1 + stim_abeta38) - kout_abeta38 * abeta38_csf
    d/dt(abeta37_csf) <- kin_abeta37 * (1 + stim_abeta37) - kout_abeta37 * abeta37_csf
    abeta42_csf(0) <- rbase_abeta42
    abeta40_csf(0) <- rbase_abeta40
    abeta38_csf(0) <- rbase_abeta38
    abeta37_csf(0) <- rbase_abeta37

    # Eq. (1): empirical placebo drift, serial-CSF study only; time is
    # hours after the single dose (the 2-h predose sample sees no drift).
    tpos <- max(t, 0)
    plac_abeta42 <- STUDY_B7991003 * (pmax_abeta42 + etapmax_abeta42) * tpos^pow_abeta42
    plac_abeta40 <- STUDY_B7991003 * (pmax_abeta40 + etapmax_abeta40) * tpos^pow_abeta40
    plac_abeta38 <- STUDY_B7991003 * (pmax_abeta38 + etapmax_abeta38) * tpos^pow_abeta38
    plac_abeta37 <- STUDY_B7991003 * (pmax_abeta37 + etapmax_abeta37) * tpos^pow_abeta37

    # Observed CSF concentrations = baseline + treatment (state) + placebo
    abeta42 <- abeta42_csf + plac_abeta42
    abeta40 <- abeta40_csf + plac_abeta40
    abeta38 <- abeta38_csf + plac_abeta38
    abeta37 <- abeta37_csf + plac_abeta37

    Cc ~ add(addSd) + prop(propSd)
    abeta42 ~ prop(propSd_abeta42)
    abeta40 ~ prop(propSd_abeta40)
    abeta38 ~ prop(propSd_abeta38)
    abeta37 ~ prop(propSd_abeta37)
  })
}
