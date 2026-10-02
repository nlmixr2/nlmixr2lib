Morse_2022_acetaminophen <- function() {
  description <- "Two-compartment population PK model for acetaminophen (paracetamol) in healthy adult volunteers given intravenous, tablet, oral-suspension and sachet formulations of an acetaminophen + ibuprofen combination under fasted and fed conditions (Morse 2022). First-order absorption with a lag time from a depot; absorption half-life and lag time carry formulation-specific factors when fasted and formulation-specific factors when fed, all relative to the fasted tablet. Clearances scale allometrically (exponent 3/4) with normal fat mass (Ffat = 0.816) and volumes linearly with total body weight (Ffat fixed to 1), standardised to a 70 kg, 1.76 m male; fat-free mass is predicted from weight, height and sex (Janmahasatian). Combined additive + proportional residual error with between-subject variability on the residual magnitude."
  reference <- "Morse JD, Stanescu I, Atkinson HC, Anderson BJ. Population Pharmacokinetic Modelling of Acetaminophen and Ibuprofen: the Influence of Body Composition, Formulation and Feeding in Healthy Adult Volunteers. Eur J Drug Metab Pharmacokinet. 2022;47(4):497-507. doi:10.1007/s13318-022-00766-9"
  vignette <- "Morse_2022_acetaminophen_ibuprofen"
  paper_specific_etas <- c("etaltabs", "etaRUV")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight (TBM / TBW in the source)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the Janmahasatian fat-free-mass prediction (Eq. 2), normal fat mass (Eqs. 3-4) for clearances, and directly (Ffat fixed to 1, so NFM = WT) for the volumes. Standard 70 kg.",
      source_name = "TBM / TBW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Converted to metres inside model() for the Janmahasatian fat-free-mass equation (Morse 2022 Eq. 2 uses HT in m).",
      source_name = "HT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the sex-specific WHSmax / WHS50 constants of the Janmahasatian fat-free-mass equation (men 42.92 / 30.93 kg/m^2, women 37.99 / 35.98 kg/m^2). The standard individual (FFM 56.1 kg) is a 70 kg, 1.76 m male.",
      source_name = "SEX"
    ),
    FED = list(
      description = "Fed-state indicator for the oral dose, 1 = fed, 0 = fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted; at least 10 h fast)",
      notes = "Study MX-14b (fed) versus MX-14a / MXIV-01 / MXIV-06 (fasted). When FED = 1 the formulation-specific fed factors REPLACE the fasted formulation factors -- every factor is relative to the fasted tablet (Morse 2022 Table 2 footnote). Ignored for intravenous doses (no absorption).",
      source_name = "fed / fasted study"
    ),
    FORM_APAPIBU_SUSP = list(
      description = "Acetaminophen + ibuprofen ready-to-use oral suspension (Maxigesic Oral Suspension, 160 mg acetaminophen + 48 mg ibuprofen per 5 mL) indicator, 1 = suspension, 0 = other formulation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (film-coated tablet when FORM_POWDER is also 0)",
      notes = "Mutually exclusive with FORM_POWDER; both 0 selects the film-coated tablet (Maxigesic 500/150 or Maxigesic 325/97.5, pooled as one tablet formulation).",
      source_name = "formulation"
    ),
    FORM_POWDER = list(
      description = "Acetaminophen + ibuprofen sachet (powder dissolved in 200 mL water; Maxigesic Sachet) indicator, 1 = sachet, 0 = other formulation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (film-coated tablet when FORM_APAPIBU_SUSP is also 0)",
      notes = "Mutually exclusive with FORM_APAPIBU_SUSP. Comparator is the film-coated tablet.",
      source_name = "formulation"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "acetaminophen", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 116L,
    n_studies = 4L,
    n_observations = 6095L,
    age_range = "18-49 years",
    age_median = "24 years",
    weight_range = "49-116 kg",
    weight_median = "70.8 kg",
    height_range = "156-199 cm",
    height_median = "173 cm",
    bmi_range = "18.6-31.4 kg/m^2",
    sex_female_pct = 19.8,
    disease_state = "healthy adult volunteers",
    dose_range = "Single doses: acetaminophen 500 or 1000 mg IV (15 min infusion); 975 or 1000 mg orally as film-coated tablets, oral suspension or sachet, fasted or fed. Given in fixed combination with ibuprofen in every oral arm and in some IV arms.",
    regions = "Jordan",
    notes = "Pooled from four phase I crossover studies (Morse 2022 Section 2.1): MXIV-01 (n = 30, IV and tablet, fasted), MXIV-06 (n = 30, IV and tablet, fasted), MX-14a (n = 28, tablet / suspension / sachet, fasted) and MX-14b (n = 28, the same treatments fed). Demographics from Morse 2022 Table 1 (93 male / 23 female)."
  )

  ini({
    # Structural parameters standardised to a 70 kg, 1.76 m male (NFM_STD).
    lcl <- log(24.0); label("Clearance CL for the standard individual (L/h/70 kg)") # Table 2 CL = 24.0
    lvc <- log(43.7); label("Central volume V1 for the standard individual (L/70 kg)") # Table 2 V1 = 43.7 (the abstract's 43.5 is the Q2 value)
    lq <- log(43.5); label("Intercompartmental clearance Q2 for the standard individual (L/h/70 kg)") # Table 2 Q2 = 43.5
    lvp <- log(29.7); label("Peripheral volume V2 for the standard individual (L/70 kg)") # Table 2 V2 = 29.7
    lfdepot <- log(0.859); label("Oral bioavailability FPARA (fraction)") # Table 2 FPARA = 0.859
    ltabs <- log(11.5); label("Absorption half-life for the fasted tablet (min)") # Table 2 T1/2ABS = 11.5 min
    ltlag <- log(5.30); label("Absorption lag time for the fasted tablet (min)") # Table 2 TLAG = 5.30 min

    # Body-composition (normal fat mass) factors, Eqs. 3-5
    ffat_cl <- 0.816; label("Fat-mass fraction Ffat in the normal fat mass for clearances (unitless)") # Table 2 FFATCL = 0.816
    ffat_v <- fixed(1); label("Fat-mass fraction Ffat in the normal fat mass for volumes (unitless; 1 = total body weight)") # Table 2 FFATV = 1 FIX
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Section 2.3, theory-based exponent fixed at 3/4
    e_wt_vc <- fixed(1); label("Allometric exponent on V1 and V2 (unitless)") # Section 2.3, theory-based exponent fixed at 1

    # Formulation factors on T1/2ABS and TLAG when fasted (relative to the fasted tablet)
    e_form_apapibu_susp_tabs <- 0.394; label("Factor on T1/2ABS for the suspension, fasted (unitless)") # Table 2 F_FAST_TABS (Suspension) = 0.394
    e_form_powder_tabs <- 0.462; label("Factor on T1/2ABS for the sachet, fasted (unitless)") # Table 2 F_FAST_TABS (Sachet) = 0.462
    e_form_apapibu_susp_tlag <- 0.743; label("Factor on TLAG for the suspension, fasted (unitless)") # Table 2 F_FAST_LAG (Suspension) = 0.743
    e_form_powder_tlag <- 0.845; label("Factor on TLAG for the sachet, fasted (unitless)") # Table 2 F_FST_LAG (Sachet) = 0.845

    # Formulation factors on T1/2ABS and TLAG when fed (relative to the FASTED tablet)
    e_fed_tabs_tablet <- 1.87; label("Factor on T1/2ABS for the tablet, fed (unitless)") # Table 2 F_FED_TABS (Tablet) = 1.87
    e_fed_tabs_susp <- 2.52; label("Factor on T1/2ABS for the suspension, fed (unitless)") # Table 2 F_FED_TABS (Suspension) = 2.52
    e_fed_tabs_powder <- 2.30; label("Factor on T1/2ABS for the sachet, fed (unitless)") # Table 2 F_FED_TABS (Sachet) = 2.30
    e_fed_tlag_tablet <- 4.63; label("Factor on TLAG for the tablet, fed (unitless)") # Table 2 F_FED_LAG (Tablet) = 4.63
    e_fed_tlag_susp <- 2.93; label("Factor on TLAG for the suspension, fed (unitless)") # Table 2 F_FED_LAG (Suspension) = 2.93
    e_fed_tlag_powder <- 2.10; label("Factor on TLAG for the sachet, fed (unitless)") # Table 2 F_FED_LAG (Sachet) = 2.10

    # IIV. Table 2 PPV% = sqrt(omega^2) * 100, so omega^2 = (PPV/100)^2.
    # Block order CL, V1, Q2, V2; covariances = r * sd_i * sd_j with the
    # correlations of Supplementary Table S2 (acetaminophen).
    # SDs 0.173 (CL), 0.537 (V1), 0.615 (Q2), 0.466 (V2); r(V1,CL) = -0.004,
    # r(Q2,CL) = 0.005, r(Q2,V1) = -0.912, r(V2,CL) = 0.137, r(V2,V1) = -0.928,
    # r(V2,Q2) = 0.985.
    etalcl + etalvc + etalq + etalvp ~ c(
      0.029929,
      -0.000371604, 0.288369,
      0.000531975, -0.301193, 0.378225,
      0.0110447, -0.232225, 0.282291, 0.217156
    )
    etalfdepot ~ 0.021025 # Table 2 FPARA PPV 14.5% -> 0.145^2
    etaltabs ~ 0.732736 # Table 2 T1/2ABS PPV 85.6% -> 0.856^2
    etaltlag ~ 0.9216 # Table 2 TLAG PPV 96.0% -> 0.960^2
    etaRUV ~ 1.1025 # Table 2 RUV PPV 105% -> 1.05^2 (Eq. 8 eta on the residual SD)

    # Residual error (Eq. 8)
    addSd <- 0.064; label("Additive residual SD (mg/L)") # Table 2 RUV ADD = 0.064 mg/L
    propSd <- 0.070; label("Proportional residual SD (fraction)") # Table 2 RUV PROP = 7.0%
  })

  model({
    # Fat-free mass, Janmahasatian (Eq. 2); HT converted from cm to m
    ht_m <- HT / 100
    whs_max <- 42.92 * (1 - SEXF) + 37.99 * SEXF
    whs_50 <- 30.93 * (1 - SEXF) + 35.98 * SEXF
    ffm <- whs_max * ht_m^2 * WT / (whs_50 * ht_m^2 + WT)

    # Normal fat mass (Eqs. 3-4) and its standard value for a 70 kg male
    # with FFM 56.1 kg (Section 2.3)
    nfm_cl <- ffm + ffat_cl * (WT - ffm)
    nfm_v <- ffm + ffat_v * (WT - ffm)
    nfm_std_cl <- 56.1 + ffat_cl * (70 - 56.1)
    nfm_std_v <- 56.1 + ffat_v * (70 - 56.1)
    fsize_cl <- (nfm_cl / nfm_std_cl)^e_wt_cl
    fsize_v <- (nfm_v / nfm_std_v)^e_wt_vc

    cl <- exp(lcl + etalcl) * fsize_cl
    q <- exp(lq + etalq) * fsize_cl
    vc <- exp(lvc + etalvc) * fsize_v
    vp <- exp(lvp + etalvp) * fsize_v

    # Absorption: fasted factors (tablet = 1) or, when fed, the fed factors,
    # all relative to the fasted tablet
    form_tab <- 1 - FORM_APAPIBU_SUSP - FORM_POWDER
    f_tabs_fast <- form_tab + e_form_apapibu_susp_tabs * FORM_APAPIBU_SUSP + e_form_powder_tabs * FORM_POWDER
    f_tlag_fast <- form_tab + e_form_apapibu_susp_tlag * FORM_APAPIBU_SUSP + e_form_powder_tlag * FORM_POWDER
    f_tabs_fed <- e_fed_tabs_tablet * form_tab + e_fed_tabs_susp * FORM_APAPIBU_SUSP + e_fed_tabs_powder * FORM_POWDER
    f_tlag_fed <- e_fed_tlag_tablet * form_tab + e_fed_tlag_susp * FORM_APAPIBU_SUSP + e_fed_tlag_powder * FORM_POWDER
    tabs <- exp(ltabs + etaltabs) * (f_tabs_fast * (1 - FED) + f_tabs_fed * FED)
    tlag <- exp(ltlag + etaltlag) * (f_tlag_fast * (1 - FED) + f_tlag_fed * FED)
    ka <- log(2) / (tabs / 60)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag / 60

    Cc <- central / vc

    # Eq. 8: SD = sqrt((Cc * propSd)^2 + addSd^2) * exp(etaRUV)
    ruv_scale <- exp(etaRUV)
    addSd_i <- addSd * ruv_scale
    propSd_i <- propSd * ruv_scale
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
