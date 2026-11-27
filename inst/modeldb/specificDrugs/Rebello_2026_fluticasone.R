Rebello_2026_fluticasone <- function() {
  description <- paste(
    "Three-compartment population PK model for inhaled fluticasone propionate",
    "(1000 ug by metered-dose inhaler, two formulations) in 96 healthy adult",
    "men (Rebello 2026). Three parallel lagged first-order absorption depots",
    "(central lung, peripheral lung, and a delayed component the authors",
    "attribute to the gut) feed the central compartment, each with its own",
    "dose fraction, absorption rate and cumulative lag, estimated separately",
    "for each formulation; formulation 2 carries a relative bioavailability.",
    "Body weight scales the apparent volumes (exponent 1) and clearances",
    "(exponent 0.75) around 62 kg. CAUTION: the published final estimates do",
    "not reproduce the paper's own NCA or VPC -- at 1000 ug the typical",
    "subject peaks near 6 h at about 240 pg/mL with AUC0-36 about 5900",
    "pg*h/mL, against an observed Tmax of about 1.5-1.8 h, Cmax about 140",
    "pg/mL and AUC0-t about 1430 pg*h/mL; the excess comes from the delayed",
    "component. Use for relative formulation and absorption-route comparisons,",
    "not for absolute exposure."
  )
  reference <- paste(
    "Rebello J, Shukla S, Naidu KR, Patil J, Morde N, Ducharme MP.",
    "Understanding Pulmonary and Systemic Pharmacokinetics of Inhaled",
    "Fluticasone Propionate in Healthy Subjects Using Population PK Analysis.",
    "Eur J Drug Metab Pharmacokinet. 2026.",
    "doi:10.1007/s13318-026-00989-0"
  )
  vignette <- "Rebello_2026_fluticasone"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of all apparent volumes (exponent 1) and",
        "clearances (exponent 0.75) around 62 kg, the median subject weight",
        "(Table 5 footnote; supplement DIFFEQ block, (R(5)/62)**1 and",
        "(R(5)/62)**0.75). The paper calls the covariate 'average body",
        "weight (ABW)' in Table 3 and Fig. 1 and 'total body weight' in the",
        "covariate-testing text; it is the subject's measured weight."
      ),
      source_name = "ABW / R(5)"
    ),
    FORM_FP_F2 = list(
      description = paste(
        "1 = the dose is fluticasone propionate HFA MDI formulation 2,",
        "0 = formulation 1"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (formulation 1, the reference product; F = 1)",
      notes = paste(
        "Per-dose-record indicator in the two-way crossover: each subject",
        "carries 0 on one period's dose and 1 on the other's. It selects the",
        "formulation-specific absorption parameters (Table 4) and applies",
        "Frel to every formulation-2 input. In the supplement code",
        "formulation 1 is the 'R' input R(1) and formulation 2 is the 'T'",
        "input R(2), multiplied by Frel."
      ),
      source_name = "Formulation 1 / Formulation 2 (supplement R(1) / R(2), R(3))"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluticasone propionate", units = "ug", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "fluticasone propionate", units = "ug", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "fluticasone propionate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 96L,
    n_studies = 4L,
    age_range = "18-45 years (inclusion criterion)",
    age_median = "27 years (mean; SD 5.34)",
    weight_range = "> 45 kg (inclusion criterion)",
    weight_median = "62 kg (reference weight); mean 62.45 kg (SD 7.42)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state = "healthy non-smoking adult men",
    dose_range = paste(
      "Single 1000 ug dose of fluticasone propionate HFA by metered-dose",
      "inhaler (four consecutive 250 ug actuations within 2 min, no spacer),",
      "for each of two formulations in a two-way crossover with a 7-14 day",
      "washout"
    ),
    regions = "India (single centre, Mumbai)",
    notes = paste(
      "Rebello 2026 Table 1 (four randomised crossover studies of 24, 22, 22",
      "and 28 subjects, 20 PK samples per period) and Table 2 (demographics:",
      "height 1.68 m, BMI 21.90 kg/m2, FEV1 91.89% predicted). 3840",
      "concentrations; LLOQ 9.9 pg/mL. The two formulations had comparable",
      "fine-particle mass and differed mainly in the source of the active",
      "ingredient."
    )
  )

  ini({
    # All parameters were ln-normally distributed in the ADAPT 5 MLEM fit
    # (Methods). Typical values are the 'Population predicted' geometric
    # means of Table 4 (absorption) and Table 5 (disposition, for a 62 kg
    # subject); IIV variances are log(GeoCV^2 + 1) of the matching
    # 'Population predicted' GeoCV column. No covariances are reported, so
    # the IIV is diagonal.
    #
    # Suffixes _f1 / _f2 mark the formulation-1 / formulation-2 estimates of
    # the absorption block, which the authors estimated separately for each
    # formulation (supplement parameters ...R... / ...T...).

    # ---- Disposition (shared by both formulations; Table 5) ----
    lcl <- log(145); label("Apparent clearance CL/F for a 62 kg subject (L/h)") # Table 5 'CL/F' population predicted 145.0 L/h
    lvc <- log(111); label("Apparent central volume Vc/F for a 62 kg subject (L)") # Table 5 'Vc/F' population predicted 111.0 L
    lq <- log(248); label("Apparent intercompartmental clearance CLD/F to peripheral1 for a 62 kg subject (L/h)") # Table 5 'CLD/F' population predicted 248.0 L/h
    lvp <- log(80.3); label("Apparent peripheral1 volume Vp/F for a 62 kg subject (L)") # Table 5 'Vp/F' population predicted 80.3 L
    lq2 <- log(2030); label("Apparent intercompartmental clearance CLD2/F to peripheral2 for a 62 kg subject (L/h)") # Table 5 'CLD2/F' population predicted 2030.0 L/h
    lvp2 <- log(1070); label("Apparent peripheral2 volume Vp2/F for a 62 kg subject (L)") # Table 5 'Vp2/F' population predicted 1070.0 L

    # ---- Allometric exponents (fixed; Results and supplement DIFFEQ) ----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on (WT/62) for CL/F, CLD/F and CLD2/F (unitless)") # Results 'expected allometric exponents of 1 and 0.75'; supplement (R(5)/62)**0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on (WT/62) for Vc/F, Vp/F and Vp2/F (unitless)") # Results 'expected allometric exponents of 1 and 0.75'; supplement (R(5)/62)**1

    # ---- Absorption, formulation 1 (Table 4, supplement 'R' parameters) ----
    lka_f1 <- log(1.15); label("Central-lung absorption rate constant kaINHD1c, formulation 1 (1/h)") # Table 4 'kaINHD1c' population predicted 1.15 1/h
    lka2_f1 <- log(0.77); label("Peripheral-lung absorption rate constant kaINHD1p, formulation 1 (1/h)") # Table 4 'kaINHD1p' population predicted 0.77 1/h
    lka3_f1 <- log(0.07); label("Delayed (gut) absorption rate constant kaINHGD1, formulation 1 (1/h)") # Table 4 'kaINHGD1' population predicted 0.07 1/h
    lfdepot_f1 <- log(0.134); label("Fraction of the dose absorbed from the central lung, formulation 1 (fraction)") # Table 4 formulation 1 'central lung absorption' population predicted 13.4 %
    lfdepot2_f1 <- log(0.119); label("Fraction of the dose absorbed from the peripheral lung, formulation 1 (fraction)") # Table 4 formulation 1 'peripheral lung absorption' population predicted 11.9 %
    ltlag_f1 <- log(0.04); label("Central-lung absorption lag TlagD1c, formulation 1 (h)") # Table 4 'TlagD1c' population predicted 0.04 h
    ltlag2_f1 <- log(0.17); label("Peripheral-lung lag increment over the central-lung lag TlagD1p, formulation 1 (h)") # Table 4 'TlagD1p' population predicted 0.17 h; cumulative per supplement Tlg2 = P(13) + P(14)
    ltlag3_f1 <- log(0.60); label("Delayed (gut) lag increment over the peripheral-lung lag TlagGD1, formulation 1 (h)") # Table 4 'TlagGD1' population predicted 0.60 h; cumulative per supplement Tlg3 = P(13) + P(14) + P(20)

    # ---- Absorption, formulation 2 (Table 4, supplement 'T' parameters) ----
    lka_f2 <- log(0.76); label("Central-lung absorption rate constant kaINHD2c, formulation 2 (1/h)") # Table 4 'kaINHD2c' population predicted 0.76 1/h
    lka2_f2 <- log(0.59); label("Peripheral-lung absorption rate constant kaINHD2p, formulation 2 (1/h)") # Table 4 'kaINHD2p' population predicted 0.59 1/h
    lka3_f2 <- log(0.07); label("Delayed (gut) absorption rate constant kaINHGD2, formulation 2 (1/h)") # Table 4 'kaINHGD2' population predicted 0.07 1/h
    lfdepot_f2 <- log(0.200); label("Fraction of the dose absorbed from the central lung, formulation 2 (fraction)") # Table 4 formulation 2 'central lung absorption' population predicted 20.0 %
    lfdepot2_f2 <- log(0.147); label("Fraction of the dose absorbed from the peripheral lung, formulation 2 (fraction)") # Table 4 formulation 2 'peripheral lung absorption' population predicted 14.7 %
    ltlag_f2 <- log(0.04); label("Central-lung absorption lag TlagD2c, formulation 2 (h)") # Table 4 'TlagD2c' population predicted 0.04 h
    ltlag2_f2 <- log(0.14); label("Peripheral-lung lag increment over the central-lung lag TlagD2p, formulation 2 (h)") # Table 4 'TlagD2p' population predicted 0.14 h; cumulative per supplement Tlg2 = P(1) + P(12)
    ltlag3_f2 <- log(0.37); label("Delayed (gut) lag increment over the peripheral-lung lag TlagGD2, formulation 2 (h)") # Table 4 'TlagGD2' population predicted 0.37 h; cumulative per supplement Tlg3 = P(1) + P(12) + P(21)
    lfrel <- log(1.02); label("Bioavailability of formulation 2 relative to formulation 1 Frel (unitless)") # Table 4 'Frel' population predicted 1.02

    # ---- IIV: log(GeoCV^2 + 1) of the Table 4 / Table 5 population GeoCV ----
    etalcl ~ 0.070365 # Table 5 'CL/F' GeoCV 27.0 %
    etalvc ~ 0.401597 # Table 5 'Vc/F' GeoCV 70.3 %
    etalq ~ 0.309251 # Table 5 'CLD/F' GeoCV 60.2 %
    etalvp ~ 0.537859 # Table 5 'Vp/F' GeoCV 84.4 %
    etalq2 ~ 0.324371 # Table 5 'CLD2/F' GeoCV 61.9 %
    etalvp2 ~ 0.193463 # Table 5 'Vp2/F' GeoCV 46.2 %
    etalka_f1 ~ 0.198051 # Table 4 'kaINHD1c' GeoCV 46.8 %
    etalka2_f1 ~ 0.232812 # Table 4 'kaINHD1p' GeoCV 51.2 %
    etalka3_f1 ~ 0.009174 # Table 4 'kaINHGD1' GeoCV 9.6 %
    etalfdepot_f1 ~ 0.267670 # Table 4 formulation 1 central lung fraction GeoCV 55.4 %
    etalfdepot2_f1 ~ 0.288245 # Table 4 formulation 1 peripheral lung fraction GeoCV 57.8 %
    etaltlag_f1 ~ 0.095179 # Table 4 'TlagD1c' GeoCV 31.6 %
    etaltlag2_f1 ~ 0.040773 # Table 4 'TlagD1p' GeoCV 20.4 %
    etaltlag3_f1 ~ 0.203451 # Table 4 'TlagGD1' GeoCV 47.5 %
    etalka_f2 ~ 0.182162 # Table 4 'kaINHD2c' GeoCV 44.7 %
    etalka2_f2 ~ 0.279618 # Table 4 'kaINHD2p' GeoCV 56.8 %
    etalka3_f2 ~ 0.031538 # Table 4 'kaINHGD2' GeoCV 17.9 %
    etalfdepot_f2 ~ 0.112458 # Table 4 formulation 2 central lung fraction GeoCV 34.5 %
    etalfdepot2_f2 ~ 0.098071 # Table 4 formulation 2 peripheral lung fraction GeoCV 32.1 %
    etaltlag_f2 ~ 0.056913 # Table 4 'TlagD2c' GeoCV 24.2 %
    etaltlag2_f2 ~ 0.228767 # Table 4 'TlagD2p' GeoCV 50.7 %
    etaltlag3_f2 ~ 0.337868 # Table 4 'TlagGD2' GeoCV 63.4 %
    etalfrel ~ 0.110614 # Table 4 'Frel' GeoCV 34.2 %

    # ---- Residual error (supplement VARMOD: SD = ADD + PROP * Y) ----
    propSd <- 0.0897; label("Proportional residual SD (fraction)") # Results 'residual variability in fitted concentrations ... was 8.97%'
    addSd <- fixed(0); label("Additive residual SD (pg/mL; 0 -- not reported in the source)") # supplement VARMOD has an additive term PV(1); its estimate is not published
  })
  model({
    # Formulation selector: 1 = formulation 2, 0 = formulation 1.
    f2 <- FORM_FP_F2
    f1 <- 1 - FORM_FP_F2

    # Disposition with allometric scaling around 62 kg (supplement DIFFEQ).
    cl <- exp(lcl + etalcl) * (WT / 62)^e_wt_cl_q
    q <- exp(lq + etalq) * (WT / 62)^e_wt_cl_q
    q2 <- exp(lq2 + etalq2) * (WT / 62)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 62)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 62)^e_wt_vc_vp
    vp2 <- exp(lvp2 + etalvp2) * (WT / 62)^e_wt_vc_vp

    # Formulation-specific absorption parameters. Each subject carries one
    # set per formulation (the crossover gives every subject both).
    ka_f1 <- exp(lka_f1 + etalka_f1)
    ka2_f1 <- exp(lka2_f1 + etalka2_f1)
    ka3_f1 <- exp(lka3_f1 + etalka3_f1)
    fc_f1 <- exp(lfdepot_f1 + etalfdepot_f1)
    fp_f1 <- exp(lfdepot2_f1 + etalfdepot2_f1)
    tlag_f1 <- exp(ltlag_f1 + etaltlag_f1)
    tlag2_f1 <- exp(ltlag2_f1 + etaltlag2_f1)
    tlag3_f1 <- exp(ltlag3_f1 + etaltlag3_f1)

    ka_f2 <- exp(lka_f2 + etalka_f2)
    ka2_f2 <- exp(lka2_f2 + etalka2_f2)
    ka3_f2 <- exp(lka3_f2 + etalka3_f2)
    fc_f2 <- exp(lfdepot_f2 + etalfdepot_f2)
    fp_f2 <- exp(lfdepot2_f2 + etalfdepot2_f2)
    tlag_f2 <- exp(ltlag_f2 + etaltlag_f2)
    tlag2_f2 <- exp(ltlag2_f2 + etaltlag2_f2)
    tlag3_f2 <- exp(ltlag3_f2 + etaltlag3_f2)
    frel <- exp(lfrel + etalfrel)

    ka <- ka_f1 * f1 + ka_f2 * f2
    ka2 <- ka2_f1 * f1 + ka2_f2 * f2
    ka3 <- ka3_f1 * f1 + ka3_f2 * f2
    fc <- fc_f1 * f1 + fc_f2 * f2
    fp_raw <- fp_f1 * f1 + fp_f2 * f2

    # Dose fractions (supplement DIFFEQ): central lung = alpha, peripheral
    # lung = beta, delayed = 1 - alpha - beta; when alpha + beta > 1 the
    # code sets peripheral = 1 - alpha and delayed = 0.
    fp <- fp_raw
    if (fc + fp_raw > 1) {
      fp <- 1 - fc
    }
    fg <- 1 - fc - fp

    # Lags are cumulative (supplement DIFFEQ: Tlg2 = Tlg1 + TlagP,
    # Tlg3 = Tlg2 + TlagG).
    tlag <- tlag_f1 * f1 + tlag_f2 * f2
    tlag2 <- tlag + tlag2_f1 * f1 + tlag2_f2 * f2
    tlag3 <- tlag2 + tlag3_f1 * f1 + tlag3_f2 * f2

    # Frel multiplies every formulation-2 input (supplement R(2)*P(9)).
    fform <- f1 + frel * f2

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(depot3) <- -ka3 * depot3
    d/dt(central) <- ka * depot + ka2 * depot2 + ka3 * depot3 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Each administration is three dose records of the full nominal dose,
    # one per depot; f() splits it.
    f(depot) <- fform * fc
    f(depot2) <- fform * fp
    f(depot3) <- fform * fg
    alag(depot) <- tlag
    alag(depot2) <- tlag2
    alag(depot3) <- tlag3

    # Dose in ug, volume in L: central / vc is ug/L = ng/mL; x 1000 -> pg/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
