Leung_2023_tenofovir <- function() {
  description <- "Semi-physiological 16-compartment population PK model for tenofovir (TFV) after a single oral dose of tenofovir disoproxil fumarate (TDF) in healthy premenopausal women (Leung 2023). The oral depot feeds a two-compartment plasma model by first-order absorption (ka) and, in parallel, a seven-stage gut-transit chain (ktr) whose last stage delivers drug to rectal tissue (ka_rectum) or out of the body (kel_gut). The total plasma clearance CLtt is routed as fixed fractions into vaginal, cervical and rectal tissue (fcl_<tissue>); each tissue holds parent TFV, cleared at cl_<tissue>, of which the fraction fm_<tissue> forms the active intracellular metabolite tenofovir diphosphate (TFV-DP) in a separate tissue state cleared at cl_<tissue>_tfvdp. Tissue volumes are fixed to physiological tissue weights (density 1 g/mL). Amounts in nmol, concentrations in nmol/L; dose the TFV molar equivalent of the TDF dose (mg TDF / 635.52 x 1e6 nmol) into depot. Seven outputs: plasma TFV and TFV / TFV-DP in each of the three tissues."
  reference <- "Leung E, Cottrell ML, Sykes C, White N, Kashuba ADM, Dumond JB. A multicompartment population PK model to predict tenofovir and emtricitabine mucosal tissue concentrations for HIV prevention. CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1922. doi:10.1002/psp4.13042"
  vignette <- "Leung_2023_tenofovir_emtricitabine"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tenofovir", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tenofovir", units = "nmol", specimen = "plasma", verified = TRUE),
    vagina = list(analyte = "tenofovir", units = "nmol", specimen = "tissue", verified = TRUE),
    cervix = list(analyte = "tenofovir", units = "nmol", specimen = "tissue", verified = TRUE),
    rectum = list(analyte = "tenofovir", units = "nmol", specimen = "tissue", verified = TRUE),
    vagina_tfvdp = list(analyte = "tenofovir diphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    cervix_tfvdp = list(analyte = "tenofovir diphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    rectum_tfvdp = list(analyte = "tenofovir diphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    transit1 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "tenofovir", units = "nmol", specimen = "administration site", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "IQR 22.75-31.25 years (TDF arm, N = 25 enrolled)",
    age_median = "27 years",
    weight_range = "IQR 60.4-76.8 kg",
    weight_median = "67.4 kg",
    bmi_median = "24.1 kg/m^2 (IQR 21.6-26.9)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 64, Black = 32, Asian = 4),
    disease_state = "Healthy premenopausal HIV-negative women (phase I dose-ranging tissue PK study, NCT01330199)",
    dose_range = "Single oral dose of TDF 150, 300 or 600 mg (7, 8 and 8 modelled participants)",
    regions = "USA (University of North Carolina)",
    notes = "Demographics from Leung 2023 Table S1 (25 enrolled in the TDF arm); 23 participants with 413 observations were modelled (275 plasma TFV, 23 TFV and 23 TFV-DP in each of vaginal, cervical and rectal tissue; Results). Each participant gave one vaginal, cervical and rectal biopsy at 6, 12, 24 or 48 h post-dose. No covariates were tested because of the homogeneous population."
  )

  ini({
    # Plasma two-compartment model. Leung 2023 Table 1 'Final estimate TFV'
    # column; the supplement 2 TFV control stream (run17) fixes the structure.
    lka <- log(0.863); label("First-order absorption rate from depot into plasma, Ka (1/h)") # Table 1 Ka = 0.863
    lvc <- log(331); label("Apparent central volume Vc (L)") # Table 1 Vc = 331
    lvp <- log(843); label("Apparent peripheral volume Vp (L)") # Table 1 Vp = 843
    lq <- log(142); label("Apparent intercompartmental clearance Q (L/h)") # Table 1 Q = 142
    lcl <- log(58.7); label("Total apparent plasma clearance CLtt, including the routed tissue fractions (L/h)") # Table 1 CLtt = 58.7

    # Fractions of CLtt routed into each tissue (control stream THETA(6)-(8)).
    lfcl_vagina <- log(7.90e-5); label("Fraction of CLtt routed into vaginal tissue, Fv (fraction)") # Table 1 Fv = 7.90 x 1e-5
    lfcl_cervix <- log(1.50e-5); label("Fraction of CLtt routed into cervical tissue, Fe (fraction)") # Table 1 Fe = 1.50 x 1e-5
    lfcl_rectum <- fixed(log(7.00e-5)); label("Fraction of CLtt routed into rectal tissue, Fr (fraction)") # Table 1 Fr = 7.00 x 1e-5; control stream THETA(8) '(0.00007 FIX)'

    # Vaginal tissue (parent TFV and TFV-DP)
    fm_vagina <- 0.243; label("Fraction of vaginal-tissue TFV clearance forming TFV-DP, Fvt (fraction)") # Table 1 Fvt = 0.243
    lcl_vagina <- log(0.0111); label("Total clearance of TFV out of vaginal tissue, CLttvvtp (L/h)") # Table 1 CLttvvtp = 0.0111
    lcl_vagina_tfvdp <- log(0.0410); label("Clearance of TFV-DP from vaginal tissue, CLvvtp (L/h)") # Table 1 CLvvtp = 0.0410
    lv_vagina <- fixed(log(0.09)); label("Vaginal tissue volume for parent TFV, Vv (L)") # Results: Vv fixed to 0.09 L; control stream THETA(11) '(0.09 FIX)'
    lv_vagina_tfvdp <- fixed(log(0.09)); label("Vaginal tissue volume for TFV-DP, Vvtp (L)") # Results: Vvtp fixed to 0.09 L; control stream THETA(12) '(0.09 FIX)'

    # Cervical tissue (parent TFV and TFV-DP)
    fm_cervix <- 0.0292; label("Fraction of cervical-tissue TFV clearance forming TFV-DP, Fet (fraction)") # Table 1 Fet = 0.0292
    lcl_cervix <- log(1.83e-3); label("Total clearance of TFV out of cervical tissue, CLttve (L/h)") # Table 1 CLttve = 1.83 x 1e-3
    lcl_cervix_tfvdp <- log(2.07e-3); label("Clearance of TFV-DP from cervical tissue, CLvetp (L/h)") # Table 1 CLvetp = 2.07 x 1e-3
    lv_cervix <- fixed(log(0.01)); label("Cervical tissue volume for parent TFV, Ve (L)") # Results: Ve fixed to 0.01 L; control stream THETA(16) '(0.01 FIX)'
    lv_cervix_tfvdp <- fixed(log(0.01)); label("Cervical tissue volume for TFV-DP, Vetp (L)") # Results: Vetp fixed to 0.01 L; control stream THETA(17) '(0.01 FIX)'

    # Rectal tissue (parent TFV and TFV-DP)
    fm_rectum <- fixed(1); label("Fraction of rectal-tissue TFV clearance forming TFV-DP, Frt (fraction)") # Table 1 Frt = 1.00 (footnote a, fixed); control stream THETA(19) '(1 FIX)'
    lcl_rectum <- log(4.77e-3); label("Total clearance of TFV out of rectal tissue, CLttvr (L/h)") # Table 1 CLttvr = 4.77 x 1e-3
    lcl_rectum_tfvdp <- log(0.140); label("Clearance of TFV-DP from rectal tissue, CLvrtp (L/h)") # Table 1 CLvrtp = 0.140
    lv_rectum <- fixed(log(0.17)); label("Rectal tissue volume for parent TFV, Vr (L)") # Results: Vr fixed to 0.17 L; control stream THETA(21) '(0.17 FIX)'
    lv_rectum_tfvdp <- fixed(log(0.17)); label("Rectal tissue volume for TFV-DP, Vrtp (L)") # Results: Vrtp fixed to 0.17 L; control stream THETA(22) '(0.17 FIX)'

    # Gut-transit chain to rectal tissue
    lktr <- log(0.0752); label("Gut-transit rate constant from depot and between the seven transit stages, Kg (1/h)") # Table 1 Kg = 0.0752
    lka_rectum <- log(0.0589); label("Transfer rate from the last gut-transit stage into rectal tissue, Kga (1/h)") # Table 1 Kga = 0.0589
    lkel_gut <- fixed(log(1)); label("Loss rate from the last gut-transit stage out of the body, Kgr (1/h)") # Table 1 Kgr = 1.00 (footnote a, fixed); control stream THETA(26) '(1 FIX)'

    # IIV: omega^2 = log(1 + CV^2) from the Table 1 %CV column. Vc and CLtt
    # carry the paper's 100% correlation (Table 1 off-diagonal 'Vc 100';
    # control stream $OMEGA BLOCK(2) on ETA(2)-ETA(3)).
    etalka ~ 0.132894 # Table 1 omega on Ka = 37.7 %CV
    etalvc + etalcl ~ c(0.141586, 0.0843447, 0.0502453) # Table 1 omega on Vc = 39.0 %CV, on CLtt = 22.7 %CV, correlation 100%
    etalcl_vagina_tfvdp ~ 0.66915 # Table 1 omega on CLvvtp = 97.6 %CV
    etalcl_cervix_tfvdp ~ 1.37584 # Table 1 omega on CLvetp = 172 %CV
    etalktr ~ 0.571479 # Table 1 omega on Kg = 87.8 %CV

    # Proportional residual error; Table 1 sigma rows in control-stream
    # EPS(1)-EPS(7) order. The last row, printed 'Metabolite in plasma', is
    # EPS(7) on the rectal TFV-DP observations (CMT 9).
    propSd <- 0.284; label("Proportional residual error, plasma TFV (fraction)") # Table 1 sigma parent drug in plasma = 28.4 %CV
    propSd_Cvagina <- 0.391; label("Proportional residual error, vaginal tissue TFV (fraction)") # Table 1 sigma parent drug in vaginal compartment = 39.1 %CV
    propSd_Cvagina_tfvdp <- 0.423; label("Proportional residual error, vaginal tissue TFV-DP (fraction)") # Table 1 sigma metabolite in vaginal compartment = 42.3 %CV
    propSd_Ccervix <- 0.420; label("Proportional residual error, cervical tissue TFV (fraction)") # Table 1 sigma parent drug in cervical compartment = 42.0 %CV
    propSd_Ccervix_tfvdp <- 0.311; label("Proportional residual error, cervical tissue TFV-DP (fraction)") # Table 1 sigma metabolite in cervical compartment = 31.1 %CV
    propSd_Crectum <- 0.513; label("Proportional residual error, rectal tissue TFV (fraction)") # Table 1 sigma parent drug in rectal compartment = 51.3 %CV
    propSd_Crectum_tfvdp <- 0.565; label("Proportional residual error, rectal tissue TFV-DP (fraction)") # Table 1 sigma 'metabolite in plasma' row = 56.5 %CV (EPS(7), rectal metabolite)
  })

  model({
    # Individual parameters (exponential IIV, Leung 2023 Results)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    cl <- exp(lcl + etalcl)

    fcl_vagina <- exp(lfcl_vagina)
    fcl_cervix <- exp(lfcl_cervix)
    fcl_rectum <- exp(lfcl_rectum)

    cl_vagina <- exp(lcl_vagina)
    cl_vagina_tfvdp <- exp(lcl_vagina_tfvdp + etalcl_vagina_tfvdp)
    v_vagina <- exp(lv_vagina)
    v_vagina_tfvdp <- exp(lv_vagina_tfvdp)

    cl_cervix <- exp(lcl_cervix)
    cl_cervix_tfvdp <- exp(lcl_cervix_tfvdp + etalcl_cervix_tfvdp)
    v_cervix <- exp(lv_cervix)
    v_cervix_tfvdp <- exp(lv_cervix_tfvdp)

    cl_rectum <- exp(lcl_rectum)
    cl_rectum_tfvdp <- exp(lcl_rectum_tfvdp)
    v_rectum <- exp(lv_rectum)
    v_rectum_tfvdp <- exp(lv_rectum_tfvdp)

    ktr <- exp(lktr + etalktr)
    ka_rectum <- exp(lka_rectum)
    kel_gut <- exp(lkel_gut)

    # Micro-constants. CLtt leaves plasma in full: (1 - Fv - Fe - Fr) * CLtt
    # out of the body and Fv / Fe / Fr * CLtt into the three tissues
    # (Equation 1; control stream K2T0, K2T4, K2T5, K2T6).
    kel <- (1 - fcl_vagina - fcl_cervix - fcl_rectum) * cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_vagina <- fcl_vagina * cl / vc
    k_cervix <- fcl_cervix * cl / vc
    k_rectum <- fcl_rectum * cl / vc

    # Depot empties into plasma (Ka) and into the gut-transit chain (Kg) in
    # parallel (control stream K1T2 = Ka, K1T10 = Kg).
    d/dt(depot) <- -(ka + ktr) * depot
    d/dt(central) <- ka * depot - (kel + k12 + k_vagina + k_cervix + k_rectum) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Tissue parent loses CLtt<tissue> / V<tissue> in total, of which the
    # fraction fm_<tissue> forms the metabolite (Equations 2-3).
    d/dt(vagina) <- k_vagina * central - cl_vagina / v_vagina * vagina
    d/dt(cervix) <- k_cervix * central - cl_cervix / v_cervix * cervix
    d/dt(rectum) <- k_rectum * central - cl_rectum / v_rectum * rectum + ka_rectum * transit7
    d/dt(vagina_tfvdp) <- fm_vagina * cl_vagina / v_vagina * vagina - cl_vagina_tfvdp / v_vagina_tfvdp * vagina_tfvdp
    d/dt(cervix_tfvdp) <- fm_cervix * cl_cervix / v_cervix * cervix - cl_cervix_tfvdp / v_cervix_tfvdp * cervix_tfvdp
    d/dt(rectum_tfvdp) <- fm_rectum * cl_rectum / v_rectum * rectum - cl_rectum_tfvdp / v_rectum_tfvdp * rectum_tfvdp

    # Seven-stage gut transit; the last stage drains to rectal tissue (Kga)
    # and out of the body (Kgr).
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - (ka_rectum + kel_gut) * transit7

    # Concentrations (nmol/L; tissue density 1 g/mL)
    Cc <- central / vc
    Cvagina <- vagina / v_vagina
    Cvagina_tfvdp <- vagina_tfvdp / v_vagina_tfvdp
    Ccervix <- cervix / v_cervix
    Ccervix_tfvdp <- cervix_tfvdp / v_cervix_tfvdp
    Crectum <- rectum / v_rectum
    Crectum_tfvdp <- rectum_tfvdp / v_rectum_tfvdp

    Cc ~ prop(propSd)
    Cvagina ~ prop(propSd_Cvagina)
    Cvagina_tfvdp ~ prop(propSd_Cvagina_tfvdp)
    Ccervix ~ prop(propSd_Ccervix)
    Ccervix_tfvdp ~ prop(propSd_Ccervix_tfvdp)
    Crectum ~ prop(propSd_Crectum)
    Crectum_tfvdp ~ prop(propSd_Crectum_tfvdp)
  })
}
