Leung_2023_emtricitabine <- function() {
  description <- "Semi-physiological 16-compartment population PK model for emtricitabine (FTC) after a single oral dose in healthy premenopausal women (Leung 2023). The oral depot feeds a two-compartment plasma model by first-order absorption (ka) and, in parallel, a seven-stage gut-transit chain (ktr) whose last stage delivers drug to rectal tissue (ka_rectum) or out of the body (kel_gut). The total plasma clearance CLtt is routed as fixed fractions into vaginal, cervical and rectal tissue (fcl_<tissue>); each tissue holds parent FTC, cleared at cl_<tissue>, of which the fraction fm_<tissue> forms the active intracellular metabolite emtricitabine triphosphate (FTC-TP) in a separate tissue state cleared at cl_<tissue>_ftctp. Tissue volumes are fixed to physiological tissue weights (density 1 g/mL). Amounts in nmol, concentrations in nmol/L; dose the FTC dose in nmol (mg / 247.24 x 1e6) into depot. Seven outputs: plasma FTC and FTC / FTC-TP in each of the three tissues."
  reference <- "Leung E, Cottrell ML, Sykes C, White N, Kashuba ADM, Dumond JB. A multicompartment population PK model to predict tenofovir and emtricitabine mucosal tissue concentrations for HIV prevention. CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1922. doi:10.1002/psp4.13042"
  vignette <- "Leung_2023_tenofovir_emtricitabine"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "emtricitabine", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "emtricitabine", units = "nmol", specimen = "plasma", verified = TRUE),
    vagina = list(analyte = "emtricitabine", units = "nmol", specimen = "tissue", verified = TRUE),
    cervix = list(analyte = "emtricitabine", units = "nmol", specimen = "tissue", verified = TRUE),
    rectum = list(analyte = "emtricitabine", units = "nmol", specimen = "tissue", verified = TRUE),
    vagina_ftctp = list(analyte = "emtricitabine triphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    cervix_ftctp = list(analyte = "emtricitabine triphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    rectum_ftctp = list(analyte = "emtricitabine triphosphate", units = "nmol", specimen = "tissue", verified = TRUE),
    transit1 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "emtricitabine", units = "nmol", specimen = "administration site", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "IQR 21-27 years (FTC arm, N = 24 enrolled)",
    age_median = "22 years",
    weight_range = "IQR 57.7-72.3 kg",
    weight_median = "62.8 kg",
    bmi_median = "22.5 kg/m^2 (IQR 20.8-26.5)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 75, Black = 16.7, Asian = 4.2, AmericanIndian = 4.2),
    disease_state = "Healthy premenopausal HIV-negative women (phase I dose-ranging tissue PK study, NCT01330199)",
    dose_range = "Single oral dose of FTC 100, 200 or 400 mg (8, 7 and 8 modelled participants)",
    regions = "USA (University of North Carolina)",
    notes = "Demographics from Leung 2023 Table S1 (24 enrolled in the FTC arm); 23 participants with 573 observations were modelled (276 plasma FTC; 22 FTC and 22 FTC-TP in vaginal tissue; 23 FTC and 23 FTC-TP in each of cervical and rectal tissue; Results). Each participant gave one vaginal, cervical and rectal biopsy at 6, 12, 24 or 48 h post-dose. No covariates were tested because of the homogeneous population."
  )

  ini({
    # Plasma two-compartment model. Leung 2023 Table 1 'Final estimate FTC'
    # column; the supplement 2 FTC control stream (run15) fixes the structure.
    lka <- log(0.649); label("First-order absorption rate from depot into plasma, Ka (1/h)") # Table 1 Ka = 0.649
    lvc <- log(72.3); label("Apparent central volume Vc (L)") # Table 1 Vc = 72.3
    lvp <- log(122); label("Apparent peripheral volume Vp (L)") # Table 1 Vp = 122
    lq <- log(6.06); label("Apparent intercompartmental clearance Q (L/h)") # Table 1 Q = 6.06
    lcl <- log(18.9); label("Total apparent plasma clearance CLtt, including the routed tissue fractions (L/h)") # Table 1 CLtt = 18.9

    # Fractions of CLtt routed into each tissue (control stream THETA(6)-(8)).
    lfcl_vagina <- log(131e-5); label("Fraction of CLtt routed into vaginal tissue, Fv (fraction)") # Table 1 Fv = 131 x 1e-5
    lfcl_cervix <- log(6.90e-5); label("Fraction of CLtt routed into cervical tissue, Fe (fraction)") # Table 1 Fe = 6.90 x 1e-5
    lfcl_rectum <- log(401e-5); label("Fraction of CLtt routed into rectal tissue, Fr (fraction)") # Table 1 Fr = 401 x 1e-5

    # Vaginal tissue (parent FTC and FTC-TP)
    fm_vagina <- 0.325; label("Fraction of vaginal-tissue FTC clearance forming FTC-TP, Fvt (fraction)") # Table 1 Fvt = 0.325
    lcl_vagina <- log(0.0101); label("Total clearance of FTC out of vaginal tissue, CLttvvtp (L/h)") # Table 1 CLttvvtp printed '0101'; read as 0.0101 (control stream initial 0.0123; see vignette)
    lcl_vagina_ftctp <- log(0.0399); label("Clearance of FTC-TP from vaginal tissue, CLvvtp (L/h)") # Table 1 CLvvtp = 0.0399
    lv_vagina <- fixed(log(0.09)); label("Vaginal tissue volume for parent FTC, Vv (L)") # Results: Vv fixed to 0.09 L; control stream THETA(11) '(0.09 FIX)'
    lv_vagina_ftctp <- fixed(log(0.09)); label("Vaginal tissue volume for FTC-TP, Vvtp (L)") # Results: Vvtp fixed to 0.09 L; control stream THETA(12) '(0.09 FIX)'

    # Cervical tissue (parent FTC and FTC-TP)
    fm_cervix <- fixed(1); label("Fraction of cervical-tissue FTC clearance forming FTC-TP, Fet (fraction)") # Table 1 Fet = 1.00 (footnote a, fixed); control stream THETA(14) '(1 FIX)'
    lcl_cervix <- log(0.947e-3); label("Total clearance of FTC out of cervical tissue, CLttve (L/h)") # Table 1 CLttve = 0.947 x 1e-3
    lcl_cervix_ftctp <- log(8.2e-3); label("Clearance of FTC-TP from cervical tissue, CLvetp (L/h)") # Table 1 CLvetp = 8.2 x 1e-3
    lv_cervix <- fixed(log(0.01)); label("Cervical tissue volume for parent FTC, Ve (L)") # Results: Ve fixed to 0.01 L; control stream THETA(16) '(0.01 FIX)'
    lv_cervix_ftctp <- fixed(log(0.01)); label("Cervical tissue volume for FTC-TP, Vetp (L)") # Results: Vetp fixed to 0.01 L; control stream THETA(17) '(0.01 FIX)'

    # Rectal tissue (parent FTC and FTC-TP)
    fm_rectum <- 0.0107; label("Fraction of rectal-tissue FTC clearance forming FTC-TP, Frt (fraction)") # Table 1 Frt = 0.0107
    lcl_rectum <- log(22.3e-3); label("Total clearance of FTC out of rectal tissue, CLttvr (L/h)") # Table 1 CLttvr = 22.3 x 1e-3
    lcl_rectum_ftctp <- log(0.647); label("Clearance of FTC-TP from rectal tissue, CLvrtp (L/h)") # Table 1 CLvrtp = 0.647
    lv_rectum <- fixed(log(0.17)); label("Rectal tissue volume for parent FTC, Vr (L)") # Results: Vr fixed to 0.17 L; control stream THETA(21) '(0.17 FIX)'
    lv_rectum_ftctp <- fixed(log(0.17)); label("Rectal tissue volume for FTC-TP, Vrtp (L)") # Results: Vrtp fixed to 0.17 L; control stream THETA(22) '(0.17 FIX)'

    # Gut-transit chain to rectal tissue
    lktr <- log(0.0724); label("Gut-transit rate constant from depot and between the seven transit stages, Kg (1/h)") # Table 1 Kg = 0.0724
    lka_rectum <- fixed(log(1)); label("Transfer rate from the last gut-transit stage into rectal tissue, Kga (1/h)") # Table 1 Kga = 1.00 (footnote a, fixed); control stream THETA(25) '(1 FIX)'
    lkel_gut <- fixed(log(1)); label("Loss rate from the last gut-transit stage out of the body, Kgr (1/h)") # Table 1 Kgr = 1.00 (footnote a, fixed); control stream THETA(26) '(1 FIX)'

    # IIV: omega^2 = log(1 + CV^2) from the Table 1 %CV column; control
    # stream ETA(1)-ETA(8), all diagonal.
    etalka ~ 0.23688 # Table 1 omega on Ka = 51.7 %CV
    etalvc ~ 0.0411659 # Table 1 omega on Vc = 20.5 %CV
    etalcl_vagina ~ 0.136879 # Table 1 omega on CLttvvtp = 38.3 %CV
    etalcl_vagina_ftctp ~ 0.147731 # Table 1 omega on CLvvtp = 39.9 %CV
    etalcl_cervix ~ 0.111842 # Table 1 omega on CLttve = 34.4 %CV
    etalcl_cervix_ftctp ~ 1.31351 # Table 1 omega on CLvetp = 164.9 %CV
    etalcl_rectum ~ 0.563551 # Table 1 omega on CLttvr = 87 %CV
    etalktr ~ 0.399716 # Table 1 omega on Kg = 70.1 %CV

    # Proportional residual error; Table 1 sigma rows in control-stream
    # EPS(1)-EPS(7) order. The last row, printed 'Metabolite in plasma', is
    # EPS(7) on the rectal FTC-TP observations (CMT 9).
    propSd <- 0.316; label("Proportional residual error, plasma FTC (fraction)") # Table 1 sigma parent drug in plasma = 31.6 %CV
    propSd_Cvagina <- 0.242; label("Proportional residual error, vaginal tissue FTC (fraction)") # Table 1 sigma parent drug in vaginal compartment = 24.2 %CV
    propSd_Cvagina_ftctp <- 0.831; label("Proportional residual error, vaginal tissue FTC-TP (fraction)") # Table 1 sigma metabolite in vaginal compartment = 83.1 %CV
    propSd_Ccervix <- 0.669; label("Proportional residual error, cervical tissue FTC (fraction)") # Table 1 sigma parent drug in cervical compartment = 66.9 %CV
    propSd_Ccervix_ftctp <- 0.582; label("Proportional residual error, cervical tissue FTC-TP (fraction)") # Table 1 sigma metabolite in cervical compartment = 58.2 %CV
    propSd_Crectum <- 0.0151; label("Proportional residual error, rectal tissue FTC (fraction)") # Table 1 sigma parent drug in rectal compartment = 1.51 %CV
    propSd_Crectum_ftctp <- 0.819; label("Proportional residual error, rectal tissue FTC-TP (fraction)") # Table 1 sigma 'metabolite in plasma' row = 81.9 %CV (EPS(7), rectal metabolite)
  })

  model({
    # Individual parameters (exponential IIV, Leung 2023 Results)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq)
    cl <- exp(lcl)

    fcl_vagina <- exp(lfcl_vagina)
    fcl_cervix <- exp(lfcl_cervix)
    fcl_rectum <- exp(lfcl_rectum)

    cl_vagina <- exp(lcl_vagina + etalcl_vagina)
    cl_vagina_ftctp <- exp(lcl_vagina_ftctp + etalcl_vagina_ftctp)
    v_vagina <- exp(lv_vagina)
    v_vagina_ftctp <- exp(lv_vagina_ftctp)

    cl_cervix <- exp(lcl_cervix + etalcl_cervix)
    cl_cervix_ftctp <- exp(lcl_cervix_ftctp + etalcl_cervix_ftctp)
    v_cervix <- exp(lv_cervix)
    v_cervix_ftctp <- exp(lv_cervix_ftctp)

    cl_rectum <- exp(lcl_rectum + etalcl_rectum)
    cl_rectum_ftctp <- exp(lcl_rectum_ftctp)
    v_rectum <- exp(lv_rectum)
    v_rectum_ftctp <- exp(lv_rectum_ftctp)

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
    d/dt(vagina_ftctp) <- fm_vagina * cl_vagina / v_vagina * vagina - cl_vagina_ftctp / v_vagina_ftctp * vagina_ftctp
    d/dt(cervix_ftctp) <- fm_cervix * cl_cervix / v_cervix * cervix - cl_cervix_ftctp / v_cervix_ftctp * cervix_ftctp
    d/dt(rectum_ftctp) <- fm_rectum * cl_rectum / v_rectum * rectum - cl_rectum_ftctp / v_rectum_ftctp * rectum_ftctp

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
    Cvagina_ftctp <- vagina_ftctp / v_vagina_ftctp
    Ccervix <- cervix / v_cervix
    Ccervix_ftctp <- cervix_ftctp / v_cervix_ftctp
    Crectum <- rectum / v_rectum
    Crectum_ftctp <- rectum_ftctp / v_rectum_ftctp

    Cc ~ prop(propSd)
    Cvagina ~ prop(propSd_Cvagina)
    Cvagina_ftctp ~ prop(propSd_Cvagina_ftctp)
    Ccervix ~ prop(propSd_Ccervix)
    Ccervix_ftctp ~ prop(propSd_Ccervix_ftctp)
    Crectum ~ prop(propSd_Crectum)
    Crectum_ftctp ~ prop(propSd_Crectum_ftctp)
  })
}
