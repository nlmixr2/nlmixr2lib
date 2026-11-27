Leil_2014_4betaHydroxycholesterol <- function() {
  description <- "Bayesian mechanism-based PK/PD model of CYP3A4 activity in healthy adults (Leil 2014), fit jointly to midazolam plasma concentrations and to plasma total cholesterol, 4beta-hydroxycholesterol (4betaHC) and 4alpha-hydroxycholesterol (4alphaHC) during 14 days of ketoconazole 400 mg or rifampin 600 mg once daily. Midazolam: oral dose through a three-state gut chain whose middle state carries the gut-wall CYP3A4 first-pass loss, two-compartment disposition, well-stirred hepatic clearance from a Michaelis-Menten intrinsic clearance. Ketoconazole (not measured, highly informative priors): oral two-compartment PK with saturable well-stirred hepatic clearance; its unbound concentration competitively inhibits midazolam intrinsic clearance in gut and liver (Ki 0.125 nmol/L) and reversibly inhibits 4betaHC formation (Ki 50 nmol/L). Rifampin: a step-on induction of liver (3.43-fold) and gut (29.6-fold) CYP3A4 synthesis gated by the CONMED_RIFAMPICIN indicator, with CYP3A4 turning over at a 70 h half-life. Endogenous: total cholesterol turnover feeds CYP3A4-dependent 4betaHC formation and CYP3A4-independent 4alphaHC formation."
  reference <- paste(
    "Leil TA, Kasichayanula S, Boulton DW, LaCreta F. Evaluation of",
    "4beta-Hydroxycholesterol as a Clinical Biomarker of CYP3A4 Drug",
    "Interactions Using a Bayesian Mechanism-Based Pharmacometric Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2014;3(6):e120.",
    "doi:10.1038/psp.2014.18.",
    "Structure from the NONMEM control stream in the Supplementary Data",
    "(MCMC chain 1); posterior parameter values from Table 1.",
    sep = " "
  )
  vignette <- "Leil_2014_4betaHydroxycholesterol"
  dosing <- c("depot", "central", "depot_ketoconazole")
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "midazolam", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "midazolam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "mg", specimen = "tissue", verified = TRUE),
    depot_ketoconazole = list(
      analyte = "ketoconazole",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ketoconazole = list(analyte = "ketoconazole", units = "mg", specimen = "tissue", verified = TRUE),
    total_cholesterol = list(analyte = "total cholesterol", units = "umol/L", specimen = "plasma", verified = TRUE),
    hydroxycholesterol_4beta = list(
      analyte = "4beta-hydroxycholesterol",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    hydroxycholesterol_4alpha = list(
      analyte = "4alpha-hydroxycholesterol",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    enzyme_liver = list(
      analyte = "CYP3A4 (hepatic), relative to baseline",
      units = "(unitless)",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_gut = list(
      analyte = "CYP3A4 (gut wall), relative to baseline",
      units = "(unitless)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONMED_RIFAMPICIN = list(
      description = "Rifampin (rifampicin) co-administration indicator: 1 while the subject is inside the rifampin dosing window, 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampin)",
      notes = "Time-varying within subject. Rifampin PK is not modelled; the indicator switches on the CYP3A4 induction drive (1 + foldind * (1 - exp(-kdeg_enzyme * tind))) on the liver and gut enzyme synthesis rates. In the source NONMEM code the drive is on for the rifampin arm (TRT = 2) between INDSTART = 48 h and INDEND = 360 h (600 mg once daily for 2 weeks, starting on day 3); set CONMED_RIFAMPICIN = 1 over the same window. Induction then builds up and washes out through the 70 h CYP3A4 turnover."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "18-45 years",
    weight_range = "(not reported)",
    sex_female_pct = NA_real_,
    race_ethnicity = "(not reported)",
    disease_state = "Healthy adult volunteers.",
    dose_range = "Midazolam 2 mg oral and 0.4 mg i.v. on days 1-2 and again on days 15-16, after 14 days of ketoconazole 400 mg once daily, rifampin 600 mg once daily, or placebo.",
    regions = "(not reported)",
    notes = "Randomized three-arm clinical study (Leil 2014 Methods, ref. 23). Analysis data: 1,600 midazolam, 606 total-cholesterol, 606 4betaHC and 605 4alphaHC plasma concentrations; predose 4betaHC / 4alphaHC / cholesterol samples daily from day 1 to day 30. Ketoconazole concentrations were not measured. Fit by NONMEM VII MCMC Bayesian estimation with informative priors from published midazolam and ketoconazole meta-analyses."
  )

  ini({
    # ---------------------------------------------------------------
    # Midazolam PK -- Table 1 posterior means; structure from the
    # Supplementary NONMEM $PK / $DES blocks.
    # ---------------------------------------------------------------
    lfpe <- log(1.0)
    label("Midazolam gut first-pass scaling factor FPE (unitless)")
    # Table 1 'F_M MDZ bioavailability' posterior mean 1.0 (SD 0.1). In the
    # control stream this parameter (FPE) multiplies the gut-wall
    # clearance in the first-pass transit state; it is not a bioavailability.
    lka <- log(3.91)
    label("Midazolam absorption rate constant KAm (1/h)") # Table 1 K_Am 3.91 (0.42)
    lq <- log(13.4)
    label("Midazolam intercompartmental clearance Qm (L/h)") # Table 1 Q_m 13.4 (1.43)
    lvc <- log(45.3)
    label("Midazolam central volume VCm (L)") # Table 1 V_Cm 45.3 (6.4)
    lvp <- log(114)
    label("Midazolam peripheral volume VPm (L)") # Table 1 V_Pm 114 (18.4)
    lvmax <- log(1100)
    label("Midazolam Vmax at baseline CYP3A4 (umol/h)") # Table 1 V_max_m 1,100 (122)
    lkm <- log(2.10)
    label("Midazolam Km (umol/L)") # Table 1 K_m_m 2.10 (0.064)
    lfu <- log(0.040)
    label("Midazolam fraction unbound in plasma (unitless)") # Table 1 F_us 0.040 (0.0009)
    lqh <- log(80)
    label("Hepatic blood flow QH (L/h)") # Table 1 'Q Hepatic blood flow' 80 (0.078)
    lmtt <- log(0.544)
    label("Midazolam mean gut transit time MTT (h)") # Table 1 MTT 0.544 (0.057)

    # ---------------------------------------------------------------
    # Ketoconazole PK (perpetrator, not measured) -- Table 1.
    # ---------------------------------------------------------------
    logitfdepot_ketoconazole <- log(0.70 / (1 - 0.70))
    label("Ketoconazole oral bioavailability FK, logit scale (unitless)")
    # Table 1 F_K 0.70 (0.0002). The control stream estimates FK on the logit
    # scale as F4 = 1 - expit(THETA(10) + ETA(10)); the equivalent form
    # expit(logit(0.70) + eta) is used here.
    lka_ketoconazole <- log(0.61)
    label("Ketoconazole absorption rate constant KA (1/h)") # Table 1 K_A 0.61 (0.005)
    lq_ketoconazole <- log(3.22)
    label("Ketoconazole intercompartmental clearance Q (L/h)") # Table 1 Q 3.22 (0.1)
    lvc_ketoconazole <- log(20.4)
    label("Ketoconazole central volume VC (L)") # Table 1 V_C 20.4 (0.7)
    lvp_ketoconazole <- log(24.5)
    label("Ketoconazole peripheral volume VP (L)") # Table 1 V_P 24.5 (0.79)
    lvmax_ketoconazole <- log(23.0)
    label("Ketoconazole Vmax (umol/h)") # Table 1 V_max 23.0 (0.70)
    lkm_ketoconazole <- log(0.48)
    label("Ketoconazole Km (umol/L)") # Table 1 K_m 0.48 (0.004)
    lfu_ketoconazole <- log(0.029)
    label("Ketoconazole fraction unbound in plasma Fuk (unitless)") # Table 1 F_uk 0.029 (0.0007)

    # ---------------------------------------------------------------
    # Endogenous cholesterol / hydroxycholesterol turnover -- Table 1.
    # ---------------------------------------------------------------
    lrbase_hydroxycholesterol_4beta <- log(81.2)
    label("Baseline 4beta-hydroxycholesterol B_BHCH (nmol/L)") # Table 1 B_BHCH 81.2 (4.7)
    lthalf_hydroxycholesterol_4beta <- log(124)
    label("Half-life of 4beta-hydroxycholesterol BHC_HL (h)") # Table 1 BHC_HL 124 (7.6)
    lthalf_total_cholesterol <- fixed(log(192))
    label("Half-life of total cholesterol CH_HL (h)") # Table 1 CH_HL 192 (0); THETA 192 FIX in the control stream
    lrbase_total_cholesterol <- log(4244)
    label("Baseline total cholesterol BCH (umol/L)") # Table 1 BCH 4,244 (127)
    lrbase_hydroxycholesterol_4alpha <- log(10.6)
    label("Baseline 4alpha-hydroxycholesterol B_AHCH (nmol/L)") # Table 1 B_AHCH 10.6 (0.53)

    # ---------------------------------------------------------------
    # CYP3A4 turnover, inhibition and induction -- Table 1.
    # ---------------------------------------------------------------
    lthalf_enzyme <- fixed(log(70))
    label("CYP3A4 half-life CYP_HL (h)") # Table 1 CYP_HL 70 (0); THETA 70 FIX in the control stream
    lki <- log(0.125)
    label("Ketoconazole unbound Ki for midazolam intrinsic clearance (nmol/L)")
    # Table 1 Ki 0.125 (0.226) nmol/L; Results text 'Ki = 0.125 nmol/l'.
    lki_hydroxycholesterol_4beta <- log(50)
    label("Ketoconazole unbound Ki for 4beta-hydroxycholesterol formation (nmol/L)")
    # Table 1 K_i4B 50 (17) nmol/L; Results text 'Ki = 50 nmol/l'.
    lfoldind_liver <- log(3.43)
    label("Rifampin fold induction of hepatic (systemic) CYP3A4 FOLDINDS (unitless)") # Table 1 FOLDINDS 3.43 (0.46)
    lfoldind_gut <- log(29.6)
    label("Rifampin fold induction of gut CYP3A4 FOLDINDG (unitless)") # Table 1 FOLDINDG 29.6 (1.5)
    tind <- fixed(312)
    label("Duration of the rifampin induction window TIND (h)")
    # Control stream: TIND = INDEND - INDSTART = 360 - 48 = 312 h (fixed design constant).

    # ---------------------------------------------------------------
    # Between-subject variability. Table 1 'Posterior ISV' is read as the
    # SD of the eta (omega); squared here to give the variance. All etas
    # are log-normal except FK, which is additive on the logit scale.
    # ---------------------------------------------------------------
    etalfpe ~ 0.1764 # Table 1 posterior ISV F_M 0.42 -> 0.42^2
    etalka ~ 0.3721 # Table 1 posterior ISV K_Am 0.61 -> 0.61^2
    etalq ~ 0.2704 # Table 1 posterior ISV Q_m 0.52 -> 0.52^2
    etalvc ~ 0.7569 # Table 1 posterior ISV V_Cm 0.87 -> 0.87^2
    etalvp ~ 2.56 # Table 1 posterior ISV V_Pm 1.6 -> 1.6^2
    etalvmax ~ 0.2304 # Table 1 posterior ISV V_max_m 0.48 -> 0.48^2
    etalkm ~ 0.001369 # Table 1 posterior ISV K_m_m 0.037 -> 0.037^2
    etalfu ~ 0.001369 # Table 1 posterior ISV F_us 0.037 -> 0.037^2
    etalqh ~ 0.001521 # Table 1 posterior ISV Q 0.039 -> 0.039^2
    etalmtt ~ 0.2401 # Table 1 posterior ISV MTT 0.49 -> 0.49^2
    etalogitfdepot_ketoconazole ~ 0.001369 # Table 1 posterior ISV F_K 0.037 -> 0.037^2
    etalka_ketoconazole ~ 0.001369 # Table 1 posterior ISV K_A 0.037 -> 0.037^2
    etalq_ketoconazole ~ 0.001369 # Table 1 posterior ISV Q (ketoconazole) 0.037 -> 0.037^2
    etalvc_ketoconazole ~ 0.001369 # Table 1 posterior ISV V_C 0.037 -> 0.037^2
    etalvp_ketoconazole ~ 0.001369 # Table 1 posterior ISV V_P 0.037 -> 0.037^2
    etalvmax_ketoconazole ~ 0.001369 # Table 1 posterior ISV V_max 0.037 -> 0.037^2
    etalkm_ketoconazole ~ 0.001369 # Table 1 posterior ISV K_m 0.037 -> 0.037^2
    etalfu_ketoconazole ~ 0.001024 # Table 1 posterior ISV F_uk 0.032 -> 0.032^2
    etalrbase_hydroxycholesterol_4beta ~ 0.1024 # Table 1 posterior ISV B_BHCH 0.32 -> 0.32^2
    etalthalf_hydroxycholesterol_4beta ~ 0.04 # Table 1 posterior ISV BHC_HL 0.20 -> 0.20^2
    etalthalf_total_cholesterol ~ fixed(0.001) # $OMEGA 0.001 FIX; Table 1 posterior ISV CH_HL 0.031
    etalrbase_total_cholesterol ~ 0.0225 # Table 1 posterior ISV BCH 0.15 -> 0.15^2
    etalrbase_hydroxycholesterol_4alpha ~ 0.0484 # Table 1 posterior ISV B_AHCH 0.22 -> 0.22^2
    etalthalf_enzyme ~ fixed(0.001) # $OMEGA 0.001 FIX; Table 1 posterior ISV CYP_HL 0.032
    etalki ~ 4.84 # Table 1 posterior ISV Ki 2.2 -> 2.2^2
    etalki_hydroxycholesterol_4beta ~ 0.49 # Table 1 posterior ISV K_i4B 0.70 -> 0.70^2
    etalfoldind_liver ~ 0.1764 # Table 1 posterior ISV FOLDINDS 0.42 -> 0.42^2
    etalfoldind_gut ~ 0.001024 # Table 1 posterior ISV FOLDINDG 0.032 -> 0.032^2

    # ---------------------------------------------------------------
    # Residual error -- Table 1 sigma rows (SDs).
    # ---------------------------------------------------------------
    propSd <- 0.36
    label("Proportional residual SD, midazolam (fraction)") # Table 1 sigma_1 0.36 (0.077)
    addSd <- 0.046
    label("Additive residual SD, midazolam (nmol/L)") # Table 1 sigma_2 0.046 (0.037)
    addSd_total_cholesterol <- 1170
    label("Additive residual SD, total cholesterol (umol/L)")
    # Table 1 sigma_3 1,170 (284); printed unit 'nmol/l' must be umol/L -- see vignette.
    addSd_hydroxycholesterol_4beta <- 28.3
    label("Additive residual SD, 4beta-hydroxycholesterol (nmol/L)") # Table 1 sigma_4 28.3 (7.2)
    addSd_hydroxycholesterol_4alpha <- 4.19
    label("Additive residual SD, 4alpha-hydroxycholesterol (nmol/L)") # Table 1 sigma_5 4.19 (1.02)
  })
  model({
    # Individual parameters ------------------------------------------------
    fpe <- exp(lfpe + etalfpe)
    ka <- exp(lka + etalka)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm + etalkm)
    fu <- exp(lfu + etalfu)
    qh <- exp(lqh + etalqh)
    mtt <- exp(lmtt + etalmtt)
    ktr <- 3 / mtt

    fdepot_ketoconazole <- expit(logitfdepot_ketoconazole + etalogitfdepot_ketoconazole)
    ka_ketoconazole <- exp(lka_ketoconazole + etalka_ketoconazole)
    q_ketoconazole <- exp(lq_ketoconazole + etalq_ketoconazole)
    vc_ketoconazole <- exp(lvc_ketoconazole + etalvc_ketoconazole)
    vp_ketoconazole <- exp(lvp_ketoconazole + etalvp_ketoconazole)
    vmax_ketoconazole <- exp(lvmax_ketoconazole + etalvmax_ketoconazole)
    km_ketoconazole <- exp(lkm_ketoconazole + etalkm_ketoconazole)
    fu_ketoconazole <- exp(lfu_ketoconazole + etalfu_ketoconazole)

    rbase_hydroxycholesterol_4beta <- exp(lrbase_hydroxycholesterol_4beta + etalrbase_hydroxycholesterol_4beta)
    thalf_hydroxycholesterol_4beta <- exp(lthalf_hydroxycholesterol_4beta + etalthalf_hydroxycholesterol_4beta)
    thalf_total_cholesterol <- exp(lthalf_total_cholesterol + etalthalf_total_cholesterol)
    rbase_total_cholesterol <- exp(lrbase_total_cholesterol + etalrbase_total_cholesterol)
    rbase_hydroxycholesterol_4alpha <- exp(lrbase_hydroxycholesterol_4alpha + etalrbase_hydroxycholesterol_4alpha)

    thalf_enzyme <- exp(lthalf_enzyme + etalthalf_enzyme)
    ki <- exp(lki + etalki)
    ki_hydroxycholesterol_4beta <- exp(lki_hydroxycholesterol_4beta + etalki_hydroxycholesterol_4beta)
    foldind_liver <- exp(lfoldind_liver + etalfoldind_liver)
    foldind_gut <- exp(lfoldind_gut + etalfoldind_gut)

    # Turnover rates (the control stream uses 0.693 for log(2)) --------------
    kout_hydroxycholesterol_4beta <- 0.693 / thalf_hydroxycholesterol_4beta
    kout_total_cholesterol <- 0.693 / thalf_total_cholesterol
    kin_total_cholesterol <- rbase_total_cholesterol * kout_total_cholesterol
    # The control stream sets the 4alphaHC elimination rate from the
    # cholesterol half-life (KELA = 0.693/CH_HL), not from AHC_HL.
    kout_hydroxycholesterol_4alpha <- 0.693 / thalf_total_cholesterol
    kform_hydroxycholesterol_4beta <- kout_hydroxycholesterol_4beta * rbase_hydroxycholesterol_4beta / rbase_total_cholesterol
    kform_hydroxycholesterol_4alpha <- kout_hydroxycholesterol_4alpha * rbase_hydroxycholesterol_4alpha / rbase_total_cholesterol
    kdeg_enzyme <- 0.693 / thalf_enzyme

    # Ketoconazole concentrations and CYP3A4 inhibition ----------------------
    Cc_ketoconazole <- central_ketoconazole / vc_ketoconazole # ug/mL
    # unbound ketoconazole in nmol/L (MW 531 g/mol)
    cu_ketoconazole <- Cc_ketoconazole * 1000 / 531 * fu_ketoconazole * 1000
    inh <- 1 + cu_ketoconazole / ki
    inh_hydroxycholesterol_4beta <- 1 - cu_ketoconazole / (cu_ketoconazole + ki_hydroxycholesterol_4beta)

    # Rifampin induction drive on CYP3A4 synthesis --------------------------
    ind_liver <- 1 + foldind_liver * CONMED_RIFAMPICIN * (1 - exp(-kdeg_enzyme * tind))
    ind_gut <- 1 + foldind_gut * CONMED_RIFAMPICIN * (1 - exp(-kdeg_enzyme * tind))

    # Midazolam clearance (well-stirred, Michaelis-Menten intrinsic) ---------
    cmdz_um <- central / vc * 1000 / 326 # midazolam in umol/L (MW 326 g/mol)
    clint_liver <- fu * vmax * enzyme_liver / (km * inh + fu * cmdz_um)
    clint_gut <- fu * vmax * enzyme_gut / (km * inh + fu * cmdz_um)
    cl_liver <- qh * clint_liver / (qh + clint_liver)
    cl_gut <- qh * clint_gut / (qh + clint_gut)

    # Ketoconazole saturable clearance ---------------------------------------
    clint_ketoconazole <- vmax_ketoconazole / (km_ketoconazole + Cc_ketoconazole)
    cl_ketoconazole <- qh * clint_ketoconazole / (qh + clint_ketoconazole)

    # ODEs -------------------------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1 - cl_gut * fpe * transit1
    d/dt(transit2) <- ktr * transit1 - ka * transit2
    d/dt(central) <- ka * transit2 + q / vp * peripheral1 - q / vc * central - cl_liver / vc * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    d/dt(depot_ketoconazole) <- -ka_ketoconazole * depot_ketoconazole
    d/dt(central_ketoconazole) <- ka_ketoconazole * depot_ketoconazole -
      q_ketoconazole / vc_ketoconazole * central_ketoconazole +
      q_ketoconazole / vp_ketoconazole * peripheral1_ketoconazole -
      cl_ketoconazole / vc_ketoconazole * central_ketoconazole
    # As published (DADT(6) = Q/VC*A(4) - Q/VP*A(6)): the peripheral input is
    # driven by the depot amount, not the central amount.
    d/dt(peripheral1_ketoconazole) <- q_ketoconazole / vc_ketoconazole * depot_ketoconazole -
      q_ketoconazole / vp_ketoconazole * peripheral1_ketoconazole

    d/dt(total_cholesterol) <- kin_total_cholesterol - kout_total_cholesterol * total_cholesterol
    d/dt(hydroxycholesterol_4beta) <- kform_hydroxycholesterol_4beta * enzyme_liver * inh_hydroxycholesterol_4beta * total_cholesterol - kout_hydroxycholesterol_4beta * hydroxycholesterol_4beta
    d/dt(hydroxycholesterol_4alpha) <- kform_hydroxycholesterol_4alpha * total_cholesterol - kout_hydroxycholesterol_4alpha * hydroxycholesterol_4alpha
    d/dt(enzyme_liver) <- kdeg_enzyme * ind_liver - kdeg_enzyme * enzyme_liver
    d/dt(enzyme_gut) <- kdeg_enzyme * ind_gut - kdeg_enzyme * enzyme_gut

    f(depot_ketoconazole) <- fdepot_ketoconazole

    total_cholesterol(0) <- rbase_total_cholesterol
    hydroxycholesterol_4beta(0) <- rbase_hydroxycholesterol_4beta
    hydroxycholesterol_4alpha(0) <- rbase_hydroxycholesterol_4alpha
    enzyme_liver(0) <- 1
    enzyme_gut(0) <- 1

    # Observations -----------------------------------------------------------
    Cc <- central / vc * 1000 / 0.326 # midazolam, nmol/L
    Cc ~ add(addSd) + prop(propSd)
    total_cholesterol ~ add(addSd_total_cholesterol)
    hydroxycholesterol_4beta ~ add(addSd_hydroxycholesterol_4beta)
    hydroxycholesterol_4alpha ~ add(addSd_hydroxycholesterol_4alpha)
  })
}
