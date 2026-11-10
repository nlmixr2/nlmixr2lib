Wan_2025_sunitinib_qsp <- function() {
  description <- "QSP (human virtual population). QSP-PK-TD model of sunitinib-induced left ventricular systolic dysfunction (LVEF decline) from the Wan 2025 in vitro-to-in vivo extrapolation of hiPSC-cardiomyocyte toxicity for 21 tyrosine kinase inhibitors. Population two-compartment plasma PK with lagged first-order absorption (Wan 2025 Table S2 typical values and CVs, taken from Wang 2020, Cancer Chemother Pharmacol 86:181). The plasma concentration is converted to an effective myocardial concentration Cin = Cp x fu x Kp x H (Eqs. 17, 18, 21; fu Table S5, Kp Table S7, H Table S8). Cin drives an in vitro toxicodynamic model fitted to hiPSC-CM data (Table 1): threshold-gated cell injury through three transit stages (Eqs. 1-8), inhibition of mitochondrial membrane potential production driving ATP (Eqs. 9-11), and a drug-raised ATP50 for average contractile force (Eqs. 12-14). Surviving fraction x average contractile force replaces the drug effect on stroke-volume production in the Sang 2021 human cardiovascular QSP model (SV, LVEDV, HR and TPR turnover with MAP and LVESV feedback), with population-specific baselines and FB_LVESV for cardiovascular-healthy, hypertensive and dilated-cardiomyopathy patients (Table 2). Wan 2025 classifies sunitinib as cardiotoxic (Figures 1-3)."
  reference <- "Wan Z, Wang C, Luo S, Zhu J, He H, Hao K (2025). Bridging the Gap Between hiPSC-CMs Cardiotoxicity Assessment and Clinical LVEF Decline Risk: A Case Study of 21 Tyrosine Kinase Inhibitors. Pharmaceuticals (Basel) 18(4):450. doi:10.3390/ph18040450. Cardiovascular QSP structure and dissipation rate constants from Sang L et al. (2021) CPT Pharmacometrics Syst Pharmacol 10(12):1512-1524 (doi:10.1002/psp4.12719; Sang_2021_doxorubicin_human_qsp); tau, kin_ATP, kin_MMP and n from Sang L, Zhou Z, Luo S, Zhang Y, Qian H, Zhou Y, He H, Hao K (2024) Pharm Res 41(2):247-262 (doi:10.1007/s11095-023-03644-4), Table I."
  vignette <- "Wan_2025_tki_cardiotoxicity"

  paper_specific_compartments <- c(
    "fnor",
    "finj1",
    "finj2",
    "finj3",
    "mmp",
    "atp",
    "sv",
    "edv",
    "tpr",
    "aemmp",
    "aeaac"
  )

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L (plasma Cc); the effective myocardial concentration driving the toxicodynamics is in umol/L"
  )

  covariateData <- list(
    LVEDV_BL = list(
      description = "Baseline (pre-treatment) left ventricular end-diastolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Wan 2025 Table 2 typical values 113 mL (cardiovascular healthy), 100.95 mL (hypertension) and 146 mL (dilated cardiomyopathy), CV 10% / 20% / 20%. Initial condition of the LVEDV state.",
      source_name = "LVEDV0"
    ),
    LVESV_BL = list(
      description = "Baseline (pre-treatment) left ventricular end-systolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Wan 2025 parameterises the baseline as SV0 = Index_SV0 x LVEDV0 (Table 2: Index_SV0 0.6 healthy and hypertension, 0.55 dilated cardiomyopathy; CV 5%), so LVESV_BL = LVEDV_BL x (1 - Index_SV0). Baseline SV = LVEDV_BL - LVESV_BL.",
      source_name = "LVEDV0 x (1 - Index_SV0)"
    ),
    HR_BL = list(
      description = "Baseline (pre-treatment) heart rate",
      units = "beats/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Wan 2025 Table 2 HR0 70 (healthy and hypertension) and 57 (dilated cardiomyopathy); no CV given. Initial condition of the HR state.",
      source_name = "HR0"
    ),
    MAP_BL = list(
      description = "Baseline (pre-treatment) mean arterial pressure",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Wan 2025 Table 2 MAP0 94.92 (healthy), 106 (hypertension; at or above 106 mmHg) and 91.54 mmHg (dilated cardiomyopathy), CV 10% / 20% / 30%. Baseline TPR = MAP_BL / (HR_BL x SV_BL), which reproduces Table 2 TPR0.",
      source_name = "MAP0"
    ),
    DIS_HYPERT = list(
      description = "Hypertension comorbidity indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (cardiovascular healthy)",
      notes = "Selects the hypertension FB_LVESV of Wan 2025 Table 2 (0.0045 vs 0.0042 1/mL). Wan defines the hypertensive virtual population by MAP0 at or above 106 mmHg (35% of the virtual cohort). Mutually exclusive with DIS_DCM; if both are 1, DIS_DCM takes precedence.",
      source_name = "Hypertension population (Table 2 column)"
    ),
    DIS_DCM = list(
      description = "Dilated cardiomyopathy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no dilated cardiomyopathy)",
      notes = "Selects the dilated-cardiomyopathy FB_LVESV of Wan 2025 Table 2 (0.0028 vs 0.0042 1/mL). Wan defines the population by LVEDV0 at or above 146 mL (0.03% of the virtual cohort).",
      source_name = "Dilated Cardiomyopathy population (Table 2 column)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib", units = "mg", specimen = "tissue", verified = TRUE),
    fnor = list(
      analyte = "normal cardiomyocyte fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    finj1 = list(
      analyte = "early-phase injured cardiomyocyte fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    finj2 = list(
      analyte = "middle-phase injured cardiomyocyte fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    finj3 = list(
      analyte = "late-phase injured cardiomyocyte fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    mmp = list(
      analyte = "mitochondrial membrane potential",
      units = "fraction of baseline",
      specimen = "not applicable",
      verified = TRUE
    ),
    atp = list(
      analyte = "cardiomyocyte ATP",
      units = "fraction of baseline",
      specimen = "not applicable",
      verified = TRUE
    ),
    sv = list(analyte = "stroke volume", units = "mL", specimen = "not applicable", verified = TRUE),
    edv = list(
      analyte = "left ventricular end-diastolic volume",
      units = "mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    hr = list(analyte = "heart rate", units = "beats/min", specimen = "not applicable", verified = TRUE),
    tpr = list(
      analyte = "total peripheral resistance",
      units = "mmHg*min/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aemmp = list(
      analyte = "area under the EMMP effect curve",
      units = "h",
      specimen = "not applicable",
      verified = TRUE
    ),
    aeaac = list(
      analyte = "area under the EAAC effect curve",
      units = "h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15000,
    disease_state = "Virtual adult cancer patients treated with sunitinib (50 mg, QD), in three cardiovascular populations: cardiovascular healthy, hypertension (MAP0 >= 106 mmHg) and dilated cardiomyopathy (LVEDV0 >= 146 mL)",
    dose_range = "50 mg, QD for one year (Wan 2025 Table S4)",
    notes = "Virtual population only; no human data were fitted. Wan 2025 generated 5000 virtual patients per cardiovascular population by Monte Carlo sampling (Table 2 typical values and CVs) and weighted the incidence 64.97% healthy, 35% hypertension, 0.03% dilated cardiomyopathy (Eq. 16). The PK model is a published population model condensed to typical values and CVs in Table S2; the toxicodynamic parameters were fitted to digitised high-throughput hiPSC-CM viability, ATP and contractility data (Sharma 2017) after 72 h exposure."
  )

  ini({
    # PK (Wan 2025 Table S2; source model: Wang 2020, Cancer Chemother Pharmacol 86:181). Wan used the typical values
    # as fixed inputs, so every structural value is fixed().
    lka <- fixed(log(0.38)); label("Absorption rate constant, Ka (1/h)") # Table S2 Sunitinib Ka = 0.38 1/h
    ltlag <- fixed(log(0.64)); label("Absorption lag time, Tlag (h)") # Table S2 Sunitinib Tlag = 0.64 h
    lcl <- fixed(log(24.1)); label("Clearance, CL (L/h)") # Table S2 Sunitinib CL = 24.1 L/h
    lvc <- fixed(log(1070)); label("Central volume, Vc (L)") # Table S2 Sunitinib Vc = 1070 L
    lq <- fixed(log(0.28)); label("Intercompartmental clearance, Q (L/h)") # Table S2 Sunitinib Q = 0.28 L/h
    lvp <- fixed(log(63.8)); label("Peripheral volume, Vp (L)") # Table S2 Sunitinib Vp = 63.8 L
    lfdepot <- fixed(log(1)); label("Bioavailability, F (unitless)") # Table S2 Sunitinib F = 1

    # Plasma-to-myocardium exposure factors (Wan 2025 Section 4.5)
    mw <- fixed(398); label("Molecular weight (g/mol)") # Table S6 Sunitinib molecular weight = 398 g/mol
    fu <- fixed(0.05); label("Fraction unbound in plasma, fu (unitless)") # Table S5 Sunitinib fu,plasma = 5%
    lkp_heart <- fixed(log(2.2)); label("Heart-to-plasma partition coefficient, Kp (unitless)") # Table S7 Sunitinib Kp_observation = 2.2 (Kp_prediction 3.85 not used: an observed value is reported)
    h_het <- fixed(0.234); label("Heterogeneous-distribution factor H = IC50,2D / IC50,3D (unitless)") # Table S8 Sunitinib IC50,2D 3.8 / IC50,3D 16.24 uM = 0.234 (Eq. 21)

    # In vitro toxicodynamic parameters for sunitinib (Wan 2025 Table 1)
    lkinj <- fixed(log(0.1)); label("Rate constant of conversion of normal to injured cardiomyocytes, kinj (1/h)") # Table 1 Sunitinib kinj = 10 x 1e-2 (header unit printed as uM; Eq. 2 requires 1/h)
    lec50_inj <- fixed(log(22.02)); label("Cin above threshold giving half-maximal injury, EC50inj (umol/L)") # Table 1 Sunitinib EC50inj = 22.02 uM
    hill_inj <- fixed(3); label("Hill coefficient of the injury effect, hill_inj (unitless)") # Table 1 Sunitinib hill_inj = 3
    lcth <- fixed(log(2.2)); label("Threshold concentration for cell injury, Cth (umol/L)") # Table 1 Sunitinib Cth = 2.2 uM
    lec50_mmp <- fixed(log(5.01)); label("Concentration halving MMP production, EC50MMP (umol/L)") # Table 1 Sunitinib EC50MMP = 5.01 uM
    leec <- fixed(log(1.42)); label("Concentration giving half-maximal effect on average contractile force, EEC (umol/L)") # Table 1 Sunitinib EEC = 1.42 uM

    # System-dependent toxicodynamic parameters shared by all TKIs
    emax_inj <- fixed(1); label("Maximal injury effect, Emax_inj (unitless)") # Section 4.2.1: Emax_inj fixed at 1
    ltau_inj <- fixed(log(4.9)); label("Mean transit time of each injured-cell stage, tau (h)") # Sang 2024 Table I DOX-TRZ tau = 4.9 h (adopted per Wan Section 4.2.4)
    lkin_mmp <- fixed(log(0.96)); label("Zero-order MMP production rate constant, kin_MMP (1/h)") # Sang 2024 Table I kin_MMP = 0.96 1/h (adopted per Wan Section 4.2.4)
    lkin_atp <- fixed(log(0.032)); label("ATP production rate constant, kin_ATP (1/h)") # Sang 2024 Table I kin_ATP = 0.032 1/h (adopted per Wan Section 4.2.4)
    n_atp <- fixed(0.8); label("Power of MMP on ATP production, n (unitless)") # Sang 2024 Table I n = 0.8 (adopted per Wan Section 4.2.4)
    latp50_ref <- fixed(log(1.35)); label("Relative ATP giving half-maximal contractile force before treatment, ATP50_ref (unitless)") # not reported in Wan 2025 or Sang 2024; back-solved by the maintainers from the Wan Figure 3 model curves of all 21 TKIs (drug-bootstrap 95% interval 0.75-2.31)

    # Cardiovascular QSP (Sang 2021 human model, Table 1; Wan 2025 Section 4.3)
    lkout_sv <- fixed(log(0.0308)); label("Dissipation rate constant of SV, kout_SV (1/h)") # Sang 2021 Table 1 human 0.0308
    lkout_edv <- fixed(log(0.0308)); label("Dissipation rate constant of LVEDV, kout_LVEDV (1/h)") # Sang 2021 Table 1 human 0.0308
    lkout_hr <- fixed(log(2.83)); label("Dissipation rate constant of HR, kout_HR (1/h)") # Sang 2021 Table 1 human 2.83
    lkout_tpr <- fixed(log(0.875)); label("Dissipation rate constant of TPR, kout_TPR (1/h)") # Sang 2021 Table 1 human 0.875
    lfb0 <- fixed(log(0.0029)); label("MAP feedback constant at the reference MAP, FB_MAP_0 (1/mmHg)") # Sang 2021 Table 1 FB_MAP_0 2.9e-3
    e_bslmap_fb <- fixed(-1.98); label("Power exponent of baseline MAP on the MAP feedback constant (unitless)") # Sang 2021 Eq. 8 exponent -1.98
    lfb_lvesv <- fixed(log(0.0042)); label("Feedback of LVESV on LVEDV dissipation, FB_LVESV, cardiovascular healthy (1/mL)") # Wan 2025 Table 2 FB_LVESV Healthy 0.0042
    lfb_lvesv_hypert <- fixed(log(0.0045)); label("FB_LVESV in hypertension (1/mL)") # Wan 2025 Table 2 FB_LVESV Hypertension 0.0045
    lfb_lvesv_dcm <- fixed(log(0.0028)); label("FB_LVESV in dilated cardiomyopathy (1/mL)") # Wan 2025 Table 2 FB_LVESV Dilated Cardiomyopathy 0.0028

    # Between-patient variability: Table S2 CVs (log-normal, omega^2 = log(1 + CV^2))
    # and Table 2 CV_FB_LVESV.
    etalka ~ fixed(0.573462) # Table S2 Sunitinib CVka = 0.88 -> log(1 + 0.88^2)
    etalcl ~ fixed(0.109392) # Table S2 Sunitinib CVCL = 0.34 -> log(1 + 0.34^2)
    etalvc ~ fixed(0.0560022) # Table S2 Sunitinib CVVc = 0.24 -> log(1 + 0.24^2)
    etalfb_lvesv ~ fixed(0.015259) # Table 2 CV_FB_LVESV 0.124 -> log(1 + 0.124^2)

    # Simulation-only model: no residual error is reported.
    propSd <- fixed(0); label("Proportional residual SD on Cc (fraction; not reported)") # simulation-only model, no residual error reported
    addSd_LVEF <- fixed(0); label("Additive residual SD on LVEF (percent; not reported)") # simulation-only model, no residual error reported
    propSd_LVEDV <- fixed(0); label("Proportional residual SD on LVEDV (fraction; not reported)") # simulation-only model, no residual error reported
    propSd_LVESV <- fixed(0); label("Proportional residual SD on LVESV (fraction; not reported)") # simulation-only model, no residual error reported
    addSd_MAP <- fixed(0); label("Additive residual SD on MAP (mmHg; not reported)") # simulation-only model, no residual error reported
  })

  model({
    # ---- PK ----
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc

    # ---- Effective myocardial concentration (umol/L), Eqs. 17, 18, 21 ----
    kp_heart <- exp(lkp_heart)
    cin <- Cc / mw * 1000 * fu * kp_heart * h_het

    # ---- Cardiomyocyte viability (Eqs. 1-8) ----
    kinj <- exp(lkinj)
    ec50_inj <- exp(lec50_inj)
    cth <- exp(lcth)
    tau_inj <- exp(ltau_inj)
    cin_th <- cin - cth
    if (cin_th < 0) {
      cin_th <- 0
    }
    einj <- emax_inj * cin_th^hill_inj / (cin_th^hill_inj + ec50_inj^hill_inj)
    d/dt(fnor) <- -einj * fnor * kinj
    d/dt(finj1) <- einj * fnor * kinj - finj1 / tau_inj
    d/dt(finj2) <- finj1 / tau_inj - finj2 / tau_inj
    d/dt(finj3) <- finj2 / tau_inj - finj3 / tau_inj
    fnor(0) <- 1
    survival <- fnor + finj1 + finj2 + finj3

    # ---- MMP and ATP (Eqs. 9-11); MMP and ATP are normalised to a baseline of 1,
    # so kout = kin keeps the untreated steady state at 1 ----
    kin_mmp <- exp(lkin_mmp)
    kout_mmp <- kin_mmp
    kin_atp <- exp(lkin_atp)
    kout_atp <- kin_atp
    emmp <- cin / (cin + exp(lec50_mmp))
    d/dt(mmp) <- kin_mmp * (1 - emmp) - kout_mmp * mmp
    d/dt(atp) <- kin_atp * mmp^n_atp - kout_atp * atp
    mmp(0) <- 1
    atp(0) <- 1

    # ---- Average contractile force (Eqs. 12-14) ----
    atp50_ref <- exp(latp50_ref)
    eaac <- cin / (cin + exp(leec))
    # Eq. 13 is applied as ATP50 = ATP50_ref / (1 - EAAC): the printed
    # ATP50_ref x (1 - EAAC) would make the drug raise contractility. Eq. 14 is
    # multiplied through by (1 - EAAC) so it stays finite as EAAC -> 1.
    aac <- atp * (1 + atp50_ref) * (1 - eaac) / (atp * (1 - eaac) + atp50_ref)
    etki <- survival * aac

    # Areas under the EMMP / EAAC effect curves (Section 4.4.2); difference
    # these over a steady-state dosing interval to obtain AEMMP / AEAAC
    d/dt(aemmp) <- emmp
    d/dt(aeaac) <- eaac

    # ---- Cardiovascular QSP (Sang 2021 human model) ----
    kout_sv <- exp(lkout_sv)
    kout_edv <- exp(lkout_edv)
    kout_hr <- exp(lkout_hr)
    kout_tpr <- exp(lkout_tpr)
    lfb_pop <- lfb_lvesv
    if (DIS_HYPERT == 1) {
      lfb_pop <- lfb_lvesv_hypert
    }
    if (DIS_DCM == 1) {
      lfb_pop <- lfb_lvesv_dcm
    }
    fb_lvesv <- exp(lfb_pop + etalfb_lvesv)

    # 106.596 mmHg is the reference MAP of the Sang 2021 MAP-feedback relation
    sv_bl <- LVEDV_BL - LVESV_BL
    tpr_bl <- MAP_BL / (HR_BL * sv_bl)
    fb_map <- exp(lfb0) * (MAP_BL / 106.596)^e_bslmap_fb

    kin_sv <- kout_sv * sv_bl / (1 - fb_map * MAP_BL)
    kin_edv <- kout_edv * (1 - fb_lvesv * LVESV_BL) * LVEDV_BL / (1 + fb_map * MAP_BL)
    kin_hr <- kout_hr * HR_BL / (1 - fb_map * MAP_BL)
    kin_tpr <- kout_tpr * tpr_bl / (1 - fb_map * MAP_BL)

    LVESV <- edv - sv
    MAP <- sv * hr * tpr
    LVEF <- sv / edv * 100
    LVEDV <- edv

    # Numerical guard of the Sang 2021 code
    flag <- 0
    if (sv > 0 & LVESV > 0 & LVEF > 20) {
      flag <- 1
    }
    fbmap_clamped <- fb_map * MAP
    if (fbmap_clamped > 1) {
      fbmap_clamped <- 1
    }

    # ETKI replaces the Sang 2021 drug effect on SV production (Eq. 15)
    d/dt(sv) <- (kin_sv * (1 - fb_map * MAP) * etki - kout_sv * sv) * flag
    d/dt(edv) <- (kin_edv * (1 + fbmap_clamped) - kout_edv * (1 - fb_lvesv * LVESV) * edv) * flag
    d/dt(hr) <- (kin_hr * (1 - fbmap_clamped) - kout_hr * hr) * flag
    d/dt(tpr) <- (kin_tpr * (1 - fbmap_clamped) - kout_tpr * tpr) * flag

    sv(0) <- sv_bl
    edv(0) <- LVEDV_BL
    hr(0) <- HR_BL
    tpr(0) <- tpr_bl

    Cc ~ prop(propSd)
    LVEF ~ add(addSd_LVEF)
    LVEDV ~ prop(propSd_LVEDV)
    LVESV ~ prop(propSd_LVESV)
    MAP ~ add(addSd_MAP)
  })
}
