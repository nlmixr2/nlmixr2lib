Sang_2021_doxorubicin_rat_qsp <- function() {
  description <- paste(
    "QSP (rat). Translational QSP-PK-PD model of doxorubicin-induced",
    "cardiac dysfunction in healthy, isoproterenol-hypertrophied and",
    "spontaneously hypertensive rats (Sang 2021). A five-compartment",
    "doxorubicin PK model (intraperitoneal depot, plasma, heart and two",
    "peripheral compartments) drives the cumulative heart-tissue AUC. Four",
    "indirect-response turnover states describe stroke volume (SV), left",
    "ventricular end-diastolic volume (LVEDV), heart rate (HR) and total",
    "peripheral resistance (TPR); MAP = HR * TPR * SV exerts a negative",
    "feedback on the production of SV, HR and TPR and a positive feedback",
    "on LVEDV production, and LVESV = LVEDV - SV feeds back on LVEDV",
    "dissipation. The afterload structure and its dissipation-rate constants",
    "are adapted from Snelder 2014. Heart AUC impairs bioenergy production",
    "(sigmoid, Hill 3, entering SV production as exp(-E)) and myocardial",
    "compliance (Emax, Hill 1); the compliance effect lowers SV production",
    "directly and lowers LVEDV production through a three-compartment",
    "transit delay. The compliance arm is switched off for Sprague-Dawley",
    "rats (systolic-dysfunction phenotype) and on for Wistar-Kyoto and",
    "spontaneously hypertensive rats (diastolic-dysfunction phenotype).",
    "Doses are per-animal amounts in ug; LVEF is reported in percent.",
    sep = " "
  )
  reference <- paste(
    "Sang L, Yuan Y, Zhou Y, Zhou Z, Jiang M, Liu X, Hao K, He H (2021).",
    "A quantitative systems pharmacology approach to predict the",
    "safe-equivalent dose of doxorubicin in patients with cardiovascular",
    "comorbidity. CPT Pharmacometrics Syst Pharmacol. 10(12):1512-1524.",
    "doi:10.1002/psp4.12719. Structure and constants from the deposited",
    "Mlxtran code (Supporting Information, 'Mlxtran code for rat QSP/PK/PD",
    "model') and Tables S4-S5. Afterload model from Snelder N et al. (2014)",
    "Br J Pharmacol 171:5076-5092.",
    sep = " "
  )
  vignette <- "Sang_2021_doxorubicin_cardiotoxicity"

  paper_specific_compartments <- c("sv", "edv", "tpr", "auc_heart")

  units <- list(
    time = "h",
    dosing = "ug (per-animal doxorubicin amount; a mg/kg dose is multiplied by body weight in kg and by 1000)",
    concentration = "ug/mL (plasma Cc and heart-tissue Cheart; the heart AUC driving the PD is in ug*h/mL)"
  )

  covariateData <- list(
    LVEDV_BL = list(
      description = "Baseline (pre-dose) left ventricular end-diastolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Initial condition of the LVEDV state and anchor of its steady-state production rate. Typical rat value 0.385 mL (Sang 2021 Table S5, from Snelder 2014).",
      source_name = "LVEDVbase"
    ),
    LVESV_BL = list(
      description = "Baseline (pre-dose) left ventricular end-systolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Baseline SV = LVEDV_BL - LVESV_BL; also anchors the LVESV feedback on LVEDV dissipation. Typical rat value 0.085 mL ('LVESV0 = 0.085' in the deposited Mlxtran code).",
      source_name = "LVESVbase"
    ),
    MAP_BL = list(
      description = "Baseline (pre-dose) mean arterial pressure",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Baseline TPR = MAP_BL / (HR0 * SV_BL), and the MAP feedback constant scales as (MAP_BL / 106.596)^-1.98. 106.596 mmHg is the typical rat MAP implied by Table S5 (0.3 mL * 423 beats/min * 0.84 mmHg*min/mL).",
      source_name = "MAPbase"
    ),
    STRAIN_SD = list(
      description = "Sprague-Dawley rat strain indicator (1 = Sprague-Dawley, 0 = Wistar-Kyoto lineage incl. SHR)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Wistar-Kyoto / spontaneously hypertensive rat)",
      notes = "The Mlxtran regressor 'Species' multiplies the myocardial-compliance drug effect by (1 - Species) under the comment 'Drug Effect on MC (Wistar)'. The literature studies 11-13 (Table S1) used Sprague-Dawley rats and are the systolic-dysfunction panels of Figure 2a; the in-house study 14A-C used WKY and SHR rats and is the diastolic-dysfunction panel of Figure 2b.",
      source_name = "Species"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "doxorubicin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "doxorubicin", units = "ug", specimen = "plasma", verified = TRUE),
    heart = list(analyte = "doxorubicin", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral1 = list(analyte = "doxorubicin", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "doxorubicin", units = "ug", specimen = "tissue", verified = TRUE),
    auc_heart = list(analyte = "doxorubicin", units = "ug*h/mL", specimen = "tissue", verified = TRUE),
    effect1 = list(
      analyte = "myocardial compliance effect",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect2 = list(
      analyte = "myocardial compliance effect",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect3 = list(
      analyte = "myocardial compliance effect",
      units = "(unitless)",
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
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley, Wistar-Kyoto, spontaneously hypertensive)",
    n_subjects = 114,
    n_studies = 6,
    weight_range = "250-329 g (study means, PD studies; Table S1)",
    disease_state = "Doxorubicin-induced cardiac dysfunction in healthy rats, isoproterenol-induced myocardial hypertrophy (WKY) and spontaneous hypertension (SHR)",
    dose_range = "1.25-3.75 mg/kg intraperitoneal, 4-16 doses (PD studies 11-14); PK from 2-6 mg/kg i.v. or i.p. single doses (studies 1-10)",
    notes = "QSP and PD parameters estimated jointly (Monolix 2018R1 SAEM) on three literature studies (Chang 2015, Kim 2012, Lee 2014; 90 SD rats) and the in-house study 14A-C (8 WKY, 8 hypertrophied WKY, 8 SHR per the Table S1 counts). The PK model was estimated first on ten literature PK studies (Table S1) and fixed. Table S1 and Supplementary Methods."
  )

  ini({
    # Doxorubicin PK (Sang 2021 Table S4; values as coded in the deposited
    # Mlxtran, which fixed the PK while the QSP / PD parameters were
    # estimated). Mlxtran compartment map: Cc (i.p. abdominal depot, V = 1),
    # Ap = central, Ah = heart, A3 = peripheral1, A4 = peripheral2.
    lka <- log(4.532); label("First-order absorption rate from the intraperitoneal depot, ka (1/h)") # Mlxtran 'k=4.532'; Table S4 ka = 4.53 (RSE 16.07%)
    lkel <- log(1.07); label("First-order elimination rate from plasma, ke (1/h)") # Mlxtran 'ke=1.07'; Table S4 ke = 1.07 (RSE 14%)
    lvc <- log(442); label("Apparent plasma volume, V1 (mL per animal)") # Mlxtran 'V1=442'; Table S4 V1 = 442 (RSE 17.2%, unit printed as L/kg -- see vignette)
    lv_heart <- log(60.5); label("Apparent heart volume, V2 (mL per animal)") # Mlxtran 'V2=60.5' (C = Ah/V2); Table S4 V2 = 60.5 (RSE 23.9%)
    lkin_heart <- log(6.75); label("Plasma-to-heart transfer rate, kin_heart (1/h)") # Mlxtran 'K12=6.75'; Table S4 k_in_heart = 6.75 (RSE 26.1%)
    lkout_heart <- log(0.605); label("Heart-to-plasma return rate, kout_heart (1/h)") # Mlxtran 'K21=0.605'; Table S4 k_out_heart = 0.605 (RSE 19.2%)
    lk12 <- log(6.74); label("Plasma-to-peripheral1 rate, k12 (1/h)") # Mlxtran 'K13=6.74'; Table S4 k12 = 6.74 (RSE 95.4%)
    lk21 <- log(13.4); label("Peripheral1-to-plasma rate, k21 (1/h)") # Mlxtran 'K31=13.4'; Table S4 k21 = 13.4 (RSE 68.2%)
    lk13 <- log(3.07); label("Plasma-to-peripheral2 rate, k13 (1/h)") # Mlxtran 'K14=3.07'; Table S4 k13 = 3.07 (RSE 26.1%)
    lk31 <- log(0.0568); label("Peripheral2-to-plasma rate, k31 (1/h)") # Mlxtran 'K41=0.0568'; Table S4 prints k31 = 0.0586 (RSE 19.1%) -- as-run code value kept

    # Cardiovascular system (QSP) parameters fixed from Snelder 2014
    # (Sang 2021 Table S5 'Source (1)' and Table 1 rat column).
    lrbase_hr <- fixed(log(423)); label("Baseline heart rate, HR0 (beats/min)") # Mlxtran 'HR0 = 423'; Table S5 HR0 = 423 (unit printed as beats/h)
    lkout_hr <- fixed(log(11.58)); label("Dissipation rate constant of HR, kout_HR (1/h)") # Mlxtran 'KoutHR =11.58'; Table 1 rat 11.6
    lkout_tpr <- fixed(log(3.58)); label("Dissipation rate constant of TPR, kout_TPR (1/h)") # Mlxtran 'KoutTPR =3.58'; Table 1 rat 3.58
    lkout_edv <- fixed(log(0.126)); label("Dissipation rate constant of LVEDV, kout_LVEDV (1/h)") # Mlxtran 'KoutLVEDV = 0.126'; Table 1 rat 0.126
    lkout_sv <- fixed(log(0.126)); label("Dissipation rate constant of SV, kout_SV (1/h)") # Table S5 k_out_SV = 0.126 (Source (1), no RSE); Table 1 rat 0.126
    lfb0 <- fixed(log(0.0029)); label("MAP feedback constant at the reference MAP, FB_MAP_0 (1/mmHg)") # Mlxtran 'FB_MAP = 0.0029*...'; Table 1 FB_MAP_0 = 2.9e-3
    e_bslmap_fb <- fixed(-1.98); label("Power exponent of baseline MAP on the MAP feedback constant (unitless)") # Eq. 8 and Methods 'a fixed exponent of -1.98' (Snelder 2014)

    # Estimated QSP / PD parameters (Table S5)
    lfb_lvesv <- log(1.43); label("Feedback of LVESV on LVEDV dissipation, FB_LVESV (1/mL)") # Table S5 FB_LVESV = 1.43 (RSE 12.4%)
    lktr <- log(0.021); label("Transit rate constant of the myocardial-compliance effect, kt (1/h)") # Table S5 kt = 0.021 (RSE 10.0%)
    lauc50_ep <- log(1390); label("Heart AUC at half-maximal bioenergy-production impairment, AUC50_EP (ug*h/mL)") # Table S5 AUC50_EP = 1390 (RSE 7.06%; unit printed as mg*h/mL)
    lauc50_mc <- log(1704); label("Heart AUC at half-maximal myocardial-compliance impairment, AUC50_MC (ug*h/mL)") # Table S5 AUC50_MC = 1704 (RSE 14.5%; unit printed as mg*h/mL)
    hill_ep <- fixed(3); label("Hill coefficient of the bioenergy-production effect (unitless)") # Results 'manually tuned ... bioenergy production (h = 3)'; Mlxtran AUC^3
    hill_mc <- fixed(1); label("Hill coefficient of the myocardial-compliance effect (unitless)") # Results 'myocardial compliance (h = 1)'; Mlxtran AUC/(AUC + AUC50MC)

    # Residual error models are named in Supplementary Methods Eqs. 1-4 but
    # their estimates are not reported anywhere in the paper or supplement.
    addSd_LVEF <- fixed(0); label("Additive residual SD on LVEF (percent; not reported)") # Suppl. Methods Eq. 1 constant error, value not reported
    propSd_LVEDV <- fixed(0); label("Proportional residual SD on LVEDV (fraction; not reported)") # Suppl. Methods Eq. 2 proportional error, value not reported
    propSd_LVESV <- fixed(0); label("Proportional residual SD on LVESV (fraction; not reported)") # Suppl. Methods Eq. 3 proportional error, value not reported
    addSd_MAP <- fixed(0); label("Additive residual SD on MAP (mmHg; not reported)") # Suppl. Methods Eq. 4 constant error, value not reported
  })

  model({
    # Individual parameters (no between-animal variability is reported)
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)
    v_heart <- exp(lv_heart)
    kin_heart <- exp(lkin_heart)
    kout_heart <- exp(lkout_heart)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    hr0 <- exp(lrbase_hr)
    kout_hr <- exp(lkout_hr)
    kout_tpr <- exp(lkout_tpr)
    kout_edv <- exp(lkout_edv)
    kout_sv <- exp(lkout_sv)
    fb_lvesv <- exp(lfb_lvesv)
    ktr <- exp(lktr)
    auc50_ep <- exp(lauc50_ep)
    auc50_mc <- exp(lauc50_mc)

    # Baselines (Mlxtran: SVbase = LVEDVbase - LVESVbase,
    # TPRbase = MAPbase/(HR0*SVbase)); 106.596 mmHg is the reference MAP
    # of Eq. 8 as coded.
    sv_bl <- LVEDV_BL - LVESV_BL
    tpr_bl <- MAP_BL / (hr0 * sv_bl)
    fb_map <- exp(lfb0) * (MAP_BL / 106.596)^e_bslmap_fb

    # Zero-order production rates from the pre-dose steady state
    kin_sv <- kout_sv * sv_bl / (1 - fb_map * MAP_BL)
    kin_edv <- kout_edv * (1 - fb_lvesv * LVESV_BL) * LVEDV_BL / (1 + fb_map * MAP_BL)
    kin_hr <- kout_hr * hr0 / (1 - fb_map * MAP_BL)
    kin_tpr <- kout_tpr * tpr_bl / (1 - fb_map * MAP_BL)

    # Doxorubicin PK (per-animal amounts in ug, volumes in mL)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + kout_heart * heart + k21 * peripheral1 + k31 * peripheral2 -
      (kin_heart + k12 + kel + k13) * central
    d/dt(heart) <- kin_heart * central - kout_heart * heart
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    Cheart <- heart / v_heart

    # Cumulative heart AUC drives the PD (Eqs. 9-10; Mlxtran ddt_AUC = C)
    d/dt(auc_heart) <- Cheart

    # Drug effects (Eqs. 9, 10, 13; Mlxtran EdrugEP / EdrugMC)
    e_drug_ep <- auc_heart^hill_ep / (auc_heart^hill_ep + auc50_ep^hill_ep)
    e_drug_mc <- auc_heart^hill_mc / (auc_heart^hill_mc + auc50_mc^hill_mc) * (1 - STRAIN_SD)
    e_ep_sv <- exp(-e_drug_ep)
    e_mc_sv <- 1 - e_drug_mc

    # Myocardial-compliance transit delay onto LVEDV production (Eqs. 10-12)
    d/dt(effect1) <- (e_drug_mc - effect1) * ktr
    d/dt(effect2) <- (effect1 - effect2) * ktr
    d/dt(effect3) <- (effect2 - effect3) * ktr
    e_mc_edv <- 1 - effect3

    # Algebraic haemodynamics (Eqs. 2, 3, 7)
    LVESV <- edv - sv
    MAP <- sv * hr * tpr
    LVEF <- sv / edv * 100
    LVEDV <- edv

    # Numerical guard in the deposited code: all turnover states freeze once
    # SV or LVESV is non-positive or LVEF falls to 20% or below.
    flag <- 0
    if (sv > 0 & LVESV > 0 & LVEF > 20) {
      flag <- 1
    }

    # MAP feedback as coded: clamped at 1 for HR, TPR and LVEDV, unclamped
    # for SV (Mlxtran 'FBMAP' vs 'FB_MAP*MAP').
    fbmap_clamped <- fb_map * MAP
    if (fbmap_clamped > 1) {
      fbmap_clamped <- 1
    }

    # Turnover states (Eqs. 1, 4, 5, 6)
    d/dt(sv) <- (kin_sv * (1 - fb_map * MAP) * e_ep_sv * e_mc_sv - kout_sv * sv) * flag
    d/dt(edv) <- (kin_edv * (1 + fbmap_clamped) * e_mc_edv - kout_edv * (1 - fb_lvesv * LVESV) * edv) * flag
    d/dt(hr) <- (kin_hr * (1 - fbmap_clamped) - kout_hr * hr) * flag
    d/dt(tpr) <- (kin_tpr * (1 - fbmap_clamped) - kout_tpr * tpr) * flag

    sv(0) <- sv_bl
    edv(0) <- LVEDV_BL
    hr(0) <- hr0
    tpr(0) <- tpr_bl

    LVEF ~ add(addSd_LVEF)
    LVEDV ~ prop(propSd_LVEDV)
    LVESV ~ prop(propSd_LVESV)
    MAP ~ add(addSd_MAP)
  })
}
