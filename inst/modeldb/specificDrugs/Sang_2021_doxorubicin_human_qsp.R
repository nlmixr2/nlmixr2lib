Sang_2021_doxorubicin_human_qsp <- function() {
  description <- paste(
    "QSP (human, translated from rat). Human QSP-PD model of",
    "doxorubicin-induced systolic dysfunction (Sang 2021), obtained by",
    "scaling the rat QSP-PK-PD model (Sang_2021_doxorubicin_rat_qsp) to",
    "humans. Four indirect-response turnover states describe stroke volume",
    "(SV), left ventricular end-diastolic volume (LVEDV), heart rate (HR)",
    "and total peripheral resistance (TPR), coupled by MAP = HR * TPR * SV",
    "feedback and by the LVESV feedback on LVEDV dissipation. Dissipation",
    "rate constants are allometrically scaled from rat (250 g) to human",
    "(70 kg) with exponent -0.25; FB_LVESV is scaled by the ratio of",
    "baseline LVESV; FB_MAP_0 and AUC50_EP are kept at their rat values.",
    "The cumulative heart-tissue doxorubicin AUC impairs bioenergy",
    "production (sigmoid, Hill 3, entering SV production as exp(-E)); the",
    "rat myocardial-compliance arm is not part of the human model. The",
    "paper predicted human heart exposure with the separately published He",
    "2018 multiscale doxorubicin PBPK model, which is not in this file: the",
    "heart-tissue doxorubicin concentration must be supplied as the",
    "time-varying covariate CEFFECT (ug/mL). Baselines (LVEDV, LVESV, HR,",
    "MAP) are per-patient covariates; LVEF is reported in percent.",
    sep = " "
  )
  reference <- paste(
    "Sang L, Yuan Y, Zhou Y, Zhou Z, Jiang M, Liu X, Hao K, He H (2021).",
    "A quantitative systems pharmacology approach to predict the",
    "safe-equivalent dose of doxorubicin in patients with cardiovascular",
    "comorbidity. CPT Pharmacometrics Syst Pharmacol. 10(12):1512-1524.",
    "doi:10.1002/psp4.12719. Human parameters from Table 1; equations",
    "from the deposited rat Mlxtran code. Heart exposure in the source",
    "came from He H et al. (2018) Pharm Res 35:174",
    "(doi:10.1007/s11095-018-2456-8), not included here.",
    sep = " "
  )
  vignette <- "Sang_2021_doxorubicin_cardiotoxicity"

  paper_specific_compartments <- c("sv", "edv", "tpr", "auc_heart")

  units <- list(
    time = "h",
    dosing = "none (no dosing compartment; doxorubicin exposure enters as the heart-tissue concentration covariate CEFFECT)",
    concentration = "ug/mL (heart-tissue doxorubicin supplied as CEFFECT; the cumulative heart AUC driving the PD is in ug*h/mL)"
  )

  covariateData <- list(
    CEFFECT = list(
      description = "Time-varying doxorubicin concentration in heart tissue",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Heart-tissue doxorubicin concentration (ug/mL, equivalent to ug/g tissue) predicted by an external PK model. Sang 2021 used the He 2018 multiscale PBPK model (Pharm Res 35:174). Integrated inside the model to the cumulative heart AUC that drives the bioenergy-production effect; set to 0 before the first dose. Supply it on a dense time grid (e.g. hourly to daily): rxode2 does not restart the integrator where a covariate changes, so a profile given only at its change points can be stepped over and its AUC lost.",
      source_name = "C (heart concentration in the rat Mlxtran code, 'C =Ah/V2')"
    ),
    LVEDV_BL = list(
      description = "Baseline (pre-treatment) left ventricular end-diastolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Table 1 typical human values 113 mL (cardiovascular healthy) and 141 mL (diseased), CV 30%.",
      source_name = "LVEDV0"
    ),
    LVESV_BL = list(
      description = "Baseline (pre-treatment) left ventricular end-systolic volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Baseline SV = LVEDV_BL - LVESV_BL. Table 1 typical healthy value is LVEDV0 - SV0 = 113 - 65 = 48 mL (the value behind the Table 1 FB_LVESV scaling).",
      source_name = "LVEDV0 - SV0"
    ),
    HR_BL = list(
      description = "Baseline (pre-treatment) heart rate",
      units = "beats/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Table 1 human value 70 (CV 30%; unit printed as beats/h). Initial condition of the HR state.",
      source_name = "HR0"
    ),
    MAP_BL = list(
      description = "Baseline (pre-treatment) mean arterial pressure",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. Baseline TPR = MAP_BL / (HR_BL * SV_BL). Table 1 gives TPR0 0.02 (healthy) / 0.025 (diseased) mmHg*min/mL; the healthy typical MAP is 70 * 0.02 * 65 = 91 mmHg, as stated in the Methods.",
      source_name = "MAPbase (derived from HR0 * TPR0 * SV0)"
    )
  )

  compartmentData <- list(
    auc_heart = list(analyte = "doxorubicin", units = "ug*h/mL", specimen = "tissue", verified = TRUE),
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
    species = "human",
    n_subjects = 13994,
    disease_state = "Virtual adult cancer patients treated with doxorubicin, with (n = 5395) or without (n = 8599) cardiovascular comorbidity (LVEDV index >= 81.5 mL/m^2 and/or MAP >= 115 mmHg)",
    weight_range = "70 kg typical (CV 30%)",
    dose_range = "Cumulative doxorubicin 120-900 mg/m^2, infused every 3 weeks over 1 year",
    notes = "Virtual population only; no human data were fitted. Physiological baselines drawn by Monte Carlo with 30% variation around the Table 1 typical values (Drafts 2013 and Snelder 2014). The model was checked against the reported incidence of doxorubicin-induced cardiac dysfunction of Von Hoff 1979."
  )

  ini({
    # Rat-to-human dissipation rate constants, allometrically scaled from
    # 250 g to 70 kg with exponent -0.25 (Eq. 14): (0.25/70)^0.25 = 0.2445.
    lkout_sv <- fixed(log(0.0308)); label("Dissipation rate constant of SV, kout_SV (1/h)") # Table 1 human 0.0308 (= 0.126 * 0.2445)
    lkout_edv <- fixed(log(0.0308)); label("Dissipation rate constant of LVEDV, kout_LVEDV (1/h)") # Table 1 human 0.0308 (= 0.126 * 0.2445)
    lkout_hr <- fixed(log(2.83)); label("Dissipation rate constant of HR, kout_HR (1/h)") # Table 1 human 2.83 (= 11.58 * 0.2445)
    lkout_tpr <- fixed(log(0.875)); label("Dissipation rate constant of TPR, kout_TPR (1/h)") # Table 1 human 0.875 (= 3.58 * 0.2445)

    # Parameters kept constant across species
    lfb0 <- fixed(log(0.0029)); label("MAP feedback constant at the reference MAP, FB_MAP_0 (1/mmHg)") # Table 1 FB_MAP_0 human 2.9e-3 ('Constant across species')
    e_bslmap_fb <- fixed(-1.98); label("Power exponent of baseline MAP on the MAP feedback constant (unitless)") # Eq. 8, fixed exponent -1.98 (Snelder 2014)
    lauc50_ep <- fixed(log(1390)); label("Heart AUC at half-maximal bioenergy-production impairment, AUC50_EP (ug*h/mL)") # Table 1 AUC50_EP human 1390 ('Constant across species'; unit printed as mg*h/mL)
    hill_ep <- fixed(3); label("Hill coefficient of the bioenergy-production effect (unitless)") # Results 'bioenergy production (h = 3)'

    # FB_LVESV scaled by the baseline-LVESV ratio (Eq. 15): 1.43 * 0.085 / 48
    lfb_lvesv <- fixed(log(2.532e-3)); label("Feedback of LVESV on LVEDV dissipation, FB_LVESV (1/mL)") # Table 1 FB_LVESV human 2.532e-3

    # Between-patient variability used in the virtual trials (Table 1 CV
    # column); log-normal variance = log(1 + CV^2).
    etalfb_lvesv ~ fixed(0.015259) # Table 1 FB_LVESV CV 12.4% -> log(1 + 0.124^2)
    etalauc50_ep ~ fixed(0.004972) # Table 1 AUC50_EP CV 7.06% -> log(1 + 0.0706^2)

    # No residual error was used or reported for the human simulations.
    addSd_LVEF <- fixed(0); label("Additive residual SD on LVEF (percent; not reported)") # simulation-only model, no residual error reported
    propSd_LVEDV <- fixed(0); label("Proportional residual SD on LVEDV (fraction; not reported)") # simulation-only model, no residual error reported
    propSd_LVESV <- fixed(0); label("Proportional residual SD on LVESV (fraction; not reported)") # simulation-only model, no residual error reported
    addSd_MAP <- fixed(0); label("Additive residual SD on MAP (mmHg; not reported)") # simulation-only model, no residual error reported
  })

  model({
    kout_sv <- exp(lkout_sv)
    kout_edv <- exp(lkout_edv)
    kout_hr <- exp(lkout_hr)
    kout_tpr <- exp(lkout_tpr)
    fb_lvesv <- exp(lfb_lvesv + etalfb_lvesv)
    auc50_ep <- exp(lauc50_ep + etalauc50_ep)

    # Baselines; 106.596 mmHg is the Eq. 8 reference MAP of the rat code,
    # which Table 1 does not list among the rescaled parameters.
    sv_bl <- LVEDV_BL - LVESV_BL
    tpr_bl <- MAP_BL / (HR_BL * sv_bl)
    fb_map <- exp(lfb0) * (MAP_BL / 106.596)^e_bslmap_fb

    kin_sv <- kout_sv * sv_bl / (1 - fb_map * MAP_BL)
    kin_edv <- kout_edv * (1 - fb_lvesv * LVESV_BL) * LVEDV_BL / (1 + fb_map * MAP_BL)
    kin_hr <- kout_hr * HR_BL / (1 - fb_map * MAP_BL)
    kin_tpr <- kout_tpr * tpr_bl / (1 - fb_map * MAP_BL)

    # Cumulative heart AUC from the supplied heart concentration
    d/dt(auc_heart) <- CEFFECT

    e_drug_ep <- auc_heart^hill_ep / (auc_heart^hill_ep + auc50_ep^hill_ep)
    e_ep_sv <- exp(-e_drug_ep)

    LVESV <- edv - sv
    MAP <- sv * hr * tpr
    LVEF <- sv / edv * 100
    LVEDV <- edv

    # Numerical guard carried over from the deposited rat code
    flag <- 0
    if (sv > 0 & LVESV > 0 & LVEF > 20) {
      flag <- 1
    }
    fbmap_clamped <- fb_map * MAP
    if (fbmap_clamped > 1) {
      fbmap_clamped <- 1
    }

    d/dt(sv) <- (kin_sv * (1 - fb_map * MAP) * e_ep_sv - kout_sv * sv) * flag
    d/dt(edv) <- (kin_edv * (1 + fbmap_clamped) - kout_edv * (1 - fb_lvesv * LVESV) * edv) * flag
    d/dt(hr) <- (kin_hr * (1 - fbmap_clamped) - kout_hr * hr) * flag
    d/dt(tpr) <- (kin_tpr * (1 - fbmap_clamped) - kout_tpr * tpr) * flag

    sv(0) <- sv_bl
    edv(0) <- LVEDV_BL
    hr(0) <- HR_BL
    tpr(0) <- tpr_bl

    LVEF ~ add(addSd_LVEF)
    LVEDV ~ prop(propSd_LVEDV)
    LVESV ~ prop(propSd_LVESV)
    MAP ~ add(addSd_MAP)
  })
}
