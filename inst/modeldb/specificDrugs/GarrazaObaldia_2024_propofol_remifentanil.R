GarrazaObaldia_2024_propofol_remifentanil <- function() {
  description <- "Semi-mechanistic turnover PD model of the haemodynamic effects of propofol and remifentanil during general anaesthesia in adult women undergoing gynaecological surgery (Garraza-Obaldia 2024). Mean arterial pressure (MAP) and heart rate (HR) are two coupled indirect-response states: propofol inhibits MAP production (Emax in the predicted plasma propofol concentration), remifentanil inhibits HR production (sigmoid Emax in the predicted plasma remifentanil concentration), relative HR feeds back on MAP production as a power function, and the surgical stimulus (hysteroscopy) transiently raises MAP production with a monoexponentially dissipating effect. The predicted plasma concentrations from the Schnider (propofol) and Minto (remifentanil) target-controlled-infusion models are supplied as time-varying covariates; the model has no PK states."
  reference <- paste(
    "Garraza-Obaldia M, Jaramillo S, Parra-Guillen ZP, Valencia JF,",
    "Gambus PL, Troconiz IF.",
    "Pharmacodynamic Model of the Hemodynamic Effects of Propofol and",
    "Remifentanil and Their Interaction with Noxious Stimulation.",
    "Pharmaceutics. 2024;16(12):1615.",
    "doi:10.3390/pharmaceutics16121615.",
    sep = " "
  )
  vignette <- "GarrazaObaldia_2024_propofol_remifentanil"
  units <- list(
    time = "min",
    dosing = "none (PD-only model; drug exposure enters as the CP_PROPOFOL_UGML and CP_REMIFENTANIL_NGML covariates)",
    concentration = "mmHg (map), beats/min (hr)"
  )

  covariateData <- list(
    CP_PROPOFOL_UGML = list(
      description = "Predicted plasma propofol concentration at the record time",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying; 0 when no propofol is on board. Garraza-Obaldia 2024 used the plasma concentration predicted by the Schnider propofol PK model that drove the target-controlled infusion (Section 2.1); concentrations were not measured. Effect-site propofol was tested as the driver and rejected (Section 3.1).",
      source_name = "CPropo"
    ),
    CP_REMIFENTANIL_NGML = list(
      description = "Predicted plasma remifentanil concentration at the record time",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying; 0 when no remifentanil is on board. Garraza-Obaldia 2024 used the plasma concentration predicted by the Minto remifentanil PK model that drove the target-controlled infusion (Section 2.1); concentrations were not measured. Effect-site remifentanil was tested as the driver and rejected (Section 3.1).",
      source_name = "CRemi"
    ),
    T_SURG_START = list(
      description = "Time of the start of the surgical procedure (onset of the surgical noxious stimulus) relative to model time zero",
      units = "min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject; tStim of Garraza-Obaldia 2024 Equation (7). The surgical-stimulus multiplier is 1 before this time and 1 + emax_stim * exp(-kstim * (t - T_SURG_START)) from it onwards, with no end-of-surgery term. Set to a time beyond the simulation horizon for a subject who has no surgery. In the founding cohort the stimulus is hysteroscopy.",
      source_name = "tStim"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained (Section 3.1: none of the covariates analysed was significant, p > 0.05)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained (Section 3.1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling; not retained (Section 3.1)."
    ),
    T_SURG = list(
      description = "Duration of the surgical procedure",
      units = "h",
      type = "continuous",
      notes = "Screened by stepwise covariate modelling, including against the stimulus dissipation rate; not retained (Section 3.1)."
    )
  )

  compartmentData <- list(
    map = list(analyte = "mean arterial pressure", units = "mmHg", specimen = "not applicable", verified = TRUE),
    hr = list(analyte = "heart rate", units = "beats/min", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 56L,
    n_studies = 1L,
    age_range = "mean 42 (SD 11) years in the 56-patient training group; mean 43 (SD 12) years over all 77 patients",
    weight_range = "mean 62 (SD 9) kg in the training group; mean 62 (SD 10) kg over all 77 patients",
    sex_female_pct = 100,
    race_ethnicity = NA_character_,
    disease_state = "ASA I-II adult women undergoing outpatient gynaecological surgery under propofol and remifentanil target-controlled-infusion general anaesthesia; patients with chronic hypertension excluded",
    dose_range = "Propofol (Schnider model, effect-site targets about 0.8-5.1 ug/mL) and remifentanil (Minto model, effect-site targets about 0.5-4 ng/mL) given alone or together in four step-up / step-down groups (Table S1)",
    regions = "Spain (Hospital Clinic de Barcelona)",
    notes = "Final model fitted to the 56-patient training group, who underwent hysteroscopy; the 21-patient test group (other procedures) was used for external validation (Section 2.2, Results, Table 1). MAP measured non-invasively every 1 min; HR from continuous ECG."
  )

  ini({
    # ---- Propofol effect on MAP production (Equation 3; gamma removed) ----
    lec50_propofol <- log(8.79); label("Plasma propofol concentration giving half of Emax on MAP production, C50Propo (ug/mL)") # Table 2, C50Propo = 8.79 ug/mL
    logitemax_propofol <- log(0.88 / (1 - 0.88)); label("Maximum fractional inhibition of MAP production by propofol, EmaxPropo (logit scale)") # Table 2, EmaxPropo = 0.88; logistic transform per Section 2.4

    # ---- Remifentanil effect on HR production (Equation 4) ----
    lec50_remifentanil <- log(4.57); label("Plasma remifentanil concentration giving half of Emax on HR production, C50Remi (ng/mL)") # Table 2, C50Remi = 4.57 ng/mL
    logitemax_remifentanil <- log(0.69 / (1 - 0.69)); label("Maximum fractional inhibition of HR production by remifentanil, EmaxRemi,HR (logit scale)") # Table 2, EmaxRemi,HR = 0.69
    lhill_remifentanil <- log(1.20); label("Hill coefficient of the remifentanil effect on HR production, gammaRemi (unitless)") # Table 2, gammaRemi = 1.20

    # ---- Cardiovascular system ----
    lrbase_map <- log(88.2); label("Baseline mean arterial pressure before drug administration, MAP0 (mmHg)") # Table 2, MAP0 = 88.2 mmHg
    lrbase_hr <- log(70.8); label("Baseline heart rate before drug administration, HR0 (beats/min)") # Table 2, HR0 = 70.8 bpm
    lkout <- log(0.123); label("First-order dissipation rate constant shared by MAP and HR, kout (1/min)") # Table 2, kout = 0.123 1/min (one typical value for both responses, Section 3.1)
    lfb_hr <- log(0.806); label("Exponent of the relative-HR feedback on MAP production, FBHR (unitless)") # Table 2, FBHR = 0.806

    # ---- Surgical stimulus (hysteroscopy) on MAP production (Equation 7) ----
    lkstim <- log(0.021); label("First-order dissipation rate constant of the surgical stimulus effect, kStim (1/min)") # Table 2, kStim = 0.021 1/min
    lemax_stim <- log(0.284); label("Fractional increase in MAP production at the onset of the surgical stimulus, thetaStim (unitless)") # Table 2, thetaStim = 0.284

    # ---- IIV (log-normal; omega^2 = log(1 + CV^2), Table 2 footnote CV = sqrt(exp(omega^2) - 1)) ----
    # Block on C50Propo and kout,MAP: Table 2 'Covariance between random effects'
    # C50Propo, kout,MAP = 0.76 (a covariance: its bootstrap 95% CI reaches 1.43).
    # Variances: C50Propo CV 104.0% -> 0.733137; kout,MAP CV 142.6% -> 1.109710.
    etalec50_propofol + etalkout_map ~ c(0.733137, 0.76, 1.109710)
    etalec50_remifentanil ~ 0.730139 # Table 2, IIV C50Remi 103.7% CV
    etalrbase_map ~ 0.00742009 # Table 2, IIV MAP0 8.63% CV
    etalrbase_hr ~ 0.0205239 # Table 2, IIV HR0 14.4% CV
    etalkout_hr ~ 0.991474 # Table 2, IIV kout,HR 130.2% CV
    etalkstim ~ fixed(0.693147) # Section 3.1, variance held at the equivalent of a 100% CV: log(1 + 1^2) = log(2)

    # ---- Residual error (additive, Section 2.4) ----
    addSd_map <- 5.53; label("Additive residual error on MAP (mmHg)") # Table 2, residual error in MAP = 5.53 mmHg
    addSd_hr <- 4.30; label("Additive residual error on HR (beats/min)") # Table 2, residual error in HR = 4.30 bpm
  })

  model({
    # Individual parameters
    ec50_propofol <- exp(lec50_propofol + etalec50_propofol)
    emax_propofol <- expit(logitemax_propofol)
    ec50_remifentanil <- exp(lec50_remifentanil + etalec50_remifentanil)
    emax_remifentanil <- expit(logitemax_remifentanil)
    hill_remifentanil <- exp(lhill_remifentanil)
    rbase_map <- exp(lrbase_map + etalrbase_map)
    rbase_hr <- exp(lrbase_hr + etalrbase_hr)
    kout_map <- exp(lkout + etalkout_map)
    kout_hr <- exp(lkout + etalkout_hr)
    fb_hr <- exp(lfb_hr)
    kstim <- exp(lkstim + etalkstim)
    emax_stim <- exp(lemax_stim)

    # Production rates from the drug-free steady state (Section 2.4)
    kin_map <- kout_map * rbase_map
    kin_hr <- kout_hr * rbase_hr

    # Drug effects on production (Equations 3 and 4)
    f_propofol <- emax_propofol * CP_PROPOFOL_UGML / (CP_PROPOFOL_UGML + ec50_propofol)
    g_remifentanil <- emax_remifentanil * CP_REMIFENTANIL_NGML^hill_remifentanil /
      (CP_REMIFENTANIL_NGML^hill_remifentanil + ec50_remifentanil^hill_remifentanil)

    # Surgical stimulus (Equation 7): 1 before T_SURG_START. The elapsed time is
    # zeroed before the stimulus so exp() cannot overflow for a far-future start.
    stim_on <- t >= T_SURG_START
    t_stim <- (t - T_SURG_START) * stim_on
    h_stim <- 1 + emax_stim * exp(-kstim * t_stim) * stim_on

    # Turnover ODEs (Equations 1, 2 and 5; final model per Figure 2B)
    map(0) <- rbase_map
    hr(0) <- rbase_hr
    d/dt(map) <- kin_map * (1 - f_propofol) * (hr / rbase_hr)^fb_hr * h_stim - kout_map * map
    d/dt(hr) <- kin_hr * (1 - g_remifentanil) - kout_hr * hr

    map ~ add(addSd_map)
    hr ~ add(addSd_hr)
  })
}
