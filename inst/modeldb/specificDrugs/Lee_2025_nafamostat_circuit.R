Lee_2025_nafamostat_circuit <- function() {
  description <- paste(
    "Two-compartment population PK model with an indirect-response (turnover)",
    "aPTT model for continuously infused nafamostat mesylate in adults on",
    "veno-venous or veno-arterial ECMO, fitted to samples drawn from the",
    "ECMO circuit at the oxygenator outflow (the 'ECMO model' of Lee 2025).",
    "Nafamostat is infused into the ECMO drainage line upstream of the pump,",
    "so circuit concentrations exceed systemic ones. The ECMO sweep-gas flow",
    "rate is an exponential covariate on clearance. Nafamostat inhibits the",
    "first-order loss of aPTT (Imax model). The sibling model fitted to",
    "central-venous samples from the same patients is",
    "Lee_2025_nafamostat_patient.",
    sep = " "
  )
  reference <- paste(
    "Lee DH, Lee JH, Jang JH, Kim YK, Kang G, Jung SY, Her M, Jang HJ.",
    "(2025).",
    "Pharmacokinetic and pharmacodynamic analyses of nafamostat in ECMO",
    "patients: comparing central vein and ECMO machine samples.",
    "Front Pharmacol 16:1541131. doi:10.3389/fphar.2025.1541131.",
    sep = " "
  )
  vignette <- "Lee_2025_nafamostat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    Q_ECMO_GAS = list(
      description = "ECMO sweep-gas (oxygenator gas) flow rate",
      units = "L/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters clearance as exp(e_qecmogas_cl * (Q_ECMO_GAS - 3.75)),",
        "centred on the cohort median of 3.75 L/min (Table 1, 'Gas flow rate,",
        "L/min' median 3.75, IQR 2.5-5.5). The paper does not say whether the",
        "covariate was time-varying within a subject; treat it as the value",
        "prevailing during the sampling window."
      ),
      source_name = "Gas flow rate"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "nafamostat",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    aPTT = list(
      analyte = "aPTT",
      units = "s",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "adults older than 18 years; mean 60.3 (SD 12.1) years",
    age_median = "61 years (IQR 57.8-67.5)",
    weight_range = "mean 70.3 (SD 15.7) kg",
    weight_median = "68.2 kg (IQR 62.9-76.8)",
    sex_female_pct = 29.2,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Critically ill adults on VA-ECMO (n = 13) or VV-ECMO (n = 11) for",
      "respiratory and/or cardiac failure (pneumonia 10, cardiogenic shock",
      "or ventricular fibrillation 4, trauma 4, gastrointestinal infection 3,",
      "interstitial lung disease 1, pulmonary thromboembolism 1, aortic",
      "dissection 1); 7 also on CRRT"
    ),
    dose_range = paste(
      "Continuous IV infusion of nafamostat mesylate into the ECMO drainage",
      "line upstream of the pump, started at 15 mg/h without a bolus and",
      "titrated (15-45 mg/h observed) to an aPTT of 40-80 s"
    ),
    regions = "Republic of Korea (Haeundae Paik Hospital, Busan)",
    n_observations = paste(
      "162 ECMO-circuit nafamostat concentrations; 95 ECMO-circuit aPTT",
      "values"
    ),
    notes = paste(
      "Demographics from Lee 2025 Table 1. Samples for this model were drawn",
      "from the ECMO circuit line returning oxygenated blood to the patient.",
      "Sweep-gas flow rate",
      "median 3.75 L/min (IQR 2.5-5.5)."
    )
  )

  ini({
    # PK: Table 2, ECMO model
    lcl <- log(85.2)
    label("Clearance (CL, L/h)") # Table 2 ECMO model: thetaCL = 85.2 L/h (RSE 7.82%)
    lvc <- log(3.83)
    label("Central volume of distribution (V1, L)") # Table 2 ECMO model: thetaV1 = 3.83 L (RSE 17.9%)
    lq <- log(46.7)
    label("Intercompartmental clearance (Q, L/h)") # Table 2 ECMO model: thetaQ = 46.7 L/h (RSE 46.9%)
    lvp <- log(36.8)
    label("Peripheral volume of distribution (V2, L)") # Table 2 ECMO model: thetaV2 = 36.8 L (RSE 36.4%)
    e_qecmogas_cl <- 0.0999
    label("Exponential effect of ECMO sweep-gas flow rate on CL (per L/min)") # Table 2 ECMO model: thetaGas flow rate_CL = 0.0999 (RSE 16%); Results equation CL = thetaCL * EXP(theta * (Gas flow rate - 3.75))

    # PD: Table 3, ECMO model; dAPTT/dt = Kin - Kout * (1 - Imax*Cp/(IC50 + Cp)) * aPTT
    limax <- log(0.436)
    label("Maximum fractional inhibition of aPTT loss (Imax, fraction)") # Table 3 ECMO model: thetaImax = 0.436 (RSE 9.25%)
    lic50 <- fixed(log(581))
    label("Nafamostat concentration for half-maximal inhibition (IC50, ug/L)") # Table 3 ECMO model: thetaIC50 = 581 ug/L, footnote a (fixed)
    lkin <- log(50.6)
    label("Zero-order aPTT production rate (Kin, s/h)") # Table 3 ECMO model: thetaKin = 50.6 s/h (RSE 26.9%)
    lrbase <- log(33.0)
    label("Baseline aPTT in the absence of nafamostat (Base, s)") # Table 3 ECMO model: thetaBase = 33.0 s (RSE 4.43%)

    # IIV: the table percentages are omega x 100 (footnote b gives RSE 'for
    # standard deviation'), so variance = (percent / 100)^2.
    etalcl ~ 0.295^2 # Table 2 ECMO model: omegaCL = 29.5%
    etalvc ~ fixed(0.628^2) # Table 2 ECMO model: omegaV1 = 62.8%, footnote a
    etalq ~ fixed(1.34^2) # Table 2 ECMO model: omegaQ = 134%, footnote a
    etalvp ~ fixed(0) # Table 2 ECMO model: omegaV2 = 0.000%, footnote a
    etalimax ~ fixed(0.273^2) # Table 3 ECMO model: omegaImax = 27.3%, footnote a
    etalic50 ~ fixed(0) # Table 3 ECMO model: row 'omegaEC50' (IIV of IC50) = 0.000%, footnote a
    etalkin ~ 0.758^2 # Table 3 ECMO model: omegaKin = 75.8% (no RSE, bootstrap CI or footnote a printed)
    etalrbase ~ 0.212^2 # Table 3 ECMO model: omegaBase = 21.2%

    # Residual error
    propSd <- 0.309
    label("Proportional residual error on nafamostat concentration (fraction)") # Table 2 ECMO model: sigma proportional = 30.9%
    addSd_aPTT <- 3.43
    label("Additive residual error on aPTT (s)") # Table 3 ECMO model: sigma additive = 3.43 s
  })
  model({
    cl <- exp(lcl + etalcl) * exp(e_qecmogas_cl * (Q_ECMO_GAS - 3.75))
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    imax <- exp(limax + etalimax)
    ic50 <- exp(lic50 + etalic50)
    kin <- exp(lkin + etalkin)
    rbase <- exp(lrbase + etalrbase)
    kout <- kin / rbase

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L -> ug/L
    Cc <- 1000 * central / vc

    aPTT(0) <- rbase
    d/dt(aPTT) <- kin - kout * (1 - imax * Cc / (ic50 + Cc)) * aPTT

    Cc ~ prop(propSd)
    aPTT ~ add(addSd_aPTT)
  })
}
