Lee_2025_nafamostat_patient <- function() {
  description <- paste(
    "Two-compartment population PK model with an indirect-response (turnover)",
    "aPTT model for continuously infused nafamostat mesylate in adults on",
    "veno-venous or veno-arterial ECMO, fitted to samples drawn from the",
    "patient's central venous catheter (the 'patient model' of Lee 2025).",
    "Nafamostat is infused into the ECMO drainage line upstream of the pump.",
    "The ECMO sweep-gas flow rate is an exponential covariate on the",
    "peripheral volume. Nafamostat inhibits the first-order loss of aPTT",
    "(Imax model). The sibling model fitted to ECMO-circuit samples from the",
    "same patients is Lee_2025_nafamostat_circuit.",
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
        "Enters the peripheral volume as exp(e_qecmogas_vp * (Q_ECMO_GAS - 3.75)),",
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
      "162 central-venous nafamostat concentrations; 95 central-venous aPTT",
      "values"
    ),
    notes = paste(
      "Demographics from Lee 2025 Table 1. Samples for this model were drawn",
      "from the patient's central venous catheter. Sweep-gas flow rate",
      "median 3.75 L/min (IQR 2.5-5.5)."
    )
  )

  ini({
    # PK: Table 2, patient model
    lcl <- log(189)
    label("Clearance (CL, L/h)") # Table 2 patient model: thetaCL = 189 L/h (RSE 14.9%)
    lvc <- log(7.01)
    label("Central volume of distribution (V1, L)") # Table 2 patient model: thetaV1 = 7.01 L (RSE 42.5%)
    lq <- fixed(log(350))
    label("Intercompartmental clearance (Q, L/h)") # Table 2 patient model: thetaQ = 350 L/h, footnote a (fixed)
    lvp <- log(55.0)
    label("Peripheral volume of distribution (V2, L)") # Table 2 patient model: thetaV2 = 55.0 L (RSE 19.2%)
    e_qecmogas_vp <- -0.852
    label("Exponential effect of ECMO sweep-gas flow rate on V2 (per L/min)") # Table 2 patient model: thetaGas flow rate_V2 = -0.852 (RSE 31.1%); Results equation V2 = thetaV2 * EXP(theta * (Gas flow rate - 3.75))

    # PD: Table 3, patient model; dAPTT/dt = Kin - Kout * (1 - Imax*Cp/(IC50 + Cp)) * aPTT
    limax <- log(0.355)
    label("Maximum fractional inhibition of aPTT loss (Imax, fraction)") # Table 3 patient model: thetaImax = 0.355 (RSE 21.8%)
    lic50 <- fixed(log(350))
    label("Nafamostat concentration for half-maximal inhibition (IC50, ug/L)") # Table 3 patient model: thetaIC50 = 350 ug/L, footnote a (fixed)
    lkin <- log(47.5)
    label("Zero-order aPTT production rate (Kin, s/h)") # Table 3 patient model: thetaKin = 47.5 s/h (RSE 32.5%)
    lrbase <- log(34.5)
    label("Baseline aPTT in the absence of nafamostat (Base, s)") # Table 3 patient model: thetaBase = 34.5 s (RSE 4.47%)

    # IIV: the table percentages are omega x 100 (footnote b gives RSE 'for
    # standard deviation'), so variance = (percent / 100)^2.
    etalcl ~ 0.613^2 # Table 2 patient model: omegaCL = 61.3%
    etalvc ~ fixed(0.695^2) # Table 2 patient model: omegaV1 = 69.5%, footnote a
    etalq ~ fixed(2.43^2) # Table 2 patient model: omegaQ = 243%, footnote a
    etalvp ~ fixed(0.492^2) # Table 2 patient model: omegaV2 = 49.2%, footnote a
    etalimax ~ fixed(0.507^2) # Table 3 patient model: omegaImax = 50.7%, footnote a
    etalic50 ~ fixed(0.706^2) # Table 3 patient model: row 'omegaEC50' (IIV of IC50) = 70.6%, footnote a
    etalkin ~ 1.05^2 # Table 3 patient model: omegaKin = 105%
    etalrbase ~ 0.207^2 # Table 3 patient model: omegaBase = 20.7%

    # Residual error
    propSd <- 0.287
    label("Proportional residual error on nafamostat concentration (fraction)") # Table 2 patient model: sigma proportional = 28.7%
    addSd_aPTT <- 2.83
    label("Additive residual error on aPTT (s)") # Table 3 patient model: sigma additive = 2.83 s
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * exp(e_qecmogas_vp * (Q_ECMO_GAS - 3.75))

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
