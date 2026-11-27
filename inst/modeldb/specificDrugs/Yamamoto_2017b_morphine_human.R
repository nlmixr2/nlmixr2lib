Yamamoto_2017b_morphine_human <- function() {
  description <- paste(
    "Translational semi-physiological multi-compartmental brain distribution",
    "model for morphine in humans: a two-compartment plasma model estimated",
    "from human plasma data, linked to a brain extracellular-fluid",
    "compartment and four cerebrospinal-fluid compartments (lateral",
    "ventricle, third + fourth ventricle, cisterna magna, subarachnoid space)",
    "of fixed human physiological volume. The plasma-to-ECF clearance and the",
    "drug-dispersion flow are the rat estimates scaled allometrically to a",
    "70 kg human, so brain ECF and CSF concentrations are predictions, not",
    "fits. All concentrations are unbound."
  )
  reference <- paste(
    "Yamamoto Y, Valitalo PA, van den Berg DJ, Hartman R, van den Brink W,",
    "Wong YC, Huntjens DR, Proost JH, Vermeulen A, Krauwinkel W, Bakshi S,",
    "Aranzana-Climent V, Marchand S, Dahyot-Fizelier C, Couet W, Danhof M,",
    "van Hasselt JGC, de Lange ECM. A Generic Multi-Compartmental CNS",
    "Distribution Model Structure for 9 Drugs Allows Prediction of Human Brain",
    "Target Site Concentrations. Pharm Res. 2017;34(2):333-351.",
    "doi:10.1007/s11095-016-2065-3."
  )
  vignette <- "Yamamoto_2017b_cns_distribution"

  units <- list(time = "min", dosing = "ng", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "morphine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "morphine", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 3,
    n_studies = 2,
    disease_state = "traumatic brain injury (neurointensive care, brain microdialysis)",
    dose_range = "10 mg morphine IV over 10 min",
    notes = paste(
      "Table II: two published studies of morphine in traumatic-brain-injury",
      "patients with brain microdialysis -- 2 patients (Table II reference 34)",
      "and 1 patient (reference 33) -- with unbound plasma 0-3 h and brain ECF",
      "0-3 h from both the physiologically 'normal' and the 'injured' side of",
      "the brain. Unbound plasma concentrations were reported by the original",
      "studies."
    )
  )

  ini({
    # ---- Plasma, estimated from human plasma data (Table IV) ----
    lcl <- log(3070); label("Plasma clearance CL_PL (mL/min)") # Table IV CLPL 3070 (RSE 15.8%)
    lvc <- log(16000); label("Central volume V_PL (mL)") # Table IV VPL 16000 (RSE 35.3%)
    lq <- log(3030); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table IV QPL_PER1 3030 (RSE 0.60%)
    lvp <- log(95400); label("Peripheral volume V_PER1 (mL)") # Table IV VPER1 95400 (RSE 2.50%)

    # ---- Brain drug-specific parameters, allometrically scaled (Table IV, Eq. 9) ----
    # P_human = P_rat * (70 / 0.25)^0.75 = P_rat * 68.4. Q_PL_ECF scales the
    # rat 4 mg/kg value (Table III 0.00750 mL/min); Q_DIFF scales 0.0200 mL/min.
    lclin <- fixed(log(0.513)); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table IV QPL_ECF 0.513 FIX
    lq_diff <- fixed(log(1.37)); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table IV QDIFF 1.37 FIX

    # ---- Human physiological volumes, fixed (Table IV footnote a) ----
    v_brain_ecf <- fixed(240); label("Brain ECF volume V_ECF (mL)") # Table IV VECF 240 FIX
    v_brain_csf_lv <- fixed(22.5); label("Lateral-ventricle CSF volume V_LV (mL)") # Table IV VLV 22.5 FIX
    v_brain_csf_tfv <- fixed(22.5); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table IV VTFV 22.5 FIX
    v_brain_csf_cm <- fixed(7.5); label("Cisterna-magna CSF volume V_CM (mL)") # Table IV VCM 7.5 FIX
    v_brain_csf_sas <- fixed(90); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table IV VSAS 90 FIX

    # ---- Between-subject variability (Table IV; exponential, Eq. 6) ----
    # Table IV gives omega as standard deviations; variance = SD^2.
    etalcl ~ 0.073441 # Table IV omega_CLPL 0.271 (RSE 19.9%)
    etalvc ~ 0.355216 # Table IV omega_VPL 0.596 (RSE 20.0%)

    # ---- Residual error (Table IV; proportional, Eq. 7) ----
    propSd <- 0.0960; label("Proportional residual SD, plasma (fraction)") # Table IV sigma_plasma 0.0960 (RSE 22.9%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    clin <- exp(lclin)
    q_diff <- exp(lq_diff)

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 -
      clin * C_PL + q_diff * C_SAS
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1

    # Brain ECF and CSF compartments in series (Fig. 2)
    d/dt(brain_ecf) <- clin * C_PL - q_diff * C_ECF
    d/dt(brain_csf_lv) <- q_diff * C_ECF - q_diff * C_LV
    d/dt(brain_csf_tfv) <- q_diff * C_LV - q_diff * C_TFV
    d/dt(brain_csf_cm) <- q_diff * C_TFV - q_diff * C_CM
    d/dt(brain_csf_sas) <- q_diff * C_CM - q_diff * C_SAS

    Cc <- C_PL
    Cbrain_ecf <- C_ECF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_tfv <- C_TFV
    Cbrain_csf_cm <- C_CM
    Cbrain_csf_sas <- C_SAS

    Cc ~ prop(propSd)
  })
}
