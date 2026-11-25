Yamamoto_2017b_acetaminophen_human <- function() {
  description <- paste(
    "Translational semi-physiological multi-compartmental brain distribution",
    "model for acetaminophen in humans: a two-compartment plasma model",
    "estimated from human plasma data, linked to a brain extracellular-fluid",
    "compartment and four cerebrospinal-fluid compartments (lateral",
    "ventricle, third + fourth ventricle, cisterna magna, subarachnoid space)",
    "of fixed human physiological volume. The plasma-to-ECF clearance and the",
    "drug-dispersion flow are the rat estimates scaled allometrically to a",
    "70 kg human, so brain ECF and CSF concentrations are predictions, not",
    "fits. Optional external ventricular drainage removes lateral-ventricle",
    "CSF into a collection state. All concentrations are unbound."
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

  # evd is the cumulative amount of drug drained through the external
  # ventricular drain; the concentration of a collected EVD sample is the
  # increase of evd over the collection interval divided by the volume
  # drained in that interval.
  paper_specific_compartments <- c("evd")

  covariateData <- list(
    CSF_DRAIN_VOL_24H = list(
      description = "Cerebrospinal-fluid drainage through an external ventricular drain, as volume per 24 h",
      units = "mL/24h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "REQUIRED model input; 0 for a patient without an EVD (for example the",
        "Study 2 nerve-root-compression patients). The paper's quantity is the",
        "flow from the lateral ventricle into the EVD, Q_LV_EVD, recorded per",
        "collection interval in mL/h (Table SIV, 1.5-30 mL/h); supply it as",
        "Q_LV_EVD x 24. Enters model() directly as the lateral-ventricle",
        "drainage flow q_evd = CSF_DRAIN_VOL_24H / 24 / 60 mL/min.",
        "Time-varying within subject."
      ),
      source_name = "QLV_EVD"
    )
  )

  compartmentData <- list(
    central = list(analyte = "acetaminophen", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acetaminophen", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "acetaminophen", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    evd = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 7,
    n_studies = 2,
    age_range = "42-79 years (Study 1)",
    weight_range = "56-74 kg (Study 1)",
    sex_female_pct = 71,
    disease_state = paste(
      "Study 1: neurointensive-care patients with an external ventricular",
      "drain (Table II 'traumatic brain injury'; Table SIII lists subarachnoid",
      "haemorrhage as the reason for admission for all seven). Study 2:",
      "patients with nerve-root compression pain (published mean profile)."
    ),
    dose_range = "1 g acetaminophen IV over 30 min (Study 1); 2 g propacetamol IV short infusion (Study 2)",
    regions = "France (Poitiers University Hospital, Study 1)",
    notes = paste(
      "Table II and Table SIII. Study 1: 7 patients (5 female), plasma 0-6 h",
      "and EVD-collected CSF 0-5.5 h; most had received acetaminophen before",
      "inclusion, so the pre-study plasma concentration was used as the",
      "initial plasma amount. Study 2: one mean plasma and subarachnoid CSF",
      "profile from the literature (Table II reference 31). Total plasma",
      "concentrations were converted to unbound with a plasma unbound",
      "fraction of 85 percent."
    )
  )

  ini({
    # ---- Plasma, estimated from human plasma data (Table IV) ----
    lcl <- log(562); label("Plasma clearance CL_PL (mL/min)") # Table IV CLPL 562 (RSE 20.1%)
    lvc <- log(9880); label("Central volume V_PL (mL)") # Table IV VPL 9880 (RSE 41.1%)
    lq <- log(2060); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table IV QPL_PER1 2060 (RSE 31.1%)
    lvp <- log(51900); label("Peripheral volume V_PER1 (mL)") # Table IV VPER1 51900 (RSE 18.3%)

    # ---- Brain drug-specific parameters, allometrically scaled (Table IV, Eq. 9) ----
    # P_human = P_rat * (70 / 0.25)^0.75 = P_rat * 68.4 from the rat
    # acetaminophen estimates (Table III: Q_PL_ECF 0.0281, Q_DIFF 0.0556 mL/min).
    lclin <- fixed(log(1.92)); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table IV QPL_ECF 1.92 FIX
    lq_diff <- fixed(log(3.81)); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table IV QDIFF 3.81 FIX

    # ---- Human physiological volumes, fixed (Table IV footnote a) ----
    v_brain_ecf <- fixed(240); label("Brain ECF volume V_ECF (mL)") # Table IV VECF 240 FIX
    v_brain_csf_lv <- fixed(22.5); label("Lateral-ventricle CSF volume V_LV (mL)") # Table IV VLV 22.5 FIX
    v_brain_csf_tfv <- fixed(22.5); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table IV VTFV 22.5 FIX
    v_brain_csf_cm <- fixed(7.5); label("Cisterna-magna CSF volume V_CM (mL)") # Table IV VCM 7.5 FIX
    v_brain_csf_sas <- fixed(90); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table IV VSAS 90 FIX

    # ---- Between-subject variability (Table IV; exponential, Eq. 6) ----
    # Table IV gives omega as standard deviations; variance = SD^2.
    etalcl ~ 0.2401 # Table IV omega_CLPL 0.490 (RSE 30.2%)
    etalvp ~ 0.055225 # Table IV omega_VPER1 0.235 (RSE 22.5%)

    # ---- Residual error (Table IV; proportional, Eq. 7) ----
    propSd <- 0.250; label("Proportional residual SD, plasma (fraction)") # Table IV sigma_plasma 0.250 (RSE 8.20%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    clin <- exp(lclin)
    q_diff <- exp(lq_diff)
    # EVD flow out of the lateral ventricle, mL/min
    q_evd <- CSF_DRAIN_VOL_24H / 24 / 60

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 -
      clin * C_PL + q_diff * C_SAS
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1

    # Brain ECF and CSF compartments in series (Fig. 2), with drainage of
    # lateral-ventricle CSF into the EVD collection (Fig. 2 grey box)
    d/dt(brain_ecf) <- clin * C_PL - q_diff * C_ECF
    d/dt(brain_csf_lv) <- q_diff * C_ECF - q_diff * C_LV - q_evd * C_LV
    d/dt(brain_csf_tfv) <- q_diff * C_LV - q_diff * C_TFV
    d/dt(brain_csf_cm) <- q_diff * C_TFV - q_diff * C_CM
    d/dt(brain_csf_sas) <- q_diff * C_CM - q_diff * C_SAS
    d/dt(evd) <- q_evd * C_LV

    Cc <- C_PL
    Cbrain_ecf <- C_ECF
    Cbrain_csf_lv <- C_LV
    Cbrain_csf_tfv <- C_TFV
    Cbrain_csf_cm <- C_CM
    Cbrain_csf_sas <- C_SAS

    Cc ~ prop(propSd)
  })
}
