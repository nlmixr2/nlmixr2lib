Yamamoto_2017b_remoxipride_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for remoxipride: an empirical two-compartment",
    "plasma model linked to a brain extracellular-fluid compartment and four",
    "cerebrospinal-fluid compartments (lateral ventricle, third + fourth",
    "ventricle, cisterna magna, subarachnoid space) of fixed rat physiological",
    "volume. Drug enters brain ECF from plasma by a one-way clearance and is",
    "carried through ECF and the CSF compartments in series, and back to",
    "plasma, by a single estimated drug-dispersion flow. Plasma and brain ECF",
    "concentrations are unbound. Naive-pooled fit: no between-animal",
    "variability."
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
    central = list(analyte = "remoxipride", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remoxipride", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "remoxipride", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "remoxipride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "remoxipride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "remoxipride", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "remoxipride", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 29,
    n_studies = 1,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "4, 8 or 16 mg/kg intravenous infusion over 30 min",
    regions = "Netherlands (Leiden University)",
    notes = paste(
      "Table I: 29 rats, 4, 8 and 16 mg/kg over 30 min, plasma sampled 0-240",
      "min and brain-ECF microdialysate 0-240 min; data previously published",
      "by Stevens et al. (Table I reference 26). A separate external",
      "validation set (65 rats, 0.7, 5.2 and 14 mg/kg over 10 min, with CSF",
      "lateral-ventricle and cisterna-magna sampling) was not used for",
      "estimation. Unbound plasma fraction 74 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, remoxipride column) ----
    lcl <- log(114); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 114 (RSE 2.70%)
    lvc <- log(286); label("Central volume V_PL (mL)") # Table III VPL 286 (RSE 16.6%)
    lq <- log(105); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table III QPL_PER1 105 (RSE 10.7%)
    lvp <- log(2310); label("Peripheral volume V_PER1 (mL)") # Table III VPER1 2310 (RSE 6.40%)

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.0141); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table III QPL_ECF 0.0141 (RSE 8.90%)
    lq_diff <- log(0.0176); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0176 (RSE 8.20%)

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 0.348; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 0.348 (RSE 6.80%)
    propSd_Cbrain_ecf <- 0.673; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 0.673 (RSE 6.80%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
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
    Cbrain_ecf ~ prop(propSd_Cbrain_ecf)
  })
}
