Yamamoto_2017b_methotrexate_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for methotrexate: an empirical three-compartment",
    "plasma model linked to a brain extracellular-fluid compartment and four",
    "cerebrospinal-fluid compartments (lateral ventricle, third + fourth",
    "ventricle, cisterna magna, subarachnoid space) of fixed rat physiological",
    "volume. Drug enters brain ECF from plasma by a one-way clearance and is",
    "carried through ECF and the CSF compartments in series, and back to",
    "plasma, by a single estimated drug-dispersion flow; methotrexate alone of",
    "the nine drugs also has a one-way efflux clearance from lateral-ventricle",
    "CSF to plasma across the blood-CSF barrier. Probenecid (MRP / OAT / OATP",
    "inhibitor) co-administration raises both the plasma-to-ECF clearance and",
    "the CSF-to-plasma efflux. Plasma, brain ECF and CSF concentrations are",
    "unbound. Naive-pooled fit: no between-animal variability."
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

  covariateData <- list(
    CONMED_PROBENECID = list(
      description = "Probenecid co-administration indicator (1 = probenecid pre-treatment, 0 = methotrexate alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no probenecid",
      notes = paste(
        "Table I: probenecid 150 mg/kg intravenously over 10 min, given before",
        "the methotrexate infusion, as an inhibitor of MRPs, OATs and OATPs.",
        "Acts by Eq. 3, P = P_PAT * (1 + theta_cov * Cov), on the plasma-to-ECF",
        "clearance (theta 4.09) and on the CSF-to-plasma efflux clearance",
        "(theta 0.410). No effect on plasma clearance was identified.",
        "Time-fixed per animal."
      ),
      source_name = "probenecid"
    )
  )

  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "ng", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "methotrexate", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "methotrexate", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "methotrexate", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "methotrexate", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "methotrexate", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 23,
    n_studies = 1,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "40 or 80 mg/kg intravenous infusion over 10 min, with or without probenecid 150 mg/kg",
    regions = "Netherlands (Leiden University)",
    co_medication = "probenecid 150 mg/kg IV over 10 min in part of the animals",
    notes = paste(
      "Table I: 23 rats, 40 and 80 mg/kg over 10 min, plasma sampled 0-300",
      "min and microdialysate (brain ECF, CSF lateral ventricle, CSF cisterna",
      "magna) 0-300 min; data previously published by Westerhout et al.",
      "(Table I reference 25). Unbound plasma fraction 45 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, methotrexate column) ----
    lcl <- log(8.12); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 8.12 (RSE 3.90%)
    lvc <- log(51.2); label("Central volume V_PL (mL)") # Table III VPL 51.2 (RSE 20.3%)
    lq <- log(28.1); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table III QPL_PER1 28.1 (RSE 18.0%)
    lvp <- log(210); label("Peripheral 1 volume V_PER1 (mL)") # Table III VPER1 210 (RSE 7.20%)
    lq2 <- log(1.50); label("Intercompartmental clearance Q_PL_PER2 (mL/min)") # Table III QPL_PER2 1.50 (RSE 14.1%)
    lvp2 <- log(114); label("Peripheral 2 volume V_PER2 (mL)") # Table III VPER2 114 (RSE 6.60%)

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.00109); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table III QPL_ECF 0.00109 (RSE 10.6%)
    lq_diff <- log(0.0598); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0598 (RSE 9.30%)
    lcl_csf <- log(0.105); label("One-way efflux clearance from lateral-ventricle CSF to plasma Q_LV_PL (mL/min)") # Table III QLV_PL 0.105 (RSE 10.6%)

    # Probenecid effects, Eq. 3: P = P_PAT * (1 + theta * CONMED_PROBENECID)
    e_conmed_probenecid_clin <- 4.09; label("Fractional increase in Q_PL_ECF with probenecid (fraction)") # Table III impact on QPL_ECF 4.09 (RSE 3.7%), footnote a probenecid
    e_conmed_probenecid_cl_csf <- 0.410; label("Fractional increase in Q_LV_PL with probenecid (fraction)") # Table III impact on QLV_PL 0.410 (RSE 16.0%), footnote a probenecid

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 0.522; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 0.522 (RSE 7.00%)
    propSd_Cbrain_ecf <- 0.529; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 0.529 (RSE 6.20%)
    propSd_Cbrain_csf_lv <- 0.663; label("Proportional residual SD, CSF lateral ventricle (fraction)") # Table III sigma_CSFLV 0.663 (RSE 5.70%)
    propSd_Cbrain_csf_cm <- 1.00; label("Proportional residual SD, CSF cisterna magna (fraction)") # Table III sigma_CSFCM 1.00 (RSE 9.30%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    clin <- exp(lclin) * (1 + e_conmed_probenecid_clin * CONMED_PROBENECID)
    cl_csf <- exp(lcl_csf) * (1 + e_conmed_probenecid_cl_csf * CONMED_PROBENECID)
    q_diff <- exp(lq_diff)

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 -
      q2 * C_PL + q2 / vp2 * peripheral2 -
      clin * C_PL + q_diff * C_SAS + cl_csf * C_LV
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1
    d/dt(peripheral2) <- q2 * C_PL - q2 / vp2 * peripheral2

    # Brain ECF and CSF compartments in series (Fig. 2), with the
    # lateral-ventricle efflux to plasma (dashed arrow Q_LV_PL)
    d/dt(brain_ecf) <- clin * C_PL - q_diff * C_ECF
    d/dt(brain_csf_lv) <- q_diff * C_ECF - q_diff * C_LV - cl_csf * C_LV
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
    Cbrain_csf_lv ~ prop(propSd_Cbrain_csf_lv)
    Cbrain_csf_cm ~ prop(propSd_Cbrain_csf_cm)
  })
}
