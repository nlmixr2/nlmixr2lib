Yamamoto_2017b_quinidine_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for quinidine: an empirical two-compartment",
    "plasma model linked to a brain extracellular-fluid compartment, a brain",
    "intracellular-fluid compartment, and four cerebrospinal-fluid",
    "compartments (lateral ventricle, third + fourth ventricle, cisterna",
    "magna, subarachnoid space) of fixed rat physiological volume. Drug enters",
    "brain ECF from plasma by a one-way clearance, exchanges with brain ICF,",
    "and is carried through ECF and the CSF compartments in series, and back",
    "to plasma, by a single estimated drug-dispersion flow. P-glycoprotein",
    "inhibition with tariquidar raises the plasma-to-ECF clearance. Plasma,",
    "brain and CSF concentrations are unbound. Naive-pooled fit: no",
    "between-animal variability."
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
    CONMED_TARIQUIDAR = list(
      description = "Tariquidar co-administration indicator (1 = P-gp inhibitor pre-treatment, 0 = quinidine alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no tariquidar",
      notes = paste(
        "Table I: tariquidar 15 mg/kg intravenously over 10 min, given before",
        "the quinidine infusion. Acts by Eq. 3, P = P_PAT * (1 + theta_cov *",
        "Cov), on the plasma-to-ECF clearance (theta 4.43, a 443 percent",
        "increase). Time-fixed per animal."
      ),
      source_name = "tariquidar"
    )
  )

  compartmentData <- list(
    central = list(analyte = "quinidine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "quinidine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "quinidine", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_icf = list(analyte = "quinidine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_csf_lv = list(analyte = "quinidine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "quinidine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "quinidine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "quinidine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 41,
    n_studies = 1,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "10 or 20 mg/kg intravenous infusion over 10 min, with or without tariquidar 15 mg/kg",
    regions = "Netherlands (Leiden University)",
    co_medication = "tariquidar 15 mg/kg IV over 10 min in part of the animals",
    notes = paste(
      "Table I: 41 rats, 10 and 20 mg/kg over 10 min, plasma sampled 0-360",
      "min and microdialysate (brain ECF, CSF lateral ventricle, CSF cisterna",
      "magna) 0-360 min; data previously published by Westerhout et al.",
      "(Table I reference 24). Unbound plasma fraction 14 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, quinidine column) ----
    lcl <- log(152); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 152 (RSE 2.40%)
    lvc <- log(194); label("Central volume V_PL (mL)") # Table III VPL 194 (RSE 37.5%)
    lq <- log(1070); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table III QPL_PER1 1070 (RSE 5.80%)
    lvp <- log(13300); label("Peripheral volume V_PER1 (mL)") # Table III VPER1 13300 (RSE 3.00%)

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.0354); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table III QPL_ECF 0.0354 (RSE 3.90%)
    lq_brain_ecf_brain_icf <- log(0.0250); label("Intercompartmental clearance between brain ECF and brain ICF Q_ECF_ICF (mL/min)") # Table III QECF_ICF 0.0250 (RSE 6.70%)
    lq_diff <- log(0.0237); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0237 (RSE 2.40%)

    # P-gp inhibition, Eq. 3: P = P_PAT * (1 + theta * CONMED_TARIQUIDAR)
    e_conmed_tariquidar_clin <- 4.43; label("Fractional increase in Q_PL_ECF with tariquidar (fraction)") # Table III impact on QPL_ECF 4.43 (RSE 2.80%), footnote b GF120918 or tariquidar

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_icf <- fixed(1.44); label("Brain ICF volume V_ICF (mL)") # Table III VICF 1.44 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 0.418; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 0.418 (RSE 4.90%)
    propSd_Cbrain_ecf <- 0.628; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 0.628 (RSE 3.30%)
    propSd_Cbrain_csf_lv <- 0.629; label("Proportional residual SD, CSF lateral ventricle (fraction)") # Table III sigma_CSFLV 0.629 (RSE 4.10%)
    propSd_Cbrain_csf_cm <- 0.466; label("Proportional residual SD, CSF cisterna magna (fraction)") # Table III sigma_CSFCM 0.466 (RSE 4.10%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    clin <- exp(lclin) * (1 + e_conmed_tariquidar_clin * CONMED_TARIQUIDAR)
    q_ecf_icf <- exp(lq_brain_ecf_brain_icf)
    q_diff <- exp(lq_diff)

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_ICF <- brain_icf / v_brain_icf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 -
      clin * C_PL + q_diff * C_SAS
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1

    # Brain ECF (exchanging with ICF) and CSF compartments in series (Fig. 2)
    d/dt(brain_ecf) <- clin * C_PL - q_diff * C_ECF -
      q_ecf_icf * C_ECF + q_ecf_icf * C_ICF
    d/dt(brain_icf) <- q_ecf_icf * C_ECF - q_ecf_icf * C_ICF
    d/dt(brain_csf_lv) <- q_diff * C_ECF - q_diff * C_LV
    d/dt(brain_csf_tfv) <- q_diff * C_LV - q_diff * C_TFV
    d/dt(brain_csf_cm) <- q_diff * C_TFV - q_diff * C_CM
    d/dt(brain_csf_sas) <- q_diff * C_CM - q_diff * C_SAS

    Cc <- C_PL
    Cbrain_ecf <- C_ECF
    Cbrain_icf <- C_ICF
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
