Yamamoto_2017b_morphine_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for morphine: an empirical three-compartment",
    "plasma model linked to a brain extracellular-fluid compartment and four",
    "cerebrospinal-fluid compartments (lateral ventricle, third + fourth",
    "ventricle, cisterna magna, subarachnoid space) of fixed rat physiological",
    "volume. Drug enters brain ECF from plasma by a one-way clearance and is",
    "carried through ECF and the CSF compartments in series, and back to",
    "plasma, by a single estimated drug-dispersion flow. The plasma-to-ECF",
    "clearance is higher at the 4 mg/kg dose than at 10 and 40 mg/kg",
    "(categorical dose effect) and is raised by P-glycoprotein inhibition with",
    "elacridar (GF120918). Plasma and brain ECF concentrations are unbound.",
    "Naive-pooled fit of two pooled studies: no between-animal variability."
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
    DOSE_LOW = list(
      description = "Low-dose indicator (1 = 4 mg/kg morphine, 0 = 10 or 40 mg/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = 10 or 40 mg/kg",
      notes = paste(
        "Brain-ECF morphine concentrations were not proportional to dose or",
        "plasma concentration, so the paper estimated a categorical dose effect:",
        "Table III gives Q_PL_ECF = 0.00458 mL/min for the 10 and 40 mg/kg",
        "groups and 0.00750 mL/min for the 4 mg/kg group (footnote d). Set to 1",
        "for the 4 mg/kg dose level only; time-fixed per animal (single dose)."
      ),
      source_name = "dosage of 4 mg/kg"
    ),
    CONMED_ELACRIDAR = list(
      description = "Elacridar (GF120918) co-administration indicator (1 = P-gp inhibitor infused, 0 = morphine alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no elacridar",
      notes = paste(
        "Table I: GF120918 6 mg/kg as a continuous intravenous infusion in the",
        "first morphine study (4, 10, 40 mg/kg). Acts by Eq. 3,",
        "P = P_PAT * (1 + theta_cov * Cov), on the plasma-to-ECF clearance",
        "(theta 1.62, a 162 percent increase). Time-fixed per animal."
      ),
      source_name = "GF120918"
    )
  )

  compartmentData <- list(
    central = list(analyte = "morphine", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "morphine", units = "ng", specimen = "tissue", verified = TRUE),
    brain_ecf = list(analyte = "morphine", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "morphine", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 83,
    n_studies = 2,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "4, 10 or 40 mg/kg intravenous infusion over 10 min, with or without elacridar (GF120918)",
    regions = "Netherlands (Leiden University)",
    co_medication = "GF120918 (elacridar) 6 mg/kg continuous IV infusion in part of the first study",
    notes = paste(
      "Table I: two published morphine studies pooled in one fit -- 65 rats at",
      "4, 10 and 40 mg/kg with or without GF120918 (Table I reference 70) and",
      "18 rats at 10 and 40 mg/kg (Table I reference 71); plasma sampled",
      "0-360 min and brain-ECF microdialysate 0-360 min. Unbound plasma",
      "fraction 83 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, morphine column) ----
    lcl <- log(21.6); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 21.6 (RSE 2.60%)
    lvc <- log(118); label("Central volume V_PL (mL)") # Table III VPL 118 (RSE 6.40%)
    lq <- log(8.72); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table III QPL_PER1 8.72 (RSE 3.80%)
    lvp <- log(1210); label("Peripheral 1 volume V_PER1 (mL)") # Table III VPER1 1210 (RSE 7.80%)
    lq2 <- log(53.3); label("Intercompartmental clearance Q_PL_PER2 (mL/min)") # Table III QPL_PER2 53.3 (RSE 5.80%)
    lvp2 <- log(570); label("Peripheral 2 volume V_PER2 (mL)") # Table III VPER2 570 (RSE 4.40%)

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.00458); label("One-way clearance from plasma to brain ECF Q_PL_ECF at 10 and 40 mg/kg (mL/min)") # Table III QPL_ECF 0.00458 (RSE 7.40%)
    # Table III reports a separate Q_PL_ECF of 0.00750 mL/min (RSE 8.90%) for
    # the 4 mg/kg group (footnote d); carried as the log-ratio to the
    # 10/40 mg/kg value so exp(lclin + e_dose_low_clin) returns 0.00750.
    e_dose_low_clin <- log(0.00750 / 0.00458); label("Log-ratio of Q_PL_ECF at 4 mg/kg to Q_PL_ECF at 10 and 40 mg/kg (unitless)") # Table III QPL_ECF 0.00750 at 4 mg/kg vs 0.00458
    lq_diff <- log(0.0200); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0200 (RSE 4.30%)

    # P-gp inhibition, Eq. 3: P = P_PAT * (1 + theta * CONMED_ELACRIDAR)
    e_conmed_elacridar_clin <- 1.62; label("Fractional increase in Q_PL_ECF with elacridar (fraction)") # Table III impact on QPL_ECF 1.62 (RSE 11.4%), footnote b GF120918 or tariquidar

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 0.647; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 0.647 (RSE 4.00%)
    propSd_Cbrain_ecf <- 0.779; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 0.779 (RSE 7.40%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    clin <- exp(lclin + e_dose_low_clin * DOSE_LOW) *
      (1 + e_conmed_elacridar_clin * CONMED_ELACRIDAR)
    q_diff <- exp(lq_diff)

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 -
      q2 * C_PL + q2 / vp2 * peripheral2 -
      clin * C_PL + q_diff * C_SAS
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1
    d/dt(peripheral2) <- q2 * C_PL - q2 / vp2 * peripheral2

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
