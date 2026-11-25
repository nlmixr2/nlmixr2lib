Yamamoto_2017b_risperidone_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for risperidone: an empirical one-compartment",
    "plasma model linked to a brain extracellular-fluid compartment and four",
    "cerebrospinal-fluid compartments (lateral ventricle, third + fourth",
    "ventricle, cisterna magna, subarachnoid space) of fixed rat physiological",
    "volume. Drug enters brain ECF from plasma by a one-way clearance and is",
    "carried through ECF and the CSF compartments in series, and back to",
    "plasma, by a single estimated drug-dispersion flow. P-glycoprotein",
    "inhibition with tariquidar raises the plasma-to-ECF clearance. Plasma,",
    "brain ECF and CSF concentrations are unbound. Naive-pooled fit: no",
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
      description = "Tariquidar co-administration indicator (1 = P-gp inhibitor pre-treatment, 0 = risperidone alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no tariquidar",
      notes = paste(
        "Table I: tariquidar 15 mg/kg intravenously over 10 min, given before",
        "the risperidone infusion. Acts by Eq. 3, P = P_PAT * (1 + theta_cov *",
        "Cov), on the plasma-to-ECF clearance (theta 1.24, a 124 percent",
        "increase). Time-fixed per animal."
      ),
      source_name = "tariquidar"
    )
  )

  compartmentData <- list(
    central = list(analyte = "risperidone", units = "ng", specimen = "plasma", verified = TRUE),
    brain_ecf = list(analyte = "risperidone", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "risperidone", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "risperidone", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "risperidone", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "risperidone", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 16,
    n_studies = 1,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "2 mg/kg intravenous infusion over 20 min, with or without tariquidar 15 mg/kg",
    regions = "Netherlands (Leiden University)",
    co_medication = "tariquidar 15 mg/kg IV over 10 min in part of the animals",
    notes = paste(
      "Table I: 16 rats, 2 mg/kg over 20 min (newly produced data), plasma",
      "sampled 0-360 min and microdialysate (brain ECF, CSF cisterna magna)",
      "0-240 min. Unbound plasma fraction 7 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, risperidone column; no peripheral compartment) ----
    lcl <- log(465); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 465 (RSE 13.0%)
    lvc <- log(60000); label("Central volume V_PL (mL)") # Table III VPL 60000 (RSE 13.0%)

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.0247); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table III QPL_ECF 0.0247 (RSE 18.6%)
    lq_diff <- log(0.0254); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0254 (RSE 15.0%)

    # P-gp inhibition, Eq. 3: P = P_PAT * (1 + theta * CONMED_TARIQUIDAR)
    e_conmed_tariquidar_clin <- 1.24; label("Fractional increase in Q_PL_ECF with tariquidar (fraction)") # Table III impact on QPL_ECF 1.24 (RSE 16.3%), footnote b GF120918 or tariquidar

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 1.44; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 1.44 (RSE 5.90%)
    propSd_Cbrain_ecf <- 0.911; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 0.911 (RSE 6.80%)
    propSd_Cbrain_csf_cm <- 0.827; label("Proportional residual SD, CSF cisterna magna (fraction)") # Table III sigma_CSFCM 0.827 (RSE 9.40%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    clin <- exp(lclin) * (1 + e_conmed_tariquidar_clin * CONMED_TARIQUIDAR)
    q_diff <- exp(lq_diff)

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    d/dt(central) <- -cl * C_PL - clin * C_PL + q_diff * C_SAS

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
    Cbrain_csf_cm ~ prop(propSd_Cbrain_csf_cm)
  })
}
