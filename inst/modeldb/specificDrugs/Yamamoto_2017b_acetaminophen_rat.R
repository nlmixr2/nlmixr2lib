Yamamoto_2017b_acetaminophen_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Semi-physiological multi-compartmental",
    "brain distribution model for acetaminophen: an empirical two-compartment",
    "plasma model with enterohepatic recirculation, linked to a brain",
    "extracellular-fluid compartment and four cerebrospinal-fluid compartments",
    "(lateral ventricle, third + fourth ventricle, cisterna magna,",
    "subarachnoid space) of fixed rat physiological volume. Drug enters brain",
    "ECF from plasma by a one-way clearance and is carried through ECF and the",
    "CSF compartments in series, and back to plasma, by a single estimated",
    "drug-dispersion flow. Plasma, brain ECF and CSF concentrations are unbound.",
    "Naive-pooled fit: no between-animal variability."
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
    central = list(analyte = "acetaminophen", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acetaminophen", units = "ng", specimen = "tissue", verified = TRUE),
    gallbladder = list(analyte = "acetaminophen", units = "ng", specimen = "bile", verified = TRUE),
    brain_ecf = list(analyte = "acetaminophen", units = "ng", specimen = "brain ISF", verified = TRUE),
    brain_csf_lv = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_tfv = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_cm = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE),
    brain_csf_sas = list(analyte = "acetaminophen", units = "ng", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 16,
    n_studies = 1,
    weight_range = "225-275 g",
    disease_state = "healthy",
    dose_range = "16 mg/kg intravenous infusion over 10 min",
    regions = "Netherlands (Leiden University)",
    notes = paste(
      "Table I: 16 rats, 16 mg/kg over 10 min, plasma sampled 0-240 min and",
      "microdialysate (brain ECF, CSF lateral ventricle, CSF cisterna magna)",
      "0-240 min; data previously published by Westerhout et al. (Table I",
      "reference 6). Unbound plasma fraction 81 percent (Table SII)."
    )
  )

  ini({
    # ---- Plasma (Table III, acetaminophen column) ----
    lcl <- log(15.9); label("Plasma clearance CL_PL (mL/min)") # Table III CLPL 15.9 (RSE 4.80%)
    lvc <- log(65.7); label("Central volume V_PL (mL)") # Table III VPL 65.7 (RSE 25.4%)
    lq <- log(29.2); label("Intercompartmental clearance Q_PL_PER1 (mL/min)") # Table III QPL_PER1 29.2 (RSE 19.9%)
    lvp <- log(219); label("Peripheral volume V_PER1 (mL)") # Table III VPER1 219 (RSE 8.90%)

    # Enterohepatic recirculation. Table III 'fraction' is the percentage of
    # drug reabsorbed by enterohepatic circulation. The paper does not write
    # the recirculation equations; the gallbladder structure and the 120-min
    # switch-on time follow the same group's control stream for the same rat
    # acetaminophen data (Yamamoto 2017 CPT Pharmacometrics Syst Pharmacol,
    # Supplementary Material S3).
    f_ehc <- 0.933; label("Fraction of eliminated drug reabsorbed by enterohepatic circulation (fraction)") # Table III fraction 93.3% (RSE 1.20%)
    t_ehc <- fixed(120); label("Time after dose at which enterohepatic recirculation starts (min)") # not in this paper; sibling control stream S3

    # ---- Brain distribution (Table III) ----
    lclin <- log(0.0281); label("One-way clearance from plasma to brain ECF Q_PL_ECF (mL/min)") # Table III QPL_ECF 0.0281 (RSE 12.1%)
    lq_diff <- log(0.0556); label("Drug dispersion flow through brain ECF and CSF Q_DIFF (mL/min)") # Table III QDIFF 0.0556 (RSE 10.3%)

    # Rat physiological volumes, fixed (Table III footnote e)
    v_brain_ecf <- fixed(0.29); label("Brain ECF volume V_ECF (mL)") # Table III VECF 0.29 FIX
    v_brain_csf_lv <- fixed(0.05); label("Lateral-ventricle CSF volume V_LV (mL)") # Table III VLV 0.05 FIX
    v_brain_csf_tfv <- fixed(0.05); label("Third + fourth ventricle CSF volume V_TFV (mL)") # Table III VTFV 0.05 FIX
    v_brain_csf_cm <- fixed(0.017); label("Cisterna-magna CSF volume V_CM (mL)") # Table III VCM 0.017 FIX
    v_brain_csf_sas <- fixed(0.18); label("Subarachnoid-space CSF volume V_SAS (mL)") # Table III VSAS 0.18 FIX

    # ---- Residual error (Table III; proportional, Eq. 7) ----
    propSd <- 0.341; label("Proportional residual SD, plasma (fraction)") # Table III sigma_plasma 0.341 (RSE 8.60%)
    propSd_Cbrain_ecf <- 1.88; label("Proportional residual SD, brain ECF (fraction)") # Table III sigma_brainECF 1.88 (RSE 6.70%)
    propSd_Cbrain_csf_lv <- 0.607; label("Proportional residual SD, CSF lateral ventricle (fraction)") # Table III sigma_CSFLV 0.607 (RSE 6.40%)
    propSd_Cbrain_csf_cm <- 0.640; label("Proportional residual SD, CSF cisterna magna (fraction)") # Table III sigma_CSFCM 0.640 (RSE 8.20%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    clin <- exp(lclin)
    q_diff <- exp(lq_diff)

    frac_ehc <- 0
    if (t >= t_ehc) frac_ehc <- f_ehc

    C_PL <- central / vc
    C_ECF <- brain_ecf / v_brain_ecf
    C_LV <- brain_csf_lv / v_brain_csf_lv
    C_TFV <- brain_csf_tfv / v_brain_csf_tfv
    C_CM <- brain_csf_cm / v_brain_csf_cm
    C_SAS <- brain_csf_sas / v_brain_csf_sas

    # Plasma: elimination feeds the gallbladder, which empties at the
    # elimination rate constant; the fraction frac_ehc of the emptied amount
    # returns to plasma. The brain loop leaves plasma at clin and returns from
    # the subarachnoid space at q_diff.
    d/dt(central) <- -cl * C_PL - q * C_PL + q / vp * peripheral1 +
      (cl / vc) * frac_ehc * gallbladder -
      clin * C_PL + q_diff * C_SAS
    d/dt(peripheral1) <- q * C_PL - q / vp * peripheral1
    d/dt(gallbladder) <- cl * C_PL - (cl / vc) * gallbladder

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
    Cbrain_csf_lv ~ prop(propSd_Cbrain_csf_lv)
    Cbrain_csf_cm ~ prop(propSd_Cbrain_csf_cm)
  })
}
