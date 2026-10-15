Chauzy_2022_metronidazole_pbpk <- function() {
  description <- "PBPK (minimal, CNS). Unbound metronidazole in plasma, brain extracellular fluid (ECF) and cerebrospinal fluid (CSF) of brain-injured neuro-ICU adults (Chauzy 2022): a blood compartment with a linear total blood clearance, exchanging by perfusion with one lumped non-CNS tissue compartment and with a three-compartment CNS (brain vasculature, brain ECF, cranial CSF). The CNS compartments are linked by passive permeability-surface-area products across the blood-brain and blood-CSF barriers, ECF-to-CSF bulk flow, CSF sink flow back to the brain vasculature, and drainage of CSF through an external ventricular drain (EVD). Volumes and flows are fixed physiological values; only the tissue partition coefficient, clearance and residual errors were estimated."
  reference <- "Chauzy A, Bouchene S, Aranzana-Climent V, Clarhaut J, Adier C, Gregoire N, Couet W, Dahyot-Fizelier C, Marchand S. A Minimal Physiologically Based Pharmacokinetic Model to Characterize CNS Distribution of Metronidazole in Neuro Care ICU Patients. Antibiotics (Basel). 2022;11(10):1293. doi:10.3390/antibiotics11101293"
  vignette <- "Chauzy_2022_metronidazole_pbpk"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  # Two bookkeeping states, not biological compartments. auc_brain_ecf
  # integrates the brain ECF concentration so the interval-averaged
  # microdialysate concentration of Eq. 5 can be read off the solve as
  # (auc_brain_ecf(t2) - auc_brain_ecf(t1)) / (t2 - t1). evd is the cumulative
  # amount of drug drained into the EVD collection bag (Eq. 8); the bag
  # concentration for a collection interval is the amount drained over that
  # interval divided by the volume collected.
  paper_specific_compartments <- c("auc_brain_ecf", "evd")

  compartmentData <- list(
    central = list(analyte = "metronidazole", units = "mg", specimen = "whole blood", verified = TRUE),
    res_tis = list(analyte = "metronidazole", units = "mg", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "metronidazole", units = "mg", specimen = "whole blood", verified = TRUE),
    brain_ecf = list(analyte = "metronidazole", units = "mg", specimen = "brain ISF", verified = TRUE),
    brain_csf = list(analyte = "metronidazole", units = "mg", specimen = "CSF", verified = TRUE),
    auc_brain_ecf = list(analyte = "metronidazole", units = "mg*h/L", specimen = "brain ISF", verified = TRUE),
    evd = list(analyte = "metronidazole", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "REQUIRED model input. Not a covariate effect: the non-CNS tissue volume is",
        "the physiological remainder V_tissue = TBW - V_blood - V_brain,vasc - V_ECF - V_CSF",
        "(Methods text under Eq. 2, citing Cao and Jusko 2012), with body weight in kg read",
        "as litres (density 1 kg/L). Cohort 75-115 kg (Table 3)."
      ),
      source_name = "TBW"
    ),
    CSF_DRAIN_VOL_24H = list(
      description = "Volume of cerebrospinal fluid removed through the external ventricular drain, expressed per 24 hours",
      units = "mL/24h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "REQUIRED model input (0 for a patient without an EVD). The paper's covariate is",
        "the EVD flow QEVD in L/h, measured per collection interval (0.5-1 h) and",
        "time-varying (Table S3: 0 to 0.040 L/h). The model uses",
        "QEVD = CSF_DRAIN_VOL_24H / 1000 / 24, so a source QEVD in L/h is supplied as",
        "QEVD x 24000 (e.g. 0.010 L/h = 240 mL/24h). QEVD removes drug from the CSF",
        "compartment (Eq. 6) and lowers the CSF sink flow, Qsink = Qsink,physio - QEVD",
        "(Eq. 7), floored at 0 when QEVD exceeds Qsink,physio = 0.024 L/h (Figure 4",
        "caption and Discussion)."
      ),
      source_name = "QEVD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "34-73 years",
    weight_range = "75-115 kg",
    height_range = "172-180 cm",
    sex_female_pct = 0,
    disease_state = "Brain-injured adults in a neurointensive care unit (4 traumatic brain injury, 3 subarachnoid haemorrhage, 1 ventricular haemorrhage), treated with metronidazole for a lung infection. Four patients had a brain microdialysis probe in a frontal lobe (plasma + brain ECF sampling) and four had an external ventricular drain in a lateral ventricle (plasma + EVD-collected CSF sampling).",
    renal_function = "Creatinine clearance 84-306 mL/min (Table 3) -- normal to augmented.",
    dose_range = "Metronidazole 500 mg every 8 h as a 30-min intravenous infusion; pharmacokinetic sampling at steady state after at least 2 days of treatment.",
    regions = "Single centre: neurointensive care unit, University Hospital of Poitiers, France.",
    notes = "Plasma unbound concentrations by ultrafiltration (7-13 blood samples per patient over one dosing interval); brain dialysates (corrected by the per-patient in vivo probe recovery) and EVD-collected CSF over 8 h at 0.5 h intervals for the first 2-4 h and 1 h intervals thereafter. Measured unbound plasma concentrations were converted to unbound blood concentrations with the PK-Sim-predicted blood-to-plasma ratio for fitting. Estimation by FOCE-I in NONMEM 7.4; 95% CIs by sampling importance resampling. Data previously reported by Frasca 2014 (Antimicrob Agents Chemother 58:1019-1023 and 58:1024-1027)."
  )

  ini({
    # ---- Estimated drug-specific parameters (Table 1) ----
    lcl <- log(7.28); label("Total blood clearance of unbound drug (L/h)") # Table 1: CL 7.28 L/h (95% CI 5.77-9.53)
    lkp_rest <- log(0.796); label("Unbound tissue-to-blood partition coefficient, non-CNS tissue (unitless)") # Table 1: Kp 0.796 (0.693-0.923)
    fd <- fixed(0.86); label("Fraction of cardiac output perfusing the non-CNS tissue compartment (unitless)") # Table 1: fd 0.86, footnote b: fixed to the maximum allowed value (identifiability)
    ps_ecf <- fixed(6.4); label("Passive permeability-surface-area product across the blood-brain barrier (L/h)") # Table 1: PSECF 6.4, footnote c: Simcyp prediction from BBB surface area, log P and MW
    ps_csf <- fixed(3.2); label("Passive permeability-surface-area product across the blood-CSF barrier (L/h)") # Table 1: PSCSF 3.2, footnote d: assumed half of PSECF

    # ---- System-specific physiological parameters (Table 4) ----
    co <- fixed(312); label("Cardiac output (L/h)") # Table 4: CO 312 L/h
    q_brain <- fixed(42); label("Brain blood flow (L/h)") # Table 4: Qbrain 42 L/h
    q_bulk <- fixed(0.0105); label("Bulk flow from brain ECF to CSF (L/h)") # Table 4: Qbulk 0.0105 L/h
    q_sink_physio <- fixed(0.024); label("Physiological CSF sink (absorption) flow (L/h)") # Table 4: Qsink,physio 0.024 L/h
    v_blood <- fixed(5.85); label("Blood volume (L)") # Table 4: Vblood 5.85 L
    v_brain_vasc <- fixed(0.0637); label("Brain vascular volume (L)") # Table 4: Vbrain,vasc 0.0637 L
    v_ecf <- fixed(0.24); label("Brain ECF volume (L)") # Table 4: VECF 0.24 L
    v_csf <- fixed(0.130); label("Cranial CSF volume (L)") # Table 4: VCSF 0.130 L

    # ---- Drug-specific input (Table 4) ----
    bpr <- fixed(0.82); label("Blood-to-plasma concentration ratio (unitless)") # Table 4: BP 0.82 (PK-Sim prediction)

    # ---- Between-subject variability (Table 1) ----
    # IIV on CL only; Table 1 prints 35.2 %CV (95% CI 24.0-57.8), converted as
    # omega^2 = log(0.352^2 + 1).
    etalcl ~ 0.11681

    # ---- Residual error (Table 1) ----
    propSd <- 0.144; label("Proportional residual error, unbound plasma (fraction)") # Table 1: sigma prop,plasma 14.4% (9.74-19.1)
    addSd <- 1.18; label("Additive residual error, unbound plasma (mg/L)") # Table 1: sigma add,plasma 1.18 ug/mL (0.320-2.90)
    propSd_Cecf <- 0.228; label("Proportional residual error, brain ECF dialysate (fraction)") # Table 1: sigma prop,ECF 22.8% (17.5-29.1)
    propSd_Ccsf <- 0.282; label("Proportional residual error, EVD-collected CSF (fraction)") # Table 1: sigma prop,CSF 28.2% (22.0-36.5)
  })
  model({
    cl <- exp(lcl + etalcl)
    kp <- exp(lkp_rest)

    # Non-CNS tissue volume closes the body-weight balance (Methods, under Eq. 2):
    # V_blood + V_tissue + V_brain,vasc + V_ECF + V_CSF = TBW.
    v_tissue <- WT - v_blood - v_brain_vasc - v_ecf - v_csf

    # EVD flow (L/h) and the CSF sink flow it displaces (Eq. 7), floored at 0
    # when the drain exceeds the physiological sink flow (Figure 4 caption).
    q_evd <- CSF_DRAIN_VOL_24H / 1000 / 24
    q_sink <- q_sink_physio - q_evd
    if (q_sink < 0) q_sink <- 0

    c_blood <- central / v_blood
    c_tissue <- res_tis / v_tissue
    c_brain_vasc <- brain_vascular / v_brain_vasc
    c_ecf <- brain_ecf / v_ecf
    c_csf <- brain_csf / v_csf

    # Eq. 1: blood; the IV infusion enters here.
    d/dt(central) <- fd * co * c_tissue / kp + q_brain * c_brain_vasc - c_blood * (fd * co + q_brain + cl)
    # Eq. 2: perfusion-limited, well-stirred non-CNS tissue.
    d/dt(res_tis) <- fd * co * (c_blood - c_tissue / kp)
    # Eq. 3: brain vasculature.
    d/dt(brain_vascular) <- q_brain * c_blood + ps_ecf * c_ecf + c_csf * (ps_csf + q_sink) - c_brain_vasc * (q_brain + ps_ecf + ps_csf)
    # Eq. 4: brain ECF.
    d/dt(brain_ecf) <- c_brain_vasc * ps_ecf - c_ecf * (ps_ecf + q_bulk)
    # Eq. 6: cranial CSF.
    d/dt(brain_csf) <- c_brain_vasc * ps_csf + c_ecf * q_bulk - c_csf * (ps_csf + q_sink + q_evd)
    # Eq. 5 integrand: cumulative ECF AUC for interval-averaged dialysate.
    d/dt(auc_brain_ecf) <- c_ecf
    # Eq. 8: cumulative drug amount drained into the EVD collection bag.
    d/dt(evd) <- c_csf * q_evd

    # Observations. The model state is the unbound blood concentration; the
    # unbound plasma concentration is that divided by the blood-to-plasma ratio
    # (Methods: B/P 'used to convert measured unbound plasma concentrations
    # into unbound blood concentrations in the model').
    Cc <- c_blood / bpr
    Cecf <- c_ecf
    Ccsf <- c_csf

    Cc ~ add(addSd) + prop(propSd)
    Cecf ~ prop(propSd_Cecf)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
