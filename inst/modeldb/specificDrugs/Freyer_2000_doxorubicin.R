Freyer_2000_doxorubicin <- function() {
  description <- "One-compartment IV population PK model for doxorubicin in 24 small cell lung cancer patients on the AVI regimen (Freyer 2000). Structural reduction from the paper's three-compartment model: only central CL and V were reported in Table 2, so the peripheral compartments are omitted."
  reference <- paste(
    "Freyer G, Tranchand B, Ligneau B, Ardiet C, Souquet P-J,",
    "Court-Fortune I, Riou R, Rebattu P, Boissel J-P,",
    "Trillet-Lenoir V, Girard P.",
    "Population pharmacokinetics of doxorubicin, etoposide and ifosfamide",
    "in small cell lung cancer patients: results of a multicentre study.",
    "Br J Clin Pharmacol. 2000;50(4):315-324.",
    "doi:10.1046/j.1365-2125.2000.00269.x.",
    sep = " "
  )
  vignette <- "Freyer_2000_AVI_regimen"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 24L,
    n_studies        = 1L,
    n_courses        = 47L,
    age_range        = "45-70 years",
    age_median       = "57 years",
    weight_range     = "54-93 kg",
    weight_median    = "69 kg",
    height_range     = "159-179 cm",
    height_median    = "170 cm",
    bsa_range        = "1.54-2.10 m^2",
    bsa_median       = "1.79 m^2",
    sex_female_pct   = "not reported",
    race_ethnicity   = "not reported (French multicentre study)",
    disease_state    = paste(
      "Small cell lung cancer (SCLC), either limited to the thorax or",
      "extensive; no severe hepatic or renal impairment (patients with",
      "ASAT/ALAT/alkaline phosphatase > 2x upper limit of normal and/or",
      "serum creatinine > 1.5x upper limit of normal were excluded from",
      "the primary therapeutic trials)."
    ),
    dose_range       = paste(
      "AVI regimen: doxorubicin 50 mg/m^2 IV over 15 min on day 1;",
      "etoposide 120 mg/m^2 IV over 30 min on days 1, 2, and 3;",
      "ifosfamide 2000 mg/m^2 IV over 2 h on days 1 and 2.",
      "Doxorubicin doses at 1.79 m^2 median BSA correspond to ~90 mg per",
      "day-1 infusion. Antiemetic support: ondansetron 8 mg/day and",
      "methylprednisolone 120 mg/day, days 1-3."
    ),
    regions          = "France (Lyon Saint-Etienne Thoracic Oncology Group, GLOT; multicentre)",
    covariates_tested = paste(
      "Screened but not retained in the final doxorubicin model: age, sex,",
      "stage of disease, body weight, height, body surface area, serum",
      "creatinine, ASAT, ALAT, alkaline phosphatase, gamma-GT, total",
      "bilirubin, LDH, total protein (Freyer 2000 Methods Data analysis;",
      "Results Doxorubicin). No covariate demonstrated a significant impact",
      "on the doxorubicin CL, V, or AUC in this cohort."
    ),
    notes            = paste(
      "Baseline characteristics from Freyer 2000 Table 1. Serum creatinine",
      "median 82 umol/L (range 44-144), ASAT median 29.5 IU/L (range 9-152),",
      "ALAT median 28.5 IU/L (range 3-222), ALKP median 102.5 IU/L (range",
      "28-869), GGT median 67 IU/L (range 21-374), bilirubin median 8.5",
      "umol/L (range 5-16), LDH median 557.5 IU/L (range 277-2524), total",
      "protein median 70 g/L (range 56-79). Sampling strategy: 19 samples",
      "per course in the first 7 patients (extensive sampling), reduced to",
      "10 samples per course thereafter (limited sampling, D-optimal). Data",
      "were analysed with NONMEM V (FO method)."
    )
  )

  ini({
    # Structural PK (Freyer 2000 Table 2, doxorubicin row). Peripheral
    # compartments dropped because k12, k21, k13, k31 (or Q, Vp) are not
    # reported anywhere in the paper. The 1-compartment reduction preserves
    # the reported CL and central V and the AUC = Dose/CL relationship
    # (matches Table 3 systemic exposures); early-time concentration
    # profiles and terminal half-life are NOT recovered.
    lcl <- log(54.0); label("Typical clearance (L/h)")                     # Freyer 2000 Table 2: CL = 54.0 L/h (s.d. 4.96); abstract typo 32.0 superseded by Table 2 / Results / Discussion
    lvc <- log(9.3);  label("Typical central volume of distribution (L)")  # Freyer 2000 Table 2: V_d = 9.3 L (s.d. 0.970) -- interpreted as V_central (V1) of the 3-cmt model, not V_ss

    # IIV -- Freyer 2000 Table 2 reports NONMEM omega variances directly;
    # the in-text %CV matches sqrt(exp(omega^2) - 1) exactly, confirming
    # the tabulated values are variances on the log-normal scale.
    etalcl ~ 0.0296  # Freyer 2000 Table 2: omega^2_CL = 0.0296 (s.d. 0.0012); sqrt(exp(0.0296)-1) = 17.3 pct CV, matches text 17.2 pct for CL
    etalvc ~ 0.0369  # Freyer 2000 Table 2: omega^2_Vd = 0.0369 (s.d. 0.0020); sqrt(exp(0.0369)-1) = 19.4 pct CV, matches text 19.2 pct for Vd

    # Residual error -- combined proportional + additive (Freyer 2000
    # Methods Data analysis: "A mixed (proportional and additive) model was
    # used to estimate the residual error for doxorubicin and ifosfamide").
    # Table 2 values are NONMEM sigma variances; SD = sqrt(sigma^2).
    propSd <- 0.1183; label("Proportional residual error (fraction)")  # Freyer 2000 Table 2: sigma^2_prop = 0.0140 -> SD = sqrt(0.0140) = 0.1183 (11.8%)
    addSd  <- 0.0387; label("Additive residual error (mg/L)")           # Freyer 2000 Table 2: sigma^2_add = 0.0015 -> SD = sqrt(0.0015) = 0.0387 mg/L
  })

  model({
    # Individual PK parameters (log-normal IIV, no covariates retained).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # ODE system -- 1-compartment IV (doxorubicin 15-min IV infusion, day 1
    # only). Dose is delivered directly into the central compartment.
    d/dt(central) <- -kel * central

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
