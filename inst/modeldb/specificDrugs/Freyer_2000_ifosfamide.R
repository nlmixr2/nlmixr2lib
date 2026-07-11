Freyer_2000_ifosfamide <- function() {
  description <- "One-compartment IV population PK model for ifosfamide in 24 small cell lung cancer patients on the AVI regimen (Freyer 2000), with day-1/day-2 clearance collapsed to the mean per operator instruction. Structural reduction from the paper's two-compartment model: only central CL and V were reported in Table 2, so the peripheral compartment is omitted; the paper's day-1 to day-2 autoinduction of clearance is also not represented in this extraction."
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
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

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
      "AVI regimen: ifosfamide 2000 mg/m^2 IV over 2 h on days 1 and 2.",
      "Total ifosfamide dose per course = 4000 mg/m^2. At the 1.79 m^2",
      "median BSA that is ~3580 mg per day-1 infusion."
    ),
    regions          = "France (Lyon Saint-Etienne Thoracic Oncology Group, GLOT; multicentre)",
    covariates_tested = paste(
      "Screened but not retained in the final ifosfamide model: age, sex,",
      "stage of disease, body weight, height, body surface area, serum",
      "creatinine, ASAT, ALAT, alkaline phosphatase, gamma-GT, total",
      "bilirubin, LDH, total protein (Freyer 2000 Methods; Results",
      "Ifosfamide). No covariate reached significance."
    ),
    notes            = paste(
      "Baseline characteristics from Freyer 2000 Table 1. Ifosfamide",
      "assayed by GC-NPD with trofosfamide internal standard (LOQ 0.5",
      "ug/mL, linear 0.5-30 ug/mL, CV 3.1-10%). Freyer 2000 fitted",
      "two separate clearance parameters -- CL_day1 = 5.6 L/h and CL_day2",
      "= 7.95 L/h (a 42% autoinduction increase) -- with no structural",
      "induction mechanism; each day's samples were fit under its own",
      "CL. Inter-occasion (course-to-course) variability was not",
      "estimable for ifosfamide."
    )
  )

  ini({
    # Structural PK (Freyer 2000 Table 2, ifosfamide row).
    #
    # The paper estimated two separate clearances -- CL_day1 = 5.6 L/h and
    # CL_day2 = 7.95 L/h -- describing an unmodelled 42% autoinduction of
    # ifosfamide metabolism between the two dosing days. This extraction
    # collapses the two values to their arithmetic mean per operator
    # instruction (sidecar 001 response 2026-06-21, q2 = BB): the
    # day-specific structure is dropped in favour of a single clearance.
    #
    # The 1-compartment reduction (from the paper's 2-cmt structure) is a
    # separate simplification -- required because the transfer rate
    # constants / peripheral volume are not reported anywhere in the paper.
    lcl <- log(6.775); label("Typical clearance, mean of day-1 and day-2 values (L/h)")  # (5.6 + 7.95)/2 = 6.775 L/h; Freyer 2000 Table 2 CL_day1 = 5.6 (s.d. 0.330) and CL_day2 = 7.95 (s.d. 0.450); mean per sidecar 001 q2 = BB
    lvc <- log(26.0);  label("Typical central volume of distribution (L)")               # Freyer 2000 Table 2: V_d = 26.0 L (s.d. 4.49) -- interpreted as V_central (V1) of the 2-cmt model, not V_ss

    # IIV -- Freyer 2000 Table 2 reports NONMEM omega variances.
    etalcl ~ 0.0100  # Freyer 2000 Table 2: omega^2_CL = 0.0100 (s.d. 0.0044); sqrt(exp(0.0100)-1) = 10.0 pct CV, matches text 10.1 pct for CL
    etalvc ~ 0.0296  # Freyer 2000 Table 2: omega^2_Vd = 0.0296 (s.d. 0.0084); sqrt(exp(0.0296)-1) = 17.3 pct CV, matches text 17.1 pct for Vd

    # Residual error -- combined proportional + additive (Freyer 2000
    # Methods Data analysis: "A mixed (proportional and additive) model was
    # used to estimate the residual error for doxorubicin and ifosfamide").
    # Table 2 values are NONMEM sigma variances; SD = sqrt(sigma^2).
    propSd <- 0.0648; label("Proportional residual error (fraction)")  # Freyer 2000 Table 2: sigma^2_prop = 0.0042 (s.d. 0.0043) -> SD = sqrt(0.0042) = 0.0648 (6.5%)
    addSd  <- 0.4583; label("Additive residual error (mg/L)")           # Freyer 2000 Table 2: sigma^2_add = 0.2100 -> SD = sqrt(0.2100) = 0.4583 mg/L
  })

  model({
    # Individual PK parameters (log-normal IIV, no covariates retained).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # ODE system -- 1-compartment IV (ifosfamide 2-h IV infusion on days 1
    # and 2). Dose is delivered directly into the central compartment.
    d/dt(central) <- -kel * central

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
