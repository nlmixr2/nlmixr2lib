Murinova_2022_oxacillin <- function() {
  description <- "One-compartment IV-infusion population PK model for oxacillin in 24 adult patients with staphylococcal infections, fitted to routine therapeutic-drug-monitoring total plasma concentrations (Murinova 2022). The model is parameterized by the volume of distribution (11.2 L, no covariate) and the elimination rate constant Ke, which carries an UNCENTERED exponential effect of the creatinine-based CKD-EPI eGFR expressed in mL/s/1.73 m2: Ke = 0.73 * exp(0.3 * eGFR) 1/h, i.e. 1.15 1/h (t1/2 0.6 h) at the cohort median eGFR of 1.53 mL/s/1.73 m2. Log-normal IIV on V and Ke; proportional residual error. An unbound concentration Cu = fu * Cc uses the cohort median measured protein binding of 86% (fu = 0.14), which the paper used to evaluate fT > MIC targets for intermittent, extended and continuous infusion."
  reference <- "Murinova I, Svidrnoch M, Gucky T, Hlavac J, Michalek P, Slanar O, Sima M. Population Pharmacokinetic Analysis Proves Superiority of Continuous Infusion in PK/PD Target Attainment with Oxacillin in Staphylococcal Infections. Antibiotics (Basel). 2022;11(12):1736. doi:10.3390/antibiotics11121736"
  vignette <- "Murinova_2022_oxacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "oxacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine-based CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Murinova 2022 Methods 4.2: 'the glomerular filtration rate (eGFR) was calculated according to",
        "the chronic kidney disease epidemiology collaboration formula.' UNITS: the paper reports eGFR in",
        "mL/s/1.73 m^2 (SI, the Czech clinical convention), NOT in the canonical mL/min/1.73 m^2 of this",
        "register. Table 1 gives median 1.53, IQR 0.98-1.71, range 0.59-2.10 mL/s/1.73 m^2, i.e. median",
        "91.8, IQR 58.8-102.6, range 35.4-126.0 mL/min/1.73 m^2. This column carries the CANONICAL",
        "mL/min/1.73 m^2 value and model() divides it by 60 before applying the published coefficient",
        "0.3, which is per mL/s/1.73 m^2.",
        "The effect is EXPONENTIAL and UNCENTERED on Ke (Results 2.1 equation:",
        "Log(Ke) = log(Ke_pop) + beta_Ke_eGFR x eGFR + eta_Ke), so it is 1 at eGFR = 0. Patients with",
        "severe renal impairment, augmented renal clearance, renal replacement therapy or extracorporeal",
        "life support were not in the cohort, so the model should not be extrapolated outside",
        "roughly 35-126 mL/min/1.73 m^2.",
        collapse = " "
      ),
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as a continuous covariate (Methods 4.5 step 2; Table 1 median 55, IQR 45-72, range 26-84 years) and not retained. No point estimate published."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested as a continuous covariate (Table 1 median 84, IQR 74-96, range 57-145 kg) and not retained; Vd carries no body-size scaling. No point estimate published."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Tested (Table 1 median 174, range 153-195 cm) and not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested (Table 1 median 28, range 20-39 kg/m^2) and not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested (Table 1 median 82, range 50-151 umol/L) and not retained in favour of the CKD-EPI eGFR derived from it."
    ),
    BUN = list(
      description = "Serum urea (the canonical BUN column accepts mmol/L urea)",
      units = "mmol/L",
      type = "continuous",
      notes = "Tested (Table 1 median 4.5, range 1.8-11.8 mmol/L urea) and not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested on the PK parameters (Table 1 median 32, range 22-44 g/L) and not retained; also not significantly associated with oxacillin plasma protein binding (Results 2, linear regression)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Tested as a categorical covariate (10 of 24 patients female) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 32L,
    age_range = "26-84 years",
    age_median = "55 years",
    weight_range = "57-145 kg",
    weight_median = "84 kg",
    sex_female_pct = 41.7,
    disease_state = "Staphylococcal infection treated with intravenous oxacillin: central nervous system (n = 7), orthopedic (n = 10), sepsis (n = 3), other (n = 4; e.g. bacteriuria, bacteremia, endocarditis). Methicillin-susceptible Staphylococcus aureus in 21 patients, S. hominis in 2, S. warneri in 1; MIC median 0.25 mg/L (IQR 0.25-0.5).",
    renal_function = "eGFR (CKD-EPI, creatinine) median 1.53 mL/s/1.73 m2 (IQR 0.98-1.71, range 0.59-2.10), i.e. median 91.8 mL/min/1.73 m2 (range 35.4-126.0). Patients on renal replacement therapy or extracorporeal life support were excluded.",
    dose_range = "1-3 g every 4 h as a 0.5-h intravenous infusion; one patient received a 3-h extended infusion.",
    regions = "Czech Republic (single centre: Military University Hospital Prague)",
    notes = "Retrospective TDM study, June 2021 to June 2022 (Methods 4.1). 32 total plasma concentrations (24 troughs, 8 peaks after the end of infusion; 1.33 per patient). Total and ultrafiltrate-unbound oxacillin measured by LC-MS/MS; median protein binding 86% (IQR 83-88%, range 74-97%), not associated with any tested covariate. Estimation in Monolix 2021R1 (SAEM); Monte Carlo simulations in Simulx 2021. Demographics in Table 1."
  )

  ini({
    # Structural parameters -- Murinova 2022 Table 2 "Fixed effects". The model
    # is parameterized in V and Ke (Results 2.1), not in CL.
    lvc <- log(11.2); label("Volume of distribution (L)")                          # Table 2, Vd_pop = 11.2 L (R.S.E. 30.4%)
    # lkel is Ke extrapolated to eGFR = 0 (uncentered exponential covariate).
    lkel <- log(0.73); label("Elimination rate constant at eGFR = 0 (1/h)")       # Table 2, Ke_pop = 0.73 1/h (R.S.E. 26.1%)

    # Covariate effect on Ke. Results 2.1 equation:
    #   Log(Ke) = log(Ke_pop) + beta_Ke_eGFR x eGFR + eta_Ke
    # so the effect is exponential in eGFR (mL/s/1.73 m^2). Table 2 and the
    # prose describe it as an additive 0.3 1/h per unit eGFR; the equation (the
    # Monolix covariate form) is encoded here -- see the vignette.
    e_crcl_kel <- 0.3; label("Exponential eGFR effect on Ke (per mL/s/1.73 m^2)") # Table 2, beta_Ke_eGFR = 0.3 (R.S.E. 50.8%)

    # IIV. Table 2 heads these rows "Standard deviation of the random effects",
    # so the published values are SDs on the log scale and are squared here.
    etalvc ~ 0.49   # Table 2, omega_Vd = 0.7 SD -> 0.7^2 = 0.49 (R.S.E. 22.5%)
    etalkel ~ 0.0121 # Table 2, omega_Ke = 0.11 SD -> 0.11^2 = 0.0121 (R.S.E. 37.5%)

    # Protein binding. Not a fitted model parameter: the cohort median of the
    # measured binding (Results 2: 'The median (IQR) value of oxacillin binding
    # to plasma proteins was 86% (83-88%)'), used to express the unbound
    # concentration on which the paper's fT > MIC targets are defined.
    fu <- fixed(0.14); label("Fraction of oxacillin unbound in plasma (unitless)") # Results 2, median binding 86% -> fu = 1 - 0.86

    # Residual error. Results 2.1: 'A proportional error model was the most
    # accurate for residual ... variability'.
    propSd <- 0.4; label("Proportional residual error (fraction)")                # Table 2, Error model parameter b = 0.4 (R.S.E. 19.8%)
  })

  model({
    # CRCL is supplied in the canonical mL/min/1.73 m^2; the published
    # coefficient 0.3 is per mL/s/1.73 m^2, so convert first.
    eGFRsi <- CRCL / 60

    # Murinova 2022 Results 2.1:
    #   Log(Vd) = log(Vd_pop) + eta_Vd
    #   Log(Ke) = log(Ke_pop) + beta_Ke_eGFR x eGFR + eta_Ke
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + e_crcl_kel * eGFRsi + etalkel)
    cl <- kel * vc

    # One compartment, linear elimination, IV infusion only (Results 2.1).
    d/dt(central) <- -kel * central

    # The assay measured TOTAL plasma oxacillin, so Cc is the fitted output;
    # Cu is the unbound concentration used for the fT > MIC targets.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
