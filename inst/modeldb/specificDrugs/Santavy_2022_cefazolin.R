Santavy_2022_cefazolin <- function() {
  description <- "Two-compartment population PK model for a prophylactic 2 g intravenous bolus of cefazolin in 65 adults undergoing cardiac surgery without cardiopulmonary bypass, with standard extracorporeal circulation (ECC) or with minimally invasive extracorporeal circulation (MiECC) (Santavy 2022). All covariate effects are exponential in the RAW, UNCENTERED covariate (Monolix log-linear form): CKD-EPI eGFR (per mL/s/1.73 m2) raises clearance, body surface area raises the central volume, and either type of bypass circuit lowers the peripheral volume relative to off-pump surgery. The model was fitted on a minute time scale; the published clearances (printed as L/h) are per minute and are converted here to L/h. Proportional residual error."
  reference <- "Santavy P, Sima M, Zuscich O, Kubickova V, Michalickova D, Slanar O, Urbanek K. Population Pharmacokinetics of Prophylactic Cefazolin in Cardiac Surgery with Standard and Minimally Invasive Extracorporeal Circulation. Antibiotics (Basel). 2022;11(11):1582. doi:10.3390/antibiotics11111582"
  vignette <- "Santavy_2022_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine-based CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Santavy 2022 Methods 3.5.2: 'eGFR according to the CKD-EPI formula' was tested as a continuous",
        "covariate. UNITS: the paper reports eGFR in mL/s/1.73 m^2 (Table 1, the Czech SI convention);",
        "the abstract's 'mL/min/1.73 m2' is a misprint, since the Table 1 medians (1.36-1.52) and the",
        "Results 2.3 augmented-renal-clearance boundary (2.17 mL/s) are only physiologic in mL/s. This",
        "column carries the CANONICAL mL/min/1.73 m^2 value and model() divides it by 60 before applying",
        "the published coefficient 0.49, which is per (mL/s/1.73 m^2). A user who supplies mL/s here will",
        "understate clearance by a large factor. The effect is EXPONENTIAL and UNCENTERED: exp(0.49 x",
        "eGFR_SI) is 1 at eGFR = 0 and 2.0 at the cohort median of about 1.42 mL/s/1.73 m^2. Table 1",
        "medians (range): off-pump 1.36 (0.48-1.60), MiECC 1.52 (0.65-2.52), ECC 1.36 (range misprinted",
        "as 1.68-1.85); Results 2.3 gives the study range as 0.67-2.8 mL/s.",
        collapse = " "
      ),
      source_name = "eGFR"
    ),
    BSA = list(
      description = "Body surface area by the DuBois formula",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1 footnote: 'BSA is body surface area estimated according to DuBois formula'. Medians",
        "(range): off-pump 1.98 (1.62-2.46), MiECC 1.99 (1.58-2.37), ECC 2.01 (1.81-2.28) m^2. The effect",
        "on the central volume is EXPONENTIAL and UNCENTERED: exp(0.51 x BSA) is 1 at BSA = 0 and 2.76 at",
        "BSA 1.99 m^2.",
        collapse = " "
      ),
      source_name = "BSA"
    ),
    SURG_CPB_STD = list(
      description = "Cardiac surgery performed on standard (conventional) cardiopulmonary bypass / extracorporeal circulation",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (off-pump surgery, or surgery on minimally invasive extracorporeal circulation)",
      notes = paste(
        "Subject-level, time-fixed indicator of the type of circuit used for the operation, not an",
        "intra-operative phase window (see CPB_ON for that). Santavy 2022 Methods 3.2.1: Stockert S5",
        "roller pump, Terumo FX 25 hard-shell oxygenator, hematocrit close to 25%, hypothermia to 34-35",
        "degrees C for difficult procedures; median CPB duration 124 (70-184) min (Results 2.1). 20 of 65",
        "patients. Mutually exclusive with SURG_CPB_MINI; both 0 = off-pump (20 patients). Enters the",
        "peripheral volume as exp(-0.77 x SURG_CPB_STD), a 54% lower Vd2 than off-pump surgery.",
        collapse = " "
      ),
      source_name = "ECC"
    ),
    SURG_CPB_MINI = list(
      description = "Cardiac surgery performed on minimally invasive extracorporeal circulation (MiECC)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (off-pump surgery, or surgery on standard extracorporeal circulation)",
      notes = paste(
        "Subject-level, time-fixed indicator of the type of circuit used for the operation. Santavy 2022",
        "Methods 3.2.2: Stockert S5 centrifugal pump (Sorin Revolution), Terumo FX 15 oxygenator without",
        "cardiotomy suction, closed system, hematocrit near 28%, normothermia 36-37 degrees C; median CPB",
        "duration 71 (45-110) min (Results 2.1). 25 of 65 patients. Mutually exclusive with SURG_CPB_STD.",
        "Enters the peripheral volume as exp(-0.79 x SURG_CPB_MINI), a 55% lower Vd2 than off-pump",
        "surgery. The Table 2 footnote names this indicator 'MECC'.",
        collapse = " "
      ),
      source_name = "MiECC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight. Tested as a continuous covariate; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 86 (51-117), MiECC 89 (56-117), ECC 85 (71-107) kg. Listed among the continuous covariates in Methods 3.5.2; no point estimate published."
    ),
    HT = list(
      description = "Height. Tested as a continuous covariate; not retained.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 171 (158-194), MiECC 172 (156-196), ECC 176 (168-185) cm. No point estimate published."
    ),
    AGE = list(
      description = "Age. Tested as a continuous covariate; not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 72 (53-86), MiECC 73 (53-85), ECC 68 (44-82) years. No point estimate published."
    ),
    SEXF = list(
      description = "Female sex indicator. Tested as a categorical covariate; not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Results 2.1: 10 of 65 patients female (off-pump 3/20, MiECC 6/25, ECC 1/20). No point estimate published."
    ),
    TBILI = list(
      description = "Serum total bilirubin. Tested as a continuous covariate; not retained.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 12.5 (5.0-50.0), MiECC 11.0 (6.0-51.0), ECC 12.0 (4.0-43.0) umol/L. No point estimate published."
    ),
    ALB = list(
      description = "Serum albumin. Tested as a continuous covariate; not retained.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 45.3 (35.2-50.7), MiECC 45.0 (39.2-51.0), ECC 44.0 (35.1-50.8) g/L. No point estimate published."
    ),
    TPRO = list(
      description = "Total serum protein. Tested as a continuous covariate; not retained.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 70.6 (60.3-81.0), MiECC 73.5 (59.4-84.3), ECC 70.2 (60.7-80.1) g/L. No point estimate published."
    ),
    BUN = list(
      description = "Serum urea. Tested as a continuous covariate; not retained.",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported as urea in mmol/L (Table 1 medians (range): off-pump 4.9 (1.5-12.2), MiECC 4.8 (3.0-8.9), ECC 5.3 (4.1-9.4)). No point estimate published."
    ),
    CREAT = list(
      description = "Serum creatinine. Tested as a continuous covariate; not retained (eGFR was retained instead).",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 medians (range): off-pump 84 (50-188), MiECC 83 (57-130), ECC 90 (64-145) umol/L. No point estimate published."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 1L,
    n_observations = 461L,
    age_range = "44-86 years (group medians 68-73)",
    weight_range = "51-117 kg (group medians 85-89)",
    sex_female_pct = 15.4,
    race_ethnicity = NULL,
    disease_state = "Adults undergoing elective cardiac surgery: 20 off-pump, 25 on minimally invasive extracorporeal circulation (MiECC), 20 on standard extracorporeal circulation (ECC)",
    dose_range = "Single 2 g intravenous bolus of cefazolin 30-60 min before the start of surgery",
    regions = "Czech Republic (single centre: University Hospital Olomouc)",
    renal_function = "CKD-EPI eGFR group medians 1.36-1.52 mL/s/1.73 m^2 (about 82-91 mL/min/1.73 m^2); study range 0.67-2.8 mL/s (Results 2.3). Serum creatinine > 200 umol/L and chronic dialysis were exclusion criteria.",
    notes = paste(
      "Prospective observational open-label PK study (Czech Health Research Council grant 17-31540A).",
      "Plasma sampled at 15, 30, 45, 60, 120 and 180 min after the dose and at the end of surgery; 461",
      "samples, HPLC-UV (Methods 3.3-3.4). Median surgery duration 135 min off-pump, 198 min ECC and",
      "191 min MiECC; median CPB duration 124 min (ECC) and 71 min (MiECC). Exclusion: BMI < 17 or > 35",
      "kg/m^2, liver enzymes > 3 x ULN, cefazolin in the preceding 3 days, extracorporeal elimination",
      "during surgery. Estimation by SAEM in Monolix 2021R1 with time in minutes.",
      collapse = " "
    )
  )

  ini({
    # Structural parameters -- Santavy 2022 Table 2 'Fixed effects'. Table 2
    # footnote (Monolix log-linear covariate model):
    #   Log(CL)  = log(CL_pop)  + beta_CL_eGFR * eGFR + eta_CL
    #   Log(Vd1) = log(Vd1_pop) + beta_Vd1_BSA * BSA  + eta_Vd1
    #   Log(Q)   = log(Q_pop)   + eta_Q
    #   Log(Vd2) = log(Vd2_pop) + beta_Vd2_MECC * MECC + beta_Vd2_ECC * ECC + eta_Vd2
    # The covariates are NOT centred, so CL_pop and Vd1_pop are the values
    # extrapolated to eGFR 0 / BSA 0, not typical-patient values.
    # TIME UNIT: Table 2 labels CL_pop and Q_pop 'L/h', but the fit was in
    # minutes (sampling 15-180 min, PTA horizon 407 min, Figures 1-4 on a
    # minute axis). Read per minute, the model reproduces Table 3's PTA and
    # Figure 4's medians; read per hour, every PTA is 100% and the 407-min
    # median is about 120 mg/L against Figure 4's 7-23 mg/L. Both are
    # therefore multiplied by 60 to give L/h.
    lcl <- log(0.045 * 60); label("Clearance at eGFR 0 (L/h)") # Table 2, CL_pop = 0.045 (R.S.E. 21.7%), per minute -> x 60 = 2.7 L/h
    lvc <- log(4.91); label("Central volume of distribution at BSA 0 (L)") # Table 2, Vd1_pop = 4.91 L (R.S.E. 47.6%)
    lq <- log(0.28 * 60); label("Intercompartmental clearance (L/h)") # Table 2, Q_pop = 0.28 (R.S.E. 17.7%), per minute -> x 60 = 16.8 L/h
    lvp <- log(22.07); label("Peripheral volume of distribution, off-pump surgery (L)") # Table 2, Vd2_pop = 22.07 L (R.S.E. 18.2%)

    # Covariate effects -- exponential in the raw covariate.
    e_crcl_cl <- 0.49; label("Exponential eGFR effect on CL (per mL/s/1.73 m^2)") # Table 2, beta_CL_eGFR = 0.49 (R.S.E. 24.5%)
    e_bsa_vc <- 0.51; label("Exponential BSA effect on Vd1 (per m^2)") # Table 2, beta_Vd1_BSA = 0.51 (R.S.E. 44.8%)
    e_cpbstd_vp <- -0.77; label("Log-scale shift in Vd2 for standard ECC vs off-pump (unitless)") # Table 2, beta_Vd2_ECC = -0.77 (R.S.E. 32.2%)
    e_cpbmini_vp <- -0.79; label("Log-scale shift in Vd2 for MiECC vs off-pump (unitless)") # Table 2, beta_Vd2_MiECC = -0.79 (R.S.E. 29.8%)

    # IIV. Table 2 reports these under 'Standard deviation of the random
    # effects', so each published value is an SD on the log scale and is
    # squared to the variance nlmixr2 expects. No correlations were reported.
    etalcl ~ 0.1225 # Table 2, omega_CL = 0.35 SD -> 0.35^2 = 0.1225 (R.S.E. 16.8%)
    etalvc ~ 0.01 # Table 2, omega_Vd1 = 0.10 SD -> 0.10^2 = 0.01 (R.S.E. 35.7%)
    etalq ~ 0.5929 # Table 2, omega_Q = 0.77 SD -> 0.77^2 = 0.5929 (R.S.E. 15.6%)
    etalvp ~ 0.3364 # Table 2, omega_Vd2 = 0.58 SD -> 0.58^2 = 0.3364 (R.S.E. 15.3%)

    # Residual error. Results 2.2: 'A proportional error model was the most
    # accurate'; Monolix's proportional error parameter b is the SD of the
    # proportional error.
    propSd <- 0.16; label("Proportional residual error (fraction)") # Table 2, 'Proportional' = 0.16 (R.S.E. 4.26%)
  })

  model({
    # CRCL is supplied in mL/min/1.73 m^2; the published coefficient is per
    # mL/s/1.73 m^2.
    crclSi <- CRCL / 60

    cl <- exp(lcl + e_crcl_cl * crclSi + etalcl)
    vc <- exp(lvc + e_bsa_vc * BSA + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + e_cpbstd_vp * SURG_CPB_STD + e_cpbmini_vp * SURG_CPB_MINI + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Results 2.2: two compartments, first-order elimination from the central
    # compartment. The 2 g dose is an IV bolus into `central` (Methods 3.3).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
