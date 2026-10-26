MartinezCasanova_2023_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous piperacillin",
    "(given as piperacillin-tazobactam by continuous infusion after a",
    "4 g loading dose) in critically ill adults with preserved renal",
    "function (CKD-EPI eGFR >= 60 mL/min/1.73 m^2) (Martinez-Casanova",
    "2023). First-order elimination from the central compartment, with",
    "clearance proportional to CKD-EPI estimated creatinine clearance",
    "centered at 99.24 mL/min/1.73 m^2. Between-subject variability on",
    "clearance only and a proportional residual error. Predicts TOTAL",
    "plasma piperacillin; tazobactam was not measured."
  )
  reference <- paste(
    "Martinez-Casanova J, Esteve-Pitarch E, Colom-Codina H,",
    "Gumucio-Sanguino VD, Cobo-Sacristan S, Shaw E, Maisterra-Santos K,",
    "Sabater-Riera J, Perez-Fernandez XL, Rigo-Bonnin R, Tubau-Quintano F,",
    "Carratala J, Padulles-Zamora A (2023).",
    "Predictive Factors of Piperacillin Exposure and the Impact on Target",
    "Attainment after Continuous Infusion Administration to Critically Ill",
    "Patients.",
    "Antibiotics (Basel) 12(3):531.",
    "doi:10.3390/antibiotics12030531.",
    sep = " "
  )
  vignette <- "MartinezCasanova_2023_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "piperacillin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "piperacillin",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CKD-EPI estimated creatinine clearance (glomerular filtration",
        "rate), mL/min/1.73 m^2 -- the canonical BSA-normalized units of",
        "CRCL, so no value transformation is needed. Martinez-Casanova 2023",
        "Results 2.2 compared three estimating equations on CL; CKD-EPI gave",
        "the largest drop in objective function (dMOFV -29.88) ahead of",
        "Cockcroft-Gault (-27.27) and MDRD-4 (-18.62), and only CKD-EPI was",
        "retained. Supply the CKD-EPI value, not Cockcroft-Gault or MDRD-4.",
        "Enters clearance as the through-origin linear ratio CRCL / 99.24",
        "(Table 2 footnote: TVCL = theta x (CLCR_CKD-EPI / 99.24), theta",
        "being CL of the typical patient at 99.24 mL/min/1.73 m^2). Table 1",
        "reports a cohort median of 97.1 (interquartile range 86-114); the",
        "paper does not say which records 99.24 is the median of, and the",
        "printed equation is followed. Covariates were recorded at baseline",
        "and on each sampling occasion (Methods 4.3); the paper does not",
        "state whether CRCL entered as time-varying. Inclusion required",
        "CKD-EPI >= 60 mL/min/1.73 m^2 (Methods 4.2), and the dosing",
        "simulations span 60-200 mL/min/1.73 m^2 (Methods 4.6), so the model",
        "does not describe renal impairment."
      ),
      source_name = "CLCRCKD-EPI"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 4.5.2, power relationship with estimated or allometrically fixed exponent) but not retained; the Discussion states that no size or body-composition metric was significant on any parameter. Cohort median 72 kg (IQR 65-84, Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 26.2 kg/m^2 (IQR 23-29, Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'gender' but not retained. 71 of 106 patients (66.98%) were male (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened but not retained. Cohort median 29 g/L (IQR 26-33, Table 1)."
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained; 63 of 106 patients (59.43%) were ventilated (Table 1). The paper's other screened clinical indicators -- neurocritical status (36, 33.96%), drainage carriage (37, 34.91%) and vasoactive drugs (45, 42.45%) -- have no canonical covariate column and are listed in population$notes. Vasoactive drugs were significant on their own (dMOFV -8.83) but not once a renal-function covariate was in the model (Results 2.2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 106L,
    n_studies = 1L,
    n_observations = 218L,
    age_range = "Adults (>= 18 years by inclusion criterion); median 65 years, IQR 50-72",
    age_median = "65 years",
    weight_range = "65-84 kg (interquartile range)",
    weight_median = "72 kg",
    sex_female_pct = 33.0,
    race_ethnicity = "Caucasian 96.23% (102 of 106); remainder not reported.",
    disease_state = paste(
      "Critically ill adults in intensive care treated with",
      "piperacillin-tazobactam by continuous infusion, with preserved renal",
      "function (CKD-EPI estimated creatinine clearance >= 60",
      "mL/min/1.73 m^2). Patients on renal replacement therapy or",
      "extracorporeal membrane oxygenation, pregnant women and patients",
      "under 18 were excluded. Median CKD-EPI 97.1 mL/min/1.73 m^2 (IQR",
      "86-114); augmented renal clearance (CKD-EPI >= 130) in 11 patients",
      "(10.38%). Median albumin 29 g/L. Neurocritical status 33.96%,",
      "drainage carriers 34.91%, mechanically ventilated 59.43%, on",
      "vasoactive drugs 42.45% (Table 1)."
    ),
    renal_function = "CKD-EPI 97.1 (86-114), Cockcroft-Gault 109.4 mL/min (84-143), MDRD-4 116.3 mL/min/1.73 m^2 (89-146); median (IQR), Table 1.",
    dose_range = paste(
      "Loading dose of piperacillin-tazobactam 4/0.5 g over 30 min on day",
      "1, then a continuous 24 h infusion of piperacillin 500 mg/h (12 g",
      "piperacillin / 1.5 g tazobactam per day). Patients already on",
      "intermittent piperacillin-tazobactam for >= 24 h received no loading",
      "dose. Doses were adjusted by therapeutic drug monitoring, but only on",
      "four occasions in three patients (Methods 4.3, Results 2.1)."
    ),
    regions = "Spain (Bellvitge University Hospital ICU, Barcelona).",
    notes = paste(
      "Observational prospective single-centre study, August 2015 to June",
      "2019. 218 total plasma piperacillin concentrations (2.1 per patient)",
      "drawn 30-60 min after the loading dose and at steady state (from 6 h",
      "after the start of the continuous infusion); geometric mean Cmax 112",
      "mg/L and Css 39.5 mg/L (Results 2.1). UHPLC-MS/MS, LLOQ 0.54 mg/L.",
      "Fit in NONMEM 7.4 with FOCE-I; 200-sample bootstrap (Table 2). Screened",
      "covariates without a canonical column: neurocritical status, drainage",
      "carriage and vasoactive drugs (none retained). The paper converts",
      "total to unbound concentration for its target-attainment simulations",
      "with 30% protein binding, i.e. fu = 0.7 (Methods 4.6)."
    )
  )

  ini({
    # Martinez-Casanova 2023 Table 2, 'Final Model Parameter Estimates
    # (RSE%) [SHR%]'.
    lcl <- log(12.0)
    label("Clearance at a CKD-EPI creatinine clearance of 99.24 mL/min/1.73 m^2 (L/h)")
    # Table 2 row 'CL (L/h) = 12.0 (6.03) x (CLCRCKD-EPI/99.24)'; bootstrap
    # median 11.95 (95% CI 10.36-13.42)

    lvc <- log(20.7)
    label("Central volume of distribution (L)")
    # Table 2 row 'Vc (L) = 20.7 (8.94)'; bootstrap 20.29 (9.20-25.05)

    lvp <- log(62.4)
    label("Peripheral volume of distribution (L)")
    # Table 2 row 'Vp (L) = 62.4 (50.80)'; bootstrap 66.32 (15.45-181.76)

    lq <- log(4.77)
    label("Intercompartmental (distributional) clearance (L/h)")
    # Table 2 row 'CLD (L/h) = 4.77 (26.4)'; the abstract rounds it to 4.8.
    # Bootstrap 5.16 (3.19-52.05)

    # Between-patient variability on CL only (Results 2.2), modelled
    # exponentially (Methods 4.5.1). Table 2 reports the variance directly.
    etalcl ~ 0.190
    # Table 2 row 'omega2 CL = 0.190 (24.42) [17.16]'; bootstrap 0.181

    # Proportional residual error (Results 2.2). Table 2 reports the NONMEM
    # variance sigma^2 = 0.140, so the SD is sqrt(0.140) = 0.3742.
    propSd <- 0.3742
    label("Proportional residual error (fraction)")
    # Table 2 row 'sigma2 = 0.140 (15.7) [13.45]'; propSd = sqrt(0.140)
  })

  model({
    # Table 2 footnote: TVCL = theta x (CLCR_CKD-EPI / 99.24), a linear
    # through-origin relationship with no estimated exponent.
    cl <- exp(lcl + etalcl) * (CRCL / 99.24)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two compartments with first-order elimination from central (Results
    # 2.2). Doses are intravenous infusions into central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L. TOTAL plasma piperacillin; the paper's unbound
    # concentration is 0.7 * Cc (Methods 4.6).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
