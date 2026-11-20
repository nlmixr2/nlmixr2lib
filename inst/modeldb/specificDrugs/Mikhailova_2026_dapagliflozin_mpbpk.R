Mikhailova_2026_dapagliflozin_mpbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based). Dapagliflozin minimal PBPK model",
    "in adults (healthy volunteers and patients with type 2 diabetes)",
    "calibrated by Bayesian inference (Stan/Torsten) on digitised aggregate",
    "plasma profiles from 18 published clinical studies. Oral absorption",
    "through a chain of six transit compartments (rate ktr between transit",
    "states, ka from the last transit state into plasma), two-compartment",
    "disposition with linear clearance, and a satellite renal-filtration",
    "chain (proximal tubule S1/S2 -> S3 -> bladder -> urine) fed by the",
    "unbound filtered load GFR x fu x Cp that predicts urinary recovery.",
    "The renal chain does not deplete the plasma compartment (the total",
    "clearance Cl already includes renal elimination). Typical-value model:",
    "no between-subject variability; exponential residual error. Parameter",
    "values are the Torsten posterior medians under the extended",
    "inverse-variance weighted prior.",
    sep = " "
  )
  reference <- paste(
    "Mikhailova A, Peskov K, Helmlinger G, Sokolov V.",
    "Bayesian Workflow for Minimal PBPK Models: Case Study of Dapagliflozin.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(8):e70307.",
    "doi:10.1002/psp4.70307.",
    "Renal-filtration constants inherited from Sokolov V, Yakovleva T,",
    "Stolbov L, et al. A mechanistic modeling platform of SGLT2 inhibition:",
    "implications for type 1 diabetes. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(6):831-841 (doi:10.1002/psp4.12956), Table S1.",
    sep = " "
  )
  vignette <- "Mikhailova_2026_dapagliflozin_mpbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  dosing <- "transit1"

  paper_specific_compartments <- c("pct", "pst", "bladder")

  covariateData <- list(
    CRCL = list(
      description = "Glomerular filtration rate, BSA-normalised, used as the absolute filtration rate of the renal chain",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converted inside model() to the model's GFR in L/h as",
        "CRCL * 60 / 1000, so the BSA-normalised value is used directly as",
        "the absolute filtration rate (as in the Sokolov SGLT2 renal model",
        "the renal chain is inherited from). Mikhailova 2026 took GFR from",
        "each validation study when reported, otherwise as 0.9 x reported",
        "creatinine clearance, and otherwise used 6.696 L/h/1.73 m^2 (the",
        "mean of the healthy-volunteer GFR range), which is 111.6",
        "mL/min/1.73 m^2 on this column's scale. GFR enters only the",
        "renal-filtration chain; plasma concentrations do not depend on it."
      ),
      source_name = "GFR"
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dapagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dapagliflozin", units = "mg", specimen = "tissue", verified = TRUE),
    pct = list(analyte = "dapagliflozin", units = "mg", specimen = "urine", verified = TRUE),
    pst = list(analyte = "dapagliflozin", units = "mg", specimen = "urine", verified = TRUE),
    bladder = list(analyte = "dapagliflozin", units = "mg", specimen = "urine", verified = TRUE),
    urine = list(analyte = "dapagliflozin", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 578,
    n_studies = 18,
    disease_state = paste(
      "Healthy adult volunteers and adults with type 2 diabetes receiving",
      "dapagliflozin monotherapy; paediatric studies and studies in",
      "populations with additional chronic conditions (e.g. inadequate",
      "glycaemic control, obesity) were excluded."
    ),
    dose_range = "1, 2.5, 5, 10, 20, 25, 50 and 100 mg oral, single dose and once daily (38 study arms)",
    data_type = paste(
      "Aggregate (study-arm mean) plasma concentration-time profiles digitised",
      "from 18 publications identified by a PRISMA systematic search",
      "(PubMed and ClinicalTrials.gov, 2024-11-01); no individual patient data."
    ),
    calibration = paste(
      "Posterior inference on the 9 parallel-design studies; the 9 crossover",
      "studies were held out for external validation of the plasma model,",
      "and urinary recovery from 5 studies validated the renal chain."
    ),
    notes = paste(
      "Mikhailova 2026 Results 'Systematic Review of Dapagliflozin PK Data'.",
      "Age, weight, sex and race distributions of the pooled studies are not",
      "reported."
    )
  )

  ini({
    # Structural PK: Torsten posterior medians, extended inverse-variance
    # weighted prior (Supplement Table S6; the prior used for all reported
    # downstream analyses, Results 'Generalized Model Parameter Inference').
    lktr <- log(11.38); label("Transit rate constant between transit compartments ktr (1/h)") # Table S6 ktr median 11.38
    lka <- log(9.99); label("Absorption rate constant from the sixth transit compartment ka (1/h)") # Table S6 ka median 9.99
    lcl <- log(17.59); label("Apparent clearance from the central compartment Cl (L/h)") # Table S6 Cl median 17.59
    lvc <- log(77.52); label("Apparent central volume of distribution Vd (L)") # Table S6 Vd median 77.52
    lq <- log(8.40); label("Apparent intercompartmental clearance Q (L/h)") # Table S6 Q median 8.40
    lvp <- log(143.25); label("Apparent peripheral volume of distribution Vp (L)") # Table S6 Vp median 143.25

    # Renal-filtration chain: literature constants held constant in
    # Mikhailova 2026 (Figure 2 'Literature' parameter source), inherited
    # from the cited Sokolov 2023 SGLT2 platform, Supplement Table S1.
    fu <- fixed(0.086); label("Dapagliflozin fraction unbound in plasma available for filtration Fupd (fraction)") # Sokolov 2023 Table S1 FUPD = 0.086
    v_pct <- fixed(0.045); label("Proximal tubule S1/S2 segment volume Vs12 (L)") # Sokolov 2023 Table S1 VS12 = 0.045
    v_pst <- fixed(0.019); label("Proximal tubule S3 segment volume Vs3 (L)") # Sokolov 2023 Table S1 VS3 = 0.019
    v_bladder <- fixed(0.2); label("Bladder volume Vbladder (L)") # Sokolov 2023 Table S1 Vbladder = 0.2
    q_pct <- fixed(2.7); label("Tubular flow from S1/S2 to S3 Qlumen (L/h)") # Sokolov 2023 Table S1 Qlumen = 2.7
    q_pst <- fixed(0.72); label("Tubular flow from S3 to the bladder Qbladder (L/h)") # Sokolov 2023 Table S1 Qbladder = 0.72
    q_urine <- fixed(0.055); label("Urine flow from the bladder Qurine (L/h)") # Sokolov 2023 Table S1 Qurine = 0.055

    # Residual error: 'exponential additive' model (Methods), i.e. additive
    # on the log scale.
    expSd <- 0.41; label("Exponential (log-scale additive) residual error sigma (SD, log scale)") # Table S6 sigma median 0.41
  })

  model({
    ktr <- exp(lktr)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # GFR in L/h from the BSA-normalised CRCL column (mL/min/1.73 m^2).
    gfr <- CRCL * 60 / 1000

    # Data S1 equations. The printed dTr5/dt drains Tr5 at ka while dTr6/dt
    # fills Tr6 at ktr * Tr5; Figure 2 and the dTr6/dt input show Tr5 -> Tr6
    # at ktr, with ka acting only from the sixth transit compartment.
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ka * transit6
    d/dt(central) <- ka * transit6 - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # Satellite renal chain (Data S1): the filtered load is not removed
    # from central; Cl is the total clearance.
    d/dt(pct) <- fu * gfr * central / vc - q_pct * pct / v_pct
    d/dt(pst) <- q_pct * pct / v_pct - q_pst * pst / v_pst
    d/dt(bladder) <- q_pst * pst / v_pst - q_urine * bladder / v_bladder
    d/dt(urine) <- q_urine * bladder / v_bladder

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
