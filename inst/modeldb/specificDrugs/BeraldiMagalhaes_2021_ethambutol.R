BeraldiMagalhaes_2021_ethambutol <- function() {
  description <- paste(
    "Two-compartment oral population PK model of ethambutol given as a crushed or whole",
    "rifampin/isoniazid/pyrazinamide/ethambutol fixed-dose-combination (FDC) tablet to 30 Brazilian",
    "adults with tuberculosis: 10 mechanically ventilated ICU patients (tablet crushed and given by",
    "nasogastric tube) and 20 outpatients (tablet swallowed whole). Fitted non-parametrically with",
    "NPAG in Pmetrics 1.5.0. Every structural parameter -- clearance, central volume, first- and",
    "second-occasion absorption rate constants, bioavailability, intercompartmental clearance and",
    "peripheral volume -- was estimated separately for the ICU and outpatient groups inside one",
    "model, so each carries an _icu / _outpt stratum suffix selected by DIS_CRITILL, with its own",
    "log-normal IIV approximated from the Table 3 %CV. Creatinine clearance normalised to",
    "101 mL/min scales clearance and total body weight normalised to 56 kg scales both volumes with",
    "a 0.25 allometric exponent. Residual error is fixed(0) because the Pmetrics error-model",
    "coefficients are not published. The published ICU parameters do not reproduce the paper's own",
    "observed ICU exposures; see the vignette.",
    sep = " "
  )
  reference <- paste(
    "Beraldi-Magalhaes F, Parker SL, Sanches C, Sousa Garcia L, Souza Carvalho BK, Fachi MM,",
    "de Liz MV, Pontarolo R, Lipman J, Cordeiro-Santos M, Roberts JA.",
    "Is Dosing of Ethambutol as Part of a Fixed-Dose Combination Product Optimal for Mechanically",
    "Ventilated ICU Patients with Tuberculosis? A Population Pharmacokinetic Study.",
    "Antibiotics (Basel). 2021;10(12):1559. doi:10.3390/antibiotics10121559. PMCID: PMC8698281.",
    sep = " "
  )
  vignette <- "BeraldiMagalhaes_2021_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DIS_CRITILL = list(
      description = paste(
        "Study group: 1 = mechanically ventilated ICU patient, 0 = outpatient. Selects the",
        "_icu or _outpt parameter set."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Results 2.2: 'All pharmacokinetic parameters for ICU patients were significantly different",
        "from outpatients and were estimated separately in the model using selective execution",
        "statements.' The ICU group differs from the outpatients in more than critical illness: all",
        "10 were mechanically ventilated and received the FDC tablet crushed, suspended in 20 mL of",
        "water and given through a nasogastric tube, whereas outpatients swallowed whole tablets",
        "(Methods 4.2). The indicator therefore carries route and formulation handling as well as",
        "illness, and the authors attribute the ka and F differences to administration (Discussion)."
      ),
      source_name = "ICU / outpatient group"
    ),
    CRCL = list(
      description = "Measured creatinine clearance, NOT BSA-normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods 4.2: 'A measured 8 h creatinine clearance was obtained.' Table 1 reports it in",
        "mL/min (median 92.3 in the ICU and 113.88 in outpatients). Results 2.2: 'creatinine",
        "clearance normalized to 101 mL/min on clearance'. The functional form is not printed; the",
        "through-origin ratio (CRCL/101) is used because it is what 'normalized to' denotes and no",
        "exponent was reported. Methods 4.6 labels the simulation grid 'mL/min/1.73 m2', which",
        "conflicts with Table 1; the Table 1 unit (raw mL/min) is used here."
      ),
      source_name = "Clcr"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results 2.2: 'total body weight normalized to 56 kg as a covariate with an allometric",
        "scaler (raised to the 25th power) on central and peripheral volume of distribution'.",
        "Read as (WT/56)^0.25; see e_wt_vc_vp. Also sets the FDC dose band (20-35 kg 550 mg,",
        "36-50 kg 825 mg, > 50 kg 1100 mg ethambutol; Methods 4.2)."
      ),
      source_name = "TBW"
    ),
    OCC = list(
      description = paste(
        "Sampled dosing occasion: 1 = first sampled dose (enrolment day 1), 2 = second sampled",
        "dose (enrolment day 3). Selects ka for occasion 1 or 2."
      ),
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Methods 4.4: 'the inclusion of occasion for the first and second dose was tested for the",
        "rate of absorption, bioavailability, lag time and clearance'; only absorption was retained",
        "(Table 3 'Ka1' and 'Ka2', 'absorption rate constant for the 1st and 2nd dose'). OCC >= 2",
        "uses the second-occasion value. Samples were taken on enrolment days 1 and 3 (Methods 4.2)."
      ),
      source_name = "occasion"
    )
  )

  covariatesDataExcluded <- list(
    HIV = list(
      description = "HIV co-infection status (and viral load, CD4 count, antiretroviral therapy)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened but NOT retained. Results 2.2: 'HIV status, HIV viral load, CD4 cell count and",
        "the antiretroviral therapy used by the subjects were tested as covariates ... but a linear",
        "regression returned a correlation coefficient < 0.2 and did not improve the",
        "pharmacokinetic model'. 9/10 ICU and 15/20 outpatients were HIV-positive (Table 1)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = paste(
        "Reported in Table 1 (median 3.1 ICU, 3.5 outpatients) and screened among the 'renal and",
        "liver function' covariates (Methods 4.4), but not retained."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "Adults >= 18 years; median 31.0 (IQR 29-40) ICU, 39.5 (IQR 32.7-46.2) outpatients",
    weight_range = "Median 51.2 kg (IQR 46.2-58.6) ICU, 58.35 kg (IQR 53.2-67) outpatients",
    sex_female_pct = 20,
    race_ethnicity = "Not reported (Amazonas State, Brazil)",
    disease_state = paste(
      "Active pulmonary or extrapulmonary tuberculosis. 10 mechanically ventilated ICU patients",
      "(median SOFA 10, APACHE II 20.5, 8/10 on vasoactive drugs) and 20 outpatients. HIV",
      "co-infection in 9/10 ICU patients and 15/20 outpatients. Measured creatinine clearance",
      "median 92.3 mL/min (ICU) and 113.88 mL/min (outpatients); patients on any form of dialysis",
      "or renal replacement therapy were excluded."
    ),
    dose_range = paste(
      "Once-daily weight-banded FDC tablets (rifampin/isoniazid/pyrazinamide/ethambutol) per the",
      "Brazilian Ministry of Health guideline: 550 mg (20-35 kg), 825 mg (36-50 kg) or 1100 mg",
      "(> 50 kg) ethambutol. ICU: crushed, suspended in 20 mL water, via nasogastric tube.",
      "Outpatients: whole tablets orally, directly observed."
    ),
    regions = "Brazil (Fundacao de Medicina Tropical Dr. Heitor Vieira Dourado, Manaus, Amazonas)",
    notes = paste(
      "Beraldi-Magalhaes 2021 Table 1. 352 plasma concentrations; samples pre-dose and at 0.5, 1,",
      "2, 4, 6, 8, 12 and 24 h on enrolment days 1 and 3, after a median 10 (ICU) or 11",
      "(outpatients) days of treatment. Total ethambutol by LC-MS/MS (0.2-5 mg/L with a dilution",
      "QC). The pre-dose concentration of each subject was used as the model initial condition",
      "during fitting; that data-fitting device is not encoded here (see the vignette)."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # All structural values are the MEAN of the NPAG non-parametric marginal
    # distribution in Beraldi-Magalhaes 2021 Table 3 ('Estimates of ethambutol
    # pharmacokinetic parameters for the final covariate model'), reported
    # separately for the ICU (n = 10) and outpatient (n = 20) groups. The
    # median is recorded in each comment. The mean is used for consistency
    # with the other Pmetrics/NPAG extractions (Setiawan_2023_sulbactam,
    # Bustinduy_2016_praziquantel, Braune_2018_meropenem).
    #
    # Table 3 footnote: 'Clearance, relative clearance; Volume, relative volume
    # of distribution of central compartment'. F was estimated alongside, so
    # CL and V are the values that the estimated F is applied against.
    # ------------------------------------------------------------------------

    # ---- ICU patients (DIS_CRITILL = 1) ----
    lcl_icu <- log(1.2); label("Clearance at CRCL = 101 mL/min, ICU patients (L/h)")
    # Table 3 ICU 'Clearance (L/h)' mean 1.2 (SD 1.5, median 0.9, CV 120.9%)
    lvc_icu <- log(64.8); label("Central volume at WT = 56 kg, ICU patients (L)")
    # Table 3 ICU 'Volume (L)' mean 64.8 (SD 11.7, median 61.1, CV 18.1%)
    lka_occ1_icu <- log(0.72); label("Absorption rate constant, first sampled dose, ICU patients (1/h)")
    # Table 3 ICU 'Ka1 (h-1)' mean 0.72 (SD 0.05, median 0.7, CV 7.4%)
    lka_occ2_icu <- log(0.75); label("Absorption rate constant, second sampled dose, ICU patients (1/h)")
    # Table 3 ICU 'Ka2 (h-1)' mean 0.75 (SD 0.10, median 0.8, CV 13.8%)
    lfdepot_icu <- log(0.80); label("Bioavailability, ICU patients (crushed FDC via nasogastric tube) (fraction)")
    # Table 3 ICU 'F' mean 0.80 (SD 0.06, median 0.8, CV 7.9%)
    lq_icu <- log(7.3); label("Intercompartmental clearance, ICU patients (L/h)")
    # Table 3 ICU 'Q(L/h)' mean 7.3 (SD 3.5, median 6.6, CV 48.6%)
    lvp_icu <- log(348.6); label("Peripheral volume at WT = 56 kg, ICU patients (L)")
    # Table 3 ICU 'Vp (L)' mean 348.6 (SD 30.1, median 361.6, CV 8.6%)

    # ---- Outpatients (DIS_CRITILL = 0) ----
    lcl_outpt <- log(17.5); label("Clearance at CRCL = 101 mL/min, outpatients (L/h)")
    # Table 3 outpatients 'Clearance (L/h)' mean 17.5 (SD 13.3, median 11.1, CV 75.8%)
    lvc_outpt <- log(137.2); label("Central volume at WT = 56 kg, outpatients (L)")
    # Table 3 outpatients 'Volume (L)' mean 137.2 (SD 55.1, median 170.4, CV 40.1%)
    lka_occ1_outpt <- log(0.35); label("Absorption rate constant, first sampled dose, outpatients (1/h)")
    # Table 3 outpatients 'Ka1 (h-1)' mean 0.35 (SD 0.12, median 0.3, CV 35.5%)
    lka_occ2_outpt <- log(0.39); label("Absorption rate constant, second sampled dose, outpatients (1/h)")
    # Table 3 outpatients 'Ka2 (h-1)' mean 0.39 (SD 0.18, median 0.4, CV 44.9%)
    lfdepot_outpt <- log(0.14); label("Bioavailability, outpatients (whole FDC tablet orally) (fraction)")
    # Table 3 outpatients 'F' mean 0.14 (SD 0.13, median 0.1, CV 87.1%)
    lq_outpt <- log(2.66); label("Intercompartmental clearance, outpatients (L/h)")
    # Table 3 outpatients 'Q(L/h)' mean 2.66 (SD 2.02, median 3.8, CV 75.9%)
    lvp_outpt <- log(343.3); label("Peripheral volume at WT = 56 kg, outpatients (L)")
    # Table 3 outpatients 'Vp (L)' mean 343.3 (SD 78.2, median 400.0, CV 22.8%)

    # ---- Covariate effect shared by both groups ----
    e_wt_vc_vp <- fixed(0.25); label("Allometric exponent of (WT/56) on central and peripheral volume (unitless)")
    # Results 2.2: 'total body weight normalized to 56 kg as a covariate with an
    # allometric scaler (raised to the 25th power) on central and peripheral
    # volume of distribution'. '25th power' read as 0.25 (a literal exponent of
    # 25 is not a plausible allometric scaler); no estimate or uncertainty is
    # reported, so it is carried as fixed.

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution, not an omega. Table 3 prints a %CV per parameter per group
    # that equals SD/mean on the linear scale (all 14 rows reproduce within
    # rounding), carried here as a log-normal approximation
    # omega^2 = log(CV^2 + 1). Each group has its own variance because the two
    # groups were estimated separately. No covariances are reported.
    # ------------------------------------------------------------------------
    etalcl_icu ~ 0.900844 # Table 3 ICU Clearance CV 120.9%
    etalvc_icu ~ 0.032236 # Table 3 ICU Volume CV 18.1%
    etalka_occ1_icu ~ 0.005461 # Table 3 ICU Ka1 CV 7.4%
    etalka_occ2_icu ~ 0.018865 # Table 3 ICU Ka2 CV 13.8%
    etalfdepot_icu ~ 0.006222 # Table 3 ICU F CV 7.9%
    etalq_icu ~ 0.212039 # Table 3 ICU Q CV 48.6%
    etalvp_icu ~ 0.007369 # Table 3 ICU Vp CV 8.6%

    etalcl_outpt ~ 0.453978 # Table 3 outpatients Clearance CV 75.8%
    etalvc_outpt ~ 0.149110 # Table 3 outpatients Volume CV 40.1%
    etalka_occ1_outpt ~ 0.118694 # Table 3 outpatients Ka1 CV 35.5%
    etalka_occ2_outpt ~ 0.183655 # Table 3 outpatients Ka2 CV 44.9%
    etalfdepot_outpt ~ 0.564541 # Table 3 outpatients F CV 87.1%
    etalq_outpt ~ 0.454941 # Table 3 outpatients Q CV 75.9%
    etalvp_outpt ~ 0.050678 # Table 3 outpatients Vp CV 22.8%

    # ------------------------------------------------------------------------
    # Residual error. Methods 4.4 describes the Pmetrics assay-error polynomial
    # SD = C0 + C1*Y with an additive (lambda) or multiplicative (gamma) term,
    # but neither the C0/C1 coefficients nor the chosen lambda/gamma are
    # reported anywhere in the paper, and the supplementary files are the
    # figure images only. Carried as fixed(0) rather than invented.
    # ------------------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
    addSd <- fixed(0); label("Additive residual SD (mg/L; not reported in the source)")
  })

  model({
    # 1. Group and occasion selectors. DIS_CRITILL = 1 selects the ICU
    #    parameter set; OCC = 1 selects the first-occasion absorption rate.
    icu <- DIS_CRITILL
    occ1 <- (OCC <= 1)

    # 2. Group-specific individual parameters (Table 3), each with its own IIV.
    cl_icu <- exp(lcl_icu + etalcl_icu)
    vc_icu <- exp(lvc_icu + etalvc_icu)
    ka_occ1_icu <- exp(lka_occ1_icu + etalka_occ1_icu)
    ka_occ2_icu <- exp(lka_occ2_icu + etalka_occ2_icu)
    fdepot_icu <- exp(lfdepot_icu + etalfdepot_icu)
    q_icu <- exp(lq_icu + etalq_icu)
    vp_icu <- exp(lvp_icu + etalvp_icu)

    cl_outpt <- exp(lcl_outpt + etalcl_outpt)
    vc_outpt <- exp(lvc_outpt + etalvc_outpt)
    ka_occ1_outpt <- exp(lka_occ1_outpt + etalka_occ1_outpt)
    ka_occ2_outpt <- exp(lka_occ2_outpt + etalka_occ2_outpt)
    fdepot_outpt <- exp(lfdepot_outpt + etalfdepot_outpt)
    q_outpt <- exp(lq_outpt + etalq_outpt)
    vp_outpt <- exp(lvp_outpt + etalvp_outpt)

    # 3. Covariates (Results 2.2): creatinine clearance normalised to 101 mL/min
    #    on clearance; total body weight normalised to 56 kg with a 0.25
    #    allometric exponent on both volumes.
    cl <- (icu * cl_icu + (1 - icu) * cl_outpt) * (CRCL / 101)
    vc <- (icu * vc_icu + (1 - icu) * vc_outpt) * (WT / 56)^e_wt_vc_vp
    vp <- (icu * vp_icu + (1 - icu) * vp_outpt) * (WT / 56)^e_wt_vc_vp
    q <- icu * q_icu + (1 - icu) * q_outpt
    ka_icu <- occ1 * ka_occ1_icu + (1 - occ1) * ka_occ2_icu
    ka_outpt <- occ1 * ka_occ1_outpt + (1 - occ1) * ka_occ2_outpt
    ka <- icu * ka_icu + (1 - icu) * ka_outpt
    fdepot <- icu * fdepot_icu + (1 - icu) * fdepot_outpt

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. Two-compartment disposition with first-order absorption and linear
    #    elimination from the central compartment (Results 2.2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # 6. Observation: total plasma ethambutol (Methods 4.3).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
