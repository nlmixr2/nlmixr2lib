Setiawan_2022_levofloxacin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for levofloxacin in hospitalised adult",
    "patients (ICU and non-ICU wards) in Surabaya, Indonesia, most of them treated for",
    "pneumonia. Fitted non-parametrically with the NPAG algorithm in Pmetrics 1.9.",
    "Clearance is a linear function of the CKD-EPI estimated glomerular filtration rate,",
    "CL = 0.044 * eGFR + 0.358 L/h (Results); central volume, intercompartmental clearance",
    "and peripheral volume carry no covariates. Inter-individual variability is a log-normal",
    "approximation to the published NPAG CV%. The unbound concentration Cu = 0.7 * Cc is",
    "exposed using the fixed 30% protein binding the authors applied in their fAUC0-24/MIC",
    "Monte Carlo simulations. Residual unexplained variability is carried as fixed(0) because",
    "the Pmetrics assay-error polynomial and the selected lambda/gamma term were never",
    "published.",
    sep = " "
  )
  reference <- paste(
    "Setiawan E, Abdul-Aziz MH, Cotta MO, Susaniwati S, Cahjono H, Sari IY, Wibowo T,",
    "Marpaung FR, Roberts JA. Population pharmacokinetics and dose optimization of",
    "intravenous levofloxacin in hospitalized adult patients. Sci Rep. 2022;12:8930.",
    "doi:10.1038/s41598-022-12627-1. PMCID: PMC9142570.",
    sep = " "
  )
  vignette <- "Setiawan_2022_levofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI creatinine equation, BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The ONLY covariate retained in the final model. Enters clearance as the additive",
        "linear relationship CL = 0.044 * eGFR + 0.358 (Setiawan 2022 Results, 'Population PK",
        "model'), i.e. a 0.358 L/h intercept plus 0.044 L/h per mL/min/1.73 m^2; there is no",
        "centring value. eGFR was computed with the CKD-EPI equation from serum creatinine",
        "recorded on the day of recruitment (Methods). Cohort mean 52.7 +/- 33.7 mL/min/1.73 m^2",
        "(Table 1). The authors simulated eGFR 20, 50, 80 and 120 mL/min/1.73 m^2; patients on",
        "renal replacement therapy were excluded, so the relationship does not apply to them.",
        sep = " "
      ),
      source_name = "eGFR CKD-EPI"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened as a candidate covariate on Vd and CL (Methods, 'gender') but not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Screened as a candidate covariate on Vd and CL (Methods) but not retained; eGFR",
        "CKD-EPI, which is derived from it, was retained instead. Cohort mean 1.99 +/- 1.48",
        "mg/dL (Table 1).",
        sep = " "
      )
    ),
    DIS_CRITILL = list(
      description = "Hospitalisation type: intensive care unit (1) versus non-ICU ward (0)",
      units = "(binary)",
      type = "categorical",
      notes = "Screened as a candidate covariate ('hospitalisation type (ICU and non-ICU)', Methods) but not retained; 6 of 26 patients were in the ICU (Table 1)."
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Screened as a candidate covariate (Methods) but not retained; 4 of 26 patients, all in the ICU, were ventilated (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    n_concentrations = 121L,
    age_range = "Adults >= 18 years; mean 58.8 +/- 16.4 years",
    weight_range = "Mean 61.6 +/- 12.1 kg (weight recorded in 11 of 26 patients)",
    sex_female_pct = 38.5,
    race_ethnicity = "Indonesian (two hospitals in Surabaya, East Java)",
    disease_state = paste(
      "Hospitalised adults receiving intravenous levofloxacin, 6 in the intensive care unit",
      "(4 mechanically ventilated) and 20 on general wards; 88.5% were diagnosed with",
      "pneumonia. Patients planned for renal replacement therapy or ECMO at the time of",
      "sampling, and pregnant women, were excluded.",
      sep = " "
    ),
    renal_function = "Serum creatinine mean 1.99 +/- 1.48 mg/dL; eGFR CKD-EPI mean 52.7 +/- 33.7 mL/min/1.73 m^2",
    dose_range = paste(
      "500 mg or 750 mg once daily as a 30-minute intravenous infusion (500 mg: 8 patients;",
      "750 mg: 14; 500 then 750 mg: 1; 750 then 500 mg: 3). No patient received more than",
      "750 mg in 24 h.",
      sep = " "
    ),
    regions = "Indonesia (Dr. Mohamad Soewandhie Public Hospital and PHC Hospital, Surabaya)",
    notes = paste(
      "Baseline demographics from Setiawan 2022 Table 1. Up to six samples per patient within",
      "one dosing interval (2-6 per patient, mean 4.77); 5 implausible and 3 contaminated",
      "concentrations were removed, leaving 121 concentrations. Total plasma levofloxacin by",
      "UHPLC-MS/MS, linear 0.1-50 mg/L (LLOQ 0.1 mg/L). Study period November 2018 to",
      "November 2019.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Clearance. Setiawan 2022 Results ('Population PK model'): "The CL of
    # levofloxacin was best described as the following equation:
    # CL = (0.044 * eGFR CKD-EPI) + 0.358". Encoded verbatim as an additive
    # linear relationship (precedent Conil_2007_ceftazidime.R). Table 2's CL row
    # (mean 1.12, median 0.90 L/h) cannot be the clearance of a typical patient
    # under this equation (it would need eGFR ~17); see vignette Errata, where
    # the paper's own PTA curves (Figures 3-4) confirm the equation as the
    # typical clearance.
    # ------------------------------------------------------------------------
    lcl <- log(0.358)
    label("Intercept of the linear CL ~ eGFR relationship (L/h)") # Results: CL = (0.044*eGFR) + 0.358
    e_crcl_cl <- 0.044
    label("eGFR slope on CL (L/h per mL/min/1.73 m^2)") # Results: CL = (0.044*eGFR) + 0.358

    # Volumes and intercompartmental clearance: Table 2 MEAN of the NPAG
    # distribution (the values quoted in the Abstract; same convention as the
    # sibling Setiawan_2023_sulbactam.R).
    lvc <- log(27.6)
    label("Central volume of distribution (L)") # Table 2 Vc mean 27.6 L (SD 19.1, CV 69.3%, median 26.4)
    lq <- log(30.9)
    label("Intercompartmental clearance (L/h)") # Table 2 Q mean 30.9 L/h (SD 16.4, CV 53.2%, median 33.2)
    lvp <- log(28.2)
    label("Peripheral volume of distribution (L)") # Table 2 Vp mean 28.2 L (SD 16.2, CV 57.7%, median 27.9)

    # ------------------------------------------------------------------------
    # Inter-individual variability. NPAG estimates a discrete non-parametric
    # distribution; Methods states "The %CV was used to describe
    # inter-individual PK variability", so the Table 2 CV% is carried as a
    # log-normal approximation, omega^2 = log(CV^2 + 1). No covariances are
    # reported, so the etas are independent.
    # ------------------------------------------------------------------------
    etalcl ~ 0.239332 # Table 2 CL CV% = 52 -> log(0.52^2 + 1)
    etalvc ~ 0.392210 # Table 2 Vc CV% = 69.3 -> log(0.693^2 + 1)
    etalq ~ 0.249220 # Table 2 Q CV% = 53.2 -> log(0.532^2 + 1)
    etalvp ~ 0.287379 # Table 2 Vp CV% = 57.7 -> log(0.577^2 + 1)

    # Protein binding: a literature value the authors applied in the dosing
    # simulations, not an estimate.
    fu <- fixed(0.7)
    label("Fraction of levofloxacin unbound in plasma (unitless)") # Methods 'Dosing simulations': protein binding set at 30%; 1 - 0.30 = 0.70

    # ------------------------------------------------------------------------
    # Residual unexplained variability is NOT reported. Methods states only
    # that "Both lambda (ranging from 0.1 to 0.9) and gamma (ranging from 1 to
    # 9) error models were tested"; neither the selected error model nor the
    # Pmetrics assay-error polynomial appears in the paper or its supplement.
    # Carried as fixed(0) rather than invented -- see vignette Errata.
    # ------------------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })

  model({
    # 1. Individual parameters. eGFR CKD-EPI (CRCL, mL/min/1.73 m^2) acts on
    #    clearance only, additively on the linear scale (Results equation), with
    #    log-normal IIV around the covariate-predicted value.
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment ODE system; dosing is a 30-minute intravenous infusion
    #    into central (Methods).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. The assay measured TOTAL plasma levofloxacin, so Cc is the
    #    total concentration and Cu the unbound concentration behind the
    #    fAUC0-24/MIC >= 80 target.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd)
  })
}
