Royer_2021_palbociclib <- function() {
  description <- "One-compartment population PK model with first-order absorption, a fixed absorption lag time and first-order elimination for oral palbociclib in women with metastatic breast cancer followed in routine care (therapeutic drug monitoring, mostly one sample per patient). Apparent oral clearance increases with Cockcroft-Gault creatinine clearance through a power model centred on the cohort mean of 78.9 mL/min (exponent 0.419). Correlated inter-individual variability on CL/F and ka; combined additive + proportional residual error."
  reference <- paste(
    "Royer B, Kaderbhai C, Fumet JD, Hennequin A, Desmoulins I, Ladoire S, Ayati S, Mayeur D, Ilie S, Schmitt A. (2021).",
    "Population Pharmacokinetics of Palbociclib in a Real-World Situation.",
    "Pharmaceuticals 14(3):181. doi:10.3390/ph14030181.",
    sep = " "
  )
  vignette <- "Royer_2021_palbociclib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")
  # Unit note. Table 2 reports CL/F in L/h and V/F in L; doses are in mg
  # (Results: 75, 100 and 125 mg per day), so central / vc is in mg/L.
  # Concentrations are reported in ug/L throughout (LLOQ 6 ug/L, Section 4.2;
  # additive error in ug/L, Table 2), hence the 1000 ug/L per mg/L factor on Cc.

  compartmentData <- list(
    depot = list(analyte = "palbociclib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "palbociclib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault formula (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Royer 2021 Table 1 row 'Creatinine Clearance - Cockroft & Gault (CRCL - ml/min)': mean 78.9, median 72.1, range 23.4-282.3 mL/min. No BSA normalization is described, so values are raw mL/min (the CRCL register entry accepts raw Cockcroft-Gault when the source does not BSA-normalize). Section 4.3 states continuous covariates entered as a power model 'normalized using the mean value' of the population, so the centring value is the Table 1 MEAN 78.9 mL/min, not the median.",
      source_name = "CRCL"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 2: significant in the forward step on CL/F (dOFV -13.5) but removed in backward deletion once CRCL was retained. Table 1 mean 69.7 kg, range 37.0-140.0 kg.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 2: significant in the forward step on CL/F (dOFV -14.6) but removed in backward deletion once CRCL was retained. Table 1 mean 67.4 years, range 40.7-92.2.",
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Section 2: dOFV -7.74 on CL/F but only a 2.2% reduction in IIV CL/F, so not carried to backward deletion. Table 1 mean 74.6 umol/L.",
      source_name = "Serum creatinine"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 124L,
    n_studies = 1L,
    age_range = "40.7-92.2 years (mean 67.4, median 70.4)",
    weight_range = "37.0-140.0 kg (mean 69.7)",
    sex_female_pct = 100,
    disease_state = "Women receiving palbociclib for breast cancer as part of its indication (hormone-receptor-positive, HER2-negative metastatic breast cancer), combined with an aromatase inhibitor (letrozole or anastrozole) or fulvestrant.",
    dose_range = "Oral palbociclib 75 mg (n = 15), 100 mg (n = 33) or 125 mg (n = 101) once daily, 21 days on / 7 days off (28-day cycle).",
    regions = "France (Centre Georges-Francois Leclerc, Dijon)",
    observations = "151 plasma concentrations from 124 patients (routine therapeutic drug monitoring, 28 October 2018 to 3 August 2020). 27 patients had two samples; for 25 of them the samples came from different cycles and were modelled as independent individuals (no inter-occasion variability). Samples drawn 1-21 days after the start of the cycle (30 in the first 8 days) and 0.9-197.25 h after the previous dose; concentrations 6-226 ug/L (mean 81.8), none below the 6 ug/L LLOQ.",
    notes = "Baseline characteristics from Royer 2021 Table 1: serum creatinine mean 74.6 umol/L (31.0-301.8); Cockcroft-Gault CRCL mean 78.9 mL/min (median 72.1, 23.4-282.3); albumin mean 39.9 g/L. Race / ethnicity was not reported. Table 1 prints a body-weight median of 98.0 kg against a mean of 69.7 kg, which is not plausible and is probably a typesetting error; weight is not used by the model."
  )

  ini({
    # Structural parameters: Royer 2021 Table 2 (final model, NONMEM 7.4 FOCE-I).
    lka <- log(0.187); label("ka: first-order absorption rate constant (1/h)") # Table 2: Ka = 0.187 /h (RSE 19.3%); bootstrap 0.107-0.370
    lcl <- log(58.3); label("CL/F: apparent oral clearance at CRCL = 78.9 mL/min (L/h)") # Table 2: CL/F = 58.3 L/h (RSE 3.3%); bootstrap 54.2-62.8
    lvc <- log(1580); label("V/F: apparent volume of distribution (L)") # Table 2: V/F = 1580 L (RSE 16.2%); bootstrap 930-2568
    # Tlag was fixed to the value of the earlier palbociclib popPK analysis
    # (Sun and Wang 2014, reference 7) -- Section 2: 'We then decided to fix
    # the Tlag to the already published value'. Table 2 prints no RSE for it
    # and mislabels its unit as '(L)'; it is a time in hours.
    ltlag <- fixed(log(0.658)); label("Tlag: absorption lag time (h)") # Table 2: Tlag = 0.658 (RSE '-', fixed per Section 2)

    # Covariate effect, Section 4.3 power model normalized by the population
    # mean: CLi = CLpop * (CRCL / 78.9)^theta_CRCL (Table 1 mean CRCL 78.9 mL/min).
    e_crcl_cl <- 0.419; label("Exponent of the (CRCL / 78.9 mL/min) power effect on CL/F (unitless)") # Table 2: CRCL on CL/F = 0.419 (RSE 13.9%); bootstrap 0.287-0.560

    # Inter-individual variability (exponential, Section 4.3). Table 2 prints
    # IIV as a percentage; converted to log-scale variance with
    # omega^2 = log(1 + CV^2):
    #   CL/F: log(1 + 0.313^2) = 0.09346
    #   ka:   log(1 + 1.261^2) = 0.95170
    # Covariance from the printed correlation -34.2%:
    #   -0.342 * sqrt(0.09346 * 0.95170) = -0.10200
    etalcl + etalka ~ c(0.09346, -0.10200, 0.95170) # Table 2: IIV CL/F 31.3% (shrinkage 14.9%), IIV Ka 126.1% (shrinkage 63.5%), correlation -34.2%

    # Residual error: combined model (Section 2). Additive term is an SD in ug/L
    # (Table 2 unit tag); the proportional term is read as an SD (fraction).
    addSd <- 8.14; label("Additive residual error (ug/L)") # Table 2: Additional Error = 8.14 ug/L (RSE 17.6%); bootstrap 4.33-14.80
    propSd <- 0.0689; label("Proportional residual error (fraction)") # Table 2: Proportional Error = 0.0689 (RSE 31.3%); bootstrap 0.0136-0.135
  })

  model({
    # ---- 1. Individual parameters ----
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 78.9)^e_crcl_cl
    vc <- exp(lvc)
    tlag <- exp(ltlag)

    # ---- 2. Micro-constant ----
    kel <- cl / vc

    # ---- 3. ODE system: one compartment, first-order absorption ----
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    # ---- 4. Absorption lag time ----
    alag(depot) <- tlag

    # ---- 5. Observation and error ----
    # central / vc is mg/L; 1 mg/L = 1000 ug/L.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
