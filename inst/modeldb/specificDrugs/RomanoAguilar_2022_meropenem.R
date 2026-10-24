RomanoAguilar_2022_meropenem <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for",
    "intravenously infused meropenem in critically ill Mexican adults",
    "(Romano-Aguilar 2022; 78 patients, sparse sampling). Clearance",
    "(11.9 L/h at the reference) is proportional to Cockcroft-Gault",
    "creatinine clearance normalised to 102.23 mL/min/1.73 m^2; volume of",
    "distribution 25.2 L. Exponential interindividual variability on",
    "clearance and volume and additive residual error."
  )
  reference <- paste(
    "Romano-Aguilar M, Ortiz-Alvarez A, Medellin-Garibay S,",
    "Martinez-Gutierrez F, Jung-Cook H, Milan-Segovia RC, Romano-Moreno S",
    "(2022). 611. Meropenem Dosage Optimization in Critically Ill Patients",
    "Based on a Population Pharmacokinetic Approach. Open Forum Infectious",
    "Diseases 9(Suppl 2):S334-S335 (IDWeek 2022 poster abstract).",
    "doi:10.1093/ofid/ofac492.663.",
    "Parameter estimates from the abstract's Table 2 (an image panel) and",
    "the covariate equation from its Results text."
  )
  vignette <- "RomanoAguilar_2022_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated creatinine clearance by the Cockcroft-Gault formula, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source name CLCr. Table 1 reports 'Estimated creatinine clearance",
        "(mL/min/1.73 m2)', footnoted 'Estimated by the Cockcroft-Gault",
        "formula', median 103.14 (range 13.5-241.4); the Table 2 footnote",
        "repeats 'CLCr: creatinine clearance in mL/min/1.73 m2 calculated by",
        "Cockcroft-Gault formula'. Enters clearance linearly through the",
        "origin: Results 'CL (L/h) = 11.9 * (CLCr/102.23)'. Table 2 prints",
        "the same relation rounded as 'CL = theta1 x (CLCr/102)'. The",
        "abstract does not state what statistic 102.23 is (Table 1's median",
        "is 103.14), nor the body weight used in the Cockcroft-Gault formula",
        "or how it was normalised to 1.73 m^2."
      ),
      source_name = "CLCr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_range = "17-92 years (median 47)",
    weight_range = "40-98 kg (median 67)",
    height_median = "1.62 m (range 1-1.9)",
    sex_female_pct = 100 * 45 / 78,
    race_ethnicity = "Mexican",
    disease_state = paste(
      "Critically ill adults receiving meropenem for serious infection,",
      "prescribed on clinical, biochemical and microbiological findings."
    ),
    renal_function = paste(
      "Serum creatinine median 0.93 mg/dL (range 0.1-5.11); Cockcroft-Gault",
      "creatinine clearance median 103.14 mL/min/1.73 m^2 (range",
      "13.5-241.4), including patients with augmented renal clearance."
    ),
    dose_range = paste(
      "Meropenem 500-6000 mg/day as intermittent intravenous infusion.",
      "Daily dose 3000 mg in 52 of 78 patients (500 mg 1, 1000 mg 6,",
      "1500 mg 2, 2000 mg 11, 6000 mg 6); dosing interval 8 h in 60 (6 h 4,",
      "12 h 12, 24 h 2); infusion time 0.5 h in 33, 1 h in 16, 2 h in 6,",
      "3 h in 23 (Table 1)."
    ),
    regions = "Mexico (Hospital Central 'Dr. Ignacio Morones Prieto', San Luis Potosi).",
    notes = paste(
      "Prospective observational study. Blood sampled pre-dose and 1, 3 and",
      "6 h post-dose; plasma meropenem by HPLC. Population PK modelling and",
      "Monte Carlo simulation in NONMEM; internal validation by bootstrap",
      "(n = 1000) and VPC, external validation against a separate dataset",
      "(not described). Demographics from Table 1 of the abstract."
    )
  )

  ini({
    # Structural parameters -- Table 2 'Mean' column (final model; bootstrap
    # medians 11.93 and 24.92 agree).
    lcl <- log(11.9); label("Clearance at CRCL = 102.23 mL/min/1.73 m^2 (L/h)")  # Table 2 theta1 = 11.9 (RSE 8%); Results 'CL (L/h) = 11.9 * (CLCr/102.23)'
    lvc <- log(25.2); label("Volume of distribution (L)")                        # Table 2 theta2 = 25.2 (RSE 11%); Results 'V (L) = 25.2'

    # Interindividual variability. Table 2 reports both terms as CV% under
    # an omega^2 symbol and does not give the conversion; the log-normal
    # form omega^2 = log(1 + CV^2) is used (see vignette Assumptions).
    etalcl ~ 0.27446  # Table 2 'Interindividual variability on CL (CV%)' = 56.2 (RSE 8%, shrinkage 7%); log(1 + 0.562^2)
    etalvc ~ 0.19959  # Table 2 'Interindividual variability on V (CV%)' = 47 (RSE 18%, shrinkage 31%); log(1 + 0.47^2)

    # Residual error. Table 2 'Residual variability (ug/mL)' sigma = 3.53
    # (RSE 42%, shrinkage 28%); reported in concentration units, so additive.
    addSd <- 3.53; label("Additive residual error (mg/L)")  # Table 2 sigma = 3.53 ug/mL
  })

  model({
    # 1. Individual PK parameters. Clearance is proportional to
    # Cockcroft-Gault creatinine clearance (Results:
    # CL (L/h) = 11.9 * (CLCr/102.23)); volume has no covariate.
    cl <- exp(lcl + etalcl) * (CRCL / 102.23)
    vc <- exp(lvc + etalvc)

    # 2. One-compartment disposition with linear elimination. Meropenem was
    # given as intravenous infusions (0.5-3 h) into the central compartment.
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # 3. Observation. Plasma meropenem in mg/L (= ug/mL), additive error.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
