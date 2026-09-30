Jeong_2021_cefaclor <- function() {
  description <- "Population pharmacokinetic model for cefaclor after a single 250 mg oral capsule (Ceclor) in healthy adult Korean males: one compartment with first-order absorption, an absorption lag time and first-order elimination. Creatinine clearance (Cockcroft-Gault, raw mL/min) enters apparent clearance and body weight enters apparent volume, both as power functions normalized to the cohort medians (110.92 mL/min and 66.05 kg)."
  reference <- paste(
    "Jeong SH, Jang JH, Cho HY, Lee YB. (2021).",
    "Population Pharmacokinetic Analysis of Cefaclor in Healthy Korean Subjects.",
    "Pharmaceutics 13(5):754.",
    "doi:10.3390/pharmaceutics13050754.",
    sep = " "
  )
  vignette <- "Jeong_2021_cefaclor"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V/F as the power term (WT / 66.05)^e_wt_vc (Jeong 2021 Section 3.4",
        "final-model equation). 66.05 kg is the cohort MEDIAN weight from Table 1",
        "(range 50.00-88.70 kg); the paper states continuous covariates were",
        "normalized to the observed median."
      ),
      source_name = "Weight"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units = "mL/min (raw Cockcroft-Gault, NOT BSA-normalized)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Jeong 2021 Section 2.3: 'CrCl was calculated using the commonly used",
        "Cockcroft-Gault equation'. Enters CL/F as the power term",
        "(CRCL / 110.92)^e_crcl_cl (Section 3.4 final-model equation). 110.92 mL/min",
        "is the cohort MEDIAN CrCl from Table 1 (range 67.50-170.42 mL/min). Values are",
        "raw mL/min and are NOT normalized to 1.73 m^2. The cohort is healthy young men,",
        "so the model is not informed below ~67 mL/min."
      ),
      source_name = "CrCl"
    )
  )

  # Screened by Jeong 2021 (Section 2.3 candidate list, Figure 2, Figure S1 and the
  # Table 4 stepwise search) but NOT retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected (Table 1, 19-26 years); no significant correlation with any PK parameter reported. Enters the model only indirectly, as a Cockcroft-Gault input to CRCL."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Collected (Table 1, 164.0-187.3 cm); not retained."
    ),
    BSA = list(
      description = "Body surface area (Mosteller equation)",
      units = "m^2",
      type = "continuous",
      notes = "Table 4 'BSA on CL/F': dOFV -2.755 (short of -3.84); 'CrCl and BSA on CL/F': dOFV +0.119 vs CrCl alone. Not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Table 4 'BMI on V/F': dOFV -2.831 (short of -3.84). Not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL as reported by Jeong 2021 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Section 2.3 candidate covariate; not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/dL as reported by Jeong 2021 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Section 3.4 lists total protein among candidates 'not valid for model improvement'."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      notes = "Section 3.4 lists BUN among candidates 'not valid for model improvement'."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Section 2.3 candidate covariate; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Section 3.4 lists AST among candidates 'not valid for model improvement'."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Section 3.4 lists ALT among candidates 'not valid for model improvement'."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Section 3.4 lists ALP among candidates 'not valid for model improvement'."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Collected (Table 1, 0.70-1.30 mg/dL); enters the model only indirectly, as the Cockcroft-Gault input to CRCL."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "cefaclor",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "cefaclor", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 2L,
    n_observations = 521L,
    age_range = "19-26 years",
    age_median = "23 years",
    weight_range = "50.0-88.7 kg",
    weight_median = "66.05 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers",
    dose_range = "Single 250 mg oral dose of cefaclor (Ceclor capsule, reference formulation) with 240 mL of water.",
    regions = "Republic of Korea (Chonnam National University, Gwangju)",
    renal_function = "Normal; CrCl 67.50-170.42 mL/min, median 110.92 mL/min (Cockcroft-Gault, Table 1)",
    notes = paste(
      "48 healthy Korean males from the reference (Ceclor) arms of two randomized, single-dose,",
      "open-label, two-way crossover bioequivalence studies (IRB permits 112, 2004, and 121, 2006;",
      "7-day washout). Sampling at 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4 and 5 h post dose;",
      "HPLC-UV assay, LLOQ 0.1 ug/mL; BLOQ samples treated as missing. Analysis in Phoenix NLME 8.3",
      "(FOCE with extended least squares, eta-epsilon interaction)."
    )
  )

  ini({
    # Structural parameters -- Jeong 2021 Table 5 (final model). Volume and clearance
    # are printed in mL and mL/h; they are divided by 1000 here so that a dose in mg
    # gives a concentration in mg/L = ug/mL, the concentration unit the paper reports.
    lka <- log(5.203); label("First-order absorption rate constant Ka (1/h)") # Table 5: tvKa = 5.203 1/h (RSE 18.02%)
    lvc <- log(22593.260 / 1000); label("Apparent volume of distribution V/F (L)") # Table 5: tvV/F = 22,593.260 mL (RSE 3.50%)
    lcl <- log(27166.883 / 1000); label("Apparent clearance CL/F (L/h)") # Table 5: tvCL/F = 27,166.883 mL/h (RSE 3.61%)
    ltlag <- log(0.245); label("Absorption lag time Tlag (h)") # Table 5: tvTlag = 0.245 h (RSE 0.60%)

    # Covariate effects -- power exponents on median-normalized covariates
    # (Section 3.4 final-model equation).
    e_crcl_cl <- 0.436; label("Power exponent of CrCl on CL/F (unitless)") # Table 5: dCL/FdCrCl = 0.436 (RSE 42.21%)
    e_wt_vc <- 0.581; label("Power exponent of body weight on V/F (unitless)") # Table 5: dV/FdWeight = 0.581 (RSE 31.54%)

    # IIV -- Table 5 reports omega^2 (variances). The IIV (%) column equals
    # 100 * sqrt(omega^2) (e.g. sqrt(0.971) = 0.985 -> 98.534%). No IIV on Tlag:
    # Table 3 selects step 02-01-04 'Remove IIV Tlag'.
    etalvc ~ 0.011 # Table 5: omega^2 V/F = 0.011 (RSE 51.90%, IIV 10.250%)
    etalcl ~ 0.034 # Table 5: omega^2 CL/F = 0.034 (RSE 30.22%, IIV 18.388%)
    etalka ~ 0.971 # Table 5: omega^2 Ka = 0.971 (RSE 30.55%, IIV 98.534%)

    # Residual error -- proportional (Table 3 step 02-01 selected). Phoenix NLME
    # reports the residual-error parameter as a standard deviation.
    propSd <- 0.270; label("Proportional residual error (fraction)") # Table 5: sigma = 0.270 (RSE 8.56%)
  })

  model({
    # Individual parameters (Section 3.4 final-model equation)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 66.05)^e_wt_vc
    cl <- exp(lcl + etalcl) * (CRCL / 110.92)^e_crcl_cl
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
