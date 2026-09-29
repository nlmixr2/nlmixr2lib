Jang_2019_cefprozil_total <- function() {
  description <- "Population pharmacokinetic model for TOTAL cefprozil (cis- plus trans-isomer) after a single 1000 mg oral dose in healthy adult Korean males: one compartment with first-order absorption, an absorption lag time and first-order elimination. Creatinine clearance enters apparent clearance as a linear centred effect around the cohort median of 124.41 mL/min. The kinetics are flip-flop (Ka 0.432 1/h is well below Kel 1.171 1/h), so the terminal slope is set by absorption rather than by CL/V. One of three independently fitted models in Jang 2019; see also modellib('Jang_2019_cefprozil_cis') and modellib('Jang_2019_cefprozil_trans')."
  reference <- paste(
    "Jang JH, Jeong SH, Cho HY, Lee YB. (2019).",
    "Population Pharmacokinetics of Cis-, Trans-, and Total Cefprozil",
    "in Healthy Male Koreans.",
    "Pharmaceutics 11(10):531.",
    "doi:10.3390/pharmaceutics11100531.",
    sep = " "
  )
  vignette <- "Jang_2019_cefprozil"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units = "mL/min (raw Cockcroft-Gault, NOT BSA-normalized)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Jang 2019 Methods 2.3: 'CrCl was estimated using the Cockcroft-Gault formula'.",
        "Enters CL as the linear centred multiplier (1 + e_crcl_cl * (CRCL - 124.41)),",
        "per the Total cefprozil final-model equation on p. 7.",
        "The centring constant 124.41 mL/min is the cohort MEDIAN CrCl from Table 1;",
        "the cohort range is 86.57-159.05 mL/min, so this model is only informed over",
        "a narrow band of supranormal renal function in healthy young men and should not",
        "be extrapolated to renal impairment without care (Jang 2019 Discussion makes the",
        "same point). Values are raw mL/min and are NOT normalized to 1.73 m^2."
      ),
      source_name = "CrCl"
    )
  )

  # Screened by Jang 2019 but NOT retained in any final model (Table 3 stepwise
  # search plus Methods 2.3 candidate list). Documented here so the covariate
  # screen's provenance is preserved without declaring unused covariates.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; cohort 21-27 years (Table 1). Not carried into the stepwise search reported in Table 3."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Table 3 'Weight on volume' for total cefprozil: dOFV -2.904, short of the -3.84 forward-selection threshold. Not retained."
    ),
    BSA = list(
      description = "Body surface area (Mosteller equation)",
      units = "m^2",
      type = "continuous",
      notes = "Table 3 'BSA on clearance' for total cefprozil: dOFV -2.413, short of the -3.84 threshold. Not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/dL as reported by Jang 2019 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Table 3 'Total protein on clearance' for total cefprozil: dOFV -0.136. Not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL as reported by Jang 2019 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Table 3 'Albumin on clearance' for total cefprozil: dOFV -0.487. Not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Not carried into Table 3."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Not carried into Table 3."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; enters the model only indirectly, as the Cockcroft-Gault input to CRCL. Not retained as a covariate in its own right."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "total cefprozil (cis + trans)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "total cefprozil (cis + trans)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    age_range = "21-27 years",
    age_median = "24 years",
    weight_range = "53.1-91.8 kg",
    weight_median = "69.5 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers",
    dose_range = "Single 1000 mg oral dose of cefprozil with 240 mL of water after an overnight fast (Jang 2019 Methods 2.1). For THIS model the dose amount is the full 1000 mg, because the observation is total (cis + trans) cefprozil; see notes for the cis / trans dose basis.",
    regions = "Republic of Korea (Chonnam National University, Gwangju)",
    renal_function = "Normal to supranormal; CrCl 86.57-159.05 mL/min, median 124.41 mL/min (Cockcroft-Gault, Table 1)",
    notes = paste(
      "35 healthy Korean males from the reference arm of a single-dose, randomized,",
      "two-way, open-label, crossover bioequivalence study (Bioequivalence Test No. 611).",
      "420 plasma samples per analyte; sampling at 0, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 8",
      "and 12 h post dose. Cis- and trans-cefprozil were assayed separately by UPLC-ESI-MS/MS",
      "(LLOQ 5 and 15 ng/mL respectively) and total cefprozil was computed as their sum.",
      "DOSE BASIS: Jang 2019 never states the dose amount entered into each analyte's dataset.",
      "It is back-solved here from the paper's own reported noncompartmental results",
      "(Discussion, p. 13) and is the isomer content of the 1000 mg dose at the ~9:1",
      "cis:trans ratio quoted in the Introduction: 1000 mg for total, 900 mg for cis,",
      "100 mg for trans. Inverting Dose = AUC * CL on the three reported mean NCA AUCs",
      "gives 1022 mg (cis), 113 mg (trans) and 1137 mg (total) -- a back-solved cis:trans",
      "ratio of 9.05:1, and a cis + trans sum (1135 mg) that matches the total (1137 mg)",
      "to 0.2%. The competing reading that all three datasets were dosed at 1000 mg is",
      "falsified: it would predict a trans AUC of 56.5 ug*h/mL against the 6.38 ug*h/mL",
      "reported. See the vignette Errata for the full arithmetic.",
      "PPK analysis was run in Phoenix NLME 8.1 (FOCE with extended least squares),",
      "not NONMEM."
    )
  )

  ini({
    # Structural parameters -- Jang 2019 Table 4, 'Total cefprozil / Final model'.
    # Volume and clearance are printed in mL and mL/h; they are divided by 1000
    # here so that a dose in mg gives a concentration in mg/L = ug/mL, the
    # concentration unit the paper reports.
    lka <- log(0.432); label("First-order absorption rate constant Ka (1/h)") # Table 4, Total cefprozil Final model: tvKa = 0.432 1/h (RSE 1.28%)
    lvc <- log(14713.10 / 1000); label("Apparent central volume of distribution V/F (L)") # Table 4, Total cefprozil Final model: tvV = 14,713.10 mL (RSE 8.16%)
    lcl <- log(17226.20 / 1000); label("Apparent clearance CL/F (L/h)") # Table 4, Total cefprozil Final model: tvCl = 17,226.20 mL/h (RSE 2.35%)
    ltlag <- log(0.351); label("Absorption lag time (h)") # Table 4, Total cefprozil Final model: tvTlag = 0.351 h (RSE 3.74%)

    # Covariate effect. Linear and centred, NOT log-scale: the final-model
    # equation on p. 7 reads Cl = Cltv * (1 + (CrCl - 124.41) * dCldCrCl) * exp(etaCl).
    # Table 4 prints this coefficient rounded to 0.003; Table 5 prints the same
    # final-model estimate to three significant figures as 2.87e-3, which is used here.
    e_crcl_cl <- 2.87e-3; label("Linear effect of creatinine clearance on CL/F (per mL/min, centred at 124.41 mL/min)") # Table 5, Total cefprozil: dCldCrCl = 2.87 x 10^-3 (95% CI 3.00e-4 to 6.04e-3); Table 4 rounds to 0.003, RSE 56.16%

    # IIV. Jang 2019 reports omega^2 (variances) directly, so these are entered
    # as variances without transformation. No IIV was estimated on Ka for total
    # cefprozil: Table 2 selects step 02-03-03 'Remove IIV Ka' as the final IIV model.
    etalvc ~ 0.124 # Table 4, Total cefprozil Final model: omega^2 V = 0.124 (shrinkage 7.08%); Discussion quotes omega = 35.3%, i.e. sqrt(0.124)
    etalcl ~ 0.016 # Table 4, Total cefprozil Final model: omega^2 Cl = 0.016 (shrinkage 10.39%)
    etaltlag ~ 0.021 # Table 4, Total cefprozil Final model: omega^2 Tlag = 0.021 (shrinkage 21.39%)

    # Residual error. Table 2 step 02-03 'Log additive' is the SELECTED residual
    # model for total cefprozil, i.e. additive error on log-transformed data
    # (Methods Equation 2), which is lnorm() in nlmixr2. Table 4 labels sigma
    # '(ug/mL)', but a log-scale SD is dimensionless; see vignette Errata.
    expSd <- 0.189; label("Log-scale additive residual standard deviation") # Table 4, Total cefprozil Final model: sigma = 0.189 (RSE 9.02%)
  })

  model({
    # Individual parameters
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 124.41))
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
