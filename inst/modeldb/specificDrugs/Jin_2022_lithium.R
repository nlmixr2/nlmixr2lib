Jin_2022_lithium <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model with first-order",
    "absorption (ka fixed) and first-order elimination for lithium in 268",
    "Chinese patients with bipolar disorder on lithium carbonate maintenance",
    "therapy (trough therapeutic-drug-monitoring samples), from Jin 2022.",
    "Apparent clearance CL/F scales as a power of total daily lithium",
    "carbonate dose (centred on 600 mg/day), body weight (62 kg) and",
    "Cockcroft-Gault creatinine clearance (116 mL/min). Doses are in mmol of",
    "lithium ion (1 mg lithium carbonate = 2 / 73.89 mmol Li) and",
    "concentrations in mmol/L.",
    sep = " "
  )
  reference <- paste(
    "Jin Z-b, Wu Z, Cui Y-f, Liu X-p, Liang H-b, You J-y, Wang C-y.",
    "Population Pharmacokinetics and Dosing Regimen of Lithium in Chinese",
    "Patients With Bipolar Disorder.",
    "Front Pharmacol. 2022;13:913935. doi:10.3389/fphar.2022.913935.",
    "PMC9289112. Parameter estimates are in Table 2; the final covariate",
    "model is Results Equations 8-10.",
    sep = " "
  )
  vignette <- "Jin_2022_lithium"
  units <- list(time = "h", dosing = "mmol", concentration = "mmol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F centred on the cohort median of 62 kg",
        "(Equation 8, (WT/62)^0.33; Table 1 median 62.0, range 35.0-110 kg).",
        sep = " "
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = paste(
        "Creatinine clearance by the Cockcroft-Gault equation, RAW mL/min and",
        "NOT BSA-normalised. Table 1 footnote: CRCL = [(140 - Age) x weight",
        "(kg)] / [0.818 x Scr (umol/L)] x k, with k = 1 for men and 0.85 for",
        "women; the weight is total body weight.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F centred on the cohort median of 116 mL/min",
        "(Equation 8, (CRCL/116)^0.186; Table 1 median 116, range",
        "61.7-226 mL/min). Simulations in the paper extrapolate down to",
        "30 mL/min (Figure 3), below the observed range.",
        sep = " "
      ),
      source_name = "CRCL"
    ),
    DOSE_LITHIUM_CARBONATE_MGD = list(
      description = paste(
        "The patient's own total daily dose of lithium carbonate (mg of the",
        "salt per day, NOT mmol of lithium ion), summed over the daily",
        "administrations.",
        sep = " "
      ),
      units = "mg/d",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F centred on the cohort median of 600 mg/day",
        "(Equation 8, (TDD/600)^0.354; Table 1 median 600, range",
        "150-1500 mg/day). Per-dose-record covariate; it must agree with the",
        "dose events, which are given in mmol of lithium ion",
        "(DOSE_LITHIUM_CARBONATE_MGD * 2 / 73.89 mmol Li per day). The",
        "positive exponent makes CL/F rise with dose, which the authors",
        "attribute to saturable tubular reabsorption of lithium (Discussion).",
        sep = " "
      ),
      source_name = "TDD"
    )
  )

  # Covariates the source screened and did not retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The CL/F random effects were correlated with sex in the graphical",
        "screen, but adding sex in the stepwise forward inclusion 'did not",
        "meet the criteria for statistical significance (p < 0.05)'",
        "(Results 3.2). No coefficient is reported.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lithium", units = "mmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lithium", units = "mmol", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 268L,
    n_studies = 1L,
    n_observations = 476L,
    age_range = "13-77 years",
    age_median = "31.0 years",
    age_mean = "35.0 +/- 14.5 years",
    weight_range = "35.0-110 kg",
    weight_median = "62.0 kg",
    weight_mean = "63.7 +/- 11.1 kg",
    sex_female_pct = 66.8,
    race_ethnicity = c(Asian = 100),
    disease_state = "Bipolar disorder on lithium carbonate maintenance treatment.",
    dose_range = paste(
      "Oral lithium carbonate 150-1500 mg/day (median 600 mg/day); 76.1%",
      "sustained-release and 23.9% ordinary tablets.",
      sep = " "
    ),
    renal_function = "Cockcroft-Gault CRCL 61.7-226 mL/min (median 116).",
    regions = "Xuzhou, Jiangsu, China",
    notes = paste(
      "Retrospective therapeutic-drug-monitoring cohort at the Affiliated",
      "Xuzhou Eastern Hospital of Xuzhou Medical University, September",
      "2016-August 2021 (Methods 2.1). 241 adults (89.9%) and 27 children",
      "(age 16 years or younger). All serum lithium samples were troughs",
      "taken before the morning dose (ADVIA 1800, calibration range",
      "0.19-3.0 mmol/L). Patients on interacting medication (diuretics,",
      "renin-angiotensin system antagonists, serotonergic drugs) were",
      "excluded. NONMEM 7.4.2, FOCE-I. Baseline demographics: Table 1.",
      sep = " "
    )
  )

  ini({
    # All values: Jin 2022 Table 2 ('Population-pharmacokinetic parameter
    # estimates and bootstrap evaluation'), column 'Final model, Parameter
    # estimates (RSE%)', and Results Equations 8-10. All disposition
    # parameters are apparent (relative to the unestimated F).

    # ---- Absorption ------------------------------------------------
    lka <- fixed(log(0.293))
    label("First-order absorption rate constant ka (1/h)")
    # Table 2: Ka (h-1) = 0.293 [fixed]; Equation 10; Methods 2.2.1: 'fixed at
    # 0.293 h-1 based on published data because no sampling was collected
    # during the absorption phase'.

    # ---- Disposition -----------------------------------------------
    lcl <- log(0.909)
    label("Apparent clearance CL/F at TDD = 600 mg/day, WT = 62 kg, CRCL = 116 mL/min (L/h)")
    # Table 2: CL (L/h) = 0.909 (RSE 3%); Equation 8

    lvc <- log(10.9)
    label("Apparent volume of distribution V/F (L)")
    # Table 2: V (L) = 10.9 (RSE 12%); Equation 9

    # ---- Covariate effects on CL/F (power, median-centred) -----------
    e_dose_lithium_cl <- 0.354
    label("Power exponent of (DOSE_LITHIUM_CARBONATE_MGD / 600 mg/day) on CL/F (unitless)")
    # Table 2: 'TDD on CL' = 0.354 (RSE 12%); Equation 8 (TDD/600)^0.354

    e_wt_cl <- 0.33
    label("Power exponent of (WT / 62 kg) on CL/F (unitless)")
    # Table 2: 'WT on CL' = 0.33 (RSE 29%); Equation 8 (WT/62)^0.33

    e_crcl_cl <- 0.186
    label("Power exponent of (CRCL / 116 mL/min) on CL/F (unitless)")
    # Table 2: 'CRCL on CL' = 0.186 (RSE 29%); Equation 8 (CRCL/116)^0.186

    # ---- Between-subject variability -------------------------------
    # Exponential IIV (Methods 2.2.1). Equations 8 and 9 print the raw
    # omega^2 in the exponent; Table 2 reports the same values as
    # 100 * sqrt(omega^2) (sqrt(0.027) = 16.4%, sqrt(0.162) = 40.2%).
    etalcl ~ 0.027
    # Equation 8: e^0.027; Table 2: BSV CL = 16.4% (RSE 10%), shrinkage 30%
    etalvc ~ 0.162
    # Equation 9: e^0.162; Table 2: BSV V = 40.2% (RSE 20%), shrinkage 62%

    # ---- Residual error --------------------------------------------
    # Additive error only. Table 2 prints 0.0218 under the heading
    # 'Additive error (mmol/L)'. Read as the NONMEM $SIGMA variance, like the
    # raw omega^2 printed in Equations 8-9, because the Figure 1A
    # observed-vs-individual-prediction scatter has a residual SD of about
    # 0.13-0.15 mmol/L, matching sqrt(0.0218) = 0.148 at 20% shrinkage and
    # not 0.0218 (see the vignette).
    addSd <- 0.1476482
    label("Additive residual error SD for serum lithium (mmol/L)")
    # Table 2: additive error = 0.0218 (RSE 13%), shrinkage 20%; sqrt(0.0218) = 0.1476482
  })

  model({
    # ---- Individual parameters (Equations 8-10) ----------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (DOSE_LITHIUM_CARBONATE_MGD / 600)^e_dose_lithium_cl *
      (WT / 62)^e_wt_cl *
      (CRCL / 116)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # ---- ODE system --------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- Observation -------------------------------------------------
    # Doses in mmol lithium and volume in L, so central / vc is mmol/L.
    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
