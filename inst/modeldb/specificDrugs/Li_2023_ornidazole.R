Li_2023_ornidazole <- function() {
  description <- paste(
    "One-compartment intravenous population PK of ornidazole in breastfeeding",
    "women after caesarean section, with breast-milk concentration linked to",
    "plasma by an estimated milk-to-plasma concentration ratio that rises as a",
    "power of time postpartum. Apparent clearance decreases with total",
    "bilirubin."
  )
  reference <- paste(
    "Li S, Cao M, Zhou Y, Shu C, Wang Y. Ornidazole Transfer into Colostrum and",
    "Assessment of Exposure Risk for Breastfeeding Infant: A Population",
    "Pharmacokinetic Analysis. Pharmaceutics. 2023;15(11):2524.",
    "doi:10.3390/pharmaceutics15112524"
  )
  vignette <- "Li_2023_ornidazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    TBILI = list(
      description = "Total bilirubin; power effect on apparent clearance",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Population median 9.1 umol/L (Table 1) is the centering value in",
        "Equation 6. The paper reports TBIL in umol/L, which is already the",
        "canonical TBILI unit, so no conversion is applied. Simulations in the",
        "source paper stratify at 17 umol/L (normal vs abnormal liver",
        "function)."
      ),
      source_name = "TBIL"
    ),
    TPP = list(
      description = "Time postpartum; power effect on the milk-to-plasma concentration ratio",
      units = "week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column PST is the postpartum sampling time in HOURS, spanning",
        "0-120 h (Figure 4A x-axis) over the colostrum period. Canonical TPP is",
        "in weeks, so model() converts to hours before forming the ratio in",
        "Equation 5. Only the ratio enters, so the conversion is exact. The",
        "centering value of 54 h is back-solved from four paper-reported",
        "anchors, not reported by the authors; see the validation vignette",
        "Errata."
      ),
      source_name = "PST"
    )
  )

  # Covariates the authors recorded and screened (Methods 2.3, Table 2) but did
  # not retain in the final model. Documentation only; not referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = "Screened on CL (Table 2 steps 3, 12, 14); dropped in backward elimination (step 15, dOFV 5.90, p > 0.01)."
    ),
    AGE = list(
      description = "Maternal age", units = "year", type = "continuous",
      notes = "Screened on CL (Table 2 step 4); not significant (dOFV 0.25)."
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance", units = "mL/min", type = "continuous",
      notes = "Screened on CL (Table 2 steps 5, 13); dropped in favour of TBILI."
    ),
    SCR = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      notes = "Screened on CL (Table 2 step 6); not significant (dOFV 1.39)."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "IU/L", type = "continuous",
      notes = "Screened on CL (Table 2 step 8); not significant (dOFV 0.04)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "IU/L", type = "continuous",
      notes = "Screened on CL (Table 2 step 9); not significant (dOFV 0.07)."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "IU/L", type = "continuous",
      notes = "Screened on CL (Table 2 step 10); not significant (dOFV 2.04)."
    ),
    BUN = list(
      description = "Blood urea nitrogen", units = "mmol/L", type = "continuous",
      notes = "Recorded and screened as a renal marker (Methods 2.3); not retained."
    ),
    UA = list(
      description = "Uric acid", units = "umol/L", type = "continuous",
      notes = "Recorded and screened as a renal marker (Methods 2.3); not retained."
    ),
    CYSC = list(
      description = "Serum cystatin C", units = "mg/L", type = "continuous",
      notes = "Recorded and screened as a renal marker (Methods 2.3); not retained."
    ),
    DBILI = list(
      description = "Direct bilirubin", units = "umol/L", type = "continuous",
      notes = "Recorded and screened as a hepatic marker (Methods 2.3); total bilirubin retained instead."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "ornidazole", units = "mg", specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    n_observations = "87 plasma + 123 breast-milk samples",
    age_range = "18.0-41.0 years",
    age_median = "30.0 years",
    weight_range = "43.7-90.0 kg",
    weight_median = "64.0 kg",
    sex_female_pct = 100,
    disease_state = paste(
      "postpartum after caesarean delivery; prophylaxis or treatment of",
      "anaerobic infection"
    ),
    dose_range = paste(
      "1000 mg IV 1-2 h before the procedure, then 500 mg IV at 12 h and 500 mg",
      "IV at 24 h after delivery (2000 mg total)"
    ),
    regions = "China (Wuhan)",
    lactation_stage = "colostrum, days 1-4 postpartum",
    renal_function = "normal; CrCL 108.52-307.12 mL/min (median 185.51)",
    hepatic_function = paste(
      "TBIL 3.4-20.9 umol/L (median 9.1); severe hepatic impairment excluded"
    ),
    notes = paste(
      "Prospective trial at Wuhan Children's Hospital (ethics 2021R141-E01);",
      "baseline demographics in Table 1. Estimation in Phoenix NLME 8.3.4 by",
      "FOCE-ELS, so no NONMEM control stream exists. IIV on V was omitted by",
      "the authors because of high shrinkage (Table 3 footnote)."
    )
  )

  ini({
    lvc <- log(35.75)
    label("Volume of distribution (V, L)")
    # Table 3: V = 35.75 L (RSE 4.89%; bootstrap 31.99-39.81)
    lcl <- log(1.89)
    label("Apparent clearance (CL, L/h)")
    # Table 3: CL = 1.89 L/h (RSE 2.78%; bootstrap 1.78-2.00)
    lcmpr <- log(0.58)
    label("Milk-to-plasma concentration ratio (MPRcon, unitless)")
    # Table 3: MPRcon = 0.58 (RSE 8.63%; bootstrap 0.48-0.68)

    e_tbili_cl <- -0.17
    label("Total bilirubin power exponent on CL (unitless)")
    # Table 3: theta_TBIL = -0.17 (RSE -31.99%); Equation 6 exponent
    e_tpp_cmpr <- 1.37
    label("Time-postpartum power exponent on the milk-to-plasma ratio (unitless)")
    # Table 3: theta_PST = 1.37 (RSE 15.32%); Equation 5 exponent

    etalcl ~ 0.024
    # Table 3: omega^2_CL = 0.024 (variance; bootstrap 0.011-0.035)
    etalcmpr ~ 0.327
    # Table 3: omega^2_MPRcon = 0.327 (variance; bootstrap 0.189-0.449)
    # IIV on V is deliberately absent: Table 3 footnote reports it was omitted
    # from the model because of a large amount of shrinkage.

    propSd <- 0.0818
    label("Proportional residual error, plasma (fraction)")
    # Table 3: sigma1 = 8.18% proportional error on plasma concentrations
    propSd_Cmilk <- 0.3175
    label("Proportional residual error, breast milk (fraction)")
    # Table 3: sigma2 = 31.75% proportional error on milk concentrations
  })

  model({
    # Equation 6: CL = TVCL * (TBIL / MedianTBIL)^theta_TBIL.
    # MedianTBIL = 9.1 umol/L is the Table 1 population median.
    cl <- exp(lcl + etalcl) * (TBILI / 9.1)^e_tbili_cl
    vc <- exp(lvc)

    # Equation 5: MPRcon = TVMPRcon * (PST / MedianPST)^theta_PST.
    # The paper's PST is in hours; canonical TPP is in weeks, so convert first
    # (1 week = 168 h). Only the ratio enters, so this is exact.
    # MedianPST = 54 h is BACK-SOLVED -- the authors never report it. Four
    # independent paper-reported anchors agree (Table S2 MPRauc for days 1-4,
    # the Figure 4A median curve, the Figure 5A milk Cmax, and the Figure 5
    # safety gap); see the validation vignette Errata.
    tppHours <- TPP * 168
    cmpr <- exp(lcmpr + etalcmpr) * (tppHours / 54)^e_tpp_cmpr

    # Equations 2 and 4: dA/dt = -CL * C1 with C1 = A / V.
    d/dt(central) <- -(cl / vc) * central
    Cc <- central / vc
    # Equation 3: C2 = MPRcon * C1. Breast-milk concentration is an algebraic
    # observable, not a compartment -- the authors state the milk data were too
    # sparse to support a separate milk compartment (Discussion).
    Cmilk <- cmpr * Cc

    Cc ~ prop(propSd)              # Table 3 sigma1, plasma
    Cmilk ~ prop(propSd_Cmilk)     # Table 3 sigma2, breast milk
  })
}
