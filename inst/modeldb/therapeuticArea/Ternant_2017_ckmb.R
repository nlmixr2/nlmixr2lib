Ternant_2017_ckmb <- function() {
  description <- paste(
    "Population kinetic model of serum creatine kinase MB isoenzyme (CK-MB) released by",
    "the injured myocardium after ST-segment elevation myocardial",
    "infarction (STEMI) treated by primary percutaneous coronary",
    "intervention (Ternant 2017). The CK-MB level at the first blood sample",
    "(time zero, before PCI) is an estimated initial condition B0; the",
    "remainder of the total released amount (B = BTOT - B0) enters the",
    "central compartment through a single Savic transit-release",
    "function, with a zero-order baseline",
    "production and first-order elimination (one compartment). BTOT",
    "rises exponentially with the area at risk and is lower in",
    "patients randomized to a cardioprotective conditioning arm",
    "(remote ischemic perconditioning with or without local ischemic",
    "postconditioning, RIPOST-MI trial). The model works in",
    "concentration units throughout (no volume).",
    sep = " "
  )
  reference <- paste(
    "Ternant D, Ivanes F, Prunier F, Mewton N, Bejan-Angoulvant T,",
    "Paintaud G, Ovize M, Angoulvant D. Revisiting myocardial necrosis",
    "biomarkers: assessment of the effect of conditioning therapies on",
    "infarct size by kinetic modelling. Sci Rep. 2017;7:10709.",
    "doi:10.1038/s41598-017-11352-4.",
    sep = " "
  )
  vignette <- "Ternant_2017_myocardial_necrosis_biomarkers"
  units <- list(
    time = "h",
    dosing = "none (endogenous release; no dose records)",
    concentration = "IU/L"
  )

  covariateData <- list(
    AAR_PCT = list(
      description = paste(
        "Myocardial area at risk, as the circumferential extent of",
        "abnormally contracting segments (ACS) on biplane left",
        "ventriculography, in percent of the left-ventricular",
        "circumference"
      ),
      units = "% of ACS",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed, measured at the index PCI. Uncentred exponential",
        "effect on BTOT: BTOT = BTOT0 * exp(e_aar_btot * AAR_PCT), so the",
        "tabulated BTOT is the value at AAR = 0 (Methods, 'Covariates').",
        "Cohort median 37.4% [IQR 30.5-45.7] (Table 1, CK-MB column)."
      ),
      source_name = "AAR"
    ),
    TRT_CONDITIONING = list(
      description = paste(
        "Randomized to a cardioprotective conditioning arm (1) versus the",
        "control arm (0)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (control arm)",
      notes = paste(
        "For the CK-MB data (RIPOST-MI trial) 1 = remote ischemic",
        "perconditioning (RIPer) with or without local ischemic",
        "postconditioning (IPost); the two conditioned arms did not differ",
        "(p = 0.89) and were merged against control. Effect on log(BTOT):",
        "additive beta = -0.75 (Results, 'Covariate models'; Table 2 'Arm",
        "on BTOT')."
      ),
      source_name = "ARM"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "creatine kinase MB isoenzyme",
      units = "IU/L",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    n_observations = 735L,
    age_range = "59 years median [IQR 49-72]",
    weight_range = "76 kg median [IQR 70-83]",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Adults with STEMI treated by primary PCI within 6 h of symptom",
      "onset, enrolled in the RIPOST-MI trial (Angers University Hospital)."
    ),
    dose_range = paste(
      "No drug dose; arms were control (16), remote ischemic",
      "perconditioning (17) and perconditioning plus local ischemic",
      "postconditioning (16)."
    ),
    regions = "France",
    notes = paste(
      "Table 1 'CK-MB' column. Model built on a random 2/3 learning subset",
      "(32 patients) and externally validated on the remaining 1/3 (17).",
      "Table 1 prints sex as 'female/male 40/19', which sums to 59 rather",
      "than 49 and is very likely transposed, so sex_female_pct is left",
      "unset. Samples were taken before PCI (time zero), every 4 h for",
      "24 h, then every 6 h for 48 h."
    )
  )

  ini({
    # Release amount (concentration units, IU/L), with covariates
    lbtot <- log(167)
    label("Total CK-MB released by the lesion at AAR = 0, control arm, BTOT (IU/L)") # Table 2 CK-MB 'B TOT' = 167 (RSE 18%)
    e_trt_btot <- -0.75
    label("Additive effect of conditioning arm on log(BTOT) (unitless)") # Table 2 CK-MB 'Arm on BTOT' = -0.75 (RSE 29%)
    e_aar_btot <- 0.037
    label("Exponential effect of area at risk on BTOT (1/% ACS)") # Table 2 CK-MB 'AAR on BTOT' = 0.037 (RSE 24%)
    lb0 <- log(43.1)
    label("CK-MB level at the first blood sample, initial condition B0 (IU/L)") # Table 2 CK-MB 'B1' = 43.1 (RSE 36%)

    # Single-transit release (Savic transit function)
    lntr1 <- log(1.0)
    label("Number of transit compartments of the release chain, n (unitless)") # Table 2 CK-MB 'n1' = 1.0 (RSE 19%)
    lktr1 <- log(0.34)
    label("Transit rate constant of the release chain, ktr (1/h)") # Table 2 CK-MB 'ktr1' = 0.34 (RSE 12%)

    # Distribution / elimination (one compartment)
    lkin <- log(0.32)
    label("Zero-order baseline CK-MB production, kprod (IU/L/h)") # Table 2 CK-MB 'kprod' = 0.32 (RSE 8%)
    lkel <- log(0.096)
    label("First-order CK-MB elimination rate constant, ke (1/h)") # Table 2 CK-MB 'ke' = 0.096 (RSE 2%)

    # IIV: Table 2 reports omega as SDs; variance = omega^2.
    etalbtot ~ 0.2601 # Table 2 CK-MB 'omega Btot' = 0.51 (RSE 14%) -> 0.51^2
    etalb0 ~ 1.21 # Table 2 CK-MB 'omega B1' = 1.1 (RSE 29%) -> 1.1^2
    etalntr1 ~ 0.4096 # Table 2 CK-MB 'omega n1' = 0.64 (RSE 23%) -> 0.64^2
    etalktr1 ~ 0.0324 # Table 2 CK-MB 'omega Ktr1' = 0.18 (RSE 22%) -> 0.18^2
    etalkin ~ 0.1369 # Table 2 CK-MB 'omega Kprod' = 0.37 (RSE 15%) -> 0.37^2
    etalkel ~ 0.01 # Table 2 CK-MB 'omega Ke' = 0.10 (RSE 14%) -> 0.10^2

    # Combined additive + proportional residual error (Methods, 'Error model')
    addSd <- 0.19
    label("Additive residual error (IU/L)") # Table 2 CK-MB 'sigma add' = 0.19 (RSE 20%)
    propSd <- 0.12
    label("Proportional residual error (fraction)") # Table 2 CK-MB 'sigma prop' = 0.12 (RSE 5%)
  })

  model({
    # Total released amount and its pre-sampling part (Methods, 'Biomarker
    # release'; 'Covariates': ln(theta) shifted by beta_ARM, theta * exp(beta_AAR * AAR))
    btot <- exp(lbtot + etalbtot + e_trt_btot * TRT_CONDITIONING + e_aar_btot * AAR_PCT)
    b0 <- exp(lb0 + etalb0)
    brel <- btot - b0

    ntr1 <- exp(lntr1 + etalntr1)
    ktr1 <- exp(lktr1 + etalktr1)
    kin <- exp(lkin + etalkin)
    kel <- exp(lkel + etalkel)

    # Savic transit release function (paper reference 43), time from the
    # first blood sample; integrates to 1 over t > 0.
    tpos <- t + 1e-10
    rel1 <- exp(log(ktr1) + ntr1 * log(ktr1 * tpos) - ktr1 * tpos - lgamma(ntr1 + 1))
    release <- brel * rel1

    # Eq. 3: dB/dt = R(t) + kprod - ke * B, with B(0) = B0
    d/dt(central) <- release + kin - kel * central
    central(0) <- b0

    ck_mb <- central
    ck_mb ~ add(addSd) + prop(propSd)
  })
}
