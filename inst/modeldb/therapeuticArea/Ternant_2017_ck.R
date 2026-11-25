Ternant_2017_ck <- function() {
  description <- paste(
    "Population kinetic model of serum creatine kinase (CK) released by",
    "the injured myocardium after ST-segment elevation myocardial",
    "infarction (STEMI) treated by primary percutaneous coronary",
    "intervention (Ternant 2017). The CK level at the first blood sample",
    "(time zero, before PCI) is an estimated initial condition B0; the",
    "remainder of the total released amount (B = BTOT - B0) enters the",
    "central compartment through the sum of two Savic transit-release",
    "functions (fraction Fr through chain 1), with a zero-order baseline",
    "production and first-order elimination (one compartment). BTOT",
    "rises exponentially with the area at risk and is lower in",
    "patients randomized to a cardioprotective conditioning arm",
    "(ischemic or cyclosporine A postconditioning). The model works in",
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
    AREA_AT_RISK_PCT = list(
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
        "effect on BTOT: BTOT = BTOT0 * exp(e_aar_btot * AREA_AT_RISK_PCT), so the",
        "tabulated BTOT is the value at AAR = 0 (Methods, 'Covariates').",
        "Cohort median 34.8% [IQR 27.0-46.5] (Table 1)."
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
        "For the CK data (PC, PC-CsA and PCNR trials) 1 = ischemic",
        "postconditioning (PC, PCNR) or cyclosporine A pharmacological",
        "postconditioning (PC-CsA). Effect on log(BTOT): additive",
        "beta = -0.43 (Methods, 'Covariates'; Table 2 'Arm on BTOT')."
      ),
      source_name = "ARM"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "creatine kinase",
      units = "IU/L",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 132L,
    n_studies = 3L,
    n_observations = 1965L,
    age_range = "58 years median [IQR 49-68]",
    weight_range = "75 kg median [IQR 67-84]",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Adults with STEMI (left anterior descending, right coronary or",
      "dominant left circumflex occlusion) treated by primary PCI within",
      "6 h of symptom onset, enrolled in the PC (two trials), PC-CsA and",
      "PCNR conditioning trials (Hospices Civils de Lyon)."
    ),
    dose_range = paste(
      "No drug dose; arms were control (61) versus ischemic or",
      "cyclosporine A postconditioning (71)."
    ),
    regions = "France",
    notes = paste(
      "Table 1 'CK and cTnI' column. Model built on a random 2/3 learning",
      "subset (84 patients) and externally validated on the remaining",
      "1/3 (48). Table 1 prints sex as 'female/male 101/31', which is",
      "almost certainly transposed (a STEMI cohort is predominantly male),",
      "so sex_female_pct is left unset. Samples were taken before PCI",
      "(time zero), every 4 h for 24 h, then every 6 h for 48 h."
    )
  )

  ini({
    # Release amount (concentration units, IU/L), with covariates
    lbtot <- log(2120)
    label("Total CK released by the lesion at AAR = 0, control arm, BTOT (IU/L)") # Table 2 CK 'B TOT' = 2120 (RSE 23%)
    e_trt_btot <- -0.43
    label("Additive effect of conditioning arm on log(BTOT) (unitless)") # Table 2 CK 'Arm on BTOT' = -0.43 (RSE 34%)
    e_aar_btot <- 0.020
    label("Exponential effect of area at risk on BTOT (1/% ACS)") # Table 2 CK 'AAR on BTOT' = 0.020 (RSE 30%)
    lb0 <- log(173)
    label("CK level at the first blood sample, initial condition B0 (IU/L)") # Table 2 CK 'B1' = 173 (RSE 12%)

    # Two-transit release (Savic transit functions). Chain 1 has the longer
    # mean release time (n + 1) / ktr, hence fslow. The paper states log-normal
    # IIV for every parameter; Fr is encoded logit-normal so that 0 < Fr < 1
    # (log-normal Fr exceeds 1 for a sizeable share of subjects).
    logitfslow <- logit(0.15)
    label("Fraction of B released through the slower transit chain 1, Fr (fraction)") # Table 2 CK 'Fr' = 0.15 (RSE 31%)
    lntr1 <- log(2.5)
    label("Number of transit compartments of release chain 1, n1 (unitless)") # Table 2 CK 'n1' = 2.5 (RSE 7%)
    lktr1 <- log(0.68)
    label("Transit rate constant of release chain 1, ktr1 (1/h)") # Table 2 CK 'ktr1' = 0.68 (RSE 11%)
    lntr2 <- log(50.7)
    label("Number of transit compartments of release chain 2, n2 (unitless)") # Table 2 CK 'n2' = 50.7 (RSE 8%)
    lktr2 <- log(14.4)
    label("Transit rate constant of release chain 2, ktr2 (1/h)") # Table 2 CK 'ktr2' = 14.4 (RSE 7%)

    # Distribution / elimination (one compartment)
    lkin <- log(6.9)
    label("Zero-order baseline CK production, kprod (IU/L/h)") # Table 2 CK 'kprod' = 6.9 (RSE 11%)
    lkel <- log(0.062)
    label("First-order CK elimination rate constant, ke (1/h)") # Table 2 CK 'ke' = 0.062 (RSE 3%)

    # IIV: Table 2 reports omega as SDs; variance = omega^2. omega_n2 and
    # omega_ktr2 were fixed to 0 (Table 2) and are omitted.
    etalbtot ~ 0.4225 # Table 2 CK 'omega Btot' = 0.65 (RSE 8%) -> 0.65^2
    etalb0 ~ 1.21 # Table 2 CK 'omega B1' = 1.1 (RSE 8%) -> 1.1^2
    etalogitfslow ~ 0.6241 # Table 2 CK 'omega Fr' = 0.79 (RSE 20%) -> 0.79^2, logit scale
    etalntr1 ~ 0.2401 # Table 2 CK 'omega n1' = 0.49 (RSE 11%) -> 0.49^2
    etalktr1 ~ 0.6084 # Table 2 CK 'omega Ktr1' = 0.78 (RSE 10%) -> 0.78^2
    etalkin ~ 0.6889 # Table 2 CK 'omega Kprod' = 0.83 (RSE 10%) -> 0.83^2
    etalkel ~ 0.0529 # Table 2 CK 'omega Ke' = 0.23 (RSE 10%) -> 0.23^2

    # Combined additive + proportional residual error (Methods, 'Error model')
    addSd <- 9.3
    label("Additive residual error (IU/L)") # Table 2 CK 'sigma add' = 9.3 (RSE 15%)
    propSd <- 0.099
    label("Proportional residual error (fraction)") # Table 2 CK 'sigma prop' = 0.099 (RSE 3%)
  })

  model({
    # Total released amount and its pre-sampling part (Methods, 'Biomarker
    # release'; 'Covariates': ln(theta) shifted by beta_ARM, theta * exp(beta_AAR * AAR))
    btot <- exp(lbtot + etalbtot + e_trt_btot * TRT_CONDITIONING + e_aar_btot * AREA_AT_RISK_PCT)
    b0 <- exp(lb0 + etalb0)
    brel <- btot - b0

    fslow <- expit(logitfslow + etalogitfslow)
    ntr1 <- exp(lntr1 + etalntr1)
    ktr1 <- exp(lktr1 + etalktr1)
    ntr2 <- exp(lntr2)
    ktr2 <- exp(lktr2)
    kin <- exp(lkin + etalkin)
    kel <- exp(lkel + etalkel)

    # Savic transit release functions (paper reference 43), time from the
    # first blood sample; each integrates to 1 over t > 0.
    tpos <- t + 1e-10
    rel1 <- exp(log(ktr1) + ntr1 * log(ktr1 * tpos) - ktr1 * tpos - lgamma(ntr1 + 1))
    rel2 <- exp(log(ktr2) + ntr2 * log(ktr2 * tpos) - ktr2 * tpos - lgamma(ntr2 + 1))
    release <- brel * (fslow * rel1 + (1 - fslow) * rel2)

    # Eq. 3: dB/dt = R(t) + kprod - ke * B, with B(0) = B0
    d/dt(central) <- release + kin - kel * central
    central(0) <- b0

    cpk <- central
    cpk ~ add(addSd) + prop(propSd)
  })
}
