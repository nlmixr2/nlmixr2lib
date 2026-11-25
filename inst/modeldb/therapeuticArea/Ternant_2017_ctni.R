Ternant_2017_ctni <- function() {
  description <- paste(
    "Population kinetic model of serum cardiac troponin I (cTnI) released by",
    "the injured myocardium after ST-segment elevation myocardial",
    "infarction (STEMI) treated by primary percutaneous coronary",
    "intervention (Ternant 2017). The cTnI level at the first blood sample",
    "(time zero, before PCI) is an estimated initial condition B0; the",
    "remainder of the total released amount (B = BTOT - B0) enters the",
    "central compartment through the sum of two Savic transit-release",
    "functions (fraction Fr through chain 1), with first-order",
    "elimination from the central compartment and first-order exchange",
    "with a peripheral compartment (two compartments, no baseline",
    "production). BTOT",
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
    concentration = "mg/L"
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
        "For the cTnI data (PC, PC-CsA and PCNR trials) 1 = ischemic",
        "postconditioning (PC, PCNR) or cyclosporine A pharmacological",
        "postconditioning (PC-CsA). Effect on log(BTOT): additive",
        "beta = -0.56 (Methods, 'Covariates'; Table 2 'Arm on BTOT')."
      ),
      source_name = "ARM"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "cardiac troponin I",
      units = "mg/L",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "cardiac troponin I",
      units = "mg/L",
      specimen = "not applicable",
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
    # Release amount (concentration units as printed, mg/L), with covariates
    lbtot <- log(196)
    label("Total cTnI released by the lesion at AAR = 0, control arm, BTOT (mg/L)") # Table 2 cTnI 'B TOT' = 196 (RSE 29%)
    e_trt_btot <- -0.56
    label("Additive effect of conditioning arm on log(BTOT) (unitless)") # Table 2 cTnI 'Arm on BTOT' = -0.56 (RSE 38%)
    e_aar_btot <- 0.019
    label("Exponential effect of area at risk on BTOT (1/% ACS)") # Table 2 cTnI 'AAR on BTOT' = 0.019 (RSE 39%)
    lb0 <- log(0.11)
    label("Central cTnI level at the first blood sample, initial condition B1(0) (mg/L)") # Table 2 cTnI 'B1' = 0.11 (RSE 41%)

    # Two-transit release (Savic transit functions). Chain 1 has the longer
    # mean release time (n + 1) / ktr, hence fslow. The paper states log-normal
    # IIV for every parameter; Fr is encoded logit-normal so that 0 < Fr < 1
    # (log-normal Fr exceeds 1 for a sizeable share of subjects).
    logitfslow <- logit(0.53)
    label("Fraction of B released through the slower transit chain 1, Fr (fraction)") # Table 2 cTnI 'Fr' = 0.53 (RSE 14%)
    lntr1 <- log(2.0)
    label("Number of transit compartments of release chain 1, n1 (unitless)") # Table 2 cTnI 'n1' = 2.0 (RSE 10%)
    lktr1 <- log(0.29)
    label("Transit rate constant of release chain 1, ktr1 (1/h)") # Table 2 cTnI 'ktr1' = 0.29 (RSE 6%)
    lntr2 <- log(4.2)
    label("Number of transit compartments of release chain 2, n2 (unitless)") # Table 2 cTnI 'n2' = 4.2 (RSE 6%)
    lktr2 <- log(0.88)
    label("Transit rate constant of release chain 2, ktr2 (1/h)") # Table 2 cTnI 'ktr2' = 0.88 (RSE 1%)

    # Distribution / elimination (two compartments; no kprod for cTnI)
    lkel <- log(0.097)
    label("First-order cTnI elimination rate constant from central, ke (1/h)") # Table 2 cTnI 'ke' = 0.097 (RSE 6%)
    lk12 <- log(0.11)
    label("Central-to-peripheral rate constant, k12 (1/h)") # Table 2 cTnI 'k12' = 0.11 (RSE 8%)
    lk21 <- log(0.044)
    label("Peripheral-to-central rate constant, k21 (1/h)") # Table 2 cTnI 'k21' = 0.044 (RSE 12%)

    # IIV: Table 2 reports omega as SDs; variance = omega^2. No omega is
    # reported for ktr2 (Table 2 '-') and it is omitted.
    etalbtot ~ 0.64 # Table 2 cTnI 'omega Btot' = 0.80 (RSE 8%) -> 0.80^2
    etalb0 ~ 13.69 # Table 2 cTnI 'omega B1' = 3.7 (RSE 18%) -> 3.7^2
    etalogitfslow ~ 0.3481 # Table 2 cTnI 'omega Fr' = 0.59 (RSE 13%) -> 0.59^2, logit scale
    etalntr1 ~ 0.2116 # Table 2 cTnI 'omega n1' = 0.46 (RSE 11%) -> 0.46^2
    etalktr1 ~ 0.1849 # Table 2 cTnI 'omega Ktr1' = 0.43 (RSE 14%) -> 0.43^2
    etalntr2 ~ 0.0625 # Table 2 cTnI 'omega n2' = 0.25 (RSE 27%) -> 0.25^2
    etalkel ~ 0.0529 # Table 2 cTnI 'omega Ke' = 0.23 (RSE 15%) -> 0.23^2
    etalk12 ~ 0.0256 # Table 2 cTnI 'omega 12' = 0.16 (RSE 30%) -> 0.16^2
    etalk21 ~ 0.6241 # Table 2 cTnI 'omega 21' = 0.79 (RSE 14%) -> 0.79^2

    # Combined additive + proportional residual error (Methods, 'Error model')
    addSd <- 0.82
    label("Additive residual error (mg/L)") # Table 2 cTnI 'sigma add' = 0.82 (RSE 20%)
    propSd <- 0.15
    label("Proportional residual error (fraction)") # Table 2 cTnI 'sigma prop' = 0.15 (RSE 4%)
  })

  model({
    # Total released amount and its pre-sampling part (Methods, 'Biomarker
    # release'; 'Covariates': ln(theta) shifted by beta_ARM, theta * exp(beta_AAR * AAR))
    btot <- exp(lbtot + etalbtot + e_trt_btot * TRT_CONDITIONING + e_aar_btot * AAR_PCT)
    b0 <- exp(lb0 + etalb0)
    brel <- btot - b0

    fslow <- expit(logitfslow + etalogitfslow)
    ntr1 <- exp(lntr1 + etalntr1)
    ktr1 <- exp(lktr1 + etalktr1)
    ntr2 <- exp(lntr2 + etalntr2)
    ktr2 <- exp(lktr2)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # Savic transit release functions (paper reference 43), time from the
    # first blood sample; each integrates to 1 over t > 0.
    tpos <- t + 1e-10
    rel1 <- exp(log(ktr1) + ntr1 * log(ktr1 * tpos) - ktr1 * tpos - lgamma(ntr1 + 1))
    rel2 <- exp(log(ktr2) + ntr2 * log(ktr2 * tpos) - ktr2 * tpos - lgamma(ntr2 + 1))
    release <- brel * (fslow * rel1 + (1 - fslow) * rel2)

    # Eqs. 4-5 (sign of the k21 term in Eq. 5 corrected to -k21 * B2):
    # B1(0) = B0; B2(0) was not estimated (Table 2 'B2' = '-') and is 0.
    d/dt(central) <- release - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    central(0) <- b0

    troponin_i <- central
    troponin_i ~ add(addSd) + prop(propSd)
  })
}
