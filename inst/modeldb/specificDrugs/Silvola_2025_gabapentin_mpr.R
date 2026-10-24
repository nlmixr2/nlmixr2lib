Silvola_2025_gabapentin_mpr <- function() {
  description <- paste(
    "One-compartment oral population PK model of gabapentin in immediately",
    "postpartum lactating women after cesarean delivery (Model A of Silvola",
    "2025), with first-order absorption and linear elimination. Breast-milk",
    "concentration is an algebraic observable equal to an estimated",
    "milk-to-plasma concentration ratio times the plasma concentration (no milk",
    "compartment). The absorption rate constant is parameterised as a multiple",
    "of each subject's own elimination rate constant (ka = karatio * CL/V,",
    "karatio bounded >= 1 to exclude flip-flop kinetics), so ka inherits the",
    "between-subject variability of CL/F and V/F. No covariates were retained.",
    "Sibling Silvola_2025_gabapentin_milkcmt is the authors' Model B, which",
    "replaces the ratio with an explicit milk compartment.",
    sep = " "
  )
  reference <- paste(
    "Silvola R, O'Kane A, Heathman M, Marotta H, Trussel H, Ray B, Dowden S,",
    "Masters AR, Haas DM, Quinney SK. Population Pharmacokinetics and Transfer",
    "of Gabapentin When Used as a Pain Adjunct for Cesarean Deliveries. CPT",
    "Pharmacometrics Syst Pharmacol. 2025;14(3):551-560.",
    "doi:10.1002/psp4.13295. PMCID: PMC11919262.",
    sep = " "
  )
  vignette <- "Silvola_2025_gabapentin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "gabapentin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Maternal body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Explored graphically but not retained: 'the addition of covariates to",
        "the model led to model instability and reduced reliability of the",
        "model' (Methods 2.6; Discussion names body weight explicitly).",
        sep = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Explored graphically but not retained (Results 3.1, 'weight/BMI')."
    ),
    CRCL = list(
      description = "eGFR by the 2021 CKD-EPI creatinine equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Source column GFR. Explored graphically but not retained (Discussion)."
    ),
    TPP = list(
      description = "Time after delivery",
      units = "weeks",
      type = "continuous",
      notes = paste(
        "CL/F and V/F showed non-significant decreases with time after birth",
        "(linear regression p = 0.07, Figure S3) but the effect could not be",
        "formally assessed because of model instability (Results 3.1).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = "125 plasma + 35 breast-milk concentrations (15 women contributed milk)",
    age_range = "19-41 years",
    age_median = "29 years",
    weight_range = "68.5-154 kg",
    weight_median = "107.5 kg",
    sex_female_pct = 100,
    race_ethnicity = paste(
      "Caucasian 48%, Black/African 48%, Native Hawaiian/Pacific Islander 5%;",
      "Hispanic 14% (Table 1, n = 21 enrolled)",
      sep = " "
    ),
    disease_state = paste(
      "Healthy singleton pregnancies, immediately postpartum after elective or",
      "urgent cesarean delivery under an enhanced-recovery-after-surgery",
      "(ERAS) protocol, planning to breastfeed",
      sep = " "
    ),
    dose_range = paste(
      "300 mg oral gabapentin before cesarean delivery, then 300 mg every 6 h",
      "for 48 h after delivery. One woman taking 600 mg every 6 h was excluded",
      "from model development.",
      sep = " "
    ),
    regions = "United States (Indianapolis, Indiana; single centre)",
    lactation_stage = "colostrum, first 48 h postpartum",
    renal_function = "eGFR (2021 CKD-EPI) 151 (61-218) mL/min/1.73 m^2",
    notes = paste(
      "Prospective observational study NCT05099484 (February-August 2021).",
      "Table 1 demographics are for the 21 enrolled women: BMI 35.8",
      "(26.2-48.2), lean body mass 55.2 (39.8-71.8) kg, height 1.7 (1.5-1.8)",
      "m, estimated gestational age at delivery 39.1 (31.1-41.3) weeks. Plasma",
      "sampled pre-dose and 0.25, 0.5, 1, 1.5, 2, 4 and 6 h after one",
      "postpartum dose; milk at least every 2 h over the same interval. Fit in",
      "NONMEM 7.5 (ADVAN13, FOCE-I).",
      sep = " "
    )
  )

  ini({
    lkaratio <- c(0, log(7.828)); label("Ratio of the absorption to the elimination rate constant, ka/kel (unitless)")
    # Data S1 $PK: KE = CL/V2; TVKA = THETA(1)*KE; $THETA (1, 7.67) ;KA -- lower
    # bound 1 (log scale 0) keeps ka above kel (Results 3.1). The final THETA(1)
    # is not printed; Table 2 reports the derived 'KA (h-1)' = 0.43 (RSE 67%),
    # so THETA(1) = 0.43 / (15.6 / 284) = 7.828.
    lcl <- log(15.6); label("Apparent oral clearance, CL/F (L/h)")
    # Table 2 'CL/F (L/h)' = 15.6 (RSE 15%)
    lvc <- log(284); label("Apparent central volume of distribution, V2/F (L)")
    # Table 2 'V2/F (L)' = 284 (RSE 31%)
    lcmpr <- log(0.401); label("Milk-to-plasma concentration ratio (unitless)")
    # Table 2 'MPR' = 0.401 (RSE 12%)

    # Between-subject variability: exponential, Table 2 omega^2 are variances.
    # Data S1 carries '$OMEGA 0 FIX ;IIV/BSV KA'; BSV on ka is absent here
    # because the authors kept BSV on CL/F and V/F 'but not ka' (Results 3.1).
    etalcl ~ 0.0249 # Table 2 omega^2 CL/F = 0.0249 (RSE 51%, shrinkage 35.11%)
    etalvc ~ 0.163 # Table 2 omega^2 V2/F = 0.163 (RSE 49%, shrinkage 6.19%)
    etalcmpr ~ 0.153 # Table 2 omega^2 MPR = 0.153 (RSE 43%, shrinkage 25.64%)

    # Residual error: proportional on both outputs, Data S1 $ERROR
    # Y = IPRED*(1 + EPS(1)) plasma and Y = IPRED*(1 + EPS(2)) milk.
    propSd <- 0.11314; label("Proportional residual SD, plasma (fraction)")
    # Table 2 'sigma^2 plasma' = 0.0128 (variance); sqrt(0.0128) = 0.11314
    propSd_Cmilk <- 0.20100; label("Proportional residual SD, breast milk (fraction)")
    # Table 2 'sigma^2 breast milk' = 0.0404 (variance); sqrt(0.0404) = 0.20100
  })

  model({
    # Individual parameters (Data S1 $PK). No covariates in the final model.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc
    karatio <- exp(lkaratio)
    ka <- karatio * kel
    cmpr <- exp(lcmpr + etalcmpr)

    # Data S1 $DES: DADT(1) = -KA*A(1); DADT(2) = KA*A(1) - KE*A(2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Plasma concentration in mg/L (dose mg, volume L; the source reports
    # ng/mL, which is 1000 x mg/L).
    Cc <- central / vc
    # Data S1 $ERROR: IPRED(milk) = (A(2)/S2)*MPR, i.e. a concentration ratio.
    Cmilk <- cmpr * Cc

    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
