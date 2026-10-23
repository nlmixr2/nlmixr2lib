Silvola_2025_gabapentin_milkcmt <- function() {
  description <- paste(
    "One-compartment oral population PK model of gabapentin in immediately",
    "postpartum lactating women after cesarean delivery (Model B of Silvola",
    "2025), with first-order absorption, linear elimination and a reversible",
    "breast-milk compartment. Drug moves between plasma and milk with",
    "first-order rate constants and is not eliminated from milk. The milk",
    "volume (30 mL, typical colostrum volume of the first 2 days postpartum)",
    "and the central-to-milk rate constant (1/h) were held at fixed values for",
    "identifiability; the milk-to-central rate constant was estimated with",
    "between-subject variability. The absorption rate constant is",
    "parameterised as a multiple of each subject's own elimination rate",
    "constant (ka = karatio * CL/V, karatio bounded >= 1 to exclude flip-flop",
    "kinetics). No covariates were retained. Sibling",
    "Silvola_2025_gabapentin_mpr is the authors' Model A, which replaces the",
    "milk compartment with a milk-to-plasma concentration ratio.",
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
    central = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE),
    milk = list(analyte = "gabapentin", units = "mg", specimen = "milk", verified = TRUE)
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
    lkaratio <- c(0, log(9.904)); label("Ratio of the absorption to the elimination rate constant, ka/kel (unitless)")
    # Data S2 $PK: KE = CL/V2; TVKA = THETA(1)*KE (same form as Model A, whose
    # $THETA gives THETA(1) a lower bound of 1). The published Data S2 stream
    # stops inside $ERROR, so the final THETA(1) is not printed; Table 3 reports
    # the derived 'KA (h-1)' = 0.51 (RSE 37%), so THETA(1) = 0.51 / (15.5 / 301)
    # = 9.904.
    lcl <- log(15.5); label("Apparent oral clearance, CL/F (L/h)")
    # Table 3 'CL/F (L/h)' = 15.5 (RSE 11%)
    lvc <- log(301); label("Apparent central volume of distribution, V2/F (L)")
    # Table 3 'V2/F (L)' = 301 (RSE 20%)
    lvmilk <- fixed(log(0.03)); label("Breast-milk compartment volume, V3/F (L)")
    # Table 3 'V3/F (L)' = 0.03, not estimated: 30 mL, the typical milk volume
    # of the first 2 days postpartum (Methods 2.6)
    lk_central_milk <- fixed(log(1)); label("Central-to-milk transfer rate constant, K23 (1/h)")
    # Table 3 'K23 (h-1)' = 1, not estimated (Methods 2.6, to make K32 identifiable)
    lk_milk_central <- log(24800); label("Milk-to-central transfer rate constant, K32 (1/h)")
    # Table 3 'K32 (h-1)' = 24,800 (RSE 21%)

    # Between-subject variability: exponential, Table 3 omega^2 are variances.
    # No BSV on ka (Results 3.1); Data S2 places ETA(4) on K32.
    etalcl ~ 0.029 # Table 3 omega^2 CL/F = 0.029 (RSE 44%, shrinkage 28.9%)
    etalvc ~ 0.106 # Table 3 omega^2 V2/F = 0.106 (RSE 95%, shrinkage 12.6%)
    etalk_milk_central ~ 0.328 # Table 3 omega^2 K32 = 0.328 (RSE 58%, shrinkage 25.5%)

    # Residual error: proportional on both outputs (Data S2 $PROBLEM tag
    # '*prop*', same $ERROR form as the complete Model A stream in Data S1).
    propSd <- 0.11747; label("Proportional residual SD, plasma (fraction)")
    # Table 3 'sigma^2 plasma' = 0.0138 (variance); sqrt(0.0138) = 0.11747
    propSd_Cmilk <- 0.19799; label("Proportional residual SD, breast milk (fraction)")
    # Table 3 'sigma^2 breast milk' = 0.0392 (variance); sqrt(0.0392) = 0.19799
  })

  model({
    # Individual parameters (Data S2 $PK). No covariates in the final model.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc
    karatio <- exp(lkaratio)
    ka <- karatio * kel
    vmilk <- exp(lvmilk)
    k_central_milk <- exp(lk_central_milk)
    k_milk_central <- exp(lk_milk_central + etalk_milk_central)

    # Data S2 $DES:
    #   DADT(1) = -KA*A(1)
    #   DADT(2) = KA*A(1) - KE*A(2) + K32*A(3) - K23*A(2)
    #   DADT(3) = K23*A(2) - K32*A(3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central + k_milk_central * milk - k_central_milk * central
    d/dt(milk) <- k_central_milk * central - k_milk_central * milk

    # Concentrations in mg/L (dose mg, volumes L; the source reports ng/mL,
    # which is 1000 x mg/L: Data S2 S2 = V2/1000 and S3 = V3/1000).
    Cc <- central / vc
    Cmilk <- milk / vmilk

    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
