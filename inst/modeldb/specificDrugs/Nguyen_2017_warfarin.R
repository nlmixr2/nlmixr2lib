Nguyen_2017_warfarin <- function() {
  description <- "Warfarin PK/PD simulation model from the ISoP Model Evaluation Group tutorial on goodness-of-fit graphics (Nguyen 2017, pharmacokinetic-pharmacodynamic case example). Total warfarin after a single oral dose follows a one-compartment model with first-order absorption (parameterised by an absorption half-life), an absorption lag time and first-order elimination, with CL and V allometrically scaled to 70 kg (fixed exponents 0.75 and 1). Warfarin inhibits the zero-order production of prothrombin complex activity (PCA) in a turnover model with complete maximal inhibition, first-order PCA loss parameterised by its half-life Teq, and PCA at its drug-free baseline before dosing. These are the authors' data-generating (true) values from Supplementary Table S2, used to simulate the O'Reilly warfarin study design (32 subjects, 1.5 mg/kg); the misspecified models fitted to the simulated data have no printed estimates and are not packaged."
  reference <- paste(
    "Nguyen THT, Mouksassi M-S, Holford N, Al-Huniti N, Freedman I, Hooker AC, John J,",
    "Karlsson MO, Mould DR, Perez Ruixo JJ, Plan EL, Savic R, van Hasselt JGC, Weber B,",
    "Zhou C, Comets E, Mentre F; Model Evaluation Group of the International Society of",
    "Pharmacometrics (ISoP) Best Practice Committee. Model Evaluation of Continuous Data",
    "Pharmacometric Models: Metrics and Graphics. CPT Pharmacometrics Syst Pharmacol.",
    "2017;6(2):87-109. doi:10.1002/psp4.12161. PMID 27884052. PMC5321813."
  )
  vignette <- "Nguyen_2017_model_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    PCA = list(analyte = "prothrombin complex", units = "% of normal activity", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL (exponent 0.75, fixed) and V (exponent 1, fixed) with reference weight 70 kg (Nguyen 2017 Eq. 25). Simulated cohort weights 40-102 kg (O'Reilly study design).",
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    weight_range = "40-102 kg",
    sex_female_pct = 100 * 5 / 32,
    disease_state = "Simulated subjects reproducing the design of the O'Reilly warfarin studies (healthy volunteers); the data analysed in the tutorial are simulated from this model, not observed.",
    dose_range = "Single oral dose of 1.5 mg/kg racemic warfarin (60-153 mg in the deposited simulated dataset).",
    regions = "United States (design of the O'Reilly studies)",
    notes = paste(
      "Nguyen 2017 PK/PD case example, 'Study design and data simulation': 27 male and 5 female",
      "patients, body weight 40-102 kg, doses calculated at 1.5 mg/kg, total warfarin",
      "concentrations (mg/L) and PCA (PCA units) as in O'Reilly et al. 1963 and O'Reilly and",
      "Aggeler 1968. Data were simulated in NONMEM 7.3 and the PD models fitted in Phoenix NLME 7.0",
      "with the PK parameters fixed."
    )
  )

  ini({
    lcl <- fixed(log(0.135)); label("Clearance at 70 kg (L/h)") # Supp Table S2 'muCL (L/h)' = 0.135
    lvc <- fixed(log(8.03)); label("Volume of distribution at 70 kg (L)") # Supp Table S2 'muV (L)' = 8.03
    ltabs <- fixed(log(0.588)); label("Absorption half-life (h)") # Supp Table S2 'muTabs (h)' = 0.588
    ltlag <- fixed(log(0.842)); label("Absorption lag time (h)") # Supp Table S2 'muTlag (h)' = 0.842
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/70 on CL (unitless)") # Supp Table S2 'beta Weight on CL' = 0.75 (fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent of WT/70 on V (unitless)") # Supp Table S2 'beta Weight on V' = 1 (fixed)

    lrbase <- fixed(log(96.6)); label("Baseline PCA0 (% of normal activity)") # Supp Table S2 'muPCA0 (PCA units)' = 96.6
    limax <- fixed(log(1)); label("Maximal fractional inhibition of PCA production (unitless)") # Supp Table S2 'muEMAX (-)' = -1 (fixed); magnitude 1 = complete inhibition, sign per vignette
    lic50 <- fixed(log(1.17)); label("Warfarin concentration giving half-maximal inhibition C50 (mg/L)") # Supp Table S2 'muC50 (mg/L)' = 1.17
    lthalf <- fixed(log(13.0)); label("PCA elimination half-life Teq (h)") # Supp Table S2 'muTEQ (h)' = 13.0

    etalcl ~ fixed(0.0703) # Supp Table S2 omega2 CL = 0.0703
    etalvc ~ fixed(0.0195) # Supp Table S2 omega2 V = 0.0195
    etaltabs ~ fixed(0.551) # Supp Table S2 omega2 Tabs = 0.551
    etaltlag ~ fixed(0.347) # Supp Table S2 omega2 Tlag = 0.347
    etalrbase ~ fixed(0.00286) # Supp Table S2 omega2 PCA0 = 0.00286
    etalic50 ~ fixed(0.2) # Supp Table S2 omega2 C50 = 0.2
    etalthalf ~ fixed(0.0114) # Supp Table S2 omega2 TEQ = 0.0114

    addSd <- fixed(sqrt(0.094)); label("Additive residual SD for warfarin (mg/L)") # Supp Table S2 'sigma2 inter,PK (mg/L)2' = 0.094
    propSd <- fixed(sqrt(0.00677)); label("Proportional residual SD for warfarin (fraction)") # Supp Table S2 'sigma2 slope,PK' = 0.00677
    addSd_PCA <- fixed(sqrt(14.6)); label("Additive residual SD for PCA (% of normal activity)") # Supp Table S2 'sigma2 inter,PD (PCA units)' = 14.6
  })

  model({
    # Eq. 25: theta_i = mu * (Weight/70)^beta * exp(eta_i)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    tabs <- exp(ltabs + etaltabs)
    tlag <- exp(ltlag + etaltlag)
    rbase <- exp(lrbase + etalrbase)
    imax <- exp(limax)
    ic50 <- exp(lic50 + etalic50)
    thalf <- exp(lthalf + etalthalf)

    ka <- log(2) / tabs
    kel <- cl / vc
    kout <- log(2) / thalf
    kin <- rbase * kout

    # Eq. 24
    d/dt(depot) <- -ka * depot
    alag(depot) <- tlag
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    d/dt(PCA) <- kin * (1 - imax * Cc / (Cc + ic50)) - kout * PCA
    PCA(0) <- rbase

    Cc ~ add(addSd) + prop(propSd) + combined2()
    PCA ~ add(addSd_PCA)
  })
}
