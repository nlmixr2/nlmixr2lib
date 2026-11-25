Eigenmann_2017_erlotinib_gefitinib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, patient-derived NSCLC LXFA677 xenograft).",
    "Tumor-growth-inhibition model with acquired resistance for the EGFR",
    "tyrosine kinase inhibitors erlotinib and gefitinib, fitted jointly to",
    "both drugs (Eigenmann 2017 acquired resistance model, Equations 1-3 and",
    "9-17). A one-compartment oral PK model with first-order absorption and",
    "elimination (gefitinib volume decreasing linearly with the mg/kg dose)",
    "drives a Koch-type exponential-then-linear growth of sensitive cells;",
    "drug kills sensitive cells into a three-stage damaged-cell transit",
    "chain, from whose last stage damaged cells either die or convert",
    "first-order into a drug-resistant population that grows by the same",
    "law with both growth rates scaled by beta. Tumor growth, transit and",
    "conversion parameters are shared by the two drugs; the drug potency k2",
    "is drug-specific, and the 100 mg/kg dose group carries its own",
    "conversion rate (both drugs) and erlotinib potency. An optional",
    "threshold kill of resistant cells, used only for the paper's",
    "simulations, is off by default."
  )
  reference <- paste(
    "Eigenmann MJ, Frances N, Lave T, Walz AC.",
    "PKPD modeling of acquired resistance to anti-cancer drug treatment.",
    "J Pharmacokinet Pharmacodyn. 2017;44(6):617-630.",
    "doi:10.1007/s10928-017-9553-x.",
    "PK parameters from Eigenmann MJ, Frances N, Hoffmann G, Lave T, Walz AC.",
    "Mol Cancer Ther. 2016;15(12):3110-3119. doi:10.1158/1535-7163.MCT-16-0076.",
    sep = " "
  )
  vignette <- "Eigenmann_2017_erlotinib_gefitinib_mouse"
  units <- list(time = "day", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    TRT_GEFITINIB = list(
      description = "Drug indicator: 1 = gefitinib, 0 = erlotinib",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the drug-specific PK parameters (ke, V, the dose effect on V",
        "and their IIV) and the drug potency k2 of the jointly fitted model.",
        "Vehicle-control animals carry no drug, so the value is immaterial."
      ),
      source_name = "drug (erlotinib / gefitinib dataset)"
    ),
    DOSE_HIGH = list(
      description = "1 = animal in the 100 mg/kg dose group; 0 = 6.25 or 25 mg/kg (or vehicle)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Eigenmann 2017 Covariates: 'high dose treatment (100 mg/kg) was",
        "identified as covariate for the transformation rate kSR ... and, in",
        "case of erlotinib, for the efficacy parameter k2', assessed as a",
        "Monolix categorical covariate (low/mid vs. high dose). Applies to",
        "kSR for both drugs and to k2 for erlotinib only."
      ),
      source_name = "high dose group"
    ),
    COHDOSE = list(
      description = "Weight-normalised oral dose of the animal's dose group",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Eigenmann 2017 Equation 3, Cp = Ap / (V - D * I): D is the",
        "weight-normalised dose (6.25, 25 or 100 mg/kg). Enters only through",
        "the gefitinib-specific dose effect on volume; I is fixed to 0 for",
        "erlotinib."
      ),
      source_name = "D"
    ),
    TUM_VOL = list(
      description = "Tumor volume at the start of the study (initial condition of the sensitive cells)",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Eigenmann 2017 Equation 10: S(0) = TV0, 'the tumor volume at the",
        "start of the study'. TV0 is not among the estimated parameters",
        "(Table 1: seven parameters, all listed in Table 2), so it is taken",
        "per animal from the data. The authors' simulation code (S3-S5)",
        "uses TV0 = 107 mm^3."
      ),
      source_name = "TV0"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "erlotinib or gefitinib", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "erlotinib or gefitinib", units = "ug", specimen = "plasma", verified = TRUE),
    cycling_cells = list(
      analyte = "drug-sensitive proliferating tumor cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    ),
    damaged_cells1 = list(
      analyte = "drug-damaged tumor cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "drug-damaged tumor cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "drug-damaged tumor cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    ),
    resistant_cells = list(
      analyte = "drug-resistant proliferating tumor cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # The drug-resistant tumor-cell pool of a Simeoni / Koch damaged-cell chain
  # has no registered canonical; the Stein-model `growth` state is a
  # different (exponential, closed-form) construct.
  paper_specific_compartments <- c("resistant_cells")

  population <- list(
    species = "mouse (patient-derived NSCLC xenograft LXFA677)",
    n_subjects = 56L,
    n_studies = 1L,
    disease_state = "Subcutaneous patient-derived non-small cell lung cancer xenograft (LXFA677)",
    dose_range = paste(
      "Oral erlotinib or gefitinib 6.25, 25 or 100 mg/kg once daily, or vehicle;",
      "days 3-16 (erlotinib high dose: days 3-10 and 14-16 after a drug holiday",
      "for dose-limiting toxicity)"
    ),
    regions = "Roche Innovation Center Basel, Switzerland (preclinical)",
    notes = paste(
      "Seven groups of 8 mice (vehicle; erlotinib and gefitinib at low, mid and",
      "high dose; Eigenmann 2017 Methods and Supplementary Table S1). Tumor",
      "volume by caliper every 2-3 days over 30 days; sparse PK on days 10 and",
      "16. 7% of tumor observations missing (Supplementary Table S2); BLQ",
      "values were treated as missing. Fitted in Monolix 4.2.2 with the",
      "individual PK estimates of the earlier report (Eigenmann 2016) held",
      "fixed. Doses are in ug for a 25 g mouse (supplement code S3-S5:",
      "'C = 25 ; Conversion to flat dose in ug'), i.e. 1 mg/kg = 25 ug."
    )
  )

  ini({
    # --- PK: one-compartment oral model (Eigenmann 2017 Table 2) -------------
    # ka was fixed for both drugs; ke and V were estimated in the earlier PK
    # report (Eigenmann 2016) and carried into this analysis.
    lka <- fixed(log(55))
    label("First-order absorption rate constant ka, both drugs (1/day)") # Table 2: ka = 55.0 (X) erlotinib and gefitinib
    lkel <- log(7.56)
    label("Erlotinib elimination rate constant ke (1/day)") # Table 2: ke erlotinib = 7.56 (RSE 10%)
    lkel_gefitinib <- log(3.87)
    label("Gefitinib elimination rate constant ke (1/day)") # Table 2: ke gefitinib = 3.87 (RSE 10%)
    lvc <- log(0.127)
    label("Erlotinib apparent volume of distribution V (L)") # Table 2: V erlotinib = 0.127 (RSE 15%)
    lvc_gefitinib <- log(1.40)
    label("Gefitinib apparent volume of distribution V at zero dose (L)") # Table 2: V gefitinib = 1.40 (RSE 11%)
    e_cohdose_vc <- 0.00772
    label("Gefitinib linear decrease of V per mg/kg dose, I (L per mg/kg)") # Table 2: I gefitinib = 0.00772 (RSE 22%); erlotinib 0 (X)

    # --- Tumor growth (shared by both drugs; Table 2) ------------------------
    ltumorExpGrowth <- log(0.217)
    label("Exponential growth rate of sensitive cells lambda0 (1/day)") # Table 2: lambda0 = 0.217 (RSE 11%)
    ltumorLinGrowth <- log(42.8)
    label("Linear growth rate of sensitive cells lambda1 (mm^3/day)") # Table 2: lambda1 = 42.8 (RSE 10%)
    lgrowthRatio_resistant <- log(0.869)
    label("Ratio of resistant to sensitive growth rates beta (unitless)") # Table 2: beta = 0.869 (RSE 18%)

    # --- Drug effect and resistance (Table 2) --------------------------------
    ldamageTransit <- log(1.52)
    label("Damaged-cell transit rate constant k1, both drugs (1/day)") # Table 2: k1 = 1.52 (RSE 7%); tau = 3/k1 = 1.97 day
    lkmut <- log(7.4e-3)
    label("Conversion rate of damaged to resistant cells kSR, 6.25 and 25 mg/kg (1/day)") # Table 2: kSR = 7.4e-3 (RSE 8%)
    lkmut_highdose <- log(1.7e-3)
    label("Conversion rate of damaged to resistant cells kSR, 100 mg/kg (1/day)") # Table 2: kSR (high dose) = 1.7e-3 (RSE 30%), p = 2e-6
    ldrugSlope <- log(9.21e-5)
    label("Erlotinib potency k2, 6.25 and 25 mg/kg (L/ug/day)") # Table 2: k2 erlotinib = 9.21e-5 (RSE 23%)
    ldrugSlope_highdose <- log(3.15e-4)
    label("Erlotinib potency k2, 100 mg/kg (L/ug/day)") # Table 2: k2 (high dose) erlotinib = 3.15e-4 (RSE 17%), p = 1.1e-5
    ldrugSlope_gefitinib <- log(6.19e-4)
    label("Gefitinib potency k2, all doses (L/ug/day)") # Table 2: k2 gefitinib = 6.19e-4 (RSE 8%)

    # --- Threshold kill of resistant cells (simulation-only extension) -------
    # Not identifiable in the fit (Fig. 1 legend), so it is off by default.
    # The authors' simulations set k2R = 9.21e-5 L/ug/day (supplement S4/S5).
    drugSlope_resistant <- fixed(0)
    label("Potency on resistant cells above the threshold k2R (L/ug/day); 0 in the fitted model") # Fig. 1 legend; Eq. 17; S4/S5 code uses 0.0000921 for simulations
    lcthres <- fixed(log(7152))
    label("Total plasma threshold concentration for killing resistant cells Cth (ug/L)") # Eq. 18: 1 uM x 393.4 g/mol / fu 0.055; S4/S5 code 'Cth = 7152'

    # --- Between-animal variability (Table 2 'Random effect', Monolix SD) ----
    # Monolix reports omega as a standard deviation; variances are SD^2.
    etalkel ~ 0.110224 # Table 2: omega ke erlotinib = 0.332 (RSE 39%)
    etalkel_gefitinib ~ 0.123904 # Table 2: omega ke gefitinib = 0.352 (RSE 20%)
    etalvc ~ fixed(0.063001) # Table 2: omega V erlotinib = 0.251 (X)
    etalvc_gefitinib ~ fixed(0.00077284) # Table 2: omega V gefitinib = 0.0278 (X)
    etaltumorExpGrowth ~ 0.237169 # Table 2: omega lambda0 = 0.487 (RSE 18%)
    etaltumorLinGrowth ~ 0.426409 # Table 2: omega lambda1 = 0.653 (RSE 11%)
    etalgrowthRatio_resistant ~ 0.550564 # Table 2: omega beta = 0.742 (RSE 18%)
    etaldamageTransit ~ fixed(0.04) # Table 2: omega k1 = 0.2 (X)
    etalkmut ~ fixed(0.01) # Table 2: omega kSR = 0.1 (X)
    etaldrugSlope ~ 0.429025 # Table 2: omega k2 erlotinib = 0.655 (RSE 23%)
    etaldrugSlope_gefitinib ~ fixed(0.04) # Table 2: omega k2 gefitinib = 0.2 (X)

    # --- Residual error -------------------------------------------------------
    # Tumor volume: combined error (Methods); Table 1 prints a single value,
    # 0.136, which is the proportional term (an additive term in mm^3 of that
    # size would be negligible). The additive term is not reported.
    propSd_tumor_vol <- 0.136
    label("Proportional residual error on tumor volume (fraction)") # Table 1: residual error, acquired resistance model = 0.136
    addSd_tumor_vol <- fixed(0)
    label("Additive residual error on tumor volume (mm^3); not reported") # Methods: combined error model; additive value not reported
    # PK: proportional error (Methods); value from the earlier PK report, not
    # printed here.
    propSd <- fixed(0)
    label("Proportional residual error on plasma concentration (fraction); not reported") # Methods: PK residual error proportional; value not reported
  })

  model({
    # Drug selection: TRT_GEFITINIB = 1 for gefitinib, 0 for erlotinib.
    gef <- TRT_GEFITINIB
    erl <- 1 - TRT_GEFITINIB

    ka <- exp(lka)
    kel_erlotinib <- exp(lkel + etalkel)
    kel_gefitinib <- exp(lkel_gefitinib + etalkel_gefitinib)
    kel <- erl * kel_erlotinib + gef * kel_gefitinib
    vc_erlotinib <- exp(lvc + etalvc)
    vc_gefitinib <- exp(lvc_gefitinib + etalvc_gefitinib)
    # Eq. 3: Cp = Ap / (V - D * I), with I = 0 for erlotinib.
    v_eff <- erl * vc_erlotinib + gef * (vc_gefitinib - e_cohdose_vc * COHDOSE)

    tumorExpGrowth <- exp(ltumorExpGrowth + etaltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth + etaltumorLinGrowth)
    growthRatio_resistant <- exp(lgrowthRatio_resistant + etalgrowthRatio_resistant)
    # Eq. 15: beta = lambda0R / lambda0 = lambda1R / lambda1
    tumorExpGrowth_resistant <- growthRatio_resistant * tumorExpGrowth
    tumorLinGrowth_resistant <- growthRatio_resistant * tumorLinGrowth

    damageTransit <- exp(ldamageTransit + etaldamageTransit)
    # The high-dose group shares the random effect of the reference group.
    kmut <- exp(lkmut + (lkmut_highdose - lkmut) * DOSE_HIGH + etalkmut)
    drugSlope_erlotinib <- exp(ldrugSlope + (ldrugSlope_highdose - ldrugSlope) * DOSE_HIGH + etaldrugSlope)
    drugSlope_gefitinib <- exp(ldrugSlope_gefitinib + etaldrugSlope_gefitinib)
    drugSlope <- erl * drugSlope_erlotinib + gef * drugSlope_gefitinib
    cthres <- exp(lcthres)
    # Eq. 16: average time to eradication or conversion
    tau <- 3 / damageTransit

    # --- PK (Eqs. 1-3; amounts in ug) ---------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / v_eff

    # --- Tumor growth inhibition with acquired resistance (Eqs. 9-14, 17) ---
    # Eq. 17: resistant cells are killed only while Cp exceeds Cth.
    kill_resistant <- 0
    if (Cc > cthres) {
      kill_resistant <- drugSlope_resistant * Cc
    }

    d/dt(cycling_cells) <- 2 * tumorExpGrowth * tumorLinGrowth * cycling_cells /
      (2 * tumorExpGrowth * cycling_cells + tumorLinGrowth) -
      drugSlope * Cc * cycling_cells
    d/dt(damaged_cells1) <- drugSlope * cycling_cells * Cc - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3) - kmut * damaged_cells3
    d/dt(resistant_cells) <- 2 * tumorExpGrowth_resistant * tumorLinGrowth_resistant * resistant_cells /
      (2 * tumorExpGrowth_resistant * resistant_cells + tumorLinGrowth_resistant) +
      kmut * damaged_cells3 - kill_resistant * resistant_cells

    cycling_cells(0) <- TUM_VOL
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0
    resistant_cells(0) <- 0

    # Eq. 9: total tumor volume
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3 + resistant_cells
    frac_resistant <- resistant_cells / tumor_vol

    Cc ~ prop(propSd)
    tumor_vol ~ add(addSd_tumor_vol) + prop(propSd_tumor_vol)
  })
}
