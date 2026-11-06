Tan_2025_adomeglivant_igrh <- function() {
  description <- paste(
    "Population PK/PD model for adomeglivant (LY2409021, oral glucagon",
    "receptor antagonist) in adults with type 2 diabetes (Tan 2025,",
    "Integrated Glucose-RBC-HbA1c (IGRH) model). One-compartment PK with",
    "first-order absorption and body-weight effects on CL/F and V/F, driving",
    "an inhibitory Imax effect on the zero-order input of average glucose",
    "(AG, indirect response) together with a fractional placebo effect on",
    "the same input. Red blood cells age through 12 non-glycated and 12",
    "glycated transit compartments; the RBC lifespan (120 days at the",
    "149 mg/dL reference AG) shortens with a power of AG, glycation is",
    "first-order in AG, and a fraction of the hemoglobin is already",
    "glycated when RBC precursors are released. HbA1c is the glycated",
    "percentage of the 24 compartments. One of four alternative",
    "semi-mechanistic HbA1c models fit sequentially to the same 4-week",
    "Phase 1b data (siblings: Tan_2025_adomeglivant_adopt, _ffh, _fhh)."
  )
  reference <- paste(
    "Tan AYF, Schneck K, Garhyan P, Chan ECY, Tham LS.",
    "Evaluation of Four Semi-Mechanistic Models for Predicting Glycemic",
    "Control With a Glucagon Receptor Antagonist in People With Type 2",
    "Diabetes. CPT Pharmacometrics Syst Pharmacol. 2025;14:1381-1390.",
    "doi:10.1002/psp4.70058.",
    "IGRH structure: Lledo-Garcia R, Mazer NA, Karlsson MO.",
    "J Pharmacokinet Pharmacodyn. 2013;40:129-142.",
    "doi:10.1007/s10928-012-9289-6; equations as deposited in Kunina H et al.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11:1443-1457 (supplement S6).",
    "doi:10.1002/psp4.12854."
  )
  vignette <- "Tan_2025_adomeglivant_hba1c"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  paper_specific_compartments <- c(
    sprintf("hba_n%02d", 1:12),
    sprintf("hba_g%02d", 1:12)
  )

  rbcNon <- list(analyte = "non-glycated hemoglobin", units = "g/L", specimen = "blood cell", verified = TRUE)
  rbcGly <- list(analyte = "glycated hemoglobin (HbA1c)", units = "g/L", specimen = "blood cell", verified = TRUE)
  compartmentData <- c(
    list(
      depot = list(analyte = "adomeglivant", units = "mg", specimen = "administration site", verified = TRUE),
      central = list(analyte = "adomeglivant", units = "mg", specimen = "plasma", verified = TRUE),
      glucose = list(analyte = "glucose (24-h average)", units = "mmol/L", specimen = "plasma", verified = TRUE)
    ),
    stats::setNames(rep(list(rbcNon), 12), sprintf("hba_n%02d", 1:12)),
    stats::setNames(rep(list(rbcGly), 12), sprintf("hba_g%02d", 1:12))
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value; reference 70 kg (Table S2 footnotes c/d). Exponential effect on CL/F, power effect on V/F.",
      source_name = "Baseline weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47L,
    n_studies = 1L,
    age_range = "21-70 years (eligibility)",
    age_mean = "53.6 years (SD 9.40)",
    weight_mean = "82.2 kg (SD 15.2)",
    sex_female_pct = 17.0,
    race_ethnicity = c(White = 46.8),
    disease_state = paste(
      "Type 2 diabetes treated with diet and exercise alone or a stable dose",
      "of metformin (68.1% on metformin); BMI 20-40 kg/m^2; HbA1c 6.5-10.0%."
    ),
    dose_range = "Adomeglivant 5, 30, 60 or 90 mg orally once daily for 28 days, or placebo (3:1)",
    regions = "Not reported",
    notes = paste(
      "Model-building data: Phase 1b study NCT01606397 (Table 1, 47 patients;",
      "928 PK, 78 HbA1c, 138 average-glucose records). Average glucose is",
      "the AUC of glucose over the sampling day divided by its duration.",
      "External evaluation: 24-week Phase 2b study NCT01241448 (254",
      "patients; placebo, 2.5, 10 or 20 mg once daily)."
    )
  )

  ini({
    # ---- PK (Table S2) ----
    lka <- log(0.385); label("First-order absorption rate constant KA (1/h)") # Table S2 KA = 0.385
    lcl <- log(0.283); label("Apparent clearance CL/F at 70 kg (L/h)") # Table S2 CL/F = 0.283
    lvc <- log(24.3); label("Apparent volume V/F at 70 kg (L)") # Table S2 V/F = 24.3
    e_wt_cl <- 0.0116; label("Exponential body-weight coefficient on CL/F, exp(e * (WT - 70)) (1/kg)") # Table S2 footnote c = 0.0116
    e_wt_vc <- 0.956; label("Power exponent of (WT/70) on V/F (unitless)") # Table S2 footnote d = 0.956

    # ---- Average glucose turnover (Table S6) ----
    lrbase_glucose <- log(10.1); label("Baseline average glucose AG0 (mmol/L)") # Table S6 AG0 = 10.1
    lkout_glucose <- log(0.164); label("First-order loss rate constant of average glucose KOUTG (1/h)") # Table S6 KOUTG = 0.164
    plac <- 0.122; label("Fractional placebo inhibition of the glucose input PLAC (fraction)") # Table S6 PLAC = 0.122
    limax <- log(0.151); label("Maximal fractional inhibition of the glucose input IMAX (fraction)") # Table S6 IMAX = 0.151
    lic50 <- log(1450); label("Adomeglivant concentration at half-maximal inhibition IC50 (ug/L)") # Table S6 IC50 = 1450

    # ---- Red blood cell / HbA1c model (Table S6) ----
    lkin_hb <- fixed(log(1.2)); label("Zero-order release constant of RBC hemoglobin KINHB (g/L/day)") # Table S6 KINHB = 1.2 Fixed (footnote c: sensitivity analysis)
    lkg <- log(6.76e-6); label("Glycosylation rate constant KGL (1/day per mg/dL)") # Table S6 KGL = 6.76 x 10^-6
    lls <- fixed(log(120)); label("Red blood cell lifespan LS at the reference average glucose 149 mg/dL (day)") # Table S6 LS = 120 Fixed (footnote d: physiological value)
    llsp <- log(7.12); label("Red blood cell precursor lifespan LSP (day)") # Table S6 LSP = 7.12
    # Table S6 reports GAMMA = 0.427 as the exponent of (reference AG / AG);
    # it is carried here as the equivalent exponent of (AG / reference AG).
    gamma_ls <- -0.427; label("Power exponent of (AG / 149 mg/dL) on the RBC lifespan (unitless)") # Table S6 GAMMA = 0.427, sign-flipped for (AG/149)

    # ---- IIV ----
    # Table S2 footnote b: CV% = 100 * sqrt(exp(OMEGA) - 1), so OMEGA = log(1 + CV^2).
    # Table S2 footnote e: correlation = COV / sqrt(OMEGA(CL) * OMEGA(V)).
    etalka ~ log(1 + 0.262^2) # Table S2 KA IIV 26.2 %CV
    etalcl + etalvc ~ c(
      log(1 + 0.180^2),
      0.778 * sqrt(log(1 + 0.180^2) * log(1 + 0.124^2)),
      log(1 + 0.124^2)
    ) # Table S2 CL/F IIV 18.0 %CV, V/F IIV 12.4 %CV, correlation 0.778
    etalrbase_glucose ~ log(1 + 0.185^2) # Table S6 AG0 IIV 18.5 %CV (exponential)
    # Table S6 footnote b: additive IIV on PLAC, CV% = 100 * sqrt(OMEGA) / THETA,
    # so OMEGA = (CV * THETA)^2.
    etaplac ~ (0.393 * 0.122)^2 # Table S6 PLAC IIV 39.3 %CV (additive)
    etallsp ~ log(1 + 0.523^2) # Table S6 LSP IIV 52.3 %CV (exponential)

    # ---- Residual error ----
    propSd <- 0.290; label("Proportional residual error, adomeglivant concentration (fraction)") # Table S2 proportional residual error 0.290
    propSd_glucose <- 0.0804; label("Proportional residual error, average glucose (fraction)") # Table S6 glucose residual error 0.0804
    propSd_Hba1c <- 0.0409; label("Proportional residual error, HbA1c (fraction)") # Table S6 HbA1c residual error 0.0409
  })
  model({
    # ---- PK ----
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * exp(e_wt_cl * (WT - 70))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Dose in mg, volume in L: mg/L x 1000 = ug/L, the unit of IC50.
    Cc <- 1000 * central / vc

    # ---- Average glucose (indirect response, inhibited input) ----
    rbase_glucose <- exp(lrbase_glucose + etalrbase_glucose)
    kout_glucose <- exp(lkout_glucose)
    kin_glucose <- kout_glucose * rbase_glucose
    plac_i <- plac + etaplac
    imax <- exp(limax)
    ic50 <- exp(lic50)
    inh_drug <- imax * Cc / (ic50 + Cc)
    glucose(0) <- rbase_glucose
    d/dt(glucose) <- kin_glucose * (1 - plac_i) * (1 - inh_drug) - kout_glucose * glucose

    # ---- Red blood cell transit model ----
    # KGL is per mg/dL, so AG (mmol/L) is converted with 18.016 mg/dL per
    # mmol/L (glucose 180.16 g/mol). Day-based rates are divided by 24 to
    # match the hourly model time.
    mmoll_to_mgdl <- 18.016
    ag_ref <- 149
    nc <- 12
    kin_hb <- exp(lkin_hb) / 24
    kg <- exp(lkg)
    ls <- exp(lls)
    lsp <- exp(llsp + etallsp)

    cag <- glucose * mmoll_to_mgdl
    kls <- nc / (ls * (cag / ag_ref)^gamma_ls * 24)
    glyc <- kg * cag / 24
    prec <- exp(-kg * cag * lsp)

    # Steady state at the baseline average glucose.
    cag0 <- rbase_glucose * mmoll_to_mgdl
    kls0 <- nc / (ls * (cag0 / ag_ref)^gamma_ls * 24)
    glyc0 <- kg * cag0 / 24
    prec0 <- exp(-kg * cag0 * lsp)
    rho0 <- kls0 / (kls0 + glyc0)
    ic_n01 <- prec0 * kin_hb / (kls0 + glyc0)
    ic_n02 <- ic_n01 * rho0
    ic_n03 <- ic_n02 * rho0
    ic_n04 <- ic_n03 * rho0
    ic_n05 <- ic_n04 * rho0
    ic_n06 <- ic_n05 * rho0
    ic_n07 <- ic_n06 * rho0
    ic_n08 <- ic_n07 * rho0
    ic_n09 <- ic_n08 * rho0
    ic_n10 <- ic_n09 * rho0
    ic_n11 <- ic_n10 * rho0
    ic_n12 <- ic_n11 * rho0
    ic_g01 <- (ic_n01 * glyc0 + kin_hb * (1 - prec0)) / kls0
    ic_g02 <- ic_g01 + ic_n02 * glyc0 / kls0
    ic_g03 <- ic_g02 + ic_n03 * glyc0 / kls0
    ic_g04 <- ic_g03 + ic_n04 * glyc0 / kls0
    ic_g05 <- ic_g04 + ic_n05 * glyc0 / kls0
    ic_g06 <- ic_g05 + ic_n06 * glyc0 / kls0
    ic_g07 <- ic_g06 + ic_n07 * glyc0 / kls0
    ic_g08 <- ic_g07 + ic_n08 * glyc0 / kls0
    ic_g09 <- ic_g08 + ic_n09 * glyc0 / kls0
    ic_g10 <- ic_g09 + ic_n10 * glyc0 / kls0
    ic_g11 <- ic_g10 + ic_n11 * glyc0 / kls0
    ic_g12 <- ic_g11 + ic_n12 * glyc0 / kls0

    hba_n01(0) <- ic_n01
    hba_n02(0) <- ic_n02
    hba_n03(0) <- ic_n03
    hba_n04(0) <- ic_n04
    hba_n05(0) <- ic_n05
    hba_n06(0) <- ic_n06
    hba_n07(0) <- ic_n07
    hba_n08(0) <- ic_n08
    hba_n09(0) <- ic_n09
    hba_n10(0) <- ic_n10
    hba_n11(0) <- ic_n11
    hba_n12(0) <- ic_n12
    hba_g01(0) <- ic_g01
    hba_g02(0) <- ic_g02
    hba_g03(0) <- ic_g03
    hba_g04(0) <- ic_g04
    hba_g05(0) <- ic_g05
    hba_g06(0) <- ic_g06
    hba_g07(0) <- ic_g07
    hba_g08(0) <- ic_g08
    hba_g09(0) <- ic_g09
    hba_g10(0) <- ic_g10
    hba_g11(0) <- ic_g11
    hba_g12(0) <- ic_g12

    d/dt(hba_n01) <- kin_hb * prec - (kls + glyc) * hba_n01
    d/dt(hba_n02) <- kls * hba_n01 - (kls + glyc) * hba_n02
    d/dt(hba_n03) <- kls * hba_n02 - (kls + glyc) * hba_n03
    d/dt(hba_n04) <- kls * hba_n03 - (kls + glyc) * hba_n04
    d/dt(hba_n05) <- kls * hba_n04 - (kls + glyc) * hba_n05
    d/dt(hba_n06) <- kls * hba_n05 - (kls + glyc) * hba_n06
    d/dt(hba_n07) <- kls * hba_n06 - (kls + glyc) * hba_n07
    d/dt(hba_n08) <- kls * hba_n07 - (kls + glyc) * hba_n08
    d/dt(hba_n09) <- kls * hba_n08 - (kls + glyc) * hba_n09
    d/dt(hba_n10) <- kls * hba_n09 - (kls + glyc) * hba_n10
    d/dt(hba_n11) <- kls * hba_n10 - (kls + glyc) * hba_n11
    d/dt(hba_n12) <- kls * hba_n11 - (kls + glyc) * hba_n12
    d/dt(hba_g01) <- glyc * hba_n01 + kin_hb * (1 - prec) - kls * hba_g01
    d/dt(hba_g02) <- glyc * hba_n02 + kls * (hba_g01 - hba_g02)
    d/dt(hba_g03) <- glyc * hba_n03 + kls * (hba_g02 - hba_g03)
    d/dt(hba_g04) <- glyc * hba_n04 + kls * (hba_g03 - hba_g04)
    d/dt(hba_g05) <- glyc * hba_n05 + kls * (hba_g04 - hba_g05)
    d/dt(hba_g06) <- glyc * hba_n06 + kls * (hba_g05 - hba_g06)
    d/dt(hba_g07) <- glyc * hba_n07 + kls * (hba_g06 - hba_g07)
    d/dt(hba_g08) <- glyc * hba_n08 + kls * (hba_g07 - hba_g08)
    d/dt(hba_g09) <- glyc * hba_n09 + kls * (hba_g08 - hba_g09)
    d/dt(hba_g10) <- glyc * hba_n10 + kls * (hba_g09 - hba_g10)
    d/dt(hba_g11) <- glyc * hba_n11 + kls * (hba_g10 - hba_g11)
    d/dt(hba_g12) <- glyc * hba_n12 + kls * (hba_g11 - hba_g12)

    rbc_non <- hba_n01 + hba_n02 + hba_n03 + hba_n04 + hba_n05 + hba_n06 +
      hba_n07 + hba_n08 + hba_n09 + hba_n10 + hba_n11 + hba_n12
    rbc_gly <- hba_g01 + hba_g02 + hba_g03 + hba_g04 + hba_g05 + hba_g06 +
      hba_g07 + hba_g08 + hba_g09 + hba_g10 + hba_g11 + hba_g12
    # Total hemoglobin (g/L); not an observed output of this model.
    hb <- rbc_non + rbc_gly
    Hba1c <- 100 * rbc_gly / hb

    Cc ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    Hba1c ~ prop(propSd_Hba1c)
  })
}
