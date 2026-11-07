Tan_2025_adomeglivant_fhh <- function() {
  description <- paste(
    "Population PK/PD model for adomeglivant (LY2409021, oral glucagon",
    "receptor antagonist) in adults with type 2 diabetes (Tan 2025, FPG-Hb-",
    "HbA1c (FHH) model). One-compartment PK with first-order absorption and",
    "body-weight effects on CL/F and V/F, driving an inhibitory Imax effect",
    "on the zero-order input of fasting plasma glucose (FPG, indirect",
    "response) together with a fractional placebo effect on the same input.",
    "Red blood cells age through four non-glycated and four glycated transit",
    "compartments (mean residence time 120 days, sex-specific hemoglobin",
    "production); glycation is driven nonlinearly by FPG (FPG^gamma). HbA1c",
    "is the glycated percentage and hemoglobin the total of the eight",
    "compartments. One of four alternative semi-mechanistic HbA1c models fit",
    "sequentially to the same 4-week Phase 1b data (siblings:",
    "Tan_2025_adomeglivant_adopt, _ffh, _igrh)."
  )
  reference <- paste(
    "Tan AYF, Schneck K, Garhyan P, Chan ECY, Tham LS.",
    "Evaluation of Four Semi-Mechanistic Models for Predicting Glycemic",
    "Control With a Glucagon Receptor Antagonist in People With Type 2",
    "Diabetes. CPT Pharmacometrics Syst Pharmacol. 2025;14:1381-1390.",
    "doi:10.1002/psp4.70058.",
    "FHH structure: Hamren B et al. Clin Pharmacol Ther. 2008;84:228-235.",
    "doi:10.1038/clpt.2008.2; equations as deposited in Kunina H et al.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11:1443-1457 (supplement S5).",
    "doi:10.1002/psp4.12854."
  )
  vignette <- "Tan_2025_adomeglivant_hba1c"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  paper_specific_compartments <- c(
    "hba_n01",
    "hba_n02",
    "hba_n03",
    "hba_n04",
    "hba_g01",
    "hba_g02",
    "hba_g03",
    "hba_g04"
  )

  rbcNon <- list(analyte = "non-glycated hemoglobin", units = "g/L", specimen = "blood cell", verified = TRUE)
  rbcGly <- list(analyte = "glycated hemoglobin (HbA1c)", units = "g/L", specimen = "blood cell", verified = TRUE)
  compartmentData <- list(
    depot = list(analyte = "adomeglivant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "adomeglivant", units = "mg", specimen = "plasma", verified = TRUE),
    glucose = list(analyte = "glucose (fasting)", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hba_n01 = rbcNon,
    hba_n02 = rbcNon,
    hba_n03 = rbcNon,
    hba_n04 = rbcNon,
    hba_g01 = rbcGly,
    hba_g02 = rbcGly,
    hba_g03 = rbcGly,
    hba_g04 = rbcGly
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value; reference 70 kg (Table S2 footnotes c/d). Exponential effect on CL/F, power effect on V/F.",
      source_name = "Baseline weight"
    ),
    SEXF = list(
      description = "Sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the sex-specific zero-order red-blood-cell (hemoglobin) release constant KINHB (Table S5: male 1.21, female 1.12 g/L/day).",
      source_name = "Sex"
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
      "928 PK, 460 fasting serum glucose, 78 HbA1c, 316 hemoglobin records).",
      "Fasting serum glucose was used in place of fasting plasma glucose.",
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

    # ---- Fasting glucose turnover (Table S5) ----
    lrbase_glucose <- log(8.72); label("Baseline fasting glucose FG0 (mmol/L)") # Table S5 FG0 = 8.72
    lkout_glucose <- log(5.10); label("First-order loss rate constant of fasting glucose KOUTG (1/h)") # Table S5 KOUTG = 5.10
    plac <- fixed(0.055); label("Fractional placebo inhibition of the glucose input PLAC (fraction)") # Table S5 PLAC = 0.055 Fixed (footnote c: from placebo data alone)
    limax <- log(0.320); label("Maximal fractional inhibition of the glucose input IMAX (fraction)") # Table S5 IMAX = 0.320
    lic50 <- log(1900); label("Adomeglivant concentration at half-maximal inhibition IC50 (ug/L)") # Table S5 IC50 = 1900

    # ---- Red blood cell / HbA1c model (Table S5) ----
    lkin_hb_male <- log(1.21); label("Zero-order release constant of RBC hemoglobin KINHB, males (g/L/day)") # Table S5 KINHB (male) = 1.21
    lkin_hb_female <- log(1.12); label("Zero-order release constant of RBC hemoglobin KINHB, females (g/L/day)") # Table S5 KINHB (female) = 1.12
    lkgl <- log(0.000224); label("Glycosylation rate constant KGL (1/day per (mmol/L)^gamma)") # Table S5 KGL = 0.000224
    lmtt_rbc <- fixed(log(120)); label("Mean residence time of red blood cells MRTE (day)") # Table S5 MRTE = 120 Fixed (footnote d: physiological value)
    lgamma <- log(0.741); label("Power exponent of FPG in the glycation rate GAM (unitless)") # Table S5 GAM = 0.741

    # ---- IIV ----
    # Table S2 footnote b: CV% = 100 * sqrt(exp(OMEGA) - 1), so OMEGA = log(1 + CV^2).
    # Table S2 footnote e: correlation = COV / sqrt(OMEGA(CL) * OMEGA(V)).
    etalka ~ log(1 + 0.262^2) # Table S2 KA IIV 26.2 %CV
    etalcl + etalvc ~ c(
      log(1 + 0.180^2),
      0.778 * sqrt(log(1 + 0.180^2) * log(1 + 0.124^2)),
      log(1 + 0.124^2)
    ) # Table S2 CL/F IIV 18.0 %CV, V/F IIV 12.4 %CV, correlation 0.778
    etalrbase_glucose ~ log(1 + 0.183^2) # Table S5 FG0 IIV 18.3 %CV (exponential)
    # Table S5 footnote b: additive IIV on PLAC, CV% = 100 * sqrt(OMEGA) / THETA,
    # so OMEGA = (CV * THETA)^2.
    etaplac ~ (1.41 * 0.055)^2 # Table S5 PLAC IIV 141 %CV (additive)
    etalmtt_rbc ~ log(1 + 0.0553^2) # Table S5 MRTE IIV 5.53 %CV (exponential)
    etalgamma ~ log(1 + 0.0417^2) # Table S5 GAM IIV 4.17 %CV (exponential)

    # ---- Residual error ----
    propSd <- 0.290; label("Proportional residual error, adomeglivant concentration (fraction)") # Table S2 proportional residual error 0.290
    propSd_glucose <- 0.0964; label("Proportional residual error, fasting glucose (fraction)") # Table S5 glucose residual error 0.0964
    propSd_Hba1c <- 0.0420; label("Proportional residual error, HbA1c (fraction)") # Table S5 HbA1c residual error 0.0420
    propSd_hb <- 0.0407; label("Proportional residual error, hemoglobin (fraction)") # Table S5 Hb residual error 0.0407
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

    # ---- Fasting glucose (indirect response, inhibited input) ----
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
    # Day-based parameters are converted to the hourly model time scale.
    kin_hb <- (exp(lkin_hb_male) * (1 - SEXF) + exp(lkin_hb_female) * SEXF) / 24
    mtt_rbc <- exp(lmtt_rbc + etalmtt_rbc) * 24
    ktr <- 4 / mtt_rbc
    kgl <- exp(lkgl) / 24
    gam <- exp(lgamma + etalgamma)
    glyc <- kgl * glucose^gam

    # Steady state at the baseline fasting glucose.
    glyc0 <- kgl * rbase_glucose^gam
    rho0 <- ktr / (ktr + glyc0)
    hba_n01(0) <- kin_hb / (ktr + glyc0)
    hba_n02(0) <- kin_hb / (ktr + glyc0) * rho0
    hba_n03(0) <- kin_hb / (ktr + glyc0) * rho0^2
    hba_n04(0) <- kin_hb / (ktr + glyc0) * rho0^3
    hba_g01(0) <- kin_hb / (ktr + glyc0) * glyc0 / ktr
    hba_g02(0) <- kin_hb / (ktr + glyc0) * glyc0 / ktr * (1 + rho0)
    hba_g03(0) <- kin_hb / (ktr + glyc0) * glyc0 / ktr * (1 + rho0 + rho0^2)
    hba_g04(0) <- kin_hb / (ktr + glyc0) * glyc0 / ktr * (1 + rho0 + rho0^2 + rho0^3)

    d/dt(hba_n01) <- kin_hb - (ktr + glyc) * hba_n01
    d/dt(hba_n02) <- ktr * hba_n01 - (ktr + glyc) * hba_n02
    d/dt(hba_n03) <- ktr * hba_n02 - (ktr + glyc) * hba_n03
    d/dt(hba_n04) <- ktr * hba_n03 - (ktr + glyc) * hba_n04
    d/dt(hba_g01) <- glyc * hba_n01 - ktr * hba_g01
    d/dt(hba_g02) <- glyc * hba_n02 + ktr * (hba_g01 - hba_g02)
    d/dt(hba_g03) <- glyc * hba_n03 + ktr * (hba_g02 - hba_g03)
    d/dt(hba_g04) <- glyc * hba_n04 + ktr * (hba_g03 - hba_g04)

    rbc_non <- hba_n01 + hba_n02 + hba_n03 + hba_n04
    rbc_gly <- hba_g01 + hba_g02 + hba_g03 + hba_g04
    hb <- rbc_non + rbc_gly
    Hba1c <- 100 * rbc_gly / hb

    Cc ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    Hba1c ~ prop(propSd_Hba1c)
    hb ~ prop(propSd_hb)
  })
}
