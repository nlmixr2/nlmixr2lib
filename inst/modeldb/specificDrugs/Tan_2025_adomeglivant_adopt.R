Tan_2025_adomeglivant_adopt <- function() {
  description <- paste(
    "Population PK/PD model for adomeglivant (LY2409021, oral glucagon",
    "receptor antagonist) in adults with type 2 diabetes (Tan 2025, ADOPT",
    "HbA1c model). One-compartment PK with first-order absorption and",
    "body-weight effects on CL/F and V/F, driving an inhibitory Imax effect",
    "on the zero-order input of average glucose (AG, indirect response)",
    "together with a fractional placebo effect on the same input; HbA1c is",
    "a second indirect response whose input is driven by (AG + beta). One of",
    "four alternative semi-mechanistic HbA1c models fit sequentially to the",
    "same 4-week Phase 1b data (siblings: Tan_2025_adomeglivant_ffh,",
    "_fhh, _igrh)."
  )
  reference <- paste(
    "Tan AYF, Schneck K, Garhyan P, Chan ECY, Tham LS.",
    "Evaluation of Four Semi-Mechanistic Models for Predicting Glycemic",
    "Control With a Glucagon Receptor Antagonist in People With Type 2",
    "Diabetes. CPT Pharmacometrics Syst Pharmacol. 2025;14:1381-1390.",
    "doi:10.1002/psp4.70058.",
    "ADOPT structure: Moller JB et al. CPT Pharmacometrics Syst Pharmacol.",
    "2013;2:e82. doi:10.1038/psp.2013.58; equations as deposited in",
    "Kunina H et al. CPT Pharmacometrics Syst Pharmacol. 2022;11:1443-1457",
    "(supplement S2). doi:10.1002/psp4.12854."
  )
  vignette <- "Tan_2025_adomeglivant_hba1c"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "adomeglivant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "adomeglivant", units = "mg", specimen = "plasma", verified = TRUE),
    glucose = list(analyte = "glucose (24-h average)", units = "mmol/L", specimen = "plasma", verified = TRUE),
    Hba1c = list(analyte = "glycated hemoglobin (HbA1c)", units = "%", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline value. Reference 70 kg (Table S2 footnotes c/d,",
        "Table S3 footnote c). Exponential effect on CL/F, power effect on",
        "V/F, and a linear-deviation effect on the baseline average glucose."
      ),
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
      "928 PK, 460 fasting serum glucose, 78 HbA1c, 138 average-glucose",
      "records). External evaluation: 24-week Phase 2b study NCT01241448",
      "(254 patients; placebo, 2.5, 10 or 20 mg once daily). Baseline HbA1c",
      "8.00% (SD 0.986), fasting glucose 8.64 mmol/L (SD 1.98), average",
      "glucose 10.2 mmol/L (SD 2.33) in Phase 1b."
    )
  )

  ini({
    # ---- PK (Table S2) ----
    lka <- log(0.385); label("First-order absorption rate constant KA (1/h)") # Table S2 KA = 0.385
    lcl <- log(0.283); label("Apparent clearance CL/F at 70 kg (L/h)") # Table S2 CL/F = 0.283
    lvc <- log(24.3); label("Apparent volume V/F at 70 kg (L)") # Table S2 V/F = 24.3
    e_wt_cl <- 0.0116; label("Exponential body-weight coefficient on CL/F, exp(e * (WT - 70)) (1/kg)") # Table S2 footnote c = 0.0116
    e_wt_vc <- 0.956; label("Power exponent of (WT/70) on V/F (unitless)") # Table S2 footnote d = 0.956

    # ---- Average glucose turnover (Table S3) ----
    lrbase_glucose <- log(9.62); label("Baseline average glucose AG0 at 70 kg (mmol/L)") # Table S3 AG0 = 9.62
    e_wt_rbase_glucose <- 0.00507; label("Linear body-weight coefficient on AG0, (1 + e * (WT - 70)) (1/kg)") # Table S3 footnote c = 0.00507
    lkout_glucose <- log(0.429); label("First-order loss rate constant of average glucose KOUTG (1/h)") # Table S3 KOUTG = 0.429
    plac <- 0.122; label("Fractional placebo inhibition of the glucose input PLAC (fraction)") # Table S3 PLAC = 0.122
    limax <- log(0.167); label("Maximal fractional inhibition of the glucose input IMAX (fraction)") # Table S3 IMAX = 0.167
    lic50 <- log(1620); label("Adomeglivant concentration at half-maximal inhibition IC50 (ug/L)") # Table S3 IC50 = 1620

    # ---- HbA1c turnover (Table S3) ----
    lkin_hba1c <- log(0.000706); label("HbA1c input rate constant KINH per unit (AG + beta) (%/h per mmol/L)") # Table S3 KINH = 0.000706
    lkout_hba1c <- log(0.00118); label("First-order HbA1c loss rate constant KOUTH (1/h)") # Table S3 KOUTH = 0.00118
    beta_ag <- fixed(3); label("Intercept beta of the AG-HbA1c steady-state relationship (mmol/L)") # Table S3 beta = 3 Fixed (footnote d: sensitivity analysis)

    # ---- IIV ----
    # Table S2 footnote b: CV% = 100 * sqrt(exp(OMEGA) - 1), so OMEGA = log(1 + CV^2).
    # Table S2 footnote e: correlation = COV / sqrt(OMEGA(CL) * OMEGA(V)).
    etalka ~ log(1 + 0.262^2) # Table S2 KA IIV 26.2 %CV
    etalcl + etalvc ~ c(
      log(1 + 0.180^2),
      0.778 * sqrt(log(1 + 0.180^2) * log(1 + 0.124^2)),
      log(1 + 0.124^2)
    ) # Table S2 CL/F IIV 18.0 %CV, V/F IIV 12.4 %CV, correlation 0.778
    etalrbase_glucose ~ log(1 + 0.132^2) # Table S3 AG0 IIV 13.2 %CV (exponential)
    # Table S3 footnote b: additive IIV on PLAC, CV% = 100 * sqrt(OMEGA) / THETA,
    # so OMEGA = (CV * THETA)^2.
    etaplac ~ (0.597 * 0.122)^2 # Table S3 PLAC IIV 59.7 %CV (additive)

    # ---- Residual error ----
    propSd <- 0.290; label("Proportional residual error, adomeglivant concentration (fraction)") # Table S2 proportional residual error 0.290
    propSd_glucose <- 0.111; label("Proportional residual error, average glucose (fraction)") # Table S3 glucose residual error 0.111
    propSd_Hba1c <- 0.0279; label("Proportional residual error, HbA1c (fraction)") # Table S3 HbA1c residual error 0.0279
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
    rbase_glucose <- exp(lrbase_glucose + etalrbase_glucose) * (1 + e_wt_rbase_glucose * (WT - 70))
    kout_glucose <- exp(lkout_glucose)
    kin_glucose <- kout_glucose * rbase_glucose
    plac_i <- plac + etaplac
    imax <- exp(limax)
    ic50 <- exp(lic50)
    # Placebo and drug act multiplicatively on the glucose input from the
    # start of treatment (time 0); see the vignette for the form.
    inh_drug <- imax * Cc / (ic50 + Cc)
    glucose(0) <- rbase_glucose
    d/dt(glucose) <- kin_glucose * (1 - plac_i) * (1 - inh_drug) - kout_glucose * glucose

    # ---- HbA1c (indirect response driven by AG + beta) ----
    kin_hba1c <- exp(lkin_hba1c)
    kout_hba1c <- exp(lkout_hba1c)
    Hba1c(0) <- kin_hba1c * (rbase_glucose + beta_ag) / kout_hba1c
    d/dt(Hba1c) <- kin_hba1c * (glucose + beta_ag) - kout_hba1c * Hba1c

    Cc ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    Hba1c ~ prop(propSd_Hba1c)
  })
}
