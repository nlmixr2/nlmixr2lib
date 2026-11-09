Tan_2025_adomeglivant_ffh <- function() {
  description <- paste(
    "Population PK/PD model for adomeglivant (LY2409021, oral glucagon",
    "receptor antagonist) in adults with type 2 diabetes (Tan 2025, FPG-FSI-",
    "HbA1c (FFH) model). One-compartment PK with first-order absorption and",
    "body-weight effects on CL/F and V/F. Fasting plasma glucose (FPG) and",
    "fasting serum insulin (FSI) follow the steady-state solution of the",
    "HOMA glucose-insulin feedback: baseline beta-cell function and insulin",
    "sensitivity are derived from each subject's baseline FPG and FSI and",
    "progress linearly on the logit scale over time (with an additional",
    "beta-cell term in drug-treated subjects); adomeglivant (inhibitory",
    "Imax) and a fractional placebo effect reduce the glucose input. HbA1c",
    "is an indirect response driven linearly by FPG, with a baseline-FPG",
    "effect on its loss rate. One of four alternative semi-mechanistic",
    "HbA1c models fit sequentially to the same 4-week Phase 1b data",
    "(siblings: Tan_2025_adomeglivant_adopt, _fhh, _igrh)."
  )
  reference <- paste(
    "Tan AYF, Schneck K, Garhyan P, Chan ECY, Tham LS.",
    "Evaluation of Four Semi-Mechanistic Models for Predicting Glycemic",
    "Control With a Glucagon Receptor Antagonist in People With Type 2",
    "Diabetes. CPT Pharmacometrics Syst Pharmacol. 2025;14:1381-1390.",
    "doi:10.1002/psp4.70058.",
    "FFH structure: de Winter W et al. J Pharmacokinet Pharmacodyn.",
    "2006;33:313-343. doi:10.1007/s10928-006-9008-2; steady-state FPG/FSI",
    "solution and HOMA constants as deposited in Kunina H et al. CPT",
    "Pharmacometrics Syst Pharmacol. 2022;11:1443-1457 (supplements S3, S4).",
    "doi:10.1002/psp4.12854."
  )
  vignette <- "Tan_2025_adomeglivant_hba1c"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "adomeglivant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "adomeglivant", units = "mg", specimen = "plasma", verified = TRUE),
    Hba1c = list(analyte = "glycated hemoglobin (HbA1c)", units = "%", specimen = "blood cell", verified = TRUE)
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
    FPG = list(
      description = "Baseline fasting glucose (fasting serum glucose used in place of fasting plasma glucose)",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-subject observed baseline. Sets the HOMA baseline beta-cell",
        "function FSI / (5 * (FPG - 3.5)) and insulin sensitivity",
        "22.5 / (FPG * FSI), so FPG at time 0 equals this value; also a",
        "linear-deviation covariate on KOUTH centred on the population",
        "median 8.2 mmol/L (Table S4 footnote e). Must exceed 3.5 mmol/L,",
        "and together with INS_BL must give HOMA beta-cell function and",
        "insulin sensitivity below 1 (100%) for the logit progression."
      ),
      source_name = "Baseline fasting glucose (FSG)"
    ),
    INS_BL = list(
      description = "Baseline fasting serum insulin",
      units = "uU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-subject observed baseline (mU/L = uU/mL, the HOMA unit). Sets",
        "the HOMA baseline beta-cell function and insulin sensitivity",
        "together with FPG, so FSI at time 0 equals this value. Requires",
        "INS_BL < 5 * (FPG - 3.5) and INS_BL * FPG > 22.5."
      ),
      source_name = "Baseline FSI"
    ),
    PLACEBO = list(
      description = "Placebo-arm indicator: 1 = placebo, 0 = adomeglivant-treated",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adomeglivant-treated)",
      notes = paste(
        "Time-fixed. Adds the drug-treated-group beta-cell disease-progression",
        "coefficient (Table S4, -0.2 logit/year additional) when PLACEBO = 0.",
        "The concentration-driven Imax effect is modelled separately."
      ),
      source_name = "Treatment group"
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
      "928 PK, 460 fasting serum glucose, 318 fasting serum insulin, 78 HbA1c",
      "records). Fasting serum glucose was used in place of fasting plasma",
      "glucose. External evaluation: 24-week Phase 2b study NCT01241448 (254",
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

    # ---- Disease progression, logit scale per year (Table S4) ----
    # Table S4 footnote c: negative values indicate improvement, positive deterioration.
    dp_bcf <- fixed(-0.2); label("Beta-cell function disease-progression rate, placebo group (logit/year)") # Table S4 beta-cell function of placebo group = -0.2 Fixed
    dp_bcf_drug <- fixed(-0.2); label("Additional beta-cell function disease-progression rate in drug-treated subjects (logit/year)") # Table S4 beta-cell function of drug-treated group (additional) = -0.2 Fixed
    dp_is <- fixed(0.1); label("Insulin-sensitivity disease-progression rate, all subjects (logit/year)") # Table S4 insulin sensitivity = 0.1 Fixed

    # ---- Glucose input effects (Table S4) ----
    plac <- fixed(0.0556); label("Fractional placebo inhibition of the glucose input PLAC (fraction)") # Table S4 PLAC = 0.0556 Fixed (footnote d: from placebo data alone)
    limax <- log(0.698); label("Maximal fractional inhibition of the glucose input IMAX (fraction)") # Table S4 IMAX = 0.698
    lic50 <- log(2210); label("Adomeglivant concentration at half-maximal inhibition IC50 (ug/L)") # Table S4 IC50 = 2210

    # ---- HbA1c turnover (Table S4) ----
    lkin_hba1c <- log(0.000906); label("HbA1c input rate constant KINH per unit FPG (%/h per mmol/L)") # Table S4 KINH = 0.000906
    lkout_hba1c <- log(0.000955); label("First-order HbA1c loss rate constant KOUTH at baseline FPG 8.2 mmol/L (1/h)") # Table S4 KOUTH = 0.000955
    e_fpg_kout_hba1c <- 0.0563; label("Linear baseline-FPG coefficient on KOUTH, (1 + e * (FPG - 8.2)) (L/mmol)") # Table S4 footnote e = 0.0563

    # ---- IIV ----
    # Table S2 footnote b: CV% = 100 * sqrt(exp(OMEGA) - 1), so OMEGA = log(1 + CV^2).
    # Table S2 footnote e: correlation = COV / sqrt(OMEGA(CL) * OMEGA(V)).
    etalka ~ log(1 + 0.262^2) # Table S2 KA IIV 26.2 %CV
    etalcl + etalvc ~ c(
      log(1 + 0.180^2),
      0.778 * sqrt(log(1 + 0.180^2) * log(1 + 0.124^2)),
      log(1 + 0.124^2)
    ) # Table S2 CL/F IIV 18.0 %CV, V/F IIV 12.4 %CV, correlation 0.778
    # Table S4 footnote b: additive IIV on PLAC, CV% = 100 * sqrt(OMEGA) / THETA,
    # so OMEGA = (CV * THETA)^2. Footnote d: fixed to 150% by sensitivity analysis.
    etaplac ~ fixed((1.50 * 0.0556)^2) # Table S4 PLAC IIV 150 %CV (additive), set by sensitivity analysis
    etalkout_hba1c ~ log(1 + 0.0507^2) # Table S4 KOUTH IIV 5.07 %CV (exponential)

    # ---- Residual error ----
    propSd <- 0.290; label("Proportional residual error, adomeglivant concentration (fraction)") # Table S2 proportional residual error 0.290
    propSd_fpg <- 0.102; label("Proportional residual error, fasting glucose (fraction)") # Table S4 glucose residual error 0.102
    propSd_insulin <- 0.744; label("Proportional residual error, fasting insulin (fraction)") # Table S4 insulin residual error 0.744
    propSd_Hba1c <- 0.0315; label("Proportional residual error, HbA1c (fraction)") # Table S4 HbA1c residual error 0.0315
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

    # ---- HOMA steady-state glucose-insulin system ----
    # HOMA constants: insulin secretion FSI = ki * BCF * (FPG - 3.5) and
    # glucose balance FPG = kg_homa / (IS * FSI).
    ki <- 5
    kg_homa <- 22.5
    bcf0 <- INS_BL / (ki * (FPG - 3.5))
    ins_sens0 <- kg_homa / (INS_BL * FPG)
    # Progression on the logit scale with time in years: BCF(t) =
    # 1 / (1 + exp(log((1 - BCF0) / BCF0) + rate * t)), so a positive rate lowers it.
    t_yr <- t / (24 * 365.25)
    bcf <- expit(logit(bcf0) - (dp_bcf + dp_bcf_drug * (1 - PLACEBO)) * t_yr)
    ins_sens <- expit(logit(ins_sens0) - dp_is * t_yr)

    # The placebo effect starts with treatment at time 0; gating on t > 0
    # keeps the algebraic FPG and FSI at their observed baselines at t = 0.
    plac_i <- (plac + etaplac) * (t > 0)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    inh_drug <- imax * Cc / (ic50 + Cc)
    eff_glucose <- (1 - plac_i) * (1 - inh_drug)

    # Positive root of FSI^2 + 3.5 * ki * BCF * FSI - ki * BCF * kg_homa * eff / IS = 0.
    pp <- 3.5 * ki * bcf
    qq <- ki * bcf * kg_homa * eff_glucose / ins_sens
    insulin <- (-pp + sqrt(pp^2 + 4 * qq)) / 2
    fpg <- kg_homa * eff_glucose / (ins_sens * insulin)

    # ---- HbA1c (indirect response driven by FPG) ----
    kin_hba1c <- exp(lkin_hba1c)
    kout_hba1c <- exp(lkout_hba1c + etalkout_hba1c) * (1 + e_fpg_kout_hba1c * (FPG - 8.2))
    Hba1c(0) <- kin_hba1c * FPG / kout_hba1c
    d/dt(Hba1c) <- kin_hba1c * fpg - kout_hba1c * Hba1c

    Cc ~ prop(propSd)
    fpg ~ prop(propSd_fpg)
    insulin ~ prop(propSd_insulin)
    Hba1c ~ prop(propSd_Hba1c)
  })
}
