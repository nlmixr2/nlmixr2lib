Huppe_2023_fosfomycin <- function() {
  description <- "Two-compartment population PK model for multiple-dose intravenous fosfomycin in critically ill adults with renal insufficiency during continuous venovenous hemodialysis (CVVHD). Total clearance is the sum of a renal arm driven linearly by measured urinary creatinine clearance (gated off in anuric patients) and a CVVHD arm given by the Michaels hemodialyzer equation as a function of blood flow rate, dialysate flow rate, and a mass transfer-area coefficient. Central volume increases linearly with time since first dose (Huppe 2023)."
  reference <- "Huppe T, Gotz KM, Meiser A, de Faria Fernandes A, Maurer F, Groesdonk HV, Volk T, Lehr T, Kreuer S. Population pharmacokinetic modeling of multiple-dose intravenous fosfomycin in critically ill patients during continuous venovenous hemodialysis. Sci Rep. 2023;13:18132. doi:10.1038/s41598-023-45084-5"
  vignette <- "Huppe_2023_fosfomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Urinary creatinine clearance, measured from a 12-hour urine collection",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "MEASURED urinary creatinine clearance, NOT a creatinine-based estimate and NOT",
        "BSA-normalized. Huppe 2023 Methods Eq. 2 computes it as",
        "(urine creatinine [mg/dL] * urine volume [mL]) / (plasma creatinine [mg/dL] * collection time [min]),",
        "recalculated daily from urine collected over 12 h. Raw mL/min, so the canonical",
        "mL/min/1.73 m^2 normalisation is deliberately NOT applied -- the effect coefficient",
        "0.0723 per (mL/min) was calibrated against raw values. Cohort mean 20.7 +/- 44.9 mL/min",
        "(Table 1); the paper stresses (Discussion) that Cockcroft-Gault estimates would",
        "overestimate residual renal function in this overweight cohort and are NOT",
        "interchangeable with this measured value. Time-varying (recalculated daily).",
        "Enters the renal clearance arm as a linear, uncentered effect. Anuric patients have",
        "CRCL = 0 by construction (no urine collected)."
      ),
      source_name        = "CLCR"
    ),
    URINE_VOL_24H = list(
      description        = "Residual diuresis: total urine volume collected over 24 hours",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used ONLY as the preserved-diuresis gate on the renal clearance arm. Huppe 2023",
        "Methods, Population pharmacokinetic modeling: 'In patients without preserved diuresis,",
        "intrinsic fosfomycin elimination was fixed to zero. Preserved diuresis was defined as a",
        "residual diuresis exceeding 100 mL in 24 h.' The model therefore multiplies the renal arm",
        "by (URINE_VOL_24H > 100). Cohort mean 24-h urine output 480 +/- 740 mL (Table 1);",
        "6 of 15 patients were anuric. The 100 mL/24h cutoff coincides with the anuria cutoff",
        "already documented in this canonical's register entry."
      ),
      source_name        = "residual diuresis"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "CVVHD-active indicator (1 while continuous venovenous hemodialysis is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CVVHD not running)",
      notes              = paste(
        "Time-varying WITHIN subject: every patient contributed one PK profile with CVVHD",
        "running and one after CVVHD was interrupted for tubing-system replacement",
        "(Huppe 2023 Methods, Protocol and sample collection). Gates the Michaels-equation",
        "dialysis clearance arm -- 'Dialysis clearance was fixed to zero in periods without",
        "CVVHD treatment' (Results, Population pharmacokinetic modeling). The modality is",
        "CONTINUOUS venovenous hemodialysis, hence the CRRT rather than the intermittent-",
        "hemodialysis member of the RRT_<modality>_<kind> canonical family."
      ),
      source_name        = "CVVHD"
    ),
    BFR = list(
      description        = "Blood flow rate through the CVVHD extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject; an individual input variable to the Michaels equation",
        "(Huppe 2023 Results: 'Individual BFR and DFR values were input variables for Eq. 1').",
        "Protocol start value 100 mL/min, subsequently adjusted per patient (Methods, CVVHD);",
        "cohort mean 110 +/- 25 mL/min (Table 1), per-patient values 100-200 mL/min.",
        "Meaningful only while RRT_CRRT_ACTIVE = 1; the Michaels term is gated off otherwise."
      ),
      source_name        = "BFR"
    ),
    DFR = list(
      description        = "Dialysate flow rate through the CVVHD extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject; the second individual input variable to the Michaels",
        "equation. Stored in the canonical mL/min and converted inside model() to the L/h in",
        "which Huppe 2023 recorded it (protocol start 2 L/h = 33.33 mL/min; cohort mean",
        "2.3 +/- 0.5 L/h; per-patient values 2-4 L/h, Table 1). The conversion is load-bearing:",
        "the published K0A estimate was calibrated against DFR expressed in L/h -- see the",
        "lkoa source-trace comment in ini() and the vignette Errata.",
        "Meaningful only while RRT_CRRT_ACTIVE = 1."
      ),
      source_name        = "DFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a covariate on central volume and not retained (Huppe 2023 Results: 'Variations in weight, creatinine clearance and serum levels of creatinine, albumin, total protein, urea, potassium, and sodium failed to adequately account for changes in central compartment volume over time'). Cohort mean 88.5 +/- 20.5 kg. Fosfomycin dosing in this study was not weight-adjusted (Limitations)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Huppe 2023 Methods: 'Time after first dose, age, weight, creatinine clearance and serum levels of creatinine, albumin, total protein, urea, potassium, and sodium were tested as covariates') and not retained. Cohort mean 60 +/- 8 years."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 1.7 +/- 0.9 mg/dL (Table 1). Distinct from the retained CRCL, which is a measured urinary clearance rather than a serum concentration."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 24 +/- 3.5 g/L (Table 1). The Limitations note that fosfomycin albumin binding is low, so total-concentration modelling was considered acceptable."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 54 +/- 7.6 g/L (Table 1)."
    ),
    UREA = list(
      description = "Serum urea",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 44 +/- 27 mg/dL (Table 1)."
    ),
    POTASSIUM = list(
      description = "Serum potassium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 4.3 +/- 0.4 mmol/L (Table 1)."
    ),
    SODIUM = list(
      description = "Serum sodium",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened and not retained. Cohort mean 144 +/- 3.8 mmol/L (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    n_observations = 300L,
    age_range      = "49-80 years (mean 60 +/- 8)",
    weight_range   = "55-122 kg (mean 88.5 +/- 20.5)",
    sex_female_pct = 13,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill ICU adults with infection caused by fosfomycin-susceptible bacteria AND",
      "renal insufficiency requiring continuous venovenous hemodialysis. Indications for CVVHD",
      "were acute renal failure, acute-on-chronic renal failure, chronic renal failure with",
      "permanent dialysis, and hepatorenal syndrome (Table 1). 6 of 15 patients were anuric.",
      "Body mass index < 35 kg/m^2 and age >= 18 years were inclusion criteria; pregnancy and",
      "breastfeeding were exclusions. SAPS II 33-79 points. 4 of 15 received vasoactive drugs and",
      "7 of 15 were mechanically ventilated."
    ),
    renal_function = paste(
      "Serum creatinine 1.7 +/- 0.9 mg/dL; measured urinary creatinine clearance from 12-h urine",
      "20.7 +/- 44.9 mL/min; 12-h urine output 250 +/- 430 mL; 24-h urine output 480 +/- 740 mL."
    ),
    rrt_settings   = paste(
      "CVVHD with a multiFiltrate Ci-Ca system and Ultraflux AV 1000S polysulfone membrane",
      "hemofilters (Fresenius Medical Care), calcium- and phosphate-free dialysate (Ci-Ca Dialysate",
      "K4, potassium 4 mmol/L), regional citrate anticoagulation with calcium re-substitution.",
      "Dialysate flow rate 2.3 +/- 0.5 L/h, blood flow rate 110 +/- 25 mL/min, ultrafiltration",
      "88 +/- 83 mL/h. Protocol start values were BFR 100 mL/min and DFR 2 L/h, then individually",
      "adjusted."
    ),
    co_medication  = "Concomitant antimicrobials in 11 of 15 patients: meropenem, piperacillin/tazobactam, vancomycin, clindamycin, tigecycline, linezolid, cotrimoxazole (Table 1).",
    dose_range     = "5 g fosfomycin (InfectoFos) intravenously over 120 min every 8 h (multiple dosing).",
    regions        = "Single-centre observational study, Saarland University Medical Center, Homburg (Saar), Germany.",
    notes          = paste(
      "Baseline demographics from Huppe 2023 Table 1. 13 of 15 patients were male (87%).",
      "Each patient contributed two PK series, one during CVVHD and one after CVVHD was",
      "interrupted, separated on average by 41.5 +/- 40.3 h; for 60% of patients the CVVHD series",
      "came first. Samples were drawn 5 min before and 15, 30, 60, 90, 180, 240, 300 and 360 min",
      "after the start of each infusion. Fosfomycin was quantified by HPLC-MS in negative SIM mode",
      "(m/z -137) over a 25-700 ug/mL calibration range."
    )
  )

  ini({
    # Structural parameters (Huppe 2023 Table 2, "Structural model parameters").
    # Renal arm intercept: the clearance of a patient with preserved diuresis and
    # CRCL = 0. Reported as CL_Renal = 0.263 L/h (RSE 16%); the Results text
    # rounds it to 0.26 L/h.
    lcl_renal <- log(0.263); label("Renal fosfomycin clearance intercept in patients with preserved diuresis (CL_Renal, L/h)")  # Huppe 2023 Table 2
    lvc       <- log(18.20); label("Central compartment volume at time since first dose = 0 (V_C, L)")                          # Huppe 2023 Table 2
    lvp       <- log(20.80); label("Peripheral compartment volume (V_P, L)")                                                    # Huppe 2023 Table 2

    # Intercompartmental clearance was FIXED (Table 2 footnote c; RSE reported
    # as "-").
    lq        <- fixed(log(5.08)); label("Intercompartmental clearance (Q, L/h)")                                              # Huppe 2023 Table 2 (fixed, footnote c)

    # Hemodialyzer mass transfer-area coefficient for the Michaels equation
    # (Eq. 1), estimated with RSE 9%. Table 2 labels the units "mL/min", but that
    # label is not self-consistent with Eq. 1 as printed: substituting 0.0288
    # mL/min with BFR and DFR in matching mL/min yields CL_CVVHD = 0.0288 mL/min
    # = 0.0017 L/h, i.e. a negligible dialysis clearance, which contradicts the
    # paper's own Fig. 3b (dialysis supplies 100% of total clearance in anuric
    # patients) and its NCA result that CVVHD lowered AUC by 28.7%. Eq. 1
    # reproduces Fig. 3b exactly when BFR is supplied in mL/min and DFR in L/h
    # -- the raw units of the study's own data columns (Table 1 reports CVVHD
    # settings as "Bloodflow [ml/min] / Dialysate [l/h]", e.g. 100/2) -- and the
    # result is read as L/min. See the vignette Errata for the five-point
    # back-solve of Fig. 3b that pins CL_CVVHD to 1.72 L/h (CV 0.9%) against the
    # 1.715 L/h this parameterisation predicts.
    lkoa      <- log(0.0288); label("Hemodialyzer mass transfer-area coefficient (K0A) for Eq. 1, with BFR in mL/min and DFR in L/h")  # Huppe 2023 Table 2

    # Covariate effects (Huppe 2023 Table 2, "Covariates"; forms in footnotes a and b).
    # Linear, uncentered effect of measured urinary creatinine clearance on the
    # renal arm: CL_Renal = theta_CLRenal * (1 + theta_CLCR * CLCR). FIXED
    # (footnote c; RSE reported as "-").
    e_crcl_cl_renal <- fixed(0.0723); label("Linear effect of urinary creatinine clearance on renal clearance (per mL/min)")  # Huppe 2023 Table 2 (fixed, footnote c)

    # Linear increase of central volume with time since first dose:
    # V_C = theta_VC * (1 + theta_TSFD * TSFD) * exp(eta_VC), RSE 24%. TSFD is in
    # MINUTES -- Table 2's legend defines V_C as "the central compartment's volume
    # of distribution at time since first dose of 0 min", and every time in the
    # paper (120 min infusion, 15-360 min sampling, AUC in ug*min/mL, t1/2 in min)
    # is in minutes. model() converts the model's hour-based time accordingly.
    e_tsfd_vc <- 0.0008; label("Linear increase of central volume per minute since the first dose (per min)")  # Huppe 2023 Table 2

    # IIV. Table 2 reports "IIV V_C (%CV) = 62.10" for an exponential random
    # effect (Methods: "Interindividual variability was implemented using
    # exponential random effect models"), so omega^2 = log(CV^2 + 1)
    #   = log(0.6210^2 + 1) = 0.326287.
    # IIV was also tested on CL_Renal but dropped from the final model: including
    # CRCL on CL_Renal "completely characterized the IIV on CL_Renal" and the IIV
    # then "provided no longer a benefit of statistical significance".
    etalvc ~ 0.326287  # Huppe 2023 Table 2 (IIV V_C, 62.10 %CV)

    # Residual error: additive only (Results: "An additive error model was used to
    # explain residual variability"), reported directly as an SD in ug/mL.
    addSd <- 30.58; label("Additive residual error (ug/mL)")  # Huppe 2023 Table 2
  })

  model({
    # ---- 1. Derived covariate terms ------------------------------------------
    # Preserved-diuresis gate on the renal arm. Huppe 2023 Methods: intrinsic
    # elimination was fixed to zero in patients without preserved diuresis,
    # defined as residual diuresis > 100 mL/24 h. Table 2 footnote a prints the
    # ungated form, but the gate is required by both the Methods text and Fig. 3b,
    # where anuric patients (CRCL = 0) sit at exactly 100% dialysis contribution;
    # without the gate that point would be 1.715/(1.715 + 0.263) = 86.8%.
    diuresis <- (URINE_VOL_24H > 100)

    # Time since first dose, converted from the model's hours to the minutes in
    # which theta_TSFD was estimated. TSFD is taken to be the model's own time
    # variable, so an event table MUST place the first dose at time = 0 or the
    # central-volume trajectory will be wrong.
    tsfd_min <- t * 60

    # ---- 2. Individual parameters -------------------------------------------
    # Renal arm: linear, uncentered in measured urinary creatinine clearance,
    # switched off entirely in anuric patients (Huppe 2023 Table 2 footnote a).
    cl_renal <- diuresis * exp(lcl_renal) * (1 + e_crcl_cl_renal * CRCL)

    koa <- exp(lkoa)

    # CVVHD arm: Michaels hemodialyzer equation (Huppe 2023 Eq. 1),
    #   CL_CVVHD = BFR * (exp(z) - 1) / (exp(z) - BFR/DFR),  z = (K0A/BFR)(1 - BFR/DFR)
    # written here in the algebraically identical arrangement already used by
    # Liesenfeld_2013_dabigatran (same equation, same senior author) so the two
    # registered Michaels models read alike.
    #
    # BFR enters in mL/min (canonical) and DFR is converted from the canonical
    # mL/min back to the L/h of the study's own data column, because K0A = 0.0288
    # was calibrated in those mixed units; the resulting value is in L/min and is
    # scaled by 60 to the L/h used everywhere else in the model. See the lkoa
    # source-trace comment in ini() and the vignette Errata.
    dfr_lh <- DFR * 60 / 1000

    # bfr_safe / dfr_safe keep the rational expression finite when BFR and DFR
    # carry sentinel zeros outside a CVVHD period; the whole arm is then gated to
    # zero by the leading RRT_CRRT_ACTIVE multiplier. The fallbacks are the
    # study's protocol start settings (BFR 100 mL/min, DFR 2 L/h) rather than a
    # bare 1, because equal fallbacks would make the numerator AND denominator
    # vanish together (0/0 -> NaN, which then propagates through the gate and
    # silently blanks the whole solve). Note the expression is singular only when
    # bfr_safe == dfr_safe, which cannot occur for physiological settings since
    # BFR is in mL/min (order 100) and DFR in L/h (order 2).
    bfr_safe <- BFR * RRT_CRRT_ACTIVE + (1 - RRT_CRRT_ACTIVE) * 100
    dfr_safe <- dfr_lh * RRT_CRRT_ACTIVE + (1 - RRT_CRRT_ACTIVE) * 2

    michaels_exp    <- exp(-koa * (dfr_safe - bfr_safe) / (bfr_safe * dfr_safe))
    cl_cvvhd_lmin   <- bfr_safe * (1 - michaels_exp) / (1 - (bfr_safe / dfr_safe) * michaels_exp)
    cl_cvvhd        <- RRT_CRRT_ACTIVE * 60 * cl_cvvhd_lmin

    # Total clearance is the sum of the two arms (Huppe 2023 Table 2 footnote a).
    cl <- cl_renal + cl_cvvhd

    # Central volume grows linearly with time since the first dose.
    vc <- exp(lvc + etalvc) * (1 + e_tsfd_vc * tsfd_min)
    vp <- exp(lvp)
    q  <- exp(lq)

    # ---- 3. Micro-constants -------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ------------------------------------------------------
    # NOTE on the degenerate zero-clearance limit. An anuric patient (diuresis
    # gate 0) with CVVHD off has total clearance of EXACTLY zero -- a state the
    # source model genuinely predicts and simulates (Fig. 6, "0 mL/min" column,
    # black line; Discussion: anuric patients "showed hardly any elimination of
    # i.v. fosfomycin without CVVHD"). Concentrations still plateau rather than
    # diverge, because V_C grows with time since first dose. Under rxode2 5.1.3 a
    # model with a residual endpoint and exactly-zero elimination returns NaN
    # states when dosed by INFUSION (via either `dur` or `rate`); a bolus solves
    # correctly, as does the equivalent endpoint-free plain rxode2 model. To
    # simulate this limit with an infusion, set RRT_CRRT_ACTIVE to a negligible
    # non-zero value (1e-9), which reproduces the analytic cl -> 0 plateau; see
    # the vignette Errata.
    # Huppe 2023 Table 2 legend gives the concentration form
    #   dC1/dt = -(Q/VC) C1 + (Q/VP) C2 - (CL/VC) C1;  dC2/dt = (Q/VC) C1 - (Q/VP) C2
    # which is the standard two-compartment system; written here on amounts.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 5. Observation and error -------------------------------------------
    # Dose in mg and vc in L give mg/L, i.e. the ug/mL the paper reports.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
