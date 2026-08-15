Kobuchi_2025_dapagliflozin <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for oral",
    "dapagliflozin coupled to a turnover (indirect-response) Emax model",
    "for the HbA1c-lowering effect, in Japanese outpatients with type 2",
    "diabetes mellitus treated with 5 mg once daily for one year in",
    "routine clinical practice (long-term real-world data). Absorption",
    "rate constant and apparent volume of distribution were FIXED from a",
    "prior dapagliflozin population PK analysis (Melin 2024) because all",
    "study samples were trough concentrations; apparent clearance was",
    "estimated with a power effect of time-varying body weight. HbA1c is",
    "described by a zero-order production / first-order loss turnover",
    "pool whose loss rate is parameterised by the HbA1c half-life; drug",
    "effect is a fractional Emax inhibition of production scaled by the",
    "subject's own baseline HbA1c relative to a reference baseline of",
    "8.0% above a 5.0% normoglycemia floor. Two endpoints are declared:",
    "plasma dapagliflozin concentration (Cc, ng/mL) and HbA1c (%).",
    sep = " "
  )
  reference <- paste(
    "Kobuchi S, Sakai S, Terada R, Kato KI, Hayakawa T, Sakaeda T.",
    "Population Pharmacokinetic-pharmacodynamic Model Analysis of",
    "Dapagliflozin for HbA1c-lowering Effects in Japanese Patients with",
    "Type 2 Diabetes Mellitus using Long-term Real-world Data.",
    "Int J Med Sci. 2025;22(10):2333-2341. doi:10.7150/ijms.111519.",
    "The HbA1c turnover-with-Emax structure is inherited from",
    "de Winter W, Dunne A, de Trixhe XW, et al. Br J Clin Pharmacol.",
    "2017;83(5):1072-1081 (canagliflozin).",
    "ka and V/F were fixed from Melin J, Parkinson J, Hamren B, et al.",
    "Br J Clin Pharmacol. 2024;90(2):606-612.",
    sep = " "
  )
  vignette <- "Kobuchi_2025_dapagliflozin"

  # `hba1c` is not a canonical compartment name; it is the paper's HbA1c
  # turnover pool. Precedent: Baron_2016_empagliflozin.R declares the same
  # state as a paper-specific compartment.
  paper_specific_compartments <- c("hba1c")

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot   = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dapagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    hba1c   = list(analyte = "HbA1c", units = "%", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING. Kobuchi 2025 Results 'Population PK model':",
        "covariate screening first identified baseline body weight on",
        "CL/F, but because dapagliflozin lowers body weight the authors",
        "re-analysed with body weight at the blood-sampling time and",
        "found that the time-varying form was the significant covariate",
        "retained in the final model. Supply WT at every event and",
        "observation row. Normalised by 77.0 kg, the study median",
        "baseline body weight (Kobuchi 2025 Table 1; the CL/F equation",
        "printed in Table 2 uses BW/77.0). Observed sequential change",
        "over 12 months was median -1.5 kg (range -10.8 to +6.0 kg)."
      ),
      source_name        = "BW (Kobuchi 2025 Table 2 CL/F equation)"
    ),
    HBA1C = list(
      description        = "Baseline (pre-treatment) glycated hemoglobin",
      units              = "% (NGSP)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject BASELINE value, time-fixed. Kobuchi 2025 PD model",
        "section: 'The baseline HbA1c value for each patient was set to",
        "the corresponding observed value.' Used three ways inside",
        "model(): (1) as the initial condition of the hba1c turnover",
        "pool, (2) to derive the zero-order production rate at the",
        "assumed baseline steady state (kin = HBA1C * kout, from",
        "dH(0) = 0 in Kobuchi 2025 equation 2), and (3) as the baseline",
        "scaling factor (HBA1C - 5) / (8 - 5) of the drug effect in",
        "equation 3. Study baseline was 6.8% (SD 0.5, range 5.6-8.8;",
        "Kobuchi 2025 Table 1). NOT a covariate on a structural",
        "parameter -- it is a structural input of the turnover model."
      ),
      source_name        = "HbA1c at baseline (Kobuchi 2025 equations 2-3, Table 1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 85L,
    n_studies      = 1L,
    age_range      = "37-75 years",
    age_median     = "59 years",
    weight_range   = "49.0-118.0 kg",
    weight_median  = "77.0 kg",
    sex_female_pct = 28.2,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Japanese outpatients with type 2 diabetes mellitus at an early",
      "stage of diabetic nephropathy (urinary albumin-to-creatinine",
      "ratio < 30 mg/g Cr) or with inadequate glycemic control. Median",
      "duration of diabetes 10.5 years (range 1.0-33.0). Baseline HbA1c",
      "6.8% (SD 0.5, range 5.6-8.8), BMI 28.3 kg/m2 (SD 3.5), eGFR 74.3",
      "mL/min/1.73 m2 (range 39.2-137.8). Comorbidities: diabetic",
      "nephropathy 28/85, hypertension 59/85, hypercholesterolemia",
      "67/85 (Kobuchi 2025 Table 1)."
    ),
    dose_range     = paste(
      "Dapagliflozin 5 mg orally once daily after breakfast for 12",
      "months. This is the Japanese approved T2DM dose; escalation to",
      "10 mg/day is permitted in Japan but no patient in this cohort",
      "received it, so the 10 mg profiles in Kobuchi 2025 Figure 4 are",
      "model-based extrapolation, not observed data."
    ),
    regions        = "Japan (single centre, Tonami General Hospital, Toyama)",
    n_observations = "415 plasma dapagliflozin concentrations (all trough) and 508 HbA1c values",
    co_medication  = paste(
      "All patients received metformin (1000 mg/day or more) plus a",
      "dipeptidyl peptidase-4 inhibitor daily before dapagliflozin was",
      "started (Kobuchi 2025 Data sources)."
    ),
    notes          = paste(
      "Real-world observational single-centre single-arm cohort; visits",
      "at 1, 3, 6, 9 and 12 months with overnight fasting. Because every",
      "PK sample was a trough, ka and V/F are not identifiable from this",
      "dataset and were fixed from Melin 2024. The PK-PD model was fit",
      "sequentially: individual post hoc PK parameters from the PK model",
      "were fixed when the HbA1c parameters were estimated (Kobuchi 2025",
      "PD model section). Estimation used Phoenix NLME v8.4 (FOCE with",
      "the extended least squares method)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PK -- one-compartment, first-order absorption.
    # Kobuchi 2025 Table 2 (final population PK model).
    # -----------------------------------------------------------------

    lka <- fixed(log(57.4))
    label("Absorption rate constant ka (log 1/day)")
    # Kobuchi 2025 Table 2: ka = 57.4 1/day, CV% column reads "Fix".
    # PK model section: "ka and V/F were fixed based on the previous
    # report of population PK analysis [11]" (Melin 2024). Not
    # estimable here because all samples were troughs.

    lcl <- log(229.3)
    label("Apparent clearance CL/F at 77.0 kg (log L/day)")
    # Kobuchi 2025 Table 2: theta_CL/F = 229.3 L/day (CV% 1.5;
    # bootstrap median 229.6, 2.5-97.5th 222.6-236.4). Results text:
    # "The estimated population mean CL/F was 229.3 L/day (1.5% CV)".

    e_wt_cl <- 0.41
    label("Power exponent of (WT / 77.0) on CL/F (unitless)")
    # Kobuchi 2025 Table 2: theta_BW = 0.41 (CV% 21.3; bootstrap median
    # 0.40, 2.5-97.5th 0.22-0.58). Results text: "0.41 (2.5%-97.5%CI:
    # 0.24-0.58) is the exponential coefficient of body weight".
    # Table 2 prints the model as CL/F = theta_CL/F * (BW/77.0)^theta_BW.

    lvc <- fixed(log(73.9))
    label("Apparent volume of distribution V/F (log L)")
    # Kobuchi 2025 Table 2: V/F = 73.9 L, CV% column reads "Fix".
    # Fixed from Melin 2024 for the same reason as ka.

    # -----------------------------------------------------------------
    # PD -- HbA1c turnover with Emax drug effect on production.
    # Kobuchi 2025 Table 3 (final population PK-PD model), equations
    # 2-5. Structure inherited from de Winter 2017 (canagliflozin).
    # -----------------------------------------------------------------

    lthalfrec <- log(16.1)
    label("HbA1c half-life t1/2HbA1c (log day)")
    # Kobuchi 2025 Table 3: t1/2HbA1c = 16.1 day (CV% 4.1; bootstrap
    # median 16.0, 2.5-97.5th 15.3-16.5). Equation 5 defines the
    # turnover loss rate as kout = ln2 / t1/2HbA1c, so the estimated
    # parameter is the half-life and kout is derived in model().
    # Discussion notes this is shorter than de Winter 2017's 28.2 day.

    lemax <- log(0.034)
    label("Maximum HbA1c-lowering effect Emax (log HbA1c %/day)")
    # Kobuchi 2025 Table 3: Emax = 0.034 HbA1c %/day (CV% 3.1;
    # bootstrap median 0.034, 2.5-97.5th 0.031-0.035). Table 3 footnote:
    # "Emax, maximum HbA1c-lowering effect of dapagliflozin for a
    # typical patient with HbA1c at baseline of 8.0%" -- hence the
    # (HBA1C - 5) / (8 - 5) rescaling in equation 3. Emax is a RATE
    # (it is subtracted from kin in equation 2), not a concentration.

    lec50 <- log(23.7)
    label("Dapagliflozin concentration at half-maximal effect EC50 (log ng/mL)")
    # Kobuchi 2025 Table 3: EC50 = 23.7 ng/mL (CV% 5.8; bootstrap
    # median 21.9, 2.5-97.5th 13.5-24.4). Table 3 footnote: "EC50,
    # dapagliflozin plasma exposure at which half-maximal effect is
    # reached"; equation 4 uses the instantaneous C(t).

    # -----------------------------------------------------------------
    # Inter-individual variability. Kobuchi 2025 equation 1 is
    # P_i = P_pop * exp(eta_i) with eta ~ N(0, omega^2); Tables 2 and 3
    # head this block "Inter-individual variability (omega)" and label
    # the rows "omega CL/F (%)" / "omega t1/2 HbA1c (%)", i.e. the
    # tabulated percentage is omega x 100 (log-scale SD), not a
    # lognormal CV. Variance therefore = (pct / 100)^2. See the
    # vignette Errata for the reproduction check that supports this
    # reading against the paper's own simulated percentile bands.
    # -----------------------------------------------------------------

    etalcl ~ 0.019321
    # Kobuchi 2025 Table 2: omega CL/F = 13.9% (CV% 21.3; bootstrap
    # median 13.7, 2.5-97.5th 11.2-16.2). 0.139^2 = 0.019321.

    etalthalfrec ~ 1.079521
    # Kobuchi 2025 Table 3: omega t1/2HbA1c = 103.9% (CV% 11.2;
    # bootstrap median 104.1, 2.5-97.5th 101.7-106.4). 1.039^2 =
    # 1.079521. Results text: "The inter-individual variability of
    # t1/2HbA1c was estimated to be 103.9% (CV% = 11.2), which was
    # relatively high" -- the 11.2 there is the RSE column, not the IIV.
    # No IIV was reported on Emax or EC50: PD model section states
    # "We employed an inter-individual variability model for t1/2HbA1c"
    # (singular), and Table 3 lists no other omega.

    # -----------------------------------------------------------------
    # Residual variability.
    # -----------------------------------------------------------------

    propSd <- 0.398
    label("Proportional residual error on dapagliflozin concentration (fraction)")
    # Kobuchi 2025 Table 2: sigma = 39.8% (CV% 5.4; bootstrap median
    # 39.8, 2.5-97.5th 35.6-44.1). Results text: "A one-compartment
    # model with first-order absorption using multiplicative error
    # correction best described the time profiles" -- multiplicative
    # == proportional.

    addSd_hba1c <- 0.24
    label("Additive residual error on HbA1c (HbA1c %)")
    # Kobuchi 2025 Table 3: sigma = 0.24 HbA1c % (CV% 5.2; bootstrap
    # median 0.24, 2.5-97.5th 0.21-0.27). Results text: "described
    # using a previously reported turnover with the Emax model [8] with
    # an additive error correction". Units are absolute HbA1c percent.
  })

  model({
    # -----------------------------------------------------------------
    # 1. Individual PK parameters.
    # WT is time-varying (body weight at the sampling time); 77.0 kg is
    # the normalisation constant printed in the Kobuchi 2025 Table 2
    # CL/F equation and equals the study median baseline body weight
    # (Table 1).
    # -----------------------------------------------------------------
    ref_wt <- 77.0

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc)

    kel <- cl / vc

    # -----------------------------------------------------------------
    # 2. PK ODEs. One compartment, first-order absorption, no lag,
    # F fixed at 1 (CL and V are apparent, CL/F and V/F).
    # -----------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Unit conversion: amounts are mg and vc is L, so central / vc is
    # mg/L. The paper reports concentrations and EC50 in ng/mL, and
    # 1 mg/L = 1000 ng/mL. Keeping Cc in ng/mL is required for the
    # equation-4 Emax term to use the published EC50 directly.
    mgL_to_ngmL <- 1000
    Cc <- central / vc * mgL_to_ngmL

    # -----------------------------------------------------------------
    # 3. HbA1c turnover with Emax effect (Kobuchi 2025 equations 2-5).
    #
    #   eq 5:  kout = ln2 / t1/2HbA1c
    #   eq 4:  Efc  = Emax * C(t) / (C(t) + EC50)
    #   eq 3:  Ef   = Efc * (H(0) - 5) / (8 - 5)
    #   eq 2:  dH/dt = kin - Ef - kout * H(t)
    #
    # The 5.0 floor is the paper's normoglycemia correction ("the
    # HbA1c level was corrected for a lower boundary of 5.0% because
    # normoglycemia is typically associated with that value") and 8.0
    # is the reference baseline at which Emax is defined (Table 3
    # footnote). Both are structural constants of the published model,
    # not estimated parameters.
    #
    # kin: the paper states "At baseline, dH(0)=0 was assumed and kin
    # was estimated as kin = H(0)/kout". Taken literally that is
    # dimensionally impossible -- kin must be HbA1c %/day, and
    # H(0)/kout has units of % * day. Setting dH/dt = 0 at t = 0 in
    # equation 2 (where C(0) = 0, so Ef = 0) gives kin = kout * H(0),
    # which is both dimensionally correct and the standard turnover
    # baseline identity. The printed equation is followed over the
    # prose; see the vignette Errata.
    # -----------------------------------------------------------------
    hba1c_floor <- 5.0
    hba1c_ref   <- 8.0

    thalfrec <- exp(lthalfrec + etalthalfrec)
    kout     <- log(2) / thalfrec
    kin      <- HBA1C * kout

    emax <- exp(lemax)
    ec50 <- exp(lec50)

    efc <- emax * Cc / (Cc + ec50)
    ef  <- efc * (HBA1C - hba1c_floor) / (hba1c_ref - hba1c_floor)

    d/dt(hba1c) <- kin - ef - kout * hba1c
    hba1c(0)    <- HBA1C

    # -----------------------------------------------------------------
    # 4. Endpoints.
    # -----------------------------------------------------------------
    Cc    ~ prop(propSd)
    hba1c ~ add(addSd_hba1c)
  })
}
