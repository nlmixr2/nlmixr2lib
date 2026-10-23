Reed_2025_meropenem_vaborbactam <- function() {
  description <- paste(
    "Joint population PK model for intravenous meropenem and vaborbactam",
    "(1:1 fixed-ratio combination) in 18 acutely ill hospitalised adults with",
    "varying renal function, three of them on continuous renal replacement",
    "therapy (CRRT) (Reed 2025). Pmetrics NPAG non-parametric fit. Each drug",
    "runs through its own copy of a three-state system: the patient's",
    "pre-filter blood (central, the sampled matrix), the CRRT hemofilter",
    "(circuit_plasma) and the post-filter return limb (circuit_return).",
    "Off CRRT the circuit states are disconnected and each drug is",
    "one-compartment with clearance proportional to BSA-normalised",
    "Cockcroft-Gault creatinine clearance (reference 120 mL/min/1.73 m^2);",
    "on CRRT the body clearance is a separate residual clearance and the",
    "filter removes drug at a rate proportional to the effluent flow",
    "standardised to 2000 mL/h. Volume scales linearly with body weight",
    "(reference 80 kg). Meropenem uses the unsuffixed canonical names;",
    "vaborbactam carries the sibling-drug suffix _vab."
  )
  reference <- paste(
    "Reed G, Valadez A, Smith BJ, McCreary EK, Kline EG, Neely MN, Squires KM,",
    "Shields RK, Rhodes NJ. Population pharmacokinetics of",
    "meropenem-vaborbactam in acutely ill hospitalized patients with various",
    "degrees of renal dysfunction including continuous renal replacement",
    "therapy. Antimicrob Agents Chemother. 2025;69(7):e00108-25.",
    "doi:10.1128/aac.00108-25",
    sep = " "
  )
  vignette <- "Reed_2025_meropenem_vaborbactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_plasma = list(analyte = "meropenem", units = "mg", specimen = "whole blood", verified = TRUE),
    circuit_return = list(analyte = "meropenem", units = "mg", specimen = "whole blood", verified = TRUE),
    central_vab = list(analyte = "vaborbactam", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_plasma_vab = list(analyte = "vaborbactam", units = "mg", specimen = "whole blood", verified = TRUE),
    circuit_return_vab = list(analyte = "vaborbactam", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  dosing <- c("central", "central_vab")

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance normalised to 1.73 m^2 body surface area",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods, 'Design and patients': CrCL 'estimated using the Cockcroft-Gault method",
        "and normalized to a body surface area of 1.73 m2'; Results, 'Model development':",
        "the BSA-normalised form fitted better than raw mL/min (AIC 7492.841 -> 7288.027).",
        "Enters the off-CRRT clearance as the through-origin ratio CRCL / 120 (Equation 1;",
        "Table S1 model 5 'CrCL/120 mL/min/1.73m2 on CL'). Not read while",
        "RRT_CRRT_ACTIVE = 1. Baseline cohort mean 123.3 +/- 103.7 (Table 1)."
      ),
      source_name = "CrCL"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters both volumes as the linear ratio WT / 80 (Equation 2; Table S1 model 4",
        "'WT/80 kg on V'). Cohort mean 90.7 +/- 32.6 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Continuous renal replacement therapy running (1) or not (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CRRT not running)",
      notes = paste(
        "The 'CRRT' term of Equation 1, and the gate on the circuit transfer constants:",
        "Methods, 'Model development', K12, K23 and K31 move drug through the circuit",
        "'when CRRT was active. If CRRT was not active, elimination occurred only via",
        "K10'. Carried as the time-varying ACTIVE gate rather than the subject-level",
        "RRT_CRRT_STATUS flag because the switch is structural: it connects and",
        "disconnects the circuit states. For the paper's cohort, in which the three",
        "CRRT patients were on CRRT throughout sampling, the two are numerically",
        "identical. 3/18 patients (16.7%; Table 1); modality not reported beyond",
        "the dialysate plus replacement flows of Results, 'Patient characteristics'."
      ),
      source_name = "CRRT"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "CRRT effluent flow rate (sum of ultrafiltration and dialysate flow rates)",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods, 'Model development': 'The effluent flow rate (i.e., the sum of",
        "ultrafiltration and dialysate flow rates) was incorporated into the model using",
        "the term FLOW, which was standardized to the median flow rate (mL/h)'. The",
        "median is 2,000 mL/h (Results, 'Patient characteristics'), so FLOW =",
        "RRT_CRRT_EFFLUENT_FLOW / 2000. Multiplies the filter elimination constant",
        "K20 only (equation 2). Supply 0 (not NA) for subjects off CRRT. The paper's",
        "Monte Carlo simulations used 1,750 mL/h (25 mL/kg/h for a 70 kg patient)."
      ),
      source_name = "FLOW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    n_observations = "83 plasma samples, each assayed for both meropenem and vaborbactam",
    age_range = "mean 54.7 +/- 14.4 years",
    weight_range = "mean 90.7 +/- 32.6 kg",
    sex_female_pct = 55.6,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Acutely ill hospitalised adults treated with meropenem-vaborbactam for",
      "suspected or confirmed carbapenem-resistant Enterobacterales infection.",
      "Two patients (11.1%) had acute kidney injury without CRRT; three (16.7%)",
      "required CRRT. Patients on intermittent hemodialysis, peritoneal or",
      "low-efficiency dialysis were excluded."
    ),
    renal_function = paste(
      "Baseline Cockcroft-Gault CrCL 123.3 +/- 103.7 mL/min/1.73 m^2; serum",
      "creatinine 0.99 +/- 0.74 mg/dL (Table 1). CRRT patients: median effluent",
      "flow 2,000 mL/h, blood flow 250 mL/min, dialysate flow 1,500 mL/h,",
      "pre- and post-replacement flow 250 mL/h."
    ),
    dose_range = paste(
      "Meropenem-vaborbactam per institutional renal-dosing protocol, 3-hour IV",
      "infusions; at least 48 h of therapy before sampling."
    ),
    regions = "United States (University of Pittsburgh Medical Center, December 2017 - July 2022).",
    notes = paste(
      "Demographics from Table 1 (mean +/- SD). Samples pre-dose and 1.5, 3, 4, 6,",
      "8 and 12 h post-dose when feasible; pre-filter blood in CRRT patients.",
      "Total concentrations by LC-MS/MS (linear 0.1-100 ug/mL for both analytes).",
      "Fitted with NPAG in Pmetrics 2.1.1; parameters are the weighted medians of",
      "the non-parametric population distribution (Table 2)."
    )
  )

  ini({
    # Structural parameters: Reed 2025 Table 2, 'Median (weighted)' column of
    # the Pmetrics NPAG non-parametric population distribution. The 95%
    # credible intervals and shrinkage printed alongside are noted per line.
    #
    # Map of the paper's symbols (Methods, 'Model development', equations
    # 1-3; Results Equations 1-2) to the canonical names used here:
    #   CL_nonCRRT -> lcl            K12 -> lk_central_circuit
    #   CL_CRRT    -> lcl_residual   K23 -> lk_circuit_return
    #   V1 (Vd)    -> lvc            K31 -> lk_return_central
    #                                K20 -> lk_circuit_effluent
    # Table 2 labels the on-CRRT clearance 'CL residual'; Equation 1 calls
    # the same quantity CL_CRRT. It is the patient's own clearance (the K10
    # route) while the circuit runs, not the circuit's clearance, which is
    # carried by the circuit states.

    # ---- Meropenem ----
    lcl <- log(10.6)
    label("Meropenem off-CRRT clearance at CRCL = 120 mL/min/1.73 m^2 (L/h)") # Table 2 MEM 'CL non-CRRT (L/h)' 10.6 (95% CrI 9.13-12.5; shrinkage 10.4%)
    lcl_residual <- log(2.1)
    label("Meropenem residual body clearance while on CRRT (L/h)") # Table 2 MEM 'CL residual (L/h)' 2.1 (95% CrI 1.24-2.85; shrinkage 19.2%)
    lvc <- log(27.4)
    label("Meropenem volume of distribution at WT = 80 kg (L)") # Table 2 MEM 'Vd (L)' 27.4 (95% CrI 23.2-38.4; shrinkage 9.1%)
    lk_central_circuit <- log(0.77)
    label("Meropenem transfer rate constant, pre-filter blood to CRRT filter (1/h)") # Table 2 MEM 'K12 (h-1)' 0.77 (95% CrI 0.43-1.16; shrinkage 10.9%)
    lk_circuit_effluent <- log(4.44)
    label("Meropenem filter elimination rate constant at the median effluent flow (1/h)") # Table 2 MEM 'K20 (h-1)' 4.44 (95% CrI 2.21-5.50; shrinkage 11.8%)
    lk_circuit_return <- log(29.3)
    label("Meropenem transfer rate constant, CRRT filter to post-filter blood (1/h)") # Table 2 MEM 'K23 (h-1)' 29.3 (95% CrI 20.6-31.9; shrinkage 9.9%)
    lk_return_central <- log(1.67)
    label("Meropenem transfer rate constant, post-filter blood to pre-filter blood (1/h)") # Table 2 MEM 'K31 (h-1)' 1.67 (95% CrI 0.51-35.0; shrinkage 24.6%)

    # ---- Vaborbactam ----
    lcl_vab <- log(6.64)
    label("Vaborbactam off-CRRT clearance at CRCL = 120 mL/min/1.73 m^2 (L/h)") # Table 2 VAB 'CL non-CRRT (L/h)' 6.64 (95% CrI 3.99-7.84; shrinkage 10.2%)
    lcl_residual_vab <- log(1.15)
    label("Vaborbactam residual body clearance while on CRRT (L/h)") # Table 2 VAB 'CL residual (L/h)' 1.15 (95% CrI 0.90-1.42; shrinkage 12.0%)
    lvc_vab <- log(24.3)
    label("Vaborbactam volume of distribution at WT = 80 kg (L)") # Table 2 VAB 'Vd (L)' 24.3 (95% CrI 19.7-31.1; shrinkage 9.1%)
    lk_central_circuit_vab <- log(1.45)
    label("Vaborbactam transfer rate constant, pre-filter blood to CRRT filter (1/h)") # Table 2 VAB 'K12 (h-1)' 1.45 (95% CrI 0.97-2.70; shrinkage 20.0%)
    lk_circuit_effluent_vab <- log(0.59)
    label("Vaborbactam filter elimination rate constant at the median effluent flow (1/h)") # Table 2 VAB 'K20 (h-1)' 0.59 (95% CrI 0.18-0.66; shrinkage 13.5%)
    lk_circuit_return_vab <- log(22.2)
    label("Vaborbactam transfer rate constant, CRRT filter to post-filter blood (1/h)") # Table 2 VAB 'K23 (h-1)' 22.2 (95% CrI 20.3-26.2; shrinkage 22.0%)
    lk_return_central_vab <- log(4.31)
    label("Vaborbactam transfer rate constant, post-filter blood to pre-filter blood (1/h)") # Table 2 VAB 'K31 (h-1)' 4.31 (95% CrI 2.56-30.7; shrinkage 27.7%)

    # Between-subject variability. NPAG estimates a discrete non-parametric
    # joint distribution, and every one of the 14 parameters above carried
    # one (Table 2 reports a shrinkage for each). Table 2 summarises the
    # distributions only by their weighted median and a 95% credible
    # interval: no SD, CV% or covariance is reported anywhere in the paper
    # or Supplement Table S1. The intervals are strongly asymmetric (MEM K31
    # 0.51-35.0 around 1.67) and too narrow on clearance (MEM CL non-CRRT
    # 9.13-12.5) to be the spread of an ICU population, so no omega can be
    # read off them. The random effects are therefore declared and carried
    # as fixed(0) rather than invented; see the vignette 'Assumptions and
    # deviations'.
    etalcl ~ fixed(0)
    etalcl_residual ~ fixed(0)
    etalvc ~ fixed(0)
    etalk_central_circuit ~ fixed(0)
    etalk_circuit_effluent ~ fixed(0)
    etalk_circuit_return ~ fixed(0)
    etalk_return_central ~ fixed(0)
    etalcl_vab ~ fixed(0)
    etalcl_residual_vab ~ fixed(0)
    etalvc_vab ~ fixed(0)
    etalk_central_circuit_vab ~ fixed(0)
    etalk_circuit_effluent_vab ~ fixed(0)
    etalk_circuit_return_vab ~ fixed(0)
    etalk_return_central_vab ~ fixed(0)

    # Residual error. The Pmetrics assay-error polynomial and the gamma /
    # lambda terms are not reported in the paper or its supplement, so both
    # components are carried as fixed(0) for each analyte.
    propSd <- fixed(0)
    label("Meropenem proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd <- fixed(0)
    label("Meropenem additive residual SD (mg/L; 0 -- not reported in the source)")
    propSd_vab <- fixed(0)
    label("Vaborbactam proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd_vab <- fixed(0)
    label("Vaborbactam additive residual SD (mg/L; 0 -- not reported in the source)")
  })
  model({
    # Reference values printed in Results Equations 1-2 (CrCL / 120, WT / 80)
    # and the FLOW standardisation of Methods, 'Model development' (median
    # effluent flow 2,000 mL/h, Results, 'Patient characteristics').
    crcl_ref <- 120
    wt_ref <- 80
    flow_ref <- 2000
    flow_std <- RRT_CRRT_EFFLUENT_FLOW / flow_ref

    # Body clearance (the K10 route), Results Equation 1:
    #   CL = CL_nonCRRT * (CrCL/120) * (1 - CRRT) + CL_CRRT * CRRT
    # The gated sum is assigned to `cl` itself (see the RRT_CRRT_ACTIVE entry
    # of covariate-columns.md); the two arms sit on their own lines so each
    # eta stays mu-referenced.
    cl_off <- exp(lcl + etalcl) * (CRCL / crcl_ref)
    cl_residual <- exp(lcl_residual + etalcl_residual)
    cl <- cl_off * (1 - RRT_CRRT_ACTIVE) + cl_residual * RRT_CRRT_ACTIVE
    cl_off_vab <- exp(lcl_vab + etalcl_vab) * (CRCL / crcl_ref)
    cl_residual_vab <- exp(lcl_residual_vab + etalcl_residual_vab)
    cl_vab <- cl_off_vab * (1 - RRT_CRRT_ACTIVE) + cl_residual_vab * RRT_CRRT_ACTIVE

    # Volume, Results Equation 2: Vd = V1 * (WT / 80).
    vc <- exp(lvc + etalvc) * (WT / wt_ref)
    vc_vab <- exp(lvc_vab + etalvc_vab) * (WT / wt_ref)

    kel <- cl / vc
    kel_vab <- cl_vab / vc_vab

    # Circuit transfer constants. K12, K23 and K31 act only 'when CRRT was
    # active' (Methods, 'Model development'), so each is gated; K20 enters
    # the filter equation multiplied by FLOW, which is 0 off CRRT.
    k_central_circuit <- exp(lk_central_circuit + etalk_central_circuit) * RRT_CRRT_ACTIVE
    k_circuit_return <- exp(lk_circuit_return + etalk_circuit_return) * RRT_CRRT_ACTIVE
    k_return_central <- exp(lk_return_central + etalk_return_central) * RRT_CRRT_ACTIVE
    k_circuit_effluent <- exp(lk_circuit_effluent + etalk_circuit_effluent)

    k_central_circuit_vab <- exp(lk_central_circuit_vab + etalk_central_circuit_vab) * RRT_CRRT_ACTIVE
    k_circuit_return_vab <- exp(lk_circuit_return_vab + etalk_circuit_return_vab) * RRT_CRRT_ACTIVE
    k_return_central_vab <- exp(lk_return_central_vab + etalk_return_central_vab) * RRT_CRRT_ACTIVE
    k_circuit_effluent_vab <- exp(lk_circuit_effluent_vab + etalk_circuit_effluent_vab)

    # Methods, 'Model development', equations 1-3, one copy per drug:
    #   dX[1] = Rate_in - (K12 + K10) * X[1] + K31 * X[3]
    #   dX[2] = K12 * X[1] - (K23 + FLOW * K20) * X[2]
    #   dX[3] = K23 * X[2] - K31 * X[3]
    # X[1] pre-filter blood = central; X[2] CRRT filter = circuit_plasma;
    # X[3] post-filter blood = circuit_return. The circuit states hold
    # amounts only; the paper assigns them no volume.
    d/dt(central) <- -(k_central_circuit + kel) * central + k_return_central * circuit_return
    d/dt(circuit_plasma) <- k_central_circuit * central -
      (k_circuit_return + flow_std * k_circuit_effluent) * circuit_plasma
    d/dt(circuit_return) <- k_circuit_return * circuit_plasma - k_return_central * circuit_return

    d/dt(central_vab) <- -(k_central_circuit_vab + kel_vab) * central_vab +
      k_return_central_vab * circuit_return_vab
    d/dt(circuit_plasma_vab) <- k_central_circuit_vab * central_vab -
      (k_circuit_return_vab + flow_std * k_circuit_effluent_vab) * circuit_plasma_vab
    d/dt(circuit_return_vab) <- k_circuit_return_vab * circuit_plasma_vab -
      k_return_central_vab * circuit_return_vab

    # Total plasma concentrations (mg/L): dose in mg, volume in L. Pre-filter
    # blood was the sampled matrix in CRRT patients (Methods, 'Dosing and
    # sample collection').
    Cc <- central / vc
    Cc_vab <- central_vab / vc_vab

    Cc ~ add(addSd) + prop(propSd)
    Cc_vab ~ add(addSd_vab) + prop(propSd_vab)
  })
}
