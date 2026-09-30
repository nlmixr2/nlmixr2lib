Otto_2021_s_ketamine <- function() {
  description <- "Joint parent-metabolite population PK model of intravenous S-ketamine and its metabolite S-norketamine in healthy adult volunteers (Otto 2021). S-ketamine has three-compartment disposition with the population parameters fixed to the Fanta 2015 model (V1 = 133 L, CL = 95.2 L/h, Vp1 = 187 L, Q1 = 23.2 L/h, Vp2 = 98.8 L, Q2 = 157 L/h at 70 kg); inter-individual variability on V1 and CL was re-estimated. All S-ketamine clearance is assumed to form S-norketamine, which has a newly estimated two-compartment disposition (V = 98.6 L, CL = 57.7 L/h, Vp = 160 L, Q = 42.8 L/h at 70 kg) with no transit compartments. Central volumes and clearances of both analytes are scaled allometrically to 70 kg (exponent 1 for volumes, 0.75 for clearances); residual error is proportional for each analyte."
  reference <- paste(
    "Otto ME, Bergmann KR, Jacobs G, van Esdonk MJ.",
    "Predictive performance of parent-metabolite population pharmacokinetic",
    "models of (S)-ketamine in healthy volunteers.",
    "Eur J Clin Pharmacol. 2021;77(8):1181-1192. doi:10.1007/s00228-021-03104-1.",
    "S-ketamine population parameters from Fanta S, Kinnunen M, Backman JT,",
    "Kalso E. Population pharmacokinetics of S-ketamine and norketamine in",
    "healthy volunteers after intravenous and oral dosing.",
    "Eur J Clin Pharmacol. 2015;71(4):441-447. doi:10.1007/s00228-015-1826-y.",
    sep = " "
  )
  vignette <- "Otto_2021_s_ketamine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight (kg).",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to a 70 kg reference (Otto 2021 Table 3 footnote and supplementary NONMEM code): exponent 1 on the central volumes and 0.75 on the elimination clearances of both S-ketamine and S-norketamine. Peripheral volumes and intercompartmental clearances are NOT weight-scaled in the as-run code.",
      source_name = "WGHT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "S-ketamine", units = "mg", specimen = "plasma", verified = TRUE),
    central_snk = list(
      analyte = "S-norketamine",
      units = "mg (S-ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_snk = list(
      analyte = "S-norketamine",
      units = "mg (S-ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 2L,
    age_range = "adults; mean 23.0 (SD 3.6) years in CHDR1311 and 23.6 (SD 5.1) years in CHDR1016",
    weight_range = "mean 68.0 (SD 7.2) kg in CHDR1311 and 71.3 (SD 8.5) kg in CHDR1016",
    bmi = "mean 21.6 (SD 2.0) kg/m^2 in CHDR1311 and 22.4 (SD 2.0) kg/m^2 in CHDR1016",
    sex_female_pct = 45.8,
    disease_state = "Healthy volunteers",
    dose_range = "CHDR1311: 10 mg S-ketamine as a 30-min IV infusion. CHDR1016: 2-h stepped IV infusion on two occasions (low and high, high = double the low); high-occasion post-amendment regimen 0.026 mg/kg bolus then 0.425 mg/kg/h (0-14 min), 0.275 mg/kg/h (15-39 min) and 0.15 mg/kg/h (40-120 min), with female rates increased by 5-15 percent.",
    regions = "Centre for Human Drug Research, Leiden, The Netherlands",
    n_observations = "268 samples in CHDR1311 (0 percent BLQ) and 864 samples in CHDR1016 (5.9 percent BLQ); BLQ samples were removed",
    sampling = "Venous",
    notes = "Demographics from Otto 2021 Table 1 (CHDR1311 N = 17, 9 male; CHDR1016 N = 31, 17 male). The S-ketamine population parameters come from Fanta 2015 (11 healthy male volunteers, IV bolus and oral S-ketamine); only the S-ketamine IIV was re-estimated on the CHDR data. The S-norketamine structural parameters, IIV and all residual errors were estimated on CHDR1311 + CHDR1016."
  )

  ini({
    # ---- S-ketamine: population parameters fixed to Fanta 2015 ----
    # Otto 2021 Results 'Model redevelopment': 'Population parameters of the
    # (S)-ketamine model of Fanta et al. were fixed to their reported values'.
    # Table 3 prints them without RSE; the supplementary NONMEM code marks
    # THETA(1)-(6) FIX.
    lvc <- fixed(log(133)); label("S-ketamine central volume V1 at 70 kg (L)") # Table 3 (S)-ketamine Vd = 133; supplement THETA(1) 133.0 FIX
    lcl <- fixed(log(95.2)); label("S-ketamine clearance CL at 70 kg (L/h)") # Table 3 (S)-ketamine CL = 95.2; supplement THETA(2) 95.2 FIX
    lvp <- fixed(log(187)); label("S-ketamine first peripheral volume Vp1 (L)") # Table 3 (S)-ketamine Vp1 = 187; supplement THETA(3) 187.0 FIX
    lq <- fixed(log(23.2)); label("S-ketamine first intercompartmental clearance Q1 (L/h)") # Table 3 (S)-ketamine Q1 = 23.2; supplement THETA(4) 23.2 FIX
    lvp2 <- fixed(log(98.8)); label("S-ketamine second peripheral volume Vp2 (L)") # Table 3 (S)-ketamine Vp2 = 98.8; supplement THETA(5) 98.8 FIX
    lq2 <- fixed(log(157)); label("S-ketamine second intercompartmental clearance Q2 (L/h)") # Table 3 (S)-ketamine Q2 = 157; supplement THETA(6) 157.0 FIX

    # ---- S-norketamine: population parameters estimated in Otto 2021 ----
    lvc_snk <- log(98.638); label("S-norketamine central volume at 70 kg (L)") # Table 3 (S)-norketamine Vd = 98.6 (RSE 3.77%); supplement THETA(7) 98.638
    lcl_snk <- log(57.72); label("S-norketamine clearance at 70 kg (L/h)") # Table 3 (S)-norketamine CL = 57.7 (RSE 5.7%); supplement THETA(8) 57.72
    lvp_snk <- log(160.04); label("S-norketamine peripheral volume (L)") # Table 3 (S)-norketamine Vp1 = 160 (RSE 17.6%); supplement THETA(9) 160.04
    lq_snk <- log(42.778); label("S-norketamine intercompartmental clearance (L/h)") # Table 3 (S)-norketamine Q1 = 42.8 (RSE 7.51%); supplement THETA(10) 42.778

    # ---- Allometric exponents (fixed; reference weight 70 kg) ----
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT/70) on S-ketamine V1 (unitless)") # Table 3 footnote: exponent 1 for Vd; supplement (WGHT/70)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on S-ketamine CL (unitless)") # Table 3 footnote: exponent 0.75 for CL; supplement (WGHT/70)**(3/4)
    e_wt_vc_snk <- fixed(1); label("Allometric exponent of (WT/70) on S-norketamine V (unitless)") # Table 3 footnote: exponent 1 for Vd; supplement (WGHT/70)
    e_wt_cl_snk <- fixed(0.75); label("Allometric exponent of (WT/70) on S-norketamine CL (unitless)") # Table 3 footnote: exponent 0.75 for CL; supplement (WGHT/70)**(3/4)

    # ---- IIV (log-scale variances; exponential IIV in the NONMEM code) ----
    # Supplement '$OMEGA BLOCK (2)' order is Vd then CL for each analyte.
    etalvc + etalcl ~ c(0.084177, 0.02261, 0.026048) # Table 3 (S)-ketamine omega2 Vd 0.084, cov 0.023, CL 0.026; supplement OMEGA BLOCK(2)
    etalvc_snk + etalcl_snk ~ c(0.040426, 0.04456, 0.10319) # Table 3 (S)-norketamine omega2 Vd 0.040, cov 0.044, CL 0.103; supplement OMEGA BLOCK(2)

    # ---- Residual error (proportional; sigma^2 converted to SD) ----
    propSd <- 0.23609; label("S-ketamine proportional residual SD (fraction)") # Table 3 sigma2 = 0.056; supplement $SIGMA 0.055736 -> sqrt = 0.23609
    propSd_snk <- 0.14261; label("S-norketamine proportional residual SD (fraction)") # Table 3 sigma2 = 0.020; supplement $SIGMA 0.020338 -> sqrt = 0.14261
  })

  model({
    # Individual parameters (supplementary NONMEM $PK)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vp <- exp(lvp)
    q <- exp(lq)
    vp2 <- exp(lvp2)
    q2 <- exp(lq2)

    vc_snk <- exp(lvc_snk + etalvc_snk) * (WT / 70)^e_wt_vc_snk
    cl_snk <- exp(lcl_snk + etalcl_snk) * (WT / 70)^e_wt_cl_snk
    vp_snk <- exp(lvp_snk)
    q_snk <- exp(lq_snk)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_snk <- cl_snk / vc_snk
    k12_snk <- q_snk / vc_snk
    k21_snk <- q_snk / vp_snk

    # S-ketamine: IV dose into central. All S-ketamine clearance forms
    # S-norketamine (fraction metabolised assumed 1, Otto 2021 Methods
    # 'Model redevelopment'); the flux is carried in S-ketamine mass with no
    # molecular-weight conversion, exactly as in the supplementary $DES.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(central_snk) <- kel * central - kel_snk * central_snk - k12_snk * central_snk + k21_snk * peripheral1_snk
    d/dt(peripheral1_snk) <- k12_snk * central_snk - k21_snk * peripheral1_snk

    # Venous plasma concentrations (mg/L; multiply by 1000 for ug/L = ng/mL)
    Cc <- central / vc
    Cc_snk <- central_snk / vc_snk

    Cc ~ prop(propSd)
    Cc_snk ~ prop(propSd_snk)
  })
}
