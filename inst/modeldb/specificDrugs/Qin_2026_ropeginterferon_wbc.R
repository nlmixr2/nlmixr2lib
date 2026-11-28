Qin_2026_ropeginterferon_wbc <- function() {
  description <- paste0(
    "Sequential PK-PD model of white blood cell count for subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 179 patients with essential ",
    "thrombocythemia (SURPASS-ET and EXCEED-ET). The fixed ",
    "quasi-equilibrium TMDD pharmacokinetic layer of ",
    "Qin_2026_ropeginterferon (patient branch) drives a sigmoid Imax ",
    "inhibition of white cell production by total serum ropeg in an ",
    "indirect-response model; from the start of treatment production is ",
    "further scaled by (1 - EP), an empirical fractional reduction. ",
    "Baseline JAK2 V617F allele burden (power) and Asian race (1.34-fold) ",
    "raise the initial white cell count."
  )
  reference <- paste(
    "Qin A, Masarova L, Chen H, Xiao Z, Jin J, Gill H, Reeves BN,",
    "El Chaer F, Kirito K, Komatsu N, Abu-Zeinah G, Shimoda K,",
    "Zagrijtschuk O, Sato T, Wu D, Wu L, Su X, Gao Y, Cui J, Zhang L,",
    "Mesa RA.",
    "Exposure-Response Analyses of Ropeginterferon Alfa-2b in",
    "Essential Thrombocythemia.",
    "J Clin Pharmacol. 2026;66(10):e70258. doi:10.1002/jcph.70258.",
    sep = " "
  )
  vignette <- "Qin_2026_ropeginterferon"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ug/L total serum ropeg (Cc, = ng/mL); white blood cell count circ_wbc in 10^9/L"
  )

  compartmentData <- list(
    depot = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "ropeginterferon alfa-2b (total, free plus target-bound)",
      units = "ug",
      specimen = "serum",
      verified = TRUE
    ),
    total_target = list(
      analyte = "ropeg target (total, free plus drug-bound binding capacity)",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    circ_wbc = list(analyte = "white blood cells", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Inherited pharmacokinetic covariate: CL = CL_TV * (BMI / 24)^1.03 (Qin 2026 Table 2). See Qin_2026_ropeginterferon.",
      source_name = "BMI"
    ),
    TUM_JAK2_V617F_VAF = list(
      description = "Baseline JAK2 V617F variant allele frequency (allele burden).",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the initial white cell count:",
        "WBC0_i = WBC0_TV * (JAK2 / 13.7)^0.178 (Qin 2026 Table 3",
        "'JAK2 on WBC0' 0.178; control stream WBCJAK2 =",
        "THETA(9)*LOG(JAK2I/13.7)). The control stream replaces any value",
        "<= 0 with 13.7 (IF(JAK2I.LE.0) JAK2I = 13.7), so JAK2",
        "V617F-negative patients (CALR- or MPL-driven disease; the EXCEED-ET",
        "median is 0%) and patients with a missing value (6 in EXCEED-ET)",
        "sit at the reference and receive no effect. model() reproduces",
        "that substitution: supply 0 for a JAK2 V617F-negative patient."
      ),
      source_name = "JAK2"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator; 1 = Asian, 0 = Caucasian or other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or other race)",
      notes = paste(
        "Multiplies the initial white cell count by 1.34 (Qin 2026",
        "Table 3 'Asian on WBC0'). The control stream sets WBCRACE =",
        "THETA(11) for RACE = 4 and THETA(10) = 0 FIX for RACE = 2, 3",
        "and 5, inside the exponent of WBC0, so the table value is the",
        "fold change exp(THETA(11)); see the vignette for the check",
        "against the SURPASS-ET / EXCEED-ET baseline counts."
      ),
      source_name = "RACE (= 4)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 179L,
    n_studies = 2L,
    age_range = "median 57.0 (EXCEED-ET) and 61.0 (SURPASS-ET) years, range 21.0-84.0 (Qin 2026 Table 1)",
    weight_range = "median 76.9 (EXCEED-ET) and 60.0 (SURPASS-ET) kg, range 41.5-153 (Qin 2026 Table 1)",
    sex_female_pct = 57.5,
    race_ethnicity = "SURPASS-ET 94.5% Asian; EXCEED-ET 77.3% Caucasian, 8.0% Asian, 14.8% other (Qin 2026 Table 1)",
    disease_state = "Essential thrombocythemia; baseline white cell count median 11.6 (SURPASS-ET) and 7.47 (EXCEED-ET) x 10^9/L; baseline JAK2 V617F allele burden median 31.7% (SURPASS-ET) and 0% (EXCEED-ET)",
    dose_range = "HIDAT regimen: 250 ug subcutaneously at week 0, 350 ug at week 2, 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    notes = "White blood cell counts measured every 2 weeks; counts were not followed after treatment discontinuation or interruption (Qin 2026 Discussion)."
  )

  ini({
    # ==================================================================
    # PHARMACOKINETIC LAYER -- ENTIRELY FIXED (sequential PK-PD: the
    # control stream reads CL, VC, KA, ALAG1, RTOT0, RTOTSS, KINT, KDEG,
    # KD and KDEC as individual post hoc values from the PK fit). Values
    # are Qin 2026 Table 2; see Qin_2026_ropeginterferon for the full
    # source trace. All PD subjects are ET patients, so the patient
    # branch (Vc fold-change and Rtot decline) is hard-coded in model().
    # ==================================================================
    lka <- fixed(log(0.123)); label("First-order absorption rate constant (1/day; divided by 24 in model())") # Qin 2026 Table 2: ka 0.123 (RSE 3.82%)
    ltlag <- fixed(log(0.651)); label("Absorption lag time (h)") # Qin 2026 Table 2: Tlag 0.651 h
    lcl <- fixed(log(0.606)); label("Apparent linear clearance at BMI 24 kg/m^2 (L/day; divided by 24 in model())") # Qin 2026 Table 2: CL/F 0.606 L/day
    lvc <- fixed(log(1.88)); label("Apparent serum volume in healthy participants (L)") # Qin 2026 Table 2: Vc/F 1.88 L
    lrtot0 <- fixed(log(0.17)); label("Baseline maximum target binding capacity (ng/mL)") # Qin 2026 Table 2: Rtot,0 0.17 ng/mL
    lrtot_ss <- fixed(log(0.0142)); label("Steady-state maximum target binding capacity (SS, ng/mL)") # Qin 2026 Table 2: Rtot,SS 0.0142 ng/mL
    lkint <- fixed(log(0.12)); label("Elimination rate constant of the drug-target complex (1/h)") # Qin 2026 Table 2: kint 0.12 1/h
    lkdeg <- fixed(log(0.579)); label("Degradation rate constant of free target (1/h)") # Qin 2026 Table 2: kdeg 0.579 1/h
    lkd <- fixed(log(0.0455)); label("Equilibrium dissociation constant (ng/mL)") # Qin 2026 Table 2: KD 0.0455 ng/mL
    lkdecay <- fixed(log(0.0203)); label("Rate constant of the decline in binding capacity (1/day)") # Qin 2026 Table 2: kdec 0.0203 1/day
    e_bmi_cl <- fixed(1.03); label("Power exponent on (BMI / 24) for clearance (unitless)") # Qin 2026 Table 2: CL/F-BMI power 1.03
    e_et_vc <- fixed(log(2.007)); label("Log fold-change in vc for patients with essential thrombocythemia (unitless)") # Qin 2026 Table 2: Fold-increase in Vc/F for ET patient 2.007

    etalcl ~ fixed(0.097344) # Qin 2026 Table 2: omega (SD) for CL/F 0.312; 0.312^2
    etalvc ~ fixed(0.751689) # Qin 2026 Table 2: omega (SD) for Vc/F 0.867; 0.867^2
    etalka ~ fixed(0.1156) # Qin 2026 Table 2: omega (SD) for ka 0.34; 0.34^2
    etalrtot0 ~ fixed(1.8496) # Qin 2026 Table 2: omega (SD) for Rtot,0 1.36; 1.36^2
    etalkint ~ fixed(2.7889) # Qin 2026 Table 2: omega (SD) for kint 1.67; 1.67^2 (shared with Rtot,SS)

    # ==================================================================
    # WHITE BLOOD CELL LAYER -- Qin 2026 Table 3, "WBC" rows, with the
    # structure from the supplementary "NONMEM Model Code for the Final
    # WBC Model". The $THETA / $OMEGA there are initial estimates.
    #
    # TIME UNITS. Table 3 labels kout '1/h', but the control stream
    # labels THETA(1) '1/day' and divides it by 24
    # (TVKOUTH = THETA(1)/24), so the value is per day.
    # ==================================================================
    lkout <- log(0.07); label("First-order loss rate constant of circulating white blood cells (1/day; divided by 24 in model())") # Qin 2026 Table 3, WBC: kout 0.07 (RSE 1.92%), printed '1/h'; per day per the control stream
    logitep <- logit(0.147); label("Logit of EP, the empirical fractional reduction of white cell production from the start of treatment (logit)") # Qin 2026 Table 3, WBC: EPW 0.147 (RSE 1.65%); control stream TREP = LOG(TVEP/(1-TVEP))
    lrbase <- log(7.67); label("Typical initial white blood cell count at JAK2 V617F 13.7% in a non-Asian patient (10^9/L)") # Qin 2026 Table 3, WBC: WBC0 (10^9/L) 7.67 (RSE 1.65%)
    logitimax <- logit(0.493); label("Logit of the maximum fractional inhibition of white cell production by ropeg (logit)") # Qin 2026 Table 3, WBC: Imax,W 0.493 (RSE 1.76%); control stream TRIMAXH = LOG(TVIMAXH/(1-TVIMAXH))
    lic50 <- log(7.09); label("Total serum ropeg concentration giving half of Imax on white cell production (ng/mL)") # Qin 2026 Table 3, WBC: IC50,W (ng/mL) 7.09 (RSE 1.79%)
    lhill <- log(0.923); label("Hill coefficient of the Imax function (unitless)") # Qin 2026 Table 3, WBC: gamma 0.923 (RSE 2.1%)
    e_jak2_rbase <- 0.178; label("Power exponent on (JAK2 V617F / 13.7) for the initial white cell count (unitless)") # Qin 2026 Table 3, WBC: JAK2 on WBC0 0.178 (RSE 15.2%)
    e_asian_rbase <- log(1.34); label("Log fold-change in the initial white cell count for Asian patients (unitless)") # Qin 2026 Table 3, WBC: Asian on WBC0 1.34 (RSE 14%), a fold change; log(1.34) = 0.2927

    # Table 3 omega = SD and 'COV' = covariance: the control-stream
    # initial variances equal the squares of the Table 3 values
    # (kout 1.26^2 = 1.588 vs initial 1.566; WBC0 0.305^2 = 0.0930 vs
    # 0.0977). Banded BLOCK(5) in the order kout, EP, WBC0, Imax, IC50
    # with zeros off the first sub-diagonal, as in the control stream.
    etalkout + etalogitep + etalrbase + etalogitimax + etalic50 ~ c(
      1.5876,
      -0.604, 1.69,
      0, 0.163, 0.093025,
      0, 0, 0.0539, 0.692224,
      0, 0, 0, 0.376, 0.492804
    ) # Qin 2026 Table 3, WBC: omega kout 1.26, EPW 1.3, WBC0 0.305, HImax,W 0.832, IC50,W 0.702 (squared); COV kout-EPW -0.604, EPW-WBC0 0.163, WBC0-HImax,W 0.0539, IC50,W-HImax,W 0.376

    propSd_circ_wbc <- 0.162; label("Proportional residual SD on white blood cell count (fraction)") # Qin 2026 Table 3, WBC: sigma of proportional residual error 0.162 (RSE 0.552%)
    # Table 3 prints 0.0000146 with RSE 0; the control stream fixes the
    # variance at 2.12565E-012 (SD 1.458e-6). The fixed control-stream
    # value is used; either is negligible next to counts of ~5 x 10^9/L.
    addSd_circ_wbc <- fixed(1.458e-6); label("Additive residual SD on white blood cell count (10^9/L)") # control stream $SIGMA 2.12565E-012 FIX; sqrt = 1.458e-6 (Table 3 prints 0.0000146)
  })

  model({
    # ================= Inherited PK layer (fixed, patient branch) =====
    ka <- exp(lka + etalka) / 24
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) / 24 * (BMI / 24)^e_bmi_cl
    vc <- exp(lvc + etalvc) * exp(e_et_vc)
    rtot0 <- exp(lrtot0 + etalrtot0)
    rtot_ss <- exp(lrtot_ss + etalkint)
    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg)
    kd <- exp(lkd)
    kdecay <- exp(lkdecay)

    rtot_cap <- rtot_ss - (rtot_ss - rtot0) * exp(-kdecay * t / 24)
    ksyn <- kdeg * rtot_cap

    ctot <- central / vc
    disc <- ctot - total_target - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    complex <- ctot - cfree
    total_target(0) <- rtot0

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - kint * complex * vc
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * complex
    alag(depot) <- tlag

    Cc <- ctot

    # ================= White blood cell layer ========================
    # JAK2I = JAK2, replaced by 13.7 when <= 0 (control stream).
    jak2_ref <- TUM_JAK2_V617F_VAF
    if (TUM_JAK2_V617F_VAF <= 0) jak2_ref <- 13.7
    rbase <- exp(lrbase + etalrbase + e_jak2_rbase * log(jak2_ref / 13.7) + e_asian_rbase * RACE_ASIAN)
    kout <- exp(lkout + etalkout) / 24
    kin <- rbase * kout
    ep <- expit(logitep + etalogitep)
    imax <- expit(logitimax + etalogitimax)
    ic50 <- exp(lic50 + etalic50)
    hill <- exp(lhill)

    # INHH = IMAXH*CTOT**HILLH/(IC50H**HILLH + CTOT**HILLH), driven by
    # TOTAL serum ropeg. max(Cc, 0): past washout, integration noise can
    # take Cc slightly negative, where the non-integer power is NaN.
    Cc_pos <- max(Cc, 0)
    inh <- imax * Cc_pos^hill / (ic50^hill + Cc_pos^hill)

    # DADT(4) = KINH*(1 - INHH)*(1 - EP) - A(4)*KOUTH, A_0(4) = WBC0.
    # EP = EPH for TIME >= 0, i.e. from the first record, so the drug-
    # free steady state of a treated patient is WBC0*(1 - EP), not WBC0.
    circ_wbc(0) <- rbase
    d/dt(circ_wbc) <- kin * (1 - inh) * (1 - ep) - kout * circ_wbc

    circ_wbc ~ add(addSd_circ_wbc) + prop(propSd_circ_wbc)
  })
}
