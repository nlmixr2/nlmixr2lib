Qin_2026_ropeginterferon <- function() {
  description <- paste0(
    "Quasi-equilibrium target-mediated drug disposition (TMDD) ",
    "population pharmacokinetic model for subcutaneous ",
    "ropeginterferon alfa-2b (ropeg, a mono-PEGylated interferon ",
    "alfa-2b) in 273 pooled participants: 94 healthy volunteers from ",
    "three single-dose phase I studies (A09-102, A17-101, A17-102; ",
    "24-300 ug) and 179 patients with essential thrombocythemia from ",
    "SURPASS-ET and EXCEED-ET (250 ug, 350 ug at week 2, then 500 ug ",
    "every 2 weeks). First-order absorption with a lag into a single ",
    "serum compartment carrying linear clearance of free drug plus ",
    "saturable binding to a turnover target pool; the drug-target ",
    "complex is internalised at kint. Body mass index acts on ",
    "clearance as a power function and the serum volume is 2.007-fold ",
    "higher in patients with essential thrombocythemia, in whom the ",
    "target binding capacity also declines from Rtot0 to Rtot,SS ",
    "under chronic dosing. Companion platelet, white-cell and ",
    "exposure-response models in the Qin_2026_ropeginterferon_* family."
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
    concentration = "ug/L total serum ropeg (Cc), numerically identical to the ng/mL Qin 2026 reports. Doses in ug with vc in L give ug/L directly"
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
    )
  )

  covariateData <- list(
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate on clearance in the final model (Qin 2026",
        "Results, PopPK Model; Table 2 'CL/F (L/day) at median BMI of",
        "24 kg/m2'). Enters as the Equation (1) power model",
        "CL_i = CL_TV * (BMI / 24)^theta. The published control stream",
        "centres at exactly 24 (CLBMI = THETA(12)*LOG(BMII/24)); the",
        "Table 1 pooled median is 24.3. The control stream imputes a",
        "missing BMI (coded negative) as 24, i.e. the reference value;",
        "18 EXCEED-ET patients (20.5%) had no BMI."
      ),
      source_name = "BMI"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator; 1 = healthy volunteer, 0 = patient with essential thrombocythemia.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (healthy volunteer), at which vc takes its Table 2 value of 1.88 L and the binding capacity is held at Rtot0",
      notes = paste(
        "The control stream carries ET = 1 for patients with essential",
        "thrombocythemia, so DIS_HEALTHY = 1 - ET. It acts twice:",
        "(1) Vc_i = Vc_TV * exp(theta_ET * ET), Table 2 'Fold-increase",
        "in Vc/F for ET patient' 2.007; (2) the decline of the target",
        "binding capacity from Rtot0 to Rtot,SS runs only when ET = 1",
        "($DES: IF(ET.EQ.1.AND.DAY.GT.TSTART) RTOT = ...). All 94",
        "phase I participants were healthy and all 179 SURPASS-ET and",
        "EXCEED-ET participants were patients (Table 1)."
      ),
      source_name = "ET (inverted)"
    )
  )

  covariatesDataExcluded <- list(
    ALT = list(
      description = "Baseline alanine aminotransferase activity.",
      units = "U/L",
      type = "continuous",
      notes = "Present in the final control stream as CLALT = THETA(13)*LOG(ALT/19) with THETA(13) = 0 FIX, i.e. switched off. Not part of the final model (Table 2)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Present in the final control stream as VCAGE = THETA(14)*LOG(AGE/45) with THETA(14) = 0 FIX. Not part of the final model."
    ),
    PLT = list(
      description = "Baseline platelet count.",
      units = "10^9/L",
      type = "continuous",
      notes = "Present in the final control stream as VCPLT = THETA(16)*LOG(PLAT/638) with THETA(16) = 0 FIX. Not part of the final model."
    ),
    HGB = list(
      description = "Baseline haemoglobin concentration.",
      units = "g/dL",
      type = "continuous",
      notes = "Present in the final control stream as VCHGB = THETA(17)*LOG(HGB/14.2) with THETA(17) = 0 FIX. Not part of the final model."
    ),
    SEXF = list(
      description = "Female-sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Present in the final control stream as VCSEX = THETA(18)*SEXF with THETA(18) = 0 FIX. Not part of the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 273L,
    n_studies = 5L,
    age_range = "median 45.0 years, range 19.0-84.0 pooled (Qin 2026 Table 1, Overall); healthy phase I medians 27.0-31.5 years, ET patient medians 57.0 (EXCEED-ET) and 61.0 (SURPASS-ET) years",
    weight_range = "median 70.0 kg, range 41.5-153 (Qin 2026 Table 1, Overall)",
    bmi_range = "median 24.3 kg/m^2, range 16.9-49.7; 18 EXCEED-ET patients missing (Qin 2026 Table 1, Overall)",
    sex_female_pct = 42.1,
    race_ethnicity = "Asian 138 (50.5%), Caucasian 114 (41.8%), other 21 (7.7%) (Qin 2026 Table 1, Overall). SURPASS-ET was 94.5% Asian; EXCEED-ET (United States and Canada) 77.3% Caucasian",
    disease_state = "94 healthy volunteers (A09-102 n = 32 Caucasian; A17-101 n = 27 Chinese; A17-102 n = 35 Japanese and Caucasian) and 179 patients with essential thrombocythemia (SURPASS-ET n = 91, phase 3 second-line versus anagrelide; EXCEED-ET n = 88, phase 2 single arm). Baseline platelet count median 942 (SURPASS-ET) and 796 (EXCEED-ET) x 10^9/L",
    dose_range = "phase I single subcutaneous doses of 24-270 ug (A09-102), 90-270 ug (A17-101) and 100-300 ug (A17-102); in ET patients the HIDAT regimen, 250 ug at week 0, 350 ug at week 2 and 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "Europe (A09-102), China (A17-101), Japan (A17-102), China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    anti_drug_antibody = "30 of 273 (11.0%) ADA-positive (Qin 2026 Table 1); not a covariate in the final model",
    notes = paste0(
      "NONMEM 7.5.0 FOCE-I (Qin 2026 Software and Data Analysis); the ",
      "published final control stream (supplement) is the structural ",
      "source for this file. All fixed and random effects had RSE below ",
      "15% except Rtot,0 (15.9%), kint (19.8%) and KD (32.2%) ",
      "(Qin 2026 Table 2)."
    )
  )

  ini({
    # ==================================================================
    # Qin 2026 Table 2, "Parameter Estimates of Final PopPK Model", with
    # the structure taken from the supplementary final control stream
    # ("NONMEM Model Code for the Final Population Pharmacokinetic
    # Model"). The $THETA / $OMEGA / $SIGMA blocks of that stream hold
    # INITIAL estimates (e.g. CL 0.611602 vs the final 0.606), so every
    # value below is the Table 2 final estimate.
    #
    # TIME UNIT. NONMEM TIME is in hours: the control stream divides the
    # day-based thetas by 24 (MU_1 = LOG(THETA(1)/24) for CL in L/day,
    # MU_3 = LOG(THETA(3)/24) for ka in 1/day), uses kint and kdeg
    # unscaled with their '1/hr' labels, and reads the lag in 'hr'.
    # Supplementary Table S1 reports AUC in hr*ng/mL with
    # AUC/Cavg = 336 h = one 2-week interval. The model time unit is
    # therefore hours and the /24 conversions are applied in model(),
    # exactly as in the control stream, so the values below stay in the
    # Table 2 units.
    # ==================================================================

    # ----- Absorption -----
    lka <- log(0.123); label("First-order absorption rate constant from the subcutaneous depot (1/day; divided by 24 in model())") # Qin 2026 Table 2: ka (1/day) 0.123 (RSE 3.82%)
    ltlag <- log(0.651); label("Absorption lag time on the subcutaneous depot (h)") # Qin 2026 Table 2: Tlag (h) 0.651 (RSE 2.09%)

    # ----- Linear disposition -----
    lcl <- log(0.606); label("Apparent linear clearance of free drug at BMI 24 kg/m^2 (L/day; divided by 24 in model())") # Qin 2026 Table 2: CL/F (L/day) at median BMI of 24 kg/m2 0.606 (RSE 2.68%)
    lvc <- log(1.88); label("Apparent volume of the serum compartment in healthy participants (L)") # Qin 2026 Table 2: Vc/F (L) 1.88 (RSE 8.37%)

    # ----- Target-mediated disposition (quasi-equilibrium) -----
    lrtot0 <- log(0.17); label("Baseline maximum target binding capacity, the initial condition of the total-target pool (ng/mL)") # Qin 2026 Table 2: Rtot,0 (ng/mL) 0.17 (RSE 15.9%)
    lrtot_ss <- log(0.0142); label("Steady-state maximum target binding capacity reached under chronic dosing in patients (SS, ng/mL)") # Qin 2026 Table 2: Rtot,SS (ng/mL) 0.0142 (RSE 13.3%)
    lkint <- log(0.12); label("First-order elimination rate constant of the drug-target complex (1/h)") # Qin 2026 Table 2: kint (1/h) 0.12 (RSE 19.8%)
    lkdeg <- log(0.579); label("First-order degradation rate constant of free target (1/h)") # Qin 2026 Table 2: kdeg (1/h) 0.579 (RSE 7.63%)
    lkd <- log(0.0455); label("Equilibrium dissociation constant of ropeg for its target (ng/mL)") # Qin 2026 Table 2: KD (ng/mL) 0.0455 (RSE 32.2%)
    lkdecay <- log(0.0203); label("First-order rate constant of the decline in binding capacity from Rtot0 to Rtot,SS (1/day; applied to t/24 in model())") # Qin 2026 Table 2: kdec (1/day) 0.0203 (RSE 2.7%)

    # ----- Covariate effects -----
    e_bmi_cl <- 1.03; label("Power exponent on (BMI / 24) for clearance (unitless)") # Qin 2026 Table 2: CL/F-BMI power 1.03 (RSE 14.8%)
    # Table 2 prints the FOLD change; the control stream adds
    # THETA(19)*ET inside the exponent, so theta = log(2.007).
    e_et_vc <- log(2.007); label("Log fold-change in vc for patients with essential thrombocythemia versus healthy participants (unitless)") # Qin 2026 Table 2: Fold-increase in Vc/F for ET patient 2.007 (RSE 11.9%); log(2.007) = 0.6966

    # ----- Inter-individual variability -----
    # Table 2's omega column is the STANDARD DEVIATION of the log-normal
    # eta (variance = omega^2). The control-stream initial variances
    # match the squares of the Table 2 values (ka: 0.34^2 = 0.1156 vs
    # initial 0.117492; CL: 0.312^2 = 0.0973 vs 0.122131), not the
    # values themselves (ka 0.34 vs 0.117).
    etalcl ~ 0.097344 # Qin 2026 Table 2: omega for CL/F 0.312 (RSE 3.87%); 0.312^2
    etalvc ~ 0.751689 # Qin 2026 Table 2: omega for Vc/F 0.867 (RSE 5.79%); 0.867^2
    etalka ~ 0.1156 # Qin 2026 Table 2: omega for ka 0.34 (RSE 7.25%, shrinkage 15.1%); 0.34^2
    etalrtot0 ~ 1.8496 # Qin 2026 Table 2: omega for Rtot,0 1.36 (RSE 7.28%, shrinkage 10.3%); 1.36^2
    etalkint ~ 2.7889 # Qin 2026 Table 2: omega for kint 1.67 (RSE 7.26%, shrinkage 33.7%); 1.67^2. The control stream applies this SAME eta to Rtot,SS (RTOTSS = THETA(6)*EXP(ETA(5)))

    # ----- Residual unexplained variability -----
    # Control stream: Y = F + F*EPS(1) + EPS(2); Table 2 gives sigma as
    # standard deviations (initial variance 0.0443364 = 0.2106^2 vs
    # final 0.209).
    propSd <- 0.209; label("Proportional residual SD on total serum ropeg concentration (fraction)") # Qin 2026 Table 2: sigma of proportional residual error 0.209 (RSE 0.637%, shrinkage 20.6%)
    addSd <- 0.155; label("Additive residual SD on total serum ropeg concentration (ng/mL)") # Qin 2026 Table 2: sigma of additional residual error (ng/mL) 0.155 (RSE 5.18%, shrinkage 23%)
  })

  model({
    # ----- Individual parameters (control stream $PK) -----
    et <- 1 - DIS_HEALTHY
    ka <- exp(lka + etalka) / 24
    tlag <- exp(ltlag)
    # Qin 2026 Equation (1): CL_i = CL_TV * (BMI / 24)^theta
    cl <- exp(lcl + etalcl) / 24 * (BMI / 24)^e_bmi_cl
    vc <- exp(lvc + etalvc + e_et_vc * et)
    rtot0 <- exp(lrtot0 + etalrtot0)
    rtot_ss <- exp(lrtot_ss + etalkint)
    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg)
    kd <- exp(lkd)
    kdecay <- exp(lkdecay)

    # ----- Time-varying target binding capacity (control stream $DES) -----
    # RTOT = RTOTSS - (RTOTSS - RTOT0)*EXP(-KDEC*(DAY - TSTART)) for
    # ET = 1 and DAY > TSTART (TSTART = 0 FIX), else RTOT0. kdec is per
    # day, so it is applied to t/24.
    rtot_cap <- rtot0 + et * (rtot_ss - rtot0) * (1 - exp(-kdecay * t / 24))
    ksyn <- kdeg * rtot_cap

    # ----- Quasi-equilibrium binding (control stream $DES) -----
    # CPUB = 0.5*((CTOT - A(3) - KD) + sqrt((CTOT - A(3) - KD)^2 + 4*KD*CTOT))
    ctot <- central / vc
    disc <- ctot - total_target - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    complex <- ctot - cfree

    total_target(0) <- rtot0

    # ----- ODEs -----
    # DADT(2) = KA*A(1) - K20*FREE - KINT*(A(2) - FREE), FREE = VC*CPUB
    # DADT(3) = KSYN - KDEG*A(3) - (KINT - KDEG)*(CTOT - CPUB)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - kint * complex * vc
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * complex
    alag(depot) <- tlag

    # ----- Observation -----
    # IPRED = CTOT2 = A(2)/S2 with S2 = VC: total serum ropeg, the
    # concentration the PK-PD models name CS.
    Cc <- ctot
    Cc ~ add(addSd) + prop(propSd)
  })
}
