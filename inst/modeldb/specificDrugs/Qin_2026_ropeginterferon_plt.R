Qin_2026_ropeginterferon_plt <- function() {
  description <- paste0(
    "Sequential PK-PD model of platelet count for subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 179 patients with essential ",
    "thrombocythemia (SURPASS-ET and EXCEED-ET). The fixed ",
    "quasi-equilibrium TMDD pharmacokinetic layer of ",
    "Qin_2026_ropeginterferon (patient branch) drives a sigmoid Imax ",
    "inhibition of platelet production by total serum ropeg in an ",
    "indirect-response model; from the start of treatment production is ",
    "further scaled by (1 - EP), an empirical fractional reduction. ",
    "No covariates on the platelet parameters."
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
    concentration = "ug/L total serum ropeg (Cc, = ng/mL); platelet count circ_plt in 10^9/L"
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
    circ_plt = list(analyte = "platelets", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Inherited pharmacokinetic covariate: CL = CL_TV * (BMI / 24)^1.03 (Qin 2026 Table 2). See Qin_2026_ropeginterferon. The platelet layer itself has no covariates.",
      source_name = "BMI"
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
    disease_state = "Essential thrombocythemia; baseline platelet count median 942 (SURPASS-ET) and 796 (EXCEED-ET) x 10^9/L, range 397-2140",
    dose_range = "HIDAT regimen: 250 ug subcutaneously at week 0, 350 ug at week 2, 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    notes = "Platelet counts measured every 2 weeks; counts were not followed after treatment discontinuation or interruption (Qin 2026 Discussion)."
  )

  ini({
    # ==================================================================
    # PHARMACOKINETIC LAYER -- ENTIRELY FIXED (sequential PK-PD; the
    # control stream reads the PK parameters as individual post hoc
    # values). Values are Qin 2026 Table 2; see Qin_2026_ropeginterferon
    # for the full source trace. All PD subjects are ET patients, so the
    # patient branch is hard-coded in model().
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
    # PLATELET LAYER -- Qin 2026 Table 3, "PLT" rows, with the structure
    # from the supplementary "NONMEM Model Code for the Final platelet
    # Model". The $THETA / $OMEGA there are initial estimates.
    #
    # TIME UNITS. Table 3 labels kout '1/h', but the control stream
    # labels THETA(1) '1/day' and divides it by 24, so the value is per
    # day.
    # ==================================================================
    lkout <- log(0.033); label("First-order loss rate constant of circulating platelets (1/day; divided by 24 in model())") # Qin 2026 Table 3, PLT: kout 0.033 (RSE 1.73%), printed '1/h'; per day per the control stream
    lrbase <- log(884); label("Typical initial platelet count (10^9/L)") # Qin 2026 Table 3, PLT: PLT0 (10^9/L) 884 (RSE 1.67%)
    logitimax <- logit(0.598); label("Logit of the maximum fractional inhibition of platelet production by ropeg (logit)") # Qin 2026 Table 3, PLT: Imax,P 0.598 (RSE 1.26%)
    lic50 <- log(7.02); label("Total serum ropeg concentration giving half of Imax on platelet production (ng/mL)") # Qin 2026 Table 3, PLT: IC50,P (ng/mL) 7.02 (RSE 1.86%)
    logitep <- logit(0.121); label("Logit of EP, the empirical fractional reduction of platelet production from the start of treatment (logit)") # Qin 2026 Table 3, PLT: EPP 0.121 (RSE 1.86%)
    lhill <- log(1.86); label("Hill coefficient of the Imax function (unitless)") # Qin 2026 Table 3, PLT: gamma 1.86 (RSE 2.02%)

    # Table 3 omega = SD and 'COV' = covariance (control-stream initial
    # variances: kout 0.391 vs 0.632^2 = 0.399; PLT0 0.132 vs
    # 0.353^2 = 0.125; covariance initial 0.112 vs final 0.0998).
    etalkout ~ 0.399424 # Qin 2026 Table 3, PLT: omega for kout 0.632 (RSE 6.17%, shrinkage 10.1%); 0.632^2
    etalrbase + etalogitimax ~ c(
      0.124609,
      0.0998, 0.346921
    ) # Qin 2026 Table 3, PLT: omega PLT0 0.353 and HImax,P 0.589 (squared); COV PLT0-HImax,P 0.0998
    etalic50 ~ 1.0404 # Qin 2026 Table 3, PLT: omega for IC50,P 1.02 (RSE 9.3%, shrinkage 24%); 1.02^2
    etalogitep ~ 2.25 # Qin 2026 Table 3, PLT: omega for EPP 1.5 (RSE 11.5%, shrinkage 31.9%); 1.5^2

    propSd_circ_plt <- 0.116; label("Proportional residual SD on platelet count (fraction)") # Qin 2026 Table 3, PLT: sigma of proportional residual error 0.116 (RSE 0.541%)
    addSd_circ_plt <- 0.103; label("Additive residual SD on platelet count (10^9/L)") # Qin 2026 Table 3, PLT: sigma of additional residual error (10^9/L) 0.103 (RSE 3.67%)
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

    # ================= Platelet layer ================================
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout) / 24
    kin <- rbase * kout
    imax <- expit(logitimax + etalogitimax)
    ic50 <- exp(lic50 + etalic50)
    ep <- expit(logitep + etalogitep)
    hill <- exp(lhill)

    inh <- imax * Cc^hill / (ic50^hill + Cc^hill)

    # DADT(4) = KINH*(1 - INHH)*(1 - EP) - A(4)*KOUTH, A_0(4) = PLT0;
    # EP applies from TIME >= 0.
    circ_plt(0) <- rbase
    d/dt(circ_plt) <- kin * (1 - inh) * (1 - ep) - kout * circ_plt

    circ_plt ~ add(addSd_circ_plt) + prop(propSd_circ_plt)
  })
}
