Kathman_2022_bentracimab <- function() {
  description <- "Semi-mechanistic population PK-PD model of ticagrelor reversal by bentracimab (PB2452, a ticagrelor-neutralising monoclonal antibody Fab fragment) in healthy volunteers pretreated with oral ticagrelor (Kathman 2022). Ticagrelor: two transit compartments into a two-compartment disposition; its active metabolite (TAM, AR-C124910XX): two-compartment disposition fed by a metabolic flux whose fraction rises with the uncomplexed PB2452 concentration; PB2452: two-compartment disposition. PB2452 binds ticagrelor and TAM (second-order association) to form complexes cleared with PB2452; complexes also pass at a weight-dependent rate ktr into delayed complex pools from which ticagrelor and TAM return to plasma. Platelet reactivity (VerifyNow PRU) falls from each subject's observed pre-ticagrelor baseline through two additive sigmoid Emax terms on uncomplexed ticagrelor and TAM."
  reference <- "Kathman SJ, Wheeler JJ, Bhatt DL, Arnold SE, Lee JS. Population pharmacokinetic-pharmacodynamic modeling of PB2452, a monoclonal antibody fragment being developed as a ticagrelor reversal agent, in healthy volunteers. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):68-81. doi:10.1002/psp4.12734"
  vignette <- "Kathman_2022_bentracimab"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as a log-linear effect centred at log(WT) = 4.35 (WT = exp(4.35) = 77.5 kg) on ktr, on the TAM EC50 for PRU and on the TAM clearance, i.e. a power model (WT / 77.5)^theta (Kathman 2022 Table 3 and Supplementary Material S1 MU_4, MU_8, MU_10). The cohort weight distribution is not reported.",
      source_name = "WT"
    ),
    BL_PRU = list(
      description = "Subject-specific observed baseline platelet reactivity in P2Y12 reaction units, measured before the first ticagrelor dose",
      units = "PRU",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed. The PRU baseline is not estimated: 'Base = Log(BPRU) baseline PRU, fixed in model to observed baseline values' with a 10% fixed IIV to allow for measurement error (Kathman 2022 Table 3; Supplementary Material S1 'Base=EXP(LOG(BPRU) + ETA(14))'). The cohort distribution is not reported.",
      source_name = "BPRU"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (ETA-versus-covariate plots, Supplementary Material S4/S5) but not retained (Kathman 2022 Results 'Examining covariates')."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Kathman 2022 Methods 'Data analysis'; source column SEXN)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (source column BAST)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (source column BALT)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained (source column BALKP)."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened but not retained (source column BHCT)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (abbreviated MDRD)",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Screened but not retained; cohort range 74.75-162.80 mL/min/1.73 m^2 (Kathman 2022 Discussion; source column EGFR)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ticagrelor", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ticagrelor", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "ticagrelor", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ticagrelor", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ticagrelor", units = "nmol", specimen = "plasma", verified = TRUE),
    central_tam = list(analyte = "TAM (AR-C124910XX)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_tam = list(analyte = "TAM (AR-C124910XX)", units = "nmol", specimen = "plasma", verified = TRUE),
    target = list(analyte = "bentracimab", units = "nmol", specimen = "plasma", verified = TRUE),
    complex = list(analyte = "bentracimab-ticagrelor complex", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_tam = list(analyte = "bentracimab-TAM complex", units = "nmol", specimen = "plasma", verified = TRUE),
    target_peripheral1 = list(analyte = "bentracimab", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_peripheral1 = list(
      analyte = "bentracimab-ticagrelor complex",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_peripheral1_tam = list(
      analyte = "bentracimab-TAM complex",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61,
    n_studies = 1,
    age_range = "18-50 years",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = "Healthy volunteers; cohorts 4-10 pretreated with oral ticagrelor to steady state",
    dose_range = "Ticagrelor 180 mg oral loading dose then 90 mg twice daily (5 doses over 48 h). PB2452 0.1-18 g IV: single 30-minute infusions of 0.1-9 g (cohorts 1-6) or 18 g as a bolus plus one or two prolonged infusions (cohorts 7-10; e.g. 6 g over 15 min + 6 g over 4 h + 6 g over 12 h).",
    regions = "United States (single centre)",
    notes = "Phase I randomised, double-blind, placebo-controlled single-ascending-dose trial (10 cohorts; Kathman 2022 Table 1): 48 subjects received PB2452 and 13 placebo. PB2452 started immediately after the 5th ticagrelor dose (48 h; cohorts 4-6) or 2 h after it (50 h; cohorts 7-10). Time zero is the first ticagrelor dose. PB2452, ticagrelor and TAM concentrations were converted to nmol/L before modelling. Bayesian MCMC estimation in NONMEM 7.4 (5000 burn-in, 20000 samples per chain). eGFR 74.75-162.80 mL/min/1.73 m^2."
  )

  ini({
    # Ticagrelor (TICA) PK -- fixed from the published ticagrelor PK-PD model
    # (Kathman 2022 reference 14, Astrand et al.) with small fixed IIVs.
    lka <- fixed(2.3); label("TICA absorption rate constant, also the depot -> transit1 -> transit2 -> central transfer rate (1/h)") # Kathman 2022 Table 3, KA = EXP(2.3), 'Fixed in model'; Supplementary Material S1 DADT(1)-DADT(4)
    lcl <- fixed(2.81); label("TICA apparent clearance CL/F (L/h)") # Kathman 2022 Table 3, CL/F = EXP(2.81), 'Fixed in model'
    lvc <- fixed(5.04); label("TICA apparent central volume V1/F (L)") # Kathman 2022 Table 3, V1/F = EXP(5.04), 'Fixed in model'
    lvp <- fixed(4.02); label("TICA apparent peripheral volume V2/F (L)") # Kathman 2022 Table 3, V2/F = EXP(4.02), 'Fixed in model'
    lq <- fixed(2.34); label("TICA apparent intercompartmental clearance Q1/F (L/h)") # Kathman 2022 Table 3, Q1/F = EXP(2.34), 'Fixed in model'
    fm <- fixed(0.3); label("TICA -> TAM metabolic rate as a fraction of the TICA elimination rate constant CL/V1 in the absence of PB2452 (unitless)") # Kathman 2022 Table 3, 'fm = 0.3', 'Fixed in model'

    # Ticagrelor active metabolite (TAM)
    lcl_tam <- 1.93; label("TAM clearance at WT = 77.5 kg (L/h)") # Kathman 2022 Table 3, CLM = EXP(THETA11 + THETA10*(LOG(WT)-4.35)), THETA11 = 1.93
    e_wt_cl_tam <- 1.31; label("Power exponent of body weight on TAM clearance (unitless)") # Kathman 2022 Table 3, THETA10 = 1.31
    lvc_tam <- fixed(1.95); label("TAM central volume VM1 (L)") # Kathman 2022 Table 3, VM1 = EXP(1.95), 'Fixed in model'
    lvp_tam <- fixed(3.74); label("TAM peripheral volume VM2 (L)") # Kathman 2022 Table 3, VM2 = EXP(3.74), 'Fixed in model'
    lq_tam <- fixed(1.48); label("TAM intercompartmental clearance Q2M (L/h)") # Kathman 2022 Table 3, Q2M = EXP(1.48), 'Fixed in model'

    # PB2452 (bentracimab) PK -- fixed from the PB2452-alone fit of cohorts 1-3
    # (Kathman_2022_bentracimab_pk). The values printed in Table 3 and used in
    # the control stream (Q_ant -0.765, V_ant_perp 1.28) differ slightly from
    # the Table 2 estimates (-0.770, 1.24); the as-run Table 3 values are used.
    lcl_target <- fixed(0.631); label("PB2452 clearance CL_ant, shared by uncomplexed PB2452 and both PB2452 complexes (L/h)") # Kathman 2022 Table 3, CL_ant = EXP(0.631), 'Fixed in model'
    lvc_target <- fixed(1.05); label("PB2452 central volume V_ant (L)") # Kathman 2022 Table 3, V_ant = EXP(1.05), 'Fixed in model'
    lqp <- fixed(-0.765); label("PB2452 intercompartmental clearance Q_ant (L/h)") # Kathman 2022 Table 3, Q_ant = EXP(-0.765), 'Fixed in model'
    lvp_target <- fixed(1.28); label("PB2452 peripheral volume V_ant_perp (L)") # Kathman 2022 Table 3, V_ant_perp = EXP(1.28), 'Fixed in model'

    # PB2452 binding to TICA and TAM
    lk1 <- -5.56; label("Second-order association rate constant of PB2452 with TICA, Kon (1/(nmol/L)/h)") # Kathman 2022 Table 3, Kon = EXP(THETA2), THETA2 = -5.56
    lk1_tam <- -3.74; label("Second-order association rate constant of PB2452 with TAM, Kon2 (1/(nmol/L)/h)") # Kathman 2022 Table 3, Kon2 = EXP(THETA5), THETA5 = -3.74
    lkd <- fixed(-4); label("Equilibrium dissociation constant Kd setting the plasma-complex off-rate Koff = Kon * Kd (nmol/L)") # Kathman 2022 Table 3, Kd = EXP(-4), 'Fixed in model' (from the preclinical model, reference 11)
    lkd2 <- 2.04; label("Dissociation constant Kd2 setting the release rate from the delayed complex pools Koff2 = Kon * Kd2 (nmol/L)") # Kathman 2022 Table 3, Kd2 = EXP(THETA3), THETA3 = 2.04
    lktr <- -1.22; label("Transfer rate constant of the plasma complexes into the delayed complex pools at WT = 77.5 kg, Ktr (1/h)") # Kathman 2022 Table 3, Ktr = EXP(THETA4 + THETA13*(LOG(WT)-4.35)), THETA4 = -1.22
    e_wt_ktr <- 1.46; label("Power exponent of body weight on Ktr (unitless)") # Kathman 2022 Table 3, THETA13 = 1.46

    # PB2452-stimulated TICA -> TAM metabolism
    lemax_fm <- 2.98; label("Maximal fractional increase of fm driven by uncomplexed PB2452, Emaxf (unitless)") # Kathman 2022 Table 3, Emaxf = EXP(THETA6), THETA6 = 2.98
    lec50_fm <- 9.36; label("Uncomplexed PB2452 concentration at half-maximal fm stimulation, ECf (nmol/L)") # Kathman 2022 Table 3, ECf = EXP(THETA7), THETA7 = 9.36
    hill_fm <- fixed(2); label("Hill coefficient of the PB2452 effect on fm (unitless)") # Kathman 2022 Results: Hill coefficients 'were all set equal to 2'; Supplementary Material S1 '**2'

    # PRU response to uncomplexed TICA and TAM
    lec50 <- 10.6; label("Uncomplexed TICA concentration at half-maximal PRU inhibition, EC50 (nmol/L)") # Kathman 2022 Table 3, EC50 = EXP(THETA1), THETA1 = 10.6
    lemax <- fixed(-0.1); label("Maximal fractional PRU inhibition by TICA, Emax (unitless)") # Kathman 2022 Table 3, Emax = EXP(-0.1) (90%), 'Fixed in model'
    lec50_tam <- 4.59; label("Uncomplexed TAM concentration at half-maximal PRU inhibition at WT = 77.5 kg, EC502 (nmol/L)") # Kathman 2022 Table 3, EC502 = EXP(THETA8 + THETA12*(LOG(WT)-4.35)), THETA8 = 4.59 (98.5 nmol/L, Discussion)
    e_wt_ec50_tam <- -0.965; label("Power exponent of body weight on EC502 (unitless)") # Kathman 2022 Table 3, THETA12 = -0.965
    lemax_tam <- 0.0181; label("Maximal fractional PRU inhibition by TAM, Emax2 (unitless)") # Kathman 2022 Table 3, Emax2 = EXP(THETA9), THETA9 = 0.0181
    hill <- fixed(2); label("Hill coefficient of the TICA and TAM effects on PRU (unitless)") # Kathman 2022 Results: Hill coefficients 'were all set equal to 2'; Supplementary Material S1 '**2'

    # IIV: exponential; variance = CV^2 (Table 3 prints OMEGA 0.01 as '10%' and
    # 0.0025 as '5%'). The OMEGA BLOCK(10) covariances were not reported, so the
    # estimated etas are encoded as uncorrelated.
    etalec50 ~ 0.054289 # Kathman 2022 Table 3, EC50 IIV 23.3% -> 0.233^2 (ETA1)
    etalk1 ~ 0.190096 # Kathman 2022 Table 3, Kon IIV 43.6% -> 0.436^2 (ETA2)
    etalkd2 ~ 0.067081 # Kathman 2022 Table 3, Kd2 IIV 25.9% -> 0.259^2 (ETA3)
    etalktr ~ 0.0625 # Kathman 2022 Table 3, Ktr IIV 25.0% -> 0.250^2 (ETA4)
    etalk1_tam ~ 0.092416 # Kathman 2022 Table 3, Kon2 IIV 30.4% -> 0.304^2 (ETA5)
    etalemax_fm ~ 0.055696 # Kathman 2022 Table 3, Emaxf IIV 23.6% -> 0.236^2 (ETA6)
    etalec50_fm ~ 0.357604 # Kathman 2022 Table 3, ECf IIV 59.8% -> 0.598^2 (ETA7)
    etalec50_tam ~ 0.064009 # Kathman 2022 Table 3, EC502 IIV 25.3% -> 0.253^2 (ETA8)
    etalemax_tam ~ 0.042849 # Kathman 2022 Table 3, Emax2 IIV 20.7% -> 0.207^2 (ETA9)
    etalcl_tam ~ 0.057121 # Kathman 2022 Table 3, CLM IIV 23.9% -> 0.239^2 (ETA10)
    etalka ~ fixed(0.01) # Kathman 2022 Table 3, KA IIV 10%; Supplementary Material S1 OMEGA(11,11) = 0.01
    etalcl ~ fixed(0.01) # Kathman 2022 Table 3, CL/F IIV 10%; Supplementary Material S1 OMEGA(12,12) = 0.01
    etalq ~ fixed(0.01) # Kathman 2022 Table 3, Q1/F IIV 10%; Supplementary Material S1 OMEGA(13,13) = 0.01
    etalrbase ~ fixed(0.01) # Kathman 2022 Table 3, Base IIV 10%; Supplementary Material S1 OMEGA(14,14) = 0.01
    etalvp_tam ~ fixed(0.01) # Kathman 2022 Table 3, VM2 IIV 10%; Supplementary Material S1 OMEGA(15,15) = 0.01
    etalcl_target ~ fixed(0.0025) # Kathman 2022 Table 3, CL_ant IIV 5%; Supplementary Material S1 OMEGA(16,16) = 0.0025
    etalqp ~ fixed(0.0025) # Kathman 2022 Table 3, Q_ant IIV 5%; Supplementary Material S1 OMEGA(17,17) = 0.0025
    etalvc_target ~ fixed(0.0025) # Kathman 2022 Table 3, V_ant IIV 5%; Supplementary Material S1 OMEGA(18,18) = 0.0025
    etalvp_target ~ fixed(0.0025) # Kathman 2022 Table 3, V_ant_perp IIV 5%; Supplementary Material S1 OMEGA(19,19) = 0.0025

    # Residual error (Table 3 'Residual Variability')
    propSd <- 0.424; label("Proportional residual error, total TICA (fraction)") # Kathman 2022 Table 3, 'Total TICA: 42.4%'
    addSd <- 332; label("Additive residual error, total TICA (nmol/L)") # Kathman 2022 Table 3, 'Total TICA: Additive SD = 332'
    propSd_tam <- 0.231; label("Proportional residual error, total TAM (fraction)") # Kathman 2022 Table 3, 'Total TAM: 23.1%'
    propSd_Cc_target <- 0.282; label("Proportional residual error, uncomplexed PB2452 (fraction)") # Kathman 2022 Table 3, 'Uncomplexed PB2452: 28.2%'
    propSd_Ctotal_target <- 0.137; label("Proportional residual error, total PB2452 (fraction)") # Kathman 2022 Table 3, 'Total PB2452: 13.7%'
    addSd_Ctotal_target <- 504; label("Additive residual error, total PB2452 (nmol/L)") # Kathman 2022 Table 3, 'Total PB2452: Additive SD = 504'
    addSd_PRU <- 20.3; label("Additive residual error, PRU (PRU)") # Kathman 2022 Table 3, 'PRU: Additive SD = 20.3'
  })

  model({
    # Ticagrelor
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)

    # TAM; weight enters as (log(WT) - 4.35), i.e. (WT / 77.5)^theta
    cl_tam <- exp(lcl_tam + e_wt_cl_tam * (log(WT) - 4.35) + etalcl_tam)
    vc_tam <- exp(lvc_tam)
    vp_tam <- exp(lvp_tam + etalvp_tam)
    q_tam <- exp(lq_tam)

    # PB2452
    cl_target <- exp(lcl_target + etalcl_target)
    vc_target <- exp(lvc_target + etalvc_target)
    qp <- exp(lqp + etalqp)
    vp_target <- exp(lvp_target + etalvp_target)

    # Binding
    k1 <- exp(lk1 + etalk1)
    k1_tam <- exp(lk1_tam + etalk1_tam)
    kd <- exp(lkd)
    kd2 <- exp(lkd2 + etalkd2)
    ktr <- exp(lktr + e_wt_ktr * (log(WT) - 4.35) + etalktr)
    # Both off-rates use the TICA association constant Kon, including the TAM
    # complexes (Supplementary Material S1: koff = kon*Kd, koff2 = kon*Kd2).
    k2 <- k1 * kd
    k2_rec <- k1 * kd2

    # PD
    ec50 <- exp(lec50 + etalec50)
    emax <- exp(lemax)
    ec50_tam <- exp(lec50_tam + e_wt_ec50_tam * (log(WT) - 4.35) + etalec50_tam)
    emax_tam <- exp(lemax_tam + etalemax_tam)
    emax_fm <- exp(lemax_fm + etalemax_fm)
    ec50_fm <- exp(lec50_fm + etalec50_fm)
    rbase <- exp(log(BL_PRU) + etalrbase)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_tam <- cl_tam / vc_tam
    k12_tam <- q_tam / vc_tam
    k21_tam <- q_tam / vp_tam

    # Uncomplexed concentrations (nmol/L)
    ctica <- central / vc
    ctam <- central_tam / vc_tam
    ctarget <- target / vc_target

    # Fraction metabolised rises with uncomplexed PB2452 (Kathman 2022 Results
    # 'Final model structure'); metabolism is an additional TICA loss on top of
    # CL/F, as coded in Supplementary Material S1 DADT(4).
    frac_met <- fm * (1 + emax_fm * ctarget^hill_fm / (ec50_fm^hill_fm + ctarget^hill_fm))
    kmet <- frac_met * kel

    # Binding fluxes (nmol/h), Supplementary Material S1: kon * A(8) * A(4) / V1
    bind_tica <- k1 * target * ctica
    bind_tam <- k1_tam * target * ctam

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(central) <- ka * transit2 - kel * central - k12 * central + k21 * peripheral1 -
      bind_tica + k2 * complex - kmet * central + k2_rec * complex_peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_tam) <- kmet * central - k12_tam * central_tam + k21_tam * peripheral1_tam -
      kel_tam * central_tam - bind_tam + k2 * complex_tam + k2_rec * complex_peripheral1_tam
    d/dt(peripheral1_tam) <- k12_tam * central_tam - k21_tam * peripheral1_tam
    # PB2452 released by dissociation does not return to circulation (Figure 1).
    d/dt(target) <- -cl_target * ctarget - qp * (ctarget - target_peripheral1 / vp_target) -
      bind_tica - bind_tam
    d/dt(complex) <- bind_tica - k2 * complex - ktr * complex - cl_target / vc_target * complex
    d/dt(complex_tam) <- bind_tam - k2 * complex_tam - ktr * complex_tam - cl_target / vc_target * complex_tam
    # The control stream carries the PB2452 peripheral state A(11) as a
    # concentration; here it is an amount (A(11) * V_ant_perp), which is exact.
    d/dt(target_peripheral1) <- qp * (ctarget - target_peripheral1 / vp_target)
    d/dt(complex_peripheral1) <- ktr * complex - k2_rec * complex_peripheral1
    d/dt(complex_peripheral1_tam) <- ktr * complex_tam - k2_rec * complex_peripheral1_tam

    # Observations (Supplementary Material S1 $ERROR). Total TICA and total TAM
    # include the plasma complex but not the delayed complex pools.
    Cc <- ctica + complex / vc_target
    Cc_tam <- ctam + complex_tam / vc_target
    Cc_target <- ctarget
    Ctotal_target <- (target + complex + complex_tam) / vc_target
    PRU <- rbase * (1 - emax * ctica^hill / (ec50^hill + ctica^hill) -
      emax_tam * ctam^hill / (ec50_tam^hill + ctam^hill))

    Cc ~ prop(propSd) + add(addSd)
    Cc_tam ~ prop(propSd_tam)
    Cc_target ~ prop(propSd_Cc_target)
    Ctotal_target ~ prop(propSd_Ctotal_target) + add(addSd_Ctotal_target)
    PRU ~ add(addSd_PRU)
  })
}
