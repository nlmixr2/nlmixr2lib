Ruppel_2017_insulin <- function() {
  description <- "Population PK and PK-glucose infusion rate (GIR) model for inhaled Technosphere insulin (TI) and subcutaneous regular human insulin (RHI) in healthy volunteers under a euglycemic glucose clamp (Ruppel 2017). Two-compartment disposition with first-order elimination fixed from Potocka 2010; RHI is absorbed through two successive first-order steps (s.c. depot and transit compartment), TI by a single first-order step whose rate falls with body weight, and TI bioavailability is relative to RHI (so volumes and clearance are apparent, relative to RHI bioavailability). GIR follows a sigmoid Emax function of an effect compartment that equilibrates with the CENTRAL AMOUNT (the authors' code fixes the effect-compartment scale to a nominal 1 L). Inter-individual variability on the s.c. and TI absorption rates, TI bioavailability, both volumes, Emax and EC50; inter-occasion variability on TI bioavailability and on the GIR baseline; proportional error on insulin, additive error on GIR."
  reference <- paste(
    "Ruppel D, Dahmen R, Boss A, Jager R, Grant M, Baughman R, Klabunde T.",
    "A Population Dose-Response Model for Inhaled Technosphere Insulin",
    "Administered to Healthy Subjects.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(6):365-372.",
    "doi:10.1002/psp4.12189.",
    "Insulin elimination rate constant fixed from Potocka E, et al.",
    "J Diabetes Sci Technol. 2010 (Ruppel 2017 reference 8).",
    sep = " "
  )
  vignette <- "Ruppel_2017_insulin"
  units <- list(
    time = "h",
    dosing = "mU",
    concentration = "mU/L",
    GIR = "mg/kg/min"
  )
  # s.c. RHI doses go to `depot`, inhaled TI doses to `depot_lung`.
  dosing <- c("depot", "depot_lung")

  compartmentData <- list(
    depot = list(
      analyte = "insulin human (regular human insulin, s.c.)",
      units = "mU",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "insulin human (regular human insulin, s.c.)",
      units = "mU",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_lung = list(
      analyte = "insulin human (Technosphere inhalation powder)",
      units = "mU",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "insulin human", units = "mU", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "insulin human", units = "mU", specimen = "not applicable", verified = TRUE),
    effect = list(analyte = "insulin human", units = "mU", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the TI absorption rate constant, ka_TI = ka_TI,0 * (WT/70)^gamma_WT (Supplementary Table S1 footnote; PK control stream KATI = TKATI*(WT/70)**GAMWT). Range in the analysed cohort 49-100 kg.",
      source_name = "WT"
    ),
    OCC = list(
      description = "Clamp occasion index (1-5) of the five-period crossover (four TI doses and one RHI dose)",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Carries inter-occasion variability on TI bioavailability and on the GIR baseline. The PK control stream indexes the bioavailability IOV by OCC = 3-6 (the four TI occasions) and the PK-GIR stream indexes the baseline IOV by VISIT = 3-7 (all five occasions); both use $OMEGA BLOCK(1) SAME. Here both are indexed by one OCC = 1-5 column; the bioavailability IOV slot of the RHI occasion is inert because only TI doses carry a bioavailability. OCC = 0 removes both IOV contributions.",
      source_name = "OCC (PK) / VISIT (PK-GIR)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "18-53 years",
    weight_range = "49-100 kg",
    sex_female_pct = 13 / 31 * 100,
    disease_state = "Healthy volunteers",
    dose_range = "Single doses of 10, 30, 60 and 80 U inhaled Technosphere insulin and 15 IU subcutaneous regular human insulin, one per occasion of a crossover",
    regions = "United States (Profil Institute for Clinical Research, Chula Vista, CA)",
    notes = "Study NCT01490762: 32 healthy volunteers in a hyperinsulinemic euglycemic glucose clamp crossover (background intravenous insulin lispro infusion, variable 20% dextrose). One subject with occasionally very high exposure was excluded, leaving 31 (18 male, 13 female). Insulin concentrations were corrected for endogenous insulin estimated from C-peptide. Insulin sampled to 3 h after TI and GIR recorded for 4 h after TI and 10 h after RHI. Demographics in the 'Study design and population' section."
  )

  ini({
    # Absorption (Supplementary Table S1, Fit 629783)
    lka_sc <- log(0.34); label("First s.c. absorption rate constant of RHI, depot to transit, ka,sc1 (1/h)") # Table S1: ka sc1 = 0.34 (0.22-0.47) 1/h
    lktr <- log(0.86); label("Second s.c. absorption rate constant of RHI, transit to central, ka,sc2 (1/h)") # Table S1: ka sc2 = 0.86 (0.55-1.17) 1/h
    lka_lung <- log(1.64); label("TI absorption rate constant at 70 kg, ka,TI,0 (1/h)") # Table S1: ka TI,0 = 1.64 (1.49-1.78) 1/h
    e_wt_ka_lung <- -1.7; label("Power exponent of (WT/70) on the TI absorption rate constant, gamma_WT (unitless)") # Table S1: WT = -1.7 (-2.2 to -1.3)
    lfdepot_lung <- log(0.25); label("Bioavailability of TI relative to s.c. RHI, F_TI (fraction)") # Table S1: F TI = 0.25 (0.21-0.28)

    # Disposition, apparent (relative to RHI bioavailability; Table S1 footnote)
    lvc <- log(4.5); label("Apparent central volume V4 (L)") # Table S1: V4 = 4.5 (3.9-5.1) L
    lvp <- fixed(log(56.4486)); label("Apparent peripheral volume V5 (L)") # Table S1: V5 = 56.4 (-) L; PK control stream $THETA 56.448600 FIX ('V5 from fit 612372')
    lq <- fixed(log(29.5602)); label("Apparent intercompartmental clearance Q (L/h)") # Table S1: Q = 29.6 (-) L/h; PK control stream $THETA 29.560200 FIX ('Q from fit612526')
    lkel <- fixed(log(8.55)); label("Elimination rate constant from the central compartment, ke0 (1/h)") # Table S1: ke0 = 8.55 (-) 1/h; PK control stream $THETA 8.55 FIX ('Potocka 43.7/5.11'); Results: taken from Ref. 8

    # PK-GIR model (Table 1, Fit119360)
    le0 <- log(2.7); label("Baseline glucose infusion rate GIR0 (mg/kg/min)") # Table 1: GIR0 = 2.7 (2.3-3.1) mg/kg/min
    lemax <- fixed(log(11.9059)); label("Maximum insulin-induced GIR increase GIRmax (mg/kg/min)") # Table 1: GIRmax = 11.9 (-); PK-GIR control stream $THETA (11.905900 11.905900 11.905900), fixed by profiling (Supplementary Figure S2)
    lec50 <- log(152); label("Effect-compartment insulin amount producing half-maximal GIR increase, EC50 (mU, nominal 1 L effect compartment)") # Table 1: EC50 = 152 (130-174), labelled mU/L
    lhill <- log(1.8); label("Hill coefficient of the GIR Emax function, gamma (unitless)") # Table 1: gamma = 1.8 (1.5-2.1)
    lke0 <- log(1.24); label("Effect-compartment equilibration rate constant, kin = kout (1/h)") # Table 1: Kin = Kout = 1.24 (1.08-1.39) 1/h; one parameter KEFF in the PK-GIR control stream

    # Inter-individual variability. The tables print each omega as a percent
    # equal to 100 * sqrt(omega^2) (log-scale SD): the Table S1 V5 row
    # '102 (76-123)' is symmetric on the variance scale only under that
    # reading (see vignette), so omega^2 = (percent / 100)^2.
    etalka_sc ~ 0.2304 # Table S1: kasc1 = 48 %
    etalka_lung ~ 0.0441 # Table S1: kaTI = 21 %
    etalfdepot_lung ~ 0.0841 # Table S1: FTI = 29 %
    etalvc ~ 0.0484 # Table S1: V4 = 22 %
    etalvp ~ 1.0404 # Table S1: V5 = 102 %
    etalemax ~ 0.0961 # Table 1: omega Emax = 31 %
    etalec50 ~ 0.1444 # Table 1: omega EC50 = 38 %

    # Inter-occasion variability ($OMEGA BLOCK(1) SAME in both streams)
    etaiov_fdepot_lung_1 ~ 0.09 # Table S1: IOV FTI = 30 %
    etaiov_fdepot_lung_2 ~ fixed(0.09) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_fdepot_lung_3 ~ fixed(0.09) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_fdepot_lung_4 ~ fixed(0.09) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_fdepot_lung_5 ~ fixed(0.09) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_e0_1 ~ 0.3721 # Table 1: IOV GIR0 = 61 %
    etaiov_e0_2 ~ fixed(0.3721) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_e0_3 ~ fixed(0.3721) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_e0_4 ~ fixed(0.3721) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_e0_5 ~ fixed(0.3721) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)

    # Residual error. Both streams code Y = IPRED + SQRT(AERR^2 + PERR^2 * IPRED^2) * EPS
    # with $SIGMA 1 FIX, so the THETAs are standard deviations; the PK additive
    # THETA and the GIR proportional THETA are fixed to 0.
    propSd <- 0.31; label("Proportional residual error on insulin concentration (fraction)") # Table S1: p = 31 % (26-35)
    addSd_GIR <- 2.1; label("Additive residual error on GIR (mg/kg/min)") # Table 1: sigma A = 2.1 (1.9-2.4) mg/kg/min
  })
  model({
    # Occasion indicators; OCC = 0 zeroes both IOV contributions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_fdepot_lung <- oc1 * etaiov_fdepot_lung_1 + oc2 * etaiov_fdepot_lung_2 +
      oc3 * etaiov_fdepot_lung_3 + oc4 * etaiov_fdepot_lung_4 + oc5 * etaiov_fdepot_lung_5
    iov_e0 <- oc1 * etaiov_e0_1 + oc2 * etaiov_e0_2 + oc3 * etaiov_e0_3 +
      oc4 * etaiov_e0_4 + oc5 * etaiov_e0_5

    # PK
    ka_sc <- exp(lka_sc + etalka_sc)
    ktr <- exp(lktr)
    ka_lung <- exp(lka_lung + etalka_lung) * (WT / 70)^e_wt_ka_lung
    fdepot_lung <- exp(lfdepot_lung + etalfdepot_lung + iov_fdepot_lung)
    kel <- exp(lkel)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    cl <- kel * vc
    k12 <- q / vc
    k21 <- q / vp

    # PK-GIR
    e0 <- exp(le0 + iov_e0)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)
    ke0 <- exp(lke0)

    d/dt(depot) <- -ka_sc * depot
    d/dt(transit1) <- ka_sc * depot - ktr * transit1
    d/dt(depot_lung) <- -ka_lung * depot_lung
    d/dt(central) <- ktr * transit1 + ka_lung * depot_lung - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # The PK-GIR stream integrates DADT(6) = KEFF*(A(4) - A(6)) with S6 = 1:
    # the effect compartment tracks the central AMOUNT, not the concentration.
    d/dt(effect) <- ke0 * (central - effect)

    # The PK stream's logit round trip F3 = exp(PHI)/(1 + exp(PHI)) with
    # PHI = LOG(FTI/(1 - FTI)) returns FTI unchanged.
    f(depot_lung) <- fdepot_lung

    Cc <- central / vc
    Cc ~ prop(propSd)

    GIR <- e0 + emax * effect^hill / (ec50^hill + effect^hill)
    GIR ~ add(addSd_GIR)
  })
}
