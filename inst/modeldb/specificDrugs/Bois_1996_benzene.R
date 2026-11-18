Bois_1996_benzene <- function() {
  description <- paste(
    "PBPK (whole-body, human) population toxicokinetics of benzene, an",
    "environmental/occupational carcinogen (Bois et al. 1996, Environ",
    "Health Perspect). Five flow-limited tissue compartments (poorly",
    "perfused, well perfused, fat, bone marrow, liver) with instantaneous",
    "pulmonary equilibrium between alveolar air, venous and arterial",
    "blood. Michaelis-Menten metabolism in liver and bone marrow (linear",
    "over the 0.001-10 ppm range studied) feeds a phenolic-metabolite",
    "pool that, with an endogenous baseline formation term, is excreted",
    "to urine; urinary phenol concentration is a model output. The 24",
    "estimated parameters are body-weight scaling coefficients fit by",
    "Markov-chain Monte Carlo to inhalation-chamber data (blood and",
    "exhaled-air benzene, urinary phenol) from three male volunteers;",
    "their posterior population geometric means and heterogeneity (Table",
    "2) are carried as typical values and lognormal between-subject",
    "variability. Inhalation exposure enters as a time-constant air",
    "concentration (CEXP, over duration TEXP), not as a dose event."
  )
  reference <- paste(
    "Bois FY, Jackson ET, Pekari K, Smith MT. Population toxicokinetics",
    "of benzene. Environ Health Perspect. 1996;104 Suppl 6:1405-1411.",
    "doi:10.1289/ehp.961041405. Structural model and scaling functions",
    "from Table 1 and the Methods ODEs; typical values and between-",
    "subject variability from the posterior Table 2 (population geometric",
    "mean exp(mu) and population SD exp(Sigma))."
  )
  vignette <- "Bois_1996_benzene"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # Benzene tissue concentrations (mg/L) and amount pools (mg). None map to
  # a canonical compartment role, so they are whitelisted here. Venous blood
  # and exhaled (alveolar) air benzene, and urinary phenol, are algebraic
  # observations (Cc, Cexh, Cphu), not ODE states.
  paper_specific_compartments <- c(
    "Cpp",
    "Cwp",
    "Cfat",
    "Cbm",
    "Cliv",
    "Amet",
    "AmetL",
    "AmetBM",
    "Ainh"
  )

  compartmentData <- list(
    Cpp = list(analyte = "benzene", units = "mg/L", specimen = "tissue", verified = TRUE),
    Cwp = list(analyte = "benzene", units = "mg/L", specimen = "tissue", verified = TRUE),
    Cfat = list(analyte = "benzene", units = "mg/L", specimen = "tissue", verified = TRUE),
    Cbm = list(analyte = "benzene", units = "mg/L", specimen = "tissue", verified = TRUE),
    Cliv = list(analyte = "benzene", units = "mg/L", specimen = "tissue", verified = TRUE),
    Amet = list(analyte = "phenolic metabolites", units = "mg", specimen = "not applicable", verified = TRUE),
    AmetL = list(analyte = "benzene metabolized in liver", units = "mg", specimen = "not applicable", verified = TRUE),
    AmetBM = list(
      analyte = "benzene metabolized in bone marrow",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    Ainh = list(analyte = "benzene inhaled", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      source_name = "body weight (Pekari subjects 55, 73, 90 kg; Table 1 prior exp(M) = 73)",
      notes = paste(
        "Drives lean body mass LM = WT*(1 - scvf), all organ volumes and",
        "blood flows, and the allometric Vmax,liver = scvmaxl*LM^0.7. The",
        "poorly perfused volume is computed by difference so organ masses",
        "plus skeleton (17% of LM) sum to WT (Table 1 footnote e)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 3L,
    n_studies = 1L,
    age_range = "young adult males",
    weight_range = "55-90 kg (subjects 1-3: 55, 73, 90 kg)",
    sex_female_pct = 0,
    race_ethnicity = "white (3 young white male subjects)",
    disease_state = paste(
      "Healthy volunteers exposed to benzene in an inhalation chamber for",
      "4 h at 1.7 ppm (5.2 ug/L) and, one month later, 10 ppm (30 ug/L).",
      "Blood and exhaled-air benzene and urinary phenol were measured",
      "during exposure and over the following 2 days (Pekari et al. 1992)."
    ),
    dose_range = "inhalation 1.7 ppm (5.2 ug/L) and 10 ppm (30 ug/L), 4 h",
    regions = "Finland (Institute of Occupational Health, Helsinki)",
    notes = paste(
      "Bayesian (MCMC) population PBPK. Parameters are body-weight scaling",
      "coefficients, not direct physiologic quantities. Population means",
      "from Table 2 column exp(mu); between-subject variability from",
      "column exp(Sigma) carried as lognormal IIV. Primary reported",
      "results: median fraction of benzene metabolized ~52% (90% CI",
      "47-67%), independent of exposure up to 10 ppm; ~2-40 mg/day",
      "metabolized in bone marrow at continuous 1 ppm."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters are the body-weight SCALING COEFFICIENTS the
    # authors estimated (Table 1 "scaling function" column defines how
    # each maps to a physiologic quantity). Typical values are the
    # posterior population geometric means, Table 2 column exp(mu); IIV
    # variances are log(exp(Sigma))^2 from Table 2 column exp(Sigma).
    # Flows in L/min, volumes in L, Vmax in mg/min, Km in mg.
    # ------------------------------------------------------------------
    lvpr <- log(1.7)
    label("Ventilation/perfusion ratio (Qalv/Qc), unitless") # Table 2 VPR exp(mu) = 1.7

    lscfwp <- log(0.34)
    label("Blood flow per unit mass, well-perfused tissues (1/min)") # Table 2 Fwp exp(mu) = 0.34
    lscfpp <- log(0.041)
    label("Blood flow per unit mass, poorly perfused tissues (1/min)") # Table 2 Fpp exp(mu) = 0.041
    lscff <- log(0.029)
    label("Blood flow per unit mass, fat (1/min)") # Table 2 Ff exp(mu) = 0.029
    lscfbm <- log(0.19)
    label("Blood flow per unit mass, bone marrow (1/min)") # Table 2 Fbm exp(mu) = 0.19
    lscfl <- log(1.0)
    label("Blood flow per unit mass, liver (1/min)") # Table 2 Fl exp(mu) = 1.0

    lscvwp <- log(0.23)
    label("Well-perfused tissue volume as fraction of lean body mass (unitless)") # Table 2 Vwp exp(mu) = 0.23
    lscvf <- log(0.25)
    label("Fat mass as fraction of body weight (unitless)") # Table 2 Vf exp(mu) = 0.25
    lscvbm <- log(0.052)
    label("Bone marrow volume as fraction of lean body mass (unitless)") # Table 2 Vbm exp(mu) = 0.052
    lscvl <- log(0.031)
    label("Liver volume as fraction of lean body mass (unitless)") # Table 2 Vl exp(mu) = 0.031

    lpba <- log(9.65)
    label("Blood/air partition coefficient (unitless)") # Table 2 Pba exp(mu) = 9.65
    lpwp <- log(1.0)
    label("Well-perfused tissue/blood partition coefficient (unitless)") # Table 2 Pwp exp(mu) = 1.0
    lppp <- log(1.7)
    label("Poorly perfused tissue/blood partition coefficient (unitless)") # Table 2 Ppp exp(mu) = 1.7
    lpf <- log(19)
    label("Fat/blood partition coefficient (unitless)") # Table 2 Pf exp(mu) = 19
    lpbm <- log(7.67)
    label("Bone marrow/blood partition coefficient (unitless)") # Table 2 Pbm exp(mu) = 7.67
    lpl <- log(1.8)
    label("Liver/blood partition coefficient (unitless)") # Table 2 Pl exp(mu) = 1.8

    lscvmaxl <- log(5.8)
    label("Scaling coefficient for maximal metabolic rate in liver (Vmax = scvmaxl*LM^0.7, mg/min)") # Table 2 Vmaxl exp(mu) = 5.8
    lscvmaxbm <- log(0.16)
    label("Bone marrow Vmax as fraction of liver Vmax (unitless)") # Table 2 Vmaxbm exp(mu) = 0.16
    lrkml <- log(0.60)
    label("First-order metabolic rate in liver, Vmax,l/Km,l (1/min)") # Table 2 Kml scaling coefficient exp(mu) = 0.60
    lrkmbm <- log(2.8)
    label("First-order metabolic rate in bone marrow, Vmax,bm/Km,bm (1/min)") # Table 2 Kmbm scaling coefficient exp(mu) = 2.8

    lkf <- log(0.014)
    label("Endogenous phenolic-metabolite formation rate in liver (mg/min)") # Table 2 Kf exp(mu) = 0.014
    lke <- log(0.0016)
    label("Metabolite excretion rate constant to urine (1/min)") # Table 2 Ke exp(mu) = 0.0016
    lfu <- log(0.0025)
    label("Urine formation rate (L/min)") # Table 2 Fu exp(mu) = 0.0025
    lfp <- log(0.78)
    label("Phenol fraction of excreted urinary metabolites (unitless)") # Table 2 fp exp(mu) = 0.78

    # Between-subject variability: variance = log(exp(Sigma))^2 from the
    # Table 2 population-SD column (geometric SDs). Applied lognormally.
    etalvpr ~ 0.017168 # Table 2 VPR exp(Sigma) = 1.14
    etalscfwp ~ 0.003395 # Table 2 Fwp exp(Sigma) = 1.06
    etalscfpp ~ 0.002380 # Table 2 Fpp exp(Sigma) = 1.05
    etalscff ~ 0.002380 # Table 2 Ff exp(Sigma) = 1.05
    etalscfbm ~ 0.002380 # Table 2 Fbm exp(Sigma) = 1.05
    etalscfl ~ 0.002380 # Table 2 Fl exp(Sigma) = 1.05
    etalscvwp ~ 0.002380 # Table 2 Vwp exp(Sigma) = 1.05
    etalscvf ~ 0.036336 # Table 2 Vf exp(Sigma) = 1.21
    etalscvbm ~ 0.002380 # Table 2 Vbm exp(Sigma) = 1.05
    etalscvl ~ 0.002380 # Table 2 Vl exp(Sigma) = 1.05
    etalpba ~ 0.017168 # Table 2 Pba exp(Sigma) = 1.14
    etalpwp ~ 0.017168 # Table 2 Pwp exp(Sigma) = 1.14
    etalppp ~ 0.017168 # Table 2 Ppp exp(Sigma) = 1.14
    etalpf ~ 0.017168 # Table 2 Pf exp(Sigma) = 1.14
    etalpbm ~ 0.017168 # Table 2 Pbm exp(Sigma) = 1.14
    etalpl ~ 0.017168 # Table 2 Pl exp(Sigma) = 1.14
    etalscvmaxl ~ 0.118054 # Table 2 Vmaxl exp(Sigma) = 1.41
    etalscvmaxbm ~ 0.118054 # Table 2 Vmaxbm exp(Sigma) = 1.41
    etalrkml ~ 0.108441 # Table 2 Kml exp(Sigma) = 1.39
    etalrkmbm ~ 0.118054 # Table 2 Kmbm exp(Sigma) = 1.41
    etalkf ~ 0.049793 # Table 2 Kf exp(Sigma) = 1.25
    etalke ~ 0.060940 # Table 2 Ke exp(Sigma) = 1.28
    etalfu ~ 0.039542 # Table 2 Fu exp(Sigma) = 1.22
    etalfp ~ 0.017168 # Table 2 fp exp(Sigma) = 1.14

    # Inhalation exposure scenario (overridden per scenario in the
    # vignette). CEXP is the time-constant chamber air concentration while
    # t < TEXP; after TEXP exposure stops. 10 ppm = 30 ug/L = 0.030 mg/L.
    CEXP <- fixed(0.030)
    label("Inhaled benzene air concentration during exposure (mg/L)") # Methods: 10 ppm = 30 ug/L
    TEXP <- fixed(240)
    label("Exposure duration (min)") # Methods: 4-h chamber exposure

    # Lognormal (log-scale) measurement-error SDs, Results: 0.22, 0.17,
    # 0.22 in log space (CV ~ 24%, 19%, 24%) for blood, exhaled air, urine.
    expSd <- 0.22
    label("Venous blood benzene residual, lognormal log-scale SD") # Results: 0.22 +/- 0.016 (log space)
    expSd_Cexh <- 0.17
    label("Exhaled-air benzene residual, lognormal log-scale SD") # Results: 0.17 +/- 0.013 (log space)
    expSd_Cphu <- 0.22
    label("Urinary phenol residual, lognormal log-scale SD") # Results: 0.22 +/- 0.017 (log space)
  })

  model({
    # ================================================================
    # 0. Back-transform scaling coefficients (lognormal IIV)
    # ================================================================
    vpr <- exp(lvpr + etalvpr)
    scfwp <- exp(lscfwp + etalscfwp)
    scfpp <- exp(lscfpp + etalscfpp)
    scff <- exp(lscff + etalscff)
    scfbm <- exp(lscfbm + etalscfbm)
    scfl <- exp(lscfl + etalscfl)
    scvwp <- exp(lscvwp + etalscvwp)
    scvf <- exp(lscvf + etalscvf)
    scvbm <- exp(lscvbm + etalscvbm)
    scvl <- exp(lscvl + etalscvl)
    pba <- exp(lpba + etalpba)
    pwp <- exp(lpwp + etalpwp)
    ppp <- exp(lppp + etalppp)
    pf <- exp(lpf + etalpf)
    pbm <- exp(lpbm + etalpbm)
    pl <- exp(lpl + etalpl)
    scvmaxl <- exp(lscvmaxl + etalscvmaxl)
    scvmaxbm <- exp(lscvmaxbm + etalscvmaxbm)
    rkml <- exp(lrkml + etalrkml)
    rkmbm <- exp(lrkmbm + etalrkmbm)
    kf <- exp(lkf + etalkf)
    ke <- exp(lke + etalke)
    fu <- exp(lfu + etalfu)
    fp <- exp(lfp + etalfp)

    # ================================================================
    # 1. Physiology from body weight (Table 1 scaling functions).
    #    Fat mass = scvf*WT (fat volume = scvf*WT/0.92, density 0.92).
    #    Lean body mass LM = WT - fat mass. Vpp by difference so organ
    #    masses + skeleton (17% of LM) sum to WT (Table 1 footnote e).
    # ================================================================
    fatmass <- scvf * WT
    LM <- WT - fatmass
    Vfat <- fatmass / 0.92
    Vwp <- scvwp * LM
    Vbm <- scvbm * LM
    Vliv <- scvl * LM
    skel <- 0.17 * LM
    Vpp <- WT - Vwp - Vbm - Vliv - skel - fatmass

    # Blood flows per unit mass (fat flow is per unit fat MASS = Vfat*0.92)
    Fwp <- scfwp * Vwp
    Fpp <- scfpp * Vpp
    Ffat <- scff * Vfat * 0.92
    Fbm <- scfbm * Vbm
    Fliv <- scfl * Vliv
    Qc <- Fwp + Fpp + Ffat + Fbm + Fliv
    Qalv <- vpr * Qc

    # ================================================================
    # 2. Metabolic parameters. Vmax,l = scvmaxl*LM^0.7; Km,l = Vmax,l/rkml
    #    so the first-order rate Vmax,l/Km,l equals rkml. Bone marrow
    #    Vmax is a fraction scvmaxbm of liver Vmax.
    # ================================================================
    Vmaxl <- scvmaxl * LM^0.7
    Kml <- Vmaxl / rkml
    Vmaxbm <- scvmaxbm * Vmaxl
    Kmbm <- Vmaxbm / rkmbm

    # ================================================================
    # 3. Inhalation forcing and instantaneous pulmonary equilibrium.
    #    Cart = (Qc*Cven + Qalv*Cinh)/(Qc + Qalv/Pba), where
    #    Qc*Cven = sum of F_i*C_i/P_i over tissues.
    # ================================================================
    Cinh <- CEXP * (t < TEXP)
    venflux <-
      Fpp * Cpp / ppp + Fwp * Cwp / pwp + Ffat * Cfat / pf +
      Fbm * Cbm / pbm + Fliv * Cliv / pl
    Cart <- (venflux + Qalv * Cinh) / (Qc + Qalv / pba)
    Cven <- venflux / Qc

    # Michaelis-Menten metabolic loss (concentration rate, mg/L/min).
    metabL <- Vmaxl * Cliv / (Kml + Cliv * Vliv)
    metabBM <- Vmaxbm * Cbm / (Kmbm + Cbm * Vbm)

    # ================================================================
    # 4. Tissue ODEs (flow-limited); liver and bone marrow metabolize.
    # ================================================================
    d / dt(Cpp) <- Fpp / Vpp * (Cart - Cpp / ppp)
    d / dt(Cwp) <- Fwp / Vwp * (Cart - Cwp / pwp)
    d / dt(Cfat) <- Ffat / Vfat * (Cart - Cfat / pf)
    d / dt(Cbm) <- Fbm / Vbm * (Cart - Cbm / pbm) - metabBM
    d / dt(Cliv) <- Fliv / Vliv * (Cart - Cliv / pl) - metabL

    # ================================================================
    # 5. Phenolic-metabolite pool and cumulative amounts. The pool gains
    #    benzene metabolized in liver and bone marrow (amount rate =
    #    concentration rate * tissue volume) plus endogenous Kf, and is
    #    excreted first-order (Ke) to urine. AmetL/AmetBM/Ainh accumulate
    #    benzene metabolized and inhaled for the fraction-metabolized
    #    check. The pool starts at its endogenous steady state Kf/Ke.
    # ================================================================
    d / dt(Amet) <- metabL * Vliv + metabBM * Vbm + kf - ke * Amet
    d / dt(AmetL) <- metabL * Vliv
    d / dt(AmetBM) <- metabBM * Vbm
    d / dt(Ainh) <- Qalv * Cinh
    Amet(0) <- kf / ke

    # ================================================================
    # 6. Observations: venous blood benzene, exhaled (alveolar) air
    #    benzene (= Cart/Pba), urinary phenol concentration (fraction fp
    #    of the metabolite excretion rate divided by the urine rate).
    # ================================================================
    Cc <- Cven
    Cexh <- Cart / pba
    Cphu <- fp * ke * Amet / fu

    Cc ~ lnorm(expSd)
    Cexh ~ lnorm(expSd_Cexh)
    Cphu ~ lnorm(expSd_Cphu)
  })
}
