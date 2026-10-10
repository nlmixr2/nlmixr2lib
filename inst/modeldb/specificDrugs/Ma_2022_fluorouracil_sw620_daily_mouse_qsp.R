Ma_2022_fluorouracil_sw620_daily_mouse_qsp <- function() {
  description <- "QSP. Preclinical (mouse, human SW620 colon carcinoma tumour). Ma 2022 integrated mechanistic PK / cellular / tumour growth inhibition (TGI) model of 5-fluorouracil (5-FU) in colon cancer, TGI model (5) of nine (5 mg/kg daily for 9 days, human SW620 colon carcinoma in immunodeficient mice (xenograft); Ma 2022 Table 3). A two-compartment plasma PK model with Michaelis-Menten elimination (fitted to mouse plasma after a single 100 mg/kg intraperitoneal dose) drives a cellular sub-model of 5-FU in tumour interstitial fluid, saturable uptake into and efflux out of tumour cells, anabolism to fluoronucleotides (FNUC), FUTP / FdUTP incorporation into RNA and DNA with delayed excision, FdUMP binding to thymidylate synthase (TS) in competition with dUMP, the resulting dNTP-pool imbalance and nucleotide-salvage feedback, and delayed DNA double-strand-break (DSB) induction. Anabolites inhibit tumour growth (sigmoid Imax on logistic growth) and DSBs above a drug-free baseline convert proliferating cells into damaged cells that die after a constant life span (life-span / delay-differential TGI model). The PK and cellular parameters are shared by all nine TGI models; only the control-growth, efficacy and residual-error parameters differ. Time is in minutes. The plasma, interstitial and cellular states are concentrations, so a dose is entered as the instantaneous plasma-concentration increment it produces (ug/mL; 100 mg/kg -> 92.9936 ug/mL in the authors' code). The five delay terms use rxode2 delay(), which needs a dense solver: solve with method = 'dop853'."
  reference <- paste(
    "Ma C, Almasan A, Gurkan-Cavusoglu E. (2022). Computational analysis of",
    "5-fluorouracil anti-tumor activity in colon cancer using a mechanistic",
    "pharmacokinetic/pharmacodynamic model. PLoS Comput Biol 18(11):e1010685.",
    "doi:10.1371/journal.pcbi.1010685.",
    "Code and MCMC chains deposited by the authors at",
    "https://github.com/Chenhui88569/Computational-analysis-of-5-fluorouracil-anti-tumor-activity",
    "(Zenodo doi:10.5281/zenodo.7267874)."
  )
  vignette <- "Ma_2022_fluorouracil_qsp"

  # Mechanistic states of the Ma 2022 cellular sub-model with no canonical
  # analogue in inst/references/compartment-names.md (A4, A5, A6, A7, A8, A9,
  # Nucp, U and the drug-free DSB reference trajectory N_DSB,0).
  paper_specific_compartments <- c(
    "intracellular",
    "anabolites",
    "fu_rna",
    "fu_dna",
    "dump",
    "ts_fdump",
    "nucp",
    "salvage",
    "dsb",
    "dsb_baseline"
  )

  units <- list(
    time = "min",
    dosing = "ug/mL (instantaneous plasma-concentration increment; the PK states are concentrations)",
    concentration = "ug/mL for plasma Cc; pmol/mg tissue for the cellular species; cm^3 for tumour volume"
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "ug/mL", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "ug/mL", specimen = "tissue", verified = TRUE),
    isf = list(analyte = "5-fluorouracil", units = "pmol/mg tissue", specimen = "tumor", verified = TRUE),
    intracellular = list(analyte = "5-fluorouracil", units = "pmol/mg tissue", specimen = "tumor", verified = TRUE),
    anabolites = list(
      analyte = "5-FU anabolites (fluoronucleotides and fluoronucleosides, FNUC)",
      units = "pmol/mg tissue",
      specimen = "tumor",
      verified = TRUE
    ),
    fu_rna = list(
      analyte = "FUTP incorporated into RNA (FU-RNA)",
      units = "pmol/mg tissue",
      specimen = "tumor",
      verified = TRUE
    ),
    fu_dna = list(
      analyte = "FdUTP incorporated into DNA (FU-DNA)",
      units = "pmol/mg tissue",
      specimen = "tumor",
      verified = TRUE
    ),
    dump = list(
      analyte = "total dUMP (free plus TS-bound)",
      units = "pmol/mg tissue",
      specimen = "tumor",
      verified = TRUE
    ),
    ts_fdump = list(analyte = "TS-FdUMP complex", units = "pmol/mg tissue", specimen = "tumor", verified = TRUE),
    nucp = list(
      analyte = "dNTP-pool imbalance (Eq 16 perturbation measure)",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    salvage = list(
      analyte = "nucleotide-salvage pathway activity U",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    dsb = list(
      analyte = "DNA double-strand breaks",
      units = "thousands of breaks",
      specimen = "tumor",
      verified = TRUE
    ),
    dsb_baseline = list(
      analyte = "DNA double-strand breaks in the absence of 5-FU (N_DSB,0)",
      units = "thousands of breaks",
      specimen = "tumor",
      verified = TRUE
    ),
    cycling_cells = list(analyte = "proliferating tumour cells P", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "damaged quiescent tumour cells D",
      units = "cm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUM_VOL = list(
      description = "Tumour volume at the start of treatment; initial condition of the proliferating-cell state (P(0) = T0, Ma 2022 Eq 26).",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = "Ma 2022 sets P(0) to the last tumour-volume measurement before dosing starts. The values used for the nine TGI models are in the authors' deposited code (CollectDataIntoNestedStructure.m): 125, 125, 100, 100, 56.76, 56.5, 71.3, 125 and 125 mm^3 for models 1-9; this model's value is 56.76 mm^3. Divided by 1000 in model() because the TGI states are in cm^3.",
      source_name = "InitialTv"
    )
  )

  population <- list(
    species = "mouse (tumour-bearing; human SW620 colon carcinoma in immunodeficient mice (xenograft)) with cellular-layer data pooled from rat and in-vitro sources",
    n_subjects = "Literature mean data only: no individual animals were modelled. Each sub-model was fitted to published group means (Ma 2022 Table 2 and Table 3).",
    n_studies = 12L,
    disease_state = "Subcutaneous colon tumour (human SW620 colon carcinoma); TGI model (5) of Ma 2022",
    dose_range = "PK and cellular layers fitted at a single 100-115 mg/kg intraperitoneal dose (80 mg/kg for the RNA, DNA and dUMP data); this TGI model at 5 mg/kg daily for 9 days (Ma 2022 Table 3, reference [75])",
    regions = "Literature data (Ma 2022 Tables 2 and 3)",
    notes = "PK: C57BL/6 mice bearing Colon 38 tumours, 100 mg/kg IP (Ma 2022 ref [65]). Cellular layer: intracellular 5-FU and anabolites in C57BL/6 mice with C38 tumours at 115 mg/kg [66]; FU-RNA and FU-DNA in BALB/c mice with Colon 26-B tumours at 80 mg/kg [67, 68]; TS-FdUMP and total TS in Wistar rats with DMH-induced colon carcinoma at 100 mg/kg [57]; dUMP in BALB/c mice with colon tumour 06/A at 80 mg/kg [69]; dNTP pool in L5178 lymphoma cells in vitro (1 uM) [70]; DSBs in OCI-AML2 leukaemia cells in vitro (3 umol/L) [61]. The literature tumour volumes were normalised to their starting values. No inter-individual variability: the model is deterministic, with Bayesian (MCMC) posterior means as point estimates."
  )

  ini({
    # ==================================================================
    # PK model -- Ma 2022 Table 1 'PK model', posterior means. The plasma
    # states are concentrations (ug/mL), so every PK constant enters
    # model() as a rate: Q/V1, Q/V2 and Vmax/(Km + C).
    # ==================================================================
    lvmax <- log(329.539)
    label("Maximal plasma elimination rate Vmax1 (ug/mL/min)")  # Table 1, Vmax
    lkm <- log(17334.977)
    label("Michaelis-Menten constant of plasma elimination Km1 (ug/mL)")  # Table 1, Km
    lq <- log(320.071)
    label("Intercompartmental clearance Q12 (mL/mg/min)")  # Table 1, Q21
    lvc <- log(16823.472)
    label("Central volume V1 (mL/mg tissue)")  # Table 1, V1
    lvp <- log(2999.821)
    label("Peripheral volume V2 (mL/mg tissue)")  # Table 1, V2

    # ==================================================================
    # 5-FU in interstitial fluid and its anabolism -- Table 1.
    # ==================================================================
    q_isf <- 0.981
    label("Plasma-to-interstitial-fluid transfer Q31 (mL/mg/min)")  # Table 1, Q31
    vmax_influx <- 538.451
    label("Maximal uptake rate into tumour cells Vmax,influx (pmol/mg/min)")  # Table 1, Vmax,influx
    km_influx <- 171.153
    label("Half-saturation constant of uptake Km,influx (pmol/mg)")  # Table 1, Km,influx
    vmax_efflux <- 536.642
    label("Maximal efflux rate out of tumour cells Vmax,efflux (pmol/mg/min)")  # Table 1, Vmax,efflux
    km_efflux <- 6.704
    label("Half-saturation constant of efflux Km,efflux (pmol/mg)")  # Table 1, Km,efflux
    k03 <- 135.995
    label("First-order elimination from interstitial fluid k03 (1/min)")  # Table 1, k03
    vmax_anab <- 27.138
    label("Maximal rate of 5-FU anabolism Vmax,54 (pmol/mg/min)")  # Table 1, Vmax,54
    km_anab <- 2325.154
    label("Half-saturation constant of 5-FU anabolism Km,54 (pmol/mg)")  # Table 1, Km,54

    # ==================================================================
    # 5-FU incorporation into RNA and DNA -- Table 1. Td,RNA and Td,DNA
    # are reported in days and converted to minutes in model().
    # ==================================================================
    k65 <- 0.089
    label("First-order incorporation of FUTP into RNA k65 (1/min)")  # Table 1, k65
    k56 <- 0.521
    label("First-order removal of FUTP from RNA k56 (1/min)")  # Table 1, k56
    k06 <- 0.114
    label("Delayed-excision coefficient for FU-RNA k06 (mg/pmol; Table 1 prints 1/min)")  # Table 1, k06
    gamma_rna <- 0.310
    label("Shape exponent of delayed FU-RNA excision gamma_lag,RNA (unitless)")  # Table 1, gamma_lag,RNA
    vmax_dna_in <- 0.0199
    label("Maximal FdUTP incorporation rate into DNA Vmax,75 (pmol/mg/min)")  # Table 1, Vmax,75
    km_dna_in <- 0.8398
    label("Half-saturation constant of FdUTP incorporation Km,75 (pmol/mg)")  # Table 1, Km,75
    vmax_dna_out <- 5.534
    label("Maximal removal rate of FdUTP from DNA Vmax,7 (pmol/mg/min)")  # Table 1, Vmax,57
    km_dna_out <- 0.95908
    label("Half-saturation constant of FdUTP removal Km,7 (pmol/mg)")  # Table 1, Km,57
    k07 <- 0.808
    label("Delayed-excision coefficient for FU-DNA k07 (mg/pmol; Table 1 prints 1/min)")  # Table 1, k07
    gamma_dna <- 0.554
    label("Shape exponent of delayed FU-DNA excision gamma_lag,DNA (unitless)")  # Table 1, gamma_lag,DNA
    td_rna <- 0.396
    label("Delay of FU-RNA excision Td,RNA (day)")  # Table 1, Td,RNA
    td_dna <- 4.360
    label("Delay of FU-DNA excision Td,DNA (day)")  # Table 1, Td,DNA

    # ==================================================================
    # TS inhibition -- Table 1. alpha_TS, k_d and TS0 are printed without a
    # credible interval and are constants in the authors' code.
    # ==================================================================
    k95 <- 0.038
    label("TS-FdUMP complex formation rate constant k95 (mg/pmol/min; Table 1 prints 1/min)")  # Table 1, k95
    k59 <- 2.882
    label("TS-FdUMP complex dissociation rate constant k59 (1/min)")  # Table 1, k59
    k09 <- 0.152
    label("TS-FdUMP complex degradation rate constant k09 (1/min)")  # Table 1, k09
    kd_dump <- 21.935
    label("Dissociation constant of the TS-dUMP complex KdUMP (pmol/mg)")  # Table 1, KdUMP
    g0_dump <- 0.155
    label("Zero-order dUMP synthesis rate G0 (pmol/mg/min; Table 1 prints mg/min)")  # Table 1, G0
    kcat <- 17.430
    label("TS catalytic rate constant kcat (1/min)")  # Table 1, kcat
    k08 <- 0.032
    label("First-order dUMP degradation rate constant k08 (1/min)")  # Table 1, k08
    alpha_ts <- fixed(2.021)
    label("Ratio of drug-stimulated to normal TS synthesis alpha_TS (unitless)")  # Table 1, alpha_TS
    kd_ts <- fixed(1.034)
    label("Rate constant of the rise in total TS after a dose k_d (1/min)")  # Table 1, k_d
    ts0 <- fixed(0.0186)
    label("Total TS before 5-FU TS0 (pmol/mg)")  # Table 1, TS0

    # ==================================================================
    # dNTP pool imbalance -- Table 1 (all unitless).
    # ==================================================================
    k1_dntp <- 0.213
    label("Maximal TS-driven production of dNTP imbalance k1 (1/min)")  # Table 1, k1
    k2_dntp <- 32.50
    label("Half-saturation term of dNTP-imbalance production k2 (unitless)")  # Table 1, k2
    k3_dntp <- 3.713
    label("Self-amelioration rate of dNTP imbalance k3 (1/min)")  # Table 1, k3
    k4_dntp <- 73.389
    label("Saturation coefficient of dNTP self-amelioration k4 (unitless)")  # Table 1, k4
    k5_dntp <- 3.695
    label("Salvage-mediated reduction of dNTP imbalance k5 (1/min)")  # Table 1, k5
    k6_dntp <- 0.303
    label("Saturation coefficient of salvage activity k6 (unitless)")  # Table 1, k6
    k7_dntp <- 3.339
    label("TS-perturbation attenuation coefficient of salvage effect k7 (unitless)")  # Table 1, k7
    k8_dntp <- 1.382
    label("Half-saturation term of salvage-activity production k8 (unitless)")  # Table 1, k8
    k9_dntp <- 0.27672
    label("Loss rate of salvage activity k9 (1/min)")  # Table 1, k9
    kb_dntp <- 35.65
    label("Saturation coefficient of dNTP imbalance in the salvage term kB (unitless)")  # Table 1, kB
    ka_dntp <- 6.038
    label("Saturation coefficient of salvage-activity loss kA (unitless)")  # Table 1, kA
    k10_dntp <- 0.0794
    label("Maximal TS-driven production of salvage activity k10 (1/min)")  # Table 1, k10
    gamma_dntp <- fixed(0.15)
    label("Hill exponent of the TS perturbation gamma_dNTP (unitless)")  # Table 1, gamma_dNTP

    # ==================================================================
    # DSB induction -- Table 1. K_dNTP appears in Eq 21 but is missing from
    # Table 1; its value is the posterior mean in the authors' deposited
    # MCMC chain (Para_dNTP_DSB_col_2.mat, the same 20% burn-in that
    # reproduces every printed Table 1 value).
    # ==================================================================
    vmax_dntp <- 10.516
    label("Maximal fold increase in DSB generation by dNTP imbalance Vmax,dNTP (unitless)")  # Table 1, Vmax,dNTP
    km_dntp <- 151.918
    label("Half-saturation dNTP imbalance for DSB generation K_dNTP (unitless)")  # deposited MCMC chain posterior mean; absent from Table 1
    vmax_hr <- 10.082
    label("Maximal DSB repair rate by homologous recombination Vmax,HR (thousand breaks/min)")  # Table 1, Vmax,HR
    km_hr <- 194.020
    label("Half-saturation DSB count for repair Km,HR (thousand breaks)")  # Table 1, Km,HR
    ki_hr <- 0.121
    label("Inhibition coefficient of repair by dNTP imbalance ki (unitless)")  # Table 1, ki
    k0_dsb <- 2.803
    label("Zero-order DSB generation rate k0 (thousand breaks/min)")  # Table 1, k0
    gamma_dsb <- fixed(0.6)
    label("Shape exponent of the delayed DSB response gamma_DSB (unitless)")  # Table 1, gamma_DSB
    td_dsb <- 1.801
    label("Delay of DSB induction after dNTP imbalance Td,DSB (day)")  # Table 1, Td,DSB

    # ==================================================================
    # TGI model (5): 5 mg/kg daily for 9 days, human SW620 colon carcinoma -- Ma 2022 Table 4.
    # ==================================================================
    td_tv <- 1253.7
    label("Mean life span of damaged tumour cells Td,Tv (min)")  # Table 4, Model 5
    ic50 <- 0.69275
    label("Anabolite level halving the tumour growth rate IC50 (pmol/mg)")  # Table 4, Model 5
    emax_damage <- 0.78173
    label("Maximal DSB-driven fractional increase in damaged-cell production Emax,damage (unitless)")  # Table 4, Model 5
    ec50_damage <- 0.57608
    label("Relative DSB excess giving half-maximal damage EC50,damage (unitless)")  # Table 4, Model 5
    lambda_g <- 0.00022528
    label("Tumour growth rate lambda_g (1/min)")  # Table 4, Model 5
    pmax <- 2.3627
    label("Tumour carrying capacity Pmax (cm^3)")  # Table 4, Model 5
    lambda_d <- 3.1695e-06
    label("Natural tumour-cell death rate lambda_d (1/min)")  # Table 4, Model 5
    gamma_tv <- fixed(0.2)
    label("Shape exponent of the anabolite and DSB effects gamma_Tv (unitless)")  # Table 4, gamma_Tv, footnote a

    # ==================================================================
    # Residual error. Ma 2022 uses additive normal noise on every fitted
    # variable and reports the noise VARIANCE s^2; each SD below is
    # sqrt(s^2). The PK and FU-RNA variances are not printed in Table 1 and
    # come from the deposited MCMC chains (Para_PK_col.mat and
    # Para_Cellular_col.mat, posterior means over the authors' own burn-in,
    # which reproduces every printed Table 1 value to the last digit).
    # ==================================================================
    addSd <- 7.21890
    label("Additive residual SD of plasma 5-FU (ug/mL)")  # sqrt(52.112), deposited PK chain posterior mean
    addSd_intracellular <- 14.1421
    label("Additive residual SD of intracellular 5-FU (pmol/mg)")  # Table 1, s2_intra = 199.9999
    addSd_anabolites <- 14.1423
    label("Additive residual SD of 5-FU anabolites (pmol/mg)")  # Table 1, s2_FNUC = 200.004
    addSd_fu_rna <- 2.23406
    label("Additive residual SD of FU-RNA (pmol/mg)")  # sqrt(4.991), deposited cellular chain posterior mean
    addSd_fu_dna <- 1.41739
    label("Additive residual SD of FU-DNA (pmol/mg)")  # Table 1, s2_DNA = 2.009
    addSd_ts_fdump <- 0.450555
    label("Additive residual SD of the TS-FdUMP complex (pmol/mg)")  # Table 1, s2_TS = 0.203
    addSd_dump <- 1.41421
    label("Additive residual SD of dUMP (pmol/mg)")  # Table 1, s2_dUMP = 2
    addSd_nucp <- 0.114018
    label("Additive residual SD of dNTP-pool imbalance (unitless)")  # Table 1, s2_dNTP = 0.013
    addSd_dsb <- 14.1421
    label("Additive residual SD of DSB count (thousand breaks)")  # Table 1, s2_DSB = 199.9999
    addSd_tumor_vol <- 0.0250777
    label("Additive residual SD of treated-group tumour volume (cm^3)")  # Table 4, Model 5, s2_treated = 0.00062889
    addSd_tumor_vol_control <- 0.496276
    label("Additive residual SD of control-group tumour volume (cm^3)")  # Table 4, Model 5, s2_control = 0.24629
  })

  model({
    # Time is in minutes throughout (Ma 2022 Tables 1 and 4).
    vmax <- exp(lvmax)
    km <- exp(lkm)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)

    # 5-FU molar mass (g/mol) converting plasma ug/mL to pmol/mL for the
    # Q31 * Cp influx term (Ma 2022 Results, 'PK model'; 130.077 g/mol).
    mw_fu <- 130.077

    # Delays reported in days (Table 1) -> minutes.
    td_rna_min <- td_rna * 1440
    td_dna_min <- td_dna * 1440
    td_dsb_min <- td_dsb * 1440

    # Total TS rises mono-exponentially after each dose (Eq 10). The
    # authors restart the integration clock at every dose, so t in Eq 10 is
    # the time since the most recent dose; before any dose it is the time
    # since the start of the simulation.
    tsd <- tad()
    if (is.na(tsd)) tsd <- t
    ts_total <- ts0 * (alpha_ts + (1 - alpha_ts) * exp(-kd_ts * tsd))
    # Free FdUMP-binding sites (Eq 11) and the rapid-equilibrium TS-dUMP
    # complex in terms of total dUMP (Eq 13).
    ts_free <- ts_total - ts_fdump
    ts_dump <- dump * ts_free / (kd_dump + ts_free)

    # ---- PK model (Eqs 1-3), written on concentrations -----------------
    d/dt(central) <- q / vc * (peripheral1 - central) - vmax * central / (km + central)
    d/dt(peripheral1) <- q / vp * (central - peripheral1)
    Cc <- central

    # ---- 5-FU in interstitial fluid and its anabolism (Eqs 4-6) --------
    influx <- vmax_influx * isf / (km_influx + isf)
    efflux <- vmax_efflux * intracellular / (km_efflux + intracellular)
    anabolism <- vmax_anab * intracellular / (km_anab + intracellular)
    d/dt(isf) <- q_isf * central * 1e6 / mw_fu + efflux - influx - k03 * isf
    d/dt(intracellular) <- influx - efflux - anabolism

    # ---- FU-RNA / FU-DNA with delayed excision (Eqs 6-8) ---------------
    # A6(t) = A7(t) = 0 before the first dose, which is the constant
    # initial-condition history delay() uses.
    rna_out <- k56 * fu_rna / (1 + k06 * delay(fu_rna, td_rna_min))^gamma_rna
    dna_in <- vmax_dna_in * anabolites / (km_dna_in + anabolites)
    dna_out <- vmax_dna_out * fu_dna / (km_dna_out + fu_dna) /
      (1 + k07 * delay(fu_dna, td_dna_min))^gamma_dna
    ts_binding <- k95 * (ts_total - ts_dump) * anabolites
    d/dt(anabolites) <- anabolism + rna_out + k59 * ts_fdump -
      (k65 * anabolites + dna_in + ts_binding)
    d/dt(fu_rna) <- k65 * anabolites - rna_out
    d/dt(fu_dna) <- dna_in - dna_out

    # ---- TS-FdUMP complex and dUMP (Eqs 9, 14, 15) ---------------------
    d/dt(ts_fdump) <- ts_binding - k59 * ts_fdump - k09 * ts_fdump
    d/dt(dump) <- g0_dump - kcat * ts_dump - k08 * dump
    dump(0) <- 3.878

    # ---- dNTP pool imbalance (Eqs 17-19) -------------------------------
    # The perturbation enters as (1 - TSf/TStotal)^gamma_dNTP; the max()
    # keeps the base non-negative when solver round-off pushes the
    # TS-FdUMP state a hair below zero.
    ts_p <- max(1 - ts_free / ts_total, 0)^gamma_dntp
    d/dt(nucp) <- k1_dntp * ts_p / (k2_dntp^gamma_dntp + ts_p) -
      k3_dntp * nucp / (1 + k4_dntp * nucp) -
      k5_dntp * salvage * nucp / ((1 + k6_dntp * salvage) *
        (1 + k7_dntp^gamma_dntp * ts_p) * (1 + kb_dntp * nucp))
    d/dt(salvage) <- k10_dntp * ts_p / (k8_dntp^gamma_dntp + ts_p) -
      k9_dntp * salvage / (1 + ka_dntp * salvage)

    # ---- DSB induction (Eqs 20-22) -------------------------------------
    # Nucp(t) = 0 before the first dose, so for t <= Td,DSB Eq 21 reduces to
    # Eq 20. max() guards the fractional power against round-off below 0.
    nucp_lag <- max(delay(nucp, td_dsb_min), 0)
    d/dt(dsb) <- k0_dsb * (1 + vmax_dntp * nucp_lag^gamma_dsb / (km_dntp^gamma_dsb + nucp_lag^gamma_dsb)) -
      vmax_hr * dsb / (km_hr + dsb) / (nucp_lag^gamma_dsb * ki_hr^gamma_dsb + 1)
    # N_DSB,0(t): the same DSB equation without drug (Nucp = 0).
    d/dt(dsb_baseline) <- k0_dsb - vmax_hr * dsb_baseline / (km_hr + dsb_baseline)
    dsb(0) <- 42.755
    dsb_baseline(0) <- 42.755
    dsb_deviation <- max((dsb - dsb_baseline) / dsb_baseline, 0)
    dsb_deviation_lag <- max((delay(dsb, td_tv) - delay(dsb_baseline, td_tv)) / delay(dsb_baseline, td_tv), 0)

    # ---- Tumour growth inhibition (Eqs 23-28) --------------------------
    e_a <- ic50^gamma_tv / (ic50^gamma_tv + max(anabolites, 0)^gamma_tv)
    e_dsb <- emax_damage * dsb_deviation^gamma_tv / (ec50_damage^gamma_tv + dsb_deviation^gamma_tv)
    e_dsb_lag <- emax_damage * dsb_deviation_lag^gamma_tv / (ec50_damage^gamma_tv + dsb_deviation_lag^gamma_tv)
    cycling_cells(0) <- TUM_VOL / 1000
    d/dt(cycling_cells) <- e_a * lambda_g * cycling_cells * (1 - cycling_cells / pmax) -
      (1 + e_dsb) * lambda_d * cycling_cells
    # Life-span model: cells damaged at t - Td,Tv leave the tumour at t.
    # Before Td,Tv the delayed DSB excess is 0, so the delayed term vanishes
    # as Eq 26 requires (P(t) = D(t) = 0 for t < 0).
    d/dt(damaged_cells1) <- e_dsb * lambda_d * cycling_cells -
      e_dsb_lag * lambda_d * delay(cycling_cells, td_tv)
    tumor_vol <- cycling_cells + damaged_cells1
    # Ma 2022 fits the drug-free (control) and the treated tumour volumes
    # with separate noise variances; both arms follow these same equations
    # (a control animal simply receives no dose).
    tumor_vol_control <- tumor_vol

    Cc ~ add(addSd)
    intracellular ~ add(addSd_intracellular)
    anabolites ~ add(addSd_anabolites)
    fu_rna ~ add(addSd_fu_rna)
    fu_dna ~ add(addSd_fu_dna)
    ts_fdump ~ add(addSd_ts_fdump)
    dump ~ add(addSd_dump)
    nucp ~ add(addSd_nucp)
    dsb ~ add(addSd_dsb)
    tumor_vol ~ add(addSd_tumor_vol)
    tumor_vol_control ~ add(addSd_tumor_vol_control)
  })
}
