vandenBrink_2019_quinpirole_rat_plasma <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Cluster-based PK/PD model of the plasma",
    "metabolomic (biogenic amine and amino acid) response to the dopamine",
    "D2/3 agonist quinpirole. Quinpirole PK is the two-compartment",
    "plasma / striatal brain extracellular fluid model of van den Brink 2018",
    "(linear elimination from plasma, asymmetric plasma-to-brain influx and",
    "brain-to-plasma efflux clearances). Twenty-eight plasma biomarkers are",
    "grouped into seven response clusters, each driven by the unbound plasma",
    "quinpirole concentration: clusters 1 and 4 are stimulatory Emax",
    "turnover models whose drug effect is delayed by a gamma-shaped transit",
    "function of time after dose, clusters 2, 3 and 7 are inhibitory Imax",
    "turnover models, and clusters 5 and 6 are inhibitory Imax pool",
    "(precursor-release) models. Each cluster output is the biomarker level",
    "relative to its typical pre-dose baseline, so the absolute level of a",
    "metabolite is its baseline times the output of its cluster."
  )
  reference <- paste(
    "van den Brink WJ, Hartman R, van den Berg DJ, Flik G,",
    "Gonzalez-Amoros B, Koopman N, Elassais-Schaap J, van der Graaf PH,",
    "Hankemeier T, de Lange ECM. Blood-Based Biomarkers of Quinpirole",
    "Pharmacology: Cluster-Based PK/PD and Metabolomics to Unravel the",
    "Underlying Dynamics in Rat Plasma and Brain. CPT Pharmacometrics Syst",
    "Pharmacol. 2019;8(2):107-117. doi:10.1002/psp4.12370.",
    "Quinpirole PK model: van den Brink WJ, van den Berg DJ, Bonsel FEM,",
    "Hartman R, Wong YC, van der Graaf PH, de Lange ECM. Fingerprints of",
    "CNS drug effects: a plasma neuroendocrine reflection of D2 receptor",
    "activation using multi-biomarker pharmacokinetic/pharmacodynamic",
    "modelling. Br J Pharmacol. 2018;175(19):3832-3843.",
    "doi:10.1111/bph.14452."
  )
  vignette <- "vandenBrink_2019_quinpirole_rat"
  units <- list(time = "min", dosing = "ug", concentration = "ug/L")

  paper_specific_compartments <- c(
    "resp_c1",
    "resp_c2",
    "resp_c3",
    "resp_c4",
    "resp_c5",
    "resp_c6",
    "resp_c7",
    "pool_c5",
    "pool_c6"
  )

  compartmentData <- list(
    central = list(analyte = "quinpirole", units = "ug", specimen = "plasma", verified = TRUE),
    brain_striatum = list(analyte = "quinpirole", units = "ug", specimen = "brain ISF", verified = TRUE),
    resp_c1 = list(
      analyte = "plasma biomarker cluster 1 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c2 = list(
      analyte = "plasma biomarker cluster 2 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c3 = list(
      analyte = "plasma biomarker cluster 3 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c4 = list(
      analyte = "plasma biomarker cluster 4 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c5 = list(
      analyte = "plasma biomarker cluster 5 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c6 = list(
      analyte = "plasma biomarker cluster 6 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c7 = list(
      analyte = "plasma biomarker cluster 7 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    pool_c5 = list(
      analyte = "plasma biomarker cluster 5 release pool",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    pool_c6 = list(
      analyte = "plasma biomarker cluster 6 release pool",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 44,
    n_studies = 1,
    disease_state = "healthy, with femoral artery / vein cannulas and bilateral striatal microdialysis probes",
    dose_range = paste(
      "0 (n = 12), 0.17 (n = 16) or 0.86 (n = 16) mg/kg quinpirole IV on",
      "experiment days 1 and 8, with the same dose given SC daily in between"
    ),
    weight_range = "not reported",
    regions = "Leiden, The Netherlands",
    notes = paste(
      "Plasma sampled at -5, 5, 7.5, 10, 15, 25, 45, 90, 120 and 180 min.",
      "Biomarkers measured on two platforms (BrainsOnline monoamine platform",
      "A, day 1 only; Accq-tag UPLC-MS/MS biogenic amine platform B, days 1",
      "and 8). The quinpirole PK model was fitted in the parent study to",
      "0.17-2.14 mg/kg IV data (van den Brink 2018). Cluster membership of",
      "the 28 plasma biomarkers (Figure 1): cluster 1 DL-3-aminoisobutyric",
      "acid, alpha-aminobutyric acid; cluster 2 homocysteine, glycine,",
      "serine; cluster 3 methionine sulfone, citrulline, homoserine,",
      "methionine, threonine, glycine (platform A), alanine, asparagine;",
      "cluster 4 5-HIAA, taurine, S-methylcysteine; cluster 5 glutamic acid",
      "(platforms A and B), valine, isoleucine, histidine, leucine,",
      "phenylalanine; cluster 6 serotonin; cluster 7 histamine, arginine,",
      "methionine sulfoxide, gamma-glutamyl-alanine."
    )
  )

  ini({
    # Quinpirole PK (van den Brink 2018 Br J Pharmacol Table 1; structure
    # Supporting Information Equation S3). Clearances are reported in L/h and
    # converted to L/min in model() because the PD time unit is minutes.
    lcl <- log(0.71); label("Elimination clearance from plasma CL_PL,o (L/h)") # vdB 2018 Table 1: CL_PL,o = 0.71 L/h (RSE 9%)
    lqin <- log(2.5); label("Plasma-to-brain ECF distribution clearance CL_PL,ECF (L/h)") # vdB 2018 Table 1: CL_PL,ECF = 2.5 L/h (RSE 20%)
    lqout <- log(0.52); label("Brain ECF-to-plasma distribution clearance CL_ECF,PL (L/h)") # vdB 2018 Table 1: CL_ECF,PLASMA = 0.52 L/h (RSE 24%)
    lvc <- log(1.0); label("Plasma (central) volume V_CENTRAL (L)") # vdB 2018 Table 1: V_CENTRAL = 1.0 L (RSE 6%)
    lvbr <- log(0.12); label("Brain ECF volume V_ECF (L)") # vdB 2018 Table 1: V_ECF = 0.12 L (RSE 13%)

    # Cluster 1: Emax turnover with transit-delayed drug effect (Table 3, plasma)
    lemax_c1 <- log(46.5); label("Cluster 1 maximal stimulation Emax (fraction of baseline)") # Table 3 plasma cluster 1: Emax = 4650% (RSE 41.1%)
    lec50_c1 <- log(383); label("Cluster 1 EC50 of unbound plasma quinpirole (nM)") # Table 3 plasma cluster 1: EC50 = 383 nM (RSE 54.3%)
    lkout_c1 <- log(0.035); label("Cluster 1 turnover rate constant kout (1/min)") # Table 3 plasma cluster 1: kout = 0.035 1/min (RSE 42.3%)
    lktr_c1 <- log(0.044); label("Cluster 1 transit rate constant ktr (1/min)") # Table 3 plasma cluster 1: ktransit = 0.044 1/min (RSE 33.1%)
    lntr_c1 <- log(8.3); label("Cluster 1 number of transit compartments (unitless)") # Table 3 plasma cluster 1: ntransit = 8.3 (RSE 19.2%)

    # Clusters 2 and 3 share Imax; clusters 2 and 4 share the potency; clusters
    # 3 and 7 share the potency (identical estimate and RSE in Table 3).
    limax_c23 <- log(0.20); label("Clusters 2 and 3 maximal inhibition Imax (fraction)") # Table 3 plasma clusters 2 and 3: Imax = -20% (RSE 30.1%)
    lec50_c24 <- log(113); label("Clusters 2 and 4 IC50 / EC50 of unbound plasma quinpirole (nM)") # Table 3 plasma cluster 2 IC50 = cluster 4 EC50 = 113 nM (RSE 98.5%)
    lec50_c37 <- log(17.2); label("Clusters 3 and 7 IC50 of unbound plasma quinpirole (nM)") # Table 3 plasma clusters 3 and 7: IC50 = 17.2 nM (RSE 50.6%)

    # Cluster 2: Imax turnover
    lkout_c2 <- log(0.057); label("Cluster 2 turnover rate constant kout (1/min)") # Table 3 plasma cluster 2: kout = 0.057 1/min (RSE 38.3%)

    # Cluster 3: Imax turnover
    lkout_c3 <- log(0.11); label("Cluster 3 turnover rate constant kout (1/min)") # Table 3 plasma cluster 3: kout = 0.11 1/min (RSE 12.2%)

    # Cluster 4: Emax turnover with transit-delayed drug effect
    lemax_c4 <- log(3.63); label("Cluster 4 maximal stimulation Emax (fraction of baseline)") # Table 3 plasma cluster 4: Emax = 363% (RSE 67.5%)
    lkout_c4 <- log(9.58); label("Cluster 4 turnover rate constant kout (1/min)") # Table 3 plasma cluster 4: kout = 9.58 1/min (RSE 104%)
    lktr_c4 <- log(0.0052); label("Cluster 4 transit rate constant ktr (1/min)") # Table 3 plasma cluster 4: ktransit = 0.0052 1/min (RSE 46.8%)
    lntr_c4 <- log(1.79); label("Cluster 4 number of transit compartments (unitless)") # Table 3 plasma cluster 4: ntransit = 1.79 (RSE 17.9%)

    # Cluster 5: Imax pool model
    limax_c5 <- log(0.41); label("Cluster 5 maximal inhibition Imax (fraction)") # Table 3 plasma cluster 5: Imax = -41% (RSE 14.6%)
    lec50_c5 <- log(339); label("Cluster 5 IC50 of unbound plasma quinpirole (nM)") # Table 3 plasma cluster 5: IC50 = 339 nM (RSE 32.8%)
    lkout_c5 <- log(0.11); label("Cluster 5 turnover rate constant kout (1/min)") # Table 3 plasma cluster 5: kout = 0.11 1/min (RSE 12.5%)
    lkrel_c5 <- log(0.018); label("Cluster 5 release rate constant from the pool krel (1/min)") # Table 3 plasma cluster 5: krel = 0.018 1/min (RSE 27.5%)

    # Cluster 6: Imax pool model
    limax_c6 <- log(0.90); label("Cluster 6 maximal inhibition Imax (fraction)") # Table 3 plasma cluster 6: Imax = -90% (RSE 0.3%)
    lec50_c6 <- fixed(log(0.001)); label("Cluster 6 IC50 of unbound plasma quinpirole (nM)") # Table 3 plasma cluster 6: IC50 = 0.001 nM (fix)
    lkout_c6 <- log(0.10); label("Cluster 6 turnover rate constant kout (1/min)") # Table 3 plasma cluster 6: kout = 0.10 1/min (RSE 18.4%)
    lkrel_c6 <- log(0.89); label("Cluster 6 release rate constant from the pool krel (1/min)") # Table 3 plasma cluster 6: krel = 0.89 1/min (RSE 19.7%)

    # Cluster 7: Imax turnover
    limax_c7 <- log(0.41); label("Cluster 7 maximal inhibition Imax (fraction)") # Table 3 plasma cluster 7: Imax = -41% (RSE 6.4%)
    lkout_c7 <- log(0.060); label("Cluster 7 turnover rate constant kout (1/min)") # Table 3 plasma cluster 7: kout = 0.060 1/min (RSE 13.5%)

    # Each output is normalised to the typical pre-dose baseline of its
    # biomarker, so the typical baseline is 1 by construction; the per-biomarker
    # baselines are listed in the validation vignette.
    lrbase <- fixed(log(1)); label("Typical baseline of every cluster output (relative to the biomarker's typical baseline, unitless)") # normalisation; Data S2 BSL0 = exp(THETA(1..28) + ETA(1))

    etalcl ~ 0.12 # vdB 2018 Table 1: IIV CL_PL,o = 0.12 (RSE 32%, shrinkage 7%)
    etalrbase ~ 0.917136 # Data S2 OMEGA = 9.17136E-01 (exponential IIV on the biomarker baseline BSL0)

    # Log-transformed data with an additive residual on the log scale; Data S2
    # estimates one SIGMA for all seven clusters, repeated here per output.
    expSd <- sqrt(0.08); label("Residual error, plasma quinpirole (log scale SD)") # vdB 2018 Table 1: RUV C_QP,PL = 0.08 (variance)
    expSd_Cbrain_striatum <- sqrt(0.12); label("Residual error, brain ECF quinpirole (log scale SD)") # vdB 2018 Table 1: RUV C_QP,ECF = 0.12 (variance)
    expSd_resp_c1 <- sqrt(0.0756288); label("Residual error, cluster 1 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c2 <- sqrt(0.0756288); label("Residual error, cluster 2 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c3 <- sqrt(0.0756288); label("Residual error, cluster 3 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c4 <- sqrt(0.0756288); label("Residual error, cluster 4 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c5 <- sqrt(0.0756288); label("Residual error, cluster 5 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c6 <- sqrt(0.0756288); label("Residual error, cluster 6 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
    expSd_resp_c7 <- sqrt(0.0756288); label("Residual error, cluster 7 (log scale SD)") # Data S2 SIGMA = 7.56288E-02 (variance, shared by all clusters)
  })
  model({
    # Quinpirole free base molecular weight (g/mol), used to express the ug/L
    # concentration in nM as Table 3 does: Data S1 IC50 exp(THETA(26)) =
    # 26.8 ug/L corresponds to the Table 3 value of 122 nM.
    mw <- 219.33

    # Individual PK parameters; clearances converted from L/h to L/min
    cl <- exp(lcl + etalcl) / 60
    qin <- exp(lqin) / 60
    qout <- exp(lqout) / 60
    vc <- exp(lvc)
    vbr <- exp(lvbr)

    # PD parameters
    rbase <- exp(lrbase + etalrbase)
    emax_c1 <- exp(lemax_c1)
    ec50_c1 <- exp(lec50_c1)
    kout_c1 <- exp(lkout_c1)
    ktr_c1 <- exp(lktr_c1)
    ntr_c1 <- exp(lntr_c1)
    imax_c23 <- exp(limax_c23)
    ec50_c24 <- exp(lec50_c24)
    ec50_c37 <- exp(lec50_c37)
    kout_c2 <- exp(lkout_c2)
    kout_c3 <- exp(lkout_c3)
    emax_c4 <- exp(lemax_c4)
    kout_c4 <- exp(lkout_c4)
    ktr_c4 <- exp(lktr_c4)
    ntr_c4 <- exp(lntr_c4)
    imax_c5 <- exp(limax_c5)
    ec50_c5 <- exp(lec50_c5)
    kout_c5 <- exp(lkout_c5)
    krel_c5 <- exp(lkrel_c5)
    imax_c6 <- exp(limax_c6)
    ec50_c6 <- exp(lec50_c6)
    kout_c6 <- exp(lkout_c6)
    krel_c6 <- exp(lkrel_c6)
    imax_c7 <- exp(limax_c7)
    kout_c7 <- exp(lkout_c7)

    # Quinpirole concentrations: Cc and Cbrain_striatum in ug/L, cnm in nM
    Cc <- central / vc
    Cbrain_striatum <- brain_striatum / vbr
    cnm <- Cc * 1000 / mw

    # Time after the quinpirole dose of the experiment day (Data S2 TIM);
    # zero before the first dose so the transit term stays defined.
    tdose <- tad(central)
    if (is.na(tdose)) tdose <- 0

    # Gamma-shaped transit delay of the drug effect (Data S2 TRANSIT), with
    # the Stirling form of the gamma function: it rises from 0 at the dose to
    # a peak of 1 / sqrt(2 * pi * ntr) at tdose = ntr / ktr.
    transit_c1 <- exp(-ktr_c1 * tdose) * (ktr_c1 * tdose)^ntr_c1 /
      (exp(-ntr_c1) * sqrt(2 * pi) * ntr_c1^(ntr_c1 + 0.5))
    transit_c4 <- exp(-ktr_c4 * tdose) * (ktr_c4 * tdose)^ntr_c4 /
      (exp(-ntr_c4) * sqrt(2 * pi) * ntr_c4^(ntr_c4 + 0.5))

    # Drug effects on biomarker production (turnover) or release (pool)
    eff_c1 <- 1 + emax_c1 * cnm / (ec50_c1 + cnm) * transit_c1
    eff_c2 <- 1 - imax_c23 * cnm / (ec50_c24 + cnm)
    eff_c3 <- 1 - imax_c23 * cnm / (ec50_c37 + cnm)
    eff_c4 <- 1 + emax_c4 * cnm / (ec50_c24 + cnm) * transit_c4
    eff_c5 <- 1 - imax_c5 * cnm / (ec50_c5 + cnm)
    eff_c6 <- 1 - imax_c6 * cnm / (ec50_c6 + cnm)
    eff_c7 <- 1 - imax_c7 * cnm / (ec50_c37 + cnm)

    # Quinpirole PK (van den Brink 2018 Supporting Information Eq. S3a-b)
    d/dt(central) <- -cl / vc * central - qin / vc * central + qout / vbr * brain_striatum
    d/dt(brain_striatum) <- qin / vc * central - qout / vbr * brain_striatum

    # Turnover clusters: kin = rbase * kout, so the pre-dose steady state is rbase
    d/dt(resp_c1) <- rbase * kout_c1 * eff_c1 - kout_c1 * resp_c1
    d/dt(resp_c2) <- rbase * kout_c2 * eff_c2 - kout_c2 * resp_c2
    d/dt(resp_c3) <- rbase * kout_c3 * eff_c3 - kout_c3 * resp_c3
    d/dt(resp_c4) <- rbase * kout_c4 * eff_c4 - kout_c4 * resp_c4
    d/dt(resp_c7) <- rbase * kout_c7 * eff_c7 - kout_c7 * resp_c7

    # Pool clusters: zero-order input rbase * kout into the pool, drug-inhibited
    # release krel * pool into the measured compartment
    d/dt(pool_c5) <- rbase * kout_c5 - krel_c5 * pool_c5 * eff_c5
    d/dt(resp_c5) <- krel_c5 * pool_c5 * eff_c5 - kout_c5 * resp_c5
    d/dt(pool_c6) <- rbase * kout_c6 - krel_c6 * pool_c6 * eff_c6
    d/dt(resp_c6) <- krel_c6 * pool_c6 * eff_c6 - kout_c6 * resp_c6

    resp_c1(0) <- rbase
    resp_c2(0) <- rbase
    resp_c3(0) <- rbase
    resp_c4(0) <- rbase
    resp_c5(0) <- rbase
    resp_c6(0) <- rbase
    resp_c7(0) <- rbase
    pool_c5(0) <- rbase * kout_c5 / krel_c5
    pool_c6(0) <- rbase * kout_c6 / krel_c6

    Cc ~ lnorm(expSd)
    Cbrain_striatum ~ lnorm(expSd_Cbrain_striatum)
    resp_c1 ~ lnorm(expSd_resp_c1)
    resp_c2 ~ lnorm(expSd_resp_c2)
    resp_c3 ~ lnorm(expSd_resp_c3)
    resp_c4 ~ lnorm(expSd_resp_c4)
    resp_c5 ~ lnorm(expSd_resp_c5)
    resp_c6 ~ lnorm(expSd_resp_c6)
    resp_c7 ~ lnorm(expSd_resp_c7)
  })
}
