vandenBrink_2019_quinpirole_rat_brain <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Cluster-based PK/PD model of the",
    "striatal brain extracellular fluid (microdialysate) metabolomic",
    "(biogenic amine and amino acid) response to the dopamine D2/3 agonist",
    "quinpirole. Quinpirole PK is the two-compartment plasma / striatal",
    "brain extracellular fluid model of van den Brink 2018 (linear",
    "elimination from plasma, asymmetric plasma-to-brain influx and",
    "brain-to-plasma efflux clearances). Nineteen responding brain",
    "biomarkers are grouped into six response clusters (2-7; cluster 1",
    "showed no response), each an inhibitory Imax turnover model driven by",
    "the unbound brain extracellular fluid quinpirole concentration, with a",
    "high-potency group (clusters 2-4) and a 122 nM group (clusters 5-7,",
    "which hold dopamine, DOPAC and HVA). Each cluster output is the",
    "biomarker level relative to its typical pre-dose baseline."
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
    "resp_c2",
    "resp_c3",
    "resp_c4",
    "resp_c5",
    "resp_c6",
    "resp_c7"
  )

  compartmentData <- list(
    central = list(analyte = "quinpirole", units = "ug", specimen = "plasma", verified = TRUE),
    brain_striatum = list(analyte = "quinpirole", units = "ug", specimen = "brain ISF", verified = TRUE),
    resp_c2 = list(
      analyte = "brain ECF biomarker cluster 2 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c3 = list(
      analyte = "brain ECF biomarker cluster 3 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c4 = list(
      analyte = "brain ECF biomarker cluster 4 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c5 = list(
      analyte = "brain ECF biomarker cluster 5 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c6 = list(
      analyte = "brain ECF biomarker cluster 6 (relative to baseline)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    resp_c7 = list(
      analyte = "brain ECF biomarker cluster 7 (relative to baseline)",
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
      "Striatal (caudate putamen) microdialysate collected in 20-min",
      "intervals from -200 to 180 min. The quinpirole PK model was fitted in",
      "the parent study to 0.17-2.14 mg/kg IV data (van den Brink 2018).",
      "Cluster membership of the brain ECF biomarkers (Figure 1, Figure 3):",
      "cluster 1 (no response, not modelled) N6,N6,N6-trimethyllysine;",
      "cluster 2 gamma-aminobutyric acid, 5-hydroxytryptophan,",
      "2-aminoadipic acid, cystathionine, homocitrulline, 4-hydroxyproline,",
      "lysine, tryptophan, tyrosine; cluster 3 3-methoxytyramine,",
      "hydroxylysine, proline; cluster 4 sarcosine, ornithine, glutamine;",
      "cluster 5 putrescine, DOPAC; cluster 6 dopamine; cluster 7 HVA."
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

    # Imax per cluster (Table 3, brainECF; Data S1 IMAX = exp(THETA) / 10)
    limax_c2 <- log(0.20); label("Cluster 2 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 2: Imax = -20% (RSE 6.1%)
    limax_c3 <- log(0.29); label("Cluster 3 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 3: Imax = -29% (RSE 7.1%)
    limax_c4 <- log(0.15); label("Cluster 4 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 4: Imax = -15% (RSE 13.5%)
    limax_c5 <- log(0.41); label("Cluster 5 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 5: Imax = -41% (RSE 9.0%)
    limax_c6 <- log(0.67); label("Cluster 6 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 6: Imax = -67% (RSE 4.9%)
    limax_c7 <- log(0.60); label("Cluster 7 maximal inhibition Imax (fraction)") # Table 3 brainECF cluster 7: Imax = -60% (RSE 9.3%)

    # Two shared potencies (Data S1 THETA(22) for clusters 2-4, THETA(26) for 5-7)
    lec50_c234 <- fixed(log(0.001)); label("Clusters 2-4 IC50 of unbound brain ECF quinpirole (nM)") # Table 3 brainECF clusters 2-4: IC50 = 0.001 nM (fix)
    lec50_c567 <- log(122); label("Clusters 5-7 IC50 of unbound brain ECF quinpirole (nM)") # Table 3 brainECF clusters 5-7: IC50 = 122 nM (RSE 51.4%)

    # Turnover rate constants; clusters 3 and 5 share kout (Data S1 THETA(30))
    lkout_c2 <- log(0.056); label("Cluster 2 turnover rate constant kout (1/min)") # Table 3 brainECF cluster 2: kout = 0.056 1/min (RSE 27.9%)
    lkout_c35 <- log(0.13); label("Clusters 3 and 5 turnover rate constant kout (1/min)") # Table 3 brainECF clusters 3 and 5: kout = 0.13 1/min (RSE 13.3%)
    lkout_c4 <- log(0.14); label("Cluster 4 turnover rate constant kout (1/min)") # Table 3 brainECF cluster 4: kout = 0.14 1/min (RSE 32.7%)
    lkout_c6 <- log(0.44); label("Cluster 6 turnover rate constant kout (1/min)") # Table 3 brainECF cluster 6: kout = 0.44 1/min (RSE 47.9%)
    lkout_c7 <- log(0.031); label("Cluster 7 turnover rate constant kout (1/min)") # Table 3 brainECF cluster 7: kout = 0.031 1/min (RSE 28.9%)

    # Each output is normalised to the typical pre-dose baseline of its
    # biomarker, so the typical baseline is 1 by construction; the per-biomarker
    # baselines are listed in the validation vignette.
    lrbase <- fixed(log(1)); label("Typical baseline of every cluster output (relative to the biomarker's typical baseline, unitless)") # normalisation; Data S1 BSL0 = exp(THETA(1..20) + ETA(1))

    etalcl ~ 0.12 # vdB 2018 Table 1: IIV CL_PL,o = 0.12 (RSE 32%, shrinkage 7%)
    etalrbase ~ 0.175155 # Data S1 OMEGA = 1.75155E-01 (exponential IIV on the biomarker baseline BSL0)

    # Log-transformed data with an additive residual on the log scale (one
    # SIGMA per cluster in Data S1, reported there as variances)
    expSd <- sqrt(0.08); label("Residual error, plasma quinpirole (log scale SD)") # vdB 2018 Table 1: RUV C_QP,PL = 0.08 (variance)
    expSd_Cbrain_striatum <- sqrt(0.12); label("Residual error, brain ECF quinpirole (log scale SD)") # vdB 2018 Table 1: RUV C_QP,ECF = 0.12 (variance)
    expSd_resp_c2 <- sqrt(0.0904392); label("Residual error, cluster 2 (log scale SD)") # Data S1 SIGMA(1) = 9.04392E-02 (variance)
    expSd_resp_c3 <- sqrt(0.109598); label("Residual error, cluster 3 (log scale SD)") # Data S1 SIGMA(2) = 1.09598E-01 (variance)
    expSd_resp_c4 <- sqrt(0.0718387); label("Residual error, cluster 4 (log scale SD)") # Data S1 SIGMA(3) = 7.18387E-02 (variance)
    expSd_resp_c5 <- sqrt(0.0538513); label("Residual error, cluster 5 (log scale SD)") # Data S1 SIGMA(4) = 5.38513E-02 (variance)
    expSd_resp_c6 <- sqrt(0.0394161); label("Residual error, cluster 6 (log scale SD)") # Data S1 SIGMA(5) = 3.94161E-02 (variance)
    expSd_resp_c7 <- sqrt(0.0281621); label("Residual error, cluster 7 (log scale SD)") # Data S1 SIGMA(6) = 2.81621E-02 (variance)
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
    imax_c2 <- exp(limax_c2)
    imax_c3 <- exp(limax_c3)
    imax_c4 <- exp(limax_c4)
    imax_c5 <- exp(limax_c5)
    imax_c6 <- exp(limax_c6)
    imax_c7 <- exp(limax_c7)
    ec50_c234 <- exp(lec50_c234)
    ec50_c567 <- exp(lec50_c567)
    kout_c2 <- exp(lkout_c2)
    kout_c35 <- exp(lkout_c35)
    kout_c4 <- exp(lkout_c4)
    kout_c6 <- exp(lkout_c6)
    kout_c7 <- exp(lkout_c7)

    # Quinpirole concentrations: Cc and Cbrain_striatum in ug/L, cbrnm in nM
    Cc <- central / vc
    Cbrain_striatum <- brain_striatum / vbr
    cbrnm <- Cbrain_striatum * 1000 / mw

    # Drug inhibition of biomarker production by brain ECF quinpirole
    eff_c2 <- 1 - imax_c2 * cbrnm / (ec50_c234 + cbrnm)
    eff_c3 <- 1 - imax_c3 * cbrnm / (ec50_c234 + cbrnm)
    eff_c4 <- 1 - imax_c4 * cbrnm / (ec50_c234 + cbrnm)
    eff_c5 <- 1 - imax_c5 * cbrnm / (ec50_c567 + cbrnm)
    eff_c6 <- 1 - imax_c6 * cbrnm / (ec50_c567 + cbrnm)
    eff_c7 <- 1 - imax_c7 * cbrnm / (ec50_c567 + cbrnm)

    # Quinpirole PK (van den Brink 2018 Supporting Information Eq. S3a-b)
    d/dt(central) <- -cl / vc * central - qin / vc * central + qout / vbr * brain_striatum
    d/dt(brain_striatum) <- qin / vc * central - qout / vbr * brain_striatum

    # Turnover clusters: kin = rbase * kout, so the pre-dose steady state is rbase
    d/dt(resp_c2) <- rbase * kout_c2 * eff_c2 - kout_c2 * resp_c2
    d/dt(resp_c3) <- rbase * kout_c35 * eff_c3 - kout_c35 * resp_c3
    d/dt(resp_c4) <- rbase * kout_c4 * eff_c4 - kout_c4 * resp_c4
    d/dt(resp_c5) <- rbase * kout_c35 * eff_c5 - kout_c35 * resp_c5
    d/dt(resp_c6) <- rbase * kout_c6 * eff_c6 - kout_c6 * resp_c6
    d/dt(resp_c7) <- rbase * kout_c7 * eff_c7 - kout_c7 * resp_c7

    resp_c2(0) <- rbase
    resp_c3(0) <- rbase
    resp_c4(0) <- rbase
    resp_c5(0) <- rbase
    resp_c6(0) <- rbase
    resp_c7(0) <- rbase

    Cc ~ lnorm(expSd)
    Cbrain_striatum ~ lnorm(expSd_Cbrain_striatum)
    resp_c2 ~ lnorm(expSd_resp_c2)
    resp_c3 ~ lnorm(expSd_resp_c3)
    resp_c4 ~ lnorm(expSd_resp_c4)
    resp_c5 ~ lnorm(expSd_resp_c5)
    resp_c6 ~ lnorm(expSd_resp_c6)
    resp_c7 ~ lnorm(expSd_resp_c7)
  })
}
