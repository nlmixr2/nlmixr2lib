Bouhaddou_2020_iadademstat_invitro <- function() {
  description <- "QSP. In vitro (NCI-H510A small-cell lung cancer cell line). Semimechanistic PD model of the covalent LSD1 (KDM1A) inhibitor iadademstat (ORY-1001) from Bouhaddou 2020. A constant drug concentration in the culture (assumed equal to the free intracellular concentration; Eq. 6, dROc/dt = 0) irreversibly inactivates LSD1 (kinact / Ki), the drug-bound LSD1 is degraded by Michaelis-Menten kinetics, percent target engagement (TE) suppresses GRP mRNA through a sigmoidal inhibition with a self-regulating production rate, and the GRP level switches cells between a proliferating and a quiescent population (reversible cytostasis). Fitted by least squares to target-engagement, GRP mRNA, drug-free growth and 10-day viability data under continuous and pulsed (wash-off) exposure. Exposure is set by dosing the drug concentration (nM) into roc and ended by a replacement (evid = 5, amt = 0) event, reproducing the wash-off experiments. The same PD parameters drive the in vivo model modellib('Bouhaddou_2020_iadademstat_mouse'), which differs only in the growth constant kP and its PK front end. Deterministic: no between-well variability and no residual error."
  reference <- paste(
    "Bouhaddou M, Yu LJ, Lunardi S, Stamatelos SK, Mack F, Gallo JM,",
    "Birtwistle MR, Walz AC. (2020). Predicting In Vivo Efficacy from In Vitro",
    "Data: Quantitative Systems Pharmacology Modeling for an Epigenetic Modifier",
    "Drug in Cancer. Clin Transl Sci 13(2):419-429. doi:10.1111/cts.12727.",
    "Parameter values from Supplementary Table S2 and the deposited MATLAB code",
    "(Supplementary zip, vitro_model/rateconstants_vitro_BEST.txt, RunModelVitro.m)."
  )
  vignette <- "Bouhaddou_2020_iadademstat"

  # Paper-mechanistic states with no canonical analogue in
  # inst/references/compartment-names.md: the culture / intracellular drug
  # concentration ROc (a state only so that exposure and wash-off can be
  # given as events), drug-bound (inactivated) LSD1, the GRP mRNA biomarker,
  # and the proliferating / quiescent cell pools of the epigenetic switch
  # (named as in Mazzocco_2015_temozolomide.R).
  paper_specific_compartments <- c("roc", "lsd1b", "grp", "prolif", "quiesc")

  dosing <- "roc"

  units <- list(
    time = "h",
    dosing = "nM",
    concentration = "nM for roc and the LSD1 states; number of cells for cell_count"
  )

  # Issue #482. verified = TRUE -- checked against Bouhaddou 2020 Eqs. 6-20,
  # Table S1 and the deposited RunModelVitro.m. An in vitro culture has no
  # biological matrix in the specimen vocabulary, hence 'not applicable'.
  compartmentData <- list(
    roc = list(
      analyte = "iadademstat (free drug in the culture, equal to the free intracellular concentration)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    lsd1b = list(
      analyte = "LSD1 (KDM1A) covalently bound to iadademstat",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    grp = list(
      analyte = "gastrin-releasing peptide (GRP) mRNA, relative to the untreated level",
      units = "unitless (1 = untreated)",
      specimen = "not applicable",
      verified = TRUE
    ),
    prolif = list(
      analyte = "proliferating NCI-H510A cells",
      units = "number of cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    quiesc = list(
      analyte = "quiescent NCI-H510A cells",
      units = "number of cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (NCI-H510A human small-cell lung cancer cell line)",
    n_subjects = "Cell-culture experiments: target engagement (3 concentrations x 4 time points), GRP mRNA (3 concentrations x 3 time points, plus a 5 nM pulsed / continuous arm), drug-free growth (6 time points) and 10-day viability (9 concentrations, plus 3 pulsed schedules); 3 biological replicates per condition",
    n_studies = 1L,
    disease_state = "Small-cell lung cancer (SCLC) cell line responsive to LSD1 inhibition",
    dose_range = "0.1-1000 nM iadademstat in the culture medium, continuous (up to 10 days, drug replenished at day 6) or pulsed (24 h, or 3, 5 or 7 days on, then washed off)",
    regions = "Roche Innovation Center New York (USA) and Oryzon Genomics (Spain)",
    notes = "Cells seeded at 8000 per well (96-well plates) for growth and viability; alamarBlue fluorescence converted to cell number with a linear calibration (y = 6.14x - 813.5). GRP mRNA by TaqMan RT-qPCR relative to vehicle at 24 h."
  )

  ini({
    # ==================================================================
    # Target engagement -- Eqs. 6-9, Table S2 'Value in vitro'.
    # ==================================================================
    ki <- fixed(19)
    label("Inhibitory constant Ki (nM), measured") # Table S2 Ki = 19 (Experiment)
    lkinact <- log(0.87052)
    label("Maximal inactivation rate constant kinact (1/h)") # Table S2 kinact = 0.87052
    lvmax_lsd1b <- log(0.03515)
    label("Maximal degradation rate of drug-bound LSD1 Vmax (nM/h)") # Table S2 vmax = 0.03515
    lkm_lsd1b <- log(0.9028)
    label("Michaelis-Menten constant of bound-LSD1 degradation Km (nM)") # Table S2 km = 0.9028
    llsd1tot <- log(3.0833)
    label("Total LSD1 in the cell LSD1TOTAL (nM)") # Table S2 LSD1TOTAL = 3.0833

    # ==================================================================
    # GRP mRNA biomarker -- Eqs. 10-12, Table S2.
    # ==================================================================
    lte50 <- log(26.513)
    label("Target engagement giving half-maximal GRP suppression K50 (%)") # Table S2 K50 = 26.513
    hill_te <- 2.0819
    label("Hill coefficient of TE on GRP production n (unitless)") # Table S2 n = 2.0819
    lkdeg_grp <- log(0.025567)
    label("First-order GRP mRNA degradation rate constant kdeg (1/h)") # Table S2 kdeg = 0.025567
    b_grp <- 0.083048
    label("Intercept b of the GRP-dependent production rate kmax = -m*GRP + b (1/h)") # Table S2 b = 0.083048

    # ==================================================================
    # Cell growth and epigenetic switch -- Eqs. 13-20, Table S2.
    # ==================================================================
    lkp <- log(0.023)
    label("Intrinsic growth rate constant of proliferating cells kP (1/h)") # Table S2 kP in vitro = 0.023
    lk50p <- log(100000)
    label("Proliferating-cell number giving half-maximal growth k50P (cells)") # Table S2 k50P = 100000
    lkmaxpq <- log(0.5932)
    label("Maximal proliferating-to-quiescent switching rate kmaxPQ (1/h)") # Table S2 kmaxPQ = 0.5932
    lkmaxqp <- log(4.4908)
    label("Maximal quiescent-to-proliferating switching rate kmaxQP (1/h)") # Table S2 kmaxQP = 4.4908
    k50pq <- 0.8836
    label("Biomarker level giving half-maximal P-to-Q switching k50PQ (unitless)") # Table S2 k50PQ = 0.8836
    k50qp <- 0.99897
    label("Biomarker level giving half-maximal Q-to-P switching k50QP (unitless)") # Table S2 k50QP = 0.99897
    hill_pq <- 3.8313
    label("Hill coefficient of the P-to-Q switch nPQ (unitless)") # Table S2 nPQ = 3.8313
    hill_qp <- 37.455
    label("Hill coefficient of the Q-to-P switch nQP (unitless)") # Table S2 nQP = 37.455
    bs <- 0.033328
    label("Basal bias term of Q-to-P switching bs (unitless)") # Table S2 bs = 0.033328

    prolif0 <- fixed(8000)
    label("Initial number of proliferating cells P(0) (cells)") # Table S1 P in vitro = 8000 cells (seeding density)
  })

  model({
    # Eq. 6: the drug is stable over the experiment, so its concentration is
    # constant between events. A dose of amt nM into roc starts an exposure;
    # a replacement event (evid = 5, amt = 0) is the wash-off.
    d/dt(roc) <- 0

    # Target engagement, Eqs. 7-9. abs() on the unbound LSD1 mirrors
    # RunModelVitro.m (LSD1U = abs(LSD1_0 - LSD1B)); it is inactive because
    # bound LSD1 never exceeds the total.
    kinact <- exp(lkinact)
    vmax_lsd1b <- exp(lvmax_lsd1b)
    km_lsd1b <- exp(lkm_lsd1b)
    lsd1tot <- exp(llsd1tot)
    lsd1u <- abs(lsd1tot - lsd1b)
    d/dt(lsd1b) <- kinact * roc / (ki + roc) * lsd1u -
      vmax_lsd1b / (km_lsd1b + lsd1b) * lsd1b
    TE <- lsd1b / lsd1tot * 100

    # GRP mRNA biomarker, Eqs. 10-12 (steady state GRP = 1 without drug).
    te50 <- exp(lte50)
    kdeg_grp <- exp(lkdeg_grp)
    m_grp <- b_grp - kdeg_grp
    kmax_grp <- -m_grp * grp + b_grp
    d/dt(grp) <- kmax_grp * te50^hill_te / (te50^hill_te + TE^hill_te) -
      kdeg_grp * grp
    grp(0) <- 1

    # Cell growth and epigenetic switch, Eqs. 13-20. BM = 1 - GRP (Eq. 18);
    # abs() as in RunModelVitro.m ('To prevent this from going negative').
    kp <- exp(lkp)
    k50p <- exp(lk50p)
    kmaxpq <- exp(lkmaxpq)
    kmaxqp <- exp(lkmaxqp)
    bm_grp <- abs(1 - grp)
    bmpq <- bm_grp
    bmqp <- 1 - bm_grp
    vp <- kp * (1 - prolif / (k50p + prolif)) * prolif
    vpq <- kmaxpq * bmpq^hill_pq / (k50pq^hill_pq + bmpq^hill_pq) * prolif
    vqp <- kmaxqp * (bmqp^hill_qp / (k50qp^hill_qp + bmqp^hill_qp) + bs) * quiesc
    d/dt(prolif) <- vp - vpq + vqp
    d/dt(quiesc) <- vpq - vqp
    prolif(0) <- prolif0

    # Total cell number, the Figure 3 viability / growth readout
    # (Plot_CG.m, Plot_CV.m: sum(yout(:,4:5), 2)).
    cell_count <- prolif + quiesc
  })
}
